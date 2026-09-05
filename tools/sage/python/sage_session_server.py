#!/usr/bin/env python3
"""Persistent Sage session HTTP bridge for Julia (and other clients).

Run under Sage's Python so ``from sage.all import *`` works:

    sage -python sage_session_server.py --host 127.0.0.1 --port 8766

Unlike a ``sage -c``-per-request bridge, this keeps one in-process
``sage.all`` namespace across calls, so group elements and named objects
survive between Julia round-trips.
"""

from __future__ import annotations

import argparse
import contextlib
import hmac
import io
import json
import threading
import traceback
from http.server import BaseHTTPRequestHandler, HTTPServer
from typing import Any

# Requires sage -python
from sage.all import *  # noqa: F401,F403

EXEC_GLOBALS: dict[str, Any] = globals()


def _to_jsonable(value: Any) -> Any:
    """Best-effort conversion of Sage/Python values to JSON."""
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, (list, tuple)):
        return [_to_jsonable(v) for v in value]
    if isinstance(value, dict):
        return {str(k): _to_jsonable(v) for k, v in value.items()}
    # Integers larger than JSON safe range still serialize via str fallback below
    try:
        import sage.rings.integer  # type: ignore

        if isinstance(value, sage.rings.integer.Integer):
            return int(value)
    except Exception:
        pass
    try:
        import sage.rings.rational  # type: ignore

        if isinstance(value, sage.rings.rational.Rational):
            return str(value)
    except Exception:
        pass
    return repr(value)


class SageSessionHandler(BaseHTTPRequestHandler):
    server_version = "SageSession/0.1"

    def _send_json(self, status_code: int, payload: dict[str, Any]) -> None:
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        self.send_response(status_code)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, format: str, *args) -> None:  # noqa: A003
        return

    def _reject_untrusted(self) -> bool:
        """Block browser-originated (drive-by / DNS-rebinding) requests.

        This endpoint executes arbitrary Sage/Python. Even bound to loopback it
        must refuse cross-origin browser calls: a malicious page can POST here.
        Trusted local clients (Julia SageBridge) send no ``Origin``/``Sec-Fetch``
        headers and a loopback ``Host``.
        """
        allowed_hosts = {
            f"127.0.0.1:{self.server.server_address[1]}",
            f"localhost:{self.server.server_address[1]}",
            f"[::1]:{self.server.server_address[1]}",
        }
        host = self.headers.get("Host", "")
        if host and host not in allowed_hosts:
            self._send_json(403, {"ok": False, "error": "forbidden host"})
            return True
        expected = f"Bearer {self.server.auth_token}"
        supplied = self.headers.get("Authorization", "")
        if not hmac.compare_digest(supplied, expected):
            self._send_json(401, {"ok": False, "error": "unauthorized"})
            return True
        if self.headers.get("Origin"):
            self._send_json(403, {"ok": False, "error": "cross-origin requests are forbidden"})
            return True
        if self.headers.get("Sec-Fetch-Site") not in (None, "none", "same-origin"):
            self._send_json(403, {"ok": False, "error": "cross-site requests are forbidden"})
            return True
        return False

    def do_GET(self) -> None:  # noqa: N802
        if self._reject_untrusted():
            return
        if self.path == "/health":
            self._send_json(
                200,
                {
                    "ok": True,
                    "service": "sage-session",
                    "version": str(sage.version.version),
                },
            )
            return
        self._send_json(404, {"ok": False, "error": "not found"})

    def do_POST(self) -> None:  # noqa: N802
        if self._reject_untrusted():
            return
        if self.path == "/exec":
            self._handle_exec()
            return
        if self.path == "/reset":
            for key in list(EXEC_GLOBALS):
                if key not in BASE_GLOBAL_KEYS:
                    del EXEC_GLOBALS[key]
            self._send_json(200, {"ok": True, "message": "user namespace cleared"})
            return
        if self.path == "/shutdown":
            self._send_json(200, {"ok": True, "message": "shutting down"})
            threading.Thread(target=self.server.shutdown, daemon=True).start()
            return
        self._send_json(404, {"ok": False, "error": "not found"})

    def _handle_exec(self) -> None:
        content_type = self.headers.get("Content-Type", "")
        if "application/json" not in content_type.lower():
            self._send_json(415, {"ok": False, "error": "Content-Type must be application/json"})
            return
        try:
            content_len = int(self.headers.get("Content-Length", "0"))
        except ValueError:
            content_len = 0

        raw = self.rfile.read(content_len)
        try:
            payload = json.loads(raw.decode("utf-8")) if raw else {}
        except json.JSONDecodeError:
            self._send_json(400, {"ok": False, "error": "invalid JSON"})
            return

        code = payload.get("code") or ""
        expr = payload.get("expr")
        bind = payload.get("bind")  # optional name to store expr result under

        stdout_buf = io.StringIO()
        stderr_buf = io.StringIO()

        try:
            with contextlib.redirect_stdout(stdout_buf), contextlib.redirect_stderr(
                stderr_buf
            ):
                result = None
                if code:
                    exec(code, EXEC_GLOBALS, EXEC_GLOBALS)
                if expr is not None and str(expr).strip() != "":
                    result = eval(str(expr), EXEC_GLOBALS, EXEC_GLOBALS)
                    if bind:
                        EXEC_GLOBALS[str(bind)] = result
                elif bind and code:
                    # allow bind of last assigned name only if client sets it via expr
                    pass

            self._send_json(
                200,
                {
                    "ok": True,
                    "stdout": stdout_buf.getvalue(),
                    "stderr": stderr_buf.getvalue(),
                    "result": _to_jsonable(result),
                    "result_repr": repr(result) if result is not None else None,
                    "result_type": type(result).__name__ if result is not None else None,
                },
            )
        except Exception as ex:
            self._send_json(
                500,
                {
                    "ok": False,
                    "error": str(ex),
                    "stdout": stdout_buf.getvalue(),
                    "stderr": stderr_buf.getvalue(),
                    "traceback": traceback.format_exc(),
                },
            )


BASE_GLOBAL_KEYS = frozenset(EXEC_GLOBALS) | {"BASE_GLOBAL_KEYS"}


def main() -> None:
    parser = argparse.ArgumentParser(description="Persistent Sage session HTTP bridge")
    parser.add_argument(
        "--host",
        default="127.0.0.1",
        choices=("127.0.0.1", "localhost"),
        help="loopback only: this token-authenticated service executes Sage code",
    )
    parser.add_argument("--port", default=8766, type=int)
    parser.add_argument("--token", required=True, help=argparse.SUPPRESS)
    args = parser.parse_args()
    if len(args.token) != 64 or any(character not in "0123456789abcdef" for character in args.token):
        parser.error("--token must be 64 lowercase hexadecimal characters")

    # Deliberately single-threaded: exec/eval share one namespace and redirect
    # process-global stdout/stderr, so concurrent handlers would race.
    server = HTTPServer((args.host, args.port), SageSessionHandler)
    server.auth_token = args.token
    print(
        f"Sage session server listening on http://{args.host}:{args.port} "
        f"(sage {sage.version.version})",
        flush=True,
    )
    server.serve_forever()


if __name__ == "__main__":
    main()
