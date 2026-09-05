#!/usr/bin/env bash
# Start persistent Sage session server for Julia SageBridge.
# Uses setsid so the server survives the end of the launching WSL session.
set -euo pipefail

PORT="${1:-8766}"
SCRIPT="${2:?usage: start_server_wsl.sh PORT SERVER_SCRIPT TOKEN}"
TOKEN="${3:?usage: start_server_wsl.sh PORT SERVER_SCRIPT TOKEN}"
[[ "$PORT" =~ ^[0-9]+$ ]] && (( PORT >= 1 && PORT <= 65535 )) || {
  echo "invalid port: $PORT" >&2
  exit 2
}
[[ "$TOKEN" =~ ^[0-9a-f]{64}$ ]] || {
  echo "invalid authentication token" >&2
  exit 2
}
HOST="127.0.0.1"
LOG="/tmp/sage_session_${PORT}.log"

# shellcheck disable=SC1090
source "$HOME/dev/miniforge3/etc/profile.d/conda.sh"
conda activate sage

if [[ ! -f "$SCRIPT" ]]; then
  echo "missing server script: $SCRIPT" >&2
  exit 1
fi

# If already healthy, do nothing.
if curl -fsS -H "Authorization: Bearer ${TOKEN}" "http://127.0.0.1:${PORT}/health" >/dev/null 2>&1; then
  echo "already healthy on port ${PORT}"
  curl -fsS -H "Authorization: Bearer ${TOKEN}" "http://127.0.0.1:${PORT}/health"
  echo
  exit 0
fi

# Best-effort cleanup of a dead/half-open listener.
pkill -f "[s]age_session_server.py.*--port ${PORT}" >/dev/null 2>&1 || true
sleep 0.3

# Detach from this WSL invocation's process group so the server keeps running.
setsid nohup bash -c \
  "source \"$HOME/dev/miniforge3/etc/profile.d/conda.sh\" && conda activate sage && exec sage -python \"$SCRIPT\" --host \"$HOST\" --port \"$PORT\" --token \"$TOKEN\"" \
  >"$LOG" 2>&1 < /dev/null &
echo "PID=$!"
echo "LOG=$LOG"

for i in $(seq 1 90); do
  if curl -fsS -H "Authorization: Bearer ${TOKEN}" "http://127.0.0.1:${PORT}/health" >/dev/null 2>&1; then
    curl -fsS -H "Authorization: Bearer ${TOKEN}" "http://127.0.0.1:${PORT}/health"
    echo
    exit 0
  fi
  sleep 0.5
done

echo "server failed to become healthy" >&2
tail -n 80 "$LOG" >&2 || true
exit 1
