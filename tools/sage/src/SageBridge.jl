"""
Thin Julia client for a persistent Sage session.

Design intent (owner boundary):
- plot/layout kernels stay pure Julia and must not require Sage
- mature general-purpose CAS ops (esp. group theory) are reused via Sage/GAP
  through this boundary, not reimplemented in Julia

Transports on this Windows research machine:
1. `:wsl_file` (default, reliable): Julia writes a request file, WSL `curl`
   talks to the Sage session on WSL localhost, response is written to a file.
2. `:wsl_oneshot`: no daemon; each call runs `sage -python` once (slower,
   always works if conda sage is installed).
3. `:direct`: HTTP.jl → host port (only if Windows→WSL TCP works).
"""
module SageBridge

using HTTP
using JSON
using Random

export SageSession, connect!, ensure_server!, health, exec!, eval!,
       sage_version, shutdown!, default_port, default_base_url,
       symmetric_group!, group_order, group_is_abelian, perm_product_repr, element_order,
       oneshot_eval

const DEFAULT_PORT = 8766
const DEFAULT_HOST = "127.0.0.1"

default_port() = DEFAULT_PORT
default_base_url(; host=DEFAULT_HOST, port=DEFAULT_PORT) = "http://$host:$port"

function server_script_path()
    return normpath(joinpath(@__DIR__, "..", "python", "sage_session_server.py"))
end

function windows_to_wsl(path::AbstractString)
    p = replace(abspath(path), "\\" => "/")
    m = match(r"^([A-Za-z]):/(.*)$", p)
    m === nothing && return p
    drive = lowercase(m.captures[1])
    rest = m.captures[2]
    return "/mnt/$drive/$rest"
end

"""Run a bash -lc command inside WSL, capturing stdout via a temp file.

Avoids Julia `read(\`wsl ...\`)` which has crashed with ACCESS_VIOLATION on
this host when curl writes to stderr.
"""
function wsl_bash_capture(inner::AbstractString; timeout_s::Real=120)
    out_path = joinpath(tempdir(), "sagebridge_out_$(getpid())_$(time_ns()).txt")
    err_path = joinpath(tempdir(), "sagebridge_err_$(getpid())_$(time_ns()).txt")
    wsl_out = windows_to_wsl(out_path)
    wsl_err = windows_to_wsl(err_path)
    wrapped = string(
        "set +e; ",
        "(", inner, ") >'", wsl_out, "' 2>'", wsl_err, "'; ",
        "ec=\$?; echo \$ec >'", wsl_out, ".ec'; exit 0",
    )
    try
        run(`wsl.exe bash -lc $wrapped`)
        body = isfile(out_path) ? read(out_path, String) : ""
        err = isfile(err_path) ? read(err_path, String) : ""
        ec_path = out_path * ".ec"
        ec = 1
        if isfile(ec_path)
            s = strip(read(ec_path, String))
            ec = try parse(Int, s) catch; 1 end
        end
        return (ec=ec, stdout=replace(body, "\r\n" => "\n"), stderr=replace(err, "\r\n" => "\n"))
    finally
        for p in (out_path, err_path, out_path * ".ec")
            isfile(p) && rm(p; force=true)
        end
    end
end

mutable struct SageSession
    base_url::String
    host::String
    port::Int
    timeout_s::Float64
    transport::Symbol  # :wsl_file | :wsl_oneshot | :direct
    auth_token::String
end

function SageSession(; host=DEFAULT_HOST, port=DEFAULT_PORT, timeout_s=180.0, transport::Symbol=:auto,
                     auth_token::Union{Nothing,AbstractString}=nothing)
    1 <= port <= 65535 || throw(ArgumentError("port must be in 1:65535"))
    host in ("127.0.0.1", "localhost") ||
        throw(ArgumentError("SageSession host must be loopback"))
    tr = transport === :auto ? (Sys.iswindows() ? :wsl_file : :direct) : transport
    token = auth_token === nothing ? bytes2hex(rand(RandomDevice(), UInt8, 32)) : String(auth_token)
    occursin(r"^[0-9a-f]{64}$", token) ||
        throw(ArgumentError("SageSession auth_token must be 64 lowercase hexadecimal characters"))
    return SageSession(default_base_url(; host, port), String(host), Int(port), Float64(timeout_s), tr, token)
end

struct BridgeResponse
    status::Int
    body::String
end

function _request_direct(sess::SageSession, method::Symbol, path::AbstractString; body::Union{Nothing,String}=nothing)
    url = sess.base_url * path
    headers = ["Authorization" => "Bearer $(sess.auth_token)"]
    resp = if method === :GET
        HTTP.get(url, headers; readtimeout=sess.timeout_s, connect_timeout=3, status_exception=false)
    else
        push!(headers, "Content-Type" => "application/json")
        HTTP.post(url, headers, body === nothing ? "{}" : body;
                  readtimeout=sess.timeout_s, connect_timeout=3, status_exception=false)
    end
    return BridgeResponse(Int(resp.status), String(resp.body))
end

function _request_wsl_file(sess::SageSession, method::Symbol, path::AbstractString; body::Union{Nothing,String}=nothing)
    url = "http://127.0.0.1:$(sess.port)$(path)"
    timeout = string(Int(ceil(sess.timeout_s)))
    auth_header = "-H 'Authorization: Bearer $(sess.auth_token)' "
    if method === :GET
        inner = "curl -sS -m $timeout -w '\\n%{http_code}' $auth_header'$url'"
        cap = wsl_bash_capture(inner; timeout_s=sess.timeout_s)
    else
        payload = body === nothing ? "{}" : body
        req_path = joinpath(tempdir(), "sagebridge_req_$(getpid())_$(time_ns()).json")
        write(req_path, payload)
        wsl_req = windows_to_wsl(req_path)
        try
            inner = string(
                "curl -sS -m ", timeout,
                " -w '\\n%{http_code}' ",
                auth_header,
                "-H 'Content-Type: application/json' ",
                "-d @'", wsl_req, "' ",
                "'", url, "'",
            )
            cap = wsl_bash_capture(inner; timeout_s=sess.timeout_s)
        finally
            isfile(req_path) && rm(req_path; force=true)
        end
    end
    out = cap.stdout
    lines = split(out, '\n')
    isempty(lines) && return BridgeResponse(0, isempty(cap.stderr) ? "" : cap.stderr)
    status_str = strip(lines[end])
    body_txt = join(lines[1:max(0, end - 1)], '\n')
    status = try parse(Int, status_str) catch; 0 end
    if status == 0 && !isempty(cap.stderr)
        body_txt = body_txt * "\n" * cap.stderr
    end
    return BridgeResponse(status, body_txt)
end

"""One-shot evaluation: no daemon. Good fallback and for batch verification."""
function oneshot_eval(expr::AbstractString; code::AbstractString="", timeout_s::Real=180)
    py = """
import json
from sage.all import *
$(code)
_result = $(expr)
def _j(v):
    if v is None or isinstance(v, (bool, int, float, str)):
        return v
    try:
        from sage.rings.integer import Integer as Z
        if isinstance(v, Z):
            return int(v)
    except Exception:
        pass
    if isinstance(v, (list, tuple)):
        return [_j(x) for x in v]
    if isinstance(v, dict):
        return {str(k): _j(x) for k, x in v.items()}
    return repr(v)
print(json.dumps({"ok": True, "result": _j(_result), "result_repr": repr(_result)}, ensure_ascii=False))
"""
    py_path = joinpath(tempdir(), "sagebridge_oneshot_$(getpid())_$(time_ns()).py")
    write(py_path, py)
    wsl_py = windows_to_wsl(py_path)
    try
        inner = string(
            "source \"\$HOME/dev/miniforge3/etc/profile.d/conda.sh\" && ",
            "conda activate sage && ",
            "sage -python '", wsl_py, "'",
        )
        cap = wsl_bash_capture(inner; timeout_s=timeout_s)
        # Sage/conda may return a non-zero status even when stdout is valid JSON.
        payload = try
            JSON.parse(strip(cap.stdout))
        catch
            nothing
        end
        if payload isa AbstractDict && get(payload, "ok", false) == true
            return payload
        end
        cap.ec == 0 || error("oneshot sage failed (ec=$(cap.ec)):\n$(cap.stdout)\n$(cap.stderr)")
        error("oneshot sage returned unparsable output:\n$(cap.stdout)\n$(cap.stderr)")
    finally
        isfile(py_path) && rm(py_path; force=true)
    end
end

function _request(sess::SageSession, method::Symbol, path::AbstractString; body::Union{Nothing,String}=nothing)
    if sess.transport === :direct
        return _request_direct(sess, method, path; body)
    elseif sess.transport === :wsl_file
        return _request_wsl_file(sess, method, path; body)
    elseif sess.transport === :wsl_oneshot
        error("HTTP path not available for :wsl_oneshot; use oneshot_eval / eval!")
    else
        return _request_wsl_file(sess, method, path; body)
    end
end

function _parse_json(body::AbstractString)
    isempty(strip(body)) && return Dict{String,Any}()
    return JSON.parse(body)
end

function health(sess::SageSession)
    if sess.transport === :wsl_oneshot
        try
            r = oneshot_eval("1"; timeout_s=min(sess.timeout_s, 60))
            return Dict("ok" => get(r, "ok", false) == true, "service" => "sage-oneshot", "version" => "oneshot")
        catch
            return nothing
        end
    end
    try
        resp = _request(sess, :GET, "/health")
        resp.status == 200 || return nothing
        return _parse_json(resp.body)
    catch
        return nothing
    end
end

function is_alive(sess::SageSession)
    h = health(sess)
    return h !== nothing && get(h, "ok", false) == true
end

function ensure_server!(sess::SageSession; force_restart::Bool=false)
    if sess.transport === :wsl_oneshot
        return sess
    end
    if !force_restart && is_alive(sess)
        return sess
    end
    if force_restart
        try
            shutdown!(sess)
        catch
        end
        wsl_bash_capture("pkill -f '[s]age_session_server.py.*--port $(sess.port)' >/dev/null 2>&1 || true")
        sleep(0.5)
    end

    launcher = normpath(joinpath(@__DIR__, "..", "scripts", "start_server_wsl.sh"))
    isfile(launcher) || error("missing launcher: $launcher")
    wsl_launcher = windows_to_wsl(launcher)
    wsl_server = windows_to_wsl(server_script_path())
    println("starting WSL sage session server on port $(sess.port) via transport=$(sess.transport) ...")
    inner = string(
        "mkdir -p \"\$HOME/bin\"; ",
        "cp '", wsl_launcher, "' \"\$HOME/bin/start_sage_session.sh\"; ",
        "chmod +x \"\$HOME/bin/start_sage_session.sh\"; ",
        "\"\$HOME/bin/start_sage_session.sh\" ", string(sess.port), " '", wsl_server,
        "' ", sess.auth_token,
    )
    cap = wsl_bash_capture(inner; timeout_s=180)
    println(strip(cap.stdout))
    if !isempty(strip(cap.stderr))
        println(strip(cap.stderr))
    end

    if is_alive(sess)
        return sess
    end
    deadline = time() + 45.0
    while time() < deadline
        is_alive(sess) && return sess
        sleep(1.0)
    end
    # Fall back to oneshot so the client remains usable.
    @warn "persistent Sage server unhealthy; falling back to :wsl_oneshot" port=sess.port
    sess.transport = :wsl_oneshot
    return sess
end

function connect!(sess::SageSession=SageSession())
    is_alive(sess) || error("no healthy Sage session; call ensure_server! first")
    return sess
end

function sage_version(sess::SageSession)
    if sess.transport === :wsl_oneshot
        r = oneshot_eval("str(sage.version.version)")
        return string(r["result"])
    end
    h = health(sess)
    h === nothing && error("server unreachable")
    return string(get(h, "version", "unknown"))
end

function exec!(sess::SageSession, code::AbstractString)
    if sess.transport === :wsl_oneshot
        # Oneshot has no durable namespace; execute for side effects is a no-op unless
        # the caller chains code+expr in eval!. Still run to surface syntax errors.
        oneshot_eval("None"; code=String(code))
        return Dict("ok" => true)
    end
    data = JSON.json(Dict("code" => String(code)))
    resp = _request(sess, :POST, "/exec"; body=data)
    payload = _parse_json(resp.body)
    if resp.status != 200 || get(payload, "ok", false) != true
        err = get(payload, "error", resp.body)
        tb = get(payload, "traceback", "")
        error("Sage exec failed: $err\n$tb")
    end
    return payload
end

function eval!(sess::SageSession, expr::AbstractString; bind::Union{Nothing,AbstractString}=nothing, code::AbstractString="")
    if sess.transport === :wsl_oneshot
        return oneshot_eval(expr; code=code, timeout_s=sess.timeout_s)
    end
    body = Dict{String,Any}("expr" => String(expr))
    isempty(code) || (body["code"] = String(code))
    bind === nothing || (body["bind"] = String(bind))
    data = JSON.json(body)
    resp = _request(sess, :POST, "/exec"; body=data)
    payload = _parse_json(resp.body)
    if resp.status != 200 || get(payload, "ok", false) != true
        err = get(payload, "error", resp.body)
        tb = get(payload, "traceback", "")
        error("Sage eval failed: $err\n$tb")
    end
    return payload
end

function shutdown!(sess::SageSession)
    sess.transport === :wsl_oneshot && return nothing
    try
        resp = _request(sess, :POST, "/shutdown"; body="{}")
        return _parse_json(resp.body)
    catch
        return nothing
    end
end

# ---------------------------------------------------------------------------
# Group-theory convenience
# ---------------------------------------------------------------------------

function symmetric_group!(sess::SageSession, n::Integer; name::AbstractString="G")
    if sess.transport === :wsl_oneshot
        # Stash definition in a session-local dict on the Julia side is impossible
        # for oneshot; callers should pass full expressions. Keep API for session mode.
        return name
    end
    exec!(sess, "$name = SymmetricGroup(Integer($n))")
    return name
end

function group_order(sess::SageSession, group_expr::AbstractString="G")
    payload = eval!(sess, "Integer(($group_expr).order())")
    return Int(payload["result"])
end

function group_is_abelian(sess::SageSession, group_expr::AbstractString="G")
    payload = eval!(sess, "bool(($group_expr).is_abelian())")
    return Bool(payload["result"])
end

function perm_product_repr(sess::SageSession, a_cycles::AbstractString, b_cycles::AbstractString; group_expr::AbstractString="G")
    code = """
_a = $(group_expr)("$a_cycles")
_b = $(group_expr)("$b_cycles")
_prod = _a * _b
"""
    if sess.transport === :wsl_oneshot
        code = """
$(group_expr) = SymmetricGroup(max(4, 5))
""" * code
    end
    payload = eval!(sess, "str(_prod)"; code=code)
    return string(payload["result"])
end

function element_order(sess::SageSession, element_expr::AbstractString)
    payload = eval!(sess, "Integer(($element_expr).order())")
    return Int(payload["result"])
end

end # module
