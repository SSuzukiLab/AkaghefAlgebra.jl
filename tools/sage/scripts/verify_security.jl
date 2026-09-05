#!/usr/bin/env julia
# Verify the authenticated loopback boundary for the persistent Sage service.

using SageBridge
using Test

function request_status(sess::SageSession, method::String, path::String;
                        auth_token::Union{Nothing,String}=nothing,
                        content_type::Union{Nothing,String}=nothing,
                        origin::Union{Nothing,String}=nothing,
                        body="")
    url = "http://127.0.0.1:$(sess.port)$(path)"
    parts = ["curl -sS -o /dev/null -w '%{http_code}' -X '$method'"]
    auth_token === nothing || push!(parts, "-H 'Authorization: Bearer $auth_token'")
    content_type === nothing || push!(parts, "-H 'Content-Type: $content_type'")
    origin === nothing || push!(parts, "-H 'Origin: $origin'")
    isempty(body) || push!(parts, "--data '$body'")
    push!(parts, "'$url'")
    capture = SageBridge.wsl_bash_capture(join(parts, " "); timeout_s=15)
    status_text = strip(capture.stdout)
    occursin(r"^[0-9]{3}$", status_text) ||
        error("WSL curl returned no HTTP status: $(capture.stderr)")
    return parse(Int, status_text)
end

function main()
    sess = SageSession(timeout_s=60)
    ensure_server!(sess; force_restart=true)

    payload = "{\"expr\":\"1+1\"}"

    @testset "SageBridge authenticated boundary" begin
        @test sess.transport == :wsl_file
        @test health(sess)["ok"] == true
        @test eval!(sess, "1+1")["result"] == 2
        @test request_status(sess, "GET", "/health") == 401
        @test request_status(
            sess,
            "GET",
            "/health";
            auth_token=repeat("0", 64),
        ) == 401
        @test request_status(sess, "GET", "/health"; auth_token=sess.auth_token) == 200
        @test request_status(
            sess,
            "POST",
            "/exec";
            auth_token=sess.auth_token,
            content_type="application/json",
            origin="https://evil.invalid",
            body=payload,
        ) == 403
        @test request_status(
            sess,
            "POST",
            "/exec";
            auth_token=sess.auth_token,
            body=payload,
        ) == 415
        @test request_status(
            sess,
            "POST",
            "/exec";
            auth_token=sess.auth_token,
            content_type="application/json",
            body=payload,
        ) == 200
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
