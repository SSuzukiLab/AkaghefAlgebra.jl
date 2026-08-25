# SageBridge (Julia ↔ Sage)

Thin CAS boundary for **mature general-purpose operations** (especially group
theory via Sage/GAP). This is **not** a plot/layout dependency.

## Layout

```
tools/sage/
├─ Project.toml
├─ src/SageBridge.jl
├─ python/sage_session_server.py   # persistent sage.all namespace
├─ scripts/start_server_wsl.sh     # setsid launcher inside WSL
├─ scripts/verify_groups.jl
├─ artifacts/                      # verification records
└─ README.md
```

## Machine path (verified on this research PC)

```
Windows Julia
  → wsl file/curl transport
  → WSL localhost :8766
  → sage -python session (conda env `sage`, SageMath 10.7)
```

Notes:

- Direct Windows TCP to WSL `127.0.0.1` does **not** work on this host; the
  client defaults to `:wsl_file` (curl run inside WSL, stdout captured via temp files).
- The session server is started with `setsid` so it survives the launching WSL shell.
- The server launcher receives the repository's resolved server-script path from
  `SageBridge`; it does not assume a particular Windows username or checkout path.
- The server is deliberately bound to loopback only. Its `/exec` endpoint executes
  arbitrary Sage/Python code and has no authentication; do not expose it to a LAN
  or change the bind address.
- Fallback `:wsl_oneshot` runs `sage -python` per call if the daemon is down.

## Quick verify

```sh
# from lib/algebra.jl
julia --project=tools/sage -e 'using Pkg; Pkg.instantiate()'
julia --project=tools/sage tools/sage/scripts/verify_groups.jl
```

Latest machine record: `artifacts/verify_groups_result.json`.

## Minimal API

```julia
using SageBridge
sess = SageSession()                 # Windows default transport = :wsl_file
ensure_server!(sess)
exec!(sess, "G = SymmetricGroup(4)")
eval!(sess, "Integer(G.order())")    # → result 24
```

Group convenience wrappers: `symmetric_group!`, `group_order`,
`group_is_abelian`, `perm_product_repr`, `element_order`.

One-shot (no daemon):

```julia
SageBridge.oneshot_eval("Integer(SymmetricGroup(3).order())")
```
