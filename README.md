# AkaghefAlgebra (Julia port)

MATLAB `AkaghefAlgebra` の Julia 移植。AlgLibMove プロジェクトの成果物。

Status: Phase C/D 進行中 (skeleton only)

## Usage

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

## CAS boundary (Sage / groups)

Exact cochain constructors widen integer/rational coefficients (including complex
ones) to `BigInt` rationals. `mode=:exact` uses `==`; with floating input this
means machine equality, not a mathematical exact-arithmetic certificate.

Plot/layout は pure Julia を維持する。群など成熟した汎用 CAS 演算は
再実装せず、薄い境界 `tools/sage`（永続 Sage セッション HTTP）経由で流用する。

```sh
julia --project=tools/sage -e 'using Pkg; Pkg.instantiate()'
julia --project=tools/sage tools/sage/scripts/verify_groups.jl
```

See [`tools/sage/README.md`](tools/sage/README.md).

## Finite-group cochains

`Cochain1`, `Cochain2`, and `Cochain3` store multiplicative cochains over a
supplied finite multiplication table. `is_normalized` and `is_cocycle` use
exact equality by default (`mode=:approx` enables tolerance-based checks), and
`coboundary2`/`coboundary3` implement the multiplicative differentials. The
coefficient type is generic over `Number`; `Rational` and
`Complex{<:Rational}` therefore remain available for exact finite checks.

## Docs

See [`dev-docs/projects/AlgLibMove/`](../dev-docs/projects/AlgLibMove/).
