# AkaghefAlgebra (Julia port)

MATLAB `AkaghefAlgebra` の Julia 移植。AlgLibMove プロジェクトの成果物。

Status: Phase C/D 進行中 (skeleton only)

## Usage

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

## CAS boundary (Sage / groups)

Plot/layout は pure Julia を維持する。群など成熟した汎用 CAS 演算は
再実装せず、薄い境界 `tools/sage`（永続 Sage セッション HTTP）経由で流用する。

```sh
julia --project=tools/sage -e 'using Pkg; Pkg.instantiate()'
julia --project=tools/sage tools/sage/scripts/verify_groups.jl
```

See [`tools/sage/README.md`](tools/sage/README.md).

## Docs

See [`dev-docs/projects/AlgLibMove/`](../dev-docs/projects/AlgLibMove/).
