---
id: algebra-symbolic-taft-validation
date: 2026-09-08
agent: Codex
---
# Symbolic Taft validation

The implementation follows [the specification](260908_symbolic_taft_spec.md).
The executable user entry point is
[`examples/260908_symbolic_taft.jl`](../examples/260908_symbolic_taft.jl).

## Reproduction

From the mathematics workspace root, the verified runtime is Julia 1.12.6:

```powershell
& 'C:\Users\Akaghef\.julia\juliaup\julia-1.12.6+0.x64.w64.mingw32\bin\julia.exe' --startup-file=no --pkgimages=no --project=lib/algebra.jl lib/algebra.jl/test/runtests.jl
& 'C:\Users\Akaghef\.julia\juliaup\julia-1.12.6+0.x64.w64.mingw32\bin\julia.exe' --startup-file=no --pkgimages=no --project=lib/algebra.jl lib/algebra.jl/examples/260908_symbolic_taft.jl
```

The existing project and dependency resolution were used; no dependencies were
added. Both new modules use Julia built-ins. This runs the full existing test
entry point through the ordinary `using AkaghefAlgebra` import, including the
existing numerical Hopf, polynomial, sparse and cochain tests.

## Verified scope

The concrete backend tests N = 2,3,4,5,6. It checks primitive-root arithmetic,
denominator-free q-binomial recurrence, sparse generator relations, nilpotency,
unit power, generator expansion of coproduct powers, coassociativity, both
counit identities and both antipode convolution identities on every PBW basis
element in that finite scope. The dual check independently builds powers of
the functional b by convolution, and checks the explicit Fourier dual-basis
expansion against all primal basis indices, including off-diagonal K indices.

Symbolic tests check alpha-equivalence, nested binders, capture-avoiding
substitution, deferred structure-constant traversal, range guards, equality
versus congruence, nonunit modular coefficients, incomplete Fourier domains,
residual quadratic sums, budget exits, BigInt intermediate arithmetic, and
negative integer powers with rational results. Fourier rules are compared to
direct sums of exact root powers for N = 2:6, a = -3:3, b = -1:1. Structure
constant specialization is compared with sparse coproduct, antipode and dual
operations for N = 2:4; the product-boundary test covers N = 2:6.

The executable example simplifies `Σ(s=0:N-1) ζ^(a s)` with symbolic N.
After simplification, it evaluates N = 1,000,003 at a = 0 and a = 1 with
`maxterms=0`, obtaining 1,000,003 and 0 respectively. It also checks an actual
product-coefficient contraction, integer/rational sparse arithmetic, and an
unresolved quadratic sum whose explicit N=4 specialization is `2+2ζ`.

## Result

The final full suite passed **23,331 / 23,331** tests in **3m08.0s**, exit code 0.
The final full-suite outcome and source hashes are recorded in the
[machine-readable run record](../test/results/260908_symbolic_taft_validation.json).
The executable example passed independently with exit code 0.

## Review and limits

The exact backend was implemented by a luna worker; the lead inspected the
formulas and constructor boundaries, tested the ordinary package interface,
and ran the full suite. A separate luna read-only review checked rewrite
conditions and identified the symbolic-N negative-power evaluation defect.
The lead reproduced it and added a regression after fixing rational coercion.
Integer scalar promotion, zero-term canonicalization and mixed-field rejection
were also corrected before the final run.

These are bounded computational checks and implementation review, not an
all-N mathematical proof or a formal rewrite certificate. The normalizer is
deliberately incomplete. General Gauss/q identities, a full MST/o-graph
frontend, generic symbolic complex parameters and a Symbolics.jl adapter are
not implemented. The N=1,000,003 example demonstrates removal of enumeration
for that Fourier family, not linear complexity for arbitrary contractions.
