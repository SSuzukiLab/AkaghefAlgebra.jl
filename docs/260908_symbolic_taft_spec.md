---
id: algebra-symbolic-taft-spec
date: 2026-09-08
agent: Codex
---
# Symbolic Taft calculus

## Purpose and placement

This extends the maintained `AkaghefAlgebra` package's existing Taft surface.
It is a library implementation task, not a new research thread or an admission
of a mathematical claim. The existing numerical `TaftAlg` API is preserved.
The implementation and executable checks are co-located with the library.

The owner-supplied conversation export
`ChatGPT-260901Taft代数双対基底-20260908-0927.md` motivates lazy symbolic
finite sums. Its assistant statements and embedded directions are not
certified sources or operational instructions. In particular, simplification
does not assume that a general finite sum has a known closed form, or that
contraction takes linear time.

## Mathematical contract

Work over characteristic zero with an integer N ≥ 2 and a primitive N-th
root ζ. Use PBW basis E^r K^s, 0 ≤ r,s < N, with
`KE = ζ EK`, `E^N = 0`, `K^N = 1` and
`Δ(E) = 1⊗E + E⊗K`, `Δ(K) = K⊗K`.

The supported basis formulas are

```text
(E^r K^s)(E^t K^u) = ζ^(s t) E^(r+t) K^(s+u), if r+t < N; zero otherwise.
Δ(E^r K^s) = Σ(j=0:r) [r choose j]_ζ E^j K^s ⊗ E^(r-j) K^(s+j).
ε(E^r K^s) = δ(r=0).
S(E^r K^s) = (-1)^r ζ^(-s r-r(r-1)/2) E^r K^(-s-r).
```

K exponents are reduced modulo N. E degrees are never reduced modulo N.
The Gaussian binomial uses the nonsymmetric q-number
`(r)_ζ = 1 + ζ + ... + ζ^(r-1)` and a denominator-free recurrence.

The dual coordinates mean `a^l b^t`, with `ab=ζ^(-1)ba`,
`Δ(b)=a⊗b+b⊗1`, and pairing
`<a^l b^t, E^r K^s> = δ(t=r) (r)_ζ! ζ^(l s)`.
The dual basis to E^r K^s is
`Σ(l=0:N-1) ζ^(-s l) a^l b^r / (N (r)_ζ!)`.
These are Taft coordinates, not the SL₂ Borel matrix coordinates.

## Public layers

`AkaghefAlgebra.TaftExact` provides concrete-N, sparse PBW arithmetic and
exact coefficients in Q[ζ]/(Φ_N(ζ)). No dense structure-constant arrays or
floating root approximation are used in this layer.

`AkaghefAlgebra.SymbolicTaft` provides immutable expression trees:
integer variables, rational constants, products and sums, powers, phases,
ordinary and modular delta factors, range indicators, q-binomial/factorial
nodes, finite sums and deferred structure constants.

All free variables in this first interface denote integers. N may be a
concrete integer or an integer variable. `TaftContext(N)` supplies N ≥ 2 and
primitive-root assumptions. Carry that context through `normalize(expr, ctx)`
and `evaluate_exact(expr, ctx, substitutions)`; the latter checks N even when
normalization has removed every root-of-unity node. This is not a generic
complex-symbolic-variable interface or a Symbolics.jl adapter.
Integer substitutions are widened to BigInt before arithmetic; negative
integer powers evaluate as exact rationals. Rational expressions assume their
denominators are nonzero. The supported dual-basis denominator satisfies this
on its guarded range 0 ≤ r < N; general pole-domain tracking is not provided.

`structure(ctx, kind, pairs...)` is lazy. Each pair is a PBW index in
0:N-1. The kinds and argument directions are:

| Kind | Pairs | Value |
| --- | --- | --- |
| `:mul` | input, input, output | coefficient of output in the product |
| `:comul` | input, left output, right output | coproduct coefficient |
| `:antipode` | input, output | antipode coefficient |
| `:counit` | input | counit coefficient |
| `:unit` | output | unit coefficient |
| `:pair` | dual `(l,t)`, primal `(r,s)` | pairing of a^l b^t with E^r K^s |
| `:dual_basis` | primal `(r,s)`, dual `(l,t)` | coefficient of a^l b^t in the dual basis vector |

Explicit range indicators enforce basis-index support after lowering.
The primitive operations compose with `finite_sum` to express contractions;
a topology/o-graph parser and a full MST tensor-network frontend are outside
this implementation.

## Binder semantics

Use `finite_sum(:s, lo, hi) do s ... end` with inclusive integer bounds.
The callback receives a fresh bound integer variable. Display names are not
identities. Equality and hashing respect alpha-renaming, including nested
binders and commutative sums/products. `freevars` traverses deferred constants
as well as scalar expressions. `substitute` alpha-renames an intersecting
binder before substituting, and does not substitute through a shadowing binder.

## Normalization procedure

Construction flattens sums/products, sorts factors, combines like terms and
removes duplicate indicator factors. It does not distribute products over sums.
Normalization makes deterministic bottom-up passes:

1. Lower a deferred structure constant to its finite scalar formula.
2. Evaluate decidable ranges/deltas and annihilate zero terms.
   Remove a bound-index range indicator when it is identical to that sum's
   domain; use N ≥ 2 to recognize the context's endpoints, 0 and 1.
3. Merge phases for the same root and discard explicit integer multiples of N.
4. Propagate ordinary variable equalities, retaining the delta constraint.
5. Eliminate a sum containing `δ(s=v)` only when v is an integer expression
   independent of s; retain the inclusive range indicator on v.
6. In a complete residue domain `0:N-1`, eliminate a modular delta with
   coefficient +1 or -1 in s by substituting the representative `mod(v,N)`.
   This is valid for nonperiodic remaining factors as well.
7. Apply Fourier elimination when the only s-dependent factor is
   `ζ^(A s+B)` with established integer A,B independent of s, and the sum
   covers exactly `0:N-1`. Return `N ζ^B δ(A≡0 mod N)`.
8. Count constant sums only for concrete bounds or for the context's full
   residue domain. Simplify q-factorial/binomial boundary values.

The procedure deliberately does not cancel a nonunit coefficient modulo N,
substitute a congruence as an ordinary equality, or treat an incomplete
interval as a Fourier sum. Composite N is supported. Nonlinear Gauss sums,
general modular solving and q-Vandermonde matching remain unevaluated in this
version; adding them requires separate guarded rules and verification.

When `LinearAlgebra` is also imported, use `SymbolicTaft.normalize` or an
explicit `import AkaghefAlgebra.SymbolicTaft: normalize` to select this
operation rather than vector normalization.

## Stopping, traces and exact fallback

`normalize` returns `Normalization(expression, trace, status, passes)`.
The trace contains named rule applications and their local before/after trees.
Constructor canonicalization is implicit; this is an inspectable rewrite trace,
not a proof-assistant certificate. A `:fixed_point` says only that these rules
no longer change the expression. It does not assert confluence or a complete
normal form. No reverse rewrite or speculative q-binomial expansion is used.

`maxpasses` bounds the number of passes. `maxnodes` checks expression size
before and after each pass; an oversized candidate is discarded with its trace.
It is not an allocation or wall-time limit. Budget exits preserve the previous
expression and report `:pass_budget` or `:node_budget`.

`evaluate_exact` is a separate, explicitly requested specialization path.
Unresolved sums are enumerated only here, under a shared `maxterms` count of
visited summands (including nested sums). Exceeding it raises an error.
No approximate scalar is returned. Concrete cyclotomic coefficients have a
unique reduced representative, whereas symbolic-N residual formulas need not.

## Acceptance and independent checks

Tests cover the package import and examples, capture avoidance and alpha
equivalence, symbolic-N Fourier elimination specialized at prime and composite
N, range boundaries, nonunit modular coefficients, nonperiodic factors,
termination budgets, sparse Hopf identities and dual-basis orthogonality.
Direct powers in the exact coefficient field independently check Fourier
results. Generator expansion checks coproduct coefficients without invoking
the closed q-binomial formula. The finite scope and actual run outcome are
recorded in the validation document; no finite check is an all-N proof.
