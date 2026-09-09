"""Generic q-series coefficient operations.

`QSeries` keeps the coefficient layer independent of Taft's symbolic AST.  The
canonical generic coefficient field is the rational function field `Q(q)`
provided by AbstractAlgebra.  Concrete q-series operations only require the
usual `zero`, `one`, arithmetic, and integer-power methods of their coefficient
element; they deliberately do not require `Number`.
"""
module QSeries

using AbstractAlgebra
using ..TaftExact: Cyclotomic, zeta

export QQq, QQqField, QQqRing, q, qgenerator, qring, qfield,
       QSeriesContext, qcontext, qnumber, qpochhammer, qbinomial,
       qfactorial, specialize_q, specialize, evaluate_at_root,
       clear_cache!

# ---------------------------------------------------------------------------
# Canonical Q(q) backend
# ---------------------------------------------------------------------------

const QQqRing = first(polynomial_ring(QQ, "q"))
const _QQqPolynomialGenerator = gen(QQqRing)
const QQqField = fraction_field(QQqRing)
const q = QQqField(_QQqPolynomialGenerator)
const QQq = typeof(q)

"""Return the canonical `(Q(q), q)` pair.

The pair is stable for the life of the loaded module.  This matters for
specialization: elements from another, independently-created AbstractAlgebra
parent ring are not silently treated as the canonical q field.
"""
qring() = (QQqField, q)
qfield() = QQqField
qgenerator() = q

# ---------------------------------------------------------------------------
# Bounded per-context cache
# ---------------------------------------------------------------------------

struct QSeriesContext{T}
    q::T
    max_n::Int
    qnumber_cache::Dict{Int,T}
    qfactorial_cache::Dict{Int,T}
    qbinomial_rows::Dict{Int,Vector{T}}
end

function QSeriesContext(q; max_n::Integer=64)
    lim = Int(max_n)
    lim >= 0 || throw(ArgumentError("max_n must be nonnegative"))
    T = typeof(q)
    ctx = QSeriesContext{T}(q, lim, Dict{Int,T}(), Dict{Int,T}(), Dict{Int,Vector{T}}())
    ctx.qnumber_cache[0] = zero(q)
    ctx.qfactorial_cache[0] = one(q)
    ctx.qbinomial_rows[0] = [one(q)]
    ctx
end

qcontext(q; max_n::Integer=64) = QSeriesContext(q; max_n=max_n)

function clear_cache!(ctx::QSeriesContext)
    empty!(ctx.qnumber_cache)
    empty!(ctx.qfactorial_cache)
    empty!(ctx.qbinomial_rows)
    ctx.qnumber_cache[0] = zero(ctx.q)
    ctx.qfactorial_cache[0] = one(ctx.q)
    ctx.qbinomial_rows[0] = [one(ctx.q)]
    ctx
end

_cacheable(ctx::QSeriesContext, n::Integer) = 0 <= n <= ctx.max_n

function _qnumber_uncached(qv, n::Int)
    out = zero(qv)
    power = one(qv)
    for _ in 1:n
        out += power
        power *= qv
    end
    out
end

"""Gaussian q-number `[n]_q = 1 + q + ... + q^(n-1)`."""
function qnumber(ctx::QSeriesContext, n::Integer)
    n >= 0 || throw(ArgumentError("qnumber requires n ≥ 0"))
    ni = Int(n)
    if haskey(ctx.qnumber_cache, ni)
        return ctx.qnumber_cache[ni]
    end
    out = _qnumber_uncached(ctx.q, ni)
    _cacheable(ctx, ni) && (ctx.qnumber_cache[ni] = out)
    out
end

qnumber(qv, n::Integer) = qnumber(QSeriesContext(qv; max_n=0), n)

function _qbinomial_row!(ctx::QSeriesContext, n::Int)
    haskey(ctx.qbinomial_rows, n) && return ctx.qbinomial_rows[n]
    # Find the largest cached row below n.  The usual use case starts at row
    # zero; this also keeps a manually populated context useful.
    start = 0
    for m in keys(ctx.qbinomial_rows)
        m <= n && m >= start && (start = m)
    end
    row = ctx.qbinomial_rows[start]
    for m in (start + 1):n
        nxt = [zero(ctx.q) for _ in 0:m]
        for k in 0:m
            k >= 1 && (nxt[k + 1] += row[k])
            k <= m - 1 && (nxt[k + 1] += ctx.q^k * row[k + 1])
        end
        row = nxt
        _cacheable(ctx, m) && (ctx.qbinomial_rows[m] = row)
    end
    row
end

"""Gaussian q-binomial coefficient, computed by q-Pascal recurrence.

`n < 0` is unsupported.  For a nonnegative `n`, an out-of-range `k` is the
usual zero coefficient and is returned without forming a quotient.
"""
function qbinomial(ctx::QSeriesContext, n::Integer, k::Integer)
    n >= 0 || throw(ArgumentError("qbinomial requires n ≥ 0"))
    (k < 0 || k > n) && return zero(ctx.q)
    row = _qbinomial_row!(ctx, Int(n))
    row[Int(k) + 1]
end

qbinomial(qv, n::Integer, k::Integer) = qbinomial(QSeriesContext(qv; max_n=0), n, k)

"""q-factorial, defined as `∏_{j=1}^n [j]_q`."""
function qfactorial(ctx::QSeriesContext, n::Integer)
    n >= 0 || throw(ArgumentError("qfactorial requires n ≥ 0"))
    ni = Int(n)
    haskey(ctx.qfactorial_cache, ni) && return ctx.qfactorial_cache[ni]
    # Reuse the nearest cached prefix to avoid recomputing q-numbers.
    start = 0
    out = one(ctx.q)
    for m in keys(ctx.qfactorial_cache)
        m <= ni && m >= start && (start = m; out = ctx.qfactorial_cache[m])
    end
    for m in (start + 1):ni
        out *= qnumber(ctx, m)
        _cacheable(ctx, m) && (ctx.qfactorial_cache[m] = out)
    end
    out
end

qfactorial(qv, n::Integer) = qfactorial(QSeriesContext(qv; max_n=0), n)

function qpochhammer(ctx::QSeriesContext, z, n::Integer)
    n >= 0 || throw(ArgumentError("qpochhammer requires n ≥ 0"))
    out = one(ctx.q)
    power = one(ctx.q)
    for _ in 1:n
        out *= one(ctx.q) - z * power
        power *= ctx.q
    end
    out
end

qpochhammer(qv, z, n::Integer) = qpochhammer(QSeriesContext(qv; max_n=0), z, n)

# Convenient shorthand for `(q;q)_n`.
qpochhammer(qv, n::Integer) = qpochhammer(qv, qv, n)
qpochhammer(ctx::QSeriesContext, n::Integer) = qpochhammer(ctx, ctx.q, n)

# ---------------------------------------------------------------------------
# Partial Q(q) -> Q(ζ_N) specialization
# ---------------------------------------------------------------------------

function _check_order(N::Integer)
    n = Int(N)
    n > 0 || throw(ArgumentError("cyclotomic order N must be positive"))
    n
end

function _eval_polynomial(poly, root::Cyclotomic{N}) where {N}
    out = zero(root)
    # Horner evaluation avoids constructing a second symbolic expression.
    d = degree(poly)
    d < 0 && return out
    for i in d:-1:0
        out = out * root + Cyclotomic{N}(coeff(poly, i))
    end
    out
end

function _canonical_field_element(x)
    x isa QQq || throw(ArgumentError("specialize_q expects an element of the canonical Q(q) field"))
    parent(x) === QQqField ||
        throw(ArgumentError("specialize_q rejects elements from a different Q(q) parent ring"))
    x
end

"""Evaluate a canonical `Q(q)` rational function at the primitive N-th root.

The numerator and denominator are already reduced by AbstractAlgebra's
fraction-field element.  The denominator is evaluated separately in
`Cyclotomic{N}`; a remaining zero denominator raises `ArgumentError`.
"""
function specialize_q(x::QQq, N::Integer)
    n = _check_order(N)
    f = _canonical_field_element(x)
    root = zeta(n)
    den = _eval_polynomial(denominator(f), root)
    iszero(den) && throw(ArgumentError("Q(q) denominator vanishes at a primitive $(n)-th root"))
    _eval_polynomial(numerator(f), root) / den
end

specialize(x::QQq, N::Integer) = specialize_q(x, N)
evaluate_at_root(x::QQq, N::Integer) = specialize_q(x, N)

function specialize_q(x::Cyclotomic{M}, N::Integer) where {M}
    n = _check_order(N)
    Int(M) == n || throw(ArgumentError("cannot mix Cyclotomic orders $(M) and $(n)"))
    x
end

specialize(x::Cyclotomic, N::Integer) = specialize_q(x, N)
evaluate_at_root(x::Cyclotomic, N::Integer) = specialize_q(x, N)

"""Specialize at an explicitly supplied primitive root of fixed order."""
function specialize_q(x::QQq, root::Cyclotomic{N}) where {N}
    root == zeta(N) || throw(ArgumentError("specialization root must be the canonical primitive root"))
    specialize_q(x, Int(N))
end

specialize(x::QQq, root::Cyclotomic) = specialize_q(x, root)
evaluate_at_root(x::QQq, root::Cyclotomic) = specialize_q(x, root)

end # module QSeries
