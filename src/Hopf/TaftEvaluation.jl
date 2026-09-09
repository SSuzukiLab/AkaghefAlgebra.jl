"""Compiled exact evaluator for `SymbolicTaft.Term` expressions.

The ordinary evaluator is intentionally a small, transparent interpreter.  This
module keeps the same semantics while compiling the tree once, reusing the
root-of-unity and q-binomial values for a fixed `N`, and evaluating indicator
factors before expensive coefficients.  It is intended for residual finite
sums after symbolic reduction.
"""
module TaftEvaluation

using ..SymbolicTaft
using ..TaftExact
using ..QSeries

export FastProgram, compile_fast, evaluate_fast

abstract type FastNode end
struct FConst <: FastNode; value::Any; end
struct FVar <: FastNode; name::Symbol; slot::Int; end
struct FAdd <: FastNode; args::Vector{FastNode}; end
struct FProd <: FastNode; args::Vector{FastNode}; end
struct FPow <: FastNode; base::FastNode; exponent::FastNode; end
struct FSum <: FastNode; slot::Int; lo::FastNode; hi::FastNode; body::FastNode; end
struct FUnary <: FastNode; op::Symbol; arg::FastNode; end
struct FBinary <: FastNode; op::Symbol; left::FastNode; right::FastNode; end
struct FRange <: FastNode; value::FastNode; lo::FastNode; hi::FastNode; end
struct FPhase <: FastNode; op::Symbol; n::FastNode; exponent::FastNode; end
struct FQChoose <: FastNode; op::Symbol; n::FastNode; upper::FastNode; lower::FastNode; end
struct FQFact <: FastNode; n::FastNode; index::FastNode; end
struct FQProduct <: FastNode; n::FastNode; start::FastNode; length::FastNode; end
struct FCoeff <: FastNode; value::Any; end

"Compiled tree plus the free names it expects in the evaluation environment."
struct FastProgram
    root::FastNode
    freevars::Vector{Symbol}
    nodes::Int
end

mutable struct FastState
    N::Int
    q::Any
    qctx::Any
    powers::Dict{BigInt,Any}
    qproducts::Dict{Tuple{BigInt,BigInt},Any}
    env::Dict{Symbol,BigInt}
    slots::Vector{BigInt}
    maxterms::Int
    terms::Int
end

_const_value(x::BigInt) = x
_const_value(x::Integer) = BigInt(x)
_const_value(x::Rational) = denominator(x) == 1 ? BigInt(numerator(x)) : x
_compile(x::Integer, scopes::Dict{Symbol,Int}, nslots::Int) =
    (FConst(_const_value(x)), nslots)
_compile(x::Rational, scopes::Dict{Symbol,Int}, nslots::Int) =
    (FConst(_const_value(x)), nslots)

_integer_value(x::BigInt) = x
_integer_value(x::Integer) = BigInt(x)
function _integer_value(x::Rational)
    denominator(x) == 1 || throw(DomainError(x, "integer required"))
    BigInt(numerator(x))
end
_integer_value(x) = throw(DomainError(x, "integer required"))

function _compile(x::SymbolicTaft.Term, scopes::Dict{Symbol,Int}, nslots::Int)
    op, a = x.op, x.args
    if op == :num
        return FConst(_const_value(a[1])), nslots
    elseif op == :coeff
        length(a) == 1 || throw(ArgumentError("coefficient node expects one value"))
        return FCoeff(a[1]), nslots
    elseif op == :var
        return FVar(a[1], get(scopes, a[1], 0)), nslots
    elseif op == :structure
        return _compile(SymbolicTaft.lower(x), scopes, nslots)
    elseif op == :sum
        i, lo, hi, body = a
        clo, n1 = _compile(lo, scopes, nslots)
        chi, n2 = _compile(hi, scopes, n1)
        slot = n2 + 1
        nested = copy(scopes)
        nested[i.args[1]] = slot
        cbody, n3 = _compile(body, nested, slot)
        return FSum(slot, clo, chi, cbody), n3
    elseif op == :add
        out = FastNode[]; ns = nslots
        for t in a
            y, ns = _compile(t, scopes, ns); push!(out, y)
        end
        return FAdd(out), ns
    elseif op == :mul
        out = FastNode[]; ns = nslots
        for t in sort(collect(a); by=t -> t.op in (:range, :eq, :moddelta) ? 0 : 1)
            y, ns = _compile(t, scopes, ns); push!(out, y)
        end
        # A static one is safe to remove.  Keep static zeros in place: a sibling
        # finite sum still has the same summand-budget semantics as the baseline.
        filter!(t -> !(t isa FConst && isone(t.value)), out)
        isempty(out) && return FConst(BigInt(1)), ns
        return FProd(out), ns
    elseif op == :pow
        b, n1 = _compile(a[1], scopes, nslots)
        e, n2 = _compile(a[2], scopes, n1)
        return FPow(b, e), n2
    elseif op in (:sign, :tri)
        y, ns = _compile(a[1], scopes, nslots)
        return FUnary(op, y), ns
    elseif op in (:eq, :range, :mod)
        l, ns = _compile(a[1], scopes, nslots)
        r, ns = _compile(a[2], scopes, ns)
        if op == :range
            h, ns = _compile(a[3], scopes, ns)
            return FRange(l, r, h), ns
        end
        return FBinary(op, l, r), ns
    elseif op in (:phase, :moddelta)
        n, ns = _compile(a[1], scopes, nslots)
        p, ns = _compile(a[2], scopes, ns)
        return FPhase(op, n, p), ns
    elseif op in (:qchoose, :qchoose_inv)
        n, ns = _compile(a[1], scopes, nslots)
        u, ns = _compile(a[2], scopes, ns)
        l, ns = _compile(a[3], scopes, ns)
        return FQChoose(op, n, u, l), ns
    elseif op == :qfact
        n, ns = _compile(a[1], scopes, nslots)
        i, ns = _compile(a[2], scopes, ns)
        return FQFact(n, i), ns
    elseif op == :qproduct
        n, ns = _compile(a[1], scopes, nslots)
        s, ns = _compile(a[2], scopes, ns)
        l, ns = _compile(a[3], scopes, ns)
        return FQProduct(n, s, l), ns
    else
        throw(ArgumentError("unsupported expression: $op"))
    end
end

function _nodecount(x::FastNode)
    (x isa FConst || x isa FVar) && return 1
    x isa FAdd && return 1 + sum(_nodecount, x.args)
    x isa FProd && return 1 + sum(_nodecount, x.args)
    x isa FPow && return 1 + _nodecount(x.base) + _nodecount(x.exponent)
    x isa FSum && return 1 + _nodecount(x.lo) + _nodecount(x.hi) + _nodecount(x.body)
    x isa FUnary && return 1 + _nodecount(x.arg)
    x isa FBinary && return 1 + _nodecount(x.left) + _nodecount(x.right)
    x isa FRange && return 1 + _nodecount(x.value) + _nodecount(x.lo) + _nodecount(x.hi)
    x isa FPhase && return 1 + _nodecount(x.n) + _nodecount(x.exponent)
    x isa FQChoose && return 1 + _nodecount(x.n) + _nodecount(x.upper) + _nodecount(x.lower)
    x isa FQFact && return 1 + _nodecount(x.n) + _nodecount(x.index)
    x isa FQProduct && return 1 + _nodecount(x.n) + _nodecount(x.start) + _nodecount(x.length)
    x isa FCoeff && return 1
    error("unknown compiled node")
end

"Compile once; free variables are looked up by name at evaluation time."
function compile_fast(x::SymbolicTaft.Term)
    root, _ = _compile(x, Dict{Symbol,Int}(), 0)
    FastProgram(root, sort!(collect(SymbolicTaft.freevars(x))), _nodecount(root))
end

function _qpower(st::FastState, e::BigInt)
    r = mod(e, BigInt(st.N))
    get!(st.powers, r) do
        st.q^r
    end
end
_qpower(st::FastState, e::Integer) = _qpower(st, BigInt(e))

function _qbinomial(st::FastState, n::BigInt, k::BigInt)
    n >= 0 || throw(DomainError(n, "q-binomial upper index must be nonnegative"))
    n <= typemax(Int) || throw(DomainError(n, "q-binomial upper index must fit Int"))
    0 <= k <= n && return QSeries.qbinomial(st.qctx, Int(n), Int(k))
    BigInt(0)
end

function _qfactorial(st::FastState, r::BigInt)
    r >= 0 && r <= typemax(Int) || throw(DomainError(r, "q-factorial index must be a nonnegative Int"))
    QSeries.qfactorial(st.qctx, Int(r))
end

function _qproduct(st::FastState, a::BigInt, n::BigInt)
    n >= 0 || throw(DomainError(n, "q-product length must be nonnegative"))
    key = (a, n)
    get!(st.qproducts, key) do
        out = one(st.q)
        j = BigInt(0)
        while j < n
            out *= 1 - _qpower(st, a + j)
            j += 1
        end
        out
    end
end

function _eval(x::FConst, st); x.value end
function _eval(x::FVar, st)
    x.slot == 0 ? st.env[x.name] : st.slots[x.slot]
end
function _eval(x::FAdd, st)
    isempty(x.args) && return BigInt(0)
    out = _eval(x.args[1], st)
    for t in x.args[2:end]
        v = _eval(t, st)
        iszero(v) || (out += v)
    end
    out
end
function _eval(x::FProd, st)
    out = BigInt(1)
    for t in x.args
        v = _eval(t, st)
        iszero(v) && return v
        isone(v) || (out *= v)
    end
    out
end
function _eval(x::FPow, st)
    b, e = _eval(x.base, st), _integer_value(_eval(x.exponent, st))
    b isa Integer && e < 0 && (b = BigInt(b)//BigInt(1))
    b^e
end
function _eval(x::FSum, st)
    lo, hi = _integer_value(_eval(x.lo, st)), _integer_value(_eval(x.hi, st))
    out = BigInt(0)
    lo > hi && return out
    old = length(st.slots) >= x.slot ? st.slots[x.slot] : BigInt(0)
    length(st.slots) < x.slot && resize!(st.slots, x.slot)
    k = lo
    try
        while k <= hi
            st.terms >= st.maxterms && throw(ArgumentError("exact evaluation summand budget exhausted"))
            st.terms += 1
            st.slots[x.slot] = k
            out += _eval(x.body, st)
            k += 1
        end
    finally
        st.slots[x.slot] = old
    end
    out
end
function _eval(x::FUnary, st)
    v = _integer_value(_eval(x.arg, st))
    x.op == :tri && return div(v*(v-1), 2)
    x.op == :sign && return (-BigInt(1))^v
    error("unsupported unary expression: $(x.op)")
end
function _eval(x::FBinary, st)
    l, r = _eval(x.left, st), _eval(x.right, st)
    x.op == :eq && return Int(l == r)
    x.op == :mod && return mod(_integer_value(l), _integer_value(r))
    error("unsupported binary expression: $(x.op)")
end
function _eval(x::FPhase, st)
    n = _integer_value(_eval(x.n, st)); n == st.N || throw(ArgumentError("all roots must use one N"))
    p = _integer_value(_eval(x.exponent, st))
    x.op == :phase && return _qpower(st, p)
    x.op == :moddelta && return Int(mod(p, BigInt(st.N)) == 0)
    error("unsupported phase expression: $(x.op)")
end
function _eval(x::FQChoose, st)
    n = _integer_value(_eval(x.n, st)); n == st.N || throw(ArgumentError("all roots must use one N"))
    u, l = _integer_value(_eval(x.upper, st)), _integer_value(_eval(x.lower, st))
    v = _qbinomial(st, u, l)
    iszero(v) && return v
    x.op == :qchoose_inv ? _qpower(st, -l*(u-l))*v : v
end
function _eval(x::FQFact, st)
    n = _integer_value(_eval(x.n, st)); n == st.N || throw(ArgumentError("all roots must use one N"))
    _qfactorial(st, _integer_value(_eval(x.index, st)))
end
function _eval(x::FQProduct, st)
    n = _integer_value(_eval(x.n, st)); n == st.N || throw(ArgumentError("all roots must use one N"))
    _qproduct(st, _integer_value(_eval(x.start, st)), _integer_value(_eval(x.length, st)))
end
function _eval(x::FCoeff, st)
    v = x.value
    v isa Union{Integer,Rational} && return v
    # QSeries owns specialization of both Cyclotomic and canonical Q(q)
    # coefficients, including the coefficient-parent/order checks.
    QSeries.specialize_q(v, st.N)
end

function _eval(x::FRange, st)
    v = _eval(x.value, st)
    lo = _eval(x.lo, st)
    hi = _eval(x.hi, st)
    Int(lo <= v <= hi)
end

function _prepare_env(env::AbstractDict)
    all(v -> v isa Integer, values(env)) || throw(ArgumentError("free variables must specialize to integers"))
    Dict{Symbol,BigInt}(k => BigInt(v) for (k,v) in env)
end

function evaluate_fast(program::FastProgram, n::Integer, env::AbstractDict; maxterms::Int=1000000)
    n >= 2 || throw(DomainError(n, "N must be ≥ 2"))
    n <= typemax(Int) || throw(DomainError(n, "N must fit Int"))
    maxterms >= 0 || throw(ArgumentError("maxterms must be nonnegative"))
    e = _prepare_env(env)
    if haskey(e, :N)
        e[:N] == BigInt(n) || throw(ArgumentError("environment N does not match evaluator N"))
    else
        e[:N] = BigInt(n)
    end
    q = TaftExact.zeta(Int(n))
    qctx = QSeries.qcontext(q; max_n=Int(n))
    st = FastState(Int(n), q, qctx, Dict{BigInt,Any}(BigInt(0) => one(q)),
        Dict{Tuple{BigInt,BigInt},Any}(),
        e, BigInt[], maxterms, 0)
    _eval(program.root, st)
end

evaluate_fast(program::FastProgram, n::Integer; env=Dict{Symbol,Int}(), maxterms::Int=1000000) =
    evaluate_fast(program, n, env; maxterms=maxterms)

function evaluate_fast(x::SymbolicTaft.Term, c::SymbolicTaft.TaftContext,
                       env::AbstractDict; maxterms::Int=1000000)
    nvalue = SymbolicTaft.evaluate_exact(c.N, env; maxterms=0)
    n = _integer_value(nvalue)
    haskey(env, :N) && BigInt(env[:N]) != n &&
        throw(ArgumentError("environment N does not match context N"))
    evaluate_fast(compile_fast(x), n, env; maxterms=maxterms)
end

evaluate_fast(x::SymbolicTaft.Term, c::SymbolicTaft.TaftContext;
              env=Dict{Symbol,Int}(), maxterms::Int=1000000) =
    evaluate_fast(x, c, env; maxterms=maxterms)

end # module TaftEvaluation
