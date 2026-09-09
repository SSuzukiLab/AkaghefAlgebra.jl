# Lightweight root-of-unity q-binomial rewrites (QA1--QA7, QA11).
# These rules are deliberately syntactic: an unmatched sum is left intact.
_qbound(bounds,x) = x.op == :var ? get(bounds,x.args[1],nothing) : nothing
_qinteger_value(x) = isnum(x) && denominator(x.args[1]) == 1 ? numerator(x.args[1]) : nothing
_qnonnegative(x,c,bounds=nothing) = begin
    v = _qinteger_value(x)
    (v !== nothing && v >= 0) || x == c.N || x == c.N - 1 ||
        (bounds !== nothing && x.op == :var &&
         (_qbound(bounds,x) !== nothing && _qnonnegative(_qbound(bounds,x)[1],c,bounds)))
end
_qpositive(x,c) = begin
    v = _qinteger_value(x)
    (v !== nothing && v >= 1) || x == c.N || x == c.N - 1
end
_qgeq(a,b,c,bounds=nothing) = begin
    av, bv = _qinteger_value(a), _qinteger_value(b)
    (av !== nothing && bv !== nothing && av >= bv) || a == b ||
        (b == term(0) && _qnonnegative(a,c,bounds)) ||
        (bounds !== nothing && b.op == :var && _qbound(bounds,b) !== nothing &&
         a == _qbound(bounds,b)[2])
end
_qfullbound(b,c) = b !== nothing && b[1] == term(0) && b[2] == c.N-1
_qphase(x,c) = x.op == :phase && x.args[1] == c.N ? x.args[2] : nothing
_qfactors(x) = x.op == :mul ? collect(x.args) : [x]
mutable struct _QRewriteBudget
    maxnodes::Int
    maxexpansion::Int
    nodes::Int
    expansions::Int
end
function _qvisit!(budget::_QRewriteBudget)
    budget.nodes += 1
    budget.nodes <= budget.maxnodes
end
function _qexpand_charge!(budget::_QRewriteBudget,n::Integer)
    n < 0 && return false
    budget.expansions > budget.maxexpansion - n && return false
    budget.expansions += n
    true
end
function _qmul_expand(a::Term,b::Term,budget::_QRewriteBudget)
    if a.op==:add
        _qexpand_charge!(budget,length(a.args)) || return product((a,b))
        return addition((_qmul_expand(t,b,budget) for t in a.args))
    elseif b.op==:add
        _qexpand_charge!(budget,length(b.args)) || return product((a,b))
        return addition((_qmul_expand(a,t,budget) for t in b.args))
    end
    product((a,b))
end
function _qexpand(x::Term,budget::_QRewriteBudget)
    _qvisit!(budget) || return x
    x.op==:add && return addition((_qexpand(t,budget) for t in x.args))
    x.op==:mul && return foldl((a,b)->_qmul_expand(a,b,budget),
                               (_qexpand(t,budget) for t in x.args);init=term(1))
    if x.op==:pow && x.args[2] >= 0
        n = x.args[2]
        n <= budget.maxexpansion || return x
        _qexpand_charge!(budget,n) || return x
        return foldl((a,b)->_qmul_expand(a,b,budget),
                     (_qexpand(x.args[1],budget) for _ in BigInt(1):n);init=term(1))
    end
    x
end

function _qstrip_range(fs,i)
    out=Term[]; found=false; lo=hi=nothing
    for f in fs
        if f.op == :range && f.args[1] == i && !found
            found=true; lo,hi=f.args[2],f.args[3]
        elseif f.op == :range && f.args[3] == i && f.args[2] == term(0) && !found
            # A common network gate is `inrange(base,0,i)`, i.e. base<=i.
            found=true; lo,hi=f.args[1],i
        else
            push!(out,f)
        end
    end
    out, found, lo, hi
end

function _qbinomial_big(n::BigInt,k::BigInt,budget::_QRewriteBudget)
    (n < 0 || k < 0 || k > n) && return nothing
    k = min(k,n-k)
    _qexpand_charge!(budget,k) || return nothing
    out = BigInt(1)
    j = BigInt(1)
    while j <= k
        out = div(out * (n-k+j), j)
        j += 1
    end
    out
end

function _qrewrite_choose(x::Term,c::TaftContext,bounds,rule,budget::_QRewriteBudget)
    N,n,k=x.args
    N == c.N || return x
    if n == c.N-1
        b=_qbound(bounds,k)
        if _qfullbound(b,c)
            return Term(:sign,(k,))*phase(c,-Term(:tri,(k+1,)))
        end
    elseif n == c.N
        b=_qbound(bounds,k)
        if b !== nothing && isnum(b[1],1) && b[2] == c.N-1
            return term(0)
        elseif b !== nothing && isnum(b[1],0) && b[2] == c.N-1
            return eqdelta(k,0)
        elseif isnum(k,0) || k == c.N
            return term(1)
        end
    end
    # QA11, specialized to numeric n,k.  Polynomial evaluation is not used.
    if isnum(N) && isnum(n) && isnum(k)
        nn,kk,nnN = _qinteger_value(n), _qinteger_value(k), _qinteger_value(N)
        (nn === nothing || kk === nothing || nnN === nothing) && return x
        if nn < 0 || kk < 0 || kk > nn
            return term(0)
        end
        a,b=divrem(nn,nnN); d,e=divrem(kk,nnN)
        e > b && return term(0)
        coeff = _qbinomial_big(a,d,budget)
        coeff === nothing && return x
        return (e == 0 || e == b) ? term(coeff) : term(coeff)*qchoose(c,b,e)
    end
    x
end

function _qrewrite_leaf(x::Term,c::TaftContext,bounds,rule,budget::_QRewriteBudget)
    if x.op == :qchoose && rule in (:pascal,:pascal1)
        N,n,k=x.args
        N == c.N && _qpositive(n,c) || return x
        return qchoose(c,n-1,k-1)+phase(c,k)*qchoose(c,n-1,k)
    elseif x.op == :qchoose && rule == :pascal2
        N,n,k=x.args
        N == c.N && _qpositive(n,c) || return x
        return phase(c,n-k)*qchoose(c,n-1,k-1)+qchoose(c,n-1,k)
    elseif x.op == :qchoose
        # A named summation rule must see its coefficient pattern intact;
        # endpoint/Lucas leaf rules are therefore reserved for `:auto` and
        # an explicit Lucas request.
        rule in (:auto,:qlucas,:lucas) || return x
        return _qrewrite_choose(x,c,bounds,rule,budget)
    elseif x.op == :qchoose_inv
        N,n,k=x.args
        N == c.N || return x
        # QA2: [n k]_(q^-1) = q^(-k(n-k)) [n k]_q.
        return phase(c,-k*(n-k))*qchoose(c,n,k)
    elseif x.op == :qproduct
        return x
    end
    x
end

function _qfinite_binomial_sum(i,lo,hi,b,c,bounds,budget::_QRewriteBudget)
    lo == term(0) || return nothing
    fs,hasrange,rlo,rhi=_qstrip_range(_qfactors(b),i)
    # A full q-binomial theorem has exactly sign(i), phase(tri(i)+a*i), qchoose(n,i).
    signf=nothing; phasef=nothing; choosef=nothing
    for f in fs
        if f.op==:sign && f.args[1]==i; signf=f
        elseif f.op==:phase && f.args[1]==c.N; phasef=f
        elseif f.op==:qchoose && f.args[1]==c.N && f.args[3]==i; choosef=f
        else; return nothing
        end
    end
    (signf===nothing || phasef===nothing || choosef===nothing || hasrange) && return nothing
    n=choosef.args[2]
    hi == n && _qnonnegative(n,c,bounds) && _qnonnegative(hi,c,bounds) || return nothing
    p=phasef.args[2]
    # tri(i)+a*i is represented literally; this guard prevents arbitrary sums.
    tri_i=Term(:tri,(i,))
    pieces=p.op==:add ? collect(p.args) : [p]
    j=findfirst(==(tri_i),pieces); j===nothing && return nothing
    deleteat!(pieces,j)
    length(pieces)==1 || return nothing
    ab=affine(pieces[1],i); ab===nothing && return nothing
    ab[2] == term(0) && integerexpr(ab[1]) || return nothing
    qproduct(c,ab[1],n)
end

function _qhockey_sum(i,lo,hi,b,c,bounds,budget::_QRewriteBudget)
    fs,hasrange,rlo,rhi=_qstrip_range(_qfactors(b),i)
    qf=nothing; pf=nothing; rest=Term[]
    for f in fs
        if f.op==:qchoose && f.args[1]==c.N && f.args[2]==i
            qf===nothing || return nothing
            qf=f
        elseif f.op==:phase && f.args[1]==c.N; pf=f
        else; push!(rest,f)
        end
    end
    (qf===nothing || pf===nothing) && return nothing
    base=qf.args[3]; pfab=affine(pf.args[2],i)
    pfab===nothing && return nothing
    # QA2 followed by QA5: q^(b^2-a*b-a) [a b]_q.
    _qexpand(pfab[1],budget) == _qexpand(-(base+1),budget) || return nothing
    # Either the sum starts at b, or a range factor supplies a<=i<=M.
    valid = (_qnonnegative(base,c,bounds) && _qnonnegative(hi,c,bounds) &&
             _qgeq(hi,base,c,bounds)) &&
             (lo == base || (lo == term(0) && hasrange && rlo==base && rhi==i))
    valid || return nothing
    any(depends(f,i) for f in rest) && return nothing
    p0=pfab[2]-base^2
    correction=p0-base-(base+1)*(hi-base)
    # At the common root-of-unity endpoint M=N-1, QA6 turns this into
    # δ_{b,N-1} when the lower index has the certified full residue range.
    endpoint = hi == c.N-1 && _qfullbound(_qbound(bounds,base),c)
    tail = endpoint ? eqdelta(base,hi) : qchoose(c,hi+1,base+1)
    product(rest)*phase(c,correction)*tail
end

function _qvander_sum(i,lo,hi,b,c,bounds,budget::_QRewriteBudget)
    lo == term(0) || return nothing
    fs=_qfactors(b); pf=nothing; q1=nothing; q2=nothing
    for f in fs
        if f.op==:phase && f.args[1]==c.N; pf=f
        elseif f.op==:qchoose && f.args[1]==c.N && f.args[3]==i; q1=f
        end
    end
    # Match the second coefficient by its explicit r-i shape.
    for f in fs
        if f.op==:qchoose && f.args[1]==c.N && f.args[3].op==:add && depends(f.args[3],i)
            ab=affine(f.args[3],i)
            ab !== nothing && ab[1]==term(-1) && (q2=f)
        end
    end
    (pf===nothing || q1===nothing || q2===nothing) && return nothing
    length(fs)==3 || return nothing
    m=q1.args[2]
    # q2 is [n; r-i], and the phase is (m-i)(r-i).
    nterm=affine(q2.args[3],i); nterm===nothing && return nothing
    nterm[1]==term(-1) || return nothing
    r=nterm[2]; n=q2.args[2]
    pf.args[2] == (m-i)*(r-i) || return nothing
    hi==m && _qnonnegative(m,c,bounds) && _qnonnegative(n,c,bounds) &&
        _qnonnegative(hi,c,bounds) || return nothing
    qchoose(c,m+n,r)
end

function _qrewrite_sum(x::Term,c::TaftContext,bounds,rule,budget::_QRewriteBudget)
    i,lo,hi,b=x.args
    rule in (:auto,:hockey,:hockeystick) && begin
        y=_qhockey_sum(i,lo,hi,b,c,bounds,budget); y!==nothing && return y
    end
    rule in (:auto,:qbinomial,:finite_qbinomial) && begin
        y=_qfinite_binomial_sum(i,lo,hi,b,c,bounds,budget); y!==nothing && return y
    end
    rule in (:auto,:vandermonde,:vander) && begin
        y=_qvander_sum(i,lo,hi,b,c,bounds,budget); y!==nothing && return y
    end
    x
end

function _qrewrite(x::Term,c::TaftContext,bounds,rule,budget::_QRewriteBudget)
    _qvisit!(budget) || return x
    if x.op == :sum
        i,lo,hi,b=x.args
        child=copy(bounds); child[i.args[1]]=(lo,hi)
        y=Term(:sum,(i,lo,hi,_qrewrite(b,c,child,rule,budget)))
        z=_qrewrite_sum(y,c,bounds,rule,budget)
        return z == y ? y : _qrewrite_leaf(z,c,bounds,rule,budget)
    end
    y=rebuild(x.op,maparg(t->_qrewrite(t,c,bounds,rule,budget),x.args))
    _qrewrite_leaf(y,c,bounds,rule,budget)
end

"""Apply guarded root-of-unity q-binomial rules; unresolved terms are retained."""
function qrewrite(x::Term,c::TaftContext;bounds=Dict(),rule=:auto,
                  maxnodes::Int=100000,maxexpansion::Int=10000)
    maxnodes >= 1 && maxexpansion >= 0 ||
        throw(ArgumentError("invalid qrewrite budget"))
    b=Dict{Symbol,Tuple{Term,Term}}()
    for (k,v) in bounds
        length(v)==2 || throw(ArgumentError("bounds values must be (lo,hi)"))
        b[k isa Symbol ? k : k.args[1]]=(term(v[1]),term(v[2]))
    end
    _qrewrite(x,c,b,rule,_QRewriteBudget(maxnodes,maxexpansion,0,0))
end
simplify_q(x::Term,c::TaftContext;kwargs...) = qrewrite(x,c;kwargs...)
