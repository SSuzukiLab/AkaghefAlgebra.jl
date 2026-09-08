"""Lazy, exact finite-sum calculus for the right-coproduct Taft convention.

All free variables denote integers. `TaftContext(N)` asserts N ≥ 2 and a
primitive N-th root of unity; specialization checks these assumptions.
"""
module SymbolicTaft

using ..TaftExact: zeta, qbinomial, qfactorial
export Term, integer, TaftContext, phase, eqdelta, moddelta, inrange,
       qchoose, qfact, finite_sum, substitute, freevars, normalize,
       Normalization, evaluate_exact, structure, lower, nodecount

struct Term
    op::Symbol
    args::Tuple
end
term(x::Term) = x
term(x::Integer) = Term(:num, (BigInt(x)//BigInt(1),))
term(x::Rational) = Term(:num, (BigInt(numerator(x))//BigInt(denominator(x)),))
integer(name::Symbol) = Term(:var, (name,))
isnum(x::Term) = x.op == :num
isnum(x::Term, n) = isnum(x) && x.args[1] == n

# Alpha-equivalence is structural equality, including nested shadowing.
function key(x::Term, scope=Symbol[])
    if x.op == :var
        j = findlast(==(x.args[1]), scope)
        return j === nothing ? "v:" * string(x.args[1]) : "b:" * string(length(scope)-j)
    elseif x.op == :sum
        i, lo, hi, b = x.args
        return "sum(" * key(lo,scope) * "," * key(hi,scope) * "," * key(b,[scope;i.args[1]]) * ")"
    end
    ks=[argkey(a,scope) for a in x.args]
    x.op in (:add,:mul) && sort!(ks)
    return string(x.op) * "(" * join(ks, ",") * ")"
end
argkey(a,scope) = a isa Term ? key(a,scope) : a isa Tuple ? "["*join((argkey(t,scope) for t in a),",")*"]" : repr(a)
maparg(f,a) = a isa Term ? f(a) : a isa Tuple ? Tuple(maparg(f,t) for t in a) : a
Base.:(==)(x::Term,y::Term) = key(x) == key(y)
Base.isequal(x::Term,y::Term) = x == y
Base.hash(x::Term,h::UInt) = hash(key(x),h)
Base.iszero(x::Term) = isnum(x,0)
Base.isone(x::Term) = isnum(x,1)

function product(xs)
    flat = Term[]; c = BigInt(1)//BigInt(1)
    for x in xs
        for t in (x.op == :mul ? x.args : (x,))
            if isnum(t); c *= t.args[1]; else; push!(flat,t); end
        end
    end
    iszero(c) && return term(0)
    # Indicator factors are idempotent, other factors are not.
    seen = Set{Term}(); out = Term[]
    for t in flat
        if t.op in (:eq,:moddelta,:range)
            t in seen && continue
            push!(seen,t)
        end
        push!(out,t)
    end
    sort!(out;by=key)
    c != 1 && pushfirst!(out,term(c))
    isempty(out) ? term(1) : length(out)==1 ? out[1] : Term(:mul,Tuple(out))
end
function addition(xs)
    cs = Dict{Term,Rational{BigInt}}()
    for x in xs
        for t in (x.op == :add ? x.args : (x,))
            c,b = isnum(t) ? (t.args[1],term(1)) :
                t.op == :mul && isnum(t.args[1]) ? (t.args[1].args[1],product(t.args[2:end])) : (BigInt(1)//BigInt(1),t)
            cs[b] = get(cs,b,BigInt(0)//BigInt(1)) + c
        end
    end
    out = [product((term(c),b)) for (b,c) in cs if c != 0]
    sort!(out;by=key)
    isempty(out) ? term(0) : length(out)==1 ? out[1] : Term(:add,Tuple(out))
end
Base.:+(x::Term,y::Term) = addition((x,y))
Base.:*(x::Term,y::Term) = product((x,y))
Base.:-(x::Term) = term(-1)*x
Base.:-(x::Term,y::Term) = x+(-y)
for op in (:+,:-,:*)
    @eval Base.$op(x::Term,y::Union{Integer,Rational}) = Base.$op(x,term(y))
    @eval Base.$op(x::Union{Integer,Rational},y::Term) = Base.$op(term(x),y)
end
function Base.:^(x::Term,n::Integer)
    n == 0 && return term(1)
    n == 1 && return x
    isnum(x) && return term(x.args[1]^n)
    Term(:pow,(x,BigInt(n)))
end
Base.:/(x::Term,y::Term) = x * (y^(-1))
Base.inv(x::Term) = Base.:^(x,-1)
Base.:/(x::Term,y::Union{Integer,Rational}) = x/term(y)
Base.:/(x::Union{Integer,Rational},y::Term) = term(x)/y

struct TaftContext
    N::Term
    function TaftContext(n)
        t = term(n)
        (t.op == :var || (isnum(t) && denominator(t.args[1])==1 && t.args[1]>=2)) ||
            throw(ArgumentError("N must be an integer variable or an integer ≥ 2"))
        new(t)
    end
end
phase(c::TaftContext,p) = Term(:phase,(c.N,term(p)))
eqdelta(a,b=0) = Term(:eq,(term(a),term(b)))
moddelta(c::TaftContext,a) = Term(:moddelta,(c.N,term(a)))
inrange(a,lo,hi) = Term(:range,(term(a),term(lo),term(hi)))
qchoose(c::TaftContext,n,k) = Term(:qchoose,(c.N,term(n),term(k)))
qfact(c::TaftContext,n) = Term(:qfact,(c.N,term(n)))

"""Create a fresh bound integer index. Bounds are inclusive; empty sums are zero."""
function finite_sum(f::Function,name::Symbol,lo,hi)
    integerexpr(term(lo)) && integerexpr(term(hi)) || throw(ArgumentError("sum bounds must be integer expressions"))
    i = integer(gensym(name))
    Term(:sum,(i,term(lo),term(hi),term(f(i))))
end
function freevars(x::Term)
    x.op == :var && return Set([x.args[1]])
    if x.op == :sum
        i,lo,hi,b = x.args
        return union(freevars(lo),freevars(hi),setdiff(freevars(b),Set([i.args[1]])))
    end
    out = Set{Symbol}()
    maparg(t->(union!(out,freevars(t));t),x.args)
    out
end
depends(x,i) = i.args[1] in freevars(x)
"""Capture-avoiding substitution of a free integer variable."""
function substitute(x::Term,i::Term,v)
    i.op == :var || throw(ArgumentError("substitution target must be a variable"))
    y = term(v)
    x == i && return y
    if x.op == :sum
        j,lo,hi,b = x.args
        lo2,hi2 = substitute(lo,i,y),substitute(hi,i,y)
        j == i && return Term(:sum,(j,lo2,hi2,b))
        if depends(y,j)
            fresh = integer(gensym(:bound))
            b = substitute(b,j,fresh); j = fresh
        end
        return Term(:sum,(j,lo2,hi2,substitute(b,i,y)))
    end
    rebuild(x.op,maparg(t->substitute(t,i,y),x.args))
end
function rebuild(op,args)
    op == :add && return addition(args)
    op == :mul && return product(args)
    op == :pow && return args[1]^args[2]
    Term(op,args)
end
argcount(a) = a isa Term ? nodecount(a) : a isa Tuple ? sum(argcount,a;init=0) : 0
nodecount(x::Term) = 1 + sum(argcount,x.args;init=0)

# Return a,b with p=a*i+b only when established by the expression tree.
function affine(p::Term,i::Term)
    !depends(p,i) && return (term(0),p)
    p == i && return (term(1),term(0))
    if p.op == :add
        pieces = [affine(a,i) for a in p.args]
        any(isnothing,pieces) && return nothing
        return (addition(t[1] for t in pieces),addition(t[2] for t in pieces))
    elseif p.op == :mul
        dep = [a for a in p.args if depends(a,i)]
        length(dep)==1 || return nothing
        ab = affine(only(dep),i); ab === nothing && return nothing
        c = product(a for a in p.args if !depends(a,i))
        return (c*ab[1],c*ab[2])
    end
    nothing
end

"""Deferred structure constant. Arguments are PBW index pairs in 0:N-1.

Kinds: :mul(input,input,output), :comul(input,left,right),
:antipode(input,output), :counit(input), :unit(output),
:pair(dual a^l b^t,input), :dual_basis(primal,dual monomial).
"""
function structure(c::TaftContext,kind::Symbol,pairs::Tuple...)
    count = Dict(:mul=>3,:comul=>3,:antipode=>2,:counit=>1,:unit=>1,:pair=>2,:dual_basis=>2)
    haskey(count,kind) && length(pairs)==count[kind] || throw(ArgumentError("invalid structure constant signature"))
    all(p->length(p)==2,pairs) || throw(ArgumentError("each basis index needs two coordinates"))
    Term(:structure,(c.N,kind,(Tuple(term(x) for x in p) for p in pairs)...))
end
function lower(x::Term)
    x.op == :structure || return x
    N,kind,ps... = x.args; c=TaftContext(N)
    gate = product(inrange(t,0,N-1) for p in ps for t in p)
    if kind == :mul
        (r,s),(t,u),(v,w)=ps
        return gate*phase(c,s*t)*eqdelta(v,r+t)*moddelta(c,w-s-u)
    elseif kind == :comul
        (r,s),(j,u),(k,v)=ps
        return gate*qchoose(c,r,j)*eqdelta(k,r-j)*moddelta(c,u-s)*moddelta(c,v-s-j)
    elseif kind == :antipode
        (r,s),(t,u)=ps
        return gate*Term(:sign,(r,))*phase(c,-s*r-Term(:tri,(r,)))*eqdelta(t,r)*moddelta(c,u+s+r)
    elseif kind == :counit
        (r,s)=only(ps)
        return gate*eqdelta(r,0)
    elseif kind == :unit
        (r,s)=only(ps)
        return gate*eqdelta(r,0)*eqdelta(s,0)
    elseif kind == :dual_basis
        (r,s),(l,t)=ps
        return gate*eqdelta(t,r)*phase(c,-s*l)/(N*qfact(c,r))
    else
        (l,t),(r,s)=ps
        return gate*eqdelta(t,r)*qfact(c,r)*phase(c,l*s)
    end
end

struct Normalization
    expression::Term
    trace::Vector{NamedTuple{(:rule,:before,:after),Tuple{Symbol,Term,Term}}}
    status::Symbol
    passes::Int
end
function record!(trace,rule,before,after)
    before == after || push!(trace,(rule=rule,before=before,after=after))
    after
end

function localrule(x::Term,trace,context)
    op,a=x.op,x.args
    if op == :structure
        return record!(trace,:structure_constant,x,lower(x))
    elseif op == :eq
        a[1]==a[2] && return record!(trace,:equal_indices,x,term(1))
        all(isnum,a) && return record!(trace,:distinct_indices,x,term(0))
    elseif op == :range
        t,lo,hi=a
        all(isnum,a) && return record!(trace,:range_check,x,term(lo.args[1]<=t.args[1]<=hi.args[1] ? 1 : 0))
        if context !== nothing && isnum(lo,0) && hi==context.N-1 &&
                (isnum(t,0) || isnum(t,1) || t==hi)
            return record!(trace,:context_range,x,term(1))
        end
    elseif op == :phase || op == :moddelta
        N,p=a
        if isnum(p,0)
            return record!(trace,:zero_exponent,x,term(1))
        elseif isnum(N) && isnum(p) && denominator(p.args[1])==1
            r=mod(p.args[1],N.args[1])
            y=op == :phase ? (r==0 ? term(1) : Term(op,(N,term(r)))) : term(r==0 ? 1 : 0)
            return record!(trace,:mod_reduce,x,y)
        end
        # Only remove N times an established integer polynomial (not arbitrary scalars).
        pieces = p.op == :add ? p.args : (p,)
        kept = [t for t in pieces if !multiple_of_N(t,N)]
        if length(kept)<length(pieces)
            return record!(trace,:periodicity,x,Term(op,(N,addition(kept))))
        end
    elseif op == :sign && isnum(a[1])
        return record!(trace,:integer_sign,x,term((-1)^Int(a[1].args[1])))
    elseif op == :tri && isnum(a[1])
        r=a[1].args[1]
        return record!(trace,:triangular_number,x,term(r*(r-1)/2))
    elseif op == :qchoose
        N,n,k=a
        if isnum(n) && isnum(k)
            nn,kk=n.args[1],k.args[1]
            nn<0 && throw(DomainError(nn,"q-binomial upper index must be nonnegative"))
            (kk<0 || kk>nn) && return record!(trace,:qbinomial_outside,x,term(0))
            (kk==0 || kk==nn) && return record!(trace,:qbinomial_boundary,x,term(1))
        end
    elseif op == :qfact && (isnum(a[2],0) || isnum(a[2],1))
        return record!(trace,:qfactorial_boundary,x,term(1))
    elseif op == :mul
        # Merge phases only for the same primitive root.
        phases=Dict{Term,Vector{Term}}(); others=Term[]
        for t in a
            if t.op==:phase; push!(get!(phases,t.args[1],Term[]),t.args[2]); else; push!(others,t); end
        end
        if any(length(v)>1 for v in values(phases))
            append!(others,[Term(:phase,(N,addition(v))) for (N,v) in phases])
            return record!(trace,:merge_phases,x,product(others))
        end
        # Ordinary equality substitution is safe in all factors; retain its constraint.
        for (j,d) in enumerate(a)
            d.op==:eq || continue
            u,v=d.args
            if u.op != :var && v.op==:var; u,v=v,u; end
            if u.op==:var && v.op==:var && key(u)<key(v); u,v=v,u; end
            if u.op==:var && !depends(v,u)
                out=[k==j ? d : substitute(t,u,v) for (k,t) in enumerate(a)]
                y=product(out)
                y==x || return record!(trace,:equality_propagation,x,y)
            end
        end
    elseif op == :sum
        i,lo,hi,b=a
        iszero(b) && return record!(trace,:zero_sum,x,term(0))
        if isnum(lo) && isnum(hi) && lo.args[1]>hi.args[1]
            return record!(trace,:empty_sum,x,term(0))
        end
        # Full residue domains are recognized structurally, including symbolic N.
        factors=b.op==:mul ? collect(b.args) : [b]
        retained=[t for t in factors if !(t.op==:range && t.args==(i,lo,hi))]
        if length(retained)<length(factors)
            return record!(trace,:bound_range,x,Term(:sum,(i,lo,hi,product(retained))))
        end
        for (j,d) in enumerate(factors)
            if d.op==:eq
                u,v=d.args
                u==i || ((u,v)=(v,u))
                if u==i && !depends(v,i) && integerexpr(v)
                    rest=product(t for (k,t) in enumerate(factors) if k!=j)
                    return record!(trace,:delta_sum,x,inrange(v,lo,hi)*substitute(rest,i,v))
                end
            elseif d.op==:moddelta
                N,p=d.args
                if isnum(lo,0) && hi==N-1
                    ab=affine(p,i)
                    if ab !== nothing && integerexpr(ab[2]) && (isnum(ab[1],1) || isnum(ab[1],-1))
                        v=Term(:mod,( -ab[2]/ab[1],N))
                        rest=product(t for (k,t) in enumerate(factors) if k!=j)
                        return record!(trace,:residue_delta_sum,x,substitute(rest,i,v))
                    end
                end
            end
        end
        dep=[t for t in factors if depends(t,i)]
        if length(dep)==1 && dep[1].op==:phase
            N,p=dep[1].args
            ab=affine(p,i)
            if isnum(lo,0) && hi==N-1 && ab !== nothing && integerexpr(ab[1]) && integerexpr(ab[2])
                rest=product(t for t in factors if !depends(t,i))
                return record!(trace,:fourier_sum,x,rest*N*Term(:phase,(N,ab[2]))*Term(:moddelta,(N,ab[1])))
            end
        end
        if !depends(b,i) && isnum(lo) && isnum(hi)
            return record!(trace,:constant_sum,x,(hi-lo+1)*b)
        elseif !depends(b,i) && context !== nothing && isnum(lo,0) && hi==context.N-1
            return record!(trace,:full_domain_count,x,context.N*b)
        end
    elseif op == :mod && all(isnum,a)
        return record!(trace,:representative,x,term(mod(a[1].args[1],a[2].args[1])))
    end
    x
end
integerexpr(x::Term) = x.op==:var || (isnum(x) && denominator(x.args[1])==1) ||
    (x.op in (:add,:mul,:tri) && all(integerexpr,x.args)) ||
    (x.op==:pow && x.args[2]>=0 && integerexpr(x.args[1])) ||
    (x.op==:mod && all(integerexpr,x.args))
multiple_of_N(t,N) = t==N || (t.op==:mul && N in t.args && all(integerexpr,t.args))

function sweep(x::Term,trace,context)
    if x.op==:mul
        out=Term[]
        for t in sort(collect(x.args);by=t->t.op in (:range,:eq,:moddelta) ? 0 : 1)
            y=sweep(t,trace,context)
            iszero(y) && return record!(trace,:zero_product,x,term(0))
            push!(out,y)
        end
        return localrule(product(out),trace,context)
    end
    y=rebuild(x.op,maparg(t->sweep(t,trace,context),x.args))
    localrule(y,trace,context)
end
"""Deterministic bounded simplification; unresolved sums are preserved.

`status=:fixed_point` is not a claim of complete normal form. No fallback
enumeration is performed. The trace records rule names and exact before/after IR.
"""
function normalize(x::Term; maxpasses::Int=64,maxnodes::Int=100000,context::Union{Nothing,TaftContext}=nothing)
    maxpasses>=0 && maxnodes>=1 || throw(ArgumentError("invalid normalization budget"))
    trace=NamedTuple{(:rule,:before,:after),Tuple{Symbol,Term,Term}}[]
    for pass in 1:maxpasses
        nodecount(x)>maxnodes && return Normalization(x,trace,:node_budget,pass-1)
        oldlength=length(trace)
        y=sweep(x,trace,context)
        if nodecount(y)>maxnodes
            resize!(trace,oldlength)
            return Normalization(x,trace,:node_budget,pass-1)
        end
        y==x && return Normalization(y,trace,:fixed_point,pass)
        x=y
    end
    Normalization(x,trace,:pass_budget,maxpasses)
end
normalize(x::Term,c::TaftContext;kwargs...) = normalize(x;context=c,kwargs...)

"""Specialize and evaluate with exact cyclotomic coefficients.

Each visited summand consumes `maxterms`; exhausted budgets throw rather than
return an approximate or partially evaluated scalar. Free variables must be integers.
"""
function evaluate_exact(x::Term,env::AbstractDict=Dict{Symbol,Int}();maxterms::Int=1000000)
    maxterms>=0 || throw(ArgumentError("maxterms must be nonnegative"))
    all(v isa Integer for v in values(env)) || throw(ArgumentError("free variables must specialize to integers"))
    budget=Ref(maxterms)
    exacteval(x,Dict{Symbol,Any}(k=>BigInt(v) for (k,v) in env),budget)
end
function evaluate_exact(x::Term,c::TaftContext,env::AbstractDict=Dict{Symbol,Int}();maxterms::Int=1000000)
    n=evaluate_exact(c.N,env;maxterms=0)
    asint(n)>=2 || throw(DomainError(n,"N must be ≥ 2"))
    evaluate_exact(x,env;maxterms=maxterms)
end
function asint(x)
    x isa Integer && return x
    x isa Rational && denominator(x)==1 && return numerator(x)
    throw(DomainError(x,"integer required"))
end
function exacteval(x,env,budget)
    op,a=x.op,x.args
    op==:num && return a[1]
    op==:var && return env[a[1]]
    op==:structure && return exacteval(lower(x),env,budget)
    if op==:sum
        i,lo,hi,b=a
        l,h=asint(exacteval(lo,env,budget)),asint(exacteval(hi,env,budget))
        result=BigInt(0)//BigInt(1)
        localenv=copy(env)
        for k in l:h
            budget[]-=1; budget[]>=0 || throw(ArgumentError("exact evaluation summand budget exhausted"))
            localenv[i.args[1]]=k
            result+=exacteval(b,localenv,budget)
        end
        return result
    end
    # Short-circuit products: zero domain indicators can guard undefined coefficients.
    if op==:mul
        result=BigInt(1)//BigInt(1)
        for t in sort(collect(a);by=t->t.op in (:range,:eq,:moddelta) ? 0 : 1)
            v=exacteval(t,env,budget); iszero(v) && return v; result*=v
        end
        return result
    end
    v=[t isa Term ? exacteval(t,env,budget) : t for t in a]
    op==:add && return sum(v)
    if op==:pow
        base=v[1]
        base isa Integer && v[2]<0 && (base=BigInt(base)//BigInt(1))
        return base^v[2]
    end
    op==:eq && return BigInt(v[1]==v[2])//BigInt(1)
    op==:range && return BigInt(v[2]<=v[1]<=v[3])//BigInt(1)
    op==:mod && return mod(asint(v[1]),asint(v[2]))
    op==:tri && return div(asint(v[1])*(asint(v[1])-1),2)
    op==:sign && return (-BigInt(1))^asint(v[1])
    if op in (:phase,:moddelta,:qchoose,:qfact)
        N=Int(asint(v[1])); N>=2 || throw(DomainError(N,"N must be ≥ 2"))
        op==:phase && return zeta(N)^asint(v[2])
        op==:moddelta && return BigInt(mod(asint(v[2]),N)==0)//BigInt(1)
        if op==:qchoose
            nn,kk=Int(asint(v[2])),Int(asint(v[3]))
            nn>=0 || throw(DomainError(nn,"q-binomial upper index must be nonnegative"))
            return kk<0 || kk>nn ? BigInt(0)//BigInt(1) : qbinomial(N,nn,kk)
        end
        op==:qfact && return qfactorial(N,Int(asint(v[2])))
    end
    error("unsupported expression: $op")
end

function Base.show(io::IO,x::Term)
    op,a=x.op,x.args
    if op==:num; print(io,denominator(a[1])==1 ? string(numerator(a[1])) : string(a[1]))
    elseif op==:var; print(io,a[1])
    elseif op in (:add,:mul); print(io,"("); join(io,a,op==:add ? " + " : " * "); print(io,")")
    elseif op==:sum; print(io,"Sum[",a[1],"=",a[2],":",a[3],"](",a[4],")")
    elseif op==:pow; print(io,"(",a[1],")^",a[2])
    else; print(io,op,"("); join(io,a,", "); print(io,")")
    end
end
end
