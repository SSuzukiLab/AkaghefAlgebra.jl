"""Topology-independent closed Taft networks in the E^r K^s basis, style L."""
module TaftNetworks
import ..SymbolicTaft as ST
import ..TaftExact as TE
const T=ST.Term
export network, reduce_network, reduce_with_qrules, close_sums, galois

function network(ports, signs, weights)
    n=length(ports)
    length(weights)==n || throw(ArgumentError("one integer weight per edge required"))
    all(w->w isa Integer,weights) || throw(ArgumentError("weights must be integers"))
    sort(vcat(ports...))==collect(1:4length(signs)) || throw(ArgumentError("ports must cover each vertex slot exactly once"))
    all(s->s in (-1,1),signs) || throw(ArgumentError("vertex signs must be ±1"))
    N=ST.integer(:N); ctx=ST.TaftContext(N)
    rs=[ST.integer(Symbol("r",j)) for j in 1:n]
    ss=[ST.integer(Symbol("s",j)) for j in 1:n]
    edgeof=Dict(p=>j for (j,ends) in enumerate(ports) for p in ends)
    fs=T[]
    for (j,sgn) in enumerate(signs)
        push!(fs,vertex(ctx,sgn,[(rs[edgeof[k]],ss[edgeof[k]]) for k in 4(j-1).+(1:4)]))
    end
    append!(fs,[ST.phase(ctx,-weights[j]*rs[j]) for j in 1:n])
    ctx,vcat(rs,ss),ST.product(fs)
end
function vertex(c,sign,p)
    (r1,s1),(r2,s2),(r3,s3),(r4,s4)=p
    m=r4-r2
    common=ST.eqdelta(r1,r3+m)*ST.moddelta(c,s2-s4)*
           ST.inrange(r2,0,r4)*ST.qchoose(c,r4,r2)
    if sign>0
        return common*ST.moddelta(c,s1-s3-s4-r2)*ST.phase(c,s3*m)
    end
    # S_L(E^m K^t)=(-1)^m q^(-tm-m(m+1)/2) E^m K^(-t-m).
    common*ST.moddelta(c,s1-s3+s4+r4)*T(:sign,(m,))*
        ST.phase(c,(s3-s4-r2)*m-T(:tri,(m+1,)))
end

# Expand only integer-index polynomials, never scalar finite sums/q-binomials.
function indexexpand(x::T)
    y=ST.rebuild(x.op,ST.maparg(indexexpand,x.args))
    if y.op==:mul
        acc=[ST.term(1)]
        for a in y.args
            parts=a.op==:add ? a.args : (a,)
            length(acc)*length(parts)<=100_000 || throw(ArgumentError("index expansion size budget exceeded"))
            acc=[u*v for u in acc for v in parts]
        end
        return ST.addition(acc)
    end
    y
end

# Removing mod representatives is safe in integer polynomial phase exponents,
# but NOT through triangular numbers (which have denominator two).
function periodic_lift(x::T,N)
    x.op==:mod && x.args[2]==N && return periodic_lift(x.args[1],N)
    x.op in (:add,:mul) && return ST.rebuild(x.op,Tuple(periodic_lift(a,N) for a in x.args))
    x.op==:pow && x.args[2]>=0 && return periodic_lift(x.args[1],N)^x.args[2]
    x
end
function tidy(x::T,ctx,vars)
    y=ST.rebuild(x.op,ST.maparg(t->tidy(t,ctx,vars),x.args))
    if y.op in (:phase,:moddelta)
        N,p=y.args
        y=T(y.op,(N,indexexpand(periodic_lift(p,N))))
    elseif y.op==:eq
        y=T(:eq,(indexexpand(y.args[1]-y.args[2]),ST.term(0)))
    elseif y.op==:range
        a,lo,hi=indexexpand.(y.args)
        lo==hi && return ST.eqdelta(a,lo)
        y=T(:range,(a,lo,hi))
        a in vars && lo==ST.term(0) && hi==ctx.N-1 && return ST.term(1)
        a.op==:mod && a.args[2]==ctx.N && lo==ST.term(0) && hi==ctx.N-1 && return ST.term(1)
    elseif y.op==:qchoose
        N,n,k=y.args; n=indexexpand(n); k=indexexpand(k)
        # Remaining bound variables are certified in 0:N-1.
        nonnegative=n in vars || ST.isnum(n) && n.args[1]>=0 || n==ctx.N-1
        nonnegative && (ST.isnum(k,0) || k==n) && return ST.term(1)
        y=T(:qchoose,(N,n,k))
    elseif y.op in (:sign,:tri)
        y=T(y.op,(indexexpand(y.args[1]),))
    end
    y
end

# Finite q-binomial theorem, with 0<=b<N and a summed over 0:N-1:
# sum_a (-1)^(b-a) [b;a]_q q^(-tri(b-a+1)-a-ab+a^2) [a<=b]
# = q^(-b) prod_{j=0}^{b-1}(1-q^(1-b+j)) = delta(b,0).
function eliminate_qbinomial(body,vars,ctx)
    fs=factors(body)
    for a in vars, (j,qb) in enumerate(fs)
        qb.op==:qchoose || continue
        N,b,k=qb.args
        N==ctx.N && k==a && b in vars && b!=a || continue
        sg=T(:sign,(indexexpand(b-a),))
        rg=ST.inrange(a,0,b)
        si=findfirst(==(sg),fs); ri=findfirst(==(rg),fs)
        (si===nothing || ri===nothing) && continue
        pi=findfirst(t->t.op==:phase && t.args[1]==N,fs)
        pi===nothing && continue
        expected=indexexpand(-T(:tri,(indexexpand(b-a+1),))-a-a*b+a*a)
        restphase=indexexpand(fs[pi].args[2]-expected)
        ST.depends(restphase,a) && continue
        keep=[t for (h,t) in enumerate(fs) if h ∉ (j,si,ri,pi)]
        any(ST.depends(t,a) for t in keep) && continue
        push!(keep,ST.phase(ctx,restphase),ST.eqdelta(b,0))
        return (body=ST.product(keep),variable=a,value=b,rule=:qbinomial_zero_product)
    end
    nothing
end

function eliminate_boundary(body,vars,ctx)
    for d in factors(body)
        d.op==:range || continue
        a,lo,hi=d.args
        for gap in (indexexpand(a-lo),indexexpand(hi-a)), i in vars
            gap == -i || continue
            return (body=ST.substitute(body,i,0),variable=i,value=ST.term(0),rule=:nonnegative_boundary)
        end
    end
    nothing
end
factors(x)=x.op==:mul ? collect(x.args) : [x]

"""Try legal unit-coefficient delta elimination in a reordered finite sum.
The residual range condition is retained for an ordinary equality.
"""
function eliminate_constraint(body,vars,ctx,kind)
    fs=factors(body)
    order=kind==:moddelta ? sort(vars;by=i->startswith(string(i.args[1]),"s") ? 0 : 1) : vars
    for i in order, (j,d) in enumerate(fs)
        d.op==kind || continue
        p=kind==:eq ? indexexpand(d.args[1]-d.args[2]) : d.args[2]
        ab=ST.affine(p,i)
        ab===nothing && continue
        a,b=ab
        (ST.isnum(a,1)||ST.isnum(a,-1)) || continue
        v=indexexpand(-b/a)
        ST.integerexpr(v) || continue
        kind==:moddelta && (v=T(:mod,(v,ctx.N)))
        rest=ST.product(t for (k,t) in enumerate(fs) if k!=j)
        rest=ST.substitute(rest,i,v)
        kind==:eq && (rest*=ST.inrange(v,0,ctx.N-1))
        return (body=rest,variable=i,value=v,rule=kind==:eq ? :integer_delta : :modular_delta)
    end
    nothing
end

function reduce_network(body,vars,ctx;maxsteps=100,maxnodes=100_000)
    vars=copy(vars); trace=Any[]
    for step in 1:maxsteps
        ST.nodecount(body)<=maxnodes || throw(ArgumentError("network expression size budget exceeded"))
        body=ST.normalize(tidy(body,ctx,vars),ctx).expression
        changed=false
        for kind in (:eq,:moddelta)
            hit=eliminate_constraint(body,vars,ctx,kind)
            hit===nothing && continue
            push!(trace,Dict("rule"=>string(hit.rule),"index"=>string(hit.variable),"value"=>string(hit.value)))
            body=hit.body; filter!(!=(hit.variable),vars)
            changed=true; break
        end
        changed && continue
        for eliminate in (eliminate_boundary,eliminate_qbinomial)
            hit=eliminate(body,vars,ctx)
            hit===nothing && continue
            push!(trace,Dict("rule"=>string(hit.rule),"index"=>string(hit.variable),"value"=>string(hit.value)))
            body=hit.body;filter!(!=(hit.variable),vars);changed=true;break
        end
        changed && continue
        for i in vars
            attempt=ST.normalize(T(:sum,(i,ST.term(0),ctx.N-1,body)),ctx)
            attempt.expression.op==:sum && continue
            push!(trace,Dict("rule"=>"library_sum","index"=>string(i),"rules"=>[string(t.rule) for t in attempt.trace]))
            body=attempt.expression; filter!(!=(i),vars); changed=true; break
        end
        changed || return body,vars,trace
    end
    error("symbolic reduction step budget exhausted")
end

function close_sums(body,vars,ctx)
    for i in reverse(vars)
        body=T(:sum,(i,ST.term(0),ctx.N-1,body))
    end
    body
end
function galois(x,N,M)
    x isa Union{Integer,Rational} && return x
    sum(a*TE.zeta(N)^(M*(j-1)) for (j,a) in enumerate(TE.coefficients(x)))
end

hasbound(x::T,i)=(x.op==:sum && x.args[1]==i) || any(a->hasbound(a,i),x.args)
hasbound(x::Tuple,i)=any(a->hasbound(a,i),x)
hasbound(x,i)=false

function reduce_with_qrules(body,vars,ctx;baseline=false)
    body,vars,trace=reduce_network(body,vars,ctx)
    baseline && return body,vars,trace
    # The q-rule hook is deliberately bounded; unrecognized sums remain sums.
    for pass in 1:16
        bounds=Dict(i=>(ST.term(0),ctx.N-1) for i in vars)
        next=ST.qrewrite(body,ctx;bounds=bounds)
        next=next isa ST.Term ? next : next.expression
        changed=next!=body
        body=next
        for i in copy(vars)
            trial=T(:sum,(i,ST.term(0),ctx.N-1,body))
            result=ST.qrewrite(trial,ctx;bounds=bounds)
            result=result isa ST.Term ? result : result.expression
            if !hasbound(result,i) && !ST.depends(result,i)
                push!(trace,Dict("rule"=>"qrewrite_sum","index"=>string(i)))
                body=result;filter!(!=(i),vars);changed=true;break
            end
        end
        body,vars,extra=reduce_network(body,vars,ctx);append!(trace,extra)
        !changed && return body,vars,trace
    end
    return body,vars,trace
end


end
