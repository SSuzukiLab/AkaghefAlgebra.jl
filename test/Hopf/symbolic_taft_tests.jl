using AkaghefAlgebra.SymbolicTaft
import AkaghefAlgebra.SymbolicTaft: normalize
import AkaghefAlgebra.TaftExact as TE

@testset "Symbolic Taft: binders, constraints, normalization" begin
    N=integer(:N); a=integer(:a); b=integer(:b); c=TaftContext(N)
    F=finite_sum(:s,0,N-1) do s
        phase(c,a*s+b)
    end
    result=normalize(F,c)
    @test result.status==:fixed_point
    @test any(t.rule==:fourier_sum for t in result.trace)
    @test result.expression == normalize(N*phase(c,b)*moddelta(c,a),c).expression
    @test normalize(result.expression,c).expression==result.expression
    for n in 2:6, av in -3:3, bv in -1:1
        env=Dict(:N=>n,:a=>av,:b=>bv)
        # Independent Fourier oracle: build actual powers, bypass IR evaluator.
        expected=sum(TE.zeta(n)^(av*s+bv) for s in 0:n-1)
        @test evaluate_exact(F,c,env)==expected
        @test evaluate_exact(result.expression,c,env)==expected
    end
    @test normalize(phase(c,a+N*b),c).expression==phase(c,a)
    @test normalize(phase(c,a)*phase(c,-a),c).expression==SymbolicTaft.term(1)

    G=finite_sum(:t,0,N-1) do t
        phase(c,a*t+b)
    end
    @test F==G
    @test hash(F)==hash(G)
    @test freevars(F)==Set([:N,:a,:b])
    inner=finite_sum(:i,0,2) do i
        finite_sum(:i,0,2) do j
            eqdelta(i,j)
        end
    end
    @test evaluate_exact(inner)==3
    @test evaluate_exact(normalize(inner).expression)==3
    # Deliberately attempt capture by substituting a with the actual binder.
    binder=F.args[1]
    captured=substitute(F,a,binder)
    @test binder.args[1] in freevars(captured)
    @test captured.args[1] != binder
    @test evaluate_exact(captured,c,Dict(:N=>3,:b=>0,binder.args[1]=>1))==0

    delta=finite_sum(:s,0,N-1) do s
        eqdelta(s,a)*phase(c,s*b)
    end
    residue=finite_sum(:s,0,N-1) do s
        moddelta(c,s-a)*s
    end
    nonunit=finite_sum(:s,0,N-1) do s
        moddelta(c,2*s-a)*s
    end
    for n in 2:6, av in -2:n+2
        env=Dict(:N=>n,:a=>av,:b=>2)
        for f in (delta,residue,nonunit)
            @test evaluate_exact(f,c,env)==evaluate_exact(normalize(f,c).expression,c,env)
        end
        @test evaluate_exact(normalize(residue,c).expression,c,env)==mod(av,n)
    end
    @test normalize(nonunit,c).expression.op==:sum
    # Congruence must not substitute into a nonperiodic factor outside a sum.
    bad=moddelta(c,a)*a
    @test evaluate_exact(normalize(bad,c).expression,c,Dict(:N=>3,:a=>3))==3
    quadratic=finite_sum(:s,0,N-1) do s
        phase(c,s*s)
    end
    @test normalize(quadratic,c).expression.op==:sum
    restricted=finite_sum(:s,0,N-2) do s
        phase(c,a*s)
    end
    @test normalize(restricted,c).expression.op==:sum
    halfdelta=finite_sum(:s,0,2) do s
        eqdelta(s,1//2)
    end
    @test evaluate_exact(normalize(halfdelta).expression)==0
    @test normalize(halfdelta).expression.op==:sum
    @test normalize(finite_sum(s->s,:s,3,2)).expression==SymbolicTaft.term(0)
    @test normalize(finite_sum(s->SymbolicTaft.term(7),:s,0,N-1),c).expression==7*N
    @test normalize(finite_sum(s->SymbolicTaft.term(7),:s,0,N-1)).expression.op==:sum
    @test normalize(F;maxpasses=0).status==:pass_budget
    @test normalize(F;maxnodes=1).status==:node_budget
    @test_throws ArgumentError evaluate_exact(F,c,Dict(:N=>3,:a=>1,:b=>0);maxterms=2)
    @test_throws DomainError evaluate_exact(result.expression,c,Dict(:N=>1,:a=>1,:b=>0))
    @test_throws ArgumentError TaftContext(1)
    @test_throws ArgumentError evaluate_exact(a,Dict(:a=>0.5))
    @test evaluate_exact(qchoose(c,3,4),c,Dict(:N=>4))==0
    @test evaluate_exact(normalize(qchoose(c,3,4),c).expression,c,Dict(:N=>4))==0
    @test evaluate_exact(a*a,Dict(:a=>2^40))==BigInt(2)^80
    for rr in 0:1
        dual=structure(c,:dual_basis,(rr,0),(0,rr))
        @test evaluate_exact(normalize(dual,c).expression,c,Dict(:N=>3))==1//3
    end

    # Multiple equalities must terminate deterministically.
    chain=eqdelta(a,b)*eqdelta(b,N)*phase(c,a+b)
    @test normalize(chain,c).status==:fixed_point
    @test evaluate_exact(chain,c,Dict(:N=>3,:a=>3,:b=>3))==
          evaluate_exact(normalize(chain,c).expression,c,Dict(:N=>3,:a=>3,:b=>3))
end

@testset "Symbolic structure constants remain lazy and compose" begin
    N=integer(:N); r=integer(:r); s=integer(:s); c=TaftContext(N)
    st=structure(c,:mul,(r,s),(1,0),(r+1,s))
    @test st.op==:structure
    @test :r in freevars(st)
    @test :r ∉ freevars(substitute(st,r,0))
    @test nodecount(st)<25
    for n in 2:6, rv in 0:n-1, sv in 0:n-1
        env=Dict(:N=>n,:r=>rv,:s=>sv)
        expected=rv+1<n ? TE.zeta(n)^sv : 0
        @test evaluate_exact(st,c,env)==expected
        @test evaluate_exact(normalize(st,c).expression,c,env)==expected
    end
    # Actual contraction: sum over output E degree of a product constant.
    contraction=finite_sum(:v,0,N-1) do v
        structure(c,:mul,(r,s),(1,0),(v,s))
    end
    for n in 2:6, rv in 0:n-1
        env=Dict(:N=>n,:r=>rv,:s=>1)
        @test evaluate_exact(contraction,c,env)==evaluate_exact(normalize(contraction,c).expression,c,env)
    end
    @test_throws ArgumentError structure(c,:unknown,(0,0))
end

@testset "Symbolic Hopf coefficients specialize to sparse backend" begin
    for n in 2:4
        c=TaftContext(n)
        for r in 0:n-1, s in 0:n-1
            x=TE.monomial(n,r,s)
            delta=TE.tensor_terms(TE.coproduct(x))
            sx=TE.coefficients(TE.antipode(x))
            @test evaluate_exact(structure(c,:counit,(r,s)),c)==TE.counit(x)
            @test evaluate_exact(structure(c,:unit,(r,s)),c)==(r==s==0 ? 1 : 0)
            for j in 0:n-1
                left=(j,s); right=(mod(r-j,n),mod(s+j,n))
                expr=structure(c,:comul,(r,s),left,right)
                @test evaluate_exact(normalize(expr,c).expression,c)==get(delta,(left,right),0)
            end
            for t in 0:n-1, u in 0:n-1
                ex=structure(c,:antipode,(r,s),(t,u))
                @test evaluate_exact(normalize(ex,c).expression,c)==get(sx,(t,u),0)
                px=structure(c,:pair,(t,u),(r,s))
                @test evaluate_exact(px,c)==TE.pairing(TE.dual_monomial(n,t,u),x)
            end
            # Explicit Fourier expansion, not the backend's coordinate lookup.
            for sp in 0:n-1
                dualpair=finite_sum(:l,0,n-1) do l
                    structure(c,:dual_basis,(r,s),(l,r))*structure(c,:pair,(l,r),(r,sp))
                end
                @test evaluate_exact(dualpair,c)==(s==sp ? 1 : 0)
                @test evaluate_exact(normalize(dualpair,c).expression,c)==(s==sp ? 1 : 0)
            end
        end
        invalid=structure(c,:comul,(-1,0),(0,0),(0,0))
        @test evaluate_exact(normalize(invalid,c).expression,c)==0
    end
end
