using Test
import AkaghefAlgebra.TaftExact as TE

_qb(N,n,k) = k < 0 || k > n ? zero(TE.zeta(N)) : TE.qbinomial(N,n,k)

@testset "Exact cyclotomic arithmetic" begin
    for n in 2:6
        q = TE.zeta(n)
        @test q^n == one(q)
        @test q^(n-1) == one(q) / q
        @test (q + 1) - q == one(q)
        @test (2*q) / 2 == q
        @test (q^2) * (q^3) == q^5
        @test TE.qbinomial(n, 4, 0) == one(q)
        @test TE.qbinomial(n, 4, 4) == one(q)
        for k in 0:3
            @test TE.qbinomial(n, 4, k) == _qb(n,3,k-1) + q^k * _qb(n,3,k)
        end
        @test TE.qfactorial(n, 0) == one(q)
        @test TE.qfactorial(n, 1) == one(q)
    end
    @test TE.zeta(2) == -1
    @test TE.zeta(1) == 1
    @test_throws ArgumentError TE.zeta(0)
    @test TE.qbinomial(3, 2, 3)==0
    @test_throws ArgumentError TE.qbinomial(3, -1, 0)
    @test_throws ArgumentError TE.TaftElement(3, Dict((0,0)=>1.0))
    @test_throws ArgumentError TE.TaftElement(3, Dict((0,0)=>TE.zeta(4)))
end

function _tensor_mul(A::TE.TaftTensor{N}, B::TE.TaftTensor{N}) where {N}
    n = Int(N); d = Dict{Tuple{Tuple{Int,Int},Tuple{Int,Int}},Any}()
    for ((a,b),x) in A.terms, ((c,e),y) in B.terms
        p = TE.monomial(n,a[1],a[2],x) * TE.monomial(n,c[1],c[2],y)
        q = TE.monomial(n,b[1],b[2],1) * TE.monomial(n,e[1],e[2],1)
        for ((lr,ls),lx) in p.terms, ((rr,rs),rx) in q.terms
            key = ((lr,ls),(rr,rs)); v = lx*rx
            d[key] = get(d,key,zero(v)) + v
        end
    end
    TE.TaftTensor{N,Any}(Dict(k=>v for (k,v) in d if !iszero(v)))
end

function _tensor_map_mul(A::TE.TaftTensor{N}) where {N}
    n = Int(N); out = TE.zero_element(n)
    for ((left,right),c) in A.terms
        out = out + c * (TE.monomial(n,left[1],left[2],1) *
                         TE.monomial(n,right[1],right[2],1))
    end
    out
end

function _tensor3_left(x::TE.TaftElement{N}) where {N}
    n = Int(N); d = Dict{Tuple{Tuple{Int,Int},Tuple{Int,Int},Tuple{Int,Int}},Any}()
    for ((a,b),c) in TE.coproduct(x).terms, ((u,v),w) in TE.coproduct(TE.monomial(n,a[1],a[2])).terms
        key = (u,v,b); d[key] = get(d,key,zero(c*1)) + c*w
    end
    Dict(k=>v for (k,v) in d if !iszero(v))
end

function _tensor3_right(x::TE.TaftElement{N}) where {N}
    n = Int(N); d = Dict{Tuple{Tuple{Int,Int},Tuple{Int,Int},Tuple{Int,Int}},Any}()
    for ((a,b),c) in TE.coproduct(x).terms, ((u,v),w) in TE.coproduct(TE.monomial(n,b[1],b[2])).terms
        key = (a,u,v); d[key] = get(d,key,zero(c*1)) + c*w
    end
    Dict(k=>v for (k,v) in d if !iszero(v))
end

function _apply_counit_left(x::TE.TaftElement{N}) where {N}
    n = Int(N); out = TE.zero_element(n)
    for ((a,b),c) in TE.coproduct(x).terms
        out = out + c * TE.counit(TE.monomial(n,a[1],a[2])) * TE.monomial(n,b[1],b[2])
    end
    out
end

function _apply_counit_right(x::TE.TaftElement{N}) where {N}
    n = Int(N); out = TE.zero_element(n)
    for ((a,b),c) in TE.coproduct(x).terms
        out = out + c * TE.counit(TE.monomial(n,b[1],b[2])) * TE.monomial(n,a[1],a[2])
    end
    out
end

function _b_value(n, x)
    q = TE.zeta(n)
    length(x.terms) == 1 || return zero(q)
    r = first(keys(x.terms))[1]
    r == 1 ? one(q) : zero(q)
end

function _b_power(n, t, x)
    t == 0 && return TE.counit(x)
    out = zero(TE.zeta(n))
    for ((a,b),c) in TE.coproduct(x).terms
        out += c * _b_power(n,t-1,TE.monomial(n,a[1],a[2])) *
               _b_value(n,TE.monomial(n,b[1],b[2]))
    end
    out
end

@testset "Sparse Taft Hopf operations" begin
    for n in 2:6
        q = TE.zeta(n); e = TE.E(n); k = TE.K(n)
        @test typeof(2*e) == typeof(e)
        @test 2*e + e == 3*e
        @test (1//2)*e + (1//2)*e == e
        @test iszero(0*e)
        @test k*e == q*e*k
        @test e^n == TE.zero_element(n)
        @test k^n == TE.one_element(n)
        Δe = TE.coproduct(e)
        @test Δe.terms == Dict(((0,0),(1,0))=>one(q), ((1,0),(0,1))=>one(q))
        Δk = TE.coproduct(k)
        @test Δk.terms == Dict(((0,1),(0,1))=>one(q))

        # Independent generator expansion: Δ(E)^r is multiplied in the
        # tensor product and compared with the closed q-binomial expansion.
        lhs = TE.TaftTensor{n,typeof(q)}(Dict(((0,0),(0,0))=>one(q)))
        for r in 1:n-1
            lhs = _tensor_mul(lhs, Δe)
            rhs = TE.coproduct(e^r)
            @test lhs.terms == rhs.terms
        end

        # Check the Hopf identities on every PBW basis element in this bound.
        for r in 0:n-1, s in 0:n-1
            x = TE.monomial(n,r,s)
            Δx = TE.coproduct(x)
            @test TE.counit(x) == (r == 0 ? one(q) : zero(q))
            @test _tensor3_left(x) == _tensor3_right(x)
            @test _apply_counit_left(x) == x
            @test _apply_counit_right(x) == x
            # The antipode identity is evaluated directly as m(S⊗id)Δ(x).
            sminus = TE.zero_element(n)
            for ((a,b),c) in Δx.terms
                sminus = sminus + c * (TE.antipode(TE.monomial(n,a[1],a[2])) *
                                       TE.monomial(n,b[1],b[2]))
            end
            @test sminus == TE.counit(x) * TE.one_element(n)
            splus = TE.zero_element(n)
            for ((a,b),c) in Δx.terms
                splus = splus + c * (TE.monomial(n,a[1],a[2]) *
                                     TE.antipode(TE.monomial(n,b[1],b[2])))
            end
            @test splus == TE.counit(x) * TE.one_element(n)
        end
    end
end

@testset "Taft dual monomial and coordinate pairing" begin
    for n in 2:6
        q = TE.zeta(n)
        for t in 0:n-1, r in 0:n-1, s in 0:n-1, sp in 0:n-1
            x = TE.monomial(n,r,sp)
            # Build b^t by independent convolution recursion from Δ and the
            # generator functional b(E*K^s)=1, then compare with the formula.
            @test _b_power(n,t,x) == TE.pairing(TE.dual_monomial(n,0,t),x)
            for l in 0:n-1
                expected = r == t ? TE.qfactorial(n,t) * q^(l*sp) : zero(q)
                @test TE.pairing(TE.dual_monomial(n,l,t), x) == expected
            end
            f = TE.dual_basis(n,t,s)
            @test TE.pairing(f, x) == (t == r && s == sp ? 1 : 0)
            fourier = zero(q)
            for (lt, coefficient) in TE.dual_basis_expansion(f)
                fourier += coefficient * TE.pairing(TE.dual_monomial(n,lt[1],lt[2]), x)
            end
            @test fourier == (t == r && s == sp ? 1 : 0)
        end
    end
end
