using Test
using AkaghefAlgebra
using AkaghefAlgebra.SymbolicTaft
import AkaghefAlgebra.TaftExact as TE

@testset "Symbolic Taft q-rule guards and budgets" begin
    c = TaftContext(3)
    k = integer(:k)
    @test qrewrite(qchoose(c,0,k), c; rule=:pascal) == qchoose(c,0,k)
    @test qrewrite(qchoose(c,1,k), c; rule=:pascal) != qchoose(c,1,k)

    n = integer(:n)
    unknown_n = finite_sum(:i,0,n) do i
        Term(:sign,(i,))*phase(c,Term(:tri,(i,))) * qchoose(c,n,i)
    end
    @test qrewrite(unknown_n,c;rule=:qbinomial) == unknown_n

    m = finite_sum(:i,0,2) do i
        qchoose(c,2,i) * phase(c,(2-i)*(1-i)) * qchoose(c,1,1-i)
    end
    @test qrewrite(m,c;rule=:vandermonde) == qchoose(c,3,1)

    b = integer(:b)
    bad_hockey = finite_sum(:i,b,-1) do i
        phase(c,-(b+1)*i+b^2) * qchoose(c,i,b)
    end
    @test qrewrite(bad_hockey,c;rule=:hockey) == bad_hockey

    huge = BigInt(typemax(Int)) + 10
    @test qrewrite(qchoose(c,huge,0),c) == SymbolicTaft.term(1)

    product_sum = finite_sum(:i,0,3) do i
        (integer(:x)+integer(:y)) * (integer(:u)+integer(:v))
    end
    @test qrewrite(product_sum,c;maxnodes=1) == product_sum
    @test_throws ArgumentError qrewrite(product_sum,c;maxnodes=0)
    @test_throws ArgumentError qrewrite(product_sum,c;maxexpansion=-1)

    @test zero_status(integer(:x)) == :unknown
    @test zero_status(integer(:x)-integer(:x)) == :zero
    @test zero_status(phase(c,integer(:x))) == :nonzero
    @test integer(:x) / 2 isa Term
    @test_throws ArgumentError integer(:x) / integer(:y)

    c3 = Term(:coeff,(TE.zeta(3),))
    c5 = Term(:coeff,(TE.zeta(5),))
    @test c3 != c5
    @test_throws ArgumentError c3 + c5
    @test (c3 * 2).op == :coeff
    @test evaluate_exact(c3,TaftContext(3),Dict(:N=>3)) == TE.zeta(3)
    q=AkaghefAlgebra.QSeries.qgenerator()
    coeff=SymbolicTaft.term((q^2-1)/(q-1))
    @test coeff.op==:coeff
    @test (coeff-SymbolicTaft.term(q+1))==SymbolicTaft.term(0)
    @test evaluate_exact(coeff,TaftContext(3))==1+TE.zeta(3)
    @test_throws ArgumentError evaluate_exact(coeff,TaftContext(3),Dict(:N=>5))
    shadow=Term(:sum,(integer(:N),SymbolicTaft.term(0),SymbolicTaft.term(1),coeff))
    @test evaluate_exact(shadow,TaftContext(3))==2*(1+TE.zeta(3))
    pole=SymbolicTaft.term(1/(q^2+q+1))
    @test_throws ArgumentError evaluate_exact(pole,TaftContext(3))
    @test_throws ArgumentError inv(integer(:x))
end
