using Test
import AbstractAlgebra
using AkaghefAlgebra

const QS = AkaghefAlgebra.QSeries
const TE = AkaghefAlgebra.TaftExact

@testset "QSeries coefficient backend" begin
    field, q = QS.qring()
    @test field === QS.qfield()
    @test q === QS.qgenerator()
    @test q isa QS.QQq
    @test !(q isa Number)

    @test QS.qnumber(q, 0) == zero(q)
    @test QS.qnumber(q, 4) == one(q) + q + q^2 + q^3
    @test_throws ArgumentError QS.qnumber(q, -1)

    @test QS.qpochhammer(q, q, 0) == one(q)
    @test QS.qpochhammer(q, q, 3) == (one(q) - q) * (one(q) - q^2) * (one(q) - q^3)
    @test QS.qpochhammer(q, 3) == QS.qpochhammer(q, q, 3)
    @test_throws ArgumentError QS.qpochhammer(q, q, -1)

    @test QS.qbinomial(q, 0, 0) == one(q)
    @test QS.qbinomial(q, 4, 0) == one(q)
    @test QS.qbinomial(q, 4, 4) == one(q)
    @test QS.qbinomial(q, 4, 2) == one(q) + q + 2q^2 + q^3 + q^4
    @test QS.qbinomial(q, 4, -1) == zero(q)
    @test QS.qbinomial(q, 4, 5) == zero(q)
    @test_throws ArgumentError QS.qbinomial(q, -1, 0)

    @test QS.qfactorial(q, 0) == one(q)
    @test QS.qfactorial(q, 3) == QS.qnumber(q, 1) * QS.qnumber(q, 2) * QS.qnumber(q, 3)
    @test_throws ArgumentError QS.qfactorial(q, -1)

    ctx = QS.qcontext(q; max_n=3)
    @test QS.qbinomial(ctx, 4, 2) == QS.qbinomial(q, 4, 2)
    @test all(k <= 3 for k in keys(ctx.qbinomial_rows))
    QS.clear_cache!(ctx)
    @test Set(keys(ctx.qbinomial_rows)) == Set([0])

    @test QS.specialize_q(q, 3) == TE.zeta(3)
    @test QS.specialize_q((q^3 - 1) / (q - 1), 3) == zero(TE.zeta(3))
    @test QS.specialize_q(one(q) / (q + one(q)), 3) isa TE.Cyclotomic{3}
    @test_throws ArgumentError QS.specialize_q(one(q) / (q^2 + q + one(q)), 3)

    for n in 3:9
        root = TE.zeta(n)
        @test QS.qnumber(root, n) == zero(root)
        @test !iszero(QS.qnumber(root, n - 1))
        @test QS.qfactorial(root, n) == zero(root)
        @test !iszero(QS.qfactorial(root, n - 1))
        for k in 1:(n - 1)
            @test QS.qbinomial(root, n, k) == zero(root)
            @test !iszero(QS.qbinomial(root, n - 1, k))
        end
    end

    other_ring, other_q_poly = AbstractAlgebra.polynomial_ring(AbstractAlgebra.QQ, "q"; cached=false)
    other_q = AbstractAlgebra.fraction_field(other_ring)(other_q_poly)
    @test_throws ArgumentError QS.specialize_q(other_q, 3)
    @test QS.specialize_q(TE.zeta(3), 3) == TE.zeta(3)
    @test_throws ArgumentError QS.specialize_q(TE.zeta(3), 5)
end
