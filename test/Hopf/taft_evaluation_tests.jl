using Test
using AkaghefAlgebra
using AkaghefAlgebra.SymbolicTaft
using AkaghefAlgebra.TaftEvaluation
import AkaghefAlgebra.TaftExact as TE

@testset "compiled exact Taft evaluator" begin
    N = integer(:N)
    c = TaftContext(N)
    a = integer(:a)
    r = integer(:r)

    # This exercises nested binders, all indicator nodes, integer arithmetic,
    # phases, q-binomial variants, q-factorials, q-products, powers and signs.
    expr = finite_sum(:r, 0, N-1) do rr
        finite_sum(:s, 0, N-1) do ss
            inrange(rr, 0, N-1) * eqdelta(ss, Term(:mod, (rr+1, N))) *
            moddelta(c, ss-rr-1) * phase(c, rr+2ss-a) *
            qchoose(c, rr, ss) * qchoose_inv(c, rr, ss) *
            qfact(c, rr) * qproduct(c, a+rr, 0) *
            (rr-ss)^2 * Term(:sign, (rr-ss,)) *
            phase(c, Term(:tri, (rr,)))
        end
    end

    for n in 2:6, av in -2:2
        env = Dict(:N => n, :a => av)
        @test evaluate_fast(expr, c, env; maxterms=100_000) ==
              evaluate_exact(expr, c, env; maxterms=100_000)
        p = compile_fast(expr)
        @test evaluate_fast(p, n, env; maxterms=100_000) ==
              evaluate_exact(expr, c, env; maxterms=100_000)
    end

    # Structure constants are lowered before compilation and retain exactness.
    st = finite_sum(:s, 0, N-1) do s
        structure(c, :mul, (r, s), (1, 0), (r+1, s))
    end
    for n in 2:5, rv in 0:n-1, sv in 0:n-1
        env = Dict(:N => n, :r => rv)
        @test evaluate_fast(st, c, env; maxterms=20) ==
              evaluate_exact(st, c, env; maxterms=20)
    end

    @test_throws ArgumentError evaluate_fast(expr, c, Dict(:N => 3, :a => 0); maxterms=1)
    @test_throws DomainError evaluate_fast(expr, c, Dict(:N => 1, :a => 0))
    @test_throws ArgumentError evaluate_fast(expr, c, Dict(:N => 3, :a => 0.5))

    phase_expr = phase(c, a)
    @test evaluate_fast(compile_fast(phase_expr), 3, Dict(:a => 1)) == TE.zeta(3)
    @test_throws ArgumentError evaluate_fast(compile_fast(phase_expr), 3,
                                              Dict(:N => 4, :a => 1))
    huge_upper = qchoose(c, integer(:u), 0)
    @test_throws DomainError evaluate_fast(huge_upper, c,
                                            Dict(:N => 3, :u => BigInt(typemax(Int)) + 1))
    @test_throws DomainError evaluate_fast(compile_fast(phase_expr),
                                            BigInt(typemax(Int)) + 1, Dict(:a => 0))
    fixed_context = TaftContext(9)
    @test_throws ArgumentError evaluate_fast(phase(fixed_context, a), fixed_context,
                                              Dict(:N => 3, :a => 0))

    coeff = Term(:coeff, (TE.zeta(3),))
    @test evaluate_fast(coeff, TaftContext(3)) == TE.zeta(3)
end
