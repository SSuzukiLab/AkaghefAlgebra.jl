using Test
using AkaghefAlgebra

@testset "R-convention Oq(SL2) dual PBW basis" begin
    for q in (2//1, 3//1)
        @testset "old formula is rejected at q=$q" begin
            leaks = qsl2_old_formula_leaks(q)
            @test leaks.E_pairing == q
            @test leaks.E_pairing != 1
            @test leaks.EH_to_E_hbar_coefficient == q
            @test leaks.EH_to_E_hbar_coefficient != 0
        end

        @testset "corrected generators at q=$q" begin
            generators = qsl2_corrected_generator_pairings(q)
            @test generators.E_on_E == 1
            @test generators.EH_on_E == 0
            @test generators.EH_on_EH == 1
            @test generators.F_on_F == 1
            @test generators.H_on_H == 1
        end

        @testset "root normalization through exponent 4 at q=$q" begin
            for i in 0:4, j in 0:4
                raw = qsl2_root_pairing_iterative(q, i, j)
                coefficient = qsl2_root_dual_coefficient(q, i, j)
                @test coefficient * raw == 1
            end
        end

        @testset "full Kronecker duality through exponent 2 at q=$q" begin
            for i in 0:2, j in 0:2, k in 0:2,
                r in 0:2, s in 0:2, t in 0:2
                expected = (i,j,k) == (r,s,t) ? 1 : 0
                @test qsl2_dual_pairing(q, i,j,k, r,s,t) == expected
            end
        end
    end
end
