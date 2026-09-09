using Test

import AkaghefAlgebra.TaftExact as TE
import AkaghefAlgebra.CyclotomicValues as CV

@testset "Public rooted cyclotomic values" begin
    q = TE.zeta(9)
    base = q^(-64)
    r1 = CV.RootSpec(9, 1)
    r2 = CV.RootSpec(9, 2)

    @test r1.N == 9 && r1.M == 1
    @test CV.RootSpec(9, 11).M == 2
    @test CV.RootSpec(9, -7).M == 2
    @test_throws ArgumentError CV.RootSpec(1)
    @test_throws ArgumentError CV.RootSpec(9, 3)
    @test_throws ArgumentError CV.RootSpec{9,11}(9, 11)

    v1 = CV.normalize_root(base, r1)
    v2 = CV.normalize_root(base, r2)
    @test CV.canonical_value(v1) == -q^2 - q^5
    @test CV.canonical_value(v2) == -q - q^4
    @test CV.root_spec(v1) == r1
    @test CV.root_spec(v2) == r2
    @test CV.specified_root(v2) == q^2
    @test CV.normalize_root(v2, r2) === v2
    @test_throws ArgumentError CV.normalize_root(v2, r1)

    v2relabelled = CV.with_root(v2, r1)
    @test CV.canonical_value(v2relabelled) == CV.canonical_value(v2)
    @test CV.root_spec(v2relabelled) == r1
    @test hash(v2relabelled) == hash(v2)

    @test v1 == CV.normalize_root(base, CV.RootSpec(9, 10))
    @test v1 != v2
    @test CV.Cyclotomic(r1, 1) == CV.Cyclotomic(r2, 1)
    @test CV.Cyclotomic(r1, 1) == 1
    @test !isequal(CV.Cyclotomic(r1, 1), 1)
    @test TE.coefficients(v1) == CV.coefficients(v1)

    @test v1 + v1 == 2 * v1 # same N,M arithmetic
    @test CV.Cyclotomic(r1, 1) + CV.Cyclotomic(r1, 1) == CV.Cyclotomic(r1, 2)
    @test_throws ArgumentError v1 + v2
    @test_throws ArgumentError v1 * v2
    @test_throws ArgumentError v1 + CV.normalize_root(TE.zeta(7), CV.RootSpec(7, 1))

    @test q^(-64) == (q^(-32))^2
    @test (v1^-2) * (v1^2) == 1
    @test ComplexF64(v2) ≈ ComplexF64(CV.canonical_value(v2))
    @test !occursin("ζ^7", sprint(show, v2))
end
