using Test
using AkaghefAlgebra

const S3_COCHAIN_TABLE = [
    1 2 3 4 5 6
    2 3 1 6 4 5
    3 1 2 5 6 4
    4 5 6 1 2 3
    5 6 4 3 1 2
    6 4 5 2 3 1
]

@testset "finite group cochains" begin
    n = size(S3_COCHAIN_TABLE, 1)
    values = fill(BigInt(1) // BigInt(1), n, n)
    two = Cochain2(S3_COCHAIN_TABLE, values)

    @test group_order(two) == 6
    @test group_identity(two) == 1
    @test S3_COCHAIN_TABLE[2, 4] != S3_COCHAIN_TABLE[4, 2]
    @test group_table(two) == S3_COCHAIN_TABLE
    @test is_normalized(two)
    @test is_cocycle(two)

    complex_values = fill(Complex{Rational{Int}}(1 // 1, 0 // 1), n, n)
    @test is_normalized(Cochain2(S3_COCHAIN_TABLE, complex_values))
end

@testset "exact coefficient boundary" begin
    z2 = [1 2; 2 1]
    f = Cochain1(z2, [1, typemax(Int)])
    @test eltype(f.values) == Rational{BigInt}
    @test coboundary(f).values[2,2] == big(typemax(Int))^2 // 1
    h = Cochain1(z2, Complex{Int}[1, 1+im])
    @test coboundary(h).values[2,2] == Complex{Rational{BigInt}}(0,2)
    @test_throws ArgumentError Cochain1(z2, [1.0, Inf])
    @test_throws ArgumentError Cochain1(z2, [1.0, NaN])
end

@testset "normalization and cocycle failures are explicit" begin
    n = size(S3_COCHAIN_TABLE, 1)
    values = fill(BigInt(1) // BigInt(1), n, n)
    values[1, 2] = 2 // 1
    @test !is_normalized(Cochain2(S3_COCHAIN_TABLE, values))

    values = fill(BigInt(1) // BigInt(1), n, n)
    values[2, 3] = 2 // 1
    @test !is_cocycle(Cochain2(S3_COCHAIN_TABLE, values))
    @test_throws ArgumentError is_cocycle(
        Cochain2(S3_COCHAIN_TABLE, fill(BigInt(1) // BigInt(1), n, n));
        mode=:bad)
end

@testset "multiplicative coboundary squares to one" begin
    n = size(S3_COCHAIN_TABLE, 1)
    values = BigInt[1, 2, 3, 4, 5, 6]
    one = Cochain1(S3_COCHAIN_TABLE, values)
    two = coboundary2(one)
    three = coboundary3(two)

    @test is_cocycle(two)
    @test is_cocycle(three)
    @test all(value == 1 for value in three.values)
end

@testset "exact normalized Z/2 representative" begin
    z2 = [1 2; 2 1]
    values = Array{Rational{BigInt}}(undef, 2, 2, 2)
    for a in 1:2, b in 1:2, c in 1:2
        values[a, b, c] = isone((a - 1) * (b - 1) * (c - 1)) ? -1 // 1 : 1 // 1
    end
    three = Cochain3(z2, values)
    @test is_normalized(three)
    @test is_cocycle(three)
    @test three.values[2, 2, 2] == -1 // 1
    @test is_cocycle(three; mode=:approx, atol=0, rtol=0)
end
