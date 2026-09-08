using LinearAlgebra

@testset "finite-group fusion coordinates" begin
    trivial=MatrixRepresentation([ones(ComplexF64,1,1) for _ in 1:2])
    sign=MatrixRepresentation([ones(ComplexF64,1,1),fill(-1+0im,1,1)])
    @test representation_dim(tensor_representation(sign,sign))==1
    @test length(intertwiner_basis(sign,sign))==1
    @test isempty(intertwiner_basis(sign,trivial))

    # The full nullspace is retained rather than collapsed to one channel.
    double_trivial=MatrixRepresentation([Matrix{ComplexF64}(I,2,2) for _ in 1:2])
    multiplicity_two=intertwiner_basis(double_trivial,trivial)
    @test length(multiplicity_two)==2

    L=[ComplexF64[1;0;;],ComplexF64[0;1;;]]
    B=ComplexF64[1 1;0 1]
    R=[reshape(B[:,i],2,1) for i in 1:2]
    result=fusion_matrix_result(L,R)
    @test result.rank==2
    @test result.residual <= 1e-12
    @test R[1]*result.matrix[1,1]+R[2]*result.matrix[2,1] ≈ L[1]

    tetra=raw_tetrahedron(ComplexF64[1 2;3 5],[2,3],[5,7])
    @test tetra.positive ≈ ComplexF64[1/5 2/5;3/7 5/7]
    @test tetra.negative ≈ inv(ComplexF64[1 2;3 5])./ComplexF64[2,3]
end
