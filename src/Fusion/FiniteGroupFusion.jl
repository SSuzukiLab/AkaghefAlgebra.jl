module FiniteGroupFusion

using LinearAlgebra

export MatrixRepresentation, representation_dim, tensor_representation,
       intertwiner_basis, fusion_matrix, fusion_matrix_result,
       TraceDualBasis, trace_dual_basis, RawTetrahedron, raw_tetrahedron,
       contract_trace_dual

"""A finite-group representation given by one matrix for each indexed element."""
struct MatrixRepresentation{T<:Number}
    matrices::Vector{Matrix{T}}
    function MatrixRepresentation(matrices::AbstractVector{<:AbstractMatrix{T}}) where {T<:Number}
        isempty(matrices) && throw(ArgumentError("a representation needs group matrices"))
        d=size(first(matrices),1)
        d>0 && all(size(m)==(d,d) for m in matrices) ||
            throw(DimensionMismatch("representation matrices must be nonempty and square with one common size"))
        new{T}(Matrix{T}.(matrices))
    end
end

representation_dim(rep::MatrixRepresentation)=size(first(rep.matrices),1)

function tensor_representation(left::MatrixRepresentation,
                               right::MatrixRepresentation)
    length(left.matrices)==length(right.matrices) ||
        throw(DimensionMismatch("representations must use the same indexed group"))
    MatrixRepresentation([kron(a,b) for (a,b) in zip(left.matrices,right.matrices)])
end

_complex_matrix(x)=Matrix{ComplexF64}(x)

"""
Return an orthonormal numerical basis of `Hom(source,target)`.

Each returned matrix `X` satisfies `target(g)X = X source(g)`.  The SVD rank
cutoff is relative to the largest singular value, and the residual is checked
against the same tolerance.  This routine deliberately preserves all nullspace
columns, including fusion multiplicities greater than one.
"""
function intertwiner_basis(source::MatrixRepresentation,
                           target::MatrixRepresentation; rtol::Real=1e-11)
    length(source.matrices)==length(target.matrices) ||
        throw(DimensionMismatch("representations must use the same indexed group"))
    ds,dt=representation_dim(source),representation_dim(target)
    blocks=Matrix{ComplexF64}[]
    for (S,T) in zip(source.matrices,target.matrices)
        # vec(TX-XS) = (I_ds kron T - transpose(S) kron I_dt) vec(X)
        push!(blocks,kron(Matrix{ComplexF64}(I,ds,ds),_complex_matrix(T)) -
                     kron(transpose(_complex_matrix(S)),Matrix{ComplexF64}(I,dt,dt)))
    end
    equations=reduce(vcat,blocks)
    decomposition=svd(equations;full=true)
    scale=isempty(decomposition.S) ? 0.0 : maximum(decomposition.S)
    rank=count(>(rtol*max(scale,1.0)),decomposition.S)
    vectors=decomposition.V[:,rank+1:end]
    basis=[reshape(copy(vectors[:,j]),dt,ds) for j in axes(vectors,2)]
    for X in basis, (S,T) in zip(source.matrices,target.matrices)
        norm(_complex_matrix(T)*X-X*_complex_matrix(S)) <=
            rtol*max(1,norm(X),norm(S),norm(T))*20 ||
            throw(ArgumentError("intertwiner nullspace residual exceeds tolerance"))
    end
    basis
end

"""Trace-dual basis for a nondegenerate bilinear trace pairing."""
struct TraceDualBasis{T<:Number}
    basis::Vector{Matrix{T}}
    dual::Vector{Matrix{T}}
    gram::Matrix{T}
end

function trace_dual_basis(basis::AbstractVector{<:AbstractMatrix{T}}) where {T<:Number}
    isempty(basis) && return TraceDualBasis(Matrix{T}[],Matrix{T}[],zeros(T,0,0))
    all(size(x)==size(first(basis)) for x in basis) ||
        throw(DimensionMismatch("paired maps must have a common shape"))
    gram=[tr(adjoint(x)*y) for x in basis,y in basis]
    rank(gram)==length(basis) || throw(ArgumentError("trace pairing is degenerate"))
    inverse_gram=inv(gram)
    dual=[sum(inverse_gram[j,i]*basis[j] for j in eachindex(basis))
          for i in eachindex(basis)]
    TraceDualBasis(Matrix{T}.(basis),Matrix{T}.(dual),Matrix{T}(gram))
end

_columnize(maps)=hcat((vec(x) for x in maps)...)

"""
Solve the fixed convention `R * F = associator * L`.

Columns of `L` are left-associated fusion embeddings and columns of `R` are
right-associated embeddings.  `associator` maps the left ambient tensor space
to the right ambient tensor space.  Rows of `F` therefore carry right-tree
channels and columns carry left-tree channels.  This convention is also the
one used by `raw_tetrahedron`; changing it reverses non-real cocycle phases.
"""
function fusion_matrix_result(left_maps,right_maps,associator=nothing;rtol::Real=1e-10)
    length(left_maps)==length(right_maps) ||
        throw(DimensionMismatch("left and right fusion decompositions must have equal dimension"))
    isempty(left_maps) && return (;matrix=zeros(ComplexF64,0,0),residual=0.0,rank=0)
    L=_columnize(left_maps); R=_columnize(right_maps)
    A=associator===nothing ? Matrix{ComplexF64}(I,size(L,1),size(L,1)) :
      associator isa UniformScaling ? Matrix{ComplexF64}(associator,size(L,1),size(L,1)) :
                                      _complex_matrix(associator)
    size(A)==(size(R,1),size(L,1)) ||
        throw(DimensionMismatch("associator has incompatible ambient spaces"))
    rr=rank(R;rtol=rtol)
    rr==size(R,2) || throw(ArgumentError("right fusion maps are not a basis"))
    F=R\(A*L)
    residual=norm(R*F-A*L)
    residual <= rtol*max(1,norm(A*L))*20 ||
        throw(ArgumentError("fusion trees do not span the same Hom space"))
    (;matrix=F,residual,rank=rr)
end

fusion_matrix(left_maps,right_maps,associator=nothing;kwargs...)=
    fusion_matrix_result(left_maps,right_maps,associator;kwargs...).matrix

"""Raw tetrahedron map together with the trace-dual gluing data on its faces."""
struct RawTetrahedron{T<:Number}
    positive::Matrix{T}
    negative::Matrix{T}
    right_dimensions::Vector{T}
    left_dimensions::Vector{T}
end

"""
Build dimension-weighted trace-dual tetrahedron coordinates.

For `R*F=alpha*L`, rows of `F` are right channels `y` and columns are left
channels `x`.  The positive raw tetrahedron is `F[y,x]/dim(y)`.  The negative
one is `inv(F)[x,y]/dim(x)`.  Thus a raw tetrahedron is not merely an F matrix;
the trace-dual face closure supplies the displayed dimension denominator.
"""
function raw_tetrahedron(F::AbstractMatrix{T},left_dimensions,
                         right_dimensions) where {T<:Number}
    size(F,2)==length(left_dimensions) ||
        throw(DimensionMismatch("left channel dimensions do not match F columns"))
    size(F,1)==length(right_dimensions) ||
        throw(DimensionMismatch("right channel dimensions do not match F rows"))
    all(!iszero,left_dimensions) && all(!iszero,right_dimensions) ||
        throw(ArgumentError("channel dimensions must be nonzero"))
    positive=Matrix{T}(F)
    for y in axes(positive,1)
        positive[y,:]./=right_dimensions[y]
    end
    inverse_F=inv(Matrix{T}(F))
    negative=copy(inverse_F)
    for x in axes(negative,1)
        negative[x,:]./=left_dimensions[x]
    end
    RawTetrahedron(positive,negative,T.(right_dimensions),T.(left_dimensions))
end

"""Contract one output face coordinate with one input face by trace duality."""
function contract_trace_dual(left::AbstractVector,right::AbstractVector,
                             pairing::TraceDualBasis)
    length(left)==length(right)==length(pairing.basis) ||
        throw(DimensionMismatch("face coordinate dimensions disagree"))
    sum(left[i]*right[i] for i in eachindex(left))
end

end
