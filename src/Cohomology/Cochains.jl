"""
Multiplicative finite-group cochains over a supplied multiplication table.

The table is deliberately kept as input data rather than wrapped in a group
CAS object.  A table entry `table[g,h] == k` means `g*h == k`; all indices are
one based.  Coefficients are values in a multiplicative `Number` type, so exact
`Rational` and `Complex{<:Rational}` coefficients remain exact under the
cochain operations.
"""
module Cohomology

export GroupCochain, Cochain1, Cochain2, Cochain3,
       cochain1, cochain2, cochain3,
       group_order, group_identity, group_table,
       coboundary, coboundary2, coboundary3,
       is_normalized, is_cocycle

struct GroupCochain{N,T<:Number}
    table::Matrix{Int}
    values::Array{T,N}
    identity::Int
end

const Cochain1 = GroupCochain{1}
const Cochain2 = GroupCochain{2}
const Cochain3 = GroupCochain{3}

# Avoid fixed-width integer overflow and non-closed integer division in the
# exact coefficient path. Floating equality remains a machine comparison.
_coefficient_type(::Type{T}) where {T<:Union{Integer,Rational}} = Rational{BigInt}
_coefficient_type(::Type{Complex{T}}) where {T<:Union{Integer,Rational}} = Complex{Rational{BigInt}}
_coefficient_type(::Type{T}) where {T<:Number} = T

function _group_data(table::AbstractMatrix{<:Integer})
    rows, columns = size(table)
    rows > 0 && rows == columns ||
        throw(ArgumentError("group table must be nonempty and square"))
    data = Matrix{Int}(table)
    order = rows
    all(1 <= entry <= order for entry in data) ||
        throw(ArgumentError("group table entries must lie in 1:$order"))
    identity = findfirst(candidate ->
        all(data[candidate, index] == index &&
            data[index, candidate] == index for index in 1:order), 1:order)
    identity === nothing &&
        throw(ArgumentError("group table has no two-sided identity"))
    for g in 1:order, h in 1:order, k in 1:order
        data[data[g, h], k] == data[g, data[h, k]] ||
            throw(ArgumentError("group table is not associative at ($g,$h,$k)"))
    end
    inverses = Vector{Int}(undef, order)
    for g in 1:order
        candidates = [h for h in 1:order
                      if data[g, h] == identity && data[h, g] == identity]
        length(candidates) == 1 ||
            throw(ArgumentError("group element $g has no unique two-sided inverse"))
        inverses[g] = only(candidates)
    end
    data, identity, inverses
end

function GroupCochain{N}(table::AbstractMatrix{<:Integer},
                         values::AbstractArray{<:Number,N}) where N
    data, identity, _ = _group_data(table)
    size(values) == ntuple(_ -> size(data, 1), N) ||
        throw(DimensionMismatch("cochain values must have one group axis per argument"))
    T = _coefficient_type(eltype(values))
    all(value -> !iszero(value) && isfinite(value), values) ||
        throw(ArgumentError("multiplicative cochain values must be finite and nonzero"))
    GroupCochain{N,T}(data, Array{T,N}(values), identity)
end

cochain1(table, values) = GroupCochain{1}(table, values)
cochain2(table, values) = GroupCochain{2}(table, values)
cochain3(table, values) = GroupCochain{3}(table, values)

group_order(cochain::GroupCochain) = size(cochain.table, 1)
group_identity(cochain::GroupCochain) = cochain.identity
group_table(cochain::GroupCochain) = copy(cochain.table)

function _comparison_mode(mode::Symbol, atol, rtol)
    mode in (:exact, :approx) ||
        throw(ArgumentError("comparison mode must be :exact or :approx"))
    mode == :exact && return (left, right) -> left == right
    (left, right) -> isapprox(left, right; atol=atol, rtol=rtol)
end

"""Check the normalized boundary conditions for a multiplicative cochain."""
function is_normalized(cochain::GroupCochain{N}; mode::Symbol=:exact,
                       atol=0, rtol=0) where N
    equal = _comparison_mode(mode, atol, rtol)
    e = cochain.identity
    values = cochain.values
    if N == 1
        return equal(values[e], one(eltype(values)))
    elseif N == 2
        return all(equal(values[e, h], one(eltype(values))) &&
                   equal(values[g, e], one(eltype(values)))
                   for g in axes(values, 1), h in axes(values, 2))
    elseif N == 3
        return all(equal(values[e, h, k], one(eltype(values))) &&
                   equal(values[g, e, k], one(eltype(values))) &&
                   equal(values[g, h, e], one(eltype(values)))
                   for g in axes(values, 1), h in axes(values, 2),
                       k in axes(values, 3))
    end
    throw(ArgumentError("normalization is implemented for degrees 1, 2, and 3"))
end

"""Check the inhomogeneous multiplicative cocycle identity exactly by default."""
function is_cocycle(cochain::GroupCochain{2}; mode::Symbol=:exact,
                    atol=0, rtol=0)
    equal = _comparison_mode(mode, atol, rtol)
    values, table = cochain.values, cochain.table
    for g in axes(values, 1), h in axes(values, 2), k in axes(values, 2)
        lhs = values[h, k] * values[g, table[h, k]]
        rhs = values[table[g, h], k] * values[g, h]
        equal(lhs, rhs) || return false
    end
    true
end

function is_cocycle(cochain::GroupCochain{3}; mode::Symbol=:exact,
                    atol=0, rtol=0)
    equal = _comparison_mode(mode, atol, rtol)
    values, table = cochain.values, cochain.table
    for g in axes(values, 1), h in axes(values, 2),
        k in axes(values, 2), l in axes(values, 2)
        lhs = values[h, k, l] * values[g, table[h, k], l] * values[g, h, k]
        rhs = values[table[g, h], k, l] * values[g, h, table[k, l]]
        equal(lhs, rhs) || return false
    end
    true
end

is_cocycle(::GroupCochain{1}; kwargs...) =
    throw(ArgumentError("cocycle checks are implemented for degrees 2 and 3"))

_division_type(::Type{T}) where {T<:Integer} = Rational{BigInt}
_division_type(::Type{T}) where {T<:Number} = T

"""Return the multiplicative coboundary δf of a one-cochain f."""
function coboundary(cochain::GroupCochain{1})
    values, table = cochain.values, cochain.table
    n = size(table, 1)
    T = _division_type(eltype(values))
    source = T.(values)
    result = Array{T}(undef, n, n)
    for g in 1:n, h in 1:n
        result[g, h] = source[h] * source[g] / source[table[g, h]]
    end
    GroupCochain{2}(table, result)
end

"""Return the multiplicative coboundary δc of a two-cochain c."""
function coboundary(cochain::GroupCochain{2})
    values, table = cochain.values, cochain.table
    n = size(table, 1)
    T = _division_type(eltype(values))
    source = T.(values)
    result = Array{T}(undef, n, n, n)
    for g in 1:n, h in 1:n, k in 1:n
        result[g, h, k] = source[h, k] * source[g, table[h, k]] /
            (source[table[g, h], k] * source[g, h])
    end
    GroupCochain{3}(table, result)
end

coboundary2(cochain::GroupCochain{1}) = coboundary(cochain)
coboundary3(cochain::GroupCochain{2}) = coboundary(cochain)

end
