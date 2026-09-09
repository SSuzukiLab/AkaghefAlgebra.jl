"""Public rooted values over the internal canonical cyclotomic field.

`TaftExact.Cyclotomic{N}` is the calculation type: its coefficients are
always in the power basis of the canonical root `ζ_N`.  `Cyclotomic{N,M}`
adds the selected primitive root `ζ_N^M` to a value whose coefficients have
already had the Galois substitution applied exactly once.
"""
module CyclotomicValues

import ..TaftExact as TE
using Base: +, -, *, /, ^, ==, zero, one, iszero, isequal

export RootSpec, Cyclotomic, normalize_root, root_spec, specified_root,
       canonical_value, with_root, coefficients

const _Scalar = Union{Integer,Rational,BigInt}

"""A normalized primitive-root specification for `Q(ζ_N)`."""
struct RootSpec{N,M}
    N::Int
    M::Int
    function RootSpec{N,M}(n::Int, m::Int) where {N,M}
        Int(N) == n && Int(M) == m && n >= 2 && 0 <= m < n && gcd(n, m) == 1 ||
            throw(ArgumentError("invalid normalized primitive-root specification"))
        new{N,M}(n, m)
    end
end

function RootSpec(N::Integer, M::Integer=1)
    n = try
        Int(N)
    catch
        throw(ArgumentError("cyclotomic order N must fit in Int"))
    end
    n >= 2 || throw(ArgumentError("cyclotomic order N must be at least 2"))
    m = Int(mod(BigInt(M), BigInt(n)))
    gcd(n, m) == 1 ||
        throw(ArgumentError("M must be coprime to N so that ζ_N^M is primitive"))
    RootSpec{n,m}(n, m)
end

function Base.show(io::IO, r::RootSpec{N,M}) where {N,M}
    print(io, "RootSpec(", N, ", ", M, ")")
end

"""A value with coefficients in the canonical `ζ_N` power basis and root metadata."""
struct Cyclotomic{N,M} <: Number
    canonical::TE.Cyclotomic{N}
    spec::RootSpec{N,M}
    function Cyclotomic{N,M}(x::TE.Cyclotomic{N}, spec::RootSpec{N,M}) where {N,M}
        # Do not retain the caller's mutable coefficient vector.
        new{N,M}(TE.Cyclotomic{N}(TE.coefficients(x)), spec)
    end
end

function Cyclotomic{N,M}(x::TE.Cyclotomic{N}) where {N,M}
    spec = RootSpec(Int(N), Int(M))
    Cyclotomic{N,M}(x, spec)
end

Cyclotomic(x::TE.Cyclotomic{N}, spec::RootSpec{N,M}) where {N,M} =
    Cyclotomic{N,M}(x, spec)
Cyclotomic(spec::RootSpec{N,M}, x::TE.Cyclotomic{N}) where {N,M} =
    Cyclotomic{N,M}(x, spec)
Cyclotomic(spec::RootSpec{N,M}, x::_Scalar) where {N,M} =
    Cyclotomic(TE.Cyclotomic{N}(x), spec)
Cyclotomic(N::Integer, M::Integer, x::TE.Cyclotomic) =
    normalize_root(x, RootSpec(N, M))
Cyclotomic(N::Integer, M::Integer, x::_Scalar) =
    Cyclotomic(TE.Cyclotomic{Int(N)}(x), RootSpec(N, M))

root_spec(x::Cyclotomic) = x.spec
canonical_value(x::Cyclotomic{N}) where {N} = TE.Cyclotomic{N}(TE.coefficients(x.canonical))
canonical_value(x::TE.Cyclotomic) = x
coefficients(x::Cyclotomic) = TE.coefficients(x.canonical)
TE.coefficients(x::Cyclotomic) = coefficients(x)

"""The exact internal value of the selected primitive root `ζ_N^M`."""
specified_root(r::RootSpec{N,M}) where {N,M} = TE.zeta(Int(N))^Int(M)
specified_root(x::Cyclotomic) = specified_root(root_spec(x))

function Base.show(io::IO, x::Cyclotomic{N,M}) where {N,M}
    # `canonical` is already normalized; showing it must not apply M again.
    print(io, "Cyclotomic{", N, ",", M, "}(")
    show(io, x.canonical)
    print(io, ")")
end

Base.ComplexF64(x::Cyclotomic) = ComplexF64(x.canonical)
Base.iszero(x::Cyclotomic) = iszero(x.canonical)
Base.zero(::Type{Cyclotomic{N,M}}) where {N,M} =
    Cyclotomic{N,M}(TE.Cyclotomic{N}(0))
Base.one(::Type{Cyclotomic{N,M}}) where {N,M} =
    Cyclotomic{N,M}(TE.Cyclotomic{N}(1))
Base.zero(x::Cyclotomic{N,M}) where {N,M} = zero(Cyclotomic{N,M})
Base.one(x::Cyclotomic{N,M}) where {N,M} = one(Cyclotomic{N,M})

function _galois_substitute(x::TE.Cyclotomic{N}, m::Int) where {N}
    n = Int(N)
    q = TE.zeta(n)
    out = zero(q)
    for (j, c) in enumerate(TE.coefficients(x))
        iszero(c) || (out += c * q^(mod((j - 1) * m, n)))
    end
    out
end

"""Apply the root substitution `ζ_N ↦ ζ_N^M` once to a canonical value."""
function normalize_root(x::TE.Cyclotomic{N}, r::RootSpec{N,M}) where {N,M}
    Cyclotomic{N,M}(_galois_substitute(x, Int(M)))
end

function normalize_root(x::TE.Cyclotomic{N}, r::RootSpec{L,P}) where {N,L,P}
    Int(N) == Int(L) ||
        throw(ArgumentError("root order N does not match the cyclotomic value"))
    normalize_root(x, RootSpec(Int(L), Int(P)))
end

function normalize_root(x::Cyclotomic{N,M}, r::RootSpec{L,P}) where {N,M,L,P}
    Int(N) == Int(L) && Int(M) == Int(P) && return x
    throw(ArgumentError("a finalized rooted value cannot be normalized again or relabeled by a different root"))
end

"""Relabel a finalized value without changing its canonical coefficients."""
function with_root(x::Cyclotomic{N,M}, r::RootSpec{L,P}) where {N,M,L,P}
    Int(N) == Int(L) ||
        throw(ArgumentError("root order N does not match the cyclotomic value"))
    Cyclotomic{N,P}(x.canonical)
end

with_root(x::Cyclotomic, N::Integer, M::Integer=1) = with_root(x, RootSpec(N, M))

function Base.:(==)(a::Cyclotomic{N,M}, b::Cyclotomic{N,P}) where {N,M,P}
    a.canonical == b.canonical
end
Base.:(==)(a::Cyclotomic{N,M}, b::Cyclotomic{L,P}) where {N,M,L,P} = false
Base.:(==)(a::Cyclotomic{N,M}, b::TE.Cyclotomic{N}) where {N,M} = a.canonical == b
Base.:(==)(a::TE.Cyclotomic{N}, b::Cyclotomic{N,M}) where {N,M} = a == b.canonical
Base.:(==)(a::Cyclotomic{N,M}, b::TE.Cyclotomic{L}) where {N,M,L} = false
Base.:(==)(a::TE.Cyclotomic{L}, b::Cyclotomic{N,M}) where {L,N,M} = false
Base.:(==)(a::Cyclotomic{N,M}, b::_Scalar) where {N,M} = a.canonical == TE.Cyclotomic{N}(b)
Base.:(==)(a::_Scalar, b::Cyclotomic{N,M}) where {N,M} = b == a

Base.isequal(a::Cyclotomic{N,M}, b::Cyclotomic{N,P}) where {N,M,P} =
    isequal(a.canonical, b.canonical)
Base.isequal(a::Cyclotomic{N,M}, b::Cyclotomic{L,P}) where {N,M,L,P} = false
# Internal values have a pre-existing structural hash.  Keep `isequal` false
# across this public/internal boundary while retaining mathematical `==`.
Base.isequal(a::Cyclotomic{N,M}, b::TE.Cyclotomic) where {N,M} = false
Base.isequal(a::TE.Cyclotomic, b::Cyclotomic) = false
Base.isequal(a::Cyclotomic, b::_Scalar) = false
Base.isequal(a::_Scalar, b::Cyclotomic) = false

function Base.hash(x::Cyclotomic{N,M}, h::UInt) where {N,M}
    hash((Int(N), Tuple(TE.coefficients(x.canonical))), h)
end

function _binary_mismatch(op, a::Cyclotomic{N,M}, b::Cyclotomic{L,P}) where {N,M,L,P}
    throw(ArgumentError("$op requires rooted cyclotomic values with the same N and M"))
end

function Base.:+(a::Cyclotomic{N,M}, b::Cyclotomic{N,M}) where {N,M}
    Cyclotomic{N,M}(a.canonical + b.canonical)
end
Base.:+(a::Cyclotomic{N,M}, b::Cyclotomic{L,P}) where {N,M,L,P} = _binary_mismatch("addition", a, b)
Base.:+(a::Cyclotomic{N,M}, b::_Scalar) where {N,M} = Cyclotomic{N,M}(a.canonical + b)
Base.:+(a::_Scalar, b::Cyclotomic{N,M}) where {N,M} = b + a

Base.:-(a::Cyclotomic{N,M}) where {N,M} = Cyclotomic{N,M}(-a.canonical)
function Base.:-(a::Cyclotomic{N,M}, b::Cyclotomic{N,M}) where {N,M}
    Cyclotomic{N,M}(a.canonical - b.canonical)
end
Base.:-(a::Cyclotomic{N,M}, b::Cyclotomic{L,P}) where {N,M,L,P} = _binary_mismatch("subtraction", a, b)
Base.:-(a::Cyclotomic{N,M}, b::_Scalar) where {N,M} = Cyclotomic{N,M}(a.canonical - b)
Base.:-(a::_Scalar, b::Cyclotomic{N,M}) where {N,M} = Cyclotomic{N,M}(a - b.canonical)

function Base.:*(a::Cyclotomic{N,M}, b::Cyclotomic{N,M}) where {N,M}
    Cyclotomic{N,M}(a.canonical * b.canonical)
end
Base.:*(a::Cyclotomic{N,M}, b::Cyclotomic{L,P}) where {N,M,L,P} = _binary_mismatch("multiplication", a, b)
Base.:*(a::Cyclotomic{N,M}, b::_Scalar) where {N,M} = Cyclotomic{N,M}(a.canonical * b)
Base.:*(a::_Scalar, b::Cyclotomic{N,M}) where {N,M} = b * a

function Base.:/(a::Cyclotomic{N,M}, b::Cyclotomic{N,M}) where {N,M}
    Cyclotomic{N,M}(a.canonical / b.canonical)
end
Base.:/(a::Cyclotomic{N,M}, b::Cyclotomic{L,P}) where {N,M,L,P} = _binary_mismatch("division", a, b)
Base.:/(a::Cyclotomic{N,M}, b::_Scalar) where {N,M} = Cyclotomic{N,M}(a.canonical / b)
Base.:/(a::_Scalar, b::Cyclotomic{N,M}) where {N,M} = Cyclotomic{N,M}(TE.Cyclotomic{N}(a) / b.canonical)

function Base.:^(a::Cyclotomic{N,M}, n::Integer) where {N,M}
    Cyclotomic{N,M}(a.canonical ^ n)
end

end # module CyclotomicValues
