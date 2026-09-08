"""Exact sparse Taft algebra over `Q(ζ_N)`.

The algebra uses

    K*E = ζ*E*K,  E^N = 0,  K^N = 1,
    Δ(E) = 1⊗E + E⊗K,  Δ(K) = K⊗K,

and the PBW order `E^r*K^s`.  Coefficients are represented by the small
cyclotomic field implementation below; no dense structure-constant arrays are
constructed.
"""
module TaftExact

using Base: +, -, *, /, ^, ==, zero, one

export Cyclotomic, coefficients, zeta, qbinomial, qfactorial,
       TaftElement, TaftTensor, taft_element, monomial, E, K, one_element,
       zero_element, tensor_terms, coproduct, counit, antipode,
       DualMonomial, DualBasis, dual_monomial, dual_basis,
       dual_basis_expansion, pairing

# ---------------------------------------------------------------------------
# Exact cyclotomic arithmetic
# ---------------------------------------------------------------------------

const BigRat = Rational{BigInt}

_br(x::BigInt) = x // BigInt(1)
_br(x::Integer) = BigInt(x) // BigInt(1)
_br(x::Rational) = BigInt(numerator(x)) // BigInt(denominator(x))
_br(x::BigRat) = x

function _trim_int(p::Vector{BigInt})
    isempty(p) && return BigInt[0]
    q = copy(p)
    while length(q) > 1 && q[end] == 0
        pop!(q)
    end
    q
end

function _int_mul(a::Vector{BigInt}, b::Vector{BigInt})
    out = fill(BigInt(0), length(a) + length(b) - 1)
    for i in eachindex(a), j in eachindex(b)
        out[i+j-1] += a[i] * b[j]
    end
    _trim_int(out)
end

function _int_divexact(a::Vector{BigInt}, b::Vector{BigInt})
    aa = copy(_trim_int(a)); bb = _trim_int(b)
    isempty(bb) || bb[end] != 0 || throw(DivideError())
    db = length(bb) - 1
    db < 0 && throw(DivideError())
    lc = bb[end]
    lc == 1 || lc == -1 || throw(ArgumentError("cyclotomic divisor must be monic"))
    length(aa) < length(bb) && return BigInt[0]
    q = fill(BigInt(0), length(aa) - db)
    while length(aa) - 1 >= db && !(length(aa) == 1 && aa[1] == 0)
        d = length(aa) - 1 - db
        c = aa[end] * lc
        q[d+1] += c
        for j in 0:db
            aa[d+j+1] -= c * bb[j+1]
        end
        aa = _trim_int(aa)
    end
    all(iszero, aa) || throw(ArgumentError("non-exact polynomial division"))
    _trim_int(q)
end

function _divisors(n::Int)
    n > 0 || throw(ArgumentError("cyclotomic order must be positive"))
    ds = Int[]
    for d in 1:floor(Int, sqrt(n))
        n % d == 0 || continue
        push!(ds, d)
        d*d != n && push!(ds, n ÷ d)
    end
    sort!(ds)
end

const _PHI_CACHE = Dict{Int,Vector{BigInt}}()

"""Return the monic integer cyclotomic polynomial Φ_n in ascending order."""
function _cyclotomic_poly(n::Int)
    get!(_PHI_CACHE, n) do
        den = BigInt[1]
        for d in _divisors(n)
            d == n || (den = _int_mul(den, _cyclotomic_poly(d)))
        end
        xn = fill(BigInt(0), n + 1); xn[1] = -1; xn[end] = 1
        _int_divexact(xn, den)
    end
end

"""An element of the exact field `Q(ζ_N)` in the basis `1,ζ,…,ζ^(φ(N)-1)`."""
struct Cyclotomic{N} <: Number
    coeffs::Vector{BigRat}
    function Cyclotomic{N}(coeffs::AbstractVector) where {N}
        n = Int(N); n > 0 || throw(ArgumentError("N must be positive"))
        phi = _cyclotomic_poly(n)
        raw = [_br(c) for c in coeffs]
        isempty(raw) && push!(raw, _br(0))
        while length(raw) > 1 && iszero(raw[end]); pop!(raw); end
        # Rational polynomial remainder modulo monic Φ_N.
        while length(raw) >= length(phi)
            d = length(raw) - length(phi)
            c = raw[end] / _br(phi[end])
            for j in eachindex(phi)
                raw[d+j] -= c * _br(phi[j])
            end
            while length(raw) > 1 && iszero(raw[end]); pop!(raw); end
        end
        while length(raw) < length(phi)-1; push!(raw, _br(0)); end
        new{N}(raw)
    end
end

Cyclotomic(N::Integer, coeffs::AbstractVector) = Cyclotomic{Int(N)}(coeffs)
Cyclotomic{N}(x::Union{BigInt,Integer,Rational}) where {N} = Cyclotomic{N}([x])
Cyclotomic(N::Integer, x::Union{BigInt,Integer,Rational}) = Cyclotomic{Int(N)}([x])
Base.convert(::Type{Cyclotomic{N}}, x::Cyclotomic{N}) where {N} = x
Base.convert(::Type{Cyclotomic{N}}, x::Union{BigInt,Integer,Rational}) where {N} = Cyclotomic{N}(x)
Base.promote_rule(::Type{Cyclotomic{N}}, ::Type{S}) where {N,S<:Integer} = Cyclotomic{N}
Base.promote_rule(::Type{Cyclotomic{N}}, ::Type{Rational{S}}) where {N,S<:Integer} = Cyclotomic{N}
Base.promote_rule(::Type{S}, ::Type{Cyclotomic{N}}) where {N,S<:Integer} = Cyclotomic{N}
Base.promote_rule(::Type{Rational{S}}, ::Type{Cyclotomic{N}}) where {N,S<:Integer} = Cyclotomic{N}
coefficients(x::Cyclotomic) = copy(x.coeffs)

function Base.show(io::IO, x::Cyclotomic{N}) where {N}
    c = x.coeffs; terms = String[]
    for (i,a) in enumerate(c)
        iszero(a) && continue
        term = i == 1 ? string(a) : (i == 2 ? string(a, "*ζ") : string(a, "*ζ^", i-1))
        push!(terms, term)
    end
    print(io, isempty(terms) ? "0" : join(terms, " + "))
end

Base.iszero(x::Cyclotomic) = all(iszero, x.coeffs)
Base.zero(::Type{Cyclotomic{N}}) where {N} = Cyclotomic{N}([0])
Base.one(::Type{Cyclotomic{N}}) where {N} = Cyclotomic{N}([1])
Base.zero(x::Cyclotomic{N}) where {N} = zero(Cyclotomic{N})
Base.one(x::Cyclotomic{N}) where {N} = one(Cyclotomic{N})

function Base.:(==)(a::Cyclotomic{N}, b::Cyclotomic{N}) where {N}
    a.coeffs == b.coeffs
end
Base.:(==)(a::Cyclotomic{N}, b::Union{Integer,Rational,BigInt}) where {N} = a == Cyclotomic{N}(b)
Base.:(==)(a::Union{Integer,Rational,BigInt}, b::Cyclotomic{N}) where {N} = Cyclotomic{N}(a) == b

function Base.:+(a::Cyclotomic{N}, b::Cyclotomic{N}) where {N}
    Cyclotomic{N}([a.coeffs[i] + b.coeffs[i] for i in eachindex(a.coeffs)])
end
Base.:-(a::Cyclotomic{N}) where {N} = Cyclotomic{N}([-x for x in a.coeffs])
Base.:-(a::Cyclotomic{N}, b::Cyclotomic{N}) where {N} = a + (-b)
Base.:+(a::Cyclotomic{N}, b::Union{Integer,Rational,BigInt}) where {N} = a + Cyclotomic{N}(b)
Base.:+(a::Union{Integer,Rational,BigInt}, b::Cyclotomic{N}) where {N} = Cyclotomic{N}(a) + b
Base.:-(a::Cyclotomic{N}, b::Union{Integer,Rational,BigInt}) where {N} = a - Cyclotomic{N}(b)
Base.:-(a::Union{Integer,Rational,BigInt}, b::Cyclotomic{N}) where {N} = Cyclotomic{N}(a) - b

function Base.:*(a::Cyclotomic{N}, b::Cyclotomic{N}) where {N}
    raw = fill(_br(0), length(a.coeffs) + length(b.coeffs) - 1)
    for i in eachindex(a.coeffs), j in eachindex(b.coeffs)
        raw[i+j-1] += a.coeffs[i] * b.coeffs[j]
    end
    Cyclotomic{N}(raw)
end
Base.:*(a::Cyclotomic{N}, b::Union{Integer,Rational,BigInt}) where {N} = Cyclotomic{N}([x*b for x in a.coeffs])
Base.:*(a::Union{Integer,Rational,BigInt}, b::Cyclotomic{N}) where {N} = b*a

function Base.:/(a::Cyclotomic{N}, b::Cyclotomic{N}) where {N}
    iszero(b) && throw(DivideError())
    # Solve the exact multiplication matrix for b*x=a in the power basis.
    d = length(a.coeffs)
    M = Matrix{BigRat}(undef, d, d)
    for j in 1:d
        p = fill(_br(0), d); p[j] = _br(1)
        q = Cyclotomic{N}(vcat(fill(_br(0), j-1), b.coeffs)).coeffs
        M[:,j] = q
    end
    rhs = copy(a.coeffs)
    # Bare Gaussian elimination over BigRat (stdlib only).
    for k in 1:d
        piv = findfirst(i -> !iszero(M[i,k]), k:d)
        piv === nothing && throw(DivideError("nonunit cyclotomic element"))
        piv += k - 1
        if piv != k
            M[k,:], M[piv,:] = M[piv,:], M[k,:]
            rhs[k], rhs[piv] = rhs[piv], rhs[k]
        end
        scale = M[k,k]; M[k,:] ./= scale; rhs[k] /= scale
        for i in 1:d
            i == k && continue
            c = M[i,k]; iszero(c) && continue
            M[i,:] .-= c .* M[k,:]; rhs[i] -= c*rhs[k]
        end
    end
    Cyclotomic{N}(rhs)
end
Base.inv(a::Cyclotomic) = one(a) / a
Base.:/(a::Cyclotomic{N}, b::Union{Integer,Rational,BigInt}) where {N} = a * Cyclotomic{N}(b)^(-1)
Base.:/(a::Union{Integer,Rational,BigInt}, b::Cyclotomic{N}) where {N} = Cyclotomic{N}(a) / b

function Base.:^(a::Cyclotomic{N}, n::Integer) where {N}
    n < 0 && return one(a) / (a^(-n))
    out = one(a); base = a; k = n
    while k > 0
        isodd(k) && (out *= base)
        k >>= 1
        k > 0 && (base *= base)
    end
    out
end

function zeta(N::Integer)
    N > 0 || throw(ArgumentError("N must be positive"))
    Cyclotomic{Int(N)}([0, 1])
end

function qbinomial(N::Integer, n::Integer, k::Integer)
    n >= 0 && 0 <= k <= n || throw(ArgumentError("require 0 ≤ k ≤ n"))
    q = zeta(N)
    row = [zero(q) for _ in 0:n]; row[1] = one(q)
    for m in 1:n
        nxt = [zero(q) for _ in 0:m]
        for j in 0:m
            j <= m-1 && (nxt[j+1] += q^j * row[j+1])
            j >= 1 && (nxt[j+1] += row[j])
        end
        row = nxt
    end
    row[k+1]
end

function qfactorial(N::Integer, r::Integer)
    r >= 0 || throw(ArgumentError("r must be nonnegative"))
    q = zeta(N); out = one(q)
    for m in 1:r
        out *= sum((q^j for j in 0:m-1); init=zero(q))
    end
    out
end

# ---------------------------------------------------------------------------
# Sparse Taft elements and tensors
# ---------------------------------------------------------------------------

struct TaftElement{N,T}
    terms::Dict{Tuple{Int,Int},T}
    function TaftElement{N,T}(terms::Dict{Tuple{Int,Int},T}) where {N,T}
        n = Int(N); n >= 2 || throw(ArgumentError("Taft order N must be at least 2"))
        T === Any && throw(ArgumentError("TaftExact requires a concrete exact coefficient type"))
        T <: Number || throw(ArgumentError("TaftExact coefficients must be numeric"))
        (T <: AbstractFloat || T <: Complex) &&
            throw(ArgumentError("TaftExact coefficients must be exact integers, rationals, or Cyclotomic values"))
        if T isa DataType && T <: Cyclotomic
            Int(T.parameters[1]) == n ||
                throw(ArgumentError("coefficient field order must equal Taft order N"))
        end
        out = Dict{Tuple{Int,Int},T}()
        for ((r,s), c) in terms
            0 <= r < n && 0 <= s < n ||
                throw(ArgumentError("PBW exponents must lie in 0:N-1"))
            iszero(c) || (out[(r,s)] = get(out,(r,s),zero(c)) + c)
        end
        # Recheck after combining duplicate keys (useful for externally built
        # dictionaries), and keep no mutable alias to the caller's dictionary.
        filter!(p -> !iszero(p.second), out)
        new{N,T}(out)
    end
end

function TaftElement{N,T}(terms::AbstractDict) where {N,T}
    n = Int(N); out = Dict{Tuple{Int,Int},T}()
    for (key, value) in terms
        r,s = key
        0 <= r < n && 0 <= s < n || throw(ArgumentError("PBW exponents must lie in 0:N-1"))
        c = convert(T, value); iszero(c) || (out[(Int(r),Int(s))] = get(out,(Int(r),Int(s)),zero(c)) + c)
    end
    TaftElement{N,T}(out)
end

function TaftElement(N::Integer, terms::AbstractDict)
    n = Int(N); vals = collect(values(terms))
    all(v -> v isa Union{Integer,Rational,Cyclotomic}, vals) ||
        throw(ArgumentError("TaftExact coefficients must be exact integers, rationals, or Cyclotomic values"))
    cv = [v for v in vals if v isa Cyclotomic]
    all(v -> Int(typeof(v).parameters[1]) == n, cv) ||
        throw(ArgumentError("all Cyclotomic coefficients must belong to Q(ζ_N)"))
    T = if isempty(vals) || all(v -> v isa Union{Integer,Rational}, vals)
        Cyclotomic{n}
    else
        promote_type((typeof(v) for v in vals)...)
    end
    TaftElement{n,T}(terms)
end

function taft_element(N::Integer; terms=nothing)
    terms === nothing && (terms = Dict{Tuple{Int,Int},Cyclotomic{Int(N)}}())
    TaftElement(Int(N), terms)
end
function monomial(N::Integer, r::Integer, s::Integer, c=1)
    n = Int(N); T = c isa Cyclotomic ? typeof(c) : Cyclotomic{n}
    TaftElement{n,T}(Dict((Int(r),Int(s)) => convert(T,c)))
end
E(N::Integer) = monomial(N, 1, 0, 1)
K(N::Integer) = monomial(N, 0, 1, 1)
one_element(N::Integer) = monomial(N, 0, 0, 1)
zero_element(N::Integer) = TaftElement(Int(N), Dict{Tuple{Int,Int},Cyclotomic{Int(N)}}())

Base.iszero(x::TaftElement) = isempty(x.terms)
Base.:(==)(a::TaftElement{N}, b::TaftElement{N}) where {N} = a.terms == b.terms
Base.zero(x::TaftElement{N,T}) where {N,T} = TaftElement{N,T}(Dict{Tuple{Int,Int},T}())
Base.one(x::TaftElement) = one_element(typeof(x).parameters[1])
coefficients(x::TaftElement) = copy(x.terms)

function _same_add(a::TaftElement{N,T}, b::TaftElement{N,S}) where {N,T,S}
    U = promote_type(T,S); d = Dict{Tuple{Int,Int},U}()
    for (k,v) in a.terms; d[k] = convert(U,v); end
    for (k,v) in b.terms; d[k] = get(d,k,zero(U)) + convert(U,v); end
    TaftElement{N,U}(d)
end
Base.:+(a::TaftElement{N,T}, b::TaftElement{N,S}) where {N,T,S} = _same_add(a,b)
Base.:+(a::TaftElement{N,T}, c::Union{Integer,Rational,BigInt,Cyclotomic}) where {N,T} =
    a + c * one_element(N)
Base.:+(c::Union{Integer,Rational,BigInt,Cyclotomic}, a::TaftElement) = a + c
Base.:-(a::TaftElement{N,T}) where {N,T} = TaftElement{N,T}(Dict(k => -v for (k,v) in a.terms))
Base.:-(a::TaftElement, b::TaftElement) = a + (-b)
Base.:-(a::TaftElement, c::Union{Integer,Rational,BigInt,Cyclotomic}) = a + (-c)
Base.:-(c::Union{Integer,Rational,BigInt,Cyclotomic}, a::TaftElement) = c * one_element(typeof(a).parameters[1]) - a

function Base.:*(a::TaftElement{N,T}, b::TaftElement{N,S}) where {N,T,S}
    U = promote_type(T,S,Cyclotomic{N}); q = zeta(Int(N)); out = Dict{Tuple{Int,Int},U}()
    for ((r,s),x) in a.terms, ((u,v),y) in b.terms
        r+u < N || continue
        key = (r+u, mod(s+v,N)); c = convert(U,x)*convert(U,y)*convert(U,q^(s*u))
        out[key] = get(out,key,zero(U)) + c
    end
    TaftElement{N,U}(out)
end
Base.:*(a::TaftElement, c::Union{Integer,Rational,BigInt,Cyclotomic}) =
    TaftElement{typeof(a).parameters[1],promote_type(eltype(values(a.terms)),typeof(c))}(Dict(k => v*c for (k,v) in a.terms))
Base.:*(c::Union{Integer,Rational,BigInt,Cyclotomic}, a::TaftElement) = a*c

function Base.:^(a::TaftElement, n::Integer)
    n < 0 && throw(ArgumentError("negative powers require an invertible element; use K for group powers"))
    out = one(a); base = a; k = n
    while k > 0
        isodd(k) && (out = out * base)
        k >>= 1; k > 0 && (base = base * base)
    end
    out
end

struct TaftTensor{N,T}
    terms::Dict{Tuple{Tuple{Int,Int},Tuple{Int,Int}},T}
end
tensor_terms(x::TaftTensor) = copy(x.terms)
Base.iszero(x::TaftTensor) = isempty(x.terms)

function coproduct(x::TaftElement{N,T}) where {N,T}
    n = Int(N); q = zeta(n); U = promote_type(T,typeof(q)); d = Dict{Tuple{Tuple{Int,Int},Tuple{Int,Int}},U}()
    for ((r,s),c) in x.terms, k in 0:r
        left = (k,s); right = (r-k,mod(k+s,n))
        key = (left,right); val = convert(U,c)*convert(U,qbinomial(n,r,k))
        d[key] = get(d,key,zero(U)) + val
    end
    TaftTensor{N,U}(d)
end

function counit(x::TaftElement{N,T}) where {N,T}
    out = zero(T)
    for ((r,_),c) in x.terms; r == 0 && (out += c); end
    out
end

function antipode(x::TaftElement{N,T}) where {N,T}
    n = Int(N); q = zeta(n); U = promote_type(T,typeof(q)); d = Dict{Tuple{Int,Int},U}()
    for ((r,s),c) in x.terms
        key = (r,mod(-r-s,n)); val = convert(U,c) * ((-one(q))^r) * q^(-r*(r-1)÷2-r*s)
        d[key] = get(d,key,zero(U)) + val
    end
    TaftElement{N,U}(d)
end

# ---------------------------------------------------------------------------
# Dual monomials and coordinate dual basis
# ---------------------------------------------------------------------------

struct DualMonomial{N,T}
    l::Int
    t::Int
end
struct DualBasis{N,T}
    r::Int
    s::Int
end
function dual_monomial(N::Integer, l::Integer, t::Integer)
    n = Int(N); n >= 2 && 0 <= l < n && 0 <= t < n ||
        throw(ArgumentError("dual exponents must lie in 0:N-1 for N ≥ 2"))
    DualMonomial{n,Cyclotomic{n}}(Int(l),Int(t))
end
function dual_basis(N::Integer, r::Integer, s::Integer)
    n = Int(N); n >= 2 && 0 <= r < n && 0 <= s < n ||
        throw(ArgumentError("basis exponents must lie in 0:N-1 for N ≥ 2"))
    DualBasis{n,Cyclotomic{n}}(Int(r),Int(s))
end

"""Fourier expansion of the coordinate dual `δ_(r,s)` in the `a^l b^r`
dual-monomial basis under the convention documented at module level."""
function dual_basis_expansion(f::DualBasis{N}) where {N}
    n = Int(N); 0 <= f.r < n && 0 <= f.s < n ||
        throw(ArgumentError("basis exponents must lie in 0:N-1"))
    q = zeta(n); scale = one(q) / (n * qfactorial(n,f.r))
    Dict((l,f.r) => scale * q^(-l*f.s) for l in 0:n-1)
end

function pairing(f::DualMonomial{N}, x::TaftElement{N,T}) where {N,T}
    n = Int(N); 0 <= f.l < n && 0 <= f.t < n || throw(ArgumentError("dual exponents must lie in 0:N-1"))
    q = zeta(n); out = zero(q); scale = qfactorial(n,f.t)
    for ((r,s),c) in x.terms
        r == f.t && (out += convert(typeof(q),c) * scale * q^(f.l*s))
    end
    out
end
function pairing(f::DualBasis{N}, x::TaftElement{N,T}) where {N,T}
    n = Int(N); 0 <= f.r < n && 0 <= f.s < n || throw(ArgumentError("basis exponents must lie in 0:N-1"))
    get(x.terms,(f.r,f.s),zero(T))
end
pairing(f::DualMonomial, r::Integer, s::Integer) = pairing(f, monomial(typeof(f).parameters[1],r,s))
pairing(f::DualBasis, r::Integer, s::Integer) = (f.r == r && f.s == s ? 1 : 0)
pairing(N::Integer, l::Integer, t::Integer, r::Integer, s::Integer) =
    pairing(dual_monomial(N,l,t), monomial(N,r,s))

end # module TaftExact
