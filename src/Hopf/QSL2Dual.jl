module QSL2DualMod

export qsl2_qnum, qsl2_qfactorial,
       qsl2_root_pairing_iterative, qsl2_root_dual_coefficient,
       qsl2_dual_pairing,
       qsl2_old_formula_leaks, qsl2_corrected_generator_pairings

_choose2(n::Integer) = n * (n - 1) ÷ 2

"""Symmetric q-number [n]_q, with the classical q=1 value filled by continuity."""
function qsl2_qnum(q::Number, n::Integer)
    n < 0 && throw(ArgumentError("n must be nonnegative"))
    n == 0 && return zero(q)
    q == one(q) && return n * one(q)
    (q^n - q^(-n)) / (q - q^(-1))
end

"""Symmetric q-factorial [n]_q!."""
function qsl2_qfactorial(q::Number, n::Integer)
    n < 0 && throw(ArgumentError("n must be nonnegative"))
    prod((qsl2_qnum(q, r) for r in 1:n); init=one(q))
end

"""
Root-sector pairing obtained from the coregular recurrence, not from the closed
normalization formula. It computes the coefficient of
∂_E^i ∂_F^j(x^i y^j) at the counit for the R/K-right convention.
"""
function qsl2_root_pairing_iterative(q::Number, i::Integer, j::Integer)
    (i < 0 || j < 0) && throw(ArgumentError("exponents must be nonnegative"))

    coefficient = one(q)

    # ∂_F^j(y^j) = q^(binom(j,2)) [j]_q! d^(2j)
    for m in 1:j
        coefficient *= q^(m - 1) * qsl2_qnum(q, m)
    end

    # After the F steps the d exponent is N=2j. Each E step lowers N by 2.
    m = i
    N = 2j
    while m > 0
        coefficient *= q^(2 - m - N) * qsl2_qnum(q, m)
        m -= 1
        N -= 2
    end

    coefficient
end

"""Coefficient multiplying x^i y^j in the corrected R-convention dual basis."""
function qsl2_root_dual_coefficient(q::Number, i::Integer, j::Integer)
    exponent = 2i * j - i * (i + 1) ÷ 2 - _choose2(j)
    q^exponent / (qsl2_qfactorial(q, i) * qsl2_qfactorial(q, j))
end

"""
Pair the corrected dual candidate for E^i F^j H^k against E^r F^s H^t.
The Cartan polynomial is (z-i*hbar)^k/k!, whose t-th jet at z=i*hbar is
δ_{kt}. The root scalar is computed independently by the recurrence above.
"""
function qsl2_dual_pairing(q::Number,
                           i::Integer, j::Integer, k::Integer,
                           r::Integer, s::Integer, t::Integer)
    any(x -> x < 0, (i, j, k, r, s, t)) &&
        throw(ArgumentError("exponents must be nonnegative"))
    (i == r && j == s && k == t) || return zero(q)

    root = qsl2_root_pairing_iterative(q, i, j)
    coefficient = qsl2_root_dual_coefficient(q, i, j) / factorial(k)
    cartan_jet = factorial(k)
    coefficient * root * cartan_jet
end

"""Known nonzero leaks of the withdrawn formula."""
qsl2_old_formula_leaks(q::Number) = (
    E_pairing = q,                    # <b d^-1, E>
    EH_to_E_hbar_coefficient = q,    # hbar coefficient of <x z, E>
)

"""Generator and EH checks for the corrected q^-1*x and x*(z-hbar) factors."""
function qsl2_corrected_generator_pairings(q::Number)
    root_E = qsl2_root_pairing_iterative(q, 1, 0)
    coefficient_E = qsl2_root_dual_coefficient(q, 1, 0)
    (
        E_on_E = coefficient_E * root_E,
        EH_on_E = zero(q),
        EH_on_EH = coefficient_E * root_E,
        F_on_F = qsl2_root_dual_coefficient(q, 0, 1) *
                 qsl2_root_pairing_iterative(q, 0, 1),
        H_on_H = one(q),
    )
end

end # module
