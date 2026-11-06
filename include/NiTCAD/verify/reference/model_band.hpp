// References for TVB-31.01-31.05: band gap, densities of states, intrinsic density and band-gap narrowing (Unit 55,
// TVR-3 step 3; DECISIONS 15.3, 15.9). Written from the sources, not from NiTCAD's code (§2.2):
//   - Varshni: Eg(T) = Eg0 - alpha T^2 / (T + beta) (Y. Varshni, Physica 34, 149 (1967), as reproduced by
//     Palankovski's thesis eq 3.60 and Xyce's bandgapTbl.tex; PB);
//   - the power-law densities of states N(T) = N300 (T / 300)^(3/2) (Palankovski eqs 3.91-3.92 citing Sze; PB);
//   - Passler's gap as Couderc et al., J. Appl. Phys. 115, 093705 (2014) give it, eq (10), with chi = 2T / Theta:
//     Eg(T) = Eg0 - alpha Theta [ (1 - 3 Delta^2) / (e^(Theta/T) - 1) + (3 Delta^2 / 2) ( P^(1/6) - 1 ) ],
//     P = 1 + pi^2 chi^2 / (3 (1 + Delta^2)) + (3 Delta^2 - 1) chi^3 / 4 + 8 chi^4 / 3 + chi^6
//     (transcription in DECISIONS Unit 35, read from the open PDF; PA);
//   - Couderc's densities of states, eqs (4), (5), (18), (19): N = 4.83e15 M(T) T^(3/2) [cm^-3],
//     M = (m_d / m0)^(3/2) = c0 T^3 + c1 T^2 + c2 T + c3 (PA);
//   - n_i^2 = Nc Nv exp(-Eg / kT) (Couderc eq (1); PA);
//   - E_c - E_i = Eg / 2 + (kT / 2) ln(Nc / Nv), from n = p = n_i under Boltzmann statistics (PD);
//   - Slotboom's narrowing as published, dEg = E0 [ln(N / N0) + sqrt(ln^2(N / N0) + 1/2)] (Tidy3D, Genius; PB), and
//     NiTCAD's documented form with "- sqrt(1/2)" above N0 and 0 at or below (bandgap_narrowing.hpp; the departure
//     F1 of DECISIONS 15.4, checked under MI-Q4); n_ie = n_i exp(dEg / 2kT).
// k is CODATA 2018's exact k / q in eV/K (codata2018.hpp), not base's.
//
// Every function is a template on the scalar: DD for the reference values, std::complex<double> for complex-step
// derivatives (the L0 checks and the condition sums of 15.9). Branches decide on the real part only.
// The analytic derivatives (the *_slope functions) are those NiTCAD returns partials for (thermal.hpp).
#pragma once

#include <array>
#include <complex>
#include <type_traits>

#include "NiTCAD/verify/reference/codata2018.hpp"
#include "NiTCAD/verify/reference/double_double.hpp"

namespace NiTCAD::verify::reference::model_band {

namespace detail {

template <class S>
[[nodiscard]] double re(const S& x) {
    if constexpr (std::is_same_v<S, DD>) {
        return x.hi;
    } else {
        return x.real();
    }
}

// x^(num / den) with the exponent formed in the scalar (1/6 is not a double; DD(1.0 / 6.0) would err by 1e-17).
template <class S>
[[nodiscard]] S rpow(const S& x, double num, double den) {
    if constexpr (std::is_same_v<S, DD>) {
        return pow(x, DD(num) / DD(den));
    } else {
        return std::pow(x, num / den);
    }
}

template <class S>
[[nodiscard]] S expm1_s(const S& x) {
    if constexpr (std::is_same_v<S, DD>) {
        return expm1(x);
    } else {
        return std::exp(x) - 1.0;
    }
}

[[nodiscard]] inline DD k_eV_dd() {
    static const DD k = codata2018::to_dd(codata2018::boltzmann) / codata2018::to_dd(codata2018::elementary_charge);
    return k;
}

template <class S>
[[nodiscard]] S k_eV() {
    if constexpr (std::is_same_v<S, DD>) {
        return k_eV_dd();
    } else {
        return S(k_eV_dd().to_double());
    }
}

template <class S>
[[nodiscard]] S pi_s() {
    if constexpr (std::is_same_v<S, DD>) {
        return pi();
    } else {
        return S(pi().to_double());
    }
}

}  // namespace detail

// ---- Varshni and the power-law densities of states (legacy silicon) ----

template <class S>
[[nodiscard]] S varshni_gap(const S& Eg0, const S& alpha, const S& beta, const S& T) {
    return Eg0 - alpha * T * T / (T + beta);
}

template <class S>
[[nodiscard]] S varshni_gap_slope(const S& alpha, const S& beta, const S& T) {
    const S d = T + beta;
    return -(alpha * T * (T + S(2.0) * beta)) / (d * d);
}

template <class S>
[[nodiscard]] S power_dos(const S& N300, const S& T) {
    const S r = T / S(300.0);
    return N300 * r * sqrt(r);
}

template <class S>
[[nodiscard]] S power_dos_log_slope(const S& T) {
    return S(1.5) / T;
}

// ---- Passler (Couderc eq (10)) ----

template <class S>
[[nodiscard]] S passler_P(const S& Delta, const S& chi) {
    const S pi = detail::pi_s<S>();
    const S chi2 = chi * chi;
    const S D2 = Delta * Delta;
    return S(1.0) + pi * pi * chi2 / (S(3.0) * (S(1.0) + D2)) + (S(3.0) * D2 - S(1.0)) * chi2 * chi / S(4.0) +
           S(8.0) * chi2 * chi2 / S(3.0) + chi2 * chi2 * chi2;
}

template <class S>
[[nodiscard]] S passler_gap(const S& Eg0, const S& alpha, const S& Theta, const S& Delta, const S& T) {
    const S D2 = Delta * Delta;
    const S chi = S(2.0) * T / Theta;
    const S first = (S(1.0) - S(3.0) * D2) / detail::expm1_s(Theta / T);
    const S second = S(1.5) * D2 * (detail::rpow(passler_P(Delta, chi), 1.0, 6.0) - S(1.0));
    return Eg0 - alpha * Theta * (first + second);
}

template <class S>
[[nodiscard]] S passler_gap_slope(const S& alpha, const S& Theta, const S& Delta, const S& T) {
    const S pi = detail::pi_s<S>();
    const S D2 = Delta * Delta;
    const S chi = S(2.0) * T / Theta;
    const S chi2 = chi * chi;
    // d/dT of 1 / (e^x - 1), x = Theta / T, is e^x Theta / (T^2 (e^x - 1)^2) = Theta / (T^2 (e^x - 1)(1 - e^-x)):
    // written so that e^(2x) is never formed (x = 446 at 1 K).
    const S x = Theta / T;
    const S dfirst = (S(1.0) - S(3.0) * D2) * Theta / (T * T * detail::expm1_s(x) * -detail::expm1_s(-x));
    const S dP_dchi = S(2.0) * pi * pi * chi / (S(3.0) * (S(1.0) + D2)) +
                      S(3.0) * (S(3.0) * D2 - S(1.0)) * chi2 / S(4.0) + S(32.0) * chi2 * chi / S(3.0) +
                      S(6.0) * chi2 * chi2 * chi;
    const S dsecond = S(1.5) * D2 * detail::rpow(passler_P(Delta, chi), -5.0, 6.0) / S(6.0) * dP_dchi *
                      S(2.0) / Theta;
    return -(alpha * Theta * (dfirst + dsecond));
}

// ---- Couderc's densities of states (eqs (4), (5), (18), (19)) ----

template <class S>
[[nodiscard]] S mass_polynomial(const std::array<S, 4>& c, const S& T) {
    return c[0] * T * T * T + c[1] * T * T + c[2] * T + c[3];
}

template <class S>
[[nodiscard]] S mass_polynomial_slope(const std::array<S, 4>& c, const S& T) {
    return S(3.0) * c[0] * T * T + S(2.0) * c[1] * T + c[2];
}

template <class S>
[[nodiscard]] S couderc_dos(const std::array<S, 4>& c, const S& T) {
    return S(4.83e15) * mass_polynomial(c, T) * T * sqrt(T);
}

template <class S>
[[nodiscard]] S couderc_dos_log_slope(const std::array<S, 4>& c, const S& T) {
    return S(1.5) / T + mass_polynomial_slope(c, T) / mass_polynomial(c, T);
}

// ---- intrinsic quantities ----

// ln n_i, and n_i formed by one exponential of it (sqrt(Nc Nv) e^(-Eg / 2kT) would underflow in the second factor
// while n_i is still a normal double).
template <class S>
[[nodiscard]] S intrinsic_log_density(const S& Nc, const S& Nv, const S& Eg, const S& T) {
    return log(Nc * Nv) / S(2.0) - Eg / (S(2.0) * detail::k_eV<S>() * T);
}

template <class S>
[[nodiscard]] S intrinsic_density(const S& Nc, const S& Nv, const S& Eg, const S& T) {
    return exp(intrinsic_log_density(Nc, Nv, Eg, T));
}

// E_c - E_i [eV].
template <class S>
[[nodiscard]] S conduction_to_intrinsic(const S& Nc, const S& Nv, const S& Eg, const S& T) {
    return Eg / S(2.0) + detail::k_eV<S>() * T / S(2.0) * log(Nc / Nv);
}

// d ln n_i / dT from d ln Nc / dT, d ln Nv / dT, Eg, dEg / dT and the narrowing dEg_BGN (0 for n_i).
template <class S>
[[nodiscard]] S intrinsic_log_slope(const S& dlnNc, const S& dlnNv, const S& Eg, const S& dEg, const S& dEg_bgn,
                                    const S& T) {
    const S k = detail::k_eV<S>();
    return (dlnNc + dlnNv) / S(2.0) + (Eg - T * dEg) / (S(2.0) * k * T * T) - dEg_bgn / (S(2.0) * k * T * T);
}

// ---- band-gap narrowing ----

template <class S>
[[nodiscard]] S slotboom_published(const S& E0, const S& N0, const S& N) {
    const S L = log(N / N0);
    return E0 * (L + sqrt(L * L + S(0.5)));
}

template <class S>
[[nodiscard]] S slotboom_documented(const S& E0, const S& N0, const S& N) {
    if (!(detail::re(N) > detail::re(N0))) return S(0.0);
    const S L = log(N / N0);
    return E0 * (L + sqrt(L * L + S(0.5)) - sqrt(S(0.5)));
}

// n_ie = n_i e^(dEg / 2kT), from ln n_i.
template <class S>
[[nodiscard]] S effective_intrinsic(const S& ln_n_i, const S& dEg, const S& T) {
    return exp(ln_n_i + dEg / (S(2.0) * detail::k_eV<S>() * T));
}

}  // namespace NiTCAD::verify::reference::model_band
