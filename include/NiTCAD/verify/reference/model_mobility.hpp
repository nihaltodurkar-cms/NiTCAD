// References for TVB-31.06-31.10: mobility (Unit 55, TVR-3 step 4; DECISIONS 15.3, 15.9). From the sources:
//   - Caughey-Thomas: mu = mu_min + (mu_max - mu_min) / (1 + (N / N_ref)^alpha) (Caughey and Thomas, Proc. IEEE 55,
//     2192 (1967), the form as Xyce's and Genius's analytic model reproduce it; PB), with mu_max scaled by
//     (T / 300)^gamma as the header documents (an Arora exponent on the CT mu_max: F3);
//   - Masetti: mu = mu_min1 e^(-P_c / N) + (mu_max - mu_min2) / (1 + (N / C_r)^alpha) - mu_1 / (1 + (C_s / N)^beta)
//     (Masetti, Severi, Solmi 1983 eq (4) as Allpix Squared prints it; PB); the published form is at 300 K, the
//     header scales mu_max by (T / 300)^gamma (F6, tag C);
//   - Canali: mu = mu0 / (1 + (mu0 E / v_sat)^beta)^(1/beta) (Canali et al. 1975 as Allpix and Ayalew's thesis; PB);
//   - Lombardi: 1/mu = 1/mu_b + 1/mu_ac + 1/mu_sr, mu_ac = B / E + C N^lambda / (T E^(1/3)), mu_sr = delta / E^2
//     (Lombardi et al. 1988 as Xyce and Genius reproduce it; PB);
//   - the energy-balance mobility (Unit 52, ED4 and the owner's beta = 1 ruling): the uniform-field balance
//     x^2 m(x) = a, m = (1 + x^beta)^(-1/beta), a = mu0 (3/2)(k/q)(T - T_L) / (tau_w v_sat^2); mu = mu0 m, E_eff =
//     x v_sat / mu0; for beta != 2 the regularization below a0 = 1e-8 (cubic c2 a^2 + c3 a^3 matching x and dx/da at
//     a0; 0 for a <= 0), for beta = 2 the closed form y = x^2 = a (a + sqrt(a^2 + 4)) / 2 for the mobility at every a.
//     Derived here (PD); the root is the second method's (Brent, TC5, in the L0 test).
// Templates on the scalar (DD values, complex double for complex step); branches decide on real parts.
#pragma once

#include <array>
#include <complex>
#include <stdexcept>
#include <type_traits>

#include "NiTCAD/verify/reference/codata2018.hpp"
#include "NiTCAD/verify/reference/double_double.hpp"
#include "NiTCAD/verify/reference/model_band.hpp"

namespace NiTCAD::verify::reference::model_mobility {

using model_band::detail::re;

namespace detail {

template <class S>
[[nodiscard]] S spow(const S& x, const S& y) {
    if constexpr (std::is_same_v<S, DD>) {
        return pow(x, y);
    } else {
        return std::pow(x, y);
    }
}

}  // namespace detail

// ---- Caughey-Thomas (31.06) ----

template <class S>
[[nodiscard]] S ct_mobility(const S& mu_min, const S& mu_max, const S& N_ref, const S& alpha, const S& gamma,
                            const S& N, const S& T) {
    const S mx = mu_max * detail::spow(T / S(300.0), gamma);
    const S r = re(N) == 0.0 ? S(0.0) : detail::spow(N / N_ref, alpha);
    return mu_min + (mx - mu_min) / (S(1.0) + r);
}

template <class S>
[[nodiscard]] S ct_mobility_slope(const S& mu_max, const S& N_ref, const S& alpha, const S& gamma, const S& N,
                                  const S& T) {
    const S mx = mu_max * detail::spow(T / S(300.0), gamma);
    const S r = re(N) == 0.0 ? S(0.0) : detail::spow(N / N_ref, alpha);
    return gamma * mx / T / (S(1.0) + r);
}

// ---- Masetti (31.07). x = {mu_min1, mu_min2, mu_max, mu_1, P_c, C_r, C_s, alpha, beta, gamma} ----

template <class S>
[[nodiscard]] S masetti_mobility(const std::array<S, 10>& p, const S& N, const S& T) {
    const S mx = p[2] * detail::spow(T / S(300.0), p[9]);
    if (re(N) == 0.0) return (re(p[4]) == 0.0 ? p[0] : S(0.0)) + mx - p[1];
    const S first = re(p[4]) == 0.0 ? p[0] : p[0] * exp(-p[4] / N);
    return first + (mx - p[1]) / (S(1.0) + detail::spow(N / p[5], p[7])) -
           p[3] / (S(1.0) + detail::spow(p[6] / N, p[8]));
}

template <class S>
[[nodiscard]] S masetti_mobility_slope(const std::array<S, 10>& p, const S& N, const S& T) {
    const S mx = p[2] * detail::spow(T / S(300.0), p[9]);
    const S r = re(N) == 0.0 ? S(0.0) : detail::spow(N / p[5], p[7]);
    return p[9] * mx / T / (S(1.0) + r);
}

// ---- Canali (31.08) ----

template <class S>
[[nodiscard]] S canali_mobility(const S& mu0, const S& v_sat, const S& beta, const S& E) {
    if (re(E) == 0.0) return mu0;
    const S x = mu0 * E / v_sat;
    return mu0 / detail::spow(S(1.0) + detail::spow(x, beta), S(1.0) / beta);
}

// d mu / dE; the right-hand derivative at E = 0 (0 for beta > 1, -mu0^2 / v_sat for beta = 1).
template <class S>
[[nodiscard]] S canali_mobility_slope(const S& mu0, const S& v_sat, const S& beta, const S& E) {
    if (re(E) == 0.0) return re(beta) == 1.0 ? -(mu0 * mu0) / v_sat : S(0.0);
    const S x = mu0 * E / v_sat;
    const S xb = detail::spow(x, beta);
    return -(mu0 * detail::spow(S(1.0) + xb, -S(1.0) / beta - S(1.0)) * xb / x * (mu0 / v_sat));
}

// ---- Lombardi (31.09). p = {B, C, lambda, delta} ----

template <class S>
[[nodiscard]] S lombardi_phonon(const std::array<S, 4>& p, const S& N, const S& T, const S& E) {
    if (re(N) == 0.0) return S(0.0);
    return p[1] * detail::spow(N, p[2]) / (T * detail::spow(E, S(1.0) / S(3.0)));
}

template <class S>
[[nodiscard]] S lombardi_mobility(const std::array<S, 4>& p, const S& mu_b, const S& N, const S& T, const S& E) {
    const S mu_ac = p[0] / E + lombardi_phonon(p, N, T, E);
    const S mu_sr = p[3] / (E * E);
    return S(1.0) / (S(1.0) / mu_b + S(1.0) / mu_ac + S(1.0) / mu_sr);
}

// d mu / dE, d mu / dT and d mu / d mu_b.
template <class S>
[[nodiscard]] std::array<S, 3> lombardi_partials(const std::array<S, 4>& p, const S& mu_b, const S& N, const S& T,
                                                const S& E) {
    const S ph = lombardi_phonon(p, N, T, E);
    const S mu_ac = p[0] / E + ph;
    const S mu_sr = p[3] / (E * E);
    const S mu = S(1.0) / (S(1.0) / mu_b + S(1.0) / mu_ac + S(1.0) / mu_sr);
    const S wa = (mu / mu_ac) * (mu / mu_ac);
    const S ws = (mu / mu_sr) * (mu / mu_sr);
    const S dac_dE = -p[0] / (E * E) - ph / (S(3.0) * E);
    const S dsr_dE = -S(2.0) * mu_sr / E;
    return {wa * dac_dE + ws * dsr_dE, wa * (-ph / T), (mu / mu_b) * (mu / mu_b)};
}

// ---- carrier-temperature mobility (31.10) ----

inline constexpr double regularization_a0 = 1e-8;

// a [1] from mu0 [cm^2/Vs], rise [K], tau_w [s], v_sat [cm/s]; k / q exact (CODATA 2018).
template <class S>
[[nodiscard]] S balance_a(const S& mu0, const S& rise, const S& tau_w, const S& v_sat) {
    return S(1.5) * model_band::detail::k_eV<S>() * rise * mu0 / (tau_w * v_sat * v_sat);
}

// m(x) and the share s = x^beta / (1 + x^beta), x > 0.
template <class S>
[[nodiscard]] S balance_m(const S& x, const S& beta) {
    return detail::spow(S(1.0) + detail::spow(x, beta), -S(1.0) / beta);
}

// The exact root x(a) > 0 of x^2 m(x) = a (a > 0) and dx/da = 1 / (x m (2 - s)), in double-double: closed forms for
// beta = 1 and 2, otherwise Newton on ln x^2 - ln(1 + x^beta) / beta = ln a from the asymptote's estimate.
struct Root {
    DD x, d_a;
};
[[nodiscard]] inline Root balance_root(const DD& a, double beta) {
    if (!(a.hi > 0.0)) throw std::domain_error("balance_root: a <= 0");
    DD x;
    if (beta == 1.0) {
        x = (a + sqrt(a * a + DD(4.0) * a)) / DD(2.0);
    } else if (beta == 2.0) {
        x = sqrt(a * (a + sqrt(a * a + DD(4.0))) / DD(2.0));
    } else {
        const DD b(beta);
        // A double Newton to convergence on the same equation, then double-double steps (quadratic: a few).
        double xd = a.hi < 1.0 ? std::sqrt(a.hi) : a.hi;
        for (int i = 0; i < 100; ++i) {
            const double xb = std::pow(xd, beta);
            const double g = 2.0 * std::log(xd) - std::log1p(xb) / beta - std::log(a.hi);
            const double step = g / (2.0 / xd - xb / (xd * (1.0 + xb)));
            xd -= step;
            if (std::abs(step) <= 1e-15 * xd) break;
        }
        x = DD(xd);
        for (int i = 0; i < 8; ++i) {
            const DD xb = pow(x, b);
            const DD g = DD(2.0) * log(x) - log(DD(1.0) + xb) / b - log(a);
            const DD gp = DD(2.0) / x - xb / (x * (DD(1.0) + xb));
            const DD step = g / gp;
            x = x - step;
            if (std::abs(step.hi) <= 0x1p-108 * std::abs(x.hi)) break;
        }
    }
    const DD b(beta);
    const DD xb = pow(x, b);
    const DD s = xb / (DD(1.0) + xb);
    const DD m = balance_m(x, b);
    return {x, DD(1.0) / (x * m * (DD(2.0) - s))};
}

// The regularized x(a) and dx/da of the model (the owner's ruling).
[[nodiscard]] inline Root regularized_root(const DD& a, double beta) {
    if (!(a.hi > 0.0)) return {DD(0.0), DD(0.0)};
    const DD a0(regularization_a0);
    if (!(a < a0)) return balance_root(a, beta);
    const Root e = balance_root(a0, beta);
    const DD c2 = (DD(3.0) * e.x - a0 * e.d_a) / (a0 * a0);
    const DD c3 = (a0 * e.d_a - DD(2.0) * e.x) / (a0 * a0 * a0);
    return {a * a * (c2 + c3 * a), a * (DD(2.0) * c2 + DD(3.0) * c3 * a)};
}

struct BalanceResult {
    DD mobility, d_mu0, d_rise, field, field_d_mu0, field_d_rise;
};

// x = {mu0, rise, tau_w, v_sat}.
[[nodiscard]] inline BalanceResult carrier_temperature_mobility(const std::array<DD, 4>& x, double beta) {
    const DD& mu0 = x[0];
    const DD& v_sat = x[3];
    const DD a_d_rise = balance_a(mu0, DD(1.0), x[2], v_sat);
    const DD a = a_d_rise * x[1];
    const Root r = regularized_root(a, beta);
    DD m(1.0), m_a(0.0);
    if (beta == 2.0) {
        const DD q = sqrt(a * a + DD(4.0));
        const DD y = a * (a + q) / DD(2.0);
        const DD y_a = ((a + q) + a * (DD(1.0) + a / q)) / DD(2.0);
        m = DD(1.0) / sqrt(DD(1.0) + y);
        m_a = -(m * y_a) / (DD(2.0) * (DD(1.0) + y));
    } else if (r.x.hi > 0.0) {
        const DD b(beta);
        const DD xb = pow(r.x, b);
        m = balance_m(r.x, b);
        m_a = -(m * xb / r.x) / (DD(1.0) + xb) * r.d_a;
    }
    return {mu0 * m,
            m + m_a * a,
            mu0 * m_a * a_d_rise,
            r.x * v_sat / mu0,
            v_sat * (r.d_a * a - r.x) / (mu0 * mu0),
            v_sat * r.d_a * a_d_rise / mu0};
}

}  // namespace NiTCAD::verify::reference::model_mobility
