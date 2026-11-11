// References for TVB-31.26-31.29: thermionic emission velocity, MLDA, the lattice conductivity and the thermopower
// factor (Unit 55, TVR-3 step 8; DECISIONS 15.3, 15.9, 15.18), in double-double.
//   - Emission (PD, derived here): with N = 2 (2 pi m k T / h^2)^(3/2), v = sqrt(kT / (2 pi m)) = (kT/h)(2/N)^(1/3);
//     with a Richardson constant, v = A T^2 / (q N). Through N(T): dv/dT = v (1/T - s/3), or v (2/T - s), s = d ln N
//     / dT (s = 3 / (2T) on the power law, so v / (2T) for both). Second method (L0): complex step on v(T).
//   - MLDA (PB form, Paasch-Ubensee / Hansch via Ieong and Entner): lambda = hbar / sqrt(2 m m0 k T) [cm]; f(z) = 1 -
//     exp(-z^2 / lambda^2); the mean over [0, a] = 1 - sqrt(pi) / (2x) erf(x), x = a / lambda (PD), by its series
//     sum_{k>=1} (-1)^(k+1) x^(2k) / (k! (2k + 1)) below x = 1/4. Second method (L0): tanh-sinh of f over [0, a].
//   - Conductivity: kappa = kappa_300 (T / 300)^(-a), d kappa / dT = -a kappa / T (the documented power law).
//   - Fermi orders (DLMF 25.12): F_-1 = sigma(x), F_0 = ln(1 + e^x), F_1 = -Li2(-e^x), reflected above 0 by Li2's
//     inversion (25.12.4): F_1(x) = x^2 / 2 + pi^2 / 6 + Li2(-e^-x). The thermopower factor h = 2 F_1 / F_0 - 3/2 and
//     dh/dx = 2 - 2 F_1 F_-1 / F_0^2; at x <= -1 from the series ratios S_j = sum (-1)^(k+1) z^(k-1) / k^(j+1), z =
//     e^x, so that e^x cancels. Where e^x is below double-double's domain (x < -600) it is formed as e^(x + 512 ln 2)
//     scaled by 2^-512 (the F_j are then e^x to double-double accuracy). Second method (L0): F_1 and F_0 differences
//     against tanh-sinh of F_0 and F_-1; the series against the closed forms at x = -1.
#pragma once

#include <array>
#include <cmath>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/verify/reference/double_double.hpp"
#include "NiTCAD/verify/reference/special_functions.hpp"

namespace NiTCAD::verify::reference::model_thermal {

namespace detail {
[[nodiscard]] inline DD cbrt(const DD& x) { return exp(log(x) / DD(3.0)); }
[[nodiscard]] inline DD planck() { return DD(2.0) * pi() * DD(base::hbar_J_s); }
}  // namespace detail

// x = {N [cm^-3], T}: v [cm/s] = (kT / h)(2 / N)^(1/3).
[[nodiscard]] inline DD emission_dos(const std::array<DD, 2>& x) {
    return DD(base::k_B_J_per_K) * x[1] / detail::planck() * detail::cbrt(DD(2.0) / x[0]);
}

// x = {A [A/(cm^2 K^2)], T, N}: v = A T^2 / (q N).
[[nodiscard]] inline DD emission_richardson(const std::array<DD, 3>& x) {
    return x[0] * x[1] * x[1] / (DD(base::q_C) * x[2]);
}

// The velocity and dv/dT through N(T) = N, s = d ln N / dT: the density-of-states form (A = 0) or Richardson's.
[[nodiscard]] inline std::array<DD, 2> emission(const DD& A, const DD& T, const DD& N, const DD& s) {
    if (A.hi > 0.0) {
        const DD v = emission_richardson({A, T, N});
        return {v, v * (DD(2.0) / T - s)};
    }
    const DD v = emission_dos({N, T});
    return {v, v * (DD(1.0) / T - s / DD(3.0))};
}

// x = {m [m0], T}: lambda [cm].
[[nodiscard]] inline DD mlda_wavelength(const std::array<DD, 2>& x) {
    return DD(base::hbar_J_s) / sqrt(DD(2.0) * x[0] * DD(base::m0_kg) * DD(base::k_B_J_per_K) * x[1]) * DD(100.0);
}

// x = {z, lambda}.
[[nodiscard]] inline DD mlda_factor(const std::array<DD, 2>& x) {
    const DD r = x[0] / x[1];
    return -expm1(-(r * r));
}

// x = {a, lambda}.
[[nodiscard]] inline DD mlda_mean(const std::array<DD, 2>& x) {
    const DD r = x[0] / x[1];
    if (r.hi < 0.25) {
        const DD r2 = r * r;
        DD term(1.0), sum(0.0);  // term = (-1)^(k+1) r^(2k) / k!
        for (int k = 1; k < 60; ++k) {
            term = k == 1 ? r2 : -(term * r2) / DD(static_cast<double>(k));
            const DD add = term / DD(2.0 * k + 1.0);
            sum = sum + add;
            if (std::abs(add.hi) < 1e-34 * std::abs(sum.hi)) break;
        }
        return sum;
    }
    return DD(1.0) - sqrt(pi()) / (DD(2.0) * r) * erf(r);
}

// x = {kappa_300, a, T}: {kappa, d kappa / dT}.
[[nodiscard]] inline std::array<DD, 2> conductivity(const std::array<DD, 3>& x) {
    const DD k = x[0] * exp(-(x[1] * log(x[2] / DD(300.0))));
    return {k, -(x[1] * k) / x[2]};
}

// ---- Fermi orders and the thermopower factor ----

struct FermiOrders {
    DD minus_one, zero, one;
};

namespace detail {
// e^x as {e^(x + 512 ln 2), true} where e^x is below double-double's domain, else {e^x, false}.
struct Exp {
    DD value;
    bool scaled;
};
[[nodiscard]] inline Exp exp_scaled(const DD& x) {
    if (x.hi < -600.0) return {exp(x + DD(512.0) * log(DD(2.0))), true};
    return {exp(x), false};
}
[[nodiscard]] inline DD unscale(const DD& v) { return DD(std::ldexp(v.hi, -512), std::ldexp(v.lo, -512)); }
// The three series S_-1, S_0, S_1 at z in (0, e^-1].
[[nodiscard]] inline std::array<DD, 3> series(const DD& z) {
    DD p(1.0), s0(0.0), s1(0.0);  // p = (-1)^(k+1) z^(k-1)
    for (int k = 1; k < 200; ++k) {
        const DD a0 = p / DD(static_cast<double>(k)), a1 = a0 / DD(static_cast<double>(k));
        s0 = s0 + a0;
        s1 = s1 + a1;
        if (std::abs(a0.hi) < 1e-34) break;
        p = -(p * z);
    }
    return {DD(1.0) / (DD(1.0) + z), s0, s1};
}
}  // namespace detail

[[nodiscard]] inline FermiOrders fermi_orders(const DD& x) {
    if (x.hi <= -1.0) {
        const detail::Exp e = detail::exp_scaled(x);
        if (e.scaled) {  // z < 1e-260: S_j = 1 to double-double accuracy
            const DD z = detail::unscale(e.value);
            return {z, z, z};
        }
        const auto s = detail::series(e.value);
        return {e.value * s[0], e.value * s[1], e.value * s[2]};
    }
    if (x.hi <= 0.0) {
        const DD z = exp(x);
        return {z / (DD(1.0) + z), log1p_exp(x), -li2(-z)};
    }
    const DD w = exp(-x);
    const DD p2 = pi() * pi() / DD(6.0);
    const DD tail = w.hi < 1e-260 ? -w : li2(-w);  // Li2(-w) = -w + w^2 / 4 - ...
    return {DD(1.0) / (DD(1.0) + w), x + log1p_exp(-x), x * x / DD(2.0) + p2 + tail};
}

namespace detail {
// The coefficients of S_0^2 - S_1 S_-1 = sum_{n>=1} c_n (-z)^n: c_n = sum_{i=0}^{n} [1 / ((i + 1)(n - i + 1)) - 1/(i +
// 1)^2] (c_0 = 0, c_1 = -1/4), so that dh/dx = 2 (S_0^2 - S_1 S_-1) / S_0^2 is summed without the cancellation of 2 - 2
// S_1 S_-1 / S_0^2 (dh/dx -> e^x / 2 as x -> -inf).
inline constexpr int dh_terms = 200;
[[nodiscard]] inline const std::array<DD, dh_terms>& dh_coefficients() {
    static const std::array<DD, dh_terms> c = [] {
        std::array<DD, dh_terms> out{};
        for (int n = 1; n < dh_terms; ++n) {
            DD s(0.0);
            for (int i = 0; i <= n; ++i) {
                const DD a(i + 1.0), b(n - i + 1.0);
                s = s + DD(1.0) / (a * b) - DD(1.0) / (a * a);
            }
            out[static_cast<std::size_t>(n)] = s;
        }
        return out;
    }();
    return c;
}
// S_0^2 - S_1 S_-1 at z in (0, e^-1].
[[nodiscard]] inline DD dh_numerator(const DD& z) {
    const auto& c = dh_coefficients();
    DD p = -z, sum(0.0);  // p = (-z)^n
    for (int n = 1; n < dh_terms; ++n) {
        const DD add = c[static_cast<std::size_t>(n)] * p;
        sum = sum + add;
        if (std::abs(add.hi) < 1e-34 * std::abs(sum.hi)) break;
        p = -(p * z);
    }
    return sum;
}
}  // namespace detail

// x = {x}: {h, dh/dx}.
[[nodiscard]] inline std::array<DD, 2> thermal_diffusion(const std::array<DD, 1>& xa) {
    const DD& x = xa[0];
    DD ratio, product;
    if (x.hi <= -1.0) {
        const detail::Exp e = detail::exp_scaled(x);
        if (e.scaled) {  // S_j = 1, S_0^2 - S_1 S_-1 = z / 4 to double-double accuracy
            return {DD(0.5), detail::unscale(e.value / DD(2.0))};
        }
        const auto s = detail::series(e.value);
        return {DD(2.0) * s[2] / s[1] - DD(1.5), DD(2.0) * detail::dh_numerator(e.value) / (s[1] * s[1])};
    } else {
        const FermiOrders f = fermi_orders(x);
        ratio = f.one / f.zero;
        product = ratio * f.minus_one / f.zero;
    }
    return {DD(2.0) * ratio - DD(1.5), DD(2.0) - DD(2.0) * product};
}

}  // namespace NiTCAD::verify::reference::model_thermal
