// References for TVB-31.11-31.16: lifetimes, recombination and interface traps (Unit 55, TVR-3 step 5; DECISIONS
// 15.3, 15.9, 15.15).
//   - Scharfetter lifetime tau = tau0 / (1 + N / N_ref): Gritsch's thesis sec. 4.7 prints the general relation
//     tau_min + (tau_max - tau_min) / (1 + (N/N_ref)^gamma) (Fossum and Lee 1982); NiTCAD's is its tau_min = 0,
//     gamma = 1 case (PB for the form via Gritsch and Genius's NSRH use; values L).
//   - Shockley-Read-Hall, derived here (PD) from the capture-emission kinetics of one level: electrons captured at
//     c_n n (1 - f) N_t and emitted at c_n n1 f N_t, holes captured at c_p p f N_t and emitted at c_p p1 (1 - f) N_t,
//     with n1 p1 the equilibrium product E (detailed balance). Steady state (equal net electron and hole capture)
//     gives f = (c_n n + c_p p1) / (c_n (n + n1) + c_p (p + p1)) and U = N_t c_n c_p (n p - E) / (c_n (n + n1) +
//     c_p (p + p1)). With tau_n = 1 / (c_n N_t), tau_p = 1 / (c_p N_t) and n1 = p1 = n_ie (a mid-gap level),
//     R = (n p - E) / (tau_p (n + n_ie) + tau_n (p + n_ie)). Second method: the root of the occupancy balance by
//     Brent (TC5), in the L0 test. The chain through E: E's partials are inputs (dE/dn, dE/dp).
//   - Auger and radiative, derived here (PD) from detailed balance: a rate proportional to the process's carriers,
//     less its reverse at the equilibrium product, R_Auger = (Cn n + Cp p)(n p - E), R_rad = B (n p - E).
//   - The equilibrium occupancy f = 1 / (1 + e^x) and 1 - f, df/d eta = f (1 - f); the occupied density of a
//     uniform trap band, integral of D f((E - E_F) / kT) dE, by adaptive Gauss-Kronrod (TC3) as 15.9 fixes for
//     31.16, with its closed form D kT [ln(1 + e^(-u_lo)) - ln(1 + e^(-u_hi))] for the condition sum and as the L0
//     check of the quadrature.
#pragma once

#include <array>
#include <cmath>
#include <complex>

#include "NiTCAD/verify/reference/double_double.hpp"
#include "NiTCAD/verify/reference/gauss_kronrod.hpp"
#include "NiTCAD/verify/reference/model_band.hpp"
#include "NiTCAD/verify/reference/special_functions.hpp"

namespace NiTCAD::verify::reference::model_recombination {

template <class S>
[[nodiscard]] S scharfetter_lifetime(const S& tau0, const S& N_ref, const S& N) {
    return tau0 / (S(1.0) + N / N_ref);
}

// x = {n, p, E, dE/dn, dE/dp, n_ie, tau_n, tau_p}: {R, dR/dn, dR/dp}. The partials are the quotient rule expanded
// so that its two large terms cancel symbolically: (p - E_n) D - (n p - E) tau_p = tau_n p (p + n_ie) + tau_p (p n_ie
// + E) - E_n D, and the mirror for p. The bare quotient rule loses everything when one capture term dominates D
// (n = 1e-10, p = 1e36: a relative difference of 1e-26 even in double-double).
template <class S>
[[nodiscard]] std::array<S, 3> srh(const std::array<S, 8>& x) {
    const S& n = x[0];
    const S& p = x[1];
    const S& E = x[2];
    const S& ni = x[5];
    const S& tn = x[6];
    const S& tp = x[7];
    const S D = tp * (n + ni) + tn * (p + ni);
    const S D2 = D * D;
    return {(n * p - E) / D, (tn * p * (p + ni) + tp * (p * ni + E) - x[3] * D) / D2,
            (tp * n * (n + ni) + tn * (n * ni + E) - x[4] * D) / D2};
}

// x = {n, p, E, dE/dn, dE/dp, Cn, Cp}.
template <class S>
[[nodiscard]] std::array<S, 3> auger(const std::array<S, 7>& x) {
    const S excess = x[0] * x[1] - x[2];
    const S C = x[5] * x[0] + x[6] * x[1];
    return {C * excess, x[5] * excess + C * (x[1] - x[3]), x[6] * excess + C * (x[0] - x[4])};
}

// x = {n, p, E, dE/dn, dE/dp, B}.
template <class S>
[[nodiscard]] std::array<S, 3> radiative(const std::array<S, 6>& x) {
    return {x[5] * (x[0] * x[1] - x[2]), x[5] * (x[1] - x[3]), x[5] * (x[0] - x[4])};
}

// x = {n, p, n1, dn1/dn, p1, dp1/dp, cn, cp}: {f, 1 - f, df/dn, df/dp, U/N_t, dU/dn, dU/dp} / N_t. With occ = c_n n +
// c_p p1, emp = c_n n1 + c_p p and D = occ + emp, the quotient rule expanded (D - occ = emp symbolically) gives
// df/dn = c_n (emp - n1' occ) / D^2, df/dp = -c_p (occ - p1' emp) / D^2, and dU/dn = c_p (p + p1) df/dn,
// dU/dp = -c_n (n + n1) df/dp (the hole and electron capture forms of U). The bare quotient rule cancels to
// nothing at f = 1 - 1e-18.
template <class S>
[[nodiscard]] std::array<S, 7> trap_kinetics(const std::array<S, 8>& x) {
    const S& n = x[0];
    const S& p = x[1];
    const S& n1 = x[2];
    const S& p1 = x[4];
    const S& cn = x[6];
    const S& cp = x[7];
    const S occ = cn * n + cp * p1;
    const S emp = cn * n1 + cp * p;
    const S D = occ + emp;
    const S D2 = D * D;
    const S df_dn = cn * (emp - x[3] * occ) / D2;
    const S df_dp = -(cp * (occ - x[5] * emp)) / D2;
    return {occ / D,
            emp / D,
            df_dn,
            df_dp,
            cn * cp * (n * p - n1 * p1) / D,
            cp * (p + p1) * df_dn,
            -(cn * (n + n1)) * df_dp};
}

// f = 1 / (1 + e^x), 1 - f and f (1 - f), each from the exponential that cannot overflow (real part decides).
template <class S>
[[nodiscard]] std::array<S, 3> fermi_occupancy(const S& x) {
    S f, g;
    if (model_band::detail::re(x) > 0.0) {
        const S e = exp(-x);
        f = e / (S(1.0) + e);
        g = S(1.0) / (S(1.0) + e);
    } else {
        const S e = exp(x);
        f = S(1.0) / (S(1.0) + e);
        g = e / (S(1.0) + e);
    }
    return {f, g, f * g};
}

// The occupied (k = 0) or empty (k = 1) density of a band, x = {D_it, E_lo, E_hi, E_F, T}, closed form:
// occupied D kT [ln(1 + e^(-u_lo)) - ln(1 + e^(-u_hi))], empty D kT [ln(1 + e^(u_hi)) - ln(1 + e^(u_lo))],
// u = (E - E_F) / kT.
[[nodiscard]] inline DD band_density_closed(const std::array<DD, 5>& x, int k) {
    const DD kT = model_band::detail::k_eV_dd() * x[4];
    const DD ulo = (x[1] - x[3]) / kT;
    const DD uhi = (x[2] - x[3]) / kT;
    if (k == 0) return x[0] * kT * (log1p_exp(-ulo) - log1p_exp(-uhi));
    return x[0] * kT * (log1p_exp(uhi) - log1p_exp(ulo));
}

// The same by adaptive Gauss-Kronrod of D f (TC3), with its error estimate.
struct BandIntegral {
    DD value;
    double abs_error;
    bool ok;
};
[[nodiscard]] inline BandIntegral band_density_quadrature(const std::array<DD, 5>& x, int k) {
    const DD kT = model_band::detail::k_eV_dd() * x[4];
    auto f = [&](const DD& E) { return fermi_occupancy((E - x[3]) / kT)[k]; };
    const AdaptiveResult r = qag15(f, x[1], x[2], 0.0, 1e-28);
    return {x[0] * r.value, std::abs(x[0].hi) * r.abs_error, r.status == QuadratureStatus::ok};
}

}  // namespace NiTCAD::verify::reference::model_recombination
