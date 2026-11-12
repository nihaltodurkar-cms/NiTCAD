// References for TVB-31.20-31.22: impact ionization, Kane band-to-band generation and the direct-gap WKB rate (Unit
// 55, TVR-3 step 7; DECISIONS 15.3, 15.9, 15.17).
//   - van Overstraeten-de Man: alpha = gamma A exp(-gamma B / E), d alpha / dE = alpha gamma B / E^2, with the branch
//     (A, B) chosen by the caller; gamma(T) = tanh(h / 300) / tanh(h / T), h = hw / (2 k), and d gamma / dT =
//     tanh(h / 300) h / (T^2 sinh^2(h / T)) (the temperature factor is a tag-C choice, F4).
//   - Kane: G = A F^2 exp(-B / F), dG/dF = G (2 + B/F) / F.
//   - Direct-gap WKB (Esseni et al. 2017, sec. 2.1, eqs (8), (9), (11), (12); PA, read in step 1), SI units:
//     alpha(delta) = -u + 2 S, S = sqrt(u (delta - 1/2) + u^2/4 + 1/4), u = m0 / (2 m_r), formed without cancellation
//     at the band edges as 1 + alpha = 2 u delta / (S + S0), 1 - alpha = 2 u (1 - delta) / (S + S1) (S^2 - S0^2 =
//     u delta, S1^2 - S^2 = u (1 - delta)); cos(theta) = sqrt((1 + alpha)(1 - alpha)); kappa = sqrt(m_r E_g) c / hbar.
//     The antiderivatives in delta, with theta = atan2(alpha, c) and d delta = (alpha + u) / (2u) d alpha:
//     K = sqrt(m_r E_g) / (2 u hbar) (u (theta + alpha c) / 2 - c^3 / 3), Ginv = hbar / (2 u sqrt(m_r E_g)) (u theta -
//     c); second method: tanh-sinh quadrature of kappa and 1 / kappa (L0).
//   - A segment along which delta runs linearly from da to db over length L: I = (L / D) int f d delta over the part
//     inside the gap, D = |db - da|. Its partials are Leibniz's, d I / d v_max = (L / D^2)(D f(v_max) - int f), and
//     the mirror; for D < 1e-2 the bracket is integrated as int (f(v_max) - f) d delta (tanh-sinh; one
//     Kronrod-15 rule where the in-gap part is within 1e-3 of its distance to the nearer edge), so that it does not
//     cancel (the closed-form difference loses u_dd / D^2 relative).
//   - The path rate eq (11): G = |s| / (36 hbar) / I_ik (1 - e^-z) e^(-2 I_k), z = k_m^2 I_ik, and its partials;
//     d G / d I_ik = -pre e^(-2 I_k) (1 - e^-z (1 + z)) / I_ik^2, the bracket by its series for z < 1.
//   - eq (8), the uniform-field rate: q^2 F^2 sqrt(m_r) / (18 pi hbar^2 sqrt(E_g)) exp(-pi sqrt(m_r) E_g^1.5 / (2 q F
//     hbar)); second method (L0): eq (11) with the exact uniform-field integrals, in the limit k_m^2 I_ik >> 1.
// Constants: NiTCAD's base constants as double-doubles (the model's constants are its inputs here; base's against
// CODATA is TVR-1's).
#pragma once

#include <algorithm>
#include <array>
#include <cmath>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/verify/reference/double_double.hpp"
#include "NiTCAD/verify/reference/gauss_kronrod.hpp"
#include "NiTCAD/verify/reference/model_band.hpp"
#include "NiTCAD/verify/reference/special_functions.hpp"
#include "NiTCAD/verify/reference/tanh_sinh.hpp"

namespace NiTCAD::verify::reference::model_generation {

// ---- impact ionization and Kane ----

// x = {A, B, gamma, E}: {alpha, d alpha / dE}.
template <class S>
[[nodiscard]] std::array<S, 2> vodm(const std::array<S, 4>& x) {
    const S a = x[2] * x[0] * exp(-(x[2] * x[1] / x[3]));
    return {a, a * x[2] * x[1] / (x[3] * x[3])};
}

// x = {hw [eV], T}: {gamma, d gamma / dT}.
[[nodiscard]] inline std::array<DD, 2> vodm_gamma(const std::array<DD, 2>& x) {
    const DD h = x[0] / (DD(2.0) * model_band::detail::k_eV_dd());
    const DD y = h / x[1];
    const DD t300 = tanh(h / DD(300.0));
    const DD sh = (exp(y) - exp(-y)) / DD(2.0);
    return {t300 / tanh(y), t300 * h / (x[1] * x[1] * sh * sh)};
}

// x = {A, B, F}: {G, dG/dF}.
template <class S>
[[nodiscard]] std::array<S, 2> kane(const std::array<S, 3>& x) {
    const S G = x[0] * x[2] * x[2] * exp(-(x[1] / x[2]));
    return {G, G * (S(2.0) + x[1] / x[2]) / x[2]};
}

// ---- direct-gap WKB ----

struct Band {
    DD gap_J, reduced_mass_kg;
};

[[nodiscard]] inline DD hbar() { return DD(base::hbar_J_s); }
[[nodiscard]] inline DD u_of(const Band& b) { return DD(base::m0_kg) / (DD(2.0) * b.reduced_mass_kg); }
[[nodiscard]] inline DD scale_of(const Band& b) { return sqrt(b.reduced_mass_kg * b.gap_J); }

struct AlphaCos {
    DD alpha, c;
};

// d in [0, 1] with its complement w = 1 - d given by the caller (a caller that has 1 - d without cancellation, as
// TVB-24's second method near d = 1, passes it; the one-argument form below forms it as 1 - d).
[[nodiscard]] inline AlphaCos alpha_cos(const DD& d, const DD& w, const DD& u) {
    const DD S = sqrt(u * (d - DD(0.5)) + u * u / DD(4.0) + DD(0.25));
    const DD one_p = DD(2.0) * u * d / (S + (u - DD(1.0)) / DD(2.0));
    const DD one_m = DD(2.0) * u * w / (S + (u + DD(1.0)) / DD(2.0));
    const DD alpha = d.hi < 0.5 ? one_p - DD(1.0) : DD(1.0) - one_m;
    return {alpha, sqrt(one_p * one_m)};
}

[[nodiscard]] inline AlphaCos alpha_cos(const DD& delta, const DD& u) {
    const DD d = delta.hi < 0.0 ? DD(0.0) : (delta.hi > 1.0 ? DD(1.0) : delta);
    return alpha_cos(d, DD(1.0) - d, u);
}

[[nodiscard]] inline DD kappa(const DD& delta, const Band& b) {
    return scale_of(b) * alpha_cos(delta, u_of(b)).c / hbar();
}
// kappa with 1 - delta given (delta and one_minus in [0, 1], summing to 1).
[[nodiscard]] inline DD kappa(const DD& delta, const DD& one_minus, const Band& b) {
    return scale_of(b) * alpha_cos(delta, one_minus, u_of(b)).c / hbar();
}
[[nodiscard]] inline DD kappa_antiderivative(const DD& delta, const Band& b) {
    const DD u = u_of(b);
    const AlphaCos a = alpha_cos(delta, u);
    const DD theta = atan2(a.alpha, a.c);
    return scale_of(b) / (DD(2.0) * u * hbar()) * (u * (theta + a.alpha * a.c) / DD(2.0) - a.c * a.c * a.c / DD(3.0));
}
[[nodiscard]] inline DD inverse_kappa_antiderivative(const DD& delta, const Band& b) {
    const DD u = u_of(b);
    const AlphaCos a = alpha_cos(delta, u);
    return hbar() / (DD(2.0) * u * scale_of(b)) * (u * atan2(a.alpha, a.c) - a.c);
}

// {I_k, I_ik, dI_k/da, dI_k/db, dI_ik/da, dI_ik/db}; ok false if a quadrature failed.
struct Segment {
    std::array<DD, 6> q;
    bool ok = true;
};

[[nodiscard]] inline Segment segment(const DD& da, const DD& db, const DD& L, const Band& b) {
    Segment out;
    const DD D = abs(db - da);
    const bool rising = db.hi >= da.hi;
    const DD vmin = rising ? da : db, vmax = rising ? db : da;
    auto clamp01 = [](const DD& v) { return v.hi < 0.0 ? DD(0.0) : (v.hi > 1.0 ? DD(1.0) : v); };
    const DD lo = clamp01(vmin), hi = clamp01(vmax);
    const bool in_lo = vmin.hi > 0.0 && vmin.hi < 1.0, in_hi = vmax.hi > 0.0 && vmax.hi < 1.0;
    auto fk = [&](const DD& v) { return kappa(v, b); };
    auto fi = [&](const DD& v) { return DD(1.0) / kappa(v, b); };
    if (D.hi == 0.0) {
        const DD mid = da;
        if (!(mid.hi > 0.0 && mid.hi < 1.0)) return out;
        out.q = {L * fk(mid), L * fi(mid), DD(0.0), DD(0.0), DD(0.0), DD(0.0)};
        out.ok = false;  // the zero-length limit is not a judged point
        return out;
    }
    const DD w = hi - lo;  // the part inside the gap
    // Short and well inside the gap (the in-gap part within 1e-3 of its distance to the nearer edge): the integrand is
    // analytic on a disc 1e3 times the interval, and one Kronrod-15 rule (exact to degree 22) converges far below
    // u_dd; there the antiderivative difference would lose u_dd K / (kappa D) and tanh-sinh's estimate stalls.
    const bool interior = in_lo && in_hi && w.hi <= 1e-3 * std::min(lo.hi, 1.0 - hi.hi);
    const DD intK = interior ? qk15(fk, lo, hi).value : kappa_antiderivative(hi, b) - kappa_antiderivative(lo, b);
    const DD intG = interior ? qk15(fi, lo, hi).value
                             : inverse_kappa_antiderivative(hi, b) - inverse_kappa_antiderivative(lo, b);
    const DD fk_hi = in_hi ? fk(hi) : DD(0.0), fk_lo = in_lo ? fk(lo) : DD(0.0);
    const DD fi_hi = in_hi ? fi(hi) : DD(0.0), fi_lo = in_lo ? fi(lo) : DD(0.0);
    // P_max = D f(v_max) - int f and P_min = int f - D f(v_min); for D < 1e-2 the in-gap part by quadrature of the
    // difference, the clipped part (D - w) f(v) exactly.
    auto brackets = [&](auto&& f, const DD& f_hi, const DD& f_lo, const DD& integral, DD& pmax, DD& pmin) {
        if (D.hi >= 1e-2 || w.hi == 0.0) {
            pmax = D * f_hi - integral;
            pmin = integral - D * f_lo;
            return true;
        }
        if (interior) {  // D = w here
            pmax = qk15([&](const DD& v) { return f_hi - f(v); }, lo, hi).value;
            pmin = qk15([&](const DD& v) { return f(v) - f_lo; }, lo, hi).value;
            return true;
        }
        const auto a = tanh_sinh([&](const DD& v) { return f_hi - f(v); }, lo, hi, 0.0, 1e-28);
        const auto c = tanh_sinh([&](const DD& v) { return f(v) - f_lo; }, lo, hi, 0.0, 1e-28);
        if (a.status != TanhSinhStatus::ok || c.status != TanhSinhStatus::ok) return false;
        pmax = a.value + (D - w) * f_hi;
        pmin = c.value - (D - w) * f_lo;
        return true;
    };
    DD pk_max, pk_min, pi_max, pi_min;
    if (!brackets(fk, fk_hi, fk_lo, intK, pk_max, pk_min) || !brackets(fi, fi_hi, fi_lo, intG, pi_max, pi_min)) {
        out.ok = false;
        return out;
    }
    const DD LD2 = L / (D * D);
    const DD dk_max = LD2 * pk_max, dk_min = LD2 * pk_min, di_max = LD2 * pi_max, di_min = LD2 * pi_min;
    out.q = {L * intK / D, L * intG / D, rising ? dk_min : dk_max, rising ? dk_max : dk_min,
             rising ? di_min : di_max, rising ? di_max : di_min};
    return out;
}

// x = {slope [J/m], I_k, I_ik [m^2], k_m^2 [1/m^2]}: {G, dG/dslope, dG/dI_k, dG/dI_ik, dG/dk_m^2}.
[[nodiscard]] inline std::array<DD, 5> path_rate(const std::array<DD, 4>& x) {
    if (!(x[2].hi > 0.0)) return {DD(0.0), DD(0.0), DD(0.0), DD(0.0), DD(0.0)};
    const DD pre = abs(x[0]) / (DD(36.0) * hbar());
    const DD z = x[3] * x[2];
    const DD e = exp(-z);
    const DD bracket = -expm1(-z);
    const DD t = exp(-(DD(2.0) * x[1]));
    // g(z) = 1 - e^-z (1 + z) = sum_{k>=2} (-1)^k (k - 1) z^k / k!.
    DD g;
    if (std::abs(z.hi) < 1.0) {
        DD term = z;  // z^k / k! at k = 1
        g = DD(0.0);
        for (int k = 2; k < 60; ++k) {
            term = term * z / DD(static_cast<double>(k));
            const DD add = DD(static_cast<double>(k - 1)) * term;
            g = (k % 2 == 0) ? g + add : g - add;
            if (std::abs(add.hi) < 1e-34 * std::abs(g.hi)) break;
        }
    } else {
        g = bracket - z * e;
    }
    const DD sign(x[0].hi > 0.0 ? 1.0 : (x[0].hi < 0.0 ? -1.0 : 0.0));
    const DD G = pre / x[2] * bracket * t;
    return {G, sign / (DD(36.0) * hbar()) / x[2] * bracket * t, -(DD(2.0) * G), -(pre * t * g / (x[2] * x[2])),
            pre * t * e};
}

// eq (8) at field F [V/m].
[[nodiscard]] inline DD uniform_field_rate(const DD& F, const Band& b) {
    const DD q(base::q_C), h = hbar(), pi_ = pi();
    const DD sm = sqrt(b.reduced_mass_kg);
    const DD pre = q * q * F * F * sm / (DD(18.0) * pi_ * h * h * sqrt(b.gap_J));
    return pre * exp(-(pi_ * sm * b.gap_J * sqrt(b.gap_J) / (DD(2.0) * q * F * h)));
}

}  // namespace NiTCAD::verify::reference::model_generation
