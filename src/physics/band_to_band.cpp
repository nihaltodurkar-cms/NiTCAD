#include "NiTCAD/physics/band_to_band.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::physics {

namespace {

constexpr double hbar = base::hbar_J_s;

// alpha(delta) and cos(theta) = sqrt(1 - alpha^2) without cancellation at the band edges (legacy
// _kane_alpha_c): with S = sqrt(u (delta - 1/2) + u^2/4 + 1/4), S0 = (u - 1)/2 and S1 = (u + 1)/2
// its values at delta = 0 and 1,
//     1 + alpha = 2 u delta / (S + S0),   1 - alpha = 2 u (1 - delta) / (S + S1).
struct AlphaCos {
    double alpha, c, S;
};

AlphaCos alpha_cos(double delta, double u) {
    const double d = std::clamp(delta, 0.0, 1.0);
    const double S = std::sqrt(u * (d - 0.5) + 0.25 * u * u + 0.25);
    const double one_p = 2.0 * u * d / (S + 0.5 * (u - 1.0));
    const double one_m = 2.0 * u * (1.0 - d) / (S + 0.5 * (u + 1.0));
    const double alpha = d < 0.5 ? one_p - 1.0 : 1.0 - one_m;
    return {alpha, std::sqrt(one_p * one_m), S};
}

double scale(const WkbBand& b) { return std::sqrt(b.reduced_mass_kg * b.gap_J); }

// d kappa / d delta = -sqrt(m_r E_g) alpha u / (hbar c S), from d kappa / d alpha =
// -sqrt(m_r E_g) alpha / (hbar c) and d alpha / d delta = u / S. Singular at the band edges.
double kappa_slope(double delta, const WkbBand& b) {
    const double u = wkb_u(b);
    const AlphaCos a = alpha_cos(delta, u);
    return -scale(b) * a.alpha * u / (hbar * std::max(a.c, 1e-300) * a.S);
}

}  // namespace

KaneRate kane_generation(double A, double B, double F) {
    NITCAD_EXPECTS(std::isfinite(F) && F >= 0.0);
    NITCAD_EXPECTS(std::isfinite(A) && A >= 0.0 && std::isfinite(B) && (A == 0.0 || B > 0.0));
    if (A == 0.0 || F == 0.0) return {0.0, 0.0};
    const double x = B / F;
    if (x > 700.0) return {0.0, 0.0};  // exp(-x) underflows; so would the derivative
    const double G = A * F * F * std::exp(-x);
    return {G, G * (2.0 + x) / F};  // G (2 / F + B / F^2)
}

double tunnelling_reduced_mass(const BandToBandParameters& p) noexcept {
    if (p.electron_mass == 0.0 || p.hole_mass == 0.0) return 0.0;
    return p.electron_mass * p.hole_mass / (p.electron_mass + p.hole_mass);
}

double wkb_u(const WkbBand& band) {
    NITCAD_EXPECTS(std::isfinite(band.gap_J) && band.gap_J > 0.0);
    NITCAD_EXPECTS(std::isfinite(band.reduced_mass_kg) && band.reduced_mass_kg > 0.0);
    const double u = base::m0_kg / (2.0 * band.reduced_mass_kg);
    NITCAD_EXPECTS(u > 1.0);
    return u;
}

double wkb_kappa(double delta, const WkbBand& band) {
    return scale(band) * alpha_cos(delta, wkb_u(band)).c / hbar;
}

// With theta = asin(alpha) and eq. (9)'s own d(delta) = (alpha + u) / (2 u) d(alpha):
//     int kappa d(delta)    = sqrt(m_r E_g) / (2 u hbar) (u (theta + sin cos) / 2 - cos^3 / 3),
//     int d(delta) / kappa  = hbar / (2 u sqrt(m_r E_g)) (u theta - cos).
double wkb_kappa_antiderivative(double delta, const WkbBand& band) {
    const double u = wkb_u(band);
    const AlphaCos a = alpha_cos(delta, u);
    const double theta = std::atan2(a.alpha, a.c);
    return scale(band) / (2.0 * u * hbar) *
           (0.5 * u * (theta + a.alpha * a.c) - a.c * a.c * a.c / 3.0);
}

double wkb_inverse_kappa_antiderivative(double delta, const WkbBand& band) {
    const double u = wkb_u(band);
    const AlphaCos a = alpha_cos(delta, u);
    const double theta = std::atan2(a.alpha, a.c);
    return hbar / (2.0 * u * scale(band)) * (u * theta - a.c);
}

WkbSegment wkb_segment(double da, double db, double L, const WkbBand& band) {
    NITCAD_EXPECTS(std::isfinite(da) && std::isfinite(db) && std::isfinite(L) && L >= 0.0);
    const double D = std::abs(db - da);
    WkbSegment s{};
    if (D < 1e-10) {
        // L f(midpoint); the derivative split equally between the ends. A midpoint within 1e-30
        // of the start edge would overflow kappa' / kappa^2 (~ delta^-1.5): such a segment is
        // flat at the edge and suppresses its path either way.
        const double mid = 0.5 * (da + db);
        if (!(mid > 0.0 && mid < 1.0)) return s;
        const double m = std::max(mid, 1e-30);
        const double k = wkb_kappa(m, band);
        const double dk = kappa_slope(m, band);
        s.I_k = L * k;
        s.I_ik = L / k;
        s.dI_k_da = s.dI_k_db = 0.5 * L * dk;
        s.dI_ik_da = s.dI_ik_db = -0.5 * L * dk / (k * k);
        return s;
    }
    const bool rising = db >= da;
    const double vmin = rising ? da : db, vmax = rising ? db : da;
    const double lo = std::clamp(vmin, 0.0, 1.0), hi = std::clamp(vmax, 0.0, 1.0);
    const bool in_lo = vmin > 0.0 && vmin < 1.0, in_hi = vmax > 0.0 && vmax < 1.0;
    // Both differences are >= 0 analytically; clamp the rounding of two O(1) values.
    const double dK = std::max(
        wkb_kappa_antiderivative(hi, band) - wkb_kappa_antiderivative(lo, band), 0.0);
    const double dG = std::max(wkb_inverse_kappa_antiderivative(hi, band) -
                                   wkb_inverse_kappa_antiderivative(lo, band),
                               0.0);
    const double k_hi = in_hi ? wkb_kappa(hi, band) : 0.0;
    const double k_lo = in_lo ? wkb_kappa(lo, band) : 0.0;
    const double i_hi = in_hi ? 1.0 / k_hi : 0.0;  // 1 / kappa only strictly inside the gap
    const double i_lo = in_lo ? 1.0 / k_lo : 0.0;
    s.I_k = L * dK / D;
    s.I_ik = L * dG / D;
    // d/d vmax (L dK / D) = L (f(vmax) - dK / D) / D; d/d vmin = L (-f(vmin) + dK / D) / D.
    const double dk_max = L * (k_hi - dK / D) / D, dk_min = L * (-k_lo + dK / D) / D;
    const double di_max = L * (i_hi - dG / D) / D, di_min = L * (-i_lo + dG / D) / D;
    s.dI_k_da = rising ? dk_min : dk_max;
    s.dI_k_db = rising ? dk_max : dk_min;
    s.dI_ik_da = rising ? di_min : di_max;
    s.dI_ik_db = rising ? di_max : di_min;
    return s;
}

WkbPathRate wkb_path_rate(double slope, double I_k, double I_ik, double km2) {
    if (!(I_ik > 0.0)) return {};
    const double pre = std::abs(slope) / (36.0 * hbar);
    const double e = std::exp(-km2 * I_ik);
    const double bracket = -std::expm1(-km2 * I_ik);  // 1 - exp(-km2 I_ik)
    const double t = std::exp(-2.0 * I_k);
    const double G = pre / I_ik * bracket * t;
    WkbPathRate r{};
    r.rate = G;
    r.d_slope = (slope > 0.0 ? 1.0 : slope < 0.0 ? -1.0 : 0.0) / (36.0 * hbar) / I_ik * bracket * t;
    r.d_I_k = -2.0 * G;
    r.d_I_ik = pre * t * (km2 * e / I_ik - bracket / (I_ik * I_ik));
    r.d_km2 = pre * t * e;  // pre / I_ik * t * I_ik e
    return r;
}

double wkb_uniform_field_rate(double F, const WkbBand& band) {
    constexpr double q = base::q_C;
    constexpr double pi = std::numbers::pi;
    const double mr = band.reduced_mass_kg, Eg = band.gap_J;
    const double f = std::max(std::abs(F), 1e-30);
    // 18 pi, not 18 pi^2; the exponent's denominator 2 hbar q F (legacy fixes 3 and 1).
    const double pre = q * q * f * f * std::sqrt(mr) / (18.0 * pi * hbar * hbar * std::sqrt(Eg));
    return pre * std::exp(-pi * std::sqrt(mr) * std::pow(Eg, 1.5) / (2.0 * q * f * hbar));
}

}  // namespace NiTCAD::physics
