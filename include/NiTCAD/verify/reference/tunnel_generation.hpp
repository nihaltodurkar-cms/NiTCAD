// Band-to-band generation on the p-i-n reference potential: the local Kane rate's integral, the nonlocal Kane paths
// and the direct-gap WKB paths (Unit 55, TVR-5a, TVB-24; DECISIONS Unit 55 17.2 (m)-(o), 17.3, 17.13 (h)). Double-
// double; base's constants only through model_generation.hpp.
//
// The rates are TVR-3's re-coded model_generation (kane, kappa, path_rate), reused unchanged; what this file adds is
// where they are evaluated: on pin_junction's potential (QE, three regions). Every integral is taken in the level s
// (pin_junction.hpp), dx = V_T ds / E(s), so no profile is inverted. g = E_g / V_T is the gap in levels: one material,
// no band-gap narrowing, so the conduction band at level s + g lies at the valence band's energy at s.
//   - Local Kane (17.2 (n)): P = int G(E(x)) dx = V_T int A E e^(-B / E) ds over (0, B), split at s_0 and s_d (where
//     dE/ds jumps), the integrand 0 where B / E > 700 (as physics::kane_generation).
//   - Nonlocal Kane (17.2 (o); drift_diffusion.hpp's `kane`): a path starting at level s_a crosses where E_c reaches
//     the start's E_v, at s_a + g; its length l = x(s_a + g) - x(s_a) (pin_position) and mean field F = E_g / l
//     [V/cm], G = A F^2 e^(-B / F) with the start's (A, B). P = int G dx_a = V_T int G(s_a) / E(s_a) ds_a over the
//     starts with a crossing (s_a < B - g), split at s_0, s_d, s_0 - g and s_d - g.
//   - Direct WKB (drift_diffusion.hpp's `direct_wkb`; Esseni et al. 2017 eqs (9), (11), (12), PA, as TVR-3): along
//     the path delta = (s - s_a) / g runs from 0 to 1; I_k = int kappa(delta) dx and I_ik = int dx / kappa over
//     [s_a, s_a + g]; the slope |dE_v/dx| at the start is q E(s_a); k_m^2 = min(2 m_v q V_T (s_a - s_v),
//     2 m_c q V_T (s_c - s_a - g)) / hbar^2, with s_v = min(s_a, s_lo) and s_c = max(s_a + g, s_hi) (so both windows
//     are >= 0, and 0 for a start below s_lo or a crossing beyond s_hi), s_lo and s_hi the levels on the p and n sides
//     where E falls to the flat field. NiTCAD's backward polyline (E_vmax) and forward
//     one (E_cmin) stop there (tunnel_paths.hpp:31-34, `flat_field_V_cm` = 1e3 V/cm in tunnel_paths.cpp:18); the
//     value is an input here. P over the starts in (s_lo, s_hi - g), split also where the two windows of k_m^2 meet.
//
// The one-sided start slope (TVB-24, DECISIONS 17.24): NiTCAD takes |dE_v/dx| over the path's first segment, half a
// cell on a 1D grid, which is the forward difference (psi_{i+1} - psi_i) / h: the slope times 1 + (h / 2) E'/E +
// O(h^2), E' = dE/dx (pin_field_slope). wkb_slope_coefficient returns c_1 = (1/2) int G E'/E dx / P [1/cm], so that
// the pairs' relative error has the leading term c_1 h. It is a prediction stated before the run, not a tolerance.
//
// Methods (17.3): adaptive Kronrod-15 for the outer integrals, pin_position for each path's length, tanh-sinh for the
// WKB integrals (1 / kappa has an inverse square root at both band edges). Second method: tanh-sinh for the outer
// integrals; for the lengths tanh-sinh of V_T / E (in s in the i region, in the logarithm of the distance to the bulk
// level on a side, where a path may end e^-30 from it); for the WKB integrals Kronrod-15 after delta = t^2 (3 - 2 t),
// which removes the band-edge singularities, its tolerance relative to the whole integral. Each result carries both
// values.
#pragma once

#include <expected>
#include <string>

#include "NiTCAD/verify/reference/double_double.hpp"
#include "NiTCAD/verify/reference/pin_junction.hpp"

namespace NiTCAD::verify::reference {

struct KaneInput {
    double A_per_cm3_s = 0.0;
    double B_V_per_cm = 0.0;
};

struct WkbInput {
    double gap_eV = 0.0;
    double electron_mass = 0.0;  // m_c / m0
    double hole_mass = 0.0;      // m_v / m0
    double flat_field_V_per_cm = 1e3;
};

struct TunnelPairs {
    DD pairs_cm2_s;          // int G dx, the method
    DD second_cm2_s;         // the second method
    double abs_error = 0.0;  // the method's own estimate
};

struct PathReference {
    bool exists = false;      // the start has a crossing (s_a + g < B)
    DD length_cm, field_V_per_cm, rate_cm3_s;
    DD second_length_cm;      // the length by the second method
    DD second_rate_cm3_s;     // the rate from the second method's integrals (WKB) or length (Kane)
    // WKB only: the path integrals and the windows, SI as model_generation.
    DD I_k, I_ik_m2, km2_per_m2, slope_J_per_m;
};

// Errors (all): kane A or B not positive, a gap or a mass not positive (u = m0 / (2 m_r) must exceed 1 for WKB), a
// level outside (0, B), a quadrature or Brent failing.
[[nodiscard]] std::expected<TunnelPairs, std::string> local_kane_pairs(const PinJunction& j, const KaneInput& k);

[[nodiscard]] std::expected<PathReference, std::string> kane_path(const PinJunction& j, double gap_eV,
                                                                  const KaneInput& k, const DD& s_start);
[[nodiscard]] std::expected<TunnelPairs, std::string> nonlocal_kane_pairs(const PinJunction& j, double gap_eV,
                                                                          const KaneInput& k);

struct WkbLevels {
    DD s_lo, s_hi;  // where E equals the flat field on the p and n sides
};
[[nodiscard]] std::expected<WkbLevels, std::string> wkb_levels(const PinJunction& j, const WkbInput& w);

[[nodiscard]] std::expected<PathReference, std::string> wkb_path(const PinJunction& j, const WkbInput& w,
                                                                 const DD& s_start);
[[nodiscard]] std::expected<TunnelPairs, std::string> wkb_pairs(const PinJunction& j, const WkbInput& w);

// c_1 = (1/2) int G E'/E dx / P [1/cm] of the WKB pairs (the header comment).
[[nodiscard]] std::expected<DD, std::string> wkb_slope_coefficient(const PinJunction& j, const WkbInput& w);

}  // namespace NiTCAD::verify::reference
