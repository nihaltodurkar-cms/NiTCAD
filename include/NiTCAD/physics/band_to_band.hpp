// Band-to-band tunnelling (Unit 20; legacy btbt.py).
//
// Kane rate, the local and the calibrated nonlocal model:
//
//     G(F) = A F^2 exp(-B / F)            [cm^-3 s^-1], F [V/cm]
//
// with each material's (A, B) (semiconductor.hpp BandToBandParameters; silicon: Hurkx, Klaassen and
// Knuvers, IEEE Trans. Electron Devices 39, 331 (1992), Table I, as the legacy pins them). Pure
// generation: no Hurkx D factor (owner's decision, Unit 20), so G does not vanish at equilibrium.
//
// Direct-gap WKB rate (Esseni, Pala, Palestri, Alper and Rollo, Semicond. Sci. Technol. 32, 083005
// (2017), section 2.1, eqs. (8), (9), (11), (12); legacy btbt.py M34-S1, including its three
// transcription fixes). SI units throughout (m, J, kg, s): the assembler converts at one boundary.
// Along a tunnel path delta = (E - E_v(x)) / E_g runs from 0 at the start to 1 at the end; with
// u = m0 / (2 m_r) > 1,
//     alpha(delta) = -u + 2 sqrt(u (delta - 1/2) + u^2/4 + 1/4),   cos(theta) = sqrt(1 - alpha^2),
//     kappa(delta) = sqrt(m_r E_g) cos(theta) / hbar                                  (eq. 9)
// vanishing at both band edges. The path's rate is
//     G = |dE_v/dx| / (36 hbar) / I_ik (1 - exp(-k_m^2 I_ik)) exp(-2 I_k)            (eq. 11)
// with I_k = int kappa dx, I_ik = int dx / kappa over the part of the path inside the gap, and
// k_m^2 = min(2 m_v (E_vmax - E), 2 m_c (E - E_cmin)) / hbar^2 (eq. 12). In a uniform field it
// reduces to eq. (8) when k_m^2 I_ik >> 1. This is first-principles physics of a direct gap: it is
// not calibrated to silicon (indirect, phonon-assisted) and must not be used for it.
#pragma once

#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::physics {

struct KaneRate {
    double rate;  // [cm^-3 s^-1]
    double d_dF;  // [cm^-3 s^-1 per V/cm]
};

// G and dG/dF at field F >= 0. Zero at F = 0, for A = 0, and where B / F > 700 (exp underflows).
// Preconditions (NITCAD_EXPECTS): F finite and >= 0; A >= 0 and B > 0 where A > 0.
[[nodiscard]] KaneRate kane_generation(double A_per_cm3_s, double B_V_per_cm,
                                       double field_V_per_cm);

// The tunnelling reduced mass m_c m_v / (m_c + m_v) of a material, in units of m0; 0 when either
// mass is 0 (the material has no direct-gap tunnelling data).
[[nodiscard]] double tunnelling_reduced_mass(const BandToBandParameters& p) noexcept;

// The direct-gap WKB functions of one band structure: gap E_g [J] and reduced mass m_r [kg].
struct WkbBand {
    double gap_J;
    double reduced_mass_kg;
};

// u = m0 / (2 m_r). Precondition (NITCAD_EXPECTS): u > 1 (m_r < m0 / 2), which kappa's vanishing
// at both band edges needs.
[[nodiscard]] double wkb_u(const WkbBand& band);

// kappa(delta) [1/m] (eq. 9), delta clipped to [0, 1]; and an antiderivative in delta of kappa
// [1/m] and of 1 / kappa [m]: over the whole gap they integrate to pi sqrt(m_r E_g) / (4 hbar) and
// pi hbar / (2 sqrt(m_r E_g)).
[[nodiscard]] double wkb_kappa(double delta, const WkbBand& band);
[[nodiscard]] double wkb_kappa_antiderivative(double delta, const WkbBand& band);
[[nodiscard]] double wkb_inverse_kappa_antiderivative(double delta, const WkbBand& band);

// The exact WKB integrals over a straight segment of length L [m] along which delta runs linearly
// from da to db, counting only the part inside the gap: I_k = int kappa dx [-] and
// I_ik = int dx / kappa [m^2], with their partials in da and db. Continuous in (da, db); a segment
// with |db - da| < 1e-10 takes L f(midpoint), the same limit without the cancellation.
struct WkbSegment {
    double I_k, I_ik;
    double dI_k_da, dI_k_db, dI_ik_da, dI_ik_db;
};
[[nodiscard]] WkbSegment wkb_segment(double da, double db, double length_m, const WkbBand& band);

// eq. (11) for one path from its integrals: slope = |dE_v/dx| at the start [J/m], I_k, I_ik
// [m^2] and k_m^2 [1/m^2]; the rate in m^-3 s^-1 and its partials. Zero (with zero partials) when
// I_ik is not positive.
struct WkbPathRate {
    double rate;
    double d_slope, d_I_k, d_I_ik, d_km2;
};
[[nodiscard]] WkbPathRate wkb_path_rate(double slope_J_per_m, double I_k, double I_ik,
                                        double km2_per_m2);

// eq. (8): the closed-form uniform-field rate [m^-3 s^-1] at field F [V/m], for the gates only.
[[nodiscard]] double wkb_uniform_field_rate(double field_V_per_m, const WkbBand& band);

}  // namespace NiTCAD::physics
