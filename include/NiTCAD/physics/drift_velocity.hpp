// High-field drift-velocity laws in velocity form, for beta < 1 (Unit 68 slice 4c infrastructure;
// DECISIONS.md amendment 10; owner rulings 2026-11-27). Sources: F. R. Is-haqzai, M. Deniz,
// K. Kroeninger, J. Weingarten, Diamond Relat. Mater. 164, 113558 (2026), Eqs. (5) and (30) and
// Table 11 (diamond, 300 K, 0.009-12 V/um; no temperature scaling).
//
//     Caughey-Thomas:  v(E) = mu0 E / [1 + (mu0 E / v_s)^beta]^(1/beta)
//     piecewise:       v(E) = mu01 E                     for E < E_t,
//                      v(E) = Caughey-Thomas with mu0 = v_s / E_c   for E >= E_t,
//
// with E_t where the two branches meet (owner ruling: continuity in place of the paper's fixed
// 0.08 V/um, where they differ by 2 %): E_t = E_c [(v_s / (E_c mu01))^beta - 1]^(1/beta).
//
// The velocity and its derivative dv/dE = mu0 [1 + x^beta]^(-1 - 1/beta), x = mu0 E / v_s, are
// returned directly: dv/dE is mu0 at E = 0 for every beta > 0, where d mu / dE (the Canali form of
// field_mobility.hpp) is singular for beta < 1. The low-field mobility (mu0, mu01) is the caller's,
// doping dependent; v_s, E_c and beta are the published ones. With a doping-dependent mobility the
// laws are an engineering extension, not validated by Table 11 (owner ruling).
//
// Not yet used by the assembler: no registered material carries these laws (diamond waits for its
// sources).
#pragma once

#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::physics {

struct DriftVelocity {
    double velocity;  // [cm/s], E >= 0
    double d_dE;      // dv/dE [cm^2/(V s)]
};

// A published high-field fit: v_s [cm/s], E_c [V/cm] (mu0,2 = v_s / E_c) and beta > 0.
struct HighFieldLaw {
    double v_sat_cm_s;
    double E_c_V_per_cm;
    double beta;
};

// Preconditions (NITCAD_EXPECTS): mu0, v_sat and beta finite and positive; E finite and >= 0.
[[nodiscard]] DriftVelocity caughey_thomas_velocity(double mu0, double field_V_per_cm, double v_sat_cm_s,
                                                    double beta);

// The field where mu01 E meets the law's Caughey-Thomas branch. Precondition: 0 < mu01 < v_s / E_c.
[[nodiscard]] double piecewise_switch_field(double mu01, const HighFieldLaw& law);
[[nodiscard]] DriftVelocity piecewise_velocity(double mu01, double field_V_per_cm, const HighFieldLaw& law);

// S. J. Rashid et al., IEEE Trans. Electron Devices 55, 2744 (2008), as Is-haqzai et al. (2026)
// Eqs. (10)-(11) give it, on the existing Caughey-Thomas parameters (no new law):
//     mu(T, N) = mu_min(T) + (mu_max(T) - mu_min(T)) / (1 + (T / 300)^m1 (N / C0)^gamma2),
//     mu_min/max(T) = mu_min/max,0 (T / 300)^n_min/max,
// is alpha = gamma2, N_ref(T) = C0 (T / 300)^(-m1 / gamma2).
[[nodiscard]] constexpr CaugheyThomasParameters rashid_mobility(double mu_min0, double n_min, double mu_max0,
                                                                double n_max, double C0, double m1, double gamma2) {
    return {.mu_min = mu_min0, .mu_max = mu_max0, .N_ref = C0, .alpha = gamma2, .T_exponent = n_max,
            .mu_min_T_exponent = n_min, .N_ref_T_exponent = -m1 / gamma2, .alpha_T_exponent = 0.0};
}

}  // namespace NiTCAD::physics
