// Canali velocity-saturation (field-dependent) mobility (legacy materials.mobility_field):
//
//     mu(E) = mu0 / [1 + (mu0 E / v_sat)^beta]^(1 / beta)
//
// with mu0 the low-field mobility and E >= 0 the magnitude of the driving field parallel to the
// current. mu(0) = mu0; the drift velocity mu(E) E rises monotonically to v_sat. Empirical (Canali
// et al., IEEE Trans. Electron Devices 22, 1045 (1975)); silicon beta = 2 for electrons and 1 for
// holes. The legacy applies it lagged (frozen at the previous Newton iterate); NiTCAD returns the
// exact derivative in E so the assembler can differentiate through it (ARCHITECTURE.md 6.2,
// Unit 13).
//
// Unit 68 slice 4b: with CanaliParameters::gamma > 0, the transferred-electron law of Barnes et
// al. (the form of K. Ghosh, U. Singisetti, J. Appl. Phys. 122, 035702 (2017), Table I):
//     mu(E) = (mu0 + (v_sat / E) (E / F_c)^gamma) / (1 + (E / F_c)^gamma),
// mu(0) = mu0; the velocity peaks and falls toward v_sat (negative differential mobility).
#pragma once

#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::physics {

struct FieldMobility {
    double mobility;  // [cm^2/(V s)]
    double d_dE;      // d mobility / dE [cm^3/(V^2 s)], the right-hand derivative at E = 0
};

// Preconditions (NITCAD_EXPECTS): mu0 finite and positive; field_V_per_cm finite and >= 0; the
// parameters as Semiconductor::create validates them (v_sat > 0, beta >= 1, so d_dE is finite;
// at E = 0 it is 0 for beta > 1 and -mu0^2 / v_sat for beta = 1).
[[nodiscard]] FieldMobility canali_mobility(double mu0, double field_V_per_cm,
                                            const CanaliParameters& c);

// The material's Canali parameters for one carrier at temperature T [K], v_sat by its temperature
// law (CanaliParameters::temperature_A; unchanged at 300 K and for A = 0).
[[nodiscard]] inline CanaliParameters saturation(const Semiconductor& m, Carrier carrier,
                                                 double T = 300.0) noexcept {
    CanaliParameters c = carrier == Carrier::electron ? m.parameters().electron_saturation
                                                      : m.parameters().hole_saturation;
    if (c.temperature_A != 0.0) c.v_sat_cm_s /= 1.0 - c.temperature_A + c.temperature_A * T / 300.0;
    return c;
}

}  // namespace NiTCAD::physics
