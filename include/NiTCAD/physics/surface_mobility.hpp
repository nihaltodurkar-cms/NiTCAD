// Lombardi surface mobility (Unit 46; DECISIONS.md, Unit 46): C. Lombardi, S. Manzini,
// A. Saporito, M. Vanzi, IEEE Trans. CAD 7, 1164 (1988), as reproduced by the Xyce Reference
// Guide and Genius TCAD (the paper itself was not read). Two inversion-layer scattering terms in
// the field E_perp normal to the current, combined with a bulk mobility by Matthiessen's rule:
//
//     1 / mu = 1 / mu_bulk + 1 / mu_ac + 1 / mu_sr,
//     mu_ac = B / E_perp + C N^lambda / (T E_perp^(1/3))    (acoustic phonons),
//     mu_sr = delta / E_perp^2                              (surface roughness),
//
// N the total impurity density [cm^-3], T the lattice temperature [K], E_perp [V/cm]. The bulk
// mobility is the caller's (any doping-dependent model); Canali's velocity saturation, when on,
// applies to the result.
#pragma once

namespace NiTCAD::physics {

struct LombardiParameters {
    double B;       // [cm/s]
    double C;       // [cm^(5/3 + 3 lambda) K V^(-2/3) s^-1]
    double lambda;  // the doping exponent
    double delta;   // [V/s]
};

// Silicon's electrons and holes (Lombardi et al.'s values as reproduced; DECISIONS.md, Unit 46).
inline constexpr LombardiParameters lombardi_silicon_electrons{
    .B = 4.75e7, .C = 1.74e5, .lambda = 0.125, .delta = 5.82e14};
inline constexpr LombardiParameters lombardi_silicon_holes{
    .B = 9.93e6, .C = 8.84e5, .lambda = 0.0317, .delta = 2.05e14};

struct SurfaceMobility {
    double mobility;  // mu [cm^2/(V s)]
    double d_dE;      // d mu / d E_perp [cm^3/(V^2 s)]
    double d_dT;      // d mu / d T at a fixed bulk mobility [cm^2/(V s K)]
    double d_dbulk;   // d mu / d mu_bulk: (mu / mu_bulk)^2
};

// Preconditions (NITCAD_EXPECTS): bulk_mobility finite and positive; total_impurity_cm3 finite and
// >= 0; temperature_K finite and positive; field_V_per_cm finite and positive (the assembler keeps
// it at least E_0 = 1 V/cm).
[[nodiscard]] SurfaceMobility lombardi_mobility(const LombardiParameters& p, double bulk_mobility,
                                                double total_impurity_cm3, double temperature_K,
                                                double field_V_per_cm);

}  // namespace NiTCAD::physics
