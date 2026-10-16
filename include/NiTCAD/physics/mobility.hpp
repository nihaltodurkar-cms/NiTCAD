// Caughey-Thomas doping-dependent low-field mobility (legacy materials.mobility_caughey_thomas):
//
//     mu(N, T) = mu_min + (mu_max (T / 300)^T_exponent - mu_min) / (1 + (N / N_ref)^alpha)
//
// N is the TOTAL ionised impurity concentration N_A + N_D, not the net doping: in a compensated
// region the net doping overestimates the mobility. Purely empirical (Caughey and Thomas,
// Proc. IEEE 55, 2192 (1967)).
//
// The mobility depends only on doping and lattice temperature, neither of which is a Newton
// unknown, so it has no partial derivatives to return (ARCHITECTURE.md section 5, item 2).
#pragma once

#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::physics {

// [cm^2/(V s)]. Preconditions (NITCAD_EXPECTS): total_impurity_cm3 is finite and >= 0;
// temperature_K is finite and positive.
// The legacy clamps N to at least 1 cm^-3; that is not carried over. It changed mu(0) by
// (1 / N_ref)^alpha relative, about 1e-15, and hid negative input.
[[nodiscard]] double caughey_thomas_mobility(const Semiconductor& m, Carrier carrier,
                                             double total_impurity_cm3, double temperature_K);

// Masetti doping-dependent mobility (Unit 34; DECISIONS.md, Unit 34): G. Masetti, M. Severi,
// S. Solmi, IEEE Trans. Electron Devices 30, 764 (1983), eq (4) and Table I, as reproduced by
// CERN Allpix Squared, Crosslight and Flooxs (the paper itself was not read):
//
//     mu(N) = mu_min1 exp(-P_c / N) + (mu_max(T) - mu_min2) / (1 + (N / C_r)^alpha)
//             - mu_1 / (1 + (C_s / N)^beta),       mu_max(T) = mu_max (T / 300)^T_exponent
//
// N the total ionized impurity density. The fits are at 300 K; the mu_max scaling (exponent -2.2
// for holes, -2.5 for electrons) is the one those reproductions apply, not the paper's. At N = 0
// the mobility is mu_max(T) - mu_min2 (the exponential and the last term vanish), or with
// P_c = 0 (the electron sets, eq (1); Unit 38) mu_min1 + mu_max(T) - mu_min2.
struct MasettiParameters {
    double mu_min1;     // [cm^2/(V s)]
    double mu_min2;     // [cm^2/(V s)]
    double mu_max;      // [cm^2/(V s)] at 300 K
    double mu_1;        // [cm^2/(V s)]
    double P_c;         // [cm^-3]
    double C_r;         // [cm^-3]
    double C_s;         // [cm^-3]
    double alpha;
    double beta;
    double T_exponent;  // mu_max(T) = mu_max (T / 300)^T_exponent (not from the paper)
};

// Holes in boron-doped silicon (Table I as reproduced; DECISIONS.md, Unit 34).
inline constexpr MasettiParameters masetti_boron_silicon{
    .mu_min1 = 44.9, .mu_min2 = 0.0, .mu_max = 470.5, .mu_1 = 29.0, .P_c = 9.23e16,
    .C_r = 2.23e17, .C_s = 6.10e20, .alpha = 0.719, .beta = 2.0, .T_exponent = -2.2};

// Electrons in phosphorus-doped silicon (Table I as reproduced; DECISIONS.md, Unit 38): the
// paper's eq (1), the form above with P_c = 0.
inline constexpr MasettiParameters masetti_phosphorus_silicon{
    .mu_min1 = 68.5, .mu_min2 = 68.5, .mu_max = 1414.0, .mu_1 = 56.1, .P_c = 0.0,
    .C_r = 9.20e16, .C_s = 3.41e20, .alpha = 0.711, .beta = 1.98, .T_exponent = -2.5};

// [cm^2/(V s)]. Preconditions (NITCAD_EXPECTS): total_impurity_cm3 is finite and >= 0;
// temperature_K is finite and positive.
[[nodiscard]] double masetti_mobility(const MasettiParameters& p, double total_impurity_cm3,
                                      double temperature_K);
// d mu / dT [cm^2/(V s K)] (only mu_max depends on T). Preconditions as masetti_mobility.
[[nodiscard]] double masetti_mobility_slope(const MasettiParameters& p, double total_impurity_cm3,
                                            double temperature_K);

}  // namespace NiTCAD::physics
