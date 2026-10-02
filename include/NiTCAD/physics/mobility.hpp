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

}  // namespace NiTCAD::physics
