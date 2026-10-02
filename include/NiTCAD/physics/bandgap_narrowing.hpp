// Slotboom (de Graaff) heavy-doping band-gap narrowing (legacy materials.bandgap_narrowing_slotboom
// and materials.nie_effective):
//
//     dEg(N) = E0 [ ln(N / N0) + sqrt(ln^2(N / N0) + 1/2) - sqrt(1/2) ]   for N > N0, else 0
//     n_ie   = n_i exp(dEg / (2 k T))
//
// with N the total ionised impurity N_A + N_D. The "- sqrt(1/2)" is the legacy's: it makes dEg
// zero at N0, so the cut-off below N0 introduces no step in n_ie. Empirical, fitted to bipolar
// transistor data. The effective n_ie enters the carrier statistics, the SRH and Auger equilibrium
// product and the ohmic contact values; the assembler adds ln(n_ie) to the Scharfetter-Gummel
// driving term so that equilibrium carries no current where n_ie varies (assemble/drift_diffusion).
#pragma once

#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::physics {

// [eV]. Precondition (NITCAD_EXPECTS): total_impurity_cm3 is finite and >= 0. The legacy clamp of
// N to at least 1 cm^-3 is not carried over (the narrowing is zero there anyway).
[[nodiscard]] double bandgap_narrowing_eV(const Semiconductor& m, double total_impurity_cm3);

// [cm^-3]. Preconditions: as bandgap_narrowing_eV, and temperature_K finite and positive.
[[nodiscard]] double effective_intrinsic_density(const Semiconductor& m, double total_impurity_cm3,
                                                 double temperature_K);

}  // namespace NiTCAD::physics
