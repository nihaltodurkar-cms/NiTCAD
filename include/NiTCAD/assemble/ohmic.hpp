// Ohmic contact boundary values in scaled variables (ARCHITECTURE.md 6.4; legacy
// Device1D::contact_value, Boltzmann branch): the contact node is held at its charge-neutral
// equilibrium, n0 - p0 = C and n0 p0 = n_ie^2, and its potential is shifted by the applied bias:
//     psi0 = V / V_T + ln(n0 / n_ie).
// The statistics come from physics::boltzmann_neutral_equilibrium (statistics.hpp), so Fermi-Dirac
// contacts later change only that call.
#pragma once

#include "NiTCAD/physics/statistics.hpp"

namespace NiTCAD::assemble {

struct OhmicValue {
    double psi;  // scaled potential psi0
    double n;    // scaled n0
    double p;    // scaled p0
};

// net_doping_scaled = (N_D - N_A) / Ns, n_ie_scaled = n_ie / Ns, bias_scaled = V / V_T.
// Preconditions as physics::boltzmann_neutral_equilibrium (finite doping, positive n_ie).
[[nodiscard]] inline OhmicValue ohmic_contact_value(double net_doping_scaled, double n_ie_scaled,
                                                    double bias_scaled) {
    const physics::NeutralEquilibrium e =
        physics::boltzmann_neutral_equilibrium(net_doping_scaled, n_ie_scaled);
    return {bias_scaled + e.eta, e.n, e.p};
}

}  // namespace NiTCAD::assemble
