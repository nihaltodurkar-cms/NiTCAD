// Ohmic contact boundary values in scaled variables (ARCHITECTURE.md 6.4; legacy
// Device1D::contact_value): the contact node is held at its charge-neutral equilibrium,
// n0 - p0 = C, and its potential is the equilibrium one shifted by the applied bias:
//     psi0 = V / V_T + eta0,
// with eta0 the reduced potential of the neutral equilibrium (Boltzmann: n0 p0 = n_ie^2 and
// eta0 = ln(n0 / n_ie); Fermi-Dirac, Unit 14: physics::fermi_dirac_neutral_equilibrium, the
// legacy _fd_contact_values). The statistics come from statistics.hpp; only the call differs.
#pragma once

#include "NiTCAD/physics/statistics.hpp"

namespace NiTCAD::assemble {

struct OhmicValue {
    double psi;  // scaled potential psi0
    double n;    // scaled n0
    double p;    // scaled p0
};

// From a neutral equilibrium (in scaled units) under either statistics; bias_scaled = V / V_T.
[[nodiscard]] constexpr OhmicValue ohmic_contact_value(const physics::NeutralEquilibrium& e,
                                                       double bias_scaled) noexcept {
    return {bias_scaled + e.eta, e.n, e.p};
}

// Boltzmann: net_doping_scaled = (N_D - N_A) / Ns, n_ie_scaled = n_ie / Ns, bias_scaled = V / V_T.
// Preconditions as physics::boltzmann_neutral_equilibrium (finite doping, positive n_ie).
[[nodiscard]] inline OhmicValue ohmic_contact_value(double net_doping_scaled, double n_ie_scaled,
                                                    double bias_scaled) {
    return ohmic_contact_value(
        physics::boltzmann_neutral_equilibrium(net_doping_scaled, n_ie_scaled), bias_scaled);
}

}  // namespace NiTCAD::assemble
