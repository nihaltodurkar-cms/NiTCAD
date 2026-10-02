// de Mari scaling (ARCHITECTURE.md 6.1, V3): derived once per problem and passed explicitly to
// every assembler, never recomputed per slice or hidden inside a device.
//
// Legacy definitions (core/src/device1d/inputs.cpp:41-65):
//     V_T = k T / q                              [V]
//     eps = eps_r eps0                           [F/cm], of node 0's material
//     Ns  = override if given, else max(max_i |N_D - N_A|_i, n_i)   [cm^-3], n_i of node 0
//     L_D = sqrt(eps V_T / (q Ns))               [cm]
//     J0  = q D0 Ns / L_D                        [A/cm^2]
//     R0  = D0 Ns / L_D^2                        [cm^-3 s^-1]
// with D0 = 1 cm^2/s (legacy kD0Ref = 1.0, unit not stated in the file; cm^2/s is the unit that
// makes the legacy scaled edge diffusivity mu V_T / D0 dimensionless).
// Scaled quantities: potentials / V_T, concentrations / Ns, lengths / L_D (so in D dimensions a
// control volume is divided by L_D^D and a coupling area by L_D^(D-1)), fluxes / J0, rates / R0.
#pragma once

#include <expected>
#include <optional>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"

namespace NiTCAD::assemble {

struct Scaling {
    double temperature_K;
    double V_T;           // [V]
    double eps_F_per_cm;  // reference permittivity
    double n_i;           // intrinsic density of the reference material [cm^-3]
    double Ns;            // concentration scale [cm^-3]
    double L_D;           // length scale [cm]
    double D0;            // diffusivity scale [cm^2/s], 1
    double J0;            // current-density scale [A/cm^2]
    double R0;            // rate scale [cm^-3 s^-1]
};

// The reference material is node 0's, as in the legacy. Errors: invalid_input if Ns_override is
// given and not finite and positive.
[[nodiscard]] std::expected<Scaling, base::Error> make_scaling(
    const device::Device& device, std::optional<double> Ns_override = std::nullopt);

// eps_r eps0 / eps: a relative permittivity over the scaling's reference permittivity (the legacy
// et of a material's edges, eps_ox / eps of a gate oxide).
[[nodiscard]] inline double permittivity_ratio(double eps_r, const Scaling& s) noexcept {
    return eps_r * base::eps0_F_per_cm / s.eps_F_per_cm;
}

}  // namespace NiTCAD::assemble
