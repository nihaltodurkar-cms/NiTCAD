// Impact ionization coefficients, van Overstraeten-de Man (legacy ionization.py; Unit 19):
//
//     alpha(E, T) = gamma A exp(-gamma B / E),   gamma = tanh(hw / (2 k T0)) / tanh(hw / (2 k T))
//
// with (A, B) the low branch below the switch field and the high branch from it on, E >= 0 the
// driving field [V/cm], hw the optical-phonon energy and T0 = 300 K (so gamma(300 K) = 1 and the
// coefficients are the published ones there; gamma = 1 when hw is 0; above T0 gamma exceeds 1 and,
// through the exponent, the coefficients fall with temperature). The temperature factor is the
// common one of device simulators (Sentaurus' vanOverstraetendeMan model); the legacy evaluated the
// 300 K coefficients at every temperature. The branches meet in value at silicon's hole switch
// (4e5 V/cm) but not in slope: the derivative is that of the branch in use (the high one at the
// switch itself).
#pragma once

#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::physics {

struct ImpactIonizationRate {
    double alpha;  // [1/cm]
    double d_dE;   // d alpha / dE [1/V]
};

// gamma at lattice temperature T. Precondition (NITCAD_EXPECTS): T finite and positive;
// phonon_energy_eV finite and >= 0.
[[nodiscard]] double impact_ionization_temperature_factor(double phonon_energy_eV,
                                                          double temperature_K);

// alpha and its derivative at field E for one carrier, with the temperature factor gamma. Zero
// below the field at which exp(-gamma B / E) underflows (gamma B / E > 700), and for a carrier with
// A = 0. Preconditions (NITCAD_EXPECTS): E finite and >= 0; gamma finite and positive.
[[nodiscard]] ImpactIonizationRate impact_ionization_coefficient(
    const ImpactIonizationCoefficients& c, double gamma, double field_V_per_cm);

// The material's coefficients for one carrier.
[[nodiscard]] inline const ImpactIonizationCoefficients& impact_ionization(
    const Semiconductor& m, Carrier carrier) noexcept {
    return carrier == Carrier::electron ? m.parameters().impact_ionization.electron
                                        : m.parameters().impact_ionization.hole;
}

}  // namespace NiTCAD::physics
