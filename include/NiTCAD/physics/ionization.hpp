// Incomplete dopant ionization (ARCHITECTURE.md section 5; Unit 15; legacy device.py
// ionized_eta_doping and ionized_doping, M13 phase 2):
//
//     N_D+ = N_D / (1 + g_D e^(eta_c + E_D / kT)),   eta_c = (E_F - E_c) / kT,
//     N_A- = N_A / (1 + g_A e^(eta_v + E_A / kT)),   eta_v = (E_v - E_F) / kT,
//
// with E_D the donor level's depth below E_c, E_A the acceptor level's height above E_v and g the
// degeneracy factors (physics::IonizationParameters; legacy 45 meV, 2 and 4 for silicon). Out of
// equilibrium the level is in equilibrium with its band: eta_c is the electrons' (from n), eta_v
// the holes' (legacy). A level depth of 0 means complete ionization (the dopant is not modelled).
// Applicability (legacy): a single shallow level per dopant type; invalid above the Mott transition
// (about 4e18 cm^-3 in silicon), where impurity bands form and ionization is complete again.
#pragma once

#include "NiTCAD/physics/statistics.hpp"

namespace NiTCAD::physics {

struct IonizedDensity {
    double value;  // ionized concentration, in the unit of N
    double d_eta;  // d value / d eta (<= 0)
};

// N / (1 + g e^(eta + depth_over_kT)): the ionized part of a dopant concentration N at its band's
// reduced energy eta. depth_over_kT <= 0 gives N (complete ionization, d_eta 0). Overflow-free for
// every eta. A per-node kernel: no checks.
[[nodiscard]] IonizedDensity ionized_density(double N, double eta, double depth_over_kT,
                                             double degeneracy) noexcept;

// The dopant levels of a node in units of kT, as the assemblers use them.
struct DopantLevels {
    double donor_kT;     // E_D / kT, <= 0 for complete ionization
    double acceptor_kT;  // E_A / kT
    double donor_degeneracy;
    double acceptor_degeneracy;
};

// The charge-neutral equilibrium with incompletely ionized dopants:
// n(eta) - p(eta) = N_D+(eta_c) - N_A-(eta_v), with n and p by Boltzmann or Fermi-Dirac statistics
// in the n_ie gauge (statistics.hpp: n = n_ie e^eta gamma, eta_c = eta - g_n,
// eta_v = -eta - g_p), solved for eta by safeguarded Newton from the completely ionized root.
// Homogeneous in the concentrations, like statistics.hpp. Preconditions (NITCAD_EXPECTS): the
// concentrations finite and >= 0; n_ie finite and positive; log_dos finite.
[[nodiscard]] NeutralEquilibrium ionized_neutral_equilibrium(double donors, double acceptors,
                                                             double n_ie, double log_dos_n,
                                                             double log_dos_p,
                                                             const DopantLevels& levels,
                                                             bool fermi_dirac);

}  // namespace NiTCAD::physics
