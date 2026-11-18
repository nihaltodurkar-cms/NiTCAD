// The currents of a system's state (Unit 57, DECISIONS.md GE3; system.hpp's header comment has
// their definitions): per edge, per contact, after a time step and in the frequency domain.
#pragma once

#include <complex>
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/system.hpp"

namespace NiTCAD::assemble {

// Scaled (Jn, Jp) of every mesh edge, in units of J0 times the edge's scaled coupling area (so a
// 1D edge gives the current density in units of J0).
[[nodiscard]] std::vector<std::pair<double, double>> edge_currents(const System& system,
                                                                   std::span<const double> x);
// Scaled total current entering the device through each contact (conventional current from the
// contact into the device), in device.contacts() order; zero for a gate. Physical value: times
// J0 L_D^(D-1), in A / cm^(3-D) (A/cm^2 in 1D, A/cm in 2D, A in 3D).
[[nodiscard]] std::vector<double> terminal_currents(const System& system,
                                                    std::span<const double> x);
// Per ohmic contact, a bound on how far the current through any cut of the device can differ
// from terminal_currents (same scale): the sum over the nodes of the total-current residual
// |F_n + F_p| plus 8 eps times the magnitudes of the terms that cancel in it (the one-sided flux
// terms, density times its flux coefficient, and the recombination terms). The same value for
// every ohmic contact; zero for a gate. A current below it is not resolved by the (psi, n, p)
// state.
[[nodiscard]] std::vector<double> terminal_current_resolution(const System& system,
                                                              std::span<const double> x);
// Per mesh edge, the magnitudes of the one-sided terms whose difference is the edge's total
// current: |d_c1 c_a| + |d_c2 c_b| of each carrier (density times its flux coefficient), in the
// units of edge_currents; 0 for an edge without carriers. It scales the edge current's rounding
// (8 eps times it), and half of it is the edge's conductance to a quasi-Fermi potential
// difference near equilibrium (the weights of the weighted current, Unit 32).
[[nodiscard]] std::vector<double> edge_flux_terms(const System& system,
                                                  std::span<const double> x);
// Scaled conduction current entering through each contact after `step`: terminal_currents plus,
// on an interface edge whose semiconductor node is on an ohmic contact, the current that charges
// its traps from the contact.
[[nodiscard]] std::vector<double> conduction_currents(const System& system,
                                                      std::span<const double> x,
                                                      const System::TimeStep& step);

// The small-signal total current entering through a contact, conduction plus s times its charge
// (conduction_currents and contact_charges, linearized about the steady state x):
//     dI = sum_k values[k] dx[columns[k]] + bias dV,
// dV the contact's own bias in V (a gate's charge depends on it directly). Scaled as
// terminal_currents; one row per contact, in device.contacts() order. s as small_signal_matrix's.
struct CurrentRow {
    std::vector<std::size_t> columns;
    std::vector<std::complex<double>> values;
    std::complex<double> bias;
};
[[nodiscard]] std::vector<CurrentRow> small_signal_currents(const System& system,
                                                            std::span<const double> x,
                                                            std::complex<double> s);

}  // namespace NiTCAD::assemble
