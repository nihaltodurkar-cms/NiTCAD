// The charges of a system's state (Unit 57, DECISIONS.md GE3): the contacts', the gates' and the
// interface traps', and the traps' occupancies.
#pragma once

#include <span>
#include <vector>

#include "NiTCAD/assemble/system.hpp"

namespace NiTCAD::assemble {

// The occupancy f of every trap slot at state x: after `step`, or without one the steady state.
[[nodiscard]] std::vector<double> trap_occupancies(const System& system,
                                                   std::span<const double> x,
                                                   const System::TimeStep* step = nullptr);
// Scaled charge of every contact at state x, in device.contacts() order: the displacement flux
// leaving its nodes along the edges to nodes outside it (a gate: as gate_charges). Unlike
// gate_charges it is not zero for an ohmic contact. Units as gate_charges.
[[nodiscard]] std::vector<double> contact_charges(const System& system,
                                                  std::span<const double> x,
                                                  const System::TimeStep* step = nullptr);
// Scaled charge on each gate, the sum over its nodes of G_i (psi_G,i - psi_i) (gate.hpp), and on
// each electrode, the displacement flux leaving its nodes along the edges to other nodes, in
// device.contacts() order; zero for an ohmic contact. Physical value: times q Ns L_D^D, in
// C / cm^(3-D).
[[nodiscard]] std::vector<double> gate_charges(const System& system, std::span<const double> x,
                                               const System::TimeStep* step = nullptr);
// Scaled trapped charge of each declared interface at state x (interface_edges.hpp, without the
// fixed charge), in units of q Ns L_D^D. Precondition (NITCAD_EXPECTS): x has unknowns() entries.
// In thermal equilibrium (equilibrium.hpp) gate_charges and interface_trap_charges take the traps
// at Fermi occupancy and no step (precondition); the other two are drift-diffusion only.
[[nodiscard]] std::vector<double> interface_trap_charges(const System& system,
                                                         std::span<const double> x,
                                                         const System::TimeStep* step = nullptr);

}  // namespace NiTCAD::assemble
