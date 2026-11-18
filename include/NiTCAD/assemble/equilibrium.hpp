// Thermal equilibrium (Unit 57 GE5; formerly EquilibriumPoisson): a System created with
// Formulation::equilibrium solves Poisson alone for the scaled potential, one unknown per mesh
// node, with the carriers slaved to psi by Boltzmann or, with models.fermi_dirac, Fermi-Dirac
// statistics (legacy Device1D::solve_equilibrium_boltzmann and solve_equilibrium_fd, generalised
// from 1D to the mesh graph).
//
// Scaled box-method row of node i not on a contact:
//     F_i = sum over edges e = (i, j):  c_e (psi_j - psi_i)  -  V_i (n_i - p_i - C_i)
// with c_e = coupling area / length / L_D^(D-2) times the edge's permittivity over the reference
// one (the harmonic mean of its ends'), V_i = control volume / L_D^D, n = n_ie e^(psi + s),
// p = n_ie e^-(psi + s) (physics::boltzmann_density; Fermi-Dirac: physics::fermi_dirac_density with
// eta = psi + s and -(psi + s)), C = (N_D - N_A) / Ns and n_ie / Ns from the node's material (with
// band-gap narrowing when models.bgn is set), and s the node's band shift (Unit 15; 0 for a single
// material): the intrinsic-level depth of its material below the vacuum level, less node 0's, over
// V_T. psi is then the electrostatic potential, continuous across a heterointerface. In 1D this is
// the legacy row et (psi[i+1] - psi[i]) / h - ... - dV (n - p - C) (with the M33 affinity shift).
// With models.incomplete_ionization C is the ionized doping N_D+ - N_A- (physics/ionization.hpp)
// at eta_c = psi + s - g_n and eta_v = -(psi + s) - g_p, with its psi derivative.
// Contact row (ohmic, Dirichlet): F_i = psi_i - psi0_i, psi0 from ohmic_contact_value at zero bias,
// less s; a Schottky contact's (Unit 48): psi_i - psi_S at zero bias.
// A gate node keeps its box row and gains the oxide term of gate.hpp, G_i (psi_G,i - psi_i) + S_i,
// at the gate's bias (0 unless set_bias says otherwise; legacy Device2D.solve_equilibrium). A gate
// carries no current, so a biased gate keeps the device in thermal equilibrium as long as every
// ohmic contact is at 0 V: the legacy MOS-C solve (moscap.MOSCapacitor.solve_psi), whose C-V is
// the quasi-static one.
// Insulators (Unit 15b): an insulator node's row is the box row with no charge (its carrier
// densities are 0), and an electrode node's (an insulator node) is Dirichlet, F_i = psi_i - psi_E
// with psi_E the electrode potential at its bias (scaled_device.hpp). An edge of a
// semiconductor-insulator interface with charge, traps or recombination carries the half-edge
// fluxes of interface_edges.hpp instead of its own: the fixed charge and the traps' charge sit at
// the interface potential psi_I, with Fermi occupancy at eta_I = psi_I + s.
//
// The terms (src/assemble/equilibrium.cpp, in order): the slaved charge, which assigns every
// row off the Dirichlet nodes (0 on an insulator node); the gate and the Poisson flux, shared with
// the drift-diffusion system (at stride 1 they are these rows' expressions); the interface edges
// at Fermi occupancy, between them; and the Dirichlet rows last, assigned whole. The Jacobian
// pattern is built once (diagonal plus both directions of every edge; contact rows keep their
// off-diagonal entries as explicit zeros), and evaluate() rewrites only the values, so a
// LinearSolver reuses its analysis across Newton iterations (6.10, "Known limits").
//
// Its queries: band_edges (band_edges.hpp; the quasi-Fermi levels are 0, NaN on insulator nodes),
// gate_charges and interface_trap_charges (charges.hpp, the traps at Fermi occupancy) and the
// functions below. Preconditions (NITCAD_EXPECTS) of each below: a system of
// Formulation::equilibrium, and spans of node_count() entries.
#pragma once

#include <cstddef>
#include <expected>
#include <span>
#include <vector>

#include "NiTCAD/assemble/system.hpp"

namespace NiTCAD::assemble {

// System::create(device, scaling, models, Formulation::equilibrium).
[[nodiscard]] std::expected<System, base::Error> create_equilibrium(
    const device::Device& device, const Scaling& scaling, const PhysicsModels& models = {});

// The scaled carrier densities slaved to psi (n = n_ie e^(psi + s) and p = n_ie e^-(psi + s)
// under Boltzmann statistics; 0 on an insulator node).
void equilibrium_carriers(const System& system, std::span<const double> psi, std::span<double> n,
                          std::span<double> p);

// The charge-neutral potential of every node, asinh(C / 2 n_ie) - s under Boltzmann statistics,
// the root of physics::fermi_dirac_neutral_equilibrium less s under Fermi-Dirac (ohmic and
// Schottky contact nodes: their Dirichlet value, which is the same at zero bias for an ohmic
// one); the legacy initial guess. An electrode node takes its Dirichlet value at the current
// bias, another insulator node 0.
[[nodiscard]] std::vector<double> charge_neutral_potential(const System& system);

// Per node: the Dirichlet value on an ohmic or Schottky contact or electrode node (the
// electrode's at its current bias); not meaningful elsewhere.
[[nodiscard]] std::span<const double> contact_potential(const System& system);
// Per node: whether its row is a Dirichlet row (an ohmic or Schottky contact or electrode node);
// gate nodes are not.
[[nodiscard]] bool is_contact(const System& system, std::size_t node);

}  // namespace NiTCAD::assemble
