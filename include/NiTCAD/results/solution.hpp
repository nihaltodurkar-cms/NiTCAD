// Fields and terminal quantities of a solved state (ARCHITECTURE.md 6.6): plain data in the public
// unit convention (V, cm^-3, A / cm^(3-D), C / cm^(3-D)), indexed by mesh node, mesh edge or
// contact.
#pragma once

#include <cstddef>
#include <vector>

#include "NiTCAD/results/convergence.hpp"

namespace NiTCAD::results {

// Per mesh node: electrostatic potential referenced to the intrinsic level of node 0's material
// (with insulators, of the lowest semiconductor node's), carrier densities (0 on insulator
// nodes).
struct NodeFields {
    std::vector<double> potential_V;
    std::vector<double> n_cm3;
    std::vector<double> p_cm3;
    // The lattice temperature in K with the electrothermal model (Unit 23); empty without it (the
    // device's uniform temperature). In an initial state, empty starts every node at the device's
    // reference temperature.
    std::vector<double> temperature_K;
};

// Per mesh node, the band diagram in eV, measured from the equilibrium Fermi level (the Fermi level
// of a contact at 0 V; a contact at bias V holds its carriers' Fermi level at -V): the conduction
// and valence band edges (band-gap narrowing shared between them) and the electron and hole
// quasi-Fermi levels, by the selected statistics. Unlike the potential it does not depend on which
// material node 0 is in (Unit 15). NaN on insulator nodes (Unit 15b), which have no bands in the
// model.
struct BandDiagram {
    std::vector<double> conduction_eV;
    std::vector<double> valence_eV;
    std::vector<double> electron_fermi_eV;
    std::vector<double> hole_fermi_eV;
};

// One nonlocal tunnel path at a state (Unit 20): its start node; the length from the start to the
// crossing along the path [cm]; the mean field over it, the potential drop divided by the length
// [V/cm]; the path's generation rate [cm^-3 s^-1], per volume of its start node; and whether the
// crossing lies on the path's frozen segment (false only while the paths lag the state).
struct TunnelPath {
    std::size_t start_node;
    double length_cm;
    double field_V_per_cm;
    double generation_cm3_s;
    bool reached;
};

struct EquilibriumResult {
    NodeFields fields;
    BandDiagram bands;
    // Charge on each gate electrode at zero bias, as BiasPoint::gate_charge.
    std::vector<double> gate_charge;
    // Trapped charge of each declared interface, as BiasPoint::interface_trap_charge.
    std::vector<double> interface_trap_charge;
    ConvergenceRecord convergence;
};

struct BiasPoint {
    std::vector<double> bias_V;  // per contact, in the device's contact order
    NodeFields fields;
    BandDiagram bands;
    // Conventional current entering the device through each contact, in A / cm^(3-D): A/cm^2 in
    // 1D, A/cm in 2D (per unit depth), A in 3D. The currents of all contacts sum to zero; a gate
    // carries none.
    std::vector<double> terminal_current;
    // Per contact, the resolution of terminal_current (same unit): a bound on how far the current
    // through any cut of the device can differ from it, from the continuity residuals and the
    // rounding of the large terms that cancel in them (assemble::DriftDiffusion::
    // terminal_current_resolution). A current below it is not resolved by the (psi, n, p) state,
    // whose densities carry the current only as the small difference of two large fluxes;
    // currents far below it need other unknowns (quasi-Fermi potentials). Zero for a gate and in
    // the quasi-static sweep. (Unit 15.)
    std::vector<double> terminal_current_resolution;
    // Charge on each gate electrode (a lumped gate or a meshed electrode), in C / cm^(3-D)
    // (C/cm^2 in 1D); zero for an ohmic contact. It balances the semiconductor charge and the
    // fixed and trapped interface charge, so the quasi-static capacitance is its derivative with
    // respect to the gate bias.
    std::vector<double> gate_charge;
    // Charge held by the traps of each declared interface (Unit 15b), in C / cm^(3-D): q times the
    // donor-like traps' empty and minus the acceptor-like traps' occupied sheet densities, summed
    // over the interface's area; zero for an interface without traps. The fixed charge is not in
    // it.
    std::vector<double> interface_trap_charge;
    // Electron and hole current through each mesh edge, from Edge::first to Edge::second (same
    // unit).
    std::vector<double> edge_current_n;
    std::vector<double> edge_current_p;
    ConvergenceRecord convergence;
    // Nonlocal band-to-band tunnelling (Unit 20; empty with the model off): each frozen path at
    // this state, and how many times the paths were re-traced before they agreed with the state
    // (the convergence record is the last Newton run's).
    std::vector<TunnelPath> tunnel_paths;
    int path_relocations = 0;
    // Electrothermal (Unit 23; empty with the model off): the temperature of each thermal contact
    // at this point in K (isothermal: its own; R_th: its ambient), and the heat leaving the device
    // through each, in W / cm^(3-D) (W/cm^2 in 1D, W/cm in 2D, W in 3D), in the device's
    // thermal-contact order. In a steady state the heats sum to sum_c V_c I_c over the contacts.
    std::vector<double> thermal_bias_K;
    std::vector<double> thermal_contact_heat;
};

}  // namespace NiTCAD::results
