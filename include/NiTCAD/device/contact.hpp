// Contact description (ARCHITECTURE.md 6.4): a named set of boundary nodes with a type.
//
// Data only. The boundary condition a contact imposes is a set of residual and Jacobian rows, so
// it is applied in assemble; the applied bias belongs to the solve layer, not here.
// - Ohmic: the contact node is held at its charge-neutral equilibrium
//   (physics::boltzmann_neutral_equilibrium), its potential shifted by the applied bias.
// - Gate (Unit 12): a gate electrode on a lumped oxide (legacy moscap.MOSCapacitor and
//   Device2D/3D.add_gate). The oxide is not meshed: it is a parallel-plate capacitor of thickness
//   t_ox between the electrode and each contact node, so the node's Poisson row gains the oxide
//   displacement flux (Gauss's law across the interface, a Robin condition on the potential) and
//   the fixed oxide charge. No carriers cross it: the continuity rows keep zero boundary flux, and
//   a gate carries no DC current. Schottky contacts are deferred.
// - Electrode (Unit 15b): a gate electrode on a meshed insulator. It is electrostatic only: its
//   nodes, which lie in an insulator region, hold the potential the electrode's Fermi level
//   imposes (a Dirichlet condition on psi); they have no carriers and carry no current. Its charge
//   is the displacement flux leaving it (assemble).
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "NiTCAD/mesh/mesh.hpp"

namespace NiTCAD::device {

enum class ContactKind : std::uint8_t { ohmic, gate, electrode };

// The gate electrode's work function phi_m (legacy moscap.flatband_voltage): a degenerately doped
// polysilicon gate has its Fermi level at the band edge of the semiconductor under it, so
// phi_m = chi (n+) or chi + Eg(T) (p+); a metal gate gives phi_m directly.
enum class GateElectrode : std::uint8_t { n_poly, p_poly, metal };

struct GateStack {
    // The mesh boundary patch the gate lies on. Each contact node must be on it, and the patch's
    // face area of the node is the gate area that node couples through.
    std::string boundary;
    double oxide_thickness_cm = 0.0;           // t_ox, finite and positive
    double oxide_relative_permittivity = 3.9;  // SiO2 (legacy EPS_OX_R), finite and positive
    GateElectrode electrode = GateElectrode::n_poly;
    double work_function_eV = 0.0;             // phi_m of a metal electrode; finite and positive
    double fixed_charge_cm2 = 0.0;             // Q_f: fixed positive oxide charge at the interface
};

// The electrode of an `electrode` contact. A polysilicon electrode is degenerately doped silicon:
// phi_m = chi or chi + Eg(T) of physics::silicon_parameters (it has no semiconductor under it to
// take them from, unlike a lumped gate); a metal gives phi_m.
struct Electrode {
    GateElectrode kind = GateElectrode::n_poly;
    double work_function_eV = 0.0;  // phi_m of a metal; finite and positive
};

struct Contact {
    std::string name;
    ContactKind kind;
    std::vector<mesh::NodeId> nodes;  // strictly increasing; each on a mesh boundary patch
    GateStack gate{};                 // read only when kind == ContactKind::gate
    Electrode electrode{};            // read only when kind == ContactKind::electrode
};

// A thermal contact (Unit 23, the electrothermal model; legacy thermal.ThermalBC): a heat sink on
// a mesh boundary patch, separate from the electrical contacts (a patch may carry both).
// - isothermal: every node of the patch held at temperature_K (Dirichlet; legacy "isothermal");
// - resistance: the heat flux out through each node's face area A is (T - temperature_K) / R_th
//   (Robin; legacy "resistance"), R_th = resistance_K_cm2_W per unit boundary area in every
//   dimension, temperature_K the ambient behind it.
// Every boundary without a thermal contact is adiabatic (no heat flux; legacy "adiabatic"). The
// temperatures are those of a bias point unless the solve gives others (solve/bias.hpp).
enum class ThermalContactKind : std::uint8_t { isothermal, resistance };

struct ThermalContact {
    std::string name;
    std::string boundary;  // the mesh boundary patch
    ThermalContactKind kind = ThermalContactKind::isothermal;
    double temperature_K = 300.0;
    double resistance_K_cm2_W = 0.0;  // read only when kind == resistance
};

}  // namespace NiTCAD::device
