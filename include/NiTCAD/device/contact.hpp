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
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "NiTCAD/mesh/mesh.hpp"

namespace NiTCAD::device {

enum class ContactKind : std::uint8_t { ohmic, gate };

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

struct Contact {
    std::string name;
    ContactKind kind;
    std::vector<mesh::NodeId> nodes;  // strictly increasing; each on a mesh boundary patch
    GateStack gate{};                 // read only when kind == ContactKind::gate
};

}  // namespace NiTCAD::device
