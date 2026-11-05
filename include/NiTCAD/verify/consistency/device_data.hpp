// The device as plain data for the consistency checks (Unit 55, TVR-1, DECISIONS 11.2, 11.6; TVD-19, TVD-20).
//
// The consistency layer judges results with the results, the mesh and the device only, never re-evaluating
// NiTCAD's physics. device::Device's header includes physics headers (its regions hold physics::Semiconductor),
// so this layer does not include it: the benchmark layer copies what the checks read into this struct
// (verify/benchmark/device_data.hpp), and the include boundary (tests/verify/boundary.cmake and
// tests/verify/consistency/header_boundary_test.cpp) holds transitively, not only for direct includes.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "NiTCAD/mesh/mesh.hpp"

namespace NiTCAD::verify::consistency {

enum class ContactKind : std::uint8_t { ohmic, gate, electrode, schottky };

struct ContactData {
    std::string name;
    ContactKind kind = ContactKind::ohmic;
    std::vector<mesh::NodeId> nodes;
    // A lumped gate (kind gate) only: each node's face area on the gate's patch [cm^(D-1)], in parallel with
    // nodes; the oxide's capacitance per area eps_0 eps_ox / t_ox [F/cm^2]; its fixed oxide charge Q_f [cm^-2].
    std::vector<double> gate_areas;
    double oxide_capacitance_F_cm2 = 0.0;
    double gate_fixed_charge_cm2 = 0.0;
    // A gate or electrode that exchanges carriers with the semiconductor by gate tunnelling (Unit 51): its
    // conduction current is not 0 (C6).
    bool tunnelling = false;
};

// A declared interface between two regions (device order, as results::BiasPoint::interface_trap_charge).
struct InterfaceData {
    std::int32_t region_a = -1;
    std::int32_t region_b = -1;
    // Whether its insulator-semiconductor edges carry the interface model (a fixed charge, traps or surface
    // recombination); such an edge's flux is split at an interface potential the results do not store.
    bool charged = false;
    double fixed_charge_cm2 = 0.0;   // sigma_f [cm^-2], the charge is q sigma_f
    double trap_density_cm2 = 0.0;   // total trap sheet density (levels, and bands' D_it times their width)
};

struct DeviceData {
    double temperature_K = 0.0;
    std::vector<double> donors;                 // N_D per node [cm^-3]
    std::vector<double> acceptors;              // N_A per node [cm^-3]
    std::vector<double> relative_permittivity;  // per node
    std::vector<std::uint8_t> insulator;        // 1 on insulator nodes
    std::vector<std::int32_t> region;           // per node, indexing the device's regions
    std::vector<ContactData> contacts;          // device contact order
    std::vector<InterfaceData> interfaces;      // device interface order
    bool interface_traps = false;               // any declared interface carries traps
    bool interface_fixed_charge = false;        // any declared interface carries a fixed charge
    bool incomplete_ionization = false;         // the run's model (the results store no ionized densities)
    bool electrothermal = false;                // the run's model (C5's predicate)
    // The thermopower exponents r of each node's semiconductor (tau ~ E^r; 0 on insulator nodes), a material
    // constant C5's tolerance reads (the energy a carrier carries, E_c + (r + 5/2) k T); empty when not copied.
    std::vector<double> thermopower_n;
    std::vector<double> thermopower_p;
    std::size_t thermal_contacts = 0;
};

// Field problems against the mesh (sizes, node ranges, finiteness); empty when consistent.
[[nodiscard]] std::vector<std::string> validate(const DeviceData& d, const mesh::Mesh& m);

}  // namespace NiTCAD::verify::consistency
