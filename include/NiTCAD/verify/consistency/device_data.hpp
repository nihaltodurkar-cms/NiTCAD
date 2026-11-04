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
};

struct DeviceData {
    double temperature_K = 0.0;
    std::vector<double> donors;                 // N_D per node [cm^-3]
    std::vector<double> acceptors;              // N_A per node [cm^-3]
    std::vector<double> relative_permittivity;  // per node
    std::vector<std::uint8_t> insulator;        // 1 on insulator nodes
    std::vector<ContactData> contacts;          // device contact order
    bool interface_traps = false;               // any declared interface carries traps
    bool interface_fixed_charge = false;        // any declared interface carries a fixed charge
    std::size_t thermal_contacts = 0;
};

// Field problems against the mesh (sizes, node ranges, finiteness); empty when consistent.
[[nodiscard]] std::vector<std::string> validate(const DeviceData& d, const mesh::Mesh& m);

}  // namespace NiTCAD::verify::consistency
