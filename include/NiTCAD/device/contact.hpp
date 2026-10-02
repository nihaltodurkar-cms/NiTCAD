// Contact description (ARCHITECTURE.md 6.4): a named set of boundary nodes with a type.
//
// Data only. The boundary condition a contact imposes is a set of residual and Jacobian rows, so
// it is applied in assemble; the applied bias belongs to the solve layer, not here. Only ohmic
// contacts exist: their boundary value is the charge-neutral equilibrium of the contact node
// (physics::boltzmann_neutral_equilibrium) shifted by the applied bias. Schottky, gate and
// oxide-coupled contacts are deferred.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "NiTCAD/mesh/mesh.hpp"

namespace NiTCAD::device {

enum class ContactKind : std::uint8_t { ohmic };

struct Contact {
    std::string name;
    ContactKind kind;
    std::vector<mesh::NodeId> nodes;  // strictly increasing; each on a mesh boundary patch
};

}  // namespace NiTCAD::device
