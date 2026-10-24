// MLDA quantum confinement on a device (Unit 47; DECISIONS.md QD2, QD3, QD6): each semiconductor
// node's distance to the nearest insulator interface and its two confinement factors, which the
// assemblers apply to the densities of states (scaled_device.cpp, QD4).
//
// Interface points: the midpoint of every edge joining a semiconductor node and an insulator node
// (where the assembler puts that interface; interface_edges.hpp), and every lumped-gate node
// itself (the gate's oxide is on the node's face; gate.hpp). A node's distance is to the nearest
// such point, so at a corner the nearer interface decides. Its factor for carrier c is
// physics::mlda_factor(d, lambda_c); a lumped-gate node (d = 0) takes instead the mean of f over
// [0, l / 2], l its shortest edge to another semiconductor node (physics::mlda_mean_factor).
#pragma once

#include <expected>
#include <vector>

#include "NiTCAD/assemble/models.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"

namespace NiTCAD::assemble {

struct ConfinementFactors {
    double wavelength_n_cm, wavelength_p_cm;  // lambda at the device temperature
    // Per node (insulator nodes: distance 0, factors 1, unused).
    std::vector<double> distance_cm;
    std::vector<char> half_cell;  // 1 at a lumped-gate node (the factors are QD3's means)
    std::vector<double> f_n, f_p;
};

// Errors (invalid_input): a mass not finite and positive; electrothermal on; a semiconductor
// region whose electron or hole Caughey-Thomas parameters are not silicon's; a thermionic-emission
// interface (its emission term takes the materials' densities of states); no insulator node and no
// lumped gate (nothing confines); a lumped-gate node with no edge to another semiconductor node.
// models.quantum_confinement is not read.
[[nodiscard]] std::expected<ConfinementFactors, base::Error> mlda_factors(
    const device::Device& device, const PhysicsModels& models);

}  // namespace NiTCAD::assemble
