// Fields and terminal quantities of a solved state (ARCHITECTURE.md 6.6): plain data in the public
// unit convention (V, cm^-3, A / cm^(3-D)), indexed by mesh node, mesh edge or contact.
#pragma once

#include <vector>

#include "NiTCAD/results/convergence.hpp"

namespace NiTCAD::results {

// Per mesh node: electrostatic potential referenced to the intrinsic level, carrier densities.
struct NodeFields {
    std::vector<double> potential_V;
    std::vector<double> n_cm3;
    std::vector<double> p_cm3;
};

struct EquilibriumResult {
    NodeFields fields;
    ConvergenceRecord convergence;
};

struct BiasPoint {
    std::vector<double> bias_V;  // per contact, in the device's contact order
    NodeFields fields;
    // Conventional current entering the device through each contact, in A / cm^(3-D): A/cm^2 in
    // 1D, A/cm in 2D (per unit depth), A in 3D. The currents of all contacts sum to zero.
    std::vector<double> terminal_current;
    // Electron and hole current through each mesh edge, from Edge::first to Edge::second (same
    // unit).
    std::vector<double> edge_current_n;
    std::vector<double> edge_current_p;
    ConvergenceRecord convergence;
};

}  // namespace NiTCAD::results
