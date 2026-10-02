// Steady-state bias solve (ARCHITECTURE.md section 11, Unit 9): Newton on the coupled
// drift-diffusion system at given contact biases (legacy Device1D::solve_bias, baseline models).
//
// The start is either a given state (the previous bias point of a sweep) or, without one, the
// equilibrium solution; its contact nodes are then set to the biased Dirichlet values, as in the
// legacy. Results are in the public unit convention (V, cm^-3, A / cm^(3-D)).
// DeviceState and BiasSolution are plain data defined here until the results layer exists
// (Unit 10).
#pragma once

#include <expected>
#include <optional>
#include <span>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/solve/newton.hpp"

namespace NiTCAD::solve {

struct BiasOptions {
    NewtonOptions newton;
    linalg::SolverConfig linear;
    std::optional<double> Ns_override;  // see assemble::make_scaling
    assemble::DriftDiffusionModels models;
};

// Per node: electrostatic potential (referenced to the intrinsic level), carrier densities.
struct DeviceState {
    std::vector<double> potential_V;
    std::vector<double> n_cm3;
    std::vector<double> p_cm3;
};

struct BiasSolution {
    assemble::Scaling scaling;
    std::vector<double> bias_V;            // per contact, device.contacts() order
    DeviceState state;
    // Conventional current entering the device through each contact, in A / cm^(3-D):
    // A/cm^2 in 1D, A/cm in 2D (per unit depth), A in 3D. They sum to zero (Kirchhoff).
    std::vector<double> terminal_current;
    // Electron and hole current through each mesh edge (along first -> second), same unit.
    std::vector<double> edge_current_n;
    std::vector<double> edge_current_p;
    NewtonReport newton;
};

// Errors: invalid_input if bias_V does not have one finite value per contact or the initial state
// does not have one positive density pair and a finite potential per node; those of
// solve_equilibrium (when no initial state is given), assemble::DriftDiffusion::create,
// linalg::LinearSolver::create and newton_solve (non_convergence among them).
[[nodiscard]] std::expected<BiasSolution, base::Error> solve_bias(
    const device::Device& device, std::span<const double> bias_V, const BiasOptions& options = {},
    const DeviceState* initial = nullptr);

}  // namespace NiTCAD::solve
