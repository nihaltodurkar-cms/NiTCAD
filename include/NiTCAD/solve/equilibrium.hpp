// Thermal-equilibrium solve (ARCHITECTURE.md section 11, Unit 8): Newton on the equilibrium Poisson
// system from the charge-neutral initial guess (legacy Device1D::solve_equilibrium_boltzmann).
//
// The result is plain data in the public unit convention (V, cm^-3; results/solution.hpp). The
// potential is the electrostatic potential referenced to the intrinsic level, so
// n = n_ie e^(phi/V_T). Progress events have Phase::equilibrium, point 0 of 1 (control.hpp).
#pragma once

#include <expected>
#include <optional>

#include "NiTCAD/assemble/models.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/results/solution.hpp"
#include "NiTCAD/solve/control.hpp"
#include "NiTCAD/solve/newton.hpp"

namespace NiTCAD::solve {

struct EquilibriumOptions {
    NewtonOptions newton;
    linalg::SolverConfig linear;
    std::optional<double> Ns_override;  // see assemble::make_scaling
    assemble::PhysicsModels models;     // only bgn matters at equilibrium
};

// Errors: those of assemble::make_scaling, assemble::EquilibriumPoisson::create,
// linalg::LinearSolver::create and newton_solve (non_convergence and cancelled among them).
[[nodiscard]] std::expected<results::EquilibriumResult, base::Error> solve_equilibrium(
    const device::Device& device, const EquilibriumOptions& options = {},
    const RunControl& control = {});

}  // namespace NiTCAD::solve
