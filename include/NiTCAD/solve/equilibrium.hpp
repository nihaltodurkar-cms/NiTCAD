// Thermal-equilibrium solve (ARCHITECTURE.md section 11, Unit 8): Newton on the equilibrium Poisson
// system from the charge-neutral initial guess (legacy Device1D::solve_equilibrium_boltzmann).
//
// Results are returned in the public unit convention (V, cm^-3), not in scaled units (6.2). The
// potential is the electrostatic potential referenced to the intrinsic level, so
// n = n_ie e^(phi/V_T).
// EquilibriumSolution is plain data defined here until the results layer exists (Unit 10).
#pragma once

#include <expected>
#include <optional>
#include <vector>

#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/solve/newton.hpp"

namespace NiTCAD::solve {

struct EquilibriumOptions {
    NewtonOptions newton;
    linalg::SolverConfig linear;
    std::optional<double> Ns_override;  // see assemble::make_scaling
};

struct EquilibriumSolution {
    assemble::Scaling scaling;
    std::vector<double> potential_V;  // per node
    std::vector<double> n_cm3;
    std::vector<double> p_cm3;
    NewtonReport newton;
};

// Errors: those of assemble::make_scaling, assemble::EquilibriumPoisson::create,
// linalg::LinearSolver::create and newton_solve (non_convergence among them).
[[nodiscard]] std::expected<EquilibriumSolution, base::Error> solve_equilibrium(
    const device::Device& device, const EquilibriumOptions& options = {});

}  // namespace NiTCAD::solve
