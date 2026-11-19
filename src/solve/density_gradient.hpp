// Private to the solve layer: density gradient's start (Unit 62, DECISIONS.md QP8, QP9). The thermal
// equilibrium of the model is the stride-3 equilibrium system (psi, Lambda_n, Lambda_p), solved from
// the classical equilibrium with Lambda = 0; when that Newton run fails, a coefficient continuation
// solves at beta k (k = 1/8, 1/4, 1/2, 1), each from the previous solution, a failed rung halving the
// gap to the last converged k at most 4 times. Every result is then read from the drift-diffusion
// state the equilibrium slaves (assemble::System::state_from_quantum_potential).
#pragma once

#include <expected>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/assemble/models.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/assemble/system.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/results/convergence.hpp"
#include "NiTCAD/results/solution.hpp"
#include "NiTCAD/solve/newton.hpp"

namespace NiTCAD::solve::detail {

struct QuantumEquilibrium {
    std::vector<double> psi, lambda_n, lambda_p;  // scaled, per node
    std::vector<double> rungs;                    // QP9's fractions k below 1 that were solved
};

// The density-gradient equilibrium at `bias` (one entry per contact) from the classical equilibrium
// potential psi at the same bias (scaled, per node). `record` holds the last solve's convergence (the
// solve at beta itself). Errors: those of assemble::create_equilibrium and newton_solve; a failed
// ladder is non_convergence, its message listing the rungs that converged.
[[nodiscard]] std::expected<QuantumEquilibrium, base::Error> solve_quantum_equilibrium(
    const device::Device& device, const assemble::Scaling& scaling,
    const assemble::PhysicsModels& models, std::span<const double> bias, std::span<const double> psi,
    const NewtonOptions& newton, linalg::LinearSolver& solver, results::ConvergenceRecord& record,
    const std::stop_token& stop = {}, const IterationObserver& observe = {});

// The classical equilibrium potential (scaled, per node) at `bias` from the charge-neutral one: the
// models with the confinement model off. Errors: as above.
[[nodiscard]] std::expected<std::vector<double>, base::Error> classical_equilibrium(
    const device::Device& device, const assemble::Scaling& scaling,
    const assemble::PhysicsModels& models, std::span<const double> bias, const NewtonOptions& newton,
    linalg::LinearSolver& solver, results::ConvergenceRecord& record,
    const std::stop_token& stop = {}, const IterationObserver& observe = {});

// The run-record warnings (QP9, QP15). With rungs: "density gradient: the equilibrium start used the
// coefficient continuation (k = 0.125, 0.25)"; else empty.
[[nodiscard]] std::string continuation_warning(std::span<const double> rungs);

// QP15's coarse-mesh scan over a run's results: |Lambda| above 10 kT at a node away from the
// interfaces (a semiconductor node with no edge to an insulator node, not a lumped-gate node).
class CoarseScan {
public:
    explicit CoarseScan(const device::Device& device);
    void add(const results::NodeFields& fields);  // no-op without the quantum potentials
    // "density gradient: |Lambda| above 10 kT at 3 node(s) away from the interfaces in 2 of 5
    // result(s), largest 14.2 kT (node 17); the mesh may not resolve the densities there (QP15)", or
    // empty when none.
    [[nodiscard]] std::string warning() const;

private:
    std::vector<char> bulk_;
    double kT_eV_;
    std::vector<char> flagged_;
    std::size_t results_ = 0, flagged_results_ = 0;
    double largest_ = 0.0;
    std::size_t where_ = 0;
};

}  // namespace NiTCAD::solve::detail
