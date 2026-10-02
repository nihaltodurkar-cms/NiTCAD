// Steady-state bias solves (ARCHITECTURE.md section 11, Units 9 and 10): Newton on the coupled
// drift-diffusion system at given contact biases (legacy Device1D::solve_bias, baseline models),
// for one bias point or a sweep.
//
// A sweep starts from the given state or, without one, from the equilibrium solution; each point
// starts from the previous one with its contact nodes set to the new Dirichlet values (legacy). One
// system and one linear solver serve the whole sweep, so the pattern is analyzed once.
// Cancellation and progress follow control.hpp: the stop token is checked between points and before
// every Newton iteration; progress events have Phase::equilibrium for the starting equilibrium and
// Phase::bias with the point index. Results are plain data (results/), in V, cm^-3, A / cm^(3-D).
#pragma once

#include <expected>
#include <optional>
#include <span>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/results/run.hpp"
#include "NiTCAD/results/solution.hpp"
#include "NiTCAD/solve/control.hpp"
#include "NiTCAD/solve/newton.hpp"

namespace NiTCAD::solve {

struct BiasOptions {
    NewtonOptions newton;
    linalg::SolverConfig linear;
    std::optional<double> Ns_override;  // see assemble::make_scaling
    assemble::DriftDiffusionModels models;
};

// Solves each bias point in order. `points[k]` holds one bias in V per contact.
// Errors, before anything is solved (invalid_input unless noted): no points; a point without one
// finite value per contact; an initial state without a finite potential and positive densities for
// every node; those of assemble::make_scaling, assemble::DriftDiffusion::create and
// linalg::LinearSolver::create.
// Once solving has started, nothing is an error: the Sweep holds the completed points and, if the
// run stopped early, `stopped` (cancelled, non_convergence, singular_system, ...) and the stopped
// point's convergence history in `unfinished` (empty if the starting equilibrium stopped).
[[nodiscard]] std::expected<results::Sweep, base::Error> sweep_bias(
    const device::Device& device, std::span<const std::vector<double>> points,
    const BiasOptions& options = {}, const results::NodeFields* initial = nullptr,
    const RunControl& control = {});

// One bias point: a sweep of one, with a stopped run returned as its error.
[[nodiscard]] std::expected<results::BiasPoint, base::Error> solve_bias(
    const device::Device& device, std::span<const double> bias_V, const BiasOptions& options = {},
    const results::NodeFields* initial = nullptr, const RunControl& control = {});

// The run record of a sweep: the identity digest of every input (device, options, bias points,
// initial state) and the options as named settings.
[[nodiscard]] results::RunRecord make_run_record(const device::Device& device,
                                                 const BiasOptions& options,
                                                 std::span<const std::vector<double>> points,
                                                 const results::NodeFields* initial = nullptr);

}  // namespace NiTCAD::solve
