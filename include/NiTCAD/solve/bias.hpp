// Steady-state bias solves (ARCHITECTURE.md section 11, Units 9 and 10): Newton on the coupled
// drift-diffusion system at given contact biases (legacy Device1D::solve_bias, baseline models),
// for one bias point or a sweep.
//
// A sweep starts from the given state or, without one, from the equilibrium solution; each point
// starts from the previous one with its contact nodes set to the new Dirichlet values (legacy). One
// system and one linear solver serve the whole sweep, so the pattern is analyzed once.
// Cancellation and progress follow control.hpp: the stop token is checked between points and before
// every Newton iteration; progress events have Phase::equilibrium for the starting equilibrium and
// Phase::bias with the point index. Results are plain data (results/), in V, cm^-3, A / cm^(3-D),
// C / cm^(3-D).
//
// Equations::equilibrium_poisson (Unit 12) solves each point as thermal equilibrium instead:
// Poisson alone (assemble::EquilibriumPoisson), every ohmic contact at 0 V, only gates biased. That
// is the legacy MOS-C solve (moscap.MOSCapacitor.cv_sweep), the quasi-static C-V. A sweep starts
// from the given state or the charge-neutral potential (no starting equilibrium), and reports zero
// currents. Drift-diffusion cannot replace it on a MOS capacitor: an inversion layer with no ohmic
// contact of its own reaches the substrate only through the depleted region, where the minority
// density is some 1e-14 of the inversion layer's, and Newton stalls there above threshold
// (ARCHITECTURE.md 6.4, Unit 12).
#pragma once

#include <cstdint>
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

enum class Equations : std::uint8_t {
    drift_diffusion,      // Poisson and both continuity equations
    equilibrium_poisson,  // Poisson alone, carriers in thermal equilibrium (only gates biased)
};

struct BiasOptions {
    NewtonOptions newton;
    linalg::SolverConfig linear;
    std::optional<double> Ns_override;  // see assemble::make_scaling
    assemble::PhysicsModels models;     // equilibrium_poisson: only bgn matters
    Equations equations = Equations::drift_diffusion;
};

// Solves each bias point in order. `points[k]` holds one bias in V per contact.
// Errors, before anything is solved (invalid_input unless noted): no points; a point rejected by
// assemble::check_contact_bias (one finite value per contact; with equilibrium_poisson, every
// ohmic contact at 0 V), with the point as context index, the contact named in the message and its
// bias as value; an initial state without a finite potential for every node, or, for
// drift-diffusion, without positive densities (the quasi-static sweep reads only the potential);
// those of assemble::make_scaling, assemble::DriftDiffusion::create (or
// EquilibriumPoisson::create) and linalg::LinearSolver::create.
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

// The run record of a sweep: the identity digest of every input the sweep reads (device, options,
// bias points, initial state) and those options as named settings. Inputs the sweep ignores are
// left out, so they do not change the identity: with equilibrium_poisson the mobility, SRH, Auger
// and field-mobility switches and the initial densities; the work function of a polysilicon gate.
[[nodiscard]] results::RunRecord make_run_record(const device::Device& device,
                                                 const BiasOptions& options,
                                                 std::span<const std::vector<double>> points,
                                                 const results::NodeFields* initial = nullptr);

}  // namespace NiTCAD::solve
