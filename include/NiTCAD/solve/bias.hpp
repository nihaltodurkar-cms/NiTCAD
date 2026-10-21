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

// Automatic bias step cutting (Unit 28). When a point does not converge from the previous
// converged state, the bias step to it is halved, from the previous converged state, until a
// sub-step converges; each converged sub-step lets the next grow by `growth`, up to the rest of
// the step. The point fails only when a step below 2^-max_cuts of the full one fails. Every
// contact's bias (and, with the electrothermal model, every thermal contact's temperature) moves
// linearly along the step. The sub-step states are not points of the sweep: only the requested
// point is recorded, with its sub-step count in BiasPoint::bias_substeps and its last Newton run
// in its convergence.
// Off by default (max_cuts = 0): a point that does not converge stops the sweep, as before.
// What a cut retries: non_convergence, singular_system and inaccurate_solve; cancellation and
// resource_exhausted stop at once. The first point steps from the bias of the starting state:
// zero on every contact for a drift-diffusion sweep started from equilibrium, initial_bias_V for
// a sweep given an initial state; without one (initial_bias_V empty, or the quasi-static sweep
// started from the charge-neutral potential) it is not cut.
struct StepControl {
    int max_cuts = 0;     // 0: off; at most 52 (a step of 2^-52 is at the rounding of the bias)
    double growth = 2.0;  // the factor a sub-step grows by after converging; at least 1
    // The Newton iteration limit of every attempt (the full step included), when on: a diverging
    // attempt is cut after these rather than after newton.max_iterations (the smaller applies).
    int max_iterations = 25;
    // The contact biases of the initial state given to the sweep, if known (empty: unknown).
    std::vector<double> initial_bias_V;
};

struct BiasOptions {
    NewtonOptions newton;
    StepControl step_control;
    linalg::SolverConfig linear;
    std::optional<double> Ns_override;  // see assemble::make_scaling
    assemble::PhysicsModels models;     // equilibrium_poisson: bgn, fermi_dirac and
                                        // incomplete_ionization matter
    Equations equations = Equations::drift_diffusion;
    // The weighted terminal current and its resolution in every point (Unit 32;
    // results::BiasPoint::weighted_current; DECISIONS.md W1-W8). Off by default: one more
    // factorization and a linear solve on the node graph per point. Drift-diffusion only; not with
    // nonlocal tunnelling (its pairs are not carried along the edges).
    bool weighted_current = false;
    // The generation-recombination rates of every drift-diffusion point and transient snapshot
    // (Unit 45; results::RateFields; DECISIONS.md RQ1): on by default, off to save memory and file
    // size. They are read from the state after it converged and change nothing in the solve.
    bool store_rates = true;
    // Electrothermal (models.electrothermal, Unit 23): per bias point, the temperature of each
    // thermal contact in K, in the device's thermal-contact order; empty: the device's own at
    // every point. Read only with the model on.
    std::vector<std::vector<double>> thermal_bias_K;
};

// Solves each bias point in order. `points[k]` holds one bias in V per contact.
// Errors, before anything is solved (invalid_input unless noted): no points; step_control with
// max_cuts outside 0..52, growth not finite or below 1, max_iterations below 1, or
// initial_bias_V neither empty nor one finite value per contact; a point rejected by
// assemble::check_contact_bias (one finite value per contact; with equilibrium_poisson, every
// ohmic contact at 0 V), with the point as context index, the contact named in the message and its
// bias as value; an initial state without a finite potential for every node, or, for
// drift-diffusion, without positive densities (the quasi-static sweep reads only the potential);
// those of assemble::make_scaling, assemble::DriftDiffusion::create (or
// EquilibriumPoisson::create) and linalg::LinearSolver::create. With models.electrothermal
// (Unit 23; DECISIONS.md T7, T10, T12): equilibrium_poisson; a connected part of the device with
// no thermal contact (a steady state has no heat sink there); thermal_bias_K neither empty nor
// one list per point, a list without one temperature per thermal contact, or a temperature not
// finite and positive or outside the range of a semiconductor region's models
// (physics::check_temperature), with the point as context index; an initial state whose
// temperature_K is neither empty nor one finite, positive value per node.
// With the model on, a point whose converged temperature lies outside a node's material range
// stops the run with non_convergence (context: the node and its temperature).
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
// bias points, initial state), those options as named settings, and the provenance (provenance.hpp;
// not in the identity; wall_s 0, which sweep_bias sets when it returns). Inputs the sweep ignores are
// left out, so they do not change the identity: with equilibrium_poisson the mobility, SRH, Auger,
// field-mobility and radiative switches, the interfaces' transport and the initial densities; the
// work function of a polysilicon gate.
[[nodiscard]] results::RunRecord make_run_record(const device::Device& device,
                                                 const BiasOptions& options,
                                                 std::span<const std::vector<double>> points,
                                                 const results::NodeFields* initial = nullptr);

}  // namespace NiTCAD::solve
