// Transient drift-diffusion (ARCHITECTURE.md section 11, Unit 21; legacy transient.py,
// transient2d.py, transient3d.py).
//
// From the steady state at the waveforms' starting values (t = 0) the run steps the coupled
// system to t_end with BDF formulas (assemble::DriftDiffusion::TimeStep): backward Euler, or
// variable-step BDF2 (second order, L-stable, one Newton solve per step). After the start and
// after every waveform breakpoint the integrator restarts with backward Euler (the history
// before a kink or jump is not smooth); BDF2 resumes once three smooth states exist.
// Poisson is algebraic: the potential follows a bias jump at once, the carriers do not.
//
// Step control (adaptive): the local error of an accepted step is estimated from the difference
// between the solution and the polynomial through the earlier states (Milne's device: linear for
// backward Euler, quadratic for BDF2), scaled by the formula's error constant; it is measured on
// the potential in units of kT/q, on each density relative to (density + density_ref) and on each
// trap occupancy, and the step is accepted when the largest is at most rtol. The next step is
// h (0.9 / ratio)^(1 / (order + 1)), within [0.2 h, 2 h] and at most dt_max. A step whose Newton
// solve fails is retried at h / 4. A run whose step falls below dt_min stops (non_convergence).
// Without adaptive control every step is dt_initial. Steps land exactly on every waveform
// breakpoint, output time and t_end.
// OLD / NEW / REASON: the legacy grew the step by 1.5 after a Newton solve of few iterations and
// halved it on failure, with no error control; the number of Newton iterations says nothing about
// the time error, so the steps are error-controlled here.
//
// Terminal currents include the displacement current (OLD / NEW / REASON: the legacy reported
// the conduction current only, which is not conserved in a transient): the charge of each contact
// (assemble::DriftDiffusion::contact_charges) is differenced with the step's own BDF formula, and
// the total currents of all contacts then sum to zero to the solve's tolerance (the header comment
// of drift_diffusion.hpp).
//
// The first step from the start is taken at dt_initial; there is no estimate for it.
// Cancellation and progress follow control.hpp: the starting equilibrium has Phase::equilibrium,
// the starting steady state Phase::bias (point 0 of 1), and every step attempt Phase::transient
// with the attempt's index as point, point_count 0 (unknown in advance) and its end time in
// time_s.
#pragma once

#include <cstdint>
#include <expected>
#include <optional>
#include <span>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/results/run.hpp"
#include "NiTCAD/results/solution.hpp"
#include "NiTCAD/results/transient.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/control.hpp"
#include "NiTCAD/solve/waveform.hpp"

namespace NiTCAD::solve {

enum class Integrator : std::uint8_t { backward_euler, bdf2 };

struct TransientOptions {
    // Newton, linear solver, scaling and models of every solve (equations must be
    // drift_diffusion).
    BiasOptions steady;
    Integrator integrator = Integrator::bdf2;
    double t_end_s = 0.0;                // finite and positive
    std::optional<double> dt_initial_s;  // default t_end / 1e6
    std::optional<double> dt_min_s;      // default dt_initial / 1e6
    std::optional<double> dt_max_s;      // default t_end
    bool adaptive = true;
    double rtol = 1e-3;
    double density_ref_cm3 = 1e10;  // densities far below it are not error-controlled
    // Times in (0, t_end] at which to keep field snapshots, besides 0 and t_end; steps land on
    // them. With fields_every_step, every accepted step keeps one.
    std::vector<double> output_times_s;
    bool fields_every_step = false;
};

// Errors, before anything is solved (invalid_input unless noted): waveforms not one per contact; a
// waveform value at 0 or at a breakpoint rejected by assemble::check_contact_bias (context index:
// the contact); equations not drift_diffusion; t_end, dt_initial, dt_min, dt_max, rtol or
// density_ref not finite and positive (density_ref may be 0), dt_min > dt_initial or dt_initial >
// dt_max; an output time not finite or outside (0, t_end]; an initial state as for sweep_bias;
// those of assemble::make_scaling, assemble::DriftDiffusion::create and
// linalg::LinearSolver::create.
// Once solving has started nothing is an error: the Transient keeps the accepted steps and, if the
// run stopped early, `stopped` and `unfinished` (a stopped starting solve leaves no points).
// `initial` is a starting guess for the steady state at t = 0 (else the run starts from thermal
// equilibrium, as sweep_bias).
// Electrothermal (steady.models.electrothermal, Unit 23; DECISIONS.md T7, T12): the thermal
// contacts keep the temperatures of steady.thermal_bias_K (one list) or the device's for the whole
// run; the heat rows store rho c dT/dt and the carriers' stored energy, and the error estimate
// measures every node's temperature rise off the isothermal sinks (absolute, in units of the
// device's temperature T0, as the potential in V_T). The state at t = 0 is the electrothermal
// steady state, except that a given initial temperature (initial->temperature_K) is the
// temperature at t = 0, held while the electrical rows are solved, and that a device with a part
// without any thermal contact (no steady temperature) starts held at T0; a fully adiabatic device
// is well posed in time. Errors as sweep_bias's electrothermal input checks (thermal_bias_K at most
// one list); a step whose temperature leaves a material's range stops the run (non_convergence).
// Each TimePoint carries the heat leaving through each thermal contact.
[[nodiscard]] std::expected<results::Transient, base::Error> solve_transient(
    const device::Device& device, std::span<const Waveform> waveforms,
    const TransientOptions& options, const results::NodeFields* initial = nullptr,
    const RunControl& control = {});

// The run record of a transient run: the identity digest of the device, the options, the waveforms
// and the initial state, and the options as named settings (make_run_record of bias.hpp, with
// the transient options added).
[[nodiscard]] results::RunRecord make_run_record(const device::Device& device,
                                                 const TransientOptions& options,
                                                 std::span<const Waveform> waveforms,
                                                 const results::NodeFields* initial = nullptr);

}  // namespace NiTCAD::solve
