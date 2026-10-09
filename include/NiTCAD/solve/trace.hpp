// Pseudo-arclength continuation of the drift-diffusion steady state in one contact's bias
// (ARCHITECTURE.md section 11, Unit 19; legacy continuation.py arc_length_sweep, Keller 1977).
//
// A bias sweep steps the bias and solves at each value; where the current-voltage curve folds back
// (avalanche snapback, a thyristor-like turn-over) no solution lies ahead along the bias, and a
// sweep fails or jumps to another branch. The trace parameterizes the branch by its arc length s in
// (state, bias) instead: each step predicts y + ds tau along the unit tangent tau and corrects with
// Newton on the bordered system
//
//     [ J      dF/dlambda ] [dx      ]   [ -F                          ]
//     [ c^T               ] [dlambda ] = [ ds - c . (y - y_previous)  ],
//
// with c = M tau (M the arc metric below) and lambda = V / V_T the swept bias; it stays regular at
// a turning point, where J alone is singular. Its dense last row would fill the factors, so it is
// solved through a sparse B0 holding only the row's largest entry (kept while it stays within 0.3
// of the largest, so B0's pattern rarely changes) and a Sherman-Morrison correction for the rest.
// The residual is exactly linear in the bias: F(x, lambda) = F(x, lambda_0) + dF/dlambda
// (lambda - lambda_0). With the electrothermal model (Unit 23) dF/dlambda depends on the state (the
// metal's Peltier heat at the swept contact's nodes), so the rows are evaluated at lambda itself
// with dF/dlambda at x; the thermal contacts stay at the starting point's temperatures
// (steady.thermal_bias_K, one list, or the device's), and the temperature rise in units of T0
// enters the state part of the arc metric like the potential in V_T. A thermal-runaway fold of a
// device with a thermal resistance is a turning point like any other. The next tangent solves the
// same matrix with the right-hand side (0, 1), so it keeps its orientation through a fold.
// The arc length is measured in volts: the swept bias; the swept contact's current at 1 V per
// decade of |I| + I0 (I0 ten times the starting point's current resolution); and, weighted by 0.1,
// the state (the potential in V_T and each density relative to itself above density_floor, times
// V_T, as root mean squares over the nodes). Through a snapback the current grows by decades while
// the bias hardly moves, and the current term keeps such steps short.
// A step whose new tangent turns by more than acos(0.9) from the last (the corrector found another
// branch nearby) is rejected and retried at half. A step whose corrector converges in at most 7
// iterations lets the next grow by 1.5 (to max_step_V); one taking 12 or more shrinks it by 0.7; a
// failed corrector is retried at half the step, down to min_step_V, below which the trace stops
// (non_convergence) and keeps its points.
// The trace ends at end_V (the last step lands on it exactly), when the swept contact's current
// reaches current_limit in magnitude, or after max_points. The bias where a criterion current is
// reached is a curve reading; the trace defines no breakdown voltage of its own (the corrector
// failing to converge says nothing about the device).
// Cancellation and progress follow control.hpp: the starting point as sweep_bias (Phase::bias,
// point 0 of 1), then every step attempt with Phase::bias, the attempt's index + 1 as point and
// point_count 0 (unknown in advance).
#pragma once

#include <cstddef>
#include <expected>
#include <optional>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/results/run.hpp"
#include "NiTCAD/results/solution.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/control.hpp"

namespace NiTCAD::solve {

struct TraceOptions {
    BiasOptions steady;            // equations must be drift_diffusion
    std::size_t contact = 0;       // the swept contact
    std::vector<double> start_V;   // every contact's bias at the start
    double end_V = 0.0;            // the swept contact's bias at the end
    double step_V = 0.05;          // the first step (along the bias) [V]
    double max_step_V = 1.0;       // arc length [V]
    double min_step_V = 1e-7;
    std::optional<double> current_limit;  // |I| of the swept contact [A/cm^(3-D)]
    std::size_t max_points = 10000;
    double density_floor_cm3 = 1e10;
};

// The points along the branch in trace order, the starting point first, as a Sweep (each
// BiasPoint's bias the full contact vector; convergence the corrector's).
// Errors, before anything is solved (invalid_input): equations not drift_diffusion; contact out of
// range; start_V not one finite bias per contact (and those of check_contact_bias), end_V not
// finite or equal to the start; a step not finite and positive or min_step_V > step_V > max_step_V;
// current_limit not finite and positive; max_points 0; density_floor negative; an initial state as
// for sweep_bias (and its electrothermal checks for the starting point); those of make_scaling,
// DriftDiffusion::create and LinearSolver::create. With electrothermal, an accepted point whose
// temperature leaves a material's range stops the trace (non_convergence, as sweep_bias).
[[nodiscard]] std::expected<results::Sweep, base::Error> trace_bias(
    const device::Device& device, const TraceOptions& options,
    const results::NodeFields* initial = nullptr, const RunControl& control = {});

[[nodiscard]] results::RunRecord make_run_record(const device::Device& device,
                                                 const TraceOptions& options,
                                                 const results::NodeFields* initial = nullptr);

}  // namespace NiTCAD::solve
