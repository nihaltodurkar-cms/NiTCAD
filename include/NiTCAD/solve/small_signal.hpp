// Small-signal (AC) analysis (ARCHITECTURE.md section 11, Unit 22; legacy ac.py, ac2d.py,
// ac3d.py).
//
// At each bias point the run solves the operating point (sweep_bias), then at each frequency f the
// complex linear system of a perturbation e^(i omega t), omega = 2 pi f:
//     (J + i omega t0 C + T(i omega t0)) dx = -dF/dV_j        (assemble::DriftDiffusion)
// once per contact j, with one factorization per frequency. J is the steady Jacobian, C the
// storage term's and T the interface traps' terms with their occupancy's dynamics. The
// admittance Y_ij is the total current entering through contact i, conduction plus i omega times
// its charge, per volt on contact j (results/small_signal.hpp).
// OLD / NEW / REASON: the legacy read the particle current only, through a finite-difference probe
// of the edge current, and its own reciprocity check missed by 1.2% at 1e10 Hz for want of the
// displacement current; here the current includes the displacement current and its partials are
// exact, so the columns and rows of Y sum to zero.
//
// With Equations::equilibrium_poisson the operating point is the quasi-static one (thermal
// equilibrium, only gates biased; the carriers in equilibrium at its potential), and the
// small-signal system is still the drift-diffusion one: so a MOS capacitor in inversion, which the
// drift-diffusion operating point cannot reach (bias.hpp), gives the low-frequency capacitance at
// low frequency and the high-frequency one above the minority carriers' response.
//
// Cancellation and progress follow control.hpp: the operating points as sweep_bias (Phase::bias),
// then one event per frequency of each point with Phase::small_signal, the frequency's index + 1 as
// iteration, the largest backward error as residual and the frequency in frequency_Hz. The stop
// token is checked before each frequency.
#pragma once

#include <expected>
#include <span>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/results/run.hpp"
#include "NiTCAD/results/small_signal.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/control.hpp"

namespace NiTCAD::solve {

struct SmallSignalOptions {
    // The operating points (equations: drift_diffusion or equilibrium_poisson) and the small-signal
    // system's models, scaling and linear solver configuration.
    BiasOptions steady;
    std::vector<double> frequencies_Hz;        // finite and >= 0, at least one
    std::vector<double> field_frequencies_Hz;  // each one of frequencies_Hz: keep the fields there
};

// Errors, before anything is solved (invalid_input unless noted): no frequency, a frequency not
// finite or negative, a field frequency not among the frequencies; those of sweep_bias (points,
// initial state, scaling, assembler, linear solver) and of linalg::ComplexLinearSolver::create.
// Once solving has started nothing is an error: the run keeps the completed points and, if it
// stopped early (an operating point not converging, a singular small-signal system, cancelled),
// `stopped` and `unfinished`.
[[nodiscard]] std::expected<results::SmallSignal, base::Error> solve_small_signal(
    const device::Device& device, std::span<const std::vector<double>> points,
    const SmallSignalOptions& options, const results::NodeFields* initial = nullptr,
    const RunControl& control = {});

// The run record: the identity digest of the device, the options (the small-signal system reads
// the transport models even when the operating point is quasi-static), the points, the
// frequencies and the initial state.
[[nodiscard]] results::RunRecord make_run_record(const device::Device& device,
                                                 const SmallSignalOptions& options,
                                                 std::span<const std::vector<double>> points,
                                                 const results::NodeFields* initial = nullptr);

}  // namespace NiTCAD::solve
