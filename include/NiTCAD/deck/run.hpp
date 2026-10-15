// Running a plan (ARCHITECTURE.md 11.3, Unit 31; DECISIONS.md DK9, DK10): the batch runner behind
// `nitcad run`.
//
// The analysis is solve::solve_equilibrium or solve::sweep_bias with the plan's options; then every
// output the plan names is written (io/files.hpp, io/device_file.hpp), the device and mesh first.
// A sweep that stops early (non-convergence, cancellation) is still a result: its completed points
// are written, its checkpoint is the last completed point, and the report says why it stopped.
// A restart (DK9) loads the checkpoint against the device's mesh, finds the plan point whose bias
// equals the checkpoint's (to 1e-12 relative, 1e-12 V near zero), and sweeps the points after it
// from the checkpoint's state, with step_control.initial_bias_V set to the checkpoint's bias; the
// sweep file then holds those points only.
#pragma once

#include <expected>
#include <filesystem>
#include <optional>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/results/run.hpp"
#include "NiTCAD/results/solution.hpp"
#include "NiTCAD/solve/control.hpp"

namespace NiTCAD::deck {

struct RunReport {
    std::optional<results::EquilibriumResult> equilibrium;
    std::optional<results::Sweep> sweep;
    std::size_t skipped_points = 0;               // plan points before a restart's checkpoint
    std::vector<std::filesystem::path> written;   // in the order written
};

// Errors: invalid_input for a plan the solvers reject before solving, a restart checkpoint of
// another mesh or whose bias is not a plan point (or is the last one), an MKL runtime that does
// not load; non_convergence (and the rest) when the equilibrium solve fails; resource_exhausted
// when an output cannot be written (outputs already written stay).
[[nodiscard]] std::expected<RunReport, base::Error> run(const Plan& plan,
                                                        const solve::RunControl& control = {});

// The CSV table (Outputs::table): a header "point,<contact> V...,<contact> A...,substeps" (with the
// weighted current, Unit 32, "<contact> weighted A..." and "<contact> weighted resolution A..."
// before substeps) and one row per point, every number with 17 significant digits.
[[nodiscard]] std::expected<void, base::Error> write_table(const std::filesystem::path& path,
                                                           const device::Device& device,
                                                           const results::Sweep& sweep);

}  // namespace NiTCAD::deck
