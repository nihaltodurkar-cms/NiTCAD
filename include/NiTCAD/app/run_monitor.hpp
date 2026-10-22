// The run monitor's content (Unit 26, DECISIONS.md VA3, page Figure 5): a row per solver point, the live
// convergence plot and the progress of a run, from the RunController's state on the UI thread.
#pragma once

#include <expected>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/app/run_controller.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/render/scene.hpp"

namespace NiTCAD::app {

// The table's columns and a row's cells: number, phase, the point (bias point k of n, the transient time,
// the small-signal frequency), Newton iterations, last update, last residual, converged, seconds.
[[nodiscard]] std::span<const std::string_view> monitor_columns() noexcept;
[[nodiscard]] std::vector<std::string> monitor_cells(const PointRow& row, std::size_t index);

// Every Newton iteration's largest scaled update against its number on a log axis (the iterations with a
// zero update, converged exactly, left out). Errors (invalid_input): fewer than two iterations to plot.
[[nodiscard]] std::expected<render::PlotScene, base::Error> convergence_plot(std::span<const Iteration> iterations);

// The fraction of the run done, in [0, 1]: converged bias points over the planned ones; -1 when the run
// does not know its length in advance (transient, trace) or nothing is planned.
[[nodiscard]] double progress_fraction(std::span<const PointRow> rows, std::size_t planned) noexcept;

}  // namespace NiTCAD::app
