// A dimension sweep's results in the application (Unit 26n, DT7, DECISIONS.md): the runs a sweep wrote
// (deck/dimension_sweep.hpp, one folder per value), their curves overlaid, and Unit 24's figures against the swept
// parameter with their CSV. No Win32.
#pragma once

#include <cstddef>
#include <expected>
#include <filesystem>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/app/workspace.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/render/scene.hpp"

namespace NiTCAD::app {

struct SweepRun {
    std::string value;               // as the sweep wrote it: "20nm", "1e18"
    double numeric = 0.0;            // cm for a length
    std::filesystem::path folder;    // <deck stem>_sweep_<name>/<name>=<value>
    std::optional<Workspace> workspace;  // the run's results, when its files load
    std::string error;               // why not ("not run", or the loader's message)
};

struct SweepResults {
    std::filesystem::path deck;       // the deck swept
    std::filesystem::path directory;  // <deck stem>_sweep_<name>
    std::string name;                 // the parameter
    bool length = false;
    std::vector<SweepRun> runs;       // in the sweep's order
};

// The sweep of `name` over `values` of the deck `text` saved at `path`: each value's folder and, where its run has
// written them, its results: the files the deck's [output] names (sweep, device, mesh), in the folder. Errors: the
// deck does not parse or its parameters do not resolve; no parameter `name`; a value sweep_values refuses for its
// kind.
[[nodiscard]] std::expected<SweepResults, std::string> load_sweep(std::string_view text,
                                                                const std::filesystem::path& path,
                                                                std::string_view name,
                                                                const std::vector<std::string>& values);

// A run's legend: "Lg = 20 nm", "Nd = 1e18" (a length's unit spaced, um as µm).
[[nodiscard]] std::string run_label(const SweepResults& s, const SweepRun& r);

// The overlay (DT7): each loaded run's curve, the one its figures are read from (figure_curve), on shared axes; one
// colour per run in viridis order (the first darkest), each named run_label. Errors (invalid_input): no loaded run
// holds a sweep's curve.
[[nodiscard]] std::expected<render::PlotScene, base::Error> sweep_overlay(const SweepResults& s, bool log_y = false);

// Unit 24's figures (extracted_figures) of the loaded runs: the figure names in the order first seen, one row per
// loaded run (a figure a run lacks is NaN).
struct SweepFigures {
    std::vector<std::string> names;
    std::vector<std::string> units;
    std::vector<std::size_t> runs;          // the row's index in SweepResults::runs
    std::vector<std::vector<double>> rows;  // rows[i][j]: figure j of row i
};
[[nodiscard]] SweepFigures sweep_figures(const SweepResults& s);

// The swept value as the figures show it: a length in nm, a number as it is; its axis title "Lg (nm)" / "Nd".
[[nodiscard]] double shown_value(const SweepResults& s, const SweepRun& r) noexcept;
[[nodiscard]] std::string value_title(const SweepResults& s);

// The figures as CSV: a header (value_title, then "<figure> (<unit>)" or the bare name without a unit), then a line
// per row: the shown value and each figure (shortest round trip), empty where NaN.
[[nodiscard]] std::string figures_csv(const SweepResults& s, const SweepFigures& f);
// <directory>/figures.csv
[[nodiscard]] std::filesystem::path figures_csv_path(const SweepResults& s);

// Figure `j` against the swept value, markers joined in sweep order. Errors (invalid_input): j out of range; fewer
// than one row with the figure.
[[nodiscard]] std::expected<render::PlotScene, base::Error> figure_plot(const SweepResults& s, const SweepFigures& f,
                                                                       std::size_t j);

// The sweep dialog's values (DT6): `typed` read by deck::sweep_values for row `name`'s kind, a bare number typed for a
// length taking the row's unit (as the panel takes it, dims::set_dimension). Errors: no such row, a derived row (a sum
// or a reference: it follows others), sweep_values' errors.
[[nodiscard]] std::expected<std::vector<std::string>, std::string> typed_sweep_values(std::string_view deck_text,
                                                                                    std::string_view name,
                                                                                    std::string_view typed);

}  // namespace NiTCAD::app
