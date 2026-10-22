// The results the application shows (Unit 26, DECISIONS.md VA4, VA5): a device, its drawing
// geometry and the results of one run, from a run just finished, a deck's output files or a result
// file with its device (and generated mesh) file, paired exactly as `nitcad plot` pairs them (the
// mesh digest, Unit 45; the node count for older files). Also what the results browser lists: the
// run record and provenance (Unit 42), the convergence history and Unit 24's extracted figures.
//
// A Workspace owns its device and keeps it at a fixed address (the geometry refers to its mesh), so
// it may be moved but not copied.
#pragma once

#include <cstddef>
#include <expected>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/render/geometry.hpp"
#include "NiTCAD/results/run.hpp"
#include "NiTCAD/results/small_signal.hpp"
#include "NiTCAD/results/solution.hpp"
#include "NiTCAD/results/transient.hpp"

namespace NiTCAD::app {

struct Workspace {
    std::filesystem::path source;  // the deck or result file shown
    std::unique_ptr<device::Device> device;
    std::unique_ptr<mesh::GeneratedMesh> generated;  // a generated mesh's simplices, if any
    std::unique_ptr<render::Geometry> geometry;
    std::optional<deck::Analysis> analysis;  // when from a deck
    std::optional<results::Sweep> sweep;      // a sweep or a trace
    std::optional<results::Transient> transient;
    std::optional<results::SmallSignal> small_signal;
    std::optional<results::BiasPoint> equilibrium;
    // The bias points the field views read, as `nitcad plot` lists them: a sweep's points, a
    // transient's snapshots (with their currents and rates), a small-signal run's operating points,
    // or the equilibrium.
    std::vector<results::BiasPoint> points;

    [[nodiscard]] bool has_results() const noexcept { return !points.empty(); }
    // The run record of the results, if they have one (an equilibrium has none).
    [[nodiscard]] const results::RunRecord* record() const noexcept;
};

// A finished run's plan and report (nothing read from disk).
[[nodiscard]] std::expected<Workspace, base::Error> workspace_from_run(deck::Plan plan, deck::RunReport report);
// A deck's device, and the results its outputs hold if they exist (none before the deck has run).
[[nodiscard]] std::expected<Workspace, base::Error> workspace_from_deck(const std::filesystem::path& deck);
// A .nitsweep, .nittran or .nitac file with its device file and, for a generated mesh, its mesh file.
[[nodiscard]] std::expected<Workspace, base::Error> workspace_from_result(const std::filesystem::path& result,
                                                                         const std::filesystem::path& device,
                                                                         const std::filesystem::path& mesh = {});

// The results browser (VA5): a name and a value per line.
struct Line {
    std::string name;
    std::string value;
};
// The run record (identity, mesh digest, settings, how it ended) and the build and platform of the
// provenance; for an equilibrium, its convergence.
[[nodiscard]] std::vector<Line> run_record(const Workspace& w);
// Every converged point's Newton history: point, iterations, final update.
struct PointHistory {
    std::size_t point = 0;
    std::vector<results::IterationRecord> iterations;
    bool converged = false;
};
[[nodiscard]] std::vector<PointHistory> convergence_history(const Workspace& w);

// Unit 24's figures where they apply: for a sweep whose swept contact is a gate or electrode, the
// transfer curve's threshold at maximum transconductance, the subthreshold swing, the maximum
// transconductance and the on/off ratio between the first and last gate biases (the sensed
// contact the other one carrying the most current); for a sweep of an ohmic contact, the current
// at its last bias. A figure that cannot be extracted is left out.
struct Figure {
    std::string name;
    std::string curve;  // the curve it was read from, e.g. "I(drain) vs V(gate)"
    double value = 0.0;
    std::string unit;
    bool extrapolated = false;
    bool below_resolution = false;
};
[[nodiscard]] std::vector<Figure> extracted_figures(const Workspace& w);

// The contact whose bias varies most over the points (the swept one), as `nitcad plot` chooses it.
[[nodiscard]] std::size_t swept_contact(const std::vector<std::vector<double>>& biases);

}  // namespace NiTCAD::app
