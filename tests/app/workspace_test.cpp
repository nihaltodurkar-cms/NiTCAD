// The results workspace (Unit 26, VA4, VA5): results from a run, a deck's outputs and result files,
// paired with their mesh; the run record; Unit 24's figures.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <filesystem>
#include <string>

#include "NiTCAD/app/run_controller.hpp"
#include "NiTCAD/app/workspace.hpp"
#include "NiTCAD/deck/run.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
namespace fs = std::filesystem;
using Catch::Matchers::WithinRel;

namespace {

// A deck run by the controller; its plan and report.
app::RunController::Taken run(const fs::path& deck, const std::string& text) {
    app::RunController rc;
    REQUIRE(rc.start_run(text, app_test::write(deck.parent_path(), deck.filename().string(), text)));
    rc.wait();
    rc.process();
    REQUIRE(rc.state() == app::RunState::done);
    return std::move(*rc.take_results());
}

const app::Line* find(const std::vector<app::Line>& lines, const std::string& name) {
    const auto it = std::ranges::find_if(lines, [&](const app::Line& l) { return l.name == name; });
    return it == lines.end() ? nullptr : &*it;
}

const app::Figure* figure(const std::vector<app::Figure>& f, const std::string& name) {
    const auto it = std::ranges::find_if(f, [&](const app::Figure& x) { return x.name == name; });
    return it == f.end() ? nullptr : &*it;
}

}  // namespace

TEST_CASE("app workspace: from a run, from the deck's outputs, from result files", "[app]") {
    const fs::path dir = app_test::scratch("workspace_diode");
    const fs::path deck = dir / "d.deck";
    // Before the deck has run: its device, no results.
    app_test::write(dir, "d.deck", app_test::diode_sweep());
    auto before = app::workspace_from_deck(deck);
    REQUIRE(before);
    CHECK(before->device->mesh().node_count() == 201);
    CHECK_FALSE(before->has_results());
    CHECK(before->record() == nullptr);

    auto taken = run(deck, app_test::diode_sweep());
    const results::Sweep sweep = *taken.report.sweep;
    auto from_run = app::workspace_from_run(std::move(taken.plan), std::move(taken.report));
    REQUIRE(from_run);
    REQUIRE(from_run->points.size() == 5);
    CHECK(from_run->geometry->dimension() == 1);
    REQUIRE(from_run->record());
    CHECK(from_run->record()->input_identity == sweep.run.input_identity);

    // The deck's outputs, now written, give the same results.
    auto from_deck = app::workspace_from_deck(deck);
    REQUIRE(from_deck);
    REQUIRE(from_deck->points.size() == 5);
    CHECK(from_deck->points.back().fields.potential_V == sweep.points.back().fields.potential_V);
    CHECK(from_deck->source == deck);

    // A result file with its device file.
    auto from_file = app::workspace_from_result(dir / "diode.nitsweep", dir / "diode.nitdev");
    REQUIRE(from_file);
    CHECK(from_file->points.size() == 5);
    CHECK(from_file->record()->mesh_digest == sweep.run.mesh_digest);
    // Not a result file; a missing device.
    CHECK_FALSE(app::workspace_from_result(dir / "diode.csv", dir / "diode.nitdev"));
    CHECK_FALSE(app::workspace_from_result(dir / "diode.nitsweep", dir / "none.nitdev"));
}

TEST_CASE("app workspace: a result is refused with another device's mesh", "[app]") {
    const fs::path dir = app_test::scratch("workspace_pairing");
    (void)run(dir / "d.deck", app_test::diode_sweep());
    // The same diode on a finer mesh: another device.
    std::string finer = app_test::diode_sweep();
    finer.replace(finer.find("201"), 3, "301");
    finer.replace(finer.find("diode.nitsweep"), 14, "other.nitsweep");
    finer.replace(finer.find("diode.nitdev"), 12, "other.nitdev");
    (void)run(dir / "e.deck", finer);
    const auto mismatched = app::workspace_from_result(dir / "diode.nitsweep", dir / "other.nitdev");
    REQUIRE_FALSE(mismatched);
    CHECK(mismatched.error().message.find("mesh digest") != std::string::npos);
}

TEST_CASE("app workspace: the run record and the convergence history (VA5)", "[app]") {
    const fs::path dir = app_test::scratch("workspace_record");
    auto taken = run(dir / "d.deck", app_test::diode_sweep());
    const results::Sweep sweep = *taken.report.sweep;
    auto w = app::workspace_from_run(std::move(taken.plan), std::move(taken.report));
    REQUIRE(w);
    const auto lines = app::run_record(*w);
    REQUIRE(find(lines, "Ended"));
    CHECK(find(lines, "Ended")->value == "completed");
    REQUIRE(find(lines, "Points"));
    CHECK(find(lines, "Points")->value == "5");
    REQUIRE(find(lines, "Input identity"));
    CHECK(find(lines, "Input identity")->value.size() == 16);
    REQUIRE(find(lines, "Revision"));
    CHECK(find(lines, "Revision")->value.starts_with(sweep.run.provenance.revision));
    CHECK(find(lines, "Mesh")->value == "1D, 201 nodes, 200 edges");
    // Each setting on its own line.
    for (const auto& [name, value] : sweep.run.settings) CHECK(find(lines, "  " + name));
    const auto history = app::convergence_history(*w);
    REQUIRE(history.size() == 5);
    for (std::size_t k = 0; k < 5; ++k) {
        CHECK(history[k].point == k);
        CHECK(history[k].converged);
        CHECK(history[k].iterations.size() == sweep.points[k].convergence.iterations.size());
    }
}

TEST_CASE("app workspace: Unit 24's figures (VA5)", "[app]") {
    const fs::path dir = app_test::scratch("workspace_figures");
    {
        // A diode sweep: the current at the last bias.
        auto taken = run(dir / "d.deck", app_test::diode_sweep());
        const double last = taken.report.sweep->points.back().terminal_current[0];
        auto w = app::workspace_from_run(std::move(taken.plan), std::move(taken.report));
        REQUIRE(w);
        const auto f = app::extracted_figures(*w);
        REQUIRE(f.size() == 1);
        CHECK(f[0].name == "Current at the last bias");
        CHECK(f[0].curve == "I(anode) vs V(anode)");
        CHECK(f[0].value == last);
        CHECK(f[0].unit == "A/cm\xC2\xB2");
    }
    {
        // A MOSFET's gate sweep: the transfer curve's figures, the drain sensed.
        auto taken = run(dir / "m.deck", app_test::mosfet2d);
        auto w = app::workspace_from_run(std::move(taken.plan), std::move(taken.report));
        REQUIRE(w);
        const auto f = app::extracted_figures(*w);
        const app::Figure* vth = figure(f, "Threshold (max gm)");
        REQUIRE(vth);
        CHECK(vth->curve == "I(drain) vs V(gate)");
        CHECK(vth->unit == "V");
        CHECK(vth->value > 0.0);
        CHECK(vth->value < 1.0);
        const app::Figure* ss = figure(f, "Subthreshold swing");
        REQUIRE(ss);
        CHECK(ss->value > 59.0);  // above the 300 K limit of 59.5 mV/dec, to rounding
        CHECK(figure(f, "Max transconductance"));
        const app::Figure* ratio = figure(f, "On/off ratio");
        REQUIRE(ratio);
        CHECK(ratio->value > 1.0);
    }
}
