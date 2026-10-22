// Run comparison in the application (Unit 26b, XD6, XG5): the figures' shifts (ΔV_th is the difference
// of the two thresholds), the run records side by side, and whether two runs share a mesh.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <filesystem>
#include <string>

#include "NiTCAD/app/compare.hpp"
#include "NiTCAD/app/run_controller.hpp"
#include "NiTCAD/app/workspace.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
namespace fs = std::filesystem;

namespace {

app::Workspace run(const fs::path& dir, const std::string& name, const std::string& text) {
    app::RunController rc;
    REQUIRE(rc.start_run(text, app_test::write(dir, name, text)));
    rc.wait();
    rc.process();
    REQUIRE(rc.state() == app::RunState::done);
    auto taken = std::move(*rc.take_results());
    auto w = app::workspace_from_run(std::move(taken.plan), std::move(taken.report));
    REQUIRE(w);
    return std::move(*w);
}

const app::Figure* figure(const std::vector<app::Figure>& f, const std::string& name) {
    const auto it = std::ranges::find_if(f, [&](const app::Figure& x) { return x.name == name; });
    return it == f.end() ? nullptr : &*it;
}

}  // namespace

TEST_CASE("app compare 26b: the figures' shifts and the records side by side (XD6, XG5)", "[app]") {
    const fs::path dir = app_test::scratch("compare_mosfet");
    const app::Workspace a = run(dir, "a.deck", app_test::mosfet2d);
    // The same transistor with a higher gate work function: its threshold moves up.
    std::string higher = app_test::mosfet2d;
    higher.replace(higher.find("work_function = 4.1"), 19, "work_function = 4.3");
    higher.replace(higher.find("mosfet.nitsweep"), 15, "higher.nitsweep");
    higher.replace(higher.find("mosfet.nitdev"), 13, "higher.nitdev");
    const app::Workspace b = run(dir, "b.deck", higher);
    REQUIRE(app::comparable(a, b));

    // ΔV_th is the difference of the thresholds each run's figures extract.
    const auto fa = app::extracted_figures(a), fb = app::extracted_figures(b);
    const app::Figure* va = figure(fa, "Threshold (max gm)");
    const app::Figure* vb = figure(fb, "Threshold (max gm)");
    REQUIRE(va);
    REQUIRE(vb);
    const auto shifts = app::figure_shifts(a, b);
    REQUIRE(shifts.size() == fa.size());
    const auto vth = std::ranges::find_if(shifts, [](const app::FigureShift& s) { return s.name == "Threshold (max gm)"; });
    REQUIRE(vth != shifts.end());
    CHECK(vth->curve == va->curve);
    CHECK(vth->unit == "V");
    CHECK(vth->a == va->value);
    CHECK(vth->b == vb->value);
    CHECK(vth->delta() == vb->value - va->value);
    CHECK(vth->delta() > 0.05);
    for (std::size_t k = 0; k < shifts.size(); ++k) CHECK(shifts[k].name == fa[k].name);  // in A's order

    // The records: A's lines in order, the input identity marked as differing, the mesh not.
    const auto ra = app::run_record(a);
    const auto rows = app::compare_records(a, b);
    REQUIRE(rows.size() >= ra.size());
    for (std::size_t k = 0; k < ra.size(); ++k) {
        CHECK(rows[k].name == ra[k].name);
        CHECK(rows[k].a == ra[k].value);
    }
    const auto row = [&](const std::string& name) {
        const auto it = std::ranges::find_if(rows, [&](const app::RecordRow& r) { return r.name == name; });
        REQUIRE(it != rows.end());
        return *it;
    };
    CHECK(row("Input identity").differs);
    CHECK_FALSE(row("Mesh").differs);
    CHECK_FALSE(row("Points").differs);

    // A run beside itself: no shift, nothing differs.
    for (const app::FigureShift& s : app::figure_shifts(a, a)) CHECK(s.delta() == 0.0);
    for (const app::RecordRow& r : app::compare_records(a, a)) CHECK_FALSE(r.differs);
}

TEST_CASE("app compare 26b: runs on two meshes are not compared node by node (XG5)", "[app]") {
    const fs::path dir = app_test::scratch("compare_meshes");
    const app::Workspace a = run(dir, "a.deck", app_test::diode_sweep());
    std::string finer = app_test::diode_sweep();
    finer.replace(finer.find("201"), 3, "301");
    finer.replace(finer.find("diode.nitsweep"), 14, "finer.nitsweep");
    finer.replace(finer.find("diode.nitdev"), 12, "finer.nitdev");
    const app::Workspace b = run(dir, "b.deck", finer);
    CHECK_FALSE(app::comparable(a, b));
    CHECK(app::comparable(a, a));
    // The figures still compare: one current at the last bias each.
    const auto shifts = app::figure_shifts(a, b);
    REQUIRE(shifts.size() == 1);
    CHECK(shifts[0].name == "Current at the last bias");
    // The mesh line differs; lines only B has would follow A's.
    const auto rows = app::compare_records(a, b);
    const auto mesh = std::ranges::find_if(rows, [](const app::RecordRow& r) { return r.name == "Mesh"; });
    REQUIRE(mesh != rows.end());
    CHECK(mesh->differs);
    CHECK(mesh->b == "1D, 301 nodes, 300 edges");
    // A figure of the same name on another curve is not a shift: the cathode swept instead.
    std::string cathode = app_test::diode_sweep();
    cathode.replace(cathode.find("contact anode = "), 16, "contact cathode = ");
    cathode.replace(cathode.find("diode.nitsweep"), 14, "cathode.nitsweep");
    cathode.replace(cathode.find("diode.nitdev"), 12, "cathode.nitdev");
    const app::Workspace c = run(dir, "c.deck", cathode);
    REQUIRE(app::extracted_figures(c).size() == 1);
    CHECK(app::extracted_figures(c)[0].name == "Current at the last bias");
    CHECK(app::extracted_figures(c)[0].curve != app::extracted_figures(a)[0].curve);
    CHECK(app::figure_shifts(a, c).empty());
    // The lines only B has follow A's: an unsolved deck (its mesh line only) beside a run.
    const fs::path unrun_dir = app_test::scratch("compare_unrun");
    auto unrun = app::workspace_from_deck(app_test::write(unrun_dir, "d.deck", app_test::diode_sweep()));
    REQUIRE(unrun);
    REQUIRE_FALSE(unrun->has_results());
    const auto only_b = app::compare_records(*unrun, a);
    const auto lines_b = app::run_record(a);
    REQUIRE(only_b.size() == lines_b.size());
    CHECK(only_b[0].name == "Mesh");
    CHECK_FALSE(only_b[0].differs);
    for (std::size_t k = 1; k < only_b.size(); ++k) {
        CHECK(only_b[k].name == lines_b[k].name);
        CHECK(only_b[k].a.empty());
        CHECK(only_b[k].b == lines_b[k].value);
        CHECK(only_b[k].differs);
    }
    // A workspace without a device compares with nothing.
    const app::Workspace empty{};
    CHECK_FALSE(app::comparable(a, empty));
    CHECK_FALSE(app::comparable(empty, a));
}
