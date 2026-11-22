// The run controller (Unit 26, WG1, WG2, WG3): the state machine, the UI thread never blocked by a
// solve, Stop keeping the completed part, and a run identical to `nitcad run`.
#include <catch2/catch_test_macros.hpp>

#include <windows.h>

#include <algorithm>
#include <bit>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include "NiTCAD/app/run_controller.hpp"
#include "NiTCAD/io/files.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
namespace fs = std::filesystem;
using app::RunState;

namespace {

bool same_bits(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (std::bit_cast<std::uint64_t>(a[i]) != std::bit_cast<std::uint64_t>(b[i])) return false;
    }
    return true;
}

// Every field, current and record entry of two sweeps, to the bit; provenance aside (times differ).
void check_identical(const results::Sweep& a, const results::Sweep& b) {
    CHECK(a.run.input_identity == b.run.input_identity);
    CHECK(a.run.mesh_digest == b.run.mesh_digest);
    CHECK(a.run.settings == b.run.settings);
    REQUIRE(a.points.size() == b.points.size());
    for (std::size_t k = 0; k < a.points.size(); ++k) {
        CAPTURE(k);
        const auto& p = a.points[k];
        const auto& q = b.points[k];
        CHECK(same_bits(p.bias_V, q.bias_V));
        CHECK(same_bits(p.fields.potential_V, q.fields.potential_V));
        CHECK(same_bits(p.fields.n_cm3, q.fields.n_cm3));
        CHECK(same_bits(p.fields.p_cm3, q.fields.p_cm3));
        CHECK(same_bits(p.terminal_current, q.terminal_current));
        CHECK(same_bits(p.edge_current_n, q.edge_current_n));
        CHECK(same_bits(p.edge_current_p, q.edge_current_p));
        CHECK(p.convergence.iterations.size() == q.convergence.iterations.size());
    }
}

// A hidden message-only window that the worker wakes by posting it a message, as the application's
// window does, and a message loop with a 10 ms timer: the largest gap between timer messages while
// the controller is busy is how long the UI thread could not run.
struct Pump {
    static constexpr UINT wake_message = WM_APP + 7;
    HWND hwnd = CreateWindowExW(0, L"STATIC", L"", 0, 0, 0, 0, 0, HWND_MESSAGE, nullptr, nullptr, nullptr);
    ~Pump() { DestroyWindow(hwnd); }
    app::RunController::Wake wake() const {
        return [h = hwnd] { PostMessageW(h, wake_message, 0, 0); };
    }
    struct Measured {
        double longest_gap_s = 0.0;
        int timer_messages = 0;
    };
    // Pumps until the controller is no longer busy; `on_event` runs after each process().
    template <class F>
    Measured run(app::RunController& rc, F&& on_event) {
        Measured m;
        SetTimer(hwnd, 1, 10, nullptr);
        auto last = std::chrono::steady_clock::now();
        MSG msg;
        while (rc.busy() && GetMessageW(&msg, nullptr, 0, 0) > 0) {
            if (msg.message == WM_TIMER) {
                const auto now = std::chrono::steady_clock::now();
                m.longest_gap_s = std::max(m.longest_gap_s, std::chrono::duration<double>(now - last).count());
                last = now;
                ++m.timer_messages;
            } else if (msg.message == wake_message) {
                if (rc.process()) on_event();
            } else {
                DispatchMessageW(&msg);
            }
        }
        KillTimer(hwnd, 1);
        return m;
    }
};

// The diode swept in small steps: a run of many points that takes a while.
std::string long_diode() { return app_test::diode_sweep("0 : 0.8 : 0.0005"); }

}  // namespace

TEST_CASE("app run controller: checks (WG3)", "[app]") {
    const fs::path dir = app_test::scratch("controller_check");
    app::RunController rc;
    CHECK(rc.state() == RunState::idle);
    CHECK_FALSE(rc.busy());
    rc.request_stop();  // nothing to stop
    CHECK(rc.state() == RunState::idle);

    REQUIRE(rc.start_check(app_test::diode_sweep(), dir / "d.deck"));
    CHECK(rc.state() == RunState::checking);
    CHECK(rc.busy());
    rc.wait();
    CHECK(rc.process());
    CHECK(rc.state() == RunState::idle);
    CHECK(rc.diagnostics().empty());
    REQUIRE(rc.plan());
    CHECK(rc.plan()->points.size() == 5);
    CHECK_FALSE(rc.report());
    CHECK(rc.points().empty());

    // A deck that does not build: failed, the error at its line (4: the comment, [mesh], type, then
    // the unknown key).
    std::string bad = app_test::diode_sweep();
    bad.insert(bad.find("x = linspace"), "bogus = 1\n");
    REQUIRE(rc.start_check(bad, dir / "d.deck"));
    rc.wait();
    rc.process();
    CHECK(rc.state() == RunState::failed);
    REQUIRE(rc.diagnostics().size() == 1);
    CHECK(rc.diagnostics()[0].line == 4);
    REQUIRE(rc.error());
    CHECK(rc.error()->message == rc.diagnostics()[0].message);
    CHECK_FALSE(rc.plan());

    // A run of a deck that does not build fails the same way, without solving.
    REQUIRE(rc.start_run(bad, dir / "d.deck"));
    rc.wait();
    rc.process();
    CHECK(rc.state() == RunState::failed);
    CHECK(rc.diagnostics().size() == 1);
    CHECK(rc.points().empty());
}

TEST_CASE("app run controller: a run, its progress rows and results (WG3)", "[app]") {
    const fs::path dir = app_test::scratch("controller_run");
    app::RunController rc;
    REQUIRE(rc.start_run(app_test::diode_sweep(), dir / "d.deck"));
    CHECK(rc.state() == RunState::running);
    // Busy: no second start.
    CHECK_FALSE(rc.start_check(app_test::diode_sweep(), dir / "d.deck"));
    CHECK_FALSE(rc.start_run(app_test::diode_sweep(), dir / "d.deck"));
    rc.wait();
    rc.process();
    CHECK(rc.state() == RunState::done);
    CHECK_FALSE(rc.error());
    CHECK_FALSE(rc.stopped_early());
    CHECK(rc.planned_points() == 5);
    REQUIRE(rc.report());
    REQUIRE(rc.report()->sweep);
    const auto& sweep = *rc.report()->sweep;
    CHECK(sweep.points.size() == 5);
    // The rows: the equilibrium, then one per bias point, each converged with the sweep's iterations.
    const auto& rows = rc.points();
    REQUIRE(rows.size() == 6);
    CHECK(rows[0].phase == solve::Phase::equilibrium);
    for (std::size_t k = 0; k < 5; ++k) {
        CAPTURE(k);
        CHECK(rows[k + 1].phase == solve::Phase::bias);
        CHECK(rows[k + 1].point == k);
        CHECK(rows[k + 1].point_count == 5);
        CHECK(rows[k + 1].converged);
        CHECK(rows[k + 1].iterations == static_cast<int>(sweep.points[k].convergence.iterations.size()));
        CHECK(rows[k + 1].update == sweep.points[k].convergence.iterations.back().update);
        CHECK(rows[k + 1].seconds >= 0.0);
    }
    // The convergence plot's iterations: numbered from 1, each pointing at its row.
    const auto& its = rc.iterations();
    REQUIRE_FALSE(its.empty());
    for (std::size_t i = 0; i < its.size(); ++i) {
        CHECK(its[i].number == i + 1);
        CHECK(its[i].row < rows.size());
    }
    // The outputs are written as `nitcad run` writes them.
    CHECK(fs::exists(dir / "diode.nitsweep"));
    CHECK(fs::exists(dir / "diode.nitdev"));
    CHECK(fs::exists(dir / "diode.csv"));
    CHECK(rc.elapsed_seconds() > 0.0);
    // Taking the results out leaves the controller done, without them.
    auto taken = rc.take_results();
    REQUIRE(taken);
    CHECK(taken->report.sweep->points.size() == 5);
    CHECK(taken->plan.points.size() == 5);
    CHECK(rc.state() == RunState::done);
    CHECK_FALSE(rc.report());
    CHECK_FALSE(rc.take_results());
}

TEST_CASE("app run controller: the UI thread keeps running during a long run (WG1)", "[app]") {
    const fs::path dir = app_test::scratch("controller_pump");
    Pump pump;
    REQUIRE(pump.hwnd != nullptr);
    app::RunController rc(pump.wake());
    int events = 0;
    REQUIRE(rc.start_run(long_diode(), dir / "d.deck"));
    const auto m = pump.run(rc, [&] { ++events; });
    CHECK(rc.state() == RunState::done);
    REQUIRE(rc.report());
    CHECK(rc.report()->sweep->points.size() == 1601);
    // The run took many timer periods, and no gap between timer messages reached 100 ms.
    CAPTURE(m.timer_messages, m.longest_gap_s, rc.elapsed_seconds());
    CHECK(m.timer_messages >= 20);
    CHECK(1e3 * m.longest_gap_s < app_test::wg1_limit_ms());
    // Progress arrived while it ran, in batches (one wake per non-empty queue).
    CHECK(events > 1);
    CHECK(rc.points().size() == 1602);
}

TEST_CASE("app run controller: Stop keeps the completed part (WG1)", "[app]") {
    const fs::path dir = app_test::scratch("controller_stop");
    Pump pump;
    app::RunController rc(pump.wake());
    REQUIRE(rc.start_run(long_diode(), dir / "d.deck"));
    bool asked = false;
    pump.run(rc, [&] {
        // Stop once three bias points have converged.
        const auto converged = std::ranges::count_if(rc.points(), [](const app::PointRow& r) {
            return r.phase == solve::Phase::bias && r.converged;
        });
        if (!asked && converged >= 3) {
            rc.request_stop();
            CHECK(rc.state() == RunState::stopping);
            asked = true;
        }
    });
    REQUIRE(asked);
    CHECK(rc.state() == RunState::done);
    const auto stopped = rc.stopped_early();
    REQUIRE(stopped);
    CHECK(stopped->code == base::ErrorCode::cancelled);
    REQUIRE(rc.report());
    const auto& sweep = *rc.report()->sweep;
    CHECK(sweep.points.size() >= 3);
    CHECK(sweep.points.size() < 1601);
    // The completed part is written.
    auto saved = io::load_sweep(dir / "diode.nitsweep");
    REQUIRE(saved);
    CHECK(saved->points.size() == sweep.points.size());
}

TEST_CASE("app run controller: a run is bit-identical to nitcad run (WG2)", "[app]") {
    // The same deck in two directories: one run by the controller's worker thread, one by the nitcad
    // executable.
    const fs::path a = app_test::scratch("controller_identical_app");
    const fs::path b = app_test::scratch("controller_identical_cli");
    for (const std::string& text : {app_test::diode_sweep(), app_test::mosfet2d}) {
        const char* name = text == app_test::mosfet2d ? "mosfet.nitsweep" : "diode.nitsweep";
        CAPTURE(name);
        app::RunController rc;
        REQUIRE(rc.start_run(text, app_test::write(a, "d.deck", text)));
        rc.wait();
        rc.process();
        REQUIRE(rc.state() == RunState::done);
        const results::Sweep& mine = *rc.report()->sweep;
        const auto cli = app_test::nitcad("run " + app_test::quoted(app_test::write(b, "d.deck", text)));
        CAPTURE(cli.output);
        REQUIRE(cli.status == 0);
        auto theirs = io::load_sweep(b / name);
        REQUIRE(theirs);
        check_identical(mine, *theirs);
        // And the file the app's run wrote is the same sweep.
        auto written = io::load_sweep(a / name);
        REQUIRE(written);
        check_identical(*written, *theirs);
    }
}

namespace {

std::string pn_example() {
    std::ifstream f(fs::path(NITCAD_EXAMPLES_DIR) / "devices" / "pn_diode_1d.deck", std::ios::binary);
    return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

}  // namespace

TEST_CASE("app run controller: a dimension sweep runs each value into its folder, as nitcad sweep (DT6)", "[app][sweep_run]") {
    const fs::path dir = app_test::scratch("controller_sweep");
    const std::string text = pn_example();
    const fs::path deck = app_test::write(dir, "pn.deck", text);
    Pump pump;
    app::RunController rc(pump.wake());
    std::vector<std::size_t> seen;
    REQUIRE(rc.start_sweep(text, deck, "Lp", {"1um", "1.5um"}));
    REQUIRE(rc.sweep());
    CHECK(rc.sweep()->name == "Lp");
    pump.run(rc, [&] {
        if (rc.sweep() && (seen.empty() || seen.back() != rc.sweep()->current)) seen.push_back(rc.sweep()->current);
        // The rows are the current value's run's only.
        if (!rc.points().empty()) CHECK(rc.points().front().point <= rc.points().back().point);
    });
    CHECK(rc.state() == RunState::done);
    CHECK(seen == std::vector<std::size_t>{0, 1});
    CHECK(rc.sweep()->completed == 2);
    CHECK_FALSE(rc.sweep()->stopped);
    // The deck's own outputs are not written, and its file is not edited.
    CHECK_FALSE(fs::exists(dir / "pn_diode_1d.nitsweep"));

    // nitcad sweep on a copy writes the same runs.
    const fs::path cli_dir = app_test::scratch("controller_sweep_cli");
    const fs::path cli_deck = app_test::write(cli_dir, "pn.deck", text);
    const auto cli = app_test::nitcad("sweep " + app_test::quoted(cli_deck) + " Lp=1um,1.5um");
    INFO(cli.output);
    REQUIRE(cli.status == 0);
    for (const std::string v : {"Lp=1um", "Lp=1.5um"}) {
        auto a = io::load_sweep(dir / "pn_sweep_Lp" / v / "pn_diode_1d.nitsweep");
        auto b = io::load_sweep(cli_dir / "pn_sweep_Lp" / v / "pn_diode_1d.nitsweep");
        REQUIRE(a);
        REQUIRE(b);
        check_identical(*a, *b);
    }
}

TEST_CASE("app run controller: Stop ends a dimension sweep in its current value (DT6)", "[app][sweep_run]") {
    const fs::path dir = app_test::scratch("controller_sweep_stop");
    std::string text = pn_example();
    text.replace(text.find("0 : 0.7 : 0.05"), 14, "0 : 0.7 : 0.0005");
    const fs::path deck = app_test::write(dir, "pn.deck", text);
    Pump pump;
    app::RunController rc(pump.wake());
    REQUIRE(rc.start_sweep(text, deck, "Lp", {"1um", "1.5um"}));
    bool asked = false;
    pump.run(rc, [&] {
        const auto converged = std::ranges::count_if(rc.points(), [](const app::PointRow& r) {
            return r.phase == solve::Phase::bias && r.converged;
        });
        if (!asked && converged >= 3) {
            rc.request_stop();
            asked = true;
        }
    });
    REQUIRE(asked);
    CHECK(rc.state() == RunState::done);
    REQUIRE(rc.sweep());
    CHECK(rc.sweep()->stopped);
    CHECK(rc.sweep()->completed == 0);
    CHECK(rc.sweep()->current == 0);
    // The stopped value's completed part is written; the next value never ran.
    auto saved = io::load_sweep(dir / "pn_sweep_Lp" / "Lp=1um" / "pn_diode_1d.nitsweep");
    REQUIRE(saved);
    CHECK(saved->points.size() >= 3);
    CHECK_FALSE(fs::exists(dir / "pn_sweep_Lp" / "Lp=1.5um"));
}

TEST_CASE("app run controller: a sweep of an unknown parameter or a wrong value fails before any run (DT6)", "[app][sweep_run]") {
    const fs::path dir = app_test::scratch("controller_sweep_refused");
    const std::string text = pn_example();
    const fs::path deck = app_test::write(dir, "pn.deck", text);
    for (const auto& [name, value] : {std::pair<std::string, std::string>{"Wfin", "1nm"}, {"Lp", "5"}, {"Na", "1nm"}}) {
        app::RunController rc;
        REQUIRE(rc.start_sweep(text, deck, name, {value}));
        rc.wait();
        rc.process();
        CHECK(rc.state() == RunState::failed);
        REQUIRE(rc.error());
        CHECK(rc.error()->message.find(name) != std::string::npos);
        CHECK_FALSE(fs::exists(dir / ("pn_sweep_" + name)));
    }
    app::RunController rc;
    CHECK_FALSE(rc.start_sweep(text, deck, "Lp", {}));
}
