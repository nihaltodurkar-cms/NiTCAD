// Checking the deck while it is edited and the structure preview (Unit 26b, XD1, XD4, XG1, XG3).
#include <catch2/catch_test_macros.hpp>

#include <atomic>
#include <condition_variable>
#include <filesystem>
#include <mutex>
#include <string>
#include <thread>

#include "NiTCAD/app/live_check.hpp"
#include "NiTCAD/app/run_controller.hpp"
#include "NiTCAD/app/view_settings.hpp"
#include "NiTCAD/app/workspace.hpp"
#include "NiTCAD/render/image.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
namespace fs = std::filesystem;

TEST_CASE("app live check: errors at their lines, then the plan, off the UI thread (XG1)", "[app]") {
    const auto ui = std::this_thread::get_id();
    std::atomic<int> wakes = 0;
    std::atomic<bool> on_ui = false;
    app::LiveChecker live([&] { ++wakes; }, [&](const std::string& text, const fs::path& path) {
        if (std::this_thread::get_id() == ui) on_ui = true;
        return app::check_deck(text, path);
    });
    CHECK_FALSE(live.pending());
    CHECK_FALSE(live.process());

    std::string text = app_test::diode_sweep();
    const auto at = text.find("kind = ohmic");
    text.replace(at, 12, "kind = ohmc");
    CHECK(live.submit(text, {}) == 1);
    CHECK(live.pending());
    live.wait();
    CHECK(live.process());
    CHECK_FALSE(live.pending());
    CHECK(live.applied() == 1);
    REQUIRE(live.result().diagnostics.size() == 1);
    CHECK(live.result().diagnostics[0].line == 14);
    CHECK_FALSE(live.result().plan);

    CHECK(live.submit(app_test::diode_sweep(), {}) == 2);
    live.wait();
    CHECK(live.process());
    CHECK(live.result().diagnostics.empty());
    auto plan = live.take_plan();
    REQUIRE(plan);
    CHECK(plan->device.mesh().node_count() == 201);
    CHECK_FALSE(live.take_plan());
    CHECK(wakes == 2);
    CHECK_FALSE(on_ui);
    CHECK(live.dropped() == 0);
}

TEST_CASE("app live check: a stale check's result is never applied (XG1)", "[app]") {
    // Each check waits inside the worker until the test lets that call go. The test knows a result
    // is posted from the wake, which the worker calls after posting it: the check returning is
    // earlier, and process() run in between finds nothing (it failed so on CI, under load).
    std::mutex m;
    std::condition_variable cv;
    int started = 0, released = 0, posted = 0;
    app::LiveChecker live(
        [&] {
            {
                std::scoped_lock lock(m);
                ++posted;
            }
            cv.notify_all();
        },
        [&](const std::string& text, const fs::path& path) {
            int call = 0;
            {
                std::unique_lock lock(m);
                call = ++started;
                cv.notify_all();
                cv.wait(lock, [&] { return released >= call; });
            }
            return app::check_deck(text, path);
        });
    const auto until = [&](auto done) {
        std::unique_lock lock(m);
        cv.wait(lock, done);
    };
    const auto let_go = [&](int call) {
        {
            std::scoped_lock lock(m);
            released = call;
        }
        cv.notify_all();
    };
    std::string broken = app_test::diode_sweep();
    broken.replace(broken.find("material = silicon"), 18, "material = silcon");
    live.submit(broken, {});
    until([&] { return started == 1; });
    // Edited while the first check runs: the newer text is what the user sees.
    live.submit(app_test::diode_sweep(), {});
    let_go(1);
    until([&] { return posted == 1; });
    // The first check's result has arrived and is dropped: the error it found is not shown.
    CHECK_FALSE(live.process());
    CHECK(live.applied() == 0);
    CHECK(live.result().diagnostics.empty());
    CHECK(live.dropped() == 1);
    CHECK(live.pending());
    until([&] { return started == 2; });
    let_go(2);
    live.wait();
    CHECK(live.process());
    CHECK(live.applied() == 2);
    CHECK(live.result().diagnostics.empty());
    CHECK(live.result().plan);
    // A text replaced before its check starts is never checked: the worker is held in check 3 while
    // two more texts arrive; only the last of them is checked next.
    live.submit(broken, {});
    until([&] { return started == 3; });
    live.submit(broken, {});
    live.submit(app_test::diode_sweep(), {});
    let_go(3);
    until([&] { return started == 4; });
    let_go(4);
    live.wait();
    CHECK(live.process());
    CHECK(live.applied() == 5);
    CHECK(live.result().diagnostics.empty());
    CHECK(started == 4);
    CHECK(live.dropped() == 2);
}

TEST_CASE("app live check: the deck is checked while a run goes on (XG1)", "[app]") {
    const fs::path dir = app_test::scratch("live_during_run");
    const std::string long_run = app_test::diode_sweep("0 : 0.8 : 0.0005");
    const fs::path deck = app_test::write(dir, "d.deck", long_run);
    app::RunController rc;
    REQUIRE(rc.start_run(long_run, deck));
    app::LiveChecker live;
    std::string edited = long_run;
    edited.insert(edited.find("[sweep]"), "[models]\nsrh = maybe\n\n");
    live.submit(edited, deck);
    live.wait();
    REQUIRE(live.process());
    rc.process();
    // The run is still going when the edited deck's error is known.
    CHECK(rc.state() == app::RunState::running);
    REQUIRE(live.result().diagnostics.size() == 1);
    CHECK(live.result().diagnostics[0].line == 22);
    rc.request_stop();
    rc.wait();
    rc.process();
    CHECK(rc.state() == app::RunState::done);
}

namespace {

render::Renderer& warp() {
    static auto r = render::Renderer::create({.adapter = render::Adapter::warp});
    REQUIRE(r);
    return *r;
}

}  // namespace

TEST_CASE("app structure preview: the deck's device as nitcad plot draws it, nothing solved (XD4, XG3)", "[app]") {
    const fs::path out = app_test::scratch("structure_preview");
    struct Case {
        const char* name;
        std::string text;
    };
    const std::vector<Case> cases{{"mosfet2d", app_test::mosfet2d}, {"generate3d", app_test::generate3d}};
    for (const Case& c : cases) {
        CAPTURE(c.name);
        const fs::path deck = app_test::write(out, std::string(c.name) + ".deck", c.text);
        app::LiveChecker live;
        live.submit(c.text, deck);
        live.wait();
        REQUIRE(live.process());
        auto plan = live.take_plan();
        REQUIRE(plan);
        auto w = app::workspace_from_structure(std::move(*plan), deck);
        REQUIRE(w);
        CHECK_FALSE(w->has_results());
        CHECK_FALSE(fs::exists(out / "mosfet.nitsweep"));
        for (const app::ViewKind kind : {app::ViewKind::solid, app::ViewKind::mesh}) {
            app::ViewSettings s;
            s.kind = kind;
            auto scene = app::build_scene(*w, s);
            REQUIRE(scene);
            auto mine = warp().draw(std::get<render::FieldScene>(*scene), 400, 300);
            REQUIRE(mine);
            const fs::path png = out / (std::string(c.name) + "_" + std::string(app::to_string(kind)) + ".png");
            std::string args = "plot " + app_test::quoted(deck);
            for (const std::string& a : app::plot_arguments(s)) args += " \"" + a + "\"";
            args += " --size 400x300 --adapter warp --out " + app_test::quoted(png);
            const auto cli = app_test::nitcad(args);
            CAPTURE(args, cli.output);
            REQUIRE(cli.status == 0);
            auto theirs = render::load_png(png);
            REQUIRE(theirs);
            const auto d = render::compare(*mine, *theirs, 0);
            REQUIRE(d);
            CHECK(d->max_channel == 0);
        }
    }
}

TEST_CASE("app structure preview: results kept while the mesh is the same (XD4)", "[app]") {
    const auto plan_of = [](const std::string& text) {
        auto r = app::check_deck(text, {});
        REQUIRE(r.plan);
        return std::move(*r.plan);
    };
    auto w = app::workspace_from_structure(plan_of(app_test::mosfet2d), {});
    REQUIRE(w);
    // Another sweep, a model or a doping changes nothing of the mesh; its nodes do.
    std::string other = app_test::mosfet2d;
    other.replace(other.find("0 : 1.5 : 0.25"), 14, "0 : 1 : 0.5");
    other.replace(other.find("1e16"), 4, "2e16");
    CHECK(app::same_mesh(*w, plan_of(other).device.mesh()));
    std::string finer = app_test::mosfet2d;
    finer.replace(finer.find("linspace(0, 2um, 21)"), 20, "linspace(0, 2um, 41)");
    CHECK_FALSE(app::same_mesh(*w, plan_of(finer).device.mesh()));
    app::Workspace empty;
    CHECK_FALSE(app::same_mesh(empty, plan_of(finer).device.mesh()));
}
