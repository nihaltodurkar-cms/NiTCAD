// The window's shell decisions (Unit 26d, DECISIONS.md ZD4-ZD6, ZG4-ZG6): what a drop opens, the taskbar's
// progress through each run outcome, the title bar's theme.
#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <vector>

#include "NiTCAD/app/integration.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
namespace fs = std::filesystem;
using app::RunState;
using app::TaskbarFlag;
using app::TaskbarProgress;

TEST_CASE("app shell 26d: what a path and a drop open (ZD4, ZG4)", "[app]") {
    const fs::path dir = app_test::scratch("shell_drop");
    // A deck opens as a deck; a result with its device file and, when there, its mesh file.
    const auto deck = app::open_target(dir / "a.deck");
    CHECK(deck.deck);
    CHECK(deck.path == dir / "a.deck");
    const auto plain = app::open_target(dir / "r.nitsweep");
    CHECK_FALSE(plain.deck);
    CHECK(plain.device == dir / "r.nitdev");
    CHECK(plain.mesh.empty());
    app_test::write(dir, "g.nitmesh", "");
    const auto generated = app::open_target(dir / "g.nittran");
    CHECK(generated.device == dir / "g.nitdev");
    CHECK(generated.mesh == dir / "g.nitmesh");

    // One file is taken; two, none or a folder are refused with the reason.
    const std::vector<fs::path> one{dir / "a.deck"};
    REQUIRE(app::dropped_file(one));
    CHECK(*app::dropped_file(one) == dir / "a.deck");
    const std::vector<fs::path> two{dir / "a.deck", dir / "b.deck"};
    CHECK(app::dropped_file(two).error() == "Drop one deck or result file at a time.");
    CHECK(app::dropped_file(std::vector<fs::path>{}).error() == "Drop one deck or result file at a time.");
    const std::vector<fs::path> folder{dir};
    CHECK(app::dropped_file(folder).error() == "A folder does not open: drop a deck or a result file.");
}

TEST_CASE("app shell 26d: the taskbar's progress through each run outcome (ZD5, ZQ3, ZG5)", "[app]") {
    const TaskbarProgress idle;
    const auto at = [](TaskbarProgress p, RunState s, bool run, double f = -1.0) { return app::taskbar_progress(s, run, f, p); };
    // A run: its fraction while it goes on (clamped), indeterminate when its length is not known.
    const TaskbarProgress running = at(idle, RunState::running, true, 0.25);
    CHECK(running == TaskbarProgress{TaskbarFlag::normal, 0.25});
    CHECK(at(idle, RunState::running, true, 1.5) == TaskbarProgress{TaskbarFlag::normal, 1.0});
    CHECK(at(idle, RunState::running, true, -1.0) == TaskbarProgress{TaskbarFlag::indeterminate, 0.0});
    CHECK(at(running, RunState::stopping, true, 0.5) == TaskbarProgress{TaskbarFlag::normal, 0.5});
    // Done (or stopped early, or by Stop: the controller ends done): cleared.
    CHECK(at(running, RunState::done, true) == idle);
    // Failed: the error state at the fraction reached (full when unknown), kept while idle...
    const TaskbarProgress failed = at(running, RunState::failed, true);
    CHECK(failed == TaskbarProgress{TaskbarFlag::error, 0.25});
    CHECK(at(at(idle, RunState::running, true, -1.0), RunState::failed, true) == TaskbarProgress{TaskbarFlag::error, 1.0});
    CHECK(at(failed, RunState::idle, false) == failed);
    CHECK(at(failed, RunState::idle, true) == failed);
    // ...until the next check (or an open's) or the next run.
    CHECK(at(failed, RunState::checking, false) == idle);
    CHECK(at(failed, RunState::running, true, 0.0) == TaskbarProgress{TaskbarFlag::normal, 0.0});
    // A check's own outcome is not shown: it neither sets the error nor shows progress.
    CHECK(at(idle, RunState::checking, false) == idle);
    CHECK(at(idle, RunState::failed, false) == idle);
    CHECK(at(idle, RunState::done, false) == idle);
}

TEST_CASE("app shell 26d: the title bar follows the app theme (ZD6, ZQ4, ZG6)", "[app]") {
    CHECK(app::dark_title_bar(0ul));
    CHECK_FALSE(app::dark_title_bar(1ul));
    CHECK_FALSE(app::dark_title_bar(2ul));
    CHECK_FALSE(app::dark_title_bar(std::nullopt));  // no value: Windows' default, light
}
