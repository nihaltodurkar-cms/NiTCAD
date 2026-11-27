// Unit 26p (GU3, GUC3): the guide pane, by window messages and its own actions, over a plain parent window.
#include <catch2/catch_test_macros.hpp>

#include "guide_pane.hpp"

#include <windows.h>

#include <string>
#include <vector>

using namespace NiTCAD::app;

namespace {

struct Recorder final : GuidePane::Host {
    std::vector<std::size_t> opened;
    std::vector<LessonStep> acted;
    void guide_open(std::size_t lesson) override { opened.push_back(lesson); }
    void guide_act(const LessonStep& step) override { acted.push_back(step); }
};

struct Parent {
    HWND hwnd = CreateWindowExW(0, L"STATIC", L"guide test", WS_OVERLAPPEDWINDOW | WS_CLIPCHILDREN, 0, 0, 800, 600, nullptr,
                                nullptr, GetModuleHandleW(nullptr), nullptr);
    ~Parent() { DestroyWindow(hwnd); }
};

LessonStep step(std::string heading, LessonAction::Kind kind = LessonAction::Kind::none, std::string name = {},
                std::string value = {}) {
    LessonStep s;
    s.heading = std::move(heading);
    s.prose = "First paragraph.\n\nSecond paragraph.";
    s.action.kind = kind;
    s.action.name = std::move(name);
    s.action.value = std::move(value);
    return s;
}

std::vector<Lesson> two_lessons() {
    Lesson a;
    a.title = "The gate and the channel";
    a.order = 1;
    a.steps = {step("What the gate does"), step("Sweep the doping", LessonAction::Kind::sweep, "Na", "1e16, 1e17, 1e18"),
               step("Run it", LessonAction::Kind::run)};
    Lesson b;
    b.title = "Velocity saturation";
    b.track = LessonTrack::advanced;
    b.order = 12;
    b.steps = {step("Turn it off", LessonAction::Kind::edit, "field_mobility", "off"),
               step("Thinner oxide", LessonAction::Kind::set, "tox", "2nm")};
    b.steps[0].action.section = "models";
    return {a, b};
}

void click(const GuidePane& g, HWND button) {
    SendMessageW(g.hwnd(), WM_COMMAND, MAKEWPARAM(GetDlgCtrlID(button), BN_CLICKED), reinterpret_cast<LPARAM>(button));
}

}  // namespace

TEST_CASE("guide pane: shows the first step", "[guide_window]") {
    Parent p;
    Recorder r;
    GuidePane g(p.hwnd, r, 7);
    g.set_lessons(two_lessons());
    g.show_lesson(0);
    CHECK(g.lesson() == 0);
    CHECK(g.step() == 0);
    CHECK(g.heading_text() == "What the gate does");
    CHECK(g.prose_text() == "First paragraph.\r\n\r\nSecond paragraph.");
    CHECK(g.position_text() == "Step 1 of 3");
    CHECK(g.action_label().empty());
    CHECK_FALSE(IsWindowEnabled(g.action_button()));
}

TEST_CASE("guide pane: Back and Next clamp at both ends", "[guide_window]") {
    Parent p;
    Recorder r;
    GuidePane g(p.hwnd, r, 7);
    g.set_lessons(two_lessons());
    g.show_lesson(0);
    g.go(-1);
    CHECK(g.step() == 0);
    g.go(1);
    g.go(1);
    g.go(1);
    CHECK(g.step() == 2);
    CHECK(g.position_text() == "Step 3 of 3");
    g.show_lesson(1);
    CHECK(g.step() == 0);
    CHECK(g.heading_text() == "Turn it off");
}

TEST_CASE("guide pane: the action's label per kind", "[guide_window]") {
    Parent p;
    Recorder r;
    GuidePane g(p.hwnd, r, 7);
    g.set_lessons(two_lessons());
    g.show_lesson(0);
    g.go(1);
    CHECK(g.action_label() == "Sweep Na = 1e16, 1e17, 1e18");
    g.go(1);
    CHECK(g.action_label() == "Run");
    g.show_lesson(1);
    CHECK(g.action_label() == "Set field_mobility = off");
    g.go(1);
    CHECK(g.action_label() == "Set tox = 2nm");
}

TEST_CASE("guide pane: the action button calls the host with the step", "[guide_window]") {
    Parent p;
    Recorder r;
    GuidePane g(p.hwnd, r, 7);
    g.set_lessons(two_lessons());
    g.show_lesson(0);
    g.go(1);
    REQUIRE(IsWindowEnabled(g.action_button()));
    click(g, g.action_button());
    REQUIRE(r.acted.size() == 1);
    CHECK(r.acted[0].action.kind == LessonAction::Kind::sweep);
    CHECK(r.acted[0].action.name == "Na");
}

TEST_CASE("guide pane: an error keeps the step", "[guide_window]") {
    Parent p;
    Recorder r;
    GuidePane g(p.hwnd, r, 7);
    g.set_lessons(two_lessons());
    g.show_lesson(0);
    g.go(1);
    g.set_report("Threshold (max gm): 0.31 V");
    CHECK(g.report_text() == "Threshold (max gm): 0.31 V");
    g.set_error("the deck has no parameter 'Na'");
    CHECK(g.step() == 1);
    CHECK(g.report_text() == "the deck has no parameter 'Na'");
    g.go(1);  // a new step clears what the last one gave
    CHECK(g.report_text().empty());
}

TEST_CASE("guide pane: busy turns the action off, not Back and Next", "[guide_window]") {
    Parent p;
    Recorder r;
    GuidePane g(p.hwnd, r, 7);
    g.set_lessons(two_lessons());
    g.show_lesson(0);
    g.go(1);
    g.set_busy(true);
    CHECK_FALSE(IsWindowEnabled(g.action_button()));
    g.act();
    CHECK(r.acted.empty());
    g.go(1);
    CHECK(g.step() == 2);
    g.set_busy(false);
    CHECK(IsWindowEnabled(g.action_button()));
}

TEST_CASE("guide pane: picking a lesson asks the host", "[guide_window]") {
    Parent p;
    Recorder r;
    GuidePane g(p.hwnd, r, 7);
    g.set_lessons(two_lessons());
    HWND picker = FindWindowExW(g.hwnd(), nullptr, L"COMBOBOX", nullptr);
    REQUIRE(picker);
    CHECK(SendMessageW(picker, CB_GETCOUNT, 0, 0) == 2);
    SendMessageW(picker, CB_SETCURSEL, 1, 0);
    SendMessageW(g.hwnd(), WM_COMMAND, MAKEWPARAM(GetDlgCtrlID(picker), CBN_SELCHANGE), reinterpret_cast<LPARAM>(picker));
    CHECK(r.opened == std::vector<std::size_t>{1});
}

TEST_CASE("guide pane: relayout at the DPI", "[guide_window]") {
    Parent p;
    Recorder r;
    GuidePane g(p.hwnd, r, 7);
    g.set_lessons(two_lessons());
    g.show_lesson(0);
    for (const int dpi : {144, 192}) {
        const int width = 320 * dpi / 96, height = 700 * dpi / 96;
        MoveWindow(g.hwnd(), 0, 0, width, height, FALSE);
        g.set_dpi(dpi, nullptr, nullptr);
        RECT prose{};
        GetWindowRect(g.prose_edit(), &prose);
        MapWindowPoints(nullptr, g.hwnd(), reinterpret_cast<POINT*>(&prose), 2);
        const int margin = 8 * dpi / 96;
        CHECK(prose.left == margin);
        CHECK(prose.right == width - margin);
        for (HWND c = GetWindow(g.hwnd(), GW_CHILD); c; c = GetWindow(c, GW_HWNDNEXT)) {
            RECT rc{};
            GetWindowRect(c, &rc);
            MapWindowPoints(nullptr, g.hwnd(), reinterpret_cast<POINT*>(&rc), 2);
            CHECK(rc.left >= 0);
            CHECK(rc.right <= width);
            CHECK(rc.bottom <= height);
        }
    }
}

// --- the guide in the main window (Task 6: GU4, GU5) -------------------------------------------------

#include "main_window.hpp"
#include "fixtures.hpp"
#include "NiTCAD/app/commands.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iterator>

namespace {

namespace fs = std::filesystem;

void pump() {
    MSG m;
    while (PeekMessageW(&m, nullptr, 0, 0, PM_REMOVE)) {
        TranslateMessage(&m);
        DispatchMessageW(&m);
    }
}

bool pump_until(const std::function<bool()>& done, double seconds = 600.0) {
    const auto end = std::chrono::steady_clock::now() + std::chrono::duration<double>(seconds);
    MSG m;
    while (!done()) {
        if (std::chrono::steady_clock::now() > end) return false;
        if (PeekMessageW(&m, nullptr, 0, 0, PM_REMOVE)) {
            TranslateMessage(&m);
            DispatchMessageW(&m);
        } else {
            MsgWaitForMultipleObjects(0, nullptr, FALSE, 5, QS_ALLINPUT);
        }
    }
    return true;
}

std::string read_file(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

// A lesson over the planar MOSFET: set Na, sweep Na, then a step whose set cannot be taken.
fs::path lesson_examples(const std::string& name) {
    const fs::path dir = app_test::scratch(name);
    const std::string deck = read_file(fs::path(NITCAD_EXAMPLES_DIR) / "devices" / "mosfet_2d.deck");
    fs::create_directories(dir / "lessons");
    app_test::write(dir / "lessons", "test_lesson.deck",
                    "# lesson: Test lesson\n# track: basics\n# order: 1\n"
                    "# step: Set it\n# do: set Na = 3e17\n"
                    "# step: Sweep it\n# do: sweep Na = 1e16, 1e17\n# show: results\n# report: Threshold (max gm)\n"
                    "# step: A bad value\n# do: set Na = abc\n" +
                        deck);
    return dir;
}

std::unique_ptr<MainWindow> guide_window(const fs::path& examples) {
    WindowOptions o;
    o.adapter = NiTCAD::render::Adapter::warp;
    o.show = false;
    o.quit_on_close = false;
    o.examples_dir = examples;
    o.examples_copy_dir = examples / "copies";
    auto w = MainWindow::create(o);
    REQUIRE(w);
    return w;
}

}  // namespace

TEST_CASE("guide window: View > Guide toggles the pane", "[guide_window]") {
    auto w = guide_window(lesson_examples("guide_toggle"));
    CHECK(w->guide() == nullptr);
    w->execute(Command::guide);
    REQUIRE(w->guide() != nullptr);
    CHECK((GetWindowLongPtrW(w->guide()->hwnd(), GWL_STYLE) & WS_VISIBLE) != 0);  // the test window itself is hidden
    CHECK(w->guide()->lessons().size() == 1);
    w->execute(Command::guide);
    CHECK((GetWindowLongPtrW(w->guide()->hwnd(), GWL_STYLE) & WS_VISIBLE) == 0);
}

TEST_CASE("guide window: a lesson opens as a copy", "[guide_window]") {
    const fs::path ex = lesson_examples("guide_copy");
    const std::string original = read_file(ex / "lessons" / "test_lesson.deck");
    auto w = guide_window(ex);
    w->execute(Command::guide);
    w->guide_open(0);
    REQUIRE(pump_until([&] { return !w->runner().busy(); }));  // the deck's build after opening it
    CHECK(w->document().path() == ex / "copies" / "test_lesson.deck");
    CHECK(w->guide()->heading_text() == "Set it");
    CHECK(read_file(ex / "lessons" / "test_lesson.deck") == original);
}

TEST_CASE("guide window: a set edits one line, and Undo takes it back", "[guide_window]") {
    auto w = guide_window(lesson_examples("guide_set"));
    w->execute(Command::guide);
    w->guide_open(0);
    REQUIRE(pump_until([&] { return !w->runner().busy(); }));  // the deck's build after opening it
    const std::string before = w->deck().text();
    REQUIRE(before.find("Na = 1e17") != std::string::npos);
    w->guide()->act();
    const std::string after = w->deck().text();
    CHECK(after.find("Na = 3e17") != std::string::npos);
    CHECK(after.size() == before.size());  // one value of the same length changed, nothing else
    w->execute(Command::tab_geometry);
    pump();
    w->execute(Command::undo);
    CHECK(w->deck().text() == before);
}

TEST_CASE("guide window: an action's error keeps the step", "[guide_window]") {
    auto w = guide_window(lesson_examples("guide_error"));
    w->execute(Command::guide);
    w->guide_open(0);
    REQUIRE(pump_until([&] { return !w->runner().busy(); }));  // the deck's build after opening it
    const std::string before = w->deck().text();
    w->guide()->go(2);
    w->guide()->act();
    CHECK(w->guide()->step() == 2);
    CHECK(w->guide()->report_text().find("Na") != std::string::npos);
    CHECK(w->deck().text() == before);
}

TEST_CASE("guide window: a sweep runs, the action waits for it, the report follows (GU5)", "[.guide_run]") {
    auto w = guide_window(lesson_examples("guide_sweep"));
    w->execute(Command::guide);
    w->guide_open(0);
    REQUIRE(pump_until([&] { return !w->runner().busy(); }));  // the deck's build after opening it
    w->guide()->go(1);
    w->guide()->act();
    REQUIRE(w->runner().sweep());
    CHECK(w->runner().sweep()->name == "Na");
    CHECK(w->runner().busy());
    CHECK_FALSE(IsWindowEnabled(w->guide()->action_button()));
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->sweep_results() != nullptr; }));
    pump();
    CHECK(IsWindowEnabled(w->guide()->action_button()));
    CHECK(w->guide()->report_text().starts_with("Threshold (max gm): "));
    CHECK(w->tab() == MainWindow::tab_results);
}

// --- final review fixes ---------------------------------------------------------------------------

TEST_CASE("guide pane: no lesson is shown until one is picked", "[guide_window]") {
    Parent p;
    Recorder r;
    GuidePane g(p.hwnd, r, 7);
    g.set_lessons(two_lessons());
    CHECK(g.current() == nullptr);
    CHECK(g.heading_text() == "Pick a lesson");
    CHECK_FALSE(IsWindowEnabled(g.action_button()));
    g.act();
    g.go(1);
    CHECK(r.acted.empty());
    HWND picker = FindWindowExW(g.hwnd(), nullptr, L"COMBOBOX", nullptr);
    CHECK(SendMessageW(picker, CB_GETCURSEL, 0, 0) == CB_ERR);  // so picking the first lesson is a change
}

TEST_CASE("guide window: a step refuses to act on a deck that is not the lesson's", "[guide_window]") {
    const fs::path ex = lesson_examples("guide_other_deck");
    auto w = guide_window(ex);
    w->execute(Command::guide);
    w->guide_open(0);
    REQUIRE(pump_until([&] { return !w->runner().busy(); }));
    // The learner opens a deck of their own; the pane still shows the lesson's step.
    const fs::path mine = app_test::write(ex, "mine.deck", read_file(fs::path(NITCAD_EXAMPLES_DIR) / "devices" / "mosfet_2d.deck"));
    REQUIRE(w->open_deck(mine));
    REQUIRE(pump_until([&] { return !w->runner().busy(); }));
    const std::string before = w->deck().text();
    w->guide()->act();  // "Set Na = 3e17"
    CHECK(w->deck().text() == before);
    CHECK(w->guide()->report_text().find("lesson") != std::string::npos);
}

TEST_CASE("guide window: a sweep with a failed run reports the runs that finished (GU5)", "[.guide_run]") {
    const fs::path dir = app_test::scratch("guide_sweep_fail");
    fs::create_directories(dir / "lessons");
    const std::string deck = read_file(fs::path(NITCAD_EXAMPLES_DIR) / "devices" / "mosfet_2d.deck");
    // tox = 0: the oxide straddle leaves the axis and that run's deck does not build.
    app_test::write(dir / "lessons", "fail.deck",
                    "# lesson: Fails\n# track: basics\n# order: 1\n# step: Sweep\n# Sweep it.\n"
                    "# do: sweep tox = 10nm, 0nm\n# report: Threshold (max gm)\n" + deck);
    auto w = guide_window(dir);
    w->execute(Command::guide);
    w->guide_open(0);
    REQUIRE(pump_until([&] { return !w->runner().busy(); }));
    w->guide()->act();
    REQUIRE(pump_until([&] { return !w->runner().busy(); }));
    pump();
    const std::string report = w->guide()->report_text();
    INFO(report);
    CHECK(report.starts_with("Threshold (max gm): "));
    CHECK(report.find(" \xE2\x86\x92 \xE2\x80\x94") != std::string::npos);  // " → —" for the run that failed
}
