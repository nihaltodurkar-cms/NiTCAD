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
