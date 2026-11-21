// Unit 26m slice 0 (DECISIONS.md MT1-MT3): the modeller's way of working without a window: what a press does to the
// view (MC2's rule), where the floating bar goes (MC3's table: inside the viewport, off the selection, at 96, 144 and
// 192 dpi), the deck line between two texts (MC4's function), the typed value, the groups.
#include <catch2/catch_test_macros.hpp>

#include <random>
#include <string>

#include "NiTCAD/app/modeller.hpp"

using namespace NiTCAD::app;
using namespace NiTCAD::app::modeller;

TEST_CASE("modeller: what a press does to the view, camera free and fixed (MT1, MQ3)", "[app][modeller]") {
    // Fixed: nothing moves the view; the left button always goes to the model.
    for (const bool three_d : {false, true}) {
        for (const bool background : {false, true}) {
            for (const bool shift : {false, true}) {
                for (const bool select : {false, true}) {
                    CHECK(navigation(true, three_d, {Button::left, background, shift, select, false}) == Nav::model);
                    CHECK(navigation(true, three_d, {Button::left, background, shift, select, true}) == Nav::none);
                    CHECK(navigation(true, three_d, {Button::middle, background, shift, select, false}) == Nav::none);
                    CHECK(navigation(true, three_d, {Button::right, background, shift, select, false}) == Nav::none);
                }
            }
        }
    }
    CHECK_FALSE(wheel_zooms(true));
    CHECK(wheel_zooms(false));
    // Free, 3D: background orbits, Shift pans; on the model, or with a drawing tool, the model.
    CHECK(navigation(false, true, {Button::left, true, false, true, false}) == Nav::orbit);
    CHECK(navigation(false, true, {Button::left, true, true, true, false}) == Nav::pan);
    CHECK(navigation(false, true, {Button::left, false, false, true, false}) == Nav::model);
    CHECK(navigation(false, true, {Button::left, true, false, false, false}) == Nav::model);
    CHECK(navigation(false, true, {Button::left, false, false, true, true}) == Nav::pan);
    CHECK(navigation(false, true, {Button::middle}) == Nav::pan);
    CHECK(navigation(false, true, {Button::right}) == Nav::orbit);
    // Free, 2D: background pans, Shift keeps the rubber band; no right-drag orbit.
    CHECK(navigation(false, false, {Button::left, true, false, true, false}) == Nav::pan);
    CHECK(navigation(false, false, {Button::left, true, true, true, false}) == Nav::model);
    CHECK(navigation(false, false, {Button::left, false, false, true, false}) == Nav::model);
    CHECK(navigation(false, false, {Button::middle}) == Nav::pan);
    CHECK(navigation(false, false, {Button::right}) == Nav::none);
}

TEST_CASE("modeller: the floating bar lies inside the viewport, off the selection, at 96, 144 and 192 dpi (MC3)", "[app][modeller]") {
    // Nothing picked: the corner.
    const Box vp{0, 30, 800, 630};
    const Placement none = place_bar(vp, std::nullopt, 300, 60, 8);
    CHECK(none.corner);
    CHECK(none.bar == Box{8, 38, 308, 98});
    // Each side in turn, as the selection moves to the viewport's edges.
    CHECK(place_bar(vp, Box{100, 100, 200, 200}, 300, 60, 8).bar == Box{208, 100, 508, 160});  // right
    CHECK(place_bar(vp, Box{600, 100, 780, 200}, 300, 60, 8).bar == Box{492, 208, 792, 268});  // below, slid inside
    CHECK(place_bar(vp, Box{600, 100, 780, 600}, 300, 60, 8).bar == Box{292, 100, 592, 160});  // left
    CHECK(place_bar(vp, Box{100, 500, 780, 600}, 300, 60, 8).bar == Box{100, 432, 400, 492});  // above
    // Filling the view: the corner.
    CHECK(place_bar(vp, Box{-50, 0, 900, 700}, 300, 60, 8).corner);

    // The table: random selections, three dpi; every placement inside, and off the selection unless in the corner,
    // which only happens when no side has room.
    std::mt19937 rng(2611);
    int corners = 0, total = 0;
    for (const int dpi : {96, 144, 192}) {
        const double k = dpi / 96.0, w = 300 * k, h = 62 * k, gap = 8 * k;
        const Box view{0, 32 * k, 1000 * k, 700 * k};
        std::uniform_real_distribution<double> x(-0.1 * view.right, 1.1 * view.right), y(-0.1 * view.bottom, 1.1 * view.bottom);
        for (int n = 0; n < 2000; ++n) {
            double a = x(rng), b = x(rng), c = y(rng), d = y(rng);
            if (a > b) std::swap(a, b);
            if (c > d) std::swap(c, d);
            const Box s{a, c, b, d};
            const Placement p = place_bar(view, s, w, h, gap);
            ++total;
            INFO("dpi " << dpi << " selection " << a << ' ' << c << ' ' << b << ' ' << d);
            CHECK(p.bar.inside(view));
            if (p.corner) {
                ++corners;
                // No side had room: each candidate's size beside the selection does not fit inside the viewport.
                const bool rows = h <= view.height() - 2 * gap, cols = w <= view.width() - 2 * gap;  // the slide keeps it in
                const bool right = rows && s.right + gap >= view.left && s.right + gap + w <= view.right;
                const bool left = rows && s.left - gap - w >= view.left && s.left - gap <= view.right;
                const bool below = cols && s.bottom + gap >= view.top && s.bottom + gap + h <= view.bottom;
                const bool above = cols && s.top - gap - h >= view.top && s.top - gap <= view.bottom;
                CHECK_FALSE(right);
                CHECK_FALSE(left);
                CHECK_FALSE(below);
                CHECK_FALSE(above);
            } else {
                CHECK_FALSE(p.bar.overlaps(s));
            }
        }
    }
    INFO("corner placements " << corners << " of " << total);
    CHECK(corners < total / 2);
}

TEST_CASE("modeller: the deck line is what the edit writes (MC4's function)", "[app][modeller]") {
    const std::string before = "[mesh]\r\ntype = generate3d\r\n\r\n[region ox]\r\nmaterial = sio2\r\nbox = 0 1um 0 1um 0 1um\r\n\r\n[equilibrium]\r\n";
    // A changed line.
    std::string after = before;
    after.replace(after.find("box = 0 1um"), 11, "box = 0 2um");
    CHECK(deck_line(before, after) == "box = 0 2um 0 1um 0 1um");
    // A region added: its lines, blank ones dropped.
    after = before;
    after.insert(after.find("[equilibrium]"), "[region silicon_1]\r\nmaterial = silicon\r\nbox = 0 1um 0 1um 0 0.5um\r\n\r\n");
    CHECK(deck_line(before, after) == "[region silicon_1]  \xC2\xB7  material = silicon  \xC2\xB7  box = 0 1um 0 1um 0 0.5um");
    // Lines removed only.
    CHECK(deck_line(after, before) == "removed: [region silicon_1]  \xC2\xB7  material = silicon  \xC2\xB7  box = 0 1um 0 1um 0 0.5um");
    CHECK(deck_line(before, before).empty());
    CHECK(deck_line("", "a = 1\n") == "a = 1");
}

TEST_CASE("modeller: the typed value, with or without a unit", "[app][modeller]") {
    Entry e;
    CHECK(e.empty());
    CHECK_FALSE(e.value());
    CHECK_FALSE(e.type(L'n'));  // a unit before a digit
    for (const wchar_t c : std::wstring(L"0.25")) CHECK(e.type(c));
    CHECK(e.value() == 2500);  // quanta: 0.1 nm
    CHECK(e.type(L'n'));
    CHECK_FALSE(e.value());  // "0.25n": a unit half typed
    CHECK(e.type(L'm'));
    CHECK(e.value() == 3);  // 0.25 nm: 2.5 quanta, rounded half away from zero as the deck's own parse
    CHECK(e.backspace());
    CHECK(e.backspace());
    CHECK(e.value() == 2500);
    CHECK_FALSE(e.type(L'x'));
    e.clear();
    for (const wchar_t c : std::wstring(L"-20nm")) CHECK(e.type(c));
    CHECK(e.value() == -200);
    e.clear();
    CHECK_FALSE(e.backspace());
    CHECK(chip_length(0) == "0 \xC2\xB5m");
    CHECK(chip_length(1200) == "0.12 \xC2\xB5m");
}

TEST_CASE("modeller: what is picked opens its group; each group's tools", "[app][modeller]") {
    CHECK(group_for(Picked::nothing) == Group::solid);
    CHECK(group_for(Picked::solid) == Group::transform);
    CHECK(group_for(Picked::edge) == Group::edges);
    CHECK(group_for(Picked::face) == Group::faces);
    CHECK(group_for(Picked::contact) == Group::faces);
    for (int g = 0; g < group_count; ++g) {
        const auto list = tools(static_cast<Group>(g));
        CHECK(list.size() <= 5);  // at most five tools showing
        CHECK_FALSE(group_name(static_cast<Group>(g)).empty());
    }
    CHECK(tools(Group::solid).size() == 3);
    CHECK(tools(Group::edges).empty());
    CHECK_FALSE(later(Group::edges).empty());
}
