// The main window's panes (Unit 26, WG5): they tile the client area, and at twice the DPI every size
// doubles (per-monitor DPI re-layout).
#include <catch2/catch_test_macros.hpp>

#include "NiTCAD/app/layout.hpp"

using namespace NiTCAD;
using app::Panes;
using app::Rect;

namespace {

void check_tiles(const Panes& p, int w, int h) {
    CHECK(p.toolbar.x == 0);
    CHECK(p.toolbar.y == 0);
    CHECK(p.toolbar.width == w);
    CHECK(p.status.bottom() == h);
    CHECK(p.status.width == w);
    // Left to right: tree, splitter, centre, splitter, inspector; top to bottom in the centre: tab strip,
    // page, metrics.
    CHECK(p.tree.x == 0);
    CHECK(p.left_splitter.x == p.tree.right());
    CHECK(p.tabs.x == p.left_splitter.right());
    CHECK(p.right_splitter.x == p.tabs.right());
    CHECK(p.inspector.x == p.right_splitter.right());
    CHECK(p.inspector.right() == w);
    for (const Rect& r : {p.tree, p.tabs, p.inspector, p.left_splitter, p.right_splitter}) {
        CHECK(r.y == p.toolbar.bottom());
        CHECK(r.bottom() == p.status.y);
        CHECK(r.width >= 0);
    }
    CHECK(p.page.x == p.tabs.x);
    CHECK(p.page.width == p.tabs.width);
    CHECK(p.page.y > p.tabs.y);
    CHECK(p.metrics.y == p.page.bottom());
    CHECK(p.metrics.bottom() == p.status.y);
}

}  // namespace

TEST_CASE("app layout: the panes tile the window and scale with the DPI (WG5)", "[app]") {
    CHECK(app::scaled(215, 96) == 215);
    CHECK(app::scaled(215, 192) == 430);
    CHECK(app::scaled(215, 144) == 323);  // 322.5 rounded
    CHECK(app::scaled(-10, 120) == -13);

    const app::Layout defaults;
    const Panes a = app::layout_panes(1280, 800, 96, defaults, 30, 22);
    check_tiles(a, 1280, 800);
    CHECK(a.tree.width == app::tree_width_96);
    CHECK(a.inspector.width == app::inspector_width_96);
    CHECK(a.metrics.height == app::metrics_height_96);
    // The same window on a 200% monitor: every pane twice the size.
    const Panes b = app::layout_panes(2560, 1600, 192, defaults, 60, 44);
    check_tiles(b, 2560, 1600);
    for (const auto& [r, s] : {std::pair{a.tree, b.tree}, std::pair{a.tabs, b.tabs}, std::pair{a.page, b.page},
                               std::pair{a.metrics, b.metrics}, std::pair{a.inspector, b.inspector},
                               std::pair{a.status, b.status}, std::pair{a.toolbar, b.toolbar}}) {
        CHECK(s.x == 2 * r.x);
        CHECK(s.y == 2 * r.y);
        CHECK(s.width == 2 * r.width);
        CHECK(s.height == 2 * r.height);
    }
    // Saved widths (96-DPI pixels) are scaled too.
    app::Layout saved;
    saved.tree_width = 300;
    saved.inspector_width = 200;
    const Panes c = app::layout_panes(2560, 1600, 192, saved, 60, 44);
    CHECK(c.tree.width == 600);
    CHECK(c.inspector.width == 400);
    // A narrow window narrows the sides first, keeping the centre's minimum, then everything: at 600
    // pixels the default sides (215 + 260 + 2 splitters) would leave the centre 117.
    const Panes d = app::layout_panes(600, 500, 96, defaults, 30, 22);
    check_tiles(d, 600, 500);
    CHECK(d.tabs.width == app::min_centre_96);
    CHECK(d.tree.width >= app::min_side_96);
    CHECK(d.inspector.width >= app::min_side_96);
    CHECK(d.tree.width < app::tree_width_96);
    const Panes e = app::layout_panes(300, 200, 96, defaults, 30, 22);
    check_tiles(e, 300, 200);
    const Panes f = app::layout_panes(0, 0, 96, defaults, 30, 22);
    CHECK(f.page.width >= 0);
    CHECK(f.page.height >= 0);
}
