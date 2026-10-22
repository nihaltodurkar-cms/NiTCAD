// Where the main window's panes go (Unit 26, DECISIONS.md VA1, WG5): the tool bar on top, the status
// bar at the bottom, the scene tree on the left and the inspector on the right with splitters between
// them and the centre, the tab strip over the page, and the metrics strip under the device view.
// Sizes are given at 96 DPI and scaled to the monitor's DPI, so a per-monitor DPI change re-lays the
// window out at the new scale (WG5).
#pragma once

#include "NiTCAD/app/app_state.hpp"

namespace NiTCAD::app {

struct Rect {
    int x = 0, y = 0, width = 0, height = 0;
    [[nodiscard]] int right() const noexcept { return x + width; }
    [[nodiscard]] int bottom() const noexcept { return y + height; }
};

struct Panes {
    Rect toolbar, tree, tabs, page, metrics, inspector, status;
    Rect left_splitter, right_splitter;
};

// Pixels at 96 DPI to the monitor's (rounded to nearest).
[[nodiscard]] int scaled(int pixels_at_96, int dpi) noexcept;

// Default and minimum sizes at 96 DPI.
inline constexpr int tree_width_96 = 215, inspector_width_96 = 260, splitter_96 = 4, metrics_height_96 = 46,
                     tab_strip_96 = 28, min_side_96 = 120, min_centre_96 = 200;

// The panes for a client area at a DPI. toolbar_height and status_height are the controls' own heights
// (in pixels at this DPI); the tree and inspector widths come from the saved layout (96-DPI pixels, 0 the
// default), narrowed so the centre keeps its minimum width while the sides keep theirs.
[[nodiscard]] Panes layout_panes(int client_width, int client_height, int dpi, const Layout& saved,
                                 int toolbar_height, int status_height);

}  // namespace NiTCAD::app
