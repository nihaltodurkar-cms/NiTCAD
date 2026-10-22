// The main window's panes (layout.hpp).
#include "NiTCAD/app/layout.hpp"

#include <algorithm>

namespace NiTCAD::app {

int scaled(int pixels_at_96, int dpi) noexcept {
    const long long v = static_cast<long long>(pixels_at_96) * dpi;
    return static_cast<int>(v >= 0 ? (v + 48) / 96 : -((-v + 48) / 96));
}

Panes layout_panes(int client_width, int client_height, int dpi, const Layout& saved, int toolbar_height,
                   int status_height) {
    Panes p;
    const int w = std::max(client_width, 0), h = std::max(client_height, 0);
    p.toolbar = {0, 0, w, std::min(toolbar_height, h)};
    p.status = {0, std::max(h - status_height, p.toolbar.bottom()), w, 0};
    p.status.height = h - p.status.y;
    const int top = p.toolbar.bottom(), bottom = p.status.y;
    const int body = std::max(bottom - top, 0);

    const int split = scaled(splitter_96, dpi);
    const int min_side = scaled(min_side_96, dpi), min_centre = scaled(min_centre_96, dpi);
    int tree = scaled(saved.tree_width > 0 ? saved.tree_width : tree_width_96, dpi);
    int inspector = scaled(saved.inspector_width > 0 ? saved.inspector_width : inspector_width_96, dpi);
    // Narrow the sides (each down to its minimum) until the centre has its minimum; then they may
    // shrink further when the window itself is narrower than that.
    const auto centre = [&] { return w - tree - inspector - 2 * split; };
    if (centre() < min_centre) {
        const int excess = min_centre - centre();
        const int from_tree = std::min(excess / 2, std::max(tree - min_side, 0));
        tree -= from_tree;
        inspector -= std::min(excess - from_tree, std::max(inspector - min_side, 0));
        if (centre() < min_centre) tree -= std::min(min_centre - centre(), std::max(tree - min_side, 0));
    }
    if (centre() < 0) {
        const int total = std::max(w - 2 * split, 0);
        tree = total / 2;
        inspector = total - tree;
    }
    p.tree = {0, top, tree, body};
    p.left_splitter = {tree, top, split, body};
    p.inspector = {w - inspector, top, inspector, body};
    p.right_splitter = {p.inspector.x - split, top, split, body};
    const int cx = p.left_splitter.right(), cw = std::max(p.right_splitter.x - cx, 0);
    p.tabs = {cx, top, cw, body};
    const int strip = std::min(scaled(tab_strip_96, dpi), body);
    const int metrics = std::min(scaled(metrics_height_96, dpi), body - strip);
    p.page = {cx, top + strip, cw, body - strip - metrics};
    p.metrics = {cx, p.page.bottom(), cw, metrics};
    return p;
}

}  // namespace NiTCAD::app
