// Unit 26m slice 0 (DECISIONS.md MT1-MT3), driven by window messages (no screen automation): the camera lock in the 3D
// modeller and the 2D editor (MC2), the floating bar's place and the groups picks open (MC3), the deck line shown
// during a drag against the line the edit writes, and typed values (MC4), the cost of the bar and the deck line
// with 100 solids (MC9).
#include <catch2/catch_test_macros.hpp>

#include "main_window.hpp"
#include "fixtures.hpp"

#include <windows.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <format>
#include <functional>
#include <numbers>
#include <optional>
#include <tuple>
#include <fstream>
#include <string>
#include <vector>

#include "NiTCAD/app/modeller.hpp"
#include "NiTCAD/app/solid.hpp"
#include "NiTCAD/app/solid_view.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/solve/equilibrium.hpp"
#include "win32.hpp"

using namespace NiTCAD;
using namespace NiTCAD::app;

namespace {

void pump() {
    MSG m;
    while (PeekMessageW(&m, nullptr, 0, 0, PM_REMOVE)) {
        TranslateMessage(&m);
        DispatchMessageW(&m);
    }
}

std::unique_ptr<MainWindow> window() {
    app::WindowOptions o;
    o.adapter = render::Adapter::warp;
    o.show = false;
    o.quit_on_close = false;
    auto w = MainWindow::create(o);
    REQUIRE(w);
    return w;
}

std::filesystem::path write(const std::string& name, const std::string& text) {
    const auto dir = app_test::scratch("modeller_window_" + name);
    const auto p = dir / (name + ".deck");
    std::ofstream(p, std::ios::binary) << text;
    return p;
}

bool has(const std::string& text, const std::string& what) { return text.find(what) != std::string::npos; }

const std::string mos = R"([mesh]
type = generate3d
layered = on
domain = 0 1um 0 0.6um 0 0.6um
patch gate = 0.3um 0 0.6um  0.7um 0.6um 0.6um
patch body = 0 0 0  1um 0.6um 0
max_size = 0.15um
interface_resolution = 0.01um

[region si]
material = silicon

[region ox]
material = sio2
box = 0 1um 0 0.6um 0.4um 0.6um

[doping]
profile = acceptor constant 1e16

[contact body]
kind = ohmic
patch = body

[contact gate]
kind = electrode
patch = gate
electrode = metal
work_function = 4.4

[equilibrium]
)";

const std::string layered2d = R"([mesh]
type = generate2d
layered = on
domain = 0 2um 0 1um
patch source = 0 0, 0 0.6um
max_size = 0.25um
interface_resolution = 20nm

[region p]
material = silicon

[region oxide]
material = sio2
box = 0 2um 0.8um 1um

[doping]
profile = acceptor constant 1e16

[contact source]
kind = ohmic
patch = source

[equilibrium]
)";

LPARAM lp(double x, double y) { return MAKELPARAM(static_cast<int>(std::lround(x)), static_cast<int>(std::lround(y))); }
void send(HWND h, UINT m, WPARAM w, double x, double y) { SendMessageW(h, m, w, lp(x, y)); }
void wheel(HWND canvas, double x, double y, int delta) {
    POINT p{static_cast<LONG>(x), static_cast<LONG>(y)};
    ClientToScreen(canvas, &p);
    SendMessageW(canvas, WM_MOUSEWHEEL, MAKEWPARAM(0, delta), MAKELPARAM(p.x, p.y));
}
void type(HWND canvas, const std::wstring& s) {
    for (const wchar_t c : s) SendMessageW(canvas, WM_CHAR, c, 0);
}

bool same(const render::Camera& a, const render::Camera& b) {
    return std::memcmp(&a.target, &b.target, sizeof a.target) == 0 && std::memcmp(&a.azimuth_deg, &b.azimuth_deg, sizeof(float)) == 0 &&
           std::memcmp(&a.elevation_deg, &b.elevation_deg, sizeof(float)) == 0 && std::memcmp(&a.extent_um, &b.extent_um, sizeof(float)) == 0 &&
           a.projection == b.projection && std::memcmp(&a.zoom, &b.zoom, sizeof(float)) == 0 && std::memcmp(&a.pan, &b.pan, sizeof a.pan) == 0;
}

// A pixel of the canvas where the pick finds `want` (scanning a grid), or none.
std::optional<std::array<double, 2>> find_pick(const SolidPage& page, const std::function<bool(const solid::Pick&)>& want) {
    RECT c{};
    GetClientRect(page.canvas(), &c);
    for (int y = 2; y < c.bottom; y += 3) {
        for (int x = 2; x < c.right; x += 3) {
            if (want(solid::pick(page.model(), page.view(), x, y, page.selection()))) return std::array<double, 2>{double(x), double(y)};
        }
    }
    return std::nullopt;
}

struct Opened {
    std::unique_ptr<MainWindow> w;
    SolidPage* page = nullptr;
};
Opened open_3d(const std::string& name, const std::string& deck) {
    Opened o{window()};
    REQUIRE(o.w->open_deck(write(name, deck)));
    pump();
    o.w->execute(Command::tab_geometry);
    REQUIRE(o.w->solid_shown());
    o.page = &o.w->solid_page();
    REQUIRE(o.page->render() >= 0.0);
    return o;
}

// Draws the box of the 26l test (0 0.2um 0 0.3um 0.4um 0.5um) and leaves it selected.
solid::SolidRef draw_box(MainWindow& w, SolidPage& page) {
    page.set_work_plane({2, 4 * solid::per_um / 10});
    page.set_tool(SolidPage::Tool::box);
    const auto a = page.screen_of_um({0.0, 0.0, 0.4}), b = page.screen_of_um({0.2, 0.3, 0.4});
    send(page.canvas(), WM_LBUTTONDOWN, MK_LBUTTON, a[0], a[1]);
    send(page.canvas(), WM_MOUSEMOVE, MK_LBUTTON, (a[0] + b[0]) / 2, (a[1] + b[1]) / 2);
    send(page.canvas(), WM_MOUSEMOVE, MK_LBUTTON, b[0], b[1]);
    send(page.canvas(), WM_LBUTTONUP, 0, b[0], b[1]);
    const auto h = page.screen_of_um({0.2, 0.3, 0.5});
    send(page.canvas(), WM_MOUSEMOVE, 0, h[0], h[1]);
    send(page.canvas(), WM_LBUTTONDOWN, MK_LBUTTON, h[0], h[1]);
    send(page.canvas(), WM_LBUTTONUP, 0, h[0], h[1]);
    REQUIRE(has(w.deck().text(), "box = 0 0.2um 0 0.3um 0.4um 0.5um\n"));
    REQUIRE(page.selection().solids.size() == 1);
    page.set_tool(SolidPage::Tool::select);
    page.render();
    return page.selection().solids[0];
}

// The face handle's pixel (x_max: face 1) of the selected solid, checked by the pick.
std::array<double, 2> face_handle(const SolidPage& page, solid::SolidRef s, int face) {
    const auto hs = solid::handles(page.model(), page.view(), s);
    const auto h = std::ranges::find_if(hs, [&](const solid::Handle& k) { return k.kind == solid::HandleKind::face && k.index == face; });
    REQUIRE(h != hs.end());
    auto p = page.screen_of_um(h->pts[0]);
    p = {std::round(p[0]), std::round(p[1])};
    const solid::Pick k = solid::pick(page.model(), page.view(), p[0], p[1], page.selection());
    REQUIRE(k.kind == solid::PickKind::handle);
    REQUIRE(k.handle == solid::HandleKind::face);
    return p;
}

}  // namespace

TEST_CASE("app window modeller: the 3D camera free and fixed, Reset view, the state kept over a reload (MC2)", "[app][window][modeller]") {
    auto [w, pp] = open_3d("camera3d", mos);
    SolidPage& page = *pp;
    const HWND c = page.canvas();
    CHECK_FALSE(page.camera_fixed());  // MQ2: free on open
    const render::Camera start = page.view().camera;
    // A background pixel: the ray misses the domain.
    const auto bg = find_pick(page, [](const solid::Pick& p) { return p.kind == solid::PickKind::none; });
    REQUIRE(bg);
    const double x = (*bg)[0], y = (*bg)[1];

    // Free: a left drag on the background orbits; the wheel zooms; Reset view restores the start exactly.
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, x, y);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, x + 30, y + 12);
    send(c, WM_LBUTTONUP, 0, x + 30, y + 12);
    CHECK(page.view().camera.azimuth_deg != start.azimuth_deg);
    CHECK(page.view().camera.elevation_deg != start.elevation_deg);
    wheel(c, 300, 200, 120);
    CHECK(page.view().camera.zoom != start.zoom);
    // Shift+left on the background pans.
    const auto pan0 = page.view().camera.pan;
    send(c, WM_LBUTTONDOWN, MK_LBUTTON | MK_SHIFT, x, y);
    send(c, WM_MOUSEMOVE, MK_LBUTTON | MK_SHIFT, x + 20, y + 10);
    send(c, WM_LBUTTONUP, 0, x + 20, y + 10);
    CHECK(page.view().camera.pan != pan0);
    page.reset_view();
    CHECK(same(page.view().camera, start));

    // Fixed (L): left, middle and right drags and the wheel leave the camera bit-identical.
    SendMessageW(c, WM_KEYDOWN, 'L', 0);
    REQUIRE(page.camera_fixed());
    page.render();
    const render::Camera fixed = page.view().camera;
    for (const auto& [down, up, keys] : {std::tuple{WM_LBUTTONDOWN, WM_LBUTTONUP, MK_LBUTTON}, std::tuple{WM_MBUTTONDOWN, WM_MBUTTONUP, MK_MBUTTON},
                                        std::tuple{WM_RBUTTONDOWN, WM_RBUTTONUP, MK_RBUTTON}}) {
        send(c, down, keys, x, y);
        send(c, WM_MOUSEMOVE, keys, x + 40, y + 25);
        send(c, up, 0, x + 40, y + 25);
        CHECK(same(page.view().camera, fixed));
        send(c, down, keys | MK_SHIFT, x, y);
        send(c, WM_MOUSEMOVE, keys | MK_SHIFT, x + 40, y + 25);
        send(c, up, 0, x + 40, y + 25);
        CHECK(same(page.view().camera, fixed));
    }
    wheel(c, 300, 200, 120);
    wheel(c, 300, 200, -240);
    CHECK(same(page.view().camera, fixed));
    page.set_tool(SolidPage::Tool::pan);
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, x, y);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, x + 40, y + 25);
    send(c, WM_LBUTTONUP, 0, x + 40, y + 25);
    CHECK(same(page.view().camera, fixed));
    page.set_tool(SolidPage::Tool::select);

    // Fixed, a left drag on a face handle still resizes the solid.
    const solid::SolidRef box = draw_box(*w, page);
    CHECK(same(page.view().camera, fixed));
    const auto f = face_handle(page, box, 1);
    const auto to = page.screen_of_um({0.3, 0.15, 0.45});
    const auto from = page.screen_of_um({0.2, 0.15, 0.45});
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, f[0], f[1]);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, f[0] + (to[0] - from[0]) / 2, f[1] + (to[1] - from[1]) / 2);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, f[0] + to[0] - from[0], f[1] + to[1] - from[1]);
    send(c, WM_LBUTTONUP, 0, f[0] + to[0] - from[0], f[1] + to[1] - from[1]);
    CHECK(has(w->deck().text(), "box = 0 0.3um 0 0.3um 0.4um 0.5um\n"));
    CHECK(same(page.view().camera, fixed));

    // The state survives a reload of the deck and the tab shown again.
    page.reload(true);
    w->execute(Command::tab_geometry);
    CHECK(page.camera_fixed());
    SendMessageW(c, WM_KEYDOWN, 'L', 0);
    CHECK_FALSE(page.camera_fixed());
}

TEST_CASE("app window modeller: the 2D camera free and fixed; Shift keeps the rubber band (MC2, MQ6)", "[app][window][modeller]") {
    auto w = window();
    REQUIRE(w->open_deck(write("camera2d", layered2d)));
    pump();
    w->execute(Command::tab_geometry);
    GeometryPage& page = w->geometry_page();
    REQUIRE(page.geometry().editable());
    page.set_view(geo::View{200.0, 50.0, 450.0, 96});
    const HWND c = page.canvas();
    CHECK_FALSE(page.camera_fixed());
    // Outside the domain (x < 0): background.
    const double x = page.view().sx(-1000), y = page.view().sy(5000);

    // Free: a left drag on the background pans by the drag; the wheel zooms.
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, x, y);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, x + 30, y + 12);
    send(c, WM_LBUTTONUP, 0, x + 30, y + 12);
    CHECK(page.view().origin_x == 80.0);
    CHECK(page.view().origin_y == 462.0);
    wheel(c, 300, 200, 120);
    CHECK(page.view().px_per_um != 200.0);

    // Fixed: no drag or wheel moves the view; a left drag on the background is the rubber band.
    page.set_view(geo::View{200.0, 50.0, 450.0, 96});
    SendMessageW(c, WM_KEYDOWN, 'L', 0);
    REQUIRE(page.camera_fixed());
    for (const auto& [down, up, keys] : {std::tuple{WM_LBUTTONDOWN, WM_LBUTTONUP, MK_LBUTTON}, std::tuple{WM_MBUTTONDOWN, WM_MBUTTONUP, MK_MBUTTON},
                                        std::tuple{WM_RBUTTONDOWN, WM_RBUTTONUP, MK_RBUTTON}}) {
        send(c, down, keys, x, y);
        send(c, WM_MOUSEMOVE, keys, x + 40, y + 25);
        send(c, up, 0, x + 40, y + 25);
        CHECK(page.view().origin_x == 50.0);
        CHECK(page.view().origin_y == 450.0);
    }
    wheel(c, 300, 200, 120);
    CHECK(page.view().px_per_um == 200.0);
    // Fixed, a drag on the oxide still moves it.
    page.set_tool(GeometryPage::Tool::select);
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, page.view().sx(10000), page.view().sy(9000));
    send(c, WM_MOUSEMOVE, MK_LBUTTON, page.view().sx(10000), page.view().sy(8500));
    send(c, WM_MOUSEMOVE, MK_LBUTTON, page.view().sx(10000), page.view().sy(8000));
    send(c, WM_LBUTTONUP, 0, page.view().sx(10000), page.view().sy(8000));
    CHECK(has(w->deck().text(), "box = 0 2um 0.7um 0.9um\n"));
    CHECK(page.view().origin_x == 50.0);

    // Free again, Shift+drag on the background is the rubber band: the view stays, the oxide is selected.
    SendMessageW(c, WM_KEYDOWN, 'L', 0);
    REQUIRE_FALSE(page.camera_fixed());
    send(c, WM_LBUTTONDOWN, MK_LBUTTON | MK_SHIFT, page.view().sx(-1000), page.view().sy(12000));
    send(c, WM_MOUSEMOVE, MK_LBUTTON | MK_SHIFT, page.view().sx(21000), page.view().sy(6000));
    send(c, WM_LBUTTONUP, 0, page.view().sx(21000), page.view().sy(6000));
    CHECK(page.view().origin_x == 50.0);
    REQUIRE(page.selection().size() == 1);
    CHECK(page.selection()[0].region == 1);
}

TEST_CASE("app window modeller: picks open their groups; the bar lies inside the view and off the selection at three dpi and six "
          "directions (MC3)",
          "[app][window][modeller]") {
    auto [w, pp] = open_3d("bar", mos);
    SolidPage& page = *pp;
    const HWND c = page.canvas();
    CHECK(page.group() == modeller::Group::solid);
    page.set_tool(SolidPage::Tool::box);
    CHECK(page.group() == modeller::Group::solid);
    page.set_tool(SolidPage::Tool::select);
    page.render();
    const solid::SolidRef ox{1, 0};  // the oxide's box, on top: seen from the default view

    // A click on the solid: Transform; on its face handle: Faces; on its edge: Edges; on a contact: Faces; on nothing: Solid.
    const auto click = [&](std::array<double, 2> p) {
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, p[0], p[1]);
        send(c, WM_LBUTTONUP, 0, p[0], p[1]);
        page.render();
    };
    const auto on_box = find_pick(page, [&](const solid::Pick& p) { return p.kind == solid::PickKind::face && p.solid == ox; });
    REQUIRE(on_box);
    click(*on_box);
    CHECK(page.group() == modeller::Group::transform);
    CHECK(page.bar().line().starts_with("Move:"));
    click(face_handle(page, ox, 1));
    CHECK(page.group() == modeller::Group::faces);
    CHECK(page.bar().line().starts_with("Push/Pull:"));
    const auto edge = find_pick(page, [&](const solid::Pick& p) { return p.kind == solid::PickKind::edge; });
    REQUIRE(edge);
    click(*edge);
    CHECK(page.group() == modeller::Group::edges);
    CHECK_FALSE(has(page.bar().line(), "slice"));  // slice 2 gave Edges its tools
    const auto contact = find_pick(page, [&](const solid::Pick& p) { return p.kind == solid::PickKind::contact; });
    REQUIRE(contact);
    click(*contact);
    CHECK(page.group() == modeller::Group::faces);
    const auto none = find_pick(page, [](const solid::Pick& p) { return p.kind == solid::PickKind::none; });
    REQUIRE(none);
    click(*none);
    CHECK(page.group() == modeller::Group::solid);
    // A tool's key opens a group holding it; the switcher's group arms its first tool.
    SendMessageW(c, WM_KEYDOWN, 'C', 0);
    CHECK(page.group() == modeller::Group::faces);
    page.set_group(modeller::Group::solid);
    CHECK(page.tool() == SolidPage::Tool::select);

    // The bar's place for a small box (drawn, selected): three window sizes (CI's screen is 1024 x 768, its view too
    // narrow for the bar's 300 px at 200 % before), three dpi, six directions.
    const solid::SolidRef box = draw_box(*w, page);
    REQUIRE(page.select_region(page.model().regions[static_cast<std::size_t>(box.region)].name));
    int corners = 0, total = 0;
    for (const auto& [ww, wh] : {std::pair{0, 0}, std::pair{1024, 768}, std::pair{800, 600}}) {
        if (ww > 0) SetWindowPos(w->hwnd(), nullptr, 0, 0, ww, wh, SWP_NOMOVE | SWP_NOZORDER | SWP_NOACTIVATE);
        pump();
        for (const int dpi : {96, 144, 192}) {
            if (ww == 800 && dpi == 192) continue;  // a 301 x 195 px view: no room for the bar at 200 % in any layout
            HFONT font = win32::ui_font(dpi), bold = win32::ui_font(dpi, true);
            page.set_dpi(dpi, font, bold);
            for (int v = 0; v < 6; ++v) {
                page.named_view(v);
                REQUIRE(page.render() >= 0.0);
                RECT r{};
                GetClientRect(c, &r);
                const modeller::Box view{0, 0, double(r.right), double(r.bottom)};
                const modeller::Placement p = page.bar_placement();
                const auto s = page.selection_bounds();
                REQUIRE(s);
                INFO("window " << ww << 'x' << wh << " (view " << r.right << 'x' << r.bottom << "), dpi " << dpi << ", view " << v
                               << ": bar " << p.bar.left << ' ' << p.bar.top << ' ' << p.bar.right << ' ' << p.bar.bottom);
                CHECK(p.bar.inside(view));
                // The page places in its panel's coordinates (the canvas below the strip): equal to rounding.
                const modeller::Box e = modeller::place_bar(view, *s, p.bar.width(), p.bar.height(), MulDiv(8, dpi, 96)).bar;
                CHECK(std::max({std::abs(p.bar.left - e.left), std::abs(p.bar.top - e.top), std::abs(p.bar.right - e.right),
                                std::abs(p.bar.bottom - e.bottom)}) <= 1e-9);
                if (p.corner) ++corners;
                else CHECK_FALSE(p.bar.overlaps(*s));
                ++total;
            }
            page.set_dpi(96, nullptr, nullptr);
            DeleteObject(font);
            DeleteObject(bold);
        }
    }
    std::printf("MC3: the bar in the corner for %d of %d placements (the selection leaving no side room)\n", corners, total);
    CHECK(corners < total);
}

TEST_CASE("app window modeller: the deck line shown while dragging is the line written; typed values set a drag exactly (MC4, MT3)",
          "[app][window][modeller]") {
    auto [w, pp] = open_3d("line3d", mos);
    SolidPage& page = *pp;
    const HWND c = page.canvas();

    // Drawing a box: the line at the height's last move is the line the click writes.
    page.set_work_plane({2, 4 * solid::per_um / 10});
    page.set_tool(SolidPage::Tool::box);
    std::string before = w->deck().text();
    const auto a = page.screen_of_um({0.0, 0.0, 0.4}), b = page.screen_of_um({0.2, 0.3, 0.4});
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, a[0], a[1]);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, b[0], b[1]);
    CHECK(has(page.live_line(), "[region silicon_1]"));
    CHECK(has(page.chip(), "0.2 \xC2\xB5m \xC3\x97 0.3 \xC2\xB5m"));
    send(c, WM_LBUTTONUP, 0, b[0], b[1]);
    const auto h = page.screen_of_um({0.2, 0.3, 0.5});
    send(c, WM_MOUSEMOVE, 0, h[0], h[1]);
    const std::string shown = page.live_line();
    CHECK(has(page.chip(), "height 0.1 \xC2\xB5m"));
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, h[0], h[1]);
    send(c, WM_LBUTTONUP, 0, h[0], h[1]);
    std::string after = w->deck().text();
    REQUIRE(has(after, "box = 0 0.2um 0 0.3um 0.4um 0.5um\n"));
    CHECK(shown == modeller::deck_line(before, after));
    CHECK(page.live_line().empty());
    CHECK(page.chip().empty());

    // A face pushed 0.1 µm by its handle: the line while dragging; then a typed 0.05 µm, Enter applying it.
    page.set_tool(SolidPage::Tool::select);
    page.render();
    const solid::SolidRef box = page.selection().solids[0];
    auto f = face_handle(page, box, 1);
    const auto step = page.screen_of_um({0.3, 0.15, 0.45}), at = page.screen_of_um({0.2, 0.15, 0.45});
    const double dx = step[0] - at[0], dy = step[1] - at[1];
    before = after;
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, f[0], f[1]);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, f[0] + dx, f[1] + dy);
    const std::string pushing = page.live_line();
    CHECK(has(page.chip(), "push"));
    send(c, WM_LBUTTONUP, 0, f[0] + dx, f[1] + dy);
    after = w->deck().text();
    CHECK(pushing == modeller::deck_line(before, after));
    CHECK(pushing == "box = 0 0.3um 0 0.3um 0.4um 0.5um");

    page.render();
    f = face_handle(page, box, 1);
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, f[0], f[1]);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, f[0] + 0.3 * dx, f[1] + 0.3 * dy);
    SendMessageW(c, WM_KEYDOWN, 'M', 0);  // a unit letter during a typed drag: not the Measure tool
    type(c, L"0.05");
    CHECK(page.typed().text() == "0.05");
    CHECK(page.tool() == SolidPage::Tool::select);
    CHECK(page.live_line() == "box = 0 0.35um 0 0.3um 0.4um 0.5um");
    SendMessageW(c, WM_KEYDOWN, VK_RETURN, 0);
    send(c, WM_LBUTTONUP, 0, f[0], f[1]);
    CHECK(has(w->deck().text(), "box = 0 0.35um 0 0.3um 0.4um 0.5um\n"));
    // Typed with a unit, Esc cancelling: the deck stays.
    page.render();
    f = face_handle(page, box, 1);
    before = w->deck().text();
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, f[0], f[1]);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, f[0] + dx, f[1] + dy);
    type(c, L"20nm");
    CHECK(page.live_line() == "box = 0 0.37um 0 0.3um 0.4um 0.5um");
    SendMessageW(c, WM_KEYDOWN, VK_ESCAPE, 0);
    send(c, WM_LBUTTONUP, 0, f[0] + dx, f[1] + dy);
    CHECK(w->deck().text() == before);
    CHECK(page.typed().empty());

    // The box's height typed: 0.15 µm.
    page.set_tool(SolidPage::Tool::box);
    page.set_work_plane({2, 0});
    const auto p0 = page.screen_of_um({0.5, 0.1, 0.0}), p1 = page.screen_of_um({0.8, 0.4, 0.0});
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, p0[0], p0[1]);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, p1[0], p1[1]);
    send(c, WM_LBUTTONUP, 0, p1[0], p1[1]);
    const auto ph = page.screen_of_um({0.8, 0.4, 0.05});
    send(c, WM_MOUSEMOVE, 0, ph[0], ph[1]);
    type(c, L"0.15");
    CHECK(has(page.chip(), "height 0.15 \xC2\xB5m"));
    SendMessageW(c, WM_KEYDOWN, VK_RETURN, 0);
    CHECK(has(w->deck().text(), "box = 0.5um 0.8um 0.1um 0.4um 0 0.15um\n"));
}

TEST_CASE("app window modeller: the 2D deck line and typed values (MC4, MQ6)", "[app][window][modeller]") {
    auto w = window();
    REQUIRE(w->open_deck(write("line2d", layered2d)));
    pump();
    w->execute(Command::tab_geometry);
    GeometryPage& page = w->geometry_page();
    page.set_view(geo::View{200.0, 50.0, 450.0, 96});
    const HWND c = page.canvas();
    const auto at = [&](geo::Coord x, geo::Coord y) { return std::array<double, 2>{page.view().sx(x), page.view().sy(y)}; };

    // A rectangle: the line at the last move is the line written.
    page.set_material("sio2");
    page.set_tool(GeometryPage::Tool::rectangle);
    std::string before = w->deck().text();
    auto a = at(4000, 2000), b = at(8000, 5000);
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, a[0], a[1]);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, b[0], b[1]);
    const std::string shown = page.live_line();
    CHECK(has(page.chip(), "0.4 \xC2\xB5m \xC3\x97 0.3 \xC2\xB5m"));
    send(c, WM_LBUTTONUP, 0, b[0], b[1]);
    std::string after = w->deck().text();
    REQUIRE(has(after, "box = 0.4um 0.8um 0.2um 0.5um\n"));
    CHECK(shown == modeller::deck_line(before, after));

    // Typed: the larger offset (x) set to 0.25 µm.
    a = at(10000, 2000), b = at(14000, 5000);
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, a[0], a[1]);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, b[0], b[1]);
    SendMessageW(c, WM_KEYDOWN, 'R', 0);
    type(c, L"0.25");
    CHECK(page.tool() == GeometryPage::Tool::rectangle);
    CHECK(has(page.live_line(), "box = 1um 1.25um 0.2um 0.5um"));
    send(c, WM_LBUTTONUP, 0, b[0], b[1]);
    CHECK(has(w->deck().text(), "box = 1um 1.25um 0.2um 0.5um\n"));
    // Typed while dragging leftwards: the value in the drag's direction.
    a = at(19000, 2000), b = at(16000, 4000);
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, a[0], a[1]);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, b[0], b[1]);
    type(c, L"0.2");
    send(c, WM_LBUTTONUP, 0, b[0], b[1]);
    CHECK(has(w->deck().text(), "box = 1.7um 1.9um 0.2um 0.4um\n"));

    // A move: the line while dragging; typed 0.1 µm along y, Enter applying it.
    page.set_tool(GeometryPage::Tool::select);
    before = w->deck().text();
    a = at(6000, 3500);
    b = at(6000, 4500);
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, a[0], a[1]);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, b[0], b[1]);
    CHECK(page.live_line() == "box = 0.4um 0.8um 0.3um 0.6um");
    type(c, L"0.2");
    CHECK(page.live_line() == "box = 0.4um 0.8um 0.4um 0.7um");
    SendMessageW(c, WM_KEYDOWN, VK_RETURN, 0);
    send(c, WM_LBUTTONUP, 0, b[0], b[1]);
    after = w->deck().text();
    CHECK(has(after, "box = 0.4um 0.8um 0.4um 0.7um\n"));
    CHECK(modeller::deck_line(before, after) == "box = 0.4um 0.8um 0.4um 0.7um");
}

TEST_CASE("app window modeller: the bar and the deck line add at most 1 ms a frame with 100 solids (MC9)", "[app][window][modeller][perf]") {
    std::string deck = "[mesh]\ntype = generate3d\nlayered = on\ndomain = 0 10um 0 10um 0 1um\nmax_size = 0.5um\n\n[region si]\nmaterial = silicon\n";
    for (int k = 0; k < 100; ++k) {
        deck += std::format("\n[region r{}]\nmaterial = {}\nbox = {}um {}.8um {}um {}.8um 0 0.{}um\n", k, k % 2 ? "sio2" : "silicon", k % 10, k % 10,
                            k / 10, k / 10, 3 + k % 7);
    }
    deck += "\n[equilibrium]\n";
    auto [w, pp] = open_3d("hundred", deck);
    SolidPage& page = *pp;
    REQUIRE(page.model().regions.size() == 101);
    REQUIRE(page.select_region("r55"));
    page.render();
    const std::string text = w->deck().text();
    const solid::SolidRef ref = page.selection().solids[0];
    std::vector<double> added;
    for (int k = 1; k <= 31; ++k) {
        solid::Model m = page.model();
        solid::move_solid(m, ref, {k * 100, 0, 0});
        const auto t0 = std::chrono::steady_clock::now();
        const auto out = solid::write_model(text, m);
        REQUIRE(out);
        const std::string line = modeller::deck_line(text, *out);
        const auto s = page.selection_bounds();
        const auto p = modeller::place_bar({0, 0, 800, 600}, s, 300, 62, 8);
        added.push_back(1e3 * std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        REQUIRE(has(line, "box = "));
        REQUIRE(std::abs(p.bar.width() - 300) < 1e-9);
    }
    std::ranges::sort(added);
    std::printf("MC9: the deck line and the bar's place with 100 solids: median %.3f ms (max %.3f)\n", added[added.size() / 2], added.back());
#ifdef NDEBUG
    CHECK(added[added.size() / 2] <= 1.0);
#endif
}

// --- Unit 26m slice 1 (MT4, MT5): the overlap switch and the Transform modes, by window messages ------------------

namespace {

// A pixel where the pick finds face `face` of solid `s` (not a handle), and the hit point (µm).
std::pair<std::array<double, 2>, solid::V3> face_pixel(const SolidPage& page, solid::SolidRef s, int face) {
    const auto px = find_pick(page, [&](const solid::Pick& p) { return p.kind == solid::PickKind::face && p.solid == s && p.index == face; });
    if (!px) {
        std::string seen;
        for (int f = 0; f < 6; ++f) {
            if (find_pick(page, [&](const solid::Pick& p) { return p.kind == solid::PickKind::face && p.solid == s && p.index == f; })) seen += std::to_string(f);
        }
        UNSCOPED_INFO("faces of the solid seen: " << seen << "; solid " << s.region << "," << s.solid);
    }
    REQUIRE(px);
    const solid::Pick p = solid::pick(page.model(), page.view(), (*px)[0], (*px)[1], page.selection());
    return {*px, p.at};
}

// A drag from a face's pixel by `d` µm in the face's plane; the deck line shown at the last move.
std::string drag_by(SolidPage& page, std::array<double, 2> from, solid::V3 at, solid::V3 d, const std::function<void()>& during = {}) {
    const HWND c = page.canvas();
    const auto to = page.screen_of_um({at[0] + d[0], at[1] + d[1], at[2] + d[2]});
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, from[0], from[1]);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, (from[0] + to[0]) / 2, (from[1] + to[1]) / 2);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, to[0], to[1]);
    if (during) during();
    const std::string shown = page.live_line();
    send(c, WM_LBUTTONUP, 0, to[0], to[1]);
    page.render();
    return shown;
}

}  // namespace

TEST_CASE("app window modeller: Copy, Mirror, Align, Array and Fill/Replace in 3D, each one undo step (MC4, MC6, MT4, MT5)",
          "[app][window][modeller][transform]") {
    auto [w, pp] = open_3d("transform3d", mos);
    SolidPage& page = *pp;
    const HWND c = page.canvas();
    const solid::SolidRef a = draw_box(*w, page);  // 0 0.2um 0 0.3um 0.4um 0.5um, region 2
    // Its top raised to the domain's (a face under the oxide is hidden from every pick): 0.4 to 0.6um.
    win32::set_text(page.control(SolidPage::c_z1), "0.6");
    REQUIRE(page.commit(SolidPage::c_z1));
    REQUIRE(has(w->deck().text(), "box = 0 0.2um 0 0.3um 0.4um 0.6um\n"));
    page.render();
    std::vector<std::string> texts{w->deck().text()};
    const auto step = [&](const std::string& shown) {
        const std::string now = w->deck().text();
        CHECK(now != texts.back());
        if (!shown.empty()) CHECK(shown == modeller::deck_line(texts.back(), now));  // MC4
        texts.push_back(now);
    };

    // Copy: the drag moves a copy; the original stays.
    page.set_mode(modeller::Mode::copy);
    CHECK(page.group() == modeller::Group::transform);
    CHECK(page.tool() == SolidPage::Tool::select);
    auto [px, at] = face_pixel(page, a, 5);
    std::string shown = drag_by(page, px, at, {0.8, 0.0, 0.0});  // clear of the gate contact over x 0.3-0.7um
    REQUIRE(has(w->deck().text(), "box = 0 0.2um 0 0.3um 0.4um 0.6um\n"));
    REQUIRE(has(w->deck().text(), "box = 0.8um 1um 0 0.3um 0.4um 0.6um\n"));
    step(shown);
    REQUIRE(page.selection().solids.size() == 1);
    const solid::SolidRef b = page.selection().solids[0];
    CHECK(b.region == 3);
    // A copy pressed and released in place, or cancelled by Esc, is no copy.
    std::tie(px, at) = face_pixel(page, a, 5);
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, px[0], px[1]);
    send(c, WM_LBUTTONUP, 0, px[0], px[1]);
    CHECK(w->deck().text() == texts.back());
    CHECK(page.model().regions.size() == 4);
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, px[0], px[1]);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, px[0] + 30, px[1]);
    CHECK(page.model().regions.size() == 5);
    SendMessageW(c, WM_KEYDOWN, VK_ESCAPE, 0);
    send(c, WM_LBUTTONUP, 0, px[0] + 30, px[1]);
    CHECK(w->deck().text() == texts.back());
    CHECK(page.model().regions.size() == 4);
    page.render();
    // Mirror: a plain click on B's top mirrors about its centre (a box is symmetric: nothing to write); Ctrl+click,
    // across that face's own plane (z = 0.6um to 0.8um: out of the domain, refused, the deck kept).
    page.set_mode(modeller::Mode::mirror);
    std::tie(px, at) = face_pixel(page, b, 5);
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, px[0], px[1]);
    send(c, WM_LBUTTONUP, 0, px[0], px[1]);
    CHECK(w->deck().text() == texts.back());
    send(c, WM_LBUTTONDOWN, MK_LBUTTON | MK_CONTROL, px[0], px[1]);
    send(c, WM_LBUTTONUP, MK_CONTROL, px[0], px[1]);
    CHECK(w->deck().text() == texts.back());
    page.render();

    // Align: B lowered by the inspector (its top to 0.5um), then A's top clicked: B moved up, its top to 0.6um.
    REQUIRE(page.selection().solids == std::vector<solid::SolidRef>{b});
    win32::set_text(page.control(SolidPage::c_z1), "0.5");
    REQUIRE(page.commit(SolidPage::c_z1));
    REQUIRE(has(w->deck().text(), "box = 0.8um 1um 0 0.3um 0.4um 0.5um\n"));
    step("");
    page.render();
    page.set_mode(modeller::Mode::align);
    std::tie(px, at) = face_pixel(page, a, 5);
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, px[0], px[1]);
    send(c, WM_LBUTTONUP, 0, px[0], px[1]);
    REQUIRE(has(w->deck().text(), "box = 0.8um 1um 0 0.3um 0.5um 0.6um\n"));
    step("");
    page.render();

    // Array: a drag of 0.1um along y from A, the wheel up once: four of it.
    page.set_mode(modeller::Mode::array);
    std::tie(px, at) = face_pixel(page, a, 5);
    shown = drag_by(page, px, at, {0.0, 0.1, 0.0}, [&] {
        CHECK(has(page.chip(), "\xC3\x97 3"));
        wheel(c, px[0], px[1], 120);
        CHECK(page.array_count() == 4);
        CHECK(has(page.chip(), "\xC3\x97 4"));
    });
    const std::string arrayed = w->deck().text();
    CHECK(has(arrayed, "box = 0 0.2um 0.1um 0.4um 0.4um 0.6um\n"));
    CHECK(has(arrayed, "box = 0 0.2um 0.2um 0.5um 0.4um 0.6um\n"));
    CHECK(has(arrayed, "box = 0 0.2um 0.3um 0.6um 0.4um 0.6um\n"));
    CHECK(page.model().regions.size() == 7);
    step(shown);

    // Fill / Replace on A: its region's line; the inspector's switch shows it, and takes it back.
    page.set_mode(modeller::Mode::none);
    std::tie(px, at) = face_pixel(page, a, 5);
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, px[0], px[1]);
    send(c, WM_LBUTTONUP, 0, px[0], px[1]);
    REQUIRE(page.selection().solids == std::vector<solid::SolidRef>{a});
    REQUIRE(page.toggle_overlap());
    CHECK(has(w->deck().text(), "[region silicon_1]\nmaterial = silicon\noverlap = fill\nbox = 0 0.2um 0 0.3um 0.4um 0.6um\n"));
    CHECK(SendMessageW(page.control(SolidPage::c_overlap), BM_GETCHECK, 0, 0) == BST_CHECKED);
    step("");
    SendMessageW(page.control(SolidPage::c_overlap), BM_SETCHECK, BST_UNCHECKED, 0);
    REQUIRE(page.commit(SolidPage::c_overlap));
    CHECK(w->deck().text() == texts[texts.size() - 2]);
    step("");

    // Each one undo step: undone one by one, every earlier text comes back byte for byte; redone, the last.
    for (std::size_t k = texts.size() - 1; k-- > 0;) {
        REQUIRE(page.undo());
        CHECK(w->deck().text() == texts[k]);
    }
    for (std::size_t k = 1; k < texts.size(); ++k) REQUIRE(page.redo());
    CHECK(w->deck().text() == texts.back());
}

TEST_CASE("app window modeller: Copy, Mirror, Align, Array and Fill/Replace in 2D, each one undo step (MC4, MC6, MQ6)",
          "[app][window][modeller][transform]") {
    auto w = window();
    REQUIRE(w->open_deck(write("transform2d", layered2d)));
    pump();
    w->execute(Command::tab_geometry);
    GeometryPage& page = w->geometry_page();
    page.set_view(geo::View{200.0, 50.0, 450.0, 96});
    const HWND c = page.canvas();
    const auto at = [&](geo::Coord x, geo::Coord y) { return std::array<double, 2>{page.view().sx(x), page.view().sy(y)}; };
    const auto click = [&](std::array<double, 2> p, WPARAM keys = 0) {
        send(c, WM_LBUTTONDOWN, MK_LBUTTON | keys, p[0], p[1]);
        send(c, WM_LBUTTONUP, keys, p[0], p[1]);
    };
    const auto drag = [&](std::array<double, 2> a, std::array<double, 2> b, const std::function<void()>& during = {}) {
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, a[0], a[1]);
        send(c, WM_MOUSEMOVE, MK_LBUTTON, (a[0] + b[0]) / 2, (a[1] + b[1]) / 2);
        send(c, WM_MOUSEMOVE, MK_LBUTTON, b[0], b[1]);
        if (during) during();
        const std::string shown = page.live_line();
        send(c, WM_LBUTTONUP, 0, b[0], b[1]);
        return shown;
    };
    // A rectangle A: 0.4 0.8 0.2 0.5um.
    page.set_material("sio2");
    page.set_tool(GeometryPage::Tool::rectangle);
    drag(at(4000, 2000), at(8000, 5000));
    REQUIRE(has(w->deck().text(), "box = 0.4um 0.8um 0.2um 0.5um\n"));
    std::vector<std::string> texts{w->deck().text()};
    const auto step = [&](const std::string& shown) {
        const std::string now = w->deck().text();
        CHECK(now != texts.back());
        if (!shown.empty()) CHECK(shown == modeller::deck_line(texts.back(), now));  // MC4
        texts.push_back(now);
    };

    // Copy: the copy moved 0.6um right, A kept.
    page.set_mode(modeller::Mode::copy);
    CHECK(page.tool() == GeometryPage::Tool::select);
    std::string shown = drag(at(6000, 3500), at(12000, 3500));
    REQUIRE(has(w->deck().text(), "box = 0.4um 0.8um 0.2um 0.5um\n"));
    REQUIRE(has(w->deck().text(), "box = 1um 1.4um 0.2um 0.5um\n"));
    step(shown);
    REQUIRE(page.selection().size() == 1);
    const geo::ShapeRef b = page.selection()[0];
    CHECK(b.region == 3);
    click(at(6000, 3500));  // pressed and released in place: no copy
    CHECK(w->deck().text() == texts.back());
    CHECK(page.geometry().regions.size() == 4);

    // Mirror: Shift+click B, a quarter turn about its centre (1.2, 0.35): 0.3 wide, 0.4 high.
    page.set_mode(modeller::Mode::mirror);
    click(at(12000, 3500), MK_SHIFT);
    REQUIRE(has(w->deck().text(), "box = 1.05um 1.35um 0.15um 0.55um\n"));
    step("");

    // Align: B picked, a click near A's bottom side: B's bottom to y = 0.2um.
    page.set_mode(modeller::Mode::align);
    REQUIRE(page.selection() == std::vector<geo::ShapeRef>{b});
    click(at(6000, 2100));
    REQUIRE(has(w->deck().text(), "box = 1.05um 1.35um 0.2um 0.6um\n"));
    step("");

    // Array: from A, 0.1um up; the wheel down once from 3 leaves 2.
    page.set_mode(modeller::Mode::array);
    shown = drag(at(6000, 3500), at(6000, 4500), [&] {
        wheel(c, at(6000, 4500)[0], at(6000, 4500)[1], -120);
        CHECK(page.array_count() == 2);
        CHECK(has(page.chip(), "\xC3\x97 2"));
    });
    CHECK(has(w->deck().text(), "box = 0.4um 0.8um 0.3um 0.6um\n"));
    CHECK(page.geometry().regions.size() == 5);
    step(shown);

    // Fill / Replace on A, and back by the inspector's switch.
    page.set_mode(modeller::Mode::none);
    click(at(5000, 2500));
    REQUIRE(page.selection() == std::vector<geo::ShapeRef>{{2, 0}});
    REQUIRE(page.toggle_overlap());
    CHECK(has(w->deck().text(), "material = sio2\noverlap = fill\nbox = 0.4um 0.8um 0.2um 0.5um\n"));
    CHECK(SendMessageW(page.control(GeometryPage::c_overlap), BM_GETCHECK, 0, 0) == BST_CHECKED);
    step("");
    SendMessageW(page.control(GeometryPage::c_overlap), BM_SETCHECK, BST_UNCHECKED, 0);
    REQUIRE(page.commit(GeometryPage::c_overlap));
    CHECK(w->deck().text() == texts[texts.size() - 2]);
    step("");

    for (std::size_t k = texts.size() - 1; k-- > 0;) {
        REQUIRE(page.undo());
        CHECK(w->deck().text() == texts[k]);
    }
    for (std::size_t k = 1; k < texts.size(); ++k) REQUIRE(page.redo());
    CHECK(w->deck().text() == texts.back());
}

// --- Unit 26m slice 2 (MT6, MT7, MT12; MQ6): the Faces and Edges tools and Cylinder, by window messages ------------

namespace {

// A silicon block alone (no cap hiding the solids' faces), s = 0.01um (6 s = 600 quanta, 4 s = 400).
const std::string block = R"([mesh]
type = generate3d
layered = on
domain = 0 1um 0 0.6um 0 0.6um
patch body = 0 0 0  1um 0.6um 0
max_size = 0.15um
interface_resolution = 0.01um

[region si]
material = silicon

[doping]
profile = acceptor constant 1e16

[contact body]
kind = ohmic
patch = body

[equilibrium]
)";

// A pixel where the pick finds what `want` says, else the test stops.
// A pixel pressing face `face` of `s`: the face itself, or its handle when `s` is selected; and the point on the face (µm).
std::pair<std::array<double, 2>, solid::V3> face_or_handle(const SolidPage& page, solid::SolidRef s, int face) {
    const auto px = find_pick(page, [&](const solid::Pick& p) {
        return p.solid == s && p.index == face && (p.kind == solid::PickKind::face || (p.kind == solid::PickKind::handle && p.handle == solid::HandleKind::face));
    });
    REQUIRE(px);
    const auto hs = solid::handles(page.model(), page.view(), s);
    const auto h = std::ranges::find_if(hs, [&](const solid::Handle& k) { return k.kind == solid::HandleKind::face && k.index == face; });
    REQUIRE(h != hs.end());
    return {*px, h->pts[0]};
}

std::array<double, 2> must_find(const SolidPage& page, const std::function<bool(const solid::Pick&)>& want, const char* what) {
    const auto p = find_pick(page, want);
    INFO("looking for " << what);
    REQUIRE(p);
    return *p;
}

}  // namespace

TEST_CASE("app window modeller: Extrude, Chamfer, Fillet, Round vertex, Taper, Inset and Cylinder in 3D, each one undo step "
          "(MC4, MC6, MT6, MT7, MT12)",
          "[app][window][modeller][tools]") {
    auto [w, pp] = open_3d("tools3d", block);
    SolidPage& page = *pp;
    const HWND c = page.canvas();
    CHECK(page.resolution() == 100);
    page.set_material("sio2");
    const solid::SolidRef a = draw_box(*w, page);  // 0 0.2um 0 0.3um 0.4um 0.5um, region 1
    REQUIRE(a == solid::SolidRef{1, 0});
    page.set_material("silicon");  // what the tools add from here (of A's own material they would join its region)
    // The domain's own faces hide what lies inside it (26l): a section at z = 0.58um opens the top (slice 1's way).
    page.set_section({2, 5800});
    page.render();
    std::vector<std::string> texts{w->deck().text()};
    // Each step: the deck is exactly the model function's edit of the deck before (by minimal edits), the line shown
    // during the drag the line written (MC4).
    const auto step = [&](const solid::Model& expected, const std::string& shown) {
        const auto want = solid::write_model(texts.back(), expected);
        REQUIRE(want);
        const std::string now = w->deck().text();
        CHECK(now == *want);
        CHECK(now != texts.back());
        CHECK(shown == modeller::deck_line(texts.back(), now));  // every step here is a drag
        texts.push_back(now);
        page.render();
    };
    const auto press_drag = [&](std::array<double, 2> from, std::array<double, 2> to, const std::wstring& typed) {
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, from[0], from[1]);
        send(c, WM_MOUSEMOVE, MK_LBUTTON, (from[0] + to[0]) / 2, (from[1] + to[1]) / 2);
        send(c, WM_MOUSEMOVE, MK_LBUTTON, to[0], to[1]);
        type(c, typed);
        const std::string shown = page.live_line();
        send(c, WM_LBUTTONUP, 0, to[0], to[1]);
        return shown;
    };

    // Extrude A's top by 0.05um (typed): a new silicon region on it.
    page.set_mode(modeller::Mode::extrude);
    CHECK(page.group() == modeller::Group::faces);
    {
        solid::Model e = page.model();
        REQUIRE(solid::extrude_face(e, a, 5, 500, "silicon") == 2);
        auto [px, at] = face_or_handle(page, a, 5);
        // Dragged, before anything is typed: the depth is measured from the pressed point on the face (a handle's
        // press too: its own place on the face), so 0.1um up reads 0.1um (on the view's snap step).
        const auto to = page.screen_of_um({at[0], at[1], at[2] + 0.1});
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, px[0], px[1]);
        send(c, WM_MOUSEMOVE, MK_LBUTTON, (px[0] + to[0]) / 2, (px[1] + to[1]) / 2);
        send(c, WM_MOUSEMOVE, MK_LBUTTON, to[0], to[1]);
        INFO("chip: " << page.chip());
        CHECK(has(page.chip(), "depth 0.1 "));
        type(c, L"0.05");
        const std::string shown = page.live_line();
        send(c, WM_LBUTTONUP, 0, to[0], to[1]);
        step(e, shown);
    }
    const solid::SolidRef top{2, 0};
    // The largest interface_resolution (MT12): written by the inspector's button as one edit, then put back by undo
    // (the tools below keep the deck's s).
    {
        const auto s = page.automatic_resolution();
        INFO((s ? std::string() : s.error()));
        REQUIRE(s);
        REQUIRE(std::isfinite(*s));
        REQUIRE(page.use_automatic_resolution());
        const auto q = static_cast<solid::Coord>(std::floor(*s * 1e8));
        CHECK(has(w->deck().text(), "interface_resolution = " + geo::format_length(q) + "\n"));
        CHECK(page.resolution() == q);
        REQUIRE(page.undo());
        CHECK(w->deck().text() == texts.back());
        CHECK(page.resolution() == 100);
        page.render();
    }

    // Chamfer A's edge at x hi, y hi along z by 0.05um: A becomes a prism along z.
    page.set_mode(modeller::Mode::chamfer);
    CHECK(page.group() == modeller::Group::edges);
    auto [fx, fa] = face_or_handle(page, a, 1);  // a press on a face picks the solid (its edges then pickable) and moves nothing
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, fx[0], fx[1]);
    send(c, WM_LBUTTONUP, 0, fx[0], fx[1]);
    REQUIRE(page.selection().solids == std::vector<solid::SolidRef>{a});
    CHECK(w->deck().text() == texts.back());
    page.render();
    const auto edge_of = [&](int i, int j) {
        const auto list = page.model().regions[1].solids[0].edges();
        for (std::size_t k = 0; k < list.size(); ++k) {
            if (list[k] == std::pair{i, j}) return static_cast<int>(k);
        }
        return -1;
    };
    {
        const int edge = edge_of(3, 7);
        REQUIRE(edge >= 0);
        const auto px = must_find(page, [&](const solid::Pick& p) { return p.kind == solid::PickKind::edge && p.solid == a && p.index == edge; }, "the x hi, y hi edge");
        solid::Model e = page.model();
        REQUIRE(solid::chamfer_edge(e, a, edge, 500));
        const std::string shown = press_drag(px, {px[0] - 15, px[1] + 5}, L"0.05");
        step(e, shown);
        CHECK(page.model().regions[1].solids[0].kind == solid::SolidKind::prism);
    }
    // Fillet the x lo, y hi edge (the prism's vertex 4, edge 14) at 0.1um: 2 chords, as 6 s allows (MQ4, MT12).
    page.set_mode(modeller::Mode::fillet);
    page.render();
    {
        const auto px = must_find(page, [&](const solid::Pick& p) { return p.kind == solid::PickKind::edge && p.solid == a && p.index == 14; }, "edge 14");
        solid::Model e = page.model();
        REQUIRE(solid::fillet_edge(e, a, 14, 1000, 2));
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, px[0], px[1]);
        send(c, WM_MOUSEMOVE, MK_LBUTTON, px[0] + 10, px[1] + 10);
        type(c, L"0.1");
        CHECK(page.chords() == 2);
        CHECK(has(page.chip(), "2 chords"));
        const std::string shown = page.live_line();
        send(c, WM_LBUTTONUP, 0, px[0] + 10, px[1] + 10);
        step(e, shown);
    }
    // Round vertex: the prism's vertex 1 (x hi, y lo) at 0.05um from its handle.
    page.set_mode(modeller::Mode::round);
    page.render();
    {
        const auto px = must_find(page, [&](const solid::Pick& p) {
            return p.kind == solid::PickKind::handle && p.handle == solid::HandleKind::vertex && p.solid == a && p.index == 1;
        }, "vertex handle 1");
        solid::Model e = page.model();
        REQUIRE(solid::fillet_vertex(e, a, 1, 500, 1));
        const std::string shown = press_drag(px, {px[0] - 10, px[1]}, L"0.05");
        step(e, shown);
    }
    // Taper: vertex 0's edge 0.02um along +x; the wall angle shows.
    page.set_mode(modeller::Mode::taper);
    page.render();
    {
        const auto px = must_find(page, [&](const solid::Pick& p) { return p.kind == solid::PickKind::edge && p.solid == a && p.index == 2; }, "edge 2");
        solid::Model e = page.model();
        REQUIRE(solid::taper_edge(e, a, 2, {200, 0}));
        const auto to = page.screen_of_um({0.05, 0.0, 0.45});
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, px[0], px[1]);
        send(c, WM_MOUSEMOVE, MK_LBUTTON, to[0], to[1]);
        type(c, L"0.02");
        CHECK(has(page.chip(), "wall"));
        const std::string shown = page.live_line();
        send(c, WM_LBUTTONUP, 0, to[0], to[1]);
        step(e, shown);
    }
    // Inset on the extruded top: the margin 0.05um, released, then the depth 0.05um up, a click ending it. Of the top's
    // own material it joins its region.
    page.set_mode(modeller::Mode::inset);
    page.render();
    {
        auto [px, at] = face_or_handle(page, top, 5);
        solid::Model e = page.model();
        REQUIRE(solid::inset_face(e, top, 5, 500, 500, "silicon") == 2);
        const auto in = page.screen_of_um({at[0] + 0.03, at[1] + 0.03, at[2]});
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, px[0], px[1]);
        send(c, WM_MOUSEMOVE, MK_LBUTTON, in[0], in[1]);
        type(c, L"0.05");
        CHECK(has(page.chip(), "margin 0.05"));
        send(c, WM_LBUTTONUP, 0, in[0], in[1]);
        CHECK(w->deck().text() == texts.back());  // not yet written
        const auto up = page.screen_of_um({at[0], at[1], at[2] + 0.1});
        send(c, WM_MOUSEMOVE, 0, up[0], up[1]);
        // Moved, before typing: the depth from the pressed point on the face, 0.1um up.
        INFO("chip: " << page.chip());
        CHECK(has(page.chip(), "depth 0.1 "));
        type(c, L"0.05");
        const std::string shown = page.live_line();
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, up[0], up[1]);
        send(c, WM_LBUTTONUP, 0, up[0], up[1]);
        step(e, shown);
    }
    // A typed margin under 4 s is raised to 4 s (MT12): shown in the chip, then Esc leaves the deck as it was.
    {
        auto [px, at] = face_or_handle(page, top, 5);
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, px[0], px[1]);
        send(c, WM_MOUSEMOVE, MK_LBUTTON, px[0] + 5, px[1] + 5);
        type(c, L"0.01");
        CHECK(has(page.chip(), "margin 0.04"));
        SendMessageW(c, WM_KEYDOWN, VK_ESCAPE, 0);
        send(c, WM_LBUTTONUP, 0, px[0] + 5, px[1] + 5);
        CHECK(w->deck().text() == texts.back());
    }
    // Cylinder on the work plane z = 0.2um: centre (0.6, 0.3)um, radius 0.1um (10 sides: 16 would be under 6 s), height 0.1um.
    page.set_mode(modeller::Mode::none);
    page.set_work_plane({2, 2000});
    page.set_tool(SolidPage::Tool::cylinder);
    CHECK(page.group() == modeller::Group::solid);
    {
        solid::Model e = page.model();
        REQUIRE(solid::add_cylinder(e, 2, 2000, 3000, {6000, 3000}, 1000, 10, "silicon") == 3);
        const auto ctr = page.screen_of_um({0.6, 0.3, 0.2}), rim = page.screen_of_um({0.68, 0.3, 0.2});
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, ctr[0], ctr[1]);
        send(c, WM_MOUSEMOVE, MK_LBUTTON, rim[0], rim[1]);
        type(c, L"0.1");
        CHECK(has(page.chip(), "10 sides"));
        send(c, WM_LBUTTONUP, 0, rim[0], rim[1]);
        const auto up = page.screen_of_um({0.7, 0.3, 0.3});
        send(c, WM_MOUSEMOVE, 0, up[0], up[1]);
        type(c, L"0.1");
        const std::string shown = page.live_line();
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, up[0], up[1]);
        send(c, WM_LBUTTONUP, 0, up[0], up[1]);
        step(e, shown);
    }
    page.set_tool(SolidPage::Tool::select);

    // The inset left the top's face a ring: a holed piece meshes since Unit 64, but this layout's interfaces end on the
    // boundary along one line, so still no s passes; the button says so, the deck kept (MT12).
    {
        const auto s = page.automatic_resolution();
        REQUIRE_FALSE(s);
        CHECK(has(s.error(), "two or more interfaces end on the boundary along one line"));
        CHECK_FALSE(page.use_automatic_resolution());
        CHECK(w->deck().text() == texts.back());
    }

    for (std::size_t k = texts.size() - 1; k-- > 0;) {
        REQUIRE(page.undo());
        CHECK(w->deck().text() == texts[k]);
    }
    for (std::size_t k = 1; k < texts.size(); ++k) REQUIRE(page.redo());
    CHECK(w->deck().text() == texts.back());
}

TEST_CASE("app window modeller: Chamfer and Round vertex on 2D corners, each one undo step (MC4, MC6, MQ6)", "[app][window][modeller][tools]") {
    auto w = window();
    REQUIRE(w->open_deck(write("tools2d", layered2d)));
    pump();
    w->execute(Command::tab_geometry);
    GeometryPage& page = w->geometry_page();
    page.set_view(geo::View{200.0, 50.0, 450.0, 96});
    const HWND c = page.canvas();
    const auto at = [&](geo::Coord x, geo::Coord y) { return std::array<double, 2>{page.view().sx(x), page.view().sy(y)}; };
    page.set_material("sio2");
    page.set_tool(GeometryPage::Tool::rectangle);
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, at(4000, 2000)[0], at(4000, 2000)[1]);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, at(8000, 5000)[0], at(8000, 5000)[1]);
    send(c, WM_LBUTTONUP, 0, at(8000, 5000)[0], at(8000, 5000)[1]);
    REQUIRE(has(w->deck().text(), "box = 0.4um 0.8um 0.2um 0.5um\n"));
    const geo::ShapeRef a{2, 0};
    std::vector<std::string> texts{w->deck().text()};
    const auto corner = [&](std::array<double, 2> from, std::array<double, 2> to, const std::wstring& typed, const geo::Geometry& expected) {
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, from[0], from[1]);
        send(c, WM_MOUSEMOVE, MK_LBUTTON, to[0], to[1]);
        type(c, typed);
        const std::string shown = page.live_line();
        send(c, WM_LBUTTONUP, 0, to[0], to[1]);
        const auto want = geo::write_geometry(texts.back(), expected);
        REQUIRE(want);
        CHECK(w->deck().text() == *want);
        CHECK(shown == modeller::deck_line(texts.back(), w->deck().text()));
        texts.push_back(w->deck().text());
    };

    // The 3D-only tools are not on the 2D bar.
    page.set_mode(modeller::Mode::chamfer);
    CHECK(page.group() == modeller::Group::edges);
    CHECK(modeller::tools_for(modeller::Group::edges, false).size() == 2);
    CHECK(modeller::tools_for(modeller::Group::faces, false).size() == 2);
    // Chamfer: a press inside A picks it; then from its upper-right corner handle, 0.05um.
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, at(6000, 3500)[0], at(6000, 3500)[1]);
    send(c, WM_LBUTTONUP, 0, at(6000, 3500)[0], at(6000, 3500)[1]);
    REQUIRE(page.selection() == std::vector<geo::ShapeRef>{a});
    CHECK(w->deck().text() == texts.back());
    {
        geo::Geometry e = page.geometry();
        REQUIRE(geo::chamfer_shape(e, a, 2, 500));
        corner(at(8000, 5000), at(7600, 4600), L"0.05", e);
    }
    // Round vertex at the lower-left corner (vertex 0), 0.2um: 2 chords (3 would be under 6 s = 0.12um).
    page.set_mode(modeller::Mode::round);
    {
        geo::Geometry e = page.geometry();
        REQUIRE(geo::fillet_shape(e, a, 0, 2000, 2));
        corner(at(4000, 2000), at(4500, 2500), L"0.2", e);
    }
    for (std::size_t k = texts.size() - 1; k-- > 0;) {
        REQUIRE(page.undo());
        CHECK(w->deck().text() == texts[k]);
    }
    for (std::size_t k = 1; k < texts.size(); ++k) REQUIRE(page.redo());
    CHECK(w->deck().text() == texts.back());
}

TEST_CASE("app window modeller: slice 2's per-frame work adds at most 1 ms with 100 solids (MC9, MT7, MT12)", "[app][window][modeller][perf]") {
    std::string deck = "[mesh]\ntype = generate3d\nlayered = on\ndomain = 0 10um 0 10um 0 1um\nmax_size = 0.5um\n\n[region si]\nmaterial = silicon\n";
    for (int k = 0; k < 100; ++k) {
        deck += std::format("\n[region r{}]\nmaterial = {}\nbox = {}um {}.8um {}um {}.8um 0 0.{}um\n", k, k % 2 ? "sio2" : "silicon", k % 10, k % 10,
                            k / 10, k / 10, 3 + k % 7);
    }
    deck += "\n[equilibrium]\n";
    auto [w, pp] = open_3d("hundred2", deck);
    SolidPage& page = *pp;
    REQUIRE(page.model().regions.size() == 101);
    REQUIRE(page.select_region("r55"));
    const std::string text = w->deck().text();
    const solid::SolidRef ref = page.selection().solids[0];
    const solid::Solid moving = page.model().regions[static_cast<std::size_t>(ref.region)].solids[0];
    std::vector<double> copy, fillet;
    for (int k = 1; k <= 31; ++k) {
        // A Copy drag's frame: kept clear of all 100 solids, then its deck line.
        auto t0 = std::chrono::steady_clock::now();
        solid::Model m = page.model();
        const solid::Coord d = solid::keep_clear(m, {moving}, {ref}, 0, k * 130, 400);
        REQUIRE(solid::duplicate_solid(m, ref, {d, 0, 0}) >= 0);
        auto out = solid::write_model(text, m);
        REQUIRE(out);
        std::string line = modeller::deck_line(text, *out);
        copy.push_back(1e3 * std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        REQUIRE_FALSE(line.empty());
        // A Fillet drag's frame: the chord count at 6 s, the fillet, its deck line.
        t0 = std::chrono::steady_clock::now();
        m = page.model();
        const solid::Coord r = 300 + 10 * k;
        const int n = geo::arc_segments(std::numbers::pi / 2, static_cast<double>(r), 600.0);
        REQUIRE(solid::fillet_edge(m, ref, 2, r, std::max(n, 1)));
        out = solid::write_model(text, m);
        REQUIRE(out);
        line = modeller::deck_line(text, *out);
        fillet.push_back(1e3 * std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        REQUIRE_FALSE(line.empty());
    }
    std::ranges::sort(copy);
    std::ranges::sort(fillet);
    std::printf("MC9 slice 2, 100 solids: Copy frame (keep_clear + line) median %.3f ms (max %.3f); Fillet frame median %.3f ms (max %.3f)\n",
                copy[copy.size() / 2], copy.back(), fillet[fillet.size() / 2], fillet.back());
#ifdef NDEBUG
    CHECK(copy[copy.size() / 2] <= 1.0);
    CHECK(fillet[fillet.size() / 2] <= 1.0);
#endif
}

// --- Unit 26m slice 3 (MT8, MC7, MC4): parameters in the window -----------------------------------------------------

TEST_CASE("app window modeller: parameters in the inspector, a tied face dragged, N naming an extent (MC4, MC7, MT8)",
          "[app][window][modeller][parameters]") {
    const std::string tied = R"([parameters]
x0 = 30nm
Lg = 20nm

[mesh]
type = generate3d
layered = on
domain = 0 100nm 0 60nm 0 14nm
patch body = 0 0 0  100nm 60nm 0
max_size = 0.05um
interface_resolution = 2nm

[region si]
material = silicon

[region gate]
material = sio2
box = $x0 $x0 + $Lg 40nm 50nm 0 14nm

[region cap]
material = sio2
box = $x0 $x0 + $Lg 50nm 60nm 0 14nm

[doping]
profile = acceptor constant 1e16

[contact body]
kind = ohmic
patch = body

[equilibrium]
)";
    auto [w, pp] = open_3d("params3d", tied);
    SolidPage& page = *pp;
    const HWND c = page.canvas();
    const std::string original = w->deck().text();
    constexpr solid::Coord nm = 10;
    const auto model = [&] {
        const auto m = solid::read_model(w->deck().text());
        REQUIRE(m);
        return *m;
    };
    CHECK(page.parameters_text() == "x0 = 30nm\r\nLg = 20nm");

    // MC7: Lg changed in the inspector moves every tied face; one undo step.
    REQUIRE(page.apply_parameters("x0 = 30nm\r\nLg = 25nm\r\n"));
    CHECK(has(w->deck().text(), "\nLg = 25nm\n"));
    CHECK(has(w->deck().text(), "box = $x0 $x0 + $Lg 40nm 50nm 0 14nm\n"));
    CHECK(model().regions[1].solids[0].box.hi.x == 55 * nm);
    CHECK(model().regions[2].solids[0].box.hi.x == 55 * nm);
    CHECK_FALSE(page.apply_parameters("Lg = $Lg"));  // a cycle: nothing written
    CHECK(has(w->deck().text(), "\nLg = 25nm\n"));
    REQUIRE(page.undo());
    CHECK(w->deck().text() == original);

    // A tied face dragged (the gate's x_max, typed 5nm): its parameter moves, the cap's tied face follows; the line
    // shown is the first line written (MC4).
    REQUIRE(page.select_region("gate"));
    page.set_section({2, 100});  // open the top: the domain's faces hide what lies inside (26l)
    page.render();
    const solid::SolidRef gate{1, 0};
    auto f = face_handle(page, gate, 1);
    const auto step = page.screen_of_um({0.06, 0.045, 0.007}), at = page.screen_of_um({0.05, 0.045, 0.007});
    std::string before = w->deck().text();
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, f[0], f[1]);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, f[0] + 0.5 * (step[0] - at[0]), f[1] + 0.5 * (step[1] - at[1]));
    type(c, L"5nm");
    const std::string shown = page.live_line();
    SendMessageW(c, WM_KEYDOWN, VK_RETURN, 0);
    send(c, WM_LBUTTONUP, 0, f[0], f[1]);
    std::string after = w->deck().text();
    CHECK(shown == modeller::deck_line(before, after));
    CHECK(has(after, "\nLg = " + NiTCAD::app::geo::format_length(25 * nm) + "\n"));
    CHECK(has(after, "box = $x0 $x0 + $Lg 40nm 50nm 0 14nm\n"));
    CHECK(model().regions[2].solids[0].box.hi.x == 55 * nm);
    REQUIRE(page.undo());
    CHECK(w->deck().text() == original);

    // N during a face drag (the gate's z_max), nothing typed: a new parameter L1 of the gate's height, its top
    // written `0 + $L1`; the chip says so.
    page.render();
    f = face_handle(page, gate, 5);
    before = w->deck().text();
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, f[0], f[1]);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, f[0], f[1] + 1);
    type(c, L"n");
    CHECK(has(page.chip(), "named $L1 (new"));
    // N again: each existing length parameter in turn, then a new one again.
    type(c, L"n");
    CHECK(has(page.chip(), "named $x0 (N: next)"));
    type(c, L"n");
    CHECK(has(page.chip(), "named $Lg (N: next)"));
    type(c, L"n");
    CHECK(has(page.chip(), "named $L1 (new"));
    type(c, L"-0.002");  // down: the domain ends at 14nm
    const std::string named = page.live_line();
    SendMessageW(c, WM_KEYDOWN, VK_RETURN, 0);
    send(c, WM_LBUTTONUP, 0, f[0], f[1]);
    after = w->deck().text();
    CHECK(named == modeller::deck_line(before, after));
    CHECK(has(after, "\nLg = 20nm\nL1 = " + NiTCAD::app::geo::format_length(12 * nm) + "\n"));
    CHECK(has(after, "0 0 + $L1\n"));
    CHECK(model().regions[1].solids[0].box.hi.z == 12 * nm);
    REQUIRE(page.undo());
    CHECK(w->deck().text() == original);
}

TEST_CASE("app window modeller: Region profile, Implant, Junction, Max size and Refine box, each one line and one undo step "
          "(MC4, MC8, MT9, MT10)",
          "[app][window][modeller][doping]") {
    auto [w, pp] = open_3d("doping3d", block);
    SolidPage& page = *pp;
    const HWND c = page.canvas();
    page.set_material("sio2");
    const solid::SolidRef a = draw_box(*w, page);  // x 0..0.2um, y 0..0.3um, z 0.4..0.5um
    REQUIRE(a == solid::SolidRef{1, 0});
    page.set_section({2, 5800});
    page.render();
    std::vector<std::string> texts{w->deck().text()};
    const auto step = [&](const std::string& want, const std::string& shown) {
        const std::string now = w->deck().text();
        CHECK(now == want);
        CHECK(now != texts.back());
        if (!shown.empty()) CHECK(shown == modeller::deck_line(texts.back(), now));
        // One undo step, then redone.
        REQUIRE(page.undo());
        CHECK(w->deck().text() == texts.back());
        REQUIRE(page.redo());
        CHECK(w->deck().text() == now);
        texts.push_back(now);
        page.render();
    };
    using NiTCAD::device::Dopant;

    // Region profile: the wheel up a decade (1e19), D to acceptors, then a press on A dragged up two 1-2-5 steps
    // (60 px): 5e19.
    page.set_mode(modeller::Mode::profile);
    CHECK(page.group() == modeller::Group::doping);
    {
        auto [px, at] = face_or_handle(page, a, 5);
        wheel(c, px[0], px[1], 120);
        type(c, L"d");
        const auto want = solid::add_profile(texts.back(), solid::region_profile(page.model(), a, Dopant::acceptor, 5e19));
        REQUIRE(want);
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, px[0], px[1]);
        send(c, WM_MOUSEMOVE, MK_LBUTTON, px[0], px[1] - 30);
        send(c, WM_MOUSEMOVE, MK_LBUTTON, px[0], px[1] - 60);
        INFO("chip: " << page.chip());
        CHECK(has(page.chip(), "acceptor 5e+19"));
        const std::string shown = page.live_line();
        send(c, WM_LBUTTONUP, 0, px[0], px[1] - 60);
        step(*want, shown);
    }
    // Implant into A's top (z = 0.5um): a window dragged on the face (x 0.05..0.15um, y 0.05..0.25um), released; the
    // peak 30nm in (typed), a click; the lateral spread 10nm (typed), a click.
    page.set_mode(modeller::Mode::implant);
    {
        auto [px, at] = face_or_handle(page, a, 5);
        const auto w0 = page.screen_of_um({0.05, 0.05, 0.5}), w1 = page.screen_of_um({0.15, 0.25, 0.5});
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, w0[0], w0[1]);
        send(c, WM_MOUSEMOVE, MK_LBUTTON, w1[0], w1[1]);
        INFO("chip: " << page.chip());
        CHECK(has(page.chip(), "window 0.1"));
        send(c, WM_LBUTTONUP, 0, w1[0], w1[1]);
        const auto to = page.screen_of_um({0.1, 0.15, 0.45});
        send(c, WM_MOUSEMOVE, 0, to[0], to[1]);
        CHECK(has(page.chip(), "peak depth "));
        CHECK(has(page.chip(), "acceptor 5e+19"));
        type(c, L"0.03");
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, to[0], to[1]);
        send(c, WM_LBUTTONUP, 0, to[0], to[1]);
        CHECK(w->deck().text() == texts.back());  // not yet: the lateral spread
        send(c, WM_MOUSEMOVE, 0, to[0] + 5, to[1]);
        type(c, L"0.01");
        CHECK(has(page.chip(), "lateral 0.01"));
        const std::string shown = page.live_line();
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, to[0] + 5, to[1]);
        send(c, WM_LBUTTONUP, 0, to[0] + 5, to[1]);
        const solid::Box3 window{{500, 500, 5000}, {1500, 2500, 5000}};
        const auto want = solid::add_profile(texts.back(), solid::implant_profile({2, 1}, 5000, window, 300, 100, Dopant::acceptor, 5e19));
        REQUIRE(want);
        step(*want, shown);
    }
    // A press on A's edge in Implant moves nothing (the tool takes faces; a drag would move A in Select).
    {
        const auto px = must_find(page, [&](const solid::Pick& p) { return p.kind == solid::PickKind::edge && p.solid == a; }, "an edge of A");
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, px[0], px[1]);
        send(c, WM_MOUSEMOVE, MK_LBUTTON, px[0] + 30, px[1] + 20);
        send(c, WM_LBUTTONUP, 0, px[0] + 30, px[1] + 20);
        CHECK(w->deck().text() == texts.back());
        page.render();
    }
    // Junction from a face looking up: refused, the deck as it was.
    page.set_mode(modeller::Mode::junction);
    {
        auto [px, at] = face_or_handle(page, a, 5);
        const auto to = page.screen_of_um({at[0], at[1], at[2] - 0.05});
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, px[0], px[1]);
        send(c, WM_MOUSEMOVE, MK_LBUTTON, to[0], to[1]);
        CHECK(has(page.live_line(), "looking down an axis"));
        send(c, WM_LBUTTONUP, 0, to[0], to[1]);
        CHECK(w->deck().text() == texts.back());
    }
    // Max size: a drag, then 0.12um typed; the chip estimates the nodes.
    page.set_mode(modeller::Mode::max_size);
    CHECK(page.group() == modeller::Group::mesh);
    {
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, 300, 300);
        send(c, WM_MOUSEMOVE, MK_LBUTTON, 300, 200);
        INFO("chip: " << page.chip());
        CHECK(has(page.chip(), "max_size "));
        CHECK(has(page.chip(), "nodes"));
        type(c, L"0.12");
        const std::string shown = page.live_line();
        send(c, WM_LBUTTONUP, 0, 300, 200);
        const auto want = solid::set_mesh_length(texts.back(), "max_size", 1200);
        REQUIRE(want);
        step(*want, shown);
    }
    // Refine box: a rectangle on the work plane (z = 0.4um), released; its height 0.1um (typed), a click; its size
    // 0.02um (typed), a click.
    page.set_mode(modeller::Mode::refine);
    {
        const auto p0 = page.screen_of_um({0.5, 0.1, 0.4}), p1 = page.screen_of_um({0.9, 0.5, 0.4});
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, p0[0], p0[1]);
        send(c, WM_MOUSEMOVE, MK_LBUTTON, p1[0], p1[1]);
        send(c, WM_LBUTTONUP, 0, p1[0], p1[1]);
        send(c, WM_MOUSEMOVE, 0, p1[0], p1[1] - 10);
        type(c, L"0.1");
        INFO("chip: " << page.chip());
        CHECK(has(page.chip(), "height 0.1"));
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, p1[0], p1[1] - 10);
        send(c, WM_LBUTTONUP, 0, p1[0], p1[1] - 10);
        CHECK(w->deck().text() == texts.back());  // not yet: the size
        send(c, WM_MOUSEMOVE, 0, p1[0], p1[1] - 20);
        type(c, L"0.02");
        CHECK(has(page.chip(), "size 0.02"));
        const std::string shown = page.live_line();
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, p1[0], p1[1] - 20);
        send(c, WM_LBUTTONUP, 0, p1[0], p1[1] - 20);
        const auto boxes = solid::size_boxes(w->deck().text());
        REQUIRE(boxes.size() == 1);
        CHECK(boxes[0].size == 200);
        CHECK(boxes[0].box.lo.z == 4000);
        CHECK(boxes[0].box.hi.z == 5000);
        CHECK(boxes[0].box.hi.x > boxes[0].box.lo.x);
        const auto want = solid::add_size_box(texts.back(), boxes[0]);
        REQUIRE(want);
        step(*want, shown);
    }
    // The whole deck still builds: the doping and mesh lines are the build's.
    const auto parsed = deck::parse(w->deck().text());
    REQUIRE(parsed);
    CHECK(deck::build(*parsed).has_value());
    CHECK(solid::profiles(w->deck().text()).size() == 3);
}

namespace {

const std::string fin_start = R"([parameters]
Lg = 0.04um

[mesh]
type = generate3d
layered = on
domain = 0 0.12um 0 0.16um 0 0.1um
max_size = 0.03um
interface_resolution = 0.004um

[region ox]
material = sio2

[equilibrium]
)";

// A box drawn by the Box tool: the base's corners on the work plane at z, then the height typed.
void draw_box_at(SolidPage& page, solid::V3 lo, solid::V3 hi) {
    const HWND c = page.canvas();
    page.set_work_plane({2, static_cast<solid::Coord>(std::llround(lo[2] * solid::per_um))});
    page.set_tool(SolidPage::Tool::box);
    page.render();
    const auto a = page.screen_of_um(lo), b = page.screen_of_um({hi[0], hi[1], lo[2]});
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, a[0], a[1]);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, (a[0] + b[0]) / 2, (a[1] + b[1]) / 2);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, b[0], b[1]);
    send(c, WM_LBUTTONUP, 0, b[0], b[1]);
    type(c, std::format(L"{}", hi[2] - lo[2]));
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, b[0], b[1]);
    send(c, WM_LBUTTONUP, 0, b[0], b[1]);
    page.set_tool(SolidPage::Tool::select);
    page.render();
}

// A contact dragged on a face of the domain, from `a` to `b` (um, on the face).
void drag_contact(SolidPage& page, solid::V3 a, solid::V3 b) {
    const HWND c = page.canvas();
    page.set_tool(SolidPage::Tool::contact);
    page.render();
    const auto p = page.screen_of_um(a), q = page.screen_of_um(b);
    send(c, WM_LBUTTONDOWN, MK_LBUTTON, p[0], p[1]);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, (p[0] + q[0]) / 2, (p[1] + q[1]) / 2);
    send(c, WM_MOUSEMOVE, MK_LBUTTON, q[0], q[1]);
    send(c, WM_LBUTTONUP, 0, q[0], q[1]);
    page.set_tool(SolidPage::Tool::select);
    page.render();
}

}  // namespace

TEST_CASE("app window modeller: a FinFET drawn by window messages builds and runs its equilibrium (MC10)",
          "[app][window][modeller][doping][finfet]") {
    auto [w, pp] = open_3d("finfet3d", fin_start);
    SolidPage& page = *pp;
    const HWND c = page.canvas();
    const auto deck = [&] { return w->deck().text(); };
    // The first fin: silicon, x through the domain, y 0.02..0.06um, z 0.02..0.07um.
    page.set_material("silicon");
    draw_box_at(page, {0.0, 0.02, 0.02}, {0.12, 0.06, 0.07});
    INFO(deck());
    // The base as the pixels give it (on the view's snap step), the height as typed.
    REQUIRE(page.selection().solids.size() == 1);
    const solid::SolidRef fin = page.selection().solids.at(0);
    const solid::Box3 fb = page.model().regions[static_cast<std::size_t>(fin.region)].solids[static_cast<std::size_t>(fin.solid)].bounds();
    CHECK(fb.lo.x == 0);
    CHECK(fb.hi.x == 1200);
    CHECK(std::abs(fb.lo.y - 200) <= 10);
    CHECK(std::abs(fb.hi.y - 600) <= 10);
    CHECK(fb.lo.z == 200);
    CHECK(fb.hi.z == 700);
    page.set_section({2, 900});
    page.render();
    // Array: two fins at a pitch of 0.08um along y (typed; the wheel down from 3).
    page.set_mode(modeller::Mode::array);
    {
        auto [px, at] = face_pixel(page, fin, 5);
        drag_by(page, px, at, {0.0, 0.08, 0.0}, [&] {
            wheel(c, px[0], px[1], -120);
            type(c, L"0.08");
        });
    }
    REQUIRE(has(deck(), std::format("box = 0 0.12um {} {} 0.02um 0.07um\n", geo::format_length(fb.lo.y + 800), geo::format_length(fb.hi.y + 800))));
    // Fillet on each fin's top edge at y hi (vertices 6, 7: x lo and hi, y hi, z hi) at 0.02um.
    page.set_mode(modeller::Mode::fillet);
    std::vector<solid::SolidRef> fins;
    for (std::size_t r = 0; r < page.model().regions.size(); ++r) {
        for (std::size_t k = 0; k < page.model().regions[r].solids.size(); ++k) {
            if (page.model().regions[r].material == "silicon") fins.push_back({static_cast<int>(r), static_cast<int>(k)});
        }
    }
    REQUIRE(fins.size() == 2);
    for (const solid::SolidRef f : fins) {
        auto [fx, fa] = face_or_handle(page, f, 5);
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, fx[0], fx[1]);
        send(c, WM_LBUTTONUP, 0, fx[0], fx[1]);
        page.render();
        const auto list = page.model().regions[static_cast<std::size_t>(f.region)].solids[static_cast<std::size_t>(f.solid)].edges();
        int edge = -1;
        for (std::size_t k = 0; k < list.size(); ++k) {
            if (list[k] == std::pair{6, 7}) edge = static_cast<int>(k);
        }
        REQUIRE(edge >= 0);
        const auto px = must_find(page, [&](const solid::Pick& p) { return p.kind == solid::PickKind::edge && p.solid == f && p.index == edge; }, "the top y hi edge");
        const std::string before = deck();
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, px[0], px[1]);
        send(c, WM_MOUSEMOVE, MK_LBUTTON, px[0] + 10, px[1] + 10);
        type(c, L"0.02");
        INFO("chip: " << page.chip());
        send(c, WM_LBUTTONUP, 0, px[0] + 10, px[1] + 10);
        CHECK(deck() != before);
        page.render();
    }
    // The gate oxide: a sio2 box over the fins, x 0.04..0.07um, then its x_max face pushed out with N naming Lg
    // (N: a new L1, N again: Lg); the region fills only empty space (Fill), so the fins stay.
    page.set_mode(modeller::Mode::none);
    page.set_material("sio2");
    draw_box_at(page, {0.04, 0.0, 0.02}, {0.07, 0.16, 0.09});
    const solid::SolidRef gate = page.selection().solids.at(0);
    REQUIRE(page.toggle_overlap());
    {
        const auto hp = face_handle(page, gate, 1);
        const auto hs = solid::handles(page.model(), page.view(), gate);
        const auto h = std::ranges::find_if(hs, [&](const solid::Handle& k) { return k.kind == solid::HandleKind::face && k.index == 1; });
        REQUIRE(h != hs.end());
        const auto to = page.screen_of_um({h->pts[0][0] + 0.01, h->pts[0][1], h->pts[0][2]});
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, hp[0], hp[1]);
        send(c, WM_MOUSEMOVE, MK_LBUTTON, to[0], to[1]);
        type(c, L"nn");
        INFO("chip: " << page.chip());
        CHECK(has(page.chip(), "$Lg"));
        send(c, WM_LBUTTONUP, 0, to[0], to[1]);
    }
    CHECK(has(deck(), "um + $Lg 0 0.16um 0.02um 0.09um\n"));
    // Each fin's top face (looking up z): a box's face 5, a filleted prism's side.
    const auto top_of = [&](solid::SolidRef f) {
        const solid::Solid& s = page.model().regions[static_cast<std::size_t>(f.region)].solids[static_cast<std::size_t>(f.solid)];
        for (int k = 0; k < 2 + static_cast<int>(s.polygon.size()) + 6; ++k) {
            const auto n = solid::face_normal(s, k);
            if (n && n->axis == 2 && n->sign > 0 && (s.kind == solid::SolidKind::box ? k == 5 : k >= 2)) return k;
        }
        return -1;
    };
    // A Gaussian source/drain implant, the gate masked: donors (the tools' start) at 1e20, two windows on each fin's
    // top (x 0..0.04um and 0.08..0.12um, the face's width), the peak 10nm in, the lateral spread 2nm.
    page.set_mode(modeller::Mode::implant);
    wheel(c, 300, 300, 120);  // 1e19
    wheel(c, 300, 300, 120);  // 1e20
    for (const solid::SolidRef f : fins) {
        for (const auto [x0, x1] : {std::pair{0.0, 0.04}, std::pair{0.08, 0.12}}) {
            auto [px, at] = face_or_handle(page, f, top_of(f));
            const double y0 = at[1] - 0.009, y1 = at[1] + 0.009;  // inside the flat top (0.02um wide after the fillet)
            const auto p = page.screen_of_um({x0 + (x0 == 0.0 ? 0.002 : 0.0), y0, at[2]}),
                       q = page.screen_of_um({x1 - (x1 == 0.12 ? 0.002 : 0.0), y1, at[2]});
            REQUIRE(solid::pick(page.model(), page.view(), std::round(p[0]), std::round(p[1]), page.selection()).solid == f);
            const std::string before = deck();
            send(c, WM_LBUTTONDOWN, MK_LBUTTON, p[0], p[1]);
            send(c, WM_MOUSEMOVE, MK_LBUTTON, q[0], q[1]);
            send(c, WM_LBUTTONUP, 0, q[0], q[1]);
            type(c, L"0.01");
            send(c, WM_LBUTTONDOWN, MK_LBUTTON, q[0], q[1]);
            send(c, WM_LBUTTONUP, 0, q[0], q[1]);
            type(c, L"0.002");
            send(c, WM_LBUTTONDOWN, MK_LBUTTON, q[0], q[1]);
            send(c, WM_LBUTTONUP, 0, q[0], q[1]);
            CHECK(deck() != before);
            CHECK(has(deck(), "lateral=2e-07cm\n"));
        }
    }
    // The fins' background: acceptors 1e17 (Region profile).
    page.set_mode(modeller::Mode::profile);
    type(c, L"d");  // acceptors
    wheel(c, 300, 300, -120);
    wheel(c, 300, 300, -120);
    wheel(c, 300, 300, -120);  // 1e17
    for (const solid::SolidRef f : fins) {
        auto [px, at] = face_or_handle(page, f, top_of(f));
        send(c, WM_LBUTTONDOWN, MK_LBUTTON, px[0], px[1]);
        send(c, WM_LBUTTONUP, 0, px[0], px[1]);
    }
    CHECK(solid::profiles(deck()).size() == 6);
    // Contacts: source and drain on the x faces over the fins, the gate on the top over the oxide box.
    page.set_mode(modeller::Mode::none);
    page.set_section({});
    page.render();
    drag_contact(page, {0.12, 0.002, 0.002}, {0.12, 0.158, 0.098});  // the whole face: no patch edge by the fins
    drag_contact(page, {0.04, 0.002, 0.1}, {0.08, 0.158, 0.1});  // inside the face: its edges pick the sides
    REQUIRE(page.selection().contact == 1);  // the gate: an electrode, by the inspector's Kind
    SendMessageW(page.control(SolidPage::c_kind), CB_SELECTSTRING, static_cast<WPARAM>(-1), reinterpret_cast<LPARAM>(L"electrode"));
    REQUIRE(page.commit(SolidPage::c_kind));
    {  // turned half round (as an orbit would): the x_min face in view
        render::Camera cam = page.view().camera;
        cam.azimuth_deg += 180.0f;
        page.set_camera(cam);
    }
    page.render();
    drag_contact(page, {0.0, 0.002, 0.002}, {0.0, 0.158, 0.098});
    INFO(deck());
    REQUIRE(page.model().contacts.size() == 3);
    // The inspector's largest interface_resolution the layout's rules allow (MT12).
    {
        const auto s = page.automatic_resolution();
        INFO((s ? std::to_string(*s) : s.error()));
        REQUIRE(s);
        REQUIRE(page.use_automatic_resolution());
        std::printf("MC10: interface_resolution %s\n", geo::format_length(page.resolution()).c_str());
    }
    const auto parsed = deck::parse(deck());
    REQUIRE(parsed);
    const auto plan = deck::build(*parsed);
    INFO((plan ? std::string() : plan.error().message));
    REQUIRE(plan);
    const std::size_t nodes = plan->device.mesh().node_count();
    std::printf("MC10: %zu nodes\n", nodes);
    CHECK(nodes <= 20000);
    const auto t0 = std::chrono::steady_clock::now();
    const auto eq = solve::solve_equilibrium(plan->device);
    INFO((eq ? std::string() : eq.error().message));
    CHECK(eq.has_value());
    std::printf("MC10: equilibrium %.1f s\n", std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
}

TEST_CASE("app window modeller: slice 4's per-frame work adds at most 1 ms with 100 solids and 100 profiles (MC9, MT9, MT10)",
          "[app][window][modeller][perf]") {
    std::string deck = "[mesh]\ntype = generate3d\nlayered = on\ndomain = 0 10um 0 10um 0 1um\nmax_size = 0.5um\n\n[region si]\nmaterial = silicon\n";
    for (int k = 0; k < 100; ++k) {
        deck += std::format("\n[region r{}]\nmaterial = {}\nbox = {}um {}.8um {}um {}.8um 0 0.{}um\n", k, k % 2 ? "sio2" : "silicon", k % 10, k % 10,
                            k / 10, k / 10, 3 + k % 7);
    }
    deck += "\n[doping]\nprofile = acceptor constant 1e16\n";
    for (int k = 0; k < 99; ++k) {
        deck += std::format("profile = donor gaussian 1e19 axis=z center={}um width=0.02um x={}um:{}.5um y={}um:{}.5um\n", 0.1 + 0.001 * k, k % 10,
                            k % 10, k / 10, k / 10);
    }
    deck += "\n[equilibrium]\n";
    auto [w, pp] = open_3d("hundred4", deck);
    SolidPage& page = *pp;
    REQUIRE(page.model().regions.size() == 101);
    const std::string text = w->deck().text();
    const auto existing = solid::profiles(text);
    REQUIRE(existing.size() == 100);
    std::vector<solid::Volume> volumes;
    for (const auto& p : existing) volumes.push_back({solid::profile_extent(page.model().domain, p), solid::doping_color(true, p.peak_cm3)});
    std::vector<double> implant, junction, drawn;
    const solid::Box3 window{{0, 0, 0}, {4000, 4000, 0}};  // its centre in the first donor window
    for (int k = 1; k <= 31; ++k) {
        // An Implant drag's frame: the profile, its deck line.
        auto t0 = std::chrono::steady_clock::now();
        const auto ip = solid::implant_profile({2, -1}, 0, window, 100 + 10 * k, 0, NiTCAD::device::Dopant::donor, 1e19);
        auto out = solid::add_profile(text, ip);
        REQUIRE(out);
        std::string line = modeller::deck_line(text, *out);
        implant.push_back(1e3 * std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        REQUIRE_FALSE(line.empty());
        // A Junction drag's frame: the other doping at the depth (100 profiles), the erfc, its deck line.
        t0 = std::chrono::steady_clock::now();
        const auto jp = solid::junction_profile(existing, {2, -1}, 0, window, 500 + 10 * k, NiTCAD::device::Dopant::acceptor, 1e21);
        INFO((jp ? std::string() : jp.error()));
        REQUIRE(jp);
        out = solid::add_profile(text, *jp);
        REQUIRE(out);
        line = modeller::deck_line(text, *out);
        junction.push_back(1e3 * std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count());
        REQUIRE_FALSE(line.empty());
        // The 100 profiles' volumes in the scene: the scene with them less the scene without.
        t0 = std::chrono::steady_clock::now();
        const auto bare = solid::solid_scene(page.model(), page.view(), page.selection());
        const double without = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        t0 = std::chrono::steady_clock::now();
        solid::Overlay ov;
        ov.volumes = volumes;
        const auto full = solid::solid_scene(page.model(), page.view(), page.selection(), {}, ov);
        drawn.push_back(1e3 * (std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() - without));
        REQUIRE(full.triangles.size() > bare.triangles.size());
    }
    std::ranges::sort(implant);
    std::ranges::sort(junction);
    std::ranges::sort(drawn);
    std::printf("MC9 slice 4, 100 solids, 100 profiles: Implant frame median %.3f ms (max %.3f); Junction frame median %.3f ms (max %.3f); "
                "100 volumes drawn median +%.3f ms\n",
                implant[implant.size() / 2], implant.back(), junction[junction.size() / 2], junction.back(), drawn[drawn.size() / 2]);
#ifdef NDEBUG
    CHECK(implant[implant.size() / 2] <= 1.0);
    CHECK(junction[junction.size() / 2] <= 1.0);
    CHECK(drawn[drawn.size() / 2] <= 1.0);
#endif
}
