// The Geometry tab (Unit 26k slice 2; DECISIONS.md GD6-GD8, GE5, GE8): the canvas driven by its messages (no
// screen automation): a rectangle drawn, moved and resized, a polygon drawn and a vertex dragged, a contact drawn
// along the boundary, undo and redo, each landing in the deck's text; the drag's frame time with 200 shapes.
#include <catch2/catch_test_macros.hpp>

#include "main_window.hpp"
#include "fixtures.hpp"

#include <windows.h>

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

using namespace NiTCAD;
using app::Command;
using app::GeometryPage;
using app::MainWindow;
namespace geo = app::geo;

namespace {

bool pump_until(const std::function<bool()>& done, double seconds = 60.0) {
    const auto end = std::chrono::steady_clock::now() + std::chrono::duration<double>(seconds);
    MSG m;
    while (!done()) {
        if (std::chrono::steady_clock::now() > end) return false;
        if (PeekMessageW(&m, nullptr, 0, 0, PM_REMOVE)) {
            TranslateMessage(&m);
            DispatchMessageW(&m);
        } else {
            MsgWaitForMultipleObjects(0, nullptr, FALSE, 10, QS_ALLINPUT);
        }
    }
    return true;
}

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
    const auto dir = app_test::scratch("geometry_window_" + name);
    const auto p = dir / (name + ".deck");
    std::ofstream(p, std::ios::binary) << text;
    return p;
}

const std::string layered = R"([mesh]
type = generate2d
layered = on
domain = 0 2um 0 1um
patch source = 0 0, 0 0.6um
max_size = 0.25um
interface_resolution = 20nm

[region p]
material = silicon

[region n]
material = silicon
polygon = 1um 0, 2um 0, 2um 0.5um, 2um 0.8um, 1um 0.8um

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

constexpr geo::Coord um = geo::per_um;

// The canvas's pixel of a model point (µm × 10⁴) in the test's view.
struct Mouse {
    GeometryPage& page;
    [[nodiscard]] LPARAM at(geo::Coord x, geo::Coord y) const {
        const int px = static_cast<int>(std::lround(page.view().sx(x))), py = static_cast<int>(std::lround(page.view().sy(y)));
        return MAKELPARAM(px, py);
    }
    void down(geo::Coord x, geo::Coord y, WPARAM keys = 0) const { SendMessageW(page.canvas(), WM_LBUTTONDOWN, MK_LBUTTON | keys, at(x, y)); }
    void move(geo::Coord x, geo::Coord y, WPARAM keys = MK_LBUTTON) const { SendMessageW(page.canvas(), WM_MOUSEMOVE, keys, at(x, y)); }
    void up(geo::Coord x, geo::Coord y) const { SendMessageW(page.canvas(), WM_LBUTTONUP, 0, at(x, y)); }
    void drag(geo::Coord x0, geo::Coord y0, geo::Coord x1, geo::Coord y1) const {
        down(x0, y0);
        move((x0 + x1) / 2, (y0 + y1) / 2);
        move(x1, y1);
        up(x1, y1);
    }
    void click(geo::Coord x, geo::Coord y) const {
        down(x, y);
        up(x, y);
    }
    void key(UINT k) const { SendMessageW(page.canvas(), WM_KEYDOWN, k, 0); }
};

bool has(const std::string& text, const std::string& what) { return text.find(what) != std::string::npos; }

}  // namespace

TEST_CASE("app window geometry: draw, move and resize a rectangle; a polygon and its vertex; a contact; undo, redo (GE8)",
          "[app][window][geometry]") {
    auto w = window();
    REQUIRE(w->open_deck(write("draw", layered)));
    pump();
    w->execute(Command::tab_geometry);
    REQUIRE(w->tab() == MainWindow::tab_geometry);
    GeometryPage& page = w->geometry_page();
    REQUIRE(page.geometry().editable());
    REQUIRE(page.geometry().mesh == geo::MeshKind::layered);
    // 200 px/µm: the grid step is 0.1 µm, and every tenth of a µm is a whole pixel.
    page.set_view(geo::View{200.0, 50.0, 450.0, 96});
    REQUIRE(geo::grid_step(page.view()) == 1000);
    const Mouse m{page};
    const auto deck = [&] { return w->deck().text(); };

    // A rectangle (R), drawn with the material chosen for drawing.
    page.set_material("sio2");
    m.key('R');
    CHECK(page.tool() == GeometryPage::Tool::rectangle);
    m.drag(4000, 2000, 8000, 5000);
    INFO(deck());
    REQUIRE(has(deck(), "[region sio2_1]\nmaterial = sio2\nbox = 0.4um 0.8um 0.2um 0.5um\n"));
    REQUIRE(page.selection().size() == 1);
    CHECK(page.selection()[0].region == 3);
    CHECK(w->document().modified());

    // Moved by 0.1 µm (V, then a drag from inside it), then resized by its upper right handle.
    m.key('V');
    m.drag(6000, 3500, 7000, 3500);
    REQUIRE(has(deck(), "box = 0.5um 0.9um 0.2um 0.5um\n"));
    m.drag(9000, 5000, 10000, 6000);
    REQUIRE(has(deck(), "box = 0.5um 1um 0.2um 0.6um\n"));

    // Undo and redo step through the deck's texts.
    w->execute(Command::undo);
    CHECK(has(deck(), "box = 0.5um 0.9um 0.2um 0.5um\n"));
    w->execute(Command::undo);
    CHECK(has(deck(), "box = 0.4um 0.8um 0.2um 0.5um\n"));
    w->execute(Command::redo);
    w->execute(Command::redo);
    CHECK(has(deck(), "box = 0.5um 1um 0.2um 0.6um\n"));
    CHECK_FALSE(page.can_redo());
    // Down 0.1 µm, and back by undo.
    m.drag(7000, 4000, 7000, 3000);
    CHECK(has(deck(), "box = 0.5um 1um 0.1um 0.5um\n"));
    w->execute(Command::undo);
    CHECK(has(deck(), "box = 0.5um 1um 0.2um 0.6um\n"));

    // A polygon (P): its vertices clicked, Enter closes it.
    m.key('P');
    for (const auto& [x, y] : std::vector<std::pair<geo::Coord, geo::Coord>>{
             {12000, 1000}, {18000, 1000}, {18000, 3000}, {15000, 3000}, {15000, 6000}, {12000, 6000}}) {
        m.click(x, y);
    }
    m.key(VK_RETURN);
    REQUIRE(has(deck(), "polygon = 1.2um 0.1um, 1.8um 0.1um, 1.8um 0.3um, 1.5um 0.3um, 1.5um 0.6um, 1.2um 0.6um\n"));
    // Its vertex dragged up 0.1 µm: the neighbours follow along their edges.
    m.key('V');
    m.click(13000, 2000);
    REQUIRE(page.selection().size() == 1);
    m.drag(18000, 3000, 18000, 4000);
    CHECK(has(deck(), "polygon = 1.2um 0.1um, 1.8um 0.1um, 1.8um 0.4um, 1.5um 0.4um, 1.5um 0.6um, 1.2um 0.6um\n"));

    // A contact (C) dragged along the bottom of the domain (on the silicon: an ohmic contact on the oxide is the
    // build's error).
    m.key('C');
    m.drag(2000, 0, 6000, 0);
    CHECK(has(deck(), "patch contact_1 = 0.2um 0, 0.6um 0\n"));
    CHECK(has(deck(), "[contact contact_1]\nkind = ohmic\npatch = contact_1\n"));
    CHECK(page.selected_contact() == 1);

    // Delete removes the selection; Ctrl+D duplicates; the rest of the deck is as it was.
    m.key('V');
    m.click(13000, 2000);
    w->execute(Command::duplicate);
    CHECK(page.geometry().regions.size() == 6);
    w->execute(Command::delete_selection);
    CHECK(page.geometry().regions.size() == 5);
    CHECK(has(deck(), "[doping]\nprofile = acceptor constant 1e16\n"));
    CHECK(has(deck(), "patch source = 0 0, 0 0.6um\n"));

    // The canvas draws (its rulers in µm).
    CHECK(page.render() >= 0.0);
    CHECK_FALSE(page.ruler_text(true).empty());
    // The deck, as drawn, builds: the live check finds no error.
    const auto submitted = w->live_checker().submitted();
    REQUIRE(pump_until([&] { return w->live_checker().submitted() > submitted && !w->live_checker().pending(); }, 60.0));
    for (const auto& d : w->deck().diagnostics()) UNSCOPED_INFO("line " << d.line << ": " << d.message);
    CHECK(w->deck().diagnostics().empty());
}

TEST_CASE("app window geometry: the inspector's Geometry section edits the selection; the tree selects on the canvas (GD8)",
          "[app][window][geometry]") {
    auto w = window();
    REQUIRE(w->open_deck(write("inspector", layered)));
    pump();
    w->execute(Command::tab_geometry);
    GeometryPage& page = w->geometry_page();
    REQUIRE(page.select_region("oxide"));
    CHECK((GetWindowLongPtrW(page.inspector(), GWL_STYLE) & WS_VISIBLE) != 0);
    // Typed values, taken as Enter does.
    SetWindowTextW(page.control(GeometryPage::c_h), L"0.15");
    page.commit(GeometryPage::c_h);
    CHECK(has(w->deck().text(), "[region oxide]\nmaterial = sio2\nbox = 0 2um 0.8um 0.95um\n"));
    SetWindowTextW(page.control(GeometryPage::c_name), L"gate_oxide");
    page.commit(GeometryPage::c_name);
    CHECK(has(w->deck().text(), "[region gate_oxide]\n"));
    // A bad value is refused and the deck kept.
    const std::string before = w->deck().text();
    SetWindowTextW(page.control(GeometryPage::c_w), L"-1");
    page.commit(GeometryPage::c_w);
    CHECK(w->deck().text() == before);
    // The background: the domain's size.
    REQUIRE(page.select_region("p"));
    SetWindowTextW(page.control(GeometryPage::c_w), L"2.5");
    page.commit(GeometryPage::c_w);
    CHECK(has(w->deck().text(), "domain = 0 2.5um 0 1um\n"));
    // A polygon's vertex list.
    REQUIRE(page.select_region("n"));
    SetWindowTextW(page.control(GeometryPage::c_vertices), L"1 0\r\n2 0\r\n2 0.7\r\n1 0.7\r\n");
    page.commit(GeometryPage::c_apply);
    CHECK(has(w->deck().text(), "polygon = 1um 0, 2um 0, 2um 0.7um, 1um 0.7um\n"));
}

TEST_CASE("app window geometry: decks it does not draw say why; a deck without a mesh starts a 2D layout (GD6)",
          "[app][window][geometry]") {
    {
        auto w = window();
        REQUIRE(w->open_deck(std::filesystem::path(NITCAD_EXAMPLES_DIR) / "mosfet_3d.deck"));
        w->execute(Command::tab_geometry);
        CHECK(w->geometry_page().message() == "3D geometry editing is a later unit");
        CHECK(w->geometry_page().render() >= 0.0);
    }
    {
        auto w = window();
        REQUIRE(w->open_deck(write("empty", "[device]\ntemperature = 300\n")));
        w->execute(Command::tab_geometry);
        GeometryPage& page = w->geometry_page();
        CHECK_FALSE(page.geometry().editable());
        CHECK((GetWindowLongPtrW(page.control(GeometryPage::c_new_layout), GWL_STYLE) & WS_VISIBLE) != 0);
        page.commit(GeometryPage::c_new_layout);
        CHECK(page.geometry().editable());
        CHECK(page.geometry().mesh == geo::MeshKind::layered);
        CHECK(has(w->deck().text(), "layered = on\n"));
    }
    {
        // The tensor example: boxes only; its domain the built mesh's.
        auto w = window();
        REQUIRE(w->open_deck(std::filesystem::path(NITCAD_EXAMPLES_DIR) / "mosfet_2d.deck"));
        pump_until([&] { return w->workspace() != nullptr; }, 60.0);
        w->execute(Command::tab_geometry);
        GeometryPage& page = w->geometry_page();
        REQUIRE(page.geometry().mesh == geo::MeshKind::tensor);
        CHECK(page.geometry().domain == geo::Rect{{0, 0}, {2 * um, um}});
        page.set_tool(GeometryPage::Tool::polygon);
        CHECK(page.tool() == GeometryPage::Tool::polygon);
    }
}

TEST_CASE("app window geometry: a drag redraws within 16 ms a frame with 200 shapes and 5000 vertices (GE5)",
          "[app][window][geometry]") {
    std::string deck = "[mesh]\ntype = generate2d\nlayered = on\ndomain = 0 20um 0 10um\nmax_size = 1um\ninterface_resolution = 0.1um\n\n"
                       "[region bg]\nmaterial = silicon\n";
    for (int k = 0; k < 199; ++k) {
        deck += std::format("\n[region r{}]\nmaterial = sio2\nbox = {}um {}um {}um {}um\n", k, k % 20, k % 20 + 0.5, k / 20,
                            k / 20 + 0.5);
    }
    // A staircase of 4802 vertices.
    deck += "\n[region stair]\nmaterial = silicon\npolygon = 0 0";
    for (int k = 0; k < 2400; ++k) deck += std::format(", {}nm {}nm, {}nm {}nm", (k + 1) * 0.5, k * 0.5, (k + 1) * 0.5, (k + 1) * 0.5);
    deck += ", 0 1200nm\n\n[doping]\nprofile = acceptor constant 1e16\n\n[equilibrium]\n";
    auto w = window();
    REQUIRE(w->open_deck(write("speed", deck)));
    w->execute(Command::tab_geometry);
    GeometryPage& page = w->geometry_page();
    REQUIRE(page.geometry().regions.size() == 201);
    std::size_t vertices = 0;
    for (const auto& r : page.geometry().regions) {
        for (const auto& s : r.shapes) vertices += s.kind == geo::ShapeKind::box ? 4 : s.points.size();
    }
    CHECK(vertices >= 5000);
    page.fit();
    const Mouse m{page};
    // r10's box (10..10.5 µm, 0..0.5 µm) dragged 60 frames.
    m.down(102500, 2500);
    REQUIRE(page.selection().size() == 1);
    std::vector<double> ms;
    for (int k = 1; k <= 60; ++k) {
        m.move(102500 + k * 500, 2500 + k * 200);
        ms.push_back(page.render());
    }
    m.up(102500 + 60 * 500, 2500 + 60 * 200);
    std::ranges::sort(ms);
    const double median = ms[ms.size() / 2];
    WARN("GE5: a drag frame, median " << median << " ms (slowest " << ms.back() << " ms), " << vertices << " vertices, software Direct2D");
    CHECK(median <= 16.0);
}
