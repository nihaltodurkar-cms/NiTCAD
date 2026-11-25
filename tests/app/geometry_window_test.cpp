// The Geometry tab (Unit 26k slice 2; DECISIONS.md GD6-GD8, GE5, GE8): the canvas driven by its messages (no
// screen automation): a rectangle drawn, moved and resized, a polygon drawn and a vertex dragged, a contact drawn
// along the boundary, undo and redo, each landing in the deck's text; the drag's frame time with 200 shapes.
#include <catch2/catch_test_macros.hpp>

#include "main_window.hpp"
#include "fixtures.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/syntax.hpp"

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
        REQUIRE(page.can_start());
        // No button: a drag on the empty canvas draws the device, its size shown live; a click alone draws nothing.
        page.set_view(geo::View{200.0, 50.0, 450.0, 96});
        const Mouse m{page};
        m.click(um / 2, um / 2);
        CHECK(page.start_readout().empty());
        CHECK_FALSE(page.geometry().editable());
        m.down(0, 0);
        m.move(um, um / 2);
        CHECK(page.start_readout() == "1 × 0.5 µm");
        m.move(3 * um, um);
        CHECK(page.start_readout() == "3 × 1 µm");
        m.up(3 * um, um);
        CHECK(page.start_readout().empty());
        CHECK(page.geometry().editable());
        CHECK(page.geometry().mesh == geo::MeshKind::layered);
        CHECK(has(w->deck().text(), "type = generate2d\nlayered = on\ndomain = 0 3um 0 1um\nmax_size = 0.15um\n"));
        CHECK_FALSE(w->solid_shown());
    }
    {
        // Esc cancels a drag in progress.
        auto w = window();
        REQUIRE(w->open_deck(write("empty", "[device]\ntemperature = 300\n")));
        w->execute(Command::tab_geometry);
        GeometryPage& page = w->geometry_page();
        page.set_view(geo::View{200.0, 50.0, 450.0, 96});
        const Mouse m{page};
        const std::string before = w->deck().text();
        m.down(0, 0);
        m.move(um, um);
        m.key(VK_ESCAPE);
        CHECK(page.start_readout().empty());
        m.up(um, um);
        CHECK(w->deck().text() == before);
        CHECK(page.can_start());
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
    // Recorded on the CI runner (owner, 2026-11-25, as SE5): its Debug build drew 17.6 ms on WARP; this workstation is
    // the gate.
    if (!app_test::on_ci()) CHECK(median <= 16.0);
}

// --- slice 3: live validation, Run, the mesh preview, a drawn MOSFET (GD9, GE7, GE8) ------------------------------

namespace {

// Waits for the live check of the deck's text now (the page then knows whether it builds).
bool checked_now(MainWindow& w, GeometryPage& page) {
    return pump_until([&] { return !w.live_checker().pending() && page.check() != GeometryPage::Check::pending; }, 120.0);
}

bool run_enabled(MainWindow& w) {
    // The toolbar's Run button (the menus are filled as they open).
    HWND bar = GetDlgItem(w.hwnd(), MainWindow::id_toolbar);
    return SendMessageW(bar, TB_ISBUTTONENABLED, static_cast<WPARAM>(Command::run), 0) != 0;
}

void replace_text(MainWindow& w, const std::string& from, const std::string& to) {
    const std::string text = w.deck().text();
    const auto at = text.find(from);
    REQUIRE(at != std::string::npos);
    SendMessageW(w.deck().editor(), EM_SETSEL, at, at + from.size());
    const std::wstring wide(to.begin(), to.end());
    SendMessageW(w.deck().editor(), EM_REPLACESEL, TRUE, reinterpret_cast<LPARAM>(wide.c_str()));
}

}  // namespace

TEST_CASE("app window geometry: a layout that does not build is marked, Run is off, the mesh previewed (GD9, GE7)",
          "[app][window][geometry]") {
    auto w = window();
    REQUIRE(w->open_deck(write("validate", layered)));
    REQUIRE(pump_until([&] { return !w->runner().busy(); }));
    w->execute(Command::tab_geometry);
    GeometryPage& page = w->geometry_page();
    page.set_view(geo::View{200.0, 50.0, 450.0, 96});
    const Mouse m{page};
    REQUIRE(checked_now(*w, page));
    CHECK(page.check() == GeometryPage::Check::valid);
    CHECK(page.problems().empty());
    CHECK(page.run_blocked().empty());
    CHECK(page.banner().empty());
    CHECK(run_enabled(*w));

    // The mesh preview: the mesh the check built, drawn over the layout.
    page.set_mesh_preview(true);
    CHECK(page.mesh_shown());
    CHECK(page.mesh().edges.size() > 50);
    CHECK(page.render() >= 0.0);

    // A contact drawn over the source's patch: the build's error, on both contacts; Run off with the reason.
    const std::size_t errors = w->console_log().count(app::LogLevel::error);
    m.key('C');
    m.drag(0, 3000, 0, 9000);
    REQUIRE(has(w->deck().text(), "patch contact_1 = 0 0.3um, 0 0.9um\n"));
    CHECK_FALSE(page.mesh_shown());  // the mesh is of the text before
    REQUIRE(checked_now(*w, page));
    REQUIRE(page.check() == GeometryPage::Check::invalid);
    REQUIRE(page.problems().size() == 1);
    const geo::Problem& p = page.problems().front();
    const auto check = app::check_deck(w->deck().text(), w->document().path());
    REQUIRE(check.diagnostics.size() == 1);
    INFO(p.message);
    CHECK(p.message == check.diagnostics.front().message);  // the build's own message
    CHECK(p.message.find("patch") != std::string::npos);  // the build's message (the new patch splits the source's edge)
    CHECK(std::ranges::count(p.markers, geo::Marker::Kind::contact, &geo::Marker::kind) == 2);
    CHECK(has(page.banner(), p.message));
    CHECK(has(page.run_blocked(), p.message));
    CHECK_FALSE(run_enabled(*w));
    CHECK(w->console_log().count(app::LogLevel::error) > errors);  // in the console's Messages
    CHECK(page.render() >= 0.0);
    // Run asked for anyway (F5 reaches execute only when enabled; a direct call too): nothing starts.
    w->execute(Command::run);
    CHECK(w->runner().state() == app::RunState::idle);
    CHECK(w->tab() == MainWindow::tab_geometry);

    // On the deck's tab Run is as before; back on the Geometry tab it is off again.
    w->select_tab(MainWindow::tab_deck);
    CHECK(run_enabled(*w));
    w->execute(Command::tab_geometry);
    CHECK_FALSE(run_enabled(*w));

    // Undone: valid once checked, Run on, the mesh shown again.
    w->execute(Command::undo);
    REQUIRE(checked_now(*w, page));
    CHECK(page.check() == GeometryPage::Check::valid);
    CHECK(page.problems().empty());
    CHECK(run_enabled(*w));
    CHECK(page.mesh_shown());
}

TEST_CASE("app window geometry: a MOSFET drawn on a new layout runs its sweep (GE8)", "[app][window][geometry]") {
    auto w = window();
    REQUIRE(w->open_deck(write("mosfet", "# a MOSFET drawn on the Geometry tab\n")));
    REQUIRE(pump_until([&] { return !w->runner().busy(); }));
    w->execute(Command::tab_geometry);
    GeometryPage& page = w->geometry_page();
    REQUIRE_FALSE(page.geometry().editable());
    page.set_view(geo::View{200.0, 50.0, 450.0, 96});
    Mouse{page}.drag(0, 0, 2 * um, um);  // the device drawn on the empty canvas: 2 × 1 µm
    REQUIRE(page.geometry().mesh == geo::MeshKind::layered);
    page.set_view(geo::View{200.0, 50.0, 450.0, 96});
    const Mouse m{page};
    const auto deck = [&] { return w->deck().text(); };
    const auto name = [&](const wchar_t* n) {
        SetWindowTextW(page.control(GeometryPage::c_name), n);
        page.commit(GeometryPage::c_name);
    };

    // The substrate is the template's background; the oxide drawn over its top 50 nm and named.
    page.set_material("sio2");
    m.key('R');
    m.drag(0, 9500, 20000, 10000);
    INFO(deck());
    REQUIRE(has(deck(), "material = sio2\nbox = 0 2um 0.95um 1um\n"));
    name(L"oxide");
    REQUIRE(has(deck(), "[region oxide]\n"));

    // The contacts along the boundary: source and drain on the sides below the oxide, the gate over it, the body
    // along the bottom.
    m.key('C');
    m.drag(0, 6000, 0, 9000);
    name(L"source");
    m.drag(20000, 6000, 20000, 9000);
    name(L"drain");
    m.drag(5000, 10000, 15000, 10000);
    name(L"gate");
    SendMessageW(page.control(GeometryPage::c_kind), CB_SELECTSTRING, static_cast<WPARAM>(-1), reinterpret_cast<LPARAM>(L"electrode"));
    page.commit(GeometryPage::c_kind);
    m.drag(0, 0, 20000, 0);
    name(L"body");
    REQUIRE(page.geometry().contacts.size() == 4);
    REQUIRE(has(deck(), "[contact gate]\nkind = electrode\n"));

    // What the canvas does not draw (a later unit): the gate's metal, the pockets' doping, the sweep, typed.
    replace_text(*w, "[contact gate]\nkind = electrode\n", "[contact gate]\nkind = electrode\nelectrode = metal\nwork_function = 4.1\n");
    replace_text(*w, "profile = acceptor constant 1e16\n",
                 "profile = acceptor constant 1e16\n"
                 "profile = donor gaussian 2e19 axis=y center=0.95um width=0.15um x=0:0.5um lateral=20nm\n"
                 "profile = donor gaussian 2e19 axis=y center=0.95um width=0.15um x=1.5um:2um lateral=20nm\n");
    replace_text(*w, "[equilibrium]\n", "[sweep]\ncontact gate = 0 : 1 : 0.25\ncontact drain = 0.1\n");
    w->execute(Command::tab_geometry);
    REQUIRE(checked_now(*w, page));
    for (const auto& p : page.problems()) UNSCOPED_INFO(p.message);
    REQUIRE(page.check() == GeometryPage::Check::valid);
    REQUIRE(run_enabled(*w));

    w->execute(Command::run);
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace() && w->workspace()->has_results(); }, 600.0));
    CHECK(w->runner().state() == app::RunState::done);
    CHECK_FALSE(w->preview());
}

TEST_CASE("app window geometry: a slanted gate drawn, a vertex moved alone, Shift keeps an edge straight (slice 4, GD10)",
          "[app][window][geometry]") {
    auto w = window();
    REQUIRE(w->open_deck(write("slanted", layered)));
    pump();
    w->execute(Command::tab_geometry);
    GeometryPage& page = w->geometry_page();
    page.set_view(geo::View{200.0, 50.0, 450.0, 96});
    const Mouse m{page};
    const auto deck = [&] { return w->deck().text(); };
    // A trapezoid gate: its sidewalls slant (no Shift, the clicks are taken as they are).
    page.set_material("sio2");
    m.key('P');
    for (const auto& [x, y] : std::vector<std::pair<geo::Coord, geo::Coord>>{{6000, 8000}, {14000, 8000}, {13000, 10000}, {7000, 10000}}) {
        m.click(x, y);
    }
    m.key(VK_RETURN);
    INFO(deck());
    REQUIRE(has(deck(), "polygon = 0.6um 0.8um, 1.4um 0.8um, 1.3um 1um, 0.7um 1um\n"));
    // With Shift a vertex is kept vertical from the last: (1.66, 0.6) comes in as (1.6, 0.6), not as (1.7, 0.6).
    m.key(VK_ESCAPE);
    m.key('P');
    m.click(14000, 2000);
    m.click(16000, 2000);
    m.down(16600, 6000, MK_SHIFT);
    m.up(16600, 6000);
    m.click(14000, 6000);
    m.key(VK_RETURN);
    REQUIRE(has(deck(), "polygon = 1.4um 0.2um, 1.6um 0.2um, 1.6um 0.6um, 1.4um 0.6um\n"));
    // The well's corner (1, 0.8) Ctrl-dragged to (1.1, 0.8): it moves alone, its left edge slanting.
    m.key('V');
    m.click(17000, 4000);
    REQUIRE(page.selection().size() == 1);
    m.down(10000, 8000, MK_CONTROL);
    m.move(10500, 8000, MK_LBUTTON | MK_CONTROL);
    m.move(11000, 8000, MK_LBUTTON | MK_CONTROL);
    m.up(11000, 8000);
    CHECK(has(deck(), "polygon = 1um 0, 2um 0, 2um 0.5um, 2um 0.8um, 1.1um 0.8um\n"));
    // The layout builds: Run is on.
    const auto d = NiTCAD::deck::parse(deck(), "s.deck");
    REQUIRE(d);
    const auto b = NiTCAD::deck::build(*d);
    if (!b) FAIL(b.error().message);
}
