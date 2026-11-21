// Unit 26l slice 3 (SD5): the 3D modeller's live validation and mesh preview, through the main window's live check.
// SE7: for a table of invalid layouts the page's problem is the deck build's own message, its markers where expected,
// Run off with the reason; every valid layout builds. SE8: a 3D MOSFET drawn on the page (substrate, oxide, gate,
// source, drain, body; coarse) runs its equilibrium and one gate bias.
#include <catch2/catch_test_macros.hpp>

#include "main_window.hpp"
#include "fixtures.hpp"

#include <windows.h>
#include <commctrl.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

#include "NiTCAD/app/document.hpp"
#include "NiTCAD/app/solid.hpp"
#include "NiTCAD/app/solid_view.hpp"

using namespace NiTCAD;
using namespace NiTCAD::app;

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
    const auto dir = app_test::scratch("solid_check_" + name);
    const auto p = dir / (name + ".deck");
    std::ofstream(p, std::ios::binary) << text;
    return p;
}

bool has(const std::string& text, const std::string& what) { return text.find(what) != std::string::npos; }

bool checked_now(MainWindow& w, SolidPage& page) {
    return pump_until([&] { return !w.live_checker().pending() && page.check() != SolidPage::Check::pending; }, 180.0);
}

bool run_enabled(MainWindow& w) {
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

const std::string mos = R"([mesh]
type = generate3d
layered = on
domain = 0 1um 0 0.6um 0 0.6um
patch gate = 0.3um 0 0.6um  0.7um 0.6um 0.6um
patch body = 0 0 0  1um 0.6um 0
max_size = 0.15um
interface_resolution = 0.02um

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

std::string with(std::string text, const std::string& from, const std::string& to) {
    const auto at = text.find(from);
    REQUIRE(at != std::string::npos);
    text.replace(at, from.size(), to);
    return text;
}

int count_of(const solid::Problem& p, solid::Marker::Kind k) {
    return static_cast<int>(std::ranges::count(p.markers, k, &solid::Marker::kind));
}

}  // namespace

TEST_CASE("app window solid check: invalid layouts carry the build's message and markers; valid ones build (SE7)",
          "[app][window][solid][check]") {
    struct Row {
        const char* what;
        std::string deck;
        bool valid;
        std::string where;  // a word of the build's message
        solid::Marker::Kind marker = solid::Marker::Kind::domain;
        int markers = 0;  // of that kind
        bool point = false;  // a point named
    };
    const std::vector<Row> rows{
        {"the stack", mos, true, ""},
        {"a box reaching the boundary", with(mos, "[doping]", "[region well]\nmaterial = silicon\nbox = 0 0.4um 0 0.3um 0.2um 0.4um\n\n[doping]"), true, ""},
        {"a flat box", with(mos, "box = 0 1um 0 0.6um 0.4um 0.6um", "box = 0 1um 0 0.6um 0.4um 0.4um"), false, "box", solid::Marker::Kind::solid, 1},
        {"a patch on no boundary face", with(mos, "patch gate = 0.3um 0 0.6um  0.7um 0.6um 0.6um", "patch gate = 0.3um 0.1um 0.3um  0.7um 0.5um 0.3um"),
         false, "patch", solid::Marker::Kind::contact, 1},
        {"a degenerate domain", with(mos, "domain = 0 1um 0 0.6um 0 0.6um", "domain = 0 1um 0 0.6um 0 0"), false, "domain", solid::Marker::Kind::domain, 1},
        {"a box inside the oxide's interface plane (a hole in it)",
         with(mos, "[doping]", "[region well]\nmaterial = silicon\nbox = 0.3um 0.6um 0.1um 0.3um 0.4um 0.5um\n\n[doping]"), false, "interface",
         solid::Marker::Kind::point, 1, true},
    };
    auto w = window();
    for (const Row& row : rows) {
        INFO(row.what);
        REQUIRE(w->open_deck(write("table" + std::to_string(&row - rows.data()), row.deck)));
        pump_until([&] { return !w->runner().busy(); });
        w->execute(Command::tab_geometry);
        REQUIRE(w->solid_shown());
        SolidPage& page = w->solid_page();
        REQUIRE(checked_now(*w, page));
        const auto build = app::check_deck(w->deck().text(), w->document().path());
        if (row.valid) {
            for (const auto& d : build.diagnostics) UNSCOPED_INFO(d.message);
            CHECK(build.diagnostics.empty());
            CHECK(page.check() == SolidPage::Check::valid);
            CHECK(page.problems().empty());
            CHECK(page.run_blocked().empty());
            CHECK(run_enabled(*w));
            continue;
        }
        REQUIRE(build.diagnostics.size() == 1);
        REQUIRE(page.check() == SolidPage::Check::invalid);
        REQUIRE(page.problems().size() == 1);
        const solid::Problem& p = page.problems().front();
        INFO(p.message);
        CHECK(p.message == build.diagnostics.front().message);  // the build's own message
        CHECK(has(p.message, row.where));
        CHECK(count_of(p, row.marker) == row.markers);
        if (row.point) {
            // The point lies on the si/ox interface piece with the hole (z = 0.4 um): the oxide, whose face it is, marked.
            const solid::Marker& mk = *std::ranges::find(p.markers, solid::Marker::Kind::point, &solid::Marker::kind);
            CHECK(mk.at.z == 4000);
            CHECK(std::ranges::any_of(p.markers, [](const solid::Marker& k) { return k.kind == solid::Marker::Kind::solid && k.solid.region == 1; }));
        }
        CHECK(has(page.banner(), p.message));
        CHECK(has(page.run_blocked(), p.message));
        CHECK_FALSE(run_enabled(*w));
        CHECK(page.render() >= 0.0);
        CHECK(std::ranges::any_of(page.scene().lines, [](const render::LineBatch& b) { return b.name == "problems" && !b.points.empty(); }));
        const auto before = w->runner().state();  // the open's own check may have failed
        w->execute(Command::run);
        CHECK_FALSE(w->runner().busy());
        CHECK(w->runner().state() == before);
    }
}

TEST_CASE("app window solid check: the mesh preview shows the built mesh's surface, and only for the text it was built from",
          "[app][window][solid][check]") {
    auto w = window();
    REQUIRE(w->open_deck(write("preview", mos)));
    pump_until([&] { return !w->runner().busy(); });
    w->execute(Command::tab_geometry);
    REQUIRE(w->solid_shown());
    SolidPage& page = w->solid_page();
    REQUIRE(checked_now(*w, page));
    REQUIRE(page.check() == SolidPage::Check::valid);
    CHECK_FALSE(page.mesh_shown());
    page.set_mesh_preview(true);
    CHECK(page.mesh_shown());
    CHECK(page.mesh_segments() > 100);
    CHECK(page.banner().empty());
    CHECK(page.render() >= 0.0);
    CHECK(std::ranges::any_of(page.scene().lines, [](const render::LineBatch& b) { return b.name == "mesh" && b.points.size() > 200; }));
    // The si/ox interface is in it (its node pairs lie either side of z = 0.4 um, FD7), away from the boundary.
    const auto& e = page.mesh();
    bool inner = false;
    for (std::size_t k = 0; k + 1 < e.size(); k += 2) {
        const auto in = [](const solid::V3& p) { return std::abs(p[2] - 0.4) < 0.05 && p[0] > 1e-6 && p[0] < 1 - 1e-6 && p[1] > 1e-6 && p[1] < 0.6 - 1e-6; };
        inner = inner || (in(e[k]) && in(e[k + 1]));
    }
    CHECK(inner);
    // A mesh of another text is not shown.
    const auto kept = page.mesh();
    page.set_mesh("another text", kept);
    CHECK_FALSE(page.mesh_shown());
    page.set_mesh(w->deck().text(), kept);
    CHECK(page.mesh_shown());
    // An edit: the mesh is of the text before until the next check.
    REQUIRE(page.select_region("ox"));
    SetWindowTextW(page.control(SolidPage::c_z0), L"0.45");
    page.commit(SolidPage::c_z0);
    CHECK_FALSE(page.mesh_shown());
    CHECK(has(page.banner(), "Building the mesh"));
    REQUIRE(checked_now(*w, page));
    CHECK(page.mesh_shown());
}

TEST_CASE("app window solid check: a 3D MOSFET drawn on the page runs its equilibrium and one gate bias (SE8)",
          "[app][window][solid][check][run]") {
    const std::string start = R"(# a 3D MOSFET drawn on the Geometry tab
[mesh]
type = generate3d
layered = on
domain = 0 1um 0 0.6um 0 0.6um
max_size = 0.1um
interface_resolution = 0.01um

[region substrate]
material = silicon

[doping]
profile = acceptor constant 1e16

[equilibrium]
)";
    auto w = window();
    REQUIRE(w->open_deck(write("mosfet", start)));
    pump_until([&] { return !w->runner().busy(); });
    w->execute(Command::tab_geometry);
    REQUIRE(w->solid_shown());
    SolidPage& page = w->solid_page();
    REQUIRE(page.model().editable());
    REQUIRE(page.render() >= 0.0);
    const auto deck = [&] { return w->deck().text(); };
    const auto px = [&](const solid::V3& p) {
        const auto s = page.screen_of_um(p);
        return MAKELPARAM(static_cast<int>(std::lround(s[0])), static_cast<int>(std::lround(s[1])));
    };
    const auto drag = [&](const solid::V3& a, const solid::V3& b) {
        SendMessageW(page.canvas(), WM_LBUTTONDOWN, MK_LBUTTON, px(a));
        SendMessageW(page.canvas(), WM_MOUSEMOVE, MK_LBUTTON, px({(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2}));
        SendMessageW(page.canvas(), WM_MOUSEMOVE, MK_LBUTTON, px(b));
        SendMessageW(page.canvas(), WM_LBUTTONUP, 0, px(b));
    };
    const auto name = [&](const wchar_t* n) {
        SetWindowTextW(page.control(SolidPage::c_name), n);
        page.commit(SolidPage::c_name);
    };

    // The oxide: a box drawn over the whole top on the work plane z = 0.5 µm, its height clicked, then its bottom
    // typed as 0.55 µm in the Geometry section; named.
    page.set_material("sio2");
    page.set_work_plane({2, 5000});
    page.set_tool(SolidPage::Tool::box);
    drag({0.0, 0.0, 0.5}, {1.0, 0.6, 0.5});
    SendMessageW(page.canvas(), WM_MOUSEMOVE, 0, px({1.0, 0.6, 0.6}));
    SendMessageW(page.canvas(), WM_LBUTTONDOWN, MK_LBUTTON, px({1.0, 0.6, 0.6}));
    SendMessageW(page.canvas(), WM_LBUTTONUP, 0, px({1.0, 0.6, 0.6}));
    REQUIRE(has(deck(), "material = sio2\nbox = 0 1um 0 0.6um 0.5um 0.6um\n"));
    SetWindowTextW(page.control(SolidPage::c_z0), L"0.55");
    page.commit(SolidPage::c_z0);
    REQUIRE(has(deck(), "box = 0 1um 0 0.6um 0.55um 0.6um\n"));
    name(L"oxide");
    REQUIRE(has(deck(), "[region oxide]\n"));

    // The contacts, each on a face seen straight on: the source and drain on the ends below the oxide, the gate on
    // the top over the channel, the body on the bottom.
    page.set_tool(SolidPage::Tool::contact);
    const auto look = [&](float azimuth, float elevation) {
        render::Camera c = page.view().camera;
        c.azimuth_deg = azimuth;
        c.elevation_deg = elevation;
        c.zoom = 1.0f;
        c.pan = {0.0f, 0.0f};
        page.set_camera(c);
        page.render();
    };
    look(-90.0f, 0.0f);  // from -x
    drag({0.0, 0.0, 0.2}, {0.0, 0.6, 0.5});
    name(L"source");
    look(90.0f, 0.0f);  // from +x
    drag({1.0, 0.0, 0.2}, {1.0, 0.6, 0.5});
    name(L"drain");
    look(0.0f, 0.0f);  // from +z
    drag({0.3, 0.0, 0.6}, {0.7, 0.6, 0.6});
    name(L"gate");
    SendMessageW(page.control(SolidPage::c_kind), CB_SELECTSTRING, static_cast<WPARAM>(-1), reinterpret_cast<LPARAM>(L"electrode"));
    page.commit(SolidPage::c_kind);
    look(180.0f, 0.0f);  // from -z
    drag({0.0, 0.0, 0.0}, {1.0, 0.6, 0.0});
    name(L"body");
    REQUIRE(page.model().contacts.size() == 4);
    REQUIRE(has(deck(), "[contact gate]\nkind = electrode\n"));

    // What the modeller does not draw: the gate's metal, the pockets' doping, the bias, typed.
    replace_text(*w, "[contact gate]\nkind = electrode\n", "[contact gate]\nkind = electrode\nelectrode = metal\nwork_function = 4.1\n");
    replace_text(*w, "profile = acceptor constant 1e16\n",
                 "profile = acceptor constant 1e16\n"
                 "profile = donor gaussian 2e19 axis=z center=0.55um width=0.15um x=0:0.3um lateral=20nm\n"
                 "profile = donor gaussian 2e19 axis=z center=0.55um width=0.15um x=0.7um:1um lateral=20nm\n");
    replace_text(*w, "[equilibrium]\n", "[sweep]\ncontact gate = 0 : 0.5 : 0.5\n");
    w->execute(Command::tab_geometry);
    REQUIRE(checked_now(*w, page));
    for (const auto& p : page.problems()) UNSCOPED_INFO(p.message);
    REQUIRE(page.check() == SolidPage::Check::valid);
    REQUIRE(run_enabled(*w));
    REQUIRE(w->workspace());
    const std::size_t nodes = w->workspace()->device->mesh().node_count();
    UNSCOPED_INFO("nodes " << nodes);
    CHECK(nodes <= 20000);

    const auto t0 = std::chrono::steady_clock::now();
    w->execute(Command::run);
    REQUIRE(pump_until([&] { return !w->runner().busy() && w->workspace() && w->workspace()->has_results(); }, 1800.0));
    CHECK(w->runner().state() == app::RunState::done);
    CHECK_FALSE(w->preview());
    CHECK_FALSE(w->runner().stopped_early());
    // The equilibrium, then the sweep's two points (the gate at 0 and 0.5 V), each converged.
    const auto& rows = w->runner().points();
    REQUIRE(rows.size() == 3);
    CHECK(rows[0].phase == solve::Phase::equilibrium);
    CHECK(std::ranges::all_of(rows, [](const app::PointRow& r) { return r.converged; }));
    CHECK(std::ranges::count_if(rows, [](const app::PointRow& r) { return r.phase != solve::Phase::equilibrium; }) == 2);
    WARN("SE8: the drawn 3D MOSFET, " << nodes << " nodes, ran in "
                                      << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << " s");
}
