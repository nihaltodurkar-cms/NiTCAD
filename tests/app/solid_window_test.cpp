// Unit 26l slice 2 (SD4): the 3D modeller on the Geometry tab, driven by its window messages (SE8; no screen
// automation): a box drawn on the work plane and given its height, moved by the gizmo, resized by a face handle;
// a prism drawn and a vertex dragged; a contact drawn on a face of the domain; undo and redo through them all; the
// deck's text as expected after each. SE5: the frame time while orbiting and dragging with 100 solids on WARP.
#include <catch2/catch_test_macros.hpp>

#include "main_window.hpp"
#include "fixtures.hpp"

#include <windows.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <format>
#include <fstream>
#include <string>
#include <vector>

#include "NiTCAD/app/solid.hpp"
#include "NiTCAD/app/solid_view.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/syntax.hpp"
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
    const auto dir = app_test::scratch("solid_window_" + name);
    const auto p = dir / (name + ".deck");
    std::ofstream(p, std::ios::binary) << text;
    return p;
}

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

const std::string extruded = R"([mesh]
type = extrude3d
axis = z
stations = 0, 0.4um, 0.6um
cell left = 0 0, 0.5um 0, 0.5um 0.6um, 0 0.6um
cell right = 0.5um 0, 1um 0, 1um 0.6um, 0.5um 0.6um
patch body = 0 0 0  1um 0.6um 0
max_size = 0.15um

[region si]
material = silicon
extrude = left right 0:0.4um

[region ox]
material = sio2
extrude = left right 0.4um:0.6um

[contact body]
kind = ohmic
patch = body

[equilibrium]
)";

bool has(const std::string& text, const std::string& what) { return text.find(what) != std::string::npos; }

// The canvas's pixel of a model point (µm).
struct Mouse {
    SolidPage& page;
    [[nodiscard]] std::array<double, 2> px(const solid::V3& p) const { return page.screen_of_um(p); }
    [[nodiscard]] static LPARAM lp(std::array<double, 2> p) { return MAKELPARAM(static_cast<int>(std::lround(p[0])), static_cast<int>(std::lround(p[1]))); }
    void down(std::array<double, 2> p, WPARAM keys = 0) const { SendMessageW(page.canvas(), WM_LBUTTONDOWN, MK_LBUTTON | keys, lp(p)); }
    void move(std::array<double, 2> p, WPARAM keys = MK_LBUTTON) const { SendMessageW(page.canvas(), WM_MOUSEMOVE, keys, lp(p)); }
    void up(std::array<double, 2> p) const { SendMessageW(page.canvas(), WM_LBUTTONUP, 0, lp(p)); }
    void click(const solid::V3& p) const {
        down(px(p));
        up(px(p));
    }
    void drag(std::array<double, 2> a, std::array<double, 2> b) const {
        down(a);
        move({(a[0] + b[0]) / 2, (a[1] + b[1]) / 2});
        move(b);
        up(b);
    }
    void key(UINT k) const { SendMessageW(page.canvas(), WM_KEYDOWN, k, 0); }
};

std::array<double, 2> plus(std::array<double, 2> p, std::array<double, 2> a, std::array<double, 2> b) { return {p[0] + b[0] - a[0], p[1] + b[1] - a[1]}; }

const solid::Solid& only_solid(const solid::Model& m, int region) {
    REQUIRE(region >= 0);
    REQUIRE(region < static_cast<int>(m.regions.size()));
    REQUIRE(m.regions[static_cast<std::size_t>(region)].solids.size() == 1);
    return m.regions[static_cast<std::size_t>(region)].solids[0];
}

// The deck read back holds this region's solid.
const solid::Solid& deck_solid(const std::string& text, const std::string& region, solid::Model& keep) {
    auto m = solid::read_model(text);
    REQUIRE(m);
    keep = std::move(*m);
    for (std::size_t r = 0; r < keep.regions.size(); ++r) {
        if (keep.regions[r].name == region) return only_solid(keep, static_cast<int>(r));
    }
    FAIL("no region " << region << " in\n" << text);
    return keep.regions[0].solids[0];
}

constexpr solid::Coord um = solid::per_um;

}  // namespace

TEST_CASE("app window solid: draw a box and its height, move it by the gizmo, resize a face; a prism and its vertex; a contact; "
          "undo, redo (SE8)",
          "[app][window][solid]") {
    auto w = window();
    REQUIRE(w->open_deck(write("draw", mos)));
    pump();
    w->execute(Command::tab_geometry);
    REQUIRE(w->tab() == MainWindow::tab_geometry);
    REQUIRE(w->solid_shown());
    SolidPage& page = w->solid_page();
    CHECK(IsWindowVisible(page.hwnd()) == IsWindowVisible(w->hwnd()));
    CHECK(page.message().empty());
    CHECK(page.model().editable());
    REQUIRE(page.render() >= 0.0);
    const std::string original = w->deck().text();
    std::vector<std::string> texts{original};
    const Mouse mouse{page};

    // A box: its base on the work plane z = 0.4 µm, then its height 0.1 µm, a click ending it.
    page.set_work_plane({2, 4 * um / 10});
    page.set_tool(SolidPage::Tool::box);
    mouse.drag(mouse.px({0.0, 0.0, 0.4}), mouse.px({0.2, 0.3, 0.4}));  // from the ox corner: the si/ox interface keeps no hole
    mouse.move(mouse.px({0.2, 0.3, 0.5}), 0);
    mouse.click({0.2, 0.3, 0.5});
    CHECK(has(w->deck().text(), "\n[region silicon_1]\nmaterial = silicon\nbox = 0 0.2um 0 0.3um 0.4um 0.5um\n\n[doping]"));
    REQUIRE(page.selection().solids.size() == 1);
    const solid::SolidRef box = page.selection().solids[0];
    CHECK(box.region == 2);
    texts.push_back(w->deck().text());

    // Moved +0.2 µm along x by the gizmo's x arrow, grabbed where the pick finds it.
    page.set_tool(SolidPage::Tool::select);
    page.render();
    const auto hs = solid::handles(page.model(), page.view(), box);
    const auto arrow = std::ranges::find_if(hs, [](const solid::Handle& h) { return h.kind == solid::HandleKind::axis && h.index == 0; });
    REQUIRE(arrow != hs.end());
    std::optional<std::array<double, 2>> grab;
    for (double t = 0.95; t > 0.3 && !grab; t -= 0.05) {
        const auto a = mouse.px(arrow->pts[0]), b = mouse.px(arrow->pts[1]);
        const std::array<double, 2> p{std::round(a[0] + t * (b[0] - a[0])), std::round(a[1] + t * (b[1] - a[1]))};
        const solid::Pick k = solid::pick(page.model(), page.view(), p[0], p[1], page.selection());
        if (k.kind == solid::PickKind::handle && k.handle == solid::HandleKind::axis && k.index == 0) grab = p;
    }
    REQUIRE(grab);
    const solid::V3 o = arrow->pts[0];
    mouse.drag(*grab, plus(*grab, mouse.px(o), mouse.px({o[0] + 0.2, o[1], o[2]})));
    CHECK(has(w->deck().text(), "box = 0.2um 0.4um 0 0.3um 0.4um 0.5um\n"));
    texts.push_back(w->deck().text());

    // Its x_max face dragged out 0.1 µm by its handle.
    page.render();
    const solid::V3 face{0.4, 0.15, 0.45};
    const auto at_face = mouse.px(face);
    const solid::Pick fk = solid::pick(page.model(), page.view(), std::round(at_face[0]), std::round(at_face[1]), page.selection());
    CHECK(fk.kind == solid::PickKind::handle);
    CHECK(fk.handle == solid::HandleKind::face);
    CHECK(fk.index == 1);
    mouse.drag(at_face, mouse.px({0.5, 0.15, 0.45}));
    CHECK(has(w->deck().text(), "box = 0.2um 0.5um 0 0.3um 0.4um 0.5um\n"));
    texts.push_back(w->deck().text());

    // A prism: an L clicked on z = 0, closed by Enter, 0.2 µm high.
    page.set_work_plane({2, 0});
    page.set_tool(SolidPage::Tool::prism);
    for (const auto [x, y] : std::vector<std::pair<double, double>>{{0.1, 0.1}, {0.4, 0.1}, {0.4, 0.2}, {0.2, 0.2}, {0.2, 0.4}, {0.1, 0.4}}) {
        mouse.click({x, y, 0.0});
    }
    mouse.key(VK_RETURN);
    mouse.move(mouse.px({0.1, 0.4, 0.2}), 0);
    mouse.click({0.1, 0.4, 0.2});
    {
        solid::Model keep;
        const solid::Solid& p = deck_solid(w->deck().text(), "silicon_2", keep);
        CHECK(p.kind == solid::SolidKind::prism);
        CHECK(p.axis == 2);
        CHECK(p.from == 0);
        CHECK(p.to == 2 * um / 10);
        CHECK(p.polygon == std::vector<solid::Pt>{{um / 10, um / 10}, {4 * um / 10, um / 10}, {4 * um / 10, 2 * um / 10}, {2 * um / 10, 2 * um / 10},
                                                  {2 * um / 10, 4 * um / 10}, {um / 10, 4 * um / 10}});
    }
    CHECK(has(w->deck().text(), "prism = z 0 0.2um, "));
    REQUIRE(page.selection().solids.size() == 1);
    texts.push_back(w->deck().text());

    // Its third vertex dragged 0.1 µm along x on its top cap; the vertical edge's other end follows.
    page.set_tool(SolidPage::Tool::select);
    page.render();
    const auto v2 = mouse.px({0.4, 0.2, 0.2});
    const solid::Pick vk = solid::pick(page.model(), page.view(), std::round(v2[0]), std::round(v2[1]), page.selection());
    CHECK(vk.kind == solid::PickKind::handle);
    CHECK(vk.handle == solid::HandleKind::vertex);
    CHECK(vk.index == 2);
    mouse.drag(v2, mouse.px({0.5, 0.2, 0.2}));
    {
        solid::Model keep;
        const solid::Solid& p = deck_solid(w->deck().text(), "silicon_2", keep);
        CHECK(p.polygon == std::vector<solid::Pt>{{um / 10, um / 10}, {5 * um / 10, um / 10}, {5 * um / 10, 2 * um / 10}, {2 * um / 10, 2 * um / 10},
                                                  {2 * um / 10, 4 * um / 10}, {um / 10, 4 * um / 10}});
    }
    texts.push_back(w->deck().text());

    // A contact dragged on the domain's y_max face.
    page.set_tool(SolidPage::Tool::contact);
    mouse.drag(mouse.px({0.2, 0.6, 0.1}), mouse.px({0.4, 0.6, 0.3}));
    CHECK(has(w->deck().text(), "patch contact_1 = 0.2um 0.6um 0.1um  0.4um 0.6um 0.3um\n"));
    CHECK(has(w->deck().text(), "[contact contact_1]\nkind = ohmic\npatch = contact_1\n"));
    CHECK(page.selection().contact == 2);
    texts.push_back(w->deck().text());
    CHECK(page.render() >= 0.0);

    // Undo through every edit to the deck opened, redo back.
    for (std::size_t k = texts.size() - 1; k > 0; --k) {
        REQUIRE(page.undo());
        CHECK(w->deck().text() == texts[k - 1]);
    }
    CHECK_FALSE(page.can_undo());
    CHECK(w->deck().text() == original);
    for (std::size_t k = 1; k < texts.size(); ++k) {
        w->execute(Command::redo);
        CHECK(w->deck().text() == texts[k]);
    }
    CHECK_FALSE(page.can_redo());
    // The result builds as the deck says: the layered stack and its new contact.
    auto m = solid::read_model(w->deck().text());
    REQUIRE(m);
    CHECK(m->regions.size() == 4);
    CHECK(m->contacts.size() == 3);
    const auto d = deck::parse(w->deck().text(), "s.deck");
    REQUIRE(d);
    const auto b = deck::build(*d);
    if (!b) FAIL(b.error().message);
    CHECK(b->device.mesh().dimension() == 3);
    CHECK(b->device.regions().size() == 4);
    CHECK(b->device.contacts().size() == 3);
}

TEST_CASE("app window solid: the inspector's numbers, the material, the work plane, Delete (SE8)", "[app][window][solid]") {
    auto w = window();
    REQUIRE(w->open_deck(write("inspector", mos)));
    pump();
    w->execute(Command::tab_geometry);
    REQUIRE(w->solid_shown());
    SolidPage& page = w->solid_page();
    REQUIRE(page.select_region("ox"));
    CHECK(win32::text_of(page.control(SolidPage::c_z0)) == "0.4");
    SetWindowTextW(page.control(SolidPage::c_z0), L"0.45");
    page.commit(SolidPage::c_z0);
    CHECK(has(w->deck().text(), "box = 0 1um 0 0.6um 0.45um 0.6um\n"));
    SendMessageW(page.control(SolidPage::c_material), CB_SELECTSTRING, static_cast<WPARAM>(-1), reinterpret_cast<LPARAM>(L"insulator"));
    page.commit(SolidPage::c_material);
    CHECK(has(w->deck().text(), "[region ox]\nmaterial = insulator\n"));
    // A refused value leaves the deck.
    const std::string before = w->deck().text();
    SetWindowTextW(page.control(SolidPage::c_z1), L"0.2");
    page.commit(SolidPage::c_z1);
    CHECK(w->deck().text() == before);
    // The work plane from the inspector.
    SendMessageW(page.control(SolidPage::c_plane), CB_SETCURSEL, 0, 0);
    SetWindowTextW(page.control(SolidPage::c_offset), L"0.25");
    page.commit(SolidPage::c_offset);
    CHECK(page.view().plane == solid::WorkPlane{0, 2500});
    // Delete from the Edit menu.
    REQUIRE(page.select_region("ox"));
    SetFocus(page.canvas());
    w->execute(Command::delete_selection);
    CHECK_FALSE(has(w->deck().text(), "[region ox]"));
    w->execute(Command::undo);
    CHECK(w->deck().text() == before);
}

TEST_CASE("app window solid: an extrude3d deck is shown, not edited, until converted; other decks keep the 2D page", "[app][window][solid]") {
    {
        auto w = window();
        REQUIRE(w->open_deck(write("extruded", extruded)));
        pump();
        w->execute(Command::tab_geometry);
        REQUIRE(w->solid_shown());
        SolidPage& page = w->solid_page();
        CHECK(page.model().mesh == solid::MeshKind::extruded);
        CHECK(page.render() >= 0.0);
        const std::string before = w->deck().text();
        const Mouse mouse{page};
        page.set_tool(SolidPage::Tool::box);
        mouse.down(mouse.px({0.1, 0.1, 0.0}));
        CHECK(has(w->status(3), "edited once converted"));  // the cursor's part says why
        mouse.move(mouse.px({0.3, 0.3, 0.0}));
        mouse.up(mouse.px({0.3, 0.3, 0.0}));
        CHECK(w->deck().text() == before);
        CHECK((GetWindowLongPtrW(page.control(SolidPage::c_convert), GWL_STYLE) & WS_VISIBLE) != 0);
        REQUIRE(page.commit(SolidPage::c_convert));
        CHECK(has(w->deck().text(), "type = generate3d\nlayered = on\n"));
        CHECK(page.model().editable());
        REQUIRE(page.undo());
        CHECK(w->deck().text() == before);
    }
    {
        auto w = window();
        REQUIRE(w->open_deck(std::filesystem::path(NITCAD_EXAMPLES_DIR) / "mosfet_3d.deck"));
        w->execute(Command::tab_geometry);
        CHECK_FALSE(w->solid_shown());
        CHECK(IsWindowVisible(w->solid_page().hwnd()) == FALSE);
    }
}

TEST_CASE("app window solid: orbiting and dragging with 100 solids take at most 16 ms a frame on WARP (SE5)", "[app][window][solid][perf]") {
    std::string deck = "[mesh]\ntype = generate3d\nlayered = on\ndomain = 0 10um 0 10um 0 1um\nmax_size = 0.5um\n\n[region si]\nmaterial = silicon\n";
    for (int k = 0; k < 100; ++k) {
        deck += std::format("\n[region r{}]\nmaterial = {}\nbox = {}um {}.8um {}um {}.8um 0 0.{}um\n", k, k % 2 ? "sio2" : "silicon", k % 10, k % 10,
                            k / 10, k / 10, 3 + k % 7);
    }
    deck += "\n[equilibrium]\n";
    auto w = window();
    REQUIRE(w->open_deck(write("hundred", deck)));
    pump();
    w->execute(Command::tab_geometry);
    REQUIRE(w->solid_shown());
    SolidPage& page = w->solid_page();
    REQUIRE(page.model().regions.size() == 101);
    REQUIRE(page.render() >= 0.0);  // the swap chain and the pipelines
    const auto median = [](std::vector<double> v) {
        std::ranges::sort(v);
        return v[v.size() / 2];
    };
    // Orbiting: right drag, a frame per step.
    std::vector<double> orbit, build;
    SendMessageW(page.canvas(), WM_RBUTTONDOWN, MK_RBUTTON, MAKELPARAM(400, 300));
    for (int k = 1; k <= 31; ++k) {
        SendMessageW(page.canvas(), WM_MOUSEMOVE, MK_RBUTTON, MAKELPARAM(400 + 4 * k, 300 + k));
        const double ms = page.render();
        REQUIRE(ms >= 0.0);
        orbit.push_back(ms);
        build.push_back(page.scene_ms());
    }
    SendMessageW(page.canvas(), WM_RBUTTONUP, 0, MAKELPARAM(400 + 124, 331));
    // Dragging a solid by the gizmo's x arrow, a frame per step.
    REQUIRE(page.select_region("r55"));
    page.render();
    const solid::SolidRef ref = page.selection().solids[0];
    const auto hs = solid::handles(page.model(), page.view(), ref);
    const auto arrow = std::ranges::find_if(hs, [](const solid::Handle& h) { return h.kind == solid::HandleKind::axis && h.index == 0; });
    REQUIRE(arrow != hs.end());
    const auto a = page.screen_of_um(arrow->pts[0]), b = page.screen_of_um(arrow->pts[1]);
    const int gx = static_cast<int>(std::lround(a[0] + 0.85 * (b[0] - a[0]))), gy = static_cast<int>(std::lround(a[1] + 0.85 * (b[1] - a[1])));
    REQUIRE(solid::pick(page.model(), page.view(), gx, gy, page.selection()).kind == solid::PickKind::handle);
    std::vector<double> drag;
    SendMessageW(page.canvas(), WM_LBUTTONDOWN, MK_LBUTTON, MAKELPARAM(gx, gy));
    for (int k = 1; k <= 31; ++k) {
        SendMessageW(page.canvas(), WM_MOUSEMOVE, MK_LBUTTON, MAKELPARAM(gx + 2 * k, gy));
        const double ms = page.render();
        REQUIRE(ms >= 0.0);
        drag.push_back(ms);
        build.push_back(page.scene_ms());
    }
    SendMessageW(page.canvas(), WM_LBUTTONUP, 0, MAKELPARAM(gx + 62, gy));
    std::printf("SE5: 100 solids on WARP: orbit median %.2f ms (max %.2f), drag median %.2f ms (max %.2f); the scene builds' median %.2f ms\n", median(orbit),
                *std::ranges::max_element(orbit), median(drag), *std::ranges::max_element(drag), median(build));
    CHECK(median(orbit) <= 16.0);
    CHECK(median(drag) <= 16.0);
}
