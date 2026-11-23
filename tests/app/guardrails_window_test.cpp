// Unit 26n slice 4 (DT8, DC7): drawing keeps solids inside the domain, by window messages: 500 seeded random drags of
// the Box, Prism and Cylinder tools and of the move gizmo; typed out-of-domain sizes refused with the reason; typed
// references and Name this value in the inspector; the 2D editor's rectangle held at the domain.
#include <catch2/catch_test_macros.hpp>

#include "main_window.hpp"
#include "fixtures.hpp"
#include "NiTCAD/app/commands.hpp"
#include "NiTCAD/app/solid_view.hpp"
#include "win32.hpp"

#include <windows.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <optional>
#include <random>
#include <string>

using namespace NiTCAD::app;
namespace fs = std::filesystem;

namespace {

void pump() {
    MSG m;
    while (PeekMessageW(&m, nullptr, 0, 0, PM_REMOVE)) {
        TranslateMessage(&m);
        DispatchMessageW(&m);
    }
}

std::unique_ptr<MainWindow> window() {
    WindowOptions o;
    o.adapter = NiTCAD::render::Adapter::warp;
    o.show = false;
    o.quit_on_close = false;
    auto w = MainWindow::create(o);
    REQUIRE(w);
    return w;
}

constexpr solid::Coord um = solid::per_um;

const std::string deck3d = R"([parameters]
W = 0.3um  # Width

[mesh]
type = generate3d
layered = on
domain = 0 1um 0 0.6um 0 0.6um
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

[equilibrium]
)";

struct Mouse {
    SolidPage& page;
    [[nodiscard]] std::array<double, 2> px(const solid::V3& p) const { return page.screen_of_um(p); }
    [[nodiscard]] static LPARAM lp(std::array<double, 2> p) {
        const auto c = [](double v) { return static_cast<int>(std::lround(std::clamp(v, -30000.0, 30000.0))); };
        return MAKELPARAM(static_cast<short>(c(p[0])), static_cast<short>(c(p[1])));
    }
    void down(std::array<double, 2> p) const { SendMessageW(page.canvas(), WM_LBUTTONDOWN, MK_LBUTTON, lp(p)); }
    void move(std::array<double, 2> p, WPARAM keys = MK_LBUTTON) const { SendMessageW(page.canvas(), WM_MOUSEMOVE, keys, lp(p)); }
    void up(std::array<double, 2> p) const { SendMessageW(page.canvas(), WM_LBUTTONUP, 0, lp(p)); }
    void click(std::array<double, 2> p) const {
        down(p);
        up(p);
    }
    void drag(std::array<double, 2> a, std::array<double, 2> b) const {
        down(a);
        move({(a[0] + b[0]) / 2, (a[1] + b[1]) / 2});
        move(b);
        up(b);
    }
    void key(UINT k) const { SendMessageW(page.canvas(), WM_KEYDOWN, k, 0); }
};

}  // namespace

TEST_CASE("guardrails window: 500 seeded random drags never leave a solid outside the domain (DC7)", "[app][window][guardrails]") {
    const fs::path dir = app_test::scratch("guardrails_drags");
    auto w = window();
    REQUIRE(w->open_deck(app_test::write(dir, "g.deck", deck3d)));
    pump();
    w->execute(Command::tab_geometry);
    REQUIRE(w->solid_shown());
    SolidPage& page = w->solid_page();
    REQUIRE(page.model().editable());
    const Mouse mouse{page};
    std::mt19937 rng(2611u);
    const auto uni = [&](double a, double b) { return std::uniform_real_distribution<double>(a, b)(rng); };
    // A point anywhere from half a domain before it to half a domain past it.
    const auto anywhere = [&](double z) -> solid::V3 { return {uni(-0.5, 1.5), uni(-0.3, 0.9), z}; };
    int edits = 0, held = 0;
    std::string before = w->deck().text();
    for (int k = 0; k < 500; ++k) {
        const int what = std::uniform_int_distribution<int>(0, 3)(rng);
        const double z = std::round(uni(0.0, 0.6) * 10.0) / 10.0;  // a work plane on the 0.1 µm grid
        page.render();
        if (what == 0) {  // Box: a base, then a height anywhere along z, a click
            page.set_work_plane({2, static_cast<solid::Coord>(std::lround(z * um))});
            page.set_tool(SolidPage::Tool::box);
            const auto a = anywhere(z), b = anywhere(z);
            mouse.drag(mouse.px(a), mouse.px(b));
            const solid::V3 h{b[0], b[1], uni(-1.0, 1.6)};
            mouse.move(mouse.px(h), 0);
            mouse.click(mouse.px(h));
        } else if (what == 1) {  // Prism: three clicks, Enter, a height, a click
            page.set_work_plane({2, static_cast<solid::Coord>(std::lround(z * um))});
            page.set_tool(SolidPage::Tool::prism);
            for (int j = 0; j < 3; ++j) mouse.click(mouse.px(anywhere(z)));
            mouse.key(VK_RETURN);
            const solid::V3 h{uni(0.0, 1.0), uni(0.0, 0.6), uni(-1.0, 1.6)};
            mouse.move(mouse.px(h), 0);
            mouse.click(mouse.px(h));
        } else if (what == 2) {  // Cylinder: a centre, a rim anywhere, a height, a click
            page.set_work_plane({2, static_cast<solid::Coord>(std::lround(z * um))});
            page.set_tool(SolidPage::Tool::cylinder);
            const solid::V3 c{uni(0.0, 1.0), uni(0.0, 0.6), z};
            mouse.drag(mouse.px(c), mouse.px(anywhere(z)));
            const solid::V3 h{c[0], c[1], uni(-1.0, 1.6)};
            mouse.move(mouse.px(h), 0);
            mouse.click(mouse.px(h));
        } else {  // The move gizmo: a solid picked, an arrow dragged far along its axis
            page.set_tool(SolidPage::Tool::select);
            const auto& m = page.model();
            std::vector<solid::SolidRef> all;
            for (int r = 0; r < static_cast<int>(m.regions.size()); ++r) {
                for (int s = 0; s < static_cast<int>(m.regions[static_cast<std::size_t>(r)].solids.size()); ++s) all.push_back({r, s});
            }
            if (all.empty()) continue;
            const solid::SolidRef pick = all[std::uniform_int_distribution<std::size_t>(0, all.size() - 1)(rng)];
            REQUIRE(page.select_region(m.regions[static_cast<std::size_t>(pick.region)].name));
            page.render();
            const solid::SolidRef sel = page.selection().solids.at(0);
            const int axis = std::uniform_int_distribution<int>(0, 2)(rng);
            const auto hs = solid::handles(page.model(), page.view(), sel);
            const auto arrow = std::ranges::find_if(hs, [&](const solid::Handle& h) { return h.kind == solid::HandleKind::axis && h.index == axis; });
            if (arrow == hs.end()) continue;
            std::optional<std::array<double, 2>> grab;
            for (double t = 0.95; t > 0.3 && !grab; t -= 0.05) {
                const auto a = mouse.px(arrow->pts[0]), b = mouse.px(arrow->pts[1]);
                const std::array<double, 2> p{std::round(a[0] + t * (b[0] - a[0])), std::round(a[1] + t * (b[1] - a[1]))};
                const solid::Pick q = solid::pick(page.model(), page.view(), p[0], p[1], page.selection());
                if (q.kind == solid::PickKind::handle && q.handle == solid::HandleKind::axis && q.index == axis) grab = p;
            }
            if (!grab) continue;
            const solid::V3 o = arrow->pts[0];
            solid::V3 to = o;
            to[static_cast<std::size_t>(axis)] += uni(-2.0, 2.0);
            const auto a = mouse.px(o), b = mouse.px(to);
            mouse.drag(*grab, {(*grab)[0] + b[0] - a[0], (*grab)[1] + b[1] - a[1]});
        }
        mouse.key(VK_ESCAPE);  // whatever is half done
        const std::string now = w->deck().text();
        INFO("drag " << k << " (" << what << ")\n" << now);
        const auto read = solid::read_model(now);
        REQUIRE(read);
        CHECK(solid::outside_domain(*read).empty());
        CHECK(solid::outside_domain(page.model()).empty());
        if (now != before) ++edits;
        if (!page.refused().empty()) ++held;
        before = now;
    }
    std::printf("DC7: 500 drags, %d edits written, %d refusals seen\n", edits, held);
    CHECK(edits >= 100);  // the drags did draw and move things
}

TEST_CASE("guardrails window: a typed size past the domain is refused with the reason; references and Name this value (DT8)",
          "[app][window][guardrails]") {
    const fs::path dir = app_test::scratch("guardrails_typed");
    auto w = window();
    REQUIRE(w->open_deck(app_test::write(dir, "t.deck", deck3d)));
    pump();
    w->execute(Command::tab_geometry);
    REQUIRE(w->solid_shown());
    SolidPage& page = w->solid_page();
    REQUIRE(page.select_region("ox"));
    const std::string original = w->deck().text();

    // Past the domain: refused, the reason in the status, the deck unchanged, the field back to its value.
    win32::set_text(page.control(SolidPage::c_z1), "0.9");
    page.commit(SolidPage::c_z1);
    CHECK(w->deck().text() == original);
    CHECK(page.refused().find("region 'ox' reaches 0.3um past the domain's z_max face") != std::string::npos);
    CHECK(win32::text_of(page.control(SolidPage::c_z1)) == "0.6");

    // A reference typed: written as typed, shown as typed.
    win32::set_text(page.control(SolidPage::c_y1), "$W");
    page.commit(SolidPage::c_y1);
    CHECK(w->deck().text().find("box = 0 1um 0 $W 0.4um 0.6um\n") != std::string::npos);
    CHECK(win32::text_of(page.control(SolidPage::c_y1)) == "$W");
    // A number typed over it replaces it.
    win32::set_text(page.control(SolidPage::c_y1), "0.5");
    page.commit(SolidPage::c_y1);
    CHECK(w->deck().text().find("box = 0 1um 0 0.5um 0.4um 0.6um\n") != std::string::npos);
    // A reference past the domain is refused as a number is.
    win32::set_text(page.control(SolidPage::c_x1), "$W + 1um");
    page.commit(SolidPage::c_x1);
    CHECK(w->deck().text().find("box = 0 1um 0 0.5um 0.4um 0.6um\n") != std::string::npos);
    CHECK(page.refused().find("x_max") != std::string::npos);

    // Name this value: a new parameter with the field's value, the field its reference; one undo takes both back.
    const std::string named_before = w->deck().text();
    CHECK(page.can_name_field(SolidPage::c_z0));
    CHECK_FALSE(page.can_name_field(SolidPage::c_name));
    REQUIRE(page.name_field(SolidPage::c_z0));
    CHECK(w->deck().text().find("\nL1 = 0.4um\n") != std::string::npos);
    CHECK(w->deck().text().find("box = 0 1um 0 0.5um $L1 0.6um\n") != std::string::npos);
    CHECK(win32::text_of(page.control(SolidPage::c_z0)) == "$L1");
    REQUIRE(page.undo());
    CHECK(w->deck().text() == named_before);
}

TEST_CASE("guardrails window: the 2D editor's rectangle stops at the domain (DT8)", "[app][window][guardrails]") {
    const std::string deck2d = R"([mesh]
type = generate2d
layered = on
domain = 0 2um 0 1um
max_size = 0.25um
interface_resolution = 20nm

[region p]
material = silicon

[region oxide]
material = sio2
box = 0 2um 0.8um 1um

[doping]
profile = acceptor constant 1e16

[equilibrium]
)";
    const fs::path dir = app_test::scratch("guardrails_2d");
    auto w = window();
    REQUIRE(w->open_deck(app_test::write(dir, "g2.deck", deck2d)));
    pump();
    w->execute(Command::tab_geometry);
    GeometryPage& page = w->geometry_page();
    REQUIRE(page.geometry().mesh == geo::MeshKind::layered);
    page.set_view(geo::View{200.0, 50.0, 450.0, 96});
    const auto at = [&](geo::Coord x, geo::Coord y) {
        return MAKELPARAM(static_cast<int>(std::lround(page.view().sx(x))), static_cast<int>(std::lround(page.view().sy(y))));
    };
    page.set_tool(GeometryPage::Tool::rectangle);
    SendMessageW(page.canvas(), WM_LBUTTONDOWN, MK_LBUTTON, at(4000, 2000));
    SendMessageW(page.canvas(), WM_MOUSEMOVE, MK_LBUTTON, at(15000, 6000));
    SendMessageW(page.canvas(), WM_MOUSEMOVE, MK_LBUTTON, at(30000, 6000));  // 1 µm past x_max
    SendMessageW(page.canvas(), WM_LBUTTONUP, 0, at(30000, 6000));
    INFO(w->deck().text());
    CHECK(w->deck().text().find("box = 0.4um 2um 0.2um 0.6um\n") != std::string::npos);
    CHECK(geo::outside_domain(page.geometry()).empty());
}

TEST_CASE("guardrails window: a box drawn and moved past the domain lands on its faces (DT8)", "[app][window][guardrails]") {
    const fs::path dir = app_test::scratch("guardrails_faces");
    auto w = window();
    REQUIRE(w->open_deck(app_test::write(dir, "f.deck", deck3d)));
    pump();
    w->execute(Command::tab_geometry);
    REQUIRE(w->solid_shown());
    SolidPage& page = w->solid_page();
    const Mouse mouse{page};
    // The base dragged 0.5 µm past x_max, the height 1 µm past z_max: the box ends on both faces.
    page.set_work_plane({2, 0});
    page.set_tool(SolidPage::Tool::box);
    page.render();
    mouse.drag(mouse.px({0.5, 0.1, 0.0}), mouse.px({1.5, 0.3, 0.0}));
    mouse.move(mouse.px({1.5, 0.3, 1.6}), 0);
    mouse.click(mouse.px({1.5, 0.3, 1.6}));
    INFO(w->deck().text());
    REQUIRE(w->deck().text().find("box = 0.5um 1um 0.1um 0.3um 0 0.6um\n") != std::string::npos);
    // Its x arrow dragged 2 µm towards x_min: the box stops on x_min.
    page.set_tool(SolidPage::Tool::select);
    REQUIRE(page.selection().solids.size() == 1);
    const solid::SolidRef box = page.selection().solids[0];
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
    const auto a = mouse.px(o), b = mouse.px({o[0] - 2.0, o[1], o[2]});
    mouse.drag(*grab, {(*grab)[0] + b[0] - a[0], (*grab)[1] + b[1] - a[1]});
    INFO(w->deck().text());
    CHECK(w->deck().text().find("box = 0 0.5um 0.1um 0.3um 0 0.6um\n") != std::string::npos);
}
