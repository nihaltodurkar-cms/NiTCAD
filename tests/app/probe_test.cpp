// Probe, measure and the camera moves (Unit 26, VA8): picking the node drawn under a pixel, its value,
// distances; zoom and pan of the view.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/app/probe.hpp"
#include "NiTCAD/app/view_settings.hpp"
#include "NiTCAD/render/image.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
namespace fs = std::filesystem;
using Catch::Matchers::WithinAbs;

namespace {

app::Workspace unrun(const std::string& name, const std::string& text) {
    const fs::path dir = app_test::scratch(name);
    auto w = app::workspace_from_deck(app_test::write(dir, "d.deck", text));
    REQUIRE(w);
    return std::move(*w);
}

const std::string bar3d = R"(# A 3D bar
[mesh]
type = tensor
x = linspace(0, 2um, 5)
y = linspace(0, 1um, 3)
z = linspace(0, 1um, 3)
[region silicon]
material = silicon
[doping]
profile = donor constant 1e16
[contact left]
kind = ohmic
patch = x_min
[contact right]
kind = ohmic
patch = x_max
[equilibrium]
)";

}  // namespace

TEST_CASE("app probe: the node under a pixel in 2D, its value; the distance between two", "[app]") {
    const app::Workspace w = unrun("probe2d", app_test::mosfet2d);
    const app::Prober prober(w);
    CHECK(prober.candidates().size() == w.device->mesh().node_count());
    app::ViewSettings s;
    s.kind = app::ViewKind::field;
    s.field = "doping";
    auto scene = app::build_scene(w, s);
    REQUIRE(scene);
    const auto& f = std::get<render::FieldScene>(*scene);
    auto values = render::field_values(*w.device, results::BiasPoint{}, render::Field::net_doping);
    REQUIRE(values);
    constexpr int W = 800, H = 600;
    // Every node is found where it is drawn.
    for (const std::size_t n : {std::size_t{0}, std::size_t{220}, std::size_t{440}}) {
        const auto& p = w.device->mesh().points()[n];
        const auto at = render::screen_position(f, W, H, {static_cast<float>(p[0] * 1e4), static_cast<float>(p[1] * 1e4), 0.0f});
        const auto hit = prober.pick(f, W, H, at.x + 2.0f, at.y - 1.0f, *values);
        REQUIRE(hit);
        CHECK(hit->node == static_cast<mesh::NodeId>(n));
        CHECK(hit->value == (*values)[n]);
        CHECK_THAT(hit->position_um[0], WithinAbs(p[0] * 1e4, 1e-12));
        CHECK(hit->distance_px < 3.0f);
    }
    // Far from every node: nothing.
    CHECK_FALSE(prober.pick(f, W, H, -100.0f, -100.0f));
    // The distance between the first node (0, 0) and the last (2, 1) um.
    const auto& p0 = w.device->mesh().points().front();
    const auto& p1 = w.device->mesh().points().back();
    const auto a0 = render::screen_position(f, W, H, {static_cast<float>(p0[0] * 1e4), static_cast<float>(p0[1] * 1e4), 0.0f});
    const auto a1 = render::screen_position(f, W, H, {static_cast<float>(p1[0] * 1e4), static_cast<float>(p1[1] * 1e4), 0.0f});
    const auto h0 = prober.pick(f, W, H, a0.x, a0.y);
    const auto h1 = prober.pick(f, W, H, a1.x, a1.y);
    REQUIRE(h0);
    REQUIRE(h1);
    CHECK(std::isnan(h0->value));  // no values given
    CHECK_THAT(app::distance_um(*h0, *h1), WithinAbs(std::sqrt(5.0), 1e-12));
}

TEST_CASE("app probe 26b: the point of a 2D view under a pixel, from any camera (XD5)", "[app]") {
    const app::Workspace w = unrun("point2d", app_test::mosfet2d);
    const app::Prober prober(w);
    app::ViewSettings s;
    s.zoom = 1.7f;
    s.pan_x = 0.2f;
    auto scene = app::build_scene(w, s);
    REQUIRE(scene);
    auto f = std::get<render::FieldScene>(*scene);
    constexpr int W = 800, H = 600;
    // The view the window draws (looking down -z, zoomed and panned), then one turned and tilted:
    // there each screen coordinate depends on both x and y, so the whole inverse is needed.
    for (const auto [azimuth, elevation] : {std::pair{0.0f, 0.0f}, std::pair{30.0f, 40.0f}}) {
        f.camera.azimuth_deg = azimuth;
        f.camera.elevation_deg = elevation;
        for (const auto [x, y] : {std::pair{0.0f, 0.0f}, std::pair{1.3f, 0.87f}, std::pair{2.0f, 0.25f}}) {
            const auto at = render::screen_position(f, W, H, {x, y, 0.0f});
            const auto p = prober.point_at(f, W, H, at.x, at.y);
            REQUIRE(p);
            CHECK_THAT((*p)[0], WithinAbs(x, 1e-4));
            CHECK_THAT((*p)[1], WithinAbs(y, 1e-4));
            CHECK((*p)[2] == 0.0);
        }
    }
    // A 1D device has no point under a pixel.
    const app::Workspace d = unrun("point1d", app_test::diode_sweep());
    CHECK_FALSE(app::Prober(d).point_at(f, W, H, 400.0f, 300.0f));
}

TEST_CASE("app probe: 3D picks the front-most boundary node", "[app]") {
    const app::Workspace w = unrun("probe3d", bar3d);
    const app::Prober prober(w);
    // A 5 x 3 x 3 grid: every node but the one inside is on the boundary.
    CHECK(prober.candidates().size() == 45 - 3);
    app::ViewSettings s;  // the solid view, from the default camera
    auto scene = app::build_scene(w, s);
    REQUIRE(scene);
    const auto& f = std::get<render::FieldScene>(*scene);
    // The camera looks from +z (azimuth 35, elevation 25): on the pixel of the corner (0, 1, 1) um,
    // that node, not one behind it.
    const auto at = render::screen_position(f, 800, 600, {0.0f, 1.0f, 1.0f});
    const auto hit = prober.pick(f, 800, 600, at.x, at.y);
    REQUIRE(hit);
    CHECK_THAT(hit->position_um[0], WithinAbs(0.0, 1e-12));
    CHECK_THAT(hit->position_um[1], WithinAbs(1.0, 1e-12));
    CHECK_THAT(hit->position_um[2], WithinAbs(1.0, 1e-12));
    // Looking straight down -z, a front-face node hides the back-face node at the same pixel: the front
    // one (z = 1 um) is picked.
    app::ViewSettings head_on;
    head_on.azimuth = 0.0f;
    head_on.elevation = 0.0f;
    auto straight = app::build_scene(w, head_on);
    REQUIRE(straight);
    const auto& g = std::get<render::FieldScene>(*straight);
    const auto mid = render::screen_position(g, 800, 600, {1.0f, 0.5f, 1.0f});
    const auto front = prober.pick(g, 800, 600, mid.x, mid.y);
    REQUIRE(front);
    CHECK_THAT(front->position_um[0], WithinAbs(1.0, 1e-12));
    CHECK_THAT(front->position_um[2], WithinAbs(1.0, 1e-12));
    // A 1D device has nothing to pick.
    const app::Workspace d = unrun("probe1d", app_test::diode_sweep());
    CHECK(app::Prober(d).candidates().empty());
}

TEST_CASE("app view settings: orbit, zoom, pan and fit (VA8)", "[app]") {
    app::ViewSettings s;
    app::orbit(s, 10.0f, 5.0f);
    CHECK(s.azimuth == 45.0f);
    CHECK(s.elevation == 30.0f);
    app::orbit(s, 300.0f, 100.0f);
    CHECK(s.azimuth == -15.0f);  // 345 wrapped into (-180, 180]
    CHECK(s.elevation == 89.0f);
    app::zoom_by(s, 2.0f);
    app::zoom_by(s, 4.0f);
    CHECK(s.zoom == 8.0f);
    app::zoom_by(s, 1e6f);
    CHECK(s.zoom == 1024.0f);
    app::zoom_by(s, -1.0f);  // ignored
    CHECK(s.zoom == 1024.0f);
    app::pan_by(s, 0.25f, -0.5f);
    CHECK(s.pan_x == 0.25f);
    CHECK(s.pan_y == -0.5f);
    app::fit(s);
    CHECK(s.zoom == 1.0f);
    CHECK(s.pan_x == 0.0f);
    CHECK(s.azimuth == -15.0f);  // the angles stay
}

TEST_CASE("app view settings: zoom and pan move the picture about the view's centre", "[app]") {
    const app::Workspace w = unrun("zoom2d", app_test::mosfet2d);
    app::ViewSettings s;
    auto fitted = app::build_scene(w, s);
    REQUIRE(fitted);
    const auto& a = std::get<render::FieldScene>(*fitted);
    constexpr int W = 800, H = 600;
    // The centre of the drawn bounds (the contact pads reach beyond the mesh).
    const render::Float3 corner{0.0f, 0.0f, 0.0f},
        centre{0.5f * (a.lo[0] + a.hi[0]), 0.5f * (a.lo[1] + a.hi[1]), 0.0f};
    const auto c0 = render::screen_position(a, W, H, centre);
    const auto k0 = render::screen_position(a, W, H, corner);
    app::zoom_by(s, 2.0f);
    auto zoomed = app::build_scene(w, s);
    const auto& b = std::get<render::FieldScene>(*zoomed);
    const auto c1 = render::screen_position(b, W, H, centre);
    const auto k1 = render::screen_position(b, W, H, corner);
    // The centre of the bounds stays; the corner is twice as far from it.
    CHECK_THAT(c1.x, WithinAbs(c0.x, 1e-3));
    CHECK_THAT(c1.y, WithinAbs(c0.y, 1e-3));
    CHECK_THAT(k1.x - c1.x, WithinAbs(2.0f * (k0.x - c0.x), 1e-2));
    CHECK_THAT(k1.y - c1.y, WithinAbs(2.0f * (k0.y - c0.y), 1e-2));
    // Pan by half the half-width to the right: the picture moves a quarter of the view's width right.
    app::fit(s);
    app::pan_by(s, 0.5f, 0.0f);
    auto panned = app::build_scene(w, s);
    const auto c2 = render::screen_position(std::get<render::FieldScene>(*panned), W, H, centre);
    const float view_width = static_cast<float>(W) - 66.0f - 12.0f;  // the 2D map's margins (no colorbar)
    CHECK_THAT(c2.x - c0.x, WithinAbs(0.25f * view_width, 1e-2));
    CHECK_THAT(c2.y, WithinAbs(c0.y, 1e-3));
    // And up by half the half-height: a quarter of the view's height up (y down in pixels).
    app::pan_by(s, 0.0f, 0.5f);
    auto up = app::build_scene(w, s);
    const auto& u = std::get<render::FieldScene>(*up);
    const auto c3 = render::screen_position(u, W, H, centre);
    const float top = u.title.empty() ? 12.0f : 36.0f;
    const float view_height = static_cast<float>(H) - 48.0f - top;
    CHECK_THAT(c3.y - c0.y, WithinAbs(-0.25f * view_height, 1e-2));
    // A redraw without rebuilding: apply_camera gives the scene build_scene gives.
    render::FieldScene again = a;
    app::apply_camera(again, s);
    CHECK(again.camera.zoom == u.camera.zoom);
    CHECK(again.camera.pan == u.camera.pan);
}

TEST_CASE("app probe 26c: hidden regions are not picked; in 3D the faces they opened are (YD4)", "[app]") {
    for (const bool three : {false, true}) {
        CAPTURE(three);
        const app::Workspace w = unrun(three ? "hide3d" : "hide2d", three ? app_test::mosfet3d : app_test::mosfet2d);
        app::Prober prober(w);
        const std::vector<mesh::NodeId> all(prober.candidates().begin(), prober.candidates().end());
        const auto regions = w.device->regions();
        const auto oxide = static_cast<std::int32_t>(std::ranges::find(regions, std::string("oxide"), &device::Region::name) -
                                                     regions.begin());
        const std::vector<std::int32_t> hidden{oxide};
        prober.hide(hidden);
        const auto picked = [&](mesh::NodeId n) { return std::ranges::find(prober.candidates(), n) != prober.candidates().end(); };
        for (std::size_t i = 0; i < w.device->mesh().node_count(); ++i) {
            const auto n = static_cast<mesh::NodeId>(i);
            if (w.device->node_region()[i] == oxide) CHECK_FALSE(picked(n));
        }
        if (three) {
            // The silicon's top inside (x 1, y 0.9, z 0.4 um), under the oxide: picked once the oxide is hidden.
            const auto points = w.device->mesh().points();
            const auto inside = std::ranges::find_if(points, [](const mesh::Point& p) {
                return std::abs(p[0] - 1e-4) < 1e-9 && std::abs(p[1] - 0.9e-4) < 1e-9 && std::abs(p[2] - 0.4e-4) < 1e-9;
            });
            REQUIRE(inside != points.end());
            const auto n = static_cast<mesh::NodeId>(inside - points.begin());
            CHECK(std::ranges::find(all, n) == all.end());
            CHECK(picked(n));
        } else {
            std::size_t silicon = 0;
            for (const auto r : w.device->node_region()) silicon += r != oxide ? 1 : 0;
            CHECK(prober.candidates().size() == silicon);
        }
        // Shown again: every candidate as before.
        prober.hide({});
        CHECK(std::ranges::equal(prober.candidates(), all));
    }
}
