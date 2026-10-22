// A window's swap chain, the camera's zoom and pan, and screen positions (Unit 26, DECISIONS.md VA4, VA8).
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <windows.h>

#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/render/image.hpp"
#include "NiTCAD/render/scene.hpp"

using namespace NiTCAD;
using namespace NiTCAD::render;
using Catch::Matchers::WithinAbs;

namespace {

std::vector<double> uniform(double a, double b, int n) {
    std::vector<double> v;
    for (int k = 0; k < n; ++k) v.push_back(k == n - 1 ? b : a + (b - a) * k / (n - 1));
    return v;
}

// A 2D silicon square, 2 x 1 um, with a contact on each side. The geometry refers to the device's mesh,
// so both live where they are made.
struct Square {
    const device::Device& device;
    Geometry geometry;
};

const Square& square() {
    static const device::Device square_device = [] {
        const auto x = uniform(0.0, 2e-4, 9), y = uniform(0.0, 1e-4, 5);
        mesh::Mesh m = *mesh::make_tensor_grid(x, y);
        auto cells = *mesh::TensorCells::create(m, x, y);
        const std::size_t n = m.node_count();
        auto left = m.find_boundary("x_min")->nodes, right = m.find_boundary("x_max")->nodes;
        auto d = device::Device::create({.mesh = std::move(m),
                                         .temperature_K = 300.0,
                                         .regions = {{"silicon", physics::silicon()}},
                                         .node_region = std::vector<device::RegionId>(n, 0),
                                         .donors = std::vector<double>(n, 1e16),
                                         .acceptors = std::vector<double>(n, 0.0),
                                         .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                                                      {"right", device::ContactKind::ohmic, std::move(right)}},
                                         .cells = std::move(cells)});
        if (!d) FAIL(d.error().message);
        return std::move(*d);
    }();
    static const Square s{square_device, *Geometry::tensor(square_device.mesh(), *square_device.cells())};
    return s;
}

Renderer& warp() {
    static auto r = Renderer::create({.adapter = Adapter::warp});
    REQUIRE(r);
    return *r;
}

}  // namespace

TEST_CASE("render camera: zoom and pan about the view's centre; screen positions", "[render]") {
    const Square& s = square();
    auto fitted = solid_view(s.device, s.geometry);
    REQUIRE(fitted);
    constexpr int W = 640, H = 480;
    // A 2D map without a colorbar fills [66, W - 12] x [top, H - 48], top below the title band (36 with a
    // title, else 12); the fitted bounds are centred in it.
    const double top = fitted->title.empty() ? 12.0 : 36.0;
    const ScreenPoint c = screen_position(*fitted, W, H, {1.0f, 0.5f, 0.0f});
    CHECK_THAT(c.x, WithinAbs(0.5 * (66.0 + W - 12.0), 1e-3));
    CHECK_THAT(c.y, WithinAbs(0.5 * (top + H - 48.0), 1e-3));
    const ScreenPoint k = screen_position(*fitted, W, H, {0.0f, 0.0f, 0.0f});
    FieldScene zoomed = *fitted;
    zoomed.camera.zoom = 4.0f;
    const ScreenPoint kz = screen_position(zoomed, W, H, {0.0f, 0.0f, 0.0f});
    CHECK_THAT(kz.x - c.x, WithinAbs(4.0 * (k.x - c.x), 1e-2));
    CHECK_THAT(kz.y - c.y, WithinAbs(4.0 * (k.y - c.y), 1e-2));
    FieldScene panned = *fitted;
    panned.camera.pan = {-1.0f, 0.5f};  // a half width left, a quarter height up
    const ScreenPoint cp = screen_position(panned, W, H, {1.0f, 0.5f, 0.0f});
    CHECK_THAT(cp.x - c.x, WithinAbs(-0.5 * (W - 12.0 - 66.0), 1e-2));
    CHECK_THAT(cp.y - c.y, WithinAbs(-0.25 * (H - 48.0 - top), 1e-2));

    // Zoom 1 and no pan is the fitted image to the bit; a zoomed one differs and draws.
    const Image a = *warp().draw(*fitted, W, H);
    FieldScene unit = *fitted;
    unit.camera.zoom = 1.0f;
    unit.camera.pan = {0.0f, 0.0f};
    CHECK(warp().draw(unit, W, H)->rgba == a.rgba);
    const auto z = warp().draw(zoomed, W, H);
    REQUIRE(z);
    CHECK(z->rgba != a.rgba);
    // Zoomed and panned so the bounds run out of the view on the right: the axes stop at the view, so
    // the right margin (12 pixels, no colorbar) stays white as in the fitted image.
    FieldScene side = *fitted;
    side.camera.zoom = 2.0f;
    side.camera.pan = {0.6f, 0.0f};
    const auto s2 = warp().draw(side, W, H);
    REQUIRE(s2);
    const auto white_right_margin = [&](const Image& im) {
        for (int y = 0; y < H; ++y) {
            for (int x = W - 10; x < W; ++x) {
                const std::size_t i = 4 * (static_cast<std::size_t>(y) * W + static_cast<std::size_t>(x));
                if (im.rgba[i] < 250 || im.rgba[i + 1] < 250 || im.rgba[i + 2] < 250) return false;
            }
        }
        return true;
    };
    CHECK(white_right_margin(a));
    CHECK(white_right_margin(*s2));
    // Panned entirely out of view: the axes have nothing to label, and drawing still works.
    FieldScene away = *fitted;
    away.camera.pan = {10.0f, 10.0f};
    CHECK(warp().draw(away, W, H).has_value());
}

TEST_CASE("render window: a swap chain shows field and plot scenes and follows resizes", "[render]") {
    const Square& s = square();
    HWND hwnd = CreateWindowExW(0, L"STATIC", L"", WS_OVERLAPPEDWINDOW, 0, 0, 400, 300, nullptr, nullptr, nullptr, nullptr);
    REQUIRE(hwnd != nullptr);
    {
        auto target = warp().window_target(hwnd, 320, 240);
        if (!target) {
            DestroyWindow(hwnd);
            SKIP("no swap chain here: " + target.error().message);
        }
        CHECK(target->width() == 320);
        CHECK(target->height() == 240);
        auto scene = solid_view(s.device, s.geometry);
        REQUIRE(scene);
        CHECK(warp().present(*scene, *target, 320, 240).has_value());
        CHECK(warp().present(*scene, *target, 500, 260).has_value());  // resized
        CHECK(target->width() == 500);
        const analysis::Curve curve{{0.0, 1.0, 2.0}, {1.0, 4.0, 9.0}, {}};
        const std::vector<SeriesInput> series{{"y", curve}};
        auto plot = curve_plot(series, "x", "y");
        REQUIRE(plot);
        CHECK(warp().present(*plot, *target, 500, 260).has_value());
        // The image drawn afterwards is unaffected by the window.
        CHECK(warp().draw(*scene, 320, 240).has_value());
        CHECK_FALSE(warp().present(*scene, *target, 8, 8).has_value());  // below 16 pixels
    }
    DestroyWindow(hwnd);
    CHECK_FALSE(warp().window_target(nullptr, 8, 8).has_value());
}
