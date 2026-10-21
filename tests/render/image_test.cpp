// Drawing (Unit 25a, VG4, VG7): WARP renders of the views compared with committed references, the
// image pipeline's own checks, the debug layer, and the benchmark.
//
// The references are tests/render/references/*.png, rendered on WARP. A render matches its
// reference when no more than 0.5% of its pixels differ by more than 24 in any channel: rounding in
// the rasterizer and glyph rendering may differ between Windows builds, a wrong view differs over
// whole regions. To write new references, run with NITCAD_UPDATE_REFERENCES=1, look at every image,
// and commit them with the reason.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>

#include "NiTCAD/analysis/curve.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/render/image.hpp"
#include "NiTCAD/render/scene.hpp"
#include "NiTCAD/solve/bias.hpp"

using namespace NiTCAD;
using namespace NiTCAD::render;
namespace fs = std::filesystem;

namespace {

const fs::path references = NITCAD_RENDER_REFERENCES;
const fs::path output = fs::temp_directory_path() / "nitcad_render_test";

Renderer& warp() {
    static Renderer r = [] {
        auto created = Renderer::create({.adapter = Adapter::warp,
#ifndef NDEBUG
                                         .debug_layer = true
#endif
        });
        if (!created) FAIL(created.error().message);
        return std::move(*created);
    }();
    return r;
}

// Draws, saves the render next to the test's output, and compares with the reference (or writes it).
void check_image(const Image& image, const std::string& name) {
    fs::create_directories(output);
    if (auto saved = save_png(image, output / (name + ".png")); !saved) FAIL(saved.error().message);
    const fs::path ref = references / (name + ".png");
    char* update = nullptr;
    std::size_t length = 0;
    const bool updating = _dupenv_s(&update, &length, "NITCAD_UPDATE_REFERENCES") == 0 && update &&
                          std::string(update) == "1";
    std::free(update);
    if (updating) {
        fs::create_directories(references);
        REQUIRE(save_png(image, ref).has_value());
        WARN("reference written: " << ref.string());
        return;
    }
    INFO("reference " << ref.string() << ", render " << (output / (name + ".png")).string());
    const auto expected = load_png(ref);
    REQUIRE(expected.has_value());
    const auto d = compare(image, *expected, 24);
    REQUIRE(d.has_value());
    CAPTURE(d->max_channel, d->fraction_over);
    CHECK(d->fraction_over <= 0.005);
}

std::vector<double> uniform(double a, double b, int n) {
    std::vector<double> v;
    for (int k = 0; k < n; ++k) v.push_back(k == n - 1 ? b : a + (b - a) * k / (n - 1));
    return v;
}

// A small n-channel MOSFET in 2D or 3D (cm): silicon under an oxide 0.1 um thick on top
// (y in [0.9, 1] um), n+ source and drain at the top corners, a p channel (1e17) between them over a
// lightly doped substrate (1e15); an electrode on the oxide over the channel, ohmic source and drain
// on the pockets' x faces, the body on y_min. Solved at Vg = 1 V, Vd = 0.5 V.
struct Mosfet {
    device::Device device;
    results::BiasPoint point;
};

const Mosfet& mosfet(bool three_d) {
    static auto make = [](bool three) {
        const auto x = uniform(0.0, 2e-4, 21);
        const auto y = uniform(0.0, 1e-4, 21);
        const auto z = three ? uniform(0.0, 0.8e-4, 7) : std::vector<double>{};
        mesh::Mesh m = three ? *mesh::make_tensor_grid(x, y, z) : *mesh::make_tensor_grid(x, y);
        const std::size_t n = m.node_count();
        std::vector<device::RegionId> region(n);
        std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            const auto& p = m.points()[i];
            region[i] = p[1] > 0.92e-4 ? 0 : 1;
            if (region[i] == 0) continue;
            const bool top = p[1] > 0.55e-4;
            acceptors[i] = top ? 1e17 : 1e15;
            if (top && (p[0] < 0.55e-4 || p[0] > 1.45e-4)) donors[i] = 2e19;
        }
        std::vector<mesh::NodeId> gate, source, drain;
        for (const mesh::NodeId v : m.find_boundary("y_max")->nodes) {
            const double px = m.points()[static_cast<std::size_t>(v)][0];
            if (px > 0.6e-4 && px < 1.4e-4) gate.push_back(v);
        }
        for (const char* side : {"x_min", "x_max"}) {
            for (const mesh::NodeId v : m.find_boundary(side)->nodes) {
                const double py = m.points()[static_cast<std::size_t>(v)][1];
                if (py > 0.6e-4 && py < 0.92e-4) (side[2] == 'm' && side[3] == 'i' ? source : drain).push_back(v);
            }
        }
        auto body = m.find_boundary("y_min")->nodes;
        device::Contact g{"gate", device::ContactKind::electrode, std::move(gate)};
        g.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.1};
        auto cells = three ? *mesh::TensorCells::create(m, x, y, z) : *mesh::TensorCells::create(m, x, y);
        auto d = device::Device::create(
            {.mesh = std::move(m),
             .temperature_K = 300.0,
             .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
             .node_region = std::move(region),
             .donors = std::move(donors),
             .acceptors = std::move(acceptors),
             .contacts = {std::move(g),
                          {"source", device::ContactKind::ohmic, std::move(source)},
                          {"drain", device::ContactKind::ohmic, std::move(drain)},
                          {"body", device::ContactKind::ohmic, std::move(body)}},
             .cells = std::move(cells)});
        if (!d) FAIL(d.error().message);
        const std::vector<std::vector<double>> points{{0.0, 0.0, 0.0, 0.0}, {0.5, 0.0, 0.0, 0.0},
                                                      {1.0, 0.0, 0.0, 0.0}, {1.0, 0.0, 0.25, 0.0},
                                                      {1.0, 0.0, 0.5, 0.0}};
        auto sweep = solve::sweep_bias(*d, points);
        if (!sweep) FAIL(sweep.error().message);
        if (sweep->stopped) FAIL(sweep->stopped->message);
        return Mosfet{std::move(*d), sweep->points.back()};
    };
    static const Mosfet two = make(false);
    static const Mosfet three = make(true);
    return three_d ? three : two;
}

Geometry geometry_of(const device::Device& d) { return *Geometry::tensor(d.mesh(), *d.cells()); }

FieldOptions field_options(const device::Device& d, const results::BiasPoint& p, Field f,
                           std::vector<double>& values) {
    values = *field_values(d, p, f);
    FieldOptions o;
    o.title = std::string(field_info(f).title);
    o.scale = *ColorScale::fit(field_info(f).default_scale, values, 1e10);
    return o;
}

}  // namespace

TEST_CASE("render image: WARP renders the pipeline's own test patterns exactly") {
    Renderer& r = warp();
    INFO(r.adapter_name());
    CHECK(r.is_warp());
    // A unit square split in two triangles, flat red, filling a 2D view: the middle pixel is red
    // and a corner white; the image has the requested size.
    const auto line = *mesh::make_tensor_grid(uniform(0.0, 1e-4, 2), uniform(0.0, 1e-4, 2));
    FieldScene s;
    s.dimension = 2;
    TriangleBatch b{"flat", Shading::flat, false, false, {}};
    const Rgba red{1.0f, 0.0f, 0.0f, 1.0f};
    for (const auto& [x, y] : {std::pair{0.0f, 0.0f}, {1.0f, 0.0f}, {1.0f, 1.0f}, {0.0f, 0.0f}, {1.0f, 1.0f}, {0.0f, 1.0f}}) {
        b.vertices.push_back({{x, y, 0.0f}, {0.0f, 0.0f, 1.0f}, 0.0f, red});
    }
    s.triangles.push_back(b);
    s.lo = {0.0f, 0.0f, 0.0f};
    s.hi = {1.0f, 1.0f, 0.0f};
    s.camera.azimuth_deg = 0.0f;
    s.camera.elevation_deg = 0.0f;
    s.axis_triad = false;
    const auto img = r.draw(s, 200, 160);
    if (!img) FAIL(img.error().message);
    CHECK(img->width == 200);
    CHECK(img->height == 160);
    REQUIRE(img->rgba.size() == 200u * 160u * 4u);
    const auto px = [&](int x, int y) { return &img->rgba[static_cast<std::size_t>((y * 200 + x) * 4)]; };
    // The view is the image less the 2D axes' margins; its centre is red.
    const int cx = (66 + (200 - 12)) / 2, cy = (12 + (160 - 48)) / 2;
    CHECK(px(cx, cy)[0] == 255);
    CHECK(px(cx, cy)[1] == 0);
    CHECK(px(cx, cy)[2] == 0);
    CHECK(px(1, 1)[0] == 255);
    CHECK(px(1, 1)[1] == 255);  // white background
    // A field batch: the colormap at position 0.5 everywhere.
    FieldScene f = s;
    f.triangles[0].shading = Shading::field;
    for (SceneVertex& v : f.triangles[0].vertices) v.value = 0.5f;
    f.colorbar = Colorbar{Colormap::viridis, *ColorScale::create(ScaleKind::linear, 0.0, 1.0), "t"};
    const auto fi = r.draw(f, 200, 160);
    if (!fi) FAIL(fi.error().message);
    const Rgb mid = sample(Colormap::viridis, 0.5);
    const std::uint8_t* c = &fi->rgba[static_cast<std::size_t>((cy * 200 + (66 + 70) / 2) * 4)];
    CHECK(std::abs(c[0] - mid.r * 255.0f) <= 2.0f);
    CHECK(std::abs(c[1] - mid.g * 255.0f) <= 2.0f);
    CHECK(std::abs(c[2] - mid.b * 255.0f) <= 2.0f);
    // Sizes outside 16..16384 are refused.
    CHECK(r.draw(s, 8, 100).error().code == base::ErrorCode::invalid_input);
    // PNG round trip.
    fs::create_directories(output);
    if (auto saved = save_png(*fi, output / "roundtrip.png"); !saved) FAIL(saved.error().message);
    const auto back = load_png(output / "roundtrip.png");
    REQUIRE(back.has_value());
    CHECK(compare(*fi, *back, 0)->max_channel == 0);
    for (const std::string& m : r.debug_messages()) UNSCOPED_INFO(m);
    CHECK(r.debug_errors() == 0);
}

TEST_CASE("render image: order-independent transparency does not depend on the order") {
    Renderer& r = warp();
    // Two overlapping translucent squares at different depths, in either order.
    FieldScene s;
    s.dimension = 3;
    const auto square = [](float z, Rgba c) {
        TriangleBatch b{"t", Shading::flat, true, true, {}};
        for (const auto& [x, y] : {std::pair{0.0f, 0.0f}, {1.0f, 0.0f}, {1.0f, 1.0f}, {0.0f, 0.0f}, {1.0f, 1.0f}, {0.0f, 1.0f}}) {
            b.vertices.push_back({{x + z * 0.3f, y, z}, {0.0f, 0.0f, 1.0f}, 0.0f, c});
        }
        return b;
    };
    s.triangles = {square(0.0f, {0.9f, 0.2f, 0.2f, 0.5f}), square(0.5f, {0.2f, 0.3f, 0.9f, 0.5f})};
    s.lo = {0.0f, 0.0f, 0.0f};
    s.hi = {1.15f, 1.0f, 0.5f};
    s.camera.azimuth_deg = 0.0f;
    s.camera.elevation_deg = 0.0f;
    s.axis_triad = false;
    const auto a = r.draw(s, 160, 160);
    std::swap(s.triangles[0], s.triangles[1]);
    const auto b = r.draw(s, 160, 160);
    if (!a) FAIL(a.error().message);
    REQUIRE(b.has_value());
    CHECK(compare(*a, *b, 1)->fraction_over == 0.0);
    for (const std::string& m : r.debug_messages()) UNSCOPED_INFO(m);
    CHECK(r.debug_errors() == 0);
}

TEST_CASE("render image: the views match their WARP references (VG4)") {
    Renderer& r = warp();
    INFO(r.adapter_name());
    constexpr int W = 960, H = 720;
    const Mosfet& m2 = mosfet(false);
    const Mosfet& m3 = mosfet(true);
    const Geometry g2 = geometry_of(m2.device);
    const Geometry g3 = geometry_of(m3.device);
    std::vector<double> values;

    SECTION("2D potential map") {
        FieldOptions o = field_options(m2.device, m2.point, Field::potential, values);
        o.contours = true;
        const auto s = field_view(m2.device, g2, values, o);
        REQUIRE(s.has_value());
        check_image(*r.draw(*s, W, H), "mosfet2d_potential");
    }
    SECTION("2D regions and mesh") {
        const auto s = mesh_view(m2.device, g2);
        REQUIRE(s.has_value());
        check_image(*r.draw(*s, W, H), "mosfet2d_mesh");
    }
    SECTION("band diagram of a 1D cut") {
        // The 1D band diagram plot from the 2D device's column through the channel's middle.
        const auto line = *mesh::make_tensor_grid(uniform(0.0, 0.9e-4, 19));
        results::BandDiagram b;
        for (std::size_t i = 0; i < m2.device.mesh().node_count(); ++i) {
            const auto& p = m2.device.mesh().points()[i];
            if (std::abs(p[0] - 1e-4) > 1e-12 || p[1] > 0.91e-4) continue;
            b.conduction_eV.push_back(m2.point.bands.conduction_eV[i]);
            b.valence_eV.push_back(m2.point.bands.valence_eV[i]);
            b.electron_fermi_eV.push_back(m2.point.bands.electron_fermi_eV[i]);
            b.hole_fermi_eV.push_back(m2.point.bands.hole_fermi_eV[i]);
        }
        REQUIRE(b.conduction_eV.size() == 19);
        const auto s = band_plot(line, b);
        REQUIRE(s.has_value());
        check_image(*r.draw(*s, W, H), "band_diagram");
    }
    SECTION("I-V plot") {
        analysis::Curve id{{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6},
                           {2e-15, 4e-13, 9e-11, 2e-8, 4e-6, 6e-5, 2e-4},
                           {1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14, 1e-14}};
        analysis::Curve ig{{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6},
                           {1e-12, 3e-12, 8e-12, 2e-11, 5e-11, 1e-10, 2e-10}, {}};
        const std::vector<SeriesInput> in{{"drain", id}, {"leakage", ig}};
        const auto s = curve_plot(in, "VG (V)", "I (A)", AxisScale::linear, AxisScale::log10, "Transfer");
        REQUIRE(s.has_value());
        check_image(*r.draw(*s, W, H), "iv_plot");
    }
    SECTION("3D solid") {
        const auto s = solid_view(m3.device, g3);
        REQUIRE(s.has_value());
        check_image(*r.draw(*s, W, H), "mosfet3d_solid");
    }
    SECTION("3D surfaces of n, p and the potential") {
        for (const Field f : {Field::electrons, Field::holes, Field::potential}) {
            const FieldOptions o = field_options(m3.device, m3.point, f, values);
            const auto s = field_view(m3.device, g3, values, o);
            REQUIRE(s.has_value());
            check_image(*r.draw(*s, W, H), std::string("mosfet3d_surface_") + std::string(field_info(f).key));
        }
    }
    SECTION("3D mesh") {
        const auto s = mesh_view(m3.device, g3);
        REQUIRE(s.has_value());
        check_image(*r.draw(*s, W, H), "mosfet3d_mesh");
    }
    SECTION("3D slices") {
        const FieldOptions o = field_options(m3.device, m3.point, Field::electrons, values);
        const std::vector<Plane> planes{{{0.0, 0.0, 1.0}, 0.4e-4}, {{1.0, 0.0, 0.0}, 1.0e-4}};
        const auto s = slice_view(m3.device, g3, values, planes, o);
        REQUIRE(s.has_value());
        check_image(*r.draw(*s, W, H), "mosfet3d_slice");
    }
    SECTION("3D opaque oxide, perspective") {
        FieldOptions o = field_options(m3.device, m3.point, Field::potential, values);
        o.view.transparent_insulators = false;
        o.view.projection = Projection::perspective;
        const auto s = field_view(m3.device, g3, values, o);
        REQUIRE(s.has_value());
        check_image(*r.draw(*s, W, H), "mosfet3d_opaque_perspective");
    }
    for (const std::string& m : r.debug_messages()) UNSCOPED_INFO(m);
    CHECK(r.debug_errors() == 0);
}

TEST_CASE("render image: drawing benchmark (VG7)", "[.benchmark]") {
    // Scene building and drawing of about 1e5 and 1e6 triangles: a 3D tensor box's boundary field
    // view, on WARP and on the GPU when there is one. Times per draw (including readback) and the
    // implied frames per second; the GPU budget is 60 at 1e6 triangles.
    for (const Adapter adapter : {Adapter::hardware, Adapter::warp}) {
        auto r = Renderer::create({.adapter = adapter});
        if (!r) {
            WARN("no renderer for this adapter: " << r.error().message);
            continue;
        }
        for (const int n : {130, 500}) {
            using clock = std::chrono::steady_clock;
            const auto axis = uniform(0.0, 1e-4, n);
            const auto m = *mesh::make_tensor_grid(axis, axis, std::vector<double>{0.0, 1e-5});
            const auto cells = *mesh::TensorCells::create(m, axis, axis, std::vector<double>{0.0, 1e-5});
            const auto t0 = clock::now();
            auto g = *Geometry::tensor(m, cells);
            std::vector<double> v;
            for (const auto& p : m.points()) v.push_back(std::sin(4e4 * p[0]) * std::cos(3e4 * p[1]));
            // A field scene by hand (no device): the boundary faces as one field batch.
            FieldScene s;
            s.dimension = 3;
            TriangleBatch b{"field", Shading::field, true, false, {}};
            const auto scale = *ColorScale::fit(ScaleKind::linear, v);
            for (const Polygon& f : g.boundary_faces()) {
                for (const SurfaceVertex& x : f.vertices) {
                    b.vertices.push_back({{static_cast<float>(x.position[0] * 1e4), static_cast<float>(x.position[1] * 1e4),
                                           static_cast<float>(x.position[2] * 1e4)},
                                          {0.0f, 0.0f, 1.0f},
                                          static_cast<float>(scale.position(interpolate(x.at, v))),
                                          {}});
                }
            }
            const std::size_t triangles = b.vertices.size() / 3;
            s.triangles.push_back(std::move(b));
            s.lo = {0.0f, 0.0f, 0.0f};
            s.hi = {1.0f, 1.0f, 0.1f};
            const double build = std::chrono::duration<double>(clock::now() - t0).count();
            REQUIRE(r->draw(s, 1280, 960).has_value());  // warm up
            const int frames = adapter == Adapter::warp ? 2 : 10;
            const auto t1 = clock::now();
            for (int k = 0; k < frames; ++k) REQUIRE(r->draw(s, 1280, 960).has_value());
            const double per = std::chrono::duration<double>(clock::now() - t1).count() / frames;
            const auto frame = r->frame_time(s, 1280, 960, adapter == Adapter::warp ? 3 : 30);
            REQUIRE(frame.has_value());
            std::printf("render benchmark %s (%s): %zu triangles, scene %.2f s, image %.1f ms (%.1f per second), "
                        "redraw on resident geometry %.2f ms (%.0f frames per second)\n",
                        r->is_warp() ? "WARP" : "GPU", r->adapter_name().c_str(), triangles, build, 1e3 * per,
                        1.0 / per, *frame, 1e3 / *frame);
        }
    }
}
