// Drawing (Unit 25a, VG4, VG7): WARP renders of the views compared with committed references, the
// image pipeline's own checks, the debug layer, and the benchmark.
//
// The references are tests/render/references/*.png, rendered on WARP. A render matches its
// reference when no more than 0.5% of its pixels differ by more than 24 in any channel: rounding in
// the rasterizer and glyph rendering may differ between Windows builds, a wrong view differs over
// whole regions. To write new references, run with NITCAD_UPDATE_REFERENCES=1, look at every image,
// and commit them with the reason.
#include <catch2/catch_test_macros.hpp>

#include <windows.h>
#include <dxgi1_4.h>
#include <wrl/client.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <optional>
#include <utility>
#include <string>
#include <vector>

#include "NiTCAD/analysis/curve.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/render/compare.hpp"
#include "NiTCAD/render/edges.hpp"
#include "NiTCAD/render/image.hpp"
#include "NiTCAD/render/scene.hpp"
#include "NiTCAD/render/streamlines.hpp"
#include "NiTCAD/render/vectors.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "image_check.hpp"

using namespace NiTCAD;
using namespace NiTCAD::render;
namespace fs = std::filesystem;

namespace {

using render_test::check_image;
using render_test::Mosfet;
using render_test::mosfet;
using render_test::output;
using render_test::references;
using render_test::warp;

std::vector<double> uniform(double a, double b, int n) {
    std::vector<double> v;
    for (int k = 0; k < n; ++k) v.push_back(k == n - 1 ? b : a + (b - a) * k / (n - 1));
    return v;
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

// A vector magnitude with nitcad plot's floor: values below ten decades under the largest raised to it.
std::vector<double> magnitudes_floored(const device::Device& d, const Geometry& g, const results::BiasPoint& p,
                                       Field f) {
    std::vector<double> v = *field_values(d, g, p, f);
    double top = 0.0;
    for (const double x : v) if (std::isfinite(x)) top = std::max(top, x);
    const double floor = std::pow(10.0, std::ceil(std::log10(top)) - 10.0);
    for (double& x : v) if (std::isfinite(x) && x < floor) x = floor;
    return v;
}

// Streamlines of a vector field from a contact: equal current for J, evenly for E.
std::vector<Streamline> lines_of(const Mosfet& m, const Geometry& g, VectorField f, std::size_t contact,
                                 std::size_t count) {
    const auto nodal = streamline_field(m.device, g, *element_vectors(m.device, g, m.point, f), f);
    const Locator loc(g);
    const std::vector<double> currents =
        f == VectorField::electric ? std::vector<double>{} : contact_node_currents(m.device, m.point, contact);
    std::vector<Streamline> out;
    for (const Point& seed : contact_seeds(m.device, loc, contact, count, currents)) {
        out.push_back(trace_streamline(loc, nodal, seed));
    }
    return out;
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
        auto s = field_view(m2.device, g2, values, o);
        REQUIRE(s.has_value());
        CHECK(s->look == Look::classic);
        const auto classic = r.draw(*s, W, H);
        REQUIRE(classic.has_value());
        check_image(*classic, "mosfet2d_potential");
        // Unit 26f (KG1): 2D ignores the look, to the bit.
        set_look(*s, Look::enhanced);
        const auto enhanced = r.draw(*s, W, H);
        REQUIRE(enhanced.has_value());
        CHECK(enhanced->rgba == classic->rgba);
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
    SECTION("band diagrams from source to drain along cut lines (Unit 26b, XG4)") {
        // Under the oxide, through the source pocket, the channel and the drain pocket.
        const auto c2 = cut_line(g2, {0.0, 0.87e-4, 0.0}, {2e-4, 0.87e-4, 0.0});
        REQUIRE(c2.has_value());
        const auto s2 = cut_band_plot(m2.device.mesh(), *c2, m2.point.bands);
        REQUIRE(s2.has_value());
        check_image(*r.draw(*s2, W, H), "mosfet2d_cut_bands");
        const auto c3 = cut_line(g3, {0.0, 0.87e-4, 0.4e-4}, {2e-4, 0.87e-4, 0.4e-4});
        REQUIRE(c3.has_value());
        const auto s3 = cut_band_plot(m3.device.mesh(), *c3, m3.point.bands);
        REQUIRE(s3.has_value());
        check_image(*r.draw(*s3, W, H), "mosfet3d_cut_bands");
        // The line drawn through the 3D device.
        auto solid = solid_view(m3.device, g3);
        REQUIRE(solid.has_value());
        add_cut_line(*solid, *c3);
        check_image(*r.draw(*solid, W, H), "mosfet3d_cut_line");
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
    // Unit 26f (KG4): the enhanced look on the solid and the potential surface, each its own reference
    // and each unlike its classic render.
    SECTION("3D solid and potential surface, enhanced") {
        const FieldOptions o = field_options(m3.device, m3.point, Field::potential, values);
        auto solid = solid_view(m3.device, g3);
        auto surface = field_view(m3.device, g3, values, o);
        REQUIRE(solid.has_value());
        REQUIRE(surface.has_value());
        for (auto* s : {&*solid, &*surface}) {
            const auto classic = r.draw(*s, W, H);
            REQUIRE(classic.has_value());
            set_look(*s, Look::enhanced);
            const auto enhanced = r.draw(*s, W, H);
            REQUIRE(enhanced.has_value());
            check_image(*enhanced, s == &*solid ? "mosfet3d_solid_enhanced" : "mosfet3d_surface_potential_enhanced");
            const auto d = compare(*enhanced, *classic, 24);
            REQUIRE(d.has_value());
            CHECK(d->fraction_over > 0.05);
        }
    }
    // Unit 26h (JG7): the enhanced solid at a named view (Right), its axis box and view cube turned with it.
    SECTION("3D solid at the Right view, enhanced") {
        auto solid = solid_view(m3.device, g3);
        REQUIRE(solid.has_value());
        set_look(*solid, Look::enhanced);
        solid->camera.azimuth_deg = 90.0f;
        solid->camera.elevation_deg = 0.0f;
        check_image(*r.draw(*solid, W, H), "mosfet3d_solid_right_enhanced");
    }
    // Unit 26g (MG6): the enhanced 3D field with mesh edges (the mesh view's) and contours.
    SECTION("3D potential with mesh edges and contours, enhanced") {
        FieldOptions o = field_options(m3.device, m3.point, Field::potential, values);
        o.contours = true;
        auto s = field_view(m3.device, g3, values, o);
        const auto mesh = mesh_view(m3.device, g3);
        REQUIRE(s.has_value());
        REQUIRE(mesh.has_value());
        for (const LineBatch& l : mesh->lines) {
            if (l.name == "mesh") s->lines.push_back(l);
        }
        set_look(*s, Look::enhanced);
        check_image(*r.draw(*s, W, H), "mosfet3d_field_mesh_contours_enhanced");
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
    // Unit 25b (BG4): |E| and |J| with their lines, isosurfaces, 3D current lines.
    SECTION("2D |E| with field lines from the gate") {
        values = magnitudes_floored(m2.device, g2, m2.point, Field::electric_field);
        FieldOptions o;
        o.title = std::string(field_info(Field::electric_field).title);
        o.scale = *ColorScale::fit(ScaleKind::log10, values, 1e10);
        auto s = field_view(m2.device, g2, values, o);
        REQUIRE(s.has_value());
        add_streamlines(*s, lines_of(m2, g2, VectorField::electric, 0, 16), Rgba{0.05f, 0.05f, 0.05f, 0.85f});
        check_image(*r.draw(*s, W, H), "mosfet2d_efield_lines");
    }
    SECTION("2D |J| with current lines from the source") {
        values = magnitudes_floored(m2.device, g2, m2.point, Field::current_density);
        FieldOptions o;
        o.title = std::string(field_info(Field::current_density).title);
        o.scale = *ColorScale::fit(ScaleKind::log10, values, 1e10);
        auto s = field_view(m2.device, g2, values, o);
        REQUIRE(s.has_value());
        add_streamlines(*s, lines_of(m2, g2, VectorField::current, 1, 16), Rgba{0.05f, 0.05f, 0.05f, 0.85f});
        check_image(*r.draw(*s, W, H), "mosfet2d_current_lines");
    }
    SECTION("3D isosurfaces of n") {
        const FieldOptions o = field_options(m3.device, m3.point, Field::electrons, values);
        const std::vector<double> levels{1e14, 1e18};
        const auto s = iso_view(m3.device, g3, values, levels, o);
        REQUIRE(s.has_value());
        check_image(*r.draw(*s, W, H), "mosfet3d_iso_n");
    }
    SECTION("3D current lines colored by |J|") {
        values = magnitudes_floored(m3.device, g3, m3.point, Field::current_density);
        const ColorScale scale = *ColorScale::fit(ScaleKind::log10, values, 1e10);
        auto s = outline_view(m3.device, g3);
        REQUIRE(s.has_value());
        s->colorbar = Colorbar{Colormap::turbo, scale, std::string(field_info(Field::current_density).title)};
        add_streamlines(*s, lines_of(m3, g3, VectorField::current, 1, 24), *s->colorbar);
        check_image(*r.draw(*s, W, H), "mosfet3d_current_lines");
    }
    // Unit 25c (CG5): arrows, the clipping box, the exploded view, a filmstrip.
    SECTION("2D |J| with current arrows") {
        values = magnitudes_floored(m2.device, g2, m2.point, Field::current_density);
        FieldOptions o;
        o.title = std::string(field_info(Field::current_density).title);
        o.scale = *ColorScale::fit(ScaleKind::log10, values, 1e10);
        auto s = field_view(m2.device, g2, values, o);
        REQUIRE(s.has_value());
        const auto field = streamline_field(m2.device, g2, *element_vectors(m2.device, g2, m2.point, VectorField::current),
                                            VectorField::current);
        std::vector<std::uint8_t> insulator;  // no current arrows in the oxide
        for (std::size_t i = 0; i < m2.device.mesh().node_count(); ++i) {
            insulator.push_back(m2.device.is_insulator(static_cast<NodeId>(i)) ? 1 : 0);
        }
        double h = 0.0;
        const auto arrows = grid_arrows(Locator(g2), field, 250, h, insulator);
        add_arrows(*s, arrows, 0.8 * h, Rgba{0.05f, 0.05f, 0.05f, 1.0f});
        check_image(*r.draw(*s, W, H), "mosfet2d_current_arrows");
    }
    SECTION("3D E arrows on the surface") {
        values = magnitudes_floored(m3.device, g3, m3.point, Field::electric_field);
        auto s = solid_view(m3.device, g3);
        REQUIRE(s.has_value());
        s->colorbar = Colorbar{Colormap::turbo, *ColorScale::fit(ScaleKind::log10, values, 1e10),
                               std::string(field_info(Field::electric_field).title)};
        const auto field = streamline_field(m3.device, g3, *element_vectors(m3.device, g3, m3.point, VectorField::electric),
                                            VectorField::electric);
        double h = 0.0;
        const auto arrows = surface_arrows(Locator(g3), field, exposed_boundary(m3.device, g3), 300, h);
        add_arrows(*s, arrows, 0.8 * h, *s->colorbar);
        check_image(*r.draw(*s, W, H), "mosfet3d_efield_arrows");
    }
    SECTION("3D quarter cut of n") {
        FieldOptions o = field_options(m3.device, m3.point, Field::electrons, values);
        o.view.clip = Box{{1.0e-4, 0.5e-4, 0.4e-4}, {3e-4, 3e-4, 3e-4}};
        const auto s = field_view(m3.device, g3, values, o);
        REQUIRE(s.has_value());
        check_image(*r.draw(*s, W, H), "mosfet3d_clip_n");
    }
    SECTION("3D exploded") {
        const auto s = exploded_view(m3.device, g3, 1, default_gap(m3.device, 1));
        REQUIRE(s.has_value());
        check_image(*r.draw(*s, W, H), "mosfet3d_exploded");
    }
    SECTION("2D filmstrip of n over the gate sweep") {
        std::vector<FilmFrame> frames;
        for (const results::BiasPoint& p : m2.sweep) {
            char label[48];
            std::snprintf(label, sizeof label, "VG = %.2f V, VD = %.2f V", p.bias_V[0], p.bias_V[2]);
            frames.push_back({label, *field_values(m2.device, p, Field::electrons)});
        }
        FieldOptions o;
        o.title = std::string(field_info(Field::electrons).title);
        const auto scenes = filmstrip_scenes(m2.device, g2, frames, ScaleKind::log10, o);
        REQUIRE(scenes.has_value());
        std::vector<Image> tiles;
        for (const FieldScene& sc : *scenes) tiles.push_back(*r.draw(sc, 480, 360));
        check_image(*tile_images(tiles, 3), "mosfet2d_filmstrip_n");
    }
    SECTION("2D difference map: the drain at 0.5 V less the drain at 0 V (Unit 26b, XG5)") {
        const auto d = difference_field(m2.device, m2.device, g2, m2.sweep[2], m2.sweep.back(), Field::potential);
        REQUIRE(d.has_value());
        FieldOptions o;
        o.title = d->title;
        o.map = Colormap::coolwarm;
        o.scale = difference_scale(d->values);
        o.contours = true;
        const auto s = field_view(m2.device, g2, d->values, o);
        REQUIRE(s.has_value());
        check_image(*r.draw(*s, W, H), "mosfet2d_difference_potential");
    }
    SECTION("3D opaque oxide, perspective") {
        FieldOptions o = field_options(m3.device, m3.point, Field::potential, values);
        o.view.transparent_insulators = false;
        o.view.projection = Projection::perspective;
        const auto s = field_view(m3.device, g3, values, o);
        REQUIRE(s.has_value());
        check_image(*r.draw(*s, W, H), "mosfet3d_opaque_perspective");
    }
    // Unit 26c (YG2): the oxide hidden, its scale fitted over the silicon.
    SECTION("2D potential with the oxide hidden") {
        values = *field_values(m2.device, m2.point, Field::potential);
        FieldOptions o;
        o.title = std::string(field_info(Field::potential).title);
        o.view.hidden_regions = {0};
        o.view.mesh_edges = true;
        o.contours = true;
        o.scale = *ColorScale::fit(ScaleKind::linear, shown_values(m2.device, values, o.view.hidden_regions), 1e10);
        const auto s = field_view(m2.device, g2, values, o);
        REQUIRE(s.has_value());
        check_image(*r.draw(*s, W, H), "mosfet2d_hidden_oxide");
    }
    SECTION("3D n with the oxide hidden") {
        FieldOptions o = field_options(m3.device, m3.point, Field::electrons, values);
        o.view.hidden_regions = {0};
        o.scale = *ColorScale::fit(ScaleKind::log10, shown_values(m3.device, values, o.view.hidden_regions), 1e10);
        const auto s = field_view(m3.device, g3, values, o);
        REQUIRE(s.has_value());
        check_image(*r.draw(*s, W, H), "mosfet3d_hidden_oxide");
    }
    for (const std::string& m : r.debug_messages()) UNSCOPED_INFO(m);
    CHECK(r.debug_errors() == 0);
}

namespace {

// Mean luminance (0-255, Rec. 709 weights) of the pixels within `radius` of an image point.
double mean_luminance(const Image& im, float cx, float cy, int radius) {
    double sum = 0.0;
    int count = 0;
    for (int y = static_cast<int>(cy) - radius; y <= static_cast<int>(cy) + radius; ++y) {
        for (int x = static_cast<int>(cx) - radius; x <= static_cast<int>(cx) + radius; ++x) {
            if (x < 0 || y < 0 || x >= im.width || y >= im.height) continue;
            const std::size_t i = 4 * (static_cast<std::size_t>(y) * static_cast<std::size_t>(im.width) + static_cast<std::size_t>(x));
            sum += 0.2126 * im.rgba[i] + 0.7152 * im.rgba[i + 1] + 0.0722 * im.rgba[i + 2];
            ++count;
        }
    }
    return count ? sum / count : 0.0;
}

}  // namespace

// Unit 26g (MG1, MG3): classic ignores the finishes to the bit; enhanced draws the electrodes as metal.
TEST_CASE("render image: finishes in classic and enhanced (Unit 26g, MG1, MG3)") {
    Renderer& r = warp();
    constexpr int W = 960, H = 720;
    const Mosfet& m3 = mosfet(true);
    const Geometry g3 = geometry_of(m3.device);
    auto solid = solid_view(m3.device, g3);
    REQUIRE(solid.has_value());
    REQUIRE(std::ranges::any_of(solid->triangles, [](const TriangleBatch& b) { return b.finish == Finish::metal; }));
    REQUIRE(std::ranges::any_of(solid->triangles, [](const TriangleBatch& b) { return b.finish == Finish::glass; }));
    FieldScene matte = *solid;
    for (TriangleBatch& b : matte.triangles) b.finish = Finish::matte;
    // MG1: a classic scene with metal and glass batches draws exactly as with matte ones.
    const auto classic = r.draw(*solid, W, H);
    const auto classic_matte = r.draw(matte, W, H);
    REQUIRE(classic.has_value());
    REQUIRE(classic_matte.has_value());
    CHECK(classic->rgba == classic_matte->rgba);
    // MG3 (metal): the gate electrode's top face (y = 1.06 um, the slab 3% of the 2 um extent on the
    // oxide) is lighter than classic's and lighter than its front face (z = 0.8 um).
    set_look(*solid, Look::enhanced);
    const auto enhanced = r.draw(*solid, W, H);
    REQUIRE(enhanced.has_value());
    const ScreenPoint top = screen_position(*solid, W, H, {1.0f, 1.06f, 0.4f});
    const ScreenPoint front = screen_position(*solid, W, H, {1.0f, 1.03f, 0.8f});
    const double top_enhanced = mean_luminance(*enhanced, top.x, top.y, 4);
    const double top_classic = mean_luminance(*classic, top.x, top.y, 4);
    const double front_enhanced = mean_luminance(*enhanced, front.x, front.y, 3);
    CAPTURE(top.x, top.y, front.x, front.y, top_enhanced, top_classic, front_enhanced);
    CHECK(top_enhanced > top_classic);
    CHECK(top_enhanced > front_enhanced);
    // The same scene with matte electrodes in enhanced is darker on top: the finish does it.
    set_look(matte, Look::enhanced);
    const auto enhanced_matte = r.draw(matte, W, H);
    REQUIRE(enhanced_matte.has_value());
    CHECK(top_enhanced > mean_luminance(*enhanced_matte, top.x, top.y, 4));
    // MG3 (glass, as ruled 2026-11-18): within 3 px of the oxide's outline (its top edges, y = 1 um, away
    // from the gate) the glass is lighter than the same pixels with a matte oxide.
    for (const Float3 p : {Float3{0.3f, 1.0f, 0.0f}, Float3{1.7f, 1.0f, 0.0f}, Float3{0.3f, 1.0f, 0.8f},
                           Float3{1.7f, 1.0f, 0.8f}}) {
        const ScreenPoint q = screen_position(*solid, W, H, p);
        const double glass = mean_luminance(*enhanced, q.x, q.y, 3);
        const double plain = mean_luminance(*enhanced_matte, q.x, q.y, 3);
        CAPTURE(p[0], p[2], glass, plain);
        CHECK(glass > plain + 10.0);
    }
    for (const std::string& m : r.debug_messages()) UNSCOPED_INFO(m);
    CHECK(r.debug_errors() == 0);
}

// Unit 26g (MG5): the ground shadow under the device in the enhanced look; none looking up from below,
// none in classic, none when the scene says so (clipped and exploded views, scene_test).
TEST_CASE("render image: the ground shadow (Unit 26g, MG5)") {
    Renderer& r = warp();
    constexpr int W = 960, H = 720;
    const Mosfet& m3 = mosfet(true);
    const Geometry g3 = geometry_of(m3.device);
    auto s = solid_view(m3.device, g3);
    REQUIRE(s.has_value());
    REQUIRE(s->ground_shadow);
    FieldScene none = *s;
    none.ground_shadow = false;
    // Classic: never drawn.
    CHECK(r.draw(*s, W, H)->rgba == r.draw(none, W, H)->rgba);
    set_look(*s, Look::enhanced);
    set_look(none, Look::enhanced);
    const auto with = r.draw(*s, W, H);
    const auto without = r.draw(none, W, H);
    REQUIRE(with.has_value());
    REQUIRE(without.has_value());
    // Under the front bottom edge's middle: darker than the gradient alone.
    const float cx = 0.5f * (s->lo[0] + s->hi[0]);
    const ScreenPoint q = screen_position(*s, W, H, {cx, s->lo[1], s->hi[2]});
    const double a = mean_luminance(*with, q.x, q.y + 10.0f, 2);
    const double b = mean_luminance(*without, q.x, q.y + 10.0f, 2);
    CAPTURE(q.x, q.y, a, b);
    CHECK(a < b - 8.0);
    // Far from the device (the view's top left): the gradient only.
    CHECK(mean_luminance(*with, 40.0f, 60.0f, 2) == mean_luminance(*without, 40.0f, 60.0f, 2));
    // Looking up from below (elevation -10 degrees): not drawn.
    s->camera.elevation_deg = -10.0f;
    none.camera.elevation_deg = -10.0f;
    CHECK(r.draw(*s, W, H)->rgba == r.draw(none, W, H)->rgba);
    // Just above the limit (-4 degrees): still drawn.
    s->camera.elevation_deg = -4.0f;
    none.camera.elevation_deg = -4.0f;
    CHECK(r.draw(*s, W, H)->rgba != r.draw(none, W, H)->rgba);
    for (const std::string& m : r.debug_messages()) UNSCOPED_INFO(m);
    CHECK(r.debug_errors() == 0);
}

// Unit 26g (MD7): in enhanced the colourbar's frame is lighter and its ticks show inside the bar too.
TEST_CASE("render image: the enhanced colourbar (Unit 26g, MD7)") {
    Renderer& r = warp();
    constexpr int W = 960, H = 720;
    const Mosfet& m3 = mosfet(true);
    const Geometry g3 = geometry_of(m3.device);
    std::vector<double> values;
    const FieldOptions o = field_options(m3.device, m3.point, Field::potential, values);
    auto s = field_view(m3.device, g3, values, o);
    REQUIRE(s.has_value());
    const auto classic = r.draw(*s, W, H);
    set_look(*s, Look::enhanced);
    const auto enhanced = r.draw(*s, W, H);
    REQUIRE(classic.has_value());
    REQUIRE(enhanced.has_value());
    // The bar as draw() lays it out: 16 px wide from x = W - 130 + 14, from 64 px to H - 24.
    const float left = W - 130.0f + 14.0f, right = left + 16.0f, top = 64.0f, bottom = H - 24.0f;
    const float mid = 0.5f * (top + bottom);
    CHECK(mean_luminance(*enhanced, left, mid, 0) > mean_luminance(*classic, left, mid, 0) + 10.0);
    // Inside the right edge at each major tick: the frame's grey in enhanced, the colormap in classic.
    int inner = 0;
    for (const Tick& t : s->colorbar->scale.ticks()) {
        if (!t.major) continue;
        const float y = std::round(bottom - (bottom - top) * static_cast<float>(t.position));
        if (y < top + 3.0f || y > bottom - 3.0f) continue;
        const std::size_t i = 4 * (static_cast<std::size_t>(y) * W + static_cast<std::size_t>(right - 3.0f));
        int most = 0;
        for (int c = 0; c < 3; ++c) most = std::max(most, std::abs(int(enhanced->rgba[i + c]) - int(classic->rgba[i + c])));
        CAPTURE(t.position, y, most);
        CHECK(most > 20);
        ++inner;
    }
    CHECK(inner >= 2);
}

namespace {

// Unit 26i: a hidden window's swap chain on the WARP renderer, or nothing where none can be made.
struct HiddenTarget {
    HWND hwnd = nullptr;
    std::optional<WindowTarget> target;
    HiddenTarget(Renderer& r, int w, int h) {
        hwnd = CreateWindowExW(0, L"STATIC", L"", WS_OVERLAPPEDWINDOW, 0, 0, w, h, nullptr, nullptr, nullptr, nullptr);
        if (!hwnd) return;
        if (auto t = r.window_target(hwnd, w, h)) target.emplace(std::move(*t));
    }
    HiddenTarget(const HiddenTarget&) = delete;
    HiddenTarget& operator=(const HiddenTarget&) = delete;
    ~HiddenTarget() {
        target.reset();
        if (hwnd) DestroyWindow(hwnd);
    }
};

// The window's frame after the last present, against draw() of the same scene: bit for bit.
void check_same_as_draw(Renderer& r, const FieldScene& s, int w, int h) {
    const auto shown = r.last_image();
    REQUIRE(shown.has_value());
    const auto drawn = r.draw(s, w, h);
    REQUIRE(drawn.has_value());
    CHECK(shown->rgba == drawn->rgba);
}

}  // namespace

// Unit 26i (NG1): after a camera change, present_camera draws the kept buffers again, and the frame is
// the one a full draw makes, byte for byte, in each look and kind of 3D scene and at 8 camera poses.
TEST_CASE("render window: camera-only presents equal full draws (Unit 26i, NG1)") {
    Renderer& r = warp();
    constexpr int W = 480, H = 360;
    HiddenTarget t(r, W, H);
    if (!t.target) SKIP("no swap chain here");
    const Mosfet& m3 = mosfet(true);
    const Geometry g3 = geometry_of(m3.device);
    std::vector<double> values;
    std::vector<std::pair<std::string, FieldScene>> scenes;
    {
        auto solid = solid_view(m3.device, g3);
        REQUIRE(solid.has_value());
        scenes.emplace_back("classic solid", *solid);
        set_look(*solid, Look::enhanced);
        scenes.emplace_back("enhanced solid (metal and glass, shadow, axis box)", *solid);
        FieldOptions o = field_options(m3.device, m3.point, Field::potential, values);
        o.contours = true;
        auto field = field_view(m3.device, g3, values, o);
        const auto mesh = mesh_view(m3.device, g3);
        REQUIRE(field.has_value());
        REQUIRE(mesh.has_value());
        for (const LineBatch& l : mesh->lines) {
            if (l.name == "mesh") field->lines.push_back(l);
        }
        set_look(*field, Look::enhanced);
        scenes.emplace_back("enhanced potential with edges and contours", *field);
        FieldOptions c = field_options(m3.device, m3.point, Field::electrons, values);
        c.view.clip = Box{{1.0e-4, 0.5e-4, 0.4e-4}, {3e-4, 3e-4, 3e-4}};
        auto clipped = field_view(m3.device, g3, values, c);
        REQUIRE(clipped.has_value());
        set_look(*clipped, Look::enhanced);
        scenes.emplace_back("enhanced clip", *clipped);
        auto exploded = exploded_view(m3.device, g3, 1, default_gap(m3.device, 1));
        REQUIRE(exploded.has_value());
        set_look(*exploded, Look::enhanced);
        scenes.emplace_back("enhanced exploded", *exploded);
    }
    struct Pose {
        float azimuth, elevation, zoom, pan_x, pan_y;
        Projection projection;
    };
    const Pose poses[] = {{10.0f, 20.0f, 1.0f, 0.0f, 0.0f, Projection::orthographic},
                          {90.0f, 0.0f, 1.0f, 0.0f, 0.0f, Projection::orthographic},
                          {200.0f, 45.0f, 1.5f, 0.1f, -0.1f, Projection::orthographic},
                          {300.0f, -20.0f, 1.0f, 0.0f, 0.0f, Projection::orthographic},  // below: no shadow
                          {35.0f, 25.0f, 1.0f, 0.0f, 0.0f, Projection::perspective},
                          {120.0f, 60.0f, 2.0f, -0.2f, 0.2f, Projection::perspective},
                          {250.0f, -30.0f, 0.8f, 0.0f, 0.0f, Projection::perspective},
                          {0.0f, 89.0f, 1.0f, 0.0f, 0.0f, Projection::orthographic}};
    for (auto& [name, s] : scenes) {
        CAPTURE(name);
        REQUIRE(r.present(s, *t.target, W, H).has_value());
        CHECK_FALSE(r.last_present_reused());
        check_same_as_draw(r, s, W, H);
        for (const Pose& p : poses) {
            CAPTURE(p.azimuth, p.elevation, p.zoom);
            s.camera.azimuth_deg = p.azimuth;
            s.camera.elevation_deg = p.elevation;
            s.camera.zoom = p.zoom;
            s.camera.pan = {p.pan_x, p.pan_y};
            s.camera.projection = p.projection;
            REQUIRE(r.present_camera(s, *t.target, W, H).has_value());
            CHECK(r.last_present_reused());
            check_same_as_draw(r, s, W, H);
        }
    }
}

// Unit 26i (NG2): a scene edit followed by a full present, then a camera move, draws the new scene; a
// target holding another scene's buffers, or buffers of other sizes, takes the full path; a resize keeps
// them; a plot in between drops them.
TEST_CASE("render window: edits and other scenes take the full path (Unit 26i, NG2)") {
    Renderer& r = warp();
    constexpr int W = 480, H = 360;
    HiddenTarget t(r, W, H);
    if (!t.target) SKIP("no swap chain here");
    const Mosfet& m3 = mosfet(true);
    const Geometry g3 = geometry_of(m3.device);
    std::vector<double> values;
    const FieldOptions o = field_options(m3.device, m3.point, Field::potential, values);
    auto s = field_view(m3.device, g3, values, o);
    REQUIRE(s.has_value());
    set_look(*s, Look::enhanced);
    const auto move = [&](FieldScene& x, float by) {
        x.camera.azimuth_deg += by;
        REQUIRE(r.present_camera(x, *t.target, W, H).has_value());
    };
    REQUIRE(r.present(*s, *t.target, W, H).has_value());
    move(*s, 15.0f);
    CHECK(r.last_present_reused());

    // Each edit, then the full present the window makes after it, then a camera move.
    const auto edited = [&](const char* what, const std::function<void(FieldScene&)>& edit) {
        CAPTURE(what);
        edit(*s);
        REQUIRE(r.present(*s, *t.target, W, H).has_value());
        CHECK_FALSE(r.last_present_reused());
        move(*s, 10.0f);
        CHECK(r.last_present_reused());
        check_same_as_draw(r, *s, W, H);
    };
    edited("colormap", [](FieldScene& x) { x.colorbar->map = Colormap::viridis; });
    edited("classic look", [](FieldScene& x) { set_look(x, Look::classic); });
    edited("enhanced look", [](FieldScene& x) { set_look(x, Look::enhanced); });
    edited("a finish", [](FieldScene& x) { x.triangles.front().finish = Finish::metal; });
    edited("transparency", [](FieldScene& x) {
        x.triangles.front().transparent = true;
        for (SceneVertex& v : x.triangles.front().vertices) v.color.a = 0.5f;
    });
    edited("no axis box", [](FieldScene& x) { x.show_axis_box = false; });
    edited("no shadow", [](FieldScene& x) { x.ground_shadow = false; });
    edited("field values", [](FieldScene& x) {
        for (TriangleBatch& b : x.triangles) {
            for (SceneVertex& v : b.vertices) v.value = 1.0f - v.value;
        }
    });

    // Another scene object: the full path, even right after a present of the first.
    auto other = solid_view(m3.device, g3);
    REQUIRE(other.has_value());
    REQUIRE(r.present(*s, *t.target, W, H).has_value());
    REQUIRE(r.present_camera(*other, *t.target, W, H).has_value());
    CHECK_FALSE(r.last_present_reused());
    check_same_as_draw(r, *other, W, H);
    // The same object with other batch sizes (a scene rebuilt in place): the full path.
    FieldScene grown = *s;
    REQUIRE(r.present(grown, *t.target, W, H).has_value());
    LineBatch extra;
    extra.name = "extra";
    extra.color = {0.2f, 0.2f, 0.2f, 1.0f};
    extra.points = {{0.0f, 0.0f, 0.0f}, {1.0f, 1.0f, 1.0f}};
    grown.lines.push_back(extra);
    REQUIRE(r.present_camera(grown, *t.target, W, H).has_value());
    CHECK_FALSE(r.last_present_reused());
    check_same_as_draw(r, grown, W, H);
    // A resize: the kept buffers still serve, at the new size.
    REQUIRE(r.present(*s, *t.target, W, H).has_value());
    s->camera.azimuth_deg += 5.0f;
    REQUIRE(r.present_camera(*s, *t.target, W + 40, H - 20).has_value());
    CHECK(r.last_present_reused());
    check_same_as_draw(r, *s, W + 40, H - 20);
    // A plot in between drops them: the next camera-only present is a full one.
    const analysis::Curve curve{{0.0, 1.0, 2.0}, {1.0, 4.0, 9.0}, {}};
    const std::vector<SeriesInput> series{{"y", curve}};
    auto plot = curve_plot(series, "x", "y");
    REQUIRE(plot.has_value());
    REQUIRE(r.present(*plot, *t.target, W, H).has_value());
    REQUIRE(r.present_camera(*s, *t.target, W, H).has_value());
    CHECK_FALSE(r.last_present_reused());
    check_same_as_draw(r, *s, W, H);
}


// Unit 26i (NG5): the kept buffers are replaced, not added to. On WARP the GPU's memory is the process's:
// after 100 field changes, each presented in full and then turned, the private bytes are within 1% of
// those after the first.
TEST_CASE("render window: kept buffers are replaced, not added to (Unit 26i, NG5)") {
    Renderer& r = warp();
    constexpr int W = 480, H = 360;
    HiddenTarget t(r, W, H);
    if (!t.target) SKIP("no swap chain here");
    const Mosfet& m3 = mosfet(true);
    const Geometry g3 = geometry_of(m3.device);
    std::vector<double> values;
    const FieldOptions o = field_options(m3.device, m3.point, Field::potential, values);
    auto s = field_view(m3.device, g3, values, o);
    REQUIRE(s.has_value());
    set_look(*s, Look::enhanced);
    const auto change = [&](int k) {
        for (TriangleBatch& b : s->triangles) {
            for (SceneVertex& v : b.vertices) v.value = std::fmod(v.value + 0.01f * static_cast<float>(k), 1.0f);
        }
        REQUIRE(r.present(*s, *t.target, W, H).has_value());
        s->camera.azimuth_deg += 3.0f;
        REQUIRE(r.present_camera(*s, *t.target, W, H).has_value());
        REQUIRE(r.last_present_reused());
    };
    // NG5 is about GPU memory, so it is measured there: the process's D3D12 memory on the WARP adapter (DXGI's
    // QueryVideoMemoryInfo, both segment groups), which counts the committed buffers and nothing else. The
    // process's private bytes, measured until 2026-11-20, also hold WARP's and the CRT's heaps: after any
    // warm-up they still drifted 1-2.5 MB in 100 changes (ratios to 1.017 locally and 1.013 on the CI runner),
    // while the GPU memory is the same to the byte after 1, 100 and 200 changes (19,202,048 bytes here).
    const auto gpu_bytes = [] {
        Microsoft::WRL::ComPtr<IDXGIFactory4> f;
        REQUIRE(SUCCEEDED(CreateDXGIFactory2(0, IID_PPV_ARGS(&f))));
        Microsoft::WRL::ComPtr<IDXGIAdapter3> a;
        REQUIRE(SUCCEEDED(f->EnumWarpAdapter(IID_PPV_ARGS(&a))));
        double sum = 0.0;
        for (const auto group : {DXGI_MEMORY_SEGMENT_GROUP_LOCAL, DXGI_MEMORY_SEGMENT_GROUP_NON_LOCAL}) {
            DXGI_QUERY_VIDEO_MEMORY_INFO info{};
            REQUIRE(SUCCEEDED(a->QueryVideoMemoryInfo(0, group, &info)));
            sum += static_cast<double>(info.CurrentUsage);
        }
        return sum;
    };
    change(0);
    const double one = gpu_bytes();
    REQUIRE(one > 0.0);
    for (int k = 1; k <= 100; ++k) change(k);
    const double hundred = gpu_bytes();
    CAPTURE(one, hundred, hundred / one);
    CHECK(hundred <= 1.01 * one);
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
            for (const render::Polygon& f : g.boundary_faces()) {
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

TEST_CASE("render image: derived-field benchmark (Unit 25b, BG7)", "[.benchmark]") {
    // |E| and |J| reconstruction (elements and nodes), one isosurface, the streamline field and 100
    // streamlines from a plane, on the 3D MOSFET and on a 3D tensor device of about 1e6 tetrahedra
    // (smooth synthetic potential and currents). Release times in s.
    using clock = std::chrono::steady_clock;
    const auto seconds = [](clock::time_point t) { return std::chrono::duration<double>(clock::now() - t).count(); };
    const auto run = [&](const char* name, const device::Device& d, const results::BiasPoint& p) {
        const Geometry g = geometry_of(d);
        auto t = clock::now();
        const auto e = *element_vectors(d, g, p, VectorField::electric);
        const auto en = nodal_vectors(g, e);
        const double t_e = seconds(t);
        t = clock::now();
        const auto j = *element_vectors(d, g, p, VectorField::current);
        const auto jn = nodal_vectors(g, j);
        const double t_j = seconds(t);
        t = clock::now();
        const IsoSurface iso = isosurface(g, p.fields.potential_V, 0.5 * (*std::ranges::min_element(p.fields.potential_V) +
                                                                          *std::ranges::max_element(p.fields.potential_V)));
        const double t_iso = seconds(t);
        t = clock::now();
        const auto field = streamline_field(d, g, j, VectorField::current);
        const Locator loc(g);
        const double t_prep = seconds(t);
        t = clock::now();
        Plane plane;
        plane.normal = {1.0, 0.0, 0.0};
        plane.offset = 0.5 * (d.mesh().points().front()[0] + d.mesh().points().back()[0]);
        std::size_t points = 0;
        for (const Point& s : plane_seeds(g, plane, 100)) points += trace_streamline(loc, field, s).points.size();
        const double t_lines = seconds(t);
        std::printf("derived benchmark %s: %zu nodes, %zu tetrahedra; E %.3f s, J %.3f s, isosurface %.3f s "
                    "(%zu triangles), streamline field and locator %.3f s, 100 streamlines %.3f s (%zu points)\n",
                    name, d.mesh().node_count(), g.element_count(), t_e, t_j, t_iso, iso.triangles.size(), t_prep,
                    t_lines, points);
        CHECK(en.size() == jn.size());
        // Unit 25c (CG8): 400 surface arrows, a quarter-cut clipped surface, the exploded view's scene.
        t = clock::now();
        double hs = 0.0;
        const auto arrows = surface_arrows(loc, field, exposed_boundary(d, g), 400, hs);
        const double t_arrows = seconds(t);
        const auto pts = d.mesh().points();
        Point lo = pts[0], hi = pts[0];
        for (const Point& q : pts) for (std::size_t a = 0; a < 3; ++a) lo[a] = std::min(lo[a], q[a]), hi[a] = std::max(hi[a], q[a]);
        Box quarter;
        for (std::size_t a = 0; a < 3; ++a) quarter.lo[a] = 0.5 * (lo[a] + hi[a]), quarter.hi[a] = hi[a] + 1.0;
        t = clock::now();
        const auto clipped = g.clipped_surface(quarter);
        const double t_clip = seconds(t);
        t = clock::now();
        const auto exploded = exploded_view(d, g, 1, default_gap(d, 1));
        const double t_exploded = seconds(t);
        REQUIRE(exploded.has_value());
        std::printf("unit 25c benchmark %s: 400 surface arrows %.3f s (%zu drawn), quarter-cut surface %.3f s (%zu polygons), "
                    "exploded view %.3f s\n",
                    name, t_arrows, arrows.size(), t_clip, clipped.size(), t_exploded);
    };
    const Mosfet& m3 = mosfet(true);
    run("3D MOSFET", m3.device, m3.point);

    const auto axis = uniform(0.0, 1e-4, 57);  // 56^3 cells, 6 tetrahedra each
    mesh::Mesh m = *mesh::make_tensor_grid(axis, axis, axis);
    auto cells = *mesh::TensorCells::create(m, axis, axis, axis);
    const std::size_t n = m.node_count();
    auto left = m.find_boundary("x_min")->nodes;
    auto right = m.find_boundary("x_max")->nodes;
    results::BiasPoint p;
    for (const Point& x : m.points()) p.fields.potential_V.push_back(std::sin(3e4 * x[0]) * std::cos(2e4 * x[1]) + 1e4 * x[2]);
    for (const mesh::Edge& ed : m.edges()) {
        const Point& a = m.points()[static_cast<std::size_t>(ed.first)];
        const Point& b = m.points()[static_cast<std::size_t>(ed.second)];
        const Point mid{0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2])};
        const Point jv{1.0 + 1e4 * mid[1], std::sin(4e4 * mid[0]), 0.3};  // A/cm^2
        const double along = (jv[0] * (b[0] - a[0]) + jv[1] * (b[1] - a[1]) + jv[2] * (b[2] - a[2])) / ed.length;
        p.edge_current_n.push_back(along * ed.coupling_area);
        p.edge_current_p.push_back(0.0);
    }
    auto d = device::Device::create({.mesh = std::move(m),
                                     .temperature_K = 300.0,
                                     .regions = {{"silicon", physics::silicon()}},
                                     .node_region = std::vector<device::RegionId>(n, 0),
                                     .donors = std::vector<double>(n, 1e16),
                                     .acceptors = std::vector<double>(n, 0.0),
                                     .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                                                  {"right", device::ContactKind::ohmic, std::move(right)}},
                                     .cells = std::move(cells)});
    REQUIRE(d.has_value());
    run("1e6 tetrahedra", *d, p);
}

// Unit 26n (DT2): a scene's text labels, drawn centred just above their point, in 2D and 3D; none by default.
TEST_CASE("render image: scene labels are drawn above their point and nowhere else (Unit 26n, DT2)") {
    Renderer& r = warp();
    constexpr int W = 800, H = 600;
    for (const bool three : {false, true}) {
        CAPTURE(three);
        const Mosfet& m = mosfet(three);
        const Geometry g = geometry_of(m.device);
        auto s = solid_view(m.device, g);
        REQUIRE(s.has_value());
        REQUIRE(s->labels.empty());
        FieldScene labelled = *s;
        const Float3 at{0.5f * (s->lo[0] + s->hi[0]), 0.5f * (s->lo[1] + s->hi[1]), 0.5f * (s->lo[2] + s->hi[2])};
        labelled.labels.push_back({at, "20 nm", {0.0f, 0.0f, 0.0f, 1.0f}});
        const auto with = r.draw(labelled, W, H);
        const auto without = r.draw(*s, W, H);
        REQUIRE(with.has_value());
        REQUIRE(without.has_value());
        const ScreenPoint q = screen_position(*s, W, H, at);
        const double a = mean_luminance(*with, q.x, q.y - 7.0f, 6);
        const double b = mean_luminance(*without, q.x, q.y - 7.0f, 6);
        CAPTURE(q.x, q.y, a, b);
        CHECK(a < b - 10.0);  // dark text over the device
        CHECK(mean_luminance(*with, 30.0f, 30.0f, 4) == mean_luminance(*without, 30.0f, 30.0f, 4));
        CHECK(mean_luminance(*with, q.x, q.y + 40.0f, 4) == mean_luminance(*without, q.x, q.y + 40.0f, 4));
    }
    CHECK(r.debug_errors() == 0);
}
