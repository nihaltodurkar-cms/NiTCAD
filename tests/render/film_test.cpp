// Filmstrips (Unit 25c, DC4, CG4): images tiled in order; frames on one scale fitted over all of them,
// titled with their labels; time labels in engineering units.
#include <catch2/catch_test_macros.hpp>

#include <cstdint>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/render/image.hpp"
#include "NiTCAD/render/scene.hpp"

using namespace NiTCAD;
using namespace NiTCAD::render;

namespace {

Image solid(int w, int h, std::uint8_t v) {
    Image i;
    i.width = w;
    i.height = h;
    i.rgba.assign(static_cast<std::size_t>(w) * static_cast<std::size_t>(h) * 4, v);
    return i;
}

std::uint8_t at(const Image& i, int x, int y) {
    return i.rgba[(static_cast<std::size_t>(y) * static_cast<std::size_t>(i.width) + static_cast<std::size_t>(x)) * 4];
}

}  // namespace

TEST_CASE("render film: images tile in rows, in order, on white (CG4)") {
    const std::vector<Image> frames{solid(4, 3, 10), solid(4, 3, 20), solid(4, 3, 30), solid(4, 3, 40), solid(4, 3, 50)};
    const auto t = tile_images(frames, 2);
    REQUIRE(t.has_value());
    CHECK(t->width == 8);
    CHECK(t->height == 9);
    CHECK(at(*t, 0, 0) == 10);
    CHECK(at(*t, 7, 2) == 20);
    CHECK(at(*t, 0, 3) == 30);
    CHECK(at(*t, 5, 5) == 40);
    CHECK(at(*t, 3, 8) == 50);
    CHECK(at(*t, 4, 6) == 255);  // the empty sixth place
    // More columns than images: one row.
    CHECK(tile_images(frames, 9)->width == 20);
    CHECK_FALSE(tile_images(std::vector<Image>{}, 2).has_value());
    CHECK_FALSE(tile_images(frames, 0).has_value());
    CHECK_FALSE(tile_images(std::vector<Image>{solid(4, 3, 1), solid(3, 3, 1)}, 2).has_value());
}

TEST_CASE("render film: frames share one scale and carry their labels (CG4, CQ2)") {
    const std::vector<double> x{0.0, 1e-4, 2e-4}, y{0.0, 1e-4};
    mesh::Mesh m = *mesh::make_tensor_grid(x, y);
    auto cells = *mesh::TensorCells::create(m, x, y);
    auto left = m.find_boundary("x_min")->nodes;
    auto right = m.find_boundary("x_max")->nodes;
    const std::size_t n = m.node_count();
    const auto d = *device::Device::create({.mesh = std::move(m),
                                            .temperature_K = 300.0,
                                            .regions = {{"silicon", physics::silicon()}},
                                            .node_region = std::vector<device::RegionId>(n, 0),
                                            .donors = std::vector<double>(n, 1e16),
                                            .acceptors = std::vector<double>(n, 0.0),
                                            .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                                                         {"right", device::ContactKind::ohmic, std::move(right)}},
                                            .cells = std::move(cells)});
    const Geometry g = *Geometry::tensor(d.mesh(), *d.cells());
    // Frame 1 spans 1..2, frame 2 spans 3..5: one scale 1..5 for both.
    const std::vector<FilmFrame> frames{{time_label(0.0), std::vector<double>(n, 1.0)},
                                        {time_label(2.5e-10), std::vector<double>(n, 5.0)}};
    std::vector<FilmFrame> spread = frames;
    spread[0].values[0] = 2.0;
    spread[1].values[0] = 3.0;
    FieldOptions o;
    o.title = "psi (V)";
    const auto s = filmstrip_scenes(d, g, spread, ScaleKind::linear, o);
    REQUIRE(s.has_value());
    REQUIRE(s->size() == 2);
    for (const FieldScene& f : *s) {
        REQUIRE(f.colorbar.has_value());
        CHECK(f.colorbar->scale.lo() == 1.0);
        CHECK(f.colorbar->scale.hi() == 5.0);
    }
    CHECK((*s)[0].title == "psi (V), t = 0 s");
    CHECK((*s)[1].title == "psi (V), t = 250 ps");
    // The frame of 5s sits at the top of the shared scale, the frame of 1s at its bottom.
    for (const TriangleBatch& b : (*s)[1].triangles) {
        if (b.name != "field") continue;
        bool top = false;
        for (const SceneVertex& v : b.vertices) top = top || v.value == 1.0f;
        CHECK(top);
    }
    // A given scale is used as is.
    o.scale = *ColorScale::create(ScaleKind::linear, 0.0, 10.0);
    CHECK(filmstrip_scenes(d, g, spread, ScaleKind::linear, o)->front().colorbar->scale.hi() == 10.0);
    CHECK_FALSE(filmstrip_scenes(d, g, std::vector<FilmFrame>{}, ScaleKind::linear, o).has_value());
}

TEST_CASE("render film: time labels") {
    CHECK(time_label(0.0) == "t = 0 s");
    CHECK(time_label(1.5) == "t = 1.50 s");
    CHECK(time_label(2.5e-3) == "t = 2.50 ms");
    CHECK(time_label(1e-6) == "t = 1.00 \xC2\xB5s");
    CHECK(time_label(1.234e-9) == "t = 1.23 ns");
    CHECK(time_label(12.5e-9) == "t = 12.5 ns");
    CHECK(time_label(325e-12) == "t = 325 ps");
    CHECK(time_label(9.9996e-10) == "t = 1.00 ns");
    CHECK(time_label(3e-15) == "t = 3.00 fs");
}
