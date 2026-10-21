// The exploded view (Unit 25c, DC3, CG3): each region's surface is closed and outward and encloses its
// nodes' control volumes; at zero gap the surfaces cover the boundary and both sides of every
// interface; regions and their contacts move by their rank's offset along the axis.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/render/scene.hpp"

using namespace NiTCAD;
using namespace NiTCAD::render;

namespace {

std::vector<double> graded(double a, double b, int n) {
    std::vector<double> v;
    for (int k = 0; k < n; ++k) {
        const double s = static_cast<double>(k) / (n - 1);
        v.push_back(k == n - 1 ? b : a + (b - a) * (s + 0.15 * s * (1.0 - s)));
    }
    return v;
}

Point sub(const Point& a, const Point& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Point cross(const Point& a, const Point& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double dot(const Point& a, const Point& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

// Silicon below y = 0.8 um, oxide above, a metal electrode on the oxide's top over the middle, ohmic
// contacts on the silicon's x faces.
device::Device mos3() {
    const auto x = graded(0.0, 2e-4, 9), y = graded(0.0, 1e-4, 8), z = graded(0.0, 0.6e-4, 4);
    mesh::Mesh m = *mesh::make_tensor_grid(x, y, z);
    auto cells = *mesh::TensorCells::create(m, x, y, z);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n, 1);
    std::vector<double> donors(n, 1e16);
    for (std::size_t i = 0; i < n; ++i) {
        if (m.points()[i][1] > 0.8e-4) region[i] = 0, donors[i] = 0.0;
    }
    const auto on = [&](const char* patch, int want, double x_lo, double x_hi) {
        std::vector<mesh::NodeId> v;
        for (const mesh::NodeId k : m.find_boundary(patch)->nodes) {
            const double px = m.points()[static_cast<std::size_t>(k)][0];
            if (region[static_cast<std::size_t>(k)] == want && px >= x_lo && px <= x_hi) v.push_back(k);
        }
        return v;
    };
    device::Contact gate{"gate", device::ContactKind::electrode, on("y_max", 0, 0.6e-4, 1.4e-4)};
    auto left = on("x_min", 1, -1.0, 1.0);
    auto right = on("x_max", 1, -1.0, 1.0);  // before the mesh moves into the device
    gate.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.4};
    auto d = device::Device::create({.mesh = std::move(m),
                                     .temperature_K = 300.0,
                                     .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
                                     .node_region = std::move(region),
                                     .donors = std::move(donors),
                                     .acceptors = std::vector<double>(n, 0.0),
                                     .contacts = {std::move(gate),
                                                  {"left", device::ContactKind::ohmic, std::move(left)},
                                                  {"right", device::ContactKind::ohmic, std::move(right)}},
                                     .cells = std::move(cells)});
    if (!d) FAIL(d.error().message);
    return std::move(*d);
}

}  // namespace

TEST_CASE("render exploded: each region's surface is closed and holds its control volumes (CG3)") {
    const device::Device d = mos3();
    const Geometry g = *Geometry::tensor(d.mesh(), *d.cells());
    const std::vector<std::int32_t> regions(d.node_region().begin(), d.node_region().end());
    const auto surfaces = g.region_surfaces(regions);
    REQUIRE(surfaces.size() == 2);
    double total_area = 0.0;
    for (std::size_t r = 0; r < 2; ++r) {
        CAPTURE(r);
        double volume = 0.0, area = 0.0;
        Point vector_area{};
        for (const auto& f : surfaces[r]) {
            CHECK(regions[static_cast<std::size_t>(f.owner)] == static_cast<std::int32_t>(r));
            for (std::size_t i = 1; i + 1 < f.points.size(); ++i) {
                volume += dot(f.points[0], cross(f.points[i], f.points[i + 1])) / 6.0;
                const Point n = cross(sub(f.points[i], f.points[0]), sub(f.points[i + 1], f.points[0]));
                for (std::size_t k = 0; k < 3; ++k) vector_area[k] += 0.5 * n[k];
                area += 0.5 * std::sqrt(dot(n, n));
            }
        }
        double expected = 0.0;
        for (std::size_t i = 0; i < regions.size(); ++i) {
            if (regions[i] == static_cast<std::int32_t>(r)) expected += d.mesh().volumes()[i];
        }
        CHECK(std::abs(volume - expected) <= 1e-10 * expected);  // closed and outward: the volume, positive
        for (std::size_t k = 0; k < 3; ++k) CHECK(std::abs(vector_area[k]) <= 1e-10 * area);
        total_area += area;
    }
    // The boundary and both sides of the interface (the plane y = the oxide's lowest control-volume face).
    double boundary = 0.0;
    for (const Polygon& f : g.boundary_faces()) {
        std::vector<Point> p;
        for (const SurfaceVertex& v : f.vertices) p.push_back(v.position);
        boundary += polygon_area(p);
    }
    const double interface_area = 2e-4 * 0.6e-4;
    CHECK(std::abs(total_area - (boundary + 2.0 * interface_area)) <= 1e-10 * total_area);
}

TEST_CASE("render exploded: regions and their contacts move by rank along the axis (CG3)") {
    const device::Device d = mos3();
    const Geometry g = *Geometry::tensor(d.mesh(), *d.cells());
    const double gap = 0.3e-4;  // cm
    const auto s = exploded_view(d, g, 1, gap);
    REQUIRE(s.has_value());
    // The oxide (the upper region along y) moves up by one gap; silicon stays. Oxide batches are the
    // translucent insulator class; its gate electrode moves with it.
    const auto zero = exploded_view(d, g, 1, 0.0);
    REQUIRE(zero.has_value());
    const auto max_y = [](const FieldScene& sc, const char* name) {
        float y = -1e9f;
        for (const TriangleBatch& b : sc.triangles) {
            if (b.name != name) continue;
            for (const SceneVertex& v : b.vertices) y = std::max(y, v.position[1]);
        }
        return y;
    };
    const std::string insulator(material_name(Material::insulator));
    CHECK(std::abs(max_y(*s, insulator.c_str()) - max_y(*zero, insulator.c_str()) - 0.3f) <= 1e-5f);
    CHECK(std::abs(max_y(*s, "electrode") - max_y(*zero, "electrode") - 0.3f) <= 1e-5f);
    // Silicon (n, light doping) stays put.
    const std::string n_light(material_name(Material::n_light));
    CHECK(max_y(*s, n_light.c_str()) == max_y(*zero, n_light.c_str()));
    // With a field the silicon carries it (positions of the field's scale) and the oxide keeps its
    // translucent class color.
    std::vector<double> f;
    for (const Point& p : d.mesh().points()) f.push_back(p[0] * 1e4);  // x in um, 0..2
    FieldOptions fo;
    fo.scale = *ColorScale::create(ScaleKind::linear, 0.0, 2.0);
    const auto fs = exploded_view(d, g, 1, gap, f, fo);
    REQUIRE(fs.has_value());
    const auto field = std::ranges::find_if(fs->triangles, [](const TriangleBatch& b) { return b.name == "field"; });
    REQUIRE(field != fs->triangles.end());
    double oxide_lowest = 1.0;  // um
    for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
        if (d.is_insulator(static_cast<NodeId>(i))) oxide_lowest = std::min(oxide_lowest, d.mesh().points()[i][1] * 1e4);
    }
    for (const SceneVertex& v : field->vertices) {
        CHECK(std::abs(v.value - v.position[0] / 2.0f) <= 1e-5f);  // the field's position at the vertex
        CHECK(v.position[1] < oxide_lowest);                        // silicon only (it stays below the oxide)
    }
    CHECK(fs->colorbar.has_value());
    CHECK(std::ranges::none_of(fs->triangles, [&](const TriangleBatch& b) { return b.name == n_light; }));
    // Along x the two regions share their mean, so the order is the regions' and one moves.
    CHECK(exploded_view(d, g, 0, gap).has_value());
    // Refused: a bad axis or gap, a field not one per node, a 2D geometry.
    CHECK(exploded_view(d, g, 3, gap).error().code == base::ErrorCode::invalid_input);
    CHECK(exploded_view(d, g, 1, -1.0).error().code == base::ErrorCode::invalid_input);
    const std::vector<double> short_field(3, 1.0);
    CHECK(exploded_view(d, g, 1, gap, short_field).error().code == base::ErrorCode::invalid_input);
    CHECK(std::abs(default_gap(d, 1) - 0.12e-4) <= 1e-18);
}
