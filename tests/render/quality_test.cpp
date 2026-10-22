// Mesh quality (Unit 26b, XD7, XG6): exact angles of equilateral, right and regular elements, the
// elements of tensor meshes, the generators' own measures, the histogram counting every element once,
// the summary, the quality map and the histogram plot.
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <numbers>
#include <numeric>
#include <optional>
#include <utility>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/render/quality.hpp"
#include "NiTCAD/render/scene.hpp"

using namespace NiTCAD;
using namespace NiTCAD::render;
using Catch::Matchers::WithinAbs;

namespace {

constexpr double degrees = 180.0 / std::numbers::pi;

std::vector<double> uniform(double a, double b, int n) {
    std::vector<double> v;
    for (int k = 0; k < n; ++k) v.push_back(k == n - 1 ? b : a + (b - a) * k / (n - 1));
    return v;
}

// A silicon device on a mesh (one ohmic contact at its first node, as a device needs; the quality views draw
// none), with its geometry. The device is held at a fixed address: the geometry refers to its mesh.
struct Meshed {
    std::unique_ptr<device::Device> device;
    std::optional<Geometry> geometry;
};

Meshed silicon(mesh::Mesh m, std::optional<mesh::TensorCells> cells, const mesh::SimplexMesh* simplices) {
    const std::size_t n = m.node_count();
    auto d = device::Device::create({.mesh = std::move(m),
                                     .temperature_K = 300.0,
                                     .regions = {{"silicon", physics::silicon()}},
                                     .node_region = std::vector<device::RegionId>(n, 0),
                                     .donors = std::vector<double>(n, 1e16),
                                     .acceptors = std::vector<double>(n, 0.0),
                                     .contacts = {{"left", device::ContactKind::ohmic, {mesh::NodeId{0}}}},
                                     .cells = std::move(cells)});
    if (!d) FAIL(d.error().message);
    Meshed out{std::make_unique<device::Device>(std::move(*d)), std::nullopt};
    auto g = simplices ? Geometry::simplex(out.device->mesh(), *simplices)
                       : Geometry::tensor(out.device->mesh(), *out.device->cells());
    REQUIRE(g.has_value());
    out.geometry = std::move(*g);
    return out;
}

Meshed tensor2d() {
    const auto x = uniform(0.0, 2e-4, 5), y = uniform(0.0, 1e-4, 3);
    auto m = mesh::make_tensor_grid(x, y);
    REQUIRE(m.has_value());
    auto c = mesh::TensorCells::create(*m, x, y);
    REQUIRE(c.has_value());
    return silicon(std::move(*m), std::move(*c), nullptr);
}

Meshed tensor3d() {
    const auto x = uniform(0.0, 1e-4, 3), y = uniform(0.0, 1e-4, 3), z = uniform(0.0, 1e-4, 3);
    auto m = mesh::make_tensor_grid(x, y, z);
    REQUIRE(m.has_value());
    auto c = mesh::TensorCells::create(*m, x, y, z);
    REQUIRE(c.has_value());
    return silicon(std::move(*m), std::move(*c), nullptr);
}

mesh::GeneratedMesh generated2d(double min_angle) {
    mesh::Geometry2D geo;
    geo.polygons = {{{{0, 0, 0}, {2e-4, 0, 0}, {2e-4, 1e-4, 0}, {0.6e-4, 1.3e-4, 0}, {0, 1e-4, 0}}, 0}};
    geo.patches = {{"left", {0, 0, 0}, {0, 1e-4, 0}}, {"right", {2e-4, 0, 0}, {2e-4, 1e-4, 0}}};
    mesh::MeshOptions2D o;
    o.max_size = 0.25e-4;
    o.boxes = {{{0.5e-4, 0.2e-4, 0}, {0.9e-4, 0.5e-4, 0}, 0.04e-4}};  // grading makes refinement work
    o.min_angle_deg = min_angle;
    auto gm = mesh::generate_mesh(geo, o);
    REQUIRE(gm.has_value());
    return std::move(*gm);
}

mesh::GeneratedMesh generated3d() {
    mesh::Geometry3D geo;
    const double L = 1e-4;
    geo.domain.vertices = {{0, 0, 0}, {L, 0, 0}, {L, L, 0}, {0, L, 0}, {0, 0, L}, {L, 0, L}, {L, L, L}, {0, L, L}};
    geo.domain.faces = {{{0, 3, 2, 1}, "z0", {}}, {{4, 5, 6, 7}, "z1", {}}, {{0, 1, 5, 4}, "y0", {}},
                        {{2, 3, 7, 6}, "y1", {}}, {{0, 4, 7, 3}, "x0", {}}, {{1, 2, 6, 5}, "x1", {}}};
    mesh::MeshOptions3D o;
    o.max_size = 0.3e-4;
    auto gm = mesh::generate_mesh(geo, o);
    REQUIRE(gm.has_value());
    return std::move(*gm);
}

}  // namespace

TEST_CASE("render quality: exact angles of equilateral, right and regular elements (XG6)", "[render]") {
    const double s3 = std::sqrt(3.0);
    const std::vector<Point> equilateral{{0, 0, 0}, {1, 0, 0}, {0.5, s3 / 2, 0}};
    CHECK_THAT(min_angle_deg(equilateral), WithinAbs(60.0, 1e-12));
    const std::vector<Point> right{{0, 0, 0}, {2e-5, 0, 0}, {0, 2e-5, 0}};
    CHECK_THAT(min_angle_deg(right), WithinAbs(45.0, 1e-12));
    // 30-60-90, in any order of the vertices.
    std::vector<Point> half{{0, 0, 0}, {s3, 0, 0}, {0, 1, 0}};
    for (int k = 0; k < 3; ++k) {
        CHECK_THAT(min_angle_deg(half), WithinAbs(30.0, 1e-12));
        std::ranges::rotate(half, half.begin() + 1);
    }
    // The regular tetrahedron: arccos(1/3) at every edge; either orientation.
    std::vector<Point> regular{{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}};
    CHECK_THAT(min_angle_deg(regular), WithinAbs(std::acos(1.0 / 3.0) * degrees, 1e-12));
    CHECK_THAT(min_angle_deg(regular), WithinAbs(70.5288, 1e-4));
    std::swap(regular[0], regular[1]);
    CHECK_THAT(min_angle_deg(regular), WithinAbs(std::acos(1.0 / 3.0) * degrees, 1e-12));
    CHECK_THAT(best_quality_deg(3), WithinAbs(std::acos(1.0 / 3.0) * degrees, 1e-15));
    CHECK(best_quality_deg(2) == 60.0);
    // The corner of a cube: three right angles at the corner's edges, arccos(1/sqrt 3) at the others.
    const std::vector<Point> corner{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    CHECK_THAT(min_angle_deg(corner), WithinAbs(std::acos(1.0 / s3) * degrees, 1e-12));
    // A cube's Kuhn tetrahedron (the path simplex): dihedral angles of 90, 60 and 45 degrees.
    const std::vector<Point> kuhn{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {1, 1, 1}};
    CHECK_THAT(min_angle_deg(kuhn), WithinAbs(45.0, 1e-12));
    // Degenerate elements: 0.
    CHECK(min_angle_deg(std::vector<Point>{{0, 0, 0}, {1, 0, 0}, {2, 0, 0}}) == 0.0);
    // Collinear to rounding: 0 exactly, not an angle of 1e-14 degrees.
    CHECK(min_angle_deg(std::vector<Point>{{0, 0, 0}, {1, 0, 0}, {0.5, 1e-16, 0}}) == 0.0);
    CHECK(min_angle_deg(std::vector<Point>{{0, 0, 0}, {1, 0, 0}, {0.5, 1e-12, 0}}) > 0.0);
    CHECK(min_angle_deg(std::vector<Point>{{0, 0, 0}, {0, 0, 0}, {0, 1, 0}}) == 0.0);
    CHECK(min_angle_deg(std::vector<Point>{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}}) == 0.0);
    CHECK(min_angle_deg(std::vector<Point>{{0, 0, 0}, {1, 0, 0}}) == 0.0);
    // A sliver: four points near a square, lifted 1e-3 out of its plane, has a dihedral angle near 0.
    const double sliver = min_angle_deg(std::vector<Point>{{0, 0, 0}, {1, 0, 0}, {1, 1, 1e-3}, {0, 1, 0}});
    CHECK(sliver > 0.0);
    CHECK(sliver < 0.2);
}

TEST_CASE("render quality: tensor meshes, and the generators' own measures (XG6)", "[render]") {
    // Square cells split into right isosceles triangles: 45 degrees each, 2 per cell.
    const Meshed t2 = tensor2d();
    auto q2 = element_quality(*t2.geometry);
    REQUIRE(q2.has_value());
    CHECK(q2->size() == 2 * 4 * 2);
    for (const double a : *q2) CHECK_THAT(a, WithinAbs(45.0, 1e-12));
    // Cubes split into Kuhn tetrahedra: 45 degrees each, 6 per cell.
    const Meshed t3 = tensor3d();
    auto q3 = element_quality(*t3.geometry);
    REQUIRE(q3.has_value());
    CHECK(q3->size() == 6 * 8);
    for (const double a : *q3) CHECK_THAT(a, WithinAbs(45.0, 1e-12));

    // The 2D generator refines every triangle below its minimum-angle setting unless it is exempt
    // (in an interface band or of interface nodes alone): at most the exempt ones fall below it,
    // and its report's minimum is the minimum over the others.
    for (const double setting : {20.7, 25.0, 30.0}) {
        const mesh::GeneratedMesh gm = generated2d(setting);
        const Meshed g = silicon(gm.mesh, std::nullopt, &gm.simplices);
        auto q = element_quality(*g.geometry);
        REQUIRE(q.has_value());
        CHECK(q->size() == gm.report.triangles);
        const auto below = static_cast<std::size_t>(std::ranges::count_if(*q, [&](double a) { return a < setting; }));
        CHECK(below <= gm.report.exempt_triangles);
        if (gm.report.exempt_triangles == 0) {
            CHECK(std::ranges::min(*q) >= setting);
            CHECK_THAT(std::ranges::min(*q), WithinAbs(gm.report.min_angle_deg, 1e-9));
        }
    }
    // The 3D generator reports the smallest dihedral angle over all its tetrahedra.
    const mesh::GeneratedMesh gm3 = generated3d();
    const Meshed g3 = silicon(gm3.mesh, std::nullopt, &gm3.simplices);
    auto q = element_quality(*g3.geometry);
    REQUIRE(q.has_value());
    CHECK(q->size() == gm3.report.tetrahedra);
    CHECK_THAT(std::ranges::min(*q), WithinAbs(gm3.report.min_dihedral_deg, 1e-9));

    // A 1D mesh or a graph without elements has no quality.
    auto line = mesh::make_tensor_grid(uniform(0.0, 1e-4, 5));
    REQUIRE(line.has_value());
    CHECK(element_quality(Geometry::graph(*line)).error().message == "mesh quality needs a 2D or 3D mesh");
    CHECK(element_quality(Geometry::graph(t2.device->mesh())).error().message == "mesh quality needs a mesh with elements");
}

TEST_CASE("render quality: the histogram counts every element once; the summary", "[render]") {
    const std::vector<double> angles{0.0, 4.999, 5.0, 19.9, 20.0, 45.0, 60.0, 70.5, 89.99, 90.0, 120.0,
                                     std::numeric_limits<double>::quiet_NaN()};
    const QualityHistogram h = quality_histogram(angles);
    CHECK(std::accumulate(h.counts.begin(), h.counts.end(), std::size_t{0}) == angles.size());
    CHECK(h.counts[0] == 3);   // 0, 4.999 and the NaN
    CHECK(h.counts[1] == 1);   // 5
    CHECK(h.counts[3] == 1);   // 19.9
    CHECK(h.counts[4] == 1);   // 20
    CHECK(h.counts[9] == 1);   // 45
    CHECK(h.counts[12] == 1);  // 60
    CHECK(h.counts[14] == 1);  // 70.5
    CHECK(h.counts[17] == 3);  // 89.99, 90 and beyond
    // Every element of a generated mesh once.
    const mesh::GeneratedMesh gm = generated2d(20.7);
    const Meshed g = silicon(gm.mesh, std::nullopt, &gm.simplices);
    const auto q = element_quality(*g.geometry);
    REQUIRE(q.has_value());
    const QualityHistogram hg = quality_histogram(*q);
    CHECK(std::accumulate(hg.counts.begin(), hg.counts.end(), std::size_t{0}) == q->size());

    const QualitySummary odd = quality_summary(std::vector<double>{30.0, 12.0, 50.0, 20.0, 19.0}, 2);
    CHECK(odd.elements == 5);
    CHECK(odd.min_deg == 12.0);
    CHECK(odd.median_deg == 20.0);
    CHECK(odd.threshold_deg == 20.0);
    CHECK(odd.below == 2);  // 12 and 19: strictly below
    const QualitySummary even = quality_summary(std::vector<double>{9.0, 40.0, 10.0, 30.0}, 3);
    CHECK(even.median_deg == 20.0);
    CHECK(even.threshold_deg == 10.0);
    CHECK(even.below == 1);
    const QualitySummary none = quality_summary({}, 3);
    CHECK(none.elements == 0);
    CHECK(none.threshold_deg == 10.0);
}

TEST_CASE("render quality: the map and the histogram plot", "[render]") {
    // 2D: each element in one color, its angle's position on 0..60; the edges and the outline.
    const mesh::GeneratedMesh gm = generated2d(20.7);
    const Meshed g = silicon(gm.mesh, std::nullopt, &gm.simplices);
    const auto q = element_quality(*g.geometry);
    REQUIRE(q.has_value());
    auto s = quality_view(*g.device, *g.geometry, *q);
    REQUIRE(s.has_value());
    REQUIRE(s->colorbar.has_value());
    CHECK(s->colorbar->scale.lo() == 0.0);
    CHECK(s->colorbar->scale.hi() == 60.0);
    CHECK(s->colorbar->title == "minimum angle (\xC2\xB0)");
    REQUIRE(s->triangles.size() == 1);
    const TriangleBatch& b = s->triangles.front();
    CHECK(b.shading == Shading::field);
    REQUIRE(b.vertices.size() == 3 * q->size());
    for (std::size_t e = 0; e < q->size(); ++e) {
        const float at = static_cast<float>((*q)[e] / 60.0);
        for (std::size_t k = 0; k < 3; ++k) CHECK(b.vertices[3 * e + k].value == at);
    }
    CHECK(std::ranges::any_of(s->lines, [](const LineBatch& l) { return l.name == "mesh" && l.on_top; }));
    CHECK(std::ranges::any_of(s->lines, [](const LineBatch& l) { return l.name == "outline"; }));
    CHECK(quality_view(*g.device, *g.geometry, std::vector<double>(q->size() - 1, 30.0)).error().message ==
          "the mesh quality is not one angle per element");
    const std::vector<Plane> plane{{{1.0, 0.0, 0.0}, 0.5e-4}};
    CHECK(quality_view(*g.device, *g.geometry, *q, plane).error().message == "slices need a 3D geometry");

    // 3D: the boundary faces in their tetrahedra's colors on 0..arccos(1/3); with planes, the sections.
    const Meshed t3 = tensor3d();
    const auto q3 = element_quality(*t3.geometry);
    REQUIRE(q3.has_value());
    auto s3 = quality_view(*t3.device, *t3.geometry, *q3);
    REQUIRE(s3.has_value());
    CHECK(s3->colorbar->scale.hi() == best_quality_deg(3));
    CHECK(s3->colorbar->title == "minimum dihedral angle (\xC2\xB0)");
    REQUIRE(s3->triangles.size() == 1);
    CHECK(s3->triangles.front().vertices.size() == 3 * t3.geometry->boundary_faces().size());
    for (const SceneVertex& v : s3->triangles.front().vertices) {
        CHECK_THAT(v.value, WithinAbs(45.0 / best_quality_deg(3), 1e-6));
    }
    auto sliced = quality_view(*t3.device, *t3.geometry, *q3, plane);
    REQUIRE(sliced.has_value());
    CHECK(sliced->triangles.front().name == "slice");
    CHECK(sliced->triangles.back().name == "outline");
    const std::vector<Plane> outside{{{1.0, 0.0, 0.0}, 5e-4}};
    CHECK(quality_view(*t3.device, *t3.geometry, *q3, outside).error().message == "a slice plane misses the device");

    // The histogram: a step through every bin's count from 0 to 90 degrees, closed at 0; the threshold.
    auto p = quality_plot(*q, 2);
    REQUIRE(p.has_value());
    REQUIRE(p->series.size() == 2);
    const Series& steps = p->series[0];
    const QualityHistogram h = quality_histogram(*q);
    REQUIRE(steps.x.size() == 2 + 2 * h.counts.size());
    CHECK(steps.x.front() == 0.0);
    CHECK(steps.y.front() == 0.0);
    CHECK(steps.x.back() == 90.0);
    CHECK(steps.y.back() == 0.0);
    for (std::size_t k = 0; k < h.counts.size(); ++k) {
        CHECK(steps.x[1 + 2 * k] == 5.0 * static_cast<double>(k));
        CHECK(steps.x[2 + 2 * k] == 5.0 * static_cast<double>(k + 1));
        CHECK(steps.y[1 + 2 * k] == static_cast<double>(h.counts[k]));
        CHECK(steps.y[2 + 2 * k] == static_cast<double>(h.counts[k]));
    }
    CHECK(p->series[1].x == std::vector<double>{20.0, 20.0});
    CHECK(p->x.lo == 0.0);
    CHECK(p->x.hi == 90.0);
    CHECK(p->x.label == "minimum angle (\xC2\xB0)");
    CHECK(quality_plot(*q3, 3)->series[1].x == std::vector<double>{10.0, 10.0});
    CHECK(quality_plot({}, 2).error().message == "a mesh quality histogram needs an element");
    CHECK(quality_plot(*q, 1).error().message == "mesh quality needs a 2D or 3D mesh");
}
