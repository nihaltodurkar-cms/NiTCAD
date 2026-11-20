// Layered 2D geometry resolved into a tiling (Unit 26k, GD2; GE2, GE5): the MOSFET's layers mesh; a hole and a
// part touching itself are cut into simple pieces; 1000 seeded random layouts tile exactly (simple polygons,
// conforming edges, every cell's region the topmost shape's, the areas exact); the errors; the speed.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <format>
#include <random>
#include <vector>

#include "NiTCAD/mesh/generate.hpp"
#include "polygon2d.hpp"

using namespace NiTCAD::mesh;

namespace {

LayeredShape box(double x0, double x1, double y0, double y1, std::int32_t region) {
    return {{{x0, y0, 0}, {x1, y0, 0}, {x1, y1, 0}, {x0, y1, 0}}, region};
}

double area(const std::vector<Point>& v) {
    double a = 0.0;
    for (std::size_t k = 0; k < v.size(); ++k) {
        const Point& p = v[k];
        const Point& q = v[(k + 1) % v.size()];
        a += p[0] * q[1] - q[0] * p[1];
    }
    return 0.5 * a;
}

// A point strictly inside an edge's interior (axis-parallel edges only).
bool inside_edge(const Point& a, const Point& b, const Point& p) {
    if (a[0] == b[0]) return p[0] == a[0] && p[1] > std::min(a[1], b[1]) && p[1] < std::max(a[1], b[1]);
    return p[1] == a[1] && p[0] > std::min(a[0], b[0]) && p[0] < std::max(a[0], b[0]);
}

// The tiling's checks on an integer lattice of unit cells over [0, n)^2: every polygon simple and
// counter-clockwise; no vertex of one polygon inside another's edge; each unit cell's centre strictly inside exactly
// one polygon, of the region `expect` gives; the areas per region the cell counts.
void check_tiling(const Geometry2D& t, int n, const std::vector<std::int32_t>& expect) {
    std::vector<double> region_area(8, 0.0), counted(8, 0.0);
    for (const RegionPolygon& p : t.polygons) {
        REQUIRE(detail::simple_2d(p.vertices));
        CHECK(area(p.vertices) > 0.0);
        region_area[static_cast<std::size_t>(p.region)] += area(p.vertices);
        for (const RegionPolygon& q : t.polygons) {
            if (&p == &q) continue;
            for (std::size_t k = 0; k < q.vertices.size(); ++k) {
                for (const Point& v : p.vertices) CHECK_FALSE(inside_edge(q.vertices[k], q.vertices[(k + 1) % q.vertices.size()], v));
            }
        }
    }
    for (int j = 0; j < n; ++j) {
        for (int i = 0; i < n; ++i) {
            const Point c{i + 0.5, j + 0.5, 0};
            int hits = 0;
            std::int32_t region = -1;
            for (const RegionPolygon& p : t.polygons) {
                if (detail::strictly_inside_2d(p.vertices, c)) {
                    ++hits;
                    region = p.region;
                }
            }
            CHECK(hits == 1);
            CHECK(region == expect[static_cast<std::size_t>(j * n + i)]);
            counted[static_cast<std::size_t>(expect[static_cast<std::size_t>(j * n + i)])] += 1.0;
        }
    }
    CHECK(region_area == counted);
}

}  // namespace

TEST_CASE("mesh layered: a junction under an oxide resolves into a tiling that meshes (Unit 26k, GE2)", "[mesh][layered]") {
    LayeredGeometry2D g{.lo = {0, 0, 0}, .hi = {2e-4, 1e-4, 0}, .background = 0, .shapes = {}, .patches = {}};
    g.shapes = {box(1e-4, 2e-4, 0, 0.8e-4, 1),   // the n side over the p background
                box(0, 2e-4, 0.8e-4, 1e-4, 2)};  // an oxide along the top, over both
    g.patches = {{"anode", {0, 0, 0}, {0, 0.6e-4, 0}}, {"cathode", {2e-4, 0, 0}, {2e-4, 0.6e-4, 0}},
                 {"gate", {0.6e-4, 1e-4, 0}, {1.4e-4, 1e-4, 0}}};
    const auto t = layered_partition2d(g);
    REQUIRE(t.has_value());
    REQUIRE(t->polygons.size() == 3);
    CHECK(t->polygons[0].vertices == std::vector<Point>{{0, 0, 0}, {1e-4, 0, 0}, {1e-4, 0.8e-4, 0}, {0, 0.8e-4, 0}});
    CHECK(t->polygons[1].region == 1);
    // The oxide carries the junction's end on its lower edge (a T-junction).
    CHECK(t->polygons[2].vertices == std::vector<Point>{{0, 0.8e-4, 0}, {1e-4, 0.8e-4, 0}, {2e-4, 0.8e-4, 0},
                                                        {2e-4, 1e-4, 0}, {0, 1e-4, 0}});
    CHECK(t->patches.size() == 3);
    MeshOptions2D o;
    o.max_size = 0.25e-4;
    o.interface_resolution = 20e-7;
    const auto m = generate_mesh(*t, o);
    if (!m) FAIL(m.error().message);
    CHECK(m->simplices.points.size() > 100);
}

TEST_CASE("mesh layered: a hole and a part touching itself are cut into simple pieces", "[mesh][layered]") {
    // A ring: region 1 around a hole of the background.
    LayeredGeometry2D ring{.lo = {0, 0, 0}, .hi = {6, 6, 0}, .background = 0, .shapes = {box(1, 5, 1, 5, 1), box(2, 4, 2, 4, 0)}, .patches = {}};
    auto t = layered_partition2d(ring);
    REQUIRE(t.has_value());
    std::vector<std::int32_t> expect(36, 0);
    for (int j = 1; j < 5; ++j) {
        for (int i = 1; i < 5; ++i) expect[static_cast<std::size_t>(j * 6 + i)] = (i >= 2 && i < 4 && j >= 2 && j < 4) ? 0 : 1;
    }
    check_tiling(*t, 6, expect);
    CHECK(std::ranges::count(t->polygons, 1, &RegionPolygon::region) == 2);  // the ring cut once
    // A C-shape whose ends touch diagonally: one part, its boundary meets itself at (2, 3).
    LayeredGeometry2D c{.lo = {0, 0, 0}, .hi = {4, 4, 0}, .background = 0,
                        .shapes = {box(0, 3, 0, 1, 1), box(0, 1, 0, 3, 1), box(1, 2, 2, 3, 1), box(2, 3, 1, 2, 1)}, .patches = {}};
    // ... ends at (1..2, 2..3) and (2..3, 1..2): touching at (2, 2), connected round through the C.
    t = layered_partition2d(c);
    REQUIRE(t.has_value());
    expect.assign(16, 0);
    for (const auto& [i, j] : {std::pair{0, 0}, {1, 0}, {2, 0}, {0, 1}, {0, 2}, {1, 2}, {2, 1}}) expect[static_cast<std::size_t>(j * 4 + i)] = 1;
    check_tiling(*t, 4, expect);
}

TEST_CASE("mesh layered: 1000 seeded random layouts tile exactly (GE2)", "[mesh][layered]") {
    std::mt19937 rng(26011);
    for (int trial = 0; trial < 1000; ++trial) {
        CAPTURE(trial);
        const int n = 12;
        LayeredGeometry2D g{.lo = {0, 0, 0}, .hi = {n, n, 0}, .background = 0, .shapes = {}, .patches = {}};
        std::vector<std::int32_t> expect(static_cast<std::size_t>(n * n), 0);
        const int count = std::uniform_int_distribution<int>(1, 30)(rng);
        for (int s = 0; s < count; ++s) {
            std::uniform_int_distribution<int> at(-2, n + 2);
            int x0 = at(rng), x1 = at(rng), y0 = at(rng), y1 = at(rng);
            if (x0 == x1 || y0 == y1) continue;
            if (x0 > x1) std::swap(x0, x1);
            if (y0 > y1) std::swap(y0, y1);
            const std::int32_t region = std::uniform_int_distribution<std::int32_t>(0, 4)(rng);
            g.shapes.push_back(box(x0, x1, y0, y1, region));
            for (int j = std::max(y0, 0); j < std::min(y1, n); ++j) {
                for (int i = std::max(x0, 0); i < std::min(x1, n); ++i) expect[static_cast<std::size_t>(j * n + i)] = region;
            }
        }
        const auto t = layered_partition2d(g);
        REQUIRE(t.has_value());
        check_tiling(*t, n, expect);
    }
}

TEST_CASE("mesh layered: the errors", "[mesh][layered]") {
    LayeredGeometry2D g{.lo = {0, 0, 0}, .hi = {4, 4, 0}, .background = -1, .shapes = {box(0, 2, 0, 4, 0)}, .patches = {}};
    auto t = layered_partition2d(g);
    REQUIRE_FALSE(t.has_value());
    CHECK(t.error().message.starts_with("the layers leave part of the domain uncovered"));
    g.background = 0;
    g.shapes = {{{{0, 0, 0}, {2, 0, 0}, {1, 2, 0}, {0, 2, 0}}, 1}};
    t = layered_partition2d(g);
    REQUIRE_FALSE(t.has_value());
    CHECK(t.error().message.find("horizontal or vertical") != std::string::npos);
    g.shapes = {box(0, 2, 0, 2, -1)};
    CHECK_FALSE(layered_partition2d(g).has_value());
    g.hi = {0, 4, 0};
    g.shapes.clear();
    CHECK_FALSE(layered_partition2d(g).has_value());
}

TEST_CASE("mesh layered: 200 rectangles resolve in at most 5 ms (GE5)", "[mesh][layered]") {
    std::mt19937 rng(7);
    LayeredGeometry2D g{.lo = {0, 0, 0}, .hi = {1000, 1000, 0}, .background = 0, .shapes = {}, .patches = {}};
    std::uniform_int_distribution<int> at(0, 1000);
    while (g.shapes.size() < 200) {
        int x0 = at(rng), x1 = at(rng), y0 = at(rng), y1 = at(rng);
        if (x0 == x1 || y0 == y1) continue;
        g.shapes.push_back(box(std::min(x0, x1), std::max(x0, x1), std::min(y0, y1), std::max(y0, y1),
                               static_cast<std::int32_t>(g.shapes.size() % 5)));
    }
    std::vector<double> ms;
    std::size_t polygons = 0;
    for (int k = 0; k < 5; ++k) {
        const auto t0 = std::chrono::steady_clock::now();
        const auto t = layered_partition2d(g);
        ms.push_back(std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
        REQUIRE(t.has_value());
        polygons = t->polygons.size();
    }
    std::ranges::sort(ms);
    WARN(std::format("GE5: 200 rectangles resolved in {:.2f} ms (median of 5), {} polygons", ms[2], polygons));
    CHECK(ms[2] <= 5.0);
}
