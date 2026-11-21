// Layered 2D geometry resolved into a tiling (Unit 26k, GD2; GE2, GE5): the MOSFET's layers mesh; a hole and a
// part touching itself are cut into simple pieces; 1000 seeded random layouts tile exactly (simple polygons,
// conforming edges, every cell's region the topmost shape's, the areas exact); the errors; the speed.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <format>
#include <random>
#include <string>
#include <vector>

#include "NiTCAD/mesh/generate.hpp"
#include "polygon2d.hpp"
#include "predicates.hpp"

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

// The painting by runs repaints a row only where a shape's crossings change (GE5): L shapes, whose crossings change
// part way up, mixed with rectangles, and no background, so that a gap left anywhere must be found.
TEST_CASE("mesh layered: 1000 seeded random L-shaped layouts tile exactly (GE2, GE5)", "[mesh][layered]") {
    std::mt19937 rng(26050);
    int refused = 0;
    for (int trial = 0; trial < 1000; ++trial) {
        CAPTURE(trial);
        const int n = 12;
        LayeredGeometry2D g{.lo = {0, 0, 0}, .hi = {n, n, 0}, .background = -1, .shapes = {}, .patches = {}};
        std::vector<std::int32_t> expect(static_cast<std::size_t>(n * n), -1);
        const int count = std::uniform_int_distribution<int>(1, 30)(rng);
        for (int s = 0; s < count; ++s) {
            std::uniform_int_distribution<int> at(-2, n + 2);
            int x0 = at(rng), x1 = at(rng), y0 = at(rng), y1 = at(rng);
            if (x0 > x1) std::swap(x0, x1);
            if (y0 > y1) std::swap(y0, y1);
            if (x1 - x0 < 2 || y1 - y0 < 2) continue;
            const std::int32_t region = std::uniform_int_distribution<std::int32_t>(0, 4)(rng);
            const int xm = std::uniform_int_distribution<int>(x0 + 1, x1 - 1)(rng);
            const int ym = std::uniform_int_distribution<int>(y0 + 1, y1 - 1)(rng);
            const bool l = std::uniform_int_distribution<int>(0, 1)(rng) == 1;
            const auto d = [](int c) { return static_cast<double>(c); };
            if (l) {  // the rectangle without its top right corner [xm, x1) x [ym, y1)
                g.shapes.push_back({{{d(x0), d(y0), 0}, {d(x1), d(y0), 0}, {d(x1), d(ym), 0}, {d(xm), d(ym), 0},
                                     {d(xm), d(y1), 0}, {d(x0), d(y1), 0}},
                                    region});
            } else {
                g.shapes.push_back(box(x0, x1, y0, y1, region));
            }
            for (int j = std::max(y0, 0); j < std::min(y1, n); ++j) {
                for (int i = std::max(x0, 0); i < std::min(x1, n); ++i) {
                    if (!(l && i >= xm && j >= ym)) expect[static_cast<std::size_t>(j * n + i)] = region;
                }
            }
        }
        const auto t = layered_partition2d(g);
        const auto gap = std::ranges::find(expect, -1);
        if (gap != expect.end()) {
            const auto c = static_cast<int>(gap - expect.begin());
            REQUIRE_FALSE(t.has_value());
            CHECK(t.error().message == "the layers leave part of the domain uncovered (near (" + std::to_string(static_cast<double>(c % n)) +
                                           ", " + std::to_string(static_cast<double>(c / n)) + ")) and no region is the background");
            ++refused;
            continue;
        }
        REQUIRE(t.has_value());
        check_tiling(*t, n, expect);
    }
    // Both outcomes are exercised.
    CHECK(refused > 50);
    CHECK(refused < 950);
}

TEST_CASE("mesh layered: the errors", "[mesh][layered]") {
    LayeredGeometry2D g{.lo = {0, 0, 0}, .hi = {4, 4, 0}, .background = -1, .shapes = {box(0, 2, 0, 4, 0)}, .patches = {}};
    auto t = layered_partition2d(g);
    REQUIRE_FALSE(t.has_value());
    CHECK(t.error().message.starts_with("the layers leave part of the domain uncovered"));
    g.background = 0;
    // A slanted edge is taken since slice 4 (GD10); a repeated vertex is not.
    g.shapes = {{{{0, 0, 0}, {2, 0, 0}, {1, 2, 0}, {0, 2, 0}}, 1}};
    CHECK(layered_partition2d(g).has_value());
    g.shapes = {{{{0, 0, 0}, {2, 0, 0}, {2, 0, 0}, {0, 2, 0}}, 1}};
    t = layered_partition2d(g);
    REQUIRE_FALSE(t.has_value());
    CHECK(t.error().message == "a layered shape repeats a vertex");
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
    // Release only, as KG6 (owner, 2026-11-19): the unoptimised Debug build takes about 20x longer; there it is recorded.
#ifdef NDEBUG
    CHECK(ms[2] <= 5.0);
#endif
}
// --- slice 4: slanted edges (GD10) ----------------------------------------------------------------------------------

namespace {

// A point strictly inside a segment (any direction), exactly.
bool inside_segment(const Point& a, const Point& b, const Point& p) {
    return detail::on_segment_2d(a, b, p) && !(p[0] == a[0] && p[1] == a[1]) && !(p[0] == b[0] && p[1] == b[1]);
}

// The tiling's checks for any layout: polygons simple and counter-clockwise; no vertex of one inside another's edge;
// the areas adding up to the domain's; each sample point (off every input edge) strictly inside exactly one polygon,
// of the topmost shape's region containing it (exact), else the background.
void check_general(const LayeredGeometry2D& g, const Geometry2D& t, const std::vector<Point>& samples) {
    double total = 0.0;
    for (const RegionPolygon& p : t.polygons) {
        REQUIRE(detail::simple_2d(p.vertices));
        CHECK(detail::polygon_turn(p.vertices) > 0);
        total += area(p.vertices);
        for (const RegionPolygon& q : t.polygons) {
            if (&p == &q) continue;
            for (std::size_t k = 0; k < q.vertices.size(); ++k) {
                for (const Point& v : p.vertices) CHECK_FALSE(inside_segment(q.vertices[k], q.vertices[(k + 1) % q.vertices.size()], v));
            }
        }
    }
    const double domain = (g.hi[0] - g.lo[0]) * (g.hi[1] - g.lo[1]);
    CHECK(std::abs(total - domain) <= 1e-9 * domain);
    // A vertex on a horizontal or vertical input edge (the domain's sides too) keeps that edge's constant coordinate
    // exactly: a computed crossing there is not rounded off the line.
    std::vector<std::array<Point, 2>> lines{{g.lo, Point{g.hi[0], g.lo[1], 0}}, {Point{g.hi[0], g.lo[1], 0}, g.hi},
                                            {g.hi, Point{g.lo[0], g.hi[1], 0}}, {Point{g.lo[0], g.hi[1], 0}, g.lo}};
    for (const LayeredShape& s : g.shapes) {
        for (std::size_t k = 0; k < s.vertices.size(); ++k) lines.push_back({s.vertices[k], s.vertices[(k + 1) % s.vertices.size()]});
    }
    const double near = 1e-9 * std::max(g.hi[0] - g.lo[0], g.hi[1] - g.lo[1]);
    for (const RegionPolygon& p : t.polygons) {
        for (const Point& v : p.vertices) {
            for (const auto& [a, b] : lines) {
                for (std::size_t k = 0; k < 2; ++k) {
                    const std::size_t m = 1 - k;
                    if (a[k] != b[k] || std::abs(v[k] - a[k]) > near) continue;
                    if (v[m] <= std::min(a[m], b[m]) || v[m] >= std::max(a[m], b[m])) continue;
                    CAPTURE(v[0], v[1], a[k]);
                    CHECK(v[k] == a[k]);
                }
            }
        }
    }
    for (const Point& c : samples) {
        std::int32_t expect = g.background;
        for (std::size_t s = g.shapes.size(); s-- > 0;) {
            if (detail::strictly_inside_2d(g.shapes[s].vertices, c)) {
                expect = g.shapes[s].region;
                break;
            }
        }
        int hits = 0;
        std::int32_t region = -1;
        for (const RegionPolygon& p : t.polygons) {
            if (detail::strictly_inside_2d(p.vertices, c)) {
                ++hits;
                region = p.region;
            }
        }
        CAPTURE(c[0], c[1]);
        CHECK(hits == 1);
        CHECK(region == expect);
    }
}

// Sample points off every input edge and off the domain's boundary.
std::vector<Point> samples_for(const LayeredGeometry2D& g, std::mt19937& rng, int count) {
    std::vector<Point> out;
    std::uniform_real_distribution<double> ux(g.lo[0], g.hi[0]), uy(g.lo[1], g.hi[1]);
    while (static_cast<int>(out.size()) < count) {
        const Point p{ux(rng), uy(rng), 0};
        bool on = false;
        for (const LayeredShape& s : g.shapes) {
            for (std::size_t k = 0; k < s.vertices.size() && !on; ++k) on = detail::on_segment_2d(s.vertices[k], s.vertices[(k + 1) % s.vertices.size()], p);
        }
        if (!on) out.push_back(p);
    }
    return out;
}

}  // namespace

TEST_CASE("mesh layered: a trapezoid gate over a slanted-wall trench resolves and meshes (GD10)", "[mesh][layered]") {
    // A silicon body; a shallow trench (STI) with slanted walls at the left, filled with oxide; a gate oxide; a
    // poly gate with slanted sidewalls on it. Coordinates in cm (µm × 1e-4).
    const double u = 1e-4;
    LayeredGeometry2D g{.lo = {0, 0, 0}, .hi = {2 * u, 1.2 * u, 0}, .background = 0, .shapes = {}, .patches = {}};
    g.shapes.push_back(box(0, 2 * u, 1.0 * u, 1.2 * u, 3));                                                            // a cap over all
    g.shapes.push_back(box(0, 2 * u, 0.98 * u, 1.0 * u, 1));                                                           // gate oxide
    g.shapes.push_back({{{0, 0.7 * u, 0}, {0.3 * u, 0.7 * u, 0}, {0.4 * u, 1.0 * u, 0}, {0, 1.0 * u, 0}}, 1});          // STI
    g.shapes.push_back({{{0.8 * u, 1.0 * u, 0}, {1.4 * u, 1.0 * u, 0}, {1.3 * u, 1.2 * u, 0}, {0.9 * u, 1.2 * u, 0}}, 2});  // gate
    g.patches = {{"body", {0, 0, 0}, {2 * u, 0, 0}}};
    const auto t = layered_partition2d(g);
    if (!t) FAIL(t.error().message);
    std::mt19937 rng(4);
    check_general(g, *t, samples_for(g, rng, 2000));
    // The STI's slanted wall crosses the gate oxide's foot at a computed point (its y kept exactly), shared.
    int sharing = 0;
    for (const RegionPolygon& p : t->polygons) {
        for (const Point& v : p.vertices) sharing += v[1] == 0.98 * u && std::abs(v[0] - (0.3 + 0.1 * 0.28 / 0.3) * u) < 1e-12 * u;
    }
    CHECK(sharing >= 2);
    MeshOptions2D o;
    o.max_size = 0.1 * u;
    o.interface_resolution = 5e-7;
    const auto m = generate_mesh(*t, o);
    if (!m) FAIL(m.error().message);
    CHECK(m->simplices.points.size() > 100);
}

TEST_CASE("mesh layered: a slanted ring is cut into simple pieces (GD10)", "[mesh][layered]") {
    // A diamond of region 1 with a diamond hole of the background, and a triangle crossing both.
    LayeredGeometry2D g{.lo = {0, 0, 0}, .hi = {8, 8, 0}, .background = 0, .shapes = {}, .patches = {}};
    g.shapes.push_back({{{4, 0.5, 0}, {7.5, 4, 0}, {4, 7.5, 0}, {0.5, 4, 0}}, 1});
    g.shapes.push_back({{{4, 2, 0}, {6, 4, 0}, {4, 6, 0}, {2, 4, 0}}, 0});
    const auto t = layered_partition2d(g);
    if (!t) FAIL(t.error().message);
    std::mt19937 rng(5);
    check_general(g, *t, samples_for(g, rng, 2000));
    CHECK(std::ranges::count(t->polygons, 1, &RegionPolygon::region) >= 2);
    g.shapes.push_back({{{1, 1, 0}, {7, 3, 0}, {3, 7, 0}}, 2});
    const auto t2 = layered_partition2d(g);
    if (!t2) FAIL(t2.error().message);
    check_general(g, *t2, samples_for(g, rng, 2000));
}

TEST_CASE("mesh layered: 1000 seeded random slanted layouts tile (GD10, GE2)", "[mesh][layered]") {
    std::mt19937 rng(26014);
    double worst_ms = 0.0;
    for (int trial = 0; trial < 1000; ++trial) {
        CAPTURE(trial);
        const int n = 12;
        LayeredGeometry2D g{.lo = {0, 0, 0}, .hi = {n, n, 0}, .background = 0, .shapes = {}, .patches = {}};
        const int count = std::uniform_int_distribution<int>(1, 20)(rng);
        std::uniform_int_distribution<int> at(-2, n + 2), kind(0, 2);
        while (static_cast<int>(g.shapes.size()) < count) {
            const std::int32_t region = std::uniform_int_distribution<std::int32_t>(0, 4)(rng);
            int x0 = at(rng), x1 = at(rng), y0 = at(rng), y1 = at(rng);
            if (x0 == x1 || y0 == y1) continue;
            if (x0 > x1) std::swap(x0, x1);
            if (y0 > y1) std::swap(y0, y1);
            const int k = kind(rng);
            if (k == 0) {
                g.shapes.push_back(box(x0, x1, y0, y1, region));
            } else if (k == 1) {  // a trapezoid: the top inset by a random amount on each side (a triangle if they meet)
                const int w = x1 - x0, a = std::uniform_int_distribution<int>(0, w)(rng), b = std::uniform_int_distribution<int>(0, w - a)(rng);
                std::vector<Point> v{{double(x0), double(y0), 0}, {double(x1), double(y0), 0}};
                if (x1 - b != x0 + a) v.push_back({double(x1 - b), double(y1), 0});
                v.push_back({double(x0 + a), double(y1), 0});
                if (v.size() == 3 && detail::polygon_turn(v) == 0) continue;
                g.shapes.push_back({v, region});
            } else {  // a random triangle
                const Point p{double(at(rng)), double(at(rng)), 0}, q{double(at(rng)), double(at(rng)), 0}, r{double(at(rng)), double(at(rng)), 0};
                std::vector<Point> v{p, q, r};
                if (detail::polygon_turn(v) == 0) continue;
                g.shapes.push_back({v, region});
            }
        }
        if (trial >= 700) {  // the last 300 in decimal cm (0.1 µm steps): no coordinate exact in binary
            const auto scale = [](Point& p) { p = {p[0] * 1e-5, p[1] * 1e-5, 0}; };
            scale(g.hi);
            for (LayeredShape& s : g.shapes) std::ranges::for_each(s.vertices, scale);
        }
        const auto t0 = std::chrono::steady_clock::now();
        const auto t = layered_partition2d(g);
        worst_ms = std::max(worst_ms, std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
        if (!t) {
            std::string d;
            for (const auto& s : g.shapes) {
                d += std::to_string(s.region) + ":";
                for (const auto& p : s.vertices) d += std::format(" {} {},", p[0], p[1]);
                d += " | ";
            }
            UNSCOPED_INFO(d);
        }
        if (!t) FAIL(t.error().message);
        check_general(g, *t, samples_for(g, rng, 200));
    }
    WARN(std::format("GD10: 1000 random slanted layouts, slowest {:.2f} ms", worst_ms));
}

TEST_CASE("mesh layered: 200 trapezoids resolve (GD10, timed)", "[mesh][layered]") {
    std::mt19937 rng(8);
    LayeredGeometry2D g{.lo = {0, 0, 0}, .hi = {1000, 1000, 0}, .background = 0, .shapes = {}, .patches = {}};
    std::uniform_int_distribution<int> at(0, 1000);
    while (g.shapes.size() < 200) {
        int x0 = at(rng), x1 = at(rng), y0 = at(rng), y1 = at(rng);
        if (std::abs(x0 - x1) < 4 || y0 == y1) continue;
        if (x0 > x1) std::swap(x0, x1);
        if (y0 > y1) std::swap(y0, y1);
        const int inset = (x1 - x0) / 4;
        g.shapes.push_back({{{double(x0), double(y0), 0}, {double(x1), double(y0), 0}, {double(x1 - inset), double(y1), 0},
                             {double(x0 + inset), double(y1), 0}},
                            static_cast<std::int32_t>(g.shapes.size() % 5)});
    }
    std::vector<double> ms;
    std::size_t polygons = 0;
    for (int k = 0; k < 3; ++k) {
        const auto t0 = std::chrono::steady_clock::now();
        const auto t = layered_partition2d(g);
        ms.push_back(std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
        if (!t) FAIL(t.error().message);
        polygons = t->polygons.size();
        if (k == 0) {
            std::mt19937 r2(9);
            check_general(g, *t, samples_for(g, r2, 500));
        }
    }
    std::ranges::sort(ms);
    WARN(std::format("GD10: 200 trapezoids resolved in {:.1f} ms (median of 3), {} polygons", ms[1], polygons));
}

TEST_CASE("mesh layered: a vertex meant to lie on a slanted edge, rounded off it, leaves no sliver (GD10)", "[mesh][layered]") {
    // Decimal µm in cm: a vertex on the line from (0, 0.7) to (2, 0.3) at x = 0.1 ... 1.9 is in binary a little above
    // or below it; both sides are taken (above, the sliver between the two shapes would be left uncovered).
    const double u = 1e-4;
    int above = 0, below = 0;
    for (int k = 1; k < 20; ++k) {
        const Point v{0.1 * k * u, (0.7 - 0.02 * k) * u, 0};
        const int side = predicates::orient2d(Point{2 * u, 0.3 * u, 0}, Point{0, 0.7 * u, 0}, v);
        if (side == 0) continue;
        (side > 0 ? above : below) += 1;
        CAPTURE(k, side);
        LayeredGeometry2D g{.lo = {0, 0, 0}, .hi = {2 * u, 1 * u, 0}, .background = 0, .shapes = {}, .patches = {}};
        g.shapes.push_back({{{0, 0, 0}, {2 * u, 0, 0}, {2 * u, 0.3 * u, 0}, {0, 0.7 * u, 0}}, 1});
        g.shapes.push_back({{{0, 0.7 * u, 0}, v, {2 * u, 0.3 * u, 0}, {2 * u, 1 * u, 0}, {0, 1 * u, 0}}, 2});
        const auto t = layered_partition2d(g);
        if (!t) FAIL(t.error().message);
        CHECK(t->polygons.size() == 2);
        std::mt19937 rng(11);
        check_general(g, *t, samples_for(g, rng, 200));
    }
    CHECK(above > 0);
    CHECK(below > 0);
}
