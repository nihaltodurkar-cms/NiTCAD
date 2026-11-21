// Unit 26l slice 0 (SD2, SE2, SE5): layered 3D solids resolved into a polyhedral complex.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <format>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "NiTCAD/mesh/complex.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "polygon2d.hpp"

using namespace NiTCAD::mesh;

namespace {

// The unoptimised Debug build runs every 20th seeded layout (both halves of each sweep); Release runs them all.
#ifdef NDEBUG
constexpr int trial_stride = 1;
#else
constexpr int trial_stride = 20;
#endif

LayeredSolid box3(double x0, double x1, double y0, double y1, double z0, double z1, std::int32_t region) {
    LayeredSolid s;
    s.kind = LayeredSolid::Kind::box;
    s.lo = {x0, y0, z0};
    s.hi = {x1, y1, z1};
    s.region = region;
    return s;
}

LayeredSolid prism(int axis, double from, double to, std::vector<Point> polygon, std::int32_t region) {
    LayeredSolid s;
    s.kind = LayeredSolid::Kind::prism;
    s.axis = axis;
    s.from = from;
    s.to = to;
    s.polygon = std::move(polygon);
    s.region = region;
    return s;
}

// Whether `p` is strictly inside solid `s` (exact).
bool inside(const LayeredSolid& s, const Point& p) {
    if (s.kind == LayeredSolid::Kind::box) {
        for (std::size_t k = 0; k < 3; ++k) {
            if (!(p[k] > s.lo[k] && p[k] < s.hi[k])) return false;
        }
        return true;
    }
    const auto a = static_cast<std::size_t>(s.axis), ua = (a + 1) % 3, va = (a + 2) % 3;
    if (!(p[a] > s.from && p[a] < s.to)) return false;
    return detail::strictly_inside_2d(s.polygon, Point{p[ua], p[va], 0});
}

std::int32_t expected_region(const LayeredGeometry3D& g, const Point& p) {
    for (std::size_t s = g.solids.size(); s-- > 0;) {
        if (inside(g.solids[s], p)) return g.solids[s].region;
    }
    return g.background;
}

struct Bounds {
    Point lo, hi;
};

// Each volume's bounding box (the kernel's volumes are boxes).
std::vector<Bounds> volume_boxes(const PolyhedralGeometry& c) {
    std::vector<Bounds> out;
    for (const ComplexVolume& v : c.volumes) {
        Bounds b{{1e300, 1e300, 1e300}, {-1e300, -1e300, -1e300}};
        for (const std::size_t f : v.facets) {
            for (const std::size_t i : c.facets[f].vertices) {
                for (std::size_t k = 0; k < 3; ++k) {
                    b.lo[k] = std::min(b.lo[k], c.vertices[i][k]);
                    b.hi[k] = std::max(b.hi[k], c.vertices[i][k]);
                }
            }
        }
        out.push_back(b);
    }
    return out;
}

// SE2's checks: the complex is valid (the kernel ran check_complex; run again here); each volume is a box exactly
// covered by its facets; every sample point lies in exactly one volume, of the topmost solid's region; the volumes per
// region add up to `volumes` when given.
void check_layout(const LayeredGeometry3D& g, const PolyhedralGeometry& c, const std::vector<Point>& samples,
                  const std::map<std::int32_t, double>* volumes = nullptr) {
    REQUIRE(check_complex(c).has_value());
    const auto boxes = volume_boxes(c);
    std::map<std::int32_t, double> got;
    double total = 0.0;
    for (std::size_t v = 0; v < boxes.size(); ++v) {
        const double w = (boxes[v].hi[0] - boxes[v].lo[0]) * (boxes[v].hi[1] - boxes[v].lo[1]) * (boxes[v].hi[2] - boxes[v].lo[2]);
        got[c.volumes[v].region] += w;
        total += w;
    }
    const double domain = (g.hi[0] - g.lo[0]) * (g.hi[1] - g.lo[1]) * (g.hi[2] - g.lo[2]);
    CHECK(std::abs(total - domain) <= 1e-12 * domain);
    if (volumes != nullptr) CHECK(got == *volumes);
    for (const Point& p : samples) {
        int hits = 0;
        std::int32_t region = -1;
        for (std::size_t v = 0; v < boxes.size(); ++v) {
            bool in = true;
            for (std::size_t k = 0; k < 3 && in; ++k) in = p[k] > boxes[v].lo[k] && p[k] < boxes[v].hi[k];
            if (in) {
                ++hits;
                region = c.volumes[v].region;
            }
        }
        CAPTURE(p[0], p[1], p[2]);
        CHECK(hits == 1);
        CHECK(region == expected_region(g, p));
    }
}

double facet_area(const PolyhedralGeometry& c, const ComplexFacet& f) {
    Point lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (const std::size_t i : f.vertices) {
        for (std::size_t k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], c.vertices[i][k]);
            hi[k] = std::max(hi[k], c.vertices[i][k]);
        }
    }
    double a = 1.0;
    for (std::size_t k = 0; k < 3; ++k) {
        if (hi[k] > lo[k]) a *= hi[k] - lo[k];
    }
    return a;
}

}  // namespace

TEST_CASE("mesh layered3d: a domain of the background is one box volume with its six faces named (SD2)", "[mesh][layered3d]") {
    LayeredGeometry3D g{.lo = {0, 0, 0}, .hi = {2, 3, 4}, .background = 0, .solids = {}, .patches = {}};
    const auto c = layered_partition3d(g);
    if (!c) FAIL(c.error().message);
    CHECK(c->volumes.size() == 1);
    CHECK(c->facets.size() == 6);
    CHECK(c->vertices.size() == 8);
    std::vector<std::string> names;
    for (const ComplexFacet& f : c->facets) names.push_back(f.patch);
    std::ranges::sort(names);
    CHECK(names == std::vector<std::string>{"x_max", "x_min", "y_max", "y_min", "z_max", "z_min"});
}

TEST_CASE("mesh layered3d: buried boxes, a prism and T-junctions resolve exactly (SD2)", "[mesh][layered3d]") {
    LayeredGeometry3D g{.lo = {0, 0, 0}, .hi = {6, 6, 6}, .background = 0, .solids = {}, .patches = {}};
    g.solids.push_back(box3(1, 5, 1, 5, 1, 5, 1));
    g.solids.push_back(box3(2, 3, 2, 3, 2, 3, 2));                                         // buried inside the first
    g.solids.push_back(prism(0, -1, 7, {{1, 4, 0}, {3, 4, 0}, {3, 5, 0}, {2, 5, 0}, {2, 6, 0}, {1, 6, 0}}, 3));  // an L along x, through
    g.solids.push_back(box3(4, 8, 0, 1, 0, 2, 0));                                         // the background again, past the domain
    const auto c = layered_partition3d(g);
    if (!c) FAIL(c.error().message);
    std::mt19937 rng(1);
    std::uniform_real_distribution<double> u(0.0, 6.0);
    std::vector<Point> samples;
    while (samples.size() < 3000) {
        const Point p{u(rng), u(rng), u(rng)};
        if (std::ranges::any_of(p, [](double x) { return x == std::floor(x); })) continue;
        samples.push_back(p);
    }
    // Exact volumes by unit cells (every coordinate an integer).
    std::map<std::int32_t, double> volumes;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            for (int k = 0; k < 6; ++k) volumes[expected_region(g, {i + 0.5, j + 0.5, k + 0.5})] += 1.0;
        }
    }
    check_layout(g, *c, samples, &volumes);
}

TEST_CASE("mesh layered3d: patch boxes split the boundary faces and name exactly what they hold (SD1)", "[mesh][layered3d]") {
    const double u = 1e-4;
    LayeredGeometry3D g{.lo = {0, 0, 0}, .hi = {2 * u, 1 * u, 1 * u}, .background = 0, .solids = {}, .patches = {}};
    g.solids.push_back(box3(0.5 * u, 1.5 * u, 0, 1 * u, 0.9 * u, 1 * u, 1));
    g.patches.push_back({"gate", {0.6 * u, 0.2 * u, 1 * u}, {1.4 * u, 0.8 * u, 1 * u}});
    g.patches.push_back({"source", {0, 0, 0}, {0, 1 * u, 0.5 * u}});
    const auto c = layered_partition3d(g);
    if (!c) FAIL(c.error().message);
    std::map<std::string, double> area;
    for (const ComplexFacet& f : c->facets) {
        if (!f.patch.empty()) area[f.patch] += facet_area(*c, f);
    }
    CHECK(std::abs(area["gate"] - 0.8 * u * 0.6 * u) <= 1e-12 * u * u);
    CHECK(std::abs(area["source"] - 1 * u * 0.5 * u) <= 1e-12 * u * u);
    CHECK(std::abs(area["x_min"] - 1 * u * 0.5 * u) <= 1e-12 * u * u);
    CHECK(std::abs(area["z_max"] - (2 * u * u - 0.48 * u * u)) <= 1e-12 * u * u);
    // 20nm and 0.02um are one plane: a patch written with other doubles still selects.
    g.patches[0].lo[2] = 1e-4 * (1.0 - 1e-15);
    CHECK(layered_partition3d(g).has_value());
    g.patches.push_back({"nothing", {0.1 * u, 0.1 * u, 0.5 * u}, {0.2 * u, 0.2 * u, 0.6 * u}});  // inside: no facet
    const auto bad = layered_partition3d(g);
    REQUIRE_FALSE(bad.has_value());
    CHECK(bad.error().message.starts_with("a patch box selects no boundary facet"));
    CHECK(bad.error().context->index == 2);
}

TEST_CASE("mesh layered3d: refusals name the solid (SD2)", "[mesh][layered3d]") {
    LayeredGeometry3D g{.lo = {0, 0, 0}, .hi = {4, 4, 4}, .background = 0, .solids = {}, .patches = {}};
    const auto refused = [&](LayeredSolid s, const std::string& message) {
        LayeredGeometry3D h = g;
        h.solids = {prism(2, 0, 1, {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}, 1), std::move(s)};
        const auto c = layered_partition3d(h);
        REQUIRE_FALSE(c.has_value());
        CHECK(c.error().message == message);
        CHECK(c.error().context->index == 1);
    };
    refused(box3(1, 1, 0, 1, 0, 1, 1), "a layered box needs lo < hi on each axis");
    refused(box3(0, 1, 0, 1, 0, 1, -1), "a layered solid's region is >= 0");
    refused(prism(0, 0, 1, {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}, 1), "slanted prisms along two axes (every slanted face must be parallel to one axis, the one all slanted prisms share)");
    refused(prism(2, 0, 1, {{0, 0, 0}, {1, 0, 0}}, 1), "a prism's polygon needs at least three vertices");
    refused(prism(2, 1, 1, {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, 1), "a prism needs from < to along its axis");
    refused(prism(3, 0, 1, {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, 1), "a prism's axis is x, y or z");
    refused(prism(2, 0, 1, {{0, 0, 0}, {2, 0, 0}, {2, 1, 0}, {0, 1, 0}, {0, 0, 0}}, 1), "a prism's polygon repeats a vertex");
    refused(prism(2, 0, 1, {{0, 0, 0}, {2, 0, 0}, {2, 1, 0}, {1, 1, 0}, {1, -1, 0}, {0, -1, 0}}, 1), "a prism's polygon is not simple");
    g.background = -1;
    const auto c = layered_partition3d(LayeredGeometry3D{g.lo, g.hi, -1, {box3(0, 4, 0, 4, 0, 3, 1)}, {}});
    REQUIRE_FALSE(c.has_value());
    CHECK(c.error().message == "the layers leave part of the domain uncovered (near (2.000000, 2.000000, 3.500000)) and no region is the background");
    CHECK(layered_partition3d(LayeredGeometry3D{g.lo, g.hi, -1, {box3(0, 4, 0, 4, 0, 3, 1), box3(-1, 5, -1, 5, 3, 5, 2)}, {}}).has_value());
}

TEST_CASE("mesh layered3d: 1000 seeded random layouts resolve exactly (SE2)", "[mesh][layered3d]") {
    std::mt19937 rng(26015);
    for (int trial = 0; trial < 1000; trial += trial_stride) {
        CAPTURE(trial);
        const int n = 8;
        LayeredGeometry3D g{.lo = {0, 0, 0}, .hi = {n, n, n}, .background = 0, .solids = {}, .patches = {}};
        const int count = std::uniform_int_distribution<int>(1, 20)(rng);
        std::uniform_int_distribution<int> at(-2, n + 2), kind(0, 2), axis(0, 2), region(0, 4);
        while (static_cast<int>(g.solids.size()) < count) {
            std::array<int, 6> c{};
            for (int& x : c) x = at(rng);
            if (c[0] == c[1] || c[2] == c[3] || c[4] == c[5]) continue;
            for (int k = 0; k < 6; k += 2) {
                if (c[static_cast<std::size_t>(k)] > c[static_cast<std::size_t>(k + 1)]) std::swap(c[static_cast<std::size_t>(k)], c[static_cast<std::size_t>(k + 1)]);
            }
            if (kind(rng) != 0) {
                g.solids.push_back(box3(c[0], c[1], c[2], c[3], c[4], c[5], region(rng)));
                continue;
            }
            // A rectilinear prism: an L, or a U when wide enough, in (u, v).
            const double u0 = c[0], u1 = c[1], v0 = c[2], v1 = c[3];
            if (u1 - u0 < 2 || v1 - v0 < 2) continue;
            const double um = u0 + 1, vm = v0 + 1;
            std::vector<Point> l{{u0, v0, 0}, {u1, v0, 0}, {u1, vm, 0}, {um, vm, 0}, {um, v1, 0}, {u0, v1, 0}};
            if (u1 - u0 >= 3 && (trial % 2) == 0) {
                l = {{u0, v0, 0}, {u1, v0, 0}, {u1, v1, 0}, {u1 - 1, v1, 0}, {u1 - 1, vm, 0}, {um, vm, 0}, {um, v1, 0}, {u0, v1, 0}};
            }
            g.solids.push_back(prism(axis(rng), c[4], c[5], std::move(l), region(rng)));
        }
        if (trial % 3 == 0) g.patches = {{"top", {1, 1, n}, {5, 4, n}}, {"left", {0, 0, 0}, {0, n, 3}}};
        const auto c = layered_partition3d(g);
        if (!c) FAIL(c.error().message);
        std::map<std::int32_t, double> volumes;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                for (int k = 0; k < n; ++k) volumes[expected_region(g, {i + 0.5, j + 0.5, k + 0.5})] += 1.0;
            }
        }
        std::uniform_real_distribution<double> u(0.0, n);
        std::vector<Point> samples;
        while (samples.size() < 200) {
            const Point p{u(rng), u(rng), u(rng)};
            if (std::ranges::any_of(p, [](double x) { return x == std::floor(x); })) continue;
            samples.push_back(p);
        }
        check_layout(g, *c, samples, &volumes);
    }
}

TEST_CASE("mesh layered3d: 100 random boxes resolve within 50 ms (SE5)", "[mesh][layered3d]") {
    std::mt19937 rng(7);
    LayeredGeometry3D g{.lo = {0, 0, 0}, .hi = {1000, 1000, 1000}, .background = 0, .solids = {}, .patches = {}};
    std::uniform_int_distribution<int> at(0, 1000);
    while (g.solids.size() < 100) {
        std::array<int, 6> c{};
        for (int& x : c) x = at(rng);
        if (c[0] == c[1] || c[2] == c[3] || c[4] == c[5]) continue;
        for (int k = 0; k < 6; k += 2) {
            if (c[static_cast<std::size_t>(k)] > c[static_cast<std::size_t>(k + 1)]) std::swap(c[static_cast<std::size_t>(k)], c[static_cast<std::size_t>(k + 1)]);
        }
        g.solids.push_back(box3(c[0], c[1], c[2], c[3], c[4], c[5], static_cast<std::int32_t>(g.solids.size() % 5)));
    }
    std::vector<double> ms;
    std::size_t volumes = 0, facets = 0;
    for (int k = 0; k < 3; ++k) {
        const auto t0 = std::chrono::steady_clock::now();
        const auto c = layered_partition3d(g);
        ms.push_back(std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
        if (!c) FAIL(c.error().message);
        volumes = c->volumes.size();
        facets = c->facets.size();
    }
    std::ranges::sort(ms);
    WARN(std::format("SE5: 100 random boxes resolved in {:.1f} ms (median of 3), {} volumes, {} facets", ms[1], volumes, facets));
    // Release only (owner, 2026-11-19, as KG6): the unoptimised Debug build records the time.
#ifdef NDEBUG
    CHECK(ms[1] <= 50.0);
#endif
}

TEST_CASE("mesh layered3d: a layered MOS stack meshes through the polyhedral generator (SD1)", "[mesh][layered3d]") {
    // Silicon, an oxide over its middle, a gate over the oxide: in µm, as cm.
    const double u = 1e-4;
    LayeredGeometry3D g{.lo = {0, 0, 0}, .hi = {1 * u, 0.6 * u, 0.6 * u}, .background = 0, .solids = {}, .patches = {}};
    g.solids.push_back(box3(0, 1 * u, 0, 0.6 * u, 0.4 * u, 0.6 * u, 1));          // oxide over the top
    g.patches.push_back({"gate", {0.3 * u, 0, 0.6 * u}, {0.7 * u, 0.6 * u, 0.6 * u}});
    g.patches.push_back({"body", {0, 0, 0}, {1 * u, 0.6 * u, 0}});
    const auto c = layered_partition3d(g);
    if (!c) FAIL(c.error().message);
    MeshOptions3D o;
    o.max_size = 0.15 * u;
    o.interface_resolution = 0.02 * u;
    const auto m = generate_mesh(*c, o);
    if (!m) FAIL(m.error().message);
    CHECK(m->simplices.points.size() > 100);
}

namespace {

// A planar polygon's area (Newell's normal).
double polygon_area(const PolyhedralGeometry& c, const ComplexFacet& f) {
    Point n{0, 0, 0};
    const std::size_t m = f.vertices.size();
    for (std::size_t k = 0; k < m; ++k) {
        const Point& p = c.vertices[f.vertices[k]];
        const Point& q = c.vertices[f.vertices[(k + 1) % m]];
        n[0] += (p[1] - q[1]) * (p[2] + q[2]);
        n[1] += (p[2] - q[2]) * (p[0] + q[0]);
        n[2] += (p[0] - q[0]) * (p[1] + q[1]);
    }
    return 0.5 * std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
}

struct Measures {
    std::map<std::int32_t, double> volume;
    std::map<std::string, double> area;  // boundary facets per patch
};

Measures measures(const PolyhedralGeometry& c) {
    const auto r = check_complex(c);
    REQUIRE(r.has_value());
    Measures m;
    for (std::size_t v = 0; v < c.volumes.size(); ++v) m.volume[c.volumes[v].region] += r->volume[v];
    std::vector<int> on(c.facets.size(), 0);
    for (const ComplexVolume& v : c.volumes) {
        for (const std::size_t f : v.facets) ++on[f];
    }
    for (std::size_t f = 0; f < c.facets.size(); ++f) {
        if (on[f] == 1) m.area[c.facets[f].patch] += polygon_area(c, c.facets[f]);
    }
    return m;
}

}  // namespace

TEST_CASE("mesh layered3d: the reduced FinFET drawn as layers equals its extrusion's volumes and patch areas (SE2)",
          "[mesh][layered3d]") {
    // decks/finfet_3d_reduced.deck: axis x, stations 0, 30, 55, 85 nm; cells in (y, z).
    const double nm = 1e-7;
    ExtrudedGeometry x;
    x.axis = 0;
    x.stations = {0, 30 * nm, 55 * nm, 85 * nm};
    const auto pts = [&](std::vector<std::pair<double, double>> v) {
        std::vector<Point> p;
        for (const auto& [a, b] : v) p.push_back({a * nm, b * nm, 0});
        return p;
    };
    x.cells = {{"box", pts({{0, 0}, {20, 0}, {20, 1.5}, {20, 2}, {20, 12}, {20, 12.5}, {20, 14}, {0, 14}})},
               {"fin", pts({{20, 2}, {50, 2}, {50, 12}, {20, 12}})},
               {"il", pts({{20, 1.5}, {50.5, 1.5}, {50.5, 12.5}, {20, 12.5}, {20, 12}, {50, 12}, {50, 2}, {20, 2}})},
               {"hk", pts({{20, 0}, {52, 0}, {52, 14}, {20, 14}, {20, 12.5}, {50.5, 12.5}, {50.5, 1.5}, {20, 1.5}})}};
    x.regions = {{0, 0, 0}, {4, 4, 4}, {1, 3, 1}, {1, 2, 1}};  // buried oxide, spacer, hfo2, interfacial, silicon
    const std::vector<PatchBox> patches{{"source", {0, 20 * nm, 2 * nm}, {0, 50 * nm, 12 * nm}},
                                        {"drain", {85 * nm, 20 * nm, 2 * nm}, {85 * nm, 50 * nm, 12 * nm}},
                                        {"gate", {30 * nm, 20 * nm, 0}, {55 * nm, 52 * nm, 14 * nm}},
                                        {"substrate", {0, 0, 0}, {85 * nm, 0, 14 * nm}}};
    x.patches = patches;
    const auto e = extrude(x);
    if (!e) FAIL(e.error().message);
    LayeredGeometry3D g{.lo = {0, 0, 0}, .hi = {85 * nm, 52 * nm, 14 * nm}, .background = 0, .solids = {}, .patches = patches};
    g.solids = {box3(0, 30 * nm, 20 * nm, 52 * nm, 0, 14 * nm, 1),       box3(55 * nm, 85 * nm, 20 * nm, 52 * nm, 0, 14 * nm, 1),
                box3(30 * nm, 55 * nm, 20 * nm, 52 * nm, 0, 14 * nm, 2), box3(30 * nm, 55 * nm, 20 * nm, 50.5 * nm, 1.5 * nm, 12.5 * nm, 3),
                box3(0, 85 * nm, 20 * nm, 50 * nm, 2 * nm, 12 * nm, 4)};
    const auto l = layered_partition3d(g);
    if (!l) FAIL(l.error().message);
    const Measures a = measures(*e), b = measures(*l);
    REQUIRE(a.volume.size() == b.volume.size());
    for (const auto& [r, v] : a.volume) {
        CAPTURE(r, v, b.volume.at(r));
        CHECK(std::abs(v - b.volume.at(r)) <= 1e-12 * v);
    }
    // The extrusion names the rest "boundary", the layers by face: compare the named patches and the total.
    double rest_a = 0.0, rest_b = 0.0;
    for (const auto& [n, v] : a.area) {
        if (b.area.contains(n)) {
            CAPTURE(n, v, b.area.at(n));
            CHECK(std::abs(v - b.area.at(n)) <= 1e-12 * v);
        } else {
            rest_a += v;
        }
    }
    for (const auto& [n, v] : b.area) {
        if (!a.area.contains(n)) rest_b += v;
    }
    CHECK(std::abs(rest_a - rest_b) <= 1e-12 * rest_a);
}

TEST_CASE("mesh layered3d: one plane written two ways (20nm, 0.02um) is one plane, and a patch box so written selects (SD2)",
          "[mesh][layered3d]") {
    // A length in nm and the same in µm, as the deck reads them: n * 1e-7 and (n / 1000) * 1e-4 cm, which can differ in
    // the last place. Take the first n where they do.
    double a = 0.0, b = 0.0;
    for (int n = 1; n < 1000 && a == b; ++n) {
        a = n * 1e-7;
        b = (n / 1000.0) * 1e-4;
    }
    REQUIRE(a != b);
    const double h = 4 * std::max(a, b);
    LayeredGeometry3D g{.lo = {0, 0, 0}, .hi = {h, h, h}, .background = 0, .solids = {}, .patches = {}};
    g.solids.push_back(box3(0, a, 0, h, 0, h, 1));  // its face at a
    g.solids.push_back(box3(b, h, 0, h, 0, h, 2));  // the next one's at b
    g.patches.push_back({"left", {0, 0, h}, {b, h, h}});  // on the top, up to b: holds the top of the first box
    g.patches.push_back({"right", {std::max(a, b), 0, h}, {h, h, h}});  // from the larger: the plane kept is the smaller
    const auto c = layered_partition3d(g);
    if (!c) FAIL(c.error().message);
    std::vector<double> xs;
    for (const Point& p : c->vertices) xs.push_back(p[0]);
    std::ranges::sort(xs);
    xs.erase(std::unique(xs.begin(), xs.end()), xs.end());
    CHECK(xs.size() == 3);  // 0, the shared plane, h: no sliver between a and b
    CHECK(c->volumes.size() == 2);
    double left = 0.0;
    for (const ComplexFacet& f : c->facets) {
        if (f.patch == "left") left += facet_area(*c, f);
    }
    CHECK(std::abs(left - std::min(a, b) * h) <= 1e-12 * h * h);
    double right = 0.0;
    for (const ComplexFacet& f : c->facets) {
        if (f.patch == "right") right += facet_area(*c, f);
    }
    CHECK(std::abs(right - (h - std::min(a, b)) * h) <= 1e-12 * h * h);
}
