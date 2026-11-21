// Unit 26l slice 4 (SD6): layered 3D solids with slanted prism faces, all parallel to one axis (SQ4).

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

// A volume of the slanted path is a polygon swept along `axis` between two planes: its cross-section (a facet lying in
// one plane across the axis) and its range along the axis.
struct Swept {
    std::vector<Point> section;  // (u, v, 0)
    double lo = 0, hi = 0;
};

std::vector<Swept> swept_volumes(const PolyhedralGeometry& c, int axis) {
    const auto a = static_cast<std::size_t>(axis), ua = (a + 1) % 3, va = (a + 2) % 3;
    std::vector<Swept> out;
    for (const ComplexVolume& v : c.volumes) {
        Swept s{{}, 1e300, -1e300};
        for (const std::size_t f : v.facets) {
            const auto& w = c.facets[f].vertices;
            for (const std::size_t i : w) {
                s.lo = std::min(s.lo, c.vertices[i][a]);
                s.hi = std::max(s.hi, c.vertices[i][a]);
            }
            const bool cap = std::ranges::all_of(w, [&](std::size_t i) { return c.vertices[i][a] == c.vertices[w[0]][a]; });
            if (cap && s.section.empty()) {
                for (const std::size_t i : w) s.section.push_back({c.vertices[i][ua], c.vertices[i][va], 0});
            }
        }
        out.push_back(std::move(s));
    }
    return out;
}

// The checks: the complex valid; its volumes fill the domain; each sample point in exactly one volume, of the topmost
// solid's region. Returns the volume per region.
std::map<std::int32_t, double> check_slanted(const LayeredGeometry3D& g, const PolyhedralGeometry& c, int axis, const std::vector<Point>& samples) {
    const auto report = check_complex(c);
    if (!report) FAIL(report.error().message);
    std::map<std::int32_t, double> got;
    double total = 0.0;
    for (std::size_t v = 0; v < c.volumes.size(); ++v) {
        got[c.volumes[v].region] += report->volume[v];
        total += report->volume[v];
    }
    const double domain = (g.hi[0] - g.lo[0]) * (g.hi[1] - g.lo[1]) * (g.hi[2] - g.lo[2]);
    CHECK(std::abs(total - domain) <= 1e-12 * domain);
    const auto vols = swept_volumes(c, axis);
    const auto a = static_cast<std::size_t>(axis), ua = (a + 1) % 3, va = (a + 2) % 3;
    for (const Point& p : samples) {
        int hits = 0;
        std::int32_t region = -1;
        for (std::size_t v = 0; v < vols.size(); ++v) {
            if (!(p[a] > vols[v].lo && p[a] < vols[v].hi)) continue;
            if (!detail::strictly_inside_2d(vols[v].section, Point{p[ua], p[va], 0})) continue;
            ++hits;
            region = c.volumes[v].region;
        }
        CAPTURE(p[0], p[1], p[2]);
        CHECK(hits == 1);
        CHECK(region == expected_region(g, p));
    }
    return got;
}

// Newell's area of a planar facet.
double area(const PolyhedralGeometry& c, const ComplexFacet& f) {
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

std::map<std::string, double> patch_areas(const PolyhedralGeometry& c) {
    std::vector<int> owners(c.facets.size(), 0);
    for (const ComplexVolume& v : c.volumes) {
        for (const std::size_t f : v.facets) ++owners[f];
    }
    std::map<std::string, double> out;
    for (std::size_t f = 0; f < c.facets.size(); ++f) {
        if (owners[f] == 1) out[c.facets[f].patch] += area(c, c.facets[f]);
    }
    return out;
}

std::vector<Point> random_points(std::mt19937& rng, const LayeredGeometry3D& g, std::size_t n) {
    std::vector<Point> out;
    std::array<std::uniform_real_distribution<double>, 3> d{std::uniform_real_distribution<double>(g.lo[0], g.hi[0]),
                                                            std::uniform_real_distribution<double>(g.lo[1], g.hi[1]),
                                                            std::uniform_real_distribution<double>(g.lo[2], g.hi[2])};
    while (out.size() < n) out.push_back({d[0](rng), d[1](rng), d[2](rng)});
    return out;
}

// A trapezoidal fin along y over a substrate, a gate box wrapping it, oxide around: (u, v) = (z, x) along y.
LayeredGeometry3D finfet() {
    LayeredGeometry3D g{.lo = {0, 0, 0}, .hi = {1, 2, 1}, .background = 0, .solids = {}, .patches = {}};
    g.solids.push_back(box3(0, 1, 0, 2, 0, 0.3, 1));                                                // substrate
    g.solids.push_back(box3(0, 1, 0.5, 1.5, 0.5, 1, 2));                                            // gate
    g.solids.push_back(prism(1, 0, 2, {{0.3, 0.3, 0}, {0.3, 0.7, 0}, {0.8, 0.6, 0}, {0.8, 0.4, 0}}, 1));  // the fin
    g.patches = {{"source", {0, 0, 0}, {1, 0, 1}}, {"gate", {0, 0.5, 1}, {1, 1.5, 1}}, {"body", {0, 0, 0}, {1, 2, 0}}};
    return g;
}

}  // namespace

TEST_CASE("mesh layered3d slanted: a trapezoidal fin under a wrapping gate resolves exactly (SD6)", "[mesh][layered3d]") {
    const LayeredGeometry3D g = finfet();
    const auto c = layered_partition3d(g);
    if (!c) FAIL(c.error().message);
    std::mt19937 rng(64);
    const auto got = check_slanted(g, *c, 1, random_points(rng, g, 500));
    // The fin's cross-section: width 0.4 at z = 0.3 to 0.2 at z = 0.8; 0.078 of it lies above z = 0.5 (in the gate).
    const std::map<std::int32_t, double> want{{0, 2.0 - 0.6 - 0.3 - 0.422}, {1, 0.6 + 0.3}, {2, 0.5 - 0.078}};
    REQUIRE(got.size() == want.size());
    for (const auto& [r, v] : want) {
        CAPTURE(r);
        CHECK(std::abs(got.at(r) - v) <= 1e-12);
    }
    const auto areas = patch_areas(*c);
    CHECK(std::abs(areas.at("source") - 1.0) <= 1e-12);
    CHECK(std::abs(areas.at("gate") - 1.0) <= 1e-12);
    CHECK(std::abs(areas.at("body") - 2.0) <= 1e-12);
    CHECK(std::abs(areas.at("y_max") - 1.0) <= 1e-12);
    CHECK(std::abs(areas.at("z_max") - 1.0) <= 1e-12);
    CHECK(std::abs(areas.at("x_min") - 2.0) <= 1e-12);
    // The slanted sidewalls are facets parallel to y and to neither other axis.
    int slanted = 0;
    for (const ComplexFacet& f : c->facets) {
        const Point& p = c->vertices[f.vertices[0]];
        const bool any_x = std::ranges::any_of(f.vertices, [&](std::size_t i) { return c->vertices[i][0] != p[0]; });
        const bool any_z = std::ranges::any_of(f.vertices, [&](std::size_t i) { return c->vertices[i][2] != p[2]; });
        const bool any_y = std::ranges::any_of(f.vertices, [&](std::size_t i) { return c->vertices[i][1] != p[1]; });
        if (any_x && any_z && any_y) ++slanted;
    }
    CHECK(slanted >= 6);  // two sidewalls through three slabs at least
}

TEST_CASE("mesh layered3d slanted: an end-plane patch box over the fin and a side patch box split exactly (SD6)", "[mesh][layered3d]") {
    LayeredGeometry3D g = finfet();
    // The drain on y_max over the fin's upper part only (the box crosses the slanted sidewalls), and a strip on x_max.
    g.patches = {{"drain", {0.2, 2, 0.5}, {0.8, 2, 0.8}}, {"strip", {1, 0.25, 0.2}, {1, 0.75, 0.6}}};
    const auto c = layered_partition3d(g);
    if (!c) FAIL(c.error().message);
    std::mt19937 rng(65);
    (void)check_slanted(g, *c, 1, random_points(rng, g, 200));
    const auto areas = patch_areas(*c);
    CHECK(std::abs(areas.at("drain") - 0.3 * 0.6) <= 1e-12);
    CHECK(std::abs(areas.at("strip") - 0.2) <= 1e-12);
    CHECK(std::abs(areas.at("x_max") - 1.8) <= 1e-12);
}

TEST_CASE("mesh layered3d slanted: a trapezoid, a triangle and decimal coordinates touching slanted edges (SD6)", "[mesh][layered3d]") {
    // Along z, in (x, y). A box whose corner lies on the trapezoid's slanted edge (decimal, rounded off it), a triangle
    // crossing it, and rectilinear prisms along x and y crossing the slabs.
    LayeredGeometry3D g{.lo = {0, 0, 0}, .hi = {1, 1, 1}, .background = 0, .solids = {}, .patches = {}};
    g.solids.push_back(prism(2, 0.1, 0.9, {{0.1, 0.1, 0}, {0.9, 0.1, 0}, {0.7, 0.7, 0}, {0.3, 0.7, 0}}, 1));
    g.solids.push_back(box3(0.4, 0.6, 0.4, 0.8, 0.2, 0.6, 2));
    g.solids.push_back(box3(0.2, 0.3, 0.4, 0.55, 0.0, 1.0, 3));  // its corner (0.2, 0.4) is on the edge x = 0.1 + y / 3
    g.solids.push_back(prism(2, 0.3, 0.7, {{0.0, 0.0, 0}, {0.5, 0.3, 0}, {0.2, 0.9, 0}}, 4));
    g.solids.push_back(prism(0, 0.0, 0.5, {{0.2, 0.2, 0}, {0.6, 0.2, 0}, {0.6, 0.4, 0}, {0.4, 0.4, 0}, {0.4, 0.8, 0}, {0.2, 0.8, 0}}, 2));
    g.solids.push_back(prism(1, 0.5, 1.0, {{0.7, 0.6, 0}, {0.9, 0.6, 0}, {0.9, 0.9, 0}, {0.7, 0.9, 0}}, 3));
    const auto c = layered_partition3d(g);
    if (!c) FAIL(c.error().message);
    std::mt19937 rng(66);
    (void)check_slanted(g, *c, 2, random_points(rng, g, 1000));
}

TEST_CASE("mesh layered3d slanted: refusals and the uncovered message (SD6)", "[mesh][layered3d]") {
    LayeredGeometry3D g{.lo = {0, 0, 0}, .hi = {4, 4, 4}, .background = -1, .solids = {}, .patches = {}};
    g.solids = {box3(0, 4, 0, 4, 0, 4, 0), prism(2, 0, 1, {{0, 0, 0}, {2, 0, 0}, {0, 2, 0}}, 1), prism(1, 0, 1, {{0, 0, 0}, {2, 0, 0}, {0, 2, 0}}, 1)};
    auto c = layered_partition3d(g);
    REQUIRE_FALSE(c.has_value());
    CHECK(c.error().message ==
          "slanted prisms along two axes (every slanted face must be parallel to one axis, the one all slanted prisms share)");
    CHECK(c.error().context->index == 2);
    g.solids = {box3(0, 4, 0, 4, 0, 3, 0), prism(2, 0, 4, {{0, 0, 0}, {2, 0, 0}, {0, 2, 0}}, 1)};
    c = layered_partition3d(g);
    REQUIRE_FALSE(c.has_value());
    CHECK(c.error().message.starts_with("the layers leave part of the domain uncovered (near ("));
    CHECK(c.error().message.ends_with(", 3.500000)) and no region is the background"));
    // A patch box selecting nothing.
    g.background = 0;
    g.patches = {{"lost", {1, 1, 1}, {2, 2, 2}}};
    c = layered_partition3d(g);
    REQUIRE_FALSE(c.has_value());
    CHECK(c.error().message == "a patch box selects no boundary facet (none inside it, or all taken by an earlier box)");
    CHECK(c.error().context->index == 0);
}

TEST_CASE("mesh layered3d slanted: 1000 seeded random layouts resolve exactly (SD6)", "[mesh][layered3d]") {
    std::mt19937 rng(26046);
    double slowest = 0.0;
    for (int trial = 0; trial < 1000; trial += trial_stride) {
        CAPTURE(trial);
        const int n = 8;
        const bool decimal = trial >= 700;  // the last 300 in tenths, so no coordinate but the integers is exact
        const double q = decimal ? 0.1 : 1.0;
        LayeredGeometry3D g{.lo = {0, 0, 0}, .hi = {n * q, n * q, n * q}, .background = 0, .solids = {}, .patches = {}};
        const int slant = trial % 3;  // the slanted prisms' axis
        const int count = std::uniform_int_distribution<int>(1, 12)(rng);
        std::uniform_int_distribution<int> at(-2, n + 2), kind(0, 3), axis(0, 2), region(0, 4);
        bool any_slanted = false;
        while (static_cast<int>(g.solids.size()) < count || !any_slanted) {
            std::array<int, 6> c{};
            for (int& x : c) x = at(rng);
            if (c[0] == c[1] || c[2] == c[3] || c[4] == c[5]) continue;
            for (int k = 0; k < 6; k += 2) {
                if (c[static_cast<std::size_t>(k)] > c[static_cast<std::size_t>(k + 1)]) std::swap(c[static_cast<std::size_t>(k)], c[static_cast<std::size_t>(k + 1)]);
            }
            const auto s = [&](int i) { return c[static_cast<std::size_t>(i)] * q; };
            const int k = kind(rng);
            if (k == 0) {
                g.solids.push_back(box3(s(0), s(1), s(2), s(3), s(4), s(5), region(rng)));
                continue;
            }
            if (k == 1) {  // a rectilinear L on any axis
                if (c[1] - c[0] < 2 || c[3] - c[2] < 2) continue;
                const double um = s(0) + q, vm = s(2) + q;
                g.solids.push_back(prism(axis(rng), s(4), s(5), {{s(0), s(2), 0}, {s(1), s(2), 0}, {s(1), vm, 0}, {um, vm, 0}, {um, s(3), 0}, {s(0), s(3), 0}}, region(rng)));
                continue;
            }
            // A trapezoid (top inset by one step each side when wide enough, else a triangle) along the slant axis.
            std::vector<Point> p;
            if (c[1] - c[0] >= 3 && k == 2) p = {{s(0), s(2), 0}, {s(1), s(2), 0}, {s(1) - q, s(3), 0}, {s(0) + q, s(3), 0}};
            else p = {{s(0), s(2), 0}, {s(1), s(2), 0}, {s(0), s(3), 0}};
            g.solids.push_back(prism(slant, s(4), s(5), std::move(p), region(rng)));
            any_slanted = true;
        }
        if (trial % 3 == 0) {
            const double top = n * q;
            g.patches = {{"end", {q, q, top}, {5 * q, 4 * q, top}}, {"side", {0, 0, 0}, {0, top, 3 * q}}};
            if (slant == 0) g.patches[0] = {"end", {top, q, q}, {top, 5 * q, 4 * q}};
        }
        const auto t0 = std::chrono::steady_clock::now();
        const auto c = layered_partition3d(g);
        slowest = std::max(slowest, std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count());
        if (!c) FAIL(c.error().message);
        (void)check_slanted(g, *c, slant, random_points(rng, g, 200));
    }
    WARN(std::format("SD6: the slowest of 1000 random slanted layouts resolved in {:.1f} ms", slowest));
}

namespace {

// A layout in µm, as cm.
LayeredGeometry3D in_um(LayeredGeometry3D g) {
    const double u = 1e-4;
    for (std::size_t k = 0; k < 3; ++k) {
        g.lo[k] *= u;
        g.hi[k] *= u;
    }
    for (LayeredSolid& s : g.solids) {
        for (std::size_t k = 0; k < 3; ++k) {
            s.lo[k] *= u;
            s.hi[k] *= u;
        }
        s.from *= u;
        s.to *= u;
        for (Point& p : s.polygon) p = {p[0] * u, p[1] * u, 0};
    }
    for (PatchBox& p : g.patches) {
        for (std::size_t k = 0; k < 3; ++k) {
            p.lo[k] *= u;
            p.hi[k] *= u;
        }
    }
    return g;
}

MeshOptions3D fin_options() {
    MeshOptions3D o;
    o.max_size = 0.15e-4;
    o.interface_resolution = 0.02e-4;
    return o;
}

}  // namespace

TEST_CASE("mesh layered3d slanted: a fin whose sidewalls no plane crosses meshes through the generator (SD6, SQ4)", "[mesh][layered3d]") {
    // The gate reaches down to the fin's base, so nothing crosses the sidewalls: every vertex on them is an input point.
    LayeredGeometry3D g = finfet();
    g.solids[1] = box3(0, 1, 0.5, 1.5, 0.3, 1, 2);
    const auto c = layered_partition3d(in_um(g));
    if (!c) FAIL(c.error().message);
    const auto m = generate_mesh(*c, fin_options());
    if (!m) FAIL(m.error().message);
    CHECK(m->simplices.points.size() > 100);
    WARN(std::format("SD6: the uncrossed trapezoidal fin meshed: {} nodes", m->simplices.points.size()));
}

TEST_CASE("mesh layered3d slanted: a sidewall crossed by a plane folds at a rounded point, and the generator refuses it (SQ4, recorded)",
          "[mesh][layered3d]") {
    // The gate's bottom (z = 0.5) crosses the sidewalls at a computed point: the sidewall's pieces above and below it
    // meet at a fold of about an ulp, parallel to rounding but not exactly coplanar, so the generator's separation rule
    // (SQ6: kept) measures them as two interfaces at distance zero.
    const auto c = layered_partition3d(in_um(finfet()));
    if (!c) FAIL(c.error().message);
    const auto m = generate_mesh(*c, fin_options());
    REQUIRE_FALSE(m.has_value());
    CHECK(m.error().message == "two interfaces within 4 s (a layer thinner than 4 s, or unrelated interfaces)");
}

TEST_CASE("mesh layered3d slanted: one coordinate written two ways is one, on u, v and the patch boxes (SD6)", "[mesh][layered3d]") {
    // 0.3 and 0.1 * 3 differ in the last place: a box's face and a slanted prism's vertex meant to meet, and a side patch
    // box's edge, are one line, not a sliver.
    const double a = 0.3, b = 0.1 * 3;
    REQUIRE(a != b);
    LayeredGeometry3D g{.lo = {0, 0, 0}, .hi = {1, 1, 1}, .background = 0, .solids = {}, .patches = {}};
    g.solids.push_back(box3(0, a, 0, 1, 0, 1, 1));
    g.solids.push_back(prism(2, 0, 1, {{b, 0.2, 0}, {0.8, 0.2, 0}, {0.6, 0.7, 0}}, 2));
    g.patches = {{"strip", {0, 0, 0}, {b, 0, 1}}};
    const auto c = layered_partition3d(g);
    if (!c) FAIL(c.error().message);
    std::mt19937 rng(67);
    (void)check_slanted(g, *c, 2, random_points(rng, g, 300));
    for (const Point& p : c->vertices) CHECK(p[0] != b);  // every such vertex at 0.3 itself
    CHECK(std::abs(patch_areas(*c).at("strip") - 0.3) <= 1e-12);
}
