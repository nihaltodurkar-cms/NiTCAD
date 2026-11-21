// Unit 64 slice 0 (IP1): the interface pieces the 3D generator refuses, reproduced from layered 3D layouts and pinned
// with the message each gets today; the layouts the proposal expected refused that mesh, pinned as meshing.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <numbers>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/mesh/complex.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "predicates.hpp"

using namespace NiTCAD::mesh;

namespace {

constexpr double um = 1e-4;

LayeredSolid box3(double x0, double x1, double y0, double y1, double z0, double z1, std::int32_t region) {
    LayeredSolid s;
    s.kind = LayeredSolid::Kind::box;
    s.lo = {x0 * um, y0 * um, z0 * um};
    s.hi = {x1 * um, y1 * um, z1 * um};
    s.region = region;
    return s;
}

// A prism along y: the polygon in (z, x), µm.
LayeredSolid prism_y(double from, double to, const std::vector<std::pair<double, double>>& zx, std::int32_t region) {
    LayeredSolid s;
    s.kind = LayeredSolid::Kind::prism;
    s.axis = 1;
    s.from = from * um;
    s.to = to * um;
    for (const auto& [z, x] : zx) s.polygon.push_back({z * um, x * um, 0.0});
    s.region = region;
    return s;
}

// Oxide (0) fills the 1 x 2 x 1 µm domain under the solids.
LayeredGeometry3D domain() { return {.lo = {0, 0, 0}, .hi = {1 * um, 2 * um, 1 * um}, .background = 0, .solids = {}, .patches = {}}; }

MeshOptions3D options(double s_um) {
    MeshOptions3D o;
    o.max_size = 0.15 * um;
    o.interface_resolution = s_um * um;
    return o;
}

// The generator's answer: "" when it meshes, else its message.
std::string outcome(const LayeredGeometry3D& g, double s_um) {
    const auto c = layered_partition3d(g);
    REQUIRE(c.has_value());
    const auto m = generate_mesh(*c, options(s_um));
    return m ? std::string{} : m.error().message;
}

const std::string holes = "an interface piece has holes (unsupported)";
const std::string on_domain_edge = "an interface ends on a domain edge (the boundary must be planar across a trace; unsupported)";
const std::string within_4s = "two interfaces within 4 s (a layer thinner than 4 s, or unrelated interfaces)";
const std::string short_trace = "an interface trace on the boundary is shorter than 6 s";

// A gate (2) and a wire (1) along y through it, ending on the domain's y faces: the gate's faces at y = 0.5 and 1.5
// are gate | oxide pieces with the wire's cross-section as a hole.
LayeredGeometry3D square_wire() {
    LayeredGeometry3D g = domain();
    g.solids.push_back(box3(0.1, 0.9, 0.5, 1.5, 0.1, 0.9, 2));
    g.solids.push_back(box3(0.4, 0.6, 0, 2, 0.4, 0.6, 1));
    return g;
}

LayeredGeometry3D polygon_wire(int sides, double radius) {
    LayeredGeometry3D g = domain();
    g.solids.push_back(box3(0.1, 0.9, 0.5, 1.5, 0.1, 0.9, 2));
    std::vector<std::pair<double, double>> zx;
    for (int k = 0; k < sides; ++k) {
        const double t = 2.0 * std::numbers::pi * k / sides;
        zx.emplace_back(0.5 + radius * std::cos(t), 0.5 + radius * std::sin(t));
    }
    g.solids.push_back(prism_y(0, 2, zx, 1));
    return g;
}

// Silicon (1): a substrate and a fin along y flush with the domain's x = 0 face, so the fin's top meets the domain
// edges x = 0, y = 0 and x = 0, y = 2.
LayeredGeometry3D flush_fin() {
    LayeredGeometry3D g = domain();
    g.solids.push_back(box3(0, 1, 0, 2, 0, 0.3, 1));
    g.solids.push_back(box3(0, 0.3, 0, 2, 0.3, 0.8, 1));
    return g;
}

// The flush fin with its outer top corner rounded: a fillet of radius 0.2 in four segments (each 0.078 long).
LayeredGeometry3D filleted_flush_fin() {
    LayeredGeometry3D g = domain();
    g.solids.push_back(box3(0, 1, 0, 2, 0, 0.3, 1));
    std::vector<std::pair<double, double>> zx{{0.3, 0}, {0.3, 0.3}};
    for (int k = 0; k <= 4; ++k) {
        const double t = 0.5 * std::numbers::pi * k / 4;
        zx.emplace_back(0.6 + 0.2 * std::sin(t), 0.1 + 0.2 * std::cos(t));
    }
    zx.emplace_back(0.8, 0);
    g.solids.push_back(prism_y(0, 2, zx, 1));
    return g;
}

// A gate (2) over the flush fin, flush with the domain's top and both x faces.
LayeredGeometry3D flush_fin_gate() {
    LayeredGeometry3D g = flush_fin();
    g.solids.push_back(box3(0, 1, 0.5, 1.5, 0.5, 1, 2));
    return g;
}

// A block (2) in the domain's corner (x, y low; z high), touching three faces.
LayeredGeometry3D corner_block() {
    LayeredGeometry3D g = domain();
    g.solids.push_back(box3(0, 1, 0, 2, 0, 0.3, 1));
    g.solids.push_back(box3(0, 0.3, 0, 0.3, 0.7, 1, 2));
    return g;
}

// A slanted wall of a prism (2) along y ending on the domain edge x = 0, z = 1: from (z, x) = (1, 0) to (0.4, 0.2),
// the oxide's wedge between it and the x = 0 face is atan(0.2 / 0.6) = 18.4 degrees. Through y = 0 to `to`: at
// to = 1 the wall also reaches the domain corner (0, 0, 1).
LayeredGeometry3D wall_on_edge(double to) {
    LayeredGeometry3D g = domain();
    g.solids.push_back(box3(0, 1, 0, 2, 0, 0.3, 1));
    g.solids.push_back(prism_y(0, to, {{1, 0}, {1, 0.6}, {0.4, 0.6}, {0.4, 0.2}}, 2));
    return g;
}

// Unit 26l SD6's FinFET: a trapezoidal fin (1) along y, the gate (2) from z = 0.5 crossing its sidewalls.
LayeredGeometry3D slanted_fin() {
    LayeredGeometry3D g = domain();
    g.solids.push_back(box3(0, 1, 0, 2, 0, 0.3, 1));
    g.solids.push_back(box3(0, 1, 0.5, 1.5, 0.5, 1, 2));
    g.solids.push_back(prism_y(0, 2, {{0.3, 0.3}, {0.3, 0.7}, {0.8, 0.6}, {0.8, 0.4}}, 1));
    return g;
}

Point newell_unit(const PolyhedralGeometry& c, std::size_t f) {
    const auto& vs = c.facets[f].vertices;
    Point n{0, 0, 0};
    for (std::size_t k = 0; k < vs.size(); ++k) {
        const Point& p = c.vertices[vs[k]];
        const Point& q = c.vertices[vs[(k + 1) % vs.size()]];
        n[0] += (p[1] - q[1]) * (p[2] + q[2]);
        n[1] += (p[2] - q[2]) * (p[0] + q[0]);
        n[2] += (p[0] - q[0]) * (p[1] + q[1]);
    }
    const double l = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    return {n[0] / l, n[1] / l, n[2] / l};
}

}  // namespace

TEST_CASE("mesh interface pieces: a wire through a gate makes holed pieces, refused (IP1)", "[mesh][interface_pieces]") {
    CHECK(outcome(square_wire(), 0.02) == holes);
    // 16 sides of 0.078 µm: s = 0.01 keeps every trace on the y faces over 6 s, so the hole is what is refused.
    CHECK(outcome(polygon_wire(16, 0.2), 0.01) == holes);
}

TEST_CASE("mesh interface pieces: fins and blocks flush with the domain's sides and corners already mesh (IP1, record)",
          "[mesh][interface_pieces]") {
    // An axis-aligned piece in a box domain never has an edge on a domain edge or a corner at a domain corner (its
    // plane would have to hold the edge): where it reaches a domain edge, it crosses it, and that meshes today.
    CHECK(outcome(flush_fin(), 0.02).empty());
    CHECK(outcome(flush_fin_gate(), 0.02).empty());
    CHECK(outcome(corner_block(), 0.02).empty());
    // A fillet meshes once s is at most a sixth of its segments; at s = 0.02 its traces on the y faces are too short.
    CHECK(outcome(filleted_flush_fin(), 0.01).empty());
    CHECK(outcome(filleted_flush_fin(), 0.02) == short_trace);
}

TEST_CASE("mesh interface pieces: a slanted wall ending on a domain edge is refused, and its wedge is under 60 degrees (IP1)",
          "[mesh][interface_pieces]") {
    CHECK(outcome(wall_on_edge(2), 0.02) == on_domain_edge);
    CHECK(outcome(wall_on_edge(1), 0.02) == on_domain_edge);  // reaching the domain corner too: the edge is met first
    // The two faces meet at 90 degrees, so a wall through their edge leaves a wedge of at most 45 degrees on one side:
    // the kept 60-degree rule (U12) refuses it if the domain-edge rule did not.
    CHECK(std::atan2(0.2, 0.6) * 180.0 / std::numbers::pi < 60.0);
}

TEST_CASE("mesh interface pieces: a sidewall crossed by a plane folds into pieces touching, parallel, not coplanar (IP1)",
          "[mesh][interface_pieces]") {
    CHECK(outcome(slanted_fin(), 0.02) == within_4s);
    const auto c = layered_partition3d(slanted_fin());
    REQUIRE(c.has_value());
    std::vector<std::vector<std::int32_t>> regions(c->facets.size());
    for (const ComplexVolume& v : c->volumes) {
        for (const std::size_t f : v.facets) regions[f].push_back(v.region);
    }
    std::vector<std::size_t> walls;  // the interface facets off the axes
    for (std::size_t f = 0; f < c->facets.size(); ++f) {
        if (regions[f].size() != 2 || regions[f][0] == regions[f][1]) continue;
        std::ranges::sort(regions[f]);
        const Point n = newell_unit(*c, f);
        if (std::abs(n[0]) != 1 && std::abs(n[1]) != 1 && std::abs(n[2]) != 1) walls.push_back(f);
    }
    // Every touching pair of wall facets not exactly coplanar: parallel to 1e-15 (the generator's rule takes 1e-12).
    int same_pair = 0, other_pair = 0;
    double worst = 0;
    for (std::size_t i = 0; i < walls.size(); ++i) {
        for (std::size_t j = i + 1; j < walls.size(); ++j) {
            const auto& A = c->facets[walls[i]].vertices;
            const auto& B = c->facets[walls[j]].vertices;
            const std::set<std::size_t> sa(A.begin(), A.end());
            if (std::ranges::none_of(B, [&](std::size_t v) { return sa.contains(v); })) continue;
            const bool coplanar = std::ranges::all_of(B, [&](std::size_t v) {
                return predicates::orient3d(c->vertices[A[0]], c->vertices[A[1]], c->vertices[A[2]], c->vertices[v]) == 0;
            });
            if (coplanar) continue;
            const Point a = newell_unit(*c, walls[i]), b = newell_unit(*c, walls[j]);
            const Point x{a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
            worst = std::max(worst, std::sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]));
            ++(regions[walls[i]] == regions[walls[j]] ? same_pair : other_pair);
        }
    }
    // Both kinds: one sidewall interface split into pieces of its own pair, and the fin | oxide and fin | gate
    // pieces meeting where the gate's bottom crosses (the proposal's "same two regions" was too narrow).
    CHECK(same_pair == 8);
    CHECK(other_pair == 6);
    CHECK(worst <= 1e-15);
}
