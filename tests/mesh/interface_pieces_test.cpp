// Unit 64: the interface pieces the 3D generator refused (slice 0, IP1: reproduced from layered 3D layouts and
// pinned with the message each got; the layouts the proposal expected refused that mesh, pinned as meshing), and
// slice 1 (IP2): holed pieces meshed natively, certified, the rest unchanged to the bit.

#include <catch2/catch_message.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <expected>
#include <format>
#include <map>
#include <numbers>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "NiTCAD/mesh/complex.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/interfaces.hpp"
#include "predicates.hpp"

using namespace NiTCAD::mesh;

namespace {

constexpr double um = 1e-4;

// The unoptimised Debug build runs every 20th seeded layout; Release runs them all.
#ifdef NDEBUG
constexpr int trial_stride = 1;
#else
constexpr int trial_stride = 20;
#endif

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

std::expected<GeneratedMesh, NiTCAD::base::Error> mesh_of(const LayeredGeometry3D& g, double s_um) {
    const auto c = layered_partition3d(g);
    REQUIRE(c.has_value());
    return generate_mesh(*c, options(s_um));
}

// The certificate's checks (as tests/mesh/generate_complex_interfaces_test.cpp's FC6) and the volume [cm^3].
void check_certified(const GeneratedMesh& m, double s, double volume) {
    REQUIRE(m.certificate.has_value());
    const InterfaceReport& r = m.certificate->report;
    CHECK(r.cross_faces > 0);
    CHECK(r.max_flat_deviation <= 1.0);
    CHECK(r.max_corner_deviation <= s);
    CHECK(certifies(*m.certificate, m.mesh, m.node_region));
    double total = 0.0;
    for (const double v : m.mesh.volumes()) total += v;
    CHECK(std::abs(total - volume) <= 1e-12 * volume);
}

// FNV-1a over the points' and cells' bytes.
std::uint64_t digest(const GeneratedMesh& m) {
    std::uint64_t h = 1469598103934665603ULL;
    const auto add = [&](const void* p, std::size_t n) {
        const auto* b = static_cast<const unsigned char*>(p);
        for (std::size_t i = 0; i < n; ++i) h = (h ^ b[i]) * 1099511628211ULL;
    };
    for (const Point& p : m.simplices.points) add(p.data(), sizeof(double) * 3);
    add(m.simplices.cells.data(), sizeof(m.simplices.cells[0]) * m.simplices.cells.size());
    add(m.node_region.data(), sizeof(std::int32_t) * m.node_region.size());
    return h;
}

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

// The slanted fin with the gate down to its base, so no plane crosses the sidewalls (26l SD6: meshes).
LayeredGeometry3D slanted_fin_uncrossed() {
    LayeredGeometry3D g = slanted_fin();
    g.solids[1] = box3(0, 1, 0.5, 1.5, 0.3, 1, 2);
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

TEST_CASE("mesh interface pieces: a wire through a gate makes holed pieces, which mesh and certify (IP1, IP2; IC2)",
          "[mesh][interface_pieces]") {
    // The square wire; 16 sides of 0.078 µm (s = 0.01 keeps every trace on the y faces over 6 s). The digests pin
    // the meshes (slice 1's own; deterministic): the holes' share of the band and of the clearances moves them.
    const std::vector<std::tuple<LayeredGeometry3D, double, std::uint64_t>> cases{
        {square_wire(), 0.02, 0xe556d52d1346d228ULL},
        {polygon_wire(16, 0.2), 0.01, 0x4943cd862869cdd9ULL},
    };
    for (const auto& [g, s, want] : cases) {
        CAPTURE(s);
        const auto m = mesh_of(g, s);
        if (!m) FAIL(m.error().message);
        check_certified(*m, s * um, 2.0 * um * um * um);
        const std::uint64_t got = digest(*m);
        UNSCOPED_INFO(std::format("nodes {}, digest 0x{:016x}", m->simplices.points.size(), got));
        CHECK(got == want);
    }
}

TEST_CASE("mesh interface pieces: layouts without holed pieces mesh to the same bits as before slice 1 (IC1)",
          "[mesh][interface_pieces]") {
    // Digests of points and cells, taken on the generator before Unit 64 slice 1 (be45c2d).
    const std::vector<std::tuple<std::string, LayeredGeometry3D, double, std::uint64_t>> cases{
        {"flush fin", flush_fin(), 0.02, 0x20c9b13ec0b7d2dcULL},
        {"flush fin gate", flush_fin_gate(), 0.02, 0x0dbb2f6452c8b994ULL},
        {"corner block", corner_block(), 0.02, 0x8adc4447c9a5b5aeULL},
        {"filleted flush fin", filleted_flush_fin(), 0.01, 0x9d00c7e7fa76f751ULL},
        {"slanted fin uncrossed", slanted_fin_uncrossed(), 0.02, 0xf0d248aeb847d2b4ULL},
    };
    for (const auto& [name, g, s, want] : cases) {
        CAPTURE(name);
        const auto m = mesh_of(g, s);
        if (!m) FAIL(m.error().message);
        const std::uint64_t got = digest(*m);
        UNSCOPED_INFO(std::format("{}: 0x{:016x}", name, got));
        CHECK(got == want);
    }
}

TEST_CASE("mesh interface pieces: holes too close, or pinched, are still refused (IC4)", "[mesh][interface_pieces]") {
    // Two square wires 0.05 µm (2.5 s) apart: their hole loops' sides are unrelated features under 4 s.
    LayeredGeometry3D two = domain();
    two.solids.push_back(box3(0.1, 0.9, 0.5, 1.5, 0.1, 0.9, 2));
    two.solids.push_back(box3(0.25, 0.45, 0, 2, 0.4, 0.6, 1));
    two.solids.push_back(box3(0.5, 0.7, 0, 2, 0.4, 0.6, 1));
    CHECK(outcome(two, 0.02) == "two unrelated interface features within 4 s");
    // A wire 0.05 µm from the gate's side: its hole within 4 s of the outer loop.
    LayeredGeometry3D near = domain();
    near.solids.push_back(box3(0.1, 0.9, 0.5, 1.5, 0.1, 0.9, 2));
    near.solids.push_back(box3(0.15, 0.35, 0, 2, 0.4, 0.6, 1));
    CHECK(outcome(near, 0.02) == "two unrelated interface features within 4 s");
    // Two wires touching along an edge: the gate's face has two holes meeting at a vertex.
    LayeredGeometry3D pinched = domain();
    pinched.solids.push_back(box3(0.1, 0.9, 0.5, 1.5, 0.1, 0.9, 2));
    pinched.solids.push_back(box3(0.3, 0.5, 0, 2, 0.3, 0.5, 1));
    pinched.solids.push_back(box3(0.5, 0.7, 0, 2, 0.5, 0.7, 1));
    CHECK(outcome(pinched, 0.02) == "an interface piece is not a simple polygon (a pinched vertex; unsupported)");
    // A diamond (square turned 45 degrees) wire whose corner is 30 degrees: refused by U12 as before.
    LayeredGeometry3D sharp = domain();
    sharp.solids.push_back(box3(0.1, 0.9, 0.5, 1.5, 0.1, 0.9, 2));
    sharp.solids.push_back(prism_y(0, 2, {{0.5, 0.3}, {0.5 + 0.05359, 0.5}, {0.5, 0.7}, {0.5 - 0.05359, 0.5}}, 1));
    CHECK(outcome(sharp, 0.01) == "an interface sector below 60 degrees (refused until justified, DECISIONS.md U12, FD7)");
}

namespace {

// IC5's seeded layouts, the first `count`: a gate slab (2) across a 1 x 1 x 0.6 µm oxide domain at z 0.2 to 0.4,
// pierced top to bottom by silicon (1): boxes 0.105 to 0.25 wide, or prisms along z of 4 to 6 sides of radius 0.1
// to 0.15, rotated; footprints kept 0.1 µm (> 4 s + 2 * 1.5 s) apart and from the sides. s = 0.015 µm.
void random_holes(int count) {
    const double s = 0.015;
    std::mt19937 rng(64);
    std::uniform_real_distribution<double> U(0.0, 1.0);
    int meshed = 0, refused = 0, holes_total = 0;
    std::map<std::string, int> why;
    const auto t0 = std::chrono::steady_clock::now();
    for (int trial = 0; trial < count; ++trial) {
        LayeredGeometry3D g{.lo = {0, 0, 0}, .hi = {1 * um, 1 * um, 0.6 * um}, .background = 0, .solids = {}, .patches = {}};
        g.solids.push_back(box3(0, 1, 0, 1, 0.2, 0.4, 2));
        const int want = 1 + static_cast<int>(U(rng) * 4.0);
        std::vector<std::array<double, 3>> placed;  // centre x, y and radius of each footprint
        for (int tries = 0; static_cast<int>(placed.size()) < want && tries < 400; ++tries) {
            const bool is_box = U(rng) < 0.5;
            const double r = is_box ? 0.075 + 0.05 * U(rng) : 0.1 + 0.05 * U(rng);
            const double reach = is_box ? r * std::sqrt(2.0) : r;
            const double cx = 0.1 + reach + (0.8 - 2 * reach) * U(rng);
            const double cy = 0.1 + reach + (0.8 - 2 * reach) * U(rng);
            if (std::ranges::any_of(placed, [&](const auto& p) { return std::hypot(p[0] - cx, p[1] - cy) < p[2] + reach + 0.1; })) continue;
            placed.push_back({cx, cy, reach});
            if (is_box) {
                const double rx = r, ry = r * (0.7 + 0.3 * U(rng));
                g.solids.push_back(box3(cx - rx, cx + rx, cy - ry, cy + ry, 0, 0.6, 1));
            } else {
                const int n = 4 + static_cast<int>(U(rng) * 3.0);
                const double phase = U(rng) * 2.0 * std::numbers::pi;
                LayeredSolid p;
                p.kind = LayeredSolid::Kind::prism;
                p.axis = 2;
                p.from = 0;
                p.to = 0.6 * um;
                for (int k = 0; k < n; ++k) {
                    const double t = phase + 2.0 * std::numbers::pi * k / n;
                    p.polygon.push_back({(cx + r * std::cos(t)) * um, (cy + r * std::sin(t)) * um, 0.0});
                }
                p.region = 1;
                g.solids.push_back(std::move(p));
            }
        }
        holes_total += static_cast<int>(placed.size());
        const auto c = layered_partition3d(g);
        REQUIRE(c.has_value());
        const auto m = generate_mesh(*c, options(s));
        if (m) {
            ++meshed;
            check_certified(*m, s * um, 0.6 * um * um * um);
        } else {
            ++refused;
            ++why[m.error().message];
        }
    }
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    WARN(std::format("IC5: {} meshed, {} refused, {} holes, {:.1f} s", meshed, refused, holes_total, secs));
    for (const auto& [msg, n] : why) WARN(std::format("IC5 refused {}: {}", n, msg));
    // Only the rules IC4 keeps may refuse.
    for (const auto& [msg, n] : why) {
        CAPTURE(msg);
        CHECK((msg.find("4 s") != std::string::npos || msg.find("60 degrees") != std::string::npos || msg.find("pinched") != std::string::npos ||
               msg.find("6 s") != std::string::npos));
    }
    CHECK(meshed > 0);
}

}  // namespace

TEST_CASE("mesh interface pieces: 20 seeded slabs pierced by 1 to 4 boxes and prisms mesh and certify (IC5, part)",
          "[mesh][interface_pieces]") {
    random_holes(20 / trial_stride);
}

// All 200 (Release about 8 minutes): run explicitly, as the record's IC5.
TEST_CASE("mesh interface pieces: 200 seeded slabs pierced by 1 to 4 boxes and prisms mesh and certify (IC5)",
          "[.interface_pieces_random][mesh][interface_pieces]") {
    random_holes(200);
}

TEST_CASE("mesh interface pieces: a holed piece's cost against its unholed parts (IC6, recorded)", "[mesh][interface_pieces]") {
    // The reference: the gate alone, plus what the wire alone adds to the empty domain.
    const auto count = [](const LayeredGeometry3D& g) {
        const auto m = mesh_of(g, 0.02);
        REQUIRE(m.has_value());
        return static_cast<double>(m->simplices.points.size());
    };
    LayeredGeometry3D gate = domain();
    gate.solids.push_back(box3(0.1, 0.9, 0.5, 1.5, 0.1, 0.9, 2));
    LayeredGeometry3D wire = domain();
    wire.solids.push_back(box3(0.4, 0.6, 0, 2, 0.4, 0.6, 1));
    LayeredGeometry3D empty = domain();
    const auto t0 = std::chrono::steady_clock::now();
    const double holed = count(square_wire());
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const double ref = count(gate) + count(wire) - count(empty);
    WARN(std::format("IC6: holed {} nodes, reference {} (ratio {:.3f}); generator {:.2f} s", holed, ref, holed / ref, secs));
    CHECK(holed > 0);
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
