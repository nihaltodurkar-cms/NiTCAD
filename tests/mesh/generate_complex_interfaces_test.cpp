// Gates for interfaces in the complex generator (Unit 53, slice 3; DECISIONS.md FD7, FD9; FC6, and
// FC2 for the FD7 rules): bends, junction lines and points, traces, non-convex pieces, P4.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "NiTCAD/mesh/complex.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/interfaces.hpp"
#include "NiTCAD/mesh/mesh.hpp"

using namespace NiTCAD::mesh;
using NiTCAD::base::ErrorCode;

namespace {

constexpr std::int32_t V = ExtrudedGeometry::void_region;
constexpr std::int32_t si = 0;
constexpr std::int32_t ox = 1;
constexpr std::int32_t fill = 2;
constexpr std::int32_t gate = 3;

Point P2(double u, double v) { return {u, v, 0.0}; }

ExtrusionCell quad(std::string name, Point a, Point b, Point c, Point d) { return {std::move(name), {a, b, c, d}}; }

// A conforming grid of rectangle cells over the breaks, along x, the region of each cell and
// section from `region(u_mid, v_mid, section)`.
ExtrudedGeometry grid(std::vector<double> ub, std::vector<double> vb, std::vector<double> stations,
                      const std::function<std::int32_t(double, double, std::size_t)>& region) {
    ExtrudedGeometry x;
    x.axis = 0;
    x.stations = std::move(stations);
    for (std::size_t i = 0; i + 1 < ub.size(); ++i) {
        for (std::size_t j = 0; j + 1 < vb.size(); ++j) {
            x.cells.push_back(quad("c" + std::to_string(i) + "_" + std::to_string(j), P2(ub[i], vb[j]), P2(ub[i + 1], vb[j]),
                                   P2(ub[i + 1], vb[j + 1]), P2(ub[i], vb[j + 1])));
            std::vector<std::int32_t> row;
            for (std::size_t k = 0; k + 1 < x.stations.size(); ++k) row.push_back(region(0.5 * (ub[i] + ub[i + 1]), 0.5 * (vb[j] + vb[j + 1]), k));
            x.regions.push_back(row);
        }
    }
    return x;
}

bool in(double a, double lo, double hi) { return a > lo && a < hi; }

// The liner-fin cross-section (y, z): a slab of height 0.8, a fin of width 1 and height 1.6 on it
// (both silicon), an oxide liner of thickness t on the fin's sides and top with its feet on the
// slab, a fill 0.8 wide beside and above. `top_inset` tapers the fin (0: rectangular).
struct Section {
    double t = 0.6;
    double top_inset = 0.0;
    [[nodiscard]] double y0() const { return 0.8; }             // the liner's outer left side
    [[nodiscard]] double y1() const { return 0.8 + t; }         // the fin's left foot
    [[nodiscard]] double y2() const { return 1.8 + t; }         // the fin's right foot
    [[nodiscard]] double y3() const { return 1.8 + 2.0 * t; }   // the liner's outer right side
    [[nodiscard]] double width() const { return 2.6 + 2.0 * t; }
    [[nodiscard]] double height() const { return 3.2 + t; }
};

// A fin with a liner wrapped on its sides and top, its feet on the slab, in a fill, along x
// (length 1): 90-degree bends along the fin's corners, three-region junction lines at the liner's
// feet, traces on the bottom, side, top and end faces, crossings of the domain edges. With
// `gated`, three sections (0.8 each) and the fill in the middle one a gate: an inverted-U
// (non-convex) piece at each gate edge, with junction points where the liner's lines meet it.
ExtrudedGeometry liner_fin(double t, bool gated) {
    const Section c{t, 0.0};
    const std::vector<double> ub{0.0, c.y0(), c.y1(), c.y2(), c.y3(), c.width()};
    const std::vector<double> vb{0.0, 0.8, 2.4, 2.4 + t, c.height()};
    const std::vector<double> stations = gated ? std::vector<double>{0.0, 0.8, 1.6, 2.4} : std::vector<double>{0.0, 1.0};
    return grid(ub, vb, stations, [&](double u, double v, std::size_t k) {
        if (v < 0.8) return si;
        if (in(u, c.y1(), c.y2()) && v < 2.4) return si;
        if (in(u, c.y0(), c.y3()) && v < 2.4 + t) return ox;
        return gated && k == 1 ? gate : fill;
    });
}

// The same with the fin tapered: its sides at 85 degrees to the slab (95 at its top corners).
ExtrudedGeometry tapered_liner_fin(double t) {
    const Section c{t, 1.6 / std::tan(85.0 * 3.14159265358979323846 / 180.0)};
    const double a = c.y1() + c.top_inset;
    const double b = c.y2() - c.top_inset;
    const double y0 = c.y0();
    const double y1 = c.y1();
    const double y2 = c.y2();
    const double y3 = c.y3();
    const double w = c.width();
    const double z1 = 2.4;
    const double z2 = 2.4 + t;
    const double z3 = c.height();
    ExtrudedGeometry x;
    x.axis = 0;
    x.stations = {0.0, 1.0};
    x.cells = {
        quad("slab0", P2(0, 0), P2(y0, 0), P2(y0, 0.8), P2(0, 0.8)),
        quad("slab1", P2(y0, 0), P2(y1, 0), P2(y1, 0.8), P2(y0, 0.8)),
        quad("slab2", P2(y1, 0), P2(y2, 0), P2(y2, 0.8), P2(y1, 0.8)),
        quad("slab3", P2(y2, 0), P2(y3, 0), P2(y3, 0.8), P2(y2, 0.8)),
        quad("slab4", P2(y3, 0), P2(w, 0), P2(w, 0.8), P2(y3, 0.8)),
        quad("fin", P2(y1, 0.8), P2(y2, 0.8), P2(b, z1), P2(a, z1)),
        quad("liner_l", P2(y0, 0.8), P2(y1, 0.8), P2(a, z1), P2(y0, z1)),
        quad("liner_r", P2(y2, 0.8), P2(y3, 0.8), P2(y3, z1), P2(b, z1)),
        quad("liner_tl", P2(y0, z1), P2(a, z1), P2(a, z2), P2(y0, z2)),
        quad("liner_t", P2(a, z1), P2(b, z1), P2(b, z2), P2(a, z2)),
        quad("liner_tr", P2(b, z1), P2(y3, z1), P2(y3, z2), P2(b, z2)),
        quad("fill_l", P2(0, 0.8), P2(y0, 0.8), P2(y0, z1), P2(0, z1)),
        quad("fill_r", P2(y3, 0.8), P2(w, 0.8), P2(w, z1), P2(y3, z1)),
        quad("fill_lm", P2(0, z1), P2(y0, z1), P2(y0, z2), P2(0, z2)),
        quad("fill_rm", P2(y3, z1), P2(w, z1), P2(w, z2), P2(y3, z2)),
        quad("top0", P2(0, z2), P2(y0, z2), P2(y0, z3), P2(0, z3)),
        quad("top1", P2(y0, z2), P2(a, z2), P2(a, z3), P2(y0, z3)),
        quad("top2", P2(a, z2), P2(b, z2), P2(b, z3), P2(a, z3)),
        quad("top3", P2(b, z2), P2(y3, z2), P2(y3, z3), P2(b, z3)),
        quad("top4", P2(y3, z2), P2(w, z2), P2(w, z3), P2(y3, z3)),
    };
    x.regions = {{si}, {si}, {si}, {si}, {si}, {si}, {ox}, {ox}, {ox}, {ox}, {ox},
                 {fill}, {fill}, {fill}, {fill}, {fill}, {fill}, {fill}, {fill}, {fill}};
    return x;
}

// Two regions split by an oblique interface meeting the bottom face at `deg` degrees (in the
// right-hand region) and the top at 180 - deg.
ExtrudedGeometry oblique_split(double deg) {
    const double h = 2.0;
    const double d = h / std::tan(deg * 3.14159265358979323846 / 180.0);
    ExtrudedGeometry x;
    x.axis = 0;
    x.stations = {0.0, 1.5};
    x.cells = {quad("a", P2(0, 0), P2(1.5, 0), P2(1.5 + d, h), P2(0, h)), quad("b", P2(1.5, 0), P2(4, 0), P2(4, h), P2(1.5 + d, h))};
    x.regions = {{0}, {1}};
    return x;
}

MeshOptions3D options(double h, double s) {
    MeshOptions3D o;
    o.max_size = h;
    o.interface_resolution = s;
    return o;
}

GeneratedMesh mesh_of(const ExtrudedGeometry& x, const MeshOptions3D& o) {
    const auto g = extrude(x);
    if (!g) FAIL(g.error().message);
    auto m = generate_mesh(*g, o);
    if (!m) FAIL(m.error().message);
    return std::move(*m);
}

void refused(const ExtrudedGeometry& x, const MeshOptions3D& o, const std::string& text) {
    const auto g = extrude(x);
    REQUIRE(g);
    const auto m = generate_mesh(*g, o);
    REQUIRE_FALSE(m);
    CHECK(m.error().code == ErrorCode::invalid_input);
    INFO(m.error().message);
    CHECK(m.error().message.find(text) != std::string::npos);
    CHECK((m.error().context.has_value() && m.error().context->index.has_value()));
}

// The FC6 checks on one mesh; the report for the record.
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
    UNSCOPED_INFO("nodes " << m.report.nodes << ", interface nodes " << m.report.interface_nodes << ", cross faces " << r.cross_faces
                           << ", corner faces " << r.corner_faces << ", flat " << r.max_flat_deviation << ", corner/s "
                           << r.max_corner_deviation / s << ", exempt " << m.report.exempt_tetrahedra << ", slivers " << m.report.slivers
                           << ", band facet points " << m.report.band_facet_points);
}

}  // namespace

TEST_CASE("complex interfaces: a liner wrapped on a fin, its feet on the slab (FC6)", "[mesh][complex_interfaces]") {
    const double s = 0.1;
    const double t = 0.7;
    const auto x = liner_fin(t, false);
    const auto m = mesh_of(x, options(0.5, s));
    check_certified(m, s, (2.6 + 2.0 * t) * (3.2 + t) * 1.0);
    // Measured 0.286 s (Release, 2026-10-30): at the end face, where the liner's outer corner line
    // meets it (a facet point in the band next to the trace junction).
    CHECK(m.certificate->report.max_corner_deviation <= 0.32 * s);
    CHECK(m.report.nodes <= 2000);
    const auto again = mesh_of(x, options(0.5, s));
    CHECK(std::ranges::equal(again.mesh.points(), m.mesh.points()));
    CHECK(again.node_region == m.node_region);
    WARN("liner fin");
}

TEST_CASE("complex interfaces: the tapered fin, 85 and 95 degree bends (FC6)", "[mesh][complex_interfaces]") {
    const double s = 0.1;
    const double t = 0.7;
    const auto m = mesh_of(tapered_liner_fin(t), options(0.5, s));
    check_certified(m, s, (2.6 + 2.0 * t) * (3.2 + t) * 1.0);
    CHECK(m.certificate->report.max_corner_deviation <= 0.16 * s);  // measured 0.141 s
    CHECK(m.report.nodes <= 2000);
    WARN("tapered liner fin");
}

TEST_CASE("complex interfaces: a gate section, junction points and non-convex pieces (FC6, FD9)", "[mesh][complex_interfaces]") {
    const double s = 0.1;
    const double t = 0.7;
    const auto m = mesh_of(liner_fin(t, true), options(0.5, s));
    check_certified(m, s, (2.6 + 2.0 * t) * (3.2 + t) * 2.4);
    // Measured 0.411 s: at the junction point where the liner's outer corner line meets a gate edge.
    CHECK(m.certificate->report.max_corner_deviation <= 0.45 * s);
    WARN("gated liner fin");
}

TEST_CASE("complex interfaces: a 60.5 degree sector at a trace (FC6)", "[mesh][complex_interfaces]") {
    const double s = 0.1;
    const auto m = mesh_of(oblique_split(60.5), options(0.5, s));
    check_certified(m, s, 4.0 * 2.0 * 1.5);
    // Measured 0.446 s: beside the trace, in the 60.5-degree sector (U15 measured 0.36 s at 60.3 in 2D).
    CHECK(m.certificate->report.max_corner_deviation <= 0.49 * s);
    WARN("60.5 degrees");
}

// The shaped FinFET deck's gate edge in small (Unit 53, owner ruling on slice 5): a tapered fin wrapped
// by a liner and a fill, two sections; in the second the liner and the fill are other regions,
// so at the station between them the liner's U-shaped ring is a non-convex piece beside the fill's coplanar
// piece, both against the first section's region, with oblique (about 104 and 256 degree) corners.
// Face vertices there lie a few hundredths of s outside the ring's outer edge; the decomposition
// bound alone fails P4b here (mutant checked), as in decks/finfet_3d_shaped.deck (1.1 s in its
// miniature); the grown polygons bound it (DECISIONS.md, Unit 53, the coarser offsets).
TEST_CASE("complex interfaces: a tapered liner ring at a gate edge (P4b, the shaped deck's corner)", "[mesh][complex_interfaces]") {
    const double s = 0.16;
    ExtrudedGeometry x;
    x.axis = 0;
    x.stations = {0.0, 1.0, 2.0};
    x.cells = {
        {"base", {P2(3.3, 0), P2(4, 0), P2(4, 1), P2(4, 2), P2(4, 4), P2(4, 5), P2(4, 6), P2(3.3, 6)}},
        {"fin", {P2(4, 2), P2(6, 2.5), P2(6, 3.5), P2(4, 4)}},
        {"liner", {P2(4, 1), P2(7, 1.5), P2(7, 4.5), P2(4, 5), P2(4, 4), P2(6, 3.5), P2(6, 2.5), P2(4, 2)}},
        {"fill", {P2(4, 0), P2(7.7, 0), P2(7.7, 6), P2(4, 6), P2(4, 5), P2(7, 4.5), P2(7, 1.5), P2(4, 1)}},
    };
    // Regions: 0 the base, 1 the first section's liner and fill, 2 the gate fill, 3 the gate liner, 4 the fin.
    x.regions = {{0, 0}, {4, 4}, {1, 3}, {1, 2}};
    const auto m = mesh_of(x, options(2.0, s));
    check_certified(m, s, 2.0 * (0.7 * 6.0 + 3.7 * 6.0));
    WARN("tapered liner ring");
}

TEST_CASE("complex interfaces: the FD7 rules are refused (FC2, FC6)", "[mesh][complex_interfaces]") {
    const double s = 0.1;
    const MeshOptions3D o = options(0.5, s);
    SECTION("no resolution") { refused(liner_fin(0.7, false), options(0.5, 0.0), "interface_resolution"); }
    SECTION("a 59.5 degree sector") { refused(oblique_split(59.5), o, "sector below 60 degrees"); }
    SECTION("an interface ending on a domain edge") {
        ExtrudedGeometry x;
        x.axis = 0;
        x.stations = {0.0, 1.0};
        x.cells = {quad("a", P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1)), quad("b", P2(1, 0), P2(2, 0), P2(2, 1), P2(1, 1)),
                   quad("c", P2(1, 1), P2(2, 1), P2(2, 2), P2(1, 2))};
        x.regions = {{0}, {1}, {1}};
        refused(x, o, "ends on a domain edge");
    }
    SECTION("two interfaces ending on the boundary along one line") {
        ExtrudedGeometry x;
        x.axis = 0;
        x.stations = {0.0, 1.0};
        // The flat bottom split three ways at (1, 0): 45, 90 and 45 degrees.
        x.cells = {{"a", {P2(0, 0), P2(1, 0), P2(0.4, 0.6), P2(0.4, 1), P2(0, 1)}},
                   {"m", {P2(1, 0), P2(1.6, 0.6), P2(1.6, 1), P2(0.4, 1), P2(0.4, 0.6)}},
                   {"b", {P2(1, 0), P2(2, 0), P2(2, 1), P2(1.6, 1), P2(1.6, 0.6)}}};
        x.regions = {{0}, {1}, {2}};
        refused(x, o, "two or more interfaces end on the boundary");
    }
    const auto slabs = [](std::vector<double> zb, double length) {
        return grid({0.0, 2.0}, std::move(zb), {0.0, length}, [](double, double v, std::size_t) { return v < 1.0 ? 0 : 1; });
    };
    // A layer ending on the boundary: its two crossings of each side edge are within 4 s, refused by the crossing
    // rule before the layer rule is reached (FC13 found this; the closed shell below reaches the layer rule).
    SECTION("a layer thinner than 4 s, ending on the boundary") {
        const auto x = grid({0.0, 2.0}, {0.0, 1.0, 1.3, 2.3}, {0.0, 1.0}, [](double, double v, std::size_t) { return v < 1.0 ? 0 : v < 1.3 ? 1 : 2; });
        refused(x, o, "a domain-edge vertex within 4 s of an interface crossing");
    }
    SECTION("a closed shell 3.5 s thick") {
        // A layer that ends anywhere meets another rule first (a domain-edge crossing, a short trace or feature
        // line). A closed shell around a core, inside the domain, does not: its inner and outer bend lines are
        // t sqrt 2 = 4.95 s apart, so only its parallel faces, 3.5 s apart, are refused (FC13: the rule loosened
        // to 3 s accepts it).
        const std::vector<double> b{0.0, 1.0, 1.35, 3.35, 3.7, 4.7};
        const auto x = grid(b, b, b, [](double u, double v, std::size_t k) {
            const bool core = in(u, 1.35, 3.35) && in(v, 1.35, 3.35) && k == 2;
            const bool shell = in(u, 1.0, 3.7) && in(v, 1.0, 3.7) && k >= 1 && k <= 3;
            return core ? 2 : shell ? 1 : 0;
        });
        refused(x, o, "two interfaces within 4 s");
    }
    SECTION("an interface within 4 s of a parallel boundary face") {
        // An enclosed box whose lower face is 3.5 s above the boundary face v = 0, which it does not meet. (A
        // slab reaching the side faces is refused by the crossing rule first.)
        const auto x = grid({0.0, 1.0, 2.0, 3.0}, {0.0, 0.35, 1.35, 2.35}, {0.0, 1.0, 2.0, 3.0},
                            [](double u, double v, std::size_t k) { return in(u, 1.0, 2.0) && in(v, 0.35, 1.35) && k == 1 ? 1 : 0; });
        refused(x, o, "an interface within 4 s of a boundary face it does not meet");
    }
    // The full message: the feature-line rule refuses a short trace too when the trace rule is loosened (FC13).
    SECTION("a trace shorter than 6 s") { refused(slabs({0.0, 1.0, 2.0}, 0.5), o, "an interface trace on the boundary is shorter than 6 s"); }
    SECTION("a feature line shorter than 6 s") {
        const auto x = grid({0.0, 1.0, 1.5, 2.5}, {0.0, 1.0, 1.5, 2.5}, {0.0, 1.0, 1.5, 2.5},
                            [](double u, double v, std::size_t k) { return in(u, 1.0, 1.5) && in(v, 1.0, 1.5) && k == 1 ? 1 : 0; });
        refused(x, o, "feature line is shorter than 6 s");
    }
    SECTION("an interface piece with a hole") {
        const auto x = grid({0.0, 1.0, 2.0, 3.0}, {0.0, 1.0, 2.0, 3.0}, {0.0, 1.0, 2.0},
                            [](double u, double v, std::size_t k) { return k == 1 && !(in(u, 1.0, 2.0) && in(v, 1.0, 2.0)) ? 1 : 0; });
        refused(x, o, "has holes");
    }
    SECTION("two unrelated features within 4 s") {
        // Two blocks inside region 0, 0.3 apart: their side bends are unrelated lines.
        const auto x = grid({0.0, 1.0, 2.0, 2.3, 3.3, 4.3}, {0.0, 1.0, 2.0, 3.0}, {0.0, 1.0},
                            [](double u, double v, std::size_t) { return in(v, 1.0, 2.0) && (in(u, 1.0, 2.0) || in(u, 2.3, 3.3)) ? 1 : 0; });
        refused(x, o, "two unrelated interface features within 4 s");
    }
    SECTION("a feature within 4 s of a domain edge") {
        // A block 0.3 from the side face: its traces on the end faces run 0.3 from their edges.
        const auto x = grid({0.0, 0.3, 1.3, 3.0}, {0.0, 1.0, 2.0, 3.0}, {0.0, 1.0},
                            [](double u, double v, std::size_t) { return in(v, 1.0, 2.0) && in(u, 0.3, 1.3) ? 1 : 0; });
        refused(x, o, "feature within 4 s of a domain edge");
    }
    SECTION("a domain-edge vertex within 4 s of a crossing") { refused(slabs({0.0, 1.0, 1.2, 2.0}, 1.0), o, "within 4 s of an interface crossing"); }
    // A patch edge on a trace, refused here in slice 3, is meshed since slice 4
    // (generate_complex_patches_test.cpp).
}

// FD11 (FC13): a region's own h is honoured: the slab of region 1 at half the global h holds more
// nodes than without it, and the other slab's count is about unchanged.
TEST_CASE("complex interfaces: a per-region size refines that region only (FD11)", "[mesh][complex_interfaces]") {
    const auto x = grid({0.0, 2.0}, {0.0, 1.0, 2.0}, {0.0, 1.0}, [](double, double v, std::size_t) { return v < 1.0 ? 0 : 1; });
    MeshOptions3D o = options(0.5, 0.1);
    const auto count = [](const GeneratedMesh& m, std::int32_t r) {
        return static_cast<std::size_t>(std::ranges::count(m.node_region, r));
    };
    const GeneratedMesh plain = mesh_of(x, o);
    o.region_max_size = {0.0, 0.25};
    const GeneratedMesh fine = mesh_of(x, o);
    REQUIRE(fine.certificate.has_value());
    UNSCOPED_INFO("region 1 nodes " << count(plain, 1) << " -> " << count(fine, 1) << ", region 0 " << count(plain, 0)
                                    << " -> " << count(fine, 0));
    CHECK(count(fine, 1) > 2 * count(plain, 1));
    CHECK(count(fine, 0) < 2 * count(plain, 0));
}
