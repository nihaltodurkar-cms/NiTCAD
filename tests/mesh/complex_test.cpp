// The polyhedral complex, its validation and exact point location, and the sectioned extrusion
// (Unit 53, slice 1; DECISIONS.md FD1, FD2, FD6; FC2 for these rules).
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "NiTCAD/mesh/complex.hpp"
#include "complex_index.hpp"

using namespace NiTCAD::mesh;
using NiTCAD::base::ErrorCode;
using Kind = detail::ComplexIndex::Location::Kind;

namespace {

// An axis-aligned box [lo, hi] as eight vertices from `base` and six facets, appended to g as one volume.
void add_box(PolyhedralGeometry& g, Point lo, Point hi, std::int32_t region) {
    const std::size_t b = g.vertices.size();
    for (int k = 0; k < 8; ++k) g.vertices.push_back({(k & 1) ? hi[0] : lo[0], (k & 2) ? hi[1] : lo[1], (k & 4) ? hi[2] : lo[2]});
    const std::size_t f = g.facets.size();
    const std::vector<std::vector<std::size_t>> faces{{0, 4, 6, 2}, {1, 3, 7, 5}, {0, 1, 5, 4}, {2, 6, 7, 3}, {0, 2, 3, 1}, {4, 5, 7, 6}};
    ComplexVolume v{{}, region};
    for (const auto& q : faces) {
        std::vector<std::size_t> idx;
        for (const std::size_t i : q) idx.push_back(b + i);
        g.facets.push_back({idx, "outer"});
        v.facets.push_back(g.facets.size() - 1);
    }
    (void)f;
    g.volumes.push_back(v);
}

PolyhedralGeometry cube() {
    PolyhedralGeometry g;
    add_box(g, {0, 0, 0}, {1, 1, 1}, 0);
    return g;
}

// Two unit cubes side by side along x sharing the facet x = 1 (the second cube's x_min facet removed).
PolyhedralGeometry two_cubes(std::int32_t r0, std::int32_t r1) {
    PolyhedralGeometry g;
    add_box(g, {0, 0, 0}, {1, 1, 1}, r0);
    // The second cube reuses vertices 1, 3, 5, 7 of the first (x = 1).
    const std::size_t b = g.vertices.size();
    for (const Point& p : std::vector<Point>{{2, 0, 0}, {2, 1, 0}, {2, 0, 1}, {2, 1, 1}}) g.vertices.push_back(p);
    const std::size_t x0 = 1, x2 = 3, x4 = 5, x6 = 7;  // the shared corners at x = 1: (y,z) = 00, 10, 01, 11
    const std::size_t y0 = b, y2 = b + 1, y4 = b + 2, y6 = b + 3;
    ComplexVolume v{{1}, r1};  // facet 1 is the first cube's x_max facet
    for (const auto& q : std::vector<std::vector<std::size_t>>{{y0, y2, y6, y4}, {x0, y0, y4, x4}, {x2, x6, y6, y2}, {x0, x2, y2, y0}, {x4, y4, y6, x6}}) {
        g.facets.push_back({q, "outer"});
        v.facets.push_back(g.facets.size() - 1);
    }
    g.volumes.push_back(v);
    return g;
}

ErrorCode code(const auto& r) { return r.error().code; }

// A FinFET-like cross-section in (y, z) extruded along x: a buried oxide slab, a tapered fin on it
// with a liner on its three sides, a high-k shell around the liner; the gate metal (void) outside
// the shell under the gate, spacers in its place beside it, the raised source/drain cells above.
ExtrudedGeometry finfet_like() {
    ExtrudedGeometry x;
    x.axis = 0;
    x.stations = {0.0, 3.0, 4.0, 7.0, 8.0, 11.0};
    // y up (u), z across (v). Buried oxide y in [0, 2], z in [0, 10].
    // Fin: base z in [4, 6] at y = 2, top z in [4.2, 5.8] at y = 5 (tapered).
    // Liner 0.5 thick around the fin (as a polygon with straight offsets, kept simple): left, top, right cells.
    // High-k: one shell polygon around the liner; outside it, metal or spacer.
    const auto P = [](double y, double z) { return Point{y, z, 0.0}; };
    x.cells = {
        {"box_l", {P(0, 0), P(0, 3), P(2, 3), P(2, 0)}},
        {"box_c", {P(0, 3), P(0, 7), P(2, 7), P(2, 6.5), P(2, 6), P(2, 4), P(2, 3.5), P(2, 3)}},
        {"box_r", {P(0, 7), P(0, 10), P(2, 10), P(2, 7)}},
        {"fin", {P(2, 4), P(2, 6), P(5, 5.8), P(5, 4.2)}},
        {"liner", {P(2, 3.5), P(2, 4), P(5, 4.2), P(5, 5.8), P(2, 6), P(2, 6.5), P(5.5, 6.3), P(5.5, 3.7)}},
        {"shell", {P(2, 3), P(2, 3.5), P(5.5, 3.7), P(5.5, 6.3), P(2, 6.5), P(2, 7), P(6, 7), P(6, 3)}},
        {"side_l", {P(2, 0), P(2, 3), P(6, 3), P(6, 0)}},
        {"side_r", {P(2, 7), P(2, 10), P(6, 10), P(6, 7)}},
    };
    constexpr std::int32_t V = ExtrudedGeometry::void_region;
    // Regions: 0 oxide, 1 silicon, 2 liner oxide (a separate region index), 3 high-k, 4 spacer.
    // Sections: source (0-3), spacer (3-4), gate (4-7), spacer (7-8), drain (8-11).
    x.regions = {
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0},
        {1, 1, 1, 1, 1},
        {1, 4, 2, 4, 1},  // the liner: raised source/drain silicon, spacer, liner under the gate
        {1, 4, 3, 4, 1},
        {V, 4, V, 4, V},
        {V, 4, V, 4, V},
    };
    x.patches = {{"source", {0, 0, 0}, {0, 6, 10}}, {"drain", {11, 0, 0}, {11, 6, 10}}, {"gate", {4, 2, 3}, {7, 6, 7}},
                 {"substrate", {0, 0, 0}, {11, 0, 10}}};
    return x;
}

}  // namespace

TEST_CASE("complex: a cube, two cubes with an interface or an internal facet", "[complex]") {
    const auto c = check_complex(cube());
    REQUIRE(c);
    CHECK(c->boundary_facets == 6);
    CHECK(c->interface_facets == 0);
    REQUIRE(c->volume.size() == 1);
    CHECK(std::abs(c->volume[0] - 1.0) < 1e-15);
    CHECK(c->patches == std::vector<std::string>{"outer"});

    const auto i = check_complex(two_cubes(0, 1));
    REQUIRE(i);
    CHECK(i->boundary_facets == 10);
    CHECK(i->interface_facets == 1);
    CHECK(i->internal_facets == 0);
    const auto s = check_complex(two_cubes(2, 2));
    REQUIRE(s);
    CHECK(s->interface_facets == 0);
    CHECK(s->internal_facets == 1);
}

TEST_CASE("complex: exact location inside, outside, on facets, at tie positions", "[complex]") {
    const auto idx = detail::ComplexIndex::create(two_cubes(0, 1));
    REQUIRE(idx);
    // Points on a lattice of quarter steps: every tie with the facets' edges and planes occurs.
    for (int i = -2; i <= 10; ++i) {
        for (int j = -2; j <= 6; ++j) {
            for (int k = -2; k <= 6; ++k) {
                const Point p{0.25 * i, 0.25 * j, 0.25 * k};
                const auto l = idx->locate(p);
                const bool in_x = p[0] > 0 && p[0] < 2;
                const bool in_y = p[1] > 0 && p[1] < 1;
                const bool in_z = p[2] > 0 && p[2] < 1;
                const bool closed = p[0] >= 0 && p[0] <= 2 && p[1] >= 0 && p[1] <= 1 && p[2] >= 0 && p[2] <= 1;
                if (in_x && in_y && in_z && p[0] != 1.0) {
                    REQUIRE(l.kind == Kind::inside);
                    CHECK(l.index == (p[0] < 1 ? 0 : 1));
                } else if (closed) {
                    CHECK(l.kind == Kind::on_facet);
                } else {
                    CHECK(l.kind == Kind::outside);
                }
            }
        }
    }
}

TEST_CASE("complex: refusals (FD1)", "[complex]") {
    const auto refused = [](const PolyhedralGeometry& g, std::optional<std::size_t> index) {
        const auto r = check_complex(g);
        REQUIRE_FALSE(r);
        CHECK(r.error().code == ErrorCode::invalid_input);
        REQUIRE(r.error().context);
        if (index) CHECK(r.error().context->index == index);
        return r.error().message;
    };
    SECTION("empty") { CHECK_FALSE(check_complex(PolyhedralGeometry{})); }
    SECTION("vertex not finite") {
        auto g = cube();
        g.vertices[3][1] = std::numeric_limits<double>::infinity();
        refused(g, 3);
    }
    SECTION("vertices coincide") {
        auto g = cube();
        g.vertices.push_back(g.vertices[2]);
        refused(g, 8);
    }
    SECTION("facet too small, out of range, repeated, collinear") {
        auto g = cube();
        g.facets[2].vertices = {0, 1};
        refused(g, 2);
        g = cube();
        g.facets[2].vertices[1] = 99;
        refused(g, 2);
        g = cube();
        g.facets[2].vertices = {0, 1, 1, 4};
        refused(g, 2);
        g = cube();
        g.vertices.push_back({0.5, 0, 0});
        g.facets[2].vertices = {0, 8, 1};
        CHECK(refused(g, 2).find("collinear") != std::string::npos);
    }
    SECTION("facet not planar") {
        auto g = cube();
        g.vertices[7][2] = 1.0 + 1e-15;
        const auto r = check_complex(g);
        REQUIRE_FALSE(r);
        CHECK(r.error().message.find("planar") != std::string::npos);
    }
    SECTION("facet not simple (a bow tie)") {
        auto g = cube();
        g.facets[4].vertices = {0, 2, 1, 3};
        CHECK(refused(g, 4).find("simple") != std::string::npos);
    }
    SECTION("negative region, empty volume, facet listed twice") {
        auto g = cube();
        g.volumes[0].region = -1;
        refused(g, 0);
        g = cube();
        g.volumes[0].facets.push_back(0);
        refused(g, 0);
        g = cube();
        g.volumes.push_back({{}, 0});
        refused(g, 1);
    }
    SECTION("volume not closed") {
        auto g = cube();
        g.volumes[0].facets.pop_back();
        CHECK(refused(g, 0).find("closed") != std::string::npos);
    }
    SECTION("facet on no volume, on three volumes") {
        auto g = cube();
        g.vertices.push_back({5, 5, 5});
        g.facets.push_back({{0, 1, 8}, "x"});
        CHECK(refused(g, 6).find("no volume") != std::string::npos);
        auto h = cube();
        h.volumes.push_back(h.volumes[0]);
        h.volumes.push_back(h.volumes[0]);
        CHECK(refused(h, 0).find("more than two") != std::string::npos);
    }
    SECTION("boundary facet without a patch") {
        auto g = cube();
        g.facets[5].patch.clear();
        refused(g, 5);
    }
    SECTION("two volumes on the same side of their facets") {
        auto g = cube();
        g.volumes.push_back(g.volumes[0]);
        CHECK(refused(g, 0).find("same side") != std::string::npos);
    }
    SECTION("facets intersect: overlapping cubes") {
        PolyhedralGeometry g;
        add_box(g, {0, 0, 0}, {1, 1, 1}, 0);
        add_box(g, {0.5, 0.5, 0.5}, {1.5, 1.5, 1.5}, 0);
        CHECK(refused(g, std::nullopt).find("intersect") != std::string::npos);
    }
    SECTION("facets touch improperly: a cube on another's face, off its corners (T-junctions)") {
        PolyhedralGeometry g;
        add_box(g, {0, 0, 0}, {2, 2, 1}, 0);
        add_box(g, {0.5, 0.5, 1}, {1.5, 1.5, 2}, 1);
        CHECK(refused(g, std::nullopt).find("intersect") != std::string::npos);
    }
    SECTION("nested volumes") {
        PolyhedralGeometry g;
        add_box(g, {0, 0, 0}, {3, 3, 3}, 0);
        add_box(g, {1, 1, 1}, {2, 2, 2}, 1);
        CHECK(refused(g, 1).find("overlap") != std::string::npos);
    }
    SECTION("volumes not connected") {
        PolyhedralGeometry g;
        add_box(g, {0, 0, 0}, {1, 1, 1}, 0);
        add_box(g, {2, 0, 0}, {3, 1, 1}, 0);
        CHECK(refused(g, 1).find("connected") != std::string::npos);
    }
    SECTION("two cubes touching at an edge only (not connected through a facet)") {
        PolyhedralGeometry g;
        add_box(g, {0, 0, 0}, {1, 1, 1}, 0);
        add_box(g, {1, 1, 0}, {2, 2, 1}, 0);
        CHECK_FALSE(check_complex(g));
    }
}

TEST_CASE("extrusion: an L-prism and the FinFET-like section", "[complex]") {
    ExtrudedGeometry l;
    l.axis = 2;
    l.stations = {0.0, 1.0, 2.0};
    const auto P = [](double u, double v) { return Point{u, v, 0.0}; };
    l.cells = {{"a", {P(0, 0), P(1, 0), P(1, 1), P(0, 1)}}, {"b", {P(1, 0), P(2, 0), P(2, 1), P(1, 1)}},
               {"c", {P(0, 1), P(1, 1), P(1, 2), P(0, 2)}}};
    l.regions = {{0, 0}, {0, ExtrudedGeometry::void_region}, {0, 0}};
    const auto g = extrude(l);
    REQUIRE(g);
    const auto r = check_complex(*g);
    REQUIRE(r);
    CHECK(r->volume.size() == 5);
    double total = 0.0;
    for (const double v : r->volume) total += v;
    CHECK(std::abs(total - 5.0) < 1e-14);
    CHECK(r->interface_facets == 0);
    // Internal facets: a|b, a|c in both sections, plus the station facets of a and c at z = 1.
    CHECK(r->internal_facets == 5);
    for (const Point& p : g->vertices) CHECK(p[2] >= 0.0);  // axis z: (u, v) are (x, y)

    const auto f = extrude(finfet_like());
    REQUIRE(f);
    const auto idx = detail::ComplexIndex::create(*f);
    REQUIRE(idx);
    CHECK(idx->report().interface_facets > 0);
    auto names = idx->report().patches;
    std::ranges::sort(names);
    CHECK(names == std::vector<std::string>{"boundary", "drain", "gate", "source", "substrate"});
    const auto region_at = [&](const Point& p) {
        const auto loc = idx->locate(p);
        REQUIRE(loc.kind == Kind::inside);
        return f->volumes[static_cast<std::size_t>(loc.index)].region;
    };
    // Axis x: points are (x, y, z) with (u, v) = (y, z).
    CHECK(region_at({5.5, 3.0, 5.0}) == 1);    // the fin under the gate
    CHECK(region_at({5.5, 3.0, 4.05}) == 2);   // the liner beside the tapered fin, under the gate
    CHECK(region_at({1.0, 3.0, 4.05}) == 1);   // the same cell in the source: raised silicon
    CHECK(region_at({3.5, 5.9, 8.0}) == 4);    // spacer
    CHECK(region_at({5.5, 1.0, 1.0}) == 0);    // buried oxide
    CHECK(idx->locate({5.5, 4.0, 1.0}).kind == Kind::outside);  // the gate metal (void)
    CHECK(idx->locate({1.0, 4.0, 1.0}).kind == Kind::outside);  // above the oxide beside the source
    CHECK(idx->locate({5.5, 2.0, 5.0}).kind == Kind::on_facet);  // the fin's foot
    // The gate patch: the shell's outer facets under the gate and the oxide top beside it.
    std::size_t gate = 0;
    for (const auto& F : f->facets) gate += F.patch == "gate" ? 1 : 0;
    CHECK(gate > 0);
    // Deterministic.
    const auto again = extrude(finfet_like());
    REQUIRE(again);
    CHECK(again->vertices == f->vertices);
}

TEST_CASE("extrusion: refusals (FD2)", "[complex]") {
    const auto refused = [](const ExtrudedGeometry& x, const char* what) {
        const auto r = extrude(x);
        REQUIRE_FALSE(r);
        CHECK(r.error().code == ErrorCode::invalid_input);
        INFO(r.error().message);
        CHECK(r.error().message.find(what) != std::string::npos);
    };
    const auto P = [](double u, double v) { return Point{u, v, 0.0}; };
    ExtrudedGeometry ok;
    ok.stations = {0.0, 1.0};
    ok.cells = {{"a", {P(0, 0), P(1, 0), P(1, 1), P(0, 1)}}, {"b", {P(1, 0), P(2, 0), P(2, 1), P(1, 1)}}};
    ok.regions = {{0}, {1}};
    REQUIRE(extrude(ok));
    auto x = ok;
    x.axis = 3;
    refused(x, "axis");
    x = ok;
    x.stations = {0.0};
    refused(x, "two stations");
    x = ok;
    x.stations = {0.0, 0.0};
    refused(x, "increase");
    x = ok;
    x.stations = {0.0, std::numeric_limits<double>::quiet_NaN()};
    refused(x, "finite");
    x = ok;
    x.regions.pop_back();
    refused(x, "one row per cell");
    x = ok;
    x.cells[1].name = "a";
    refused(x, "unique");
    x = ok;
    x.cells[1].polygon.pop_back();
    x.cells[1].polygon.pop_back();
    refused(x, "fewer than three");
    x = ok;
    x.cells[1].polygon[2] = P(2, 0);
    refused(x, "repeats");
    x = ok;
    x.cells[1].polygon = {P(1, 0), P(2, 1), P(2, 0), P(1, 1)};
    refused(x, "simple");
    x = ok;
    x.regions[0] = {-2};
    refused(x, "void");
    x = ok;
    x.regions = {{ExtrudedGeometry::void_region}, {ExtrudedGeometry::void_region}};
    refused(x, "every cell is void");
    x = ok;
    x.regions[0] = {0, 0};
    refused(x, "one entry per section");
    x = ok;
    x.cells[1].polygon = {P(0.5, 0), P(2, 0), P(2, 1), P(0.5, 1)};
    refused(x, "conforming");
    x = ok;
    x.cells[1].polygon = {P(1, 0.5), P(2, 0.5), P(2, 1.5), P(1, 1.5)};  // a T-junction on a's edge
    refused(x, "not conforming");
    x = ok;
    x.cells.push_back({"c", {P(0.25, 0.25), P(0.75, 0.25), P(0.75, 0.75)}});
    x.regions.push_back({0});
    refused(x, "overlap");
    x = ok;
    x.cells[1].polygon = {P(1, -1), P(2, -1), P(2, 2), P(1, 2)};  // a's edge inside b's longer edge
    refused(x, "not conforming");
    x = ok;
    x.patches = {{"boundary", {0, 0, 0}, {1, 1, 1}}};
    refused(x, "default patch");
    x = ok;
    x.patches = {{"p", {0, 0, 0}, {-1, 1, 1}}};
    refused(x, "lo above hi");
    x = ok;
    x.patches = {{"p", {5, 5, 5}, {6, 6, 6}}};
    refused(x, "selects no boundary facet");
    x = ok;
    x.cells[1].polygon = {P(3, 0), P(4, 0), P(4, 1), P(3, 1)};  // apart: two disconnected volumes
    refused(x, "connected");
}
