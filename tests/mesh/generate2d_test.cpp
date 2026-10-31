// Gates for the 2D generator (Unit 16; DECISIONS.md M3, M5, M7, M8, M11, M13, U2, U12): validity
// (make_dual_mesh accepts every generated mesh), the quality and size bounds, P4 and its certificate
// for every interface feature type, the corner rules, the refusals, determinism, contact areas,
// and speed.
#include <catch2/catch_test_macros.hpp>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <algorithm>
#include <array>
#include <string>
#include <vector>

#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/interfaces.hpp"
#include "NiTCAD/mesh/mesh.hpp"

using namespace NiTCAD::mesh;
using NiTCAD::base::ErrorCode;

namespace {

RegionPolygon rect(double x0, double y0, double x1, double y1, std::int32_t region) {
    return {{{x0, y0, 0}, {x1, y0, 0}, {x1, y1, 0}, {x0, y1, 0}}, region};
}

double total(std::span<const double> v) {
    double s = 0.0;
    for (const double x : v) s += x;
    return s;
}

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::max(std::abs(a), std::abs(b)); }

// The circumradius and smallest angle of every generated triangle, against the options.
void check_quality(const GeneratedMesh& m, const MeshOptions2D& o) {
    CHECK(m.report.min_angle_deg >= o.min_angle_deg - 1e-9);
    const auto& p = m.simplices.points;
    double worst = 0.0;
    for (std::size_t c = 0; c < m.simplices.cells.size() / 3; ++c) {
        const Point& a = p[static_cast<std::size_t>(m.simplices.cells[3 * c])];
        const Point& b = p[static_cast<std::size_t>(m.simplices.cells[3 * c + 1])];
        const Point& d = p[static_cast<std::size_t>(m.simplices.cells[3 * c + 2])];
        const double ab = std::hypot(b[0] - a[0], b[1] - a[1]);
        const double bd = std::hypot(d[0] - b[0], d[1] - b[1]);
        const double da = std::hypot(a[0] - d[0], a[1] - d[1]);
        const double area = 0.5 * std::abs((b[0] - a[0]) * (d[1] - a[1]) - (b[1] - a[1]) * (d[0] - a[0]));
        worst = std::max(worst, ab * bd * da / (4.0 * area));
    }
    // Exempt triangles (interface bands) may be larger than h / sqrt 3, but never beyond h.
    CHECK(worst <= o.max_size * (1.0 + 1e-9));
}

}  // namespace

TEST_CASE("generator: a single-region rectangle", "[mesh][generate]") {
    Geometry2D g{{rect(0, 0, 1e-4, 5e-5, 0)}, {{"left", {0, 0, 0}, {0, 5e-5, 0}}, {"right", {1e-4, 0, 0}, {1e-4, 5e-5, 0}}}};
    MeshOptions2D o;
    o.max_size = 1e-5;
    const auto m = generate_mesh(g, o);
    REQUIRE(m.has_value());
    CHECK_FALSE(m->certificate.has_value());
    CHECK(close(total(m->mesh.volumes()), 5e-9, 1e-12));
    REQUIRE(m->mesh.boundary().size() == 2);
    CHECK(close(total(m->mesh.boundary()[0].areas), 5e-5, 1e-12));
    CHECK(close(total(m->mesh.boundary()[1].areas), 5e-5, 1e-12));
    CHECK(m->report.exempt_triangles == 0);
    check_quality(*m, o);
}

TEST_CASE("generator: the same input gives the same mesh", "[mesh][generate]") {
    Geometry2D g{{rect(0, 0, 3, 1, 0), rect(0, 1, 3, 2, 1)}, {}};
    MeshOptions2D o;
    o.max_size = 0.4;
    o.interface_resolution = 0.05;
    const auto a = generate_mesh(g, o);
    const auto b = generate_mesh(g, o);
    REQUIRE(a.has_value());
    REQUIRE(b.has_value());
    CHECK(graph_digest(a->mesh) == graph_digest(b->mesh));
    CHECK(a->simplices.cells == b->simplices.cells);
    CHECK(a->node_region == b->node_region);
}

TEST_CASE("generator: non-convex domains and holes", "[mesh][generate]") {
    MeshOptions2D o;
    o.max_size = 0.1;
    SECTION("an L-shaped domain") {
        Geometry2D g{{{{{0, 0, 0}, {1, 0, 0}, {1, 0.4, 0}, {0.4, 0.4, 0}, {0.4, 1, 0}, {0, 1, 0}}, 0}}, {}};
        const auto m = generate_mesh(g, o);
        REQUIRE(m.has_value());
        CHECK(close(total(m->mesh.volumes()), 1.0 - 0.36, 1e-12));
        check_quality(*m, o);
    }
    SECTION("a frame of four rectangles of one region: a hole, and internal edges") {
        // Conforming: the long rectangles list the corners of the short ones on their edges.
        Geometry2D g{{{{{0, 0, 0}, {1, 0, 0}, {1, 0.3, 0}, {0.7, 0.3, 0}, {0.3, 0.3, 0}, {0, 0.3, 0}}, 0},
                      {{{0, 0.7, 0}, {0.3, 0.7, 0}, {0.7, 0.7, 0}, {1, 0.7, 0}, {1, 1, 0}, {0, 1, 0}}, 0},
                      {{{0, 0.3, 0}, {0.3, 0.3, 0}, {0.3, 0.7, 0}, {0, 0.7, 0}}, 0},
                      {{{0.7, 0.3, 0}, {1, 0.3, 0}, {1, 0.7, 0}, {0.7, 0.7, 0}}, 0}},
                     {}};
        const auto m = generate_mesh(g, o);
        REQUIRE(m.has_value());
        CHECK(close(total(m->mesh.volumes()), 1.0 - 0.16, 1e-12));
        CHECK_FALSE(m->certificate.has_value());
        check_quality(*m, o);
    }
}

TEST_CASE("generator: a flat interface, certified; the device needs the certificate", "[mesh][generate]") {
    Geometry2D g{{rect(0, 0, 1, 0.5, 0), rect(0, 0.5, 1, 1, 1)}, {{"bottom", {0, 0, 0}, {1, 0, 0}}, {"top", {0, 1, 0}, {1, 1, 0}}}};
    MeshOptions2D o;
    o.max_size = 0.1;
    o.interface_resolution = 1.0 / 64.0;
    const auto m = generate_mesh(g, o);
    REQUIRE(m.has_value());
    REQUIRE(m->certificate.has_value());
    CHECK(m->certificate->report.max_flat_deviation <= 1.0);
    CHECK(certifies(*m->certificate, m->mesh, m->node_region));
    for (std::size_t i = 0; i < m->node_region.size(); ++i) {
        CHECK(m->node_region[i] == (m->simplices.points[i][1] < 0.5 ? 0 : 1));
    }
    check_quality(*m, o);
}

namespace {

// No node within one resolution of a feature point: the generator places no corner nodes (U15);
// the nearest are the graded pairs, sqrt((1.5 s)^2 + (s/2)^2) = 1.58 s away.
void check_no_corner_nodes(const GeneratedMesh& m, const Point& v, double s) {
    for (const Point& p : m.simplices.points) CHECK(std::hypot(p[0] - v[0], p[1] - v[1]) > s);
}

}  // namespace

TEST_CASE("generator: every measured corner configuration (U12, U15 regression)", "[mesh][generate]") {
    // Regression bounds from the measurements of 2026-10-14 (DECISIONS.md U15), per configuration;
    // P4b (<= s) still decides acceptance, these catch a regression below it.
    MeshOptions2D o;
    o.max_size = 0.1;
    o.interface_resolution = 1.0 / 64.0;
    const double s = o.interface_resolution;
    const auto certified = [&](const Geometry2D& g) {
        auto m = generate_mesh(g, o);
        REQUIRE(m.has_value());
        REQUIRE(m->certificate.has_value());
        CHECK(m->certificate->report.max_flat_deviation <= 1.0);
        return std::move(*m);
    };
    SECTION("a flat interface: exact") {
        const auto m = certified({{rect(0, 0, 1, 0.5, 0), rect(0, 0.5, 1, 1, 1)}, {}});
        CHECK(m.certificate->report.max_corner_deviation <= 1e-12 * s);
    }
    SECTION("a 90-degree bend: exact, no corner nodes") {
        const auto m = certified({{{{{0, 0, 0}, {1, 0, 0}, {1, 0.5, 0}, {0.5, 0.5, 0}, {0.5, 1, 0}, {0, 1, 0}}, 0},
                                   {{{0.5, 0.5, 0}, {1, 0.5, 0}, {1, 1, 0}, {0.5, 1, 0}}, 1}},
                                  {}});
        CHECK(m.certificate->report.max_corner_deviation <= 1e-12 * s);
        check_no_corner_nodes(m, {0.5, 0.5, 0}, s);
    }
    SECTION("a T-junction of three regions: exact (0.31 s with bisector nodes)") {
        const auto m = certified({{{{{0, 0, 0}, {1, 0, 0}, {1, 0.5, 0}, {0.5, 0.5, 0}, {0, 0.5, 0}}, 0},
                                   rect(0, 0.5, 0.5, 1, 1), rect(0.5, 0.5, 1, 1, 2)},
                                  {}});
        CHECK(m.certificate->report.max_corner_deviation <= 1e-12 * s);
        check_no_corner_nodes(m, {0.5, 0.5, 0}, s);
    }
    SECTION("a 60.3-degree bend with an oblique boundary junction: 0.359 s (0.95 s with bisector nodes)") {
        const double top = 0.97;  // the arm meets y = 1.5 here; exactly 60 cannot be represented
        const auto m = certified({{{{{0, 0, 0}, {1.5, 0, 0}, {1.5, 0.5, 0}, {0.4, 0.5, 0}, {top, 1.5, 0}, {0, 1.5, 0}}, 0},
                                   {{{0.4, 0.5, 0}, {1.5, 0.5, 0}, {1.5, 1.5, 0}, {top, 1.5, 0}}, 1}},
                                  {}});
        CHECK(m.certificate->report.max_corner_deviation <= 0.37 * s);
        check_no_corner_nodes(m, {0.4, 0.5, 0}, s);
    }
    SECTION("a 120-degree bend: 0.321 s (0.73 s with bisector nodes)") {
        const auto m = certified({{{{{0, 0, 0}, {1.5, 0, 0}, {1.5, 0.5, 0}, {0.6, 0.5, 0}, {0.1, 1.5, 0}, {0, 1.5, 0}}, 0},
                                   {{{0.6, 0.5, 0}, {1.5, 0.5, 0}, {1.5, 1.5, 0}, {0.1, 1.5, 0}}, 1}},
                                  {}});
        CHECK(m.certificate->report.max_corner_deviation <= 0.33 * s);
        check_no_corner_nodes(m, {0.6, 0.5, 0}, s);
    }
    SECTION("a 45-degree sector is still refused (U12)") {
        const auto m = generate_mesh({{{{{0, 0, 0}, {1, 0, 0}, {1, 0.5, 0}, {0.4, 0.5, 0}, {0.9, 1.0, 0}, {0, 1, 0}}, 0},
                                       {{{0.4, 0.5, 0}, {1, 0.5, 0}, {1, 1, 0}, {0.9, 1.0, 0}}, 1}},
                                      {}},
                                     o);
        REQUIRE_FALSE(m.has_value());
        CHECK(m.error().code == ErrorCode::invalid_input);
        CHECK(m.error().message.find("U12") != std::string::npos);
    }
}

TEST_CASE("generator: thin layers", "[mesh][generate]") {
    MeshOptions2D o;
    o.max_size = 0.1;
    o.interface_resolution = 0.02;
    SECTION("a layer at least 4 s thick is meshed and certified") {
        Geometry2D g{{rect(0, 0, 1, 0.4, 0), rect(0, 0.4, 1, 0.5, 1), rect(0, 0.5, 1, 1, 2)}, {}};
        const auto m = generate_mesh(g, o);
        REQUIRE(m.has_value());
        REQUIRE(m->certificate.has_value());
        CHECK(m->certificate->report.max_flat_deviation <= 1.0);
    }
    SECTION("a thinner one is refused with a message") {
        Geometry2D g{{rect(0, 0, 1, 0.4, 0), rect(0, 0.4, 1, 0.45, 1), rect(0, 0.45, 1, 1, 2)}, {}};
        const auto m = generate_mesh(g, o);
        REQUIRE_FALSE(m.has_value());
        CHECK(m.error().code == ErrorCode::invalid_input);
        CHECK(m.error().message.find("4 s") != std::string::npos);
    }
}

TEST_CASE("generator: contacts on part of a boundary edge have their exact length", "[mesh][generate]") {
    Geometry2D g{{rect(0, 0, 1, 0.5, 0), rect(0, 0.5, 1, 1, 1)},
                 {{"source", {0.1, 1, 0}, {0.3, 1, 0}}, {"drain", {0.7, 1, 0}, {0.9, 1, 0}}, {"substrate", {0, 0, 0}, {1, 0, 0}}}};
    MeshOptions2D o;
    o.max_size = 0.1;
    o.interface_resolution = 1.0 / 64.0;
    const auto m = generate_mesh(g, o);
    REQUIRE(m.has_value());
    REQUIRE(m->mesh.boundary().size() == 3);
    CHECK(close(total(m->mesh.boundary()[0].areas), 0.2, 1e-12));
    CHECK(close(total(m->mesh.boundary()[1].areas), 0.2, 1e-12));
    CHECK(close(total(m->mesh.boundary()[2].areas), 1.0, 1e-12));
}

TEST_CASE("generator: refusals", "[mesh][generate]") {
    MeshOptions2D o;
    o.max_size = 0.1;
    o.interface_resolution = 1.0 / 64.0;
    const auto refused = [&](const Geometry2D& g, const MeshOptions2D& opt, const std::string& text) {
        const auto m = generate_mesh(g, opt);
        REQUIRE_FALSE(m.has_value());
        CHECK(m.error().code == ErrorCode::invalid_input);
        INFO(m.error().message);
        CHECK(m.error().message.find(text) != std::string::npos);
    };
    refused({{{{{0, 0, 0}, {1, 0, 0}, {2, 0, 0}}, 0}}, {}}, o, "no area");
    refused({{rect(0, 0, 1, 1, 0), rect(0.5, 0, 1.5, 1, 1)}, {}}, o, "overlap");
    refused({{rect(0, 0, 1, 1, 0), rect(1, 0, 2, 0.5, 1)}, {}}, o, "sharing whole edges");  // a T-vertex
    refused({{{{{0, 0, 0}, {2, 0, 0}, {0, 1, 0}, {1, 1, 0}}, 0}}, {}}, o, "simple");          // a bow tie
    refused({{rect(0, 0, 1, 1, 0)}, {{"p", {0.2, 0.5, 0}, {0.4, 0.5, 0}}}}, o, "boundary edge");
    refused({{{{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}}, 0}, {{{0, 0, 0}, {1, 1, 0}, {0, 1, 0}}, 1}}, {}}, o, "domain corner");
    MeshOptions2D no_s = o;
    no_s.interface_resolution = 0.0;
    refused({{rect(0, 0, 1, 0.5, 0), rect(0, 0.5, 1, 1, 1)}, {}}, no_s, "interface_resolution");
    MeshOptions2D bad = o;
    bad.max_size = -1;
    refused({{rect(0, 0, 1, 1, 0)}, {}}, bad, "max_size");
    MeshOptions2D tiny = o;
    tiny.max_nodes = 50;
    const auto capped = generate_mesh({{rect(0, 0, 1, 1, 0)}, {}}, tiny);
    REQUIRE_FALSE(capped.has_value());
    CHECK(capped.error().code == ErrorCode::resource_exhausted);
}

TEST_CASE("generator: speed, 2D (Release budget 1e5 nodes in 2 s)", "[.benchmark][mesh][generate]") {
    Geometry2D g{{rect(0, 0, 1, 0.5, 0), rect(0, 0.5, 1, 1, 1)}, {}};
    MeshOptions2D o;
    o.max_size = 0.0045;
    o.interface_resolution = 0.001;
    const auto t0 = std::chrono::steady_clock::now();
    const auto m = generate_mesh(g, o);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    REQUIRE(m.has_value());
    WARN("2D generator: " << m->report.nodes << " nodes, " << m->report.triangles << " triangles in " << seconds
                          << " s (" << m->report.refinement_points << " refinement points, "
                          << m->report.exempt_triangles << " exempt)");
    CHECK(m->report.nodes > 50000);
}

TEST_CASE("generator: refinement next to an interface respects the protected band", "[mesh][generate]") {
    // A size box hugging the interface forces circumcentres into its band: they must be refused
    // (the triangles exempt), or P4 would fail.
    Geometry2D g{{rect(0, 0, 1, 0.5, 0), rect(0, 0.5, 1, 1, 1)}, {}};
    MeshOptions2D o;
    o.max_size = 0.1;
    o.interface_resolution = 1.0 / 64.0;
    o.boxes = {{{0, 0.4, 0}, {1, 0.6, 0}, 1.0 / 32.0}};
    const auto m = generate_mesh(g, o);
    REQUIRE(m.has_value());
    REQUIRE(m->certificate.has_value());
    CHECK(m->certificate->report.max_flat_deviation <= 1.0);
    CHECK(m->report.exempt_triangles > 0);
}

TEST_CASE("generator: interface spacing rules", "[mesh][generate]") {
    MeshOptions2D o;
    o.max_size = 0.1;
    const auto refused = [&](const Geometry2D& g, double s, const std::string& text) {
        MeshOptions2D opt = o;
        opt.interface_resolution = s;
        const auto m = generate_mesh(g, opt);
        REQUIRE_FALSE(m.has_value());
        CHECK(m.error().code == ErrorCode::invalid_input);
        INFO(m.error().message);
        CHECK(m.error().message.find(text) != std::string::npos);
    };
    // A 0.05 interface segment between a boundary junction and a 90-degree bend: shorter than
    // 6 s = 0.094.
    refused({{{{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0.05, 1, 0}, {0.05, 0.5, 0}, {0, 0.5, 0}}, 0},
              {{{0, 0.5, 0}, {0.05, 0.5, 0}, {0.05, 1, 0}, {0, 1, 0}}, 1}},
             {}},
            1.0 / 64.0, "6 s");
    // Two interfaces 0.03 apart: closer than 4 s = 0.04.
    refused({{rect(0, 0, 1, 0.4, 0), rect(0, 0.4, 1, 0.43, 1), rect(0, 0.43, 1, 1, 2)}, {}}, 0.01,
            "two interfaces closer than 4 s");
    // An interface 0.05 from the bottom edge: closer than 4 s = 0.08 to the boundary.
    refused({{rect(0, 0, 1, 0.05, 0), rect(0, 0.05, 1, 1, 1)}, {}}, 0.02, "closer than 4 s to the boundary");
    // A polygon inside another (no edges cross).
    refused({{rect(0, 0, 1, 1, 0), rect(0.3, 0.3, 0.6, 0.6, 1)}, {}}, 0.01, "inside another");
}

TEST_CASE("generator: patch ends are nodes exactly, on a slanted edge too", "[mesh][generate]") {
    // The top edge runs from (1, 1) to (0, 0.5); the patch's ends are dyadic points on it.
    Geometry2D g{{{{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 0.5, 0}}, 0}},
                 {{"slant", {0.75, 0.875, 0}, {0.25, 0.625, 0}}, {"a", {0.9375, 0.96875, 0}, {0.8125, 0.90625, 0}},
                  {"b", {0.1875, 0.59375, 0}, {0.0625, 0.53125, 0}}}};
    MeshOptions2D o;
    o.max_size = 0.1;
    const auto m = generate_mesh(g, o);
    REQUIRE(m.has_value());
    for (const Point& end : {Point{0.75, 0.875, 0}, Point{0.25, 0.625, 0}, Point{0.9375, 0.96875, 0},
                             Point{0.8125, 0.90625, 0}, Point{0.1875, 0.59375, 0}, Point{0.0625, 0.53125, 0}}) {
        bool found = false;
        for (const Point& p : m->simplices.points) found = found || (p[0] == end[0] && p[1] == end[1]);
        CHECK(found);
    }
    CHECK(close(total(m->mesh.boundary()[0].areas), std::hypot(0.5, 0.25), 1e-12));
}

TEST_CASE("generator: regions that do not touch are certified (no interfaces)", "[mesh][generate]") {
    Geometry2D g{{rect(0, 0, 1, 1, 0), rect(2, 0, 3, 1, 1)}, {}};
    MeshOptions2D o;
    o.max_size = 0.2;  // no interface_resolution: there are no interfaces
    const auto m = generate_mesh(g, o);
    REQUIRE(m.has_value());
    REQUIRE(m->certificate.has_value());
    CHECK(m->certificate->report.cross_faces == 0);
}

TEST_CASE("generator: a structured boundary layer is a tensor lattice in the dual (U14)", "[mesh][generate]") {
    // The unit square with a layer on its bottom edge: columns at 0, 0.25, ..., 1 and rows at
    // 0.05, 0.1, 0.2. Inside the layer the dual is the tensor grid's: a horizontal edge in the
    // row y = 0.05 couples (0.05 + 0.05) / 2, a vertical edge between y = 0 and 0.05 at x = 0.5
    // couples (0.25 + 0.25) / 2; no other node lies in the layer's rectangle.
    Geometry2D g{{rect(0, 0, 1, 1, 0)}, {{"bottom", {0, 0, 0}, {1, 0, 0}}}};
    MeshOptions2D o;
    o.max_size = 0.1;
    o.layers = {{{0, 0, 0}, {1, 0, 0}, {0, 0.25, 0.5, 0.75, 1}, {0.05, 0.1, 0.2}}};
    const auto m = generate_mesh(g, o);
    REQUIRE(m.has_value());
    const auto node_at = [&](double x, double y) {
        for (std::size_t i = 0; i < m->mesh.node_count(); ++i) {
            if (m->mesh.points()[i][0] == x && m->mesh.points()[i][1] == y) return static_cast<NodeId>(i);
        }
        return NodeId{-1};
    };
    const auto coupling = [&](NodeId a, NodeId b) {
        for (const Edge& e : m->mesh.edges()) {
            if (e.first == std::min(a, b) && e.second == std::max(a, b)) return e.coupling_area;
        }
        return -1.0;
    };
    for (const double x : {0.0, 0.25, 0.5, 0.75, 1.0}) {
        for (const double y : {0.0, 0.05, 0.1, 0.2}) CHECK(node_at(x, y) >= 0);
    }
    INFO("couplings: " << coupling(node_at(0.25, 0.05), node_at(0.5, 0.05)) << " " << coupling(node_at(0.5, 0.0), node_at(0.5, 0.05)) << " " << coupling(node_at(0.5, 0.05), node_at(0.5, 0.1)));
    CHECK(close(coupling(node_at(0.25, 0.05), node_at(0.5, 0.05)), 0.05, 1e-12));
    CHECK(close(coupling(node_at(0.5, 0.0), node_at(0.5, 0.05)), 0.25, 1e-12));
    CHECK(close(coupling(node_at(0.5, 0.05), node_at(0.5, 0.1)), 0.25, 1e-12));
    std::size_t inside = 0;
    for (const Point& p : m->simplices.points) inside += (p[1] > 0.0 && p[1] < 0.2 && p[0] > 0.0 && p[0] < 1.0) ? 1 : 0;
    CHECK(inside == 6);  // the lattice's interior nodes only: rows 0.05 and 0.1, three inner columns each
    CHECK(close(total(m->mesh.boundary()[0].areas), 1.0, 1e-12));
}

TEST_CASE("generator: boundary layers take patch ends as columns, and refuse bad input", "[mesh][generate]") {
    MeshOptions2D o;
    o.max_size = 0.1;
    SECTION("a patch end becomes a lattice column") {
        Geometry2D g{{rect(0, 0, 1, 1, 0)}, {{"left", {0, 0, 0}, {0.375, 0, 0}}}};
        o.layers = {{{0, 0, 0}, {1, 0, 0}, {0, 0.5, 1}, {0.05}}};
        const auto m = generate_mesh(g, o);
        REQUIRE(m.has_value());
        bool found = false;
        for (const Point& p : m->simplices.points) found = found || (p[0] == 0.375 && p[1] == 0.05);
        CHECK(found);
        CHECK(close(total(m->mesh.boundary()[0].areas), 0.375, 1e-12));
    }
    const auto refused = [&](const Geometry2D& g, const BoundaryLayer& layer, const std::string& text) {
        MeshOptions2D opt = o;
        opt.layers = {layer};
        opt.interface_resolution = 1.0 / 64.0;
        const auto m = generate_mesh(g, opt);
        REQUIRE_FALSE(m.has_value());
        CHECK(m.error().code == ErrorCode::invalid_input);
        INFO(m.error().message);
        CHECK(m.error().message.find(text) != std::string::npos);
    };
    const Geometry2D square{{rect(0, 0, 1, 1, 0)}, {}};
    refused(square, {{0, 0, 0}, {1, 0, 0}, {0.1, 1}, {0.05}}, "along must run");
    refused(square, {{0, 0, 0}, {1, 0, 0}, {0, 0.6, 0.5, 1}, {0.05}}, "along not strictly increasing");
    refused(square, {{0, 0, 0}, {1, 0, 0}, {0, 1}, {0.1, 0.05}}, "offsets");
    refused(square, {{0, 0.5, 0}, {1, 0.5, 0}, {0, 1}, {0.05}}, "not inside one boundary edge");
    refused(square, {{0, 0, 0}, {1, 0, 0}, {0, 1}, {1.5}}, "leaves the polygon");
    // An interface 0.1 above the layer's top, closer than its band.
    refused({{rect(0, 0, 1, 0.3, 0), rect(0, 0.3, 1, 1, 1)}, {}}, {{0, 0, 0}, {1, 0, 0}, {0, 0.5, 1}, {0.29}}, "band");
    MeshOptions2D two = o;
    two.layers = {{{0, 0, 0}, {1, 0, 0}, {0, 1}, {0.3}}, {{0, 0, 0}, {0, 1, 0}, {0, 1}, {0.3}}};
    const auto overlap = generate_mesh(square, two);
    REQUIRE_FALSE(overlap.has_value());
    CHECK(overlap.error().message.find("overlap") != std::string::npos);
    // A small layer on the left edge lying wholly inside the bottom one: only its corners are in
    // the other.
    MeshOptions2D nested = o;
    nested.layers = {{{0, 0, 0}, {1, 0, 0}, {0, 1}, {0.3}}, {{0, 0.1, 0}, {0, 0.2, 0}, {0, 0.1}, {0.05}}};
    const auto inner = generate_mesh(square, nested);
    REQUIRE_FALSE(inner.has_value());
    CHECK(inner.error().message.find("overlap") != std::string::npos);
    std::swap(nested.layers[0], nested.layers[1]);  // the same, in the other order
    const auto outer = generate_mesh(square, nested);
    REQUIRE_FALSE(outer.has_value());
    CHECK(outer.error().message.find("overlap") != std::string::npos);
}

// Unit 49 (BC1, BC2): a layer on a whole boundary edge with its columns coarser than h. Each end
// column coincides with a domain corner and a patch end; the boundary fill once filled the gap
// from such an end to the next column at h, and the fill nodes joined the first row with slivers.
TEST_CASE("generator: no boundary fill between a layer's columns coarser than h (Unit 49)", "[mesh][generate]") {
    const double W = 1e-4, length = 2e-4, row = 1e-7;
    std::vector<double> offsets;
    for (int i = 1; i <= 200; ++i) offsets.push_back(i * row);
    for (int i = 1; i <= 80; ++i) offsets.push_back(2e-5 + i * 1e-6);
    for (const int columns : {3, 5}) {
        std::vector<double> along;
        for (int i = 0; i < columns; ++i) along.push_back(W * i / (columns - 1));
        const double dt = W / (columns - 1);
        Geometry2D g{{rect(0, 0, length, W, 0)}, {{"x_min", {0, 0, 0}, {0, W, 0}}, {"x_max", {length, 0, 0}, {length, W, 0}}}};
        MeshOptions2D o;
        o.max_size = 1e-5;
        o.layers = {{{0, 0, 0}, {0, W, 0}, along, offsets}};
        const auto m = generate_mesh(g, o);
        REQUIRE(m.has_value());
        INFO(columns << " columns");
        std::size_t on_edge = 0;
        for (const Point& p : m->simplices.points) {
            if (p[0] > offsets.back() * (1 + 1e-12)) continue;
            // Every node in the layer's rectangle is a lattice node (BC2), those on the edge only its columns (BC1).
            const bool in_row = p[0] == 0.0 || std::ranges::any_of(offsets, [&](double v) { return std::abs(p[0] - v) <= 1e-12 * W; });
            const bool in_column = std::ranges::any_of(along, [&](double v) { return std::abs(p[1] - v) <= 1e-12 * W; });
            CHECK((in_row && in_column));
            on_edge += p[0] == 0.0 ? 1 : 0;
        }
        CHECK(on_edge == along.size());
        for (std::size_t i = 0; i < m->mesh.node_count(); ++i) {
            const Point& p = m->mesh.points()[i];
            if (p[0] != 0.0) continue;
            const bool end = p[1] == 0.0 || p[1] == W;
            CHECK(close(m->mesh.volumes()[i], (end ? 0.5 : 1.0) * dt * 0.5 * row, 1e-12));
        }
        CHECK(close(total(m->mesh.boundary()[0].areas), W, 1e-12));
    }
    // The U14 square: a layer on the whole bottom edge, columns 0.25 apart, h 0.1.
    Geometry2D g{{rect(0, 0, 1, 1, 0)}, {{"bottom", {0, 0, 0}, {1, 0, 0}}}};
    MeshOptions2D o;
    o.max_size = 0.1;
    o.layers = {{{0, 0, 0}, {1, 0, 0}, {0, 0.25, 0.5, 0.75, 1}, {0.05, 0.1, 0.2}}};
    const auto m = generate_mesh(g, o);
    REQUIRE(m.has_value());
    std::size_t on_edge = 0;
    for (const Point& p : m->simplices.points) on_edge += p[1] == 0.0 ? 1 : 0;
    CHECK(on_edge == 5);
}

// Unit 49 (BC3): a column's node on the edge is the given point, the layer's end and a patch end,
// although a + t u differs from it in the last bit; no second node next to it.
TEST_CASE("generator: a layer's end and patch-end columns are the given points (Unit 49)", "[mesh][generate]") {
    const Point a{0.2, 0, 0}, b{0.9, 0, 0}, q{0.4, 0, 0};
    // The case: the lattice's own construction, a + t u, misses b in the last bit.
    const double len = std::hypot(b[0] - a[0], b[1] - a[1]);
    REQUIRE(a[0] + len * ((b[0] - a[0]) / len) != b[0]);
    Geometry2D g{{rect(0, 0, 1, 1, 0)}, {{"mid", q, b}}};
    MeshOptions2D o;
    o.max_size = 0.1;
    o.layers = {{a, b, {0, 0.7}, {0.05}}};
    const auto m = generate_mesh(g, o);
    REQUIRE(m.has_value());
    for (const Point& given : {a, b, q}) {
        INFO(given[0]);
        std::size_t exact = 0, near = 0;
        for (const Point& p : m->simplices.points) {
            exact += (p[0] == given[0] && p[1] == given[1]) ? 1 : 0;
            near += std::hypot(p[0] - given[0], p[1] - given[1]) <= 1e-9 * 0.7 ? 1 : 0;
        }
        CHECK(exact == 1);
        CHECK(near == 1);
    }
    CHECK(close(total(m->mesh.boundary()[0].areas), 0.5, 1e-12));
}

// UC1 (Unit 54): two U14 layers near the square's corner, on the bottom and left edges. Delaunay joins
// their nodes across the corner by triangles of layer nodes lying in neither lattice, which U14 exempts
// from refinement; signed layers exempt only triangles inside a lattice (Unit 54 slice 3, addition (d)),
// and U14's layers keep U14's rule. Pinned: the digest measured before Unit 53 (fea05e5), unchanged since.
TEST_CASE("generator: two boundary layers near a corner keep U14's mesh (Unit 54, UC1)", "[mesh][generate]") {
    Geometry2D g{{rect(0, 0, 1, 1, 0)}, {}};
    MeshOptions2D o;
    o.max_size = 0.1;
    o.layers = {{{0.3, 0, 0}, {1, 0, 0}, {0, 0.35, 0.7}, {0.05, 0.1}}, {{0, 0.3, 0}, {0, 1, 0}, {0, 0.35, 0.7}, {0.05, 0.1}}};
    const auto m = generate_mesh(g, o);
    REQUIRE(m.has_value());
    CHECK(m->mesh.node_count() == 215);
    CHECK(graph_digest(m->mesh) == 0x1c2949ce70fe860aULL);
}
