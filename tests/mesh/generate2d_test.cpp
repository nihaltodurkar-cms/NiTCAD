// Gates for the 2D generator (Unit 16; DECISIONS.md M3, M5, M7, M8, M11, M13, U2, U12): validity
// (make_dual_mesh accepts every generated mesh), the quality and size bounds, P4 and its certificate
// for every interface feature type, the corner rules, the refusals, determinism, contact areas,
// and speed.
#include <catch2/catch_test_macros.hpp>

#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
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

TEST_CASE("generator: interface features under the corner rules (U12)", "[mesh][generate]") {
    MeshOptions2D o;
    o.max_size = 0.1;
    o.interface_resolution = 1.0 / 64.0;
    SECTION("a 90-degree bend: no corner nodes, exact") {
        Geometry2D g{{{{{0, 0, 0}, {1, 0, 0}, {1, 0.5, 0}, {0.5, 0.5, 0}, {0.5, 1, 0}, {0, 1, 0}}, 0},
                      {{{0.5, 0.5, 0}, {1, 0.5, 0}, {1, 1, 0}, {0.5, 1, 0}}, 1}},
                     {}};
        const auto m = generate_mesh(g, o);
        REQUIRE(m.has_value());
        REQUIRE(m->certificate.has_value());
        CHECK(m->certificate->report.max_corner_deviation <= 1e-12 * o.interface_resolution);  // exact
        CHECK(m->certificate->report.max_flat_deviation <= 1.0);
    }
    SECTION("a T-junction of three regions") {
        Geometry2D g{{{{{0, 0, 0}, {1, 0, 0}, {1, 0.5, 0}, {0.5, 0.5, 0}, {0, 0.5, 0}}, 0}, rect(0, 0.5, 0.5, 1, 1),
                      rect(0.5, 0.5, 1, 1, 2)},
                     {}};
        const auto m = generate_mesh(g, o);
        REQUIRE(m.has_value());
        REQUIRE(m->certificate.has_value());
        CHECK(m->certificate->report.max_corner_deviation <= o.interface_resolution);
    }
    SECTION("a 60.3-degree bend with an oblique boundary junction: accepted only through P4b") {
        // The arm meets y = 1.5 at x = 0.97: about 60.3 degrees (exactly 60 cannot be represented, and
        // a hair below is refused).
        const double top = 0.97;
        Geometry2D g{{{{{0, 0, 0}, {1.5, 0, 0}, {1.5, 0.5, 0}, {0.4, 0.5, 0}, {top, 1.5, 0}, {0, 1.5, 0}}, 0},
                      {{{0.4, 0.5, 0}, {1.5, 0.5, 0}, {1.5, 1.5, 0}, {top, 1.5, 0}}, 1}},
                     {}};
        const auto m = generate_mesh(g, o);
        if (m.has_value()) {
            REQUIRE(m->certificate.has_value());
            CHECK(m->certificate->report.max_corner_deviation <= o.interface_resolution);
        } else {
            // The ruling's other outcome: refused because P4b failed, never silently accepted.
            CHECK(m.error().code == ErrorCode::degenerate_mesh);
            CHECK(m.error().message.find("U12") != std::string::npos);
        }
        INFO((m.has_value() ? "accepted, corner deviation / s = " +
                                  std::to_string(m->certificate->report.max_corner_deviation / o.interface_resolution)
                            : "refused: " + m.error().message));
        CHECK(m.has_value());  // measured: the generator's 60-degree construction meets P4b
    }
    SECTION("a 120-degree bend (measured 0.73 s)") {
        Geometry2D g{{{{{0, 0, 0}, {1.5, 0, 0}, {1.5, 0.5, 0}, {0.6, 0.5, 0}, {0.1, 1.5, 0}, {0, 1.5, 0}}, 0},
                      {{{0.6, 0.5, 0}, {1.5, 0.5, 0}, {1.5, 1.5, 0}, {0.1, 1.5, 0}}, 1}},
                     {}};
        const auto m = generate_mesh(g, o);
        REQUIRE(m.has_value());
        REQUIRE(m->certificate.has_value());
        CHECK(m->certificate->report.max_corner_deviation <= o.interface_resolution);
        CHECK(m->certificate->report.max_corner_deviation > 0.0);  // bisector nodes: not exact
        // The ruled rule (U12): one node per sector on its bisector at s / 2 from the bend at
        // (0.6, 0.5); the arms run toward +x and toward (0.1, 1.5).
        const double s = o.interface_resolution;
        const double ax = 1.0, ay = 0.0;
        const double bl = std::hypot(-0.5, 1.0), bx = -0.5 / bl, by = 1.0 / bl;
        const double dl = std::hypot(ax + bx, ay + by);
        const std::array<Point, 2> expected{Point{0.6 + 0.5 * s * (ax + bx) / dl, 0.5 + 0.5 * s * (ay + by) / dl, 0},
                                            Point{0.6 - 0.5 * s * (ax + bx) / dl, 0.5 - 0.5 * s * (ay + by) / dl, 0}};
        for (const Point& e : expected) {
            bool found = false;
            for (const Point& p : m->simplices.points) found = found || std::hypot(p[0] - e[0], p[1] - e[1]) <= 1e-14;
            CHECK(found);
        }
    }
    SECTION("a 45-degree sector is refused (U12)") {
        Geometry2D g{{{{{0, 0, 0}, {1, 0, 0}, {1, 0.5, 0}, {0.4, 0.5, 0}, {0.9, 1.0, 0}, {0, 1, 0}}, 0},
                      {{{0.4, 0.5, 0}, {1, 0.5, 0}, {1, 1, 0}, {0.9, 1.0, 0}}, 1}},
                     {}};
        const auto m = generate_mesh(g, o);
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
