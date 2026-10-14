// Gates for the 3D generator (Unit 16, commit 6b; DECISIONS.md M3, M6, M10, M11, M13): validity
// (make_dual_mesh accepts every generated mesh), volumes and face areas, the quality report against
// a recomputation, determinism, the refusals, and speed.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/mesh.hpp"

using namespace NiTCAD::mesh;
using NiTCAD::base::ErrorCode;

namespace {

double total(std::span<const double> v) {
    double s = 0.0;
    for (const double x : v) s += x;
    return s;
}

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::max(std::abs(a), std::abs(b)); }

Geometry3D box(double x, double y, double z) {
    Geometry3D g;
    g.domain.vertices = {{0, 0, 0}, {x, 0, 0}, {x, y, 0}, {0, y, 0}, {0, 0, z}, {x, 0, z}, {x, y, z}, {0, y, z}};
    g.domain.faces = {{{0, 3, 7, 4}, "x_min"}, {{1, 2, 6, 5}, "x_max"}, {{0, 1, 5, 4}, "y_min"},
                      {{3, 2, 6, 7}, "y_max"}, {{0, 1, 2, 3}, "z_min"}, {{4, 5, 6, 7}, "z_max"}};
    return g;
}

// A regular tetrahedron on alternate cube corners: every face oblique, dihedral angles 70.53 degrees.
Geometry3D regular_tetrahedron() {
    Geometry3D g;
    g.domain.vertices = {{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}};
    g.domain.faces = {{{1, 2, 3}, "f0"}, {{0, 2, 3}, "f1"}, {{0, 1, 3}, "f2"}, {{0, 1, 2}, "f3"}};
    return g;
}

// A triangular prism, z in [0, 1], over a nearly equilateral triangle: 60-degree vertical edges, two
// oblique faces.
Geometry3D prism() {
    Geometry3D g;
    const double a = 1.75;
    g.domain.vertices = {{0, 0, 0}, {2, 0, 0}, {1, a, 0}, {0, 0, 1}, {2, 0, 1}, {1, a, 1}};
    g.domain.faces = {{{0, 1, 2}, "bottom"}, {{3, 4, 5}, "top"}, {{0, 1, 4, 3}, "south"},
                      {{1, 2, 5, 4}, "east"}, {{2, 0, 3, 5}, "west"}};
    return g;
}

double tet_volume(const SimplexMesh& s, std::size_t c) {
    const auto P = [&](std::size_t k) -> const Point& { return s.points[static_cast<std::size_t>(s.cells[4 * c + k])]; };
    const Point u{P(1)[0] - P(0)[0], P(1)[1] - P(0)[1], P(1)[2] - P(0)[2]};
    const Point v{P(2)[0] - P(0)[0], P(2)[1] - P(0)[1], P(2)[2] - P(0)[2]};
    const Point w{P(3)[0] - P(0)[0], P(3)[1] - P(0)[1], P(3)[2] - P(0)[2]};
    return std::abs(u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) +
                    u[2] * (v[0] * w[1] - v[1] * w[0])) / 6.0;
}

// The generated mesh's volumes, face areas and report against the geometry and options.
void check_mesh(const GeneratedMesh& m, double volume, const std::vector<double>& face_areas, const MeshOptions3D& o) {
    CHECK(m.mesh.dimension() == 3);
    CHECK(m.report.nodes == m.mesh.node_count());
    CHECK(m.report.tetrahedra == m.simplices.cells.size() / 4);
    CHECK(close(total(m.mesh.volumes()), volume, 1e-12));
    double cells = 0.0;
    for (std::size_t c = 0; c < m.report.tetrahedra; ++c) cells += tet_volume(m.simplices, c);
    CHECK(close(cells, volume, 1e-12));
    REQUIRE(m.mesh.boundary().size() == face_areas.size());
    for (std::size_t f = 0; f < face_areas.size(); ++f) CHECK(close(total(m.mesh.boundary()[f].areas), face_areas[f], 1e-12));
    CHECK(m.report.min_dihedral_deg > 0.0);
    CHECK(m.report.min_dihedral_deg < 70.6);
    // Radius-edge: the bound holds unless some tetrahedra were exempt.
    if (m.report.exempt_tetrahedra == 0) CHECK(m.report.max_radius_edge <= o.radius_edge * (1.0 + 1e-9));
    for (const auto r : m.node_region) CHECK(r == 0);
}

// Every polyhedron edge is split into subsegments that are edges of the tetrahedralization, with no
// node strictly inside their diametral balls (the encroachment rule's postcondition).
void check_edges(const GeneratedMesh& m, const Geometry3D& g) {
    const auto& pts = m.simplices.points;
    std::set<std::pair<NodeId, NodeId>> cell_edges;
    for (std::size_t c = 0; c < m.simplices.cells.size(); c += 4) {
        for (std::size_t i = 0; i < 4; ++i) {
            for (std::size_t j = i + 1; j < 4; ++j) {
                const NodeId a = m.simplices.cells[c + i];
                const NodeId b = m.simplices.cells[c + j];
                cell_edges.insert({std::min(a, b), std::max(a, b)});
            }
        }
    }
    std::set<std::pair<std::size_t, std::size_t>> edges;
    for (const PolyhedronFace& f : g.domain.faces) {
        for (std::size_t k = 0; k < f.vertices.size(); ++k) {
            const std::size_t a = f.vertices[k];
            const std::size_t b = f.vertices[(k + 1) % f.vertices.size()];
            edges.insert({std::min(a, b), std::max(a, b)});
        }
    }
    for (const auto& [ia, ib] : edges) {
        const Point& A = g.domain.vertices[ia];
        const Point& B = g.domain.vertices[ib];
        const Point d{B[0] - A[0], B[1] - A[1], B[2] - A[2]};
        const double L2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        std::vector<std::pair<double, NodeId>> on;
        for (std::size_t i = 0; i < pts.size(); ++i) {
            const Point r{pts[i][0] - A[0], pts[i][1] - A[1], pts[i][2] - A[2]};
            const double t = (r[0] * d[0] + r[1] * d[1] + r[2] * d[2]) / L2;
            const Point off{r[0] - t * d[0], r[1] - t * d[1], r[2] - t * d[2]};
            if (t >= -1e-12 && t <= 1 + 1e-12 && off[0] * off[0] + off[1] * off[1] + off[2] * off[2] <= 1e-24 * L2) {
                on.emplace_back(t, static_cast<NodeId>(i));
            }
        }
        std::ranges::sort(on);
        REQUIRE(on.size() >= 2);
        for (std::size_t k = 0; k + 1 < on.size(); ++k) {
            const NodeId a = on[k].second;
            const NodeId b = on[k + 1].second;
            CHECK(cell_edges.contains({std::min(a, b), std::max(a, b)}));
            const Point& P = pts[static_cast<std::size_t>(a)];
            const Point& Q = pts[static_cast<std::size_t>(b)];
            double least = 0.0;  // (P - X).(Q - X) < 0: X strictly inside the diametral ball
            for (std::size_t i = 0; i < pts.size(); ++i) {
                if (static_cast<NodeId>(i) == a || static_cast<NodeId>(i) == b) continue;
                const Point& X = pts[i];
                least = std::min(least, (P[0] - X[0]) * (Q[0] - X[0]) + (P[1] - X[1]) * (Q[1] - X[1]) + (P[2] - X[2]) * (Q[2] - X[2]));
            }
            CHECK(least >= 0.0);
        }
    }
}

}  // namespace

TEST_CASE("generator 3D: a box", "[mesh][generate3d]") {
    MeshOptions3D o;
    o.max_size = 0.25;
    const auto m = generate_mesh(box(1, 1, 1), o);
    REQUIRE(m.has_value());
    check_mesh(*m, 1.0, {1, 1, 1, 1, 1, 1}, o);
    check_edges(*m, box(1, 1, 1));
    CHECK(m->report.nodes > 64);
    CHECK(m->report.nodes < 2000);
    INFO("nodes " << m->report.nodes << " tets " << m->report.tetrahedra << " slivers " << m->report.slivers
                  << " min dihedral " << m->report.min_dihedral_deg << " max ratio " << m->report.max_radius_edge);
    CHECK(m->report.exempt_tetrahedra == 0);
}

TEST_CASE("generator 3D: a flat box and a size box", "[mesh][generate3d]") {
    MeshOptions3D o;
    o.max_size = 0.5;
    o.boxes = {{{0, 0, 0}, {0.5, 1, 0.25}, 0.1}};
    const auto m = generate_mesh(box(2, 1, 0.25), o);
    REQUIRE(m.has_value());
    check_mesh(*m, 0.5, {0.25, 0.25, 0.5, 0.5, 2, 2}, o);
    check_edges(*m, box(2, 1, 0.25));
    // Finer in the size box: more nodes per unit volume there than outside.
    std::size_t in = 0;
    for (const Point& p : m->mesh.points()) in += p[0] < 0.5 ? 1 : 0;
    CHECK(4 * in > m->mesh.node_count());
}

TEST_CASE("generator 3D: a regular tetrahedron (oblique faces, 70.53-degree dihedral angles)", "[mesh][generate3d]") {
    MeshOptions3D o;
    o.max_size = 0.6;
    const auto m = generate_mesh(regular_tetrahedron(), o);
    REQUIRE(m.has_value());
    // Edge 2 sqrt 2: volume 8/3, face area 2 sqrt 3.
    const double face = 2.0 * std::sqrt(3.0);
    check_mesh(*m, 8.0 / 3.0, {face, face, face, face}, o);
    check_edges(*m, regular_tetrahedron());
}

TEST_CASE("generator 3D: a prism with 60-degree edges (heuristic, capped)", "[mesh][generate3d]") {
    MeshOptions3D o;
    o.max_size = 0.3;
    const auto m = generate_mesh(prism(), o);
    REQUIRE(m.has_value());
    const double side = std::hypot(1.0, 1.75);
    check_mesh(*m, 1.75, {1.75, 1.75, 2.0, side, side}, o);
    check_edges(*m, prism());
}

TEST_CASE("generator 3D: oblique faces at finer sizes (nodes rounded off their planes)", "[mesh][generate3d]") {
    // An octahedron: eight oblique faces, 109.47-degree dihedral angles.
    Geometry3D octa;
    octa.domain.vertices = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
    octa.domain.faces = {{{0, 2, 4}, "a"}, {{2, 1, 4}, "b"}, {{1, 3, 4}, "c"}, {{3, 0, 4}, "d"},
                         {{2, 0, 5}, "e"}, {{1, 2, 5}, "f"}, {{3, 1, 5}, "g"}, {{0, 3, 5}, "h"}};
    const double tri_area = std::sqrt(3.0) / 2.0;
    struct Case {
        Geometry3D g;
        double h;
        double volume;
        std::vector<double> areas;
    };
    // An irregular tetrahedron at device scale: its edge midpoints are not exact binary fractions,
    // so split points on edges and oblique faces round off their lines and planes.
    Geometry3D small;
    small.domain.vertices = {{0, 0, 0}, {3e-5, 0, 0}, {1e-5, 2.7e-5, 0}, {1.3e-5, 0.9e-5, 2.3e-5}};
    small.domain.faces = {{{1, 2, 3}, "f0"}, {{0, 2, 3}, "f1"}, {{0, 1, 3}, "f2"}, {{0, 1, 2}, "f3"}};
    const auto tri_area_of = [&](std::size_t i, std::size_t j, std::size_t k) {
        const Point& a = small.domain.vertices[i];
        const Point& b = small.domain.vertices[j];
        const Point& c = small.domain.vertices[k];
        const Point u{b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        const Point v{c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        return 0.5 * std::hypot(u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]);
    };
    const double small_volume = 3e-5 * 2.7e-5 / 2.0 * 2.3e-5 / 3.0;
    const double side = std::hypot(1.0, 1.75);
    const double face = 2.0 * std::sqrt(3.0);
    for (const Case& c : {Case{small, 0.4e-5, small_volume, {tri_area_of(1, 2, 3), tri_area_of(0, 2, 3), tri_area_of(0, 1, 3), tri_area_of(0, 1, 2)}},
                          Case{octa, 0.25, 4.0 / 3.0, std::vector<double>(8, tri_area)},
                          Case{regular_tetrahedron(), 0.3, 8.0 / 3.0, {face, face, face, face}},
                          Case{prism(), 0.15, 1.75, {1.75, 1.75, 2.0, side, side}}}) {
        MeshOptions3D o;
        o.max_size = c.h;
        const auto m = generate_mesh(c.g, o);
        REQUIRE(m.has_value());
        std::printf("oblique: %zu nodes, %zu tetrahedra, %zu flat dropped, %zu slivers, min dihedral %.2f\n",
                    m->report.nodes, m->report.tetrahedra, m->report.flat_dropped, m->report.slivers,
                    m->report.min_dihedral_deg);
        check_mesh(*m, c.volume, c.areas, o);
        check_edges(*m, c.g);
    }
}

TEST_CASE("generator 3D: the same input gives the same mesh", "[mesh][generate3d]") {
    MeshOptions3D o;
    o.max_size = 0.3;
    const auto a = generate_mesh(prism(), o);
    const auto b = generate_mesh(prism(), o);
    REQUIRE(a.has_value());
    REQUIRE(b.has_value());
    CHECK(a->simplices.points == b->simplices.points);
    CHECK(a->simplices.cells == b->simplices.cells);
    CHECK(a->simplices.facets == b->simplices.facets);
    CHECK(graph_digest(a->mesh) == graph_digest(b->mesh));
}

TEST_CASE("generator 3D: refusals", "[mesh][generate3d]") {
    MeshOptions3D o;
    o.max_size = 0.5;
    const auto refused = [&](const Geometry3D& g, const MeshOptions3D& opt, ErrorCode code, const std::string& text) {
        const auto m = generate_mesh(g, opt);
        REQUIRE_FALSE(m.has_value());
        INFO(m.error().message);
        CHECK(m.error().code == code);
        CHECK(m.error().message.find(text) != std::string::npos);
    };
    SECTION("options") {
        MeshOptions3D bad = o;
        bad.max_size = 0.0;
        refused(box(1, 1, 1), bad, ErrorCode::invalid_input, "max_size");
        bad = o;
        bad.radius_edge = 1.2;
        refused(box(1, 1, 1), bad, ErrorCode::invalid_input, "radius_edge");
        bad = o;
        bad.boxes = {{{1, 0, 0}, {0, 1, 1}, 0.1}};
        refused(box(1, 1, 1), bad, ErrorCode::invalid_input, "size box");
        Geometry3D g = box(1, 1, 1);
        g.region = -1;
        refused(g, o, ErrorCode::invalid_input, "region");
    }
    SECTION("geometry") {
        Geometry3D g = box(1, 1, 1);
        g.domain.faces.pop_back();
        refused(g, o, ErrorCode::invalid_input, "not closed");
        g = box(1, 1, 1);
        g.domain.faces[1].name = "x_min";
        refused(g, o, ErrorCode::invalid_input, "unique");
        g = box(1, 1, 1);
        g.domain.vertices[6] = {1, 1, 1.25};  // the top and two sides are no longer planar
        refused(g, o, ErrorCode::invalid_input, "planar");
        g = box(1, 1, 1);
        g.domain.faces[0].vertices = {0, 7, 3, 4};  // self-intersecting quadrilateral
        refused(g, o, ErrorCode::invalid_input, "convex and simple");
        g = box(1, 1, 1);
        g.domain.vertices[0][0] = std::nan("");
        refused(g, o, ErrorCode::invalid_input, "not finite");
        g = box(1, 1, 1);
        g.domain.faces[2].vertices = {0, 1, 5, 99};
        refused(g, o, ErrorCode::invalid_input, "out of range");
    }
    SECTION("non-convex domains (U7)") {
        // An L-shaped prism: planar faces, a re-entrant edge.
        Geometry3D g;
        g.domain.vertices = {{0, 0, 0}, {2, 0, 0}, {2, 1, 0}, {1, 1, 0}, {1, 2, 0}, {0, 2, 0},
                             {0, 0, 1}, {2, 0, 1}, {2, 1, 1}, {1, 1, 1}, {1, 2, 1}, {0, 2, 1}};
        g.domain.faces = {{{0, 1, 2, 3, 4, 5}, "bottom"}, {{6, 7, 8, 9, 10, 11}, "top"}, {{0, 1, 7, 6}, "s"},
                          {{1, 2, 8, 7}, "e1"}, {{2, 3, 9, 8}, "n1"}, {{3, 4, 10, 9}, "e2"},
                          {{4, 5, 11, 10}, "n2"}, {{5, 0, 6, 11}, "w"}};
        refused(g, o, ErrorCode::invalid_input, "");
    }
    SECTION("the node cap") {
        MeshOptions3D tight = o;
        tight.max_size = 0.05;
        tight.max_nodes = 200;
        refused(box(1, 1, 1), tight, ErrorCode::resource_exhausted, "max_nodes");
    }
}

TEST_CASE("generator 3D: speed (Release budget 1e5 nodes in 15 s)", "[.benchmark][mesh][generate3d]") {
    MeshOptions3D o;
    o.max_size = 0.03;
    o.max_nodes = 400'000;
    const auto start = std::chrono::steady_clock::now();
    const auto m = generate_mesh(box(1, 1, 1), o);
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    REQUIRE(m.has_value());
    std::printf("generator 3D: %zu nodes, %zu tetrahedra in %.2f s; slivers %zu, min dihedral %.2f, max ratio %.3f, "
                "exempt %zu\n",
                m->report.nodes, m->report.tetrahedra, seconds, m->report.slivers, m->report.min_dihedral_deg,
                m->report.max_radius_edge, m->report.exempt_tetrahedra);
    // The import checks alone (make_dual_mesh: P0-P3, P5), part of the time above.
    const auto t0 = std::chrono::steady_clock::now();
    const auto again = make_dual_mesh(m->simplices);
    const double dual_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    REQUIRE(again.has_value());
    std::printf("generator 3D: of which make_dual_mesh %.2f s\n", dual_seconds);
}
