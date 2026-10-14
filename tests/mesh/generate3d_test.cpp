// Gates for the 3D generator (Unit 16, commit 6b; DECISIONS.md M3, M6, M10, M11, M13): validity
// (make_dual_mesh accepts every generated mesh), volumes and face areas, the quality report against
// a recomputation, determinism, the refusals, and speed.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <map>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/interfaces.hpp"
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

// The total volume only (meshes with interfaces have regions other than 0).
void check_mesh_volume(const GeneratedMesh& m, double volume) {
    CHECK(close(total(m.mesh.volumes()), volume, 1e-12));
    double cells = 0.0;
    for (std::size_t c = 0; c < m.report.tetrahedra; ++c) cells += tet_volume(m.simplices, c);
    CHECK(close(cells, volume, 1e-12));
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

namespace {

// The M9 gate: each node's area in each patch of an axis-aligned face (the constant coordinate
// `axis`) against an independent 2D computation, make_dual_mesh on that patch's boundary
// triangles in face coordinates; and the patch's total against its exact area.
void check_patch_duals(const GeneratedMesh& m, int axis, const std::vector<std::pair<std::string, double>>& patches) {
    const auto& s = m.simplices;
    for (const auto& [name, area] : patches) {
        std::size_t p = 0;
        while (p < s.patch_names.size() && s.patch_names[p] != name) ++p;
        REQUIRE(p < s.patch_names.size());
        SimplexMesh flat;
        flat.dimension = 2;
        std::map<NodeId, NodeId> local;
        std::vector<NodeId> global;
        for (std::size_t k = 0; k < s.facet_patch.size(); ++k) {
            if (s.facet_patch[k] != static_cast<std::int32_t>(p)) continue;
            for (std::size_t j = 0; j < 3; ++j) {
                const NodeId g = s.facets[3 * k + j];
                auto [it, fresh] = local.try_emplace(g, static_cast<NodeId>(global.size()));
                if (fresh) {
                    global.push_back(g);
                    const Point& q = s.points[static_cast<std::size_t>(g)];
                    flat.points.push_back({q[static_cast<std::size_t>((axis + 1) % 3)], q[static_cast<std::size_t>((axis + 2) % 3)], 0.0});
                }
                flat.cells.push_back(it->second);
            }
        }
        const auto dual = make_dual_mesh(flat);
        REQUIRE(dual.has_value());
        const BoundaryPatch* bp = m.mesh.find_boundary(name);
        REQUIRE(bp != nullptr);
        CHECK(close(total(bp->areas), area, 1e-12));
        CHECK(close(total(dual->volumes()), area, 1e-12));
        // Node by node (nodes whose area in the patch is tau-zero are dropped from the 3D patch).
        std::map<NodeId, double> by_node;
        for (std::size_t k = 0; k < bp->nodes.size(); ++k) by_node[bp->nodes[k]] = bp->areas[k];
        double worst = 0.0;
        for (std::size_t i = 0; i < global.size(); ++i) {
            const auto it = by_node.find(global[i]);
            const double three_d = it == by_node.end() ? 0.0 : it->second;
            worst = std::max(worst, std::abs(three_d - dual->volumes()[i]));
        }
        CHECK(worst <= 1e-12 * area);
    }
}

Geometry3D box_with_top(std::vector<FacePatch> top) {
    Geometry3D g = box(3, 1, 1);
    g.domain.faces[5].patches = std::move(top);
    return g;
}

}  // namespace

TEST_CASE("generator 3D: contacts on one face are exact patches (M9)", "[mesh][generate3d]") {
    MeshOptions3D o;
    o.max_size = 0.35;
    SECTION("source, gate and drain strips across the top face") {
        const auto m = generate_mesh(box_with_top({{"source", {{0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}}},
                                                   {"gate", {{1, 0, 1}, {2, 0, 1}, {2, 1, 1}, {1, 1, 1}}},
                                                   {"drain", {{2, 0, 1}, {3, 0, 1}, {3, 1, 1}, {2, 1, 1}}}}),
                                     o);
        REQUIRE(m.has_value());
        check_mesh(*m, 3.0, {1, 1, 3, 3, 3, 1, 1, 1}, o);
        check_patch_duals(*m, 2, {{"source", 1.0}, {"gate", 1.0}, {"drain", 1.0}});
    }
    SECTION("a contact inside the face, the rest tiled around it") {
        const auto m = generate_mesh(
            box_with_top({{"contact", {{1, 0.25, 1}, {2, 0.25, 1}, {2, 0.75, 1}, {1, 0.75, 1}}},
                          {"left", {{0, 0, 1}, {1, 0, 1}, {1, 0.25, 1}, {1, 0.75, 1}, {1, 1, 1}, {0, 1, 1}}},
                          {"right", {{2, 0, 1}, {3, 0, 1}, {3, 1, 1}, {2, 1, 1}, {2, 0.75, 1}, {2, 0.25, 1}}},
                          {"front", {{1, 0, 1}, {2, 0, 1}, {2, 0.25, 1}, {1, 0.25, 1}}},
                          {"back", {{1, 0.75, 1}, {2, 0.75, 1}, {2, 1, 1}, {1, 1, 1}}}}),
            o);
        REQUIRE(m.has_value());
        check_mesh(*m, 3.0, {1, 1, 3, 3, 3, 0.5, 1, 1, 0.25, 0.25}, o);
        check_patch_duals(*m, 2, {{"contact", 0.5}, {"left", 1.0}, {"right", 1.0}, {"front", 0.25}, {"back", 0.25}});
    }
    SECTION("a slanted patch edge (its split points round off the line)") {
        const auto m = generate_mesh(box_with_top({{"a", {{0, 0, 1}, {0.3, 0, 1}, {0.7, 1, 1}, {0, 1, 1}}},
                                                   {"b", {{0.3, 0, 1}, {3, 0, 1}, {3, 1, 1}, {0.7, 1, 1}}}}),
                                     o);
        REQUIRE(m.has_value());
        check_mesh(*m, 3.0, {1, 1, 3, 3, 3, 0.5, 2.5}, o);
        check_patch_duals(*m, 2, {{"a", 0.5}, {"b", 2.5}});
    }
}

TEST_CASE("generator 3D: face tilings are refused unless exact and conforming", "[mesh][generate3d]") {
    MeshOptions3D o;
    o.max_size = 0.5;
    const auto refused = [&](std::vector<FacePatch> top, const std::string& text) {
        const auto m = generate_mesh(box_with_top(std::move(top)), o);
        REQUIRE_FALSE(m.has_value());
        INFO(m.error().message);
        CHECK(m.error().code == ErrorCode::invalid_input);
        CHECK(m.error().message.find(text) != std::string::npos);
    };
    const FacePatch whole{"whole", {{0, 0, 1}, {3, 0, 1}, {3, 1, 1}, {0, 1, 1}}};
    refused({{"low", {{0, 0, 1}, {3, 0, 1}, {3, 1, 0.5}, {0, 1, 1}}}}, "plane");
    refused({{"out", {{0, 0, 1}, {4, 0, 1}, {4, 1, 1}, {0, 1, 1}}}}, "outside");
    refused({{"half", {{0, 0, 1}, {1.5, 0, 1}, {1.5, 1, 1}, {0, 1, 1}}}}, "tile");
    // A T-vertex: the right polygon's left edge is split where the two left polygons meet.
    refused({{"l1", {{0, 0, 1}, {1, 0, 1}, {1, 0.5, 1}, {0, 0.5, 1}}},
             {"l2", {{0, 0.5, 1}, {1, 0.5, 1}, {1, 1, 1}, {0, 1, 1}}},
             {"r", {{1, 0, 1}, {3, 0, 1}, {3, 1, 1}, {1, 1, 1}}}},
            "shared whole");
    refused({whole, {"whole", {{0, 0, 1}, {1, 0, 1}, {1, 1, 1}}}}, "unique");
    refused({{"x_min", whole.vertices}}, "unique");
    refused({{"thin", {{0, 0, 1}, {3, 0, 1}, {1.5, 0, 1}}}}, "no area");
}

namespace {

InterfacePlane plane_through(Point a, Point b, Point c) { return InterfacePlane{{a, b, c}}; }

// The unit cube cut by planes, with seeds at the given points.
Geometry3D cut_cube(std::vector<InterfacePlane> planes, std::vector<RegionSeed> seeds) {
    Geometry3D g = box(1, 1, 1);
    g.interfaces = std::move(planes);
    g.regions = std::move(seeds);
    return g;
}

// The generated mesh with interfaces: certified (P4), each node in its seed's slab, the
// certificate matching; deviations printed for the record.
// `corner_bound`: the regression bound for this configuration (DECISIONS.md U15 style), from the
// measured corner deviation.
void check_interfaces_3d(const GeneratedMesh& m, double s, const char* label, double corner_bound) {
    REQUIRE(m.certificate.has_value());
    CHECK(certifies(*m.certificate, m.mesh, m.node_region));
    const InterfaceReport& r = m.certificate->report;
    CHECK(r.cross_faces > 0);
    CHECK(r.max_flat_deviation <= 1.0);
    CHECK(r.max_corner_deviation <= s);
    CHECK(r.max_corner_deviation <= corner_bound * s);
    std::printf("P4 3D %s: %zu nodes, %zu cross faces (%zu in corner zones), flat %.3g of tolerance, corner %.3f s, "
                "band facet points %zu, exempt %zu\n",
                label, m.report.nodes, r.cross_faces, r.corner_faces, r.max_flat_deviation, r.max_corner_deviation / s,
                m.report.band_facet_points, m.report.exempt_tetrahedra);
}

}  // namespace

TEST_CASE("generator 3D: planar interfaces are certified (M7, P4)", "[mesh][generate3d]") {
    MeshOptions3D o;
    o.max_size = 0.25;
    o.interface_resolution = 0.02;
    SECTION("one horizontal plane") {
        const auto m = generate_mesh(cut_cube({plane_through({0, 0, 0.5}, {1, 0, 0.5}, {0, 1, 0.5})},
                                              {{{0.5, 0.5, 0.2}, 0}, {{0.5, 0.5, 0.8}, 1}}),
                                     o);
        REQUIRE(m.has_value());
        check_mesh_volume(*m, 1.0);
        check_interfaces_3d(*m, o.interface_resolution, "horizontal", 1e-12);  // measured 0
        for (std::size_t i = 0; i < m->mesh.node_count(); ++i) CHECK(m->node_region[i] == (m->mesh.points()[i][2] < 0.5 ? 0 : 1));
    }
    SECTION("three slabs: a thin middle layer (5 s)") {
        const auto m = generate_mesh(cut_cube({plane_through({0, 0, 0.4}, {1, 0, 0.4}, {0, 1, 0.4}),
                                               plane_through({0, 0, 0.5}, {1, 0, 0.5}, {0, 1, 0.5})},
                                              {{{0.5, 0.5, 0.2}, 0}, {{0.5, 0.5, 0.45}, 1}, {{0.5, 0.5, 0.8}, 2}}),
                                     o);
        REQUIRE(m.has_value());
        check_mesh_volume(*m, 1.0);
        check_interfaces_3d(*m, o.interface_resolution, "thin layer", 1e-12);  // measured 0
    }
    SECTION("an oblique plane") {
        const auto m = generate_mesh(cut_cube({plane_through({0, 0, 0.4}, {1, 0, 0.5}, {0, 1, 0.45})},
                                              {{{0.5, 0.5, 0.2}, 3}, {{0.5, 0.5, 0.8}, 1}}),
                                     o);
        REQUIRE(m.has_value());
        check_mesh_volume(*m, 1.0);
        check_interfaces_3d(*m, o.interface_resolution, "oblique", 0.19);  // measured 0.180
    }
}

TEST_CASE("generator 3D: interface refusals", "[mesh][generate3d]") {
    MeshOptions3D o;
    o.max_size = 0.25;
    o.interface_resolution = 0.02;
    const auto refused = [&](const Geometry3D& g, const MeshOptions3D& opt, const std::string& text) {
        const auto m = generate_mesh(g, opt);
        REQUIRE_FALSE(m.has_value());
        INFO(m.error().message);
        CHECK(m.error().code == ErrorCode::invalid_input);
        CHECK(m.error().message.find(text) != std::string::npos);
    };
    const auto flat = [](double z) { return plane_through({0, 0, z}, {1, 0, z}, {0, 1, z}); };
    const std::vector<RegionSeed> two{{{0.5, 0.5, 0.2}, 0}, {{0.5, 0.5, 0.8}, 1}};
    MeshOptions3D no_s = o;
    no_s.interface_resolution = 0.0;
    refused(cut_cube({flat(0.5)}, two), no_s, "interface_resolution");
    refused(cut_cube({flat(0.05)}, {{{0.5, 0.5, 0.02}, 0}, {{0.5, 0.5, 0.8}, 1}}), o, "within 4 s of a polyhedron vertex");
    refused(cut_cube({flat(1.5)}, two), o, "does not cut");
    // Steeper than 60 degrees against the side faces x = 0 and x = 1.
    refused(cut_cube({plane_through({0.2, 0, 0}, {0.2, 1, 0}, {0.8, 0, 1})}, two), o, "meets a face at under 60 degrees");
    // Two planes crossing inside the cube (a junction line).
    refused(cut_cube({flat(0.5), plane_through({0, 0, 0.3}, {1, 0, 0.7}, {0, 1, 0.3})},
                     {{{0.5, 0.5, 0.1}, 0}, {{0.5, 0.5, 0.9}, 1}, {{0.9, 0.5, 0.45}, 2}}),
            o, "");
    refused(cut_cube({flat(0.5), flat(0.55)}, {{{0.5, 0.5, 0.2}, 0}, {{0.5, 0.5, 0.52}, 1}, {{0.5, 0.5, 0.8}, 2}}), o,
            "within 4 s");
    refused(cut_cube({flat(0.5)}, {{{0.5, 0.5, 0.2}, 0}}), o, "every slab");
    refused(cut_cube({flat(0.5)}, {{{0.5, 0.5, 0.2}, 0}, {{0.5, 0.5, 0.3}, 1}}), o, "two region seeds in one slab");
    refused(cut_cube({flat(0.5)}, {{{0.5, 0.5, 0.2}, 1}, {{0.5, 0.5, 0.8}, 1}}), o, "must differ");
    refused(cut_cube({}, two), o, "only for domains with interfaces");
    Geometry3D patched = cut_cube({plane_through({0, 0, 0.5}, {0, 1, 0.5}, {1, 0, 0.5})}, two);
    patched.domain.faces[0].patches = {{"lo", {{0, 0, 0}, {0, 1, 0}, {0, 1, 0.5}, {0, 0, 0.5}}},
                                       {"hi", {{0, 0, 0.5}, {0, 1, 0.5}, {0, 1, 1}, {0, 0, 1}}}};
    refused(patched, o, "carries patches");
}

TEST_CASE("generator 3D: a linear potential is exact on the generated dual (linear-field gate)", "[mesh][generate3d]") {
    // Laplace with V = 0 on x = 0, V = 1 on x = 1.3e-4 and no flux elsewhere: the box method gives
    // the linear potential exactly when every dual face is orthogonal to its edge and every control
    // volume closes (P5), whatever the tetrahedra. Solved by conjugate gradients on the free nodes.
    const double L = 1.3e-4;
    Geometry3D g = box(L, 5e-5, 4e-5);
    MeshOptions3D o;
    o.max_size = 1.1e-5;  // not commensurate with the box: genuinely unstructured
    const auto m = generate_mesh(g, o);
    REQUIRE(m.has_value());
    const Mesh& mesh = m->mesh;
    const std::size_t n = mesh.node_count();
    std::vector<int> fixed(n, -1);
    for (const NodeId v : mesh.find_boundary("x_min")->nodes) fixed[static_cast<std::size_t>(v)] = 0;
    for (const NodeId v : mesh.find_boundary("x_max")->nodes) fixed[static_cast<std::size_t>(v)] = 1;
    std::vector<double> V(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) V[i] = fixed[i] == 1 ? 1.0 : 0.0;
    // A V on the free nodes: sum over edges of w (V_i - V_j), w = coupling / length.
    const auto apply = [&](const std::vector<double>& x) {
        std::vector<double> y(n, 0.0);
        for (const Edge& e : mesh.edges()) {
            const auto a = static_cast<std::size_t>(e.first);
            const auto b = static_cast<std::size_t>(e.second);
            const double w = e.coupling_area / e.length;
            y[a] += w * (x[a] - x[b]);
            y[b] += w * (x[b] - x[a]);
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (fixed[i] >= 0) y[i] = 0.0;
        }
        return y;
    };
    const auto dot = [&](const std::vector<double>& a, const std::vector<double>& b) {
        double s = 0.0;
        for (std::size_t i = 0; i < n; ++i) s += a[i] * b[i];
        return s;
    };
    std::vector<double> r = apply(V);
    for (double& x : r) x = -x;
    std::vector<double> p = r;
    double rr = dot(r, r);
    const double r0 = std::sqrt(rr);
    for (int it = 0; it < 10000 && std::sqrt(rr) > 1e-15 * r0; ++it) {
        const std::vector<double> Ap = apply(p);
        const double alpha = rr / dot(p, Ap);
        for (std::size_t i = 0; i < n; ++i) {
            V[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }
        const double next = dot(r, r);
        for (std::size_t i = 0; i < n; ++i) p[i] = r[i] + next / rr * p[i];
        rr = next;
    }
    double worst = 0.0;
    for (std::size_t i = 0; i < n; ++i) worst = std::max(worst, std::abs(V[i] - mesh.points()[i][0] / L));
    INFO(n << " nodes; worst deviation " << worst);
    CHECK(worst <= 1e-12);
}
