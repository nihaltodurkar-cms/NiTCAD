// Gates for the simplex import path, make_dual_mesh (Unit 16; DECISIONS.md M1-M4).
// - Tensor equivalence: split rectangles (2D) and Kuhn-split boxes (3D) on graded axes give
//   make_tensor_grid's volumes, couplings and boundary areas, with every diagonal dropped as a zero
//   coupling. So does any Delaunay triangulation of the lattice, which is all ties: the triangulations
//   here come from a brute-force Delaunay search with the perturbed predicates, so they also test
//   that the tie-break yields a valid triangulation.
// - Unstructured and analytic cases: closure, total measures, exact values for regular simplices.
// - Every refusal, with its code and index.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/mesh/simplex_mesh.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "dual_tolerance.hpp"
#include "predicates.hpp"
#include "simplex_fixtures.hpp"

using namespace NiTCAD::mesh;
using NiTCAD::base::ErrorCode;
namespace pr = NiTCAD::mesh::predicates;
using namespace mesh_fixtures;

namespace {

// Graded axes (a 1e2 spacing ratio on x), in cm.
const std::vector<double> ax{0.0, 1e-7, 3e-7, 7e-7, 1.5e-6, 1.0e-5};
const std::vector<double> ay{0.0, 2.0e-6, 2.5e-6, 9.0e-6};
const std::vector<double> az{0.0, 5.0e-6, 6.0e-6};

SimplexMesh split_rectangles(const std::vector<double>& x, const std::vector<double>& y,
                             bool checker) {
    SimplexMesh s;
    s.dimension = 2;
    s.points = grid_points(x, y);
    const auto nx = static_cast<NodeId>(x.size());
    const auto v = [&](std::size_t i, std::size_t j) {
        return static_cast<NodeId>(i) + nx * static_cast<NodeId>(j);
    };
    for (std::size_t j = 0; j + 1 < y.size(); ++j) {
        for (std::size_t i = 0; i + 1 < x.size(); ++i) {
            const NodeId a = v(i, j), b = v(i + 1, j), c = v(i, j + 1), d = v(i + 1, j + 1);
            if (checker && (i + j) % 2 == 1) {
                s.cells.insert(s.cells.end(), {a, b, c, b, d, c});
            } else {
                s.cells.insert(s.cells.end(), {a, b, d, a, d, c});
            }
        }
    }
    tag_box_boundary(s, {x.front(), y.front(), 0}, {x.back(), y.back(), 0});
    return s;
}

// Each box split into the six Kuhn (Freudenthal) tetrahedra along the paths from its lowest to its
// highest corner; conforming across boxes.
SimplexMesh kuhn_boxes(const std::vector<double>& x, const std::vector<double>& y,
                       const std::vector<double>& z) {
    SimplexMesh s;
    s.dimension = 3;
    s.points = grid_points(x, y, z);
    const std::array<std::size_t, 3> n{x.size(), y.size(), z.size()};
    const auto v = [&](const std::array<std::size_t, 3>& at) {
        return static_cast<NodeId>(at[0] + n[0] * (at[1] + n[1] * at[2]));
    };
    for (std::size_t k = 0; k + 1 < n[2]; ++k) {
        for (std::size_t j = 0; j + 1 < n[1]; ++j) {
            for (std::size_t i = 0; i + 1 < n[0]; ++i) {
                std::array<std::size_t, 3> perm{0, 1, 2};
                do {
                    std::array<std::size_t, 3> at{i, j, k};
                    s.cells.push_back(v(at));
                    for (const std::size_t a : perm) {
                        ++at[a];
                        s.cells.push_back(v(at));
                    }
                } while (std::ranges::next_permutation(perm).found);
            }
        }
    }
    tag_box_boundary(s, {x.front(), y.front(), z.front()}, {x.back(), y.back(), z.back()});
    return s;
}

// The Delaunay triangulation of the points under the perturbed predicates, by brute force: every
// non-degenerate simplex whose perturbed circumsphere holds no other point.
void expect_same_graph(const Mesh& got, const Mesh& want, double rel) {
    REQUIRE(got.dimension() == want.dimension());
    REQUIRE(got.node_count() == want.node_count());
    for (std::size_t i = 0; i < want.node_count(); ++i) {
        REQUIRE(got.points()[i] == want.points()[i]);
        INFO("node " << i);
        CHECK(close(got.volumes()[i], want.volumes()[i], rel));
    }
    std::map<std::pair<NodeId, NodeId>, const Edge*> expected;
    for (const Edge& e : want.edges()) expected[{e.first, e.second}] = &e;
    CHECK(got.edges().size() == want.edges().size());
    for (const Edge& e : got.edges()) {
        INFO("edge " << e.first << "-" << e.second);
        const auto it = expected.find({e.first, e.second});
        REQUIRE(it != expected.end());
        CHECK(close(e.length, it->second->length, 1e-15));
        CHECK(close(e.coupling_area, it->second->coupling_area, rel));
    }
    REQUIRE(got.boundary().size() == want.boundary().size());
    for (std::size_t p = 0; p < want.boundary().size(); ++p) {
        const BoundaryPatch& g = got.boundary()[p];
        const BoundaryPatch& w = want.boundary()[p];
        INFO("patch " << w.name);
        CHECK(g.name == w.name);
        REQUIRE(g.nodes == w.nodes);
        for (std::size_t k = 0; k < w.areas.size(); ++k) CHECK(close(g.areas[k], w.areas[k], rel));
    }
}

std::size_t unique_edge_count(int dim, const SimplexMesh& s) {
    const std::size_t per = static_cast<std::size_t>(dim) + 1;
    std::vector<std::pair<NodeId, NodeId>> e;
    for (std::size_t c = 0; c < s.cells.size() / per; ++c) {
        for (std::size_t a = 0; a < per; ++a) {
            for (std::size_t b = a + 1; b < per; ++b) {
                const NodeId u = s.cells[c * per + a], v = s.cells[c * per + b];
                e.emplace_back(std::min(u, v), std::max(u, v));
            }
        }
    }
    std::ranges::sort(e);
    return static_cast<std::size_t>(std::ranges::unique(e).begin() - e.begin());
}

ErrorCode code_of(const SimplexMesh& s) {
    const auto m = make_dual_mesh(s);
    REQUIRE_FALSE(m.has_value());
    return m.error().code;
}

std::size_t index_of(const SimplexMesh& s) {
    const auto m = make_dual_mesh(s);
    REQUIRE_FALSE(m.has_value());
    REQUIRE(m.error().context.has_value());
    REQUIRE(m.error().context->index.has_value());
    return *m.error().context->index;
}

}  // namespace

TEST_CASE("split rectangles reproduce the tensor grid; diagonals are dropped", "[mesh][simplex]") {
    const auto tensor = make_tensor_grid(ax, ay);
    REQUIRE(tensor.has_value());
    for (const bool checker : {false, true}) {
        const SimplexMesh s = split_rectangles(ax, ay, checker);
        DualMeshReport report;
        const auto dual = make_dual_mesh(s, &report);
        REQUIRE(dual.has_value());
        expect_same_graph(*dual, *tensor, 1e-13);
        CHECK(report.dropped_edges == (ax.size() - 1) * (ay.size() - 1));
        CHECK(report.dropped_patch_nodes == 0);
        CHECK(report.max_closure_residual < 1e-14);
    }
}

TEST_CASE("Kuhn-split boxes reproduce the 3D tensor grid", "[mesh][simplex]") {
    const auto tensor = make_tensor_grid(ax, ay, az);
    REQUIRE(tensor.has_value());
    const SimplexMesh s = kuhn_boxes(ax, ay, az);
    DualMeshReport report;
    const auto dual = make_dual_mesh(s, &report);
    REQUIRE(dual.has_value());
    expect_same_graph(*dual, *tensor, 1e-13);
    CHECK(report.dropped_edges == unique_edge_count(3, s) - tensor->edges().size());
    CHECK(report.max_closure_residual < 1e-13);
}

TEST_CASE("a box in five tetrahedra (all cospherical) reproduces the tensor box", "[mesh][simplex]") {
    const std::vector<double> x{0.0, 1e-6}, y{0.0, 2e-6}, z{0.0, 3e-6};
    const auto tensor = make_tensor_grid(x, y, z);
    REQUIRE(tensor.has_value());
    SimplexMesh s;
    s.dimension = 3;
    s.points = grid_points(x, y, z);  // id = i + 2 j + 4 k
    s.cells = {0, 3, 5, 6,  1, 0, 3, 5,  2, 0, 3, 6,  4, 0, 5, 6,  7, 3, 5, 6};
    tag_box_boundary(s, {0, 0, 0}, {1e-6, 2e-6, 3e-6});
    DualMeshReport report;
    const auto dual = make_dual_mesh(s, &report);
    REQUIRE(dual.has_value());
    expect_same_graph(*dual, *tensor, 1e-13);
    CHECK(report.dropped_edges == 6);  // the six face diagonals
}

TEST_CASE("the perturbed Delaunay triangulation of a lattice is valid and tensor-equivalent",
          "[mesh][simplex][predicates]") {
    // Every lattice cell is cocircular or cospherical, so every choice is a tie the perturbation
    // decides; the result must still be a triangulation, and its dual the tensor grid.
    SECTION("2D") {
        SimplexMesh s;
        s.dimension = 2;
        s.points = grid_points(ax, ay);
        s.cells = brute_delaunay(2, s.points);
        CHECK(s.cells.size() / 3 == 2 * (ax.size() - 1) * (ay.size() - 1));
        CHECK(close(cell_measure(2, s), (ax.back() - ax.front()) * (ay.back() - ay.front()), 1e-12));
        tag_box_boundary(s, {ax.front(), ay.front(), 0}, {ax.back(), ay.back(), 0});
        const auto dual = make_dual_mesh(s);
        REQUIRE(dual.has_value());
        expect_same_graph(*dual, *make_tensor_grid(ax, ay), 1e-13);
    }
    SECTION("3D") {
        const std::vector<double> x{0.0, 1e-6, 1.5e-6}, y{0.0, 2e-6, 5e-6}, z{0.0, 1e-6, 4e-6};
        SimplexMesh s;
        s.dimension = 3;
        s.points = grid_points(x, y, z);
        s.cells = brute_delaunay(3, s.points);
        CHECK(close(cell_measure(3, s), 1.5e-6 * 5e-6 * 4e-6, 1e-12));
        tag_box_boundary(s, {0, 0, 0}, {1.5e-6, 5e-6, 4e-6});
        const auto dual = make_dual_mesh(s);
        REQUIRE(dual.has_value());
        expect_same_graph(*dual, *make_tensor_grid(x, y, z), 1e-13);
    }
}

TEST_CASE("cell orientation and order do not change the graph", "[mesh][simplex]") {
    SimplexMesh s = split_rectangles(ax, ay, true);
    const auto reference = make_dual_mesh(s);
    REQUIRE(reference.has_value());
    // Reverse every other triangle and reverse the cell order.
    for (std::size_t c = 0; c < s.cells.size() / 3; c += 2) std::swap(s.cells[3 * c], s.cells[3 * c + 1]);
    std::vector<NodeId> reversed;
    for (std::size_t c = s.cells.size() / 3; c-- > 0;) {
        reversed.insert(reversed.end(), s.cells.begin() + static_cast<std::ptrdiff_t>(3 * c),
                        s.cells.begin() + static_cast<std::ptrdiff_t>(3 * c + 3));
    }
    s.cells = reversed;
    const auto again = make_dual_mesh(s);
    REQUIRE(again.has_value());
    expect_same_graph(*again, *reference, 1e-14);
}

TEST_CASE("the same input gives a bit-identical mesh", "[mesh][simplex]") {
    const SimplexMesh s = kuhn_boxes(ax, ay, az);
    const auto a = make_dual_mesh(s);
    const auto b = make_dual_mesh(s);
    REQUIRE(a.has_value());
    REQUIRE(b.has_value());
    CHECK(std::ranges::equal(a->volumes(), b->volumes()));
    REQUIRE(a->edges().size() == b->edges().size());
    for (std::size_t e = 0; e < a->edges().size(); ++e) {
        CHECK(a->edges()[e].first == b->edges()[e].first);
        CHECK(a->edges()[e].second == b->edges()[e].second);
        CHECK(a->edges()[e].coupling_area == b->edges()[e].coupling_area);
    }
    for (std::size_t p = 0; p < a->boundary().size(); ++p) {
        CHECK(a->boundary()[p].areas == b->boundary()[p].areas);
    }
}

TEST_CASE("analytic duals: a regular hexagon and a regular tetrahedron", "[mesh][simplex]") {
    SECTION("hexagon of unit side around its centre") {
        std::vector<Point> p{{0, 0, 0}};
        for (int k = 0; k < 6; ++k) {
            const double t = k * 3.141592653589793 / 3.0;
            p.push_back({std::cos(t), std::sin(t), 0});
        }
        SimplexMesh s = two_d(p, {});
        for (NodeId k = 1; k <= 6; ++k) {
            s.cells.insert(s.cells.end(), {0, k, static_cast<NodeId>(k % 6 + 1)});
            s.facets.insert(s.facets.end(), {k, static_cast<NodeId>(k % 6 + 1)});
            s.facet_patch.push_back(0);
        }
        s.patch_names = {"rim"};
        const auto m = make_dual_mesh(s);
        REQUIRE(m.has_value());
        const double r3 = std::sqrt(3.0);
        CHECK(close(m->volumes()[0], r3 / 2.0, 1e-13));
        for (std::size_t k = 1; k <= 6; ++k) CHECK(close(m->volumes()[k], r3 / 6.0, 1e-13));
        for (const Edge& e : m->edges()) {
            CHECK(close(e.coupling_area, e.first == 0 ? 1.0 / r3 : 1.0 / (2.0 * r3), 1e-13));
        }
        for (const double a : m->boundary()[0].areas) CHECK(close(a, 1.0, 1e-13));
    }
    SECTION("regular tetrahedron of edge 2 sqrt 2") {
        SimplexMesh s;
        s.dimension = 3;
        s.points = {{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}};
        s.cells = {0, 1, 2, 3};
        s.facets = {0, 1, 2, 0, 1, 3, 0, 2, 3, 1, 2, 3};
        s.facet_patch = {0, 0, 0, 0};
        s.patch_names = {"surface"};
        const auto m = make_dual_mesh(s);
        REQUIRE(m.has_value());
        for (const double v : m->volumes()) CHECK(close(v, 2.0 / 3.0, 1e-13));
        // Dual face: (h d + h d) / 2 with h = a / (2 sqrt 3), d = a / (2 sqrt 6): a^2 / (12 sqrt 2).
        for (const Edge& e : m->edges()) CHECK(close(e.coupling_area, 8.0 / (12.0 * std::sqrt(2.0)), 1e-13));
        // Each node has a third of each of its three faces: one face area, 2 sqrt 3.
        for (const double a : m->boundary()[0].areas) CHECK(close(a, 2.0 * std::sqrt(3.0), 1e-13));
    }
}

TEST_CASE("unstructured Delaunay meshes: totals, closure, obtuse cells accepted", "[mesh][simplex]") {
    SECTION("2D: a square with uneven boundary spacing and scattered interior points") {
        std::vector<Point> p;
        for (int k = 0; k <= 5; ++k) p.push_back({k / 5.0, 0, 0});
        for (int k = 0; k <= 6; ++k) p.push_back({k / 6.0, 1, 0});
        for (int k = 1; k < 7; ++k) p.push_back({0, k / 7.0, 0});
        for (int k = 1; k < 4; ++k) p.push_back({1, k / 4.0, 0});
        SplitMix rng{7};
        for (int k = 0; k < 14; ++k) p.push_back({0.25 + 0.5 * rng.uniform(), 0.25 + 0.5 * rng.uniform(), 0});
        SimplexMesh s = two_d(p, brute_delaunay(2, p));
        REQUIRE(close(cell_measure(2, s), 1.0, 1e-12));
        int obtuse = 0;
        for (std::size_t c = 0; c < s.cells.size() / 3; ++c) {
            for (std::size_t k = 0; k < 3; ++k) {
                const Point& a = p[static_cast<std::size_t>(s.cells[3 * c + k])];
                const Point& b = p[static_cast<std::size_t>(s.cells[3 * c + (k + 1) % 3])];
                const Point& d = p[static_cast<std::size_t>(s.cells[3 * c + (k + 2) % 3])];
                if ((b[0] - a[0]) * (d[0] - a[0]) + (b[1] - a[1]) * (d[1] - a[1]) < 0) ++obtuse;
            }
        }
        CHECK(obtuse > 0);  // negative pieces occur and are accepted: only totals are checked (M2)
        tag_box_boundary(s, {0, 0, 0}, {1, 1, 0});
        DualMeshReport report;
        const auto m = make_dual_mesh(s, &report);
        REQUIRE(m.has_value());
        double volume = 0.0;
        for (const double v : m->volumes()) volume += v;
        CHECK(close(volume, 1.0, 1e-12));
        for (const BoundaryPatch& b : m->boundary()) {
            double area = 0.0;
            for (const double a : b.areas) area += a;
            CHECK(close(area, 1.0, 1e-12));
        }
        CHECK(report.max_closure_residual < 1e-13);
        CHECK(report.boundary_facets == 5 + 6 + 7 + 4);
    }
    SECTION("3D: a cube with a boundary lattice and scattered interior points") {
        std::vector<Point> p;
        for (int k = 0; k < 3; ++k) {
            for (int j = 0; j < 3; ++j) {
                for (int i = 0; i < 3; ++i) {
                    if (i == 1 && j == 1 && k == 1) continue;
                    p.push_back({0.5 * i, 0.5 * j, 0.5 * k});
                }
            }
        }
        SplitMix rng{11};
        for (int k = 0; k < 3; ++k) {
            p.push_back({0.4 + 0.2 * rng.uniform(), 0.4 + 0.2 * rng.uniform(), 0.4 + 0.2 * rng.uniform()});
        }
        SimplexMesh s;
        s.dimension = 3;
        s.points = p;
        s.cells = brute_delaunay(3, p);
        REQUIRE(close(cell_measure(3, s), 1.0, 1e-12));
        tag_box_boundary(s, {0, 0, 0}, {1, 1, 1});
        DualMeshReport report;
        const auto m = make_dual_mesh(s, &report);
        REQUIRE(m.has_value());
        double volume = 0.0;
        for (const double v : m->volumes()) volume += v;
        CHECK(close(volume, 1.0, 1e-12));
        for (const BoundaryPatch& b : m->boundary()) {
            double area = 0.0;
            for (const double a : b.areas) area += a;
            CHECK(close(area, 1.0, 1e-12));
        }
        CHECK(report.max_closure_residual < 1e-12);
        CHECK(report.boundary_facets == 6 * 8);
    }
}

TEST_CASE("a square with either diagonal: the cocircular diagonal is dropped", "[mesh][simplex]") {
    const std::vector<Point> p{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
    for (const auto& cells : {std::vector<NodeId>{0, 1, 2, 0, 2, 3}, std::vector<NodeId>{0, 1, 3, 1, 2, 3}}) {
        DualMeshReport report;
        const auto m = make_dual_mesh(two_d(p, cells), &report);
        REQUIRE(m.has_value());
        CHECK(m->edges().size() == 4);
        CHECK(report.dropped_edges == 1);
        for (const double v : m->volumes()) CHECK(close(v, 0.25, 1e-15));
        for (const Edge& e : m->edges()) CHECK(close(e.coupling_area, 0.5, 1e-15));
    }
}

TEST_CASE("malformed input is refused with invalid_input and the entry's index", "[mesh][simplex]") {
    const std::vector<Point> sq{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
    const std::vector<NodeId> two{0, 1, 2, 0, 2, 3};
    SimplexMesh s = two_d(sq, two);
    s.dimension = 1;
    CHECK(code_of(s) == ErrorCode::invalid_input);
    CHECK(code_of(two_d({}, {})) == ErrorCode::invalid_input);
    CHECK(code_of(two_d(sq, {0, 1, 2, 0})) == ErrorCode::invalid_input);  // not a multiple of 3

    SimplexMesh nan = two_d(sq, two);
    nan.points[2][0] = std::nan("");
    CHECK(code_of(nan) == ErrorCode::invalid_input);
    CHECK(index_of(nan) == 2);
    SimplexMesh off_plane = two_d(sq, two);
    off_plane.points[3][2] = 1e-9;
    CHECK(code_of(off_plane) == ErrorCode::invalid_input);
    CHECK(index_of(off_plane) == 3);

    CHECK(index_of(two_d(sq, {0, 1, 2, 0, 2, 4})) == 1);  // node out of range
    CHECK(code_of(two_d(sq, {0, 1, 2, 0, 2, 2})) == ErrorCode::invalid_input);
    CHECK(index_of(two_d(sq, {0, 1, 2, 0, 2, 2})) == 1);  // repeated node
    CHECK(code_of(two_d(sq, {0, 1, 2, 2, 0, 1})) == ErrorCode::invalid_input);
    CHECK(index_of(two_d(sq, {0, 1, 2, 2, 0, 1})) == 1);  // the same cell twice
    CHECK(code_of(two_d(sq, {0, 1, 2})) == ErrorCode::invalid_input);
    CHECK(index_of(two_d(sq, {0, 1, 2})) == 3);  // node 3 in no cell

    SimplexMesh names = two_d(sq, two);
    names.patch_names = {"a", "a"};
    CHECK(code_of(names) == ErrorCode::invalid_input);
    CHECK(index_of(names) == 1);
    SimplexMesh bad_patch = two_d(sq, two);
    bad_patch.patch_names = {"a"};
    bad_patch.facets = {0, 1};
    bad_patch.facet_patch = {1};
    CHECK(code_of(bad_patch) == ErrorCode::invalid_input);
    SimplexMesh twice = two_d(sq, two);
    twice.patch_names = {"a"};
    twice.facets = {0, 1, 1, 0};
    twice.facet_patch = {0, 0};
    CHECK(code_of(twice) == ErrorCode::invalid_input);
    CHECK(index_of(twice) == 1);
    SimplexMesh interior = two_d(sq, two);
    interior.patch_names = {"a"};
    interior.facets = {0, 1, 2, 0};  // edge 0-2 is shared by both triangles
    interior.facet_patch = {0, 0};
    CHECK(code_of(interior) == ErrorCode::invalid_input);
    CHECK(index_of(interior) == 1);
    SimplexMesh sized = two_d(sq, two);
    sized.facets = {0, 1};
    CHECK(code_of(sized) == ErrorCode::invalid_input);  // facet_patch missing
}

TEST_CASE("degenerate meshes are refused with degenerate_mesh and the cell", "[mesh][simplex]") {
    SECTION("P0: a flat triangle and a flat tetrahedron") {
        const auto m = two_d({{0, 0, 0}, {1, 0, 0}, {0.5, 0, 0}, {0, 1, 0}}, {0, 1, 3, 0, 1, 2});
        CHECK(code_of(m) == ErrorCode::degenerate_mesh);
        CHECK(index_of(m) == 1);
        SimplexMesh t;
        t.dimension = 3;
        t.points = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0.25, 0.25, 0}};
        t.cells = {0, 1, 2, 3};
        CHECK(code_of(t) == ErrorCode::degenerate_mesh);
    }
    SECTION("a facet in three cells, and two cells on one side of their facet") {
        const std::vector<Point> p{{0, 0, 0}, {1, 0, 0}, {0.5, 1, 0}, {0.5, -1, 0}, {0.5, 2, 0}};
        CHECK(code_of(two_d(p, {0, 1, 2, 0, 3, 1, 0, 1, 4})) == ErrorCode::degenerate_mesh);
        const auto fold = two_d({{0, 0, 0}, {1, 0, 0}, {0.5, 1, 0}, {0.5, 2, 0}}, {0, 1, 2, 0, 1, 3});
        CHECK(code_of(fold) == ErrorCode::degenerate_mesh);
        CHECK(index_of(fold) == 1);
    }
    SECTION("P1: not locally Delaunay") {
        // (1, -1) is inside the circumcircle of (0,0), (3,0), (1,1) (centre (1.5, -0.5), r^2 2.5).
        const auto m = two_d({{0, 0, 0}, {3, 0, 0}, {1, 1, 0}, {1, -1, 0}}, {0, 1, 2, 0, 3, 1});
        CHECK(code_of(m) == ErrorCode::degenerate_mesh);
        const auto r = make_dual_mesh(m);
        CHECK(r.error().message.find("P1") != std::string::npos);
    }
    SECTION("P2: a boundary facet that is not Gabriel") {
        // The obtuse apex (1, 0.2) is inside the diametral circle of the edge (0,0)-(2,0).
        const auto m = two_d({{0, 0, 0}, {2, 0, 0}, {1, 0.2, 0}}, {0, 1, 2});
        CHECK(code_of(m) == ErrorCode::degenerate_mesh);
        const auto r = make_dual_mesh(m);
        CHECK(r.error().message.find("P2") != std::string::npos);
        CHECK(r.error().context->value == 2.0);
        // A tetrahedron with an obtuse face on the boundary: (0.5, 0.5, 0.05) sits inside the
        // diametral sphere of the face (0,0,0), (1,0,0), (0,1,0) (centre (0.5, 0.5, 0)).
        SimplexMesh t;
        t.dimension = 3;
        t.points = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0.5, 0.5, 0.05}};
        t.cells = {0, 1, 2, 3};
        CHECK(code_of(t) == ErrorCode::degenerate_mesh);
    }
}

TEST_CASE("M2 rule: exact boundaries of keep, zero and refuse", "[mesh][simplex][tolerance]") {
    using NiTCAD::mesh::detail::accept_volume;
    using NiTCAD::mesh::detail::classify_total;
    using NiTCAD::mesh::detail::dual_tau;
    using NiTCAD::mesh::detail::TotalClass;
    constexpr double inf = std::numeric_limits<double>::infinity();
    for (const double scale : {1.0, 3.0, 1e-14, 1e-30, 2.5e8}) {
        const double t = dual_tau * scale;
        INFO("scale " << scale);
        CHECK(classify_total(0.0, scale) == TotalClass::zero);
        CHECK(classify_total(t, scale) == TotalClass::zero);  // both ends included
        CHECK(classify_total(std::nextafter(t, inf), scale) == TotalClass::keep);
        CHECK(classify_total(-t, scale) == TotalClass::zero);
        CHECK(classify_total(std::nextafter(-t, -inf), scale) == TotalClass::refuse);
        CHECK_FALSE(accept_volume(t, scale));  // a zero volume is refused, not dropped
        CHECK(accept_volume(std::nextafter(t, inf), scale));
    }
    CHECK(classify_total(std::nan(""), 1.0) == TotalClass::refuse);
    CHECK(classify_total(inf, 1.0) == TotalClass::refuse);
    CHECK_FALSE(accept_volume(std::nan(""), 1.0));
}

TEST_CASE("M2 on assembled totals: the threshold, under changes of scale", "[mesh][simplex][tolerance]") {
    // The unit square with (0, 1) moved to (0, 1 + eps): the diagonal (0,0)-(1,1) is Delaunay; its
    // total coupling is eps / sqrt 2 (the right angle at (1,0) contributes 0, the angle at
    // (0, 1 + eps) has cotangent eps), its scale max(L, |pieces|) = sqrt 2. It is dropped below
    // eps = 2e-12 (total 1.414e-12 = tau * scale) and kept above, at every coordinate scale.
    for (const double lambda : {1.0, 1e-7, 0x1p-40, 1e5}) {
        for (const auto& [eps, kept] : {std::pair{1e-12, false}, std::pair{4e-12, true}, std::pair{1e-9, true}}) {
            INFO("lambda " << lambda << ", eps " << eps);
            const std::vector<Point> p{{0, 0, 0}, {lambda, 0, 0}, {lambda, lambda, 0}, {0, lambda * (1 + eps), 0}};
            DualMeshReport report;
            const auto m = make_dual_mesh(two_d(p, {0, 1, 2, 0, 2, 3}), &report);
            REQUIRE(m.has_value());
            CHECK(report.dropped_edges == (kept ? 0u : 1u));
            if (kept) {
                const auto it = std::ranges::find_if(m->edges(), [](const Edge& e) { return e.first == 0 && e.second == 2; });
                REQUIRE(it != m->edges().end());
                CHECK(close(it->coupling_area, lambda * eps / std::sqrt(2.0), 1e-3));
            }
        }
    }
    // The same in 3D: the Kuhn-split tensor grid scaled by lambda drops the same diagonals.
    std::size_t reference = 0;
    for (const double lambda : {1.0, 1e-9, 1e3}) {
        std::vector<double> x = ax, y = ay, z = az;
        for (double& v : x) v *= lambda;
        for (double& v : y) v *= lambda;
        for (double& v : z) v *= lambda;
        DualMeshReport report;
        REQUIRE(make_dual_mesh(kuhn_boxes(x, y, z), &report).has_value());
        if (lambda == 1.0) reference = report.dropped_edges;
        CHECK(report.dropped_edges == reference);
    }
}

TEST_CASE("M2 on assembled totals: cancellation of large pieces, and validity before totals",
          "[mesh][simplex][tolerance]") {
    // (-3,4), (3,4), (0,5), (0,-5) lie on the circle of radius 5. The shared edge (-3,4)-(3,4) has
    // pieces -4 (the obtuse angle at (0,5)) and +4 (at (0,-5)), both exact: the total is exactly 0,
    // its scale max(L, |pieces|) = 8, and it is dropped, not refused.
    const std::vector<Point> p{{-3, 4, 0}, {3, 4, 0}, {0, 5, 0}, {0, -5, 0}};
    DualMeshReport report;
    REQUIRE(make_dual_mesh(two_d(p, {0, 1, 2, 1, 0, 3}), &report).has_value());
    CHECK(report.dropped_edges == 1);
    // With (0,-5) moved out by one ulp the quad is still Delaunay; the total is a tiny positive
    // rounding-level value, still dropped.
    std::vector<Point> out = p;
    out[3][1] = std::nextafter(-5.0, -10.0);
    REQUIRE(make_dual_mesh(two_d(out, {0, 1, 2, 1, 0, 3}), &report).has_value());
    CHECK(report.dropped_edges == 1);
    // Moved in by one ulp it is not Delaunay: refused by the exact P1 test, although its total
    // (a tiny negative within tau * scale) would only have been dropped by M2.
    std::vector<Point> in = p;
    in[3][1] = std::nextafter(-5.0, 0.0);
    const auto refused = make_dual_mesh(two_d(in, {0, 1, 2, 1, 0, 3}));
    REQUIRE_FALSE(refused.has_value());
    CHECK(refused.error().message.find("P1") != std::string::npos);
    // (0,-10): one negative piece (-4) and one positive (6 * 187/84 / 2): only the total, 3 (187/84
    // - 4/3) = 2.678571..., decides, and it is kept.
    std::vector<Point> far = p;
    far[3][1] = -10.0;
    const auto m = make_dual_mesh(two_d(far, {0, 1, 2, 1, 0, 3}), &report);
    REQUIRE(m.has_value());
    CHECK(report.dropped_edges == 0);
    const auto it = std::ranges::find_if(m->edges(), [](const Edge& e) { return e.first == 0 && e.second == 1; });
    REQUIRE(it != m->edges().end());
    CHECK(close(it->coupling_area, 3.0 * (187.0 / 84.0 - 4.0 / 3.0), 1e-14));
}

TEST_CASE("M2: the cancellation term decides when pieces dwarf the edge", "[mesh][simplex][tolerance]") {
    // (-101, 5100), (101, 5100), (0, 5101), (0, -5101) lie on the circle of radius 5101. The chord
    // (L = 202) has pieces -5100.5 and +5100.5, 25 times L. Moving (0, -5101) out by eta keeps the
    // quad Delaunay and makes the total about eta / 2. Between tau L (2.0e-10) and tau sum |pieces|
    // (1.02e-8) only the cancellation term of the scale treats it as zero: rounding in such a
    // total scales with the pieces, not with L.
    const auto edge_kept = [](double eta) {
        const std::vector<Point> p{{-101, 5100, 0}, {101, 5100, 0}, {0, 5101, 0}, {0, -(5101 + eta), 0}};
        const auto m = make_dual_mesh(two_d(p, {0, 1, 2, 1, 0, 3}));
        REQUIRE(m.has_value());
        return std::ranges::any_of(m->edges(), [](const Edge& e) { return e.first == 0 && e.second == 1; });
    };
    CHECK_FALSE(edge_kept(0.0));
    CHECK_FALSE(edge_kept(2e-9));  // total ~1e-9: above tau L, within tau sum |pieces|
    CHECK(edge_kept(2e-7));        // total ~1e-7: kept
}
