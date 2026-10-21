// The views' geometry (Unit 25a, VG2): elements of tensor and generated meshes cover the domain
// exactly once; 3D boundary surfaces are closed and outward; Kuhn tetrahedra fill each cell; slices
// match analytic sections; control-volume pieces of a tensor mesh are exactly its control volumes;
// interfaces lie halfway between nodes; contours of a linear field lie on the level.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <numbers>
#include <utility>
#include <vector>

#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/render/geometry.hpp"

using namespace NiTCAD;
using namespace NiTCAD::render;

namespace {

struct Tensor {
    mesh::Mesh mesh;
    mesh::TensorCells cells;
};

Tensor tensor(std::vector<double> x, std::vector<double> y = {}, std::vector<double> z = {}) {
    auto m = y.empty() ? mesh::make_tensor_grid(x)
             : z.empty() ? mesh::make_tensor_grid(x, y)
                         : mesh::make_tensor_grid(x, y, z);
    REQUIRE(m.has_value());
    auto c = mesh::TensorCells::create(*m, x, y, z);
    REQUIRE(c.has_value());
    return {std::move(*m), std::move(*c)};
}

double signed_area(const Point& a, const Point& b, const Point& c) {
    return 0.5 * ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]));
}

double tet_volume(const Point& a, const Point& b, const Point& c, const Point& d) {
    const Point u{b[0] - a[0], b[1] - a[1], b[2] - a[2]}, v{c[0] - a[0], c[1] - a[1], c[2] - a[2]},
        w{d[0] - a[0], d[1] - a[1], d[2] - a[2]};
    return (u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0]) +
            u[2] * (v[0] * w[1] - v[1] * w[0])) /
           6.0;
}

// Every directed edge of the polygons appears once, with its reverse also once (a closed, consistently
// oriented surface), or for an open 2D cover: interior edges twice (once each way), boundary once.
std::map<std::pair<NodeId, NodeId>, int> directed_edges(const std::vector<std::vector<NodeId>>& polys) {
    std::map<std::pair<NodeId, NodeId>, int> e;
    for (const auto& p : polys) {
        for (std::size_t i = 0; i < p.size(); ++i) ++e[{p[i], p[(i + 1) % p.size()]}];
    }
    return e;
}

void check_cover_2d(const Geometry& g, double domain_area) {
    double area = 0.0;
    std::vector<std::vector<NodeId>> tris;
    const auto pts = g.mesh().points();
    for (std::size_t e = 0; e < g.element_count(); ++e) {
        const auto v = g.element(e);
        const double a = signed_area(pts[v[0]], pts[v[1]], pts[v[2]]);
        CHECK(a > 0.0);
        area += a;
        tris.push_back({v[0], v[1], v[2]});
    }
    CHECK(std::abs(area - domain_area) <= 1e-12 * domain_area);
    // No edge twice in the same direction: no overlap and consistent orientation.
    for (const auto& [edge, count] : directed_edges(tris)) {
        CAPTURE(edge.first, edge.second);
        CHECK(count == 1);
    }
}

void check_closed(const std::vector<Polygon>& faces, double volume, double area) {
    std::vector<std::vector<NodeId>> polys;
    double total = 0.0;
    for (const Polygon& f : faces) {
        std::vector<NodeId> p;
        std::vector<Point> xs;
        for (const SurfaceVertex& v : f.vertices) {
            p.push_back(v.at.nodes[0]);
            xs.push_back(v.position);
        }
        polys.push_back(p);
        total += polygon_area(xs);
    }
    const auto edges = directed_edges(polys);
    for (const auto& [edge, count] : edges) {
        CAPTURE(edge.first, edge.second);
        CHECK(count == 1);
        CHECK(edges.count({edge.second, edge.first}) == 1);
    }
    CHECK(std::abs(enclosed_volume(faces) - volume) <= 1e-12 * volume);
    if (area > 0.0) CHECK(std::abs(total - area) <= 1e-12 * area);
}

}  // namespace

TEST_CASE("render geometry: tensor 2D triangles cover the domain once") {
    const Tensor t = tensor({0.0, 0.1, 0.35, 0.4, 1.0}, {0.0, 0.2, 0.25, 0.7});
    const auto g = Geometry::tensor(t.mesh, t.cells);
    REQUIRE(g.has_value());
    CHECK(g->kind() == GeometryKind::tensor);
    CHECK(g->element_count() == 2 * 4 * 3);
    CHECK(g->cell_count() == 12);
    check_cover_2d(*g, 1.0 * 0.7);
}

TEST_CASE("render geometry: generated 2D triangles cover the domain once") {
    mesh::Geometry2D geo;
    geo.polygons = {{{{0, 0, 0}, {1e-4, 0, 0}, {1e-4, 1e-4, 0}, {0, 1e-4, 0}}, 0},
                    {{{1e-4, 0, 0}, {2e-4, 0, 0}, {2e-4, 1e-4, 0}, {1e-4, 1e-4, 0}}, 1}};
    geo.patches = {{"left", {0, 0, 0}, {0, 1e-4, 0}}, {"right", {2e-4, 0, 0}, {2e-4, 1e-4, 0}}};
    mesh::MeshOptions2D o;
    o.max_size = 0.2e-4;
    o.interface_resolution = 2e-6;
    const auto gm = mesh::generate_mesh(geo, o);
    REQUIRE(gm.has_value());
    const auto g = Geometry::simplex(gm->mesh, gm->simplices);
    REQUIRE(g.has_value());
    check_cover_2d(*g, 2e-4 * 1e-4);

    // Owner pieces split every triangle without loss, each piece nearest its owner.
    std::vector<double> by_region(2, 0.0);
    double total = 0.0;
    for (std::size_t e = 0; e < g->element_count(); ++e) {
        for (const OwnedPolygon& p : g->owner_pieces(g->element_polygon(e))) {
            const double a = polygon_area(p.points);
            total += a;
            by_region[static_cast<std::size_t>(gm->node_region[static_cast<std::size_t>(p.owner)])] += a;
            const auto v = g->element(e);
            CHECK(std::find(v.begin(), v.end(), p.owner) != v.end());
        }
    }
    CHECK(std::abs(total - 2e-8) <= 1e-12 * 2e-8);
    // Each piece is nearest its owner: at the piece's centroid the owner weighs most.
    const auto pts = gm->mesh.points();
    for (std::size_t e = 0; e < g->element_count(); e += 7) {
        const auto v = g->element(e);
        for (const OwnedPolygon& p : g->owner_pieces(g->element_polygon(e))) {
            Point c{};
            for (const Point& x : p.points) {
                for (std::size_t a = 0; a < 3; ++a) c[a] += x[a] / static_cast<double>(p.points.size());
            }
            const Interpolant at = g->interpolant(e, c);
            double owner = -1.0, best = -1.0;
            for (std::size_t k = 0; k < 3; ++k) {
                best = std::max(best, at.weights[k]);
                if (at.nodes[k] == p.owner) owner = at.weights[k];
            }
            CHECK(owner == best);
            (void)pts;
            (void)v;
        }
    }
    // The interfaces separate the two regions: they lie within an interface resolution of x = 1 um.
    const auto f = g->interfaces(gm->node_region);
    REQUIRE(!f.empty());
    for (const InterfacePiece& p : f) {
        CHECK(p.region_a == 0);
        CHECK(p.region_b == 1);
        for (const Point& x : p.points) CHECK(std::abs(x[0] - 1e-4) <= 2e-6);
    }
    CHECK(std::abs(by_region[0] - 1e-8) <= 0.03e-8);
}

TEST_CASE("render geometry: imported clockwise simplices are reoriented") {
    // Two clockwise triangles of the unit square: the elements come out counterclockwise.
    mesh::SimplexMesh s;
    s.dimension = 2;
    s.points = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}};
    s.cells = {0, 2, 1, 0, 3, 2};
    const auto m = mesh::make_dual_mesh(s);
    REQUIRE(m.has_value());
    const auto g = Geometry::simplex(*m, s);
    REQUIRE(g.has_value());
    check_cover_2d(*g, 1.0);
}

TEST_CASE("render geometry: a mismatched simplex mesh is refused") {
    const Tensor t = tensor({0.0, 1.0}, {0.0, 1.0});
    mesh::SimplexMesh s;
    s.dimension = 2;
    s.points = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    s.cells = {0, 1, 2};
    CHECK(Geometry::simplex(t.mesh, s).error().code == base::ErrorCode::invalid_input);
    const Tensor other = tensor({0.0, 1.0, 2.0}, {0.0, 1.0});
    CHECK(Geometry::tensor(t.mesh, other.cells).error().code == base::ErrorCode::invalid_input);
}

TEST_CASE("render geometry: Kuhn tetrahedra fill each cell; the 3D boundary is closed") {
    const Tensor t = tensor({0.0, 0.3, 1.0, 1.2}, {0.0, 0.5, 0.6}, {0.0, 0.25, 0.4});
    const auto g = Geometry::tensor(t.mesh, t.cells);
    REQUIRE(g.has_value());
    CHECK(g->element_count() == 6 * 3 * 2 * 2);
    const auto pts = t.mesh.points();
    std::vector<double> cell(g->cell_count(), 0.0);
    for (std::size_t e = 0; e < g->element_count(); ++e) {
        const auto v = g->element(e);
        const double vol = tet_volume(pts[v[0]], pts[v[1]], pts[v[2]], pts[v[3]]);
        CHECK(vol > 0.0);
        cell[g->element_cell(e)] += vol;
    }
    for (std::size_t c = 0; c < g->cell_count(); ++c) {
        const auto [lo, hi] = g->cell_box(c);
        const double box = (hi[0] - lo[0]) * (hi[1] - lo[1]) * (hi[2] - lo[2]);
        CHECK(std::abs(cell[c] - box) <= 1e-14 * box);
    }
    const double V = 1.2 * 0.6 * 0.4;
    const double A = 2.0 * (1.2 * 0.6 + 1.2 * 0.4 + 0.6 * 0.4);
    check_closed(g->boundary_faces(), V, A);
}

TEST_CASE("render geometry: a generated 3D mesh's boundary is closed") {
    mesh::Geometry3D geo;
    const double L = 1e-4;
    geo.domain.vertices = {{0, 0, 0}, {L, 0, 0}, {L, L, 0}, {0, L, 0},
                           {0, 0, L}, {L, 0, L}, {L, L, L}, {0, L, L}};
    geo.domain.faces = {{{0, 3, 2, 1}, "z0", {}}, {{4, 5, 6, 7}, "z1", {}},
                        {{0, 1, 5, 4}, "y0", {}}, {{2, 3, 7, 6}, "y1", {}},
                        {{0, 4, 7, 3}, "x0", {}}, {{1, 2, 6, 5}, "x1", {}}};
    mesh::MeshOptions3D o;
    o.max_size = 0.4e-4;
    const auto gm = mesh::generate_mesh(geo, o);
    REQUIRE(gm.has_value());
    const auto g = Geometry::simplex(gm->mesh, gm->simplices);
    REQUIRE(g.has_value());
    check_closed(g->boundary_faces(), L * L * L, 6.0 * L * L);
}

TEST_CASE("render geometry: tensor control-volume pieces are the control volumes") {
    // 2D: the area of each node's pieces over all triangles is its control volume (per unit depth).
    const Tensor t = tensor({0.0, 0.1, 0.35, 0.4, 1.0}, {0.0, 0.2, 0.25, 0.7});
    const auto g = Geometry::tensor(t.mesh, t.cells);
    REQUIRE(g.has_value());
    std::vector<double> area(t.mesh.node_count(), 0.0);
    for (std::size_t e = 0; e < g->element_count(); ++e) {
        for (const OwnedPolygon& p : g->owner_pieces(g->element_polygon(e))) {
            area[static_cast<std::size_t>(p.owner)] += polygon_area(p.points);
        }
    }
    const auto volumes = t.mesh.volumes();
    for (std::size_t i = 0; i < area.size(); ++i) {
        CAPTURE(i);
        CHECK(std::abs(area[i] - volumes[i]) <= 1e-14);
    }
}

TEST_CASE("render geometry: slices match analytic sections") {
    // The unit cube on a 5-node lattice: a section x = 0.73 (area 1), one through grid nodes
    // z = 0.5 (area 1, each face once), and the hexagon x + y + z = 1.5 (area 3 sqrt 3 / 4).
    std::vector<double> a{0.0, 0.25, 0.5, 0.75, 1.0};
    const Tensor t = tensor(a, a, a);
    const auto g = Geometry::tensor(t.mesh, t.cells);
    REQUIRE(g.has_value());
    std::vector<double> f;
    for (const Point& p : t.mesh.points()) f.push_back(p[0] + 2.0 * p[1] + 3.0 * p[2]);
    const auto check = [&](const Plane& plane, double area) {
        double total = 0.0;
        for (const Polygon& p : g->slice(plane)) {
            std::vector<Point> xs;
            for (const SurfaceVertex& v : p.vertices) {
                xs.push_back(v.position);
                const double exact = v.position[0] + 2.0 * v.position[1] + 3.0 * v.position[2];
                CHECK(std::abs(interpolate(v.at, f) - exact) <= 1e-14);
                CHECK(std::abs(v.position[0] * plane.normal[0] + v.position[1] * plane.normal[1] +
                               v.position[2] * plane.normal[2] - plane.offset) <= 1e-14);
            }
            total += polygon_area(xs);
        }
        CAPTURE(plane.offset, total, area);
        CHECK(std::abs(total - area) <= 1e-13);
    };
    check({{1.0, 0.0, 0.0}, 0.73}, 1.0);
    check({{0.0, 0.0, 1.0}, 0.5}, 1.0);
    check({{0.0, 0.0, 1.0}, 0.0}, 1.0);
    check({{1.0, 1.0, 1.0}, 1.5}, 3.0 * std::sqrt(3.0) / 4.0);

    // 2D: a line y = 0.3 through a tensor rectangle gives segments of total length its width.
    const Tensor r = tensor({0.0, 0.1, 0.35, 0.4, 1.0}, {0.0, 0.2, 0.25, 0.7});
    const auto g2 = Geometry::tensor(r.mesh, r.cells);
    double length = 0.0;
    for (const Polygon& p : g2->slice({{0.0, 1.0, 0.0}, 0.3})) {
        REQUIRE(p.vertices.size() == 2);
        length += std::hypot(p.vertices[1].position[0] - p.vertices[0].position[0],
                             p.vertices[1].position[1] - p.vertices[0].position[1]);
    }
    CHECK(std::abs(length - 1.0) <= 1e-14);
}

TEST_CASE("render geometry: tensor interfaces lie halfway between the nodes") {
    // Regions split between the node columns x = 0.35 and x = 0.4: the interface is x = 0.375,
    // over the whole height.
    const Tensor t = tensor({0.0, 0.1, 0.35, 0.4, 1.0}, {0.0, 0.2, 0.25, 0.7});
    const auto g = Geometry::tensor(t.mesh, t.cells);
    std::vector<std::int32_t> region;
    for (const Point& p : t.mesh.points()) region.push_back(p[0] < 0.37 ? 0 : 1);
    double length = 0.0;
    for (const InterfacePiece& p : g->interfaces(region)) {
        REQUIRE(p.points.size() == 2);
        CHECK(p.points[0][0] == 0.375);
        CHECK(p.points[1][0] == 0.375);
        length += std::abs(p.points[1][1] - p.points[0][1]);
    }
    CHECK(std::abs(length - 0.7) <= 1e-15);

    // 3D: a plane z = mid of a layer, of area the cross-section.
    const Tensor c = tensor({0.0, 0.5, 1.0}, {0.0, 0.4, 1.0}, {0.0, 0.2, 0.6});
    const auto g3 = Geometry::tensor(c.mesh, c.cells);
    std::vector<std::int32_t> r3;
    for (const Point& p : c.mesh.points()) r3.push_back(p[2] < 0.3 ? 0 : 1);
    double area = 0.0;
    for (const InterfacePiece& p : g3->interfaces(r3)) {
        for (const Point& x : p.points) CHECK(x[2] == 0.4);
        area += polygon_area(p.points);
    }
    CHECK(std::abs(area - 1.0) <= 1e-15);
}

TEST_CASE("render geometry: contours of a linear field lie on the level") {
    const std::array<Point, 3> p{{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}}};
    const auto f = [](const Point& x) { return 2.0 * x[0] + 3.0 * x[1] - 0.25; };
    const std::array<double, 3> v{f(p[0]), f(p[1]), f(p[2])};
    for (const double level : {0.0, 0.5, 1.7}) {
        const auto s = contour_segment(p, v, level);
        REQUIRE(s.has_value());
        CHECK(std::abs(f((*s)[0]) - level) <= 1e-15);
        CHECK(std::abs(f((*s)[1]) - level) <= 1e-15);
    }
    CHECK_FALSE(contour_segment(p, v, 5.0).has_value());
    CHECK_FALSE(contour_segment(p, v, -1.0).has_value());
    // A level exactly at a vertex: the vertex counts as above, so one segment from it or none.
    const std::array<double, 3> w{0.0, 1.0, -1.0};
    const auto at = contour_segment(p, w, 0.0);
    REQUIRE(at.has_value());
}
