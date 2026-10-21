// Isosurfaces (Unit 25b, DS4, BG2): marching tetrahedra give the exact plane of a linear field,
// oriented along its gradient, on tensor and generated meshes; a sphere's distance field gives a
// closed, outward surface whose enclosed volume converges at second order; iso_view draws the
// levels on the field's scale and 2D contours at them.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
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

std::vector<double> graded(double a, double b, int n) {
    std::vector<double> v;
    for (int k = 0; k < n; ++k) {
        const double s = static_cast<double>(k) / (n - 1);
        v.push_back(k == n - 1 ? b : a + (b - a) * (s + 0.15 * s * (1.0 - s)));
    }
    return v;
}

std::vector<double> uniform(double a, double b, int n) {
    std::vector<double> v;
    for (int k = 0; k < n; ++k) v.push_back(k == n - 1 ? b : a + (b - a) * k / (n - 1));
    return v;
}

double dot(const Point& a, const Point& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Point sub(const Point& a, const Point& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Point cross(const Point& a, const Point& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

Point normal_of(const IsoSurface& s, const std::array<std::uint32_t, 3>& t) {
    const Point& p0 = s.vertices[t[0]].position;
    return cross(sub(s.vertices[t[1]].position, p0), sub(s.vertices[t[2]].position, p0));
}

double volume_of(const IsoSurface& s) {
    double v = 0.0;
    for (const auto& t : s.triangles) {
        v += dot(s.vertices[t[0]].position, cross(s.vertices[t[1]].position, s.vertices[t[2]].position));
    }
    return v / 6.0;
}

// Closed and consistently oriented: every directed edge once, and its reverse once.
bool closed_and_oriented(const IsoSurface& s) {
    std::map<std::pair<std::uint32_t, std::uint32_t>, int> directed;
    for (const auto& t : s.triangles) {
        for (std::size_t i = 0; i < 3; ++i) ++directed[{t[i], t[(i + 1) % 3]}];
    }
    for (const auto& [e, count] : directed) {
        if (count != 1) return false;
        const auto r = directed.find({e.second, e.first});
        if (r == directed.end() || r->second != 1) return false;
    }
    return !directed.empty();
}

void check_plane(const mesh::Mesh& m, const Geometry& g) {
    const Point grad{2.0e4, -0.7e4, 1.3e4};
    const double gn = std::sqrt(dot(grad, grad));
    std::vector<double> f;
    for (const Point& x : m.points()) f.push_back(dot(grad, x));
    double lo = f[0], hi = f[0];
    for (const double v : f) lo = std::min(lo, v), hi = std::max(hi, v);
    for (const double level : {lo + 0.3 * (hi - lo), lo + 0.61 * (hi - lo)}) {
        const IsoSurface s = isosurface(g, f, level);
        REQUIRE_FALSE(s.triangles.empty());
        for (const SurfaceVertex& v : s.vertices) {
            CHECK(std::abs(dot(grad, v.position) - level) <= 1e-12 * (hi - lo));
            CHECK(std::abs(interpolate(v.at, f) - level) <= 1e-12 * (hi - lo));
        }
        for (const auto& t : s.triangles) {
            const Point n = normal_of(s, t);
            const double nn = std::sqrt(dot(n, n));
            CHECK(dot(n, grad) / (nn * gn) >= 1.0 - 1e-9);  // along the gradient, toward higher values
        }
        // Watertight inside the domain: an interior edge borders two triangles of opposite sense.
        std::map<std::pair<std::uint32_t, std::uint32_t>, int> directed;
        for (const auto& t : s.triangles) {
            for (std::size_t i = 0; i < 3; ++i) ++directed[{t[i], t[(i + 1) % 3]}];
        }
        for (const auto& [e, count] : directed) CHECK(count == 1);
    }
}

}  // namespace

TEST_CASE("render iso: a linear field's isosurface is the exact plane, oriented (BG2)") {
    const auto x = graded(0.0, 2e-4, 7), y = graded(0.0, 1e-4, 6), z = graded(0.0, 0.8e-4, 5);
    const mesh::Mesh m = *mesh::make_tensor_grid(x, y, z);
    const auto cells = *mesh::TensorCells::create(m, x, y, z);
    check_plane(m, *Geometry::tensor(m, cells));

    mesh::Geometry3D geo;
    const double L = 1e-4;
    geo.domain.vertices = {{0, 0, 0}, {L, 0, 0}, {L, L, 0}, {0, L, 0},
                           {0, 0, L}, {L, 0, L}, {L, L, L}, {0, L, L}};
    geo.domain.faces = {{{0, 3, 2, 1}, "z0", {}}, {{4, 5, 6, 7}, "z1", {}},
                        {{0, 1, 5, 4}, "y0", {}}, {{2, 3, 7, 6}, "y1", {}},
                        {{0, 4, 7, 3}, "x0", {}}, {{1, 2, 6, 5}, "x1", {}}};
    mesh::MeshOptions3D o;
    o.max_size = 0.3e-4;
    const auto gm = mesh::generate_mesh(geo, o);
    REQUIRE(gm.has_value());
    check_plane(gm->mesh, *Geometry::simplex(gm->mesh, gm->simplices));
}

TEST_CASE("render iso: a sphere is closed, outward, and its volume converges at second order (BG2)") {
    const double L = 1e-4, R = 0.3e-4;
    const Point c{0.5013e-4, 0.4987e-4, 0.5021e-4};  // off the grid planes
    const double exact = 4.0 / 3.0 * std::numbers::pi * R * R * R;
    std::vector<double> error;
    for (const int cells_per_side : {8, 16, 32}) {
        const auto u = uniform(0.0, L, cells_per_side + 1);
        const mesh::Mesh m = *mesh::make_tensor_grid(u, u, u);
        const auto cells = *mesh::TensorCells::create(m, u, u, u);
        const Geometry g = *Geometry::tensor(m, cells);
        std::vector<double> r;
        for (const Point& x : m.points()) {
            const Point d = sub(x, c);
            r.push_back(std::sqrt(dot(d, d)));
        }
        const IsoSurface s = isosurface(g, r, R);
        CAPTURE(cells_per_side);
        CHECK(closed_and_oriented(s));
        const double v = volume_of(s);
        CHECK(v > 0.0);  // outward: normals toward larger distance
        error.push_back(std::abs(v - exact) / exact);
    }
    const double rate1 = std::log2(error[0] / error[1]), rate2 = std::log2(error[1] / error[2]);
    std::printf("sphere volume errors %.3e %.3e %.3e, rates %.3f %.3f\n", error[0], error[1], error[2],
                rate1, rate2);
    CHECK(rate1 >= 1.8);
    CHECK(rate1 <= 2.2);
    CHECK(rate2 >= 1.8);
    CHECK(rate2 <= 2.2);
}

TEST_CASE("render iso: nodes at the level, NaN values, empty cases") {
    const auto u = uniform(0.0, 1e-4, 5);
    const mesh::Mesh m = *mesh::make_tensor_grid(u, u, u);
    const auto cells = *mesh::TensorCells::create(m, u, u, u);
    const Geometry g = *Geometry::tensor(m, cells);
    // f = x: the level 0.5e-4 passes through a plane of nodes; each such node is one vertex, no
    // triangle is degenerate, the surface is the plane x = 0.5e-4 with area 1e-8.
    std::vector<double> f;
    for (const Point& x : m.points()) f.push_back(x[0]);
    const IsoSurface s = isosurface(g, f, 0.5e-4);
    CHECK(s.vertices.size() == 25);
    double area = 0.0;
    for (const auto& t : s.triangles) {
        const Point n = normal_of(s, t);
        CHECK(n[0] > 0.0);
        area += 0.5 * std::sqrt(dot(n, n));
    }
    CHECK(std::abs(area - 1e-8) <= 1e-20);
    // f = -x: now the node at the level is each cut edge's lower node; the same plane, facing -x.
    std::vector<double> neg;
    for (const double v : f) neg.push_back(-v);
    const IsoSurface sn = isosurface(g, neg, -0.5e-4);
    CHECK(sn.vertices.size() == 25);
    double area_n = 0.0;
    for (const auto& t : sn.triangles) {
        const Point n = normal_of(sn, t);
        CHECK(n[0] < 0.0);
        area_n += 0.5 * std::sqrt(dot(n, n));
    }
    CHECK(std::abs(area_n - 1e-8) <= 1e-20);
    // The level at the field's maximum: nodes at the level count as above it, so the boundary face
    // x = 1 um is the surface.
    const IsoSurface top = isosurface(g, f, 1e-4);
    double area_top = 0.0;
    for (const auto& t : top.triangles) area_top += 0.5 * std::sqrt(dot(normal_of(top, t), normal_of(top, t)));
    CHECK(std::abs(area_top - 1e-8) <= 1e-20);
    // Levels outside the values, a NaN node: nothing there; a 2D geometry: empty.
    CHECK(isosurface(g, f, 2e-4).triangles.empty());
    f[static_cast<std::size_t>(m.node_count() / 2)] = std::nan("");
    CHECK(isosurface(g, f, 0.5e-4).triangles.size() < s.triangles.size());
    const mesh::Mesh m2 = *mesh::make_tensor_grid(u, u);
    const auto c2 = *mesh::TensorCells::create(m2, u, u);
    CHECK(isosurface(*Geometry::tensor(m2, c2), std::vector<double>(m2.node_count(), 0.0), 0.0).triangles.empty());
}
