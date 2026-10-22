// The clipping box (Unit 25c, DC2, CG2): a polygon's parts outside a box; the clipped surface of tensor
// and generated 3D meshes is closed and outward, encloses the device's volume less the box's, and
// carries the interpolated field on its cut faces; views draw nothing inside the box.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/render/geometry.hpp"
#include "NiTCAD/render/scene.hpp"

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

Point sub(const Point& a, const Point& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
Point cross(const Point& a, const Point& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
double dot(const Point& a, const Point& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

// The enclosed volume (divergence theorem), the vector area (zero for a closed surface) and the area.
struct Measures {
    double volume = 0.0;
    Point vector_area{};
    double area = 0.0;
};
Measures measures(const std::vector<Polygon>& surface) {
    Measures m;
    for (const Polygon& p : surface) {
        for (std::size_t i = 1; i + 1 < p.vertices.size(); ++i) {
            const Point& a = p.vertices[0].position;
            const Point& b = p.vertices[i].position;
            const Point& c = p.vertices[i + 1].position;
            m.volume += dot(a, cross(b, c)) / 6.0;
            const Point n = cross(sub(b, a), sub(c, a));
            for (std::size_t k = 0; k < 3; ++k) m.vector_area[k] += 0.5 * n[k];
            m.area += 0.5 * std::sqrt(dot(n, n));
        }
    }
    return m;
}

double overlap(const Box& a, const Box& b) {
    double v = 1.0;
    for (std::size_t k = 0; k < 3; ++k) v *= std::max(0.0, std::min(a.hi[k], b.hi[k]) - std::max(a.lo[k], b.lo[k]));
    return v;
}

void check_closed(const std::vector<Polygon>& surface, double volume, double area_scale) {
    const Measures m = measures(surface);
    CHECK(std::abs(m.volume - volume) <= 1e-10 * volume);
    for (std::size_t k = 0; k < 3; ++k) CHECK(std::abs(m.vector_area[k]) <= 1e-10 * area_scale);
}

}  // namespace

TEST_CASE("render clip: a polygon's parts outside a box (CG2)") {
    const Box box{{1.0, 1.0, -1.0}, {3.0, 3.0, 1.0}};
    // A 4 x 4 square at z = 0 overlapping the box's 2 x 2 corner: 12 left in convex parts.
    const std::vector<Point> square{{0, 0, 0}, {4, 0, 0}, {4, 4, 0}, {0, 4, 0}};
    const auto parts = outside_box(square, box);
    double area = 0.0;
    for (const auto& p : parts) {
        area += polygon_area(p);
        // Convex: every turn the same way (the square's normal +z).
        for (std::size_t i = 0; i < p.size(); ++i) {
            const Point t = cross(sub(p[(i + 1) % p.size()], p[i]), sub(p[(i + 2) % p.size()], p[(i + 1) % p.size()]));
            CHECK(t[2] >= -1e-12);
        }
        for (const Point& x : p) CHECK_FALSE((x[0] > 1.0 && x[0] < 3.0 && x[1] > 1.0 && x[1] < 3.0));
    }
    CHECK(std::abs(area - 12.0) <= 1e-12);
    // Inside: nothing; outside: itself.
    CHECK(outside_box(std::vector<Point>{{1.5, 1.5, 0}, {2.5, 1.5, 0}, {2.5, 2.5, 0}}, box).empty());
    const auto whole = outside_box(std::vector<Point>{{5, 5, 0}, {6, 5, 0}, {6, 6, 0}}, box);
    REQUIRE(whole.size() == 1);
    CHECK(polygon_area(whole[0]) == 0.5);
}

TEST_CASE("render clip: the clipped surface is closed, outward, the device less the box (CG2)") {
    const auto x = graded(0.0, 2e-4, 9), y = graded(0.0, 1e-4, 6), z = graded(0.0, 0.6e-4, 5);
    const mesh::Mesh m = *mesh::make_tensor_grid(x, y, z);
    const auto cells = *mesh::TensorCells::create(m, x, y, z);
    const Geometry g = *Geometry::tensor(m, cells);
    const Box device{{0, 0, 0}, {2e-4, 1e-4, 0.6e-4}};
    const double V = 2e-4 * 1e-4 * 0.6e-4, A = 2e-8;
    for (const Box& box : {Box{{0.5e-4, 0.3e-4, 0.2e-4}, {1.3e-4, 0.7e-4, 0.45e-4}},   // inside
                           Box{{1.0e-4, 0.5e-4, 0.3e-4}, {5e-4, 5e-4, 5e-4}},          // a corner: the quarter cut
                           Box{{-1e-4, 0.4e-4, -1e-4}, {5e-4, 0.6e-4, 5e-4}},          // a slab through: two pieces
                           Box{{1e-4, 0.3e-4, 0.2e-4}, {2e-4, 0.7e-4, 0.45e-4}},       // a face on the boundary
                           Box{{3e-4, 3e-4, 3e-4}, {4e-4, 4e-4, 4e-4}}}) {             // missing the device
        CAPTURE(box.lo[0], box.lo[1], box.hi[0]);
        const auto surface = g.clipped_surface(box);
        check_closed(surface, V - overlap(device, box), A);
    }
    // Missing the device: the boundary itself.
    CHECK(g.clipped_surface(Box{{3e-4, 3e-4, 3e-4}, {4e-4, 4e-4, 4e-4}}).size() == g.boundary_faces().size());
    // The field on the cut faces is the interpolated one (exact for a linear field).
    std::vector<double> f;
    for (const Point& p : m.points()) f.push_back(3.0 * p[0] - 2.0 * p[1] + 5.0 * p[2]);
    for (const Polygon& p : g.clipped_surface(Box{{0.5e-4, 0.3e-4, 0.2e-4}, {1.3e-4, 0.7e-4, 0.45e-4}})) {
        for (const SurfaceVertex& v : p.vertices) {
            CHECK(std::abs(interpolate(v.at, f) - (3.0 * v.position[0] - 2.0 * v.position[1] + 5.0 * v.position[2])) <= 1e-15);
        }
    }

    // A generated mesh: the same, with the box cutting its tetrahedra anywhere.
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
    const Geometry gs = *Geometry::simplex(gm->mesh, gm->simplices);
    const Box cube{{0, 0, 0}, {L, L, L}};
    for (const Box& box : {Box{{0.31e-4, 0.27e-4, 0.36e-4}, {0.66e-4, 0.71e-4, 0.58e-4}},
                           Box{{0.47e-4, 0.53e-4, 0.41e-4}, {2e-4, 2e-4, 2e-4}}}) {
        check_closed(gs.clipped_surface(box), L * L * L - overlap(cube, box), 6.0 * L * L);
    }
}

TEST_CASE("render clip: the views draw nothing inside the box (CG2, CQ4)") {
    const auto x = graded(0.0, 2e-4, 9), y = graded(0.0, 1e-4, 6), z = graded(0.0, 0.6e-4, 5);
    mesh::Mesh m = *mesh::make_tensor_grid(x, y, z);
    auto cells = *mesh::TensorCells::create(m, x, y, z);
    // Oxide below and above the silicon (y < 0.15 um, y > 0.85 um), so interfaces are drawn.
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n, 1);
    std::vector<double> donors(n, 1e16);
    for (std::size_t i = 0; i < n; ++i) {
        const double py = m.points()[i][1];
        if (py < 0.15e-4 || py > 0.85e-4) region[i] = 0;
    }
    for (std::size_t i = 0; i < n; ++i) if (region[i] == 0) donors[i] = 0.0;
    const auto silicon_nodes = [&](const char* patch) {
        std::vector<mesh::NodeId> v;
        for (const mesh::NodeId k : m.find_boundary(patch)->nodes) {
            if (region[static_cast<std::size_t>(k)] == 1) v.push_back(k);
        }
        return v;
    };
    auto left = silicon_nodes("x_min");
    auto right = silicon_nodes("x_max");
    auto created = device::Device::create({.mesh = std::move(m),
                                            .temperature_K = 300.0,
                                            .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
                                            .node_region = std::move(region),
                                            .donors = std::move(donors),
                                            .acceptors = std::vector<double>(n, 0.0),
                                            .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                                                         {"right", device::ContactKind::ohmic, std::move(right)}},
                                            .cells = std::move(cells)});
    if (!created) FAIL(created.error().message);
    const device::Device& d = *created;
    const Geometry g = *Geometry::tensor(d.mesh(), *d.cells());
    const Box box{{1.0e-4, 0.5e-4, 0.3e-4}, {5e-4, 5e-4, 5e-4}};
    FieldOptions o;
    o.view.clip = box;
    std::vector<double> f(n);
    for (std::size_t i = 0; i < n; ++i) f[i] = d.mesh().points()[i][0];
    for (const auto& s : {field_view(d, g, f, o), solid_view(d, g, o.view)}) {
        REQUIRE(s.has_value());
        bool on_cut = false;
        for (const TriangleBatch& b : s->triangles) {
            // Contact slabs (outside the device, cut by their faces' centroids) are not the surface.
            if (b.name == "electrode" || b.name == "gate oxide") continue;
            for (const SceneVertex& v : b.vertices) {
                // um; strictly inside the box is drawn nothing.
                const bool inside = v.position[0] > 1.0f + 1e-4f && v.position[1] > 0.5f + 1e-4f &&
                                    v.position[2] > 0.3f + 1e-4f;
                CHECK_FALSE(inside);
                on_cut = on_cut || (std::abs(v.position[0] - 1.0f) < 1e-5f && v.position[1] > 0.5f);
            }
        }
        CHECK(on_cut);  // the cut face x = 1 um is drawn
    }
    // A box over the whole right end takes the right contact's pad with it.
    ViewOptions v;
    v.clip = Box{{1.5e-4, -1e-4, -1e-4}, {5e-4, 5e-4, 5e-4}};
    const auto s = solid_view(d, g, v);
    REQUIRE(s.has_value());
    for (const TriangleBatch& b : s->triangles) {
        if (b.name != "electrode") continue;
        for (const SceneVertex& vx : b.vertices) CHECK(vx.position[0] < 1.5f);
    }
}
