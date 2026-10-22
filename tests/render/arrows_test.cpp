// Arrows (Unit 25c, DC1, CG1): grid and surface sites, directions exact in a uniform field, tangential
// on no-flux faces, none where the field is undefined; the drawn arrows have the set length, are
// lifted off their surface, and carry the colormap position of |v| or the given color.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <set>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/render/scene.hpp"
#include "NiTCAD/render/streamlines.hpp"
#include "NiTCAD/render/vectors.hpp"

using namespace NiTCAD;
using namespace NiTCAD::render;

namespace {

std::vector<double> uniform(double a, double b, int n) {
    std::vector<double> v;
    for (int k = 0; k < n; ++k) v.push_back(k == n - 1 ? b : a + (b - a) * k / (n - 1));
    return v;
}

double dot(const Point& a, const Point& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

// A silicon bar, 2 um by 1 um (by 0.6 um in 3D), ohmic contacts on x_min and x_max.
device::Device bar(int D) {
    const auto x = uniform(0.0, 2e-4, 11), y = uniform(0.0, 1e-4, 6), z = uniform(0.0, 0.6e-4, 4);
    mesh::Mesh m = D == 2 ? *mesh::make_tensor_grid(x, y) : *mesh::make_tensor_grid(x, y, z);
    auto cells = D == 2 ? *mesh::TensorCells::create(m, x, y) : *mesh::TensorCells::create(m, x, y, z);
    auto left = m.find_boundary("x_min")->nodes;
    auto right = m.find_boundary("x_max")->nodes;
    const std::size_t n = m.node_count();
    auto d = device::Device::create({.mesh = std::move(m),
                                     .temperature_K = 300.0,
                                     .regions = {{"silicon", physics::silicon()}},
                                     .node_region = std::vector<device::RegionId>(n, 0),
                                     .donors = std::vector<double>(n, 1e16),
                                     .acceptors = std::vector<double>(n, 0.0),
                                     .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                                                  {"right", device::ContactKind::ohmic, std::move(right)}},
                                     .cells = std::move(cells)});
    if (!d) FAIL(d.error().message);
    return std::move(*d);
}

}  // namespace

TEST_CASE("render arrows: a 2D grid of arrows in a uniform field (CG1)") {
    const device::Device d = bar(2);
    const Geometry g = *Geometry::tensor(d.mesh(), *d.cells());
    const Locator loc(g);
    const Point v{3.0, 4.0, 0.0};
    std::vector<Point> field(d.mesh().node_count(), v);
    for (const std::size_t count : {50u, 200u, 800u}) {
        CAPTURE(count);
        double h = 0.0;
        const auto arrows = grid_arrows(loc, field, count, h);
        // The bar fills its box: about `count` arrows, on a square grid of spacing h.
        CHECK(std::abs(static_cast<double>(arrows.size()) - static_cast<double>(count)) <= 0.1 * static_cast<double>(count));
        CHECK(std::abs(h - std::sqrt(2e-8 / static_cast<double>(count))) <= 1e-12 * h);
        std::set<long long> columns;
        for (const Arrow& a : arrows) {
            CHECK(std::abs(a.direction[0] - 0.6) <= 1e-12);
            CHECK(std::abs(a.direction[1] - 0.8) <= 1e-12);
            CHECK(std::abs(a.magnitude - 5.0) <= 1e-12);
            CHECK(a.position[0] > 0.0);
            CHECK(a.position[0] < 2e-4);
            columns.insert(std::llround(a.position[0] / h * 1e6));
        }
        // Centred: the outer columns and rows are symmetric about the box's centre (1, 0.5) um.
        double x_lo = 1.0, x_hi = 0.0, y_lo = 1.0, y_hi = 0.0;
        for (const Arrow& a : arrows) {
            x_lo = std::min(x_lo, a.position[0]), x_hi = std::max(x_hi, a.position[0]);
            y_lo = std::min(y_lo, a.position[1]), y_hi = std::max(y_hi, a.position[1]);
        }
        CHECK(std::abs(0.5 * (x_lo + x_hi) - 1e-4) <= 1e-12 * 1e-4);
        CHECK(std::abs(0.5 * (y_lo + y_hi) - 0.5e-4) <= 1e-12 * 1e-4);
        // Neighbouring columns one spacing apart.
        std::vector<long long> c(columns.begin(), columns.end());
        for (std::size_t i = 1; i < c.size(); ++i) CHECK(std::llabs(c[i] - c[i - 1] - 1000000) <= 1);
    }
    // No arrow where the field is undefined (NaN) or zero.
    for (std::size_t i = 0; i < field.size(); ++i) {
        const double x = d.mesh().points()[i][0];
        if (x >= 1.4e-4) field[i] = {std::nan(""), 0.0, 0.0};
        else if (x <= 0.6e-4) field[i] = {0.0, 0.0, 0.0};
    }
    double h = 0.0;
    for (const Arrow& a : grid_arrows(loc, field, 200, h)) {
        CHECK(a.position[0] < 1.4e-4);
        CHECK(a.position[0] > 0.4e-4);  // within one cell of the zero region the field is nonzero
    }
}

TEST_CASE("render arrows: excluded nodes' control volumes get no arrow (CG1)") {
    // The upper half of the bar's nodes flagged (as an insulator's are for currents): no site whose
    // nearest node is flagged; the lower half keeps its arrows.
    const device::Device d = bar(2);
    const Geometry g = *Geometry::tensor(d.mesh(), *d.cells());
    const Locator loc(g);
    const std::vector<Point> field(d.mesh().node_count(), Point{1.0, 0.0, 0.0});
    std::vector<std::uint8_t> excluded;
    for (const Point& p : d.mesh().points()) excluded.push_back(p[1] > 0.5e-4 ? 1 : 0);
    double h = 0.0;
    const auto all = grid_arrows(loc, field, 400, h);
    const auto kept = grid_arrows(loc, field, 400, h, excluded);
    CHECK(kept.size() < all.size());
    CHECK(kept.size() > all.size() / 3);
    for (const Arrow& a : kept) {
        std::size_t nearest = 0;
        for (std::size_t k = 1; k < 3; ++k) {
            if (a.at.weights[k] > a.at.weights[nearest]) nearest = k;
        }
        CHECK(excluded[static_cast<std::size_t>(a.at.nodes[nearest])] == 0);
    }
}

TEST_CASE("render arrows: surface arrows are tangential on no-flux faces (CG1)") {
    const device::Device d = bar(3);
    const Geometry g = *Geometry::tensor(d.mesh(), *d.cells());
    const Locator loc(g);
    results::BiasPoint p;
    for (const Point& x : d.mesh().points()) p.fields.potential_V.push_back(-3e3 * x[0] - 4e3 * x[1] - 2e3 * x[2]);
    const auto field = streamline_field(d, g, *element_vectors(d, g, p, VectorField::electric), VectorField::electric);
    // The boundary except the contacts' faces (x = 0 and x = 2 um): 2 of the 6 sides, 5 x 3 cells of
    // 2 triangles each.
    const std::vector<Polygon> faces = exposed_boundary(d, g);
    CHECK(faces.size() == g.boundary_faces().size() - 2 * 5 * 3 * 2);
    for (const Polygon& f : faces) {
        bool on_contact = true;
        for (const SurfaceVertex& v : f.vertices) on_contact = on_contact && (v.position[0] == 0.0 || v.position[0] == 2e-4);
        CHECK_FALSE(on_contact);
    }
    double h = 0.0;
    const auto arrows = surface_arrows(loc, field, faces, 300, h);
    CHECK(arrows.size() >= 270);
    CHECK(arrows.size() <= 300);
    for (const Arrow& a : arrows) {
        CHECK(std::abs(dot(a.normal, a.normal) - 1.0) <= 1e-12);
        // On the y and z faces away from the contacts (whose nodes keep their full field, so the cells
        // next to them interpolate a normal part): no normal component.
        if (a.position[0] < 0.2e-4 || a.position[0] > 1.8e-4) continue;
        CHECK(std::abs(dot(a.direction, a.normal)) <= 1e-9);
    }
}

TEST_CASE("render arrows: drawn with the set length, lifted, colored (CG1)") {
    Arrow a;
    a.position = {1e-4, 0.5e-4, 0.0};
    a.direction = {0.6, 0.8, 0.0};
    a.magnitude = 1e3;
    const double length = 0.2e-4;  // cm
    // 2D: flat, one color, spanning the length along the direction, just in front of z = 0.
    FieldScene s2;
    s2.dimension = 2;
    add_arrows(s2, std::vector<Arrow>{a}, length, Rgba{0.1f, 0.2f, 0.3f, 1.0f});
    REQUIRE(s2.triangles.size() == 1);
    const TriangleBatch& b2 = s2.triangles[0];
    CHECK(b2.name == "arrows");
    CHECK_FALSE(b2.lit);
    double lo = 1e9, hi = -1e9;
    for (const SceneVertex& v : b2.vertices) {
        const double along = (v.position[0] - 1.0) * 0.6 + (v.position[1] - 0.5) * 0.8;  // um
        lo = std::min(lo, along);
        hi = std::max(hi, along);
        CHECK(v.position[2] > 0.0f);
        CHECK(v.color.r == 0.1f);
    }
    CHECK(std::abs(lo + 0.1) <= 1e-5);
    CHECK(std::abs(hi - 0.1) <= 1e-5);
    // 3D: lit, lifted off the surface along its normal, colored by |v|'s position.
    FieldScene s3;
    s3.dimension = 3;
    a.normal = {0.0, 0.0, 1.0};
    const Colorbar bar{Colormap::turbo, *ColorScale::create(ScaleKind::log10, 1e0, 1e6), "|E|"};
    add_arrows(s3, std::vector<Arrow>{a}, length, bar);
    REQUIRE(s3.triangles.size() == 1);
    const TriangleBatch& b3 = s3.triangles[0];
    CHECK(b3.lit);
    CHECK(b3.shading == Shading::field);
    double hi3 = -1e9, lo3 = 1e9;
    for (const SceneVertex& v : b3.vertices) {
        CHECK(v.value == 0.5f);
        CHECK(v.position[2] > 0.0f);  // above the surface z = 0
        const double along = (v.position[0] - 1.0) * 0.6 + (v.position[1] - 0.5) * 0.8;
        lo3 = std::min(lo3, along);
        hi3 = std::max(hi3, along);
    }
    CHECK(std::abs(lo3 + 0.1) <= 1e-5);
    CHECK(std::abs(hi3 - 0.1) <= 1e-5);
}
