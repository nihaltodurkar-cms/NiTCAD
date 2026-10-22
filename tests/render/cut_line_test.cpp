// Cut lines (Unit 26b, XD5, XG4): a linear field read exactly along the line; the samples are the
// crossings with element faces (tensor and generated meshes, 2D and 3D), found here independently by
// intersecting the segment with every face; along a 1D device the cut band diagram is band_plot's;
// a line leaving the mesh has a gap there.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <utility>
#include <vector>

#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/render/cut_line.hpp"
#include "NiTCAD/render/scene.hpp"

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

std::vector<double> uniform(double a, double b, int n) {
    std::vector<double> v;
    for (int k = 0; k < n; ++k) v.push_back(k == n - 1 ? b : a + (b - a) * k / (n - 1));
    return v;
}

mesh::GeneratedMesh generated2d() {
    mesh::Geometry2D geo;
    geo.polygons = {{{{0, 0, 0}, {1e-4, 0, 0}, {1e-4, 1e-4, 0}, {0, 1e-4, 0}}, 0},
                    {{{1e-4, 0, 0}, {2e-4, 0, 0}, {2e-4, 1e-4, 0}, {1e-4, 1e-4, 0}}, 1}};
    geo.patches = {{"left", {0, 0, 0}, {0, 1e-4, 0}}, {"right", {2e-4, 0, 0}, {2e-4, 1e-4, 0}}};
    mesh::MeshOptions2D o;
    o.max_size = 0.2e-4;
    o.interface_resolution = 2e-6;
    auto gm = mesh::generate_mesh(geo, o);
    REQUIRE(gm.has_value());
    return std::move(*gm);
}

mesh::GeneratedMesh generated3d() {
    mesh::Geometry3D geo;
    const double L = 1e-4;
    geo.domain.vertices = {{0, 0, 0}, {L, 0, 0}, {L, L, 0}, {0, L, 0},
                           {0, 0, L}, {L, 0, L}, {L, L, L}, {0, L, L}};
    geo.domain.faces = {{{0, 3, 2, 1}, "z0", {}}, {{4, 5, 6, 7}, "z1", {}},
                        {{0, 1, 5, 4}, "y0", {}}, {{2, 3, 7, 6}, "y1", {}},
                        {{0, 4, 7, 3}, "x0", {}}, {{1, 2, 6, 5}, "x1", {}}};
    mesh::MeshOptions3D o;
    o.max_size = 0.3e-4;
    auto gm = mesh::generate_mesh(geo, o);
    REQUIRE(gm.has_value());
    return std::move(*gm);
}

double linear(const Point& p) { return 0.3 + 2e4 * p[0] - 3e4 * p[1] + 5e4 * p[2]; }

std::vector<double> nodal(const mesh::Mesh& m) {
    std::vector<double> v;
    for (const Point& p : m.points()) v.push_back(linear(p));
    return v;
}

Point at(const CutLine& c, double distance) {
    const double L = std::sqrt((c.b[0] - c.a[0]) * (c.b[0] - c.a[0]) + (c.b[1] - c.a[1]) * (c.b[1] - c.a[1]) +
                               (c.b[2] - c.a[2]) * (c.b[2] - c.a[2]));
    const double t = distance / L;
    return {c.a[0] + t * (c.b[0] - c.a[0]), c.a[1] + t * (c.b[1] - c.a[1]), c.a[2] + t * (c.b[2] - c.a[2])};
}

// A linear field is read exactly at every sample, and each sample's weights place it on its line.
void check_linear(const Geometry& g, const CutLine& c) {
    const auto field = nodal(g.mesh());
    const auto values = cut_values(c, field);
    REQUIRE(values.size() == c.samples.size());
    const auto pts = g.mesh().points();
    for (std::size_t i = 0; i < values.size(); ++i) {
        const Point x = at(c, c.samples[i].distance);
        CAPTURE(i, c.samples[i].distance);
        CHECK(std::abs(values[i] - linear(x)) <= 1e-12 * 10.0);
        Point y{};
        double sum = 0.0;
        for (std::size_t k = 0; k < 4; ++k) {
            const double w = c.samples[i].at.weights[k];
            CHECK(w >= 0.0);
            sum += w;
            if (w > 0.0) {
                for (std::size_t a = 0; a < 3; ++a) y[a] += w * pts[c.samples[i].at.nodes[k]][a];
            }
        }
        CHECK(std::abs(sum - 1.0) <= 1e-15);
        for (std::size_t a = 0; a < 3; ++a) CHECK(std::abs(y[a] - x[a]) <= 1e-12 * 1e-4);
    }
}

Point sub(const Point& a, const Point& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
double dot(const Point& a, const Point& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Point cross(const Point& a, const Point& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

// Where the segment a + t (b - a) meets the face (2D: an edge p q; 3D: a triangle p q r), if it
// crosses it at one point.
std::optional<double> meet(const Point& a, const Point& b, std::span<const Point> face) {
    const Point d = sub(b, a);
    if (face.size() == 2) {
        const Point e = sub(face[1], face[0]);
        const double det = d[0] * -e[1] - d[1] * -e[0];
        if (std::abs(det) < 1e-30) return std::nullopt;
        const Point r = sub(face[0], a);
        const double t = (r[0] * -e[1] - r[1] * -e[0]) / det;
        const double s = (d[0] * r[1] - d[1] * r[0]) / det;
        if (s < -1e-12 || s > 1 + 1e-12) return std::nullopt;
        return t;
    }
    // Moeller-Trumbore.
    const Point e1 = sub(face[1], face[0]), e2 = sub(face[2], face[0]);
    const Point h = cross(d, e2);
    const double det = dot(e1, h);
    if (std::abs(det) < 1e-30) return std::nullopt;
    const Point s = sub(a, face[0]);
    const double u = dot(s, h) / det;
    const Point q = cross(s, e1);
    const double v = dot(d, q) / det;
    if (u < -1e-12 || v < -1e-12 || u + v > 1 + 1e-12) return std::nullopt;
    return dot(e2, q) / det;
}

// The line parameters where the segment crosses element faces, with its ends when inside the mesh
// (here always: the lines tested start and end inside).
std::vector<double> crossings(const Geometry& g, const Point& a, const Point& b) {
    std::vector<double> ts{0.0, 1.0};
    const auto pts = g.mesh().points();
    for (std::size_t e = 0; e < g.element_count(); ++e) {
        const auto n = g.element(e);
        for (std::size_t skip = 0; skip < n.size(); ++skip) {
            std::vector<Point> face;
            for (std::size_t k = 0; k < n.size(); ++k) {
                if (k != skip) face.push_back(pts[n[k]]);
            }
            if (const auto t = meet(a, b, face); t && *t > 0.0 && *t < 1.0) ts.push_back(*t);
        }
    }
    std::ranges::sort(ts);
    std::vector<double> out;
    for (const double t : ts) {
        if (out.empty() || t - out.back() > 1e-9) out.push_back(t);
    }
    return out;
}

void check_crossings(const Geometry& g, const Point& a, const Point& b) {
    const auto c = cut_line(g, a, b);
    REQUIRE(c.has_value());
    const double L = std::sqrt(dot(sub(c->b, c->a), sub(c->b, c->a)));
    const auto expected = crossings(g, c->a, c->b);
    REQUIRE(c->samples.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        CAPTURE(i);
        CHECK(std::abs(c->samples[i].distance / L - expected[i]) <= 1e-9);
        CHECK_FALSE(c->samples[i].gap_before);
    }
    check_linear(g, *c);
}

}  // namespace

TEST_CASE("render cut line: a tensor 2D line crosses every cell face and diagonal") {
    const Tensor t = tensor(uniform(0, 2e-4, 11), uniform(0, 1e-4, 6));
    const auto g = Geometry::tensor(t.mesh, t.cells);
    REQUIRE(g.has_value());
    // Across the middle of a row of cells: 11 cell faces and 10 diagonals.
    const auto c = cut_line(*g, {0, 0.33e-4, 0}, {2e-4, 0.33e-4, 0});
    REQUIRE(c.has_value());
    CHECK(c->samples.size() == 21);
    CHECK(c->samples.front().distance == 0.0);
    CHECK(std::abs(c->samples.back().distance - 2e-4) <= 1e-18);
    for (std::size_t i = 1; i < c->samples.size(); ++i) CHECK(c->samples[i].distance > c->samples[i - 1].distance);
    check_linear(*g, *c);
    // A slanted line, and one along grid lines (on faces shared by two elements, through nodes).
    check_crossings(*g, {0.01e-4, 0.07e-4, 0}, {1.93e-4, 0.91e-4, 0});
    const auto along = cut_line(*g, {0, 0.4e-4, 0}, {2e-4, 0.4e-4, 0});
    REQUIRE(along.has_value());
    CHECK(along->samples.size() == 11);  // the nodes: diagonals meet the row there
    check_linear(*g, *along);
    const auto diagonal = cut_line(*g, {0, 0, 0}, {0.4e-4, 0.4e-4, 0});
    REQUIRE(diagonal.has_value());
    CHECK(diagonal->samples.size() == 3);  // through the nodes (0,0), (0.2,0.2), (0.4,0.4) um
    check_linear(*g, *diagonal);
}

TEST_CASE("render cut line: generated 2D and 3D meshes, samples at the element crossings") {
    const auto m2 = generated2d();
    const auto g2 = Geometry::simplex(m2.mesh, m2.simplices);
    REQUIRE(g2.has_value());
    check_crossings(*g2, {0.013e-4, 0.071e-4, 0}, {1.987e-4, 0.913e-4, 0});
    check_crossings(*g2, {0.5e-4, 0.97e-4, 0}, {1.7e-4, 0.11e-4, 0});

    const auto m3 = generated3d();
    const auto g3 = Geometry::simplex(m3.mesh, m3.simplices);
    REQUIRE(g3.has_value());
    check_crossings(*g3, {0.011e-4, 0.13e-4, 0.07e-4}, {0.97e-4, 0.89e-4, 0.93e-4});
}

TEST_CASE("render cut line: tensor 3D, through the Kuhn tetrahedra") {
    const Tensor t = tensor(uniform(0, 1e-4, 5), uniform(0, 1e-4, 4), uniform(0, 0.6e-4, 4));
    const auto g = Geometry::tensor(t.mesh, t.cells);
    REQUIRE(g.has_value());
    check_crossings(*g, {0.03e-4, 0.11e-4, 0.05e-4}, {0.95e-4, 0.87e-4, 0.53e-4});
    // Along a grid line: the nodes on it.
    const auto c = cut_line(*g, {0, 1.0e-4 / 3, 0.2e-4}, {1e-4, 1.0e-4 / 3, 0.2e-4});
    REQUIRE(c.has_value());
    CHECK(c->samples.size() == 5);
    check_linear(*g, *c);
}

TEST_CASE("render cut line: along a 1D device, band_plot's points") {
    const Tensor t = tensor(uniform(0, 1e-4, 11));
    const auto g = Geometry::tensor(t.mesh, t.cells);
    REQUIRE(g.has_value());
    const auto c = cut_line(*g, {0, 0, 0}, {1e-4, 0, 0});
    REQUIRE(c.has_value());
    REQUIRE(c->samples.size() == 11);
    results::BandDiagram b;
    for (std::size_t i = 0; i < 11; ++i) {
        const double x = static_cast<double>(i);
        b.conduction_eV.push_back(0.56 - 0.01 * x * x);
        b.valence_eV.push_back(-0.56 - 0.01 * x * x);
        b.electron_fermi_eV.push_back(i < 4 ? std::nan("") : 0.1 * x);  // NaN as on an insulator
        b.hole_fermi_eV.push_back(-0.1 * x);
    }
    // Every sample on a node: the values are the nodes' own.
    const auto v = cut_values(*c, b.conduction_eV);
    for (std::size_t i = 0; i < 11; ++i) CHECK(v[i] == b.conduction_eV[i]);
    const auto expected = band_plot(t.mesh, b);
    const auto cut = cut_band_plot(t.mesh, *c, b);
    REQUIRE(expected.has_value());
    REQUIRE(cut.has_value());
    CHECK(cut->title == expected->title);
    CHECK(cut->x.label == "distance (\xC2\xB5m)");
    CHECK(cut->y.label == expected->y.label);
    CHECK(cut->y.lo == expected->y.lo);
    CHECK(cut->y.hi == expected->y.hi);
    REQUIRE(cut->series.size() == 4);
    for (std::size_t s = 0; s < 4; ++s) {
        CHECK(cut->series[s].name == expected->series[s].name);
        CHECK(cut->series[s].color.r == expected->series[s].color.r);
        CHECK(cut->series[s].color.g == expected->series[s].color.g);
        CHECK(cut->series[s].color.b == expected->series[s].color.b);
        REQUIRE(cut->series[s].x.size() == 11);
        for (std::size_t i = 0; i < 11; ++i) {
            CHECK(std::abs(cut->series[s].x[i] - expected->series[s].x[i]) <= 1e-12);
            const double y0 = expected->series[s].y[i], y1 = cut->series[s].y[i];
            CHECK((std::isnan(y0) ? std::isnan(y1) : y1 == y0));
        }
    }
}

TEST_CASE("render cut line: a line leaving the mesh has a gap there") {
    // A U: two arms joined at the bottom; a line across the arms runs outside between them, and
    // starts and ends outside.
    mesh::Geometry2D geo;
    geo.polygons = {{{{0, 0, 0}, {3e-4, 0, 0}, {3e-4, 2e-4, 0}, {2e-4, 2e-4, 0}, {2e-4, 1e-4, 0},
                      {1e-4, 1e-4, 0}, {1e-4, 2e-4, 0}, {0, 2e-4, 0}},
                     0}};
    mesh::MeshOptions2D o;
    o.max_size = 0.3e-4;
    const auto gm = mesh::generate_mesh(geo, o);
    REQUIRE(gm.has_value());
    const auto g = Geometry::simplex(gm->mesh, gm->simplices);
    REQUIRE(g.has_value());
    const auto c = cut_line(*g, {-0.5e-4, 1.5e-4, 0}, {3.5e-4, 1.5e-4, 0});
    REQUIRE(c.has_value());
    REQUIRE(c->samples.size() >= 4);
    CHECK(std::abs(c->samples.front().distance - 0.5e-4) <= 1e-15);  // enters at x = 0
    CHECK_FALSE(c->samples.front().gap_before);
    CHECK(std::abs(c->samples.back().distance - 3.5e-4) <= 1e-15);  // leaves at x = 3 um
    std::size_t gaps = 0;
    for (std::size_t i = 0; i < c->samples.size(); ++i) {
        if (!c->samples[i].gap_before) continue;
        ++gaps;
        CHECK(std::abs(c->samples[i - 1].distance - 1.5e-4) <= 1e-15);  // left at x = 1 um
        CHECK(std::abs(c->samples[i].distance - 2.5e-4) <= 1e-15);      // back at x = 2 um
    }
    CHECK(gaps == 1);
    check_linear(*g, *c);
    // Drawn on a scene without the gap: a segment between each two samples inside, black under white, on top.
    FieldScene scene;
    add_cut_line(scene, *c);
    REQUIRE(scene.lines.size() == 2);
    CHECK(scene.lines[0].points.size() == 2 * (c->samples.size() - 2));
    CHECK(scene.lines[1].points == scene.lines[0].points);
    CHECK(scene.lines[0].width_px > scene.lines[1].width_px);
    CHECK((scene.lines[0].on_top && scene.lines[1].on_top));
    // The plot breaks there: a NaN point between the arms.
    const auto field = nodal(gm->mesh);
    const auto p = cut_plot(gm->mesh, *c, field, "f");
    REQUIRE(p.has_value());
    REQUIRE(p->series.size() == 1);
    const auto& s = p->series[0];
    CHECK(s.x.size() == c->samples.size() + 1);
    CHECK(std::ranges::count_if(s.y, [](double y) { return std::isnan(y); }) == 1);
    CHECK(s.x.front() == c->samples.front().distance * 1e4);
    CHECK(p->x.label == "distance (\xC2\xB5m)");
    CHECK(p->y.label == "f");
    // Wholly outside: no samples, and nothing to plot.
    const auto miss = cut_line(*g, {1.2e-4, 1.5e-4, 0}, {1.8e-4, 1.9e-4, 0});
    REQUIRE(miss.has_value());
    CHECK(miss->samples.empty());
    CHECK_FALSE(cut_plot(gm->mesh, *miss, field, "f").has_value());
}

TEST_CASE("render cut line: refused") {
    const Tensor t = tensor(uniform(0, 1e-4, 3), uniform(0, 1e-4, 3));
    const auto g = Geometry::tensor(t.mesh, t.cells);
    REQUIRE(g.has_value());
    CHECK_FALSE(cut_line(*g, {0.5e-4, 0.5e-4, 0}, {0.5e-4, 0.5e-4, 1e-4}).has_value());  // z unused in 2D
    CHECK_FALSE(cut_line(*g, {std::nan(""), 0, 0}, {1e-4, 0, 0}).has_value());
    CHECK_FALSE(cut_line(Geometry::graph(t.mesh), {0, 0, 0}, {1e-4, 0, 0}).has_value());
    const auto c = cut_line(*g, {0, 0, 0}, {1e-4, 0, 0});
    REQUIRE(c.has_value());
    const std::vector<double> short_field(3, 0.0);
    CHECK_FALSE(cut_plot(t.mesh, *c, short_field, "f").has_value());
    results::BandDiagram b;
    CHECK_FALSE(cut_band_plot(t.mesh, *c, b).has_value());
}
