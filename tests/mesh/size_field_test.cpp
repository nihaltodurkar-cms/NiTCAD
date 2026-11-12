// Gates for size fields (Unit 17, slice 1; DECISIONS.md Unit 17 AD4, AD10; AC1 part, AC2): exact
// location and linear interpolation, local sizes, grading, axis refinement, and the size field in
// the 2D generator and the swept generator's cross-section (honoured, bit-identical when it does not
// refine, digest unchanged when absent, certificates kept).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <set>
#include <utility>
#include <vector>

#include "NiTCAD/mesh/complex.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/mesh/size_field.hpp"

using namespace NiTCAD::mesh;
using NiTCAD::base::ErrorCode;

namespace {

RegionPolygon rect(double x0, double y0, double x1, double y1, std::int32_t region) {
    return {{{x0, y0, 0}, {x1, y0, 0}, {x1, y1, 0}, {x0, y1, 0}}, region};
}

double circumradius(const Point& a, const Point& b, const Point& c) {
    const double ab = std::hypot(b[0] - a[0], b[1] - a[1]);
    const double bc = std::hypot(c[0] - b[0], c[1] - b[1]);
    const double ca = std::hypot(a[0] - c[0], a[1] - c[1]);
    const double area2 = std::abs((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]));
    return ab * bc * ca / (2.0 * area2);
}

// A size field on a generated mesh's own triangulation.
std::shared_ptr<const SizeField> field_on(const GeneratedMesh& m, std::vector<double> h) {
    auto f = SizeField::create(m.simplices.points, m.simplices.cells, std::move(h));
    REQUIRE(f.has_value());
    return *f;
}

// The unit square as two triangles, with sizes given at its corners.
std::shared_ptr<const SizeField> square(std::vector<double> h) {
    auto f = SizeField::create({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, {0, 1, 2, 0, 2, 3}, std::move(h));
    REQUIRE(f.has_value());
    return *f;
}

}  // namespace

TEST_CASE("size field: a linear size is reproduced exactly inside, nearest outside", "[size_field]") {
    // h = 1 + 2x + 3y at the corners: linear, so every interior point gives it to rounding.
    const auto f = square({1.0, 3.0, 6.0, 4.0});
    for (int i = 0; i <= 20; ++i) {
        for (int j = 0; j <= 20; ++j) {
            const double x = i / 20.0, y = j / 20.0;
            CHECK(std::abs(f->at({x, y, 0}) - (1.0 + 2.0 * x + 3.0 * y)) <= 1e-12 * 6.0);
        }
    }
    // Outside: the nearest point of the nearest triangle.
    CHECK(std::abs(f->at({2.0, 0.5, 0}) - f->at({1.0, 0.5, 0})) <= 1e-12);
    CHECK(std::abs(f->at({-1.0, -1.0, 0}) - 1.0) <= 1e-12);
    CHECK(f->min_size() == 1.0);
}

TEST_CASE("size field: refusals", "[size_field]") {
    const std::vector<Point> p{{0, 0, 0}, {1, 0, 0}, {0, 1, 0}};
    CHECK(SizeField::create({}, {0, 1, 2}, {}).error().code == ErrorCode::invalid_input);
    CHECK(SizeField::create(p, {}, {1, 1, 1}).error().code == ErrorCode::invalid_input);
    CHECK(SizeField::create(p, {0, 1, 2}, {1, 1}).error().code == ErrorCode::invalid_input);
    CHECK(SizeField::create(p, {0, 1, 3}, {1, 1, 1}).error().code == ErrorCode::invalid_input);
    CHECK(SizeField::create(p, {0, 1, 1}, {1, 1, 1}).error().code == ErrorCode::invalid_input);
    CHECK(SizeField::create(p, {0, 1, 2}, {1, 0, 1}).error().code == ErrorCode::invalid_input);
    CHECK(SizeField::create(p, {0, 1, 2}, {1, NAN, 1}).error().code == ErrorCode::invalid_input);
    CHECK(SizeField::create({{0, 0, 0}, {1, 0, 0}, {2, 0, 0}}, {0, 1, 2}, {1, 1, 1}).error().code ==
          ErrorCode::invalid_input);  // flat
    CHECK(SizeField::create({{0, 0, 0}, {1, 0, 1}, {0, 1, 0}}, {0, 1, 2}, {1, 1, 1}).error().code ==
          ErrorCode::invalid_input);  // z != 0
    // A clockwise triangle is accepted (reoriented), and the digest names the content.
    auto cw = SizeField::create(p, {0, 2, 1}, {1, 2, 3});
    REQUIRE(cw.has_value());
    CHECK(std::abs((*cw)->at({0.25, 0.25, 0}) - (1.0 + 0.25 + 0.5)) <= 1e-12);
    CHECK((*cw)->digest() != (*SizeField::create(p, {0, 2, 1}, {1, 2, 4}))->digest());
    CHECK((*cw)->digest() == (*SizeField::create(p, {0, 2, 1}, {1, 2, 3}))->digest());
}

TEST_CASE("size field: local sizes and grading", "[size_field]") {
    // A strip of unit squares 0..4 along x, two triangles each.
    SimplexMesh sm;
    sm.dimension = 2;
    for (int i = 0; i <= 4; ++i) {
        sm.points.push_back({static_cast<double>(i), 0, 0});
        sm.points.push_back({static_cast<double>(i), 1, 0});
    }
    for (NodeId i = 0; i < 4; ++i) {
        const NodeId a = 2 * i, b = 2 * i + 2, c = 2 * i + 3, d = 2 * i + 1;
        sm.cells.insert(sm.cells.end(), {a, b, c, a, c, d});
    }
    const auto h = local_sizes(sm);
    REQUIRE(h.size() == 10);
    CHECK(h[0] == std::sqrt(2.0));  // the diagonal at node 0
    CHECK(h[1] == 1.0);             // node (0, 1): its edges are unit

    // Grading from all nodes: every edge's ratio is at most g afterwards.
    std::vector<double> s(10, 8.0);
    s[0] = 0.5;
    REQUIRE(grade_sizes(sm, s, 1.5).has_value());
    for (std::size_t c = 0; c < sm.cells.size(); c += 3) {
        for (int i = 0; i < 3; ++i) {
            const double a = s[static_cast<std::size_t>(sm.cells[c + static_cast<std::size_t>(i)])];
            const double b = s[static_cast<std::size_t>(sm.cells[c + static_cast<std::size_t>((i + 1) % 3)])];
            CHECK(std::max(a, b) <= 1.5 * std::min(a, b) * (1 + 1e-15));
        }
    }
    CHECK(s[0] == 0.5);
    CHECK(s[2] == 0.75);    // one edge away
    CHECK(s[9] == 0.5 * 1.5 * 1.5 * 1.5 * 1.5);  // node (4, 1): four edges away along the diagonals
    // Seeds: a node not seeded keeps a small size without spreading it.
    std::vector<double> t(10, 8.0);
    t[0] = 0.5;  // a pre-existing small size, not a seed
    t[8] = 1.0;  // the refined node
    std::vector<char> seeds(10, 0);
    seeds[8] = 1;
    REQUIRE(grade_sizes(sm, t, 2.0, seeds).has_value());
    CHECK(t[0] == 0.5);
    CHECK(t[2] == 8.0);  // three edges from node 8 (8 -> 6 -> 4 -> 2): 1 * 2^3 = 8, not lower
    CHECK(t[6] == 2.0);
    CHECK(t[4] == 4.0);
    CHECK(!grade_sizes(sm, t, 0.9).has_value());
    CHECK(!grade_sizes(sm, t, 1.5, std::vector<char>(3, 1)).has_value());
}

TEST_CASE("size field: axis refinement, grading against new intervals, protected ranges", "[size_field]") {
    const std::vector<double> axis{0, 1, 2, 3, 4, 5, 6, 7, 8};
    std::vector<char> m(8, 0);
    m[3] = 1;  // [3, 4]
    auto r = refine_axis(axis, m, 1.5, 0.0);
    REQUIRE(r.has_value());
    CHECK(r->split == 1);
    CHECK(r->graded == 0);  // [3, 3.5, 4]: a ratio of 2 beside it is accepted (factor max(g, 2))
    CHECK(r->lines == std::vector<double>{0, 1, 2, 3, 3.5, 4, 5, 6, 7, 8});
    // A second level, [3, 3.25]: the given [2, 3] is 4 times it and is split, [1, 2] then 2 times.
    std::vector<char> m2(9, 0);
    m2[3] = 1;
    auto r2 = refine_axis(r->lines, m2, 1.5, 0.0);
    REQUIRE(r2.has_value());
    CHECK(r2->graded == 1);
    CHECK(r2->lines == std::vector<double>{0, 1, 2, 2.5, 3, 3.25, 3.5, 4, 5, 6, 7, 8});
    const auto& l = r2->lines;
    CHECK(std::ranges::includes(l, r->lines));  // every given line kept

    // Given ratios are kept: a straddle (protected) beside unit intervals is not graded.
    const std::vector<double> s{0, 1, 1.99, 2.01, 3, 4};
    auto keep = refine_axis(s, std::vector<char>(5, 0), 1.5, 0.0, {{1.99, 2.01}});
    REQUIRE(keep.has_value());
    CHECK(keep->lines == s);
    // A marked protected interval is skipped and counted.
    std::vector<char> ms(5, 0);
    ms[2] = 1;
    auto skip = refine_axis(s, ms, 1.5, 0.0, {{1.99, 2.01}});
    REQUIRE(skip.has_value());
    CHECK(skip->lines == s);
    CHECK(skip->protected_skips == 1);
    // min_interval: a split leaving halves below it is not made.
    std::vector<char> all(8, 1);
    auto floor = refine_axis(axis, all, 1.5, 0.6);
    REQUIRE(floor.has_value());
    CHECK(floor->lines == axis);
    auto ok = refine_axis(axis, all, 1.5, 0.5);
    REQUIRE(ok.has_value());
    CHECK(ok->lines.size() == 17);
    // Refusals.
    CHECK(!refine_axis({0.0}, {}, 1.5, 0.0).has_value());
    CHECK(!refine_axis({0.0, 0.0}, {0}, 1.5, 0.0).has_value());
    CHECK(!refine_axis(axis, {0}, 1.5, 0.0).has_value());
    CHECK(!refine_axis(axis, m, 0.5, 0.0).has_value());
    CHECK(!refine_axis(axis, m, 1.5, -1.0).has_value());
    // Deterministic.
    CHECK(refine_axis(axis, m, 1.5, 0.0)->lines == r->lines);
}

TEST_CASE("size field: the 2D generator honours a field, and one that does not refine is bit-identical",
          "[size_field]") {
    Geometry2D g;
    g.polygons = {rect(0, 0, 2, 1, 0), rect(2, 0, 4, 1, 1)};
    g.patches = {{"left", {0, 0, 0}, {0, 1, 0}}, {"right", {4, 0, 0}, {4, 1, 0}}};
    MeshOptions2D o;
    o.max_size = 0.25;
    o.interface_resolution = 0.02;
    auto base = generate_mesh(g, o);
    REQUIRE(base.has_value());
    REQUIRE(base->certificate.has_value());

    // A constant field equal to h: the same mesh (AC2).
    MeshOptions2D same = o;
    same.size_field = field_on(*base, std::vector<double>(base->simplices.points.size(), 0.25));
    auto again = generate_mesh(g, same);
    REQUIRE(again.has_value());
    CHECK(graph_digest(again->mesh) == graph_digest(base->mesh));
    CHECK(again->simplices.points == base->simplices.points);

    // Refine the left half: h / 4 for x < 1, graded (g = 1.5) from those nodes.
    std::vector<double> h = local_sizes(base->simplices);
    std::vector<char> seeds(h.size(), 0);
    for (std::size_t i = 0; i < h.size(); ++i) {
        if (base->simplices.points[i][0] < 1.0) {
            h[i] = 0.0625;
            seeds[i] = 1;
        }
    }
    REQUIRE(grade_sizes(base->simplices, h, 1.5, seeds).has_value());
    // The field's grading is at most g along every background edge from a seed outward.
    MeshOptions2D fine = o;
    fine.size_field = field_on(*base, h);
    auto refined = generate_mesh(g, fine);
    REQUIRE(refined.has_value());
    REQUIRE(refined->certificate.has_value());  // P4 kept
    CHECK(refined->report.nodes > base->report.nodes);
    // Every triangle the generator judges meets its own bound against the field at its centroid.
    const auto& p = refined->simplices.points;
    const auto& c = refined->simplices.cells;
    std::size_t checked = 0, within = 0;
    for (std::size_t t = 0; t < c.size(); t += 3) {
        const Point& a = p[static_cast<std::size_t>(c[t])];
        const Point& b = p[static_cast<std::size_t>(c[t + 1])];
        const Point& d = p[static_cast<std::size_t>(c[t + 2])];
        const Point centroid{(a[0] + b[0] + d[0]) / 3, (a[1] + b[1] + d[1]) / 3, 0};
        const double bound = std::min(0.25, fine.size_field->at(centroid)) / std::sqrt(3.0);
        ++checked;
        within += circumradius(a, b, d) <= bound * (1 + 1e-12) ? 1 : 0;
    }
    // Exempt triangles (interface pairs, their bands) may exceed it; the rest may not.
    CHECK(checked - within <= refined->report.exempt_triangles);
    // Nodes concentrate where the field is small.
    std::size_t left = 0, right = 0;
    for (const Point& q : p) {
        if (q[0] < 1.0) ++left;
        if (q[0] > 3.0) ++right;
    }
    CHECK(left > 4 * right);
    // Deterministic.
    auto twice = generate_mesh(g, fine);
    REQUIRE(twice.has_value());
    CHECK(graph_digest(twice->mesh) == graph_digest(refined->mesh));
}

TEST_CASE("size field: the swept cross-section takes the field; the digest is unchanged without one",
          "[size_field]") {
    ExtrudedGeometry x;
    x.axis = 2;
    x.stations = {0.0, 3.0};
    x.cells = {{"r", {{0, 0, 0}, {2, 0, 0}, {2, 1, 0}, {0, 1, 0}}}};
    x.regions = {{0}};
    SweptOptions o;
    o.section.max_size = 0.25;
    o.lines = {0.0, 1.0, 2.0, 3.0};
    auto base = generate_swept(x, o);
    REQUIRE(base.has_value());
    const std::uint64_t d0 = generation_digest(x, o);
    CHECK(base->generation_digest == d0);

    // The cross-section mesh, to carry a field: the 2D generator on the same polygon.
    Geometry2D g2;
    g2.polygons = {rect(0, 0, 2, 1, 0)};
    MeshOptions2D o2 = o.section;
    auto m2 = generate_mesh(g2, o2);
    REQUIRE(m2.has_value());
    SweptOptions same = o;
    same.section.size_field = field_on(*m2, std::vector<double>(m2->simplices.points.size(), 0.25));
    auto again = generate_swept(x, same);
    REQUIRE(again.has_value());
    CHECK(graph_digest(again->mesh) == graph_digest(base->mesh));
    CHECK(again->generation_digest != d0);  // the field is named in the digest

    std::vector<double> h(m2->simplices.points.size(), 0.25);
    for (std::size_t i = 0; i < h.size(); ++i) {
        if (m2->simplices.points[i][0] < 0.5) h[i] = 0.05;
    }
    SweptOptions fine = o;
    fine.section.size_field = field_on(*m2, h);
    auto refined = generate_swept(x, fine);
    REQUIRE(refined.has_value());
    CHECK(refined->report.nodes > base->report.nodes);
    CHECK(refined->report.nodes % 4 == 0);  // four axis lines
}
