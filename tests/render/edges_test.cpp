// Unit 26f (KG2, KG3): feature edges, set_look and the enhanced look's lights.
#include "NiTCAD/render/edges.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <random>
#include <string>
#include <vector>

#include "NiTCAD/render/image.hpp"

using namespace NiTCAD::render;

namespace {

void triangle(TriangleBatch& b, const Float3& p, const Float3& q, const Float3& r) {
    for (const Float3& x : {p, q, r}) b.vertices.push_back({x, {0.0f, 0.0f, 1.0f}, 0.0f, {0.5f, 0.5f, 0.5f, 1.0f}});
}

// A quad p q r s (counterclockwise) as two triangles along p-r; `flip` winds the second the other way.
void quad(TriangleBatch& b, const Float3& p, const Float3& q, const Float3& r, const Float3& s, bool flip = false) {
    triangle(b, p, q, r);
    if (flip) {
        triangle(b, p, s, r);
    } else {
        triangle(b, p, r, s);
    }
}

TriangleBatch cube(std::string name = "silicon") {
    TriangleBatch b{std::move(name), Shading::flat, true, false, {}};
    const auto c = [](int k) { return Float3{float(k & 1), float((k >> 1) & 1), float((k >> 2) & 1)}; };
    quad(b, c(0), c(2), c(3), c(1));  // z = 0
    quad(b, c(4), c(5), c(7), c(6));  // z = 1
    quad(b, c(0), c(1), c(5), c(4));  // y = 0
    quad(b, c(2), c(6), c(7), c(3));  // y = 1
    quad(b, c(0), c(4), c(6), c(2));  // x = 0
    quad(b, c(1), c(3), c(7), c(5));  // x = 1
    return b;
}

// Two unit squares meeting along the y axis at z = 0, the second turned up by `deg` about it.
TriangleBatch fold(double deg) {
    TriangleBatch b{"silicon", Shading::flat, true, false, {}};
    const double t = deg * std::numbers::pi / 180.0;
    const Float3 o{0, 0, 0}, top{0, 1, 0}, a0{-1, 0, 0}, a1{-1, 1, 0};
    const Float3 b0{float(std::cos(t)), 0.0f, float(std::sin(t))};
    const Float3 b1{b0[0], 1.0f, b0[2]};
    quad(b, a0, o, top, a1);
    quad(b, o, b0, b1, top);
    return b;
}

bool has(const std::vector<Segment>& edges, const Float3& p, const Float3& q) {
    return std::ranges::any_of(edges, [&](const Segment& s) { return (s[0] == p && s[1] == q) || (s[0] == q && s[1] == p); });
}

double length3(const std::array<float, 3>& v) { return std::sqrt(double(v[0]) * v[0] + double(v[1]) * v[1] + double(v[2]) * v[2]); }

double dot3(const std::array<float, 3>& a, const std::array<double, 3>& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

}  // namespace

TEST_CASE("render edges: a closed cube gives exactly its 12 edges") {
    const auto e = feature_edges(cube());
    REQUIRE(e.size() == 12);
    // Each a unit-length axis edge: none of the faces' diagonals.
    for (const Segment& s : e) {
        const float d = std::abs(s[0][0] - s[1][0]) + std::abs(s[0][1] - s[1][1]) + std::abs(s[0][2] - s[1][2]);
        CHECK(d == 1.0f);
    }
}

TEST_CASE("render edges: a flat square gives only its boundary") {
    TriangleBatch b{"slice", Shading::field, true, false, {}};
    // 2 x 2 cells of the unit square, 8 triangles.
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            const float x0 = 0.5f * i, x1 = 0.5f * (i + 1), y0 = 0.5f * j, y1 = 0.5f * (j + 1);
            quad(b, {x0, y0, 0}, {x1, y0, 0}, {x1, y1, 0}, {x0, y1, 0});
        }
    }
    const auto e = feature_edges(b);
    REQUIRE(e.size() == 8);
    for (const Segment& s : e) {
        const bool on_side = (s[0][0] == s[1][0] && (s[0][0] == 0.0f || s[0][0] == 1.0f)) ||
                             (s[0][1] == s[1][1] && (s[0][1] == 0.0f || s[0][1] == 1.0f));
        CHECK(on_side);
    }
}

TEST_CASE("render edges: a fold is an edge only beyond the angle") {
    const Float3 o{0, 0, 0}, top{0, 1, 0};
    CHECK_FALSE(has(feature_edges(fold(20.0)), o, top));
    CHECK(has(feature_edges(fold(40.0)), o, top));
    CHECK_FALSE(has(feature_edges(fold(40.0), 60.0f), o, top));
    // The boundary of the folded strip is there either way: 6 sides, plus the fold at 40 degrees.
    CHECK(feature_edges(fold(20.0)).size() == 6);
    CHECK(feature_edges(fold(40.0)).size() == 7);
}

TEST_CASE("render edges: winding does not matter, degenerate triangles are ignored") {
    TriangleBatch b{"slice", Shading::field, true, false, {}};
    quad(b, {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, true);  // the second triangle wound the other way
    quad(b, {1, 0, 0}, {2, 0, 0}, {2, 1, 0}, {1, 1, 0});
    triangle(b, {0, 0, 0}, {1, 0, 0}, {2, 0, 0});  // zero area, on the boundary
    triangle(b, {5, 5, 5}, {5, 5, 5}, {6, 6, 6});  // a point and a line
    const auto e = feature_edges(b);
    CHECK(e.size() == 6);  // the 2 x 1 strip's boundary
    CHECK_FALSE(has(e, {0, 0, 0}, {1, 1, 0}));
    CHECK_FALSE(has(e, {1, 0, 0}, {1, 1, 0}));
    CHECK_FALSE(has(e, {5, 5, 5}, {6, 6, 6}));
}

TEST_CASE("render edges: the result is sorted and the same whatever the triangle order") {
    TriangleBatch b = cube();
    const auto first = feature_edges(b);
    CHECK(feature_edges(b) == first);
    CHECK(std::ranges::is_sorted(first, [](const Segment& x, const Segment& y) {
        return std::ranges::lexicographical_compare(std::array{x[0], x[1]}, std::array{y[0], y[1]});
    }));
    // Shuffle whole triangles, and rotate each triangle's vertices.
    std::vector<std::array<SceneVertex, 3>> tris;
    for (std::size_t t = 0; t < b.vertices.size(); t += 3) tris.push_back({b.vertices[t], b.vertices[t + 1], b.vertices[t + 2]});
    std::mt19937 rng(7);
    std::ranges::shuffle(tris, rng);
    b.vertices.clear();
    for (std::size_t k = 0; k < tris.size(); ++k) {
        for (std::size_t i = 0; i < 3; ++i) b.vertices.push_back(tris[k][(i + k) % 3]);
    }
    CHECK(feature_edges(b) == first);
    // -0 and +0 are one position.
    TriangleBatch z = cube();
    for (SceneVertex& v : z.vertices) {
        for (float& c : v.position) {
            if (c == 0.0f) c = -0.0f;
        }
    }
    z.vertices.resize(18);  // three faces with -0 ...
    const TriangleBatch whole = cube();
    z.vertices.insert(z.vertices.end(), whole.vertices.begin() + 18, whole.vertices.end());  // ... three with +0
    CHECK(feature_edges(z).size() == 12);
}

TEST_CASE("render edges: set_look adds one feature-edge batch to enhanced 3D scenes only") {
    FieldScene s;
    s.dimension = 3;
    s.triangles.push_back(cube());
    TriangleBatch arrows = cube("arrows");
    TriangleBatch iso = cube("isosurface");
    TriangleBatch shell = cube("outline");
    TriangleBatch unlit = cube("field");
    unlit.lit = false;
    for (const TriangleBatch& b : {arrows, iso, shell, unlit}) s.triangles.push_back(b);
    // Move the excluded batches away so their edges would be told apart.
    for (std::size_t k = 1; k < s.triangles.size(); ++k) {
        for (auto& v : s.triangles[k].vertices) v.position[1] += 10.0f * static_cast<float>(k);
    }
    s.lines.push_back({"mesh", {}, 1.0f, false, {{0, 0, 0}, {1, 1, 1}}});
    CHECK(s.look == Look::classic);

    set_look(s, Look::enhanced);
    CHECK(s.look == Look::enhanced);
    REQUIRE(s.lines.size() == 2);
    CHECK(s.lines[0].name == "mesh");
    const LineBatch& e = s.lines[1];
    CHECK(e.name == "feature edges");
    CHECK_FALSE(e.on_top);
    CHECK(e.points.size() == 24);  // the solid cube's 12 edges only
    for (const Float3& p : e.points) CHECK(p[1] <= 1.0f);

    set_look(s, Look::enhanced);  // again: still one batch
    CHECK(s.lines.size() == 2);
    set_look(s, Look::classic);
    CHECK(s.look == Look::classic);
    REQUIRE(s.lines.size() == 1);
    CHECK(s.lines[0].name == "mesh");

    FieldScene flat;
    flat.dimension = 2;
    flat.triangles.push_back(cube());
    set_look(flat, Look::enhanced);
    CHECK(flat.look == Look::enhanced);
    CHECK(flat.lines.empty());
}

TEST_CASE("render edges: the enhanced look's lights (KG3)") {
    FieldScene s;
    s.dimension = 3;
    s.lo = {0.0f, 0.0f, 0.0f};
    s.hi = {2.0f, 1.0f, 0.8f};
    s.camera.azimuth_deg = 35.0f;
    s.camera.elevation_deg = 25.0f;
    const double az = 35.0 * std::numbers::pi / 180.0, el = 25.0 * std::numbers::pi / 180.0;
    const std::array<double, 3> toward{std::cos(el) * std::sin(az), std::sin(el), std::cos(el) * std::cos(az)};
    // right = forward x +y, forward = -toward.
    const std::array<double, 3> right{toward[2], 0.0, -toward[0]};
    const double rn = std::hypot(right[0], right[2]);
    const std::array<double, 3> r{right[0] / rn, 0.0, right[2] / rn};

    SECTION("classic and 2D leave them zero") {
        const FrameLighting c = frame_lighting(s, 4.0 / 3.0);
        CHECK_FALSE(c.enhanced);
        CHECK_FALSE(c.perspective);
        for (const auto* v : {&c.key, &c.fill, &c.up, &c.eye}) CHECK(*v == std::array<float, 3>{});
        s.look = Look::enhanced;
        s.dimension = 2;
        CHECK_FALSE(frame_lighting(s, 4.0 / 3.0).enhanced);
    }
    SECTION("orthographic: the eye is the direction toward the viewer") {
        s.look = Look::enhanced;
        const FrameLighting l = frame_lighting(s, 4.0 / 3.0);
        CHECK(l.enhanced);
        CHECK_FALSE(l.perspective);
        for (std::size_t i = 0; i < 3; ++i) CHECK(l.eye[i] == Catch::Approx(toward[i]).margin(1e-6));
        for (const auto* v : {&l.key, &l.fill, &l.up, &l.eye}) CHECK(length3(*v) == Catch::Approx(1.0).margin(1e-6));
        CHECK(dot3(l.key, r) < 0.0);  // the key from the left
        CHECK(dot3(l.fill, r) > 0.0);  // the fill from the right
        CHECK(dot3(l.key, toward) > 0.0);
        CHECK(dot3(l.fill, toward) > 0.0);
        CHECK(l.up[1] > 0.0f);
    }
    SECTION("perspective: the eye is the camera's position, along the view direction") {
        s.look = Look::enhanced;
        s.camera.projection = Projection::perspective;
        const FrameLighting l = frame_lighting(s, 4.0 / 3.0);
        CHECK(l.perspective);
        const std::array<double, 3> centre{1.0, 0.5, 0.4};
        std::array<double, 3> d{l.eye[0] - centre[0], l.eye[1] - centre[1], l.eye[2] - centre[2]};
        const double n = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
        CHECK(n > 1.2);  // outside the bounds
        for (std::size_t i = 0; i < 3; ++i) CHECK(d[i] / n == Catch::Approx(toward[i]).margin(1e-5));
    }
}
