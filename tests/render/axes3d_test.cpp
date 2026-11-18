// The axis box (Unit 26h, JD1, JQ1, JQ2; DECISIONS.md JG2): which faces, which edges, the ticks and labels.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <vector>

#include "NiTCAD/render/axes3d.hpp"
#include "NiTCAD/render/image.hpp"
#include "NiTCAD/render/scene.hpp"

using namespace NiTCAD::render;

namespace {

constexpr int W = 960, H = 720;

// A 3D scene that is only its bounds: 2 x 1 x 0.8 um at the default camera.
FieldScene box_scene(float az = 35.0f, float el = 25.0f) {
    FieldScene s;
    s.dimension = 3;
    s.title = "Regions";
    s.lo = {0.0f, 0.0f, 0.0f};
    s.hi = {2.0f, 1.0f, 0.8f};
    s.camera.target = {1.0f, 0.5f, 0.4f};
    s.camera.azimuth_deg = az;
    s.camera.elevation_deg = el;
    s.look = Look::enhanced;
    return s;
}

// Toward the eye, as make_view sets it from the camera's angles.
std::array<double, 3> eye_direction(const FieldScene& s) {
    const double az = s.camera.azimuth_deg * std::numbers::pi / 180.0;
    const double el = s.camera.elevation_deg * std::numbers::pi / 180.0;
    return {std::cos(el) * std::sin(az), std::sin(el), std::cos(el) * std::cos(az)};
}

double dot(const std::array<double, 3>& a, const Float3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

bool same(const Float3& a, const Float3& b) { return a[0] == b[0] && a[1] == b[1] && a[2] == b[2]; }

// The convex hull of points (monotone chain), counterclockwise in pixel coordinates.
std::vector<std::array<float, 2>> hull(std::vector<std::array<float, 2>> p) {
    std::ranges::sort(p);
    const auto cross = [](const auto& o, const auto& a, const auto& b) {
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
    };
    std::vector<std::array<float, 2>> h(2 * p.size());
    std::size_t k = 0;
    for (std::size_t i = 0; i < p.size(); ++i) {
        while (k >= 2 && cross(h[k - 2], h[k - 1], p[i]) <= 0.0f) --k;
        h[k++] = p[i];
    }
    for (std::size_t i = p.size() - 1, t = k + 1; i-- > 0;) {
        while (k >= t && cross(h[k - 2], h[k - 1], p[i]) <= 0.0f) --k;
        h[k++] = p[i];
    }
    h.resize(k - 1);
    return h;
}

bool inside_hull(const std::vector<std::array<float, 2>>& h, float x, float y) {
    for (std::size_t i = 0; i < h.size(); ++i) {
        const auto& a = h[i];
        const auto& b = h[(i + 1) % h.size()];
        if ((b[0] - a[0]) * (y - a[1]) - (b[1] - a[1]) * (x - a[0]) < 0.0f) return false;
    }
    return true;
}

// Every JG2 property of one camera.
void check_box(const FieldScene& s) {
    const AxisBox box = axis_box(s, W, H);
    const auto e = eye_direction(s);
    CAPTURE(s.camera.azimuth_deg, s.camera.elevation_deg);
    // The faces turned away: their outward normal points away from the eye (ties to the lo side).
    for (int a = 0; a < 3; ++a) {
        const double toward = e[static_cast<std::size_t>(a)];
        if (std::abs(toward) > 1e-6) CHECK(box.away_hi[static_cast<std::size_t>(a)] == (toward < 0.0));
        else CHECK_FALSE(box.away_hi[static_cast<std::size_t>(a)]);
    }
    // Their edges: nine, each on an away face.
    CHECK(box.edges.size() == 9);
    const auto on_away = [&](const Float3& p, int axis) {
        const std::size_t a = static_cast<std::size_t>(axis);
        return p[a] == (box.away_hi[a] ? s.hi[a] : s.lo[a]);
    };
    for (const auto& [p, q] : box.edges) {
        bool any = false;
        for (int a = 0; a < 3; ++a) any = any || (on_away(p, a) && on_away(q, a));
        CHECK(any);
    }
    // The silhouette in pixels, for "outside the device".
    std::vector<std::array<float, 2>> corners;
    for (int k = 0; k < 8; ++k) {
        const Float3 c{(k & 1) ? s.hi[0] : s.lo[0], (k & 2) ? s.hi[1] : s.lo[1], (k & 4) ? s.hi[2] : s.lo[2]};
        const ScreenPoint p = screen_position(s, W, H, c);
        corners.push_back({p.x, p.y});
    }
    const auto outline = hull(corners);
    struct Rect {
        float x, y, hw, hh;
    };
    std::vector<Rect> labels;
    for (const BoxAxis& ax : box.axes) {
        const std::size_t a = static_cast<std::size_t>(ax.axis);
        const std::size_t b = (a + 1) % 3, c = (a + 2) % 3;
        INFO("axis " << ax.axis);
        // The edge: parallel to the axis, on exactly one away face (an outline edge, JQ1).
        CHECK(ax.edge[0][b] == ax.edge[1][b]);
        CHECK(ax.edge[0][c] == ax.edge[1][c]);
        CHECK(on_away(ax.edge[0], static_cast<int>(b)) != on_away(ax.edge[0], static_cast<int>(c)));
        // Of the two outline edges, x and z take the nearer, y the one further left.
        Float3 other = ax.edge[0];
        other[b] = other[b] == s.lo[b] ? s.hi[b] : s.lo[b];
        other[c] = other[c] == s.lo[c] ? s.hi[c] : s.lo[c];
        Float3 mid = ax.edge[0], other_mid = other;
        mid[a] = other_mid[a] = 0.5f * (s.lo[a] + s.hi[a]);
        if (a == 1) CHECK(screen_position(s, W, H, mid).x <= screen_position(s, W, H, other_mid).x);
        else CHECK(dot(e, mid) >= dot(e, other_mid));
        // Ticks: 1-2-5 steps in um covering the bounds, 3 to 8 of them.
        REQUIRE(ax.ticks.size() >= 3);
        CHECK(ax.ticks.size() <= 8);
        const double step = ax.ticks[1].value - ax.ticks[0].value;
        const double mantissa = step / std::pow(10.0, std::floor(std::log10(step) + 1e-9));
        CHECK((std::abs(mantissa - 1.0) < 1e-6 || std::abs(mantissa - 2.0) < 1e-6 || std::abs(mantissa - 5.0) < 1e-6));
        CHECK(ax.ticks.front().value - s.lo[a] < step + 1e-6);
        CHECK(s.hi[a] - ax.ticks.back().value < step + 1e-6);
        for (std::size_t k = 0; k < ax.ticks.size(); ++k) {
            const AxisLabel& t = ax.ticks[k];
            if (k > 0) CHECK(std::abs(t.value - ax.ticks[k - 1].value - step) < 1e-6 * step + 1e-9);
            CHECK(t.value >= s.lo[a] - 1e-6);
            CHECK(t.value <= s.hi[a] + 1e-6);
            // On its edge, at its value, at its pixel.
            CHECK(t.position[b] == ax.edge[0][b]);
            CHECK(t.position[c] == ax.edge[0][c]);
            CHECK(std::abs(t.position[a] - static_cast<float>(t.value)) < 1e-6f);
            const ScreenPoint p = screen_position(s, W, H, t.position);
            CHECK(std::abs(p.x - t.x) < 0.01f);
            CHECK(std::abs(p.y - t.y) < 0.01f);
            if (!t.kept) continue;
            // The label outside the device's silhouette.
            CHECK_FALSE(inside_hull(outline, t.label_x, t.label_y));
            labels.push_back({t.label_x, t.label_y, 3.25f * static_cast<float>(t.text.size()), 7.0f});
        }
        CHECK(ax.title == (a == 0 ? "x (\xC2\xB5m)" : a == 1 ? "y (\xC2\xB5m)" : "z (\xC2\xB5m)"));
        CHECK_FALSE(inside_hull(outline, ax.title_x, ax.title_y));
    }
    // No two kept labels overlap, across all the axes.
    for (std::size_t i = 0; i < labels.size(); ++i) {
        for (std::size_t j = i + 1; j < labels.size(); ++j) {
            const bool overlap = std::abs(labels[i].x - labels[j].x) < labels[i].hw + labels[j].hw &&
                                 std::abs(labels[i].y - labels[j].y) < labels[i].hh + labels[j].hh;
            CHECK_FALSE(overlap);
        }
    }
}

}  // namespace

TEST_CASE("render axes3d: the default camera's box, edges and ticks (Unit 26h, JG2)") {
    const FieldScene s = box_scene();
    const AxisBox box = axis_box(s, W, H);
    // Seen from +x +y +z: the lo faces are turned away.
    CHECK(box.away_hi == std::array<bool, 3>{false, false, false});
    REQUIRE(box.axes.size() == 3);
    // x along the lower front edge (y lo, z hi), y up the front left edge (x lo, z hi), z along the lower
    // right side (x hi, y lo): JD1's description of the default view.
    const auto edge_of = [&](int a) { return std::ranges::find_if(box.axes, [&](const BoxAxis& x) { return x.axis == a; }); };
    CHECK(edge_of(0)->edge[0][1] == 0.0f);
    CHECK(edge_of(0)->edge[0][2] == 0.8f);
    CHECK(edge_of(1)->edge[0][0] == 0.0f);
    CHECK(edge_of(1)->edge[0][2] == 0.8f);
    CHECK(edge_of(2)->edge[0][0] == 2.0f);
    CHECK(edge_of(2)->edge[0][1] == 0.0f);
    // x: 0, 0.5, ..., 2 um.
    std::vector<double> xs;
    for (const AxisLabel& t : edge_of(0)->ticks) xs.push_back(t.value);
    CHECK(xs == std::vector<double>{0.0, 0.5, 1.0, 1.5, 2.0});
    check_box(s);
}

TEST_CASE("render axes3d: the edges follow the camera (Unit 26h, JG2, JQ1)") {
    for (const float az : {35.0f, 125.0f, 215.0f, 305.0f}) {
        for (const float el : {25.0f, 60.0f, -60.0f}) check_box(box_scene(az, el));
    }
    // The named views' angles too (26h·3): looking along an axis drops that axis's ticks.
    for (const auto& [az, el] : {std::pair{0.0f, 0.0f}, std::pair{90.0f, 0.0f}, std::pair{0.0f, 90.0f}}) {
        const FieldScene s = box_scene(az, el);
        const AxisBox box = axis_box(s, W, H);
        CAPTURE(az, el);
        CHECK(box.axes.size() == 2);
        check_box(s);
    }
    // Perspective: the faces turned away from the camera's position.
    FieldScene p = box_scene();
    p.camera.projection = Projection::perspective;
    check_box(p);
}

TEST_CASE("render axes3d: no z axis for a box of zero thickness; none in 2D (Unit 26h, JG2)") {
    FieldScene s = box_scene();
    s.hi[2] = s.lo[2];
    const AxisBox box = axis_box(s, W, H);
    CHECK(std::ranges::none_of(box.axes, [](const BoxAxis& a) { return a.axis == 2; }));
    CHECK_FALSE(box.axes.empty());
    FieldScene flat = box_scene();
    flat.dimension = 2;
    const AxisBox none = axis_box(flat, W, H);
    CHECK(none.edges.empty());
    CHECK(none.axes.empty());
}
