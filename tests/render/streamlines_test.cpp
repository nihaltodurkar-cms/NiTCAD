// Streamlines (Unit 25b, DS5, DS6, BG3): point location; a uniform field's lines are straight and end
// on the boundary; the rotational field (-y, x) closes after one turn within the integrator's
// tolerance; lines stop where the field is undefined; equal-current seeds are evenly spaced for a
// uniform contact current and spaced by its integral for a nonuniform one.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numbers>
#include <optional>
#include <utility>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/render/streamlines.hpp"

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

double dot(const Point& a, const Point& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
Point sub(const Point& a, const Point& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
double norm(const Point& a) { return std::sqrt(dot(a, a)); }
Point cross(const Point& a, const Point& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}

// A silicon bar on a 2D tensor mesh, 2 um by 1 um, ohmic contacts on x_min and x_max.
device::Device bar(const std::vector<double>& x, const std::vector<double>& y) {
    mesh::Mesh m = *mesh::make_tensor_grid(x, y);
    auto cells = *mesh::TensorCells::create(m, x, y);
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

// Edge currents of a current density J(x) = (j(y), 0): I = J . t A on each edge (A/cm^2 times cm).
template <typename F>
results::BiasPoint current_along_x(const device::Device& d, F&& j) {
    results::BiasPoint p;
    const auto pts = d.mesh().points();
    for (const mesh::Edge& e : d.mesh().edges()) {
        const Point& a = pts[static_cast<std::size_t>(e.first)];
        const Point& b = pts[static_cast<std::size_t>(e.second)];
        const double tx = (b[0] - a[0]) / e.length;
        p.edge_current_n.push_back(j(0.5 * (a[1] + b[1])) * tx * e.coupling_area);
        p.edge_current_p.push_back(0.0);
    }
    return p;
}

}  // namespace

TEST_CASE("render streamlines: point location in tensor and generated meshes") {
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
    const Geometry g = *Geometry::simplex(gm->mesh, gm->simplices);
    const Locator loc(g);
    CHECK(std::abs(loc.size() - L) <= 1e-18);
    const auto pts = gm->mesh.points();
    std::size_t found = 0;
    std::optional<std::size_t> hint;
    for (int i = 0; i < 2000; ++i) {
        // A deterministic scatter over the cube.
        const Point x{L * std::fmod(i * 0.7548776662, 1.0), L * std::fmod(i * 0.5698402910, 1.0),
                      L * std::fmod(i * 0.3477499440 + 0.01, 1.0)};
        const auto h = loc.locate(x, i % 2 ? hint : std::nullopt);
        REQUIRE(h.has_value());
        hint = h->element;
        Point back{};
        for (std::size_t k = 0; k < 4; ++k) {
            const Point& p = pts[static_cast<std::size_t>(g.element(h->element)[k])];
            for (std::size_t a = 0; a < 3; ++a) back[a] += h->weights[k] * p[a];
            CHECK(h->weights[k] >= -1e-10);
        }
        CHECK(norm(sub(back, x)) <= 1e-12 * L);
        ++found;
    }
    CHECK(found == 2000);
    CHECK_FALSE(loc.locate({1.5 * L, 0.5 * L, 0.5 * L}).has_value());
    CHECK_FALSE(loc.locate({0.5 * L, 0.5 * L, -1e-3 * L}).has_value());
}

TEST_CASE("render streamlines: a uniform field's lines are straight and end on the boundary (BG3)") {
    const auto x = graded(0.0, 2e-4, 9), y = graded(0.0, 1e-4, 6);
    const device::Device d = bar(x, y);
    const Geometry g = *Geometry::tensor(d.mesh(), *d.cells());
    const Locator loc(g);
    const Point dir{0.8, 0.6, 0.0};
    const std::vector<Point> field(d.mesh().node_count(), Point{3.2e3, 2.4e3, 0.0});  // along dir
    const Point seed{0.7e-4, 0.3e-4, 0.0};
    const Streamline s = trace_streamline(loc, field, seed);
    REQUIRE(s.points.size() >= 3);
    CHECK(s.forward == StreamlineStop::boundary);
    CHECK(s.backward == StreamlineStop::boundary);
    for (const Point& p : s.points) CHECK(std::abs(cross(sub(p, seed), dir)[2]) <= 1e-12 * 2e-4);
    // Forward it leaves through y = 1 um, backward through y = 0 (the line x - 0.7 = (y - 0.3) 4/3).
    CHECK(std::abs(s.points.back()[1] - 1e-4) <= 1e-6 * 2e-4);
    CHECK(std::abs(s.points.front()[1]) <= 1e-6 * 2e-4);
    // In order along the field.
    for (std::size_t i = 1; i < s.points.size(); ++i) CHECK(dot(sub(s.points[i], s.points[i - 1]), dir) > 0.0);
}

TEST_CASE("render streamlines: the rotational field closes after one turn (BG3)") {
    const auto u = graded(0.0, 2e-4, 21);
    const device::Device d = bar(u, u);
    const Geometry g = *Geometry::tensor(d.mesh(), *d.cells());
    const Locator loc(g);
    const Point c{1.03e-4, 0.98e-4, 0.0};
    std::vector<Point> field;
    for (const Point& p : d.mesh().points()) field.push_back({-(p[1] - c[1]), p[0] - c[0], 0.0});
    const double r = 0.6e-4;
    const Point seed{c[0] + r, c[1], 0.0};
    StreamlineOptions o;
    o.backward = false;
    o.max_length = 2.0 * std::numbers::pi * r;
    for (const double tol : {1e-5, 1e-7}) {
        o.tolerance = tol;
        const Streamline s = trace_streamline(loc, field, seed, o);
        CHECK(s.forward == StreamlineStop::length);
        CHECK(s.points.size() >= 50);  // steps of at most half an element: a smooth drawn circle
        double off = 0.0;
        for (const Point& p : s.points) off = std::max(off, std::abs(norm(sub(p, c)) - r));
        const double closure = norm(sub(s.points.back(), seed));
        std::printf("rotation: tolerance %.0e, %zu points, closure %.3e cm, radius error %.3e cm\n", tol,
                    s.points.size(), closure, off);
        // Within the integrator's tolerance: the error per step times the steps taken.
        const double budget = tol * loc.size() * static_cast<double>(s.points.size());
        CHECK(closure <= budget);
        CHECK(off <= budget);
    }
    // On a coarse mesh (the field is linear, so still exact) half an element is a large step, and only
    // the error control keeps the circle closed.
    const auto coarse = graded(0.0, 2e-4, 5);
    const device::Device dc = bar(coarse, coarse);
    const Geometry gc = *Geometry::tensor(dc.mesh(), *dc.cells());
    const Locator lc(gc);
    std::vector<Point> fc;
    for (const Point& p : dc.mesh().points()) fc.push_back({-(p[1] - c[1]), p[0] - c[0], 0.0});
    o.tolerance = 1e-12;
    const Streamline s = trace_streamline(lc, fc, seed, o);
    const double closure = norm(sub(s.points.back(), seed));
    std::printf("rotation, coarse mesh: %zu points, closure %.3e cm\n", s.points.size(), closure);
    CHECK(closure <= 1e-12 * lc.size() * static_cast<double>(s.points.size()));
}

TEST_CASE("render streamlines: a kinked field's line is exact within the tolerance (DS5, the step control)") {
    // v = (1, k |x - xm|) with xm on a node column: linear on every element (so represented exactly),
    // but its gradient jumps at xm. The line from (x0, y0) is y = y0 + k [(xm - x0)^2 - (xm - x)^2] / 2
    // before xm and y0 + k [(xm - x0)^2 + (x - xm)^2] / 2 after. Only rejecting the steps that cross
    // the kink keeps it within the tolerance there.
    std::vector<double> u;
    for (int i = 0; i < 9; ++i) u.push_back(2e-4 * i / 8);
    const device::Device d = bar(u, u);
    const Geometry g = *Geometry::tensor(d.mesh(), *d.cells());
    const Locator loc(g);
    const double xm = 1e-4, k = 1e4;
    std::vector<Point> field;
    for (const Point& p : d.mesh().points()) field.push_back({1.0, k * std::abs(p[0] - xm), 0.0});
    const double x0 = 0.1e-4, y0 = 0.2e-4;
    const auto exact = [&](double x) {
        return x <= xm ? y0 + k * ((xm - x0) * (xm - x0) - (xm - x) * (xm - x)) / 2.0
                       : y0 + k * ((xm - x0) * (xm - x0) + (x - xm) * (x - xm)) / 2.0;
    };
    StreamlineOptions o;
    o.backward = false;
    o.tolerance = 1e-11;
    const Streamline s = trace_streamline(loc, field, {x0, y0, 0.0}, o);
    CHECK(s.forward == StreamlineStop::boundary);
    double worst = 0.0;
    for (const Point& p : s.points) worst = std::max(worst, std::abs(p[1] - exact(p[0])));
    std::printf("kinked field: %zu points, worst deviation %.3e cm\n", s.points.size(), worst);
    CHECK(worst <= o.tolerance * loc.size() * static_cast<double>(s.points.size()));
    CHECK(std::abs(s.points.back()[0] - 2e-4) <= 1e-6 * 2e-4);
}

TEST_CASE("render streamlines: lines stop where the field is undefined or vanishes (DS5)") {
    const auto x = graded(0.0, 2e-4, 9), y = graded(0.0, 1e-4, 6);
    const device::Device d = bar(x, y);
    const Geometry g = *Geometry::tensor(d.mesh(), *d.cells());
    const Locator loc(g);
    std::vector<Point> field(d.mesh().node_count(), Point{1.0, 0.0, 0.0});
    const double nan = std::nan("");
    for (std::size_t i = 0; i < field.size(); ++i) {
        if (d.mesh().points()[i][0] >= 1.5e-4) field[i] = {nan, nan, nan};
    }
    const Streamline s = trace_streamline(loc, field, {0.2e-4, 0.5e-4, 0.0});
    CHECK(s.forward == StreamlineStop::no_field);
    CHECK(s.backward == StreamlineStop::boundary);
    CHECK(s.points.back()[0] < 1.5e-4);
    CHECK(s.points.back()[0] > 1.2e-4);  // through the last cell with a field
    // A sink: v = (1 um - x, 0) points at the plane x = 1 um from both sides; the line ends there in a
    // few steps instead of oscillating to its step limit.
    std::vector<Point> sink;
    for (const Point& q : d.mesh().points()) sink.push_back({1e-4 - q[0], 0.0, 0.0});
    StreamlineOptions forward;
    forward.backward = false;
    const Streamline k = trace_streamline(loc, sink, {0.3e-4, 0.5e-4, 0.0}, forward);
    CHECK(k.forward == StreamlineStop::no_field);
    CHECK(std::abs(k.points.back()[0] - 1e-4) <= 0.1e-4);
    CHECK(k.points.size() < 200);
    // A zero field: the seed alone.
    const std::vector<Point> zero(d.mesh().node_count(), Point{});
    const Streamline z = trace_streamline(loc, zero, {0.2e-4, 0.5e-4, 0.0});
    CHECK(z.points.size() == 1);
    CHECK(z.forward == StreamlineStop::no_field);
}

TEST_CASE("render streamlines: equal-current seeds on a contact (BG3, BQ2)") {
    const auto x = graded(0.0, 2e-4, 9);
    const auto y = graded(0.0, 1e-4, 11);
    const device::Device d = bar(x, y);
    const Geometry g = *Geometry::tensor(d.mesh(), *d.cells());
    const Locator loc(g);
    const std::size_t N = 8;
    // A uniform current: the contact's current per node is j times its part of the contact, and the
    // seeds lie at the midpoints of N equal parts of it.
    const auto uniform = current_along_x(d, [](double) { return 5.0; });
    const auto iu = contact_node_currents(d, uniform, 0);
    REQUIRE(iu.size() == d.contacts()[0].nodes.size());
    const auto su = contact_seeds(d, loc, 0, N, iu);
    REQUIRE(su.size() == N);
    for (std::size_t i = 0; i < N; ++i) {
        CHECK(std::abs(su[i][0]) <= 1e-18);
        CHECK(std::abs(su[i][1] - (i + 0.5) / N * 1e-4) <= 1e-12 * 1e-4);
    }
    // Without currents: evenly by length, the same points.
    const auto sl = contact_seeds(d, loc, 0, N);
    for (std::size_t i = 0; i < N; ++i) CHECK(std::abs(sl[i][1] - su[i][1]) <= 1e-12 * 1e-4);

    // j(y) = 1 + 3 y / H: each node's current is spread evenly over its part of the contact (the
    // solver's resolution), so the seeds split the integral of that piecewise-constant density into
    // N equal shares.
    const double H = 1e-4;
    const auto ramp = current_along_x(d, [&](double yy) { return 1.0 + 3.0 * yy / H; });
    const auto ir = contact_node_currents(d, ramp, 0);
    const auto sr = contact_seeds(d, loc, 0, N, ir);
    REQUIRE(sr.size() == N);
    // The density per node, over [y_k - h_below / 2, y_k + h_above / 2].
    const auto nodes = d.contacts()[0].nodes;
    std::vector<std::pair<double, double>> parts;  // (y, current) per contact node
    for (std::size_t k = 0; k < nodes.size(); ++k) parts.emplace_back(d.mesh().points()[static_cast<std::size_t>(nodes[k])][1], ir[k]);
    std::ranges::sort(parts);
    const auto integral = [&](double to) {
        double sum = 0.0;
        for (std::size_t k = 0; k < parts.size(); ++k) {
            const double lo = k == 0 ? parts[0].first : 0.5 * (parts[k - 1].first + parts[k].first);
            const double hi = k + 1 == parts.size() ? parts[k].first : 0.5 * (parts[k].first + parts[k + 1].first);
            const double density = parts[k].second / (hi - lo);
            sum += density * std::clamp(to - lo, 0.0, hi - lo);
        }
        return sum;
    };
    const double total = integral(H);
    for (std::size_t i = 0; i < N; ++i) {
        CHECK(std::abs(integral(sr[i][1]) - (i + 0.5) / N * total) <= 1e-9 * total);
    }
    // Denser where the current is: the last gap is smaller than the first.
    CHECK(sr[N - 1][1] - sr[N - 2][1] < sr[1][1] - sr[0][1]);
    // The current into the device through the left contact equals the right's out of it.
    double in = 0.0, out = 0.0;
    for (const double i : ir) in += i;
    for (const double i : contact_node_currents(d, ramp, 1)) out += i;
    CHECK(std::abs(in + out) <= 1e-12 * std::abs(in));
}

TEST_CASE("render streamlines: segment and plane seeds") {
    const auto s = segment_seeds({0, 0, 0}, {4e-4, 0, 0}, 4);
    REQUIRE(s.size() == 4);
    CHECK(s[0][0] == 0.5e-4);
    CHECK(s[3][0] == 3.5e-4);
    const auto u = graded(0.0, 1e-4, 6);
    const mesh::Mesh m = *mesh::make_tensor_grid(u, u, u);
    const auto cells = *mesh::TensorCells::create(m, u, u, u);
    const Geometry g = *Geometry::tensor(m, cells);
    Plane p;
    p.normal = {0.0, 0.0, 1.0};
    p.offset = 0.37e-4;
    const auto ps = plane_seeds(g, p, 50);
    REQUIRE(ps.size() == 50);
    double mx = 0.0, my = 0.0;
    for (const Point& q : ps) {
        CHECK(std::abs(q[2] - 0.37e-4) <= 1e-15);
        CHECK(q[0] >= 0.0);
        CHECK(q[0] <= 1e-4);
        CHECK(q[1] >= 0.0);
        CHECK(q[1] <= 1e-4);
        mx += q[0] / 50.0;
        my += q[1] / 50.0;
    }
    // Spread over the square: the mean near its centre.
    CHECK(std::abs(mx - 0.5e-4) <= 0.1e-4);
    CHECK(std::abs(my - 0.5e-4) <= 0.1e-4);
}
