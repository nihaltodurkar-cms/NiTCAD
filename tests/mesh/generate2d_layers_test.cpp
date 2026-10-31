// Unit 54 slice 2 (DECISIONS.md UD5, UD7; UC6 part, UC2 part): 2D structured layers with signed
// lines on interfaces. The lattice's nodes are bit-identical to their formula, its dual is the
// tensor box method, the result is certified (P4), and every UD5 rule is refused with its message.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <initializer_list>
#include <map>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/interfaces.hpp"
#include "NiTCAD/mesh/mesh.hpp"

using namespace NiTCAD::mesh;
using NiTCAD::base::ErrorCode;

namespace {

RegionPolygon rect(double x0, double y0, double x1, double y1, std::int32_t region) {
    return {{{x0, y0, 0}, {x1, y0, 0}, {x1, y1, 0}, {x0, y1, 0}}, region};
}

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::max(std::abs(a), std::abs(b)); }

double total(std::span<const double> v) {
    double s = 0.0;
    for (const double x : v) s += x;
    return s;
}

// A layer's lattice as the generator builds it: a + t u + o n, n the left normal of a -> b.
struct Lattice {
    Point a;
    Point u;
    Point n;
    std::vector<double> along;
    std::vector<double> lines;  // with the straddles

    Lattice(const Point& from, const Point& to, std::vector<double> t, std::vector<double> o)
        : a(from), along(std::move(t)), lines(std::move(o)) {
        const double len = std::hypot(to[0] - from[0], to[1] - from[1]);
        u = {(to[0] - from[0]) / len, (to[1] - from[1]) / len, 0.0};
        n = {-u[1], u[0], 0.0};
        std::ranges::sort(along);
        std::ranges::sort(lines);
    }
    [[nodiscard]] Point at(std::size_t i, std::size_t j) const {
        const double t = along[i], o = lines[j];
        return {a[0] + t * u[0] + o * n[0], a[1] + t * u[1] + o * n[1], 0.0};
    }
    [[nodiscard]] Point at_t(double t, double o) const { return {a[0] + t * u[0] + o * n[0], a[1] + t * u[1] + o * n[1], 0.0}; }
    // The tensor dual length of entry i of v, cut at the lattice's ends (or not, where the lattice
    // continues past them in the mesh).
    static double dual(const std::vector<double>& v, std::size_t i) {
        const double lo = i == 0 ? 0.0 : 0.5 * (v[i] - v[i - 1]);
        const double hi = i + 1 == v.size() ? 0.0 : 0.5 * (v[i + 1] - v[i]);
        return lo + hi;
    }
};

std::optional<NodeId> node_at(const GeneratedMesh& m, const Point& p) {
    for (std::size_t i = 0; i < m.mesh.node_count(); ++i) {
        if (m.mesh.points()[i][0] == p[0] && m.mesh.points()[i][1] == p[1]) return static_cast<NodeId>(i);
    }
    return std::nullopt;
}

std::map<std::pair<NodeId, NodeId>, double> couplings(const GeneratedMesh& m) {
    std::map<std::pair<NodeId, NodeId>, double> c;
    for (const Edge& e : m.mesh.edges()) c[{e.first, e.second}] = e.coupling_area;
    return c;
}

double coupling(const std::map<std::pair<NodeId, NodeId>, double>& c, NodeId a, NodeId b) {
    const auto it = c.find({std::min(a, b), std::max(a, b)});
    return it == c.end() ? 0.0 : it->second;
}

// Whether every lattice rectangle next to row j has its circumdisc within the lattice's offsets
// (or the lattice is closed by the boundary on that side), and row j is not an open edge.
bool row_inside(const Lattice& L, std::size_t j, bool closed_below, bool closed_above) {
    const std::size_t last = L.lines.size() - 1;
    if ((j == 0 && !closed_below) || (j == last && !closed_above)) return false;
    for (std::size_t r = j == 0 ? 0 : j - 1; r <= std::min(j, last - 1); ++r) {
        const double centre = 0.5 * (L.lines[r] + L.lines[r + 1]);
        for (std::size_t i = 0; i + 1 < L.along.size(); ++i) {
            const double radius = std::hypot(0.5 * (L.along[i + 1] - L.along[i]), 0.5 * (L.lines[r + 1] - L.lines[r]));
            if (!closed_below && centre - radius < L.lines.front()) return false;
            if (!closed_above && centre + radius > L.lines.back()) return false;
        }
    }
    return true;
}

// Every lattice node is a mesh node, bit for bit. A row is the tensor box method's (node volumes,
// couplings to lattice neighbours, none across the rectangles' diagonals) where every lattice
// rectangle next to it is Delaunay whatever lies outside: its circumdisc inside the lattice's
// offsets, or the lattice ending on the boundary on that side (`closed_below`, `closed_above`; the
// tests' lattices end on the boundary at both ends along the base). Returns the rows checked.
std::size_t check_tensor(const GeneratedMesh& m, const Lattice& L, bool closed_below, bool closed_above) {
    std::vector<std::vector<NodeId>> id(L.along.size(), std::vector<NodeId>(L.lines.size()));
    for (std::size_t i = 0; i < L.along.size(); ++i) {
        for (std::size_t j = 0; j < L.lines.size(); ++j) {
            const auto v = node_at(m, L.at(i, j));
            INFO("lattice node " << i << ", " << j << " at " << L.at(i, j)[0] << ", " << L.at(i, j)[1]);
            REQUIRE(v.has_value());
            id[i][j] = *v;
        }
    }
    const auto c = couplings(m);
    std::size_t checked = 0;
    for (std::size_t i = 0; i < L.along.size(); ++i) {
        for (std::size_t j = 0; j < L.lines.size(); ++j) {
            if (!row_inside(L, j, closed_below, closed_above)) continue;
            checked += i == 0 ? 1 : 0;
            INFO("node " << i << ", " << j);
            const double dt = Lattice::dual(L.along, i), dn = Lattice::dual(L.lines, j);
            INFO("volume " << m.mesh.volumes()[static_cast<std::size_t>(id[i][j])] << ", tensor " << dt * dn);
            CHECK(close(m.mesh.volumes()[static_cast<std::size_t>(id[i][j])], dt * dn, 1e-12));
            if (i + 1 < L.along.size()) CHECK(close(coupling(c, id[i][j], id[i + 1][j]), dn, 1e-12));
            if (j + 1 < L.lines.size()) CHECK(close(coupling(c, id[i][j], id[i][j + 1]), dt, 1e-12));
            // The diagonals of the lattice's rectangles are ties: no coupling.
            if (i + 1 < L.along.size() && j + 1 < L.lines.size()) {
                CHECK(std::abs(coupling(c, id[i][j], id[i + 1][j + 1])) <= 1e-12 * dt);
                CHECK(std::abs(coupling(c, id[i + 1][j], id[i][j + 1])) <= 1e-12 * dt);
            }
        }
    }
    return checked;
}

// The node at exactly one of the given points (a corner node is computed by one of its two layers'
// formulas, whichever owns it).
std::optional<NodeId> node_at_either(const GeneratedMesh& m, const Point& p, const Point& q) {
    const auto a = node_at(m, p);
    return a ? a : node_at(m, q);
}

// A tensor block X x Y of mesh nodes (ids[i][j] at (X[i], Y[j]), found by the caller): every node
// not on the block's edge whose four rectangles have their circumdiscs inside the block has the
// tensor box method's volume and couplings, and no coupling across the rectangles' diagonals.
// Returns the nodes checked.
std::size_t check_block(const GeneratedMesh& m, const std::vector<double>& X, const std::vector<double>& Y,
                        const std::vector<std::vector<NodeId>>& ids) {
    const auto c = couplings(m);
    const auto disc_inside = [&](std::size_t i, std::size_t j) {  // rectangle [i, i+1] x [j, j+1]
        const double cx = 0.5 * (X[i] + X[i + 1]), cy = 0.5 * (Y[j] + Y[j + 1]);
        const double r = std::hypot(0.5 * (X[i + 1] - X[i]), 0.5 * (Y[j + 1] - Y[j]));
        return cx - r >= X.front() && cx + r <= X.back() && cy - r >= Y.front() && cy + r <= Y.back();
    };
    std::size_t checked = 0;
    for (std::size_t i = 1; i + 1 < X.size(); ++i) {
        for (std::size_t j = 1; j + 1 < Y.size(); ++j) {
            if (!(disc_inside(i - 1, j - 1) && disc_inside(i, j - 1) && disc_inside(i - 1, j) && disc_inside(i, j))) continue;
            ++checked;
            INFO("block node " << i << ", " << j << " at " << X[i] << ", " << Y[j]);
            const NodeId v = ids[i][j];
            const double dx = Lattice::dual(X, i), dy = Lattice::dual(Y, j);
            CHECK(close(m.mesh.volumes()[static_cast<std::size_t>(v)], dx * dy, 1e-12));
            CHECK(close(coupling(c, v, ids[i + 1][j]), dy, 1e-12));
            CHECK(close(coupling(c, v, ids[i][j + 1]), dx, 1e-12));
            CHECK(std::abs(coupling(c, v, ids[i + 1][j + 1])) <= 1e-12 * dx);
            CHECK(std::abs(coupling(c, ids[i + 1][j], ids[i][j + 1])) <= 1e-12 * dx);
        }
    }
    return checked;
}

// The lines with the straddles the generator adds about the given interface offsets.
std::vector<double> with_straddles(std::vector<double> lines, std::initializer_list<double> interfaces, double s) {
    for (const double d : interfaces) {
        lines.push_back(d - 0.5 * s);
        lines.push_back(d + 0.5 * s);
    }
    std::ranges::sort(lines);
    return lines;
}

BoundaryLayer signed_layer(Point a, Point b, std::vector<double> along, std::vector<double> lines) {
    BoundaryLayer l;
    l.a = a;
    l.b = b;
    l.along = std::move(along);
    l.lines = std::move(lines);
    return l;
}

// The intersection of the lines p + x u and q + y v.
Point cross_lines(const Point& p, const Point& u, const Point& q, const Point& v) {
    const double det = u[0] * v[1] - u[1] * v[0];
    const double x = ((q[0] - p[0]) * v[1] - (q[1] - p[1]) * v[0]) / det;
    return {p[0] + x * u[0], p[1] + x * u[1], 0.0};
}

}  // namespace

TEST_CASE("corner blocks: a rectangular fin wrapped by three layers is one tensor lattice (UC6, UD6)",
          "[mesh][generate][layers]") {
    // Silicon fin x in [6, 14], y in [0, 10] on the bottom boundary, wrapped by 0.5 of oxide and
    // 1.5 of high-k, a fill to the box [0, 20] x [0, 14]. Layers on the fin's left side (upward),
    // its top, and its right side (downward); right-angle corners at the fin's top corners. Each
    // corner square is the product of the two layers' lines, so the frame is one tensor lattice.
    const double s = 0.1;
    const Geometry2D g{{rect(6, 0, 14, 10, 0),
                        {{{5.5, 0, 0}, {6, 0, 0}, {6, 10, 0}, {14, 10, 0}, {14, 0, 0}, {14.5, 0, 0}, {14.5, 10.5, 0}, {5.5, 10.5, 0}}, 1},
                        {{{4, 0, 0}, {5.5, 0, 0}, {5.5, 10.5, 0}, {14.5, 10.5, 0}, {14.5, 0, 0}, {16, 0, 0}, {16, 12, 0}, {4, 12, 0}}, 2},
                        {{{0, 0, 0}, {4, 0, 0}, {4, 12, 0}, {16, 12, 0}, {16, 0, 0}, {20, 0, 0}, {20, 14, 0}, {0, 14, 0}}, 3}},
                       {{"gate", {0, 14, 0}, {20, 14, 0}}}};
    MeshOptions2D o;
    o.max_size = 1.0;
    o.interface_resolution = s;
    const std::vector<double> top_lines{-2, -1.2, -0.6, -0.3, -0.15, 0.25, 1.0, 1.5, 2.4, 3.0};
    const std::vector<double> side_lines{-2, -1, -0.5, -0.2, 0.3, 0.8, 1.4, 2.5, 3.0};
    const BoundaryLayer left = signed_layer({6, 0, 0}, {6, 10, 0}, {0, 1, 2, 4, 6, 8, 10}, side_lines);
    const BoundaryLayer top = signed_layer({6, 10, 0}, {14, 10, 0}, {0, 2, 3, 4, 5, 6, 8}, top_lines);
    const BoundaryLayer right = signed_layer({14, 10, 0}, {14, 0, 0}, {0, 2, 4, 6, 8, 9, 10}, side_lines);
    o.layers = {top, left, right};
    const auto m = generate_mesh(g, o);
    INFO((m ? std::string("ok") : m.error().message));
    REQUIRE(m.has_value());

    const std::vector<double> A = with_straddles(top_lines, {0.0, 0.5, 2.0}, s);   // the top's lines, y - 10
    const std::vector<double> B = with_straddles(side_lines, {0.0, 0.5, 2.0}, s);  // the sides' lines, outward
    const Lattice LA({6, 10, 0}, {14, 10, 0}, {}, {});
    const Lattice LL({6, 0, 0}, {6, 10, 0}, {}, {});
    const Lattice LR({14, 10, 0}, {14, 0, 0}, {}, {});
    // The top band: x from the left square (6 - b), the top's own columns (6 + 3, 4, 5), the right
    // square (14 + b); y = 10 + a. A corner node is the top's (t = -b or 8 + b, line a) or the
    // side's (t = 10 + a or -a, line b).
    struct Column {
        double x;
        int corner;  // -1 left square, 0 the top's own, 1 right square
        double b;    // the side's line
    };
    std::vector<Column> cols;
    for (const double b : B) cols.push_back({6.0 - b, -1, b});
    for (const double t : {3.0, 4.0, 5.0}) cols.push_back({6.0 + t, 0, 0.0});
    for (const double b : B) cols.push_back({14.0 + b, 1, b});
    std::ranges::sort(cols, [](const Column& p, const Column& q) { return p.x < q.x; });
    std::vector<double> X, Y;
    for (const Column& c : cols) X.push_back(c.x);
    for (const double a : A) Y.push_back(10.0 + a);
    std::vector<std::vector<NodeId>> ids(X.size(), std::vector<NodeId>(Y.size()));
    for (std::size_t i = 0; i < cols.size(); ++i) {
        for (std::size_t j = 0; j < A.size(); ++j) {
            const Column& c = cols[i];
            const double a = A[j];
            std::optional<NodeId> v;
            if (c.corner == 0) {
                v = node_at(*m, LA.at_t(c.x - 6.0, a));
            } else if (c.corner < 0) {
                v = node_at_either(*m, LA.at_t(0.0 + -1.0 * c.b, a), LL.at_t(10.0 + 1.0 * a, c.b));
            } else {
                v = node_at_either(*m, LA.at_t(8.0 + 1.0 * c.b, a), LR.at_t(0.0 + -1.0 * a, c.b));
            }
            INFO("top band node at " << c.x << ", " << 10.0 + a);
            REQUIRE(v.has_value());
            ids[i][j] = *v;
        }
    }
    const std::size_t top_checked = check_block(*m, X, Y, ids);
    // The left band: x = 6 - b, y from the side's own columns (0, 1, 2, 4, 6) and the corner (10 + a).
    std::vector<double> XL, YL;
    for (const double b : B) XL.push_back(6.0 - b);
    std::ranges::sort(XL);
    std::vector<std::pair<double, bool>> rows;  // y, in the corner
    for (const double t : {0.0, 1.0, 2.0, 4.0, 6.0}) rows.push_back({t, false});
    for (const double a : A) rows.push_back({10.0 + a, true});
    std::vector<std::vector<NodeId>> idl(XL.size(), std::vector<NodeId>(rows.size()));
    for (std::size_t i = 0; i < XL.size(); ++i) {
        const double b = 6.0 - XL[i];
        for (std::size_t j = 0; j < rows.size(); ++j) {
            const auto [y, corner] = rows[j];
            const auto v = corner ? node_at_either(*m, LA.at_t(0.0 + -1.0 * b, y - 10.0), LL.at_t(10.0 + 1.0 * (y - 10.0), b))
                                  : node_at(*m, LL.at_t(y, b));
            INFO("left band node at " << XL[i] << ", " << y);
            REQUIRE(v.has_value());
            idl[i][j] = *v;
        }
    }
    for (const auto& [y, corner] : rows) YL.push_back(y);
    const std::size_t left_checked = check_block(*m, XL, YL, idl);
    std::printf("UC6 rectangular fin: %zu nodes; tensor-checked nodes: top band %zu, left band %zu\n",
                m->mesh.node_count(), top_checked, left_checked);
    CHECK(top_checked >= 60);
    CHECK(left_checked >= 40);
    // No node inside the frame but the lattice's.
    std::size_t inside = 0;
    for (const Point& p : m->simplices.points) {
        const bool in_top = p[0] > 3.0 && p[0] < 17.0 && p[1] > 8.0 && p[1] < 13.0;
        const bool in_sides = (p[0] > 3.0 && p[0] < 8.0 && p[1] < 13.0) || (p[0] > 12.0 && p[0] < 17.0 && p[1] < 13.0);
        inside += in_top || in_sides ? 1 : 0;
    }
    // The top band's interior, and each side band's interior columns on its own five columns and the
    // top band's lowest line (y = 8).
    const std::size_t lattice_inside = (X.size() - 2) * (Y.size() - 2) + 2 * (XL.size() - 2) * 6;
    CHECK(inside == lattice_inside);
    REQUIRE(m->certificate.has_value());
    INFO("flat " << m->certificate->report.max_flat_deviation << " corner " << m->certificate->report.max_corner_deviation);
    CHECK(m->certificate->report.max_flat_deviation <= 0.01);  // of P4a's tolerance: rounding (measured 2.0e-3)
    CHECK(m->certificate->report.max_corner_deviation <= 1e-12 * 20);  // the tensor's bends are exact (measured 1.8e-15)
    // Determinism and the order of the layers (the corner owner is the lower index).
    const auto again = generate_mesh(g, o);
    REQUIRE(again.has_value());
    CHECK(again->simplices.points == m->simplices.points);
    CHECK(again->simplices.cells == m->simplices.cells);
}

TEST_CASE("corner blocks: a tapered fin with 95-degree corners at its top and foot (UC6, UD6)", "[mesh][generate][layers]") {
    // A substrate y in [0, 4] under a fin tapering from x in [6, 14] at y = 4 to [6.875, 13.125] at
    // y = 14 (sides 5.0 degrees from vertical), wrapped with the floor by 0.5 of oxide and 1.0 of
    // high-k (uniform normal thickness), a fill to y = 19. Five layers: the floor left and right,
    // the two sides, the top; four corners at 95 degrees (two at the top, two at the foot). P4:
    // exact away from the bends (mirrored columns), the bends judged as corner zones.
    const double s = 0.1;
    const Point foot_l{6, 4, 0}, top_l{6.875, 14, 0}, top_r{13.125, 14, 0}, foot_r{14, 4, 0};
    const auto unit = [](const Point& p, const Point& q) {
        const double l = std::hypot(q[0] - p[0], q[1] - p[1]);
        return Point{(q[0] - p[0]) / l, (q[1] - p[1]) / l, 0.0};
    };
    const Point up_l = unit(foot_l, top_l), down_r = unit(top_r, foot_r);
    const Point out_l{-up_l[1], up_l[0], 0.0}, out_r{-down_r[1], down_r[0], 0.0};  // left normals: outward
    // The surface offset by d: the floor at 4 + d, the sides along their normals, the top at 14 + d.
    const auto offset = [&](double d) {
        const Point fl = cross_lines({0, 4 + d, 0}, {1, 0, 0}, {foot_l[0] + d * out_l[0], foot_l[1] + d * out_l[1], 0}, up_l);
        const Point tl = cross_lines({0, 14 + d, 0}, {1, 0, 0}, {foot_l[0] + d * out_l[0], foot_l[1] + d * out_l[1], 0}, up_l);
        const Point tr = cross_lines({0, 14 + d, 0}, {1, 0, 0}, {top_r[0] + d * out_r[0], top_r[1] + d * out_r[1], 0}, down_r);
        const Point fr = cross_lines({0, 4 + d, 0}, {1, 0, 0}, {top_r[0] + d * out_r[0], top_r[1] + d * out_r[1], 0}, down_r);
        return std::vector<Point>{{0, 4 + d, 0}, fl, tl, tr, fr, {20, 4 + d, 0}};
    };
    const auto band = [](std::vector<Point> lower, std::vector<Point> upper) {  // between two surfaces, counterclockwise
        std::ranges::reverse(upper);
        lower.insert(lower.end(), upper.begin(), upper.end());
        return lower;
    };
    const std::vector<Point> si_top{{0, 4, 0}, foot_l, top_l, top_r, foot_r, {20, 4, 0}};
    const std::vector<Point> ox = offset(0.5), hk = offset(1.5);
    std::vector<Point> silicon{{0, 0, 0}, {20, 0, 0}};
    for (auto it = si_top.rbegin(); it != si_top.rend(); ++it) silicon.push_back(*it);
    Geometry2D g{{{silicon, 0}, {band(si_top, ox), 1}, {band(ox, hk), 2}, {band(hk, {{0, 19, 0}, {20, 19, 0}}), 3}},
                 {{"body", {0, 0, 0}, {20, 0, 0}}, {"gate", {20, 19, 0}, {0, 19, 0}}}};
    MeshOptions2D o;
    o.max_size = 1.0;
    o.interface_resolution = s;
    const std::vector<double> lines{-1.5, -1.0, -0.6, -0.3, -0.15, 0.2, 0.35, 0.7, 1.1, 1.3, 1.75, 2.1, 2.5};
    const double side_len = std::hypot(0.875, 10.0);
    o.layers = {signed_layer({0, 4, 0}, foot_l, {0, 1, 2, 3, 6}, lines),
                signed_layer(foot_l, top_l, {0, 2.5, 4, 5.5, 7, 8.5, side_len}, lines),
                signed_layer(top_l, top_r, {0, 2, 3.125, 4.25, 6.25}, lines),
                signed_layer(top_r, foot_r, {0, side_len - 8.5, side_len - 7, side_len - 5.5, side_len - 4, side_len - 2.5, side_len}, lines),
                signed_layer(foot_r, {20, 4, 0}, {0, 3, 4, 5, 6}, lines)};
    const auto m = generate_mesh(g, o);
    INFO((m ? std::string("ok") : m.error().message));
    REQUIRE(m.has_value());
    REQUIRE(m->certificate.has_value());
    const auto& r = m->certificate->report;
    std::printf("UC6 tapered fin: %zu nodes, %zu triangles; P4 flat %.3e of the tolerance, corner %.3e (s = %.3g), "
                "%zu corner faces of %zu\n",
                m->mesh.node_count(), m->report.triangles, r.max_flat_deviation, r.max_corner_deviation, s, r.corner_faces,
                r.cross_faces);
    CHECK(m->report.interface_nodes == 0);  // every interface piece is under a lattice
    CHECK(r.max_flat_deviation <= 0.05);  // of P4a's tolerance: rounding (measured 6.6e-3)
    CHECK(r.max_corner_deviation <= 1e-12 * 20);  // measured 2.7e-15: one set of lines, each corner symmetric about its bisector
    CHECK(close(total(m->mesh.boundary()[0].areas), 20.0, 1e-12));
    CHECK(close(total(m->mesh.boundary()[1].areas), 20.0, 1e-12));
    // The layers' own columns, away from the corners, are exact lattice nodes: the top's at t = 2,
    // 3.125, 4.25 on every line.
    const Lattice LA(top_l, top_r, {}, {});
    for (const double t : {2.0, 3.125, 4.25}) {
        for (const double a : with_straddles(lines, {0.0, 0.5, 1.5}, s)) CHECK(node_at(*m, LA.at_t(t, a)).has_value());
    }
    // The general case: the sides on other lines (the same first and last), so no corner is
    // symmetric and its bends are judged as corner zones (P4b).
    const std::vector<double> side{-1.5, -0.8, -0.4, -0.2, 0.15, 0.3, 0.8, 1.0, 1.25, 1.9, 2.5};
    o.layers[1].lines = side;
    o.layers[3].lines = side;
    const auto u = generate_mesh(g, o);
    INFO((u ? std::string("ok") : u.error().message));
    REQUIRE(u.has_value());
    REQUIRE(u->certificate.has_value());
    const auto& ru = u->certificate->report;
    std::printf("UC6 tapered fin, unequal lines: %zu nodes; P4 flat %.3e of the tolerance, corner %.3e (s = %.3g), "
                "%zu corner faces of %zu\n",
                u->mesh.node_count(), ru.max_flat_deviation, ru.max_corner_deviation, s, ru.corner_faces, ru.cross_faces);
    CHECK(ru.max_flat_deviation <= 0.05);
    CHECK(ru.max_corner_deviation <= s);
}

TEST_CASE("corner blocks: refusals (UC2, UD6)", "[mesh][generate][layers]") {
    const double s = 0.1;
    const auto options = [&](std::vector<BoundaryLayer> layers) {
        MeshOptions2D o;
        o.max_size = 1.0;
        o.interface_resolution = s;
        o.layers = std::move(layers);
        return o;
    };
    const auto refused = [&](const Geometry2D& geo, std::vector<BoundaryLayer> layers, const std::string& text, std::size_t index) {
        const auto m = generate_mesh(geo, options(std::move(layers)));
        INFO(text);
        REQUIRE_FALSE(m.has_value());
        INFO(m.error().message);
        CHECK(m.error().code == ErrorCode::invalid_input);
        REQUIRE(m.error().context.has_value());
        CHECK(m.error().context->index == std::optional<std::size_t>{index});
        CHECK(m.error().message.find(text) != std::string::npos);
    };
    // A silicon block [2, 6] x [0, 4] on the bottom of an oxide box [0, 8] x [0, 8]; layers on its
    // top (0), left side (1, upward) and right side (2, downward), right-angle corners at (2, 4)
    // and (6, 4).
    const Geometry2D block{{rect(2, 0, 6, 4, 0), {{{0, 0, 0}, {2, 0, 0}, {2, 4, 0}, {6, 4, 0}, {6, 0, 0}, {8, 0, 0}, {8, 8, 0}, {0, 8, 0}}, 1}}, {}};
    const std::vector<double> lines{-1, -0.5, 0.2, 1};
    const auto top = [&](std::vector<double> along) { return signed_layer({2, 4, 0}, {6, 4, 0}, std::move(along), lines); };
    const auto left = [&](std::vector<double> along) { return signed_layer({2, 0, 0}, {2, 4, 0}, std::move(along), lines); };
    const auto right = [&](std::vector<double> along, std::vector<double> l) {
        return signed_layer({6, 4, 0}, {6, 0, 0}, std::move(along), std::move(l));
    };
    {
        const auto m = generate_mesh(block, options({top({0, 1, 2, 3, 4}), left({0, 2, 4}), right({0, 2, 4}, lines)}));
        INFO((m ? std::string("ok") : m.error().message));
        CHECK(m.has_value());  // the valid configuration the refusals vary
    }
    refused(block, {top({0, 1, 2, 3, 4}), left({0, 2, 4}), right({0, 2, 4}, {-1, -0.5, 0.2, 1.2})}, "do not span the same offsets", 2);
    refused(block, {top({0, 1, 2, 3.5, 4}), left({0, 2, 4}), right({0, 2, 4}, lines)}, "a given column inside a corner", 0);
    refused(block, {top({0, 1, 2, 3, 4}), left({0, 2, 4}), right({0, 0.5, 2, 4}, lines)}, "a given column inside a corner", 2);
    // A fourth layer on the right side, upward: three layers end at (6, 4).
    refused(block, {top({0, 1, 2, 3, 4}), left({0, 2, 4}), right({0, 2, 4}, lines), signed_layer({6, 0, 0}, {6, 4, 0}, {0, 2, 4}, lines)},
            "three or more layers meet at one end", 3);
    // A block 1.5 wide under layers 1 deep: the top's corners overlap.
    {
        const Geometry2D narrow{{rect(2, 0, 3.5, 4, 0), {{{0, 0, 0}, {2, 0, 0}, {2, 4, 0}, {3.5, 4, 0}, {3.5, 0, 0}, {8, 0, 0}, {8, 8, 0}, {0, 8, 0}}, 1}}, {}};
        refused(narrow,
                {left({0, 2, 4}), signed_layer({2, 4, 0}, {3.5, 4, 0}, {0, 1.5}, lines), signed_layer({3.5, 4, 0}, {3.5, 0, 0}, {0, 2, 4}, lines)},
                "the corners at its two ends overlap", 1);
    }
    // 135 degrees: the top meeting a chamfer from (5, 4) to (6, 3).
    {
        const Geometry2D chamfer{{{{{2, 0, 0}, {6, 0, 0}, {6, 3, 0}, {5, 4, 0}, {2, 4, 0}}, 0},
                                  {{{0, 0, 0}, {2, 0, 0}, {2, 4, 0}, {5, 4, 0}, {6, 3, 0}, {6, 0, 0}, {8, 0, 0}, {8, 8, 0}, {0, 8, 0}}, 1}},
                                 {}};
        refused(chamfer, {signed_layer({2, 4, 0}, {5, 4, 0}, {0, 1, 3}, lines), signed_layer({5, 4, 0}, {6, 3, 0}, {0, std::sqrt(2.0)}, lines)},
                "a corner's angle outside [60, 120] degrees", 1);
    }
    // An interface across the top inside its oblique (84 degree) corner: the oxide split at x = 5.2
    // above a block whose right side leans out to (6.4, 0).
    {
        const Geometry2D leaning{{{{{2, 0, 0}, {6.4, 0, 0}, {6, 4, 0}, {5.2, 4, 0}, {2, 4, 0}}, 0},
                                  {{{0, 0, 0}, {2, 0, 0}, {2, 4, 0}, {5.2, 4, 0}, {5.2, 8, 0}, {0, 8, 0}}, 1},
                                  {{{5.2, 4, 0}, {6, 4, 0}, {6.4, 0, 0}, {8, 0, 0}, {8, 8, 0}, {5.2, 8, 0}}, 2}},
                                 {}};
        refused(leaning,
                {top({0, 1, 2, 4}), left({0, 2, 4}), signed_layer({6, 4, 0}, {6.4, 0, 0}, {0, 2, std::hypot(0.4, 4.0)}, lines)},
                "an interface across the layer inside a corner", 0);
    }
    // A patch end inside a corner: layers on a square's bottom and right boundaries (lines from the
    // boundary inward), a patch on the right side from inside the corner.
    {
        const Geometry2D square{{rect(0, 0, 4, 4, 0)}, {{"right", {4, 0.5, 0}, {4, 4, 0}}}};
        refused(square,
                {signed_layer({0, 0, 0}, {4, 0, 0}, {0, 2, 4}, {0, 0.5, 1}), signed_layer({4, 0, 0}, {4, 4, 0}, {0, 2, 4}, {0, 0.5, 1})},
                "a patch end inside a corner", 1);
    }
}

TEST_CASE("signed layer: a strip with a layer straddling its interface is the tensor mesh (UC6)", "[mesh][generate][layers]") {
    // Two regions of the unit square meeting at y = 0.5; one layer on the interface covering the
    // whole square: its lines reach both boundaries, its end columns are the sides. The mesh is the
    // lattice alone, its dual the tensor box method, the interface exact.
    const double s = 1.0 / 64.0;
    const Geometry2D g{{rect(0, 0, 1, 0.5, 0), rect(0, 0.5, 1, 1, 1)}, {{"bottom", {0, 0, 0}, {1, 0, 0}}}};
    MeshOptions2D o;
    o.max_size = 0.1;
    o.interface_resolution = s;
    BoundaryLayer layer;
    layer.a = {0, 0.5, 0};
    layer.b = {1, 0.5, 0};
    layer.along = {0, 0.25, 0.5, 0.75, 1};
    layer.lines = {-0.5, -0.25, -0.1, 0.1, 0.25, 0.5};
    o.layers = {layer};
    const auto m = generate_mesh(g, o);
    INFO((m ? std::string("ok") : m.error().message));
    REQUIRE(m.has_value());
    const Lattice L(layer.a, layer.b, layer.along, {-0.5, -0.25, -0.1, -0.5 * s, 0.5 * s, 0.1, 0.25, 0.5});
    CHECK(m->mesh.node_count() == L.along.size() * L.lines.size());
    CHECK(m->report.interface_nodes == 0);  // no mirrored pairs: the lattice straddles the interface
    CHECK(check_tensor(*m, L, true, true) == L.lines.size());
    REQUIRE(m->certificate.has_value());
    INFO("flat " << m->certificate->report.max_flat_deviation << " corner " << m->certificate->report.max_corner_deviation);
    CHECK(m->certificate->report.max_flat_deviation <= 1e-9);
    CHECK(close(total(m->mesh.boundary()[0].areas), 1.0, 1e-12));
    // Determinism.
    const auto again = generate_mesh(g, o);
    REQUIRE(again.has_value());
    CHECK(again->simplices.points == m->simplices.points);
    CHECK(again->simplices.cells == m->simplices.cells);
}

TEST_CASE("signed layer: collinear interfaces on one line share a straddle (UC6, slice 4)", "[mesh][generate][layers]") {
    // Three regions meet at (0.5, 0.5): two below y = 0.5, side by side, and one above. The line y = 0.5
    // is two interface chains (0 | 2 and 1 | 2), one straddle pair serving both; the interface
    // x = 0.5 below it crosses the layer and is straddled by columns. The mesh is the lattice alone
    // (until slice 4 the second chain's straddle was refused as overlapping the first's).
    const double s = 1.0 / 64.0;
    const Geometry2D g{{rect(0, 0, 0.5, 0.5, 0), rect(0.5, 0, 1, 0.5, 1), {{{0, 0.5, 0}, {0.5, 0.5, 0}, {1, 0.5, 0}, {1, 1, 0}, {0, 1, 0}}, 2}},
                       {{"top", {0, 1, 0}, {1, 1, 0}}}};
    MeshOptions2D o;
    o.max_size = 0.1;
    o.interface_resolution = s;
    BoundaryLayer layer;
    layer.a = {0, 0.5, 0};
    layer.b = {1, 0.5, 0};
    layer.along = {0, 0.25, 0.75, 1};
    layer.lines = {-0.5, -0.25, -0.1, 0.1, 0.25, 0.5};
    o.layers = {layer};
    const auto m = generate_mesh(g, o);
    INFO((m ? std::string("ok") : m.error().message));
    REQUIRE(m.has_value());
    const Lattice L(layer.a, layer.b, {0, 0.25, 0.5 - 0.5 * s, 0.5 + 0.5 * s, 0.75, 1},
                    {-0.5, -0.25, -0.1, -0.5 * s, 0.5 * s, 0.1, 0.25, 0.5});
    CHECK(m->mesh.node_count() == L.along.size() * L.lines.size());
    CHECK(check_tensor(*m, L, true, true) == L.lines.size());
    REQUIRE(m->certificate.has_value());
    INFO("flat " << m->certificate->report.max_flat_deviation << " corner " << m->certificate->report.max_corner_deviation);
    CHECK(m->certificate->report.max_flat_deviation <= 1e-9);
    CHECK(m->certificate->report.max_corner_deviation <= 1e-12);  // the junction's faces lie on the lines
}

TEST_CASE("signed layer: a three-layer gate stack under one lattice (UC6)", "[mesh][generate][layers]") {
    // Si (y < 10), a 0.5-thick oxide, a 1.5-thick high-k layer to the top boundary at y = 12; the
    // layer's base the Si interface, its lines from 3 below it to the boundary. Both interfaces are
    // straddled at +-s/2, the top boundary is a line, the sides are columns; below the lattice the
    // refinement fills Si.
    const double s = 0.1;
    const Geometry2D g{{rect(0, 0, 20, 10, 0), rect(0, 10, 20, 10.5, 1), rect(0, 10.5, 20, 12, 2)},
                       {{"gate", {0, 12, 0}, {20, 12, 0}}, {"body", {0, 0, 0}, {20, 0, 0}}}};
    MeshOptions2D o;
    o.max_size = 2.0;
    o.interface_resolution = s;
    BoundaryLayer layer;
    layer.a = {0, 10, 0};
    layer.b = {20, 10, 0};
    for (int i = 0; i <= 8; ++i) layer.along.push_back(2.5 * i);
    layer.lines = {-3, -2, -1.2, -0.6, -0.3, -0.15, 0.25, 0.8, 1.2, 2.0};
    o.layers = {layer};
    const auto m = generate_mesh(g, o);
    INFO((m ? std::string("ok") : m.error().message));
    REQUIRE(m.has_value());
    std::vector<double> lines = layer.lines;
    for (const double x : {-0.5 * s, 0.5 * s, 0.5 - 0.5 * s, 0.5 + 0.5 * s}) lines.push_back(x);
    const Lattice L(layer.a, layer.b, layer.along, lines);
    // Every row but the lowest (which borders the refined Si) is the tensor box method.
    CHECK(check_tensor(*m, L, false, true) == L.lines.size() - 2);  // the two lowest rows border the refinement
    std::size_t inside = 0;
    for (const Point& p : m->simplices.points) inside += p[1] > 7.0 ? 1 : 0;
    CHECK(inside == L.along.size() * (L.lines.size() - 1));  // the lattice's nodes only above its lowest row
    REQUIRE(m->certificate.has_value());
    INFO("flat " << m->certificate->report.max_flat_deviation << " corner " << m->certificate->report.max_corner_deviation);
    CHECK(m->certificate->report.max_flat_deviation <= 1e-9);
    CHECK(m->certificate->report.max_corner_deviation <= s);
    CHECK(close(total(m->mesh.boundary()[0].areas), 20.0, 1e-12));
    CHECK(close(total(m->mesh.boundary()[1].areas), 20.0, 1e-12));
}

TEST_CASE("signed layer: an interface across a boundary layer is straddled by columns (UC6)", "[mesh][generate][layers]") {
    // Two regions side by side meeting at x = 0.5; a layer on the bottom boundary with lines 0 to
    // 0.2. The interface crosses it: columns at 0.5 +- s/2, no junction straddle, and above the
    // lattice the mirrored pairs resume 1.5 s from its top line. A patch end inside the lattice
    // (the left side's) is a lattice node exactly.
    const double s = 1.0 / 64.0;
    const Geometry2D g{{rect(0, 0, 0.5, 1, 0), rect(0.5, 0, 1, 1, 1)}, {{"left", {0, 0, 0}, {0, 1, 0}}}};
    MeshOptions2D o;
    o.max_size = 0.1;
    o.interface_resolution = s;
    BoundaryLayer layer;
    layer.a = {0, 0, 0};
    layer.b = {1, 0, 0};
    layer.along = {0, 0.125, 0.25, 0.375, 0.625, 0.75, 0.875, 1};
    layer.lines = {0, 0.05, 0.1, 0.2};
    o.layers = {layer};
    const auto m = generate_mesh(g, o);
    INFO((m ? std::string("ok") : m.error().message));
    REQUIRE(m.has_value());
    const Lattice L(layer.a, layer.b, {0, 0.125, 0.25, 0.375, 0.5 - 0.5 * s, 0.5 + 0.5 * s, 0.625, 0.75, 0.875, 1}, layer.lines);
    CHECK(check_tensor(*m, L, true, false) == 2);  // the top two rows border the refinement
    for (const Point& p : m->simplices.points) {
        if (p[1] >= 0.2) continue;
        INFO(p[0] << ", " << p[1]);
        CHECK((std::ranges::find(L.lines, p[1]) != L.lines.end() && std::ranges::find(L.along, p[0]) != L.along.end()));
    }
    // The first pair above the lattice, 1.5 s from its top line.
    CHECK(node_at(*m, {0.5 - 0.5 * s, 0.2 + 1.5 * s, 0}).has_value());
    CHECK(node_at(*m, {0.5 + 0.5 * s, 0.2 + 1.5 * s, 0}).has_value());
    REQUIRE(m->certificate.has_value());
    INFO("flat " << m->certificate->report.max_flat_deviation << " corner " << m->certificate->report.max_corner_deviation);
    CHECK(m->certificate->report.max_flat_deviation <= 1e-9);  // measured 0: the pairs resume mirrored
    CHECK(m->certificate->report.max_corner_deviation <= s);
    CHECK(close(total(m->mesh.boundary()[0].areas), 1.0, 1e-12));
}

TEST_CASE("signed layer: refusals (UC2, UD5)", "[mesh][generate][layers]") {
    const double s = 1.0 / 64.0;
    const Geometry2D strip{{rect(0, 0, 1, 0.5, 0), rect(0, 0.5, 1, 1, 1)}, {}};
    const Geometry2D halves{{rect(0, 0, 0.5, 1, 0), rect(0.5, 0, 1, 1, 1)}, {}};
    const auto refused = [&](const Geometry2D& g, const BoundaryLayer& layer, const std::string& text) {
        MeshOptions2D o;
        o.max_size = 0.1;
        o.interface_resolution = s;
        o.layers = {layer};
        const auto m = generate_mesh(g, o);
        INFO(text);
        REQUIRE_FALSE(m.has_value());
        CHECK(m.error().code == ErrorCode::invalid_input);
        if (text != "too short for its pairs") {  // a chain's rule, not the layer's: no layer index
            REQUIRE(m.error().context.has_value());
            CHECK(m.error().context->index == std::optional<std::size_t>{0});
        }
        INFO(m.error().message);
        CHECK(m.error().message.find(text) != std::string::npos);
    };
    const auto strip_layer = [](std::vector<double> lines, std::vector<double> along = {0, 0.5, 1}) {
        BoundaryLayer l;
        l.a = {0, 0.5, 0};
        l.b = {1, 0.5, 0};
        l.along = std::move(along);
        l.lines = std::move(lines);
        return l;
    };
    const auto bottom_layer = [](Point a, Point b, std::vector<double> along, std::vector<double> lines) {
        BoundaryLayer l;
        l.a = a;
        l.b = b;
        l.along = std::move(along);
        l.lines = std::move(lines);
        return l;
    };
    {
        BoundaryLayer both = strip_layer({-0.1, 0.1});
        both.offsets = {0.1};
        refused(strip, both, "not both");
    }
    {
        BoundaryLayer neg = strip_layer({-0.1, 0.1});
        neg.straddle = -1.0;
        refused(strip, neg, "straddle not finite");
    }
    {
        BoundaryLayer off = strip_layer({-0.1, 0.1});
        off.a = {0, 0.3, 0};
        off.b = {1, 0.3, 0};
        refused(strip, off, "the base is not covered");
    }
    refused(strip, strip_layer({0.1}), "at least two lines");
    refused(strip, strip_layer({0.1, -0.1}), "strictly increasing");
    refused(strip, strip_layer({-0.1, 0.1}, {0, 0.6, 0.5, 1}), "along not strictly increasing");
    refused(strip, strip_layer({-0.1, 0.1}, {0, 0.9}), "along must run");
    refused(strip, strip_layer({-0.1, 0.0, 0.1}), "a line inside the straddle");
    refused(strip, strip_layer({-0.005, 0.1}), "needs a line beyond its straddle");
    // A layer on the bottom boundary reaching to 0.47: the interface at 0.5 lies inside its band.
    refused(strip, bottom_layer({0, 0, 0}, {1, 0, 0}, {0, 1}, {0, 0.47}), "protected band");
    // Lines past the top boundary.
    refused(strip, strip_layer({-0.1, 0.5, 0.6}), "the lattice leaves the domain");
    // A boundary inside the lattice that is not one of its lines: the bottom at -0.5, lines to -0.6.
    refused(strip, strip_layer({-0.6, 0.1}), "not one of its lines");
    // A crossing interface: a column inside its straddle, or too near the first column.
    refused(halves, bottom_layer({0, 0, 0}, {1, 0, 0}, {0, 0.5, 1}, {0, 0.1}), "a column inside the straddle");
    refused(halves, bottom_layer({0.495, 0, 0}, {1, 0, 0}, {0, 0.505}, {0, 0.1}), "first or last column");
    // The uncovered part of the crossing interface too short: the lattice stops 2 s below its end.
    refused(halves, bottom_layer({0, 0, 0}, {1, 0, 0}, {0, 0.25, 0.75, 1}, {0, 1 - 2 * s}), "too short for its pairs");
    // A domain vertex on a line between columns: the polygons split at x = 0.3, one region.
    {
        const Geometry2D split{{rect(0, 0, 0.3, 1, 0), rect(0.3, 0, 1, 1, 0)}, {}};
        MeshOptions2D o;
        o.max_size = 0.1;
        o.layers = {bottom_layer({0, 0, 0}, {1, 0, 0}, {0, 0.5, 1}, {0, 0.1})};
        const auto m = generate_mesh(split, o);
        REQUIRE_FALSE(m.has_value());
        INFO(m.error().message);
        CHECK(m.error().message.find("not on one of its columns") != std::string::npos);
        o.layers = {bottom_layer({0, 0, 0}, {1, 0, 0}, {0, 0.3, 1}, {0, 0.1})};
        CHECK(generate_mesh(split, o).has_value());
    }
    // An oblique interface through the lattice.
    {
        const Geometry2D slant{{{{{0, 0, 0}, {0.4, 0, 0}, {0.6, 1, 0}, {0, 1, 0}}, 0}, {{{0.4, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0.6, 1, 0}}, 1}}, {}};
        refused(slant, bottom_layer({0, 0, 0}, {1, 0, 0}, {0, 0.25, 0.75, 1}, {0, 0.1}), "neither parallel nor perpendicular");
    }
    // A patch end inside the lattice between its lines (on the left side, at y = 0.05).
    {
        const Geometry2D g{{rect(0, 0, 1, 1, 0)}, {{"left", {0, 0.05, 0}, {0, 1, 0}}}};
        MeshOptions2D o;
        o.max_size = 0.1;
        o.layers = {bottom_layer({0, 0, 0}, {1, 0, 0}, {0, 1}, {0, 0.1})};
        const auto m = generate_mesh(g, o);
        REQUIRE_FALSE(m.has_value());
        INFO(m.error().message);
        CHECK(m.error().message.find("patch end inside the lattice") != std::string::npos);
    }
}
