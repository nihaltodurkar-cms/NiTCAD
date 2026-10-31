// Unit 54 slice 2 (DECISIONS.md UD5, UD7; UC6 part, UC2 part): 2D structured layers with signed
// lines on interfaces. The lattice's nodes are bit-identical to their formula, its dual is the
// tensor box method, the result is certified (P4), and every UD5 rule is refused with its message.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
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

}  // namespace

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
