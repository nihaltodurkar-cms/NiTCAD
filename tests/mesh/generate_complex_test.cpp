// Gates for the 3D generator over a polyhedral complex (Unit 53, slice 2; DECISIONS.md FD3-FD5;
// FC3 non-convex single-region meshes, FC4 multi-facet patches, FC5 the linear field).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <set>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/mesh/complex.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/mesh/simplex_mesh.hpp"

using namespace NiTCAD::mesh;
using NiTCAD::base::ErrorCode;

namespace {

constexpr std::int32_t V = ExtrudedGeometry::void_region;

Point P2(double u, double v) { return {u, v, 0.0}; }

// A rectangle cell [u0, u1] x [v0, v1].
ExtrusionCell rect(std::string name, double u0, double v0, double u1, double v1) {
    return {std::move(name), {P2(u0, v0), P2(u1, v0), P2(u1, v1), P2(u0, v1)}};
}

// An L-prism along z: three unit squares in (x, y), the corner cell full height, the arms half.
ExtrudedGeometry l_prism() {
    ExtrudedGeometry x;
    x.axis = 2;
    x.stations = {0.0, 1.0, 2.0};
    x.cells = {rect("a", 0, 0, 1, 1), rect("b", 1, 0, 2, 1), rect("c", 0, 1, 1, 2)};
    x.regions = {{0, 0}, {0, 0}, {0, 0}};
    return x;
}

// A U-channel along x: a slab with two walls, (u, v) = (y, z).
ExtrudedGeometry u_channel() {
    ExtrudedGeometry x;
    x.axis = 0;
    x.stations = {0.0, 2.5};
    x.cells = {rect("floor_l", 0, 0, 1, 1), rect("floor_c", 0, 1, 1, 2), rect("floor_r", 0, 2, 1, 3),
               rect("wall_l", 1, 0, 2, 1), rect("gap", 1, 1, 2, 2), rect("wall_r", 1, 2, 2, 3)};
    x.regions = {{0}, {0}, {0}, {0}, {V}, {0}};
    return x;
}

// A fin on a slab (a T), along x in four sections; `gate` selects the fin's top and sides in the
// second section (FC4).
ExtrudedGeometry t_fin() {
    ExtrudedGeometry x;
    x.axis = 0;
    x.stations = {0.0, 1.0, 2.0, 3.0};
    x.cells = {rect("slab_l", 0, 0, 1, 4), rect("slab_c", 0, 4, 1, 6), rect("slab_r", 0, 6, 1, 10), rect("fin", 1, 4, 4, 6)};
    x.regions = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    x.patches = {{"gate", {1, 1, 4}, {2, 4, 6}}};
    return x;
}

// A slab 3 x 2 cells in cross-section and three sections, the top middle cell void in the middle
// section: a pocket with reflex edges along all three axes.
ExtrudedGeometry pocket() {
    ExtrudedGeometry x;
    x.axis = 1;
    x.stations = {0.0, 1.0, 2.0, 3.0};
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 2; ++j) x.cells.push_back(rect("c" + std::to_string(i) + std::to_string(j), i, j, i + 1, j + 1));
    }
    x.regions.assign(6, {0, 0, 0});
    x.regions[3][1] = V;  // cell c11: u in [1, 2], v in [1, 2]
    return x;
}

// A tapered fin (sidewalls at about 86 degrees to the slab) on a slab, along x.
ExtrudedGeometry tapered_fin() {
    ExtrudedGeometry x;
    x.axis = 0;
    x.stations = {0.0, 3.0};
    x.cells = {rect("slab_l", 0, 0, 1, 2), rect("slab_c", 0, 2, 1, 3), rect("slab_r", 0, 3, 1, 5),
               {"fin", {P2(1, 2), P2(1, 3), P2(3, 2.86), P2(3, 2.14)}}};
    x.regions = {{0}, {0}, {0}, {0}};
    return x;
}

// An L-prism given directly, with L-shaped (non-convex) top and bottom facets.
PolyhedralGeometry l_facets() {
    PolyhedralGeometry g;
    const std::vector<Point> ring{{0, 0, 0}, {2, 0, 0}, {2, 1, 0}, {1, 1, 0}, {1, 2, 0}, {0, 2, 0}};
    for (const Point& p : ring) g.vertices.push_back(p);
    for (const Point& p : ring) g.vertices.push_back({p[0], p[1], 1.0});
    g.facets.push_back({{0, 1, 2, 3, 4, 5}, "bottom"});
    g.facets.push_back({{6, 7, 8, 9, 10, 11}, "top"});
    for (std::size_t k = 0; k < 6; ++k) g.facets.push_back({{k, (k + 1) % 6, (k + 1) % 6 + 6, k + 6}, "side"});
    g.volumes.push_back({{0, 1, 2, 3, 4, 5, 6, 7}, 0});
    return g;
}

double volume_of(const SimplexMesh& s) {
    double v = 0.0;
    for (std::size_t c = 0; c < s.cells.size(); c += 4) {
        const auto& a = s.points[static_cast<std::size_t>(s.cells[c])];
        const auto& b = s.points[static_cast<std::size_t>(s.cells[c + 1])];
        const auto& d = s.points[static_cast<std::size_t>(s.cells[c + 2])];
        const auto& e = s.points[static_cast<std::size_t>(s.cells[c + 3])];
        const Point u{b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        const Point w{d[0] - a[0], d[1] - a[1], d[2] - a[2]};
        const Point z{e[0] - a[0], e[1] - a[1], e[2] - a[2]};
        v += u[0] * (w[1] * z[2] - w[2] * z[1]) - u[1] * (w[0] * z[2] - w[2] * z[0]) + u[2] * (w[0] * z[1] - w[1] * z[0]);
    }
    return v / 6.0;
}

double total(std::span<const double> v) {
    double s = 0.0;
    for (const double x : v) s += x;
    return s;
}

// FC3's checks on one mesh.
void check_valid(const GeneratedMesh& m, double volume, const MeshOptions3D& o) {
    CHECK(m.mesh.dimension() == 3);
    CHECK(m.report.nodes == m.mesh.node_count());
    CHECK(std::abs(volume_of(m.simplices) - volume) <= 1e-12 * volume);
    CHECK(std::abs(total(m.mesh.volumes()) - volume) <= 1e-12 * volume);
    CHECK(m.report.min_dihedral_deg > 0.0);
    if (m.report.exempt_tetrahedra == 0) CHECK(m.report.max_radius_edge <= o.radius_edge * (1.0 + 1e-9));
    // Every boundary triangle Gabriel: no node strictly inside its smallest sphere (P2, rechecked
    // here in double with a margin; make_dual_mesh checks it exactly).
    const auto& pts = m.simplices.points;
    for (std::size_t f = 0; f < m.simplices.facets.size(); f += 3) {
        const Point& a = pts[static_cast<std::size_t>(m.simplices.facets[f])];
        const Point& b = pts[static_cast<std::size_t>(m.simplices.facets[f + 1])];
        const Point& c = pts[static_cast<std::size_t>(m.simplices.facets[f + 2])];
        const Point ab{b[0] - a[0], b[1] - a[1], b[2] - a[2]};
        const Point ac{c[0] - a[0], c[1] - a[1], c[2] - a[2]};
        const Point n{ab[1] * ac[2] - ab[2] * ac[1], ab[2] * ac[0] - ab[0] * ac[2], ab[0] * ac[1] - ab[1] * ac[0]};
        const double nn = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
        const double ab2 = ab[0] * ab[0] + ab[1] * ab[1] + ab[2] * ab[2];
        const double ac2 = ac[0] * ac[0] + ac[1] * ac[1] + ac[2] * ac[2];
        // Circumcentre a + (|ac|^2 (n x ab) + |ab|^2 (ac x n)) / (2 |n|^2).
        const Point nab{n[1] * ab[2] - n[2] * ab[1], n[2] * ab[0] - n[0] * ab[2], n[0] * ab[1] - n[1] * ab[0]};
        const Point acn{ac[1] * n[2] - ac[2] * n[1], ac[2] * n[0] - ac[0] * n[2], ac[0] * n[1] - ac[1] * n[0]};
        Point cc{};
        for (std::size_t k = 0; k < 3; ++k) cc[k] = a[k] + (ac2 * nab[k] + ab2 * acn[k]) / (2.0 * nn);
        const double r2 = (cc[0] - a[0]) * (cc[0] - a[0]) + (cc[1] - a[1]) * (cc[1] - a[1]) + (cc[2] - a[2]) * (cc[2] - a[2]);
        for (const Point& x : pts) {
            const double d2 = (x[0] - cc[0]) * (x[0] - cc[0]) + (x[1] - cc[1]) * (x[1] - cc[1]) + (x[2] - cc[2]) * (x[2] - cc[2]);
            if (d2 < r2 * (1.0 - 1e-9)) {
                FAIL("a boundary triangle is not Gabriel");
            }
        }
    }
}

// FC5: Laplace with V = x_k on every boundary facet normal to axis k and no flux elsewhere, by
// conjugate gradients on the generated dual; the largest deviation from x_k over its range.
double linear_field_error(const GeneratedMesh& m, std::size_t k) {
    const Mesh& mesh = m.mesh;
    const std::size_t n = mesh.node_count();
    const auto& pts = mesh.points();
    std::vector<char> fixed(n, 0);
    for (std::size_t f = 0; f < m.simplices.facets.size(); f += 3) {
        const auto a = static_cast<std::size_t>(m.simplices.facets[f]);
        const auto b = static_cast<std::size_t>(m.simplices.facets[f + 1]);
        const auto c = static_cast<std::size_t>(m.simplices.facets[f + 2]);
        if (pts[a][k] == pts[b][k] && pts[b][k] == pts[c][k]) fixed[a] = fixed[b] = fixed[c] = 1;
    }
    double lo = pts[0][k];
    double hi = pts[0][k];
    for (const Point& p : pts) {
        lo = std::min(lo, p[k]);
        hi = std::max(hi, p[k]);
    }
    std::vector<double> x(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) x[i] = fixed[i] != 0 ? pts[i][k] : 0.5 * (lo + hi);
    const auto apply = [&](const std::vector<double>& y) {
        std::vector<double> out(n, 0.0);
        for (const Edge& e : mesh.edges()) {
            const auto a = static_cast<std::size_t>(e.first);
            const auto b = static_cast<std::size_t>(e.second);
            const double w = e.coupling_area / e.length;
            out[a] += w * (y[a] - y[b]);
            out[b] += w * (y[b] - y[a]);
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (fixed[i] != 0) out[i] = 0.0;
        }
        return out;
    };
    const auto dot = [&](const std::vector<double>& a, const std::vector<double>& b) {
        double s = 0.0;
        for (std::size_t i = 0; i < n; ++i) s += a[i] * b[i];
        return s;
    };
    std::vector<double> r = apply(x);
    for (double& v : r) v = -v;
    std::vector<double> p = r;
    double rr = dot(r, r);
    const double r0 = std::sqrt(rr);
    for (int it = 0; it < 20000 && std::sqrt(rr) > 1e-16 * r0; ++it) {
        const std::vector<double> Ap = apply(p);
        const double alpha = rr / dot(p, Ap);
        for (std::size_t i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }
        const double next = dot(r, r);
        for (std::size_t i = 0; i < n; ++i) p[i] = r[i] + next / rr * p[i];
        rr = next;
    }
    double worst = 0.0;
    for (std::size_t i = 0; i < n; ++i) worst = std::max(worst, std::abs(x[i] - pts[i][k]));
    return worst / (hi - lo);
}

GeneratedMesh mesh_of(const ExtrudedGeometry& x, const MeshOptions3D& o) {
    const auto g = extrude(x);
    REQUIRE(g);
    auto m = generate_mesh(*g, o);
    if (!m) FAIL(m.error().message);
    return std::move(*m);
}

}  // namespace

TEST_CASE("complex generator: non-convex single-region meshes are valid (FC3)", "[mesh][complex_generate]") {
    MeshOptions3D o;
    o.max_size = 0.45;
    SECTION("L-prism") {
        const auto m = mesh_of(l_prism(), o);
        check_valid(m, 6.0, o);
        CHECK(m.mesh.boundary().size() == 1);  // the default patch over every facet
    }
    SECTION("U-channel") { check_valid(mesh_of(u_channel(), o), 2.5 * 5.0, o); }
    SECTION("fin on a slab") {
        o.max_size = 0.8;
        check_valid(mesh_of(t_fin(), o), 3.0 * (10.0 + 6.0), o);
    }
    SECTION("pocket") { check_valid(mesh_of(pocket(), o), 18.0 - 1.0, o); }
    SECTION("tapered fin (86 degrees, heuristic)") {
        o.max_size = 0.5;
        check_valid(mesh_of(tapered_fin(), o), 3.0 * (5.0 + 2.0 * (1.0 + 0.72) / 2.0), o);
    }
    SECTION("a non-convex facet") {
        const auto m = generate_mesh(l_facets(), o);
        REQUIRE(m);
        check_valid(*m, 3.0, o);
        REQUIRE(m->mesh.find_boundary("top") != nullptr);
        CHECK(std::abs(total(m->mesh.find_boundary("top")->areas) - 3.0) <= 3e-12);
        CHECK(std::abs(total(m->mesh.find_boundary("side")->areas) - 8.0) <= 8e-12);
    }
}

TEST_CASE("complex generator: the same input gives the same mesh", "[mesh][complex_generate]") {
    MeshOptions3D o;
    o.max_size = 0.45;
    const auto a = mesh_of(pocket(), o);
    const auto b = mesh_of(pocket(), o);
    CHECK(a.simplices.points == b.simplices.points);
    CHECK(a.simplices.cells == b.simplices.cells);
    CHECK(graph_digest(a.mesh) == graph_digest(b.mesh));
    // Pinned (FC13, FD15): a change of insertion order or of any refinement rule changes it.
    CHECK(graph_digest(a.mesh) == 0x6c05e45d2f27a5dfull);
}

TEST_CASE("complex generator: a patch over three facets has their summed areas (FC4)", "[mesh][complex_generate]") {
    MeshOptions3D o;
    o.max_size = 0.6;
    const auto m = mesh_of(t_fin(), o);
    const BoundaryPatch* gate = m.mesh.find_boundary("gate");
    REQUIRE(gate != nullptr);
    CHECK(std::abs(total(gate->areas) - 8.0) <= 8e-12);
    // Each node's gate area against the 2D duals of the gate's triangles on each of its planes.
    const auto& s = m.simplices;
    const auto gate_index = static_cast<std::int32_t>(std::ranges::find(s.patch_names, "gate") - s.patch_names.begin());
    std::map<NodeId, double> expected;
    for (const auto& [axis, value] : std::vector<std::pair<std::size_t, double>>{{1, 4.0}, {2, 4.0}, {2, 6.0}}) {
        SimplexMesh flat;
        flat.dimension = 2;
        std::map<NodeId, NodeId> local;
        std::vector<NodeId> global;
        for (std::size_t f = 0; f < s.facet_patch.size(); ++f) {
            if (s.facet_patch[f] != gate_index) continue;
            const std::array<NodeId, 3> v{s.facets[3 * f], s.facets[3 * f + 1], s.facets[3 * f + 2]};
            if (!std::ranges::all_of(v, [&](NodeId x) { return s.points[static_cast<std::size_t>(x)][axis] == value; })) continue;
            for (const NodeId x : v) {
                const auto [it, fresh] = local.emplace(x, static_cast<NodeId>(global.size()));
                if (fresh) {
                    global.push_back(x);
                    const Point& p = s.points[static_cast<std::size_t>(x)];
                    const std::size_t b = axis == 1 ? 2 : 1;
                    flat.points.push_back({p[0], p[b], 0.0});
                }
                flat.cells.push_back(it->second);
            }
        }
        // Each triangle counterclockwise in the plane.
        for (std::size_t c = 0; c < flat.cells.size(); c += 3) {
            const Point& a = flat.points[static_cast<std::size_t>(flat.cells[c])];
            const Point& b = flat.points[static_cast<std::size_t>(flat.cells[c + 1])];
            const Point& d = flat.points[static_cast<std::size_t>(flat.cells[c + 2])];
            if ((b[0] - a[0]) * (d[1] - a[1]) - (b[1] - a[1]) * (d[0] - a[0]) < 0.0) std::swap(flat.cells[c + 1], flat.cells[c + 2]);
        }
        REQUIRE_FALSE(flat.cells.empty());
        const auto dual = make_dual_mesh(flat);
        if (!dual) FAIL(dual.error().message);
        for (std::size_t i = 0; i < global.size(); ++i) expected[global[i]] += dual->volumes()[i];
    }
    REQUIRE(expected.size() == gate->nodes.size());
    for (std::size_t i = 0; i < gate->nodes.size(); ++i) {
        const double want = expected.at(gate->nodes[i]);
        CHECK(std::abs(gate->areas[i] - want) <= 1e-12 * 8.0);
    }
}

TEST_CASE("complex generator: a linear potential is exact on non-convex domains (FC5)", "[mesh][complex_generate]") {
    MeshOptions3D o;
    o.max_size = 0.45;
    for (const auto& [name, shape] : std::vector<std::pair<std::string, ExtrudedGeometry>>{{"L", l_prism()}, {"U", u_channel()}, {"pocket", pocket()}}) {
        const auto m = mesh_of(shape, o);
        for (std::size_t k = 0; k < 3; ++k) {
            const double e = linear_field_error(m, k);
            INFO(name << " axis " << k << ": " << m.mesh.node_count() << " nodes, deviation " << e);
            CHECK(e <= 1e-12);
        }
    }
    {
        o.max_size = 0.8;
        const auto m = mesh_of(t_fin(), o);
        for (std::size_t k = 0; k < 3; ++k) {
            const double e = linear_field_error(m, k);
            INFO("T axis " << k << ": deviation " << e);
            CHECK(e <= 1e-12);
        }
    }
    {
        o.max_size = 0.5;
        const auto m = mesh_of(tapered_fin(), o);
        const double e = linear_field_error(m, 0);  // the fin's axis: every other facet is parallel to it
        INFO("tapered fin axis x: deviation " << e);
        CHECK(e <= 1e-12);
    }
}

TEST_CASE("complex generator: refusals", "[mesh][complex_generate]") {
    MeshOptions3D o;
    o.max_size = 0.5;
    SECTION("interfaces (slice 3)") {
        auto x = l_prism();
        x.regions[1] = {1, 1};
        const auto g = extrude(x);
        REQUIRE(g);
        const auto m = generate_mesh(*g, o);
        REQUIRE_FALSE(m);
        CHECK(m.error().code == ErrorCode::invalid_input);
    }
    SECTION("a boundary edge in four boundary facets") {
        ExtrudedGeometry x;
        x.axis = 0;
        x.stations = {0.0, 1.0, 2.0};
        x.cells = {rect("a", 0, 0, 1, 1), rect("b", 1, 0, 2, 1), rect("c", 0, 1, 1, 2), rect("d", 1, 1, 2, 2)};
        x.regions = {{0, 0}, {V, 0}, {V, 0}, {0, 0}};
        const auto g = extrude(x);
        REQUIRE(g);
        const auto m = generate_mesh(*g, o);
        REQUIRE_FALSE(m);
        CHECK(m.error().code == ErrorCode::invalid_input);
        CHECK(m.error().message.find("manifold") != std::string::npos);
    }
    SECTION("options") {
        const auto g = extrude(l_prism());
        REQUIRE(g);
        o.max_size = 0.0;
        CHECK(generate_mesh(*g, o).error().code == ErrorCode::invalid_input);
        o.max_size = 0.5;
        o.radius_edge = 1.0;
        CHECK(generate_mesh(*g, o).error().code == ErrorCode::invalid_input);
        o.radius_edge = 2.0;
        o.max_nodes = 20;
        const auto capped = generate_mesh(*g, o);
        CHECK(capped.error().code == ErrorCode::resource_exhausted);
        CHECK(capped.error().message == "max_nodes reached");  // not the iteration cap (FC13)
    }
}
