// Gates for patches at interface traces in the complex generator (Unit 53, slice 4; DECISIONS.md
// FD8, FD10; FC7, and FC2 for the FD8 rules): a patch edge on a trace, crossing a trace and ending
// on one; a patch spanning two regions; P6.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "NiTCAD/mesh/complex.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/interfaces.hpp"
#include "NiTCAD/mesh/mesh.hpp"

using namespace NiTCAD::mesh;
using NiTCAD::base::ErrorCode;

namespace {

Point P2(double u, double v) { return {u, v, 0.0}; }

// A conforming grid of rectangle cells over the breaks, along x, the region of each cell and
// section from `region(u_mid, v_mid, section)`.
ExtrudedGeometry grid(std::vector<double> ub, std::vector<double> vb, std::vector<double> stations,
                      const std::function<std::int32_t(double, double, std::size_t)>& region) {
    ExtrudedGeometry x;
    x.axis = 0;
    x.stations = std::move(stations);
    for (std::size_t i = 0; i + 1 < ub.size(); ++i) {
        for (std::size_t j = 0; j + 1 < vb.size(); ++j) {
            x.cells.push_back({"c" + std::to_string(i) + "_" + std::to_string(j),
                               {P2(ub[i], vb[j]), P2(ub[i + 1], vb[j]), P2(ub[i + 1], vb[j + 1]), P2(ub[i], vb[j + 1])}});
            std::vector<std::int32_t> row;
            for (std::size_t k = 0; k + 1 < x.stations.size(); ++k) row.push_back(region(0.5 * (ub[i] + ub[i + 1]), 0.5 * (vb[j] + vb[j + 1]), k));
            x.regions.push_back(row);
        }
    }
    return x;
}

// Two slabs along x, 2 wide (y) and 2 high (z): region 0 below z = 1, region 1 above. The
// interface's traces run along z = 1 on the side faces y = 0 and y = 2 and on the end faces.
ExtrudedGeometry slabs(std::vector<double> stations) {
    return grid({0.0, 2.0}, {0.0, 1.0, 2.0}, std::move(stations), [](double, double v, std::size_t) { return v < 1.0 ? 0 : 1; });
}

MeshOptions3D options(double h, double s) {
    MeshOptions3D o;
    o.max_size = h;
    o.interface_resolution = s;
    return o;
}

GeneratedMesh mesh_of(const ExtrudedGeometry& x, const MeshOptions3D& o) {
    const auto g = extrude(x);
    if (!g) FAIL(g.error().message);
    auto m = generate_mesh(*g, o);
    if (!m) FAIL(m.error().message);
    return std::move(*m);
}

void refused(const ExtrudedGeometry& x, const MeshOptions3D& o, const std::string& text) {
    const auto g = extrude(x);
    REQUIRE(g);
    const auto m = generate_mesh(*g, o);
    REQUIRE_FALSE(m);
    CHECK(m.error().code == ErrorCode::invalid_input);
    INFO(m.error().message);
    CHECK(m.error().message.find(text) != std::string::npos);
    CHECK((m.error().context.has_value() && m.error().context->index.has_value()));
}

double patch_area(const Mesh& m, const std::string& name) {
    const BoundaryPatch* p = m.find_boundary(name);
    REQUIRE(p != nullptr);
    double a = 0.0;
    for (const double x : p->areas) a += x;
    return a;
}

std::size_t patch_index(const GeneratedMesh& m, const std::string& name) {
    const auto it = std::ranges::find(m.simplices.patch_names, name);
    REQUIRE(it != m.simplices.patch_names.end());
    return static_cast<std::size_t>(it - m.simplices.patch_names.begin());
}

// The FC7 checks on one mesh (P4 and P6 certified; each patch's area error within s x its zone
// length, and the mesh's area equal to the reported one); the report for the record.
void check_certified(const GeneratedMesh& m, double s, const std::vector<std::pair<std::string, double>>& areas, bool straddled) {
    REQUIRE(m.certificate.has_value());
    CHECK(certifies(*m.certificate, m.mesh, m.node_region));
    const InterfaceReport& r = m.certificate->report;
    CHECK(r.max_flat_deviation <= 1.0);
    CHECK(r.max_corner_deviation <= s);
    REQUIRE(m.certificate->patches.has_value());
    const PatchBoundaryReport& p = *m.certificate->patches;
    CHECK((p.cross_pieces > 0) == straddled);
    CHECK(m.simplices.facet_corner_patch.empty() == !straddled);
    CHECK(p.max_flat_deviation <= 1.0);
    CHECK(p.max_corner_deviation <= s);
    REQUIRE(p.area_error.size() == m.simplices.patch_names.size());
    REQUIRE(p.zone_length.size() == m.simplices.patch_names.size());
    for (const auto& [name, geometric] : areas) {
        const std::size_t k = patch_index(m, name);
        const double mesh_area = patch_area(m.mesh, name);
        CHECK(std::abs(mesh_area - geometric - p.area_error[k]) <= 1e-12 * geometric);
        CHECK(std::abs(p.area_error[k]) <= s * p.zone_length[k] + 1e-12 * geometric);
        UNSCOPED_INFO("patch " << name << ": area " << mesh_area << " of " << geometric << ", error " << p.area_error[k]
                               << ", zone length " << p.zone_length[k]);
    }
    UNSCOPED_INFO("nodes " << m.report.nodes << ", P4 flat " << r.max_flat_deviation << ", P4 corner/s " << r.max_corner_deviation / s
                           << ", P6 pieces " << p.cross_pieces << " (" << p.corner_pieces << " in zones), P6 flat " << p.max_flat_deviation
                           << ", P6 corner/s " << p.max_corner_deviation / s << ", slivers " << m.report.slivers);
}

}  // namespace

TEST_CASE("complex patches: a patch edge on an interface trace (FC7, FD8 a)", "[mesh][complex_patches]") {
    const double s = 0.1;
    auto x = slabs({0.0, 1.0});
    // The side face y = 0 below the interface: its edge along z = 1 is the trace.
    x.patches = {{"low", {-1.0, 0.0, -1.0}, {3.0, 0.0, 1.0}}};
    const auto m = mesh_of(x, options(0.5, s));
    check_certified(m, s, {{"low", 1.0}}, true);
    // Measured 2026-10-30 (Release): 350 nodes, 28 pieces (13 in zones), every deviation 0 and
    // the area exact: the trace pairs are mirror images up to the ends. Regression bounds.
    CHECK(m.certificate->patches->max_corner_deviation <= 0.05 * s);
    CHECK(std::abs(m.certificate->patches->area_error[patch_index(m, "low")]) <= 1e-12);
    const auto again = mesh_of(x, options(0.5, s));
    CHECK(std::ranges::equal(again.mesh.points(), m.mesh.points()));
    CHECK(again.simplices.facet_corner_patch == m.simplices.facet_corner_patch);
    WARN("patch edge on a trace");
}

TEST_CASE("complex patches: a patch edge crossing a trace, the patch over two regions (FC7, FD8 b, c)", "[mesh][complex_patches]") {
    const double s = 0.1;
    auto x = slabs({0.0, 1.0, 2.0});
    // The side face y = 0 for x in [0, 1], both regions: its edge x = 1 crosses the trace at (1, 0, 1).
    x.patches = {{"src", {-1.0, 0.0, -1.0}, {1.0, 0.0, 3.0}}};
    const auto m = mesh_of(x, options(0.5, s));
    check_certified(m, s, {{"src", 2.0}}, false);
    // Measured: 536 nodes; the patch edge stays constrained, so no triangle straddles a patch
    // boundary and the area is exact; P4 corner deviation 0 (the crossing nodes are a pair).
    CHECK(std::abs(m.certificate->patches->area_error[patch_index(m, "src")]) <= 1e-12);
    CHECK(m.certificate->report.max_corner_deviation <= 0.05 * s);
    WARN("patch edge crossing a trace");
}

TEST_CASE("complex patches: a patch edge ending on a trace (FC7, FD8 a, b)", "[mesh][complex_patches]") {
    const double s = 0.1;
    auto x = slabs({0.0, 1.0, 2.0});
    // The side face y = 0 for x in [0, 1] below the interface: its edge x = 1 ends on the trace at
    // (1, 0, 1), and its edge along z = 1 lies on the trace.
    x.patches = {{"corner", {-1.0, 0.0, -1.0}, {1.0, 0.0, 1.0}}};
    const auto m = mesh_of(x, options(0.5, s));
    check_certified(m, s, {{"corner", 1.0}}, true);
    // Measured: 539 nodes, 26 pieces (11 in zones), P6 deviations 0; P4 corner deviation 0.25 s
    // and the area 0.00125 short (zone length 0.6); likely the patch edge's end node, s / 2 from
    // the trace with no mirror (not traced). Regression bounds.
    CHECK(m.certificate->report.max_corner_deviation <= 0.28 * s);
    CHECK(m.certificate->patches->max_corner_deviation <= 0.05 * s);
    CHECK(std::abs(m.certificate->patches->area_error[patch_index(m, "corner")]) <= 0.0014);
    WARN("patch edge ending on a trace");
}

TEST_CASE("complex patches: a patch corner where two traces cross (Unit 53 slice 5)", "[mesh][complex_patches]") {
    const double s = 0.1;
    // A block (y in [2, 5]) on a slab (y in [0, 2]); the block is region 1 except for x in [1, 2], where it is
    // region 2. On the side face z = 0 the patch is the rectangle x in [1, 2], y in [2, 5]: its edges x = 1 and
    // x = 2 lie on the traces of the block's two interfaces and its edge y = 2 on the trace of the slab's, so it
    // turns a corner at each point where two traces cross. A corner of a straddled run is a P6 zone (it
    // was not, and the corner pieces failed P6 by s / 4).
    const auto x = [&] {
        auto g = grid({0.0, 2.0, 5.0}, {0.0, 3.0}, {0.0, 1.0, 2.0, 3.0},
                      [](double u, double, std::size_t k) { return u < 2.0 ? 0 : k == 1 ? 2 : 1; });
        g.patches = {{"gate", {1.0, 2.0, 0.0}, {2.0, 5.0, 0.0}}};
        return g;
    }();
    const auto m = mesh_of(x, options(0.5, s));
    check_certified(m, s, {{"gate", 3.0}}, true);
    CHECK(m.certificate->patches->corner_pieces > 0);
    // The patch's four corners: two on the slab's trace, two where the block's traces meet the domain edge.
    CHECK(m.certificate->patches->max_corner_deviation <= s);
    WARN("patch corner at a trace crossing");
}

TEST_CASE("complex patches: without a patch edge at an interface, no P6 (FC1)", "[mesh][complex_patches]") {
    const auto m = mesh_of(slabs({0.0, 1.0}), options(0.5, 0.1));
    REQUIRE(m.certificate.has_value());
    CHECK_FALSE(m.certificate->patches.has_value());
    CHECK(m.simplices.facet_corner_patch.empty());
}

TEST_CASE("complex patches: the FD8 rules are refused (FC2, FC7)", "[mesh][complex_patches]") {
    const double s = 0.1;
    const MeshOptions3D o = options(0.5, s);
    SECTION("a patch edge meeting a trace at under 60 degrees") {
        // On the end face x = 0: a cell edge from (0.25, 0) to (1.75, 2) meets the trace z = 1 at
        // (1, 1), 53 degrees from it (coordinates exact, so the edge is one chain); the cells left of it
        // take a patch.
        ExtrudedGeometry x;
        x.axis = 0;
        x.stations = {0.0, 1.0};
        x.cells = {{"ll", {P2(0, 0), P2(0.25, 0), P2(1, 1), P2(0, 1)}},
                   {"lr", {P2(0.25, 0), P2(2, 0), P2(2, 1), P2(1, 1)}},
                   {"ul", {P2(0, 1), P2(1, 1), P2(1.75, 2), P2(0, 2)}},
                   {"ur", {P2(1, 1), P2(2, 1), P2(2, 2), P2(1.75, 2)}}};
        x.regions = {{0}, {0}, {1}, {1}};
        x.patches = {{"left", {-1.0, -1.0, -1.0}, {0.0, 1.75, 3.0}}};
        refused(x, o, "a patch edge meets an interface trace at under 60 degrees");
    }
    SECTION("a patch-edge vertex within 4 s of where it meets a trace") {
        // On the end face x = 0: the patch edge y = 1 has a vertex at z = 1.2, 0.2 from the trace.
        ExtrudedGeometry x;
        x.axis = 0;
        x.stations = {0.0, 1.0};
        x.cells = {{"ll", {P2(0, 0), P2(1, 0), P2(1, 1), P2(0, 1)}},
                   {"lr", {P2(1, 0), P2(2, 0), P2(2, 1), P2(1, 1)}},
                   {"ul", {P2(0, 1), P2(1, 1), P2(1, 1.2), P2(1, 2), P2(0, 2)}},
                   {"ur", {P2(1, 1), P2(2, 1), P2(2, 2), P2(1, 2), P2(1, 1.2)}}};
        x.regions = {{0}, {0}, {1}, {1}};
        x.patches = {{"left", {-1.0, -1.0, -1.0}, {0.0, 1.0, 3.0}}};
        refused(x, o, "within 4 s of where a patch edge meets an interface trace");
    }
    SECTION("a patch edge within 4 s of a trace it does not meet") {
        // On the end face x = 0: a cell 0.2 above the trace, away from the domain edges, takes a
        // patch; its lower edge runs parallel to the trace.
        ExtrudedGeometry x;
        x.axis = 0;
        x.stations = {0.0, 1.0};
        x.cells = {{"lo", {P2(0, 0), P2(2, 0), P2(2, 1), P2(1.5, 1), P2(0.5, 1), P2(0, 1)}},
                   {"left", {P2(0, 1), P2(0.5, 1), P2(0.5, 1.2), P2(0.5, 1.6), P2(0.5, 2), P2(0, 2)}},
                   {"right", {P2(1.5, 1), P2(2, 1), P2(2, 2), P2(1.5, 2), P2(1.5, 1.6), P2(1.5, 1.2)}},
                   {"below", {P2(0.5, 1), P2(1.5, 1), P2(1.5, 1.2), P2(0.5, 1.2)}},
                   {"pad", {P2(0.5, 1.2), P2(1.5, 1.2), P2(1.5, 1.6), P2(0.5, 1.6)}},
                   {"above", {P2(0.5, 1.6), P2(1.5, 1.6), P2(1.5, 2), P2(0.5, 2)}}};
        x.regions = {{0}, {1}, {1}, {1}, {1}, {1}};
        x.patches = {{"pad", {-1.0, 0.5, 1.2}, {0.0, 1.5, 1.6}}};
        refused(x, o, "a patch edge within 4 s of an interface feature it does not meet");
    }
}
