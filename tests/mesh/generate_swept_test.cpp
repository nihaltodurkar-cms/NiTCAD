// Gates for the swept mesh of an extrusion (Unit 54, slice 1; DECISIONS.md UD1-UD4, UD8; UC2 part,
// UC3 product exactness, UC4 the linear field, the certificate on swept interfaces).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <map>
#include <string>
#include <utility>
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

ExtrusionCell rect(std::string name, double u0, double v0, double u1, double v1) {
    return {std::move(name), {P2(u0, v0), P2(u1, v0), P2(u1, v1), P2(u0, v1)}};
}

// One region: a rectangle 2 x 1 swept along z over [0, 3].
ExtrudedGeometry box() {
    ExtrudedGeometry x;
    x.axis = 2;
    x.stations = {0.0, 3.0};
    x.cells = {rect("r", 0, 0, 2, 1)};
    x.regions = {{0}};
    return x;
}

// One region: an L of three unit squares swept along y.
ExtrudedGeometry l_shape() {
    ExtrudedGeometry x;
    x.axis = 1;
    x.stations = {0.0, 1.0, 2.5};
    x.cells = {rect("a", 0, 0, 1, 1), rect("b", 1, 0, 2, 1), rect("c", 0, 1, 1, 2)};
    x.regions = {{0, 0}, {0, 0}, {0, 0}};
    return x;
}

// A tapered fin (region 1, sidewalls at about 86 degrees) on a slab (region 0), wrapped by an oxide
// (region 2) out to a rectangular cross-section, swept along x. Every interface meets the outer
// boundary where it is straight (the 2D generator refuses one at a domain corner).
ExtrudedGeometry tapered_fin() {
    ExtrudedGeometry x;
    x.axis = 0;
    x.stations = {0.0, 3.0};
    x.cells = {rect("slab_l", 0, 0, 1, 2), rect("slab_c", 0, 2, 1, 3), rect("slab_r", 0, 3, 1, 5),
               {"fin", {P2(1, 2), P2(1, 3), P2(3, 2.86), P2(3, 2.14)}},
               {"ox", {P2(1, 0), P2(3.5, 0), P2(3.5, 5), P2(1, 5), P2(1, 3), P2(3, 2.86), P2(3, 2.14), P2(1, 2)}}};
    x.regions = {{0}, {0}, {0}, {1}, {2}};
    return x;
}

// A gated fin: a fin (region 1) on a slab (region 0) under an overlayer that is a spacer (2) outside
// the gate section and an oxide (3) inside it; the gate patch on the overlayer's outer faces in the
// middle section only. Stations 0, 1, 2, 3 along x; the two gate edges straddled.
ExtrudedGeometry gated_fin() {
    ExtrudedGeometry x;
    x.axis = 0;
    x.stations = {0.0, 1.0, 2.0, 3.0};
    // (u, v) = (y, z): slab y in [0, 1]; fin y in [1, 2], z in [1, 2]; the overlayer a U over it to y = 2.5.
    x.cells = {{"slab", {P2(0, 0), P2(1, 0), P2(1, 1), P2(1, 2), P2(1, 3), P2(0, 3)}},
               rect("fin", 1, 1, 2, 2),
               {"over", {P2(1, 0), P2(2.5, 0), P2(2.5, 3), P2(1, 3), P2(1, 2), P2(2, 2), P2(2, 1), P2(1, 1)}}};
    x.regions = {{0, 0, 0}, {1, 1, 1}, {2, 3, 2}};
    x.patches = {{"gate", {1.0, 1.0, 0.0}, {2.0, 2.5, 3.0}}};
    return x;
}

SweptOptions options_for(const ExtrudedGeometry& x, double h, std::size_t per_section) {
    SweptOptions o;
    o.section.max_size = h;
    for (std::size_t i = 0; i + 1 < x.stations.size(); ++i) {
        for (std::size_t k = 0; k < per_section; ++k) {
            o.lines.push_back(x.stations[i] + (x.stations[i + 1] - x.stations[i]) * static_cast<double>(k) / static_cast<double>(per_section));
        }
    }
    o.lines.push_back(x.stations.back());
    return o;
}

GeneratedMesh swept(const ExtrudedGeometry& x, const SweptOptions& o) {
    auto m = generate_swept(x, o);
    if (!m) FAIL(m.error().message);
    return std::move(*m);
}

// UC4: Laplace with V = x_k on every boundary facet normal to axis k and no flux elsewhere, by
// conjugate gradients on the dual; the largest deviation from x_k over its range.
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

// UC3: every volume and coupling of the swept dual against the product of the cross-section's 2D
// dual (generated here from the same polygons and options) and the lines' 1D dual.
void check_product(const ExtrudedGeometry& x, const SweptOptions& o, const GeneratedMesh& m, const std::vector<std::int32_t>& column) {
    Geometry2D g2;
    for (std::size_t c = 0; c < x.cells.size(); ++c) g2.polygons.push_back({x.cells[c].polygon, column[c]});
    MeshOptions2D o2 = o.section;
    const auto m2 = generate_mesh(g2, o2);
    REQUIRE(m2);
    const std::size_t n2 = m2->mesh.node_count();
    const std::vector<double>& lines = o.lines;  // no straddles in these cases
    const std::size_t nl = lines.size();
    REQUIRE(m.mesh.node_count() == n2 * nl);
    std::vector<double> d1(nl, 0.0);
    for (std::size_t j = 0; j < nl; ++j) {
        if (j > 0) d1[j] += 0.5 * (lines[j] - lines[j - 1]);
        if (j + 1 < nl) d1[j] += 0.5 * (lines[j + 1] - lines[j]);
    }
    const auto& vol2 = m2->mesh.volumes();
    const auto& vol3 = m.mesh.volumes();
    for (std::size_t j = 0; j < nl; ++j) {
        for (std::size_t p = 0; p < n2; ++p) {
            const double want = vol2[p] * d1[j];
            if (!(std::abs(vol3[j * n2 + p] - want) <= 1e-12 * want)) FAIL("node volume " << vol3[j * n2 + p] << " against " << want);
        }
    }
    std::map<std::pair<NodeId, NodeId>, double> c2;
    for (const Edge& e : m2->mesh.edges()) c2[{std::min(e.first, e.second), std::max(e.first, e.second)}] = e.coupling_area;
    double scale = 0.0;
    for (const Edge& e : m.mesh.edges()) scale = std::max(scale, e.coupling_area);
    std::size_t in_plane = 0;
    std::size_t axial = 0;
    for (const Edge& e : m.mesh.edges()) {
        const auto a = static_cast<std::size_t>(e.first);
        const auto b = static_cast<std::size_t>(e.second);
        const std::size_t ja = a / n2, jb = b / n2, pa = a % n2, pb = b % n2;
        if (ja == jb) {
            const auto it = c2.find({static_cast<NodeId>(std::min(pa, pb)), static_cast<NodeId>(std::max(pa, pb))});
            const double want = it == c2.end() ? 0.0 : it->second * d1[ja];
            if (!(std::abs(e.coupling_area - want) <= 1e-12 * std::max(want, scale * 1e-3))) FAIL("in-plane coupling " << e.coupling_area << " against " << want);
            ++in_plane;
        } else if (pa == pb) {
            const double want = vol2[pa];
            if (!(std::abs(e.coupling_area - want) <= 1e-12 * want)) FAIL("axial coupling " << e.coupling_area << " against " << want);
            ++axial;
        } else {
            // A diagonal inside a prism: a tie, zero coupling (dropped or within rounding).
            if (!(std::abs(e.coupling_area) <= 1e-12 * scale)) FAIL("a prism diagonal has coupling " << e.coupling_area);
        }
    }
    CHECK(axial == n2 * (nl - 1));
    CHECK(in_plane >= c2.size() * nl);
}

}  // namespace

TEST_CASE("swept mesh: the dual is the product of the 2D and 1D duals (UC3)", "[mesh][swept]") {
    SECTION("rectangle") {
        const auto x = box();
        const auto o = options_for(x, 0.3, 7);
        const auto m = swept(x, o);
        CHECK_FALSE(m.certificate.has_value());
        check_product(x, o, m, {0});
    }
    SECTION("L-shape") {
        const auto x = l_shape();
        auto o = options_for(x, 0.35, 3);
        o.lines = {0.0, 0.25, 1.0, 1.1, 2.0, 2.5};  // uneven, station 1 a line
        const auto m = swept(x, o);
        check_product(x, o, m, {0, 0, 0});
    }
    SECTION("tapered fin on a slab (two regions)") {
        const auto x = tapered_fin();
        auto o = options_for(x, 0.4, 5);
        o.section.interface_resolution = 0.05;
        const auto m = swept(x, o);
        check_product(x, o, m, {0, 0, 0, 1, 2});
        REQUIRE(m.certificate.has_value());
        CHECK(m.certificate->report.max_flat_deviation <= 1.0);
        CHECK(m.certificate->report.max_corner_deviation <= 0.05);
    }
    SECTION("two runs are bit-identical") {
        const auto x = tapered_fin();
        auto o = options_for(x, 0.4, 5);
        o.section.interface_resolution = 0.05;
        const auto a = swept(x, o);
        const auto b = swept(x, o);
        CHECK(graph_digest(a.mesh) == graph_digest(b.mesh));
        CHECK(region_digest(a.node_region) == region_digest(b.node_region));
        CHECK(a.generation_digest == b.generation_digest);
        CHECK(a.generation_digest != 0);
    }
}

TEST_CASE("swept mesh: a linear potential is exact (UC4)", "[mesh][swept]") {
    for (const auto& [name, x, h, s] : std::vector<std::tuple<std::string, ExtrudedGeometry, double, double>>{
             {"rectangle", box(), 0.3, 0.0}, {"L", l_shape(), 0.35, 0.0}, {"tapered fin", tapered_fin(), 0.4, 0.05}}) {
        auto o = options_for(x, h, 4);
        o.section.interface_resolution = s;
        const auto m = swept(x, o);
        for (std::size_t k = 0; k < 3; ++k) {
            if (name == "tapered fin" && k != 0) continue;  // the fin's axis: every other facet is parallel to it
            const double e = linear_field_error(m, k);
            INFO(name << " axis " << k << ": " << m.mesh.node_count() << " nodes, deviation " << e);
            CHECK(e <= 1e-12);
        }
    }
}

TEST_CASE("swept mesh: straddled stations, interfaces and a partial gate patch are certified (UD3, UD8)", "[mesh][swept]") {
    const auto x = gated_fin();
    SweptOptions o;
    o.section.max_size = 0.3;
    o.section.interface_resolution = 0.04;
    o.lines = {0.0, 0.3, 0.6, 0.85, 1.5, 2.15, 2.4, 2.7, 3.0};
    o.straddles = {{1.0, 0.1}, {2.0, 0.1}};
    const auto m = swept(x, o);
    REQUIRE(m.certificate.has_value());
    const auto& c = *m.certificate;
    INFO("cross faces " << c.report.cross_faces << ", corner " << c.report.corner_faces << ", flat " << c.report.max_flat_deviation
                        << ", corner deviation " << c.report.max_corner_deviation);
    CHECK(c.report.max_flat_deviation <= 1.0);
    CHECK(c.report.max_corner_deviation <= 0.04);
    // The gate patch ends at the straddled stations: per-corner patches and P6.
    CHECK_FALSE(m.simplices.facet_corner_patch.empty());
    REQUIRE(c.patches.has_value());
    CHECK(c.patches->max_flat_deviation <= 1.0);
    const BoundaryPatch* gate = m.mesh.find_boundary("gate");
    REQUIRE(gate != nullptr);
    double area = 0.0;
    for (const double a : gate->areas) area += a;
    // The overlayer's outer faces in the gate section: z = 0 and 3 (y 1 to 2.5) and y = 2.5 (z 0 to 3), 1 long.
    CHECK(std::abs(area - (1.5 + 1.5 + 3.0)) <= 1e-12 * 6.0);
    // Regions: a node on the line just inside the gate section is oxide (3) where the overlayer is.
    std::size_t oxide = 0;
    std::size_t spacer = 0;
    for (std::size_t i = 0; i < m.node_region.size(); ++i) {
        oxide += m.node_region[i] == 3 ? 1 : 0;
        spacer += m.node_region[i] == 2 ? 1 : 0;
        if (m.node_region[i] == 3) CHECK((m.mesh.points()[i][0] > 1.0 && m.mesh.points()[i][0] < 2.0));
    }
    CHECK(oxide > 0);
    CHECK(spacer > 0);
}

TEST_CASE("swept mesh: refusals (UC2, slice 1)", "[mesh][swept]") {
    const auto refused = [](const ExtrudedGeometry& x, const SweptOptions& o, ErrorCode code, const std::string& text) {
        const auto m = generate_swept(x, o);
        REQUIRE_FALSE(m);
        CHECK(m.error().code == code);
        INFO(m.error().message);
        CHECK(m.error().message.find(text) != std::string::npos);
    };
    SECTION("a void next to a filled section (UD4)") {
        auto x = l_shape();
        x.regions[2] = {0, V};
        refused(x, options_for(x, 0.4, 2), ErrorCode::invalid_input, "void in every section or in none");
    }
    SECTION("lines") {
        const auto x = box();
        auto o = options_for(x, 0.4, 3);
        o.lines.front() = 0.1;
        refused(x, o, ErrorCode::invalid_input, "first and last axis lines");
        o = options_for(x, 0.4, 3);
        std::swap(o.lines[1], o.lines[2]);
        refused(x, o, ErrorCode::invalid_input, "increase strictly");
        o.lines = {0.0};
        refused(x, o, ErrorCode::invalid_input, "at least two axis lines");
    }
    SECTION("straddles") {
        const auto x = gated_fin();
        SweptOptions o;
        o.section.max_size = 0.4;
        o.section.interface_resolution = 0.04;
        o.lines = {0.0, 0.5, 1.5, 2.5, 3.0};
        refused(x, o, ErrorCode::invalid_input, "needs a straddle");
        o.straddles = {{1.0, 0.1}, {2.0, 0.1}, {1.0, 0.2}};
        refused(x, o, ErrorCode::invalid_input, "straddled twice");
        o.straddles = {{1.0, 0.1}, {0.0, 0.1}};
        refused(x, o, ErrorCode::invalid_input, "interior station");
        o.straddles = {{1.0, 0.1}, {2.0, -0.1}};
        refused(x, o, ErrorCode::invalid_input, "finite and positive");
        o.straddles = {{1.0, 0.1}, {2.0, 2.5}};
        refused(x, o, ErrorCode::invalid_input, "inside the neighbouring sections");
        o.straddles = {{1.0, 1.1}, {2.0, 0.1}};
        refused(x, o, ErrorCode::invalid_input, "within a straddle");
    }
    SECTION("a lateral patch changing at a station that is neither straddled nor a line") {
        auto x = box();
        x.axis = 0;
        x.stations = {0.0, 1.0, 2.0};
        x.regions = {{0, 0}};
        x.patches = {{"top", {0.0, 0.0, 1.0}, {1.0, 2.0, 1.0}}};
        auto o = options_for(x, 0.4, 2);
        o.lines = {0.0, 0.5, 1.5, 2.0};
        refused(x, o, ErrorCode::invalid_input, "lateral boundary patch changes");
        o.lines = {0.0, 0.5, 1.0, 1.5, 2.0};  // a line at the station: accepted, the patch edge a mesh line
        const auto m = swept(x, o);
        CHECK(m.simplices.facet_corner_patch.empty());
        const BoundaryPatch* top = m.mesh.find_boundary("top");
        REQUIRE(top != nullptr);
        double area = 0.0;
        for (const double a : top->areas) area += a;
        CHECK(std::abs(area - 2.0) <= 2e-12);
    }
    SECTION("a lateral patch changing where the boundary is straight inside one 2D region") {
        auto x = box();
        x.cells = {rect("a", 0, 0, 1, 1), rect("b", 1, 0, 2, 1)};
        x.regions = {{0}, {0}};
        x.patches = {{"half", {0.0, 0.0, 0.0}, {1.0, 0.0, 3.0}}};  // the y = 0 face of cell a only
        // a and b differ in their lateral patch sets, so they are two 2D regions and the vertex is a
        // junction: accepted, straddled.
        auto o = options_for(x, 0.3, 3);
        o.section.interface_resolution = 0.05;
        const auto m = generate_swept(x, o);
        if (!m) FAIL(m.error().message);
    }
    SECTION("interfaces without a resolution") {
        const auto x = tapered_fin();
        refused(x, options_for(x, 0.4, 3), ErrorCode::invalid_input, "positive interface_resolution");
    }
    SECTION("max_nodes") {
        const auto x = box();
        auto o = options_for(x, 0.3, 7);
        o.max_nodes = 50;
        refused(x, o, ErrorCode::resource_exhausted, "max_nodes");
    }
    SECTION("the cross-section's own errors are named") {
        const auto x = tapered_fin();
        auto o = options_for(x, 0.4, 3);
        o.section.interface_resolution = 0.4;  // features shorter than 6 s
        refused(x, o, ErrorCode::invalid_input, "the cross-section: ");
    }
}
