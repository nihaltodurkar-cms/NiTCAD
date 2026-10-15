// The interface rule P4 (Unit 16; DECISIONS.md M7, M8) on hand-built layouts: mirrored node pairs
// at +-s/2 across each interface, the protected band, corner and junction placements, and the
// failures the check must catch. Each layout is triangulated by brute force (perturbed Delaunay),
// imported with make_dual_mesh, and checked with check_interfaces.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/mesh/interfaces.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/mesh/simplex_mesh.hpp"
#include "simplex_fixtures.hpp"

using namespace NiTCAD::mesh;
using namespace mesh_fixtures;
using NiTCAD::base::ErrorCode;

namespace {

constexpr double s = 1.0 / 32.0;  // pair spacing and resolution
const double band = std::hypot(1.0 / 16.0, s / 2.0);  // sqrt(R^2 + s^2/4), R = 1/16 (pair columns 1/8 apart)

struct Layout {
    int dim = 2;
    std::vector<Point> p;
    std::vector<std::int32_t> region;
    void add(double x, double y, std::int32_t r) {
        p.push_back({x, y, 0});
        region.push_back(r);
    }
    void add3(double x, double y, double z, std::int32_t r) {
        p.push_back({x, y, z});
        region.push_back(r);
    }
    void pair(double x, double y, double nx, double ny, std::int32_t minus, std::int32_t plus) {
        add(x - 0.5 * s * nx, y - 0.5 * s * ny, minus);
        add(x + 0.5 * s * nx, y + 0.5 * s * ny, plus);
    }
};

InterfacePiece segment(Point a, Point b, std::int32_t ra, std::int32_t rb) { return {{a, b}, ra, rb}; }

std::expected<InterfaceReport, NiTCAD::base::Error> check(const Layout& L, const std::vector<InterfacePiece>& pieces) {
    SimplexMesh m;
    m.dimension = L.dim;
    m.points = L.p;
    m.cells = brute_delaunay(L.dim, L.p);
    tag_box_boundary(m, {0, 0, 0}, {1, 1, L.dim == 3 ? 1.0 : 0.0});
    const auto dual = make_dual_mesh(m);
    REQUIRE(dual.has_value());  // every layout here is a valid import; only P4 is in question
    return check_interfaces(m, *dual, L.region, InterfaceSet{pieces, s});
}

// The unit square split at y = 0.5 (region 0 below, 1 above): boundary nodes, the junction pairs
// where the interface meets x = 0 and x = 1, a pair at each x in pair_x, and two interior rows.
Layout flat(const std::vector<double>& pair_x) {
    Layout L;
    for (const double x : {0.0, 0.25, 0.5, 0.75, 1.0}) {
        L.add(x, 0, 0);
        L.add(x, 1, 1);
    }
    for (const double x : {0.0, 1.0}) {
        L.add(x, 0.25, 0);
        L.add(x, 0.75, 1);
        L.pair(x, 0.5, 0, 1, 0, 1);  // junction with the outer boundary: both nodes on it
    }
    for (const double x : pair_x) L.pair(x, 0.5, 0, 1, 0, 1);
    for (const double x : {0.25, 0.5, 0.75}) {
        L.add(x, 0.25, 0);
        L.add(x, 0.75, 1);
    }
    return L;
}

const std::vector<double> columns{0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875};
const std::vector<InterfacePiece> flat_piece{segment({0, 0.5, 0}, {1, 0.5, 0}, 0, 1)};

}  // namespace

TEST_CASE("P4: a flat interface with mirrored pairs is established", "[mesh][interfaces]") {
    const auto r = check(flat(columns), flat_piece);
    REQUIRE(r.has_value());
    CHECK(r->cross_faces > 0);
    CHECK(r->corner_faces > 0);                   // where the interface meets the boundary
    CHECK(r->cross_faces > r->corner_faces);      // and exact faces away from it
    CHECK(r->max_flat_deviation <= 1.0);
    CHECK(r->max_corner_deviation <= s);
}

TEST_CASE("P4: graded pair spacing along the interface stays exact", "[mesh][interfaces]") {
    // Isosceles trapezoids between neighbouring pairs are cocircular: the cross diagonals have
    // zero coupling, whatever the spacing.
    const auto r = check(flat({0.03, 0.07, 0.13, 0.22, 0.36, 0.55, 0.78}), flat_piece);
    REQUIRE(r.has_value());
    CHECK(r->max_flat_deviation <= 1.0);
}

TEST_CASE("P4: the protected band", "[mesh][interfaces]") {
    SECTION("a node inside the band breaks the interface") {
        Layout L = flat(columns);
        L.add(0.3, 0.5 + 0.3 * s, 1);
        const auto r = check(L, flat_piece);
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().code == ErrorCode::degenerate_mesh);
        CHECK(r.error().message.find("P4") != std::string::npos);
    }
    SECTION("a node just outside the band does not (the band condition is sufficient)") {
        Layout L = flat(columns);
        L.add(0.3125, 0.5 + 1.05 * band, 1);  // midway between two pair columns, the worst place
        const auto r = check(L, flat_piece);
        REQUIRE(r.has_value());
        CHECK(r->max_flat_deviation <= 1.0);
    }
}

TEST_CASE("P4: thin layers", "[mesh][interfaces]") {
    // Region 1 is a layer between y = y0 and y = y1; region 0 below, 2 above.
    const auto layer = [](double y0, double y1, double stagger) {
        Layout L;
        for (const double x : {0.0, 0.25, 0.5, 0.75, 1.0}) {
            L.add(x, 0, 0);
            L.add(x, 1, 2);
        }
        for (const double x : {0.0, 1.0}) {
            L.add(x, 0.25, 0);
            L.add(x, 0.75, 2);
            L.pair(x, y0, 0, 1, 0, 1);
            L.pair(x, y1, 0, 1, 1, 2);
        }
        for (const double x : columns) {
            L.pair(x, y0, 0, 1, 0, 1);
            if (x + stagger < 1.0) L.pair(x + stagger, y1, 0, 1, 1, 2);
        }
        if (stagger > 0.0) {
            L.pair(stagger, y1, 0, 1, 1, 2);  // the staggered columns start at x = stagger: R stays 1/16
        }
        for (const double x : {0.25, 0.5, 0.75}) {
            L.add(x, 0.25, 0);
            L.add(x, 0.75, 2);
        }
        return std::pair{L, std::vector<InterfacePiece>{segment({0, y0, 0}, {1, y0, 0}, 0, 1),
                                                         segment({0, y1, 0}, {1, y1, 0}, 1, 2)}};
    };
    SECTION("a layer wider than both bands is established") {
        const auto [L, pieces] = layer(0.45, 0.55, 1.0 / 16.0);
        const auto r = check(L, pieces);
        REQUIRE(r.has_value());
        CHECK(r->max_flat_deviation <= 1.0);
    }
    SECTION("a layer inside one interface's band, with staggered pairs, is refused") {
        const auto [L, pieces] = layer(0.48, 0.52, 1.0 / 16.0);
        const auto r = check(L, pieces);
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().code == ErrorCode::degenerate_mesh);
    }
}

TEST_CASE("P4: a corner, by the bisector rule (M8)", "[mesh][interfaces]") {
    // Region 1 is the quadrant [0.5, 1]^2; the interface bends by 90 degrees at (0.5, 0.5).
    // Pairs approach the corner at 1.5 s and 3 s; at the corner one node per region sits on its
    // sector's bisector at s / 2.
    Layout L;
    for (const double x : {0.0, 0.25, 0.5, 0.75, 1.0}) L.add(x, 0, 0);
    for (const double y : {0.25, 0.5, 0.75}) L.add(0, y, 0);
    L.add(0, 1, 0);
    L.add(0.25, 1, 0);
    L.pair(0.5, 1, 1, 0, 0, 1);  // the vertical interface meets y = 1
    L.add(0.75, 1, 1);
    L.add(1, 1, 1);
    L.add(1, 0.25, 0);
    L.pair(1, 0.5, 0, 1, 0, 1);  // the horizontal interface meets x = 1
    L.add(1, 0.75, 1);
    for (const double t : {1.5 * s, 3 * s, 0.125, 0.25, 0.375}) {
        L.pair(0.5, 0.5 + t, 1, 0, 0, 1);
        L.pair(0.5 + t, 0.5, 0, 1, 0, 1);
    }
    const double c = 0.5 * s / std::sqrt(2.0);
    L.add(0.5 - c, 0.5 - c, 0);
    L.add(0.5 + c, 0.5 + c, 1);
    L.add(0.25, 0.25, 0);
    L.add(0.25, 0.75, 0);
    L.add(0.75, 0.25, 0);
    L.add(0.75, 0.75, 1);
    const std::vector<InterfacePiece> pieces{segment({0.5, 0.5, 0}, {0.5, 1, 0}, 0, 1),
                                             segment({0.5, 0.5, 0}, {1, 0.5, 0}, 0, 1)};
    const auto r = check(L, pieces);
    REQUIRE(r.has_value());
    CHECK(r->max_corner_deviation <= s);
    CHECK(r->max_flat_deviation <= 1.0);
    INFO("corner deviation / s = " << r->max_corner_deviation / s);
    SUCCEED();
}

TEST_CASE("P4: a T-junction of three regions (M8)", "[mesh][interfaces]") {
    // Region 0 below y = 0.5; above it region 1 left of x = 0.5 and region 2 right of it.
    Layout L;
    for (const double x : {0.0, 0.25, 0.5, 0.75, 1.0}) L.add(x, 0, 0);
    L.add(0, 1, 1);
    L.add(0.25, 1, 1);
    L.pair(0.5, 1, 1, 0, 1, 2);
    L.add(0.75, 1, 2);
    L.add(1, 1, 2);
    L.add(0, 0.25, 0);
    L.pair(0, 0.5, 0, 1, 0, 1);
    L.add(0, 0.75, 1);
    L.add(1, 0.25, 0);
    L.pair(1, 0.5, 0, 1, 0, 2);
    L.add(1, 0.75, 2);
    for (const double t : {1.5 * s, 3 * s, 0.125, 0.25, 0.375}) {
        L.pair(0.5 - t, 0.5, 0, 1, 0, 1);
        L.pair(0.5 + t, 0.5, 0, 1, 0, 2);
        L.pair(0.5, 0.5 + t, 1, 0, 1, 2);
    }
    const double c = 0.5 * s / std::sqrt(2.0);
    L.add(0.5, 0.5 - 0.5 * s, 0);
    L.add(0.5 - c, 0.5 + c, 1);
    L.add(0.5 + c, 0.5 + c, 2);
    L.add(0.25, 0.25, 0);
    L.add(0.75, 0.25, 0);
    L.add(0.25, 0.75, 1);
    L.add(0.75, 0.75, 2);
    const std::vector<InterfacePiece> pieces{segment({0, 0.5, 0}, {0.5, 0.5, 0}, 0, 1),
                                             segment({0.5, 0.5, 0}, {1, 0.5, 0}, 0, 2),
                                             segment({0.5, 0.5, 0}, {0.5, 1, 0}, 1, 2)};
    const auto r = check(L, pieces);
    REQUIRE(r.has_value());
    CHECK(r->max_corner_deviation <= s);
    CHECK(r->max_flat_deviation <= 1.0);
}

// The right-angle quadrant layout (region 1 = [0.5, 1]^2) with the pairs stopping `gap` short of
// the corner, with or without the bisector nodes, and optional extra nodes.
Layout quadrant(double gap, bool bisector_nodes, double bisector_distance = 0.5 * s) {
    Layout L;
    for (const double x : {0.0, 0.25, 0.5, 0.75, 1.0}) L.add(x, 0, 0);
    for (const double y : {0.25, 0.5, 0.75}) L.add(0, y, 0);
    L.add(0, 1, 0);
    L.add(0.25, 1, 0);
    L.pair(0.5, 1, 1, 0, 0, 1);
    L.add(0.75, 1, 1);
    L.add(1, 1, 1);
    L.add(1, 0.25, 0);
    L.pair(1, 0.5, 0, 1, 0, 1);
    L.add(1, 0.75, 1);
    for (const double t : {1.5 * s, 3 * s, 0.125, 0.25, 0.375}) {
        if (t < gap) continue;
        L.pair(0.5, 0.5 + t, 1, 0, 0, 1);
        L.pair(0.5 + t, 0.5, 0, 1, 0, 1);
    }
    if (bisector_nodes) {
        const double c = bisector_distance / std::sqrt(2.0);
        L.add(0.5 - c, 0.5 - c, 0);
        L.add(0.5 + c, 0.5 + c, 1);
    }
    L.add(0.25, 0.25, 0);
    L.add(0.25, 0.75, 0);
    L.add(0.75, 0.25, 0);
    L.add(0.75, 0.75, 1);
    return L;
}

const std::vector<InterfacePiece> quadrant_pieces{segment({0.5, 0.5, 0}, {0.5, 1, 0}, 0, 1),
                                                  segment({0.5, 0.5, 0}, {1, 0.5, 0}, 0, 1)};

TEST_CASE("P4: a right-angle corner needs no corner nodes (mirror symmetry)", "[mesh][interfaces]") {
    // With mirrored pairs on both arms the corner point is equidistant from the nearest pair of
    // each arm, so the dual interface follows the bend exactly even with the pairs 8 s away.
    const auto r = check(quadrant(0.25, false), quadrant_pieces);
    REQUIRE(r.has_value());
    CHECK(r->max_corner_deviation <= 1e-12);
}

TEST_CASE("P4b: a node on the wrong side of a corner is refused", "[mesh][interfaces]") {
    // A region-1 node at (0.47, 0.47), inside region 0 and within the corner zone: its cell is a
    // region-1 cell in region-0 territory, deviating beyond s.
    Layout L = quadrant(0.125, true);
    L.add(0.47, 0.47, 1);
    const auto r = check(L, quadrant_pieces);
    REQUIRE_FALSE(r.has_value());
    CHECK(r.error().message.find("P4") != std::string::npos);  // P4a or P4b: its faces leave the zone
}

TEST_CASE("P4: a 60-degree corner by the bisector rule (M8), measured", "[mesh][interfaces]") {
    // Region 1 is the wedge between the arm along +x from (0.5, 0.5) and the arm at 60 degrees,
    // which meets y = 1 at x = 0.5 + 0.5 / tan 60. Pairs at s / 2 along both arms (normals into
    // and out of the wedge), one bisector node per sector at s / 2 from the apex.
    const double pi = 3.141592653589793;
    const double cs = std::cos(pi / 3), sn = std::sin(pi / 3);
    const double top = 0.5 + 0.5 * cs / sn;
    Layout L;
    for (const double x : {0.0, 0.25, 0.5, 0.75, 1.0}) L.add(x, 0, 0);
    for (const double y : {0.25, 0.5, 0.75, 1.0}) L.add(0, y, 0);
    L.add(0.25, 1, 0);
    L.add(0.5, 1, 0);
    L.add(top - 0.5 * s, 1, 0);  // the arm meets y = 1: a pair straddling it on the boundary
    L.add(top + 0.5 * s, 1, 1);
    for (const double t : {1.5 * s, 3 * s, 6 * s}) L.add(top - t, 1, 0);  // boundary graded toward it
    for (const double t : {1.5 * s, 3 * s}) L.add(top + t, 1, 1);
    L.add(1, 1, 1);
    L.add(1, 0.25, 0);
    L.pair(1, 0.5, 0, 1, 0, 1);
    L.add(1, 0.75, 1);
    // Pairs graded toward both ends of each arm (the apex and the boundary junction).
    const double arm2 = 0.5 / sn;  // the 60-degree arm's length to y = 1
    for (const double t : {1.5 * s, 3 * s, 0.125, 0.25, 0.375}) L.pair(0.5 + t, 0.5, 0, 1, 0, 1);
    for (const double t : {1.5 * s, 3 * s, 0.125, 0.25, 0.375, arm2 - 0.125, arm2 - 3 * s, arm2 - 1.5 * s}) {
        L.pair(0.5 + t * cs, 0.5 + t * sn, sn, -cs, 0, 1);  // region 1 below-right of the arm
    }
    const double c = 0.5 * s;
    L.add(0.5 + c * std::cos(pi / 6), 0.5 + c * std::sin(pi / 6), 1);           // wedge bisector, 30 deg
    L.add(0.5 + c * std::cos(7 * pi / 6), 0.5 + c * std::sin(7 * pi / 6), 0);   // the rest, 210 deg
    L.add(0.25, 0.25, 0);
    L.add(0.25, 0.75, 0);
    L.add(0.75, 0.25, 0);
    L.add(0.9, 0.75, 1);
    const std::vector<InterfacePiece> pieces{segment({0.5, 0.5, 0}, {1, 0.5, 0}, 0, 1),
                                             segment({0.5, 0.5, 0}, {top, 1, 0}, 0, 1)};
    const auto r = check(L, pieces);
    INFO((r.has_value() ? "established; corner deviation / s = " + std::to_string(r->max_corner_deviation / s)
                        : "refused: " + r.error().message));
    REQUIRE(r.has_value());
    CHECK(r->max_corner_deviation <= s);
    CHECK(r->max_flat_deviation <= 1.0);
}

TEST_CASE("P4: undeclared region pairs and imported meshes are not interface-valid", "[mesh][interfaces]") {
    SECTION("a region pair with no declared piece") {
        Layout L = flat(columns);
        L.region.back() = 2;  // the interior node (0.75, 0.75) now borders regions 1 and 2
        const auto r = check(L, flat_piece);
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().message.find("no declared interface") != std::string::npos);
    }
    SECTION("an imported mesh with regions assigned by position") {
        // Scattered points, a valid import, regions by y < 0.5: P4 cannot be established.
        Layout L;
        for (int k = 0; k <= 5; ++k) L.add(k / 5.0, 0, 0);
        for (int k = 0; k <= 6; ++k) L.add(k / 6.0, 1, 1);
        for (int k = 1; k < 7; ++k) L.add(0, k / 7.0, k / 7.0 < 0.5 ? 0 : 1);
        for (int k = 1; k < 4; ++k) L.add(1, k / 4.0, k / 4.0 < 0.5 ? 0 : 1);
        SplitMix rng{7};
        for (int k = 0; k < 14; ++k) {
            const double x = 0.25 + 0.5 * rng.uniform();
            const double y = 0.25 + 0.5 * rng.uniform();
            L.add(x, y, y < 0.5 ? 0 : 1);
        }
        const auto r = check(L, flat_piece);
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().code == ErrorCode::degenerate_mesh);
    }
    SECTION("malformed interface input") {
        const Layout L = flat(columns);
        SimplexMesh m;
        m.points = L.p;
        m.cells = brute_delaunay(2, L.p);
        tag_box_boundary(m, {0, 0, 0}, {1, 1, 0});
        const auto dual = make_dual_mesh(m);
        REQUIRE(dual.has_value());
        const auto invalid = [](const std::expected<InterfaceReport, NiTCAD::base::Error>& r) {
            return !r.has_value() && r.error().code == ErrorCode::invalid_input;
        };
        CHECK(invalid(check_interfaces(m, *dual, L.region, InterfaceSet{flat_piece, 0.0})));
        CHECK(invalid(check_interfaces(m, *dual, std::span(L.region).first(3), InterfaceSet{flat_piece, s})));
        CHECK(invalid(check_interfaces(m, *dual, L.region,
                                       InterfaceSet{{segment({0, 0.5, 0}, {1, 0.5, 0}, 1, 1)}, s})));
        SimplexMesh other = m;
        other.points[0][0] = 1e-3;
        CHECK(invalid(check_interfaces(other, *dual, L.region, InterfaceSet{flat_piece, s})));
    }
}

TEST_CASE("P4 in 3D: a planar interface, and a broken mirror", "[mesh][interfaces]") {
    // The unit cube split at z = 0.5: a graded lattice in x and y, the pair planes at 0.5 -+ s/2.
    const auto layout = [](double shift) {
        Layout L;
        L.dim = 3;
        for (const double z : {0.0, 0.5 - 0.5 * s, 0.5 + 0.5 * s, 1.0}) {
            for (const double y : {0.0, 0.3, 1.0}) {
                for (const double x : {0.0, 0.3, 1.0}) {
                    const bool interior = x == 0.3 && y == 0.3 && z > 0.5 && z < 1.0;
                    L.add3(interior ? x + shift : x, y, z, z < 0.5 ? 0 : 1);
                }
            }
        }
        return L;
    };
    const std::vector<InterfacePiece> plane{{{{0, 0, 0.5}, {1, 0, 0.5}, {1, 1, 0.5}, {0, 1, 0.5}}, 0, 1}};
    SECTION("mirrored") {
        const auto r = check(layout(0.0), plane);
        REQUIRE(r.has_value());
        CHECK(r->cross_faces > r->corner_faces);
        CHECK(r->max_flat_deviation <= 1.0);
    }
    SECTION("one pair node moved off its mirror") {
        const auto r = check(layout(0.05), plane);
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().code == ErrorCode::degenerate_mesh);
    }
}

TEST_CASE("P4b: the corner threshold is the resolution", "[mesh][interfaces]") {
    // The bisector-rule corner (its largest corner deviation is 0.0249, 0.80 s) checked against
    // a smaller declared resolution: accepted at 0.026, refused by P4b at 0.024, the deviating
    // faces staying inside the (smaller) corner zone.
    const Layout L = quadrant(1.5 * s, true);
    SimplexMesh m;
    m.points = L.p;
    m.cells = brute_delaunay(2, L.p);
    tag_box_boundary(m, {0, 0, 0}, {1, 1, 0});
    const auto dual = make_dual_mesh(m);
    REQUIRE(dual.has_value());
    const auto accepted = check_interfaces(m, *dual, L.region, InterfaceSet{quadrant_pieces, 0.026});
    REQUIRE(accepted.has_value());
    const auto refused = check_interfaces(m, *dual, L.region, InterfaceSet{quadrant_pieces, 0.024});
    REQUIRE_FALSE(refused.has_value());
    CHECK(refused.error().message.find("P4b") != std::string::npos);
}

TEST_CASE("P4: every vertex of a boundary dual face counts", "[mesh][interfaces]") {
    // Hand-built meshes with pieces through some vertices of each cross face and far from the
    // others: the check must see the far ones.
    SECTION("2D: the boundary end of the face (the edge midpoint)") {
        // The unit square in two triangles, node (0,0) alone in region 0. The cross faces run from
        // the circumcentre (0.5, 0.5) to the midpoints (0.5, 0) and (0, 0.5); the pieces touch
        // only the circumcentre end.
        const SimplexMesh m = two_d({{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}, {0, 1, 2, 0, 2, 3});
        const auto dual = make_dual_mesh(m);
        REQUIRE(dual.has_value());
        const std::vector<std::int32_t> region{0, 1, 1, 1};
        const std::vector<InterfacePiece> pieces{segment({0.5, 0.5, 0}, {0.6, 0.5, 0}, 0, 1)};
        const auto r = check_interfaces(m, *dual, region, InterfaceSet{pieces, 0.2});
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().message.find("P4") != std::string::npos);
    }
    SECTION("3D: the circumcentres of the boundary facets around the edge") {
        // The unit cube in five tetrahedra, node (0,0,0) alone in region 0. The face of each cube
        // edge at that node is the square (centre, two face centres, edge midpoint); thin pieces
        // join each edge midpoint to the centre, so only the face centres deviate (by 0.35).
        SimplexMesh m;
        m.dimension = 3;
        m.points = grid_points({0, 1}, {0, 1}, {0, 1});
        m.cells = {0, 3, 5, 6, 1, 0, 3, 5, 2, 0, 3, 6, 4, 0, 5, 6, 7, 3, 5, 6};
        const auto dual = make_dual_mesh(m);
        REQUIRE(dual.has_value());
        const std::vector<std::int32_t> region{0, 1, 1, 1, 1, 1, 1, 1};
        const auto thin = [](Point mid, Point centre, Point side) {
            return InterfacePiece{{mid, centre, side}, 0, 1};
        };
        const std::vector<InterfacePiece> pieces{
            thin({0.5, 0, 0}, {0.5, 0.5, 0.5}, {0.5, 0.26, 0.24}),
            thin({0, 0.5, 0}, {0.5, 0.5, 0.5}, {0.24, 0.5, 0.26}),
            thin({0, 0, 0.5}, {0.5, 0.5, 0.5}, {0.26, 0.24, 0.5})};
        const auto r = check_interfaces(m, *dual, region, InterfaceSet{pieces, 0.2});
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().message.find("P4") != std::string::npos);
    }
}
