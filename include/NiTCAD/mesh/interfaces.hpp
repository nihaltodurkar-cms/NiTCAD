// The interface rule P4 (Unit 16; DECISIONS.md M7, M8): whether the control volumes of a dual mesh
// meet the declared material interfaces.
//
// Under M7 every node belongs to one region and the interface between two regions is represented
// by the dual faces between their nodes (the faces of the edges joining a node of each). The
// assembler treats each such edge as crossing the interface at its midpoint
// (assemble/interface_edges.hpp). check_interfaces is the only way to establish that this holds
// for a mesh: make_dual_mesh does not (an imported mesh carries no interfaces), and a mesh that
// has not passed check_interfaces must not be treated as interface-valid.
//
// Metric. For each cross-region dual face F (of an edge e = ij with regions A != B, kept in the
// dual graph), with vertex set V(F) (the circumcentres of the simplices around e, and on the
// boundary the circumcentres of the boundary facets around e and e's midpoint):
//   deviation(F) = min over the declared pieces P separating A and B of max over v in V(F) of
//   dist(v, P).
// For a convex piece (a segment, or a planar convex polygon) dist(., P) is convex, and F lies in
// the convex hull of V(F): max over V(F) of dist(., P) bounds it over all of F. A non-convex
// planar polygon (3D; Unit 53, FD9) is measured by the plane distance when the hull of V(F)'s
// projection lies in the polygon (tested exactly; axis-aligned pieces, whose projection is exact),
// by the plane distance plus the grown reach when it lies in the polygon grown outward by 1e-6 s
// (mitred; the reach its largest vertex displacement), for face vertices rounded just across a
// feature line; and otherwise by the min over the triangles of its ear-clip decomposition of max
// over V(F) (an upper bound; always for oblique pieces).
// The minimum over pieces makes deviation(F) an upper bound on the largest distance from F to the
// interface (exact where one piece is nearest everywhere on F; an over-estimate near a bend, the
// safe direction).
// Corner zones: the points within 2 s of a piece's relative boundary (a segment's end points; a
// polygon's edges), s the declared resolution. A face is in a zone when its midpoint or any vertex
// is. Acceptance:
// - outside the zones (P4a, the exact part of M7): deviation(F) <= 1e-12 L_e + 1e-13 X, X the
//   largest coordinate magnitude among F's vertices and e's end points (the first term relative to
//   the edge, the second the rounding of coordinates far from the origin, about 450 ulps). This
//   holds exactly in real arithmetic for a mirrored pair set with its protected band (M7), and
//   then e is a mirrored pair and its midpoint lies on the interface;
// - inside the zones (P4b, M8): deviation(F) <= s;
// - every pair of regions joined by a kept edge must have a declared piece.
// Not checked here: that each node lies inside its region's geometry (the geometry of whole regions
// is not an input); a misplaced node shows up as a deviating face unless it is surrounded by its
// own region's nodes.
#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <optional>
#include <span>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/mesh/simplex_mesh.hpp"

namespace NiTCAD::mesh {

struct InterfacePiece {
    // 2D: a segment, two points. 3D: a planar simple polygon (convex, or not: FD9), at least three
    // points in order.
    std::vector<Point> vertices;
    std::int32_t region_a = 0;
    std::int32_t region_b = 0;
};

struct InterfaceSet {
    std::vector<InterfacePiece> pieces;
    double resolution = 0.0;  // s [cm]: corner zones have radius 2 s; faces in them may deviate s
};

struct InterfaceReport {
    std::size_t cross_faces = 0;        // dual faces between nodes of different regions
    std::size_t corner_faces = 0;       // of those, the faces in a corner zone
    double max_flat_deviation = 0.0;    // outside the zones, deviation / its tolerance (<= 1)
    double max_corner_deviation = 0.0;  // inside the zones [cm] (<= s)
};

// `dual` must be make_dual_mesh(simplices) (same dimension, node count and points).
// Errors:
// - invalid_input: node_region not one per node, or a region negative; dual not of these
//   simplices; resolution not finite and positive; a piece with the wrong vertex count, a
//   non-finite vertex, equal or negative regions, or (3D) not planar, simple and non-degenerate
//   (the index names the piece);
// - degenerate_mesh, the rule not established (the index is the edge's first node, both nodes in
//   the message): regions joined with no declared piece between them; a face outside the corner
//   zones deviating beyond its tolerance (P4a); a face in a corner zone deviating beyond s (P4b).
[[nodiscard]] std::expected<InterfaceReport, base::Error> check_interfaces(
    const SimplexMesh& simplices, const Mesh& dual, std::span<const std::int32_t> node_region,
    const InterfaceSet& interfaces);


// P6, patch boundaries (Unit 53; DECISIONS.md FD8, FD10): whether the boundary areas of a 3D mesh
// whose patch edges are not all mesh edges are split along the declared patch boundaries. Under
// FD8 (a) a patch edge lies on an interface trace and boundary triangles straddle it; each corner's
// boundary piece goes to its own patch (SimplexMesh::facet_corner_patch). On every listed boundary
// triangle, each edge whose two corners carry different patches has a boundary dual piece in the
// triangle: the segment from the edge's midpoint to the triangle's circumcentre, the boundary
// between the corners' pieces. Its deviation is min over the declared boundaries of its two
// patches of max over its two end points of the distance (an upper bound, as P4's). Zones: the
// points within 2 s of a declared zone centre; a piece is in a zone when either end point is.
// Acceptance: outside the zones deviation <= 1e-12 L_e + 1e-13 X (P4a's tolerance, L_e the edge,
// X the largest coordinate magnitude of the piece's ends and the edge's nodes); inside, <= s; two
// patches meeting on a triangle need a declared boundary between them. Reported, not judged:
// each patch's area (the dual's) against its geometric area, and the length of its declared
// boundary inside the zones.
struct PatchBoundary {
    Point a{};
    Point b{};
    std::int32_t patch_a = 0;  // indices into SimplexMesh::patch_names, patch_a != patch_b
    std::int32_t patch_b = 0;
};

struct PatchBoundarySet {
    std::vector<PatchBoundary> boundaries;
    std::vector<Point> zone_centres;  // where patch boundaries meet traces, domain edges or each other, or turn on a trace
    std::vector<double> patch_area;   // the geometric area of each patch [cm^2], one per patch name
    double resolution = 0.0;          // s [cm]
};

struct PatchBoundaryReport {
    std::size_t cross_pieces = 0;       // boundary dual pieces between corners of different patches
    std::size_t corner_pieces = 0;      // of those, the pieces in a zone
    double max_flat_deviation = 0.0;    // outside the zones, deviation / its tolerance (<= 1)
    double max_corner_deviation = 0.0;  // inside the zones [cm] (<= s)
    std::vector<double> area_error;     // per patch: the dual's area minus the geometric [cm^2]
    std::vector<double> zone_length;    // per patch: its declared boundary inside the zones [cm]
};

// `dual` must be make_dual_mesh(simplices) of a 3D mesh. Errors:
// - invalid_input: not 3D or dual not of these simplices; resolution not finite and positive;
//   patch_area not one finite non-negative entry per patch; a boundary or zone centre not finite,
//   a boundary degenerate or its patches out of range or equal (the index names it);
// - degenerate_mesh, the rule not established (the index is the listed facet): two patches meet
//   with no declared boundary between them; a piece outside the zones deviating beyond its
//   tolerance; a piece in a zone deviating beyond s.
[[nodiscard]] std::expected<PatchBoundaryReport, base::Error> check_patch_boundaries(const SimplexMesh& simplices,
                                                                                     const Mesh& dual,
                                                                                     const PatchBoundarySet& patches);

// A record that check_interfaces established P4 for one mesh, one node-region assignment and one
// set of declared interfaces (DECISIONS.md U2). Device::create requires a matching one for a
// multi-region device on a mesh that is not a tensor grid; a certificate whose digests do not
// match the device's mesh and node regions is stale and refused. With P6 checked as well (Unit 53,
// FD10), its report; the mesh digest covers the boundary patches.
struct InterfaceCertificate {
    std::uint64_t mesh_digest = 0;        // graph_digest of the dual mesh
    std::uint64_t region_digest = 0;      // region_digest of the node regions
    std::uint64_t interfaces_digest = 0;  // the declared pieces and resolution (and patch set, with P6)
    InterfaceReport report;
    std::optional<PatchBoundaryReport> patches;  // P6 (FD10)
};

// FNV-1a 64 of a node-region assignment (count, then each region as a 64-bit word).
[[nodiscard]] std::uint64_t region_digest(std::span<const std::int32_t> node_region);

// check_interfaces (and check_patch_boundaries when `patches` is given), and on success the
// certificate of what it checked. Errors: those of the checks; invalid_input when the simplices
// carry per-corner patches and `patches` is null (FD10: such a mesh is certified with P6 or not
// at all, DECISIONS.md SQ3).
[[nodiscard]] std::expected<InterfaceCertificate, base::Error> certify_interfaces(
    const SimplexMesh& simplices, const Mesh& dual, std::span<const std::int32_t> node_region,
    const InterfaceSet& interfaces, const PatchBoundarySet* patches = nullptr);

// Whether the certificate was issued for this mesh and these node regions.
[[nodiscard]] bool certifies(const InterfaceCertificate& c, const Mesh& mesh,
                             std::span<const std::int32_t> node_region);


// P4 in 1D (DECISIONS.md U16): an interface is a point between two regions. In 1D the dual face of
// an edge is its midpoint, so a mesh establishes P4 when the midpoint of every edge joining two
// regions lies on a declared point of that pair, within 1e-12 L + 1e-13 X (as P4a), and every
// pair of regions that meets is declared. A verified 1D tensor grid needs no certificate; any
// other 1D mesh with more than one region needs this one.
struct InterfacePoint1D {
    double x = 0.0;  // [cm]
    std::int32_t region_a = 0;
    std::int32_t region_b = 0;
};

// Errors: invalid_input (the mesh not 1D, node_region not one per node or negative, a point not
// finite or with equal or negative regions); degenerate_mesh, the rule not established (an edge
// between regions with no declared pair, or its midpoint off every declared point of its pair;
// the index is the edge's first node).
[[nodiscard]] std::expected<InterfaceCertificate, base::Error> certify_interfaces_1d(
    const Mesh& mesh, std::span<const std::int32_t> node_region,
    std::span<const InterfacePoint1D> interfaces);

}  // namespace NiTCAD::mesh
