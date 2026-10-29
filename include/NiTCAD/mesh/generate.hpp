// The 2D mesh generator (Unit 16; DECISIONS.md M1-M15, U2, U12).
//
// Input: a partition of the domain into region polygons, and named boundary patches. Output: a
// boundary-conforming Delaunay triangulation, its box-method graph (make_dual_mesh: P0-P3, P5),
// the region of every node, and, with two or more regions, the P4 certificate (check_interfaces)
// that Device::create requires.
//
// Geometry rules (invalid_input otherwise):
// - each polygon: at least three distinct vertices, finite, z = 0, simple, of non-zero area; its
//   region index >= 0; either orientation;
// - the polygons tile the domain: they meet only along whole shared edges or at shared vertices
//   (bit-identical coordinates), never overlap, never nest. Edges shared by polygons of different
//   regions are the material interfaces; edges shared by the same region are internal;
// - each patch is a segment inside one boundary edge, patches do not overlap, names are unique;
// - interfaces (DECISIONS.md M7, U12): every interface segment between two features is at least
//   6 s long; features are kept 4 s apart from unrelated interfaces and boundary edges; every
//   sector at an interface feature (a bend, a junction of three or more regions, a junction with
//   the outer boundary) is at least 60 degrees. Unsupported, and refused: an interface meeting the
//   boundary at a domain corner, two interfaces meeting the boundary at one point, a patch end
//   within 3 s of an interface junction.
//
// Node placement:
// - along each interface, mirrored pairs at +-s/2, graded toward each end (1.5 s, 3 s) and at
//   most min(h, 4 s) apart elsewhere;
// - no corner nodes at any interior feature (DECISIONS.md U15: the graded pairs alone measured
//   better than bisector nodes everywhere; a 90-degree bend is exact by mirror symmetry); P4b
//   judges every feature; at a junction with the boundary a node at s/2 on either side along it;
// - boundary nodes at domain corners and at patch ends (exactly at the given points), then at most h
//   apart.
// Refinement (Ruppert): encroached boundary segments are split (at midpoints; next to a domain
// corner at a power of two from it); triangles with a smallest angle below min_angle_deg or a
// circumradius above h / sqrt 3 get their circumcentre, unless it lies in an interface's
// protected band (M7) or every vertex is an interface node, when the triangle is exempt and
// counted. Determinism: the same input gives the same mesh in the same build.
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <optional>
#include <string>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/mesh/complex.hpp"
#include "NiTCAD/mesh/interfaces.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/mesh/simplex_mesh.hpp"

namespace NiTCAD::mesh {

struct RegionPolygon {
    std::vector<Point> vertices;  // cm, z = 0
    std::int32_t region = 0;
};

struct PatchSegment {
    std::string name;
    Point a;
    Point b;
};

struct Geometry2D {
    std::vector<RegionPolygon> polygons;
    std::vector<PatchSegment> patches;  // the result's boundary patches, in this order
};

// A structured boundary layer (DECISIONS.md U14): a tensor lattice along a boundary segment, for
// thin features the isotropic refinement cannot resolve within a node budget (a MOSFET's
// inversion layer). Its nodes are a + t u + o n for t in `along` and o in {0} and `offsets`, u the
// segment's direction and n the inward normal of its boundary edge. The lattice's rectangles are
// Delaunay (all ties), so its dual is exactly the tensor grid's box method; refinement never
// inserts into the layer's rectangle and does not judge triangles of layer nodes alone. Patch
// ends on the segment are added to `along`. The edge carries only the columns, however coarse
// against h; the columns at `a`, `b` and patch ends are those points exactly (Unit 49).
// Rules (invalid_input otherwise): the segment inside one boundary edge; `along` strictly
// increasing from 0 to the segment's length (to 1e-12 relative); `offsets` strictly increasing and
// positive; the lattice inside the polygon owning the edge; layers not overlapping; no interface
// within its protected band of a layer.
struct BoundaryLayer {
    Point a;
    Point b;
    std::vector<double> along;    // [cm] positions from a toward b
    std::vector<double> offsets;  // [cm] distances into the domain
};

struct SizeBox {
    Point lo;
    Point hi;
    double max_size = 0.0;  // h inside the box [cm]
};

struct MeshOptions2D {
    double max_size = 0.0;                // h [cm], finite and positive
    std::vector<double> region_max_size;  // per region index; 0 or absent: max_size
    std::vector<SizeBox> boxes;           // smaller h where given
    std::vector<BoundaryLayer> layers;    // structured layers (U14)
    double interface_resolution = 0.0;    // s [cm]; required when there are interfaces
    double min_angle_deg = 20.7;          // in (0, 33.8]
    std::size_t max_nodes = 2'000'000;
};

struct MeshReport {
    std::size_t nodes = 0;
    std::size_t triangles = 0;
    std::size_t interface_nodes = 0;   // pair, corner and junction nodes
    std::size_t refinement_points = 0;  // circumcentres inserted
    std::size_t segment_splits = 0;
    std::size_t exempt_triangles = 0;  // skinny or large, but in a band or all interface nodes
    double min_angle_deg = 0.0;        // over the triangles not exempt
    // 3D only (DECISIONS.md M10, M11).
    std::size_t tetrahedra = 0;
    std::size_t facet_splits = 0;       // boundary subfacets split at their circumcentres
    std::size_t exempt_tetrahedra = 0;  // poor or large, but not refinable (see generate_mesh 3D)
    double max_radius_edge = 0.0;       // over the tetrahedra not exempt
    double min_dihedral_deg = 0.0;      // over all tetrahedra
    std::size_t slivers = 0;            // tetrahedra with a dihedral angle below 5 degrees
    std::size_t flat_dropped = 0;       // tetrahedra dropped with all four vertices on one face
    std::size_t band_facet_points = 0;  // subfacet circumcentres inserted within an interface band (judged by P4)
};

struct GeneratedMesh {
    Mesh mesh;
    SimplexMesh simplices;
    std::vector<std::int32_t> node_region;
    std::optional<InterfaceCertificate> certificate;  // present with two or more regions
    MeshReport report;
};

// Errors: invalid_input (the geometry or option rules above); resource_exhausted (max_nodes
// reached); degenerate_mesh (the generated mesh failed an import check or P4: for an interface
// feature, the corner rule did not meet P4b, DECISIONS.md U12).
[[nodiscard]] std::expected<GeneratedMesh, base::Error> generate_mesh(const Geometry2D& geometry,
                                                                     const MeshOptions2D& options);


// The 3D generator (Unit 16; DECISIONS.md M6, M9-M13): convex polyhedra.
//
// Input: a convex polyhedron with named faces. Output: a boundary-conforming Delaunay
// tetrahedralization (every boundary triangle Gabriel, P2), its box-method graph (make_dual_mesh),
// and the region of every node.
//
// Geometry rules (invalid_input otherwise, the index naming the face or vertex):
// - vertices finite; each face at least three distinct vertex indices in range, a non-empty name
//   unique among the faces;
// - each face exactly planar (every vertex on the plane of its first three non-collinear ones, by
//   exact orientation), convex and simple (every vertex on the same side of every edge's line,
//   exactly);
// - the polyhedron closed (every edge in exactly two faces, V - E + F = 2) and strictly convex
//   (every vertex not in a face strictly on one side of its plane, the same side for all), so of
//   non-zero volume. Non-convex domains are refused (DECISIONS.md U7).
//
// Node placement and refinement (Shewchuk's 3D Delaunay refinement, DECISIONS.md M11): the
// polyhedron's vertices; then, in this priority, encroached boundary edge segments are split (at
// midpoints; next to a polyhedron vertex at a power of two from it), encroached boundary
// subfacets (the triangles of each face's 2D Delaunay triangulation) are split at their
// circumcentres, and tetrahedra with a circumradius-to-shortest-edge ratio above radius_edge or a
// circumradius above sqrt(3/8) h (a regular tetrahedron of edge h) get their circumcentre. A
// subsegment or subfacet is encroached when it is missing from the tetrahedralization or a vertex
// lies strictly inside its diametral ball. A candidate circumcentre that encroaches subsegments
// or subfacets is not inserted: they are split instead. A candidate outside the domain that
// encroaches nothing is not inserted, and its tetrahedron is exempt and counted. Termination at
// radius_edge = 2 for input angles of at least 90 degrees is cited (Shewchuk 1998, checked in
// DECISIONS.md U6); for smaller input angles, patches and interfaces it is heuristic, capped by
// max_nodes.
// Points on axis-aligned faces and edges keep their constant coordinates exactly, so they are
// exactly coplanar or collinear. Points on oblique faces are rounded off the plane by about an
// ulp; the tetrahedra whose four vertices lie on one face (by construction) are dropped, and the
// result is checked: every boundary triangle on one face, the volume equal to the polyhedron's to
// 1e-12, then make_dual_mesh's P0-P3 and P5. Determinism: the same input gives the same mesh in
// the same build.
// A contact or other boundary patch on part of a face (DECISIONS.md M9): a simple polygon of a
// tiling of the face. The polygons of a face's tiling: vertices exactly on the face's plane and in
// the face; conforming (every edge on the face's boundary or shared whole, in the opposite
// direction, by exactly one other polygon of the face); areas summing to the face's (1e-12). Their
// edges inside the face become edges of the mesh (split when encroached, as the domain edges), so
// every boundary triangle lies in one patch and each patch's area is exact. A patch vertex inside
// a domain edge is a node of both faces of the edge.
struct FacePatch {
    std::string name;  // a boundary patch of the result; unique among all patches
    std::vector<Point> vertices;
};

struct PolyhedronFace {
    std::vector<std::size_t> vertices;  // indices into ConvexPolyhedron::vertices, in order around the face
    std::string name;                   // the boundary patch of the result when `patches` is empty
    std::vector<FacePatch> patches;     // a tiling of the face into patches, or empty
};

struct ConvexPolyhedron {
    std::vector<Point> vertices;  // cm
    std::vector<PolyhedronFace> faces;
};

// A planar material interface (DECISIONS.md M7, U12; owner choice 2026-10-15: planar, no junction
// lines). The plane through three points (exact side tests); it must cut the whole domain. Rules
// (invalid_input otherwise, the index naming the interface): the points finite and not collinear;
// every polyhedron vertex at least 4 s from the plane, with vertices on both sides; every face it
// cuts met at 60 to 120 degrees, and every domain edge it crosses met by the face's trace at 60
// to 120 degrees within that face; each trace at least 6 s long; no cut face carrying patches
// (deferred); two planes never meeting in the domain and at least 4 s apart (so no line where three
// regions meet). Each slab between planes is one region, named by one seed point strictly inside
// it; the regions on the two sides of a plane differ.
// Nodes (M7): on the plane's polygon (the plane within the domain), inset by 1.5 s, a pair set P
// from the 2D generator at spacing min(h, 4 s), and for each p in P the nodes p +- (s/2) n; along
// each trace in a face, pairs at +- s/2 within the face, graded from each end (1.5 s, 3 s); on each
// crossed domain edge, nodes at +- s/2 from the crossing. The protected band b = sqrt(R^2 + s^2/4),
// R = h_P / sqrt 3: refinement inserts no circumcentre within b of the polygon, an encroached
// subfacet within b of a trace gets a new trace pair instead, and tetrahedra of interface nodes
// only are exempt from the quality bound. The result carries the P4 certificate (check_interfaces
// on the polygons, resolution s).
struct InterfacePlane {
    std::array<Point, 3> points;
};

struct RegionSeed {
    Point p;
    std::int32_t region = 0;  // >= 0
};

struct Geometry3D {
    ConvexPolyhedron domain;
    std::int32_t region = 0;                // with no interfaces: the region of every node, >= 0
    std::vector<InterfacePlane> interfaces;  // planar interfaces (6d)
    std::vector<RegionSeed> regions;         // with interfaces: one seed per slab
};

struct MeshOptions3D {
    double max_size = 0.0;                // h [cm], finite and positive
    std::vector<double> region_max_size;  // per region index; 0 or absent: max_size
    double interface_resolution = 0.0;    // s [cm]; required with interfaces
    std::vector<SizeBox> boxes;           // smaller h where given (lo and hi in all three coordinates)
    double radius_edge = 2.0;             // in [1.414, infinity)
    std::size_t max_nodes = 2'000'000;
};

// Errors: invalid_input (the geometry or option rules above); resource_exhausted (max_nodes
// reached, or refinement did not settle within its iteration cap); degenerate_mesh (the generated
// mesh failed the boundary check or an import check).
[[nodiscard]] std::expected<GeneratedMesh, base::Error> generate_mesh(const Geometry3D& geometry,
                                                                     const MeshOptions3D& options);


// The 3D generator over a polyhedral complex (Unit 53; DECISIONS.md FD3-FD5): non-convex domains.
//
// Input: a complex (complex.hpp, check_complex's rules). Output as generate_mesh(Geometry3D): a
// boundary-conforming Delaunay tetrahedralization, every boundary triangle Gabriel (P2) and inside
// one boundary facet, its box-method graph and the region of every node. Boundary patches are the
// facets' patch names, a patch over any number of facets (FD3); a node on an edge between two
// facets of one patch has the sum of its areas on both.
//
// Additional rules (invalid_input): the domain's boundary a manifold surface (every boundary edge
// in exactly two boundary facets).
//
// Interfaces (FD7-FD9): facets between volumes of different regions. A piece is a maximal
// connected set of coplanar interface facets of one region pair: a simple polygon, possibly
// non-convex, without holes. Its edges are feature lines (bends, junctions of three or more
// regions, traces where it meets the outer boundary) and their meeting points feature points.
// With interfaces, interface_resolution s is required, and (invalid_input, the index an interface
// facet): every sector about a feature line at least 60 degrees, also within a boundary face
// about a feature point; feature lines and traces at least 6 s long; unrelated features, pieces
// (thin layers), boundary faces a piece does not meet, and domain edges at least 4 s apart; every
// piece wide enough for its 1.5 s inset. Unsupported, and refused: holes, pinched pieces, a piece
// ending on a domain edge, two pieces ending on the boundary along one line, a trace at a domain
// corner, a feature line on a domain edge other than a single crossing, a patch edge meeting an
// interface (slice 4). Placement: each piece's pair set (the 2D generator's nodes in its 1.5 s
// inset at min(h of both regions, 4 s)) mirrored at +-s/2; pairs along each trace in the boundary
// face; a pair at each crossing of a domain edge, the domain-edge segment across it never split.
// Rows along every feature line sit on one lattice (multiples of 1.5 s along the line), so the
// pieces meeting there align. No node is inserted in a protected band; tetrahedra of interface
// nodes only, or with the circumcentre in a band, are exempt from the quality bound. The result
// carries the P4 certificate (check_interfaces; degenerate_mesh if it fails).
//
// Refinement as the convex generator, over the boundary's faces (maximal connected sets of
// coplanar boundary facets, possibly non-convex, with holes) and its domain edges (chains of facet
// edges where two faces meet) and patch edges (facet edges inside a face between different
// patches). A subsegment is also encroached when it is missing from the 2D triangulation of a
// face it bounds, so a non-convex face's boundary is recovered by splitting. Only tetrahedra
// inside the domain are refined. Termination at radius_edge = 2 for input angles of at least 90
// degrees (reflex included) is cited (Shewchuk 1998, FD5, to verify for non-convex input);
// otherwise heuristic, capped by max_nodes. Extraction: the tetrahedra a flood fill from the
// bounding tetrahedra cannot reach across subfacets; then the checks: every boundary triangle in
// one face, the volume and every face's area equal to the complex's within 1e-12, P0-P3 and P5.
//
// Errors: invalid_input (check_complex's rules, the options, the rules above); resource_exhausted
// (max_nodes, or the iteration cap); degenerate_mesh (a postcondition failed).
[[nodiscard]] std::expected<GeneratedMesh, base::Error> generate_mesh(const PolyhedralGeometry& geometry,
                                                                     const MeshOptions3D& options);

}  // namespace NiTCAD::mesh
