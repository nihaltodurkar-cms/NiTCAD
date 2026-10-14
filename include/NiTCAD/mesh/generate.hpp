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

#include <cstddef>
#include <cstdint>
#include <expected>
#include <optional>
#include <string>
#include <vector>

#include "NiTCAD/base/error.hpp"
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
// ends on the segment are added to `along`.
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

}  // namespace NiTCAD::mesh
