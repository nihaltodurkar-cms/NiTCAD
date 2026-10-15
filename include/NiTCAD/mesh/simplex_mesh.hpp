// Simplicial meshes and their box-method graph (Unit 16; ARCHITECTURE.md 6.11; DECISIONS.md M1-M4).
//
// A SimplexMesh is triangles (2D) or tetrahedra (3D) on points, with named boundary facets.
// make_dual_mesh turns it into the node/edge/control-volume graph (mesh.hpp) on the circumcentric
// dual, after checking that the dual is the Voronoi diagram clipped to the domain:
// - P0: no flat element (exact orientation; a negatively oriented element is reoriented);
// - P1: locally Delaunay: across every interior facet the opposite vertex is not strictly inside
//   the circumcircle or circumsphere (exact; a tie, a cocircular or cospherical vertex, is allowed);
// - P2: every boundary facet, listed or not, is Gabriel: no node strictly inside its diametral
//   circle (2D) or the smallest sphere through it (3D) (exact);
// - P3: totals non-negative (M2), with tolerance tau = 1e-12 on a scale of the same dimension as
//   the total: an edge's coupling area with s = max(L^(D-1), sum of |pieces|); a node's volume
//   with s = max(l^D, sum of |pieces|), l the mean length of the node's edges; a node's area in a
//   boundary patch with s = max(l^(D-1), sum of |pieces|). The sum of |pieces| covers cancellation
//   among signed pieces. A coupling or patch area within tau s of zero is zero: the edge, or the
//   node's membership of the patch, is dropped (a cocircular diagonal, for example). Below -tau s
//   it is refused, as is a volume below tau s;
// - P5: every control volume closes: sum over its edges of A_e t_e (outward) plus sum over its
//   boundary pieces of a n = 0 to 1e-11 of the sum of the magnitudes.
// P1 and P2 are what make every total non-negative (M3), so P3 and P5 refuse only meshes that P1
// and P2 accept through rounding: they guard the computation, not the input.
//
// Boundary facets not listed in any patch are natural (zero-flux) boundaries, still checked by P2.
// The result is deterministic for a given build: the same input gives the same mesh. Every
// topological decision is an exact predicate, so it does not depend on the platform; the measures
// are floating-point and bit identity across compilers or machines is not claimed.
//
// Not checked: global overlap of a mesh that is locally valid everywhere (a fold onto itself far
// from any shared facet), and the interface rule P4 (an imported mesh carries no interfaces).
#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <string>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/mesh/mesh.hpp"

namespace NiTCAD::mesh {

struct SimplexMesh {
    int dimension = 2;                      // 2 or 3
    std::vector<Point> points;              // cm; coordinates beyond the dimension are 0
    std::vector<NodeId> cells;              // dimension + 1 node ids per cell, flattened
    std::vector<NodeId> facets;             // dimension node ids per listed boundary facet
    std::vector<std::int32_t> facet_patch;  // per listed facet, an index into patch_names
    std::vector<std::string> patch_names;   // the boundary patches of the result, in this order
};

// What make_dual_mesh did besides building the graph.
struct DualMeshReport {
    std::size_t dropped_edges = 0;        // couplings within tau of zero
    std::size_t dropped_patch_nodes = 0;  // patch areas within tau of zero
    std::size_t boundary_facets = 0;      // topological boundary facets, listed or not
    double max_closure_residual = 0.0;    // largest P5 residual, relative
};

// Errors (the context index names the entry the message names):
// - invalid_input: dimension not 2 or 3; no points or no cells; cells, facets or facet_patch sized
//   inconsistently; a coordinate not finite or nonzero beyond the dimension (node); a cell node out
//   of range or repeated in the cell (cell); the same cell listed twice (cell); a node in no cell
//   (node); a patch name empty or repeated (patch); a facet node out of range or repeated, a patch
//   index out of range, a facet listed twice, or a listed facet that is not a boundary facet (a
//   face of exactly one cell) (facet);
// - degenerate_mesh: a flat cell (P0, cell); a facet shared by more than two cells, or two cells
//   on the same side of their shared facet (cell); not locally Delaunay (P1, cell); a boundary
//   facet with a node strictly inside its diametral circle or sphere (P2, the cell owning the
//   facet; the node in the message); a negative coupling (P3, the edge's first node; both in the
//   message), a volume not positive (P3, node), a negative patch area (P3, node), a control volume
//   that does not close (P5, node);
// - resource_exhausted: more points than NodeId can index;
// and those of Mesh::from_parts.
[[nodiscard]] std::expected<Mesh, base::Error> make_dual_mesh(const SimplexMesh& simplices,
                                                              DualMeshReport* report = nullptr);

}  // namespace NiTCAD::mesh
