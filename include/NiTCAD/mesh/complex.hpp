// The polyhedral complex and the sectioned extrusion (Unit 53; DECISIONS.md FD1, FD2, FD6): the
// geometry input of the non-convex 3D generator.
//
// A complex is vertices, planar polygonal facets over them, and volumes, each a closed polyhedron
// given by its facets. Every volume has a region index; several volumes may share one. A facet
// lies on one volume (an outer boundary facet, which names its boundary patch) or on two (an
// internal facet): between volumes of different regions it is a material interface, between
// volumes of one region it is not a boundary of anything and plays no further part.
//
// Rules (invalid_input; the index names the vertex, facet or volume), all decided exactly:
// - vertices finite and pairwise distinct;
// - each facet: at least three vertex indices in range, none repeated; planar (every vertex on the
//   plane of its first three non-collinear ones, by orient3d); simple in its plane (no two edges
//   meeting except consecutive ones at their shared vertex; consecutive collinear vertices are
//   allowed, a reversal is not);
// - each facet on one or two volumes; a boundary facet's patch name non-empty;
// - each volume: its facets distinct, every edge of its facets in exactly two of them (so facets
//   meet edge to edge, never at a vertex inside another facet's edge), consistently orientable, the
//   facets around each vertex one cycle (no pinched vertices), of non-zero volume; region >= 0;
// - facets meet only along shared edges or at shared vertices (every pair of facet triangles
//   intersects in exactly their shared vertices or edge);
// - the volumes connected through shared facets.
#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <string>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/mesh/mesh.hpp"

namespace NiTCAD::mesh {

struct ComplexFacet {
    std::vector<std::size_t> vertices;  // indices into PolyhedralGeometry::vertices, in order around the facet
    std::string patch;                  // a boundary facet's patch; ignored for internal facets
};

struct ComplexVolume {
    std::vector<std::size_t> facets;  // indices into PolyhedralGeometry::facets: a closed polyhedron
    std::int32_t region = 0;          // >= 0
};

struct PolyhedralGeometry {
    std::vector<Point> vertices;  // cm
    std::vector<ComplexFacet> facets;
    std::vector<ComplexVolume> volumes;
};

struct ComplexReport {
    std::size_t boundary_facets = 0;
    std::size_t interface_facets = 0;  // internal facets between different regions
    std::size_t internal_facets = 0;   // internal facets within one region
    std::vector<double> volume;        // per volume [cm^3] (double, for reporting)
    std::vector<std::string> patches;  // the boundary patches, in order of first appearance
};

// Validates the complex against the rules above. Errors: invalid_input.
[[nodiscard]] std::expected<ComplexReport, base::Error> check_complex(const PolyhedralGeometry& geometry);


// The sectioned extrusion (DECISIONS.md FD2; owner ruling 2026-10-29): a cross-section tiled into
// named cells, swept along a coordinate axis through stations t_0 < ... < t_k; in each section
// (t_i, t_i+1) each cell holds a region or is void (outside the domain). Cross-section coordinates
// (u, v) follow the axis cyclically: axis x gives (y, z), axis y (z, x), axis z (x, y); a cell
// polygon's points carry u and v in their first two components (the third is ignored).
//
// The complex built: one volume per non-void (cell, section), its region the assignment; the
// facets are cell edges x sections (lateral) and cells at stations, each emitted once with the one
// or two volumes beside it. Every vertex is a cell vertex at a station, so nothing is computed by
// intersection. A boundary facet takes the patch of the first patch box containing it (closed, all
// three coordinates), else `default_patch`.
//
// Rules (invalid_input; the index names the cell, station or patch box): axis 0, 1 or 2; at least
// two stations, finite and strictly increasing; each cell a simple polygon of non-zero area with a
// unique non-empty name; the cells a conforming tiling (they meet only along whole shared edges or
// at shared vertices, bit-identical coordinates; no overlap, no nesting); `regions` one row per
// cell with one entry per section, each >= 0 or `void_region`; at least one non-void entry; patch
// box names non-empty and not `default_patch`. The complex then meets check_complex's rules
// (connectedness, manifold volumes) or the error is check_complex's.
struct ExtrusionCell {
    std::string name;
    std::vector<Point> polygon;  // (u, v, ignored), either orientation
};

struct PatchBox {
    std::string name;
    Point lo;
    Point hi;
};

struct ExtrudedGeometry {
    static constexpr std::int32_t void_region = -1;
    int axis = 0;
    std::vector<double> stations;  // cm
    std::vector<ExtrusionCell> cells;
    std::vector<std::vector<std::int32_t>> regions;  // [cell][section]
    std::vector<PatchBox> patches;                   // boundary facets inside a box take its name
    std::string default_patch = "boundary";
};

// Builds and checks the complex. Errors: invalid_input.
[[nodiscard]] std::expected<PolyhedralGeometry, base::Error> extrude(const ExtrudedGeometry& geometry);

}  // namespace NiTCAD::mesh
