// Material interfaces of a polyhedral complex (Unit 53; DECISIONS.md FD7): their pieces, feature
// lines and points, the rules, and the geometry the generator places pair nodes on. Internal to
// the mesh layer.
//
// A piece is a maximal connected union of coplanar interface facets of one region pair: a planar
// simple polygon (holes refused). Its edges are feature lines: a bend (two pieces of one pair), a
// junction (three or more regions), or a trace (the piece meets the outer boundary, which must be
// planar across it). Feature points: where feature lines meet, and where a trace crosses a domain
// edge.
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "complex_index.hpp"

namespace NiTCAD::mesh::detail {

// A domain edge (two boundary faces meeting) or a patch edge (two patches meeting inside a face):
// collinear complex vertices in order.
struct BoundaryChain {
    std::vector<std::size_t> chain;
    std::array<int, 2> faces{-1, -1};
    int face = -1;  // a patch edge: its face
};

struct InterfaceInput {
    const ComplexIndex* complex = nullptr;
    std::vector<int> face_of;                      // per complex facet: its boundary face, or -1
    std::vector<std::array<Point, 3>> face_plane;  // per boundary face: three non-collinear points
    std::vector<std::vector<std::array<Point, 3>>> face_triangles;  // per boundary face, in 3D
    std::vector<BoundaryChain> chains;
    std::vector<double> region_size;  // h per region index
    double resolution = 0.0;          // s
};

struct InterfacePieceLayout {
    std::vector<std::size_t> facets;
    std::int32_t region_a = 0;  // region_a < region_b
    std::int32_t region_b = 0;
    std::vector<std::size_t> loop;  // complex vertices around the piece, collinear ones included
    std::vector<Point> polygon;     // its corners (collinear vertices dropped), counterclockwise in the frame
    int axis = -1;                  // the coordinate constant on the piece, or -1 (oblique)
    Point origin{}, u{}, w{};       // oblique: an orthonormal frame in the plane
    Point normal{};                 // unit; (u, w, normal) right-handed for oblique pieces
    std::vector<Point> local;       // the corners in the frame
    std::vector<Point> inset;       // the polygon inset by 1.5 s, in the frame
    double spacing = 0.0;           // pair spacing: min(h of both regions, 4 s)

    [[nodiscard]] Point to_local(const Point& p) const;
    // The point at frame coordinates l, offset along the normal (exact on the axis for axis pieces).
    [[nodiscard]] Point to_world(const Point& l, double offset) const;
    // Distance to the polygon (closed), in double.
    [[nodiscard]] double distance(const Point& x) const;
};

struct InterfaceTrace {
    int piece = 0;
    int face = 0;
    Point start{};
    Point dir{};     // unit, along the trace
    Point normal{};  // unit, in the face, across the trace
    double length = 0.0;
};

// A trace crossing a domain edge at chains[chain].chain[position] (never inserted as a node).
struct InterfaceCrossing {
    int chain = 0;
    std::size_t position = 0;
};

struct InterfaceLayout {
    std::vector<InterfacePieceLayout> pieces;
    std::vector<InterfaceTrace> traces;
    std::vector<InterfaceCrossing> crossings;
    std::vector<char> on_interface;  // per complex vertex: on an interface facet
    double band = 0.0;               // the protected band b = sqrt(R^2 + s^2 / 4), largest over pieces
};

// The layout and the FD7 rules (invalid_input, the index a facet of the piece concerned). Without
// interface facets the layout is empty and the resolution is not read.
[[nodiscard]] std::expected<InterfaceLayout, base::Error> read_interfaces(const InterfaceInput& in);

// The pair set P of a piece in its frame: the 2D generator's nodes on the inset polygon at the
// piece's spacing, with boundary nodes on each inset edge at the lattice of the feature line beside it.
[[nodiscard]] std::expected<std::vector<Point>, base::Error> pair_set(const InterfacePieceLayout& piece, double s,
                                                                      std::size_t max_nodes);

// The distance from x to the closed triangle abc, in double.
[[nodiscard]] double triangle_distance(const Point& x, const Point& a, const Point& b, const Point& c);

// Pair positions along a trace (from start, along unit dir): 1.5 s from each end and at the trace's
// line lattice (multiples of 1.5 s in the line's coordinate) between, as the pieces' rows beside it.
[[nodiscard]] std::vector<double> trace_positions(const Point& start, const Point& dir, double length, double s);

}  // namespace NiTCAD::mesh::detail
