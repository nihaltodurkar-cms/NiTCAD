// What the views draw on: a mesh's linear elements, its boundary surface, slices through it, the
// control-volume pieces that carry each node's region, and contour lines (Unit 25a; DECISIONS.md VR2,
// VR3, VG2).
//
// A Mesh (mesh.hpp) is a node and edge graph without cells, so the elements come from the producer:
// - a tensor mesh's TensorCells: each cell is split into 2 triangles (2D) or 6 tetrahedra (3D, the
//   Kuhn split along the cell's main diagonal, conforming between neighbouring cells);
// - a generated or imported mesh's SimplexMesh: its triangles or tetrahedra;
// - otherwise the graph alone: a 1D chain's segments, or (2D, 3D) nodes and edges only.
// Fields are interpolated linearly on these simplices.
//
// Regions, materials and doping are per node (device.hpp); in the box method a node stands for its
// control volume, and the interface between two regions runs between their nodes. The views draw a
// node's region on the part of each element nearest to it:
// - in a tensor cell, the sub-box of the corner (the cell split at its mid-planes): exactly the
//   control volumes of the box method;
// - in a simplex, the points whose largest barycentric weight is the node's (the barycentric dual,
//   bounded by the planes of equal weight).
// Pieces are found by clipping convex polygons with these planes, so they are exact up to rounding.
//
// Coordinates are the mesh's (cm). A Geometry refers to its mesh, which must outlive it.
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/mesh/simplex_mesh.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"

namespace NiTCAD::render {

using mesh::NodeId;
using mesh::Point;

enum class GeometryKind : std::uint8_t {
    graph,    // nodes and edges (a 1D chain has segments as elements)
    tensor,   // TensorCells: boxes split into simplices
    simplex,  // a SimplexMesh's cells
};

// A point as a linear combination of up to four mesh nodes (unused slots weigh 0), so any nodal field
// can be read at it: value = sum of weights[k] * field[nodes[k]].
struct Interpolant {
    std::array<NodeId, 4> nodes{};
    std::array<double, 4> weights{};
};

// A vertex of a drawn polygon: its position and how to read fields there.
struct SurfaceVertex {
    Point position{};
    Interpolant at;
};

// A convex planar polygon (3 or more vertices, in order) inside one element.
struct Polygon {
    std::vector<SurfaceVertex> vertices;
    std::size_t element = 0;  // the simplex it lies in
};

// A polygon in the control volume of one node.
struct OwnedPolygon {
    std::vector<Point> points;
    NodeId owner = 0;
};

// A segment (2D) or polygon (3D) on which the owners of different regions meet.
struct InterfacePiece {
    std::vector<Point> points;  // 2 in 2D, 3 or more in 3D
    std::int32_t region_a = 0;  // region_a < region_b
    std::int32_t region_b = 0;
    NodeId node_a = 0;          // the nodes on either side, in region_a and region_b
    NodeId node_b = 0;
    std::size_t element = 0;    // an element it lies in (tensor: the first of its cell's)
};

// A plane n . x = d (n need not be unit length).
struct Plane {
    Point normal{};
    double offset = 0.0;
};

// An axis-aligned box, lo <= hi on each axis (cm): the part a clipping box removes (Unit 25c, DC2).
struct Box {
    Point lo{};
    Point hi{};
};

// The parts of a convex planar polygon outside a box, each convex: the polygon is cut by the box's six
// planes in turn, the part beyond each kept and the rest passed on (what remains is inside and
// dropped). Parts of zero area are dropped.
[[nodiscard]] std::vector<std::vector<Point>> outside_box(std::span<const Point> polygon, const Box& box);

class Geometry {
public:
    // Errors (invalid_input): cells that do not match the mesh (TensorCells::matches).
    [[nodiscard]] static std::expected<Geometry, base::Error> tensor(const mesh::Mesh& mesh,
                                                                    const mesh::TensorCells& cells);
    // Errors (invalid_input): a dimension, point count or point differing from the mesh's; a cell
    // node out of range.
    [[nodiscard]] static std::expected<Geometry, base::Error> simplex(
        const mesh::Mesh& mesh, const mesh::SimplexMesh& simplices);
    // The graph alone. In 1D the chain's edges are its elements.
    [[nodiscard]] static Geometry graph(const mesh::Mesh& mesh);

    [[nodiscard]] GeometryKind kind() const noexcept { return kind_; }
    [[nodiscard]] int dimension() const noexcept { return mesh_->dimension(); }
    [[nodiscard]] const mesh::Mesh& mesh() const noexcept { return *mesh_; }

    // The linear elements: segments (1D), triangles (2D) or tetrahedra (3D), dimension + 1 nodes
    // each, positively oriented in 2D and 3D. Empty for a 2D or 3D graph.
    [[nodiscard]] std::size_t element_count() const noexcept;
    [[nodiscard]] std::span<const NodeId> element(std::size_t e) const noexcept;
    // The tensor cell an element was split from (tensor geometry only).
    [[nodiscard]] std::size_t element_cell(std::size_t e) const noexcept { return element_cell_[e]; }
    // A tensor cell's lower and upper corner.
    [[nodiscard]] std::size_t cell_count() const noexcept { return cell_lo_.size(); }
    [[nodiscard]] std::pair<Point, Point> cell_box(std::size_t c) const noexcept;
    // A tensor cell's 2^D corner nodes; corner k has bit a set for the upper node along axis a.
    [[nodiscard]] std::span<const NodeId> cell_corners(std::size_t c) const noexcept;

    // The element's polygon as drawn: 2D, the triangle itself.
    [[nodiscard]] Polygon element_polygon(std::size_t e) const;
    // 3D: the faces of the tetrahedra on the domain boundary, oriented outward (counterclockwise
    // seen from outside); 2D: every triangle. Each as a triangle Polygon.
    [[nodiscard]] std::vector<Polygon> boundary_faces() const;
    // The section of every element by a plane (2D and 3D; 3D gives triangles and quadrilaterals,
    // 2D segments as two-vertex polygons). Elements touching the plane only at a vertex, an edge
    // or (3D) lying in it on a face give no piece of zero area.
    [[nodiscard]] std::vector<Polygon> slice(const Plane& plane) const;
    // 3D, the surface of the device with the box removed (Unit 25c, DC2): the boundary faces' parts
    // outside the box, and the box's faces inside the device (sections of the elements by each face's
    // plane, clipped to the face), all oriented outward from what is left (the cut faces toward the
    // box), each a convex polygon in its element. A box face on or beyond the mesh's bounding box makes
    // no cut face. Empty for a geometry that is not 3D.
    [[nodiscard]] std::vector<Polygon> clipped_surface(const Box& box) const;
    // The parts of a polygon inside element `p.element` nearest to each of the element's nodes, by
    // the rule above. Pieces of zero area are dropped.
    [[nodiscard]] std::vector<OwnedPolygon> owner_pieces(const Polygon& p) const;
    // How to read fields at a point of element e: its barycentric weights there. For a tensor
    // geometry the point may lie anywhere in e's cell, and the cell's element holding it is used
    // (the nearest, when rounding puts the point just outside all of them).
    [[nodiscard]] Interpolant interpolant(std::size_t e, const Point& x) const;
    // Where the control volumes of nodes in different regions meet inside the elements (2D segments,
    // 3D polygons). node_region has one entry per mesh node.
    [[nodiscard]] std::vector<InterfacePiece> interfaces(std::span<const std::int32_t> node_region) const;

private:
    Geometry() = default;
    void split_tensor(const mesh::TensorCells& cells);

    const mesh::Mesh* mesh_ = nullptr;
    GeometryKind kind_ = GeometryKind::graph;
    int per_element_ = 0;                 // nodes per element, 0 when none
    std::vector<NodeId> elements_;
    std::vector<std::size_t> element_cell_;
    std::vector<Point> cell_lo_;
    std::vector<Point> cell_hi_;
    std::vector<NodeId> cell_corners_;    // 2^D per cell
};

// The value of a nodal field at an interpolant.
[[nodiscard]] double interpolate(const Interpolant& at, std::span<const double> field) noexcept;

// Marching triangles: where a field linear on the triangle (values at its three vertices) equals
// `level`, as a segment, or nothing. A vertex exactly at the level counts as above it, so a contour
// through a shared vertex is drawn once.
[[nodiscard]] std::optional<std::array<Point, 2>> contour_segment(const std::array<Point, 3>& p,
                                                                  const std::array<double, 3>& v,
                                                                  double level) noexcept;

// An isosurface as an indexed triangle mesh (Unit 25b, DS4): each mesh edge the surface crosses is cut
// once, at one vertex shared by every element around the edge (a crossing exactly at a node is that
// node's vertex), so the surface is watertight where it does not leave the domain. Triangles are
// counterclockwise seen from the side of higher values.
struct IsoSurface {
    std::vector<SurfaceVertex> vertices;
    std::vector<std::array<std::uint32_t, 3>> triangles;
};

// Marching tetrahedra: where the field, linear on each tetrahedron, equals `level`. As in
// contour_segment a node at the level counts as above it; elements with a NaN value are skipped,
// and triangles of zero area (two vertices at one node) dropped. Empty for a geometry that is not
// 3D or has no elements. `values` has one entry per mesh node.
[[nodiscard]] IsoSurface isosurface(const Geometry& geometry, std::span<const double> values,
                                    double level);

// Area of a planar polygon (any orientation) and the signed volume enclosed by triangles (positive
// for an outward-oriented closed surface).
[[nodiscard]] double polygon_area(std::span<const Point> points) noexcept;
[[nodiscard]] double enclosed_volume(std::span<const Polygon> triangles) noexcept;

}  // namespace NiTCAD::render
