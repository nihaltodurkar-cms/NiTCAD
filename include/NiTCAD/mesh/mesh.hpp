// Node/edge/control-volume mesh: the box-method graph every layer above uses (ARCHITECTURE.md
// principle 4 and 6.11).
//
// The public interface is the graph only:
// - nodes: coordinates and control volume;
// - edges: two end nodes, the edge length, and the coupling area (the measure of the dual
//   control-volume face the edge crosses), so a flux term is area * (u_second - u_first) / length;
// - boundary patches: named node sets with the boundary face area belonging to each node.
// There is no (i, j, k) indexing and no stride arithmetic: a tensor grid is one producer of this
// graph (tensor_grid.hpp), and an unstructured producer can be added without changing it.
// Edge directions come from the node coordinates.
//
// Units follow base/constants.hpp (cm). In D dimensions a volume is in cm^D and an area in
// cm^(D-1): a 1D mesh is per unit cross-section (area 1), a 2D mesh per unit depth.
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/base/error.hpp"

namespace NiTCAD::mesh {

using NodeId = std::int32_t;

// x, y, z in cm. Coordinates beyond the mesh dimension are 0.
using Point = std::array<double, 3>;

struct Edge {
    NodeId first;   // first < second
    NodeId second;
    double length;  // distance between the end nodes
    double coupling_area;
};

struct BoundaryPatch {
    std::string name;
    std::vector<NodeId> nodes;  // strictly increasing
    std::vector<double> areas;  // boundary face area per node, in parallel with nodes
};

class Mesh {
public:
    // Builds a mesh from its graph and validates it; every producer goes through here.
    // Errors:
    // - invalid_input: dimension not 1, 2 or 3; sizes that do not match; a non-finite coordinate
    //   or a nonzero coordinate beyond the dimension; an edge endpoint out of range, with
    //   first >= second, or listed twice; a boundary patch with an empty or repeated name, nodes
    //   out of range or not strictly increasing, or a size mismatch;
    // - degenerate_mesh: a control volume, edge length, coupling area or boundary area that is not
    //   positive and finite, or an edge length that differs from its endpoints' distance;
    // - resource_exhausted: more nodes than NodeId can index.
    // The context index names the offending node, edge or patch entry.
    [[nodiscard]] static std::expected<Mesh, base::Error> from_parts(
        int dimension, std::vector<Point> points, std::vector<double> volumes,
        std::vector<Edge> edges, std::vector<BoundaryPatch> boundary);

    [[nodiscard]] int dimension() const noexcept { return dimension_; }
    [[nodiscard]] std::size_t node_count() const noexcept { return points_.size(); }
    [[nodiscard]] std::span<const Point> points() const noexcept { return points_; }
    [[nodiscard]] std::span<const double> volumes() const noexcept { return volumes_; }
    [[nodiscard]] std::span<const Edge> edges() const noexcept { return edges_; }
    [[nodiscard]] std::span<const BoundaryPatch> boundary() const noexcept { return boundary_; }
    // The boundary patch with this name, or nullptr.
    [[nodiscard]] const BoundaryPatch* find_boundary(std::string_view name) const noexcept;

private:
    Mesh() = default;

    int dimension_ = 1;
    std::vector<Point> points_;
    std::vector<double> volumes_;
    std::vector<Edge> edges_;
    std::vector<BoundaryPatch> boundary_;
};

}  // namespace NiTCAD::mesh
