// Point location on a mesh for state transfer (Unit 17, slice 3; DECISIONS.md Unit 17 AD5).
//
// A PointLocator answers, for a point, the nodes of the old mesh whose linear (simplex) or
// multilinear (tensor) interpolant holds the point, and their weights (summing to 1). Simplex
// meshes (2D triangles, 3D tetrahedra, as generate_mesh and generate_swept give them): exact
// location by orientation predicates through a bucket grid, the lowest-indexed containing cell on a
// shared facet; tensor meshes: the cell by binary search on each axis and TensorCells::stencil. A
// point outside every cell (rounding at an oblique boundary) takes its nearest node (lowest index on
// ties) with weight 1, and `inside` is false. Linear fields are reproduced exactly (to rounding)
// at every inside point of a simplex mesh, multilinear ones on a tensor mesh.
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <optional>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/mesh/simplex_mesh.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"

namespace NiTCAD::mesh {

struct Interpolant {
    std::array<NodeId, 8> nodes{};
    std::array<double, 8> weights{};
    int count = 0;
    bool inside = true;
};

class PointLocator {
public:
    // Errors (invalid_input): dimension not 2 or 3, no cells, a cell node out of range, a flat cell.
    [[nodiscard]] static std::expected<PointLocator, base::Error> create(const SimplexMesh& simplices);
    [[nodiscard]] static PointLocator create(const TensorCells& cells);

    [[nodiscard]] Interpolant locate(const Point& p) const;
    [[nodiscard]] std::size_t node_count() const noexcept { return nodes_; }

private:
    PointLocator() = default;
    [[nodiscard]] bool contains(std::size_t cell, const Point& p) const;
    [[nodiscard]] Interpolant weights(std::size_t cell, const Point& p) const;
    [[nodiscard]] Interpolant nearest(const Point& p) const;

    std::optional<TensorCells> cells_;  // a tensor mesh
    int dim_ = 2;
    std::size_t nodes_ = 0;
    std::vector<Point> points_;
    std::vector<NodeId> simplices_;  // positively oriented
    Point lo_{};
    double cell_ = 1.0;
    std::array<std::size_t, 3> n_{1, 1, 1};
    std::vector<std::uint32_t> start_;
    std::vector<std::uint32_t> items_;
};

}  // namespace NiTCAD::mesh
