// The cells of a tensor-product mesh (Unit 20): what a path traced through the device needs and the
// node/edge graph (mesh.hpp) does not carry. A cell is the box between consecutive nodes along each
// axis, indexed by its lowest corner; a point in a cell has the multilinear stencil of the cell's
// 2^D corners (weights summing to 1, linear in the nodal values it interpolates). Node numbering
// is make_tensor_grid's: x fastest, then y, then z.
//
// An unstructured locator (Unit 16) is to provide the same operations for its own cells.
#pragma once

#include <array>
#include <cstddef>
#include <expected>
#include <optional>
#include <span>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/mesh/mesh.hpp"

namespace NiTCAD::mesh {

class TensorCells {
public:
    using Index = std::array<std::size_t, 3>;  // per axis; 0 beyond the dimension

    // The multilinear stencil of a point: count = 2^D corner nodes and their weights.
    struct Stencil {
        std::array<NodeId, 8> nodes{};
        std::array<double, 8> weights{};
        int count = 0;
    };

    // The cells of `mesh`, which must be make_tensor_grid of these axes (one per dimension).
    // Errors (invalid_input): the number of axes differs from the mesh's dimension; an axis with
    // fewer than two nodes, a non-finite value or not strictly increasing; the node count differs
    // from the product of the axis sizes; a node whose coordinates are not the axes' at its index
    // (the context index names the node).
    [[nodiscard]] static std::expected<TensorCells, base::Error> create(
        const Mesh& mesh, std::span<const double> x, std::span<const double> y = {},
        std::span<const double> z = {});

    // Whether these are the cells of `mesh`: its dimension, node count and every node's coordinates
    // at its index.
    [[nodiscard]] bool matches(const Mesh& mesh) const noexcept;

    [[nodiscard]] int dimension() const noexcept { return dimension_; }
    [[nodiscard]] std::span<const double> axis(int a) const noexcept { return axes_[a]; }
    // Nodes along axis a.
    [[nodiscard]] std::size_t size(int a) const noexcept { return axes_[a].size(); }
    [[nodiscard]] std::size_t stride(int a) const noexcept { return stride_[a]; }
    [[nodiscard]] NodeId node(const Index& index) const noexcept;
    // The index of a node. Precondition (NITCAD_EXPECTS): node in range.
    [[nodiscard]] Index index(NodeId node) const;
    // Whether a cell index is inside the mesh (each component below size - 1 on the dimension's
    // axes).
    [[nodiscard]] bool valid_cell(const Index& cell) const noexcept;
    // The multilinear stencil of point p in `cell` (p clamped into the cell). Corner k has bit a
    // of k set for the upper node along axis a. Precondition (NITCAD_EXPECTS): a valid cell.
    [[nodiscard]] Stencil stencil(const Index& cell, const Point& p) const;
    // The smallest node spacing of a cell along the dimension's axes [cm]; the largest cell
    // diagonal of the mesh [cm].
    [[nodiscard]] double smallest_width(const Index& cell) const;
    [[nodiscard]] double largest_diagonal() const noexcept { return largest_diagonal_; }

private:
    TensorCells() = default;

    // The first node whose coordinates are not the axes' at its index (node counts equal).
    [[nodiscard]] std::optional<std::size_t> mismatch(const Mesh& mesh) const noexcept;

    int dimension_ = 1;
    std::array<std::vector<double>, 3> axes_;
    std::array<std::size_t, 3> stride_{1, 1, 1};
    double largest_diagonal_ = 0.0;
};

}  // namespace NiTCAD::mesh
