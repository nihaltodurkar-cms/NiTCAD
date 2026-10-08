#include "NiTCAD/mesh/tensor_cells.hpp"

#include <algorithm>
#include <cmath>
#include <optional>
#include <string>
#include <utility>

#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::mesh {

namespace {

base::Error invalid(std::string message, std::optional<std::size_t> index = std::nullopt) {
    return {base::ErrorCode::invalid_input, std::move(message),
            base::ErrorContext{.index = index, .value = std::nullopt}};
}

}  // namespace

std::expected<TensorCells, base::Error> TensorCells::create(const Mesh& mesh,
                                                            std::span<const double> x,
                                                            std::span<const double> y,
                                                            std::span<const double> z) {
    const std::span<const double> given[3] = {x, y, z};
    int dimension = 0;
    while (dimension < 3 && !given[dimension].empty()) ++dimension;
    for (int a = dimension; a < 3; ++a) {
        if (!given[a].empty()) return std::unexpected(invalid("an axis follows an empty one"));
    }
    if (dimension != mesh.dimension()) {
        return std::unexpected(invalid("the number of axes differs from the mesh's dimension"));
    }
    TensorCells c;
    c.dimension_ = dimension;
    std::size_t count = 1;
    for (int a = 0; a < dimension; ++a) {
        const std::span<const double> axis = given[a];
        if (axis.size() < 2) return std::unexpected(invalid("an axis needs at least two nodes"));
        for (std::size_t k = 0; k < axis.size(); ++k) {
            if (!std::isfinite(axis[k]) || (k > 0 && !(axis[k] > axis[k - 1]))) {
                return std::unexpected(invalid("an axis is not finite and strictly increasing", k));
            }
        }
        c.axes_[a].assign(axis.begin(), axis.end());
        c.stride_[a] = count;
        count *= axis.size();
    }
    for (int a = dimension; a < 3; ++a) {
        c.axes_[a] = {0.0};
        c.stride_[a] = count;
    }
    if (count != mesh.node_count()) {
        return std::unexpected(invalid("the node count differs from the axes' product"));
    }
    if (const auto bad = c.mismatch(mesh)) {
        return std::unexpected(
            invalid("a node's coordinates are not the axes' at its index", *bad));
    }
    double diagonal = 0.0;
    for (int a = 0; a < dimension; ++a) {
        double widest = 0.0;
        for (std::size_t k = 0; k + 1 < c.axes_[a].size(); ++k) {
            widest = std::max(widest, c.axes_[a][k + 1] - c.axes_[a][k]);
        }
        diagonal += widest * widest;
    }
    c.largest_diagonal_ = std::sqrt(diagonal);
    return c;
}

bool TensorCells::matches(const Mesh& mesh) const noexcept {
    std::size_t count = 1;
    for (int a = 0; a < dimension_; ++a) count *= axes_[a].size();
    return mesh.dimension() == dimension_ && mesh.node_count() == count && !mismatch(mesh);
}

std::optional<std::size_t> TensorCells::mismatch(const Mesh& mesh) const noexcept {
    for (std::size_t node = 0; node < mesh.node_count(); ++node) {
        const Point& p = mesh.points()[node];
        for (int a = 0; a < dimension_; ++a) {
            if (p[a] != axes_[a][(node / stride_[a]) % axes_[a].size()]) return node;
        }
    }
    return std::nullopt;
}

NodeId TensorCells::node(const Index& index) const noexcept {
    std::size_t n = 0;
    for (int a = 0; a < dimension_; ++a) n += index[a] * stride_[a];
    return static_cast<NodeId>(n);
}

TensorCells::Index TensorCells::index(NodeId node) const {
    const auto n = static_cast<std::size_t>(node);
    NITCAD_EXPECTS(node >= 0 && n < stride_[dimension_ - 1] * axes_[dimension_ - 1].size());
    Index at{};
    for (int a = 0; a < dimension_; ++a) at[a] = (n / stride_[a]) % axes_[a].size();
    return at;
}

bool TensorCells::valid_cell(const Index& cell) const noexcept {
    for (int a = 0; a < dimension_; ++a) {
        if (cell[a] + 1 >= axes_[a].size()) return false;
    }
    return true;
}

TensorCells::Stencil TensorCells::stencil(const Index& cell, const Point& p) const {
    NITCAD_EXPECTS(valid_cell(cell));
    double t[3] = {};  // the point's fraction along each axis of the cell
    for (int a = 0; a < dimension_; ++a) {
        const double lo = axes_[a][cell[a]], hi = axes_[a][cell[a] + 1];
        t[a] = std::clamp((p[a] - lo) / (hi - lo), 0.0, 1.0);
    }
    Stencil s;
    s.count = 1 << dimension_;
    for (int k = 0; k < s.count; ++k) {
        Index corner = cell;
        double w = 1.0;
        for (int a = 0; a < dimension_; ++a) {
            const bool upper = ((k >> a) & 1) != 0;
            corner[a] += upper ? 1 : 0;
            w *= upper ? t[a] : 1.0 - t[a];
        }
        s.nodes[static_cast<std::size_t>(k)] = node(corner);
        s.weights[static_cast<std::size_t>(k)] = w;
    }
    return s;
}

double TensorCells::smallest_width(const Index& cell) const {
    NITCAD_EXPECTS(valid_cell(cell));
    double w = axes_[0][cell[0] + 1] - axes_[0][cell[0]];
    for (int a = 1; a < dimension_; ++a) {
        w = std::min(w, axes_[a][cell[a] + 1] - axes_[a][cell[a]]);
    }
    return w;
}

}  // namespace NiTCAD::mesh
