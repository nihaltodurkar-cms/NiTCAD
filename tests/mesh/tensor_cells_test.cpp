// Tensor-grid cells (Unit 20): matching a mesh to its axes, node indices, and the multilinear
// stencil (corners, weights summing to 1, exact for a linear field).
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"

using namespace NiTCAD::mesh;
namespace base = NiTCAD::base;

TEST_CASE("tensor cells: matching the mesh") {
    const std::vector<double> x{0.0, 1.0, 3.0}, y{0.0, 0.5, 2.0, 2.5};
    const Mesh m = *make_tensor_grid(x, y);
    const auto c = TensorCells::create(m, x, y);
    REQUIRE(c.has_value());
    REQUIRE(c->dimension() == 2);
    REQUIRE(c->size(0) == 3);
    REQUIRE(c->size(1) == 4);
    REQUIRE(c->stride(1) == 3);
    for (NodeId n = 0; n < static_cast<NodeId>(m.node_count()); ++n) {
        const TensorCells::Index at = c->index(n);
        REQUIRE(c->node(at) == n);
        REQUIRE(m.points()[static_cast<std::size_t>(n)][0] == x[at[0]]);
        REQUIRE(m.points()[static_cast<std::size_t>(n)][1] == y[at[1]]);
    }
    REQUIRE(c->valid_cell({1, 2, 0}));
    REQUIRE_FALSE(c->valid_cell({2, 0, 0}));
    REQUIRE(c->smallest_width({1, 0, 0}) == 0.5);
    REQUIRE(c->largest_diagonal() == std::sqrt(2.0 * 2.0 + 1.5 * 1.5));

    // Mismatches.
    const auto code = [](const auto& r) { return r.error().code; };
    REQUIRE(code(TensorCells::create(m, x)) == base::ErrorCode::invalid_input);  // too few axes
    const std::vector<double> other{0.0, 1.0, 3.5};
    const auto moved = TensorCells::create(m, other, y);
    REQUIRE_FALSE(moved.has_value());
    REQUIRE(moved.error().context->index == 2u);  // the first node at x = 3
    const std::vector<double> shorter{0.0, 1.0};
    REQUIRE(code(TensorCells::create(m, shorter, y)) == base::ErrorCode::invalid_input);
    const std::vector<double> unsorted{0.0, 3.0, 1.0};
    REQUIRE(code(TensorCells::create(m, unsorted, y)) == base::ErrorCode::invalid_input);
}

TEST_CASE("tensor cells: the multilinear stencil") {
    const std::vector<double> x{0.0, 1.0, 3.0}, y{0.0, 0.5, 2.0}, z{0.0, 4.0};
    const Mesh m = *make_tensor_grid(x, y, z);
    const TensorCells c = *TensorCells::create(m, x, y, z);
    // A linear field a + b . r is reproduced exactly; the weights sum to 1.
    const auto field = [](const Point& p) { return 0.3 + 2.0 * p[0] - 1.5 * p[1] + 0.25 * p[2]; };
    const Point p{2.2, 1.1, 1.3};
    const TensorCells::Stencil s = c.stencil({1, 1, 0}, p);
    REQUIRE(s.count == 8);
    double sum = 0.0, value = 0.0;
    for (int k = 0; k < s.count; ++k) {
        const auto i = static_cast<std::size_t>(k);
        sum += s.weights[i];
        value += s.weights[i] * field(m.points()[static_cast<std::size_t>(s.nodes[i])]);
    }
    REQUIRE(std::abs(sum - 1.0) < 1e-15);
    REQUIRE(std::abs(value - field(p)) < 1e-14);
    // At a corner the stencil is that node alone.
    const TensorCells::Stencil corner = c.stencil({1, 1, 0}, {3.0, 0.5, 4.0});
    for (int k = 0; k < corner.count; ++k) {
        const auto i = static_cast<std::size_t>(k);
        const bool is = corner.nodes[i] == c.node({2, 1, 1});
        REQUIRE(corner.weights[i] == (is ? 1.0 : 0.0));
    }
}
