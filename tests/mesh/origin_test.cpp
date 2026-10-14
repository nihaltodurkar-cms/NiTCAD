// The mesh origin tag (DECISIONS.md U2): set by make_tensor_grid, verified by Mesh::from_parts,
// not part of the graph digest.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"

using namespace NiTCAD::mesh;

TEST_CASE("mesh origin: the tensor-grid tag is verified and kept out of the digest", "[mesh]") {
    const std::vector<double> x{0.0, 1.0, 3.0};
    const std::vector<double> y{0.0, 2.0};
    const auto grid = make_tensor_grid(x, y);
    REQUIRE(grid.has_value());
    CHECK(grid->origin() == MeshOrigin::tensor_grid);
    CHECK(make_tensor_grid(x)->origin() == MeshOrigin::tensor_grid);
    CHECK(make_tensor_grid(x, y, x)->origin() == MeshOrigin::tensor_grid);
    const auto copy = [&](MeshOrigin o) {
        return Mesh::from_parts(2, {grid->points().begin(), grid->points().end()},
                                {grid->volumes().begin(), grid->volumes().end()},
                                {grid->edges().begin(), grid->edges().end()},
                                {grid->boundary().begin(), grid->boundary().end()}, o);
    };
    const auto graph = copy(MeshOrigin::graph);
    REQUIRE(graph.has_value());
    CHECK(graph->origin() == MeshOrigin::graph);
    CHECK(graph_digest(*graph) == graph_digest(*grid));  // the origin is not the identity
    CHECK(copy(MeshOrigin::tensor_grid).has_value());
    // A diagonal edge under the tensor tag is refused; as a general graph it is fine.
    const std::vector<Point> p{{0, 0, 0}, {1, 1, 0}};
    const std::vector<Edge> diagonal{{0, 1, std::sqrt(2.0), 1.0}};
    const auto tagged = Mesh::from_parts(2, p, {0.5, 0.5}, diagonal, {}, MeshOrigin::tensor_grid);
    REQUIRE_FALSE(tagged.has_value());
    CHECK(tagged.error().code == NiTCAD::base::ErrorCode::invalid_input);
    CHECK(Mesh::from_parts(2, p, {0.5, 0.5}, diagonal, {}).has_value());
}
