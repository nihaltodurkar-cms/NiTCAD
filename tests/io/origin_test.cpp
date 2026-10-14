// The mesh origin in files (DECISIONS.md U2), and the io digest against the mesh layer's.
#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <vector>

#include "NiTCAD/io/files.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"

using namespace NiTCAD;
namespace fs = std::filesystem;

TEST_CASE("io: the mesh origin round-trips; the digest equals mesh::graph_digest", "[io]") {
    const fs::path dir = fs::temp_directory_path() / "nitcad_io_origin_test";
    fs::create_directories(dir);
    const std::vector<double> x{0.0, 1e-6, 3e-6};
    const auto grid = mesh::make_tensor_grid(x, x);
    REQUIRE(grid.has_value());
    REQUIRE(io::save_mesh(dir / "grid.ntc", *grid).has_value());
    const auto back = io::load_mesh(dir / "grid.ntc");
    REQUIRE(back.has_value());
    CHECK(back->origin() == mesh::MeshOrigin::tensor_grid);
    const auto graph = mesh::Mesh::from_parts(2, {grid->points().begin(), grid->points().end()},
                                              {grid->volumes().begin(), grid->volumes().end()},
                                              {grid->edges().begin(), grid->edges().end()},
                                              {grid->boundary().begin(), grid->boundary().end()});
    REQUIRE(graph.has_value());
    REQUIRE(io::save_mesh(dir / "graph.ntc", *graph).has_value());
    const auto graph_back = io::load_mesh(dir / "graph.ntc");
    REQUIRE(graph_back.has_value());
    CHECK(graph_back->origin() == mesh::MeshOrigin::graph);
    // The io digest (unchanged since Unit 29) and the mesh layer's compute the same function.
    for (const mesh::Mesh& m : {*grid, *graph, *mesh::make_tensor_grid(x), *mesh::make_tensor_grid(x, x, x)}) {
        CHECK(io::mesh_digest(m) == mesh::graph_digest(m));
    }
    fs::remove_all(dir);
}
