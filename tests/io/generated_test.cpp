// Generated meshes on disk (Unit 16; DECISIONS.md M14): simplices, regions, report and P4
// certificate round-trip; a plain mesh, a stale certificate and inconsistent simplices are refused.
#include <catch2/catch_test_macros.hpp>

#include <filesystem>
#include <vector>

#include "NiTCAD/io/files.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"

using namespace NiTCAD;
using base::ErrorCode;
namespace fs = std::filesystem;

namespace {

mesh::GeneratedMesh two_regions() {
    mesh::Geometry2D g{{{{{0, 0, 0}, {1e-4, 0, 0}, {1e-4, 5e-5, 0}, {0, 5e-5, 0}}, 0},
                        {{{0, 5e-5, 0}, {1e-4, 5e-5, 0}, {1e-4, 1e-4, 0}, {0, 1e-4, 0}}, 1}},
                       {{"bottom", {0, 0, 0}, {1e-4, 0, 0}}, {"top", {0, 1e-4, 0}, {1e-4, 1e-4, 0}}}};
    mesh::MeshOptions2D o;
    o.max_size = 1e-5;
    o.interface_resolution = 1e-6;
    auto m = mesh::generate_mesh(g, o);
    REQUIRE(m.has_value());
    return std::move(*m);
}

}  // namespace

TEST_CASE("io: a generated mesh round-trips with its simplices, regions and certificate", "[io][generate]") {
    const fs::path dir = fs::temp_directory_path() / "nitcad_io_generated_test";
    fs::create_directories(dir);
    const mesh::GeneratedMesh m = two_regions();
    REQUIRE(m.certificate.has_value());
    REQUIRE(io::save_generated_mesh(dir / "g.ntc", m).has_value());
    const auto back = io::load_generated_mesh(dir / "g.ntc");
    REQUIRE(back.has_value());
    CHECK(mesh::graph_digest(back->mesh) == mesh::graph_digest(m.mesh));
    CHECK(back->simplices.cells == m.simplices.cells);
    CHECK(back->simplices.facets == m.simplices.facets);
    CHECK(back->simplices.facet_patch == m.simplices.facet_patch);
    CHECK(back->node_region == m.node_region);
    REQUIRE(back->certificate.has_value());
    CHECK(back->certificate->mesh_digest == m.certificate->mesh_digest);
    CHECK(back->certificate->report.max_corner_deviation == m.certificate->report.max_corner_deviation);
    CHECK(back->report.nodes == m.report.nodes);
    CHECK(back->report.min_angle_deg == m.report.min_angle_deg);
    // The same file read as a plain mesh.
    const auto plain = io::load_mesh(dir / "g.ntc");
    REQUIRE(plain.has_value());
    CHECK(mesh::graph_digest(*plain) == mesh::graph_digest(m.mesh));
    fs::remove_all(dir);
}

TEST_CASE("io: generated-mesh refusals", "[io][generate]") {
    const fs::path dir = fs::temp_directory_path() / "nitcad_io_generated_refusals";
    fs::create_directories(dir);
    SECTION("a plain mesh file") {
        REQUIRE(io::save_mesh(dir / "plain.ntc", *mesh::make_tensor_grid(std::vector<double>{0, 1})).has_value());
        const auto r = io::load_generated_mesh(dir / "plain.ntc");
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().code == ErrorCode::invalid_input);
    }
    SECTION("a certificate stale for the stored regions") {
        mesh::GeneratedMesh m = two_regions();
        m.node_region[0] = 1 - m.node_region[0];
        REQUIRE(io::save_generated_mesh(dir / "stale.ntc", m).has_value());
        const auto r = io::load_generated_mesh(dir / "stale.ntc");
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().message.find("certificate") != std::string::npos);
    }
    SECTION("simplices that do not give the stored mesh") {
        mesh::GeneratedMesh m = two_regions();
        std::swap(m.simplices.cells[0], m.simplices.cells[1]);  // a cell reoriented: same mesh
        std::swap(m.simplices.cells[3], m.simplices.cells[6]);  // two cells' vertices exchanged
        REQUIRE(io::save_generated_mesh(dir / "bad.ntc", m).has_value());
        const auto r = io::load_generated_mesh(dir / "bad.ntc");
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().code == ErrorCode::invalid_input);
    }
    SECTION("valid simplices stored with a different mesh") {
        mesh::GeneratedMesh m = two_regions();
        std::vector<double> volumes(m.mesh.volumes().begin(), m.mesh.volumes().end());
        volumes[0] *= 2.0;  // a valid graph, but not the dual of the stored simplices
        auto other = mesh::Mesh::from_parts(2, {m.mesh.points().begin(), m.mesh.points().end()}, std::move(volumes),
                                            {m.mesh.edges().begin(), m.mesh.edges().end()},
                                            {m.mesh.boundary().begin(), m.mesh.boundary().end()});
        REQUIRE(other.has_value());
        const mesh::GeneratedMesh tampered{std::move(*other), m.simplices, m.node_region, std::nullopt, m.report};
        REQUIRE(io::save_generated_mesh(dir / "other.ntc", tampered).has_value());
        const auto r = io::load_generated_mesh(dir / "other.ntc");
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().message.find("simplices") != std::string::npos);
    }
    fs::remove_all(dir);
}
