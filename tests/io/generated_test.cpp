// Generated meshes on disk (Unit 16; DECISIONS.md M14): simplices, regions, report and P4
// certificate round-trip; a plain mesh, a stale certificate and inconsistent simplices are refused.
#include <catch2/catch_test_macros.hpp>

#include <cstdint>
#include <filesystem>
#include <vector>

#include "NiTCAD/io/files.hpp"
#include "NiTCAD/mesh/complex.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/interfaces.hpp"
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

TEST_CASE("io: a complex-generated mesh round-trips with its corner patches and P6 report (Unit 53, FD14)",
          "[io][generate]") {
    const fs::path dir = fs::temp_directory_path() / "nitcad_io_generated_complex_test";
    fs::create_directories(dir);
    // Two slabs along x, 2 wide and 2 high, the interface at z = 1; the patch "low" is the side face below it,
    // so its edge lies on the interface trace and boundary triangles straddle it (FD8 a).
    mesh::ExtrudedGeometry x;
    x.axis = 0;
    x.stations = {0.0, 1.0};
    const double v[3] = {0.0, 1.0, 2.0};
    for (int j = 0; j < 2; ++j) {
        x.cells.push_back({"c" + std::to_string(j),
                           {{0.0, v[j], 0.0}, {2.0, v[j], 0.0}, {2.0, v[j + 1], 0.0}, {0.0, v[j + 1], 0.0}}});
        x.regions.push_back({j});
    }
    x.patches = {{"low", {-1.0, 0.0, -1.0}, {3.0, 0.0, 1.0}}};
    const auto geometry = mesh::extrude(x);
    REQUIRE(geometry.has_value());
    mesh::MeshOptions3D o;
    o.max_size = 0.5;
    o.interface_resolution = 0.1;
    auto generated = mesh::generate_mesh(*geometry, o);
    REQUIRE(generated.has_value());
    const mesh::GeneratedMesh& m = *generated;
    REQUIRE_FALSE(m.simplices.facet_corner_patch.empty());
    REQUIRE(m.certificate.has_value());
    REQUIRE(m.certificate->patches.has_value());

    REQUIRE(io::save_generated_mesh(dir / "g.ntc", m).has_value());
    const auto back = io::load_generated_mesh(dir / "g.ntc");
    if (!back) FAIL(back.error().message);
    CHECK(back->simplices.facet_corner_patch == m.simplices.facet_corner_patch);
    CHECK(mesh::graph_digest(back->mesh) == mesh::graph_digest(m.mesh));
    REQUIRE(back->certificate.has_value());
    REQUIRE(back->certificate->patches.has_value());
    const mesh::PatchBoundaryReport& a = *m.certificate->patches;
    const mesh::PatchBoundaryReport& b = *back->certificate->patches;
    CHECK(a.cross_pieces == b.cross_pieces);
    CHECK(a.corner_pieces == b.corner_pieces);
    CHECK(a.max_flat_deviation == b.max_flat_deviation);
    CHECK(a.max_corner_deviation == b.max_corner_deviation);
    CHECK(a.area_error == b.area_error);
    CHECK(a.zone_length == b.zone_length);
    CHECK(mesh::certifies(*back->certificate, back->mesh, back->node_region));
    // The generation digest (FD14) is the complex's and options', and round-trips.
    CHECK(m.generation_digest == mesh::generation_digest(*geometry, o));
    CHECK(m.generation_digest != 0);
    CHECK(back->generation_digest == m.generation_digest);
    fs::remove_all(dir);
}

TEST_CASE("io: the generation digest names the complex and every option (Unit 53, FD14)", "[io][generate]") {
    mesh::ExtrudedGeometry x;
    x.axis = 0;
    x.stations = {0.0, 1.0};
    x.cells.push_back({"c", {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}}});
    x.regions.push_back({0});
    x.patches = {{"low", {-1.0, -1.0, -1.0}, {2.0, 2.0, 0.0}}};
    const auto g = mesh::extrude(x);
    REQUIRE(g.has_value());
    mesh::MeshOptions3D o;
    o.max_size = 0.5;
    const std::uint64_t d = mesh::generation_digest(*g, o);
    CHECK(d == mesh::generation_digest(*g, o));
    // Each input changes it, also those that do not change the mesh (max_nodes, an unused region size).
    const auto with = [&](auto change) {
        mesh::PolyhedralGeometry h = *g;
        mesh::MeshOptions3D p = o;
        change(h, p);
        return mesh::generation_digest(h, p);
    };
    CHECK(with([](mesh::PolyhedralGeometry& h, mesh::MeshOptions3D&) { h.vertices[0][2] = 1e-12; }) != d);
    CHECK(with([](mesh::PolyhedralGeometry& h, mesh::MeshOptions3D&) { h.facets[0].patch += "x"; }) != d);
    CHECK(with([](mesh::PolyhedralGeometry& h, mesh::MeshOptions3D&) { h.volumes[0].region = 1; }) != d);
    CHECK(with([](mesh::PolyhedralGeometry&, mesh::MeshOptions3D& p) { p.max_size = 0.25; }) != d);
    CHECK(with([](mesh::PolyhedralGeometry&, mesh::MeshOptions3D& p) { p.region_max_size = {0.5}; }) != d);
    CHECK(with([](mesh::PolyhedralGeometry&, mesh::MeshOptions3D& p) { p.interface_resolution = 0.1; }) != d);
    CHECK(with([](mesh::PolyhedralGeometry&, mesh::MeshOptions3D& p) { p.boxes = {{{0, 0, 0}, {1, 1, 1}, 0.5}}; }) != d);
    CHECK(with([](mesh::PolyhedralGeometry&, mesh::MeshOptions3D& p) { p.radius_edge = 2.5; }) != d);
    CHECK(with([](mesh::PolyhedralGeometry&, mesh::MeshOptions3D& p) { p.max_nodes = 1000; }) != d);
}

TEST_CASE("io: a generated 3D mesh round-trips with its patches, regions, certificate and 3D report",
          "[io][generate]") {
    const fs::path dir = fs::temp_directory_path() / "nitcad_io_generated_3d_test";
    fs::create_directories(dir);
    mesh::Geometry3D g;
    g.domain.vertices = {{0, 0, 0}, {1e-4, 0, 0}, {1e-4, 1e-4, 0}, {0, 1e-4, 0},
                         {0, 0, 1e-4}, {1e-4, 0, 1e-4}, {1e-4, 1e-4, 1e-4}, {0, 1e-4, 1e-4}};
    g.domain.faces = {{{0, 3, 7, 4}, "x_min"}, {{1, 2, 6, 5}, "x_max"}, {{0, 1, 5, 4}, "y_min"},
                      {{3, 2, 6, 7}, "y_max"}, {{0, 1, 2, 3}, "bottom"}, {{4, 5, 6, 7}, "top"}};
    g.domain.faces[5].patches = {{"gate", {{0, 0, 1e-4}, {5e-5, 0, 1e-4}, {5e-5, 1e-4, 1e-4}, {0, 1e-4, 1e-4}}},
                                 {"rest", {{5e-5, 0, 1e-4}, {1e-4, 0, 1e-4}, {1e-4, 1e-4, 1e-4}, {5e-5, 1e-4, 1e-4}}}};
    g.interfaces = {mesh::InterfacePlane{{mesh::Point{0, 0, 6e-5}, mesh::Point{1e-4, 0, 6e-5}, mesh::Point{0, 1e-4, 6e-5}}}};
    g.regions = {{{5e-5, 5e-5, 2e-5}, 0}, {{5e-5, 5e-5, 9e-5}, 1}};
    mesh::MeshOptions3D o;
    o.max_size = 3e-5;
    o.interface_resolution = 3e-6;
    const auto m = mesh::generate_mesh(g, o);
    REQUIRE(m.has_value());
    REQUIRE(m->certificate.has_value());
    CHECK(m->generation_digest == 0);  // FD14: the convex generator's files gain no record
    REQUIRE(io::save_generated_mesh(dir / "g3.ntc", *m).has_value());
    const auto back = io::load_generated_mesh(dir / "g3.ntc");
    REQUIRE(back.has_value());
    CHECK(back->generation_digest == 0);
    CHECK(back->mesh.dimension() == 3);
    CHECK(mesh::graph_digest(back->mesh) == mesh::graph_digest(m->mesh));
    CHECK(back->simplices.cells == m->simplices.cells);
    CHECK(back->simplices.facets == m->simplices.facets);
    CHECK(back->simplices.facet_patch == m->simplices.facet_patch);
    CHECK(back->simplices.patch_names == m->simplices.patch_names);
    CHECK(back->node_region == m->node_region);
    REQUIRE(back->certificate.has_value());
    CHECK(back->certificate->interfaces_digest == m->certificate->interfaces_digest);
    CHECK(back->report.tetrahedra == m->report.tetrahedra);
    CHECK(back->report.slivers == m->report.slivers);
    CHECK(back->report.band_facet_points == m->report.band_facet_points);
    CHECK(back->report.min_dihedral_deg == m->report.min_dihedral_deg);
    CHECK(back->report.max_radius_edge == m->report.max_radius_edge);
    fs::remove_all(dir);
}
