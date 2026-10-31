// Generated meshes on disk (Unit 16; DECISIONS.md M14): simplices, regions, report and P4
// certificate round-trip; a plain mesh, a stale certificate and inconsistent simplices are refused.
#include <catch2/catch_test_macros.hpp>

#include <cstdint>
#include <filesystem>
#include <utility>
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

// Unit 54 (UD11, UC10): the gated fin of the swept-mesh tests (a slab, a fin and an overlayer whose region
// changes under the gate, x in [1, 2]), its gate patch ending at the straddled stations.
namespace {

std::pair<mesh::ExtrudedGeometry, mesh::SweptOptions> swept_gated_fin() {
    const auto P = [](double u, double v) { return mesh::Point{u, v, 0.0}; };
    mesh::ExtrudedGeometry x;
    x.axis = 0;
    x.stations = {0.0, 1.0, 2.0, 3.0};
    x.cells = {{"slab", {P(0, 0), P(1, 0), P(1, 1), P(1, 2), P(1, 3), P(0, 3)}},
               {"fin", {P(1, 1), P(2, 1), P(2, 2), P(1, 2)}},
               {"over", {P(1, 0), P(2.5, 0), P(2.5, 3), P(1, 3), P(1, 2), P(2, 2), P(2, 1), P(1, 1)}}};
    x.regions = {{0, 0, 0}, {1, 1, 1}, {2, 3, 2}};
    x.patches = {{"gate", {1.0, 1.0, 0.0}, {2.0, 2.5, 3.0}}};
    mesh::SweptOptions o;
    o.section.max_size = 0.3;
    o.section.interface_resolution = 0.04;
    o.lines = {0.0, 0.3, 0.6, 0.85, 1.5, 2.15, 2.4, 2.7, 3.0};
    o.straddles = {{1.0, 0.1}, {2.0, 0.1}};
    return {x, o};
}

}  // namespace

TEST_CASE("io: a swept mesh round-trips with its corner patches, P6 report, swept counts and digest (Unit 54, UD11)",
          "[io][generate]") {
    const fs::path dir = fs::temp_directory_path() / "nitcad_io_generated_swept_test";
    fs::create_directories(dir);
    const auto [x, o] = swept_gated_fin();
    auto generated = mesh::generate_swept(x, o);
    if (!generated) FAIL(generated.error().message);
    const mesh::GeneratedMesh& m = *generated;
    REQUIRE_FALSE(m.simplices.facet_corner_patch.empty());
    REQUIRE(m.certificate.has_value());
    REQUIRE(m.certificate->patches.has_value());
    CHECK(m.report.axis_lines == o.lines.size() + 2 * o.straddles.size());
    CHECK(m.report.prisms * 3 == m.report.tetrahedra);

    REQUIRE(io::save_generated_mesh(dir / "s.ntc", m).has_value());
    const auto back = io::load_generated_mesh(dir / "s.ntc");
    if (!back) FAIL(back.error().message);
    CHECK(mesh::graph_digest(back->mesh) == mesh::graph_digest(m.mesh));
    CHECK(back->simplices.cells == m.simplices.cells);
    CHECK(back->simplices.facet_corner_patch == m.simplices.facet_corner_patch);
    CHECK(back->node_region == m.node_region);
    REQUIRE(back->certificate.has_value());
    REQUIRE(back->certificate->patches.has_value());
    CHECK(back->certificate->patches->cross_pieces == m.certificate->patches->cross_pieces);
    CHECK(mesh::certifies(*back->certificate, back->mesh, back->node_region));
    CHECK(back->report.prisms == m.report.prisms);
    CHECK(back->report.axis_lines == m.report.axis_lines);
    CHECK(back->report.layers == m.report.layers);
    CHECK(back->report.corner_blocks == m.report.corner_blocks);
    CHECK(back->report.tetrahedra == m.report.tetrahedra);
    CHECK(back->generation_digest == mesh::generation_digest(x, o));
    // Another generator's file carries no swept counts (the record is absent, its bytes unchanged).
    REQUIRE(io::save_generated_mesh(dir / "t.ntc", mesh::GeneratedMesh{m.mesh, m.simplices, m.node_region, m.certificate,
                                                                         mesh::MeshReport{.nodes = m.report.nodes}, 0})
                .has_value());
    const auto plain = io::load_generated_mesh(dir / "t.ntc");
    if (!plain) FAIL(plain.error().message);
    CHECK(plain->report.prisms == 0);
    CHECK(fs::file_size(dir / "t.ntc") < fs::file_size(dir / "s.ntc"));
    // The four swept counts are distinct fields (this mesh has no layers): distinct values read back in place.
    mesh::GeneratedMesh counted = m;
    counted.report.layers = 3;
    counted.report.corner_blocks = 2;
    REQUIRE(io::save_generated_mesh(dir / "c.ntc", counted).has_value());
    const auto counts = io::load_generated_mesh(dir / "c.ntc");
    if (!counts) FAIL(counts.error().message);
    CHECK(counts->report.prisms == m.report.prisms);
    CHECK(counts->report.axis_lines == m.report.axis_lines);
    CHECK(counts->report.layers == 3);
    CHECK(counts->report.corner_blocks == 2);
    fs::remove_all(dir);
}

TEST_CASE("io: the swept generation digest names the extrusion and every swept option (Unit 54, UD11)", "[io][generate]") {
    const auto [x, o] = swept_gated_fin();
    const std::uint64_t d = mesh::generation_digest(x, o);
    CHECK(d == mesh::generation_digest(x, o));
    const auto with = [&, &x = x, &o = o](auto change) {
        mesh::ExtrudedGeometry h = x;
        mesh::SweptOptions p = o;
        change(h, p);
        return mesh::generation_digest(h, p);
    };
    using G = mesh::ExtrudedGeometry;
    using S = mesh::SweptOptions;
    CHECK(with([](G& h, S&) { h.axis = 1; }) != d);
    CHECK(with([](G& h, S&) { h.stations[1] = 1.01; }) != d);
    CHECK(with([](G& h, S&) { h.cells[1].name = "core"; }) != d);
    CHECK(with([](G& h, S&) { h.cells[1].polygon[0][0] = 1.01; }) != d);
    CHECK(with([](G& h, S&) { h.regions[2][1] = 2; }) != d);
    CHECK(with([](G& h, S&) { h.patches[0].hi[0] = 2.1; }) != d);
    CHECK(with([](G& h, S&) { h.default_patch = "side"; }) != d);
    CHECK(with([](G&, S& p) { p.lines[1] = 0.31; }) != d);
    CHECK(with([](G&, S& p) { p.straddles[0].width = 0.12; }) != d);
    CHECK(with([](G&, S& p) { p.straddles[0].station = 1.1; }) != d);
    CHECK(with([](G&, S& p) { p.max_nodes = 1000; }) != d);
    CHECK(with([](G&, S& p) { p.section.max_size = 0.25; }) != d);
    CHECK(with([](G&, S& p) { p.section.region_max_size = {0.3}; }) != d);
    CHECK(with([](G&, S& p) { p.section.interface_resolution = 0.05; }) != d);
    CHECK(with([](G&, S& p) { p.section.min_angle_deg = 25; }) != d);
    CHECK(with([](G&, S& p) { p.section.max_nodes = 1000; }) != d);
    CHECK(with([](G&, S& p) { p.section.boxes = {{{0, 0, 0}, {1, 1, 0}, 0.1}}; }) != d);
    mesh::BoundaryLayer l;
    l.a = {1, 1, 0};
    l.b = {2, 1, 0};
    l.along = {0, 1};
    l.lines = {-0.1, 0.1};
    const std::uint64_t layered = with([&](G&, S& p) { p.section.layers = {l}; });
    CHECK(layered != d);
    const auto layer_with = [&](auto change) {
        mesh::BoundaryLayer k = l;
        change(k);
        return with([&](G&, S& p) { p.section.layers = {k}; });
    };
    using L = mesh::BoundaryLayer;
    CHECK(layer_with([](L& k) { k.a[0] = 1.01; }) != layered);
    CHECK(layer_with([](L& k) { k.b[1] = 1.01; }) != layered);
    CHECK(layer_with([](L& k) { k.along = {0, 0.5, 1}; }) != layered);
    CHECK(layer_with([](L& k) { k.lines = {-0.1, 0.12}; }) != layered);
    CHECK(layer_with([](L& k) { k.offsets = {0.1}; }) != layered);
    CHECK(layer_with([](L& k) { k.straddle = 0.02; }) != layered);
    // The swept digest is not the complex's (FD14) for the same extrusion.
    const auto complex = mesh::extrude(x);
    REQUIRE(complex.has_value());
    CHECK(mesh::generation_digest(*complex, mesh::MeshOptions3D{}) != d);
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
