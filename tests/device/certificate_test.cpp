// P4 certification of multi-region devices (Unit 16; DECISIONS.md U2): a device with more than one
// region on a mesh that is not a tensor grid needs a certificate from mesh::certify_interfaces,
// issued for its mesh and node regions.
#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/interfaces.hpp"
#include "NiTCAD/mesh/simplex_mesh.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD::device;
using NiTCAD::base::ErrorCode;
using NiTCAD::mesh::Mesh;
using NiTCAD::mesh::NodeId;
using NiTCAD::mesh::Point;
using NiTCAD::physics::silicon;

namespace {

// A unit square in triangles on a lattice with rows at 0.5 -+ s/2: region 0 below y = 0.5 and 1
// above, the bottom and top rows as patches. The lattice makes the interface exact.
struct TwoRegion {
    NiTCAD::mesh::SimplexMesh simplices;
    Mesh dual;
    std::vector<RegionId> region;
};

TwoRegion two_region_square() {
    const double s = 1.0 / 32.0;
    const std::vector<double> xs{0.0, 0.25, 0.5, 0.75, 1.0};
    const std::vector<double> ys{0.0, 0.25, 0.5 - 0.5 * s, 0.5 + 0.5 * s, 0.75, 1.0};
    NiTCAD::mesh::SimplexMesh m;
    m.dimension = 2;
    for (const double y : ys) {
        for (const double x : xs) m.points.push_back(Point{x, y, 0});
    }
    const auto v = [&](std::size_t i, std::size_t j) {
        return static_cast<NodeId>(i + xs.size() * j);
    };
    for (std::size_t j = 0; j + 1 < ys.size(); ++j) {
        for (std::size_t i = 0; i + 1 < xs.size(); ++i) {
            m.cells.insert(m.cells.end(), {v(i, j), v(i + 1, j), v(i + 1, j + 1), v(i, j),
                                           v(i + 1, j + 1), v(i, j + 1)});
        }
    }
    m.patch_names = {"bottom", "top"};
    const std::size_t last = ys.size() - 1;
    for (std::size_t i = 0; i + 1 < xs.size(); ++i) {
        m.facets.insert(m.facets.end(), {v(i, 0), v(i + 1, 0), v(i, last), v(i + 1, last)});
        m.facet_patch.insert(m.facet_patch.end(), {0, 1});
    }
    auto dual = NiTCAD::mesh::make_dual_mesh(m);
    REQUIRE(dual.has_value());
    std::vector<RegionId> region;
    for (const Point& p : m.points) region.push_back(p[1] < 0.5 ? 0 : 1);
    return {std::move(m), std::move(*dual), std::move(region)};
}

DeviceDescription two_region_device(const TwoRegion& c) {
    const std::size_t n = c.dual.node_count();
    std::vector<NodeId> bottom, top;
    for (std::size_t i = 0; i < n; ++i) {
        if (c.dual.points()[i][1] == 0.0) bottom.push_back(static_cast<NodeId>(i));
        if (c.dual.points()[i][1] == 1.0) top.push_back(static_cast<NodeId>(i));
    }
    return DeviceDescription{.mesh = c.dual,
                             .temperature_K = 300.0,
                             .regions = {{"lower", silicon()}, {"upper", silicon()}},
                             .node_region = c.region,
                             .donors = std::vector<double>(n, 1e16),
                             .acceptors = std::vector<double>(n, 0.0),
                             .contacts = {{"bottom", ContactKind::ohmic, bottom},
                                          {"top", ContactKind::ohmic, top}}};
}

const std::vector<NiTCAD::mesh::InterfacePiece> mid_line{{{{0, 0.5, 0}, {1, 0.5, 0}}, 0, 1}};

}  // namespace

TEST_CASE("device: a multi-region device on a general mesh needs a valid P4 certificate") {
    const TwoRegion c = two_region_square();
    REQUIRE(c.dual.origin() == NiTCAD::mesh::MeshOrigin::graph);
    const auto cert = NiTCAD::mesh::certify_interfaces(
        c.simplices, c.dual, c.region, NiTCAD::mesh::InterfaceSet{mid_line, 1.0 / 32.0});
    REQUIRE(cert.has_value());
    SECTION("missing: refused") {
        const auto r = Device::create(two_region_device(c));
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().code == ErrorCode::invalid_input);
        CHECK(r.error().message.find("needs a P4 certificate") != std::string::npos);
    }
    SECTION("valid: accepted") {
        auto d = two_region_device(c);
        d.interface_certificate = *cert;
        CHECK(Device::create(std::move(d)).has_value());
    }
    SECTION("stale, for other node regions: refused") {
        auto d = two_region_device(c);
        d.interface_certificate = *cert;
        d.node_region[6] = 1;  // a node of the second row moved to the upper region
        const auto r = Device::create(std::move(d));
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().message.find("stale") != std::string::npos);
    }
    SECTION("stale, for another mesh: refused") {
        auto d = two_region_device(c);
        d.interface_certificate = *cert;
        d.interface_certificate->mesh_digest ^= 1;  // as if issued for a different graph
        const auto r = Device::create(std::move(d));
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().message.find("stale") != std::string::npos);
    }
    SECTION("a certificate that failed P4 does not exist") {
        // With the node regions split elsewhere than the declared interface, certification
        // itself fails, so there is no certificate to give the device.
        std::vector<RegionId> off(c.region);
        for (std::size_t i = 0; i < off.size(); ++i) off[i] = c.simplices.points[i][1] < 0.3 ? 0 : 1;
        CHECK_FALSE(NiTCAD::mesh::certify_interfaces(c.simplices, c.dual, off,
                                                     NiTCAD::mesh::InterfaceSet{mid_line, 1.0 / 32.0})
                        .has_value());
    }
    SECTION("one region needs none") {
        auto d = two_region_device(c);
        d.regions = {{"silicon", silicon()}};
        d.node_region.assign(d.node_region.size(), 0);
        CHECK(Device::create(std::move(d)).has_value());
    }
}

TEST_CASE("device: tensor-grid meshes are exempt from P4 certification by provenance") {
    const std::vector<double> x{0.0, 0.5, 1.0};
    const std::vector<double> y{0.0, 0.4, 0.6, 1.0};
    auto grid = NiTCAD::mesh::make_tensor_grid(x, y);
    REQUIRE(grid.has_value());
    CHECK(grid->origin() == NiTCAD::mesh::MeshOrigin::tensor_grid);
    std::vector<RegionId> region;
    for (const auto& p : grid->points()) region.push_back(p[1] < 0.5 ? 0 : 1);
    const std::size_t n = grid->node_count();
    DeviceDescription d{.mesh = *grid,
                        .temperature_K = 300.0,
                        .regions = {{"lower", silicon()}, {"upper", silicon()}},
                        .node_region = region,
                        .donors = std::vector<double>(n, 1e16),
                        .acceptors = std::vector<double>(n, 0.0),
                        .contacts = {{"bottom", ContactKind::ohmic, {0, 1, 2}},
                                     {"top", ContactKind::ohmic, {9, 10, 11}}}};
    CHECK(Device::create(std::move(d)).has_value());
}

// ---------------------------------------------------------------------------------------------
// 1D (DECISIONS.md U16): only a verified tensor grid is exempt; any other 1D mesh with two regions
// needs certify_interfaces_1d.

namespace {

// Nodes at x = 0, 1, 2, 3, 4 (um) joined in order, regions 0 0 1 1 1: the transition at 1.5 um.
Mesh chain_1d(NiTCAD::mesh::MeshOrigin origin) {
    std::vector<Point> p;
    std::vector<NiTCAD::mesh::Edge> e;
    for (int i = 0; i < 5; ++i) p.push_back({1e-4 * i, 0, 0});
    for (int i = 0; i < 4; ++i) e.push_back({i, i + 1, 1e-4, 1.0});
    return *Mesh::from_parts(1, p, std::vector<double>(5, 1e-4), e,
                             {{"x_min", {0}, {1.0}}, {"x_max", {4}, {1.0}}}, origin);
}

DeviceDescription two_region_1d(Mesh m) {
    return DeviceDescription{.mesh = std::move(m),
                             .temperature_K = 300.0,
                             .regions = {{"left", silicon()}, {"right", silicon()}},
                             .node_region = {0, 0, 1, 1, 1},
                             .donors = std::vector<double>(5, 1e16),
                             .acceptors = std::vector<double>(5, 0.0),
                             .contacts = {{"a", ContactKind::ohmic, {0}}, {"b", ContactKind::ohmic, {4}}}};
}

}  // namespace

TEST_CASE("device: a general 1D mesh with two regions needs a 1D P4 certificate", "[device]") {
    SECTION("missing: refused") {
        const auto r = Device::create(two_region_1d(chain_1d(NiTCAD::mesh::MeshOrigin::graph)));
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().message.find("needs a P4 certificate") != std::string::npos);
    }
    SECTION("the transition on a declared point: certified and accepted") {
        auto d = two_region_1d(chain_1d(NiTCAD::mesh::MeshOrigin::graph));
        const std::vector<NiTCAD::mesh::InterfacePoint1D> at{{1.5e-4, 0, 1}};
        const auto cert = NiTCAD::mesh::certify_interfaces_1d(d.mesh, d.node_region, at);
        REQUIRE(cert.has_value());
        CHECK(cert->report.cross_faces == 1);
        d.interface_certificate = *cert;
        CHECK(Device::create(std::move(d)).has_value());
    }
    SECTION("the declared point off the dual face, or undeclared: not certified") {
        const auto d = two_region_1d(chain_1d(NiTCAD::mesh::MeshOrigin::graph));
        const std::vector<NiTCAD::mesh::InterfacePoint1D> off{{1.4e-4, 0, 1}};
        const auto a = NiTCAD::mesh::certify_interfaces_1d(d.mesh, d.node_region, off);
        REQUIRE_FALSE(a.has_value());
        CHECK(a.error().code == ErrorCode::degenerate_mesh);
        const auto b = NiTCAD::mesh::certify_interfaces_1d(d.mesh, d.node_region, {});
        REQUIRE_FALSE(b.has_value());
        CHECK(b.error().message.find("no declared interface") != std::string::npos);
    }
    SECTION("stale for other regions: refused") {
        auto d = two_region_1d(chain_1d(NiTCAD::mesh::MeshOrigin::graph));
        const std::vector<NiTCAD::mesh::InterfacePoint1D> at{{1.5e-4, 0, 1}};
        d.interface_certificate = *NiTCAD::mesh::certify_interfaces_1d(d.mesh, d.node_region, at);
        d.node_region = {0, 0, 0, 1, 1};
        const auto r = Device::create(std::move(d));
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().message.find("stale") != std::string::npos);
    }
    SECTION("a verified 1D tensor grid is exempt") {
        CHECK(Device::create(two_region_1d(chain_1d(NiTCAD::mesh::MeshOrigin::tensor_grid))).has_value());
        const auto grid = NiTCAD::mesh::make_tensor_grid(std::vector<double>{0, 1e-4, 2e-4, 3e-4, 4e-4});
        REQUIRE(grid.has_value());
        CHECK(Device::create(two_region_1d(*grid)).has_value());
    }
}

TEST_CASE("mesh: the 1D tensor-grid tag needs the chain structure (U16)", "[device]") {
    const std::vector<Point> p{{0, 0, 0}, {1, 0, 0}, {2, 0, 0}};
    const auto tagged = [&](std::vector<Point> pts, std::vector<NiTCAD::mesh::Edge> e) {
        return Mesh::from_parts(1, std::move(pts), std::vector<double>(3, 1.0), std::move(e), {},
                                NiTCAD::mesh::MeshOrigin::tensor_grid);
    };
    CHECK(tagged(p, {{0, 1, 1, 1}, {1, 2, 1, 1}}).has_value());
    CHECK_FALSE(tagged(p, {{0, 2, 2, 1}, {1, 2, 1, 1}}).has_value());  // an edge skips node 1
    CHECK_FALSE(tagged(p, {{0, 1, 1, 1}}).has_value());                // not every neighbour joined
    CHECK_FALSE(tagged({{0, 0, 0}, {2, 0, 0}, {1, 0, 0}}, {{0, 1, 2, 1}, {1, 2, 1, 1}}).has_value());  // 0-2 skips x = 1
    CHECK_FALSE(tagged({{0, 0, 0}, {1, 0, 0}, {1, 0, 0}}, {{0, 1, 1, 1}, {1, 2, 1, 1}}).has_value());  // two nodes at x = 1
    // Untagged, the same graphs are general meshes, which need a certificate for two regions.
    CHECK(Mesh::from_parts(1, p, std::vector<double>(3, 1.0), {{0, 2, 2, 1}, {1, 2, 1, 1}}, {}).has_value());
}

TEST_CASE("device: a generated 3D mesh with a planar interface carries its certificate (6d)", "[device]") {
    using namespace NiTCAD::mesh;
    Geometry3D g;
    g.domain.vertices = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    g.domain.faces = {{{0, 3, 7, 4}, "x_min"}, {{1, 2, 6, 5}, "x_max"}, {{0, 1, 5, 4}, "y_min"},
                      {{3, 2, 6, 7}, "y_max"}, {{0, 1, 2, 3}, "bottom"}, {{4, 5, 6, 7}, "top"}};
    g.interfaces = {InterfacePlane{{Point{0, 0, 0.5}, Point{1, 0, 0.5}, Point{0, 1, 0.5}}}};
    g.regions = {{{0.5, 0.5, 0.25}, 0}, {{0.5, 0.5, 0.75}, 1}};
    MeshOptions3D o;
    o.max_size = 0.3;
    o.interface_resolution = 0.03;
    const auto m = generate_mesh(g, o);
    REQUIRE(m.has_value());
    REQUIRE(m->certificate.has_value());
    const std::size_t n = m->mesh.node_count();
    const auto contact = [&](const char* name) {
        const BoundaryPatch* p = m->mesh.find_boundary(name);
        REQUIRE(p != nullptr);
        return std::vector<NodeId>(p->nodes.begin(), p->nodes.end());
    };
    const auto description = [&]() {
        return DeviceDescription{.mesh = m->mesh,
                                 .temperature_K = 300.0,
                                 .regions = {{"lower", silicon()}, {"upper", silicon()}},
                                 .node_region = std::vector<RegionId>(m->node_region.begin(), m->node_region.end()),
                                 .donors = std::vector<double>(n, 1e16),
                                 .acceptors = std::vector<double>(n, 0.0),
                                 .contacts = {{"bottom", ContactKind::ohmic, contact("bottom")},
                                              {"top", ContactKind::ohmic, contact("top")}}};
    };
    const auto missing = Device::create(description());
    REQUIRE_FALSE(missing.has_value());
    CHECK(missing.error().message.find("needs a P4 certificate") != std::string::npos);
    auto d = description();
    d.interface_certificate = m->certificate;
    CHECK(Device::create(std::move(d)).has_value());
    auto stale = description();
    stale.interface_certificate = m->certificate;
    stale.node_region[0] = 1 - stale.node_region[0];
    CHECK_FALSE(Device::create(std::move(stale)).has_value());
}
