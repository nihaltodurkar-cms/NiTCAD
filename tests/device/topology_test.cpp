// The floating-region check (ARCHITECTURE.md 6.10, the topology gate assigned to Unit 6): every
// connected part of the mesh graph must contain an ohmic contact node, or its potential or carrier
// densities are undetermined. A gate (Unit 12) does not anchor a part.
#include <catch2/catch_test_macros.hpp>

#include <cstddef>
#include <expected>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD::device;
using NiTCAD::base::ErrorCode;
using NiTCAD::mesh::BoundaryPatch;
using NiTCAD::mesh::Edge;
using NiTCAD::mesh::Mesh;
using NiTCAD::mesh::NodeId;
using NiTCAD::mesh::Point;

namespace {

// A 1D mesh of unit-length edges between the given node pairs (x coordinates chosen so each edge
// has length 1), with every node on a boundary patch of its own so any node can be a contact.
Mesh graph(std::vector<double> x, const std::vector<std::pair<NodeId, NodeId>>& pairs) {
    std::vector<Point> points;
    std::vector<BoundaryPatch> boundary;
    for (std::size_t i = 0; i < x.size(); ++i) {
        points.push_back({x[i], 0.0, 0.0});
        boundary.push_back({"n" + std::to_string(i), {static_cast<NodeId>(i)}, {1.0}});
    }
    std::vector<Edge> edges;
    for (const auto& [a, b] : pairs) edges.push_back({a, b, 1.0, 1.0});
    return *Mesh::from_parts(1, std::move(points), std::vector<double>(x.size(), 1.0),
                             std::move(edges), std::move(boundary));
}

// Ohmic contacts on the given node sets, and a gate on each of the gate nodes (on its own patch).
std::expected<Device, NiTCAD::base::Error> build(Mesh m,
                                                 std::vector<std::vector<NodeId>> contacts,
                                                 const std::vector<NodeId>& gates = {}) {
    const std::size_t n = m.node_count();
    std::vector<Contact> list;
    for (std::size_t c = 0; c < contacts.size(); ++c) {
        list.push_back({"c" + std::to_string(c), ContactKind::ohmic, std::move(contacts[c])});
    }
    for (const NodeId g : gates) {
        Contact gate{"g" + std::to_string(g), ContactKind::gate, {g}};
        gate.gate = {.boundary = "n" + std::to_string(g), .oxide_thickness_cm = 1e-6};
        list.push_back(std::move(gate));
    }
    return Device::create({.mesh = std::move(m),
                           .temperature_K = 300.0,
                           .regions = {{"silicon", NiTCAD::physics::silicon()}},
                           .node_region = std::vector<RegionId>(n, 0),
                           .donors = std::vector<double>(n, 1e16),
                           .acceptors = std::vector<double>(n, 0.0),
                           .contacts = std::move(list)});
}

// Two interleaved chains: 0-2-4 (x = 0, 1, 2) and 1-3-5 (x = 10, 11, 12), not joined.
Mesh two_chains() {
    return graph({0.0, 10.0, 1.0, 11.0, 2.0, 12.0}, {{0, 2}, {2, 4}, {1, 3}, {3, 5}});
}

}  // namespace

TEST_CASE("topology: one contact anchors a connected device") {
    // A single contact is enough: the Poisson and continuity rows of the connected graph then
    // each have a Dirichlet row.
    const auto x = std::vector<double>{0.0, 1e-5, 2e-5, 3e-5};
    REQUIRE(build(*NiTCAD::mesh::make_tensor_grid(x), {{0}}).has_value());
    REQUIRE(build(*NiTCAD::mesh::make_tensor_grid(x, x), {{0}}).has_value());
}

TEST_CASE("topology: a device with no contact is floating") {
    const auto e = build(*NiTCAD::mesh::make_tensor_grid(std::vector<double>{0.0, 1.0}), {});
    REQUIRE_FALSE(e.has_value());
    REQUIRE(e.error().code == ErrorCode::invalid_input);
    REQUIRE(e.error().message.find("floating") != std::string::npos);
    REQUIRE(e.error().context->index == 0);
}

TEST_CASE("topology: each disconnected part needs its own contact") {
    // Contact on the first chain only: the second chain floats; its lowest node is 1.
    const auto first_only = build(two_chains(), {{4}});
    REQUIRE_FALSE(first_only.has_value());
    REQUIRE(first_only.error().context->index == 1);
    // Contact on the second chain only: the first floats; its lowest node is 0.
    const auto second_only = build(two_chains(), {{5}});
    REQUIRE_FALSE(second_only.has_value());
    REQUIRE(second_only.error().context->index == 0);
    // One contact on each chain anchors both.
    REQUIRE(build(two_chains(), {{0}, {3}}).has_value());
    // So does one contact holding a node of each chain.
    REQUIRE(build(two_chains(), {{2, 5}}).has_value());
}

TEST_CASE("topology: an isolated node is a part of its own") {
    // Nodes 0-1-2 joined, node 3 alone.
    const auto m = [] { return graph({0.0, 1.0, 2.0, 50.0}, {{0, 1}, {1, 2}}); };
    const auto alone = build(m(), {{0}});
    REQUIRE_FALSE(alone.has_value());
    REQUIRE(alone.error().context->index == 3);
    REQUIRE(build(m(), {{0}, {3}}).has_value());
}

TEST_CASE("topology: a gate does not anchor a part") {
    // A gate fixes no carrier density: chain 1-3-5 with only a gate on node 3 floats.
    const auto gated = build(two_chains(), {{4}}, {3});
    REQUIRE_FALSE(gated.has_value());
    REQUIRE(gated.error().message.find("ohmic") != std::string::npos);
    REQUIRE(gated.error().context->index == 1);
    // A gate alongside an ohmic contact on each chain is fine.
    REQUIRE(build(two_chains(), {{4}, {5}}, {0, 3}).has_value());
}
