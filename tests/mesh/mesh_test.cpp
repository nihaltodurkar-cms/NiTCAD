// Validation of the generic mesh graph (Mesh::from_parts), the entry point for every producer.
#include <catch2/catch_test_macros.hpp>

#include <limits>
#include <vector>

#include "NiTCAD/mesh/mesh.hpp"

using namespace NiTCAD::mesh;
using NiTCAD::base::ErrorCode;

namespace {

// A valid two-cell 1D mesh on [0, 2]: nodes 0, 1, 2 with volumes 0.5, 1, 0.5.
struct Parts {
    int dimension = 1;
    std::vector<Point> points{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}};
    std::vector<double> volumes{0.5, 1.0, 0.5};
    std::vector<Edge> edges{{0, 1, 1.0, 1.0}, {1, 2, 1.0, 1.0}};
    std::vector<BoundaryPatch> boundary{{"left", {0}, {1.0}}, {"right", {2}, {1.0}}};
};

auto build(Parts p) {
    return Mesh::from_parts(p.dimension, std::move(p.points), std::move(p.volumes),
                            std::move(p.edges), std::move(p.boundary));
}

}  // namespace

TEST_CASE("mesh: a valid graph is accepted and read back unchanged") {
    const auto m = build({});
    REQUIRE(m.has_value());
    REQUIRE(m->dimension() == 1);
    REQUIRE(m->node_count() == 3);
    REQUIRE(m->volumes()[1] == 1.0);
    REQUIRE(m->edges()[1].second == 2);
    REQUIRE(m->find_boundary("right")->nodes[0] == 2);
    REQUIRE(m->find_boundary("top") == nullptr);
}

TEST_CASE("mesh: dimension and sizes are checked") {
    for (const int d : {0, 4}) {
        Parts p;
        p.dimension = d;
        REQUIRE(build(p).error().code == ErrorCode::invalid_input);
    }
    Parts sizes;
    sizes.volumes.pop_back();
    REQUIRE(build(sizes).error().code == ErrorCode::invalid_input);
    Parts empty;
    empty.points.clear();
    empty.volumes.clear();
    empty.edges.clear();
    empty.boundary.clear();
    REQUIRE(build(empty).error().code == ErrorCode::invalid_input);
}

TEST_CASE("mesh: coordinates must be finite and zero beyond the dimension") {
    Parts beyond;
    beyond.points[1][1] = 0.5;  // y in a 1D mesh
    const auto b = build(beyond);
    REQUIRE(b.error().code == ErrorCode::invalid_input);
    REQUIRE(b.error().context->index == 1);
    Parts nan;
    nan.points[2][0] = std::numeric_limits<double>::quiet_NaN();
    REQUIRE(build(nan).error().context->index == 2);
}

TEST_CASE("mesh: non-positive control volumes are degenerate_mesh") {
    for (const double v : {0.0, -0.5, std::numeric_limits<double>::infinity()}) {
        Parts p;
        p.volumes[2] = v;
        const auto r = build(p);
        REQUIRE(r.error().code == ErrorCode::degenerate_mesh);
        REQUIRE(r.error().context->index == 2);
    }
}

TEST_CASE("mesh: edge endpoints must be in range, ordered and unique") {
    Parts range;
    range.edges[1].second = 3;
    REQUIRE(build(range).error().code == ErrorCode::invalid_input);
    Parts order;
    order.edges[0] = {1, 0, 1.0, 1.0};
    REQUIRE(build(order).error().code == ErrorCode::invalid_input);
    Parts self;
    self.edges[0] = {1, 1, 1.0, 1.0};
    REQUIRE(build(self).error().code == ErrorCode::invalid_input);
    Parts twice;
    twice.edges.push_back({0, 1, 1.0, 1.0});
    const auto t = build(twice);
    REQUIRE(t.error().code == ErrorCode::invalid_input);
    REQUIRE(t.error().context->index == 2);  // the second listing
}

TEST_CASE("mesh: edge length and coupling area are validated") {
    Parts zero_area;
    zero_area.edges[0].coupling_area = 0.0;
    REQUIRE(build(zero_area).error().code == ErrorCode::degenerate_mesh);
    // A negative dual area, the legacy unstructured hazard (clockwise triangles), is rejected.
    Parts negative_area;
    negative_area.edges[1].coupling_area = -0.25;
    const auto n = build(negative_area);
    REQUIRE(n.error().code == ErrorCode::degenerate_mesh);
    REQUIRE(n.error().context->index == 1);
    Parts wrong_length;
    wrong_length.edges[0].length = 1.001;
    REQUIRE(build(wrong_length).error().code == ErrorCode::degenerate_mesh);
}

TEST_CASE("mesh: boundary patches are validated") {
    Parts repeated;
    repeated.boundary[1].name = "left";
    REQUIRE(build(repeated).error().code == ErrorCode::invalid_input);
    Parts unnamed;
    unnamed.boundary[0].name.clear();
    REQUIRE(build(unnamed).error().code == ErrorCode::invalid_input);
    Parts unsorted;
    unsorted.boundary[0] = {"left", {2, 0}, {1.0, 1.0}};
    REQUIRE(build(unsorted).error().code == ErrorCode::invalid_input);
    Parts mismatch;
    mismatch.boundary[0].areas.push_back(1.0);
    REQUIRE(build(mismatch).error().code == ErrorCode::invalid_input);
    Parts bad_area;
    bad_area.boundary[1].areas[0] = 0.0;
    REQUIRE(build(bad_area).error().code == ErrorCode::degenerate_mesh);
}
