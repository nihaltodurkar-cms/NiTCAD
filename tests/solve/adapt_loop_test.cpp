// Gates for the adaptive loop (Unit 17, slice 3; DECISIONS.md Unit 17 AD1, AD4-AD6, AD9, AD11, AD13;
// AC6, AC10 part, AC20 part): the legacy diode from a coarse mesh reaches the reference currents
// with fewer nodes than its hand-made mesh, stopping by the goal rule and certified; a budget stop
// is never certified; two runs are bit-identical; refusals.
#include <catch2/catch_test_macros.hpp>

#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/adapt.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "legacy_graded_mesh.hpp"

using namespace NiTCAD;

namespace {

// The legacy diode (bias_test.cpp's): 2 um, abrupt at 1 um, 1e17 each side, legacy silicon.
device::Device legacy_diode(const std::vector<double>& x) {
    auto m = mesh::make_tensor_grid(x);
    REQUIRE(m.has_value());
    auto cells = mesh::TensorCells::create(*m, x);
    REQUIRE(cells.has_value());
    const std::size_t n = m->node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) (m->points()[i][0] < 1e-4 ? acceptors[i] : donors[i]) = 1e17;
    auto anode = m->find_boundary("x_min")->nodes;
    auto cathode = m->find_boundary("x_max")->nodes;
    auto d = device::Device::create({.mesh = std::move(*m),
                                     .temperature_K = 300.0,
                                     .regions = {{"silicon", physics::legacy_silicon()}},
                                     .node_region = std::vector<device::RegionId>(n, 0),
                                     .donors = std::move(donors),
                                     .acceptors = std::move(acceptors),
                                     .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                                                  {"cathode", device::ContactKind::ohmic, std::move(cathode)}},
                                     .cells = std::move(*cells)});
    REQUIRE(d.has_value());
    return std::move(*d);
}

std::vector<double> uniform(double L, int nodes) {
    std::vector<double> a(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) a[static_cast<std::size_t>(i)] = i == nodes - 1 ? L : L * i / (nodes - 1);
    return a;
}

solve::DeviceBuilder diode_builder(int initial_nodes) {
    return [initial_nodes](const solve::MeshRefinement& r) -> std::expected<solve::AdaptableMesh, base::Error> {
        const std::vector<double> x = r.axes[0].empty() ? uniform(2e-4, initial_nodes) : r.axes[0];
        return solve::AdaptableMesh{.device = legacy_diode(x), .kind = solve::AdaptMeshKind::tensor};
    };
}

// The legacy options: Caughey-Thomas holes (bias_test.cpp), step cutting for the reverse point, and
// the weighted current (the goal the loop reads when it is on).
solve::BiasOptions legacy_options() {
    solve::BiasOptions o{.models = {.hole_mobility = assemble::HoleMobility::caughey_thomas}};
    o.step_control.max_cuts = 20;
    o.weighted_current = true;  // Unit 32: the reverse current is below the terminal current's resolution
    return o;
}

// Forward and reverse: the sweep runs from equilibrium down to -1 V and up to 0.6 V.
const std::vector<std::vector<double>> points{{-1.0, 0.0}, {-0.5, 0.0}, {0.0, 0.0}, {0.3, 0.0}, {0.6, 0.0}};

double anode(const results::Sweep& s, std::size_t k) { return s.points[k].weighted_current[0]; }

}  // namespace

TEST_CASE("adapt: the legacy diode reaches its 2x reference with fewer nodes than the hand-made mesh (AC6)",
          "[adapt_loop]") {
    // The hand-made mesh (the Unit 9 gate's) and the reference on it refined twice (every interval split).
    const std::vector<double> hand = legacy_graded_mesh(2e-4, 1e-4, 1e-8, 1e-6);
    std::vector<double> ref_axis;
    for (std::size_t i = 0; i + 1 < hand.size(); ++i) {
        ref_axis.push_back(hand[i]);
        ref_axis.push_back(0.5 * (hand[i] + hand[i + 1]));
    }
    ref_axis.push_back(hand.back());
    const device::Device ref_device = legacy_diode(ref_axis);
    auto ref = solve::sweep_bias(ref_device, points, legacy_options());
    REQUIRE(ref.has_value());
    REQUIRE(!ref->stopped);

    solve::AdaptOptions o;
    o.goals = {{0}};
    o.tolerance = 1e-3;
    o.at = {0, 4};  // -1 V and 0.6 V
    o.max_passes = 12;
    auto r = solve::adapt(diode_builder(21), points, legacy_options(), o);
    REQUIRE(r.has_value());
    for (std::size_t k = 0; k < r->passes.size(); ++k) {
        const auto& p = r->passes[k];
        WARN("AC6 pass " << k << ": nodes " << p.nodes << ", I(-1 V) " << p.goal[0] << ", I(0.6 V) " << p.goal[1]
                         << ", marked " << p.marking.count << ", Newton " << p.newton_iterations << ", transfer "
                         << p.seconds_transfer << " s, solve " << p.seconds_solve << " s");
    }
    WARN("AC6 stop " << solve::to_string(r->stop) << "; verification nodes " << r->verification_nodes << " in "
                     << r->verification_seconds << " s; hand mesh " << hand.size() << " nodes");
    CHECK(r->stop == solve::AdaptStop::converged);
    const std::size_t nodes = r->mesh->device.mesh().node_count();
    CHECK(nodes < hand.size());
    for (const std::size_t k : o.at) {
        const double I = anode(r->sweep, k), I_ref = anode(*ref, k);
        WARN("AC6 point " << k << ": adapted " << I << ", reference " << I_ref << ", relative "
                          << std::abs(I - I_ref) / std::abs(I_ref));
        CHECK(std::abs(I - I_ref) <= 1e-3 * std::abs(I_ref));
    }
    REQUIRE(r->certificates.size() == 2);
    for (const auto& c : r->certificates) {
        WARN("AC6 certificate point " << c.point << ": " << solve::to_string(c.status) << ", estimate " << c.estimate
                                      << ", verified " << c.verified_error << " (p " << c.order << "), allowed "
                                      << c.allowed << ", actual " << std::abs(c.value - anode(*ref, c.point)));
        // AC20: no certified result has an actual error above its tolerance.
        if (c.status == solve::CertificateStatus::certified) {
            CHECK(std::abs(c.value - anode(*ref, c.point)) <= c.allowed);
        }
    }
}

TEST_CASE("adapt: a budget stop is unverified, never certified (AC20 part)", "[adapt_loop]") {
    solve::AdaptOptions o;
    o.goals = {{0}};
    o.tolerance = 1e-6;
    o.at = {4};
    o.max_nodes = 40;
    auto r = solve::adapt(diode_builder(21), points, legacy_options(), o);
    REQUIRE(r.has_value());
    CHECK(r->stop == solve::AdaptStop::max_nodes);
    CHECK(r->mesh->device.mesh().node_count() <= 40);
    CHECK(r->verification_nodes == 0);
    REQUIRE(r->certificates.size() == 1);
    CHECK(r->certificates[0].status == solve::CertificateStatus::unverified);
    CHECK(std::isnan(r->certificates[0].verified_error));

    o.max_nodes = 0;
    o.max_passes = 2;
    auto p = solve::adapt(diode_builder(21), points, legacy_options(), o);
    REQUIRE(p.has_value());
    CHECK(p->stop == solve::AdaptStop::max_passes);
    CHECK(p->passes.size() == 2);
    CHECK(p->certificates[0].status == solve::CertificateStatus::unverified);
    // verify = off after a goal stop: unverified too.
    o.max_passes = 12;
    o.tolerance = 1e-2;
    o.verify = false;
    auto q = solve::adapt(diode_builder(21), points, legacy_options(), o);
    REQUIRE(q.has_value());
    CHECK(q->stop == solve::AdaptStop::converged);
    CHECK(q->certificates[0].status == solve::CertificateStatus::unverified);
}

TEST_CASE("adapt: two runs are bit-identical (AD9)", "[adapt_loop]") {
    solve::AdaptOptions o;
    o.goals = {{0}};
    o.tolerance = 1e-2;
    o.at = {4};
    auto a = solve::adapt(diode_builder(21), points, legacy_options(), o);
    auto b = solve::adapt(diode_builder(21), points, legacy_options(), o);
    REQUIRE(a.has_value());
    REQUIRE(b.has_value());
    REQUIRE(a->passes.size() == b->passes.size());
    for (std::size_t k = 0; k < a->passes.size(); ++k) {
        CHECK(a->passes[k].nodes == b->passes[k].nodes);
        for (std::size_t j = 0; j < a->passes[k].goal.size(); ++j) {
            CHECK(std::bit_cast<std::uint64_t>(a->passes[k].goal[j]) == std::bit_cast<std::uint64_t>(b->passes[k].goal[j]));
        }
    }
    const auto pa = a->mesh->device.mesh().points(), pb = b->mesh->device.mesh().points();
    REQUIRE(pa.size() == pb.size());
    for (std::size_t i = 0; i < pa.size(); ++i) CHECK(pa[i] == pb[i]);
    REQUIRE(a->certificates.size() == b->certificates.size());
    CHECK(std::bit_cast<std::uint64_t>(a->certificates[0].verified_error) ==
          std::bit_cast<std::uint64_t>(b->certificates[0].verified_error));
    CHECK(a->certificates[0].status == b->certificates[0].status);
}

TEST_CASE("adapt: refusals (AC10 part)", "[adapt_loop]") {
    const auto build = diode_builder(11);
    solve::AdaptOptions o;
    o.goals = {{0}};
    CHECK(!solve::adapt(build, {}, legacy_options(), o).has_value());
    solve::AdaptOptions none = o;
    none.goals.clear();
    CHECK(!solve::adapt(build, points, legacy_options(), none).has_value());
    solve::AdaptOptions bad = o;
    bad.goals = {{5}};
    CHECK(solve::adapt(build, points, legacy_options(), bad).error().code == base::ErrorCode::invalid_input);
    bad = o;
    bad.tolerance = 0.0;
    CHECK(!solve::adapt(build, points, legacy_options(), bad).has_value());
    bad = o;
    bad.at = {9};
    CHECK(!solve::adapt(build, points, legacy_options(), bad).has_value());
    bad = o;
    bad.theta = 0.0;
    CHECK(!solve::adapt(build, points, legacy_options(), bad).has_value());
    bad = o;
    bad.max_passes = 0;
    CHECK(!solve::adapt(build, points, legacy_options(), bad).has_value());
    bad = o;
    bad.weighting = solve::Weighting::adjoint;
    CHECK(!solve::adapt(build, points, legacy_options(), bad).has_value());
}
