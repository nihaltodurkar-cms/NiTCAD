// Gates for the decision layer (Unit 17, slice 5; DECISIONS.md Unit 17 ID1-ID3, AD12, AR5; AC17,
// AC21): the a priori prediction against hand-computed scales, the choice of adaptation points, a
// deliberately wrong prediction guarded and handed to basic, each AD12 trigger forced once, and the
// certificate independent of the policy.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/adapt.hpp"
#include "NiTCAD/solve/adapt_policy.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "legacy_graded_mesh.hpp"

using namespace NiTCAD;

namespace {

constexpr double inf = std::numeric_limits<double>::infinity();

// The legacy diode of adapt_loop_test.cpp: 2 um, abrupt at 1 um, 1e17 each side, legacy silicon.
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

// A builder that keeps the 21-node start unless asked for at least `honour` lines (AD11's M' asks
// 41): the loop's own refinements change nothing, to force AD12's triggers.
solve::DeviceBuilder stuck_builder(std::size_t honour) {
    return [honour](const solve::MeshRefinement& r) -> std::expected<solve::AdaptableMesh, base::Error> {
        const std::vector<double> x = r.axes[0].size() >= honour ? r.axes[0] : uniform(2e-4, 21);
        return solve::AdaptableMesh{.device = legacy_diode(x), .kind = solve::AdaptMeshKind::tensor};
    };
}

solve::BiasOptions legacy_options() {
    solve::BiasOptions o{.models = {.hole_mobility = assemble::HoleMobility::caughey_thomas}};
    o.step_control.max_cuts = 20;
    o.weighted_current = true;
    return o;
}

const std::vector<std::vector<double>> points{{-1.0, 0.0}, {-0.5, 0.0}, {0.0, 0.0}, {0.3, 0.0}, {0.6, 0.0}};

// No prediction anywhere: every target infinite.
solve::SizeTargets no_targets(const device::Device& d) {
    const std::size_t n = d.mesh().node_count();
    return {std::vector<double>(n, inf), std::vector<std::array<double, 3>>(n, std::array<double, 3>{})};
}

// AC17's wrong prediction: 5 nm across x in the anode's bulk (x < 0.3 um) only, nowhere near the
// junction.
solve::SizeTargets substrate_only(const device::Device& d) {
    solve::SizeTargets t = no_targets(d);
    const auto pts = d.mesh().points();
    for (std::size_t i = 0; i < pts.size(); ++i) {
        if (pts[i][0] < 0.3e-4) {
            t.size[i] = 5e-7;
            t.normal[i] = {1.0, 0.0, 0.0};
        }
    }
    return t;
}

deck::Plan plan_of(const std::string& text) {
    auto d = deck::parse(text, "policy.deck", {});
    if (!d) FAIL(d.error().message);
    auto p = deck::build(*d);
    if (!p) FAIL(p.error().message);
    return std::move(*p);
}

// The MOS-C of adapt_indicators_test.cpp: p 1e17 under 1.5 nm of oxide, n+ poly, the gate at 1 V.
const std::string mosc = R"(
[device]
temperature = 300
[mesh]
type = tensor
x = linspace(0, 200nm, 101) + linspace(200nm, 1um, 41)
y = 0, 0.5um, 1um
[region silicon]
material = silicon
[doping]
profile = acceptor constant 1e17
[contact gate]
kind = gate
patch = x_min
oxide_thickness = 1.5nm
electrode = n_poly
[contact substrate]
kind = ohmic
patch = x_max
[sweep]
contact gate = 1.0
[solver]
equations = equilibrium_poisson
)";

}  // namespace

TEST_CASE("policy: the a priori prediction at a junction matches the depletion width and Debye length (ID2)",
          "[adapt_policy]") {
    const device::Device d = legacy_diode(uniform(2e-4, 41));
    auto p = solve::predict_sizes(d, points, legacy_options());
    REQUIRE(p.has_value());
    // One junction edge: 0.95 um (p) to 1.0 um (n), the junction point at 0.975 um.
    CHECK(p->junction_points == 1);
    CHECK(p->interface_points == 0);
    const double T = 300.0, vt = base::thermal_voltage(T), q = base::q_C;
    const double eps = base::eps0_F_per_cm * d.relative_permittivity(0);
    const double ni = physics::intrinsic_density(d.material(0), T);
    const double neff = 0.5e17;
    const double W = std::sqrt(2.0 * eps * (vt * std::log(1e34 / (ni * ni)) + 1.0) / (q * neff));
    const double LD = std::sqrt(eps * vt / (q * neff));
    const auto pts = d.mesh().points();
    std::size_t inside = 0;
    for (std::size_t i = 0; i < pts.size(); ++i) {
        const bool in = std::abs(pts[i][0] - 0.975e-4) <= W;
        if (in) {
            ++inside;
            CHECK(std::abs(p->targets.size[i] - 0.5 * LD) <= 1e-12 * LD);
            CHECK(p->targets.normal[i] == std::array<double, 3>{1.0, 0.0, 0.0});
        } else {
            CHECK(std::isinf(p->targets.size[i]));
        }
    }
    CHECK(p->predicted_nodes == inside);
    CHECK(inside > 2);
    // Every edge within W is longer than L_D / 2: all exceed; edges outside do not.
    for (std::size_t e = 0; e < d.mesh().edges().size(); ++e) {
        const auto& ed = d.mesh().edges()[e];
        const bool near = std::isfinite(p->targets.size[static_cast<std::size_t>(ed.first)]) ||
                          std::isfinite(p->targets.size[static_cast<std::size_t>(ed.second)]);
        CHECK(solve::edge_exceeds(d, p->targets, e) == near);
    }
    // A bias point of the wrong size is refused.
    const std::vector<std::vector<double>> bad{{0.0}};
    CHECK(solve::predict_sizes(d, bad, legacy_options()).error().code == base::ErrorCode::invalid_input);
}

TEST_CASE("policy: the a priori prediction under a gate is the classical inversion length (ID2)", "[adapt_policy]") {
    const deck::Plan p = plan_of(mosc);
    const device::Device& d = p.device;
    auto pr = solve::predict_sizes(d, p.points, p.options);
    REQUIRE(pr.has_value());
    CHECK(pr->junction_points == 0);
    const auto gate = d.find_contact("gate");
    REQUIRE(gate != nullptr);
    CHECK(pr->interface_points == gate->nodes.size());
    const double vt = base::thermal_voltage(300.0);
    const double Es = 3.9 * 1.0 / (d.relative_permittivity(gate->nodes[0]) * 1.5e-7);
    const double lambda = vt / Es;
    for (const mesh::NodeId v : gate->nodes) {
        const auto i = static_cast<std::size_t>(v);
        CHECK(std::abs(pr->targets.size[i] - lambda) <= 1e-12 * lambda);
        CHECK(std::abs(pr->targets.normal[i][0] - 1.0) <= 1e-12);  // into the silicon, along x
    }
    // The depletion region below (AR5): L_D / 2 at 1e17 within W_d = sqrt(2 eps 2 phi_F / (q N)).
    const double eps = base::eps0_F_per_cm * d.relative_permittivity(gate->nodes[0]);
    const double ni = physics::intrinsic_density(d.material(gate->nodes[0]), 300.0);
    const double Wd = std::sqrt(2.0 * eps * 2.0 * vt * std::log(1e17 / ni) / (base::q_C * 1e17));
    const double LD = std::sqrt(eps * vt / (base::q_C * 1e17));
    const auto pts = d.mesh().points();
    for (std::size_t i = 0; i < pts.size(); ++i) {
        if (pts[i][0] == 0.0) continue;
        if (pts[i][0] <= Wd * (1 - 1e-12)) {
            CHECK(std::abs(pr->targets.size[i] - 0.5 * LD) <= 1e-12 * LD);
        } else if (pts[i][0] > Wd * (1 + 1e-12)) {
            CHECK(std::isinf(pr->targets.size[i]));
        }
    }
    // Edges from the gate across the interface exceed the target; edges along it (y) do not.
    std::size_t across = 0, along = 0;
    const auto is_gate = [&](mesh::NodeId v) { return std::ranges::binary_search(gate->nodes, v); };
    for (std::size_t e = 0; e < d.mesh().edges().size(); ++e) {
        const auto& ed = d.mesh().edges()[e];
        const auto a = d.mesh().points()[static_cast<std::size_t>(ed.first)], b = d.mesh().points()[static_cast<std::size_t>(ed.second)];
        if (!is_gate(ed.first) && !is_gate(ed.second)) continue;
        const bool exceeds = solve::edge_exceeds(d, pr->targets, e);
        if (a[1] == b[1]) {
            CHECK(exceeds);
            ++across;
        } else {
            CHECK_FALSE(exceeds);
            ++along;
        }
    }
    CHECK(across > 0);
    CHECK(along > 0);
}

TEST_CASE("policy: the adaptation points are the last and the largest changes (ID3)", "[adapt_policy]") {
    const auto field = [](std::vector<double> e) {
        solve::IndicatorField f;
        f.dominant.assign(e.size(), solve::Indicator::potential);
        f.excess = std::move(e);
        return f;
    };
    const std::vector<solve::IndicatorField> f{field({0, 0, 0}), field({1, 0, 0}), field({1, 0, 0}), field({5, 0, 0}),
                                               field({5, 0, 2})};
    // Changes from the previous point: 1 (point 1), 0 (2), 4 (3); the last (4) always.
    CHECK(*solve::choose_points(f, 1) == std::vector<std::size_t>{4});
    CHECK(*solve::choose_points(f, 2) == std::vector<std::size_t>{3, 4});
    CHECK(*solve::choose_points(f, 3) == std::vector<std::size_t>{1, 3, 4});
    CHECK(*solve::choose_points(f, 9) == std::vector<std::size_t>{1, 2, 3, 4});
    // Only eligible points (a goal above its resolution): the last eligible is 3, then point 1.
    const std::vector<char> eligible{1, 1, 0, 1, 0};
    CHECK(*solve::choose_points(f, 2, eligible) == std::vector<std::size_t>{1, 3});
    CHECK(*solve::choose_points(f, 9, eligible) == std::vector<std::size_t>{1, 3});
    CHECK(!solve::choose_points(f, 2, std::vector<char>(5, 0)).has_value());
    CHECK(!solve::choose_points(f, 2, std::vector<char>(4, 1)).has_value());
    CHECK(!solve::choose_points(f, 0).has_value());
    CHECK(!solve::choose_points({}, 2).has_value());
    const std::vector<solve::IndicatorField> uneven{field({0}), field({0, 1})};
    CHECK(!solve::choose_points(uneven, 2).has_value());
}

TEST_CASE("policy: the certificate does not depend on the policy (AC21)", "[adapt_policy]") {
    // The same final mesh forced: the intelligent run's predictions are empty (every proposal a miss,
    // so AD3's marks throughout).
    solve::AdaptOptions o;
    o.goals = {{0}};
    o.tolerance = 1e-2;
    o.at = {4};
    o.max_passes = 12;
    auto b = solve::adapt(diode_builder(21), points, legacy_options(), o);
    o.policy = solve::Policy::intelligent;
    o.intelligent.prediction = no_targets;
    auto i = solve::adapt(diode_builder(21), points, legacy_options(), o);
    REQUIRE(b.has_value());
    REQUIRE(i.has_value());
    REQUIRE(b->passes.size() == i->passes.size());
    const auto pb = b->mesh->device.mesh().points(), pi = i->mesh->device.mesh().points();
    REQUIRE(pb.size() == pi.size());
    for (std::size_t k = 0; k < pb.size(); ++k) REQUIRE(pb[k] == pi[k]);
    REQUIRE(b->certificates.size() == i->certificates.size());
    const auto same = [](double x, double y) { return std::bit_cast<std::uint64_t>(x) == std::bit_cast<std::uint64_t>(y); };
    for (std::size_t c = 0; c < b->certificates.size(); ++c) {
        CHECK(same(b->certificates[c].value, i->certificates[c].value));
        CHECK(same(b->certificates[c].estimate, i->certificates[c].estimate));
        CHECK(same(b->certificates[c].verified_error, i->certificates[c].verified_error));
        CHECK(same(b->certificates[c].order, i->certificates[c].order));
        CHECK(same(b->certificates[c].allowed, i->certificates[c].allowed));
        CHECK(b->certificates[c].status == i->certificates[c].status);
    }
    // The decision layer did act: misses and a fallback recorded on the intelligent run only.
    CHECK(b->misses == 0);
    CHECK(b->fallback == solve::FallbackTrigger::none);
    CHECK(i->misses >= 2);
    CHECK(i->fallback == solve::FallbackTrigger::two_misses);
}

TEST_CASE("policy: a wrong prediction is guarded, handed to basic and still certified (AC17)", "[adapt_policy]") {
    // AC6's diode, its 2x reference on the hand-made mesh; adjoint weighting (the estimate the guard
    // compares), two adaptation points chosen from the coarse sweep (ID3).
    const std::vector<double> hand = legacy_graded_mesh(2e-4, 1e-4, 1e-8, 1e-6);
    std::vector<double> ref_axis;
    for (std::size_t k = 0; k + 1 < hand.size(); ++k) {
        ref_axis.push_back(hand[k]);
        ref_axis.push_back(0.5 * (hand[k] + hand[k + 1]));
    }
    ref_axis.push_back(hand.back());
    auto ref = solve::sweep_bias(legacy_diode(ref_axis), points, legacy_options());
    REQUIRE(ref.has_value());
    REQUIRE(!ref->stopped);
    std::vector<double> ref4;
    for (std::size_t k = 0; k + 1 < ref_axis.size(); ++k) {
        ref4.push_back(ref_axis[k]);
        ref4.push_back(0.5 * (ref_axis[k] + ref_axis[k + 1]));
    }
    ref4.push_back(ref_axis.back());
    auto r4 = solve::sweep_bias(legacy_diode(ref4), points, legacy_options());
    REQUIRE(r4.has_value());
    for (std::size_t k = 0; k < points.size(); ++k) {
        WARN("AC17 reference point " << k << ": 2x " << ref->points[k].weighted_current[0] << ", 4x "
                                     << r4->points[k].weighted_current[0] << ", relative "
                                     << std::abs(r4->points[k].weighted_current[0] - ref->points[k].weighted_current[0]) /
                                            std::abs(r4->points[k].weighted_current[0]));
    }

    solve::AdaptOptions o;
    o.goals = {{0}};
    o.tolerance = 1e-2;  // the guard's mechanics; at 1e-3 adjoint weighting needs over 20 passes here
    o.weighting = solve::Weighting::adjoint;
    o.max_passes = 20;  // the wrong field is applied over two proposals (the growth cap), then basic
    o.policy = solve::Policy::intelligent;
    o.intelligent.max_points = 2;
    o.intelligent.prediction = substrate_only;
    const auto t_i = std::chrono::steady_clock::now();
    auto r = solve::adapt(diode_builder(21), points, legacy_options(), o);
    const double seconds_i = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_i).count();
    REQUIRE(r.has_value());
    solve::AdaptOptions ob = o;
    ob.policy = solve::Policy::basic;
    ob.at = r->points;
    const auto t_b = std::chrono::steady_clock::now();
    auto b = solve::adapt(diode_builder(21), points, legacy_options(), ob);
    const double seconds_b = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_b).count();
    REQUIRE(b.has_value());
    for (std::size_t k = 0; k < r->passes.size(); ++k) {
        const auto& p = r->passes[k];
        WARN("AC17 pass " << k << ": " << solve::to_string(p.source) << ", nodes " << p.nodes << ", relative estimate "
                          << p.relative_estimate << (p.miss ? " (miss)" : "") << (p.empty_proposal ? " (empty proposal)" : ""));
    }
    WARN("AC17: points " << r->points.size() << " (last " << r->points.back() << "), coarse " << r->coarse_nodes
                         << " nodes, prediction rounds " << r->prediction_rounds << ", misses " << r->misses
                         << ", fallback " << solve::to_string(r->fallback) << " at pass " << r->fallback_pass << "; stop "
                         << solve::to_string(r->stop) << ", nodes " << r->mesh->device.mesh().node_count() << " (basic "
                         << b->mesh->device.mesh().node_count() << "), time " << seconds_i << " s (basic " << seconds_b
                         << " s)");
    CHECK(r->stop == solve::AdaptStop::converged);
    CHECK(r->misses >= 2);
    CHECK(r->fallback == solve::FallbackTrigger::two_misses);
    CHECK(r->fallback_pass >= 0);
    // AC17's time bound (total <= 2x basic) is measured on the 2D MOSFET (adapt_intelligent_test.cpp): here
    // the times are milliseconds, and the wrong field's nodes rival the whole adapted mesh (AR5).
    for (const auto& c : r->certificates) {
        // Against the 4x reference: the 2x moves 4.7e-4 from it at -0.5 V, half the tolerance (AR5).
        const double I_ref = r4->points[c.point].weighted_current[0];
        WARN("AC17 point " << c.point << ": " << solve::to_string(c.status) << ", value " << c.value << ", reference "
                           << I_ref << ", relative " << std::abs(c.value - I_ref) / std::abs(I_ref));
        CHECK(c.status == solve::CertificateStatus::certified);
        CHECK(std::abs(c.value - I_ref) <= c.allowed);
    }
}

TEST_CASE("policy: each AD12 trigger switches to basic at the pass it names (AC17)", "[adapt_policy]") {
    solve::AdaptOptions o;
    o.goals = {{0}};
    o.at = {4};
    o.policy = solve::Policy::intelligent;
    SECTION("b: a pass lowers neither the estimate nor the goal change") {
        // Only M' (41 lines, the adjoint estimate's mesh) is built: every pass solves the same mesh.
        // The proposals are not empty (two intervals split, which the builder ignores) and, at margin 0,
        // an equal estimate is no miss.
        o.tolerance = 1e-3;
        o.weighting = solve::Weighting::adjoint;
        o.max_passes = 4;
        o.intelligent.prediction = [](const device::Device& d) {
            solve::SizeTargets t = no_targets(d);
            const auto pts = d.mesh().points();
            for (std::size_t i = 0; i < pts.size(); ++i) {
                if (pts[i][0] < 0.2e-4) {
                    t.size[i] = 0.06e-4;
                    t.normal[i] = {1.0, 0.0, 0.0};
                }
            }
            return t;
        };
        o.intelligent.miss_margin = 0.0;
        auto r = solve::adapt(stuck_builder(41), points, legacy_options(), o);
        REQUIRE(r.has_value());
        CHECK(r->fallback == solve::FallbackTrigger::no_progress);
        CHECK(r->fallback_pass == 2);
        CHECK(r->stop == solve::AdaptStop::max_passes);
    }
    SECTION("c: AD11 rejects a proposed stop") {
        // The same mesh twice is a zero change, a proposed stop; M' (41 lines) is built and rejects it.
        o.tolerance = 1e-3;
        o.theta = 1.0;
        o.max_passes = 3;
        o.intelligent.prediction = no_targets;
        auto r = solve::adapt(stuck_builder(41), points, legacy_options(), o);
        REQUIRE(r.has_value());
        REQUIRE(r->passes.size() >= 2);
        CHECK(r->passes[1].rejected);
        CHECK(r->fallback == solve::FallbackTrigger::rejected_stop);
        CHECK(r->fallback_pass == 1);
    }
    SECTION("a: two consecutive proposals miss") {
        o.tolerance = 1e-2;
        o.max_passes = 6;
        o.intelligent.prediction = no_targets;
        auto r = solve::adapt(diode_builder(21), points, legacy_options(), o);
        REQUIRE(r.has_value());
        CHECK(r->passes[0].empty_proposal);
        CHECK(r->passes[1].empty_proposal);
        CHECK(r->fallback == solve::FallbackTrigger::two_misses);
        CHECK(r->fallback_pass == 1);
    }
    SECTION("basic never falls back") {
        o.policy = solve::Policy::basic;
        o.tolerance = 1e-3;
        o.max_passes = 3;
        auto r = solve::adapt(stuck_builder(41), points, legacy_options(), o);
        REQUIRE(r.has_value());
        CHECK(r->fallback == solve::FallbackTrigger::none);
        CHECK(r->misses == 0);
    }
}

TEST_CASE("policy: the cost model fits its passes (ID5, AR7)", "[adapt_policy]") {
    // Passes timed by a known model: 10 Newton iterations of 1e-9 nnz^1.3 s, 2e-6 s per node to mesh,
    // 1e-6 to transfer.
    const auto pass = [](double nnz, bool transfer) {
        const double nodes = nnz / 30.0;
        return solve::CostSample{.newton_iterations = 10.0,
                                 .matrix_nonzeros = nnz,
                                 .fill = 4.0 * nnz,
                                 .dofs = 3.0 * nodes,
                                 .nodes = nodes,
                                 .seconds_solve = 10.0 * 1e-9 * std::pow(nnz, 1.3),
                                 .seconds_mesh = 2e-6 * nodes,
                                 .seconds_transfer = transfer ? 1e-6 * nodes : 0.0};
    };
    std::vector<solve::CostSample> s;
    for (int k = 0; k < 5; ++k) s.push_back(pass(1e4 * std::exp2(k), k > 0));
    SECTION("no sample: no prediction") {
        const solve::CostModel m = solve::fit_cost_model({});
        CHECK(m.samples == 0);
        CHECK(std::isnan(m.pass_seconds(10, 1e4, 4e4, 1e3, 333)));
        CHECK(std::isnan(m.pass_seconds_nodes(333)));
    }
    SECTION("one sample: a alone at b = 1.5, exact on it") {
        const solve::CostModel m = solve::fit_cost_model(std::span(s).first(1));
        CHECK(m.samples == 1);
        CHECK(m.b == 1.5);
        CHECK(m.c == 0.0);
        CHECK(m.d == 0.0);
        CHECK(std::abs(m.per_iteration(1e4, 4e4, 1e3) / (s[0].seconds_solve / 10.0) - 1.0) < 1e-12);
        CHECK(m.transfer_per_node == 0.0);  // no pass transferred
        CHECK(m.beta == 1.0);
    }
    SECTION("two samples: b through both") {
        const solve::CostModel m = solve::fit_cost_model(std::span(s).first(2));
        CHECK(std::abs(m.b - 1.3) < 1e-12);
        CHECK(std::abs(m.transfer_per_node - 1e-6) < 1e-15);
        CHECK(std::abs(m.mesh_per_node - 2e-6) < 1e-15);
    }
    SECTION("five samples: every pass that transferred within 1%, and the next predicted within 5%") {
        const solve::CostModel m = solve::fit_cost_model(s);
        CHECK(m.samples == 5);
        CHECK(m.b >= 1.0);
        CHECK(m.b <= 2.0);
        CHECK(m.residual < 0.01);
        for (const solve::CostSample& x : std::span(s).subspan(1)) {  // a predicted pass transfers; the first did not
            const double t = x.seconds_solve + x.seconds_mesh + x.seconds_transfer;
            CHECK(std::abs(m.pass_seconds(x.newton_iterations, x.matrix_nonzeros, x.fill, x.dofs, x.nodes) / t - 1.0) < 0.01);
        }
        const solve::CostSample next = pass(1e4 * 32.0, true);
        const double t = next.seconds_solve + next.seconds_mesh + next.seconds_transfer;
        CHECK(std::abs(m.pass_seconds(10, next.matrix_nonzeros, next.fill, next.dofs, next.nodes) / t - 1.0) < 0.05);
    }
    SECTION("the node-only model is a power law in the nodes") {
        std::vector<solve::CostSample> p;
        for (int k = 0; k < 4; ++k) {
            const double nodes = 1e3 * std::exp2(k);
            p.push_back({.newton_iterations = 5.0, .matrix_nonzeros = 30.0 * nodes, .dofs = 3.0 * nodes, .nodes = nodes,
                         .seconds_solve = 3e-6 * std::pow(nodes, 1.2)});
        }
        const solve::CostModel m = solve::fit_cost_model(p);
        CHECK(std::abs(m.beta - 1.2) < 1e-9);
        CHECK(std::abs(m.pass_seconds_nodes(16e3) / (3e-6 * std::pow(16e3, 1.2)) - 1.0) < 1e-9);
        CHECK(m.residual_nodes < 1e-9);
    }
    SECTION("samples without a solve are ignored") {
        std::vector<solve::CostSample> p{solve::CostSample{}, s[0]};
        CHECK(solve::fit_cost_model(p).samples == 1);
    }
}

TEST_CASE("policy: the local orders fit the realized ratios (ID5, AR7)", "[adapt_policy]") {
    using solve::Indicator;
    const std::vector<solve::OrderSample> s{
        {Indicator::electron_density, 1, 0.25}, {Indicator::potential, 0, 0.5}, {Indicator::electron_density, 1, 0.25},
        {Indicator::potential, 0, 0.5},         {Indicator::potential, 0, 0.5}, {Indicator::electron_density, 1, 0.25},
        {Indicator::hole_density, 0, 1e-6},     {Indicator::hole_density, 0, 0.0}, {Indicator::doping, 2, inf}};
    const std::vector<solve::OrderFit> f = solve::fit_orders(s);
    REQUIRE(f.size() == 3);  // doping's only sample is not finite: no fit
    CHECK(f[0].kind == Indicator::potential);
    CHECK(f[0].samples == 3);
    CHECK(f[0].order == 1.0);  // (3 + 1) / 4: the ridge toward 1 changes nothing here
    CHECK(f[1].kind == Indicator::electron_density);
    CHECK(f[1].region == 1);
    CHECK(f[1].order == 1.75);  // (3 x 2 + 1) / 4
    CHECK(f[2].kind == Indicator::hole_density);
    CHECK(f[2].samples == 1);  // the zero ratio skipped
    CHECK(f[2].order == 3.0);  // clamped
    CHECK(solve::order_of(f, Indicator::electron_density, 1) == 1.75);
    CHECK(solve::order_of(f, Indicator::electron_density, 0) == 1.0);
    CHECK(solve::order_of({}, Indicator::potential, 0) == 1.0);
}

TEST_CASE("policy: the allocation is deterministic and its report recomputes from the passes (AC18)", "[adapt_policy]") {
    solve::AdaptOptions o;
    o.goals = {{0}};
    o.tolerance = 1e-3;
    o.at = {0, 4};
    o.max_passes = 14;
    o.weighting = solve::Weighting::adjoint;
    o.policy = solve::Policy::intelligent;
    auto a = solve::adapt(diode_builder(21), points, legacy_options(), o);
    auto b = solve::adapt(diode_builder(21), points, legacy_options(), o);
    REQUIRE(a.has_value());
    REQUIRE(b.has_value());
    // Bit-identical, timings aside (and the timing-fitted model's seconds, which only reports and
    // budgets read; AR7).
    REQUIRE(a->passes.size() == b->passes.size());
    std::size_t allocated = 0, measured = 0;
    for (std::size_t k = 0; k < a->passes.size(); ++k) {
        const auto &p = a->passes[k], &q = b->passes[k];
        CHECK(p.nodes == q.nodes);
        CHECK(p.goal == q.goal);
        CHECK(p.estimate == q.estimate);
        CHECK(p.relative_estimate == q.relative_estimate);
        CHECK(p.source == q.source);
        CHECK(p.dofs == q.dofs);
        CHECK(p.matrix_nonzeros == q.matrix_nonzeros);
        CHECK(p.factor_fill == q.factor_fill);
        CHECK(p.candidate_count == q.candidate_count);
        CHECK(p.candidates_taken == q.candidates_taken);
        CHECK(std::bit_cast<std::uint64_t>(p.predicted_relative) == std::bit_cast<std::uint64_t>(q.predicted_relative));
        CHECK(p.predicted_dofs == q.predicted_dofs);
        REQUIRE(p.orders.size() == q.orders.size());
        for (std::size_t j = 0; j < p.orders.size(); ++j) {
            CHECK(p.orders[j].order == q.orders[j].order);
            CHECK(p.orders[j].samples == q.orders[j].samples);
        }
        REQUIRE(p.candidates.size() == q.candidates.size());
        for (std::size_t j = 0; j < p.candidates.size(); ++j) {
            const auto &x = p.candidates[j], &y = q.candidates[j];
            CHECK(x.axis == y.axis);
            CHECK(x.index == y.index);
            CHECK(x.benefit == y.benefit);
            CHECK(x.work == y.work);
            CHECK(x.order == y.order);
            CHECK(std::bit_cast<std::uint64_t>(x.realized) == std::bit_cast<std::uint64_t>(y.realized));
        }
        // Recomputed from the saved passes: the orders used are the pass's fits; the predicted estimate
        // is the estimate less the taken benefits (all listed when at most 20 were taken); a
        // candidate's benefit over l levels is its share times 1 - 2^-lp, and its realized reduction
        // its share times one less its ratio.
        if (!p.candidates.empty()) ++allocated;
        double predicted = p.relative_estimate;
        for (const auto& x : p.candidates) {
            CHECK(x.order == solve::order_of(p.orders, x.indicator, x.region));
            CHECK(x.levels >= 1);
            CHECK(x.levels <= 3);
            predicted -= x.benefit;
            if (std::isfinite(x.realized)) {
                ++measured;
                const double share = x.benefit / (1.0 - std::exp2(-x.levels * x.order));
                CHECK(std::abs(x.realized - share * (1.0 - x.realized_ratio)) <= 1e-12 * std::abs(share));
            }
        }
        if (!p.candidates.empty() && p.candidates_taken == p.candidates.size()) {
            CHECK(std::abs(predicted - p.predicted_relative) <= 1e-12 * p.relative_estimate);
        }
        // The cost model: fitted on the passes up to this one; this pass predicted from the earlier.
        std::vector<solve::CostSample> s;
        for (std::size_t j = 0; j <= k; ++j) {
            const auto& r = a->passes[j];
            s.push_back({.newton_iterations = static_cast<double>(r.newton_iterations),
                         .matrix_nonzeros = static_cast<double>(r.matrix_nonzeros),
                         .fill = static_cast<double>(r.factor_fill),
                         .dofs = static_cast<double>(r.dofs),
                         .nodes = static_cast<double>(r.nodes),
                         .seconds_solve = r.seconds_solve,
                         .seconds_mesh = r.seconds_mesh,
                         .seconds_transfer = r.seconds_transfer});
        }
        const solve::CostModel all = solve::fit_cost_model(s);
        CHECK(all.a == p.cost_model.a);
        CHECK(all.b == p.cost_model.b);
        CHECK(all.alpha == p.cost_model.alpha);
        if (k > 0) {
            const solve::CostModel before = solve::fit_cost_model(std::span(s).first(k));
            CHECK(before.pass_seconds(s[k - 1].newton_iterations, s[k].matrix_nonzeros, s[k].fill, s[k].dofs, s[k].nodes) ==
                  p.predicted_seconds);
            CHECK(before.pass_seconds_nodes(s[k].nodes) == p.predicted_seconds_nodes);
        } else {
            CHECK(std::isnan(p.predicted_seconds));
        }
        CHECK(p.matrix_nonzeros > 0);
        CHECK(p.factor_fill >= p.matrix_nonzeros / 2);  // SparseLU reports its factors
    }
    CHECK(allocated >= 1);
    CHECK(measured >= 1);
    REQUIRE(a->certificates.size() == b->certificates.size());
    for (std::size_t j = 0; j < a->certificates.size(); ++j) {
        CHECK(a->certificates[j].status == b->certificates[j].status);
        CHECK(a->certificates[j].estimate == b->certificates[j].estimate);
        CHECK(std::bit_cast<std::uint64_t>(a->certificates[j].verified_error) ==
              std::bit_cast<std::uint64_t>(b->certificates[j].verified_error));
    }
    WARN("AC18 diode: stop " << solve::to_string(a->stop) << " after " << a->passes.size() << " passes, "
                             << a->mesh->device.mesh().node_count() << " nodes, certificate "
                             << solve::to_string(a->certificates[0].status));
    const auto pa = a->mesh->device.mesh().points(), pb = b->mesh->device.mesh().points();
    REQUIRE(pa.size() == pb.size());
    for (std::size_t i = 0; i < pa.size(); ++i) CHECK(pa[i] == pb[i]);
}
