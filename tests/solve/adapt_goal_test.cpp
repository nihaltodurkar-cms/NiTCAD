// AC7 (Unit 17, slice 4; DECISIONS.md Unit 17 AD3, AD6, AR3, AR4): goal weighting earns its place. A
// MOSFET in a p-well whose deep junction with the n-substrate lies off the drain current's path:
// adjoint-weighted against feature-only adaptation to the same tolerance on I_D, against uniform
// refinements extrapolated; the weighted mesh has fewer nodes, its marks avoid the well junction, and its
// adjoint estimate is within a factor 3 of the actual error at the final pass.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <format>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/solve/adapt.hpp"
#include "NiTCAD/solve/bias.hpp"

using namespace NiTCAD;

namespace {

constexpr double width_cm = 1e-4;   // x: 1 um
constexpr double depth_cm = 1.5e-4; // y: 1.5 um, the surface at y = 0
constexpr double junction_cm = 0.6e-4;

std::string axis_text(const std::vector<double>& a) {
    std::string s;
    for (std::size_t i = 0; i < a.size(); ++i) s += (i ? ", " : "") + std::format("{:.17g}cm", a[i]);
    return s;
}

std::vector<double> linspace(double L, int n) {
    std::vector<double> a(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) a[static_cast<std::size_t>(i)] = i == n - 1 ? L : L * i / (n - 1);
    return a;
}

// The depth axis: 2 nm lines over the first 20 nm under the oxide (the inversion layer), 20 nm to
// 0.2 um, 0.1 um below.
std::vector<double> surface_axis() {
    std::vector<double> y;
    for (int i = 0; i < 10; ++i) y.push_back(2e-7 * i);
    for (int i = 0; i < 9; ++i) y.push_back(2e-6 + 2e-6 * i);
    for (int i = 0; i <= 13; ++i) y.push_back(i == 13 ? depth_cm : 2e-5 + 1e-5 * i);
    return y;
}

// Every interval split into k equal parts.
std::vector<double> split(const std::vector<double>& a, int k) {
    std::vector<double> b;
    for (std::size_t i = 0; i + 1 < a.size(); ++i) {
        for (int j = 0; j < k; ++j) b.push_back(a[i] + (a[i + 1] - a[i]) * j / k);
    }
    b.push_back(a.back());
    return b;
}

std::string deck_text(const std::vector<double>& x, const std::vector<double>& y) {
    return "[device]\ntemperature = 300\n[mesh]\ntype = tensor\nx = " + axis_text(x) + "\ny = " + axis_text(y) + R"(
[region silicon]
material = silicon
[doping]
profile = donor constant 1e16
profile = acceptor erfc 1e17 axis=y center=0.6um width=0.02um
profile = donor erfc 1e20 axis=y center=0.1um width=0.01um x=-1um:0.25um
profile = donor erfc 1e20 axis=y center=0.1um width=0.01um x=0.75um:2um
[contact source]
kind = ohmic
patch = y_min
box = 0 0.15um 0 0
[contact drain]
kind = ohmic
patch = y_min
box = 0.85um 1um 0 0
[contact gate]
kind = gate
patch = y_min
box = 0.2um 0.8um 0 0
oxide_thickness = 5nm
electrode = n_poly
[contact body]
kind = ohmic
patch = x_min
box = 0 0 0.3um 0.5um
[contact substrate]
kind = ohmic
patch = y_max
[sweep]
contact gate = 1.5
contact drain = 0.5
[solver]
max_cuts = 20
weighted_current = on
)";
}

deck::Plan plan_of(const std::vector<double>& x, const std::vector<double>& y) {
    auto d = deck::parse(deck_text(x, y), "well.deck", {});
    if (!d) FAIL(d.error().message);
    auto p = deck::build(*d);
    if (!p) FAIL(p.error().message);
    return std::move(*p);
}

solve::DeviceBuilder builder() {
    return [](const solve::MeshRefinement& r) -> std::expected<solve::AdaptableMesh, base::Error> {
        const std::vector<double> x = r.axes[0].empty() ? linspace(width_cm, 21) : r.axes[0];
        const std::vector<double> y = r.axes[1].empty() ? surface_axis() : r.axes[1];
        deck::Plan p = plan_of(x, y);
        return solve::AdaptableMesh{.device = std::move(p.device), .kind = solve::AdaptMeshKind::tensor};
    };
}

std::size_t contact_index(const device::Device& d, const std::string& name) {
    std::size_t c = 0;
    while (c < d.contacts().size() && d.contacts()[c].name != name) ++c;
    INFO("contact " << name);
    REQUIRE(c < d.contacts().size());
    return c;
}

// The y-lines of a tensor axis within 0.15 um of the well junction.
std::size_t lines_near_junction(std::span<const double> y) {
    std::size_t n = 0;
    for (const double v : y) n += std::abs(v - junction_cm) <= 0.15e-4 ? 1 : 0;
    return n;
}

}  // namespace

TEST_CASE("adapt: goal weighting avoids the well junction off the current path (AC7)", "[.adapt_goal]") {
    const deck::Plan coarse = plan_of(linspace(width_cm, 21), surface_axis());
    const std::size_t drain = contact_index(coarse.device, "drain");
    const std::size_t k = coarse.points.size() - 1;

    // The reference (DECISIONS.md Unit 17 AR4), measured once (Release, 2026-11-12): the start refined
    // 2, 4, 8 and 16 times (every interval split; 2665, 10449, 41377 and 164673 nodes) gave I_D
    // 4.77792, 5.03309, 5.12877 and 5.17144; Richardson from 4x, 8x and 16x at their observed order
    // (1.165) gives 5.20577. From 2x, 4x and 8x (order 1.415) it was 5.18618: the order still drifts,
    // so the reference is uncertain by about 0.4%. The 16x solve takes 12 minutes, so the value is
    // recorded; the 2x solve is repeated as the guard that the deck and the physics are unchanged.
    constexpr double I_ref = 5.20577;
    const deck::Plan p2 = plan_of(split(linspace(width_cm, 21), 2), split(surface_axis(), 2));
    auto s2 = solve::sweep_bias(p2.device, p2.points, p2.options);
    REQUIRE(s2.has_value());
    REQUIRE(!s2->stopped);
    const double I2 = s2->points[k].weighted_current[drain];
    WARN("AC7 reference guard: 2x " << I2 << " on " << p2.device.mesh().node_count() << " nodes (recorded 4.77792)");
    REQUIRE(std::abs(I2 - 4.77792) <= 1e-5 * 4.77792);

    solve::AdaptOptions o;
    o.goals = {{drain}};
    o.tolerance = 1e-2;
    o.at = {k};
    o.max_passes = 10;
    struct Run {
        std::size_t nodes = 0;
        double error = 0.0;
        std::size_t junction_marks = 0, marks = 0;
    };
    const auto run = [&](solve::Weighting w, const char* name) {
        solve::AdaptOptions oo = o;
        oo.weighting = w;
        auto r = solve::adapt(builder(), coarse.points, coarse.options, oo);
        REQUIRE(r.has_value());
        Run out;
        for (std::size_t p = 0; p < r->passes.size(); ++p) {
            const auto& rep = r->passes[p];
            WARN("AC7 " << name << " pass " << p << ": nodes " << rep.nodes << ", I_D " << rep.goal[0] << ", estimate "
                        << (rep.estimate.empty() ? 0.0 : rep.estimate[0]) << " on " << rep.estimate_nodes
                        << " nodes in " << rep.seconds_estimate << " s, marked " << rep.marking.count
                        << (rep.rejected ? " (stop rejected)" : "") << ", solve " << rep.seconds_solve << " s");
            out.marks += rep.marking.count;
        }
        out.nodes = r->mesh->device.mesh().node_count();
        out.error = std::abs(r->sweep.points[k].weighted_current[drain] - I_ref);
        const auto& c = r->certificates[0];
        WARN("AC7 " << name << ": stop " << solve::to_string(r->stop) << ", nodes " << out.nodes << ", error "
                    << out.error / std::abs(I_ref) << " relative; certificate " << solve::to_string(c.status)
                    << ", estimate " << c.estimate << ", verified " << c.verified_error << ", actual " << out.error
                    << ", estimate / actual " << c.estimate / out.error);
        // AC20: a certified result within its tolerance of AR4's reference.
        if (c.status == solve::CertificateStatus::certified) CHECK(out.error <= c.allowed);
        return std::make_pair(out, std::move(*r));
    };
    auto [feature, fr] = run(solve::Weighting::none, "feature-only");
    auto [weighted, wr] = run(solve::Weighting::adjoint, "adjoint");
    // Both reach the tolerance.
    CHECK(feature.error <= o.tolerance * std::abs(I_ref));
    CHECK(weighted.error <= o.tolerance * std::abs(I_ref));
    // The weighted mesh has at least 25% fewer nodes.
    WARN("AC7 nodes: adjoint " << weighted.nodes << ", feature-only " << feature.nodes << ", ratio "
                               << static_cast<double>(weighted.nodes) / static_cast<double>(feature.nodes));
    CHECK(static_cast<double>(weighted.nodes) <= 0.75 * static_cast<double>(feature.nodes));
    // Its marks avoid the well junction: the final mesh's y-lines near it, against feature-only's.
    const std::size_t lw = lines_near_junction(wr.mesh->device.cells()->axis(1));
    const std::size_t lf = lines_near_junction(fr.mesh->device.cells()->axis(1));
    const std::vector<double> y0 = surface_axis();
    WARN("AC7 y-lines within 0.15 um of the well junction: adjoint " << lw << ", feature-only " << lf << " (initial "
                                                                  << lines_near_junction(y0) << ")");
    CHECK(lw < lf);
    // The adjoint estimate within a factor 3 of the actual error at the final pass.
    const double est = wr.certificates[0].estimate;
    CHECK(est <= 3.0 * weighted.error);
    CHECK(est >= weighted.error / 3.0);
}

// Fails as frozen (AR7, owner ruling pending), measured 2026-11-13: tolerance mode 12744 DOFs against
// basic's 6642 (1.92x; 3 passes against 6, relative estimate 0.013 against 0.52 at the stop: the
// growth-capped first allocation overshoots the accuracy); budget mode passes (3255 DOFs, error 6.5%
// against basic's 7.2% on 3240).
TEST_CASE("adapt: the allocation needs fewer DOFs than basic on the well MOSFET, and does better on a budget (AC16)",
          "[.adapt_goal][!shouldfail]") {
    // Unit 17 slice 5b (ID6, AR7): adjoint-weighted, policy basic against intelligent, to AC7's
    // tolerance on I_D against AC7's recorded reference (AR4).
    constexpr double I_ref = 5.20577;
    const deck::Plan coarse = plan_of(linspace(width_cm, 21), surface_axis());
    const std::size_t drain = contact_index(coarse.device, "drain");
    const std::size_t k = coarse.points.size() - 1;
    solve::AdaptOptions o;
    o.goals = {{drain}};
    o.tolerance = 1e-2;
    o.at = {k};
    o.max_passes = 10;
    o.weighting = solve::Weighting::adjoint;
    const auto dofs = [](const solve::AdaptResult& r) { return 3 * r.mesh->device.mesh().node_count(); };  // all silicon
    const auto error = [&](const solve::AdaptResult& r) {
        return std::abs(r.sweep.points[k].weighted_current[drain] - I_ref) / std::abs(I_ref);
    };
    const auto run = [&](const solve::AdaptOptions& oo, const char* name) {
        auto r = solve::adapt(builder(), coarse.points, coarse.options, oo);
        if (!r) FAIL(r.error().message);
        for (std::size_t p = 0; p < r->passes.size(); ++p) {
            const auto& rep = r->passes[p];
            WARN("AC16 well " << name << " pass " << p << ": " << solve::to_string(rep.source) << ", " << rep.dofs
                              << " DOFs, I_D " << rep.goal[0] << ", relative estimate " << rep.relative_estimate
                              << ", predicted " << rep.predicted_relative << " on " << rep.predicted_dofs << " DOFs, "
                              << rep.candidates_taken << " of " << rep.candidate_count << " candidates"
                              << (rep.miss ? " (miss)" : ""));
        }
        WARN("AC16 well " << name << ": stop " << solve::to_string(r->stop) << ", " << r->passes.size() << " passes, "
                          << dofs(*r) << " DOFs, error " << error(*r) << ", certificate "
                          << solve::to_string(r->certificates[0].status) << ", misses " << r->misses << ", fallback "
                          << solve::to_string(r->fallback));
        return std::move(*r);
    };
    solve::AdaptOptions ob = o;
    ob.policy = solve::Policy::basic;
    const solve::AdaptResult b = run(ob, "basic");
    solve::AdaptOptions oi = o;
    oi.policy = solve::Policy::intelligent;
    const solve::AdaptResult i = run(oi, "intelligent");
    CHECK(b.certificates[0].status == solve::CertificateStatus::certified);
    CHECK(i.certificates[0].status == solve::CertificateStatus::certified);
    CHECK(error(b) <= o.tolerance);
    CHECK(error(i) <= o.tolerance);
    CHECK(static_cast<double>(dofs(i)) <= 0.8 * static_cast<double>(dofs(b)));
    // Budget mode at half basic's final DOFs; basic at the same budget.
    oi.intelligent.dof_budget = dofs(b) / 2;
    const solve::AdaptResult ib = run(oi, "intelligent, budget");
    ob.max_nodes = dofs(b) / 2 / 3;
    const solve::AdaptResult bb = run(ob, "basic, budget");
    CHECK(dofs(ib) <= dofs(b) / 2);
    CHECK(error(ib) <= error(bb));
}
