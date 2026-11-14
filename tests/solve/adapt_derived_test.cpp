// Unit 17c (DECISIONS.md DG1-DG12; DC1, DC2, DC7, DC8): derived adaptation goals.
// - DG12: solve's restated figures equal analysis/dc.hpp's, bit for bit, on n- and p-channel curves.
// - DC2: each sensitivity equals a central finite difference of the figure with its window fixed.
// - DC7: adapt()'s refusals; DC1: settings without a derived goal carry no new key.
// - A small end-to-end run on the coarse tensor MOSFET: windows per pass, certificates of their kind.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/analysis/curve.hpp"
#include "NiTCAD/analysis/dc.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/device/doping_profiles.hpp"
#include "NiTCAD/mesh/locate.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/results/learn.hpp"
#include "NiTCAD/solve/adapt.hpp"
#include "NiTCAD/solve/adapt_adjoint.hpp"
#include "NiTCAD/solve/adapt_derived.hpp"
#include "NiTCAD/solve/adapt_learn.hpp"
#include "NiTCAD/solve/adapt_policy.hpp"
#include "NiTCAD/solve/bias.hpp"

using namespace NiTCAD;

namespace {

// A transfer curve: exponential below threshold (about 77 mV/decade), then the overdrive
// u = n V_T ln(1 + exp((V - 0.35) / n V_T)) with mobility degradation, I ~ u / (1 + 3 u): g_m has
// one clear peak above threshold, so a small change of a current keeps every window.
solve::TransferCurve analytic_curve(double sign, double drain_V) {
    solve::TransferCurve c;
    c.drain_V = sign * drain_V;
    for (int k = 0; k <= 24; ++k) {
        const double v = -0.2 + 0.05 * k;
        const double nvt = 1.3 * 0.025852;
        const double u = nvt * std::log1p(std::exp((v - 0.35) / nvt));
        const double i = 2e-6 * u / (1.0 + 3.0 * u) + 1e-14 * (1.0 + 0.3 * k);
        c.gate_V.push_back(sign * v);
        c.current.push_back(sign * i);
        c.resolution.push_back(1e-15);
    }
    return c;
}

analysis::Curve as_curve(const solve::TransferCurve& c) {
    return {.x = c.gate_V, .y = c.current, .resolution = c.resolution};
}

std::vector<double> uniform(double L, int nodes) {
    std::vector<double> a(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) a[static_cast<std::size_t>(i)] = i == nodes - 1 ? L : L * i / (nodes - 1);
    return a;
}

// adapt_intelligent_test.cpp's tensor MOSFET (AR6): Lg 600 nm, Lsd 300 nm, depth 200 nm, N_A 1e17,
// n+ source and drain, t_ox 5 nm, from 25 x 11 uniform nodes; contacts source, drain, body, gate.
constexpr double Lg = 6e-5, Lsd = 3e-5, depth = 2e-5, t_ox = 5e-7;

std::expected<solve::AdaptableMesh, base::Error> tensor_mosfet_mesh(const solve::MeshRefinement& r) {
    const std::vector<double> x = r.axes[0].empty() ? uniform(2 * Lsd + Lg, 25) : r.axes[0];
    const std::vector<double> y = r.axes[1].empty() ? uniform(depth, 11) : r.axes[1];
    auto m = mesh::make_tensor_grid(x, y);
    if (!m) return std::unexpected(m.error());
    auto cells = mesh::TensorCells::create(*m, x, y);
    if (!cells) return std::unexpected(cells.error());
    std::vector<mesh::NodeId> source, drain, gate;
    for (const mesh::NodeId v : m->find_boundary("y_min")->nodes) {
        const double px = m->points()[static_cast<std::size_t>(v)][0];
        (px <= Lsd ? source : (px >= Lsd + Lg ? drain : gate)).push_back(v);
    }
    auto body = m->find_boundary("y_max")->nodes;
    const double lateral = std::sqrt(2.0) * 1e-6;
    device::DopingProfile sd{.dopant = device::Dopant::donor, .peak_cm3 = 1e19, .shape = device::DepthShape::gaussian,
                             .axis = 1, .center = 0.0, .width = 5e-6, .lateral = lateral};
    sd.hi[0] = Lsd;
    device::DopingProfile dd = sd;
    dd.hi[0] = INFINITY;
    dd.lo[0] = Lsd + Lg;
    const std::vector<device::DopingProfile> profiles{sd, dd, {.dopant = device::Dopant::acceptor, .peak_cm3 = 1e17}};
    const std::vector<device::Region> regions{{"silicon", physics::silicon()}};
    const std::vector<device::RegionId> region(m->node_count(), 0);
    auto doping = device::sample(profiles, *m, region, regions);
    if (!doping) return std::unexpected(doping.error());
    device::Contact gc{"gate", device::ContactKind::gate, std::move(gate)};
    gc.gate = {.boundary = "y_min", .oxide_thickness_cm = t_ox};
    auto d = device::Device::create({.mesh = std::move(*m),
                                     .temperature_K = 300.0,
                                     .regions = regions,
                                     .node_region = region,
                                     .donors = doping->donors,
                                     .acceptors = doping->acceptors,
                                     .contacts = {{"source", device::ContactKind::ohmic, std::move(source)},
                                                  {"drain", device::ContactKind::ohmic, std::move(drain)},
                                                  {"body", device::ContactKind::ohmic, std::move(body)},
                                                  std::move(gc)},
                                     .cells = std::move(*cells)});
    if (!d) return std::unexpected(d.error());
    return solve::AdaptableMesh{.device = std::move(*d), .kind = solve::AdaptMeshKind::tensor};
}

// A gate sweep at V_d = 0.05 V: {source, drain, body, gate}.
std::vector<std::vector<double>> gate_sweep(double from, double to, int n) {
    std::vector<std::vector<double>> p;
    for (int k = 0; k < n; ++k) p.push_back({0.0, 0.05, 0.0, from + (to - from) * k / (n - 1)});
    return p;
}

solve::BiasOptions sweep_options() {
    solve::BiasOptions b;
    b.step_control.max_cuts = 20;
    return b;
}

}  // namespace

TEST_CASE("derived goals: the restated figures equal analysis/dc.hpp's, bit for bit (DG12)", "[adapt_derived]") {
    for (const double sign : {1.0, -1.0}) {
        INFO("sign " << sign);
        const solve::TransferCurve c = analytic_curve(sign, 0.05);
        // dc.hpp reads a p-channel curve scaled by -1 (currents and gate biases; the drain's sign too).
        solve::TransferCurve read = c;
        for (double& v : read.gate_V) v *= sign;
        for (double& i : read.current) i *= sign;
        const analysis::Curve curve = as_curve(read);
        auto vt = solve::derived_goal(solve::GoalKind::vt, c);
        auto gm = solve::derived_goal(solve::GoalKind::gm, c);
        auto sw = solve::derived_goal(solve::GoalKind::swing, c);
        REQUIRE(vt.has_value());
        REQUIRE(gm.has_value());
        REQUIRE(sw.has_value());
        auto ref_vt = analysis::threshold_max_gm(curve, sign * c.drain_V);
        auto ref_gm = analysis::max_transconductance(curve);
        auto ref_sw = analysis::subthreshold_swing(curve);
        REQUIRE(ref_vt.has_value());
        REQUIRE(ref_gm.has_value());
        REQUIRE(ref_sw.has_value());
        CHECK(vt->value == sign * ref_vt->value);
        CHECK(gm->value == ref_gm->value);
        CHECK(sw->value == ref_sw->value);
        CHECK(gm->first == ref_gm->window.first);
        CHECK(gm->last == ref_gm->window.last);
        CHECK(vt->first == ref_vt->window.first);
        CHECK(sw->first == ref_sw->window.first);
        CHECK(sw->last == ref_sw->window.last);
        // A sign of the device's own: V_T positive for this n-channel curve, negative for the p-channel one.
        CHECK(vt->value * sign > 0.0);
    }
}

TEST_CASE("derived goals: each sensitivity is the figure's finite difference, the window fixed (DC2)", "[adapt_derived]") {
    for (const double sign : {1.0, -1.0}) {
        for (const solve::GoalKind kind : {solve::GoalKind::vt, solve::GoalKind::gm, solve::GoalKind::swing}) {
            INFO("sign " << sign << ", " << solve::to_string(kind));
            const solve::TransferCurve c = analytic_curve(sign, 0.05);
            auto g = solve::derived_goal(kind, c);
            REQUIRE(g.has_value());
            REQUIRE(g->sensitivity.size() == g->last - g->first + 1);
            double resolution = 0.0;
            for (std::size_t k = g->first; k <= g->last; ++k) {
                const double h = 1e-6 * std::abs(c.current[k]);
                solve::TransferCurve up = c, down = c;
                up.current[k] += h;
                down.current[k] -= h;
                auto gu = solve::derived_goal(kind, up);
                auto gd = solve::derived_goal(kind, down);
                REQUIRE(gu.has_value());
                REQUIRE(gd.has_value());
                REQUIRE(gu->first == g->first);  // the window held
                REQUIRE(gd->first == g->first);
                const double fd = (gu->value - gd->value) / (2.0 * h);
                const double c_k = g->sensitivity[k - g->first];
                INFO("point " << k << ": analytic " << c_k << ", difference " << fd);
                CHECK(std::abs(fd - c_k) <= 1e-6 * std::abs(c_k) + 1e-12 * std::abs(g->value / c.current[k]));
                resolution += std::abs(c_k) * c.resolution[k];
            }
            CHECK(g->resolution == resolution);  // DG5
        }
    }
}

TEST_CASE("derived goals: a searched range admits only windows inside it (DR5)", "[adapt_derived]") {
    // g_m largest at point 1 (25); inside [3, 6] the candidates are 4 (10), 5 (8.5), 6 (2, one-sided).
    solve::TransferCurve c;
    c.gate_V = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6};
    c.current = {0.0, 2.0, 5.0, 6.0, 6.5, 8.0, 8.2};
    c.resolution.assign(7, 0.0);
    auto all = solve::derived_goal(solve::GoalKind::gm, c);
    REQUIRE(all.has_value());
    CHECK(all->point == 1);
    auto part = solve::derived_goal(solve::GoalKind::gm, c, 3, 6);
    REQUIRE(part.has_value());
    CHECK(part->first >= 3);
    CHECK(part->point == 4);
    CHECK(!solve::derived_goal(solve::GoalKind::gm, c, 3, 4).has_value());  // no point's slope inside
    // The swing (ceiling 1e-2: points 0..4): decades per segment 1, 2, 1, 0.7; the steepest segment
    // of the whole curve is 1..2, inside [2, 6] it is 2..3.
    solve::TransferCurve e = c;
    e.current = {1e-9, 1e-8, 1e-6, 1e-5, 5e-5, 0.5, 1.0};
    auto sa = solve::derived_goal(solve::GoalKind::swing, e);
    auto sp = solve::derived_goal(solve::GoalKind::swing, e, 2, 6);
    REQUIRE(sa.has_value());
    REQUIRE(sp.has_value());
    CHECK(sa->first == 1);
    CHECK(sp->first == 2);
}

TEST_CASE("derived goals: the figures' errors", "[adapt_derived]") {
    const solve::TransferCurve c = analytic_curve(1.0, 0.05);
    CHECK(!solve::derived_goal(solve::GoalKind::current, c).has_value());
    // A falling curve (an n-channel current falling with the gate): g_m nowhere positive, no rising segment.
    solve::TransferCurve falling = c;
    for (std::size_t k = 0; k < falling.current.size(); ++k) falling.current[k] = 1e-9 * std::exp(-0.5 * static_cast<double>(k));
    CHECK(!solve::derived_goal(solve::GoalKind::gm, falling).has_value());
    CHECK(!solve::derived_goal(solve::GoalKind::swing, falling).has_value());
    solve::TransferCurve short_curve = c;
    short_curve.gate_V.resize(1);
    short_curve.current.resize(1);
    short_curve.resolution.resize(1);
    CHECK(!solve::derived_goal(solve::GoalKind::vt, short_curve).has_value());
}

TEST_CASE("derived goals: adapt() refuses what DG1 and DG2 exclude (DC7)", "[adapt_derived]") {
    const auto refused = [](const std::vector<std::vector<double>>& pts, const solve::AdaptOptions& o) {
        auto r = solve::adapt(tensor_mosfet_mesh, pts, sweep_options(), o);
        REQUIRE(!r.has_value());
        return r.error().message;
    };
    solve::AdaptOptions o;
    o.goals = {{.contact = 1, .kind = solve::GoalKind::gm}};
    const auto pts = gate_sweep(0.0, 1.0, 5);
    {
        solve::AdaptOptions a = o;
        a.at = {2};
        CHECK(refused(pts, a).find("its window decides") != std::string::npos);
    }
    CHECK(refused(gate_sweep(0.0, 1.0, 2), o).find("at least 3 points") != std::string::npos);
    {
        auto two = pts;
        two[2][1] = 0.1;  // the drain moves too
        CHECK(refused(two, o).find("exactly one contact") != std::string::npos);
    }
    {
        // The body swept: an ohmic contact, not a gate or an electrode.
        std::vector<std::vector<double>> body;
        for (int k = 0; k < 4; ++k) body.push_back({0.0, 0.05, -0.1 * k, 0.5});
        CHECK(refused(body, o).find("gate or an electrode") != std::string::npos);
    }
    {
        solve::AdaptOptions a = o;
        a.vt_tolerance = 0.0;
        a.goals = {{.contact = 1, .kind = solve::GoalKind::vt}};
        CHECK(refused(pts, a).find("vt_tolerance") != std::string::npos);
    }
}

TEST_CASE("derived goals: settings and identity carry the kinds only with a derived goal (DC1, DG9)", "[adapt_derived]") {
    solve::AdaptOptions o;
    o.goals = {{.contact = 1}};
    for (const auto& [name, value] : solve::adapt_settings(o, false)) {
        CHECK(name != "adapt.goal_kind");
        CHECK(name != "adapt.vt_tolerance");
    }
    o.goals.push_back({.contact = 1, .kind = solve::GoalKind::swing});
    std::vector<double> kinds;
    bool vt_tolerance = false;
    for (const auto& [name, value] : solve::adapt_settings(o, false)) {
        if (name == "adapt.goal_kind") kinds.push_back(value);
        vt_tolerance = vt_tolerance || name == "adapt.vt_tolerance";
    }
    CHECK(kinds == std::vector<double>{0.0, 3.0});
    CHECK(vt_tolerance);
}

TEST_CASE("derived goals: a gm and a vt goal on the coarse tensor MOSFET (DG2, DG4-DG6)", "[adapt_derived]") {
    solve::AdaptOptions o;
    o.goals = {{.contact = 1, .kind = solve::GoalKind::gm}, {.contact = 1, .kind = solve::GoalKind::vt}};
    o.tolerance = 1e-2;
    o.vt_tolerance = 1e-2;
    o.max_passes = 2;
    const auto pts = gate_sweep(0.0, 1.2, 7);
    auto r = solve::adapt(tensor_mosfet_mesh, pts, sweep_options(), o);
    if (!r) FAIL(r.error().message);
    REQUIRE(!r->passes.empty());
    for (const solve::AdaptPass& p : r->passes) {
        REQUIRE(p.window_first.size() == 2);
        REQUIRE(p.goal.size() == 2);  // no current goal: the derived entries only
        for (std::size_t i = 0; i < 2; ++i) {
            CHECK(p.window_last[i] >= p.window_first[i]);
            CHECK(p.window_last[i] < pts.size());
        }
        // g_m and V_T share the largest-g_m window.
        CHECK(p.window_first[0] == p.window_first[1]);
    }
    REQUIRE(r->certificates.size() == 2);
    CHECK(r->certificates[0].kind == solve::GoalKind::gm);
    CHECK(r->certificates[1].kind == solve::GoalKind::vt);
    CHECK(r->certificates[1].allowed == 1e-2);  // DG5: vt's allowance absolute
    CHECK(r->certificates[0].allowed >= 1e-2 * std::abs(r->certificates[0].value));  // tol max(|G|, r)
    CHECK(r->points.empty());  // no current goal: no adaptation points of its own
    // V_T of a 1e17 / 5 nm MOSFET lies inside the sweep.
    CHECK(r->certificates[1].value > 0.0);
    CHECK(r->certificates[1].value < 1.2);
}

TEST_CASE("derived goals: deterministic, and the certificate independent of the policy and predictor (DC8)",
          "[adapt_derived]") {
    solve::AdaptOptions o;
    o.goals = {{.contact = 1, .kind = solve::GoalKind::gm}, {.contact = 1, .kind = solve::GoalKind::vt}};
    o.tolerance = 1e-2;
    o.vt_tolerance = 1e-2;
    o.weighting = solve::Weighting::adjoint;
    o.max_passes = 3;
    const auto pts = gate_sweep(0.0, 1.2, 7);
    const auto same = [](double x, double y) { return std::bit_cast<std::uint64_t>(x) == std::bit_cast<std::uint64_t>(y); };
    const auto run = [&](const solve::AdaptOptions& op) {
        auto r = solve::adapt(tensor_mosfet_mesh, pts, sweep_options(), op);
        if (!r) FAIL(r.error().message);
        return std::move(*r);
    };
    const auto same_certificates = [&](const solve::AdaptResult& a, const solve::AdaptResult& b) {
        REQUIRE(a.passes.size() == b.passes.size());
        for (std::size_t p = 0; p < a.passes.size(); ++p) CHECK(a.passes[p].mesh_digest == b.passes[p].mesh_digest);
        const auto pa = a.mesh->device.mesh().points(), pb = b.mesh->device.mesh().points();
        REQUIRE(pa.size() == pb.size());
        for (std::size_t k = 0; k < pa.size(); ++k) REQUIRE(pa[k] == pb[k]);
        REQUIRE(a.certificates.size() == b.certificates.size());
        for (std::size_t c = 0; c < a.certificates.size(); ++c) {
            const auto &x = a.certificates[c], &y = b.certificates[c];
            CHECK(x.kind == y.kind);
            CHECK(x.point == y.point);
            CHECK(x.window_moved == y.window_moved);
            CHECK(same(x.value, y.value));
            CHECK(same(x.estimate, y.estimate));
            CHECK(same(x.verified_error, y.verified_error));
            CHECK(same(x.order, y.order));
            CHECK(same(x.allowed, y.allowed));
            CHECK(x.status == y.status);
        }
    };
    const solve::AdaptResult b = run(o);
    // Two runs: bit-identical.
    same_certificates(b, run(o));
    // AC21 extended: the intelligent layer with empty predictions (every proposal a miss, AD3's marks).
    solve::AdaptOptions oi = o;
    oi.policy = solve::Policy::intelligent;
    oi.intelligent.prediction = [](const device::Device& d) {
        const std::size_t n = d.mesh().node_count();
        return solve::SizeTargets{std::vector<double>(n, INFINITY),
                                  std::vector<std::array<double, 3>>(n, std::array<double, 3>{})};
    };
    const solve::AdaptResult i = run(oi);
    same_certificates(b, i);
    CHECK(i.misses >= 1);
    // And with a learned predictor through learned_targets: its one tree splits on feature 12 (the
    // goal kind, DG10), 3 levels for a current goal and none for a derived one. Given the run's goal
    // kind it predicts nothing, so the start mesh and the certificates are the basic run's.
    results::LearnedModel model;
    auto start = tensor_mosfet_mesh({});
    REQUIRE(start.has_value());
    model.materials = solve::material_keys(start->device);
    model.trees = {{results::TreeNode{.feature = 12, .threshold = 0.5, .left = 1, .right = 2},
                    results::TreeNode{.value = 3.0}, results::TreeNode{.value = 0.0}}};
    model.tree_count = 1;
    const double gm_kind = static_cast<double>(static_cast<unsigned>(solve::GoalKind::gm));
    const auto finite = [](const solve::SizeTargets& t) {
        return std::ranges::count_if(t.size, [](double s) { return std::isfinite(s); });
    };
    {
        auto as_current = solve::learned_targets(model, start->device, pts, sweep_options(), 0.0);
        auto as_gm = solve::learned_targets(model, start->device, pts, sweep_options(), gm_kind);
        REQUIRE(as_current.has_value());
        REQUIRE(as_gm.has_value());
        CHECK(finite(*as_current) == static_cast<std::ptrdiff_t>(start->device.mesh().node_count()));
        CHECK(finite(*as_gm) == 0);
    }
    solve::AdaptOptions op = oi;
    op.intelligent.first_prediction = [&](const device::Device& d) {
        auto t = solve::learned_targets(model, d, pts, sweep_options(), gm_kind);
        REQUIRE(t.has_value());
        return *t;
    };
    same_certificates(b, run(op));
}

namespace {

// DC3's device and sweep: the tensor MOSFET at V_d 0.05 V, V_g -0.2 to 1.0 V in 0.05 V steps.
const std::vector<std::vector<double>> dc3_points = gate_sweep(-0.2, 1.0, 25);

// The three figures on the tensor MOSFET refined uniformly by `factor` (each interval split evenly).
std::array<double, 3> uniform_figures(int factor, int ky = 0) {
    if (ky == 0) ky = factor;
    solve::MeshRefinement r;
    r.axes[0] = uniform(2 * Lsd + Lg, 24 * factor + 1);
    r.axes[1] = uniform(depth, 10 * ky + 1);
    auto m = tensor_mosfet_mesh(r);
    REQUIRE(m.has_value());
    auto s = solve::sweep_bias(m->device, dc3_points, sweep_options());
    if (!s) FAIL(s.error().message);
    solve::TransferCurve c;
    c.drain_V = 0.05;
    for (std::size_t k = 0; k < s->points.size(); ++k) {
        const auto& p = s->points[k];
        c.gate_V.push_back(dc3_points[k][3]);
        c.current.push_back(p.terminal_current[1]);
        c.resolution.push_back(p.terminal_current_resolution.size() > 1 ? p.terminal_current_resolution[1] : 0.0);
    }
    std::array<double, 3> out{};
    const solve::GoalKind kinds[3] = {solve::GoalKind::vt, solve::GoalKind::gm, solve::GoalKind::swing};
    for (int i = 0; i < 3; ++i) {
        auto g = solve::derived_goal(kinds[i], c);
        if (!g) FAIL(g.error().message);
        out[static_cast<std::size_t>(i)] = g->value;
        std::printf("  factor %2d, %s: %.12g (window %zu..%zu)\n", factor, solve::to_string(kinds[i]), g->value, g->first,
                    g->last);
    }
    std::printf("factor %d: %zu nodes\n", factor, m->device.mesh().node_count());
    std::fflush(stdout);
    return out;
}

}  // namespace

TEST_CASE("derived goals: probe the tensor MOSFET's uniform ladder for DC3's reference", "[.adapt_derived_ref]") {
    std::vector<std::array<double, 3>> f;
    for (const int factor : {1, 2, 4, 8, 16}) f.push_back(uniform_figures(factor));
    const char* names[3] = {"vt", "gm", "swing"};
    for (std::size_t i = 0; i < 3; ++i) {
        // Richardson on 4x / 8x / 16x: the observed order and the extrapolated value.
        const double d1 = f[3][i] - f[2][i], d2 = f[4][i] - f[3][i];
        const double p = std::log2(std::abs(d1 / d2));
        const double ext = f[4][i] + d2 / (std::exp2(p) - 1.0);
        std::printf("%s: 4x %.12g, 8x %.12g, 16x %.12g, order %.3f, extrapolated %.12g\n", names[i], f[2][i], f[3][i],
                    f[4][i], p, ext);
    }
}

// The ladder per axis (as AR9: the uniform ladder is not asymptotic on this device). Only the
// windows' neighbourhoods are solved: points 2..5 (swing, window 3..4 on 4x..16x) and 11..17 (vt
// and gm, window 13..15). Swing is taken on the fixed window 3..4 (its ceiling needs the whole
// curve); vt and gm on the sub-curve 11..17, their window required interior.
TEST_CASE("derived goals: probe the tensor MOSFET's ladder per axis for DC3's reference", "[.adapt_derived_ref2]") {
    std::vector<std::pair<int, int>> ladder{{4, 4}, {4, 32}, {8, 32}, {4, 64}, {16, 32}, {8, 64}, {4, 128}, {32, 32}};
    std::vector<std::vector<double>> pts;
    std::vector<std::size_t> index;
    for (std::size_t k = 2; k <= 5; ++k) index.push_back(k);
    for (std::size_t k = 11; k <= 17; ++k) index.push_back(k);
    for (const std::size_t k : index) pts.push_back(dc3_points[k]);
    for (const auto [kx, ky] : ladder) {
        solve::MeshRefinement r;
        r.axes[0] = uniform(2 * Lsd + Lg, 24 * kx + 1);
        r.axes[1] = uniform(depth, 10 * ky + 1);
        auto m = tensor_mosfet_mesh(r);
        REQUIRE(m.has_value());
        const auto t0 = std::chrono::steady_clock::now();
        auto s = solve::sweep_bias(m->device, pts, sweep_options());
        if (!s) FAIL(s.error().message);
        REQUIRE(s->points.size() == pts.size());
        const auto I = [&](std::size_t j) { return s->points[j].terminal_current[1]; };
        const double swing = 1000.0 * (dc3_points[4][3] - dc3_points[3][3]) / (std::log10(I(2)) - std::log10(I(1)));
        solve::TransferCurve c;
        c.drain_V = 0.05;
        for (std::size_t j = 4; j < pts.size(); ++j) {
            c.gate_V.push_back(pts[j][3]);
            c.current.push_back(I(j));
            c.resolution.push_back(0.0);
        }
        auto vt = solve::derived_goal(solve::GoalKind::vt, c);
        auto gm = solve::derived_goal(solve::GoalKind::gm, c);
        if (!vt) FAIL(vt.error().message);
        if (!gm) FAIL(gm.error().message);
        CHECK(gm->first > 0);
        CHECK(gm->last + 1 < c.gate_V.size());
        const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::printf("ladder kx=%d ky=%d nodes %zu: vt %.12g gm %.12g swing %.12g (gm window %zu..%zu of 11..17) %.0f s\n",
                    kx, ky, m->device.mesh().node_count(), vt->value, gm->value, swing, gm->first + 11, gm->last + 11,
                    sec);
        std::fflush(stdout);
    }
}

// DC3's runs, printed pass by pass (a probe: the path and the cost before the gates are written).
TEST_CASE("derived goals: probe DC3's adaptive runs", "[.adapt_derived_trial]") {
    const solve::GoalKind kinds[3] = {solve::GoalKind::vt, solve::GoalKind::gm, solve::GoalKind::swing};
    for (const auto kind : kinds) {
        solve::AdaptOptions o;
        o.goals = {{.contact = 1, .kind = kind}};
        o.tolerance = 1e-2;
        o.vt_tolerance = 1e-3;
        o.weighting = solve::Weighting::adjoint;
        o.max_passes = 12;
        const auto t0 = std::chrono::steady_clock::now();
        auto r = solve::adapt(tensor_mosfet_mesh, dc3_points, sweep_options(), o);
        if (!r) FAIL(r.error().message);
        const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        for (std::size_t p = 0; p < r->passes.size(); ++p) {
            const auto& a = r->passes[p];
            std::printf("  %s pass %zu: nodes %zu dofs %zu value %.12g window %zu..%zu estimate %.4g\n",
                        solve::to_string(kind), p, a.nodes, a.dofs, a.goal.empty() ? NAN : a.goal.back(),
                        a.window_first.empty() ? 0 : a.window_first.back(),
                        a.window_last.empty() ? 0 : a.window_last.back(), a.estimate.empty() ? NAN : a.estimate.back());
        }
        REQUIRE(r->certificates.size() == 1);
        const auto& c = r->certificates[0];
        std::printf("%s: %s value %.12g estimate %.4g verified %.4g allowed %.4g moved %d stop %d, %zu passes, "
                    "verification %zu nodes, %.0f s\n",
                    solve::to_string(kind), solve::to_string(c.status), c.value, c.estimate, c.verified_error,
                    c.allowed, c.window_moved ? 1 : 0, static_cast<int>(r->stop), r->passes.size(),
                    r->verification_nodes, sec);
        std::fflush(stdout);
    }
}

// The adjoint's D_k against the actual change, 4x -> 8x uniform, in the swing's window (a probe).
TEST_CASE("derived goals: probe the adjoint differences in subthreshold", "[.adapt_derived_dprobe]") {
    const auto make = [](int k) {
        solve::MeshRefinement r;
        r.axes[0] = uniform(2 * Lsd + Lg, 24 * k + 1);
        r.axes[1] = uniform(depth, 10 * k + 1);
        auto m = tensor_mosfet_mesh(r);
        REQUIRE(m.has_value());
        return std::move(*m);
    };
    const auto c = make(4), f = make(8);
    const std::vector<std::vector<double>> pts(dc3_points.begin(), dc3_points.begin() + 16);
    auto sc = solve::sweep_bias(c.device, pts, sweep_options());
    auto sf = solve::sweep_bias(f.device, pts, sweep_options());
    REQUIRE(sc.has_value());
    REQUIRE(sf.has_value());
    const auto loc = mesh::PointLocator::create(*c.device.cells());
    for (const std::size_t k : {std::size_t{2}, std::size_t{3}, std::size_t{4}, std::size_t{14}}) {
        const std::vector<std::size_t> goals{1};
        auto e = solve::estimate_goals(c.device, loc, sc->points[k], f.device, sweep_options(), goals);
        if (!e) FAIL(e.error().message);
        const double I4 = sc->points[k].terminal_current[1], I8 = sf->points[k].terminal_current[1];
        std::printf("point %zu: I4 %.6e I8 %.6e actual %.4e D %.4e prolongated %.4e res %.3e\n", k, I4, I8, I8 - I4,
                    e->difference[0], e->prolongated_part[0],
                    sc->points[k].terminal_current_resolution.size() > 1 ? sc->points[k].terminal_current_resolution[1]
                                                                         : -1.0);
    }
}

// The same on the swing run's adapted meshes: each pass's tensor axes, halved, as M'.
TEST_CASE("derived goals: probe the adjoint differences on adapted meshes", "[.adapt_derived_dprobe2]") {
    std::vector<std::array<std::vector<double>, 2>> axes;
    solve::AdaptOptions o;
    o.goals = {{.contact = 1, .kind = solve::GoalKind::swing}};
    o.tolerance = 1e-2;
    o.weighting = solve::Weighting::adjoint;
    o.max_passes = 8;
    o.on_mesh = [&](std::size_t, const solve::AdaptableMesh& m) -> std::expected<void, base::Error> {
        std::array<std::vector<double>, 2> a;
        for (const auto& p : m.device.mesh().points()) {
            a[0].push_back(p[0]);
            a[1].push_back(p[1]);
        }
        for (auto& v : a) {
            std::ranges::sort(v);
            v.erase(std::unique(v.begin(), v.end()), v.end());
        }
        axes.push_back(std::move(a));
        return {};
    };
    auto r = solve::adapt(tensor_mosfet_mesh, dc3_points, sweep_options(), o);
    if (!r) FAIL(r.error().message);
    const std::vector<std::vector<double>> pts(dc3_points.begin(), dc3_points.begin() + 6);
    for (std::size_t p = 0; p < axes.size(); p += 2) {
        solve::MeshRefinement rc, rf;
        rc.axes[0] = axes[p][0];
        rc.axes[1] = axes[p][1];
        for (int a = 0; a < 2; ++a) {
            const auto& v = axes[p][static_cast<std::size_t>(a)];
            for (std::size_t i = 0; i + 1 < v.size(); ++i) {
                rf.axes[static_cast<std::size_t>(a)].push_back(v[i]);
                rf.axes[static_cast<std::size_t>(a)].push_back(0.5 * (v[i] + v[i + 1]));
            }
            rf.axes[static_cast<std::size_t>(a)].push_back(v.back());
        }
        auto c = tensor_mosfet_mesh(rc);
        auto f = tensor_mosfet_mesh(rf);
        REQUIRE(c.has_value());
        REQUIRE(f.has_value());
        auto sc = solve::sweep_bias(c->device, pts, sweep_options());
        auto sf = solve::sweep_bias(f->device, pts, sweep_options());
        REQUIRE(sc.has_value());
        REQUIRE(sf.has_value());
        const auto loc = mesh::PointLocator::create(*c->device.cells());
        std::printf("pass %zu: %zu x %zu nodes\n", p, axes[p][0].size(), axes[p][1].size());
        for (const std::size_t k : {std::size_t{2}, std::size_t{3}, std::size_t{4}}) {
            const std::vector<std::size_t> goals{1};
            auto e = solve::estimate_goals(c->device, loc, sc->points[k], f->device, sweep_options(), goals);
            if (!e) FAIL(e.error().message);
            const double I = sc->points[k].terminal_current[1], If = sf->points[k].terminal_current[1];
            std::printf("  point %zu: I %.6e I' %.6e actual/I %.4e D/I %.4e prolongated %.4e\n", k, I, If, (If - I) / I,
                        e->difference[0] / I, e->prolongated_part[0]);
        }
        std::fflush(stdout);
    }
}

// DC4's baseline (DQ4 (a)): `goal = current drain` on the window points, adjoint-weighted, the
// tolerance tightened; the three figures extracted from each final sweep.
TEST_CASE("derived goals: probe DC4's current-goal baseline", "[.adapt_derived_baseline]") {
    const std::vector<std::pair<const char*, std::vector<std::size_t>>> windows{{"vt/gm", {13, 14, 15}},
                                                                                 {"swing", {3, 4}}};
    for (const auto& [name, at] : windows) {
        for (const double tol : {1e-2, 3e-3, 1e-3, 3e-4}) {
            solve::AdaptOptions o;
            o.goals = {{.contact = 1}};
            o.at = at;
            o.tolerance = tol;
            o.weighting = solve::Weighting::adjoint;
            o.max_passes = 12;
            const auto t0 = std::chrono::steady_clock::now();
            auto r = solve::adapt(tensor_mosfet_mesh, dc3_points, sweep_options(), o);
            if (!r) FAIL(r.error().message);
            const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
            solve::TransferCurve c;
            c.drain_V = 0.05;
            for (std::size_t k = 0; k < r->sweep.points.size(); ++k) {
                c.gate_V.push_back(dc3_points[k][3]);
                c.current.push_back(r->sweep.points[k].terminal_current[1]);
                c.resolution.push_back(0.0);
            }
            auto vt = solve::derived_goal(solve::GoalKind::vt, c);
            auto gm = solve::derived_goal(solve::GoalKind::gm, c);
            auto sw = solve::derived_goal(solve::GoalKind::swing, c);
            if (!vt || !gm || !sw) FAIL("extraction");
            bool certified = true;
            for (const auto& g : r->certificates) certified = certified && g.status == solve::CertificateStatus::certified;
            std::printf("baseline %s tol %.0e: dofs %zu nodes %zu passes %zu certified %d: vt %.12g gm %.12g "
                        "(window %zu..%zu) swing %.12g (window %zu..%zu) %.0f s\n",
                        name, tol, r->passes.back().dofs, r->passes.back().nodes, r->passes.size(), certified ? 1 : 0,
                        vt->value, gm->value, gm->first, gm->last, sw->value, sw->first, sw->last, sec);
            std::fflush(stdout);
        }
    }
}

TEST_CASE("derived goals: a window that moves on M' is never certified (DC6, DG6)", "[.adapt_derived_slow]") {
    // Found by a probe (2026-11-14): at a 0.02 V gate step the subthreshold swing is nearly the same on
    // every segment, so its steepest segment moves between the passes and onto M'.
    solve::AdaptOptions o;
    o.goals = {{.contact = 1, .kind = solve::GoalKind::swing}};
    o.tolerance = 1e-1;
    o.max_passes = 4;
    auto r = solve::adapt(tensor_mosfet_mesh, gate_sweep(-0.2, 0.2, 21), sweep_options(), o);
    if (!r) FAIL(r.error().message);
    REQUIRE(r->certificates.size() == 1);
    const auto& c = r->certificates[0];
    CHECK(c.window_moved);
    CHECK(c.status == solve::CertificateStatus::not_certified);
    CHECK(r->rejected_stops >= 1);  // the moved window rejected a proposed stop (AR2)
    // The windows moved between passes too: no stop was proposed until two passes agreed.
    CHECK(r->passes.front().window_first[0] != r->passes.back().window_first[0]);
}

// Per-axis differences on coarse meshes (a probe for DR6): kx x ky uniform factors.
TEST_CASE("derived goals: probe per-axis differences", "[.adapt_derived_axes]") {
    for (const auto [kx, ky] : std::vector<std::pair<int, int>>{{1, 1}, {2, 1}, {1, 2}, {2, 2}, {4, 2}, {2, 4}, {4, 4},
                                                                {8, 4}, {4, 8}}) {
        (void)uniform_figures(kx, ky);
    }
}
