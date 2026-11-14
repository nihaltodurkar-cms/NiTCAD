// Unit 17c (DECISIONS.md DG1-DG12; DC1, DC2, DC7, DC8): derived adaptation goals.
// - DG12: solve's restated figures equal analysis/dc.hpp's, bit for bit, on n- and p-channel curves.
// - DC2: each sensitivity equals a central finite difference of the figure with its window fixed.
// - DC7: adapt()'s refusals; DC1: settings without a derived goal carry no new key.
// - A small end-to-end run on the coarse tensor MOSFET: windows per pass, certificates of their kind.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "NiTCAD/analysis/curve.hpp"
#include "NiTCAD/analysis/dc.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/device/doping_profiles.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/adapt.hpp"
#include "NiTCAD/solve/adapt_derived.hpp"

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
