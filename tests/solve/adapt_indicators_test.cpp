// Gates for the refinement indicators and marking (Unit 17, slice 2; DECISIONS.md Unit 17 AD2, AD3,
// AR1; AC4): hand-computed edge values on a 4-node fixture; on an abrupt diode the doping and
// potential indicators mark the junction's depletion region only; on a MOS-C in inversion the
// electron-density indicator marks the inversion layer; marking rules and refusals.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/adapt_indicators.hpp"

using namespace NiTCAD;
using solve::Indicator;

namespace {

deck::Plan plan_of(const std::string& text) {
    auto d = deck::parse(text, "t.deck", {});
    if (!d) FAIL(d.error().message);
    auto p = deck::build(*d);
    if (!p) FAIL(p.error().message);
    return std::move(*p);
}

constexpr std::size_t K(Indicator k) { return static_cast<std::size_t>(k); }

// The signed log of net doping, AR1, written out again.
double L(double N, double ni) {
    const double m = std::log(std::max(std::abs(N), ni) / ni);
    return N < 0 ? -m : m;
}

const std::string four_nodes = R"(
[device]
temperature = 300
[mesh]
type = tensor
x = 0, 1um
y = 0, 1um
[region silicon]
material = silicon
[doping]
profile = acceptor constant 1e16 axis=y x=-1um:0.5um
profile = donor constant 1e18 axis=y x=0.5um:2um
[contact a]
kind = ohmic
patch = x_min
[equilibrium]
)";

// An abrupt diode along x, 1D in 2D; no node on the junction at 1 um.
const std::string diode = R"(
[device]
temperature = 300
[mesh]
type = tensor
x = linspace(0, 2um, 80)
y = 0, 0.1um, 0.2um
[region silicon]
material = silicon
[doping]
profile = acceptor constant 1e17 axis=y x=-1um:1um
profile = donor constant 1e17 axis=y x=1um:3um
[contact anode]
kind = ohmic
patch = x_min
[contact cathode]
kind = ohmic
patch = x_max
[equilibrium]
)";

// A MOS-C on p-type 1e17 under 1.5 nm of oxide, n+ poly, the gate at 1 V (strong inversion).
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

TEST_CASE("indicators: hand-computed edge values on four nodes", "[adapt_indicators]") {
    const deck::Plan p = plan_of(four_nodes);
    const device::Device& d = p.device;
    const auto& m = d.mesh();
    REQUIRE(m.node_count() == 4);
    const auto pts = m.points();
    // A state by hand.
    results::NodeFields s;
    s.potential_V.resize(4);
    s.n_cm3.resize(4);
    s.p_cm3.resize(4);
    for (std::size_t i = 0; i < 4; ++i) {
        const double x = pts[i][0] / 1e-4, y = pts[i][1] / 1e-4;  // in um
        s.potential_V[i] = 0.1 * x + 0.02 * y;
        s.n_cm3[i] = 1e10 * std::pow(10.0, 3.0 * x + 0.5 * y);
        s.p_cm3[i] = 1e15 * std::pow(10.0, -2.0 * x);
    }
    results::RateFields r;
    r.impact.assign(4, 0.0);
    r.impact[3] = -1e20;  // generation at one node only
    auto e = solve::edge_indicators(d, &s, &r);
    REQUIRE(e.has_value());
    const double vt = base::thermal_voltage(300.0);
    const double ni = physics::intrinsic_density(d.material(0), 300.0);
    const auto edges = m.edges();
    REQUIRE(edges.size() == 4);
    for (std::size_t k = 0; k < edges.size(); ++k) {
        const auto a = static_cast<std::size_t>(edges[k].first), b = static_cast<std::size_t>(edges[k].second);
        const bool along_x = pts[a][0] != pts[b][0];
        CHECK(std::abs(e->raw[K(Indicator::potential)][k] - (along_x ? 0.1 : 0.02) / vt) <= 1e-12 * 4.0);
        CHECK(std::abs(e->raw[K(Indicator::electron_density)][k] - (along_x ? 3.0 : 0.5)) <= 1e-12);
        CHECK(std::abs(e->raw[K(Indicator::hole_density)][k] - (along_x ? 2.0 : 0.0)) <= 1e-12);
        const double dop = along_x ? std::abs(L(-1e16, ni) - L(1e18, ni)) / std::log(10.0) : 0.0;
        CHECK(std::abs(e->raw[K(Indicator::doping)][k] - dop) <= 1e-12 * std::max(1.0, dop));
        // Generation: node 3 holds all of it; share times 4 nodes = 4 on its edges, 0 elsewhere.
        const bool touches3 = a == 3 || b == 3;
        CHECK(e->raw[K(Indicator::generation)][k] == (touches3 ? 4.0 : 0.0));
    }
    // Node excess: the largest over edges and kinds, r / t - 1.
    auto f = solve::compute_indicators(d, &s, &r);
    REQUIRE(f.has_value());
    for (std::size_t i = 0; i < 4; ++i) {
        double best = 0.0;
        for (std::size_t k = 0; k < edges.size(); ++k) {
            if (static_cast<std::size_t>(edges[k].first) != i && static_cast<std::size_t>(edges[k].second) != i) continue;
            best = std::max({best, e->raw[K(Indicator::potential)][k] - 1.0, e->raw[K(Indicator::electron_density)][k] - 1.0,
                             e->raw[K(Indicator::hole_density)][k] - 1.0, e->raw[K(Indicator::doping)][k] - 1.0,
                             e->raw[K(Indicator::generation)][k] / 10.0 - 1.0});
        }
        CHECK(std::abs(f->excess[i] - best) <= 1e-12 * std::max(1.0, best));
    }
    // Doping only, without a state.
    auto dop = solve::compute_indicators(d, nullptr, nullptr);
    REQUIRE(dop.has_value());
    for (std::size_t i = 0; i < 4; ++i) CHECK(dop->dominant[i] == Indicator::doping);
    // Thresholds move the excess.
    solve::IndicatorOptions o;
    o.potential_vt = 1e9;
    o.density_decades = 1e9;
    o.doping_decades = 1e9;
    o.generation_share = 1e9;
    auto none = solve::compute_indicators(d, &s, &r, o);
    REQUIRE(none.has_value());
    for (const double x : none->excess) CHECK(x == 0.0);
    // Refusals.
    o.potential_vt = 0.0;
    CHECK(!solve::edge_indicators(d, &s, &r, o).has_value());
    results::NodeFields bad = s;
    bad.n_cm3.pop_back();
    CHECK(!solve::edge_indicators(d, &bad, nullptr).has_value());
}

TEST_CASE("indicators: on an abrupt diode, doping and potential mark the depletion region only", "[adapt_indicators]") {
    deck::Plan p = plan_of(diode);
    auto run = deck::run(p);
    REQUIRE(run.has_value());
    REQUIRE(run->equilibrium.has_value());
    const device::Device& d = p.device;
    const auto pts = d.mesh().points();
    // The depletion width of a symmetric abrupt junction: W = sqrt(2 eps V_bi (2 / N) / q), half each side.
    const double N = 1e17;
    const double ni = physics::intrinsic_density(d.material(0), 300.0);
    const double vbi = base::thermal_voltage(300.0) * std::log(N * N / (ni * ni));
    const double eps = 11.7 * base::eps0_F_per_cm;
    const double W = std::sqrt(2.0 * eps * vbi * (2.0 / N) / base::q_C);
    const double h = 2e-4 / 79.0;
    const double xj = 1e-4;

    // Doping alone, before any solve: the two nodes either side of the junction only.
    solve::IndicatorOptions doping_only;
    doping_only.enabled = {false, false, false, true, false};
    auto fd = solve::compute_indicators(d, nullptr, nullptr, doping_only);
    REQUIRE(fd.has_value());
    auto md = solve::mark_nodes(d, *fd, {});
    REQUIRE(md.has_value());
    REQUIRE(md->count > 0);
    for (std::size_t i = 0; i < pts.size(); ++i) {
        if (md->marked[i] != 0) CHECK(std::abs(pts[i][0] - xj) < h);
    }

    // Potential alone on the equilibrium: inside the depletion region (one spacing of slack).
    solve::IndicatorOptions potential_only;
    potential_only.enabled = {true, false, false, false, false};
    auto fp = solve::compute_indicators(d, &run->equilibrium->fields, nullptr, potential_only);
    REQUIRE(fp.has_value());
    auto mp = solve::mark_nodes(d, *fp, {});
    REQUIRE(mp.has_value());
    REQUIRE(mp->count >= 6);  // the depletion region spans about 6 spacings
    std::size_t near_junction = 0;
    for (std::size_t i = 0; i < pts.size(); ++i) {
        if (mp->marked[i] == 0) continue;
        CHECK(std::abs(pts[i][0] - xj) <= 0.5 * W + h);
        near_junction += std::abs(pts[i][0] - xj) < h ? 1 : 0;
        CHECK(mp->by_kind[K(Indicator::potential)] == mp->count);
    }
    CHECK(near_junction > 0);
    // Quasi-neutral regions are resolved: zero excess beyond the depletion edges.
    for (std::size_t i = 0; i < pts.size(); ++i) {
        if (std::abs(pts[i][0] - xj) > 0.5 * W + 2.0 * h) CHECK(fp->excess[i] == 0.0);
    }
    // All kinds together keep the same support.
    auto fa = solve::compute_indicators(d, &run->equilibrium->fields, nullptr);
    REQUIRE(fa.has_value());
    auto ma = solve::mark_nodes(d, *fa, {});
    REQUIRE(ma.has_value());
    for (std::size_t i = 0; i < pts.size(); ++i) {
        if (ma->marked[i] != 0) CHECK(std::abs(pts[i][0] - xj) <= 0.5 * W + h);
    }
    CHECK(ma->by_region.size() == 1);
    CHECK(ma->by_region[0] == ma->count);
}

TEST_CASE("indicators: on a MOS-C in inversion, the electron density marks the inversion layer", "[adapt_indicators]") {
    deck::Plan p = plan_of(mosc);
    auto run = deck::run(p);
    REQUIRE(run.has_value());
    REQUIRE(run->sweep.has_value());
    REQUIRE(run->sweep->points.size() == 1);
    const results::NodeFields& s = run->sweep->points[0].fields;
    const device::Device& d = p.device;
    const auto pts = d.mesh().points();
    // The inversion layer: n above the acceptor density.
    double x_inv = 0.0;
    for (std::size_t i = 0; i < pts.size(); ++i) {
        if (s.n_cm3[i] > 1e17) x_inv = std::max(x_inv, pts[i][0]);
    }
    REQUIRE(x_inv > 0.0);
    solve::IndicatorOptions n_only;
    n_only.enabled = {false, true, false, false, false};
    auto f = solve::compute_indicators(d, &s, nullptr, n_only);
    REQUIRE(f.has_value());
    auto mk = solve::mark_nodes(d, *f, {});
    REQUIRE(mk.has_value());
    REQUIRE(mk->count > 0);
    const double h = 2e-7;  // 2 nm near the surface
    bool surface = false;
    for (std::size_t i = 0; i < pts.size(); ++i) {
        if (mk->marked[i] == 0) continue;
        CHECK(pts[i][0] <= x_inv + h * (1 + 1e-9));
        surface = surface || pts[i][0] == 0.0;
    }
    CHECK(surface);
}

TEST_CASE("indicators: marking, weights, union and refusals", "[adapt_indicators]") {
    const deck::Plan p = plan_of(diode);
    const device::Device& d = p.device;
    const std::size_t n = d.mesh().node_count();
    solve::IndicatorField f;
    f.excess.assign(n, 0.0);
    f.dominant.assign(n, Indicator::potential);
    f.excess[5] = 1.0;
    f.excess[6] = 0.3;
    f.excess[7] = 0.29;
    auto m = solve::mark_nodes(d, f, {}, 0.3);
    REQUIRE(m.has_value());
    CHECK(m->count == 2);
    CHECK(m->marked[5] == 1);
    CHECK(m->marked[6] == 1);
    CHECK(m->marked[7] == 0);
    CHECK(m->threshold == 0.3);
    // A weight moves the marks: node 7 weighted 10 is now the largest.
    std::vector<double> w(n, 1.0);
    w[7] = 10.0;
    auto mw = solve::mark_nodes(d, f, w, 0.3);
    REQUIRE(mw.has_value());
    CHECK(mw->count == 2);  // 2.9 and 1.0 >= 0.87
    CHECK(mw->marked[7] == 1);
    CHECK(mw->marked[6] == 0);
    // Zero excess is never marked, even at theta 1 with all zero.
    solve::IndicatorField z = f;
    std::ranges::fill(z.excess, 0.0);
    auto mz = solve::mark_nodes(d, z, {}, 1.0);
    REQUIRE(mz.has_value());
    CHECK(mz->count == 0);
    // Union: the larger excess per node, its dominant kind with it.
    solve::IndicatorField g = z;
    g.excess[6] = 2.0;
    g.dominant[6] = Indicator::doping;
    const std::vector<solve::IndicatorField> both{f, g};
    auto u = solve::unite(both);
    REQUIRE(u.has_value());
    CHECK(u->excess[5] == 1.0);
    CHECK(u->excess[6] == 2.0);
    CHECK(u->dominant[6] == Indicator::doping);
    CHECK(u->dominant[5] == Indicator::potential);
    // Refusals.
    CHECK(!solve::mark_nodes(d, f, {}, 0.0).has_value());
    CHECK(!solve::mark_nodes(d, f, {}, 1.5).has_value());
    CHECK(!solve::mark_nodes(d, f, std::vector<double>(3, 1.0), 0.3).has_value());
    std::vector<double> neg(n, 1.0);
    neg[0] = -1.0;
    CHECK(!solve::mark_nodes(d, f, neg, 0.3).has_value());
    CHECK(!solve::unite(std::span<const solve::IndicatorField>{}).has_value());
    solve::IndicatorField small;
    small.excess = {1.0};
    small.dominant = {Indicator::potential};
    const std::vector<solve::IndicatorField> mixed{f, small};
    CHECK(!solve::unite(mixed).has_value());
}
