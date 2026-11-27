// Unit 68 slice 4b (DECISIONS.md amendment 9): beta-Ga2O3 on a deck. The n-type-only guard (owner
// ruling 7), the refusal of field mobility without a hole saturation law (ruling 1), and MAC18 in
// the owner's form (ruling 6): a reverse-biased Ni Schottky contact (barrier 1.50 eV) on 20 um of
// n- (1e16) over 0.5 um of n+ (1e19), swept with the default solver settings until the peak field
// is the calculated breakdown field / 2 (Ghosh 2018, [010]: 4.8 / 2 MV/cm), every n and p finite
// and positive.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
using deck_test::must;
using deck_test::plan_of;

namespace {

std::string diode(const std::string& doping) {
    return R"(
[device]
temperature = 300

[mesh]
type = tensor
x = linspace(0, 2um, 401) + linspace(2um, 20um, 181) + linspace(20um, 20.5um, 51)

[region ga2o3]
material = ga2o3_beta

[doping]
)" + doping + R"(

[contact anode]
kind = schottky
patch = x_min
barrier = 1.50

[contact cathode]
kind = ohmic
patch = x_max

[sweep]
contact anode = 0
)";
}

const std::string mac18_doping = "profile = donor constant 1e16 x=-1um:20um\nprofile = donor constant 1e19 x=20um:21um";

double peak_field(const device::Device& d, const std::vector<double>& psi) {
    double peak = 0.0;
    for (std::size_t i = 0; i + 1 < psi.size(); ++i) {
        peak = std::max(peak, std::abs(psi[i + 1] - psi[i]) / (d.mesh().points()[i + 1][0] - d.mesh().points()[i][0]));
    }
    return peak;
}

}  // namespace

TEST_CASE("ga2o3_beta: a node whose acceptors exceed its donors is refused (owner ruling 7)", "[ga2o3][unit68]") {
    const auto p = plan_of(diode("profile = donor constant 1e16\nprofile = acceptor constant 2e16 x=5um:6um"));
    REQUIRE_FALSE(p.has_value());
    CHECK(p.error().message.find("n-type only") != std::string::npos);
    // Acceptors as compensation, below the donors: accepted.
    CHECK(plan_of(diode("profile = donor constant 1e16\nprofile = acceptor constant 5e15 x=5um:6um")).has_value());
}

TEST_CASE("ga2o3_beta: field mobility is refused, the holes having no saturation law (owner ruling 1)",
          "[ga2o3][unit68]") {
    deck::Plan p = must(diode(mac18_doping));
    p.options.models.field_mobility = true;
    const std::vector<std::vector<double>> zero{{0.0, 0.0}};
    const auto s = solve::sweep_bias(p.device, zero, p.options);
    REQUIRE_FALSE(s.has_value());
    CHECK(s.error().message.find("no hole saturation velocity") != std::string::npos);
    CHECK(s.error().message.find("'ga2o3'") != std::string::npos);
}

// MAC18 for beta-Ga2O3 (amendments 9, 11): with the default solver settings Newton stops near 190 V, the
// depletion-region holes (n_i about 1e-18 cm^-3) unresolved in plain densities until one underflows to
// 0. With the log-density step (owner, 2026-11-27: the opt-in step counts for MAC18) it reaches 1750 V.
TEST_CASE("ga2o3_beta MAC18: the Ni Schottky n- n+ diode reverse-biased to 2.4 MV/cm converges",
          "[ga2o3][unit68]") {
    deck::Plan p = must(diode(mac18_doping));
    p.options.log_densities = true;  // amendment 11 (owner: the opt-in step counts for MAC18)
    // The one-sided depletion estimate of the bias that puts 2.4 MV/cm at the contact, rounded up.
    const double E_target = 2.4e6;
    const double V_estimate = 10.87 * base::eps0_F_per_cm * E_target * E_target / (2.0 * base::q_C * 1e16);
    std::vector<std::vector<double>> points;
    // Small first steps (a 25 V step does not converge from equilibrium), then 10 V.
    for (double V : {0.0, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0}) points.push_back({-V, 0.0});
    for (double V = 20.0; V < V_estimate + 20.0; V += 10.0) points.push_back({-V, 0.0});
    const auto s = solve::sweep_bias(p.device, points, p.options);
    REQUIRE(s.has_value());
    if (s->stopped) std::printf("ga2o3 MAC18: stopped after %zu of %zu points: %s\n", s->points.size(), points.size(), s->stopped->message.c_str());
    REQUIRE_FALSE(s->stopped.has_value());
    REQUIRE(s->points.size() == points.size());
    for (const results::BiasPoint& b : s->points) {
        CHECK(std::ranges::all_of(b.fields.n_cm3, [](double v) { return std::isfinite(v) && v > 0.0; }));
        CHECK(std::ranges::all_of(b.fields.p_cm3, [](double v) { return std::isfinite(v) && v > 0.0; }));
    }
    const double peak = peak_field(p.device, s->points.back().fields.potential_V);
    std::printf("ga2o3 MAC18: %.0f V reverse, peak field %.3f MV/cm (target 2.4)\n", -points.back()[0], peak / 1e6);
    CHECK(peak >= E_target);
}

// The log-density step (DECISIONS.md amendment 11) on the MAC18 diode at low bias: the same
// equations, so the plain step's potential and electrons; and holes resolved to their own size,
// which the plain step does not (at 190 V its neutral n- holes sat near 1e-15 cm^-3, n_i^2 / N_D near 2e-52).
TEST_CASE("log densities: the Ga2O3 diode keeps the plain solution and resolves its holes",
          "[log_densities][unit68]") {
    std::vector<std::vector<double>> points;
    for (double V : {0.0, 0.01, 0.1, 0.5, 1.0, 2.0, 5.0, 10.0}) points.push_back({-V, 0.0});
    deck::Plan plain = must(diode(mac18_doping));
    deck::Plan log = must(diode(mac18_doping));
    log.options.log_densities = true;
    const auto a = solve::sweep_bias(plain.device, points, plain.options);
    const auto b = solve::sweep_bias(log.device, points, log.options);
    REQUIRE(a.has_value());
    REQUIRE(b.has_value());
    REQUIRE_FALSE(a->stopped.has_value());
    REQUIRE_FALSE(b->stopped.has_value());
    for (std::size_t k = 0; k < points.size(); ++k) {
        CAPTURE(k);
        const auto& f = a->points[k].fields;
        const auto& g = b->points[k].fields;
        for (std::size_t i = 0; i < f.potential_V.size(); ++i) {
            CHECK(std::abs(f.potential_V[i] - g.potential_V[i]) < 1e-9);
            if (f.n_cm3[i] > 1e6) CHECK(std::abs(g.n_cm3[i] / f.n_cm3[i] - 1.0) < 1e-6);
        }
    }
    // At 0 V, equilibrium: n p is n_i^2 at every node (Boltzmann statistics, no gap narrowing),
    // through 50 decades of holes.
    const auto& e = b->points[0].fields;
    double lo = e.n_cm3[0] * e.p_cm3[0], hi = lo;
    for (std::size_t i = 0; i < e.n_cm3.size(); ++i) {
        REQUIRE(e.p_cm3[i] > 0.0);
        lo = std::min(lo, e.n_cm3[i] * e.p_cm3[i]);
        hi = std::max(hi, e.n_cm3[i] * e.p_cm3[i]);
    }
    CHECK(hi / lo - 1.0 < 1e-6);
}

TEST_CASE("log densities: refused with the line search", "[log_densities][unit68]") {
    deck::Plan p = must(diode(mac18_doping));
    p.options.log_densities = true;
    p.options.newton.line_search.enabled = true;
    const std::vector<std::vector<double>> points{{0.0, 0.0}, {-1.0, 0.0}};
    const auto s = solve::sweep_bias(p.device, points, p.options);
    const bool refused = !s.has_value() || s->stopped.has_value();
    CHECK(refused);
}
