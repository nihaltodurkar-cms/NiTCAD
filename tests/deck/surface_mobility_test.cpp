// Lombardi surface mobility on devices (Unit 46; DECISIONS.md LS5-LS11): a 1D diode unchanged
// within 1e-4 (E_perp = E_0) and the run record; and, Release only (hidden tag [.lombardi]): the
// effective mobility of long-channel n-MOSFETs against the model at E_eff and its universality, the
// 2D MOSFET example's Id-Vg, the example reflected (x and y exchanged), the electrothermal limit and
// the small-signal conductances, with the pattern and sweep times reported; and the long channel
// past the minority rounding floor (Unit 50).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <filesystem>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/surface_mobility.hpp"
#include "NiTCAD/solve/newton.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
namespace fs = std::filesystem;

using deck_test::diode_device;
using deck_test::must;

namespace {

const auto lombardi = assemble::SurfaceMobility::lombardi;

// A plan's sweep with the surface mobility set, no files written.
results::Sweep sweep(deck::Plan plan, assemble::SurfaceMobility surface) {
    plan.options.models.surface_mobility = surface;
    plan.outputs = {};
    const auto report = deck::run(plan);
    if (!report) FAIL(report.error().message);
    REQUIRE(report->sweep.has_value());
    REQUIRE_FALSE(report->sweep->stopped.has_value());
    return std::move(*report->sweep);
}

deck::Plan example(const char* file) {
    const auto d = deck::read(fs::path(NITCAD_EXAMPLES_DIR) / file);
    if (!d) FAIL(d.error().message);
    auto p = deck::build(*d);
    if (!p) FAIL(p.error().message);
    return std::move(*p);
}

// The 2D MOSFET example with x and y exchanged (a reflection): the same device.
const std::string reflected_mosfet = R"(
[mesh]
type = tensor
x = linspace(0, 0.9um, 19) + linspace(0.9um, 0.95um, 6) + linspace(0.95um, 1um, 3)
y = linspace(0, 2um, 41)

[region silicon]
material = silicon

[region oxide]
material = sio2
box = 0.95um 1um  0 2um

[doping]
profile = acceptor constant 1e16
profile = donor gaussian 2e19 axis=x center=0.95um width=0.15um y=0:0.5um lateral=20nm
profile = donor gaussian 2e19 axis=x center=0.95um width=0.15um y=1.5um:2um lateral=20nm

[contact source]
kind = ohmic
patch = y_min
box = 0.6um 0.94um  0 0

[contact drain]
kind = ohmic
patch = y_max
box = 0.6um 0.94um  2um 2um

[contact gate]
kind = electrode
patch = x_max
box = 1um 1um  0.5um 1.5um
electrode = metal
work_function = 4.1

[contact body]
kind = ohmic
patch = x_min

[sweep]
contact gate = 0 : 1.5 : 0.1
contact drain = 0.1
)";

// A long-channel n-MOSFET: a 6 um body (p-type, N_A), n+ pockets of 0.5 um at both ends, a 10 nm
// oxide and a metal gate, 0.05 nm spacing in the top 5 nm of silicon (the classical inversion
// layer is a few tenths of a nm thick at 1e6 V/cm); Vd = 10 mV, Vg 0 .. 6.5 V.
std::string long_mosfet(const char* acceptors) {
    return std::string(R"(
[mesh]
type = tensor
x = linspace(0, 6um, 61)
y = linspace(0, 1.5um, 16) + linspace(1.5um, 1.9um, 9) + linspace(1.9um, 1.98um, 9) + linspace(1.98um, 1.995um, 16) + linspace(1.995um, 2um, 101) + linspace(2um, 2.01um, 6)

[region silicon]
material = silicon

[region oxide]
material = sio2
box = 0 6um  2um 2.01um

[doping]
profile = acceptor constant )") +
           acceptors + R"(
profile = donor gaussian 1e20 axis=y center=2um width=0.15um x=0:0.5um lateral=20nm
profile = donor gaussian 1e20 axis=y center=2um width=0.15um x=5.5um:6um lateral=20nm

[contact source]
kind = ohmic
patch = x_min
box = 0 0  1.7um 1.99um

[contact drain]
kind = ohmic
patch = x_max
box = 6um 6um  1.7um 1.99um

[contact gate]
kind = electrode
patch = y_max
box = 0.5um 5.5um  2.01um 2.01um
electrode = metal
work_function = 4.1

[contact body]
kind = ohmic
patch = y_min

[sweep]
contact gate = 0 : 6.5 : 0.5
contact drain = 0.01
)";
}

// One column of silicon nodes at x, sorted by y: (y, node).
std::vector<std::pair<double, std::size_t>> column(const device::Device& d, double x,
                                                   const results::BiasPoint& p) {
    std::vector<std::pair<double, std::size_t>> c;
    for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
        const auto& q = d.mesh().points()[i];
        if (std::abs(q[0] - x) > 1e-9 || !std::isfinite(p.bands.electron_fermi_eV[i])) continue;
        c.emplace_back(q[1], i);
    }
    std::ranges::sort(c);
    return c;
}

// The box-method integral of f(node) over an interior column [per cm^2 for densities]: each
// node's control volume over the column's 0.1 um width. The top silicon node's box reaches the
// interface, the bisector of its edge to the oxide, where the inversion electrons peak.
double integral(const device::Device& d, const std::vector<std::pair<double, std::size_t>>& c,
                const auto& f) {
    double s = 0.0;
    for (const auto& [y, i] : c) s += f(i) * d.mesh().volumes()[i] / 1e-5;
    return s;
}

struct Effective {
    double E_eff;  // [V/cm]
    double mu;     // [cm^2/(V s)]
};

// The effective mobility at mid-channel per gate bias above threshold: mu_eff = Id / (q N_inv
// dE_Fn/dx), N_inv the column's electrons, dE_Fn/dx the difference of the electron-weighted
// quasi-Fermi levels of the columns 0.1 um either side; E_eff = q (N_dep + N_inv / 2) / eps_Si.
std::vector<Effective> effective(const deck::Plan& plan, const results::Sweep& s, double NA) {
    const device::Device& d = plan.device;
    const double eps = physics::silicon_parameters.eps_r * base::eps0_F_per_cm;
    std::vector<Effective> out;
    for (const results::BiasPoint& p : s.points) {
        const auto mid = column(d, 3e-4, p);
        double fermi[2];
        for (const int side : {0, 1}) {
            const auto c = column(d, side == 0 ? 2.9e-4 : 3.1e-4, p);
            const double weight = integral(d, c, [&](std::size_t i) { return p.fields.n_cm3[i]; });
            fermi[side] = integral(d, c, [&](std::size_t i) {
                              return p.fields.n_cm3[i] * p.bands.electron_fermi_eV[i];
                          }) /
                          weight;
        }
        const double gradient = std::abs(fermi[1] - fermi[0]) / 0.2e-4;  // [V/cm]
        const double N_inv = integral(d, mid, [&](std::size_t i) { return p.fields.n_cm3[i]; });
        if (N_inv < 1e11) continue;  // below threshold
        const double N_dep = integral(d, mid, [&](std::size_t i) {
            return std::max(NA - p.fields.p_cm3[i], 0.0);
        });
        const double Id = std::abs(p.terminal_current[1]);  // A/cm
        out.push_back({base::q_C * (N_dep + 0.5 * N_inv) / eps,
                       Id / (base::q_C * N_inv * gradient)});
    }
    return out;
}

// LD3's mobility for electrons at E_perp = E and the substrate's N (Caughey-Thomas bulk).
double model(double E, double N) {
    const physics::Semiconductor si = physics::silicon();
    const double bulk = physics::caughey_thomas_mobility(si, physics::Carrier::electron, N, 300.0);
    return physics::lombardi_mobility(physics::lombardi_silicon_electrons, bulk, N, 300.0, E)
        .mobility;
}

// mu_eff at E, log-linear between the bracketing points; NaN outside.
double at(const std::vector<Effective>& v, double E) {
    for (std::size_t k = 1; k < v.size(); ++k) {
        const double a = v[k - 1].E_eff, b = v[k].E_eff;
        if ((a - E) * (b - E) > 0.0) continue;
        const double t = std::log(E / a) / std::log(b / a);
        return std::exp(std::log(v[k - 1].mu) + t * (std::log(v[k].mu) - std::log(v[k - 1].mu)));
    }
    return std::nan("");
}

// The long-channel curves with Lombardi per substrate, computed once for the LS6 and LS7 cases.
// Without the model the extraction gives the bulk mobility (checked within 1%).
const std::map<double, std::vector<Effective>>& long_channel() {
    static const std::map<double, std::vector<Effective>> curves = [] {
        std::map<double, std::vector<Effective>> out;
        for (const auto& [NA, text] : {std::pair{1e16, "1e16"}, std::pair{1e17, "1e17"}}) {
            const deck::Plan p = must(long_mosfet(text));
            const double bulk = physics::caughey_thomas_mobility(
                physics::silicon(), physics::Carrier::electron, NA, 300.0);
            for (const Effective& e : effective(p, sweep(p, assemble::SurfaceMobility::off), NA)) {
                std::printf("LS6 N_A %.0e without Lombardi: E_eff %.4e V/cm, mu_eff %.2f\n", NA,
                            e.E_eff, e.mu);
                CHECK(std::abs(e.mu / bulk - 1.0) <= 0.01);
            }
            out[NA] = effective(p, sweep(p, lombardi), NA);
        }
        return out;
    }();
    return curves;
}

// LS6 on one substrate: the worst gap to the model for E_eff 3e5 .. 1e6 V/cm.
void within_model(double NA) {
    double worst = 0.0;
    std::size_t checked = 0;
    for (const Effective& e : long_channel().at(NA)) {
        const double m = model(e.E_eff, NA);
        std::printf("LS6 N_A %.0e: E_eff %.4e V/cm, mu_eff %.2f, model %.2f, gap %+.2f%%\n", NA,
                    e.E_eff, e.mu, m, 100.0 * (e.mu / m - 1.0));
        if (e.E_eff < 3e5 || e.E_eff > 1e6) continue;
        worst = std::max(worst, std::abs(e.mu / m - 1.0));
        ++checked;
    }
    CAPTURE(NA, worst, checked);
    CHECK(checked >= 3);
    CHECK(worst <= 0.15);
}

}  // namespace

TEST_CASE("surface mobility: a 1D diode within 1e-4 (E_perp = E_0) and the run record (LS5)",
          "[deck]") {
    const deck::Plan p = must(diode_device + "[sweep]\ncontact anode = 0 : 0.6 : 0.2\n");
    const auto off = sweep(p, assemble::SurfaceMobility::off);
    const auto on = sweep(p, lombardi);
    for (std::size_t k = 1; k < off.points.size(); ++k) {
        const double a = off.points[k].terminal_current[0], b = on.points[k].terminal_current[0];
        CAPTURE(k, a, b);
        CHECK(std::abs(b / a - 1.0) <= 1e-4);
        CHECK(b != a);
    }
    const auto named = [](const results::RunRecord& r) {
        return std::ranges::any_of(r.settings,
                                   [](const auto& s) { return s.first == "models.surface_mobility"; });
    };
    CHECK(named(on.run));
    CHECK_FALSE(named(off.run));
    CHECK(on.run.input_identity != off.run.input_identity);
}

TEST_CASE("surface mobility: long-channel effective mobility, 1e17 substrate, and universality "
          "(LS6, LS7)",
          "[.lombardi]") {
    within_model(1e17);
    const auto& curves = long_channel();
    double worst = 0.0;
    for (const double E : {5e5, 6e5, 7e5, 8e5, 9e5, 1e6}) {
        const double a = at(curves.at(1e16), E), b = at(curves.at(1e17), E);
        std::printf("LS7 E_eff %.1e: mu_eff 1e16 %.2f, 1e17 %.2f\n", E, a, b);
        REQUIRE(std::isfinite(a));
        REQUIRE(std::isfinite(b));
        worst = std::max(worst, std::abs(a / b - 1.0));
    }
    CAPTURE(worst);
    CHECK(worst <= 0.15);
}

// LS6 on the 1e16 substrate is not met: within 15% to about 7.5e5 V/cm, 20.0% above the model at
// 9.2e5. The classical inversion layer's field spread (mu_eff, an average over the layer, exceeds
// mu(E_eff)); kept as written for the owner's ruling.
TEST_CASE("surface mobility: long-channel effective mobility, 1e16 substrate (LS6, not met)",
          "[.lombardi][!shouldfail]") {
    within_model(1e16);
}

// Unit 47 (DECISIONS.md QC12), reported, not gated: LS6's gap with MLDA quantum confinement on as
// well, which moves the inversion electrons off the interface and narrows the field spread. At the
// default tol_update (1e-8) the 1e17 sweep stops where Newton's correction stalls near 3e-8 with
// the residual at rounding; a sweep that stops reports why and its last iterations, and is run
// again at tol_update 1e-7 for the gaps.
TEST_CASE("surface mobility: LS6 with MLDA quantum confinement (Unit 47, QC12)", "[.mlda]") {
    for (const auto& [NA, text] : {std::pair{1e16, "1e16"}, std::pair{1e17, "1e17"}}) {
        deck::Plan p = must(long_mosfet(text));
        p.options.models.quantum_confinement = assemble::QuantumConfinement::mlda;
        p.options.models.surface_mobility = lombardi;
        p.outputs = {};
        auto report = deck::run(p);
        if (!report) FAIL(report.error().message);
        REQUIRE(report->sweep.has_value());
        if (const results::Sweep& s = *report->sweep; s.stopped) {
            std::printf("QC12 LS6 with MLDA, N_A %.0e: stopped after %zu of %zu points: %s\n", NA,
                        s.points.size(), p.points.size(), s.stopped->message.c_str());
            if (s.unfinished) {
                const auto& it = s.unfinished->iterations;
                for (std::size_t k = it.size() < 12 ? 0 : it.size() - 12; k < it.size(); ++k) {
                    std::printf("QC12 LS6 with MLDA, unfinished iteration %d: correction %.3e, "
                                "residual %.3e, step %.3f\n",
                                it[k].iteration, it[k].update, it[k].residual, it[k].step);
                }
            }
            p.options.newton.tol_update = 1e-7;
            report = deck::run(p);
            if (!report) FAIL(report.error().message);
            REQUIRE(report->sweep.has_value());
            REQUIRE_FALSE(report->sweep->stopped.has_value());
            std::printf("QC12 LS6 with MLDA, N_A %.0e: at tol_update 1e-7 all %zu points\n", NA,
                        report->sweep->points.size());
        }
        double worst = 0.0;
        for (const Effective& e : effective(p, *report->sweep, NA)) {
            const double m = model(e.E_eff, NA);
            std::printf("QC12 LS6 with MLDA, N_A %.0e: E_eff %.4e V/cm, mu_eff %.2f, model %.2f, "
                        "gap %+.2f%%\n",
                        NA, e.E_eff, e.mu, m, 100.0 * (e.mu / m - 1.0));
            if (e.E_eff >= 3e5 && e.E_eff <= 1e6) worst = std::max(worst, std::abs(e.mu / m - 1.0));
        }
        std::printf("QC12 LS6 with MLDA, N_A %.0e: worst gap in 3e5 .. 1e6 V/cm %.2f%%\n", NA,
                    100.0 * worst);
    }
}

TEST_CASE("surface mobility: the 2D MOSFET example, reflected, electrothermal, small signal "
          "(LS8-LS11)",
          "[.lombardi]") {
    const deck::Plan p = example("mosfet_2d.deck");
    const auto timed = [](const deck::Plan& plan, assemble::SurfaceMobility s) {
        const auto t0 = std::chrono::steady_clock::now();
        auto result = sweep(plan, s);
        const double seconds =
            std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        return std::pair{std::move(result), seconds};
    };
    const auto [off, t_off] = timed(p, assemble::SurfaceMobility::off);
    const auto [on, t_on] = timed(p, lombardi);
    std::printf("LS11 2D MOSFET example sweep: %.2f s off, %.2f s with Lombardi\n", t_off, t_on);
    {
        const auto scaling = assemble::make_scaling(p.device);
        REQUIRE(scaling.has_value());
        auto models = p.options.models;
        const auto base = assemble::DriftDiffusion::create(p.device, *scaling, models);
        models.surface_mobility = lombardi;
        const auto with = assemble::DriftDiffusion::create(p.device, *scaling, models);
        REQUIRE(base.has_value());
        REQUIRE(with.has_value());
        std::printf("LS11 2D MOSFET example Jacobian nonzeros: %zu off, %zu with Lombardi\n",
                    base->make_jacobian().nonzeros(), with->make_jacobian().nonzeros());
    }
    // Contacts: source, drain, gate, body; gate points 0 .. 1.5 by 0.1.
    const auto id = [](const results::Sweep& s, std::size_t k) {
        return s.points[k].terminal_current[1];
    };
    const std::size_t sub = 2, top = 15;  // Vg 0.2 V and 1.5 V
    REQUIRE(on.points.size() == 16);
    const double ratio = id(on, top) / id(off, top);
    std::printf("LS8 Id(1.5 V) %.6e A/cm with Lombardi, %.6e without, ratio %.4f; Id(0.2 V) ratio "
                "%.4f\n",
                id(on, top), id(off, top), ratio, id(on, sub) / id(off, sub));
    CHECK(std::abs(id(on, top)) < std::abs(id(off, top)));
    // The measured ratio, pinned (LS8); the subthreshold part of LS8 is its own case below.
    CHECK(std::abs(ratio / 0.694738 - 1.0) <= 1e-4);

    // LS4: Newton at Vg = 1.5 V converges quadratically over its last three corrections.
    const auto& steps = on.points[top].convergence.iterations;
    for (const auto& r : steps) std::printf("LS4 iteration %d: correction %.3e\n", r.iteration, r.update);
    REQUIRE(steps.size() >= 3);
    const double u1 = steps[steps.size() - 3].update, u2 = steps[steps.size() - 2].update,
                 u3 = steps.back().update;
    const double order = std::log(u3 / u2) / std::log(u2 / u1);
    std::printf("LS4 convergence order %.2f\n", order);
    CHECK(order >= 1.8);

    // LS9: the reflected device gives the same currents.
    const auto reflected = sweep(must(reflected_mosfet), lombardi);
    for (const std::size_t k : {std::size_t{5}, top}) {
        const double a = id(on, k), b = id(reflected, k);
        CAPTURE(k, a, b);
        CHECK(std::abs(b / a - 1.0) <= 1e-9);
    }

    // LS10: electrothermal with 300 K sinks on the body and the gate side within 1e-3 of it.
    auto thermal_text = std::string(R"(
[thermal_contact substrate]
patch = y_min
kind = isothermal
temperature = 300
)");
    const auto read = deck::read(fs::path(NITCAD_EXAMPLES_DIR) / "mosfet_2d.deck");
    REQUIRE(read.has_value());
    auto with_sink = *read;
    const auto extra = deck::parse(thermal_text, "sink.deck");
    REQUIRE(extra.has_value());
    for (const auto& s : extra->sections) with_sink.sections.push_back(s);
    auto plan_t = deck::build(with_sink);
    if (!plan_t) FAIL(plan_t.error().message);
    plan_t->options.models.electrothermal = true;
    plan_t->points = p.points;
    const auto hot = sweep(*plan_t, lombardi);
    const double thermal_ratio = hot.points.back().terminal_current[1] / id(on, top);
    std::printf("LS10 electrothermal Id(1.5 V) / isothermal %.6f\n", thermal_ratio);
    CHECK(std::abs(thermal_ratio - 1.0) <= 1e-3);

    // LS10: the small-signal drain conductance and transconductance against DC differences.
    solve::SmallSignalOptions ac;
    ac.steady = p.options;
    ac.steady.models.surface_mobility = lombardi;
    ac.frequencies_Hz = {0.0};
    std::vector<double> bias = p.points[top];
    const auto y = solve::solve_small_signal(p.device, p.points, ac);
    REQUIRE(y.has_value());
    REQUIRE_FALSE(y->stopped.has_value());
    REQUIRE(y->points.size() == p.points.size());
    const auto& Y = y->points[top].admittance[0];
    const std::size_t contacts = 4;
    const double h = 1e-4;
    for (const std::size_t j : {std::size_t{1}, std::size_t{2}}) {  // drain, gate
        deck::Plan q = p;
        auto up = bias, down = bias;
        up[j] += h;
        down[j] -= h;
        q.points = p.points;
        q.points.push_back(up);
        q.points.push_back(down);
        const auto dc = sweep(q, lombardi);
        const double g = (id(dc, top + 1) - id(dc, top + 2)) / (2.0 * h);
        const double ac_g = Y[1 * contacts + j].real();
        std::printf("LS10 dI_drain/dV_%zu: small signal %.8e, DC %.8e\n", j, ac_g, g);
        CHECK(std::abs(ac_g / g - 1.0) <= 1e-5);
    }
}

// LS8's subthreshold criterion is not met: Id(0.2 V) is 27.5% lower with Lombardi (0.718 on a mesh
// refined to 0.05 nm at the surface, so not a mesh effect). The depletion field at the surface,
// about 4e4 V/cm, already lowers mu_ac; the criterion assumed it would not. Kept as written for the
// owner's ruling.
TEST_CASE("surface mobility: the 2D MOSFET example below threshold within 10% (LS8, not met)",
          "[.lombardi][!shouldfail]") {
    const deck::Plan p = example("mosfet_2d.deck");
    const auto off = sweep(p, assemble::SurfaceMobility::off);
    const auto on = sweep(p, lombardi);
    const std::size_t sub = 2;  // Vg 0.2 V
    const double ratio = on.points[sub].terminal_current[1] / off.points[sub].terminal_current[1];
    std::printf("LS8 Id(0.2 V) ratio %.4f\n", ratio);
    CHECK(std::abs(ratio - 1.0) <= 0.10);
}

// Unit 50 (RC1, RC2; DECISIONS.md RD1, RD2): the 1e17 long channel with Lombardi, where minority
// holes in the n+ pockets (2.7e-19 of the electrons) set a rounding floor. With the minority floor
// on (opt-in) and otherwise default options each point converges without chattering above the
// tolerance, its drain current equal to the current before Unit 50 (recorded below, Release, this
// machine); with the residual line search as well every point solves (none did before).
TEST_CASE("surface mobility: the long channel converges past the minority rounding floor (Unit 50)",
          "[.lombardi]") {
    constexpr std::array<double, 14> before{
        5.1301549888304837e-10, 0.00018426809546428683, 0.001585448292749148, 0.0029335061273558266,
        0.0041126887024053004,  0.0051331463263241403,  0.0060143872798963825, 0.0067764951805582129,
        0.0074378601235738605,  0.0080144771468285518,  0.0085198801065484997, 0.0089653496180598057,
        0.0093602177692615161,  0.0097121831802286628};
    deck::Plan plan = must(long_mosfet("1e17"));
    plan.options.minority_floor = true;  // opt-in (RD2)
    const auto plain = sweep(plan, lombardi);
    REQUIRE(plain.points.size() == before.size());
    for (std::size_t k = 0; k < before.size(); ++k) {
        const auto& p = plain.points[k];
        std::size_t chatter = 0;
        for (const auto& it : p.convergence.iterations) chatter += it.update >= 1e-8 && it.update < 1e-4 ? 1 : 0;
        CAPTURE(k, chatter, p.convergence.iterations.size(), p.terminal_current[1], before[k]);
        CHECK(chatter <= 2);
        CHECK(std::abs(p.terminal_current[1] / before[k] - 1.0) <= 1e-9);
    }
    plan.options.newton.line_search.enabled = true;
    plan.options.newton.line_search.test = solve::LineSearchTest::residual;
    const auto searched = sweep(plan, lombardi);
    CHECK(searched.points.size() == before.size());
    // RC1 as fixed: the residual line search alone, the floor off (the default).
    plan.options.minority_floor = false;
    const auto searched_alone = sweep(plan, lombardi);
    CHECK(searched_alone.points.size() == before.size());
}
