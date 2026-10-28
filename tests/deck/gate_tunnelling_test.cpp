// Gate tunnelling on devices (Unit 51; DECISIONS.md GC8-GC12, GC15): a 1D MOS capacitor with a
// lumped n+ polysilicon gate over 1.5 nm of SiO2, self-consistent from -2 to +2 V, its ECB current
// against the Tsu-Esaki integral on the solution's own bands; the same capacitor with the oxide
// meshed; the processes' polarity and dominance on the nMOS and its pMOS mirror; a 2D device
// against 1D; the small-signal conductance against DC differences and a transient's conservation;
// and the deck, device file, run record and sweep file.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/gate.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/io/device_file.hpp"
#include "NiTCAD/io/files.hpp"
#include "NiTCAD/physics/gate_tunnelling.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/thermionic_emission.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "NiTCAD/solve/transient.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
using deck_test::must;

namespace {

constexpr double T = 300.0;
constexpr std::size_t gate = 0, substrate = 1;

struct Mos {
    char body = 'p';              // p: acceptors, n+ poly (nMOS); n: donors, p+ poly (pMOS)
    double doping = 5e17;         // [cm^-3]
    std::string tunnelling = "on";
    bool meshed = false;          // the oxide meshed at 0.05 nm under an electrode
    int dimension = 1;
};

// 1 um of silicon, 0.01 nm spacing over its first 2 nm; with `meshed` the oxide's 1.5 nm in front
// of it, the last oxide node and the first silicon node 5 pm either side of the interface (the
// box method's interface is their edge's midpoint).
std::string deck_of(const Mos& m) {
    const char* poly = m.body == 'p' ? "n_poly" : "p_poly";
    const char* dopant = m.body == 'p' ? "acceptor" : "donor";
    const char* across = m.dimension == 1 ? "" : "y = linspace(0, 1um, 3)\n";
    char text[4096];
    if (!m.meshed) {
        std::snprintf(text, sizeof text, R"(
[device]
temperature = 300

[mesh]
type = tensor
x = linspace(0, 2nm, 201) + linspace(2nm, 50nm, 97) + linspace(50nm, 1um, 96)
%s
[region silicon]
material = silicon

[doping]
profile = %s constant %g

[contact gate]
kind = gate
patch = x_min
oxide_thickness = 1.5nm
electrode = %s
tunnelling = %s

[contact substrate]
kind = ohmic
patch = x_max

[sweep]
contact gate = 0
)",
                      across, dopant, m.doping, poly, m.tunnelling.c_str());
    } else {
        std::snprintf(text, sizeof text, R"(
[device]
temperature = 300

[mesh]
type = tensor
x = linspace(0, 1.45nm, 30) + linspace(1.495nm, 1.505nm, 2) + linspace(1.505nm, 3.505nm, 201) + linspace(3.505nm, 51.505nm, 97) + linspace(51.505nm, 1.001505um, 96)
%s
[region oxide]
material = sio2
box = 0 1.5nm

[region silicon]
material = silicon
box = 1.5nm 2um

[doping]
profile = %s constant %g

[contact gate]
kind = electrode
patch = x_min
electrode = %s
tunnelling = %s

[contact substrate]
kind = ohmic
patch = x_max

[sweep]
contact gate = 0
)",
                      across, dopant, m.doping, poly, m.tunnelling.c_str());
    }
    return text;
}

// The gate bias of every point (the substrate at 0), from 0 V in steps toward each.
results::Sweep sweep(const deck::Plan& p, const std::vector<double>& gate_V) {
    std::vector<std::vector<double>> points;
    for (const double v : gate_V) points.push_back({v, 0.0});
    auto s = solve::sweep_bias(p.device, points, p.options);
    if (!s) FAIL(s.error().message);
    if (s->stopped) {
        std::printf("stopped after %zu points: %s\n", s->points.size(),
                    s->stopped->message.c_str());
    }
    REQUIRE_FALSE(s->stopped.has_value());
    REQUIRE(s->points.size() == gate_V.size());
    return std::move(*s);
}

std::vector<double> steps(double from, double to, double by) {
    std::vector<double> v;
    const int n = static_cast<int>(std::lround((to - from) / by));
    for (int k = 0; k <= n; ++k) v.push_back(from + by * k);
    return v;
}

double relative(double a, double b) { return std::abs(a - b) / std::abs(b); }

// The ECB current density [A/cm^2] of GD4 on a 1D lumped state's own band diagram and potential:
// the vacuum level depth - psi at the surface, the gate's E_FG + phi_m with E_FG = -V_G; the range
// from max(E_c, E_c,G) (n+ poly: E_c,G = E_FG) to the barrier top. A one-way flux (GR2, GR9) puts
// the other side's level far below the band, so only `from` supplies.
enum class Supply { both, semiconductor, electrode };
double ecb_of_state(const deck::Plan& p, const results::BiasPoint& b, Supply from = Supply::both) {
    const physics::Semiconductor& si = p.device.material(0);
    const device::GateStack& g = p.device.contacts()[gate].gate;
    const physics::InsulatorParameters& ox = physics::silicon_dioxide_tunnelling_parameters;
    const double depth = physics::intrinsic_level_depth_eV(si, T);
    const double vac_s = depth - b.fields.potential_V[0];
    const double FG = -b.bias_V[gate];
    const double vac_g = FG + assemble::gate_work_function_eV(g, si, T);
    const physics::BarrierSegment seg{vac_s - ox.electron_affinity_eV,
                                      vac_g - ox.electron_affinity_eV, g.oxide_thickness_cm,
                                      ox.tunnel_mass_n};
    const double richardson = base::q_C *
                              physics::emission_velocity_cm_s(si, physics::Carrier::electron, T) *
                              physics::conduction_band_dos(si, T) / (T * T);
    const double low = std::max(b.bands.conduction_eV[0], FG);
    const double Fn = b.bands.electron_fermi_eV[0];
    return physics::tsu_esaki_current(std::span(&seg, 1),
                                      {low, std::numeric_limits<double>::infinity()},
                                      from == Supply::electrode ? Fn - 1e3 : Fn,
                                      from == Supply::semiconductor ? FG - 1e3 : FG, T, richardson)
        .J;
}

}  // namespace

// The point at 0 V is thermal equilibrium (GC4): no current, nothing relative to compare.
TEST_CASE("Gate tunnelling: a self-consistent MOS capacitor from -2 to +2 V (GC8)",
          "[gate_tunnelling]") {
    const deck::Plan p = must(deck_of({}));
    for (const auto& range : {steps(0.0, 2.0, 0.1), steps(0.0, -2.0, -0.1)}) {
        const results::Sweep s = sweep(p, range);
        for (const results::BiasPoint& b : s.points) {
            if (b.bias_V[gate] == 0.0) continue;
            const double IG = b.terminal_current[gate];
            CAPTURE(b.bias_V[gate], IG);
            REQUIRE(b.tunnelling_ecb.size() == 2);
            CHECK(relative(b.tunnelling_ecb[gate] + b.tunnelling_hvb[gate] + b.tunnelling_evb[gate],
                           IG) <= 1e-12);
            // GR2, GR9: to 1e-12 of the larger one-way flux, the scale the inputs' rounding
            // resolves (the net current near quasi-equilibrium is their small difference).
            const double J = ecb_of_state(p, b);
            const double one_way = std::max(std::abs(ecb_of_state(p, b, Supply::semiconductor)),
                                            std::abs(ecb_of_state(p, b, Supply::electrode)));
            CAPTURE(J, one_way, b.tunnelling_ecb[gate]);
            CHECK(std::abs(b.tunnelling_ecb[gate] - J) <= 1e-12 * one_way);
        }
        const results::BiasPoint& last = s.points.back();
        std::printf("GC8 V_G %+.1f V: I_G %.4e A/cm^2 (ECB %.3e, HVB %.3e, EVB %.3e)\n",
                    last.bias_V[gate], last.terminal_current[gate], last.tunnelling_ecb[gate],
                    last.tunnelling_hvb[gate], last.tunnelling_evb[gate]);
    }
    // GR1, GR6: the conservation check with generation-recombination off, so the current is
    // tunnelling-limited and resolved at the substrate (with it on, the substrate's resolution is
    // three orders above the generation-limited current); 1e-9 A/cm^2 covers the substrate
    // current's rounding floor in accumulation (measured 7.6e-10).
    deck::Plan q = must(deck_of({}));
    q.options.models.srh = false;
    q.options.models.auger = false;
    double worst = 0.0;
    for (const auto& range : {steps(0.0, 2.0, 0.1), steps(0.0, -2.0, -0.1)}) {
        const results::Sweep s = sweep(q, range);
        for (const results::BiasPoint& b : s.points) {
            if (b.bias_V[gate] == 0.0) continue;
            const double IG = b.terminal_current[gate], IS = b.terminal_current[substrate];
            CAPTURE(b.bias_V[gate], IG, IS, b.terminal_current_resolution[substrate]);
            worst = std::max(worst, std::abs(IG + IS) / std::abs(IG));
            CHECK(std::abs(IG + IS) <= 1e-10 * std::abs(IG) + 1e-9);
        }
    }
    std::printf("GC8 conservation without generation-recombination: worst %.3e relative\n", worst);
}

// Quadratic, as SC11: in the asymptotic window each correction at most the square of the last.
TEST_CASE("Gate tunnelling: Newton converges quadratically at +1.5 V (GC7)", "[gate_tunnelling]") {
    deck::Plan p = must(deck_of({}));
    p.options.newton.tol_update = 1e-12;
    const results::Sweep s = sweep(p, [] {
        auto v = steps(0.0, 1.5, 0.1);
        v.push_back(1.501);  // a small last step: the asymptotic phase in the window
        return v;
    }());
    const auto& iterations = s.points.back().convergence.iterations;
    int pairs = 0;
    for (std::size_t k = 1; k < iterations.size(); ++k) {
        const double before = iterations[k - 1].update, after = iterations[k].update;
        std::printf("GC7 iteration %d: correction %.3e\n", iterations[k].iteration, after);
        if (before >= 1e-2 || after <= 1e-13) continue;
        CAPTURE(k, before, after);
        CHECK(after <= before * before);
        ++pairs;
    }
    CHECK(pairs >= 1);
}

TEST_CASE("Gate tunnelling: the meshed oxide against the lumped gate (GC9)", "[gate_tunnelling]") {
    const deck::Plan lumped = must(deck_of({}));
    // GR5: the capacitor is near-singular (its inversion electrons almost float); the meshed
    // path's ordering exposes that to the floating-region guard, so it is off here.
    deck::Plan meshed = must(deck_of({.meshed = true}));
    meshed.options.linear.min_pivot_ratio = 0.0;
    const results::Sweep a = sweep(lumped, steps(0.0, 2.0, 0.1));
    const results::Sweep b = sweep(meshed, steps(0.0, 2.0, 0.1));
    for (const std::size_t k : {std::size_t{10}, std::size_t{20}}) {
        const double IL = a.points[k].terminal_current[gate];
        const double IM = b.points[k].terminal_current[gate];
        std::printf("GC9 V_G %+.1f V: lumped %.5e, meshed %.5e A/cm^2 (%.3f%%)\n",
                    a.points[k].bias_V[gate], IL, IM, 100.0 * (IM - IL) / IL);
        CHECK(relative(IM, IL) <= 0.02);
    }
}

TEST_CASE("Gate tunnelling: polarity and dominance on the nMOS and its pMOS mirror (GC10)",
          "[gate_tunnelling]") {
    {  // nMOS in inversion: electrons leave the channel for the gate, current enters at the gate.
        const results::Sweep s = sweep(must(deck_of({})), steps(0.0, 1.5, 0.1));
        const results::BiasPoint& b = s.points.back();
        const double ecb = b.tunnelling_ecb[gate], hvb = b.tunnelling_hvb[gate],
                     evb = b.tunnelling_evb[gate];
        std::printf("GC10 nMOS +1.5 V: ECB %.3e, HVB %.3e, EVB %.3e\n", ecb, hvb, evb);
        CHECK(ecb > 0.0);
        CHECK(hvb >= 0.0);
        CHECK(evb >= 0.0);
        CHECK(ecb > std::abs(hvb));
        CHECK(ecb > std::abs(evb));
    }
    {  // pMOS in inversion: holes leave the channel for the gate, current leaves at the gate.
        const results::Sweep s = sweep(must(deck_of({.body = 'n'})), steps(0.0, -1.5, -0.1));
        const results::BiasPoint& b = s.points.back();
        const double ecb = b.tunnelling_ecb[gate], hvb = b.tunnelling_hvb[gate],
                     evb = b.tunnelling_evb[gate];
        std::printf("GC10 pMOS -1.5 V: ECB %.3e, HVB %.3e, EVB %.3e\n", ecb, hvb, evb);
        CHECK(hvb < 0.0);
        CHECK(ecb <= 0.0);
        CHECK(evb <= 0.0);
        CHECK(-hvb > std::abs(ecb));
        CHECK(-hvb > std::abs(evb));
    }
}

TEST_CASE("Gate tunnelling: a 2D device uniform across the gate gives the 1D current (GC11)",
          "[gate_tunnelling]") {
    const results::Sweep one = sweep(must(deck_of({})), steps(0.0, 1.0, 0.25));
    const results::Sweep two = sweep(must(deck_of({.dimension = 2})), steps(0.0, 1.0, 0.25));
    for (std::size_t k = 1; k < one.points.size(); ++k) {
        const double I1 = one.points[k].terminal_current[gate];
        const double I2 = two.points[k].terminal_current[gate] / 1e-4;  // per 1 um width
        CAPTURE(one.points[k].bias_V[gate], I1, I2);
        CHECK(relative(I2, I1) <= 1e-9);
    }
}

TEST_CASE("Gate tunnelling: the small-signal conductance and a transient (GC12)",
          "[gate_tunnelling]") {
    deck::Plan p = must(deck_of({}));
    p.options.newton.tol_update = 1e-12;
    // GR3, GR7: Y_GG at 0 Hz (the inversion layer responds in about 1e3 s, so 1 Hz is its
    // high-frequency response) against the central difference of the DC gate current, at -1 V:
    // in inversion the gate current is a 1e-5 cancellation the small-signal solve does not
    // resolve (its reported resolution covers the miss).
    const double V = -1.0, h = 1e-3;
    std::vector<double> path = steps(0.0, -0.9, -0.1);
    for (const double v : {V - h, V, V + h}) path.push_back(v);
    const results::Sweep s = sweep(p, path);
    const std::size_t at = path.size() - 3;
    const double dI =
        (s.points[at + 2].terminal_current[gate] - s.points[at].terminal_current[gate]) / (2.0 * h);
    solve::SmallSignalOptions o{.steady = p.options, .frequencies_Hz = {0.0}};
    const std::vector<std::vector<double>> point{{V, 0.0}};
    const auto y = solve::solve_small_signal(p.device, point, o);
    if (!y) FAIL(y.error().message);
    REQUIRE_FALSE(y->stopped.has_value());
    const double G = y->admittance(0, 0, gate, gate).real();
    std::printf("GC12 -1.0 V: G_GG(0 Hz) %.6e S/cm^2, DC difference %.6e\n", G, dI);
    CHECK(relative(G, dI) <= 1e-5);

    // A gate step 1.4 -> 1.5 V: the contacts' currents sum to zero at every step with a current
    // above 1e-3 A/cm^2, within 1e-7 (GR8; the transient's own floor without tunnelling is 2.5e-8).
    // GR1's setup: generation-recombination off, so the substrate's current is resolved.
    deck::Plan q = must(deck_of({}));  // the deck's own Newton tolerance, as GR8 was measured
    q.options.models.srh = false;
    q.options.models.auger = false;
    solve::TransientOptions t{.steady = q.options};
    t.t_end_s = 1e-9;
    t.dt_initial_s = 1e-13;
    const std::vector<solve::Waveform> w{*solve::Waveform::step(1.4, 1.5, 0.0),
                                         *solve::Waveform::step(0.0, 0.0, 0.0)};
    const auto r = solve::solve_transient(q.device, w, t);
    if (!r) FAIL(r.error().message);
    REQUIRE_FALSE(r->stopped.has_value());
    REQUIRE(r->points.size() > 2);
    double worst = 0.0;
    for (const auto& tp : r->points) {
        const double IG = tp.terminal_current[gate], IS = tp.terminal_current[substrate];
        const double largest = std::max(std::abs(IG), std::abs(IS));
        if (largest <= 1e-3) continue;
        CAPTURE(tp.time_s, IG, IS);
        worst = std::max(worst, std::abs(IG + IS) / largest);
        CHECK(std::abs(IG + IS) <= 1e-7 * largest);
    }
    std::printf("GC12 transient: worst relative sum %.3e\n", worst);
}

TEST_CASE("Gate tunnelling: deck, device file, run record and sweep file (GC15)",
          "[gate_tunnelling]") {
    std::string text = deck_of({.tunnelling = "ecb, evb"});
    text.insert(text.find("[contact substrate]"),
                "oxide_affinity = 1.05\noxide_band_gap = 9.0\noxide_mass_n = 0.5\n"
                "oxide_mass_p = 0.5\noxide_mass_evb = 0.42\n\n");
    const deck::Plan p = must(text);
    const device::Contact& g = p.device.contacts()[gate];
    using device::TunnelProcess;
    CHECK(g.tunnelling ==
          (static_cast<std::uint8_t>(TunnelProcess::ecb) | static_cast<std::uint8_t>(TunnelProcess::evb)));
    CHECK(g.gate.oxide_band_gap_eV == 9.0);
    CHECK(g.gate.oxide_mass_evb == 0.42);

    // An insulator region's band data, by its parameter names.
    const deck::Plan m = must(deck_of({.meshed = true}) + R"(
[region high_k]
material = insulator
eps_r = 22
electron_affinity_eV = 2.0
band_gap_eV = 5.8
tunnel_mass_n = 0.18
tunnel_mass_p = 0.2
box = 0 0.2nm
)");

    const auto dir = deck_test::scratch("gate_tunnelling_round_trip");
    for (const deck::Plan* plan : {&p, &m}) {
        const auto file = dir / "mos.nitdev";
        REQUIRE(io::save_device(file, plan->device).has_value());
        const auto back = io::load_device(file);
        REQUIRE(back.has_value());
        for (std::size_t c = 0; c < 2; ++c) {
            const device::Contact& a = plan->device.contacts()[c];
            const device::Contact& r = back->contacts()[c];
            CHECK(r.tunnelling == a.tunnelling);
            CHECK(r.gate.oxide_affinity_eV == a.gate.oxide_affinity_eV);
            CHECK(r.gate.oxide_band_gap_eV == a.gate.oxide_band_gap_eV);
            CHECK(r.gate.oxide_mass_n == a.gate.oxide_mass_n);
            CHECK(r.gate.oxide_mass_p == a.gate.oxide_mass_p);
            CHECK(r.gate.oxide_mass_evb == a.gate.oxide_mass_evb);
        }
        REQUIRE(back->regions().size() == plan->device.regions().size());
        for (std::size_t k = 0; k < back->regions().size(); ++k) {
            if (!device::is_insulator(back->regions()[k])) continue;
            CHECK(std::get<physics::Insulator>(back->regions()[k].material).parameters() ==
                  std::get<physics::Insulator>(plan->device.regions()[k].material).parameters());
        }
        // The record names the processes, and the identity follows them.
        const std::vector<std::vector<double>> points{{1.0, 0.0}};
        const results::RunRecord r = solve::make_run_record(plan->device, plan->options, points);
        const auto named = std::ranges::find(r.settings, std::string("tunnelling.gate"),
                                             &std::pair<std::string, double>::first);
        REQUIRE(named != r.settings.end());
        CHECK(named->second == static_cast<double>(plan->device.contacts()[gate].tunnelling));
        CHECK(solve::make_run_record(*back, plan->options, points).input_identity ==
              r.input_identity);
    }
    const deck::Plan off = must(deck_of({.tunnelling = "off"}));
    const deck::Plan on = must(deck_of({}));
    const std::vector<std::vector<double>> points{{1.0, 0.0}};
    CHECK(solve::make_run_record(off.device, off.options, points).input_identity !=
          solve::make_run_record(on.device, on.options, points).input_identity);
    CHECK(std::ranges::none_of(solve::make_run_record(off.device, off.options, points).settings,
                               [](const auto& s) { return s.first.starts_with("tunnelling"); }));

    // The per-process currents through the sweep file.
    const results::Sweep s = sweep(on, {0.5, 1.0});
    const auto file = dir / "mos.nitsweep";
    REQUIRE(io::save_sweep(file, s).has_value());
    const auto back = io::load_sweep(file);
    REQUIRE(back.has_value());
    for (std::size_t k = 0; k < s.points.size(); ++k) {
        CHECK(back->points[k].tunnelling_ecb == s.points[k].tunnelling_ecb);
        CHECK(back->points[k].tunnelling_hvb == s.points[k].tunnelling_hvb);
        CHECK(back->points[k].tunnelling_evb == s.points[k].tunnelling_evb);
    }
    // The CLI's table: three columns for the tunnelling gate, none for the substrate.
    const auto table = dir / "mos.csv";
    REQUIRE(deck::write_table(table, on.device, s).has_value());
    std::ifstream in(table);
    std::string header;
    std::getline(in, header);
    CHECK(header == "point,gate V,substrate V,gate A,substrate A,gate ECB A,gate HVB A,"
                    "gate EVB A,substeps");
}

