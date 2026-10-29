// Gate tunnelling on devices (Unit 51; DECISIONS.md GC8-GC12, GC15): a 1D MOS capacitor with a
// lumped n+ polysilicon gate over 1.5 nm of SiO2, self-consistent from -2 to +2 V, its ECB current
// against the Tsu-Esaki integral on the solution's own bands; the same capacitor with the oxide
// meshed; the processes' polarity and dominance on the nMOS and its pMOS mirror; a 2D device
// against 1D, and a 3D corner's two paths against two 2D devices; the small-signal conductance against DC differences and a transient's conservation;
// and the deck, device file, run record and sweep file.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <sstream>
#include <fstream>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/gate.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/io/device_file.hpp"
#include "NiTCAD/io/files.hpp"
#include "NiTCAD/linalg/mkl_runtime.hpp"
#include "NiTCAD/physics/gate_tunnelling.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/thermionic_emission.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/equilibrium.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "NiTCAD/solve/transient.hpp"
#include "../solve/assembly_timing.hpp"
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
// At another node than 0 (a gate at x_max): that node's material, the vacuum level from the
// reference node's intrinsic level.
enum class Supply { both, semiconductor, electrode };
double ecb_of_state(const deck::Plan& p, const results::BiasPoint& b, Supply from = Supply::both,
                    std::size_t node = 0) {
    const physics::Semiconductor& si = p.device.material(static_cast<mesh::NodeId>(node));
    const device::GateStack& g = p.device.contacts()[gate].gate;
    const physics::InsulatorParameters& ox = physics::silicon_dioxide_tunnelling_parameters;
    const double depth = physics::intrinsic_level_depth_eV(p.device.material(p.device.reference_node()), T);
    const double vac_s = depth - b.fields.potential_V[node];
    const double FG = -b.bias_V[gate];
    const double vac_g = FG + assemble::gate_work_function_eV(g, si, T);
    const physics::BarrierSegment seg{vac_s - ox.electron_affinity_eV,
                                      vac_g - ox.electron_affinity_eV, g.oxide_thickness_cm,
                                      ox.tunnel_mass_n};
    const double richardson = base::q_C *
                              physics::emission_velocity_cm_s(si, physics::Carrier::electron, T) *
                              physics::conduction_band_dos(si, T) / (T * T);
    const double low = std::max(b.bands.conduction_eV[node], FG);
    const double Fn = b.bands.electron_fermi_eV[node];
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

// GC8's comparison where the band shift is not zero (GD1): the gate at x_max over 100 nm of
// germanium, the reference node (node 0) silicon, so the gate node's E_c is the vacuum level less
// germanium's affinity, not silicon's.
TEST_CASE("Gate tunnelling: GC8's ECB on a germanium surface, off the reference material",
          "[gate_tunnelling]") {
    // GR5 as GC9: the near-singular capacitor (its inversion electrons almost float), here with
    // the gate at the far end of the ordering, trips the floating-region guard, so it is off.
    deck::Plan p = must(R"(
[device]
temperature = 300

[mesh]
type = tensor
x = linspace(0, 950nm, 96) + linspace(950nm, 998nm, 97) + linspace(998nm, 1um, 201)

[region silicon]
material = silicon

[region germanium]
material = germanium
box = 900nm 1um

[doping]
profile = acceptor constant 5e17

[contact gate]
kind = gate
patch = x_max
oxide_thickness = 1.5nm
electrode = n_poly
tunnelling = on

[contact substrate]
kind = ohmic
patch = x_min

[sweep]
contact gate = 0
)");
    p.options.linear.min_pivot_ratio = 0.0;
    const std::size_t top = p.device.mesh().node_count() - 1;
    REQUIRE(p.device.reference_node() == 0);
    REQUIRE(p.device.material(static_cast<mesh::NodeId>(top)).parameters().electron_affinity_eV !=
            p.device.material(0).parameters().electron_affinity_eV);
    const results::Sweep s = sweep(p, steps(0.0, 1.5, 0.25));
    for (const results::BiasPoint& b : s.points) {
        if (b.bias_V[gate] == 0.0) continue;
        const double J = ecb_of_state(p, b, Supply::both, top);
        const double one_way = std::max(std::abs(ecb_of_state(p, b, Supply::semiconductor, top)),
                                        std::abs(ecb_of_state(p, b, Supply::electrode, top)));
        CAPTURE(b.bias_V[gate], J, one_way, b.tunnelling_ecb[gate]);
        CHECK(b.tunnelling_ecb[gate] != 0.0);
        CHECK(std::abs(b.tunnelling_ecb[gate] - J) <= 1e-12 * one_way);
    }
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

namespace {

// GC11 in 3D: silicon from 1.5 nm on x and on y with 1.5 nm of SiO2 before it on both (an L around
// the corner; the oxide's corner block, 0 to 1.5 nm on both, belongs to no contact) and an n+ poly
// gate on each face; in 3D the same section along z over 1 um. a and b: the gates' tunnelling.
std::string corner_deck(int dimension, const char* a, const char* b) {
    const char* axis = "linspace(0, 1.45nm, 30) + linspace(1.495nm, 1.505nm, 2) + "
                       "linspace(1.505nm, 3.505nm, 21) + linspace(3.505nm, 51.505nm, 25)";
    const char* z_mesh = dimension == 3 ? "z = linspace(0, 1um, 3)\n" : "";
    const char* z = dimension == 3 ? "  0 1um" : "";
    char text[4096];
    std::snprintf(text, sizeof text, R"(
[device]
temperature = 300

[mesh]
type = tensor
x = %s
y = %s
%s
[region oxide]
material = sio2

[region silicon]
material = silicon
box = 1.5nm 1um  1.5nm 1um%s

[doping]
profile = acceptor constant 5e17

[contact gate_a]
kind = electrode
patch = x_min
box = 0 0  1.5nm 51.505nm%s
electrode = n_poly
tunnelling = %s

[contact gate_b]
kind = electrode
patch = y_min
box = 1.5nm 51.505nm  0 0%s
electrode = n_poly
tunnelling = %s

[contact substrate]
kind = ohmic
patch = x_max
box = 51.505nm 51.505nm  1.5nm 51.505nm%s

[sweep]
contact gate_a = 0
contact gate_b = 0
)",
                  axis, axis, z_mesh, z, z, a, z, b, z);
    return text;
}

}  // namespace

// At one fixed state, so the comparison is geometric only: the 2D section's equilibrium potential
// with both gates at 1.5 V (state_from_potential), copied along z in 3D. Each 2D device tunnels
// through one gate, so its corner node has one path; the 3D device's corner node has both, each
// with its own face's area, and each 3D gate carries 1 um times its 2D device's current.
TEST_CASE("Gate tunnelling: a 3D corner node takes one path per face (GC11)", "[gate_tunnelling]") {
    const deck::Plan a = must(corner_deck(2, "on", "off"));
    const deck::Plan b = must(corner_deck(2, "off", "on"));
    const deck::Plan c = must(corner_deck(3, "on", "on"));
    const auto eq = solve::solve_equilibrium(a.device);
    REQUIRE(eq.has_value());
    // Per contact, the tunnelling current at the potential psi_V [A in 3D, A/cm in 2D].
    const auto currents = [](const deck::Plan& p, std::vector<double> psi) {
        const auto scaling = *assemble::make_scaling(p.device);
        auto system = *assemble::DriftDiffusion::create(p.device, scaling, p.options.models);
        REQUIRE(system.set_bias(std::vector<double>{1.5, 1.5, 0.0}).has_value());
        for (double& v : psi) v /= scaling.V_T;
        const auto t = system.gate_tunnelling_currents(system.state_from_potential(psi));
        const double unit = scaling.J0 * std::pow(scaling.L_D, p.device.mesh().dimension() - 1);
        std::vector<double> I;
        for (std::size_t k = 0; k < t.ecb.size(); ++k) I.push_back((t.ecb[k] + t.hvb[k] + t.evb[k]) * unit);
        return I;
    };
    std::map<std::pair<double, double>, double> section;
    const auto points = a.device.mesh().points();
    for (std::size_t i = 0; i < points.size(); ++i) {
        section[{points[i][0], points[i][1]}] = eq->fields.potential_V[i];
    }
    std::vector<double> psi3;
    for (const mesh::Point& q : c.device.mesh().points()) psi3.push_back(section.at({q[0], q[1]}));
    const std::vector<double> Ia = currents(a, eq->fields.potential_V);
    const std::vector<double> Ib = currents(b, eq->fields.potential_V);
    const std::vector<double> Ic = currents(c, psi3);
    constexpr double width = 1e-4;  // [cm]
    std::printf("GC11 3D corner: gate_a %.6e A (2D x 1 um %.6e), gate_b %.6e A (%.6e)\n", Ic[0],
                width * Ia[0], Ic[1], width * Ib[1]);
    CHECK(Ia[0] > 0.0);
    CHECK(Ia[1] == 0.0);
    CHECK(Ib[0] == 0.0);
    CHECK(Ib[1] > 0.0);
    CHECK(relative(Ic[0], width * Ia[0]) <= 1e-9);
    CHECK(relative(Ic[1], width * Ib[1]) <= 1e-9);
    CHECK(relative(Ic[0] + Ic[1], width * (Ia[0] + Ib[1])) <= 1e-9);
}

// GD2: a path ends where it meets a semiconductor node. Along x, electrode | 1.5 nm SiO2 |
// silicon A | 1.5 nm SiO2 | silicon B (A contacted on its y_max face, so it does not float): A
// tunnels to the gate, B's path meets A, so B has none. At fixed states, moving B's potential
// leaves the gate's tunnelling current unchanged; moving A's does not.
TEST_CASE("Gate tunnelling: no path through a semiconductor (GD2)", "[gate_tunnelling]") {
    const deck::Plan p = must(R"(
[device]
temperature = 300

[mesh]
type = tensor
x = linspace(0, 20nm, 401)
y = linspace(0, 1nm, 3)

[region oxide]
material = sio2

[region silicon_a]
material = silicon
box = 1.5nm 3nm  0 1nm

[region silicon_b]
material = silicon
box = 4.5nm 20nm  0 1nm

[doping]
profile = acceptor constant 5e17

[contact gate]
kind = electrode
patch = x_min
electrode = n_poly
tunnelling = on

[contact substrate]
kind = ohmic
patch = x_max

[contact a]
kind = ohmic
patch = y_max
box = 1.5nm 3nm  1nm 1nm

[sweep]
contact gate = 0
)");
    const auto scaling = *assemble::make_scaling(p.device);
    auto system = *assemble::DriftDiffusion::create(p.device, scaling, p.options.models);
    REQUIRE(system.set_bias(std::vector<double>{1.0, 0.0, 0.0}).has_value());
    const auto points = p.device.mesh().points();
    const auto current = [&](double shift_a, double shift_b) {
        std::vector<double> psi(points.size(), 0.0);
        for (std::size_t i = 0; i < points.size(); ++i) {
            if (p.device.is_insulator(static_cast<mesh::NodeId>(i))) continue;
            psi[i] = points[i][0] < 4e-7 ? shift_a : shift_b;
        }
        const auto t = system.gate_tunnelling_currents(system.state_from_potential(psi));
        return t.ecb[gate] + t.hvb[gate] + t.evb[gate];
    };
    const double base = current(0.0, 0.0);
    CHECK(base != 0.0);
    CHECK(current(0.0, 3.0) == base);
    CHECK(current(3.0, 0.0) != base);
}

// GD5, GD6, GD10: each process tunnels with its own mass. A lumped gate's stack with unequal masses
// (deck keys), at one fixed state (the equilibrium potential under 3 V, where all three processes
// carry current); raising one mass changes its own process's current and leaves the others' exact.
TEST_CASE("Gate tunnelling: each process takes its own mass (GD5, GD6, GD10)", "[gate_tunnelling]") {
    const auto plan = [](double mn, double mp, double mevb) {
        char text[2048];
        std::snprintf(text, sizeof text, R"(
[device]
temperature = 300

[mesh]
type = tensor
x = linspace(0, 2nm, 201) + linspace(2nm, 50nm, 97) + linspace(50nm, 1um, 96)

[region silicon]
material = silicon

[doping]
profile = acceptor constant 5e17

[contact gate]
kind = gate
patch = x_min
oxide_thickness = 1.5nm
electrode = n_poly
tunnelling = on
oxide_affinity = 1.05
oxide_band_gap = 9.0
oxide_mass_n = %g
oxide_mass_p = %g
oxide_mass_evb = %g

[contact substrate]
kind = ohmic
patch = x_max

[sweep]
contact gate = 0
)",
                      mn, mp, mevb);
        return must(text);
    };
    const deck::Plan base = plan(0.40, 0.60, 0.50);
    const auto eq = solve::solve_equilibrium(base.device);
    REQUIRE(eq.has_value());
    const auto currents = [&](const deck::Plan& p) {
        const auto scaling = *assemble::make_scaling(p.device);
        auto system = *assemble::DriftDiffusion::create(p.device, scaling, p.options.models);
        REQUIRE(system.set_bias(std::vector<double>{3.0, 0.0}).has_value());
        std::vector<double> psi(eq->fields.potential_V);
        for (double& v : psi) v /= scaling.V_T;
        const auto t = system.gate_tunnelling_currents(system.state_from_potential(psi));
        return std::array<double, 3>{t.ecb[gate], t.hvb[gate], t.evb[gate]};
    };
    const std::array<double, 3> b = currents(base);
    CAPTURE(b[0], b[1], b[2]);
    for (const double v : b) REQUIRE(v != 0.0);
    const std::array<double, 3> raised[3] = {currents(plan(0.45, 0.60, 0.50)),
                                             currents(plan(0.40, 0.65, 0.50)),
                                             currents(plan(0.40, 0.60, 0.55))};
    for (std::size_t k = 0; k < 3; ++k) {      // the mass raised: ECB's, HVB's, EVB's
        for (std::size_t j = 0; j < 3; ++j) {  // the process
            CAPTURE(k, j, raised[k][j]);
            if (j == k) CHECK(std::abs(raised[k][j]) < std::abs(b[j]));
            else CHECK(raised[k][j] == b[j]);
        }
    }
}

// GD10: an insulator region's band-data keys, each into its own field (five distinct values).
TEST_CASE("Gate tunnelling: an insulator region's band data from the deck (GD10)", "[gate_tunnelling]") {
    const deck::Plan p = must(R"(
[device]
temperature = 300

[mesh]
type = tensor
x = linspace(0, 1.45nm, 30) + linspace(1.495nm, 1.505nm, 2) + linspace(1.505nm, 1um, 101)

[region hfo2]
material = insulator
eps_r = 22
electron_affinity_eV = 2.55
band_gap_eV = 5.7
tunnel_mass_n = 0.18
tunnel_mass_p = 0.21
tunnel_mass_evb = 0.24

[region silicon]
material = silicon
box = 1.5nm 2um

[doping]
profile = acceptor constant 5e17

[contact gate]
kind = electrode
patch = x_min
electrode = metal
work_function = 4.4
tunnelling = on

[contact substrate]
kind = ohmic
patch = x_max

[sweep]
contact gate = 0
)");
    const device::Region& r = p.device.regions()[0];
    REQUIRE(r.name == "hfo2");
    const physics::InsulatorParameters& b = std::get<physics::Insulator>(r.material).parameters();
    CHECK(b.eps_r == 22.0);
    CHECK(b.electron_affinity_eV == 2.55);
    CHECK(b.band_gap_eV == 5.7);
    CHECK(b.tunnel_mass_n == 0.18);
    CHECK(b.tunnel_mass_p == 0.21);
    CHECK(b.tunnel_mass_evb == 0.24);
}

// GD9: the weighted terminal current (Unit 32) has no tunnelling term, so it is refused.
TEST_CASE("Gate tunnelling: the weighted current is refused (GD9)", "[gate_tunnelling]") {
    deck::Plan p = must(deck_of({}));
    p.options.weighted_current = true;
    const std::vector<std::vector<double>> points{{0.1, 0.0}};
    const auto s = solve::sweep_bias(p.device, points, p.options);
    REQUIRE_FALSE(s.has_value());
    INFO(s.error().message);
    CHECK(s.error().code == base::ErrorCode::invalid_input);
    CHECK(s.error().message.find("gate tunnelling") != std::string::npos);
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


namespace {

// The FinFET deck of decks/ with its three gates tunnelling through the interfacial SiO2 (the
// built-in band data) and the HfO2 (GR10); with mkl, PARDISO on 6 threads.
std::string finfet_text(bool tunnelling, const char* mkl = nullptr) {
    std::ifstream in(std::filesystem::path(NITCAD_DECKS_DIR) / "finfet_3d.deck");
    std::stringstream ss;
    ss << in.rdbuf();
    std::string text = ss.str();
    const auto replace_all = [&](const std::string& from, const std::string& to) {
        int count = 0;
        for (auto at = text.find(from); at != std::string::npos; at = text.find(from, at + to.size())) {
            text.replace(at, from.size(), to);
            ++count;
        }
        return count;
    };
    REQUIRE(replace_all("eps_r = 22\n", "eps_r = 22\nelectron_affinity_eV = 2.55\nband_gap_eV = 5.7\n"
                                        "tunnel_mass_n = 0.18\ntunnel_mass_p = 0.18\n") == 1);
    if (tunnelling) REQUIRE(replace_all("work_function = 4.4\n", "work_function = 4.4\ntunnelling = on\n") == 3);
    if (mkl != nullptr) {
        REQUIRE(replace_all("# backend = pardiso\n", "backend = pardiso\nmkl = " + std::string(mkl) +
                                                         "\nthreads = 6\n") == 1);
    }
    return text;
}

deck::Plan finfet(bool tunnelling, const char* mkl = nullptr) {
    return must(finfet_text(tunnelling, mkl), NITCAD_DECKS_DIR);
}

}  // namespace

// GC14 (hidden; run pinned to the P cores as Unit 40's gates): assembly with all three processes
// at most 1.5 times without, at the equilibrium potential under the sweep's last bias.
TEST_CASE("Gate tunnelling: FinFET assembly cost (GC14)", "[.gate_tunnelling_finfet]") {
    double evaluate[2] = {}, residual[2] = {};
    for (const bool on : {false, true}) {
        const deck::Plan p = finfet(on);
        const auto eq = solve::solve_equilibrium(p.device);
        REQUIRE(eq.has_value());
        const auto scaling = *assemble::make_scaling(p.device);
        auto system = *assemble::DriftDiffusion::create(p.device, scaling, p.options.models);
        REQUIRE(system.set_bias(p.points.back()).has_value());
        std::vector<double> psi(eq->fields.potential_V);
        for (double& v : psi) v /= scaling.V_T;
        const AssemblyTiming t = time_assembly(on ? "FinFET, gate tunnelling" : "FinFET",
                                               system, system.state_from_potential(psi));
        evaluate[on ? 1 : 0] = t.evaluate_us;
        residual[on ? 1 : 0] = t.residual_us;
    }
    std::printf("GC14 evaluate ratio %.3f, residual ratio %.3f\n", evaluate[1] / evaluate[0],
                residual[1] / residual[0]);
    CHECK(evaluate[1] <= 1.5 * evaluate[0]);
}

// GC14 and GC16 (hidden; about 40 min with PARDISO, NITCAD_MKL the path of mkl_rt.3.dll): the
// sweep with tunnelling converges at every point that converges without it; the gate current
// against the drain current, reported.
TEST_CASE("Gate tunnelling: the FinFET sweep (GC14, GC16)", "[.gate_tunnelling_finfet_sweep]") {
    char* value = nullptr;
    std::size_t length = 0;
    const bool given = _dupenv_s(&value, &length, "NITCAD_MKL") == 0 && value != nullptr;
    const std::string mkl = given ? value : "";
    std::free(value);
    if (mkl.empty()) SKIP("NITCAD_MKL not set");
    if (auto loaded = linalg::load_mkl_runtime(mkl); !loaded) FAIL(loaded.error().message);
    results::Sweep runs[2];
    for (const bool on : {false, true}) {
        const deck::Plan p = finfet(on, mkl.c_str());
        auto s = solve::sweep_bias(p.device, p.points, p.options);
        if (!s) FAIL(s.error().message);
        results::Sweep& r = runs[on ? 1 : 0];
        r = std::move(*s);
        const std::string stopped = r.stopped ? ", stopped: " + r.stopped->message : "";
        std::printf("GC14 FinFET %s: %zu of %zu points%s\n", on ? "with tunnelling" : "without",
                    r.points.size(), p.points.size(), stopped.c_str());
    }
    CHECK(runs[1].points.size() >= runs[0].points.size());
    const deck::Plan p = finfet(true);
    const auto index = [&](const std::string& name) {
        std::size_t c = 0;
        while (c < p.device.contacts().size() && p.device.contacts()[c].name != name) ++c;
        INFO("contact " << name);
        REQUIRE(c < p.device.contacts().size());
        return c;
    };
    const std::size_t drain = index("drain");
    for (const results::BiasPoint& b : runs[1].points) {
        double I_G = 0.0;
        for (const char* g : {"gate_top", "gate_left", "gate_right"}) {
            I_G += b.terminal_current[index(g)];
        }
        std::printf("GC16 V_G %.2f V: I_D %.4e A, I_G %.4e A (%.2e of I_D)\n",
                    b.bias_V[index("gate_top")], b.terminal_current[drain], I_G,
                    I_G / b.terminal_current[drain]);
    }
}
