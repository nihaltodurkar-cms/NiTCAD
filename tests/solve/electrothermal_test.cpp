// Electrothermal steady solves (ARCHITECTURE.md section 11, Unit 23; DECISIONS.md T1-T14 and the
// owner-approved gates): the heat equation alone against its Kirchhoff-transform solution
// (isothermal and R_th ends); the steady energy balance (heat out of the sinks equals sum I V);
// no current or heat at a uniform temperature at equilibrium; a conductivity 1e8 times larger
// against the isothermal solve, and a uniform 350 K against the isothermal device at 350 K; the
// Seebeck open-circuit voltage against the integral of the thermopower; the run record; and the
// inputs a steady electrothermal solve refuses.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/statistics.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/thermal.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "NiTCAD/solve/trace.hpp"
#include "NiTCAD/solve/transient.hpp"
#include "legacy_graded_mesh.hpp"

using namespace NiTCAD;
using base::ErrorCode;

namespace {

std::vector<double> uniform(double a, double b, int nodes) {
    std::vector<double> v(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) v[static_cast<std::size_t>(i)] = a + (b - a) * i / (nodes - 1);
    return v;
}

struct Bar {
    int nodes = 41;
    double length_cm = 1e-4;
    double temperature_K = 300.0;
    // Doping per node from its position.
    double (*donors)(double x, double L) = [](double, double) { return 1e17; };
    double (*acceptors)(double x, double L) = [](double, double) { return 0.0; };
    physics::SemiconductorParameters material = physics::silicon_parameters;
    std::vector<device::ThermalContact> thermal = {
        {"left", "x_min", device::ThermalContactKind::isothermal, 300.0, 0.0},
        {"right", "x_max", device::ThermalContactKind::isothermal, 300.0, 0.0}};
};

// A 1D silicon bar with ohmic contacts at both ends.
device::Device bar(const Bar& b) {
    mesh::Mesh m = *mesh::make_tensor_grid(uniform(0.0, b.length_cm, b.nodes));
    const std::size_t n = m.node_count();
    std::vector<double> donors(n), acceptors(n);
    for (std::size_t i = 0; i < n; ++i) {
        const double x = m.points()[i][0];
        donors[i] = b.donors(x, b.length_cm);
        acceptors[i] = b.acceptors(x, b.length_cm);
    }
    auto left = m.find_boundary("x_min")->nodes;
    auto right = m.find_boundary("x_max")->nodes;
    auto d = device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = b.temperature_K,
         .regions = {{"silicon", *physics::Semiconductor::create(b.material)}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                      {"right", device::ContactKind::ohmic, std::move(right)}},
         .thermal_contacts = b.thermal});
    if (!d) FAIL(d.error().message);
    return std::move(*d);
}

double pn_acceptors(double x, double L) { return x < 0.5 * L ? 1e17 : 0.0; }
double pn_donors(double x, double L) { return x < 0.5 * L ? 0.0 : 1e17; }

solve::BiasOptions thermal_options() {
    solve::BiasOptions o;
    o.models.electrothermal = true;
    o.newton.tol_update = 1e-12;
    return o;
}

results::BiasPoint solved(const device::Device& d, std::vector<double> bias,
                          const solve::BiasOptions& o,
                          const results::NodeFields* initial = nullptr) {
    auto p = solve::solve_bias(d, bias, o, initial);
    if (!p) FAIL(p.error().message);
    return std::move(*p);
}

// sum_c V_c I_c of a point [W / cm^(3-D)].
double electrical_power(const results::BiasPoint& p) {
    double P = 0.0;
    for (std::size_t c = 0; c < p.bias_V.size(); ++c) P += p.bias_V[c] * p.terminal_current[c];
    return P;
}

}  // namespace

TEST_CASE("electrothermal solve: off, nothing thermal in the results") {
    const device::Device d = bar({});
    const auto p = solved(d, {0.1, 0.0}, {});
    REQUIRE(p.fields.temperature_K.empty());
    REQUIRE(p.thermal_bias_K.empty());
    REQUIRE(p.thermal_contact_heat.empty());
    const auto on = solved(d, {0.1, 0.0}, thermal_options());
    REQUIRE(on.fields.temperature_K.size() == d.mesh().node_count());
    REQUIRE(on.thermal_bias_K == std::vector<double>{300.0, 300.0});
    REQUIRE(on.thermal_contact_heat.size() == 2);
}

TEST_CASE("electrothermal solve: the heat equation alone against the Kirchhoff transform") {
    // An insulator bar with silicon's conductivity law kappa = 1.48 (T / 300)^-1.33 and two
    // silicon nodes at x_min with the only ohmic contact (no current flows): steady conduction
    // with no source. Then u(T) = integral of kappa is linear in x, which the Kirchhoff form of
    // the discrete conduction reproduces exactly at the nodes (1D, any spacing). Isothermal ends at
    // 300 and 400 K; then an R_th end to a 400 K ambient, where the flux q = (T_R - 400) / R_th
    // with u(T_R) - u(300) = -q L... solved here by bisection.
    physics::InsulatorParameters oxide = physics::silicon_dioxide_parameters;
    oxide.thermal = physics::silicon_parameters.thermal;
    const auto build = [&](device::ThermalContact right) {
        // Graded spacing: the Kirchhoff form is exact on any 1D mesh.
        std::vector<double> x;
        for (int i = 0; i <= 30; ++i) x.push_back(1e-4 * std::pow(i / 30.0, 1.4));
        mesh::Mesh m = *mesh::make_tensor_grid(x);
        const std::size_t n = m.node_count();
        std::vector<device::RegionId> region(n, 1);
        region[0] = region[1] = 0;
        std::vector<double> donors(n, 0.0);
        donors[0] = donors[1] = 1e16;
        auto left = m.find_boundary("x_min")->nodes;
        auto d = device::Device::create(
            {.mesh = std::move(m),
             .temperature_K = 300.0,
             .regions = {{"silicon", physics::silicon()},
                         {"oxide", *physics::Insulator::create(oxide)}},
             .node_region = std::move(region),
             .donors = donors,
             .acceptors = std::vector<double>(n, 0.0),
             .contacts = {{"left", device::ContactKind::ohmic, std::move(left)}},
             .thermal_contacts = {{"left", "x_min", device::ThermalContactKind::isothermal,
                                   300.0, 0.0},
                                  std::move(right)}});
        if (!d) FAIL(d.error().message);
        return std::move(*d);
    };
    // u(T) = integral from 300 of kappa: 1.48 * 300 / (1 - 1.33) ((T/300)^(1-1.33) - 1).
    const auto u = [](double T) {
        return 1.48 * 300.0 / (1.0 - 1.33) * (std::pow(T / 300.0, 1.0 - 1.33) - 1.0);
    };
    const auto inverse = [&](double target) {  // u is increasing
        double lo = 200.0, hi = 600.0;
        for (int k = 0; k < 200; ++k) {
            const double mid = 0.5 * (lo + hi);
            (u(mid) < target ? lo : hi) = mid;
        }
        return 0.5 * (lo + hi);
    };
    {
        const device::Device d =
            build({"right", "x_max", device::ThermalContactKind::isothermal, 400.0, 0.0});
        const auto p = solved(d, {0.0}, thermal_options());
        double worst = 0.0;
        for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
            const double s = d.mesh().points()[i][0] / 1e-4;
            worst = std::max(worst, std::abs(p.fields.temperature_K[i] - inverse(s * u(400.0))));
        }
        CAPTURE(worst);
        REQUIRE(worst <= 1e-9);  // K
        // Heat in at the hot end leaves at the cold one: kappa flux u(400) / L [W/cm^2].
        const double q = u(400.0) / 1e-4;
        REQUIRE(std::abs(p.thermal_contact_heat[0] - q) <= 1e-9 * q);
        REQUIRE(std::abs(p.thermal_contact_heat[1] + q) <= 1e-9 * q);
    }
    {
        // R_th = 2e-5 K cm^2/W to 400 K: the flux into the bar q = (400 - T_R) / R_th and
        // u(T_R) = q L, so T_R solves u(T_R) = (400 - T_R) L / R_th.
        const double R = 2e-5;
        const device::Device d =
            build({"right", "x_max", device::ThermalContactKind::resistance, 400.0, R});
        const auto p = solved(d, {0.0}, thermal_options());
        double lo = 300.0, hi = 400.0;
        for (int k = 0; k < 200; ++k) {
            const double T = 0.5 * (lo + hi);
            (u(T) < (400.0 - T) * 1e-4 / R ? lo : hi) = T;
        }
        const double TR = 0.5 * (lo + hi);
        REQUIRE(std::abs(p.fields.temperature_K.back() - TR) <= 1e-9);
        double worst = 0.0;
        for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
            const double s = d.mesh().points()[i][0] / 1e-4;
            worst = std::max(worst, std::abs(p.fields.temperature_K[i] - inverse(s * u(TR))));
        }
        CAPTURE(TR, worst);
        REQUIRE(worst <= 1e-9);
    }
}

TEST_CASE("electrothermal solve: no current and no heat at a uniform temperature") {
    for (const double T : {300.0, 350.0}) {
        Bar b;
        b.donors = pn_donors;
        b.acceptors = pn_acceptors;
        for (auto& t : b.thermal) t.temperature_K = T;
        const device::Device d = bar(b);
        const auto p = solved(d, {0.0, 0.0}, thermal_options());
        for (const double Ti : p.fields.temperature_K) REQUIRE(std::abs(Ti - T) <= 1e-9 * T);
        for (std::size_t c = 0; c < 2; ++c) {
            CAPTURE(T, p.terminal_current[c], p.terminal_current_resolution[c]);
            REQUIRE(std::abs(p.terminal_current[c]) <= p.terminal_current_resolution[c]);
        }
        // The heats are of the size of that current times the energies of the carriers.
        for (const double h : p.thermal_contact_heat) REQUIRE(std::abs(h) <= 1e-6);
    }
}

TEST_CASE("electrothermal solve: the steady energy balance") {
    // A forward-biased diode, isothermal at x_min and R_th at x_max: the heat out of the two
    // sinks equals sum_c V_c I_c (Joule, recombination, Peltier and Thomson heat together).
    // DECISIONS.md U8, two gates at every point, at the default Newton tolerance (1e-8) and at
    // 1e-12:
    // - absolute: |heat - sum V I| <= E_max R_I. The two sides are the same energy flows summed
    //   two ways (the identity holds at any state, "the energy identity of the heat rows"), so they
    //   differ by the rounding of those flows: currents resolved only to R_I
    //   (terminal_current_resolution, a bound for the current through any cut) times the largest
    //   energy a carrier carries or a contact delivers, E_max [V]. Measured 1e-4 to 5e-3 of it.
    // - relative: |heat - sum V I| <= 1e-9 sum V I where the heating is resolved (E_max R_I below
    //   1e-9 of the power: 0.7 and 0.8 V). Measured 3.4e-12 (0.7 V at tolerance 1e-8; 7.8e-13 at
    //   1e-12) and 1.3e-14.
    // At 0.3 and 0.5 V the imbalance is 2.4e-12 and 2.3e-12 W/cm^2 at both tolerances: a floor
    // of the flows' rounding, not of the Newton tolerance; the tolerance shows at 0.7 V (1.4e-10
    // at 1e-8, 3.2e-11 at 1e-12), well inside both gates.
    Bar b;
    b.donors = pn_donors;
    b.acceptors = pn_acceptors;
    b.thermal[1] = {"right", "x_max", device::ThermalContactKind::resistance, 300.0, 1e-3};
    const device::Device d = bar(b);
    for (const double tol : {1e-8, 1e-12}) {
        auto o = thermal_options();
        o.newton.tol_update = tol;
        std::vector<std::vector<double>> points;
        for (const double V : {0.3, 0.5, 0.7, 0.8}) points.push_back({V, 0.0});
        auto sweep = solve::sweep_bias(d, points, o);
        REQUIRE(sweep.has_value());
        REQUIRE_FALSE(sweep->stopped);
        int resolved = 0;
        for (const results::BiasPoint& p : sweep->points) {
            const double power = electrical_power(p);
            const double heat = p.thermal_contact_heat[0] + p.thermal_contact_heat[1];
            double E = std::abs(p.bias_V[0]);  // the contacts' and carriers' energies [V]
            for (std::size_t i = 0; i < p.bands.conduction_eV.size(); ++i) {
                E = std::max({E, std::abs(p.bands.conduction_eV[i]),
                              std::abs(p.bands.valence_eV[i])});
            }
            E += 3.0 * base::k_B_eV_per_K * *std::ranges::max_element(p.fields.temperature_K);
            const double absolute = E * p.terminal_current_resolution[0];
            const double error = std::abs(heat - power);
            CAPTURE(tol, p.bias_V[0], power, heat, error, absolute);
            REQUIRE(power > 0.0);
            REQUIRE(error <= absolute);
            if (absolute <= 1e-9 * power) {
                ++resolved;
                REQUIRE(error <= 1e-9 * power);
            }
        }
        REQUIRE(resolved == 2);
        // The last point heats the R_th end measurably (6.5 mK: the isothermal end, 1 um away,
        // takes most of the heat).
        REQUIRE(sweep->points.back().fields.temperature_K.back() > 300.005);
    }
}

TEST_CASE("electrothermal solve: consistency with the isothermal solve") {
    // kappa 1e8 times silicon's: the temperature stays within 1e-6 K of the sinks', and the
    // currents and densities are those of the isothermal device at the sinks' temperature: at
    // 300 K (the reference), and at 350 K with the reference still 300 K (so the temperature
    // paths, not the scaling, carry the 50 K).
    for (const double T : {300.0, 350.0}) {
        Bar b;
        b.donors = pn_donors;
        b.acceptors = pn_acceptors;
        b.material.thermal.conductivity_W_cmK *= 1e8;
        for (auto& t : b.thermal) t.temperature_K = T;
        const device::Device hot = bar(b);
        Bar iso = b;
        iso.temperature_K = T;
        iso.thermal.clear();
        const device::Device reference = bar(iso);
        for (const double V : {0.2, 0.5, 0.7}) {
            const auto p = solved(hot, {V, 0.0}, thermal_options());
            solve::BiasOptions plain;
            plain.newton.tol_update = 1e-12;
            const auto q = solved(reference, {V, 0.0}, plain);
            double rise = 0.0;
            for (const double Ti : p.fields.temperature_K) rise = std::max(rise, std::abs(Ti - T));
            const double dI = std::abs(p.terminal_current[0] - q.terminal_current[0]) /
                              std::abs(q.terminal_current[0]);
            double dn = 0.0;
            for (std::size_t i = 0; i < hot.mesh().node_count(); ++i) {
                dn = std::max({dn, std::abs(p.fields.n_cm3[i] / q.fields.n_cm3[i] - 1.0),
                               std::abs(p.fields.p_cm3[i] / q.fields.p_cm3[i] - 1.0)});
            }
            // The current to 1e-10 or the two runs' resolution (at 0.2 V the current is near its
            // rounding floor, and the two scalings round differently: 6.4e-8 at 350 K).
            const double resolution = (p.terminal_current_resolution[0] +
                                       q.terminal_current_resolution[0]) /
                                      std::abs(q.terminal_current[0]);
            CAPTURE(T, V, rise, dI, resolution, dn);
            // Measured at most elsewhere: rise 1.2e-11 K, dI 3.0e-12, dn 1.1e-13.
            REQUIRE(rise <= 1e-9);
            REQUIRE(dI <= std::max(1e-10, resolution));
            REQUIRE(dn <= 1e-10);
        }
    }
}

TEST_CASE("electrothermal solve: a p-polysilicon gate at 350 K") {
    // A p-type bar (1e17) with a p-polysilicon gate (10 nm oxide) at x_min and an ohmic contact at
    // x_max, kappa 1e8 times silicon's, both ends' sinks at 350 K: the gate's work function
    // chi + Eg(350 K) and the charges equal the isothermal device's at 350 K (reference 300 K).
    const auto build = [](double T, bool sinks) {
        physics::SemiconductorParameters si = physics::silicon_parameters;
        si.thermal.conductivity_W_cmK *= 1e8;
        mesh::Mesh m = *mesh::make_tensor_grid(uniform(0.0, 1e-4, 81));
        const std::size_t n = m.node_count();
        auto gate = m.find_boundary("x_min")->nodes;
        auto body = m.find_boundary("x_max")->nodes;
        std::vector<device::ThermalContact> thermal;
        if (sinks) {
            thermal = {{"gate", "x_min", device::ThermalContactKind::isothermal, 350.0, 0.0},
                       {"body", "x_max", device::ThermalContactKind::isothermal, 350.0, 0.0}};
        }
        auto d = device::Device::create(
            {.mesh = std::move(m),
             .temperature_K = T,
             .regions = {{"silicon", *physics::Semiconductor::create(si)}},
             .node_region = std::vector<device::RegionId>(n, 0),
             .donors = std::vector<double>(n, 0.0),
             .acceptors = std::vector<double>(n, 1e17),
             .contacts = {{"gate", device::ContactKind::gate, std::move(gate),
                           {.boundary = "x_min", .oxide_thickness_cm = 1e-6,
                            .electrode = device::GateElectrode::p_poly}},
                          {"body", device::ContactKind::ohmic, std::move(body)}},
             .thermal_contacts = std::move(thermal)});
        if (!d) FAIL(d.error().message);
        return std::move(*d);
    };
    const device::Device hot = build(300.0, true), iso = build(350.0, false);
    solve::BiasOptions plain;
    plain.newton.tol_update = 1e-12;
    for (const double V : {-1.0, 0.0, 0.4}) {
        const auto p = solved(hot, {V, 0.0}, thermal_options());
        const auto q = solved(iso, {V, 0.0}, plain);
        const double dQ =
            std::abs(p.gate_charge[0] - q.gate_charge[0]) / std::abs(q.gate_charge[0]);
        CAPTURE(V, p.gate_charge[0], q.gate_charge[0], dQ);
        REQUIRE(dQ <= 1e-9);
    }
}

TEST_CASE("electrothermal solve: the Seebeck voltage of an n-type bar") {
    // N_D = 1e17, sinks at 300 and 310 K. Open circuit, J_n = -q mu n (grad phi_n + P_n grad T)
    // = 0, so V_R - V_L = -integral of P_n dT, with P_n = -(k/q) ((5/2 + r) - ln(n / Nc(T))),
    // n = N_D, Nc = Nc300 (T/300)^1.5, r = -1/2 (holes negligible). The open-circuit voltage is
    // read where the current of a sweep of the right contact crosses zero (a resistor: linear).
    Bar b;
    b.nodes = 81;
    b.length_cm = 2e-4;
    b.thermal[1].temperature_K = 310.0;
    const device::Device d = bar(b);
    auto o = thermal_options();
    const std::vector<std::vector<double>> points{{0.0, 0.0}, {0.0, 0.01}};
    const auto sweep = solve::sweep_bias(d, points, o);
    REQUIRE(sweep.has_value());
    REQUIRE(sweep->points.size() == 2);
    const double I0 = sweep->points[0].terminal_current[1];
    const double I1 = sweep->points[1].terminal_current[1];
    const double Voc = -I0 * 0.01 / (I1 - I0);
    const double TL = 300.0, TR = 310.0;
    const auto F = [](double T) { return T * std::log(T / 300.0) - T; };  // of ln(T/300)
    const double kq = base::k_B_eV_per_K;  // k/q [V/K]
    const double expected =
        kq * (2.0 * (TR - TL) - (TR - TL) * std::log(1e17 / 2.86e19) + 1.5 * (F(TR) - F(TL)));
    CAPTURE(Voc, expected, Voc / expected - 1.0);
    REQUIRE(expected > 0.0);
    REQUIRE(std::abs(Voc / expected - 1.0) <= 1e-5);  // measured 4.4e-7 (81 nodes)
}

TEST_CASE("electrothermal solve: Peltier heat at a current-carrying contact") {
    // An n-type bar (1e17) between two isothermal 300 K sinks. Electrons crossing an ohmic
    // contact from the metal take the energy e_n - E_F = kT ((5/2 + r) - ln(n / Nc)) = -q Pi_n
    // from the lattice, so the heat leaving through the right sink is -Pi_n I_R plus the Joule
    // heat. Joule is even in the current: the odd part of the heat over the current is -Pi_n.
    const device::Device d = bar({});
    const auto plus = solved(d, {1e-3, 0.0}, thermal_options());
    const auto minus = solved(d, {-1e-3, 0.0}, thermal_options());
    const double dQ = plus.thermal_contact_heat[1] - minus.thermal_contact_heat[1];
    const double dI = plus.terminal_current[1] - minus.terminal_current[1];
    const double kT = base::k_B_eV_per_K * 300.0;
    const double expected = kT * (2.0 - std::log(1e17 / 2.86e19));
    CAPTURE(dQ / dI, expected);
    REQUIRE(std::abs(dQ / dI / expected - 1.0) <= 1e-7);  // measured 4.7e-10
    // The left sink sees the opposite.
    const double dQ_left = plus.thermal_contact_heat[0] - minus.thermal_contact_heat[0];
    const double dI_left = plus.terminal_current[0] - minus.terminal_current[0];
    REQUIRE(std::abs(dQ_left / dI_left / expected - 1.0) <= 1e-7);
}

TEST_CASE("electrothermal solve: a diode with R_th against the lumped model; the runaway fold") {
    // A 1 um pn diode with kappa 1e3 times silicon's (uniform T: the internal drop is some 1e-6 K;
    // 1e8 would leave the R_th rows 1e-13 of the conduction, singular to working precision) and
    // an R_th of 2 K cm^2/W to
    // 300 K at each end (R_eff = 1 K cm^2/W). By the energy balance every watt of sum I V leaves
    // through the two resistances, so T = 300 + R_eff V I and the current is the isothermal
    // device's at that T: the lumped electrothermal model is exact here. Each lumped point is
    // solved by bisection in T over isothermal solves. The heating raises the current, which
    // raises the heating: past the fold of the I-V curve (dV/dI = 0) the branch turns back, and
    // trace_bias follows it; the fold's bias is the largest V of the lumped V(T) along
    // T - 300 = R_eff V I_iso(V, T).
    const double R = 2.0, R_eff = 1.0;
    Bar b;
    b.donors = pn_donors;
    b.acceptors = pn_acceptors;
    b.material.thermal.conductivity_W_cmK *= 1e3;
    for (auto& t : b.thermal) t = {t.name, t.boundary, device::ThermalContactKind::resistance,
                                   300.0, R};
    const device::Device d = bar(b);
    const auto isothermal_current = [&](double V, double T) {
        Bar iso = b;
        iso.temperature_K = T;
        iso.thermal.clear();
        solve::BiasOptions plain;
        plain.newton.tol_update = 1e-12;
        return solved(bar(iso), {V, 0.0}, plain).terminal_current[0];
    };
    // Lumped T at bias V: the root of T - 300 - R_eff V I_iso(V, T) below the fold (bisection
    // from 300 K up to where the heating overtakes).
    const auto lumped_T = [&](double V) {
        double lo = 300.0, hi = 300.0;
        while (hi - 300.0 <= R_eff * V * isothermal_current(V, hi)) hi += 5.0;
        for (int k = 0; k < 50; ++k) {
            const double T = 0.5 * (lo + hi);
            (T - 300.0 <= R_eff * V * isothermal_current(V, T) ? lo : hi) = T;
        }
        return 0.5 * (lo + hi);
    };
    auto o = thermal_options();
    for (const double V : {0.55, 0.62}) {
        const auto p = solved(d, {V, 0.0}, o);
        const double T = lumped_T(V);
        const double I = isothermal_current(V, T);
        double spread = 0.0;
        for (const double Ti : p.fields.temperature_K) {
            spread = std::max(spread, std::abs(Ti - p.fields.temperature_K.front()));
        }
        CAPTURE(V, T, p.fields.temperature_K.front(), I, p.terminal_current[0], spread);
        REQUIRE(spread <= 1e-6);
        REQUIRE(std::abs(p.fields.temperature_K.front() - T) <= 1e-6 * (T - 300.0));
        REQUIRE(std::abs(p.terminal_current[0] / I - 1.0) <= 1e-6);
    }

    // The fold: the trace's largest bias against the lumped V(T)'s maximum.
    solve::TraceOptions trace;
    trace.steady = o;
    trace.contact = 0;
    trace.start_V = {0.5, 0.0};
    trace.end_V = 1.0;
    trace.step_V = 0.02;
    trace.max_step_V = 0.02;
    const auto run = solve::trace_bias(d, trace);
    REQUIRE(run.has_value());
    // V and T along the branch; it must turn back.
    std::size_t top = 0;
    for (std::size_t k = 0; k < run->points.size(); ++k) {
        if (run->points[k].bias_V[0] > run->points[top].bias_V[0]) top = k;
    }
    REQUIRE(top > 0);
    REQUIRE(top + 1 < run->points.size());
    REQUIRE(run->points.back().bias_V[0] < run->points[top].bias_V[0]);
    // The branch carries on past the fold until the temperature leaves silicon's range (near
    // 857 K, holes' mu_max below mu_min): the trace stops there (T10) with its points.
    REQUIRE(run->stopped);
    REQUIRE(run->stopped->code == ErrorCode::non_convergence);
    REQUIRE(run->stopped->message.find("outside its material's range") != std::string::npos);
    REQUIRE(run->stopped->context);
    REQUIRE(run->stopped->context->value > 850.0);  // the offending temperature [K]
    REQUIRE(run->points.back().fields.temperature_K.front() > 800.0);
    // The trace's fold: the vertex of the parabola V(T) through the top point and its neighbours.
    const auto at = [&](std::size_t k) {
        return std::pair{run->points[k].fields.temperature_K.front(), run->points[k].bias_V[0]};
    };
    const auto [T1, V1] = at(top - 1);
    const auto [T2, V2] = at(top);
    const auto [T3, V3] = at(top + 1);
    const double s12 = (V2 - V1) / (T2 - T1), s23 = (V3 - V2) / (T3 - T2);
    const double curvature = (s23 - s12) / (T3 - T1);  // V = V2 + s (T - T2) + c (T - T2)^2
    const double slope = s12 + curvature * (T2 - T1);   // at T2
    const double traced = V2 - slope * slope / (4.0 * curvature);
    // The lumped fold: the largest V on T - 300 = R_eff V I_iso(V, T), V(T) by bisection in V
    // (V I_iso rises with V), its maximum by golden section in T.
    const auto lumped_V = [&](double T) {
        double lo = 0.5, hi = 0.8;
        for (int k = 0; k < 45; ++k) {
            const double V = 0.5 * (lo + hi);
            (R_eff * V * isothermal_current(V, T) < T - 300.0 ? lo : hi) = V;
        }
        return 0.5 * (lo + hi);
    };
    double a = T2 - 3.0, c = T2 + 3.0;
    const double g = 0.5 * (std::sqrt(5.0) - 1.0);
    double x1 = c - g * (c - a), x2 = a + g * (c - a);
    double f1 = lumped_V(x1), f2 = lumped_V(x2);
    for (int k = 0; k < 30; ++k) {
        if (f1 > f2) {
            c = x2;
            x2 = x1;
            f2 = f1;
            x1 = c - g * (c - a);
            f1 = lumped_V(x1);
        } else {
            a = x1;
            x1 = x2;
            f1 = f2;
            x2 = a + g * (c - a);
            f2 = lumped_V(x2);
        }
    }
    const double lumped = std::max(f1, f2);
    CAPTURE(traced, lumped, T2, 0.5 * (a + c));
    REQUIRE(std::abs(traced - lumped) <= 1e-7);  // measured 4.6e-9 V
}

namespace {

// An n-type resistor (1e17) of length L between isothermal 300 K sinks, its mobility and
// conductivity independent of T when `constant` (the parabola's linear heat equation).
Bar resistor(int nodes, double L, bool constant) {
    Bar b;
    b.nodes = nodes;
    b.length_cm = L;
    if (constant) {
        b.material.electron_mobility.T_exponent = 0.0;
        b.material.hole_mobility.T_exponent = 0.0;
        b.material.thermal.conductivity_exponent = 0.0;
    }
    return b;
}

}  // namespace

TEST_CASE("electrothermal solve: legacy G-PARABOLA and G-BC") {
    // G-PARABOLA (legacy test_m19_thermal.py): a uniform heat source H0 with a constant kappa
    // between isothermal ends gives T = T0 + H0 / (2 kappa) x (L - x). Here the source is the
    // Joule heat of a 10 um resistor at 1 mV (mobility and kappa independent of T): uniform,
    // H0 = V I / L, its thermoelectric (Thomson) part some 1e-8 of it; the rise is about 1e-6 K.
    // The control volumes' second difference is exact for a parabola, so the nodes match it.
    // G-BC: with an R_th of 1e-3 K cm^2/W at each end instead the rod runs hotter, and as the
    // problem is linear by T_res - T_iso = R Q_L + R (Q_R - Q_L) x / L, Q the heat through each
    // end: half the Joule heat each, plus or minus each contact's Peltier heat Pi I (some 200
    // times the Joule heat at 1 mV, which isothermal ends absorb where it is released).
    const double L = 1e-3, V = 1e-3, kappa = 1.48;
    const device::Device iso = bar(resistor(101, L, true));
    const auto p = solved(iso, {V, 0.0}, thermal_options());
    const double P = V * p.terminal_current[0];
    const double H0 = P / L;
    double rise = 0.0, worst = 0.0;
    for (std::size_t i = 0; i < iso.mesh().node_count(); ++i) {
        const double x = iso.mesh().points()[i][0];
        const double T = 300.0 + H0 / (2.0 * kappa) * x * (L - x);
        rise = std::max(rise, T - 300.0);
        worst = std::max(worst, std::abs(p.fields.temperature_K[i] - T));
    }
    CAPTURE(rise, worst, worst / rise);
    REQUIRE(rise > 5e-7);
    REQUIRE(worst <= 1e-6 * rise);  // measured 1.0e-7
    Bar rb = resistor(101, L, true);
    const double R = 1e-3;
    for (auto& t : rb.thermal) t = {t.name, t.boundary, device::ThermalContactKind::resistance,
                                    300.0, R};
    const auto q = solved(bar(rb), {V, 0.0}, thermal_options());
    const double QL = q.thermal_contact_heat[0], QR = q.thermal_contact_heat[1];
    double shift = 0.0;
    for (std::size_t i = 0; i < iso.mesh().node_count(); ++i) {
        const double x = iso.mesh().points()[i][0];
        const double expected = R * QL + R * (QR - QL) * x / L;
        shift = std::max(shift, std::abs(q.fields.temperature_K[i] -
                                         p.fields.temperature_K[i] - expected));
    }
    // The R_th run's own power: its ends differ by 1.3 mK (the Peltier heats through R_th), whose
    // Seebeck voltage moves the current by 8.6e-4.
    const double Pq = V * q.terminal_current[0];
    CAPTURE(R * P / 2.0, R * QL, R * QR, QL + QR - Pq, shift);
    REQUIRE(*std::ranges::max_element(q.fields.temperature_K) >
            *std::ranges::max_element(p.fields.temperature_K));
    REQUIRE(std::abs(QL + QR - Pq) <= 1e-9 * Pq);
    REQUIRE(shift <= 1e-5 * R * std::max(std::abs(QL), std::abs(QR)));  // measured 3.1e-6
}

TEST_CASE("electrothermal solve: a self-heated resistor against a reference ODE") {
    // Legacy G-ROLLOFF's direction, quantitatively. In an n-type resistor (N_D = 1e17, 10 um,
    // isothermal 300 K ends) n = N_D throughout (its space charge is some 1e-8 of it here), so the
    // model reduces to J = q mu(T) N_D (E_c' + (k/q)(1 + r) T') and the heat to H = J e_n' with
    // e_n = E_c + 2 kT: (kappa(T) T')' = -(J^2 / (q mu(T) N_D) + (3/2) (k/q) J T'), T = 300 K at
    // both ends, and V = integral of J / (q mu(T) N_D). That ODE is solved here to 2e-9 (Kirchhoff
    // form on 20000 cells, Newton) and J found for the device's V by the secant method. At 4 V the
    // centre is about 18 K hotter, the mobility lower, and the current below the isothermal one.
    const double L = 1e-3, N = 1e17, V = 4.0;
    const Bar b = resistor(201, L, false);
    const device::Device d = bar(b);
    const physics::Semiconductor si = *physics::Semiconductor::create(b.material);
    const physics::ThermalParameters& th = b.material.thermal;
    const double q = base::q_C, kq = base::k_B_eV_per_K;
    const auto mu = [&](double T) {
        return physics::caughey_thomas_mobility(si, physics::Carrier::electron, N, T);
    };
    const auto u = [&](double T) {  // the Kirchhoff transform, integral of kappa from 300 K
        const double s = 1.0 - th.conductivity_exponent;
        return th.conductivity_W_cmK * 300.0 / s * (std::pow(T / 300.0, s) - 1.0);
    };
    const auto kappa = [&](double T) { return physics::thermal_conductivity(th, T).value; };
    constexpr int M = 20000;
    const double h = L / M;
    // T on the grid for current density J, by Newton on the tridiagonal system.
    const auto profile = [&](double J) {
        std::vector<double> T(M + 1, 300.0);
        for (int it = 0; it < 50; ++it) {
            std::vector<double> a(M + 1), bb(M + 1), c(M + 1), r(M + 1);
            double change = 0.0;
            for (int k = 1; k < M; ++k) {
                // (u_{k+1} - 2 u_k + u_{k-1}) / h^2 + J^2 / (q mu N) + 1.5 kq J T' = 0.
                const double Tp = T[k + 1], Tk = T[k], Tm = T[k - 1];
                const double src = J * J / (q * mu(Tk) * N);
                const double dmu = physics::caughey_thomas_mobility_slope(
                    si, physics::Carrier::electron, N, Tk);
                r[k] = (u(Tp) - 2.0 * u(Tk) + u(Tm)) / (h * h) + src +
                       1.5 * kq * J * (Tp - Tm) / (2.0 * h);
                a[k] = kappa(Tm) / (h * h) - 1.5 * kq * J / (2.0 * h);
                c[k] = kappa(Tp) / (h * h) + 1.5 * kq * J / (2.0 * h);
                bb[k] = -2.0 * kappa(Tk) / (h * h) - src * dmu / mu(Tk);
            }
            // Thomas algorithm on rows 1..M-1 (T_0, T_M fixed): J dT = -r.
            std::vector<double> cp(M + 1), dp(M + 1);
            for (int k = 1; k < M; ++k) {
                const double den = bb[k] - (k > 1 ? a[k] * cp[k - 1] : 0.0);
                cp[k] = c[k] / den;
                dp[k] = (-r[k] - (k > 1 ? a[k] * dp[k - 1] : 0.0)) / den;
            }
            for (int k = M - 1; k >= 1; --k) {
                const double dT = dp[k] - (k < M - 1 ? cp[k] * dp[k + 1] : 0.0);
                dp[k] = dT;
                T[k] += dT;
                change = std::max(change, std::abs(dT));
            }
            if (change < 1e-12) break;
        }
        return T;
    };
    const auto voltage = [&](double J, std::vector<double>* T_out) {
        const std::vector<double> T = profile(J);
        double v = 0.0;
        for (int k = 0; k < M; ++k) {
            v += 0.5 * h * (J / (q * mu(T[k]) * N) + J / (q * mu(T[k + 1]) * N));
        }
        if (T_out != nullptr) *T_out = T;
        return v;
    };
    double J0 = V * q * mu(300.0) * N / L, J1 = 0.9 * J0;
    double V0 = voltage(J0, nullptr), V1 = voltage(J1, nullptr);
    for (int it = 0; it < 30 && std::abs(V1 - V) > 1e-14 * V; ++it) {
        const double J2 = J1 + (V - V1) * (J1 - J0) / (V1 - V0);
        J0 = J1;
        V0 = V1;
        J1 = J2;
        V1 = voltage(J1, nullptr);
    }
    std::vector<double> T_ref;
    (void)voltage(J1, &T_ref);
    const auto p = solved(d, {V, 0.0}, thermal_options());
    solve::BiasOptions plain;
    plain.newton.tol_update = 1e-12;
    Bar iso = b;
    iso.thermal.clear();
    const double I_iso = solved(bar(iso), {V, 0.0}, plain).terminal_current[0];
    const double T_mid = p.fields.temperature_K[100];
    CAPTURE(p.terminal_current[0], J1, I_iso, T_mid, T_ref[M / 2]);
    REQUIRE(T_ref[M / 2] - 300.0 > 10.0);
    REQUIRE(p.terminal_current[0] < I_iso);  // the roll-off
    REQUIRE(std::abs(p.terminal_current[0] / J1 - 1.0) <= 1e-5);  // measured 1.7e-6
    REQUIRE(std::abs(T_mid - T_ref[M / 2]) <= 1e-5 * (T_ref[M / 2] - 300.0));  // 2.2e-6
}

TEST_CASE("electrothermal solve: the Seebeck voltage under Fermi-Dirac statistics") {
    // A degenerate n-type bar (1e20 cm^-3, about 1.5 kT above the band edge) between 300 and
    // 310 K: V_R - V_L = integral of (k/q) (2 F_1(eta) / F_0(eta) - eta) dT, eta from
    // N_D = Nc(T) F_1/2(eta) (the thermopower of r = -1/2, thermal.hpp).
    Bar b;
    b.nodes = 81;
    b.length_cm = 2e-4;
    b.donors = [](double, double) { return 1e20; };
    b.thermal[1].temperature_K = 310.0;
    const device::Device d = bar(b);
    auto o = thermal_options();
    o.models.fermi_dirac = true;
    const std::vector<std::vector<double>> points{{0.0, 0.0}, {0.0, 0.01}};
    const auto sweep = solve::sweep_bias(d, points, o);
    REQUIRE(sweep.has_value());
    REQUIRE(sweep->points.size() == 2);
    const double I0 = sweep->points[0].terminal_current[1];
    const double I1 = sweep->points[1].terminal_current[1];
    const double Voc = -I0 * 0.01 / (I1 - I0);
    const physics::Semiconductor si = physics::silicon();
    const auto integrand = [&](double T) {
        const double Nc = physics::conduction_band_dos(si, T);
        const double eta =
            std::log(1e20 / Nc) -
            physics::fermi_dirac_degeneracy(1.0, std::log(Nc), 1e20).log_gamma;
        const physics::FermiOrdersZeroOne F = physics::fermi_orders_zero_one(eta);
        return base::k_B_eV_per_K * (2.0 * F.one / F.zero - eta);
    };
    double expected = 0.0;  // Simpson, 200 intervals
    const int K = 200;
    for (int k = 0; k <= K; ++k) {
        const double w = k == 0 || k == K ? 1.0 : (k % 2 == 1 ? 4.0 : 2.0);
        expected += w * integrand(300.0 + 10.0 * k / K);
    }
    expected *= 10.0 / K / 3.0;
    CAPTURE(Voc, expected, Voc / expected - 1.0);
    REQUIRE(expected > 0.0);
    REQUIRE(std::abs(Voc / expected - 1.0) <= 1e-5);
}

TEST_CASE("electrothermal solve: Peltier heat at an n+ / n junction") {
    // A 100 um bar, 1e18 cm^-3 below 50 um and 1e16 above, kappa independent of T, isothermal
    // 300 K ends, +-1 mV. Pi = T P_n = -(kT/q) (2 - ln(n / Nc)) changes across the junction, which
    // releases (Pi_L - Pi_R) I for a current I from left to right; at the bar's middle half of it
    // reaches each sink (constant kappa). With the right contact's own -Pi_R I_R, the odd part of
    // the right sink's heat over its current is -(Pi_L + Pi_R) / 2. The junction heat is spread
    // over the n side's space charge (a Debye length of 40 nm), 4e-4 of the bar from its middle.
    Bar b = resistor(401, 1e-2, true);
    b.donors = [](double x, double L) { return x < 0.5 * L ? 1e18 : 1e16; };
    const device::Device d = bar(b);
    const auto plus = solved(d, {1e-3, 0.0}, thermal_options());
    const auto minus = solved(d, {-1e-3, 0.0}, thermal_options());
    const double dQ = plus.thermal_contact_heat[1] - minus.thermal_contact_heat[1];
    const double dI = plus.terminal_current[1] - minus.terminal_current[1];
    const double kT = base::k_B_eV_per_K * 300.0;
    const double Nc = physics::conduction_band_dos(physics::silicon(), 300.0);
    const double expected =
        kT * (2.0 - 0.5 * std::log(1e18 / Nc) - 0.5 * std::log(1e16 / Nc));
    CAPTURE(dQ / dI, expected, dQ / dI / expected - 1.0);
    REQUIRE(std::abs(dQ / dI / expected - 1.0) <= 2e-3);  // measured 7.2e-4
}

TEST_CASE("electrothermal solve: 1D extruded to 2D and 3D") {
    // The diode of the energy-balance gate (isothermal x_min, R_th x_max) on 1D, 2D (5 rows) and
    // 3D (4 x 4) tensor grids: per unit area the currents and heats, and every node's temperature
    // and potential, are the 1D ones.
    const auto build = [](int D) {
        const auto x = uniform(0.0, 1e-4, 41);
        const auto y = uniform(0.0, 0.4e-4, D == 2 ? 5 : 4);
        mesh::Mesh m = D == 1   ? *mesh::make_tensor_grid(x)
                       : D == 2 ? *mesh::make_tensor_grid(x, y)
                                : *mesh::make_tensor_grid(x, y, y);
        const std::size_t n = m.node_count();
        std::vector<double> donors(n), acceptors(n);
        for (std::size_t i = 0; i < n; ++i) {
            const double px = m.points()[i][0];
            donors[i] = px < 0.5e-4 ? 0.0 : 1e17;
            acceptors[i] = px < 0.5e-4 ? 1e17 : 0.0;
        }
        auto left = m.find_boundary("x_min")->nodes;
        auto right = m.find_boundary("x_max")->nodes;
        auto dev = device::Device::create(
            {.mesh = std::move(m),
             .temperature_K = 300.0,
             .regions = {{"silicon", physics::silicon()}},
             .node_region = std::vector<device::RegionId>(n, 0),
             .donors = std::move(donors),
             .acceptors = std::move(acceptors),
             .contacts = {{"anode", device::ContactKind::ohmic, std::move(left)},
                          {"cathode", device::ContactKind::ohmic, std::move(right)}},
             .thermal_contacts = {
                 {"left", "x_min", device::ThermalContactKind::isothermal, 300.0, 0.0},
                 {"right", "x_max", device::ThermalContactKind::resistance, 300.0, 1e-3}}});
        if (!dev) FAIL(dev.error().message);
        return std::move(*dev);
    };
    const device::Device d1 = build(1);
    const auto p1 = solved(d1, {0.75, 0.0}, thermal_options());
    REQUIRE(p1.fields.temperature_K.back() > 300.0005);  // 0.68 mK
    for (const int D : {2, 3}) {
        const device::Device d = build(D);
        const auto p = solved(d, {0.75, 0.0}, thermal_options());
        const double area = D == 2 ? 0.4e-4 : 0.4e-4 * 0.4e-4;
        double worst = 0.0;
        for (std::size_t c = 0; c < 2; ++c) {
            worst = std::max(
                worst, std::abs(p.terminal_current[c] / area / p1.terminal_current[c] - 1.0));
            worst = std::max(worst, std::abs(p.thermal_contact_heat[c] / area /
                                             p1.thermal_contact_heat[c] - 1.0));
        }
        double dT = 0.0, dpsi = 0.0;
        for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
            const std::size_t i1 = i % 41;  // x varies fastest
            dT = std::max(dT, std::abs(p.fields.temperature_K[i] - p1.fields.temperature_K[i1]));
            dpsi = std::max(dpsi, std::abs(p.fields.potential_V[i] - p1.fields.potential_V[i1]));
        }
        CAPTURE(D, worst, dT, dpsi);
        REQUIRE(worst <= 1e-9);
        REQUIRE(dT <= 1e-9);
        REQUIRE(dpsi <= 1e-9);
    }
}

TEST_CASE("electrothermal solve: MOSFET self-heating", "[.mosfet]") {
    // The legacy MOSFET of mosfet_test.cpp (Lg 600 nm, n+ source and drain, N_A 1e17, 5 nm oxide,
    // n+ polysilicon gate) on a coarser graded mesh (40 x 20), with an isothermal 300 K sink on the
    // body contact's bottom face, at Vg = 1.5 V, Vd = 1 V: the heat out of the sink equals
    // sum I V (the gate carries none); the device warms, most of all on the drain side of the
    // channel (where the field and the current density are largest); and the drain current is below
    // the isothermal one (the lattice mobility falls with T). About 18 s in Release: hidden tag
    // [.mosfet], its own ctest entry (CMakeLists.txt). Measured: balance 2.3e-13, 301.6 K at the
    // drain edge, I_d 0.75% below the isothermal 6.34 A/cm.
    constexpr double Lg = 6e-5, Lsd = 3e-5, depth = 2e-5;
    const auto build = [&](bool sink) {
        const double L = 2 * Lsd + Lg;
        const auto x = legacy_graded_mesh(L, {Lsd, Lsd + Lg}, L / 800.0, L / 40.0, 1.15);
        const auto y = legacy_graded_mesh(depth, 0.0, depth / 400.0, depth / 20.0, 1.15);
        mesh::Mesh m = *mesh::make_tensor_grid(x, y);
        const std::size_t n = m.node_count();
        const auto sd = [](double px, double py, double edge, bool source) {
            const double s = source ? edge - px : px - edge;
            return 1e19 * std::exp(-(py * py) / (2.0 * 5e-6 * 5e-6)) * 0.5 *
                   std::erfc(-s / (std::sqrt(2.0) * 1e-6));
        };
        std::vector<double> donors(n);
        for (std::size_t i = 0; i < n; ++i) {
            const auto& p = m.points()[i];
            donors[i] = sd(p[0], p[1], Lsd, true) + sd(p[0], p[1], Lsd + Lg, false);
        }
        std::vector<mesh::NodeId> source, drain, gate;
        for (const mesh::NodeId v : m.find_boundary("y_min")->nodes) {
            const double px = m.points()[static_cast<std::size_t>(v)][0];
            (px <= Lsd ? source : px >= Lsd + Lg ? drain : gate).push_back(v);
        }
        auto body = m.find_boundary("y_max")->nodes;
        device::Contact g{"gate", device::ContactKind::gate, std::move(gate)};
        g.gate = {.boundary = "y_min", .oxide_thickness_cm = 5e-7};
        std::vector<device::ThermalContact> thermal;
        if (sink) thermal = {{"heat sink", "y_max", device::ThermalContactKind::isothermal, 300.0}};
        auto dev = device::Device::create(
            {.mesh = std::move(m),
             .temperature_K = 300.0,
             .regions = {{"silicon", physics::silicon()}},
             .node_region = std::vector<device::RegionId>(n, 0),
             .donors = std::move(donors),
             .acceptors = std::vector<double>(n, 1e17),
             .contacts = {{"source", device::ContactKind::ohmic, std::move(source)},
                          {"drain", device::ContactKind::ohmic, std::move(drain)},
                          {"body", device::ContactKind::ohmic, std::move(body)},
                          std::move(g)},
             .thermal_contacts = std::move(thermal)});
        if (!dev) FAIL(dev.error().message);
        return std::move(*dev);
    };
    // Contacts: source, drain, body, gate. The gate up to 1.5 V, then the drain up to 1 V.
    std::vector<std::vector<double>> points;
    for (int k = 1; k <= 6; ++k) points.push_back({0.0, 0.0, 0.0, 0.25 * k});
    for (int k = 1; k <= 10; ++k) points.push_back({0.0, 0.1 * k, 0.0, 1.5});
    const device::Device hot = build(true);
    auto o = thermal_options();
    o.newton.tol_update = 1e-10;
    const auto t0 = std::chrono::steady_clock::now();
    const auto sweep = solve::sweep_bias(hot, points, o);
    const auto t1 = std::chrono::steady_clock::now();
    REQUIRE(sweep.has_value());
    if (sweep->stopped) FAIL(sweep->stopped->message);
    const results::BiasPoint& p = sweep->points.back();
    double power = 0.0;
    for (std::size_t c = 0; c < 4; ++c) power += p.bias_V[c] * p.terminal_current[c];
    const double heat = p.thermal_contact_heat[0];
    std::size_t hottest = 0;
    for (std::size_t i = 0; i < hot.mesh().node_count(); ++i) {
        if (p.fields.temperature_K[i] > p.fields.temperature_K[hottest]) hottest = i;
    }
    solve::BiasOptions plain = o;
    plain.models.electrothermal = false;
    const device::Device cold = build(false);
    const auto t2 = std::chrono::steady_clock::now();
    const auto iso = solve::sweep_bias(cold, points, plain);
    const auto t3 = std::chrono::steady_clock::now();
    const std::chrono::duration<double> thermal_s = t1 - t0, isothermal_s = t3 - t2;
    std::size_t thermal_it = 0, isothermal_it = 0;
    for (const auto& q : sweep->points) thermal_it += q.convergence.iterations.size();
    for (const auto& q : iso->points) isothermal_it += q.convergence.iterations.size();
    REQUIRE(iso.has_value());
    REQUIRE_FALSE(iso->stopped);
    const double Id = p.terminal_current[1], Id_iso = iso->points.back().terminal_current[1];
    const double x_hot = hot.mesh().points()[hottest][0];
    CAPTURE(power, heat, p.fields.temperature_K[hottest], x_hot, Id, Id_iso, thermal_s.count(),
            isothermal_s.count(), thermal_it, isothermal_it);
    REQUIRE(std::abs(heat - power) <= 1e-9 * power);
    REQUIRE(p.fields.temperature_K[hottest] > 300.05);
    REQUIRE(x_hot > Lsd + 0.5 * Lg);
    REQUIRE(Id < Id_iso);
}

TEST_CASE("electrothermal solve: open circuit in a temperature gradient, by refinement (U2)") {
    // An n-type bar (1e17, 1 um) with one ohmic contact (x_min) between sinks at 300 and 400 K:
    // no current can flow, so the steady state is the thermoelectric equilibrium. Every edge
    // carries no current (within the state's resolution) although T varies along it, and the
    // electron quasi-Fermi level rises by integral of P_n dT = -(k/q) integral of
    // (2 - ln(N_D / Nc(T))) dT across the bar (the Seebeck voltage) for n = N_D; the model's
    // space charge at the ends departs from that by some 4e-6 of it. The edge-mean T flux is second
    // order: the differences of successive refinements fall by about 4 per halving.
    const auto build = [](int nodes) {
        mesh::Mesh m = *mesh::make_tensor_grid(uniform(0.0, 1e-4, nodes));
        const auto n = static_cast<std::size_t>(nodes);
        auto left = m.find_boundary("x_min")->nodes;
        auto d = device::Device::create(
            {.mesh = std::move(m),
             .temperature_K = 300.0,
             .regions = {{"silicon", physics::silicon()}},
             .node_region = std::vector<device::RegionId>(n, 0),
             .donors = std::vector<double>(n, 1e17),
             .acceptors = std::vector<double>(n, 0.0),
             .contacts = {{"left", device::ContactKind::ohmic, std::move(left)}},
             .thermal_contacts = {
                 {"cold", "x_min", device::ThermalContactKind::isothermal, 300.0, 0.0},
                 {"hot", "x_max", device::ThermalContactKind::isothermal, 400.0, 0.0}}});
        if (!d) FAIL(d.error().message);
        return std::move(*d);
    };
    const double TL = 300.0, TR = 400.0;
    const auto F = [](double T) { return T * std::log(T / 300.0) - T; };
    const double Nc300 = physics::conduction_band_dos(physics::silicon(), 300.0);
    const double expected = -base::k_B_eV_per_K * (2.0 * (TR - TL) -
                                                   (TR - TL) * std::log(1e17 / Nc300) +
                                                   1.5 * (F(TR) - F(TL)));
    std::vector<double> rises;
    for (const int nodes : {21, 41, 81, 161, 321}) {
        const device::Device d = build(nodes);
        const auto p = solved(d, {0.0}, thermal_options());
        double current = 0.0;
        for (std::size_t k = 0; k < p.edge_current_n.size(); ++k) {
            current = std::max({current, std::abs(p.edge_current_n[k]),
                                std::abs(p.edge_current_p[k])});
        }
        const double rise = p.bands.electron_fermi_eV.back() - p.bands.electron_fermi_eV.front();
        CAPTURE(nodes, current, p.terminal_current_resolution[0], rise, expected);
        REQUIRE(current <= p.terminal_current_resolution[0]);
        REQUIRE(std::abs(p.fields.temperature_K.back() - TR) <= 1e-9);
        rises.push_back(rise);
    }
    for (std::size_t k = 2; k < rises.size(); ++k) {
        const double before = std::abs(rises[k - 1] - rises[k - 2]);
        const double after = std::abs(rises[k] - rises[k - 1]);
        const double order = std::log2(before / after);
        CAPTURE(k, before, after, order);
        REQUIRE(order >= 1.8);
    }
    CAPTURE(rises.back(), expected);
    REQUIRE(std::abs(rises.back() / expected - 1.0) <= 1e-5);
}

TEST_CASE("electrothermal solve: a two-layer conductor, by refinement (U4)") {
    // Steady conduction through two insulator layers in series, kappa_A = 1.48 (T/300)^-1.33
    // (silicon's law; the two silicon nodes at x_min carrying the only ohmic contact share it) and
    // kappa_B = 0.3 (T/300)^0.5 W/(cm K), between isothermal 300 K (x_min) and 500 K (x_max), the
    // interface at 0.4 L plus half a spacing (an edge's midpoint, where the model puts it). Exact:
    // the heat flux q is constant, U_A(T_i) - U_A(300) = q x_i and U_B(500) - U_B(T_i) =
    // q (L - x_i) (U the Kirchhoff transforms), and inside each layer U is linear in x. Within a
    // layer the discrete conduction is exact; the interface edge's harmonic mean of the ends'
    // kappa at their own T is not: the error falls with refinement (U4's gate).
    const double L = 1e-3;
    physics::ThermalParameters A = physics::silicon_parameters.thermal;
    physics::ThermalParameters B = A;
    B.conductivity_W_cmK = 0.3;
    B.conductivity_exponent = -0.5;
    const auto U = [](const physics::ThermalParameters& t, double T) {
        const double s = 1.0 - t.conductivity_exponent;
        return t.conductivity_W_cmK * 300.0 / s * (std::pow(T / 300.0, s) - 1.0);
    };
    const auto inverse = [&](const physics::ThermalParameters& t, double target) {
        double lo = 1.0, hi = 2000.0;
        for (int k = 0; k < 200; ++k) {
            const double mid = 0.5 * (lo + hi);
            (U(t, mid) < target ? lo : hi) = mid;
        }
        return 0.5 * (lo + hi);
    };
    std::vector<double> errors;
    for (const int cells : {10, 20, 40, 80, 160}) {
        const double h = L / cells;
        const double xi = 0.4 * L + 0.5 * h;
        mesh::Mesh m = *mesh::make_tensor_grid(uniform(0.0, L, cells + 1));
        const std::size_t n = m.node_count();
        std::vector<device::RegionId> region(n, 1);
        region[0] = region[1] = 0;
        for (std::size_t i = 2; i < n; ++i) {
            if (m.points()[i][0] > xi) region[i] = 2;
        }
        std::vector<double> donors(n, 0.0);
        donors[0] = donors[1] = 1e16;
        physics::SemiconductorParameters si = physics::silicon_parameters;
        si.thermal = A;
        physics::InsulatorParameters layer_a = physics::silicon_dioxide_parameters;
        layer_a.thermal = A;
        physics::InsulatorParameters layer_b = physics::silicon_dioxide_parameters;
        layer_b.thermal = B;
        auto left = m.find_boundary("x_min")->nodes;
        auto d = device::Device::create(
            {.mesh = std::move(m),
             .temperature_K = 300.0,
             .regions = {{"silicon", *physics::Semiconductor::create(si)},
                         {"A", *physics::Insulator::create(layer_a)},
                         {"B", *physics::Insulator::create(layer_b)}},
             .node_region = std::move(region),
             .donors = std::move(donors),
             .acceptors = std::vector<double>(n, 0.0),
             .contacts = {{"left", device::ContactKind::ohmic, std::move(left)}},
             .thermal_contacts = {
                 {"cold", "x_min", device::ThermalContactKind::isothermal, 300.0, 0.0},
                 {"hot", "x_max", device::ThermalContactKind::isothermal, 500.0, 0.0}}});
        if (!d) FAIL(d.error().message);
        const auto p = solved(*d, {0.0}, thermal_options());
        // The exact interface temperature by bisection on the flux balance.
        double lo = 300.0, hi = 500.0;
        for (int k = 0; k < 200; ++k) {
            const double Ti = 0.5 * (lo + hi);
            const double qa = (U(A, Ti) - U(A, 300.0)) / xi;
            const double qb = (U(B, 500.0) - U(B, Ti)) / (L - xi);
            (qa < qb ? lo : hi) = Ti;
        }
        const double Ti = 0.5 * (lo + hi);
        const double q = (U(A, Ti) - U(A, 300.0)) / xi;
        double worst = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            const double x = d->mesh().points()[i][0];
            const double T = x < xi ? inverse(A, U(A, 300.0) + q * x)
                                    : inverse(B, U(B, 500.0) - q * (L - x));
            worst = std::max(worst, std::abs(p.fields.temperature_K[i] - T));
        }
        const double dq = std::abs(p.thermal_contact_heat[1] + q) / q;  // heat in at the hot end
        CAPTURE(cells, Ti, q, worst, dq);
        errors.push_back(worst);
    }
    for (std::size_t k = 1; k < errors.size(); ++k) {
        const double order = std::log2(errors[k - 1] / errors[k]);
        CAPTURE(k, errors[k - 1], errors[k], order);
        REQUIRE(order >= 1.8);  // measured 1.96 to 1.99
    }
    REQUIRE(errors.back() <= 1e-3);  // K, at 160 cells (measured 7.5e-4)
}

namespace {

// A 1D chain of `nodes` silicon nodes (n-type 1e17, length 1 um) at height y in a 2D mesh's node
// and edge lists, `first` its first node id, w the depth (the coupling and boundary areas).
void chain(std::vector<mesh::Point>& points, std::vector<double>& volumes,
           std::vector<mesh::Edge>& edges, int nodes, double y, double w) {
    const auto first = static_cast<mesh::NodeId>(points.size());
    const double h = 1e-4 / (nodes - 1);
    for (int k = 0; k < nodes; ++k) {
        points.push_back({h * k, y, 0.0});
        volumes.push_back((k == 0 || k == nodes - 1 ? 0.5 : 1.0) * h * w);
    }
    for (int k = 0; k + 1 < nodes; ++k) {
        edges.push_back({static_cast<mesh::NodeId>(first + k),
                         static_cast<mesh::NodeId>(first + k + 1), h, w});
    }
}

// One or two unconnected resistor chains (a, and b when `two`), contacts at both ends of each,
// an isothermal 300 K sink at chain a's left end only.
device::Device chains(bool two) {
    constexpr int N = 21;
    constexpr double w = 1e-4;
    std::vector<mesh::Point> points;
    std::vector<double> volumes;
    std::vector<mesh::Edge> edges;
    chain(points, volumes, edges, N, 0.0, w);
    if (two) chain(points, volumes, edges, N, 2e-4, w);
    std::vector<mesh::BoundaryPatch> patches{{"a_left", {0}, {w}}, {"a_right", {N - 1}, {w}}};
    if (two) {
        patches.push_back({"b_left", {N}, {w}});
        patches.push_back({"b_right", {2 * N - 1}, {w}});
    }
    const std::size_t n = points.size();
    auto m = mesh::Mesh::from_parts(2, std::move(points), std::move(volumes), std::move(edges),
                                    std::move(patches));
    if (!m) FAIL(m.error().message);
    std::vector<device::Contact> contacts{{"a1", device::ContactKind::ohmic, {0}},
                                          {"a2", device::ContactKind::ohmic, {N - 1}}};
    if (two) {
        contacts.push_back({"b1", device::ContactKind::ohmic, {N}});
        contacts.push_back({"b2", device::ContactKind::ohmic, {2 * N - 1}});
    }
    auto d = device::Device::create(
        {.mesh = std::move(*m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::vector<double>(n, 1e17),
         .acceptors = std::vector<double>(n, 0.0),
         .contacts = std::move(contacts),
         .thermal_contacts = {
             {"sink", "a_left", device::ThermalContactKind::isothermal, 300.0, 0.0}}});
    if (!d) FAIL(d.error().message);
    return std::move(*d);
}

}  // namespace

TEST_CASE("electrothermal solve: the transient's starting temperature (U5)") {
    // U5: each connected thermal domain starts on its own terms. (1) A region without a sink of
    // its own, joined to one: a 2D resistor whose silicon carries an isothermal sink on its
    // bottom face and an oxide layer above it with none; the run starts from the electrothermal
    // steady state everywhere, the oxide warmed through the silicon (not held at T0). (2) Two
    // unconnected resistors, one with a sink and one with none: the first starts from its own
    // steady state (that of the same chain alone), the second is held at T0 = 300 K at t = 0 and
    // its heat content grows adiabatically afterwards (its contacts' Peltier heat cools one end
    // by 0.08 K and heats the other more in the first ns); a steady sweep of that device is
    // refused (T7).
    solve::TransientOptions o;
    o.steady = thermal_options();
    o.t_end_s = 1e-9;
    o.dt_initial_s = o.dt_max_s = 2.5e-10;
    o.adaptive = false;
    {
        const auto x = uniform(0.0, 1e-4, 11), y = uniform(0.0, 0.5e-4, 6);
        mesh::Mesh m = *mesh::make_tensor_grid(x, y);
        const std::size_t n = m.node_count();
        std::vector<device::RegionId> region(n, 0);
        std::vector<double> donors(n, 1e17);
        for (std::size_t i = 0; i < n; ++i) {
            if (m.points()[i][1] >= 0.3e-4 - 1e-12) {
                region[i] = 1;
                donors[i] = 0.0;
            }
        }
        const auto silicon = [&](const char* patch) {
            std::vector<mesh::NodeId> v;
            for (const mesh::NodeId k : m.find_boundary(patch)->nodes) {
                if (region[static_cast<std::size_t>(k)] == 0) v.push_back(k);
            }
            return v;
        };
        // Before the mesh moves into the description (its members are initialized in order).
        auto left = silicon("x_min"), right = silicon("x_max");
        auto d = device::Device::create(
            {.mesh = std::move(m),
             .temperature_K = 300.0,
             .regions = {{"silicon", physics::silicon()}, {"oxide", physics::silicon_dioxide()}},
             .node_region = std::move(region),
             .donors = std::move(donors),
             .acceptors = std::vector<double>(n, 0.0),
             .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                          {"right", device::ContactKind::ohmic, std::move(right)}},
             .thermal_contacts = {
                 {"sink", "y_min", device::ThermalContactKind::isothermal, 300.0, 0.0}}});
        if (!d) FAIL(d.error().message);
        const auto steady = solved(*d, {0.5, 0.0}, thermal_options());
        const std::vector<solve::Waveform> w{solve::Waveform::constant(0.5),
                                             solve::Waveform::constant(0.0)};
        const auto run = solve::solve_transient(*d, w, o);
        REQUIRE(run.has_value());
        REQUIRE_FALSE(run->stopped);
        const auto& T0 = run->snapshots.front().fields.temperature_K;
        double worst = 0.0, oxide = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            worst = std::max(worst, std::abs(T0[i] - steady.fields.temperature_K[i]));
            if (d->is_insulator(static_cast<mesh::NodeId>(i))) oxide = std::max(oxide, T0[i]);
        }
        CAPTURE(worst, oxide);
        REQUIRE(worst <= 1e-9);
        REQUIRE(oxide > 300.0 + 1e-6);
    }
    {
        const device::Device both = chains(true);
        const std::vector<std::vector<double>> point{{0.5, 0.0, 0.5, 0.0}};
        const auto refused = solve::sweep_bias(both, point, thermal_options());
        REQUIRE_FALSE(refused.has_value());
        REQUIRE(refused.error().message.find("no thermal contact") != std::string::npos);
        const auto alone = solved(chains(false), {0.5, 0.0}, thermal_options());
        const std::vector<solve::Waveform> w{
            solve::Waveform::constant(0.5), solve::Waveform::constant(0.0),
            solve::Waveform::constant(0.5), solve::Waveform::constant(0.0)};
        const auto run = solve::solve_transient(both, w, o);
        REQUIRE(run.has_value());
        REQUIRE_FALSE(run->stopped);
        const auto& start = run->snapshots.front().fields.temperature_K;
        const auto& end = run->snapshots.back().fields.temperature_K;
        double worst = 0.0, held = 0.0, heated = 0.0;  // heated: sum of V_i (T_i - 300)
        const auto volumes = both.mesh().volumes();
        for (std::size_t i = 0; i < 21; ++i) {
            worst = std::max(worst, std::abs(start[i] - alone.fields.temperature_K[i]));
            held = std::max(held, std::abs(start[21 + i] - 300.0));
            heated += volumes[21 + i] * (end[21 + i] - 300.0);
        }
        CAPTURE(worst, held, heated, *std::ranges::max_element(alone.fields.temperature_K));
        REQUIRE(*std::ranges::max_element(alone.fields.temperature_K) > 300.0 + 1e-6);
        REQUIRE(worst <= 1e-9);
        REQUIRE(held == 0.0);
        REQUIRE(heated > 0.0);
    }
}

TEST_CASE("electrothermal solve: the run record carries the thermal data") {
    const device::Device d = bar({});
    const std::vector<std::vector<double>> points{{0.1, 0.0}};
    solve::BiasOptions off;
    const auto on = thermal_options();
    const auto r_off = solve::make_run_record(d, off, points);
    const auto r_on = solve::make_run_record(d, on, points);
    REQUIRE(r_off.input_identity != r_on.input_identity);
    // Off, the thermal contacts and temperatures are not read, so they do not change the identity.
    Bar other;
    other.thermal[0].temperature_K = 320.0;
    REQUIRE(solve::make_run_record(bar(other), off, points).input_identity == r_off.input_identity);
    REQUIRE(solve::make_run_record(bar(other), on, points).input_identity != r_on.input_identity);
    auto swept = on;
    swept.thermal_bias_K = {{300.0, 305.0}};
    REQUIRE(solve::make_run_record(d, swept, points).input_identity != r_on.input_identity);
}

namespace {

// A rod of length 1e-3 cm for the heat equation alone: two silicon nodes at x_min (the only ohmic
// contact; no current flows) and an insulator, both with kappa = 1.48 W/(cm K) independent of T
// and silicon's rho c, so T_t = D T_xx with D = kappa / rho c; an isothermal sink at x_min at
// 310 K, x_max adiabatic.
device::Device rod(int nodes) {
    physics::ThermalParameters thermal = physics::silicon_parameters.thermal;
    thermal.conductivity_exponent = 0.0;
    physics::SemiconductorParameters si = physics::silicon_parameters;
    si.thermal = thermal;
    physics::InsulatorParameters oxide = physics::silicon_dioxide_parameters;
    oxide.thermal = thermal;
    mesh::Mesh m = *mesh::make_tensor_grid(uniform(0.0, 1e-3, nodes));
    const auto n = static_cast<std::size_t>(nodes);
    std::vector<device::RegionId> region(n, 1);
    region[0] = region[1] = 0;
    std::vector<double> donors(n, 0.0);
    donors[0] = donors[1] = 1e16;
    auto left = m.find_boundary("x_min")->nodes;
    auto d = device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", *physics::Semiconductor::create(si)},
                     {"oxide", *physics::Insulator::create(oxide)}},
         .node_region = std::move(region),
         .donors = std::move(donors),
         .acceptors = std::vector<double>(n, 0.0),
         .contacts = {{"left", device::ContactKind::ohmic, std::move(left)}},
         .thermal_contacts = {{"hot", "x_min", device::ThermalContactKind::isothermal, 310.0,
                               0.0}}});
    if (!d) FAIL(d.error().message);
    return std::move(*d);
}

// The rod from 300 K (held at t = 0, the sink at 310 K), fixed steps of t_end / steps.
// With rtol > 0, adaptive steps from t_end / steps instead.
results::Transient heated_rod(const device::Device& d, solve::Integrator integrator, int steps,
                              double t_end, double rtol = 0.0) {
    auto start = solved(d, {0.0}, thermal_options());
    start.fields.temperature_K.assign(d.mesh().node_count(), 300.0);
    solve::TransientOptions o;
    o.steady = thermal_options();
    o.integrator = integrator;
    o.t_end_s = t_end;
    o.dt_initial_s = o.dt_max_s = t_end / steps;
    o.adaptive = rtol > 0.0;
    if (o.adaptive) {
        o.rtol = rtol;
        o.dt_max_s = t_end;
    }
    const std::vector<solve::Waveform> w{solve::Waveform::constant(0.0)};
    auto run = solve::solve_transient(d, w, o, &start.fields);
    if (!run) FAIL(run.error().message);
    if (run->stopped) FAIL(run->stopped->message);
    return std::move(*run);
}

}  // namespace

TEST_CASE("electrothermal solve: the step-heated rod, orders 1 and 2") {
    // Against the series T = T1 + (T0 - T1) sum_k 4 / ((2k+1) pi) sin(l_k x) exp(-l_k^2 D t),
    // l_k = (2k+1) pi / (2L) (x from the sink, the far end adiabatic), at t = 2e-7 s (0.18 of
    // L^2 / D); the time error's order against a run of 1280 steps on the same mesh (40 and 80
    // steps: at 20 backward Euler is not yet asymptotic after the sink's jump at t = 0): halving
    // the step halves backward Euler's error and quarters BDF2's.
    const double t_end = 2e-7, L = 1e-3;
    const double D = 1.48 / physics::silicon_parameters.thermal.heat_capacity_J_cm3K;
    const auto series = [&](double x, double t) {
        double sum = 0.0;
        for (int k = 0; k < 200; ++k) {
            const double l = (2 * k + 1) * std::numbers::pi / (2.0 * L);
            sum += 4.0 / ((2 * k + 1) * std::numbers::pi) * std::sin(l * x) *
                   std::exp(-l * l * D * t);
        }
        return 310.0 + (300.0 - 310.0) * sum;
    };
    {
        const device::Device d = rod(101);
        const auto run = heated_rod(d, solve::Integrator::bdf2, 400, t_end);
        const auto& T = run.snapshots.back().fields.temperature_K;
        double worst = 0.0;
        for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
            worst = std::max(worst, std::abs(T[i] - series(d.mesh().points()[i][0], t_end)));
        }
        CAPTURE(worst);
        REQUIRE(worst <= 5e-4);  // K, the spatial error of 101 nodes (measured 1.2e-4)
    }
    const device::Device d = rod(41);
    const auto end = [&](solve::Integrator integrator, int steps) {
        const auto run = heated_rod(d, integrator, steps, t_end);
        return run.snapshots.back().fields.temperature_K.back();
    };
    const double reference = end(solve::Integrator::bdf2, 1280);
    for (const auto& [integrator, order] :
         {std::pair{solve::Integrator::backward_euler, 1}, std::pair{solve::Integrator::bdf2, 2}}) {
        const double e1 = std::abs(end(integrator, 40) - reference);
        const double e2 = std::abs(end(integrator, 80) - reference);
        CAPTURE(order, e1, e2, e1 / e2);
        REQUIRE(std::abs(std::log2(e1 / e2) - order) <= 0.15);
    }
    // Adaptive steps (rtol 1e-6 of T0 on the rise, 0.3 mK a step, from t_end / 1000): only the
    // temperature changes, so the error estimate must measure it (T12) for the run to stay
    // accurate. Measured: 124 steps, 1.4e-3 K.
    const auto adaptive = heated_rod(d, solve::Integrator::bdf2, 1000, t_end, 1e-6);
    const double e = std::abs(adaptive.snapshots.back().fields.temperature_K.back() - reference);
    CAPTURE(adaptive.points.size(), adaptive.rejected_steps, e);
    REQUIRE(e <= 3e-3);
}

TEST_CASE("electrothermal solve: the transient energy balance") {
    // At every step, sum_c V_c I_c (conduction) = the heat out of the sinks + the lattice's
    // storage rho c dT/dt + the carriers' stored energy e_n dn/dt - e_p dp/dt (e_n = E_c + 2 kT,
    // e_p = E_v - 2 kT for r = -1/2), each rate the step's BDF difference. An n-type bar at 0.5 V,
    // from 300 K held at t = 0: adiabatic (no sink; T7 allows it in time), and with an isothermal
    // and an R_th sink; backward Euler and BDF2, fixed steps.
    for (const bool sinks : {false, true}) {
        Bar b;
        if (!sinks) {
            b.thermal.clear();
        } else {
            b.thermal[1] = {"right", "x_max", device::ThermalContactKind::resistance, 300.0,
                            1e-4};
        }
        const device::Device d = bar(b);
        auto start = solved(bar({}), {0.5, 0.0}, thermal_options());
        start.fields.temperature_K.assign(d.mesh().node_count(), 300.0);
        for (const auto integrator :
             {solve::Integrator::backward_euler, solve::Integrator::bdf2}) {
            solve::TransientOptions o;
            o.steady = thermal_options();
            o.integrator = integrator;
            o.t_end_s = 1e-8;
            o.dt_initial_s = o.dt_max_s = 1e-9;
            o.adaptive = false;
            o.fields_every_step = true;
            const std::vector<solve::Waveform> w{solve::Waveform::constant(0.5),
                                                 solve::Waveform::constant(0.0)};
            const auto run = solve::solve_transient(d, w, o, &start.fields);
            REQUIRE(run.has_value());
            REQUIRE_FALSE(run->stopped);
            REQUIRE(run->snapshots.size() == run->points.size());
            const auto volumes = d.mesh().volumes();
            const double rho_c = physics::silicon_parameters.thermal.heat_capacity_J_cm3K;
            const double k = base::k_B_eV_per_K, q = base::q_C;
            double worst = 0.0;
            for (std::size_t s = 1; s < run->points.size(); ++s) {
                const results::TimePoint& p = run->points[s];
                const double h = p.step_s;
                // The BDF difference of a node quantity, as the run's (beta, a1, a2).
                double beta = 1.0, a1 = 1.0, a2 = 0.0;
                if (p.order == 2) {
                    const double w2 = h / run->points[s - 1].step_s, dd = 1.0 + 2.0 * w2;
                    beta = (1.0 + w2) / dd;
                    a1 = (1.0 + w2) * (1.0 + w2) / dd;
                    a2 = w2 * w2 / dd;
                }
                const auto& now = run->snapshots[s].fields;
                const auto& one = run->snapshots[s - 1].fields;
                const results::NodeFields* two =
                    s >= 2 ? &run->snapshots[s - 2].fields : nullptr;
                const auto rate = [&](const std::vector<double> results::NodeFields::*f,
                                      std::size_t i) {
                    double c = (one.*f)[i] * a1;
                    if (p.order == 2) c -= a2 * (two->*f)[i];
                    return ((now.*f)[i] - c) / (beta * h);
                };
                double lattice = 0.0, carriers = 0.0, scale = 0.0;
                const auto& bands = run->snapshots[s].bands;
                for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
                    lattice += rho_c * volumes[i] * rate(&results::NodeFields::temperature_K, i);
                    if (i == 0 || i + 1 == d.mesh().node_count()) continue;  // ohmic: no storage
                    const double T = now.temperature_K[i];
                    const double en = bands.conduction_eV[i] + 2.0 * k * T;
                    const double ep = bands.valence_eV[i] - 2.0 * k * T;
                    const double dn = rate(&results::NodeFields::n_cm3, i);
                    const double dp = rate(&results::NodeFields::p_cm3, i);
                    carriers += q * volumes[i] * (en * dn - ep * dp);
                    scale += q * volumes[i] * (std::abs(en * dn) + std::abs(ep * dp));
                }
                double power = 0.0, out = 0.0;
                for (std::size_t c = 0; c < 2; ++c) {
                    power += p.bias_V[c] * p.conduction_current[c];
                }
                for (const double hq : p.thermal_contact_heat) out += hq;
                const double imbalance = power - out - lattice - carriers;
                worst = std::max(worst, std::abs(imbalance) / (std::abs(power) + scale));
            }
            CAPTURE(sinks, integrator == solve::Integrator::bdf2, worst);
            REQUIRE(worst <= 1e-9);
            // Heated: the adiabatic bar warms everywhere, rho c L dT/dt = P on average.
            if (!sinks) REQUIRE(run->snapshots.back().fields.temperature_K[20] > 300.0);
        }
    }
}

TEST_CASE("electrothermal solve: small-signal limits and the thermal pole") {
    // The diode of the lumped gate (kappa 1e2 times silicon's: uniform T; R_th 2 K cm^2/W at each
    // end, R_eff = 1) at 0.6 V. With T uniform the lumped model is exact: T~ = Z P~,
    // Z = R_eff / (1 + i w R_eff C_th), C_th = rho c L (every box), P~ = I V~ + V I~_d and
    // I~ = Y_e V~ + I_T T~, I~_d its dissipative part (the displacement current stores, it does
    // not heat), so Y = Y_e + I_T Z (I + V Re Y_e) / (1 - Z V I_T), with Y_e the isothermal
    // device's admittance at the operating temperature (its conductance and capacitance, 1.7e-7
    // F/cm^2) and I_T its dI/dT (the carriers' stored energy, which the model also has, is
    // negligible this far below the electrical poles). Low frequency: Y equals dI/dV of the
    // electrothermal sweep; across the pole (R_eff C_th = 1.6e-4 s, 1 kHz) it follows the lumped
    // Y; well above it, the isothermal admittance at the operating temperature.
    const double R = 2.0, R_eff = 1.0, V = 0.6;
    Bar b;
    b.donors = pn_donors;
    b.acceptors = pn_acceptors;
    b.material.thermal.conductivity_W_cmK *= 1e2;  // see below
    for (auto& t : b.thermal) {
        t = {t.name, t.boundary, device::ThermalContactKind::resistance, 300.0, R};
    }
    const device::Device d = bar(b);
    solve::SmallSignalOptions ac;
    ac.steady = thermal_options();
    ac.frequencies_Hz = {1e-3, 10.0, 100.0, 1e3, 1e4, 1e5, 1e8};
    const std::vector<std::vector<double>> points{{V, 0.0}};
    const auto run = solve::solve_small_signal(d, points, ac);
    REQUIRE(run.has_value());
    REQUIRE_FALSE(run->stopped);
    const results::SmallSignalPoint& p = run->points.front();
    const double I = p.dc.terminal_current[0];
    const double T = p.dc.fields.temperature_K.front();
    // The electrothermal sweep's slope.
    const double dV = 1e-5;
    const double slope = (solved(d, {V + dV, 0.0}, ac.steady).terminal_current[0] -
                          solved(d, {V - dV, 0.0}, ac.steady).terminal_current[0]) /
                         (2.0 * dV);
    // The isothermal device at T: G, I_T and its admittance Y_e at every frequency.
    const auto isothermal = [&](double temperature) {
        Bar iso = b;
        iso.temperature_K = temperature;
        iso.thermal.clear();
        return bar(iso);
    };
    solve::BiasOptions plain;
    plain.newton.tol_update = 1e-12;
    const device::Device at_T = isothermal(T);
    const double G = (solved(at_T, {V + dV, 0.0}, plain).terminal_current[0] -
                      solved(at_T, {V - dV, 0.0}, plain).terminal_current[0]) /
                     (2.0 * dV);
    const double dT = 0.01;
    const double I_T = (solved(isothermal(T + dT), {V, 0.0}, plain).terminal_current[0] -
                        solved(isothermal(T - dT), {V, 0.0}, plain).terminal_current[0]) /
                       (2.0 * dT);
    solve::SmallSignalOptions iso_ac;
    iso_ac.steady = plain;
    iso_ac.frequencies_Hz = ac.frequencies_Hz;
    const auto iso_run = solve::solve_small_signal(at_T, points, iso_ac);
    REQUIRE(iso_run.has_value());
    const double C_th = physics::silicon_parameters.thermal.heat_capacity_J_cm3K * 1e-4;
    for (std::size_t k = 0; k < ac.frequencies_Hz.size(); ++k) {
        const double f = ac.frequencies_Hz[k];
        const std::complex<double> Y = p.admittance[k][0];
        const std::complex<double> Z =
            R_eff / (1.0 + std::complex<double>(0.0, 2.0 * std::numbers::pi * f * R_eff * C_th));
        const std::complex<double> Y_e = iso_run->points.front().admittance[k][0];
        const std::complex<double> lumped =
            Y_e + I_T * Z * (I + V * Y_e.real()) / (1.0 - Z * V * I_T);
        CAPTURE(f, Y, lumped, slope, Y_e, p.pivot_ratio[k]);
        // Measured: at most 1.4e-7 up to 1e5 Hz; at 1e8 Hz 7e-7 (lumped) and 1.0e-6 (Y_e), where
        // the bar's own diffusion time (L^2 rho c / kappa, 1e-10 s with kappa 1e2 times
        // silicon's, w tau = 0.06) and the residual thermal response show.
        const double bound = f < 1e6 ? 1e-6 : 3e-6;
        REQUIRE(std::abs(Y / lumped - 1.0) <= bound);
        if (f == 1e-3) REQUIRE(std::abs(Y.real() / slope - 1.0) <= 1e-6);
        if (f == 1e8) REQUIRE(std::abs(Y / Y_e - 1.0) <= bound);
    }
    // The heating's share at low frequency is large enough to be seen.
    REQUIRE(slope > 1.01 * G);
}

TEST_CASE("electrothermal solve: refused inputs") {
    const auto refused = [](const device::Device& d, const solve::BiasOptions& o,
                            const std::string& text, std::vector<double> bias = {0.1, 0.0}) {
        const auto p = solve::solve_bias(d, bias, o);
        REQUIRE_FALSE(p.has_value());
        CAPTURE(p.error().message);
        REQUIRE(p.error().code == ErrorCode::invalid_input);
        REQUIRE(p.error().message.find(text) != std::string::npos);
    };
    const device::Device d = bar({});
    auto o = thermal_options();
    // T7: no heat sink.
    Bar bare;
    bare.thermal.clear();
    refused(bar(bare), o, "no thermal contact");
    // T12: the quasi-static sweep.
    auto q = o;
    q.equations = solve::Equations::equilibrium_poisson;
    refused(d, q, "drift-diffusion equations", {0.0, 0.0});
    // Thermal bias: per point, per thermal contact, finite, positive, in range.
    auto t = o;
    t.thermal_bias_K = {{300.0}};
    refused(d, t, "one temperature per thermal contact");
    t.thermal_bias_K = {{300.0, 300.0}, {300.0, 300.0}};
    refused(d, t, "one list per bias point");
    t.thermal_bias_K = {{300.0, -1.0}};
    refused(d, t, "finite and positive");
    t.thermal_bias_K = {{300.0, 1000.0}};
    refused(d, t, "region 'silicon'");
    // An initial temperature of the wrong size.
    results::NodeFields initial = solved(d, {0.1, 0.0}, o).fields;
    initial.temperature_K.pop_back();
    const auto p = solve::solve_bias(d, std::vector<double>{0.1, 0.0}, o, &initial);
    REQUIRE_FALSE(p.has_value());
    REQUIRE(p.error().message.find("initial temperature") != std::string::npos);
    // T7: a small-signal operating point and a trace need a heat sink too (a transient does not;
    // see the transient energy balance).
    solve::SmallSignalOptions ac;
    ac.steady = o;
    ac.frequencies_Hz = {1e6};
    const std::vector<std::vector<double>> points{{0.1, 0.0}};
    const auto small = solve::solve_small_signal(bar(bare), points, ac);
    REQUIRE_FALSE(small.has_value());
    REQUIRE(small.error().message.find("no thermal contact") != std::string::npos);
    solve::TraceOptions trace;
    trace.steady = o;
    trace.start_V = {0.0, 0.0};
    trace.end_V = 0.5;
    const auto traced = solve::trace_bias(bar(bare), trace);
    REQUIRE_FALSE(traced.has_value());
    REQUIRE(traced.error().message.find("no thermal contact") != std::string::npos);
}
