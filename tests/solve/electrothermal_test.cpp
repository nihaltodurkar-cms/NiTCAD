// Electrothermal steady solves (ARCHITECTURE.md section 11, Unit 23; DECISIONS.md T1-T14 and the
// owner-approved gates): the heat equation alone against its Kirchhoff-transform solution
// (isothermal and R_th ends); the steady energy balance (heat out of the sinks equals sum I V);
// no current or heat at a uniform temperature at equilibrium; a conductivity 1e8 times larger
// against the isothermal solve, and a uniform 350 K against the isothermal device at 350 K; the
// Seebeck open-circuit voltage against the integral of the thermopower; the run record; and the
// inputs a steady electrothermal solve refuses.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/thermal.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "NiTCAD/solve/trace.hpp"
#include "NiTCAD/solve/transient.hpp"

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
    // sinks equals sum_c V_c I_c (Joule, recombination, Peltier and Thomson heat together), to
    // round-off where the heating is resolved. Measured relative error 1.7e-10 at 0.5 V (a rise
    // of 1.6e-7 K), 7.8e-13 at 0.7 V, 1.3e-14 at 0.8 V. Below that the balance is limited by the
    // Newton tolerance, an absolute bound on the rise (T10): at 0.3 V the rise is 1.5e-10 K, under
    // the tolerance's 3e-10 K, and the balance closed to 6.4e-7.
    Bar b;
    b.donors = pn_donors;
    b.acceptors = pn_acceptors;
    b.thermal[1] = {"right", "x_max", device::ThermalContactKind::resistance, 300.0, 1e-3};
    const device::Device d = bar(b);
    auto o = thermal_options();
    std::vector<std::vector<double>> points;
    for (const double V : {0.5, 0.7, 0.8}) points.push_back({V, 0.0});
    auto sweep = solve::sweep_bias(d, points, o);
    REQUIRE(sweep.has_value());
    REQUIRE_FALSE(sweep->stopped);
    for (const results::BiasPoint& p : sweep->points) {
        const double power = electrical_power(p);
        const double heat = p.thermal_contact_heat[0] + p.thermal_contact_heat[1];
        CAPTURE(p.bias_V[0], power, heat, p.fields.temperature_K.back());
        REQUIRE(power > 0.0);
        REQUIRE(std::abs(heat - power) <= 1e-9 * power);
    }
    // The last point heats the R_th end measurably (6.5 mK: the isothermal end, 1 um away,
    // takes most of the heat).
    REQUIRE(sweep->points.back().fields.temperature_K.back() > 300.005);
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
    // Not built yet: trace, transient, small signal.
    solve::TraceOptions trace;
    trace.steady = o;
    trace.start_V = {0.0, 0.0};
    trace.end_V = 0.5;
    const auto traced = solve::trace_bias(d, trace);
    REQUIRE_FALSE(traced.has_value());
    REQUIRE(traced.error().message.find("trace_bias") != std::string::npos);
    solve::SmallSignalOptions ac;
    ac.steady = o;
    ac.frequencies_Hz = {1e6};
    const std::vector<std::vector<double>> points{{0.1, 0.0}};
    const auto small = solve::solve_small_signal(d, points, ac);
    REQUIRE_FALSE(small.has_value());
    REQUIRE(small.error().message.find("small-signal") != std::string::npos);
}
