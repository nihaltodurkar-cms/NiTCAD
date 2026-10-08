// Transient runs (ARCHITECTURE.md section 11, Unit 21): waveforms; the steady limits (a constant
// bias stays put, a long run ends at the steady state of its final bias); the conservation of the
// total current with displacement current; the measured order of backward Euler and BDF2; the
// analytic RC and dielectric-relaxation responses of a MOS capacitor (meshed oxide and lumped
// gate); interface trap emission; 2D and 3D extrusions; errors, cancellation and progress.
// The five slowest (hidden tag [.transient], 81 to 280 s each in Debug, about 9 s together in
// Release) run in the Release configuration only, as the MOSFET gates do.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stop_token>
#include <tuple>
#include <utility>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/transient.hpp"
#include "legacy_graded_mesh.hpp"
#include "legacy_turnoff.hpp"

using namespace NiTCAD;
using base::ErrorCode;
using solve::Waveform;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

// The Unit 9 diode: legacy graded mesh, 1e17 / 1e17, junction at 1 um; anode (contact 0) on x_min.
device::Device diode(mesh::Mesh m) {
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        (m.points()[i][0] < 1e-4 ? acceptors : donors)[i] = 1e17;
    }
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
}

std::vector<double> diode_axis() { return legacy_graded_mesh(2e-4, {1e-4}, 1e-8, 1e-6); }

std::vector<Waveform> waveforms(Waveform a, Waveform b) { return {std::move(a), std::move(b)}; }

double largest_abs(const std::vector<double>& v) {
    double m = 0.0;
    for (const double x : v) m = std::max(m, std::abs(x));
    return m;
}

}  // namespace

TEST_CASE("transient: waveforms") {
    const Waveform step = *Waveform::step(0.0, 1.0, 1e-9);
    REQUIRE(step.value(0.0) == 0.0);
    REQUIRE(step.left_value(1e-9) == 0.0);
    REQUIRE(step.value(1e-9) == 1.0);
    REQUIRE(step.value(2e-9) == 1.0);
    const Waveform ramp = *Waveform::ramp(0.0, 2.0, 0.0, 1e-9);
    REQUIRE(ramp.value(-1.0) == 0.0);
    REQUIRE(ramp.value(0.5e-9) == 1.0);
    REQUIRE(ramp.left_value(0.5e-9) == 1.0);
    REQUIRE(ramp.value(2e-9) == 2.0);
    // Times that add exactly in binary, so the corners are where the test expects them.
    const double t1 = 0.25e-9 * 4, w = 0.5e-9 * 4;  // 1 ns, 2 ns
    const Waveform pulse = *Waveform::pulse(0.0, 1.0, t1, w);  // legacy: jumps
    REQUIRE(pulse.value(0.0) == 0.0);
    REQUIRE(pulse.left_value(t1) == 0.0);
    REQUIRE(pulse.value(t1) == 1.0);
    REQUIRE(pulse.left_value(t1 + w) == 1.0);
    REQUIRE(pulse.value(t1 + w) == 0.0);
    REQUIRE(pulse.value(2.0 * (t1 + w)) == 0.0);
    REQUIRE(pulse.breakpoints(0.0, 1.0) == std::vector<double>{t1, t1 + w});
    // The left value at a jump is the first corner's value exactly, also where interpolating to
    // it would round (0.7 + (0.1 - 0.7) is not 0.1).
    const Waveform down = *Waveform::piecewise_linear({{0.0, 0.7}, {t1, 0.1}, {t1, 0.5}});
    REQUIRE(down.left_value(t1) == 0.1);
    REQUIRE(down.value(t1) == 0.5);
    const Waveform trapezoid = *Waveform::pulse(0.0, 1.0, t1, w, t1, t1);
    REQUIRE(close(trapezoid.value(1.5 * t1), 0.5, 1e-12));
    REQUIRE(close(trapezoid.value(2 * t1 + w + 0.5 * t1), 0.5, 1e-12));
    REQUIRE(trapezoid.breakpoints(0.0, 1.0) ==
            std::vector<double>{t1, 2 * t1, 2 * t1 + w, 3 * t1 + w});
    REQUIRE(trapezoid.breakpoints(2 * t1, 2 * t1 + w) == std::vector<double>{2 * t1 + w});
    const Waveform sine = *Waveform::sine(0.1, 0.2, 1e9, 1e-9, 0.5);
    REQUIRE(sine.value(0.0) == 0.1 + 0.2 * std::sin(0.5));
    REQUIRE(close(sine.value(1.25e-9), 0.1 + 0.2 * std::cos(0.5), 1e-12));
    REQUIRE(sine.breakpoints(0.0, 1.0) == std::vector<double>{1e-9});
    REQUIRE(Waveform::constant(0.3).value(5.0) == 0.3);
    REQUIRE(Waveform::constant(0.3).breakpoints(-1.0, 1.0) == std::vector<double>{0.0});

    REQUIRE(Waveform::piecewise_linear({}).error().code == ErrorCode::invalid_input);
    REQUIRE(Waveform::piecewise_linear({{1.0, 0.0}, {0.5, 0.0}}).error().context->index == 1);
    REQUIRE(!Waveform::piecewise_linear({{1.0, 0.0}, {1.0, 1.0}, {1.0, 2.0}}).has_value());
    REQUIRE(!Waveform::piecewise_linear({{0.0, std::nan("")}}).has_value());
    REQUIRE(!Waveform::ramp(0.0, 1.0, 1.0, 1.0).has_value());
    REQUIRE(!Waveform::pulse(0.0, 1.0, 0.0, 0.0).has_value());
    REQUIRE(!Waveform::pulse(0.0, 1.0, 0.0, 1.0, -1.0).has_value());
    REQUIRE(!Waveform::sine(0.0, 1.0, 0.0).has_value());
}

TEST_CASE("transient: a diode switched on reaches its steady state, conserving current") {
    // The Unit 9 diode stepped from 0 to 0.5 V; the transit time is about 1 ns, so by 100 ns the
    // run is at the steady state of 0.5 V: the anode current and the fields are sweep_bias's.
    const device::Device d = diode(*mesh::make_tensor_grid(diode_axis()));
    const auto wf = waveforms(*Waveform::step(0.0, 0.5, 0.0), Waveform::constant(0.0));
    const solve::TransientOptions o{.t_end_s = 1e-7, .dt_initial_s = 1e-15};
    const auto run = *solve::solve_transient(d, wf, o);
    REQUIRE(!run.stopped);
    REQUIRE(run.points.front().time_s == 0.0);
    REQUIRE(run.points.front().bias_V == std::vector<double>{0.0, 0.0});
    REQUIRE(run.points.back().time_s == 1e-7);
    REQUIRE(run.points.back().bias_V == std::vector<double>{0.5, 0.0});
    REQUIRE(run.snapshots.size() == 2);  // 0 and t_end
    const auto dc = *solve::solve_bias(d, std::vector<double>{0.5, 0.0});
    const auto& last = run.points.back();
    CAPTURE(run.points.size(), run.rejected_steps, last.terminal_current[0],
            dc.terminal_current[0]);
    // Equal to the eleven digits printed (153 steps, 4 rejected).
    REQUIRE(close(last.terminal_current[0], dc.terminal_current[0], 1e-8));
    const auto& f = run.snapshots.back().fields;
    for (std::size_t i = 0; i < f.potential_V.size(); ++i) {
        REQUIRE(std::abs(f.potential_V[i] - dc.fields.potential_V[i]) <= 1e-9);
        REQUIRE(close(f.n_cm3[i], dc.fields.n_cm3[i], 1e-7));
        REQUIRE(close(f.p_cm3[i], dc.fields.p_cm3[i], 1e-7));
    }
    double scale = 0.0, worst = 0.0;
    for (const auto& p : run.points) scale = std::max(scale, largest_abs(p.terminal_current));
    for (const auto& p : run.points) {
        worst = std::max(worst, std::abs(p.terminal_current[0] + p.terminal_current[1]) / scale);
    }
    CAPTURE(worst);
    REQUIRE(worst <= 1e-7);  // measured 4.1e-8 (the displacement spike right after the step)
}

TEST_CASE("transient: measured order of backward Euler and BDF2", "[.transient]") {
    // Fixed steps over a smooth turn-on ramp (0 to 0.6 V in 1 ns); the error of the anode's total
    // current and charge at 1 ns against a BDF2 run of 5120 steps. Measured ratios per halving:
    // backward Euler 1.92, 1.96, 1.98, 1.99; BDF2 3.33, 3.60, 3.78, 3.90 (its first two steps are
    // backward Euler, whose O(h^2) local errors fade from the ratio as h falls).
    const device::Device d = diode(*mesh::make_tensor_grid(diode_axis()));
    const double T = 1e-9;
    const auto wf = waveforms(*Waveform::ramp(0.0, 0.6, 0.0, T), Waveform::constant(0.0));
    const auto end = [&](solve::Integrator method, int steps) {
        solve::TransientOptions o{.integrator = method, .t_end_s = T, .dt_initial_s = T / steps,
                                  .adaptive = false};
        const auto run = *solve::solve_transient(d, wf, o);
        REQUIRE(!run.stopped);
        REQUIRE(run.points.size() == static_cast<std::size_t>(steps) + 1);
        return run.points.back();
    };
    const auto reference = end(solve::Integrator::bdf2, 5120);
    for (const auto method : {solve::Integrator::backward_euler, solve::Integrator::bdf2}) {
        std::vector<double> current, charge;
        for (const int steps : {20, 40, 80, 160, 320}) {
            const auto p = end(method, steps);
            current.push_back(std::abs(p.terminal_current[0] - reference.terminal_current[0]) /
                              std::abs(reference.terminal_current[0]));
            charge.push_back(std::abs(p.contact_charge[0] - reference.contact_charge[0]) /
                             std::abs(reference.contact_charge[0]));
        }
        std::vector<double> ratio;
        for (std::size_t k = 1; k < current.size(); ++k) {
            ratio.push_back(current[k - 1] / current[k]);
        }
        CAPTURE(static_cast<int>(method), current, charge, ratio);
        if (method == solve::Integrator::backward_euler) {
            REQUIRE(ratio.back() > 1.95);
            REQUIRE(ratio.back() < 2.05);
        } else {
            for (std::size_t k = 1; k < ratio.size(); ++k) REQUIRE(ratio[k] > ratio[k - 1]);
            REQUIRE(ratio.back() > 3.8);
            REQUIRE(ratio.back() < 4.2);
            REQUIRE(current.back() < 5e-5);  // 3.9e-5; backward Euler 2.7e-3 at 320 steps
        }
        for (std::size_t k = 0; k < current.size(); ++k) {
            REQUIRE(close(charge[k], current[k], 0.1));  // the charge converges alike
        }
    }
}

namespace {

// A 1D p-type MOS capacitor: a metal gate (work function 4.1 eV) on t_ox of oxide over L of
// silicon (N_A) with the substrate contact on x_max. Meshed: the oxide is a region (uniform, 20
// cells) with the electrode on x_min and the interface straddled; lumped: the Unit 12 gate on the
// silicon's x_min. The silicon mesh is graded from 1 nm at the surface.
device::Device moscap(double t_ox, double L, double N_A, bool meshed,
                      std::vector<device::Interface> interfaces = {}) {
    std::vector<double> si = legacy_graded_mesh(L, {0.0}, 1e-7, L / 200.0);
    std::vector<double> x;
    if (meshed) {
        const double h = si[1];
        for (int k = 0; k <= 20; ++k) x.push_back(-t_ox + (t_ox - 0.5 * h) * k / 20);
        for (const double v : si) x.push_back(v + 0.5 * h);
    } else {
        x = si;
    }
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n, 0);
    std::vector<double> acceptors(n, N_A);
    if (meshed) {
        for (std::size_t i = 0; i < n; ++i) {
            if (m.points()[i][0] < 0.0) {
                region[i] = 1;
                acceptors[i] = 0.0;
            }
        }
    }
    auto top = m.find_boundary("x_min")->nodes;
    auto substrate = m.find_boundary("x_max")->nodes;
    const auto kind = meshed ? device::ContactKind::electrode : device::ContactKind::gate;
    device::Contact gate{"gate", kind, std::move(top)};
    if (meshed) {
        gate.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.1};
    } else {
        gate.gate = {.boundary = "x_min",
                     .oxide_thickness_cm = t_ox,
                     .electrode = device::GateElectrode::metal,
                     .work_function_eV = 4.1};
    }
    std::vector<device::Region> regions{{"silicon", physics::silicon()}};
    if (meshed) regions.push_back({"oxide", physics::silicon_dioxide()});
    return *device::Device::create({.mesh = std::move(m),
                                    .temperature_K = 300.0,
                                    .regions = std::move(regions),
                                    .node_region = std::move(region),
                                    .donors = std::vector<double>(n, 0.0),
                                    .acceptors = std::move(acceptors),
                                    .contacts = {std::move(gate),
                                                 {"substrate", device::ContactKind::ohmic,
                                                  std::move(substrate)}},
                                    .interfaces = std::move(interfaces)});
}

// The small-signal response of the gate charge to a gate step dV at V0: a capacitance C (the
// oxide in series with the semiconductor's, the quasi-static dQ/dV) in series with the
// substrate, a resistance R = L / (q mu_p N_A) shunted by its geometric capacitance
// C_g = eps_si / L, charges with tau = R (C + C_g). Returns (measured tau, R (C + C_g)).
std::pair<double, double> rc_response(const device::Device& d, double L, double N_A, double V0,
                                      double dV) {
    // C from the steady gate charge at V0 -+ dV.
    const auto q_at = [&](double V) {
        return solve::solve_bias(d, std::vector<double>{V, 0.0})->gate_charge[0];
    };
    const double C = (q_at(V0 + dV) - q_at(V0 - dV)) / (2.0 * dV);
    const physics::Semiconductor si = physics::silicon();
    const double mu = physics::caughey_thomas_mobility(si, physics::Carrier::hole, N_A, 300.0);
    const double R = L / (base::q_C * mu * N_A);
    const double Cg = si.parameters().eps_r * base::eps0_F_per_cm / L;
    const double tau = R * (C + Cg);
    const auto wf = waveforms(*Waveform::step(V0, V0 + dV, 0.0), Waveform::constant(0.0));
    solve::TransientOptions o{.t_end_s = 12.0 * tau,
                              .dt_initial_s = 1e-6 * tau,
                              .rtol = 1e-8,
                              .output_times_s = {tau, 4.0 * tau}};
    const auto run = *solve::solve_transient(d, wf, o, nullptr);
    REQUIRE(!run.stopped);
    const double q_end = q_at(V0 + dV);
    double q1 = 0.0, q2 = 0.0;
    for (const auto& p : run.points) {
        if (p.time_s == tau) q1 = p.contact_charge[0];
        if (p.time_s == 4.0 * tau) q2 = p.contact_charge[0];
    }
    // The decay rate between tau and 4 tau, against the steady charge at V0 + dV.
    const double measured = 3.0 * tau / std::log((q_end - q1) / (q_end - q2));
    // The run ends at the steady state.
    const double q_start = run.points.front().contact_charge[0];
    const double end_error = (run.points.back().contact_charge[0] - q_end) / (q_end - q_start);
    CAPTURE(C, Cg, R, tau, measured, run.points.size(), run.rejected_steps, end_error);
    REQUIRE(std::abs(end_error) < 1e-4);
    return {measured, tau};
}

}  // namespace

TEST_CASE("transient: MOS capacitor RC response, meshed oxide and lumped gate") {
    // 10 nm of oxide over 100 um of 1e16 p-silicon, accumulated at -2 V and stepped by 1 mV:
    // tau = R (C + C_g) = 5.0 ns. Measured 0.092% (meshed) and 0.093% (lumped) short; the local
    // rate is constant to 1e-4 from 0.05 tau to 3 tau, so the response is the single exponential
    // of the circuit, and the 0.07% left is the accumulation layer's share of the substrate (its
    // holes conduct better than the bulk's).
    for (const bool meshed : {true, false}) {
        const double L = 1e-2, N_A = 1e16;
        const device::Device d = moscap(1e-6, L, N_A, meshed);
        const auto [measured, tau] = rc_response(d, L, N_A, -2.0, 1e-3);
        CAPTURE(meshed, measured / tau - 1.0);
        REQUIRE(close(measured, tau, 3e-3));
    }
}

TEST_CASE("transient: dielectric relaxation behind a thick oxide") {
    // 100 um of oxide over 10 um of 1e16 p-silicon: C (3.5e-11 F/cm^2) is 3% of the substrate's
    // geometric C_g, so tau = R (C + C_g) is 1.03 times the dielectric relaxation time
    // eps / (q mu_p N_A) = 1.59 ps. Measured 0.34% short: on this time scale the surface layer
    // (a Debye length, 0.4% of the substrate) is no longer quasi-static.
    const double L = 1e-3, N_A = 1e16;
    const device::Device d = moscap(1e-2, L, N_A, true);
    const auto [measured, tau] = rc_response(d, L, N_A, -1.0, 1e-3);
    const physics::Semiconductor si = physics::silicon();
    const double mu = physics::caughey_thomas_mobility(si, physics::Carrier::hole, N_A, 300.0);
    const double relaxation =
        si.parameters().eps_r * base::eps0_F_per_cm / (base::q_C * mu * N_A);
    CAPTURE(measured, tau, relaxation, measured / tau - 1.0);
    REQUIRE(tau / relaxation < 1.04);
    REQUIRE(close(measured, tau, 1e-2));
}

TEST_CASE("transient: an interface trap level empties at its emission rate", "[.transient]") {
    // An acceptor level 0.3 eV below midgap (N_t = 1e9 cm^-2, sigma 1e-15 cm^2, v_th 1e7 cm/s)
    // under a 10 nm oxide on 1e16 p-silicon. Near flat band (-0.9 V) the holes keep it nearly
    // empty (f = 0.055); a gate step to 0 V (weak inversion: 7e3 holes and 3e10 electrons per cm^3
    // at the surface) leaves it without carriers to capture, so it fills by hole emission and
    // relaxes at B = cp (p_I + p1) + cn (n_I + n1), dominated by e_p = cp p1 =
    // sigma v_th n_i e^(0.3 eV / kT) = 1.17e7 /s. The carriers settle in R C = 0.5 ns, so the trap
    // transient is a single exponential long after; the trap density is small enough to leave
    // the surface potential (and B) unchanged.
    const double N_A = 1e16, L = 1e-3, E = -0.3, sigma = 1e-15, v = 1e7;
    device::Interface f{"oxide", "silicon"};
    f.traps.levels = {{.type = physics::TrapType::acceptor, .density_cm2 = 1e9, .energy_eV = E,
                       .sigma_n_cm2 = sigma, .sigma_p_cm2 = sigma}};
    f.traps.thermal_velocity_n_cm_s = f.traps.thermal_velocity_p_cm_s = v;
    const device::Device d = moscap(1e-6, L, N_A, true, {f});
    const physics::Semiconductor si = physics::silicon();
    const double V_T = base::thermal_voltage(300.0);
    const double ni = physics::intrinsic_density(si, 300.0);
    const double c = sigma * v, p1 = ni * std::exp(-E / V_T), n1 = ni * std::exp(E / V_T);
    const double tau = 1.0 / (c * p1);
    const double V1 = -0.9, V2 = 0.0;
    const auto wf = waveforms(*Waveform::step(V1, V2, 0.0), Waveform::constant(0.0));
    solve::TransientOptions o{.t_end_s = 20.0 * tau,
                              .dt_initial_s = 1e-15,
                              .rtol = 1e-6,
                              .output_times_s = {tau, 3.0 * tau}};
    const auto run = *solve::solve_transient(d, wf, o);
    REQUIRE(!run.stopped);
    double f1 = 0.0, f2 = 0.0;
    std::size_t surface = 0;  // the first silicon node
    for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
        if (d.mesh().points()[i][0] > 0.0) {
            surface = i;
            break;
        }
    }
    double ps = 0.0, ns = 0.0;
    for (const auto& s : run.snapshots) {
        REQUIRE(s.trap_occupancy.size() == 1);
        if (s.time_s == tau) f1 = s.trap_occupancy[0];
        if (s.time_s == 3.0 * tau) f2 = s.trap_occupancy[0];
    }
    const auto& first = run.snapshots.front();
    const auto& last = run.snapshots.back();
    ps = last.fields.p_cm3[surface];
    ns = last.fields.n_cm3[surface];
    const double f0 = first.trap_occupancy[0], f_end = last.trap_occupancy[0];
    const double measured = std::log((f_end - f1) / (f_end - f2)) / (2.0 * tau);
    const double B = c * (p1 + ps + n1 + ns);
    // The end is the steady state at V2.
    const auto dc = *solve::solve_bias(d, std::vector<double>{V2, 0.0});
    const double q_end = run.points.back().interface_trap_charge[0];
    CAPTURE(f0, f1, f2, f_end, ps, ns, p1, B, measured, run.points.size(), q_end,
            dc.interface_trap_charge[0]);
    REQUIRE(f0 < 0.1);
    REQUIRE(f_end > 0.999);
    // n_I and p_I differ from the node's half a cell away by a factor near 1; they are 1e-5 of
    // p1 here, so that does not show.
    REQUIRE((ps + ns) < 1e-4 * p1);
    REQUIRE(close(measured, B, 1e-3));  // measured 3.2e-4
    REQUIRE(close(q_end, dc.interface_trap_charge[0], 1e-6));  // measured 2.2e-8
}

TEST_CASE("transient: a constant bias stays at its steady state") {
    // Started at its steady state, the diode stays there: every step's current is the steady
    // current, no displacement current flows, and the error control doubles the step each time.
    const device::Device d = diode(*mesh::make_tensor_grid(diode_axis()));
    const auto wf = waveforms(Waveform::constant(0.5), Waveform::constant(0.0));
    const auto run = *solve::solve_transient(d, wf, {.t_end_s = 1e-6, .dt_initial_s = 1e-12});
    REQUIRE(!run.stopped);
    const auto dc = *solve::solve_bias(d, std::vector<double>{0.5, 0.0});
    double worst = 0.0, displacement = 0.0;
    for (const auto& p : run.points) {
        worst = std::max(worst, std::abs(p.terminal_current[0] - dc.terminal_current[0]));
        displacement = std::max(displacement, std::abs(p.displacement_current[0]));
    }
    CAPTURE(run.points.size(), worst, displacement, dc.terminal_current[0]);
    REQUIRE(worst <= 1e-12 * std::abs(dc.terminal_current[0]));  // measured 1.9e-15 relative
    REQUIRE(displacement <= 1e-12 * std::abs(dc.terminal_current[0]));  // measured 0
    // 22 steps from 1 ps to 1 us: each at most twice the one before (the growth limit, below
    // BDF2's zero-stability bound 1 + sqrt(2)), the last shortened to land on t_end.
    REQUIRE(run.points.size() >= 20);
    REQUIRE(run.points.size() <= 25);
    for (std::size_t k = 2; k < run.points.size(); ++k) {
        REQUIRE(run.points[k].step_s <= 2.0 * run.points[k - 1].step_s * (1.0 + 1e-12));
    }
}

namespace {

// A 2D MOS structure: the meshed-oxide capacitor above (10 nm oxide, 2 um of 1e16 p-silicon,
// 5 columns over 1 um) with the electrode over half of the oxide, the substrate on x_max and a
// second ohmic contact ("side") on the upper half of the silicon part of y_max, so an interface
// edge ends on an ohmic node; a fixed charge, a donor level, an acceptor band and surface
// recombination.
device::Device mos2d() {
    std::vector<double> si = legacy_graded_mesh(2e-4, {0.0}, 1e-7, 1e-5);
    std::vector<double> x;
    const double t_ox = 1e-6, h = si[1];
    for (int k = 0; k <= 10; ++k) x.push_back(-t_ox + (t_ox - 0.5 * h) * k / 10);
    for (const double v : si) x.push_back(v + 0.5 * h);
    std::vector<double> y{0.0, 2.5e-5, 5e-5, 7.5e-5, 1e-4};
    mesh::Mesh m = *mesh::make_tensor_grid(x, y);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n, 0);
    std::vector<double> acceptors(n, 1e16);
    for (std::size_t i = 0; i < n; ++i) {
        if (m.points()[i][0] < 0.0) {
            region[i] = 1;
            acceptors[i] = 0.0;
        }
    }
    std::vector<mesh::NodeId> top, side;
    for (const mesh::NodeId v : m.find_boundary("x_min")->nodes) {
        if (m.points()[static_cast<std::size_t>(v)][1] < 5e-5) top.push_back(v);
    }
    for (const mesh::NodeId v : m.find_boundary("y_max")->nodes) {
        const double px = m.points()[static_cast<std::size_t>(v)][0];
        if (px > 0.0 && px < 1e-4) side.push_back(v);  // the upper half, off the substrate
    }
    auto substrate = m.find_boundary("x_max")->nodes;
    device::Contact gate{"gate", device::ContactKind::electrode, std::move(top)};
    gate.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.1};
    device::Interface f{"oxide", "silicon"};
    f.fixed_charge_cm2 = 2e11;
    f.traps.levels = {{.type = physics::TrapType::donor, .density_cm2 = 5e11, .energy_eV = -0.25,
                       .sigma_n_cm2 = 1e-15, .sigma_p_cm2 = 1e-15}};
    f.traps.bands = {{.type = physics::TrapType::acceptor, .density_cm2_eV = 1e12,
                      .energy_low_eV = 0.0, .energy_high_eV = 0.4}};
    f.recombination_velocity_n_cm_s = 1e3;
    f.recombination_velocity_p_cm_s = 1e3;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}, {"oxide", physics::silicon_dioxide()}},
         .node_region = std::move(region),
         .donors = std::vector<double>(n, 0.0),
         .acceptors = std::move(acceptors),
         .contacts = {std::move(gate),
                      {"substrate", device::ContactKind::ohmic, std::move(substrate)},
                      {"side", device::ContactKind::ohmic, std::move(side)}},
         .interfaces = {std::move(f)}});
}

}  // namespace

TEST_CASE("transient: the total currents sum to zero and integrate to the contact charges",
          "[.transient]") {
    // A gate ramp and a pulse on the side contact over a 2D MOS structure with traps, under both
    // integrators. Every step: the total currents of the three contacts sum to zero (to the solve's
    // tolerance). Backward Euler: the sum over the steps of h times each contact's displacement
    // current is the change of its charge (the BDF difference telescopes).
    const device::Device d = mos2d();
    const auto wf = std::vector<Waveform>{*Waveform::ramp(-1.0, 0.5, 0.0, 2e-9),
                                          Waveform::constant(0.0),
                                          *Waveform::pulse(0.0, 0.2, 5e-10, 5e-10, 1e-10, 1e-10)};
    for (const auto method : {solve::Integrator::backward_euler, solve::Integrator::bdf2}) {
        const auto run = *solve::solve_transient(
            d, wf, {.integrator = method, .t_end_s = 4e-9, .dt_initial_s = 1e-14});
        REQUIRE(!run.stopped);
        // Against the largest current of the run: at the start (steady state, no current) the
        // sum is rounding of nothing.
        double scale = 0.0;
        for (const auto& p : run.points) {
            scale = std::max({scale, largest_abs(p.conduction_current),
                              largest_abs(p.displacement_current)});
        }
        double worst = 0.0;
        std::vector<double> integral(3, 0.0);
        for (std::size_t k = 0; k < run.points.size(); ++k) {
            const auto& p = run.points[k];
            double sum = 0.0;
            for (const double I : p.terminal_current) sum += I;
            worst = std::max(worst, std::abs(sum) / scale);
            for (std::size_t c = 0; c < 3; ++c) {
                integral[c] += p.step_s * p.displacement_current[c];
            }
        }
        CAPTURE(static_cast<int>(method), run.points.size(), run.rejected_steps, worst);
        REQUIRE(worst <= 1e-8);  // measured 3.1e-10 (backward Euler), 2.3e-10 (BDF2)
        if (method == solve::Integrator::backward_euler) {
            std::vector<double> dQ(3);
            for (std::size_t c = 0; c < 3; ++c) {
                dQ[c] = run.points.back().contact_charge[c] - run.points.front().contact_charge[c];
            }
            for (std::size_t c = 0; c < 3; ++c) {
                CAPTURE(c, integral, dQ);
                // Exact but for the rounding of the sum.
                REQUIRE(std::abs(integral[c] - dQ[c]) <= 1e-12 * largest_abs(dQ));
            }
        }
    }
}

TEST_CASE("transient: y-uniform 2D and 3D diodes reproduce the 1D transient", "[.transient]") {
    const auto x = diode_axis();
    const std::vector<double> y{0.0, 5e-5, 1e-4}, z{0.0, 2.5e-5, 5e-5};
    const auto wf = waveforms(*Waveform::ramp(0.0, 0.6, 0.0, 1e-9), Waveform::constant(0.0));
    const solve::TransientOptions o{.t_end_s = 2e-9, .dt_initial_s = 2e-11, .adaptive = false};
    const auto r1 = *solve::solve_transient(diode(*mesh::make_tensor_grid(x)), wf, o);
    const auto r2 = *solve::solve_transient(diode(*mesh::make_tensor_grid(x, y)), wf, o);
    const auto r3 = *solve::solve_transient(diode(*mesh::make_tensor_grid(x, y, z)), wf, o);
    REQUIRE(r1.points.size() == 101);
    REQUIRE(r2.points.size() == r1.points.size());
    REQUIRE(r3.points.size() == r1.points.size());
    double worst = 0.0, largest = 0.0;
    for (std::size_t k = 0; k < r1.points.size(); ++k) {
        for (std::size_t c = 0; c < 2; ++c) {
            // 2D per unit depth over 1 um, 3D over 1 um x 0.5 um.
            const double I1 = r1.points[k].terminal_current[c];
            const double I2 = r2.points[k].terminal_current[c] / 1e-4;
            const double I3 = r3.points[k].terminal_current[c] / 5e-9;
            largest = std::max(largest, std::abs(I1));
            worst = std::max({worst, std::abs(I2 - I1), std::abs(I3 - I1)});
        }
    }
    CAPTURE(worst, largest);
    REQUIRE(worst <= 1e-9 * largest);  // measured 2.3e-13
}

TEST_CASE("transient: input errors") {
    const device::Device d = diode(*mesh::make_tensor_grid(diode_axis()));
    const auto wf = waveforms(Waveform::constant(0.0), Waveform::constant(0.0));
    const auto code = [&](const std::vector<Waveform>& w, const solve::TransientOptions& o) {
        const auto r = solve::solve_transient(d, w, o);
        return r ? ErrorCode{} : r.error().code;
    };
    const solve::TransientOptions ok{.t_end_s = 1e-9};
    REQUIRE(solve::solve_transient(d, wf, ok).has_value());
    REQUIRE(code({Waveform::constant(0.0)}, ok) == ErrorCode::invalid_input);
    REQUIRE(code(wf, {}) == ErrorCode::invalid_input);  // t_end missing
    REQUIRE(code(wf, {.t_end_s = std::nan("")}) == ErrorCode::invalid_input);
    REQUIRE(code(wf, {.t_end_s = 1e-9, .dt_initial_s = 0.0}) == ErrorCode::invalid_input);
    REQUIRE(code(wf, {.t_end_s = 1e-9, .dt_initial_s = 1e-12, .dt_min_s = 1e-11}) ==
            ErrorCode::invalid_input);
    REQUIRE(code(wf, {.t_end_s = 1e-9, .dt_initial_s = 1e-10, .dt_max_s = 1e-11}) ==
            ErrorCode::invalid_input);
    REQUIRE(code(wf, {.t_end_s = 1e-9, .rtol = 0.0}) == ErrorCode::invalid_input);
    REQUIRE(code(wf, {.t_end_s = 1e-9, .density_ref_cm3 = -1.0}) == ErrorCode::invalid_input);
    REQUIRE(code(wf, {.t_end_s = 1e-9, .output_times_s = {2e-9}}) == ErrorCode::invalid_input);
    REQUIRE(code(wf, {.t_end_s = 1e-9, .output_times_s = {0.0}}) == ErrorCode::invalid_input);
    solve::TransientOptions quasi = ok;
    quasi.steady.equations = solve::Equations::equilibrium_poisson;
    REQUIRE(code(wf, quasi) == ErrorCode::invalid_input);
    // A waveform value the contact rule rejects (a constant is not checked when made): named by
    // contact.
    const auto r = solve::solve_transient(d, waveforms(Waveform::constant(0.0),
                                                       Waveform::constant(std::nan(""))),
                                          ok);
    REQUIRE(!r.has_value());
    REQUIRE(r.error().code == ErrorCode::invalid_input);
    REQUIRE(r.error().context->index == 1);
    // A bad initial state.
    results::NodeFields initial;
    REQUIRE(solve::solve_transient(d, wf, ok, &initial).error().code == ErrorCode::invalid_input);
}

TEST_CASE("transient: cancellation keeps the accepted steps; progress is monotonic") {
    const device::Device d = diode(*mesh::make_tensor_grid(diode_axis()));
    const auto wf = waveforms(*Waveform::step(0.0, 0.5, 0.0), Waveform::constant(0.0));
    const solve::TransientOptions o{.t_end_s = 1e-7, .dt_initial_s = 1e-15};
    std::vector<solve::Progress> events;
    const auto full = *solve::solve_transient(
        d, wf, o, nullptr, {.progress = [&](const solve::Progress& e) { events.push_back(e); }});
    REQUIRE(!full.stopped);
    REQUIRE(!events.empty());
    REQUIRE(events.front().phase == solve::Phase::equilibrium);
    std::size_t transient_events = 0;
    for (std::size_t k = 1; k < events.size(); ++k) {
        const auto& a = events[k - 1];
        const auto& b = events[k];
        const auto key = [](const solve::Progress& e) {
            return std::tuple{static_cast<int>(e.phase), e.point, e.iteration};
        };
        REQUIRE(key(a) < key(b));
        if (b.phase == solve::Phase::transient) {
            ++transient_events;
            REQUIRE(b.point_count == 0);
            REQUIRE(b.time_s > 0.0);
        }
    }
    REQUIRE(transient_events > 0);
    // Stop after the 20th step attempt converges.
    std::stop_source source;
    const auto stopped = *solve::solve_transient(
        d, wf, o, nullptr,
        {.stop = source.get_token(), .progress = [&](const solve::Progress& e) {
             if (e.phase == solve::Phase::transient && e.point == 19 && e.converged) {
                 source.request_stop();
             }
         }});
    REQUIRE(stopped.stopped);
    REQUIRE(stopped.stopped->code == ErrorCode::cancelled);
    REQUIRE(stopped.points.size() >= 2);
    REQUIRE(stopped.points.size() < full.points.size());
    for (std::size_t k = 0; k < stopped.points.size(); ++k) {
        REQUIRE(stopped.points[k].time_s == full.points[k].time_s);
        REQUIRE(stopped.points[k].terminal_current == full.points[k].terminal_current);
    }
}

TEST_CASE("transient: the run record covers the waveforms and the time options") {
    const device::Device d = diode(*mesh::make_tensor_grid(diode_axis()));
    const auto wf = waveforms(*Waveform::step(0.0, 0.5, 0.0), Waveform::constant(0.0));
    const solve::TransientOptions o{.t_end_s = 1e-7};
    const auto id = solve::make_run_record(d, o, wf).input_identity;
    REQUIRE(solve::make_run_record(d, o, wf).input_identity == id);
    REQUIRE(solve::solve_transient(d, wf, {.t_end_s = 1e-9})->run.input_identity !=
            id);  // t_end
    auto o2 = o;
    o2.integrator = solve::Integrator::backward_euler;
    REQUIRE(solve::make_run_record(d, o2, wf).input_identity != id);
    o2 = o;
    o2.rtol = 1e-4;
    REQUIRE(solve::make_run_record(d, o2, wf).input_identity != id);
    o2 = o;
    o2.output_times_s = {5e-8};
    REQUIRE(solve::make_run_record(d, o2, wf).input_identity != id);
    const auto wf2 = waveforms(*Waveform::step(0.0, 0.6, 0.0), Waveform::constant(0.0));
    REQUIRE(solve::make_run_record(d, o, wf2).input_identity != id);
    const auto wf3 = waveforms(*Waveform::sine(0.0, 0.5, 1e9), Waveform::constant(0.0));
    REQUIRE(solve::make_run_record(d, o, wf3).input_identity != id);
    // Without adaptive control rtol does not enter.
    o2 = o;
    o2.adaptive = false;
    auto o3 = o2;
    o3.rtol = 0.5;
    REQUIRE(solve::make_run_record(d, o2, wf).input_identity ==
            solve::make_run_record(d, o3, wf).input_identity);
}

TEST_CASE("transient: the legacy diode turn-off, backward Euler on the legacy steps") {
    // The legacy fixture (test_m17_transient.py test_diode_turnoff_storage_delay_reference: p+n,
    // 1e19 / 1e15, L = 10 um, junction at 3 um, graded_mesh(L, [xj], 1e-7, 1e-6, 1.2), no BGN, no
    // Auger), run by the C++ port of the legacy diode and transient loop (legacy_turnoff.hpp) and
    // by NiTCAD with the same fixed backward-Euler steps of t_t / 100 over 3 t_t. (At the legacy
    // fixture's own t_t / 50 the legacy Newton stalls on a fixed step.) The legacy reports the
    // conduction current of the anode edge; NiTCAD's conduction current of the anode contact is the
    // same quantity.
    const double L = 1e-3, xj = 3e-4;
    const std::vector<double> x = legacy_graded_mesh(L, {xj}, 1e-7, 1e-6, 1.2);
    REQUIRE(x.size() == 1016);
    std::vector<double> doping(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) doping[i] = x[i] < xj ? -1e19 : 1e15;
    LegacyDiode1D legacy(x, doping);
    const LegacyTurnoffRun ref = legacy.turnoff(0.5, L, xj, 100, 3);
    REQUIRE(ref.current.size() == 301);

    mesh::Mesh m = *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        (doping[i] < 0.0 ? acceptors : donors)[i] = std::abs(doping[i]);
    }
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    const device::Device d = *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
    const auto wf = waveforms(*Waveform::step(0.5, 0.0, 0.0), Waveform::constant(0.0));
    solve::TransientOptions o{.integrator = solve::Integrator::backward_euler,
                              .t_end_s = 300 * ref.dt,
                              .dt_initial_s = ref.dt,
                              .adaptive = false};
    o.steady.models.auger = false;
    o.steady.models.bgn = false;
    const auto run = *solve::solve_transient(d, wf, o);
    REQUIRE(!run.stopped);
    REQUIRE(run.points.size() == ref.current.size());
    // The legacy loop stalls near equilibrium: the merit 0.5 |F|^2 stops decreasing at rounding
    // level, the line search shrinks lambda towards 0, and a step with lambda = 0 is reported as
    // converged (the clipped correction is below the tolerance), so the state, and the current,
    // freeze. The step where this starts depends on the rounding of the merit sum: the legacy run
    // froze from step 184 (current -1.4174e-6 A/cm^2), this port, the same algorithm in C++, from
    // step 190 (-1.2497e-6 A/cm^2), its line search collapsing from step 186. NiTCAD keeps
    // decaying to equilibrium (1.6e-14 A/cm^2 at 3 t_t).
    std::size_t frozen = ref.current.size();  // the first step that leaves the state unchanged
    for (std::size_t k = 1; k < ref.current.size(); ++k) {
        if (ref.current[k] == ref.current[k - 1]) {
            frozen = k;
            break;
        }
    }
    REQUIRE(frozen < ref.current.size());
    // The collapse: the run of steps with lambda below 1e-3 that ends in the stall.
    std::size_t collapse = frozen;
    while (collapse > 1 && ref.step_lambda[collapse - 1] < 1e-3) --collapse;
    double worst = 0.0;
    for (std::size_t k = 0; k < collapse; ++k) {
        worst = std::max(worst, std::abs(run.points[k].conduction_current[0] - ref.current[k]) /
                                    std::abs(ref.current[k]));
    }
    CAPTURE(collapse, frozen, worst, ref.forward_current, ref.current[frozen],
            run.points.back().terminal_current[0]);
    // Before the collapse: 1e-13 on most steps, 3.7e-5 at worst, where the current has fallen to
    // 1.7e-4 of the forward one (the same as against the legacy run's own output). Both stop
    // Newton at a correction below 1e-8 of the densities, at different iterates; the current, a
    // small difference of large fluxes, amplifies that.
    REQUIRE(worst <= 5e-5);
    // The stall: near equilibrium, a lambda of 0, and the current frozen to the end.
    REQUIRE(collapse < frozen);
    REQUIRE(std::abs(ref.current[frozen]) < 1e-4 * ref.forward_current);
    REQUIRE(ref.step_lambda[frozen + 1] == 0.0);
    REQUIRE(ref.current.back() == ref.current[frozen]);
    // The physics: the stored charge is pulled out first (a reverse current 447 times the forward
    // one on the first step), then the diode relaxes to equilibrium, which NiTCAD reaches.
    REQUIRE(run.points[1].conduction_current[0] < -400.0 * ref.forward_current);
    REQUIRE(std::abs(run.points.back().terminal_current[0]) < 1e-12 * ref.forward_current);
}

TEST_CASE("transient: with incomplete ionization the bound carriers are stored too") {
    // The dopants' bound electrons and holes follow the free ones at once (instantaneous
    // ionization), so the continuity rows store n - N_D+ and p - N_A-; then the charge balance and
    // so the total current's conservation hold as without the model. (Storing n and p alone breaks
    // it by the bound carriers' rate of change.)
    const device::Device d = diode(*mesh::make_tensor_grid(diode_axis()));
    const auto wf = waveforms(*Waveform::step(0.0, 0.5, 0.0), Waveform::constant(0.0));
    solve::TransientOptions o{.t_end_s = 1e-8, .dt_initial_s = 1e-15};
    o.steady.models.incomplete_ionization = true;
    const auto run = *solve::solve_transient(d, wf, o);
    REQUIRE(!run.stopped);
    double scale = 0.0, worst = 0.0;
    for (const auto& p : run.points) scale = std::max(scale, largest_abs(p.terminal_current));
    for (const auto& p : run.points) {
        worst = std::max(worst, std::abs(p.terminal_current[0] + p.terminal_current[1]) / scale);
    }
    CAPTURE(run.points.size(), worst);
    REQUIRE(worst <= 1e-7);
    o.steady.models.incomplete_ionization = false;
    const auto complete = *solve::solve_transient(d, wf, o);
    // The model matters here: 1e17 is partly ionized at 300 K.
    REQUIRE(!close(run.points.back().terminal_current[0],
                   complete.points.back().terminal_current[0], 1e-3));
}


TEST_CASE("transient: the error control follows rtol, across kinks and a jump", "[.transient]") {
    // A trapezoid with a jump back (0 to 0.6 V over 0.2 ns, held to 1.2 ns, then a jump to 0) on
    // the diode, adaptive BDF2 at rtol 1e-3 and 1e-4 against rtol 1e-6. (Fixed steps of 1 ps
    // are no reference here: after the jump their error is 1e-3 of the charge.) The error of the
    // anode charge just before the jump and 0.1 ns after it falls with rtol, and the step count
    // grows like rtol^(-1/3), BDF2's local error being O(h^3). Every corner is a step's end, the
    // jump's both sides included, and the integrator restarts there.
    const device::Device d = diode(*mesh::make_tensor_grid(diode_axis()));
    const auto wf = waveforms(*Waveform::piecewise_linear(
                                  {{0.0, 0.0}, {2e-10, 0.6}, {1.2e-9, 0.6}, {1.2e-9, 0.0}}),
                              Waveform::constant(0.0));
    const double jump = 1.2e-9, T = 1.3e-9;
    struct Result {
        double before, after;
        std::size_t steps;
    };
    const auto run_at = [&](double rtol) {
        const auto run = *solve::solve_transient(
            d, wf, {.t_end_s = T, .dt_initial_s = 1e-15, .rtol = rtol, .output_times_s = {jump}});
        REQUIRE(!run.stopped);
        Result r{0.0, run.points.back().contact_charge[0], run.points.size() - 1};
        bool landed = false;
        for (std::size_t k = 0; k + 1 < run.points.size(); ++k) {
            const auto& p = run.points[k];
            if (p.time_s == 2e-10) REQUIRE(run.points[k + 1].order == 1);  // a kink: restarted
            if (p.time_s != jump) continue;
            landed = true;
            r.before = p.contact_charge[0];
            REQUIRE(p.bias_V[0] == 0.6);  // the left side
            REQUIRE(run.points[k + 1].bias_V[0] == 0.0);
            REQUIRE(run.points[k + 1].order == 1);
        }
        REQUIRE(landed);
        return r;
    };
    const Result reference = run_at(1e-6);
    const double scale = std::abs(reference.before);  // the charge before the jump (from 0)
    std::vector<double> error;
    std::vector<std::size_t> steps;
    for (const double rtol : {1e-3, 1e-4}) {
        const Result r = run_at(rtol);
        error.push_back(std::max(std::abs(r.before - reference.before),
                                 std::abs(r.after - reference.after)) /
                        scale);
        steps.push_back(r.steps);
    }
    CAPTURE(error, steps, reference.steps);
    // Measured 4.1e-3 and 9.9e-4 (the charge 0.1 ns after the jump); 335, 670 and 2948 steps.
    REQUIRE(error[0] < 5e-3);
    REQUIRE(error[1] < 1.5e-3);
    REQUIRE(error[0] > 2.5 * error[1]);
    const double growth = static_cast<double>(steps[1]) / static_cast<double>(steps[0]);
    REQUIRE(growth > 1.7);  // 10^(1/3) = 2.15
    REQUIRE(growth < 2.6);
}

TEST_CASE("transient: traps on an ohmic node charge from the contact") {
    // A 1D oxide with one silicon node, which is the substrate contact: the interface's traps
    // exchange carriers with the contact itself, so the substrate's conduction current is the
    // traps' charging current, and the total currents still sum to zero.
    std::vector<double> x;
    for (int k = 0; k <= 10; ++k) x.push_back(-1e-6 + 0.99e-6 * k / 10);
    x.push_back(1e-8);
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n, 1);
    region[n - 1] = 0;
    std::vector<double> acceptors(n, 0.0);
    acceptors[n - 1] = 1e16;
    device::Contact gate{"gate", device::ContactKind::electrode, m.find_boundary("x_min")->nodes};
    gate.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.1};
    device::Interface f{"oxide", "silicon"};
    f.traps.levels = {{.type = physics::TrapType::acceptor, .density_cm2 = 1e11,
                       .energy_eV = -0.3, .sigma_n_cm2 = 1e-15, .sigma_p_cm2 = 1e-15}};
    auto substrate = m.find_boundary("x_max")->nodes;
    const device::Device d = *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}, {"oxide", physics::silicon_dioxide()}},
         .node_region = std::move(region),
         .donors = std::vector<double>(n, 0.0),
         .acceptors = std::move(acceptors),
         .contacts = {std::move(gate),
                      {"substrate", device::ContactKind::ohmic, std::move(substrate)}},
         .interfaces = {std::move(f)}});
    const auto wf = waveforms(*Waveform::step(-1.0, 1.0, 0.0), Waveform::constant(0.0));
    // Backward Euler: the conduction current's sum over the steps is exactly the trapped
    // charge's change (the traps' capture imbalance is their charge's difference over the step).
    const auto run = *solve::solve_transient(
        d, wf,
        {.integrator = solve::Integrator::backward_euler, .t_end_s = 1e-6,
         .dt_initial_s = 1e-13});
    REQUIRE(!run.stopped);
    double scale = 0.0, worst = 0.0, integral = 0.0;
    for (const auto& p : run.points) scale = std::max(scale, largest_abs(p.terminal_current));
    for (const auto& p : run.points) {
        worst = std::max(worst, std::abs(p.terminal_current[0] + p.terminal_current[1]) / scale);
        integral += p.step_s * p.conduction_current[1];
    }
    const double dQ = run.points.back().interface_trap_charge[0] -
                      run.points.front().interface_trap_charge[0];
    CAPTURE(run.points.size(), worst, integral, dQ);
    REQUIRE(dQ != 0.0);
    REQUIRE(close(integral, dQ, 1e-9));
    REQUIRE(worst <= 1e-8);
}
