// Continuing a transient run from its resume state (Unit 61, RD3; gate RG4): a run stopped right
// after its k-th accepted step and continued from Transient::resume, joined to its first part,
// equals the uninterrupted run bit for bit (times, steps, orders, error ratios, currents, charges,
// convergence, snapshots), for k at the start, inside BDF2's history, at a breakpoint, just after
// one and inside a ramp; with interface traps and with electrothermal on. Another run's state,
// the final step and an initial state with a resume state are refused.
#include <catch2/catch_test_macros.hpp>

#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stop_token>
#include <utility>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/transient.hpp"
#include "legacy_graded_mesh.hpp"

using namespace NiTCAD;
using solve::Waveform;

namespace {

bool same(double a, double b) { return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b); }

bool same(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (!same(a[i], b[i])) return false;
    return true;
}

void same_point(const results::TimePoint& a, const results::TimePoint& b) {
    CHECK(same(a.time_s, b.time_s));
    CHECK(same(a.step_s, b.step_s));
    CHECK(a.order == b.order);
    CHECK(same(a.error_ratio, b.error_ratio));
    CHECK(same(a.bias_V, b.bias_V));
    CHECK(same(a.terminal_current, b.terminal_current));
    CHECK(same(a.conduction_current, b.conduction_current));
    CHECK(same(a.displacement_current, b.displacement_current));
    CHECK(same(a.contact_charge, b.contact_charge));
    CHECK(same(a.interface_trap_charge, b.interface_trap_charge));
    CHECK(same(a.thermal_contact_heat, b.thermal_contact_heat));
    REQUIRE(a.convergence.iterations.size() == b.convergence.iterations.size());
    for (std::size_t i = 0; i < a.convergence.iterations.size(); ++i) {
        CHECK(same(a.convergence.iterations[i].update, b.convergence.iterations[i].update));
        CHECK(same(a.convergence.iterations[i].residual, b.convergence.iterations[i].residual));
    }
}

void same_snapshot(const results::TransientSnapshot& a, const results::TransientSnapshot& b) {
    CHECK(same(a.time_s, b.time_s));
    CHECK(same(a.fields.potential_V, b.fields.potential_V));
    CHECK(same(a.fields.n_cm3, b.fields.n_cm3));
    CHECK(same(a.fields.p_cm3, b.fields.p_cm3));
    CHECK(same(a.fields.temperature_K, b.fields.temperature_K));
    CHECK(same(a.edge_current_n, b.edge_current_n));
    CHECK(same(a.edge_current_p, b.edge_current_p));
    CHECK(same(a.trap_occupancy, b.trap_occupancy));
}

// The Unit 9 diode (anode, contact 0, on x_min), optionally with an isothermal sink on x_max.
device::Device diode(bool sink) {
    mesh::Mesh m = *mesh::make_tensor_grid(legacy_graded_mesh(2e-4, {1e-4}, 1e-8, 1e-6));
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) (m.points()[i][0] < 1e-4 ? acceptors : donors)[i] = 1e17;
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    device::DeviceDescription d{.mesh = std::move(m),
                                .temperature_K = 300.0,
                                .regions = {{"silicon", physics::silicon()}},
                                .node_region = std::vector<device::RegionId>(n, 0),
                                .donors = std::move(donors),
                                .acceptors = std::move(acceptors),
                                .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                                             {"cathode", device::ContactKind::ohmic, std::move(cathode)}}};
    if (sink) d.thermal_contacts = {{"sink", "x_max", device::ThermalContactKind::isothermal, 300.0, 0.0}};
    auto dev = device::Device::create(std::move(d));
    if (!dev) FAIL(dev.error().message);
    return std::move(*dev);
}

// A meshed-oxide MOS capacitor with an acceptor trap level at the interface (gate contact 0).
device::Device trapped_moscap() {
    const double t_ox = 1e-6, L = 1e-4;
    const std::vector<double> si = legacy_graded_mesh(L, {0.0}, 1e-7, L / 100.0);
    const double h = si[1];
    std::vector<double> x;
    for (int k = 0; k <= 10; ++k) x.push_back(-t_ox + (t_ox - 0.5 * h) * k / 10);
    for (const double v : si) x.push_back(v + 0.5 * h);
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n, 0);
    std::vector<double> acceptors(n, 1e16);
    for (std::size_t i = 0; i < n; ++i) {
        if (m.points()[i][0] < 0.0) {
            region[i] = 1;
            acceptors[i] = 0.0;
        }
    }
    device::Interface f{"oxide", "silicon"};
    f.traps.levels = {{.type = physics::TrapType::acceptor, .density_cm2 = 1e10, .energy_eV = -0.1,
                       .sigma_n_cm2 = 1e-15, .sigma_p_cm2 = 1e-15}};
    f.traps.thermal_velocity_n_cm_s = f.traps.thermal_velocity_p_cm_s = 1e7;
    device::Contact gate{"gate", device::ContactKind::electrode, m.find_boundary("x_min")->nodes};
    gate.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.1};
    auto substrate = m.find_boundary("x_max")->nodes;
    auto dev = device::Device::create({.mesh = std::move(m),
                                       .temperature_K = 300.0,
                                       .regions = {{"silicon", physics::silicon()}, {"oxide", physics::silicon_dioxide()}},
                                       .node_region = std::move(region),
                                       .donors = std::vector<double>(n, 0.0),
                                       .acceptors = std::move(acceptors),
                                       .contacts = {std::move(gate), {"substrate", device::ContactKind::ohmic, std::move(substrate)}},
                                       .interfaces = {f}});
    if (!dev) FAIL(dev.error().message);
    return std::move(*dev);
}

struct Case {
    const device::Device& device;
    std::vector<Waveform> waveforms;
    solve::TransientOptions options;
};

// The attempt that accepted each point of the uninterrupted run: the last converged attempt that
// reached the point's time (progress events: point = attempt index, time_s = the step's end).
std::vector<std::size_t> accepting_attempts(const Case& c, results::Transient& full) {
    std::vector<std::pair<std::size_t, double>> converged;
    solve::RunControl control;
    control.progress = [&](const solve::Progress& p) {
        if (p.phase == solve::Phase::transient && p.converged) converged.emplace_back(p.point, p.time_s);
    };
    auto r = solve::solve_transient(c.device, c.waveforms, c.options, nullptr, control);
    if (!r) FAIL(r.error().message);
    if (r->stopped) FAIL(r->stopped->message);
    full = std::move(*r);
    std::vector<std::size_t> out(full.points.size(), 0);
    for (std::size_t k = 1; k < full.points.size(); ++k)
        for (const auto& [attempt, t] : converged)
            if (t == full.points[k].time_s) out[k] = attempt;
    return out;
}

// Stops the run right after its k-th accepted step, continues it, and requires the joined run to
// be the uninterrupted one.
void stop_and_continue(const Case& c, const results::Transient& full, const std::vector<std::size_t>& accepting,
                       std::size_t k) {
    INFO("stopped after accepted step " << k << " of " << full.points.size() - 1);
    std::stop_source source;
    solve::RunControl control{.stop = source.get_token()};
    control.progress = [&](const solve::Progress& p) {
        if (p.phase == solve::Phase::transient && p.converged && p.point == accepting[k]) source.request_stop();
    };
    auto first = solve::solve_transient(c.device, c.waveforms, c.options, nullptr, control);
    REQUIRE(first.has_value());
    REQUIRE(first->stopped.has_value());
    CHECK(first->stopped->code == base::ErrorCode::cancelled);
    REQUIRE(first->points.size() == k + 1);
    REQUIRE(first->resume.has_value());

    auto rest = solve::solve_transient(c.device, c.waveforms, c.options, nullptr, {}, &*first->resume);
    if (!rest) FAIL(rest.error().message);
    REQUIRE_FALSE(rest->stopped.has_value());
    CHECK(rest->run.input_identity == full.run.input_identity);
    REQUIRE(rest->points.size() == full.points.size() - k);
    for (std::size_t i = 0; i < first->points.size(); ++i) same_point(first->points[i], full.points[i]);
    for (std::size_t i = 0; i < rest->points.size(); ++i) {
        INFO("point " << k + i);
        same_point(rest->points[i], full.points[k + i]);
    }
    // Snapshots: the first part's, then the continuation's after its own first (the resumed state).
    std::vector<const results::TransientSnapshot*> joined;
    for (const auto& s : first->snapshots) joined.push_back(&s);
    REQUIRE_FALSE(rest->snapshots.empty());
    CHECK(rest->snapshots.front().point == 0);
    // A snapshot's time is the scaled state time times t0, which may differ from the point's in the last bit.
    CHECK(std::abs(rest->snapshots.front().time_s - full.points[k].time_s) <= 1e-15 * full.points[k].time_s);
    for (std::size_t i = 1; i < rest->snapshots.size(); ++i) {
        CHECK(rest->snapshots[i].point + k == [&] {
            for (std::size_t j = 0; j < full.snapshots.size(); ++j)
                if (full.snapshots[j].time_s == rest->snapshots[i].time_s) return full.snapshots[j].point;
            return std::size_t{0};
        }());
        joined.push_back(&rest->snapshots[i]);
    }
    REQUIRE(joined.size() == full.snapshots.size());
    for (std::size_t i = 0; i < joined.size(); ++i) same_snapshot(*joined[i], full.snapshots[i]);
}

std::size_t point_at(const results::Transient& full, double t) {
    std::size_t k = 0;
    while (k < full.points.size() && full.points[k].time_s != t) ++k;
    INFO("no point at that time");
    REQUIRE(k < full.points.size());
    return k;
}

}  // namespace

TEST_CASE("transient resume: a pulse on the diode, BDF2 adaptive, continues bit for bit (RG4)", "[solve][resume]") {
    const device::Device d = diode(false);
    Case c{d,
           {*Waveform::pulse(0.0, 0.5, 100e-12, 200e-12, 50e-12, 50e-12), Waveform::constant(0.0)},
           {.t_end_s = 600e-12, .dt_initial_s = 1e-13, .rtol = 1e-3, .output_times_s = {200e-12, 400e-12}}};
    results::Transient full;
    const auto accepting = accepting_attempts(c, full);
    REQUIRE(full.points.size() > 12);
    bool bdf2 = false;
    for (const auto& p : full.points) bdf2 = bdf2 || p.order == 2;
    REQUIRE(bdf2);
    const std::size_t at_break = point_at(full, 100e-12);  // the pulse's start
    const std::size_t at_top = point_at(full, 150e-12);    // the rise's end
    std::size_t mid_ramp = 0;
    for (std::size_t k = 0; k < full.points.size(); ++k)
        if (full.points[k].time_s > 110e-12 && full.points[k].time_s < 140e-12) mid_ramp = k;
    REQUIRE(mid_ramp > 0);
    for (const std::size_t k : {std::size_t{1}, std::size_t{2}, std::size_t{3}, at_break, at_break + 1, mid_ramp, at_top,
                                full.points.size() - 2}) {
        stop_and_continue(c, full, accepting, k);
    }
}

TEST_CASE("transient resume: with interface traps (RG4)", "[solve][resume]") {
    const device::Device d = trapped_moscap();
    Case c{d,
           {*Waveform::step(-0.5, 0.5, 1e-11), Waveform::constant(0.0)},
           {.t_end_s = 1e-9, .dt_initial_s = 1e-14, .rtol = 1e-3, .output_times_s = {5e-10}}};
    results::Transient full;
    const auto accepting = accepting_attempts(c, full);
    REQUIRE(full.points.size() > 6);
    REQUIRE_FALSE(full.snapshots.front().trap_occupancy.empty());
    for (const std::size_t k : {std::size_t{2}, full.points.size() / 2, full.points.size() - 2})
        stop_and_continue(c, full, accepting, k);
}

TEST_CASE("transient resume: with electrothermal on (RG4)", "[solve][resume]") {
    const device::Device d = diode(true);
    solve::TransientOptions o{.t_end_s = 1e-8, .dt_initial_s = 1e-12, .rtol = 1e-3};
    o.steady.models.electrothermal = true;
    Case c{d, {*Waveform::step(0.0, 0.8, 1e-10), Waveform::constant(0.0)}, o};
    results::Transient full;
    const auto accepting = accepting_attempts(c, full);
    REQUIRE(full.points.size() > 6);
    REQUIRE_FALSE(full.points.back().thermal_contact_heat.empty());
    for (const std::size_t k : {std::size_t{1}, full.points.size() / 2, full.points.size() - 2})
        stop_and_continue(c, full, accepting, k);
}

TEST_CASE("transient resume: refusals (RG4)", "[solve][resume]") {
    const device::Device d = diode(false);
    const std::vector<Waveform> w{*Waveform::step(0.0, 0.3, 1e-11), Waveform::constant(0.0)};
    const solve::TransientOptions o{.t_end_s = 1e-10, .dt_initial_s = 1e-12};
    const auto full = solve::solve_transient(d, w, o);
    REQUIRE(full.has_value());
    REQUIRE(full->resume.has_value());

    // The final step: nothing left.
    const auto done = solve::solve_transient(d, w, o, nullptr, {}, &*full->resume);
    REQUIRE_FALSE(done.has_value());
    CHECK(done.error().message.find("final step") != std::string::npos);

    // A state stopped part-way, given to runs with other options or waveforms.
    std::stop_source source;
    solve::RunControl control{.stop = source.get_token()};
    int accepted = 0;
    control.progress = [&](const solve::Progress& p) {
        if (p.phase == solve::Phase::transient && p.converged && ++accepted == 3) source.request_stop();
    };
    const auto part = solve::solve_transient(d, w, o, nullptr, control);
    REQUIRE(part.has_value());
    REQUIRE(part->resume.has_value());
    solve::TransientOptions longer = o;
    longer.t_end_s = 2e-10;
    const auto other = solve::solve_transient(d, w, longer, nullptr, {}, &*part->resume);
    REQUIRE_FALSE(other.has_value());
    CHECK(other.error().message.find("another run") != std::string::npos);
    const std::vector<Waveform> w2{*Waveform::step(0.0, 0.4, 1e-11), Waveform::constant(0.0)};
    CHECK_FALSE(solve::solve_transient(d, w2, o, nullptr, {}, &*part->resume).has_value());

    // An initial state together with a resume state.
    results::NodeFields guess = full->snapshots.front().fields;
    CHECK_FALSE(solve::solve_transient(d, w, o, &guess, {}, &*part->resume).has_value());

    // A damaged state.
    results::TransientResume broken = *part->resume;
    broken.history.front().x.pop_back();
    CHECK_FALSE(solve::solve_transient(d, w, o, nullptr, {}, &broken).has_value());
    broken = *part->resume;
    broken.history.clear();
    CHECK_FALSE(solve::solve_transient(d, w, o, nullptr, {}, &broken).has_value());

    // A run that stops before accepting any step has no resume state.
    std::stop_source early;
    const auto none = solve::solve_transient(d, w, o, nullptr, {.stop = early.get_token(), .progress = [&](const solve::Progress&) {
                                                 early.request_stop();
                                             }});
    REQUIRE(none.has_value());
    CHECK_FALSE(none->resume.has_value());
}

TEST_CASE("transient resume: the carrier-dependent Klaassen mobility continues bit for bit (Unit 66, KC9)",
          "[solve][resume]") {
    const device::Device d = diode(false);
    Case c{d,
           {*Waveform::pulse(0.0, 0.5, 100e-12, 200e-12, 50e-12, 50e-12), Waveform::constant(0.0)},
           {.t_end_s = 400e-12, .dt_initial_s = 1e-13, .rtol = 1e-3, .output_times_s = {200e-12}}};
    c.options.steady.models.electron_mobility = assemble::ElectronMobility::klaassen;
    c.options.steady.models.hole_mobility = assemble::HoleMobility::klaassen;
    results::Transient full;
    const auto accepting = accepting_attempts(c, full);
    REQUIRE(full.points.size() > 6);
    std::size_t mid_ramp = 0;
    for (std::size_t k = 0; k < full.points.size(); ++k)
        if (full.points[k].time_s > 110e-12 && full.points[k].time_s < 140e-12) mid_ramp = k;
    REQUIRE(mid_ramp > 0);
    for (const std::size_t k : {std::size_t{2}, point_at(full, 100e-12), mid_ramp, full.points.size() - 2})
        stop_and_continue(c, full, accepting, k);
}
