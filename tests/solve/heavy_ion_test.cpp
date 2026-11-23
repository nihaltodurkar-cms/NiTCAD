// Heavy-ion tracks in transient runs (Unit 65; DECISIONS.md HD2-HD7, HD10; gates HC1, HC3-HC6,
// HC8): the mesh deposit against LET L, the source's charge over a run, collection by a
// reverse-biased diode, the steps around the pulse, the run record, refusals and restart.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <stop_token>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/assemble/system.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/heavy_ion.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/transient.hpp"
#include "legacy_graded_mesh.hpp"

using namespace NiTCAD;
using physics::HeavyIonTrack;
using solve::Waveform;

namespace {

std::vector<double> uniform(double a, double b, int cells) {
    std::vector<double> v;
    for (int k = 0; k <= cells; ++k) v.push_back(a + (b - a) * k / cells);
    return v;
}

// A p-n diode on mesh m, junction at x = xj (acceptors below), 1e17 each side; anode (contact 0)
// on x_min, cathode on x_max.
device::Device diode(mesh::Mesh m, double xj) {
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) (m.points()[i][0] < xj ? acceptors : donors)[i] = 1e17;
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    auto d = device::Device::create({.mesh = std::move(m),
                                     .temperature_K = 300.0,
                                     .regions = {{"silicon", physics::silicon()}},
                                     .node_region = std::vector<device::RegionId>(n, 0),
                                     .donors = std::move(donors),
                                     .acceptors = std::move(acceptors),
                                     .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                                                  {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
    if (!d) FAIL(d.error().message);
    return std::move(*d);
}

device::Device diode_1d() { return diode(*mesh::make_tensor_grid(legacy_graded_mesh(2e-4, {1e-4}, 1e-8, 1e-6)), 1e-4); }

std::vector<double> mesh_deposits(const device::Device& d, const std::vector<HeavyIonTrack>& tracks) {
    assemble::PhysicsModels models;
    models.heavy_ions = &tracks;
    auto scaling = assemble::make_scaling(d);
    REQUIRE(scaling.has_value());
    auto system = assemble::System::create(d, *scaling, models);
    if (!system) FAIL(system.error().message);
    return system->heavy_ion_deposits();
}

// A 1D track through the junction of diode_1d: 0.9 to 1.1 um, 1e-3 pC/um over 1 um^2.
HeavyIonTrack junction_track_1d() {
    return {.start_cm = {0.9e-4, 0.0, 0.0}, .direction = {1.0, 0.0, 0.0}, .length_cm = 0.2e-4,
            .let_C_per_cm = 1e-11, .area_cm2 = 1e-8, .t0_s = 20e-12, .sigma_s = 2e-12};
}

// The collected charge of contact c: the trapezoid integral of I_c(t) - I_c(0) over the points.
double collected(const results::Transient& run, std::size_t c) {
    double q = 0.0;
    const double i0 = run.points.front().terminal_current[c];
    for (std::size_t k = 1; k < run.points.size(); ++k) {
        const auto& a = run.points[k - 1];
        const auto& b = run.points[k];
        q += 0.5 * (b.time_s - a.time_s) * ((a.terminal_current[c] - i0) + (b.terminal_current[c] - i0));
    }
    return q;
}

solve::TransientOptions reverse_options(std::vector<HeavyIonTrack> tracks) {
    solve::TransientOptions o{.t_end_s = 1e-9, .dt_initial_s = 1e-14};
    o.steady.models.srh = false;
    o.steady.models.auger = false;
    o.steady.models.radiative = false;
    o.heavy_ions = std::move(tracks);
    return o;
}

std::vector<Waveform> reverse_bias() { return {Waveform::constant(-5.0), Waveform::constant(0.0)}; }

}  // namespace

TEST_CASE("heavy ion: the mesh deposit is LET L on a resolving mesh (HC3)") {
    const double w = 5e-6, h = w / 4;  // 0.05 um radius, spacing w / 4
    // 1D: exact (R is constant); ends between nodes.
    {
        const device::Device d = diode(*mesh::make_tensor_grid(uniform(0.0, 2e-4, 160)), 1e-4);
        HeavyIonTrack t = junction_track_1d();
        t.start_cm[0] = 0.90625e-4;  // between nodes (spacing 1.25e-6 cm)
        const double got = mesh_deposits(d, {t}).at(0);
        const double expected = t.let_C_per_cm * t.length_cm / t.area_cm2;
        CAPTURE(got, expected);
        CHECK(std::abs(got / expected - 1.0) <= 1e-12);
    }
    // 2D: track along x at mid-height, ends between nodes.
    {
        const auto x = uniform(0.0, 4e-5, 32), y = uniform(0.0, 4e-5, 32);
        const device::Device d = diode(*mesh::make_tensor_grid(x, y), 2e-5);
        const HeavyIonTrack t{.start_cm = {1e-5 + h / 2, 2e-5, 0.0}, .direction = {1.0, 0.0, 0.0},
                              .length_cm = 2e-5, .let_C_per_cm = 1e-11, .radius_cm = w, .depth_cm = 1e-4,
                              .t0_s = 1e-11, .sigma_s = 1e-12};
        const double got = mesh_deposits(d, {t}).at(0);
        const double expected = t.let_C_per_cm * t.length_cm / t.depth_cm;
        std::printf("HC3 2D: Q_m / Q_a - 1 = %.3e\n", got / expected - 1.0);
        CHECK(std::abs(got / expected - 1.0) <= 0.01);
    }
    // 3D: track along z through the centre.
    {
        const auto x = uniform(0.0, 4e-5, 32), y = uniform(0.0, 4e-5, 32), z = uniform(0.0, 4e-5, 32);
        const device::Device d = diode(*mesh::make_tensor_grid(x, y, z), 2e-5);
        const HeavyIonTrack t{.start_cm = {2e-5, 2e-5, 1e-5 + h / 2}, .direction = {0.0, 0.0, 1.0},
                              .length_cm = 2e-5, .let_C_per_cm = 1e-11, .radius_cm = w,
                              .t0_s = 1e-11, .sigma_s = 1e-12};
        const double got = mesh_deposits(d, {t}).at(0);
        const double expected = t.let_C_per_cm * t.length_cm;
        std::printf("HC3 3D: Q_m / Q_a - 1 = %.3e\n", got / expected - 1.0);
        CHECK(std::abs(got / expected - 1.0) <= 0.01);
    }
}

TEST_CASE("heavy ion: the source's charge over a run and the steps around the pulse (HC4, HC6)") {
    const device::Device d = diode_1d();
    const HeavyIonTrack t = junction_track_1d();
    const double s = t.sigma_s, t0 = t.t0_s;
    struct Variant {
        solve::Integrator integrator;
        bool adaptive;
        double dt0;
    };
    const double t_end = 1e-10;
    for (const Variant v : {Variant{solve::Integrator::bdf2, true, 1e-14}, Variant{solve::Integrator::backward_euler, true, 1e-14},
                            Variant{solve::Integrator::bdf2, true, t_end / 10}, Variant{solve::Integrator::bdf2, false, t_end / 10},
                            Variant{solve::Integrator::backward_euler, false, t_end / 10}}) {
        INFO("integrator " << static_cast<int>(v.integrator) << ", adaptive " << v.adaptive << ", dt0 " << v.dt0);
        solve::TransientOptions o = reverse_options({t});
        o.t_end_s = t_end;
        o.integrator = v.integrator;
        o.adaptive = v.adaptive;
        o.dt_initial_s = v.dt0;
        const auto run = solve::solve_transient(d, reverse_bias(), o);
        if (!run) FAIL(run.error().message);
        REQUIRE(!run->stopped);
        // HC6: the landings, and no step inside the window longer than s / 4.
        std::size_t inside = 0;
        for (const double land : {t0 - 4 * s, t0, t0 + 4 * s}) {
            CHECK(std::ranges::any_of(run->points, [&](const auto& p) { return p.time_s == land; }));
        }
        double charge = 0.0;
        for (std::size_t k = 1; k < run->points.size(); ++k) {
            const double a = run->points[k - 1].time_s, b = run->points[k].time_s;
            if (a >= t0 - 4 * s && a < t0 + 4 * s) {
                ++inside;
                // The step asked is at most s / 4; the difference of the two rounded times can
                // exceed it by the rounding of the sum t + h (an ulp of b).
                CHECK(b - a <= 0.25 * s + 2.0 * (std::nextafter(b, 1.0) - b));
            }
            charge += physics::heavy_ion_time_fraction(t, a, b);
        }
        CHECK(inside >= 32);
        // HC4: the source's charge (fraction of the deposit) is the exact integral.
        const double exact = 1.0 - 0.5 * std::erfc(t0 / s) - 0.5 * std::erfc((t_end - t0) / s);
        CAPTURE(charge, exact, inside, run->points.size());
        CHECK(std::abs(charge / exact - 1.0) <= 1e-12);
    }
}

// Release only (hidden tag [.heavy_ion], about 80 s in Release), as the other slow gates.
TEST_CASE("heavy ion: a reverse-biased diode collects the deposit (HC5)", "[.heavy_ion]") {
    // 1D.
    {
        const device::Device d = diode_1d();
        const HeavyIonTrack t = junction_track_1d();
        const auto run = solve::solve_transient(d, reverse_bias(), reverse_options({t}));
        if (!run) FAIL(run.error().message);
        REQUIRE(!run->stopped);
        const double deposit = mesh_deposits(d, {t}).at(0) * physics::heavy_ion_time_fraction(t, 0.0, 1e-9);
        const double q_cathode = collected(*run, 1), q_anode = collected(*run, 0);
        std::printf("HC5 1D: deposit %.6e C/cm^2, cathode %.6e, anode %.6e, ratio %.6f (%zu points)\n", deposit,
                    q_cathode, q_anode, std::abs(q_cathode) / deposit, run->points.size());
        CHECK(std::abs(std::abs(q_cathode) / deposit - 1.0) <= 0.01);
        CHECK(run->run.warnings.empty());
    }
    // 2D: y-uniform diode, the track along x at mid-height.
    {
        const auto x = legacy_graded_mesh(2e-4, {1e-4}, 1e-7, 4e-6);  // coarser than 1D: run time
        const auto y = uniform(0.0, 1e-4, 20);
        const device::Device d = diode(*mesh::make_tensor_grid(x, y), 1e-4);
        const HeavyIonTrack t{.start_cm = {0.9e-4, 0.5e-4, 0.0}, .direction = {1.0, 0.0, 0.0},
                              .length_cm = 0.2e-4, .let_C_per_cm = 1e-11, .radius_cm = 1e-5,
                              .depth_cm = 1e-4, .t0_s = 20e-12, .sigma_s = 2e-12};
        const auto run = solve::solve_transient(d, reverse_bias(), reverse_options({t}));
        if (!run) FAIL(run.error().message);
        REQUIRE(!run->stopped);
        const double deposit = mesh_deposits(d, {t}).at(0) * physics::heavy_ion_time_fraction(t, 0.0, 1e-9);
        const double q_cathode = collected(*run, 1);
        std::printf("HC5 2D: deposit %.6e C/cm, cathode %.6e, ratio %.6f (%zu points)\n", deposit, q_cathode,
                    std::abs(q_cathode) / deposit, run->points.size());
        CHECK(std::abs(std::abs(q_cathode) / deposit - 1.0) <= 0.01);
    }
    // With SRH on and the track deep into the neutral p side, less is collected.
    {
        const device::Device d = diode_1d();
        HeavyIonTrack t = junction_track_1d();
        t.start_cm[0] = 0.3e-4;
        t.length_cm = 0.8e-4;
        solve::TransientOptions o = reverse_options({t});
        o.steady.models.srh = true;
        const auto run = solve::solve_transient(d, reverse_bias(), o);
        if (!run) FAIL(run.error().message);
        REQUIRE(!run->stopped);
        const double deposit = mesh_deposits(d, {t}).at(0) * physics::heavy_ion_time_fraction(t, 0.0, 1e-9);
        const double q_cathode = collected(*run, 1);
        std::printf("HC5 SRH: deposit %.6e, cathode %.6e, ratio %.6f\n", deposit, q_cathode,
                    std::abs(q_cathode) / deposit);
        CHECK(std::abs(q_cathode) < deposit);
    }
}

TEST_CASE("heavy ion: run record, warnings and refusals (HC1, HD6, HD7, HD10)") {
    const device::Device d = diode_1d();
    const auto w = reverse_bias();
    solve::TransientOptions plain{.t_end_s = 1e-10};
    solve::TransientOptions with = plain;
    with.heavy_ions = {junction_track_1d()};
    const auto r0 = solve::make_run_record(d, plain, w);
    const auto r1 = solve::make_run_record(d, with, w);
    CHECK(r0.input_identity != r1.input_identity);
    CHECK(r1.settings.size() == r0.settings.size() + 5);
    CHECK(std::ranges::none_of(r0.settings, [](const auto& s) { return s.first.starts_with("heavy_ion"); }));
    solve::TransientOptions other = with;
    other.heavy_ions[0].let_C_per_cm *= 2.0;
    CHECK(solve::make_run_record(d, other, w).input_identity != r1.input_identity);
    // Refusals.
    const auto refused = [&](solve::TransientOptions o) {
        const auto r = solve::solve_transient(d, w, o);
        if (r) return false;
        INFO(r.error().message);
        return r.error().code == base::ErrorCode::invalid_input;
    };
    {
        auto o = with;
        o.steady.models.energy_balance = true;
        CHECK(refused(o));
    }
    {
        auto o = with;
        o.heavy_ions[0].start_cm[0] = 5e-4;  // beyond the device
        CHECK(refused(o));
    }
    {
        auto o = with;
        o.heavy_ions[0].t0_s = 3.0 * o.heavy_ions[0].sigma_s;
        CHECK(refused(o));
    }
    {
        auto o = with;
        o.heavy_ions[0].area_cm2 = 0.0;
        CHECK(refused(o));
    }
    // HD7: a track reaching past the device warns (its deposit is half LET L).
    {
        auto o = with;
        o.t_end_s = 4e-11;
        o.heavy_ions[0].start_cm[0] = 1.9e-4;
        const auto r = solve::solve_transient(d, w, o);
        if (!r) FAIL(r.error().message);
        REQUIRE(r->run.warnings.size() == 1);
        CHECK(r->run.warnings[0].starts_with("heavy-ion track 0"));
    }
}

TEST_CASE("heavy ion: a run stopped mid-pulse and continued is the uninterrupted one (HC8)") {
    const device::Device d = diode_1d();
    solve::TransientOptions o = reverse_options({junction_track_1d()});
    o.t_end_s = 1e-10;
    const auto w = reverse_bias();
    const auto full = solve::solve_transient(d, w, o);
    if (!full) FAIL(full.error().message);
    // Stop after the first accepted step past t0.
    std::size_t k = 0;
    while (full->points[k].time_s <= 20e-12) ++k;
    std::stop_source source;
    std::size_t accepted = 0;
    solve::RunControl control{.stop = source.get_token()};
    control.progress = [&](const solve::Progress& p) {
        if (p.phase == solve::Phase::transient && p.converged && ++accepted == k) source.request_stop();
    };
    const auto first = solve::solve_transient(d, w, o, nullptr, control);
    REQUIRE(first.has_value());
    REQUIRE(first->stopped.has_value());
    REQUIRE(first->resume.has_value());
    const auto rest = solve::solve_transient(d, w, o, nullptr, {}, &*first->resume);
    if (!rest) FAIL(rest.error().message);
    REQUIRE(!rest->stopped);
    const std::size_t from = first->points.size() - 1;
    REQUIRE(rest->points.size() == full->points.size() - from);
    CHECK(rest->run.input_identity == full->run.input_identity);
    const auto same = [](double a, double b) { return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b); };
    bool equal = true;
    for (std::size_t i = 0; i < rest->points.size(); ++i) {
        const auto& a = rest->points[i];
        const auto& b = full->points[from + i];
        equal = equal && same(a.time_s, b.time_s);
        for (std::size_t c = 0; c < 2; ++c) equal = equal && same(a.terminal_current[c], b.terminal_current[c]);
    }
    CAPTURE(from, full->points.size());
    CHECK(equal);
}
