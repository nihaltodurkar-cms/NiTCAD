// Transient and small-signal files (Unit 42; DECISIONS.md PV5, gates PG2 and PG4): every field
// round-trips bit for bit, NaN, -0 and infinities included; damaged files and files of another
// content kind are refused.
#include <catch2/catch_test_macros.hpp>

#include <bit>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <optional>
#include <string>
#include <vector>

#include "NiTCAD/io/files.hpp"

using namespace NiTCAD;
using base::ErrorCode;
namespace fs = std::filesystem;

namespace {

fs::path scratch(const std::string& name) {
    const fs::path dir = fs::temp_directory_path() / "nitcad_io_test" / name;
    fs::remove_all(dir);
    fs::create_directories(dir);
    return dir;
}

constexpr double nan = std::numeric_limits<double>::quiet_NaN();
constexpr double inf = std::numeric_limits<double>::infinity();
constexpr double tiny = std::numeric_limits<double>::denorm_min();

bool same(double a, double b) { return std::bit_cast<std::uint64_t>(a) == std::bit_cast<std::uint64_t>(b); }

bool same(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (!same(a[i], b[i])) return false;
    }
    return true;
}

bool same(const std::vector<std::complex<double>>& a, const std::vector<std::complex<double>>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (!same(a[i].real(), b[i].real()) || !same(a[i].imag(), b[i].imag())) return false;
    }
    return true;
}

bool same(const results::NodeFields& a, const results::NodeFields& b) {
    return same(a.potential_V, b.potential_V) && same(a.n_cm3, b.n_cm3) && same(a.p_cm3, b.p_cm3) &&
           same(a.temperature_K, b.temperature_K);
}

bool same(const results::BandDiagram& a, const results::BandDiagram& b) {
    return same(a.conduction_eV, b.conduction_eV) && same(a.valence_eV, b.valence_eV) &&
           same(a.electron_fermi_eV, b.electron_fermi_eV) && same(a.hole_fermi_eV, b.hole_fermi_eV);
}

bool same(const results::ConvergenceRecord& a, const results::ConvergenceRecord& b) {
    if (a.iterations.size() != b.iterations.size() || a.converged != b.converged) return false;
    for (std::size_t i = 0; i < a.iterations.size(); ++i) {
        const auto& x = a.iterations[i];
        const auto& y = b.iterations[i];
        if (x.iteration != y.iteration || !same(x.update, y.update) ||
            !same(x.residual, y.residual) || !same(x.step, y.step)) {
            return false;
        }
    }
    if (a.smallest_pivot_ratio.has_value() != b.smallest_pivot_ratio.has_value()) return false;
    return !a.smallest_pivot_ratio || same(*a.smallest_pivot_ratio, *b.smallest_pivot_ratio);
}

bool same(const std::optional<base::Error>& a, const std::optional<base::Error>& b) {
    if (a.has_value() != b.has_value()) return false;
    if (!a) return true;
    if (a->code != b->code || a->message != b->message) return false;
    if (a->context.has_value() != b->context.has_value()) return false;
    return !a->context || (a->context->index == b->context->index &&
                           a->context->value.has_value() == b->context->value.has_value() &&
                           (!a->context->value || same(*a->context->value, *b->context->value)));
}

bool same(const std::optional<results::ConvergenceRecord>& a,
          const std::optional<results::ConvergenceRecord>& b) {
    return a.has_value() == b.has_value() && (!a || same(*a, *b));
}

results::ConvergenceRecord convergence(double scale) {
    return {.iterations = {{1, 0.5 * scale, 1e3 * scale, 1.0}, {2, -0.0, tiny, 0.25}},
            .converged = true,
            .smallest_pivot_ratio = 1e-9 * scale};
}

results::NodeFields fields(double s, bool thermal) {
    return {.potential_V = {0.0, -0.0, s},
            .n_cm3 = {1e17, nan, tiny},
            .p_cm3 = {inf, -inf, 1e3},
            .temperature_K = thermal ? std::vector<double>{300.0, 301.5, 302.0}
                                     : std::vector<double>{}};
}

results::Provenance provenance() {
    return {.version = "0.0.0", .revision = "abc", .dirty = true, .build_type = "Release",
            .compiler = "MSVC 19.51.36260.0", .eigen_version = "5.0.1-dev",
            .eigen_cache_bytes = {49152, 3145728, 20971520}, .mkl_version = "",
            .os = "Windows 10.0.26300", .cpu = "cpu", .logical_cpus = 10,
            .started_utc = "2026-10-17T09:12:03.123Z", .wall_s = 1.5};
}

results::BiasPoint bias_point() {
    results::BiasPoint p;
    p.bias_V = {0.7, 0.0};
    p.fields = fields(0.7, true);
    p.bands = {{1.1, 1.0, 0.9}, {0.0, -0.1, -0.2}, {0.5, 0.5, nan}, {0.4, -0.0, 0.3}};
    p.terminal_current = {1e-3, -1e-3};
    p.terminal_current_resolution = {1e-15, 1e-15};
    p.gate_charge = {};
    p.interface_trap_charge = {2e-9};
    p.edge_current_n = {1.0, 2.0};
    p.edge_current_p = {-1.0, -0.0};
    p.convergence = convergence(1.0);
    p.thermal_bias_K = {300.0};
    p.thermal_contact_heat = {0.25};
    p.bias_substeps = 3;
    return p;
}

results::Transient transient(bool stopped) {
    results::Transient t;
    t.run = {.input_identity = 0x0123456789abcdefull,
             .settings = {{"transient.t_end_s", 1e-9}, {"transient.rtol", 1e-3}},
             .provenance = provenance()};
    for (int k = 0; k < 3; ++k) {
        results::TimePoint p{
            .time_s = k * 1e-10,
            .step_s = k == 0 ? 0.0 : 1e-10,
            .order = k,
            .error_ratio = k == 1 ? nan : 0.5,
            .bias_V = {0.1 * k, -0.0},
            .terminal_current = {1e-3 * k, -1e-3 * k},
            .conduction_current = {1e-3, -inf},
            .displacement_current = {tiny, 0.0},
            .contact_charge = {1e-12, -1e-12},
            .interface_trap_charge = k == 2 ? std::vector<double>{3e-9} : std::vector<double>{},
            .thermal_contact_heat = {0.5 * k},
            .convergence = k == 0 ? results::ConvergenceRecord{} : convergence(k)};
        t.points.push_back(std::move(p));
    }
    for (std::size_t k : {std::size_t{0}, std::size_t{2}}) {
        t.snapshots.push_back({.time_s = t.points[k].time_s,
                               .point = k,
                               .fields = fields(static_cast<double>(k), k == 2),
                               .bands = {{1.0}, {0.0}, {nan}, {-0.0}},
                               .edge_current_n = {1.0, -2.0},
                               .edge_current_p = {inf, 0.0},
                               .trap_occupancy = {0.25, 0.75}});
    }
    t.rejected_steps = 4;
    t.retraced_steps = 1;
    t.lagged_steps = 2;
    if (stopped) {
        t.stopped = base::Error{ErrorCode::non_convergence, "step below dt_min",
                                base::ErrorContext{.index = 7, .value = 1e-20}};
        t.unfinished = convergence(2.0);
    }
    return t;
}

void check(const results::Transient& a, const results::Transient& b) {
    CHECK(a.run.input_identity == b.run.input_identity);
    CHECK(a.run.settings == b.run.settings);
    CHECK(a.run.provenance == b.run.provenance);
    REQUIRE(a.points.size() == b.points.size());
    for (std::size_t k = 0; k < a.points.size(); ++k) {
        const auto& x = a.points[k];
        const auto& y = b.points[k];
        CHECK(same(x.time_s, y.time_s));
        CHECK(same(x.step_s, y.step_s));
        CHECK(x.order == y.order);
        CHECK(same(x.error_ratio, y.error_ratio));
        CHECK(same(x.bias_V, y.bias_V));
        CHECK(same(x.terminal_current, y.terminal_current));
        CHECK(same(x.conduction_current, y.conduction_current));
        CHECK(same(x.displacement_current, y.displacement_current));
        CHECK(same(x.contact_charge, y.contact_charge));
        CHECK(same(x.interface_trap_charge, y.interface_trap_charge));
        CHECK(same(x.thermal_contact_heat, y.thermal_contact_heat));
        CHECK(same(x.convergence, y.convergence));
    }
    REQUIRE(a.snapshots.size() == b.snapshots.size());
    for (std::size_t k = 0; k < a.snapshots.size(); ++k) {
        const auto& x = a.snapshots[k];
        const auto& y = b.snapshots[k];
        CHECK(same(x.time_s, y.time_s));
        CHECK(x.point == y.point);
        CHECK(same(x.fields, y.fields));
        CHECK(same(x.bands, y.bands));
        CHECK(same(x.edge_current_n, y.edge_current_n));
        CHECK(same(x.edge_current_p, y.edge_current_p));
        CHECK(same(x.trap_occupancy, y.trap_occupancy));
    }
    CHECK(a.rejected_steps == b.rejected_steps);
    CHECK(a.retraced_steps == b.retraced_steps);
    CHECK(a.lagged_steps == b.lagged_steps);
    CHECK(same(a.stopped, b.stopped));
    CHECK(same(a.unfinished, b.unfinished));
}

results::SmallSignal small_signal(bool stopped) {
    results::SmallSignal s;
    s.run = {.input_identity = 42, .settings = {{"small_signal.operating_point_equations", 1.0}},
             .provenance = provenance()};
    s.contacts = 2;
    s.frequency_Hz = {0.0, 1e6};
    for (int k = 0; k < 2; ++k) {
        results::SmallSignalPoint p;
        p.dc = bias_point();
        p.dc.bias_V = {0.1 * k, 0.0};
        p.admittance = {{{1e-3, 0.0}, {-1e-3, -0.0}, {-1e-3, nan}, {1e-3, inf}},
                        {{1e-3, 1e-9}, {-1e-3, -1e-9}, {-1e-3, -1e-9}, {tiny, 1e-9}}};
        p.pivot_ratio = {1e-8, nan};
        p.backward_error = {1e-16, 2e-16};
        p.resolution = {{1e-18, 2e-18}, {3e-18, -0.0}};
        p.fields.push_back({.frequency_Hz = 1e6,
                            .contact = 1,
                            .potential = {{1.0, 0.0}, {0.5, -0.25}},
                            .n_cm3 = {{1e10, nan}, {-0.0, 0.0}},
                            .p_cm3 = {{inf, 1.0}, {0.0, -inf}},
                            .temperature_K = k == 1 ? std::vector<std::complex<double>>{{1e-3, 2e-3}, {0.0, 0.0}}
                                                    : std::vector<std::complex<double>>{}});
        s.points.push_back(std::move(p));
    }
    if (stopped) {
        s.stopped = base::Error{ErrorCode::cancelled, "cancelled", std::nullopt};
        s.unfinished = convergence(3.0);
    }
    return s;
}

void check(const results::SmallSignal& a, const results::SmallSignal& b) {
    CHECK(a.run.input_identity == b.run.input_identity);
    CHECK(a.run.settings == b.run.settings);
    CHECK(a.run.provenance == b.run.provenance);
    CHECK(a.contacts == b.contacts);
    CHECK(same(a.frequency_Hz, b.frequency_Hz));
    REQUIRE(a.points.size() == b.points.size());
    for (std::size_t k = 0; k < a.points.size(); ++k) {
        const auto& x = a.points[k];
        const auto& y = b.points[k];
        CHECK(same(x.dc.bias_V, y.dc.bias_V));
        CHECK(same(x.dc.fields, y.dc.fields));
        CHECK(same(x.dc.bands, y.dc.bands));
        CHECK(same(x.dc.terminal_current, y.dc.terminal_current));
        CHECK(same(x.dc.convergence, y.dc.convergence));
        CHECK(x.dc.bias_substeps == y.dc.bias_substeps);
        REQUIRE(x.admittance.size() == y.admittance.size());
        for (std::size_t f = 0; f < x.admittance.size(); ++f) {
            CHECK(same(x.admittance[f], y.admittance[f]));
        }
        CHECK(same(x.pivot_ratio, y.pivot_ratio));
        CHECK(same(x.backward_error, y.backward_error));
        REQUIRE(x.resolution.size() == y.resolution.size());
        for (std::size_t f = 0; f < x.resolution.size(); ++f) {
            CHECK(same(x.resolution[f], y.resolution[f]));
        }
        REQUIRE(x.fields.size() == y.fields.size());
        for (std::size_t f = 0; f < x.fields.size(); ++f) {
            CHECK(same(x.fields[f].frequency_Hz, y.fields[f].frequency_Hz));
            CHECK(x.fields[f].contact == y.fields[f].contact);
            CHECK(same(x.fields[f].potential, y.fields[f].potential));
            CHECK(same(x.fields[f].n_cm3, y.fields[f].n_cm3));
            CHECK(same(x.fields[f].p_cm3, y.fields[f].p_cm3));
            CHECK(same(x.fields[f].temperature_K, y.fields[f].temperature_K));
        }
    }
    CHECK(same(a.stopped, b.stopped));
    CHECK(same(a.unfinished, b.unfinished));
}

std::vector<char> bytes_of(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

void write_bytes(const fs::path& p, const std::vector<char>& b) {
    std::ofstream f(p, std::ios::binary | std::ios::trunc);
    f.write(b.data(), static_cast<std::streamsize>(b.size()));
}

// A damaged copy of `good` is refused by `load` as invalid_input: truncated, one payload byte
// flipped (checksum), one byte appended.
template <class Load>
void check_damage(const fs::path& good, Load load) {
    const std::vector<char> b = bytes_of(good);
    const fs::path bad = good.parent_path() / "damaged";
    write_bytes(bad, {b.begin(), b.end() - 9});
    CHECK(load(bad).error().code == ErrorCode::invalid_input);
    std::vector<char> flipped = b;
    flipped[40] ^= 1;
    write_bytes(bad, flipped);
    CHECK(load(bad).error().code == ErrorCode::invalid_input);
    std::vector<char> longer = b;
    longer.push_back(0);
    write_bytes(bad, longer);
    CHECK(load(bad).error().code == ErrorCode::invalid_input);
}

}  // namespace

TEST_CASE("transient files: every field round-trips bit for bit") {
    const fs::path dir = scratch("transient_round_trip");
    for (const bool stopped : {false, true}) {
        CAPTURE(stopped);
        const results::Transient t = transient(stopped);
        REQUIRE(io::save_transient(dir / "t.nittran", t).has_value());
        const auto back = io::load_transient(dir / "t.nittran");
        REQUIRE(back.has_value());
        check(t, *back);
    }
    // An empty run (no points, no snapshots) too.
    REQUIRE(io::save_transient(dir / "e.nittran", results::Transient{}).has_value());
    const auto empty = io::load_transient(dir / "e.nittran");
    REQUIRE(empty.has_value());
    check(results::Transient{}, *empty);
}

TEST_CASE("small-signal files: every field round-trips bit for bit") {
    const fs::path dir = scratch("small_signal_round_trip");
    for (const bool stopped : {false, true}) {
        CAPTURE(stopped);
        const results::SmallSignal s = small_signal(stopped);
        REQUIRE(io::save_small_signal(dir / "s.nitac", s).has_value());
        const auto back = io::load_small_signal(dir / "s.nitac");
        REQUIRE(back.has_value());
        check(s, *back);
        CHECK(back->capacitance(1, 1, 0, 0) == s.capacitance(1, 1, 0, 0));
    }
}

TEST_CASE("transient and small-signal files: damage and other kinds are refused") {
    const fs::path dir = scratch("run_files_refused");
    REQUIRE(io::save_transient(dir / "t.nittran", transient(true)).has_value());
    REQUIRE(io::save_small_signal(dir / "s.nitac", small_signal(true)).has_value());
    check_damage(dir / "t.nittran", [](const fs::path& p) { return io::load_transient(p); });
    check_damage(dir / "s.nitac", [](const fs::path& p) { return io::load_small_signal(p); });

    // Each loader refuses the other kinds, and the sweep loader refuses both.
    CHECK(io::load_transient(dir / "s.nitac").error().code == ErrorCode::invalid_input);
    CHECK(io::load_small_signal(dir / "t.nittran").error().code == ErrorCode::invalid_input);
    CHECK(io::load_sweep(dir / "t.nittran").error().code == ErrorCode::invalid_input);
    CHECK(io::load_sweep(dir / "s.nitac").error().code == ErrorCode::invalid_input);
    results::Sweep sweep;
    REQUIRE(io::save_sweep(dir / "w.nitsweep", sweep).has_value());
    CHECK(io::load_transient(dir / "w.nitsweep").error().code == ErrorCode::invalid_input);
    CHECK(io::load_small_signal(dir / "w.nitsweep").error().code == ErrorCode::invalid_input);
}

TEST_CASE("small-signal files: admittances that do not fit the contacts are refused") {
    // A loaded run's admittance(point, frequency, i, j) must stay in bounds: every point has one
    // contacts x contacts matrix per frequency.
    const fs::path dir = scratch("small_signal_shape");
    results::SmallSignal s = small_signal(false);
    s.points[1].admittance[1].pop_back();
    REQUIRE(io::save_small_signal(dir / "s.nitac", s).has_value());
    CHECK(io::load_small_signal(dir / "s.nitac").error().code == ErrorCode::invalid_input);
    s = small_signal(false);
    s.points[1].admittance[1].push_back({1.0, 0.0});  // too many as well as too few
    REQUIRE(io::save_small_signal(dir / "s.nitac", s).has_value());
    CHECK(io::load_small_signal(dir / "s.nitac").error().code == ErrorCode::invalid_input);
    s = small_signal(false);
    s.points[0].admittance.pop_back();
    REQUIRE(io::save_small_signal(dir / "s.nitac", s).has_value());
    CHECK(io::load_small_signal(dir / "s.nitac").error().code == ErrorCode::invalid_input);
    s = small_signal(false);
    s.points[0].resolution[1].pop_back();
    REQUIRE(io::save_small_signal(dir / "s.nitac", s).has_value());
    CHECK(io::load_small_signal(dir / "s.nitac").error().code == ErrorCode::invalid_input);
    s = small_signal(false);
    s.points[1].fields[0].contact = 2;  // of two contacts, 0 and 1
    REQUIRE(io::save_small_signal(dir / "s.nitac", s).has_value());
    CHECK(io::load_small_signal(dir / "s.nitac").error().code == ErrorCode::invalid_input);
}

TEST_CASE("transient files: references that do not fit are refused") {
    const fs::path dir = scratch("transient_shape");
    results::Transient t = transient(false);
    t.snapshots.back().point = t.points.size();  // one past the last time point
    REQUIRE(io::save_transient(dir / "t.nittran", t).has_value());
    CHECK(io::load_transient(dir / "t.nittran").error().code == ErrorCode::invalid_input);
    t = transient(false);
    t.points[1].order = 3;  // backward Euler (1) or BDF2 (2), 0 for the start
    REQUIRE(io::save_transient(dir / "t.nittran", t).has_value());
    CHECK(io::load_transient(dir / "t.nittran").error().code == ErrorCode::invalid_input);
}
