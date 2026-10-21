// Rates and the mesh digest in files (Unit 45, RG4): sweep, transient and small-signal files carry
// both bit for bit; a point or snapshot without rates stays without; files written before Unit 45
// load with no rates and digest 0.
#include <catch2/catch_test_macros.hpp>

#include <bit>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <string>
#include <vector>

#include "NiTCAD/io/files.hpp"

using namespace NiTCAD;
namespace fs = std::filesystem;

namespace {

fs::path scratch(const std::string& name) {
    const fs::path dir = fs::temp_directory_path() / "nitcad_io_test" / name;
    fs::remove_all(dir);
    fs::create_directories(dir);
    return dir;
}

bool same(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (std::bit_cast<std::uint64_t>(a[i]) != std::bit_cast<std::uint64_t>(b[i])) return false;
    }
    return true;
}

bool same(const results::RateFields& a, const results::RateFields& b) {
    return same(a.srh, b.srh) && same(a.auger, b.auger) && same(a.radiative, b.radiative) &&
           same(a.impact, b.impact) && same(a.interface_n, b.interface_n) &&
           same(a.interface_p, b.interface_p) && same(a.band_to_band_n, b.band_to_band_n) &&
           same(a.band_to_band_p, b.band_to_band_p) && same(a.net_n, b.net_n) && same(a.net_p, b.net_p);
}

// Every vector set, with values that test the bits (signed zero, NaN, a denormal, huge and tiny).
results::RateFields full() {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double tiny = std::numeric_limits<double>::denorm_min();
    results::RateFields r;
    r.srh = {1e20, -0.0, 3.5e9};
    r.auger = {1e3, 2e3, 3e3};
    r.radiative = {nan, 1.0, 2.0};
    r.impact = {-1e18, 0.0, -tiny};
    r.interface_n = {0.0, 7e15, 0.0};
    r.interface_p = {0.0, 6.9e15, 0.0};
    r.band_to_band_n = {-1e22, 0.0, 0.0};
    r.band_to_band_p = {0.0, 0.0, -1e22};
    r.net_n = {1, 2, 3};
    r.net_p = {4, 5, 6};
    return r;
}

results::BiasPoint point(bool rates) {
    results::BiasPoint p;
    p.bias_V = {0.1, 0.0};
    p.fields = {.potential_V = {0.0, 0.1, 0.2}, .n_cm3 = {1, 2, 3}, .p_cm3 = {4, 5, 6}};
    if (rates) p.rates = full();
    return p;
}

const fs::path pre_unit42 = fs::path(NITCAD_IO_FIXTURES) / "pre_unit42";

}  // namespace

TEST_CASE("rates files: sweeps carry the rates and the mesh digest (Unit 45, RG4)") {
    const fs::path dir = scratch("rates_sweep");
    results::Sweep s;
    s.run.input_identity = 11;
    s.run.mesh_digest = 0x0123456789abcdefULL;
    s.points = {point(true), point(false)};
    REQUIRE(io::save_sweep(dir / "s.nitsweep", s).has_value());
    const auto back = io::load_sweep(dir / "s.nitsweep");
    REQUIRE(back.has_value());
    CHECK(back->run.mesh_digest == 0x0123456789abcdefULL);
    CHECK(same(back->points[0].rates, full()));
    CHECK(same(back->points[1].rates, results::RateFields{}));
    // Digest 0 (unknown) is not written and reads as 0.
    s.run.mesh_digest = 0;
    REQUIRE(io::save_sweep(dir / "z.nitsweep", s).has_value());
    CHECK(io::load_sweep(dir / "z.nitsweep")->run.mesh_digest == 0);
}

TEST_CASE("rates files: transient snapshots and small-signal operating points (Unit 45, RG4)") {
    const fs::path dir = scratch("rates_runs");
    results::Transient t;
    t.run.mesh_digest = 42;
    t.points.push_back({.time_s = 0.0, .step_s = 0.0, .order = 0, .error_ratio = 0.0, .bias_V = {0.0, 0.0},
                        .terminal_current = {0.0, 0.0}, .conduction_current = {0.0, 0.0},
                        .displacement_current = {0.0, 0.0}, .contact_charge = {0.0, 0.0},
                        .interface_trap_charge = {}, .thermal_contact_heat = {}, .convergence = {}});
    results::TransientSnapshot snap;
    snap.point = 0;
    snap.fields = point(false).fields;
    snap.rates = full();
    t.snapshots = {snap, snap};
    t.snapshots[1].rates = {};
    REQUIRE(io::save_transient(dir / "t.nittran", t).has_value());
    const auto tb = io::load_transient(dir / "t.nittran");
    REQUIRE(tb.has_value());
    CHECK(tb->run.mesh_digest == 42);
    CHECK(same(tb->snapshots[0].rates, full()));
    CHECK(same(tb->snapshots[1].rates, results::RateFields{}));

    results::SmallSignal a;
    a.run.mesh_digest = 43;
    a.contacts = 2;
    a.frequency_Hz = {1e6};
    results::SmallSignalPoint ap;
    ap.dc = point(true);
    ap.admittance = {{{1.0, 2.0}, {0.0, 0.0}, {0.0, 0.0}, {1.0, 2.0}}};
    ap.pivot_ratio = {0.1};
    ap.backward_error = {1e-16};
    a.points = {ap};
    REQUIRE(io::save_small_signal(dir / "a.nitac", a).has_value());
    const auto ab = io::load_small_signal(dir / "a.nitac");
    REQUIRE(ab.has_value());
    CHECK(ab->run.mesh_digest == 43);
    CHECK(same(ab->points[0].dc.rates, full()));
}

TEST_CASE("rates files: files written before Unit 45 load with no rates and digest 0 (RG4)") {
    const auto sweep = io::load_sweep(pre_unit42 / "bar.nitsweep");
    REQUIRE(sweep.has_value());
    CHECK(sweep->run.mesh_digest == 0);
    for (const results::BiasPoint& p : sweep->points) CHECK(same(p.rates, results::RateFields{}));
}
