// Semiconductor parameter validation and band quantities (ARCHITECTURE.md section 11, Unit 5).
// References computed here in double-double arithmetic (references.hpp) from the legacy formulas
// and the CODATA 2018 constants (not from this code).
#include <catch2/catch_test_macros.hpp>

#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <functional>
#include <limits>
#include <new>
#include <string>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "references.hpp"

// Counts every allocation through the global operator new in this test executable, so a test can
// show that successful validation does not allocate (base/error.hpp: strings only on the error
// path).
namespace {
std::atomic<long long> allocation_count{0};
}  // namespace

void* operator new(std::size_t size) {
    ++allocation_count;
    if (void* p = std::malloc(size == 0 ? 1 : size)) return p;
    throw std::bad_alloc{};
}
void operator delete(void* p) noexcept { std::free(p); }
void operator delete(void* p, std::size_t) noexcept { std::free(p); }

using namespace NiTCAD::physics;
using NiTCAD::base::ErrorCode;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

constexpr double not_a_number = std::numeric_limits<double>::quiet_NaN();
constexpr double infinity = std::numeric_limits<double>::infinity();

}  // namespace

TEST_CASE("semiconductor: the silicon parameter set is valid") {
    REQUIRE(Semiconductor::create(silicon_parameters).has_value());
    REQUIRE(silicon().parameters().eps_r == 11.7);
}

TEST_CASE("semiconductor: n_i(300 K) of silicon is 1.0674e10 cm^-3 (section 10 gate)") {
    const Semiconductor si = silicon();
    const double ni = intrinsic_density(si, 300.0);
    REQUIRE(close(ni, 1.0674e10, 1e-4));                // proposed gate (ARCHITECTURE.md section 10)
    REQUIRE((ni > 9e9 && ni < 1.6e10));                 // legacy test_ni_300k_within_accepted_band
    REQUIRE(close(ni, reference::intrinsic_density(silicon_parameters, 300.0).value(), 1e-13));
}

TEST_CASE("semiconductor: Varshni band gap and effective densities of states") {
    const Semiconductor si = silicon();
    const SemiconductorParameters& p = silicon_parameters;
    REQUIRE(close(band_gap_eV(si, 300.0), reference::band_gap(p, 300.0).value(), 1e-15));
    REQUIRE(close(band_gap_eV(si, 400.0), reference::band_gap(p, 400.0).value(), 1e-15));
    REQUIRE(conduction_band_dos(si, 300.0) == 2.86e19);
    REQUIRE(valence_band_dos(si, 300.0) == 3.10e19);
    REQUIRE(close(conduction_band_dos(si, 400.0),
                  (reference::DD(p.Nc300) * reference::dos_factor(400.0)).value(), 1e-15));
    REQUIRE(close(intrinsic_density(si, 400.0), reference::intrinsic_density(p, 400.0).value(),
                  1e-13));
    REQUIRE(close(intrinsic_density(si, 250.0), reference::intrinsic_density(p, 250.0).value(),
                  1e-13));
}

TEST_CASE("semiconductor: every parameter is validated and the error names it") {
    using P = SemiconductorParameters;
    struct Case {
        std::string name;
        std::function<void(P&)> spoil;
    };
    const std::vector<Case> cases{
        {"eps_r", [](P& p) { p.eps_r = 0.0; }},
        {"Eg0_eV", [](P& p) { p.Eg0_eV = -1.0; }},
        {"varshni_alpha_eV_per_K", [](P& p) { p.varshni_alpha_eV_per_K = -1e-4; }},
        {"varshni_beta_K", [](P& p) { p.varshni_beta_K = not_a_number; }},
        {"Nc300", [](P& p) { p.Nc300 = 0.0; }},
        {"Nv300", [](P& p) { p.Nv300 = infinity; }},
        {"electron_affinity_eV", [](P& p) { p.electron_affinity_eV = -0.1; }},
        {"electron_saturation.v_sat_cm_s",
         [](P& p) { p.electron_saturation.v_sat_cm_s = 0.0; }},
        {"hole_saturation.beta", [](P& p) { p.hole_saturation.beta = 0.9; }},
        {"electron_mobility.mu_min", [](P& p) { p.electron_mobility.mu_min = -1.0; }},
        {"electron_mobility.mu_max", [](P& p) { p.electron_mobility.mu_max = 0.0; }},
        {"electron_mobility.N_ref", [](P& p) { p.electron_mobility.N_ref = 0.0; }},
        {"electron_mobility.alpha", [](P& p) { p.electron_mobility.alpha = 0.0; }},
        {"electron_mobility.T_exponent",
         [](P& p) { p.electron_mobility.T_exponent = not_a_number; }},
        {"hole_mobility.mu_min", [](P& p) { p.hole_mobility.mu_min = 600.0; }},  // > mu_max
        {"hole_mobility.N_ref", [](P& p) { p.hole_mobility.N_ref = -infinity; }},
        {"lifetime.tau_n0", [](P& p) { p.lifetime.tau_n0 = 0.0; }},
        {"lifetime.tau_p0", [](P& p) { p.lifetime.tau_p0 = -3e-6; }},
        {"lifetime.N_ref", [](P& p) { p.lifetime.N_ref = not_a_number; }},
    };
    for (const Case& c : cases) {
        CAPTURE(c.name);
        SemiconductorParameters p = silicon_parameters;
        c.spoil(p);
        const auto m = Semiconductor::create(p);
        REQUIRE_FALSE(m.has_value());
        REQUIRE(m.error().code == ErrorCode::invalid_input);
        REQUIRE(m.error().message.find(c.name) != std::string::npos);
        REQUIRE(m.error().context.has_value());
    }
}

TEST_CASE("semiconductor: zero values allowed by the legacy material sets are accepted") {
    SemiconductorParameters p = silicon_parameters;
    p.electron_mobility.mu_min = 0.0;  // legacy GE, GAAS
    p.hole_mobility.mu_min = 0.0;
    p.varshni_beta_K = 0.0;            // legacy SIC_4H
    REQUIRE(Semiconductor::create(p).has_value());
}

TEST_CASE("semiconductor: successful validation does not allocate") {
    const long long before = allocation_count.load();
    const auto m = Semiconductor::create(silicon_parameters);
    const bool valid_at_300 = m.has_value() && check_temperature(*m, 300.0).has_value();
    const long long after = allocation_count.load();
    REQUIRE(valid_at_300);
    REQUIRE(after == before);
    // The counter does see the library's allocations: a failure builds its message.
    SemiconductorParameters p = silicon_parameters;
    p.hole_mobility.alpha = 0.0;
    const long long before_failure = allocation_count.load();
    const auto bad = Semiconductor::create(p);
    const long long after_failure = allocation_count.load();
    REQUIRE_FALSE(bad.has_value());
    REQUIRE(bad.error().message == "semiconductor parameter hole_mobility.alpha must be finite "
                                   "and positive");
    REQUIRE(after_failure > before_failure);
}

TEST_CASE("semiconductor: n_i uses the base thermal voltage for kT") {
    // At these temperatures k_B_eV_per_K * T and thermal_voltage(T) = k_B T / q differ by one
    // rounding, and Eg / 2kT is 300 to 70, so the other spelling would move n_i by about 1e-14
    // relative. With the same kT, only the prefactor's roundings differ (a few 1e-16).
    const Semiconductor si = silicon();
    for (const double T : {22.0, 30.5, 44.0, 77.5, 95.5}) {
        CAPTURE(T);
        const double expected =
            std::sqrt(conduction_band_dos(si, T) * valence_band_dos(si, T)) *
            std::exp(-band_gap_eV(si, T) / (2.0 * NiTCAD::base::thermal_voltage(T)));
        REQUIRE(close(intrinsic_density(si, T), expected, 2e-15));
    }
}

TEST_CASE("check_temperature: silicon is accepted over its working range") {
    const Semiconductor si = silicon();
    for (const double T : {1.0, 77.0, 300.0, 500.0, 850.0}) {
        CAPTURE(T);
        REQUIRE(check_temperature(si, T).has_value());
    }
}

TEST_CASE("check_temperature: T must be finite and positive") {
    const Semiconductor si = silicon();
    for (const double T : {0.0, -0.0, -300.0, not_a_number, infinity, -infinity}) {
        CAPTURE(T);
        const auto r = check_temperature(si, T);
        REQUIRE_FALSE(r.has_value());
        REQUIRE(r.error().code == ErrorCode::invalid_input);
        REQUIRE(r.error().message == "temperature must be finite and positive");
    }
}

TEST_CASE("check_temperature: the band gap must stay positive") {
    // Silicon's Varshni gap crosses zero near 2998 K.
    const Semiconductor si = silicon();
    const auto hot = check_temperature(si, 3500.0);
    REQUIRE_FALSE(hot.has_value());
    REQUIRE(hot.error().message.find("band gap") != std::string::npos);
    REQUIRE(hot.error().context->value == band_gap_eV(si, 3500.0));
    REQUIRE(*hot.error().context->value < 0.0);
    // Exactly zero is rejected: beta = 0 (as legacy SIC_4H) and alpha = 2^-10 eV/K give
    // Eg = 1 - T / 1024 eV, exact in floating point. Temperature-independent mobility isolates
    // the band-gap rule.
    SemiconductorParameters p = silicon_parameters;
    p.Eg0_eV = 1.0;
    p.varshni_alpha_eV_per_K = 0x1p-10;
    p.varshni_beta_K = 0.0;
    p.electron_mobility.T_exponent = 0.0;
    p.hole_mobility.T_exponent = 0.0;
    const Semiconductor m = *Semiconductor::create(p);
    REQUIRE(band_gap_eV(m, 1024.0) == 0.0);
    REQUIRE_FALSE(check_temperature(m, 1024.0).has_value());
    REQUIRE(check_temperature(m, 1023.0).has_value());
}

TEST_CASE("check_temperature: mu_max(T) must not fall below mu_min") {
    const Semiconductor si = silicon();
    // Holes cross at about 857 K, electrons at about 953 K; electrons are reported first.
    REQUIRE(check_temperature(si, 850.0).has_value());
    const auto holes = check_temperature(si, 860.0);
    REQUIRE_FALSE(holes.has_value());
    REQUIRE(holes.error().message.starts_with("hole_mobility"));
    REQUIRE(*holes.error().context->value < si.parameters().hole_mobility.mu_min);
    const auto electrons = check_temperature(si, 1000.0);
    REQUIRE_FALSE(electrons.has_value());
    REQUIRE(electrons.error().message.starts_with("electron_mobility"));
    // What the rule prevents: past the crossing, mobility rises with doping.
    REQUIRE(caughey_thomas_mobility(si, Carrier::electron, 0.0, 1000.0) <
            caughey_thomas_mobility(si, Carrier::electron, 1e20, 1000.0));
    REQUIRE(caughey_thomas_mobility(si, Carrier::hole, 0.0, 850.0) >
            caughey_thomas_mobility(si, Carrier::hole, 1e20, 850.0));
    // mu_max(T) overflowing to infinity at a tiny T is rejected too.
    const auto tiny = check_temperature(si, 1e-300);
    REQUIRE_FALSE(tiny.has_value());
    REQUIRE(tiny.error().message.starts_with("electron_mobility"));
}
