// Semiconductor parameter validation and band quantities (ARCHITECTURE.md section 11, Unit 5).
// Reference values were computed independently at 40 digits from the legacy formulas and the
// CODATA 2018 constants (not from this code).
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "NiTCAD/physics/semiconductor.hpp"

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
    REQUIRE(close(ni, 10673775147.815624025, 1e-13));   // 40-digit reference
}

TEST_CASE("semiconductor: Varshni band gap and effective densities of states") {
    const Semiconductor si = silicon();
    REQUIRE(close(band_gap_eV(si, 300.0), 1.1245192307692307692, 1e-15));
    REQUIRE(close(band_gap_eV(si, 400.0), 1.0969498069498069498, 1e-15));
    REQUIRE(conduction_band_dos(si, 300.0) == 2.86e19);
    REQUIRE(valence_band_dos(si, 300.0) == 3.10e19);
    REQUIRE(close(conduction_band_dos(si, 400.0), 44032580530195458307.0, 1e-15));
    REQUIRE(close(intrinsic_density(si, 400.0), 5633682064707.5052902, 1e-13));
    REQUIRE(close(intrinsic_density(si, 250.0), 79127613.386517756152, 1e-13));
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
