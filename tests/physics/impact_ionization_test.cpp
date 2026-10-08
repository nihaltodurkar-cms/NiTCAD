// Impact ionization coefficients (ARCHITECTURE.md section 11, Unit 19): the van Overstraeten-de
// Man silicon values and branches (the legacy's, with its hole-switch fix), the temperature factor,
// the exact field derivative against finite differences, the low-field cut-off and the parameter
// validation.
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/physics/impact_ionization.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD::physics;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

}  // namespace

TEST_CASE("impact ionization: silicon's van Overstraeten-de Man coefficients") {
    const Semiconductor si = silicon();
    const auto& n = impact_ionization(si, Carrier::electron);
    const auto& p = impact_ionization(si, Carrier::hole);
    REQUIRE(n == ImpactIonizationCoefficients{7.03e5, 1.231e6, 7.03e5, 1.231e6, 5.0e5});
    REQUIRE(p == ImpactIonizationCoefficients{1.582e6, 2.036e6, 6.71e5, 1.693e6, 4.0e5});
    REQUIRE(si.parameters().impact_ionization.phonon_energy_eV == 0.063);
    // A exp(-B / E) on the branch in use: electrons one fit; holes low below 4e5 V/cm, high from
    // it on (the legacy shared electron switch 5e5 put holes in [4e5, 5e5) on the low branch).
    for (const double E : {1.75e5, 3e5, 3.99e5, 4e5, 4.5e5, 6e5}) {
        CAPTURE(E);
        REQUIRE(close(impact_ionization_coefficient(n, 1.0, E).alpha,
                      7.03e5 * std::exp(-1.231e6 / E), 1e-15));
        const double hole = E < 4e5 ? 1.582e6 * std::exp(-2.036e6 / E)
                                    : 6.71e5 * std::exp(-1.693e6 / E);
        REQUIRE(close(impact_ionization_coefficient(p, 1.0, E).alpha, hole, 1e-15));
    }
    // Magnitudes from the paper's fits: alpha_n(3e5) = 1.161e4 /cm, alpha_p(4.5e5) = 1.559e4 /cm.
    REQUIRE(close(impact_ionization_coefficient(n, 1.0, 3e5).alpha, 1.161e4, 1e-3));
    REQUIRE(close(impact_ionization_coefficient(p, 1.0, 4.5e5).alpha, 1.559e4, 1e-3));
    // The hole branches meet in value at the switch (to the fits' 0.1%), not in slope.
    const double below = impact_ionization_coefficient(p, 1.0, 4e5 * (1.0 - 1e-12)).alpha;
    const double at = impact_ionization_coefficient(p, 1.0, 4e5).alpha;
    REQUIRE(close(below, at, 2e-3));
}

TEST_CASE("impact ionization: the derivative is exact") {
    const Semiconductor si = silicon();
    for (const Carrier c : {Carrier::electron, Carrier::hole}) {
        const auto& k = impact_ionization(si, c);
        for (const double gamma : {1.0, 0.8}) {
            for (const double E : {5e4, 1.5e5, 3e5, 3.7e5, 4.4e5, 7e5, 2e6}) {
                const double h = 1e-6 * E;
                const double fd = (impact_ionization_coefficient(k, gamma, E + h).alpha -
                                   impact_ionization_coefficient(k, gamma, E - h).alpha) /
                                  (2.0 * h);
                const double d = impact_ionization_coefficient(k, gamma, E).d_dE;
                CAPTURE(E, gamma, d, fd);
                REQUIRE(close(fd, d, 1e-7));
            }
        }
    }
}

TEST_CASE("impact ionization: the temperature factor") {
    // gamma = tanh(hw / 2kT0) / tanh(hw / 2kT): 1 at 300 K, above 1 above it (in the exponent it
    // dominates: the coefficients fall with temperature), below 1 below it; 1 at every T without a
    // phonon energy.
    REQUIRE(impact_ionization_temperature_factor(0.063, 300.0) == 1.0);
    REQUIRE(impact_ionization_temperature_factor(0.0, 450.0) == 1.0);
    const double k = NiTCAD::base::k_B_eV_per_K;
    for (const double T : {200.0, 350.0, 500.0}) {
        const double g = impact_ionization_temperature_factor(0.063, T);
        CAPTURE(T, g);
        REQUIRE(close(g, std::tanh(0.063 / (2.0 * k * 300.0)) / std::tanh(0.063 / (2.0 * k * T)),
                      1e-15));
        REQUIRE((T > 300.0 ? g > 1.0 : g < 1.0));
    }
    const auto& n = impact_ionization(silicon(), Carrier::electron);
    const double g = impact_ionization_temperature_factor(0.063, 400.0);
    REQUIRE(impact_ionization_coefficient(n, g, 3e5).alpha <
            impact_ionization_coefficient(n, 1.0, 3e5).alpha);
    REQUIRE(close(impact_ionization_coefficient(n, g, 3e5).alpha,
                  g * 7.03e5 * std::exp(-g * 1.231e6 / 3e5), 1e-15));
}

TEST_CASE("impact ionization: zero at zero field, without coefficients, and below underflow") {
    const auto& n = impact_ionization(silicon(), Carrier::electron);
    for (const double E : {0.0, 1e-300, 1.0, 1.7e3}) {
        CAPTURE(E);
        const ImpactIonizationRate r = impact_ionization_coefficient(n, 1.0, E);
        REQUIRE(r.alpha == 0.0);
        REQUIRE(r.d_dE == 0.0);
    }
    REQUIRE(impact_ionization_coefficient(n, 1.0, 1.8e3).alpha > 0.0);  // just above the cut
    REQUIRE(std::isfinite(impact_ionization_coefficient(n, 1.0, 1.8e3).d_dE));
    const ImpactIonizationCoefficients none{};
    REQUIRE(impact_ionization_coefficient(none, 1.0, 1e6).alpha == 0.0);
    // Other materials carry none.
    REQUIRE(gallium_arsenide_parameters.impact_ionization == ImpactIonizationParameters{});
}

TEST_CASE("impact ionization: parameter validation") {
    using NiTCAD::base::ErrorCode;
    const auto bad = [](auto edit) {
        SemiconductorParameters p = silicon_parameters;
        edit(p);
        const auto m = Semiconductor::create(p);
        REQUIRE_FALSE(m.has_value());
        return m.error();
    };
    REQUIRE(bad([](auto& p) { p.impact_ionization.electron.A_low_per_cm = -1.0; }).code ==
            ErrorCode::invalid_input);
    REQUIRE(bad([](auto& p) { p.impact_ionization.hole.B_high_V_per_cm = 0.0; }).context->value ==
            0.0);
    REQUIRE(bad([](auto& p) { p.impact_ionization.hole.switch_V_per_cm = 0.0; }).code ==
            ErrorCode::invalid_input);
    REQUIRE(bad([](auto& p) { p.impact_ionization.phonon_energy_eV = -0.01; }).code ==
            ErrorCode::invalid_input);
    REQUIRE(bad([](auto& p) { p.impact_ionization.electron.B_low_V_per_cm = std::nan(""); })
                .code == ErrorCode::invalid_input);
    REQUIRE(bad([](auto& p) { p.impact_ionization.electron.B_low_V_per_cm = 0.0; }).code ==
            ErrorCode::invalid_input);
    REQUIRE(bad([](auto& p) { p.impact_ionization.hole.A_high_per_cm = -1.0; }).code ==
            ErrorCode::invalid_input);
    // A branch without A may leave B at zero.
    SemiconductorParameters q = silicon_parameters;
    q.impact_ionization.hole.A_low_per_cm = 0.0;
    q.impact_ionization.hole.B_low_V_per_cm = 0.0;
    REQUIRE(Semiconductor::create(q).has_value());
    // A carrier without coefficients is valid.
    SemiconductorParameters p = silicon_parameters;
    p.impact_ionization.hole = {};
    REQUIRE(Semiconductor::create(p).has_value());
}
