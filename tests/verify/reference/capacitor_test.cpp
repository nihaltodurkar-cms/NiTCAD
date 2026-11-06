// L0 for the capacitor references (DECISIONS Unit 55 16.3): each closed form against a second method, a
// Gauss-Kronrod integral of the elastance density in double-double, and the limits it must meet.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "NiTCAD/verify/reference/capacitor.hpp"
#include "NiTCAD/verify/reference/gauss_kronrod.hpp"

using namespace NiTCAD::verify::reference;

namespace {

double rel(const DD& a, const DD& b) { return std::abs(((a - b) / b).hi); }

}  // namespace

TEST_CASE("capacitor: eps0 is CODATA 2018's printed value in F/cm", "[capacitor]") {
    const DD e0 = vacuum_permittivity_F_cm();
    REQUIRE(e0.hi == 8.8541878128e-14);
    // The ten printed digits are exact in the double-double (88541878128 < 2^53; 10^22 is the largest exact power of ten).
    REQUIRE(rel(e0 * DD{1e22, 0.0}, DD{88541878128.0, 0.0} / DD{100.0, 0.0}) <= 4.0 * u_dd);
}

TEST_CASE("capacitor: the series stack against the integral of dx / eps(x)", "[capacitor]") {
    const std::vector<double> t{3e-7, 1.7e-7, 5e-7};
    const std::vector<double> e{3.9, 22.0, 7.5};
    const DD C = series_capacitance(t, e);
    // Second method: elastance 1/C = (1/eps0) int_0^T dx / eps_r(x), one Kronrod integral per layer.
    DD inverse{0.0, 0.0};
    for (std::size_t i = 0; i < t.size(); ++i) {
        const double er = e[i];
        const auto q = qag15([er](const DD&) { return DD{1.0, 0.0} / DD{er, 0.0}; }, DD{0.0, 0.0}, DD{t[i], 0.0},
                             0.0, 0.5e-28);
        REQUIRE(q.status == QuadratureStatus::ok);
        inverse += q.value;
    }
    REQUIRE(rel(C, vacuum_permittivity_F_cm() / inverse) <= 16.0 * u_dd);
    // One layer is eps0 eps_r / t; splitting a layer leaves C unchanged.
    REQUIRE(rel(series_capacitance(std::vector<double>{4e-7}, std::vector<double>{3.9}),
                vacuum_permittivity_F_cm() * DD{3.9, 0.0} / DD{4e-7, 0.0}) <= 4.0 * u_dd);
    REQUIRE(rel(series_capacitance(std::vector<double>{1e-7, 3e-7}, std::vector<double>{3.9, 3.9}),
                series_capacitance(std::vector<double>{4e-7}, std::vector<double>{3.9})) <= 8.0 * u_dd);
    REQUIRE_THROWS_AS(series_capacitance(std::vector<double>{}, std::vector<double>{}), std::invalid_argument);
    REQUIRE_THROWS_AS(series_capacitance(std::vector<double>{1e-7}, std::vector<double>{-1.0}),
                      std::invalid_argument);
}

TEST_CASE("capacitor: the coaxial form against the integral of dr / (2 pi eps r)", "[capacitor]") {
    for (const double b : {2e-5, 1.1e-5, 1e-3}) {
        const double a = 1e-5, er = 3.9;
        const DD C = coaxial_capacitance(a, b, er);
        const auto q = qag15([](const DD& r) { return DD{1.0, 0.0} / r; }, DD{a, 0.0}, DD{b, 0.0}, 0.0, 0.5e-28);
        REQUIRE(q.status == QuadratureStatus::ok);
        const DD second = ldexp(pi(), 1) * vacuum_permittivity_F_cm() * DD{er, 0.0} / q.value;
        CAPTURE(b);
        REQUIRE(rel(C, second) <= 64.0 * u_dd);
    }
    // Thin shell limit: b = a (1 + d), C -> 2 pi a eps / (a d) = the planar eps / d per unit area times 2 pi a.
    const double a = 1e-4, d = 1e-6;
    const DD thin = coaxial_capacitance(a, a * (1.0 + d), 3.9);
    const DD planar = series_capacitance(std::vector<double>{a * d}, std::vector<double>{3.9}) * ldexp(pi(), 1) *
                      DD{a, 0.0};
    REQUIRE(rel(thin, planar) <= 1e-6);  // the next term is d / 2
    REQUIRE_THROWS_AS(coaxial_capacitance(2e-5, 1e-5, 3.9), std::invalid_argument);
}
