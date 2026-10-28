// MLDA quantum confinement (Unit 47; DECISIONS.md, Unit 47, QC2): the thermal wavelength from the
// constants, the factor against a quadrature of MLDA's nondegenerate j0-kernel integral (not of
// the closed form it reduces to), and the half-cell mean against a quadrature of the factor.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <numbers>

#include "NiTCAD/physics/quantum_confinement.hpp"

using namespace NiTCAD::physics;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

// 1 - sin(y) / y, by its series where the quotient cancels.
long double one_minus_j0(long double y) {
    if (y < 0.5L) {
        const long double y2 = y * y;
        long double term = y2 / 6.0L, sum = 0.0L;  // y^(2k) / (2k + 1)!, alternating
        for (int k = 1; k <= 12; ++k) {
            sum += term;
            term *= -y2 / ((2.0L * k + 2.0L) * (2.0L * k + 3.0L));
        }
        return sum;
    }
    return 1.0L - std::sin(y) / y;
}

// int_0^inf sqrt(e) e^(-e) [1 - j0(2 a sqrt(e))] de / (sqrt(pi) / 2), with e = u^2:
// int_0^inf 2 u^2 e^(-u^2) [1 - j0(2 a u)] du, composite Simpson on [0, 12].
double mlda_integral(double a) {
    constexpr int n = 200000;
    const long double h = 12.0L / n;
    auto g = [a](long double u) {
        return 2.0L * u * u * std::exp(-u * u) * one_minus_j0(2.0L * a * u);
    };
    long double s = g(0.0L) + g(12.0L);
    for (int i = 1; i < n; ++i) s += (i % 2 ? 4.0L : 2.0L) * g(i * h);
    return static_cast<double>(s * h / 3.0L / (std::sqrt(std::numbers::pi_v<long double>) / 2.0L));
}

// (1 / a) int_0^a (1 - e^(-z^2 / lambda^2)) dz by composite Simpson.
double mean_quadrature(double a, double lambda) {
    constexpr int n = 20000;
    const long double h = static_cast<long double>(a) / n;
    auto g = [lambda](long double z) {
        const long double x = z / lambda;
        return -std::expm1(-x * x);
    };
    long double s = g(0.0L) + g(a);
    for (int i = 1; i < n; ++i) s += (i % 2 ? 4.0L : 2.0L) * g(i * h);
    return static_cast<double>(s * h / 3.0L / a);
}

}  // namespace

TEST_CASE("quantum confinement: silicon's quantization masses (Unit 47)") {
    REQUIRE(mlda_silicon_electron_mass == 0.98);
    REQUIRE(mlda_silicon_hole_mass == 0.49);
}

TEST_CASE("quantum confinement: the thermal wavelength (QC2)") {
    // hbar / sqrt(2 m k T) from the CODATA 2018 values written out here.
    for (double m : {0.49, 0.98, 1.0, 0.19}) {
        for (double T : {250.0, 300.0, 400.0}) {
            const long double ref = 1.054571817e-34L /
                                    std::sqrt(2.0L * m * 9.1093837015e-31L * 1.380649e-23L * T) *
                                    100.0L;
            CHECK(close(mlda_wavelength_cm(m, T), static_cast<double>(ref), 1e-14));
        }
    }
    CHECK(std::abs(mlda_wavelength_cm(mlda_silicon_electron_mass, 300.0) - 1.2263e-7) <= 1e-11);
    CHECK(std::abs(mlda_wavelength_cm(mlda_silicon_hole_mass, 300.0) - 1.7343e-7) <= 1e-11);
    CHECK(std::abs(mlda_wavelength_cm(1.0, 300.0) - 1.2140e-7) <= 1e-11);
}

TEST_CASE("quantum confinement: the factor is MLDA's nondegenerate integral (QC2)") {
    const double lambda = 1.2263e-7;
    for (double a : {0.01, 0.03, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0}) {
        const double f = mlda_factor(a * lambda, lambda);
        CHECK(close(f, mlda_integral(a), 1e-10));
    }
    CHECK(mlda_factor(0.0, lambda) == 0.0);
    CHECK(mlda_factor(1e-4, lambda) == 1.0);  // 1 um: exp(-(d / lambda)^2) underflows
}

TEST_CASE("quantum confinement: the half-cell mean (QC2)") {
    const double lambda = 1.7343e-7;
    for (double x : {1e-4, 1e-2, 0.1, 0.3, 0.49, 0.5, 0.51, 1.0, 2.0, 5.0, 20.0}) {
        const double a = x * lambda;
        CHECK(close(mlda_mean_factor(a, lambda), mean_quadrature(a, lambda), 1e-12));
        CHECK(mlda_mean_factor(a, lambda) > 0.0);
        // Never above the factor at the far end of the interval (f increases).
        CHECK(mlda_mean_factor(a, lambda) <= mlda_factor(a, lambda));
    }
    // Small depth: x^2 / 3; large depth: 1.
    CHECK(close(mlda_mean_factor(1e-6 * lambda, lambda), 1e-12 / 3.0, 1e-9));
    CHECK(close(mlda_mean_factor(1e3 * lambda, lambda), 1.0 - std::sqrt(std::numbers::pi) / 2e3,
                1e-14));
}
