// Masetti hole mobility (Unit 34, physics/mobility.hpp): against the formula and coefficients as
// recorded in DECISIONS.md (Unit 34), evaluated here independently (AC2), and its temperature slope
// (AC7).
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "NiTCAD/physics/mobility.hpp"

using namespace NiTCAD;

namespace {

// Masetti et al. (1983), eq (4) with Table I's boron values, as reproduced by Allpix Squared,
// Crosslight and Flooxs; mu_max scaled by (T / 300)^-2.2. Written out from DECISIONS.md, not
// from the implementation's parameter struct.
double reference(double N, double T) {
    const double mu_max = 470.5 * std::pow(T / 300.0, -2.2);
    if (N == 0.0) return mu_max;
    return 44.9 * std::exp(-9.23e16 / N) + mu_max / (1.0 + std::pow(N / 2.23e17, 0.719)) -
           29.0 / (1.0 + std::pow(6.10e20 / N, 2.0));
}

}  // namespace

TEST_CASE("masetti: the boron hole mobility as published (AC2)", "[physics][masetti]") {
    const auto& p = physics::masetti_boron_silicon;
    for (const double T : {250.0, 300.0, 400.0}) {
        CHECK(physics::masetti_mobility(p, 0.0, T) == 470.5 * std::pow(T / 300.0, -2.2));
        for (int e = 13; e <= 21; ++e) {
            const double N = std::pow(10.0, e);
            const double mu = physics::masetti_mobility(p, N, T);
            CAPTURE(T, N, mu, reference(N, T));
            CHECK(std::abs(mu / reference(N, T) - 1.0) < 1e-12);
            CHECK(mu > 0.0);
        }
    }
    // The low-doping limit is continuous: mu(N) -> mu_max as N -> 0 (at 1e6 cm^-3 the C_r term
    // is (1e6 / 2.23e17)^0.719, about 7e-9).
    CHECK(std::abs(physics::masetti_mobility(p, 1e6, 300.0) / 470.5 - 1.0) < 1e-8);
    // At the clustering density C_s the mobility stays positive.
    CHECK(physics::masetti_mobility(p, 6.10e20, 300.0) > 0.0);
}

TEST_CASE("masetti: the temperature slope (AC7)", "[physics][masetti]") {
    const auto& p = physics::masetti_boron_silicon;
    for (const double T : {250.0, 300.0, 400.0}) {
        for (const double N : {0.0, 1e15, 1e17, 1e19, 1e20}) {
            const double h = 1e-4 * T;  // truncation about 2e-8 of the slope
            const double fd = (physics::masetti_mobility(p, N, T + h) -
                               physics::masetti_mobility(p, N, T - h)) / (2.0 * h);
            const double slope = physics::masetti_mobility_slope(p, N, T);
            CAPTURE(T, N, slope, fd);
            CHECK(std::abs(slope / fd - 1.0) < 1e-6);
        }
    }
}

namespace {

// Masetti et al. (1983), eq (1) with Table I's phosphorus values (P_c = 0), as reproduced by Allpix
// Squared and Palankovski's thesis; mu_max scaled by (T / 300)^-2.5. Written out from DECISIONS.md
// (Unit 38), not from the parameter struct.
double electron_reference(double N, double T) {
    const double mu_max = 1414.0 * std::pow(T / 300.0, -2.5);
    if (N == 0.0) return mu_max;
    return 68.5 + (mu_max - 68.5) / (1.0 + std::pow(N / 9.20e16, 0.711)) -
           56.1 / (1.0 + std::pow(3.41e20 / N, 1.98));
}

}  // namespace

TEST_CASE("masetti: the phosphorus electron mobility as published (ME2)", "[physics][masetti]") {
    const auto& p = physics::masetti_phosphorus_silicon;
    for (const double T : {250.0, 300.0, 400.0}) {
        CAPTURE(T);
        CHECK(std::abs(physics::masetti_mobility(p, 0.0, T) /
                       (1414.0 * std::pow(T / 300.0, -2.5)) - 1.0) < 1e-12);
        for (int e = 13; e <= 21; ++e) {
            const double N = std::pow(10.0, e);
            CAPTURE(N);
            CHECK(std::abs(physics::masetti_mobility(p, N, T) / electron_reference(N, T) - 1.0) <
                  1e-12);
        }
        // Continuous at N -> 0 (E3: the mu_min1 term stays with P_c = 0).
        CHECK(std::abs(physics::masetti_mobility(p, 1e6, T) / physics::masetti_mobility(p, 0.0, T) -
                       1.0) < 1e-6);
        for (const double N : {0.0, 1e15, 1e17, 1e19, 1e20}) {
            const double h = 1e-4 * T;
            const double fd = (physics::masetti_mobility(p, N, T + h) -
                               physics::masetti_mobility(p, N, T - h)) / (2.0 * h);
            CHECK(std::abs(physics::masetti_mobility_slope(p, N, T) / fd - 1.0) < 1e-6);
        }
    }
}
