// Lombardi surface mobility (Unit 46; DECISIONS.md, Unit 46, LS2): the formula and coefficients
// evaluated independently here (long double, the terms summed as reciprocals in another order),
// and the slopes against central differences.
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "NiTCAD/physics/surface_mobility.hpp"

using namespace NiTCAD::physics;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

// 1 / mu = 1 / mu_b + E / (B + C N^lambda E^(2/3) / T) + E^2 / delta, the reference's own form.
double reference(const LombardiParameters& p, double mu_b, double N, double T, double E) {
    const long double e = E;
    const long double ac =
        static_cast<long double>(p.B) / e +
        static_cast<long double>(p.C) * std::pow(static_cast<long double>(N), static_cast<long double>(p.lambda)) /
            (static_cast<long double>(T) * std::pow(e, 1.0L / 3.0L));
    const long double sr = static_cast<long double>(p.delta) / (e * e);
    return static_cast<double>(1.0L / (1.0L / sr +
                                       1.0L / ac + 1.0L / static_cast<long double>(mu_b)));
}

}  // namespace

TEST_CASE("surface mobility: Lombardi's coefficients as reproduced (Unit 46)") {
    REQUIRE(lombardi_silicon_electrons.B == 4.75e7);
    REQUIRE(lombardi_silicon_electrons.C == 1.74e5);
    REQUIRE(lombardi_silicon_electrons.lambda == 0.125);
    REQUIRE(lombardi_silicon_electrons.delta == 5.82e14);
    REQUIRE(lombardi_silicon_holes.B == 9.93e6);
    REQUIRE(lombardi_silicon_holes.C == 8.84e5);
    REQUIRE(lombardi_silicon_holes.lambda == 0.0317);
    REQUIRE(lombardi_silicon_holes.delta == 2.05e14);
}

TEST_CASE("surface mobility: the formula against an independent evaluation (LS2)") {
    for (const LombardiParameters& p : {lombardi_silicon_electrons, lombardi_silicon_holes}) {
        for (const double mu_b : {1400.0, 450.0, 90.0}) {
            for (double E = 1e3; E <= 1e7; E *= 10.0) {
                for (const double N : {0.0, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20}) {
                    for (const double T : {250.0, 300.0, 400.0}) {
                        CAPTURE(p.B, mu_b, E, N, T);
                        const SurfaceMobility s = lombardi_mobility(p, mu_b, N, T, E);
                        REQUIRE(close(s.mobility, reference(p, mu_b, N, T, E), 1e-12));
                        REQUIRE(s.mobility < mu_b);
                    }
                }
            }
        }
    }
}

TEST_CASE("surface mobility: slopes against central differences (LS2)") {
    for (const LombardiParameters& p : {lombardi_silicon_electrons, lombardi_silicon_holes}) {
        for (double E = 1e3; E <= 1e7; E *= 10.0) {
            for (const double N : {0.0, 1e16, 1e18, 1e20}) {
                for (const double T : {250.0, 300.0, 400.0}) {
                    const double mu_b = 800.0;
                    CAPTURE(p.B, E, N, T);
                    const SurfaceMobility s = lombardi_mobility(p, mu_b, N, T, E);
                    const double hE = 1e-5 * E, hT = 1e-5 * T, hb = 1e-5 * mu_b;
                    const double dE = (lombardi_mobility(p, mu_b, N, T, E + hE).mobility -
                                       lombardi_mobility(p, mu_b, N, T, E - hE).mobility) /
                                      (2.0 * hE);
                    const double dT = (lombardi_mobility(p, mu_b, N, T + hT, E).mobility -
                                       lombardi_mobility(p, mu_b, N, T - hT, E).mobility) /
                                      (2.0 * hT);
                    const double db = (lombardi_mobility(p, mu_b + hb, N, T, E).mobility -
                                       lombardi_mobility(p, mu_b - hb, N, T, E).mobility) /
                                      (2.0 * hb);
                    REQUIRE(std::abs(s.d_dE - dE) <= 1e-6 * std::abs(dE));
                    REQUIRE(std::abs(s.d_dbulk - db) <= 1e-6 * std::abs(db));
                    if (N == 0.0) {
                        REQUIRE(s.d_dT == 0.0);
                    } else {
                        REQUIRE(std::abs(s.d_dT - dT) <= 1e-6 * std::abs(dT) + 1e-12 * s.mobility / T);
                    }
                }
            }
        }
    }
}

TEST_CASE("surface mobility: small at inversion fields, negligible at E_0") {
    // At 1 V/cm (the assembler's floor E_0) the surface terms change the bulk mobility by < 1e-4.
    for (const LombardiParameters& p : {lombardi_silicon_electrons, lombardi_silicon_holes}) {
        const double bulk = p.B == lombardi_silicon_electrons.B ? 1400.0 : 450.0;
        const double mu = lombardi_mobility(p, bulk, 1e17, 300.0, 1.0).mobility;
        REQUIRE(close(mu, bulk, 1e-4));
    }
    // At 1e6 V/cm surface roughness dominates: mu_sr = delta / E^2 = 582 (electrons).
    const double mu = lombardi_mobility(lombardi_silicon_electrons, 1400.0, 1e17, 300.0, 1e6).mobility;
    REQUIRE(mu < 400.0);
    REQUIRE(mu > 100.0);
}
