// Canali velocity-saturation mobility (ARCHITECTURE.md section 11, Unit 13): the legacy formula and
// silicon parameters, the limits, the exact field derivative against finite differences, and the
// parameter validation.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
#include <string>

#include "NiTCAD/physics/field_mobility.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD::physics;

namespace {

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

}  // namespace

TEST_CASE("field mobility: legacy silicon Canali parameters") {
    const Semiconductor si = silicon();
    REQUIRE(saturation(si, Carrier::electron) == CanaliParameters{1.07e7, 2.0});
    REQUIRE(saturation(si, Carrier::hole) == CanaliParameters{8.37e6, 1.0});
}

TEST_CASE("field mobility: the Canali formula for beta = 1 and 2") {
    const CanaliParameters n{1.07e7, 2.0}, p{8.37e6, 1.0};
    for (const double E : {0.0, 1e2, 1e3, 1e4, 3e4, 1e5, 1e6}) {
        CAPTURE(E);
        const double xn = 1360.0 * E / 1.07e7, xp = 495.0 * E / 8.37e6;
        REQUIRE(close(canali_mobility(1360.0, E, n).mobility, 1360.0 / std::sqrt(1.0 + xn * xn),
                      1e-15));
        REQUIRE(close(canali_mobility(495.0, E, p).mobility, 495.0 / (1.0 + xp), 1e-15));
    }
    // A general beta, against the formula written out.
    const CanaliParameters g{1e7, 1.5};
    const double x = 800.0 * 2e4 / 1e7;
    REQUIRE(close(canali_mobility(800.0, 2e4, g).mobility,
                  800.0 / std::pow(1.0 + std::pow(x, 1.5), 1.0 / 1.5), 1e-15));
}

TEST_CASE("field mobility: low-field limit, monotonic decrease, velocity saturation") {
    for (const CanaliParameters c : {CanaliParameters{1.07e7, 2.0}, CanaliParameters{8.37e6, 1.0},
                                     CanaliParameters{1e7, 1.3}}) {
        CAPTURE(c.v_sat_cm_s, c.beta);
        const double mu0 = 1000.0;
        REQUIRE(canali_mobility(mu0, 0.0, c).mobility == mu0);
        double mu_prev = mu0, v_prev = 0.0;
        for (double E = 1.0; E <= 1e9; E *= 1.7) {
            const double mu = canali_mobility(mu0, E, c).mobility;
            REQUIRE(mu < mu_prev);           // mobility falls with the field
            REQUIRE(mu * E > v_prev);        // the drift velocity rises
            REQUIRE(mu * E < c.v_sat_cm_s);  // and stays below v_sat
            mu_prev = mu;
            v_prev = mu * E;
        }
        REQUIRE(close(v_prev, c.v_sat_cm_s, 1e-3));  // at 1e9 V/cm it has saturated
        // Low field: mu = mu0 (1 - x^beta / beta + ...).
        const double E = 10.0, x = mu0 * E / c.v_sat_cm_s;
        REQUIRE(close(canali_mobility(mu0, E, c).mobility,
                      mu0 * (1.0 - std::pow(x, c.beta) / c.beta), 2.0 * std::pow(x, 2.0 * c.beta)));
    }
}

TEST_CASE("field mobility: the field derivative matches finite differences") {
    // Fourth-order central differences; the error is measured against mu / E, the natural scale of
    // dmu/dE (at low field the derivative itself is tiny: x^(beta-1) mu0 / v_sat).
    double worst = 0.0;
    for (const CanaliParameters c : {CanaliParameters{1.07e7, 2.0}, CanaliParameters{8.37e6, 1.0},
                                     CanaliParameters{1e7, 1.5}}) {
        for (const double mu0 : {95.0, 800.0, 1360.0}) {
            for (double E = 10.0; E <= 1e7; E *= 3.1) {
                const double h = 1e-3 * E;
                const auto mu = [&](double e) { return canali_mobility(mu0, e, c).mobility; };
                const double fd =
                    (-mu(E + 2 * h) + 8 * mu(E + h) - 8 * mu(E - h) + mu(E - 2 * h)) / (12 * h);
                const FieldMobility exact = canali_mobility(mu0, E, c);
                worst = std::max(worst, std::abs(fd - exact.d_dE) / (exact.mobility / E));
            }
        }
    }
    CAPTURE(worst);
    REQUIRE(worst < 1e-8);
}

TEST_CASE("field mobility: the derivative at zero field is the right-hand one") {
    // beta = 2: the mobility is even in E, so the slope at 0 is 0; beta = 1: -mu0^2 / v_sat.
    REQUIRE(canali_mobility(1360.0, 0.0, {1.07e7, 2.0}).d_dE == 0.0);
    REQUIRE(close(canali_mobility(495.0, 0.0, {8.37e6, 1.0}).d_dE, -495.0 * 495.0 / 8.37e6, 1e-15));
    const double h = 1e-3;
    const double one_sided =
        (canali_mobility(495.0, h, {8.37e6, 1.0}).mobility - 495.0) / h;
    REQUIRE(close(one_sided, -495.0 * 495.0 / 8.37e6, 1e-6));
}

TEST_CASE("field mobility: Canali parameters are validated") {
    using P = SemiconductorParameters;
    const auto rejected = [](void (*spoil)(P&), const std::string& name) {
        P p = silicon_parameters;
        spoil(p);
        const auto m = Semiconductor::create(p);
        REQUIRE_FALSE(m.has_value());
        CAPTURE(m.error().message);
        REQUIRE(m.error().message.find(name) != std::string::npos);
    };
    rejected([](P& p) { p.electron_saturation.v_sat_cm_s = 0.0; },
             "electron_saturation.v_sat_cm_s");
    rejected([](P& p) { p.hole_saturation.v_sat_cm_s = -1.0; }, "hole_saturation.v_sat_cm_s");
    rejected([](P& p) { p.electron_saturation.beta = 0.5; }, "electron_saturation.beta");
    rejected([](P& p) { p.hole_saturation.beta = std::numeric_limits<double>::quiet_NaN(); },
             "hole_saturation.beta");
    P ok = silicon_parameters;
    ok.hole_saturation.beta = 1.0;  // the boundary is allowed
    REQUIRE(Semiconductor::create(ok).has_value());
}
