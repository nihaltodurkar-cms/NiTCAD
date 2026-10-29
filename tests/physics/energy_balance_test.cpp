// The energy-balance model's carrier-temperature mobility and effective field (Unit 52;
// DECISIONS.md ED4-ED6, EC2 and the owner's ruling on beta = 1): the uniform-field balance
// inverted, the closed forms against Newton, the regularization below a0, cooled carriers, and
// every partial against central differences.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/physics/energy_balance.hpp"
#include "NiTCAD/physics/field_mobility.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD::physics;

namespace {

double relative(double a, double b) { return std::abs(a - b) / std::abs(b); }

constexpr double k_over_q = NiTCAD::base::k_B_J_per_K / NiTCAD::base::q_C;

// The heating a uniform field E gives: (3/2) k rise / tau_w = q v(E) E.
double uniform_rise(double mu0, double E, double tau_w, const CanaliParameters& c) {
    const double v = canali_mobility(mu0, E, c).mobility * E;
    return 2.0 / 3.0 * tau_w * v * E / k_over_q;
}

// a for a rise, as carrier_temperature_mobility forms it.
double reduced(double mu0, double rise, double tau_w, const CanaliParameters& c) {
    return 1.5 * k_over_q * mu0 / (tau_w * c.v_sat_cm_s * c.v_sat_cm_s) * rise;
}

}  // namespace

TEST_CASE("energy balance: silicon's relaxation times (EQ1) and the default of none") {
    REQUIRE(silicon_energy_relaxation.tau_n_s == 0.35e-12);
    REQUIRE(silicon_energy_relaxation.tau_p_s == 0.4e-12);
    // Not part of the parameter set, so silicon's device files do not change (EC1).
    REQUIRE(silicon_parameters.energy_relaxation == EnergyRelaxationParameters{});
    SemiconductorParameters p = silicon_parameters;
    p.energy_relaxation.tau_n_s = -1e-12;
    REQUIRE_FALSE(Semiconductor::create(p).has_value());
    p.energy_relaxation.tau_n_s = 1e-12;
    REQUIRE(Semiconductor::create(p).has_value());
}

TEST_CASE("energy balance: a uniform field's heating gives back its field and Canali mobility") {
    // EC2: 1e-13 relative, beta = 1 and 2 (and 1.5 through Newton), 1e2-1e6 V/cm, 0.1-1 ps.
    double worst = 0.0;
    for (const CanaliParameters c :
         {CanaliParameters{1.07e7, 2.0}, CanaliParameters{8.37e6, 1.0},
          CanaliParameters{1e7, 1.5}}) {
        for (const double mu0 : {1360.0, 495.0, 92.0}) {
            for (const double tau_w : {0.1e-12, 0.35e-12, 1e-12}) {
                for (double E = 1e2; E <= 1.0001e6; E *= std::sqrt(10.0)) {
                    CAPTURE(c.beta, mu0, tau_w, E);
                    const double rise = uniform_rise(mu0, E, tau_w, c);
                    REQUIRE(reduced(mu0, rise, tau_w, c) >= energy_balance_regularization);
                    const CarrierTemperatureMobility m =
                        carrier_temperature_mobility(mu0, rise, tau_w, c);
                    const double mu_c = canali_mobility(mu0, E, c).mobility;
                    worst = std::max({worst, relative(m.field, E), relative(m.mobility, mu_c)});
                    REQUIRE(relative(m.field, E) <= 1e-13);
                    REQUIRE(relative(m.mobility, mu_c) <= 1e-13);
                }
            }
        }
    }
    UNSCOPED_INFO("worst inversion error " << worst);
}

TEST_CASE("energy balance: the closed forms equal the Newton root") {
    // EC2: beta = 1 (and 2) within 1e-14, over a from just above a0 to 1e6.
    for (const double beta : {1.0, 2.0}) {
        for (double a = 1.0001e-8; a <= 1e6; a *= 3.7) {
            CAPTURE(beta, a);
            const ReducedField closed = reduced_field(a, beta);
            const ReducedField newton = reduced_field_newton(a, beta);
            REQUIRE(relative(closed.x, newton.x) <= 1e-14);
            REQUIRE(relative(closed.d_a, newton.d_a) <= 1e-14);
            // And it solves x^2 m(x) = a.
            const double m = std::pow(1.0 + std::pow(closed.x, beta), -1.0 / beta);
            REQUIRE(relative(closed.x * closed.x * m, a) <= 4e-15);
        }
    }
}

TEST_CASE("energy balance: the regularization below a0 (owner ruling on beta = 1)") {
    constexpr double a0 = energy_balance_regularization;
    for (const double beta : {1.0, 1.5, 2.0, 3.0}) {
        CAPTURE(beta);
        // Value and slope continuous at a0.
        const ReducedField below = reduced_field(std::nextafter(a0, 0.0), beta);
        const ReducedField at = reduced_field(a0, beta);
        REQUIRE(relative(below.x, at.x) <= 1e-14);
        REQUIRE(relative(below.d_a, at.d_a) <= 1e-13);
        // Zero with zero slope at 0 and for cooled carriers; increasing in between.
        REQUIRE(reduced_field(0.0, beta).x == 0.0);
        REQUIRE(reduced_field(0.0, beta).d_a == 0.0);
        REQUIRE(reduced_field(-0.3, beta).x == 0.0);
        REQUIRE(reduced_field(-0.3, beta).d_a == 0.0);
        double previous = 0.0;
        for (int k = 1; k <= 1000; ++k) {
            const double a = a0 * k / 1000.0;
            const ReducedField r = reduced_field(a, beta);
            CAPTURE(a);
            REQUIRE(r.x > previous);
            REQUIRE(r.d_a > 0.0);
            REQUIRE(r.x <= at.x * (1.0 + 1e-15));
            previous = r.x;
        }
    }
    // The mobility the regularization changes: below a0 only, by at most 1e-4 of mu0 (beta = 1,
    // where x ~ sqrt(a) moves most; the exact root from Newton).
    const CanaliParameters holes{8.37e6, 1.0};
    double worst = 0.0;
    for (int k = 1; k <= 1000; ++k) {
        const double a = a0 * k / 1000.0;
        const double x_exact = reduced_field_newton(a, 1.0).x;
        const double x_used = reduced_field(a, 1.0).x;
        worst = std::max(worst, std::abs(1.0 / (1.0 + x_used) - 1.0 / (1.0 + x_exact)));
    }
    UNSCOPED_INFO("largest mobility change below a0, relative to mu0: " << worst);
    REQUIRE(worst <= 1e-4);
    // Silicon holes: a0 is a few microkelvin of heating.
    const double rise_a0 = a0 / reduced(495.0, 1.0, 0.4e-12, holes);
    REQUIRE(rise_a0 < 1e-5);
}

TEST_CASE("energy balance: at the lattice temperature and below it") {
    const CanaliParameters electrons{1.07e7, 2.0}, holes{8.37e6, 1.0};
    for (const CanaliParameters& c : {electrons, holes}) {
        const CarrierTemperatureMobility m = carrier_temperature_mobility(800.0, 0.0, 0.35e-12, c);
        REQUIRE(m.mobility == 800.0);
        REQUIRE(m.field == 0.0);
        REQUIRE(m.d_mu0 == 1.0);
    }
    // beta = 1: cooled carriers keep mu0 and its slope is 0 (C1 at T_L).
    const CarrierTemperatureMobility cold = carrier_temperature_mobility(495.0, -20.0, 0.4e-12, holes);
    REQUIRE(cold.mobility == 495.0);
    REQUIRE(cold.d_rise == 0.0);
    REQUIRE(carrier_temperature_mobility(495.0, 0.0, 0.4e-12, holes).d_rise == 0.0);
    // beta = 2: analytic through T_L, cooled carriers faster; the slope at T_L is -mu0 da/drise / 2
    // from either side.
    const double a_rise = reduced(1360.0, 1.0, 0.35e-12, electrons);
    const double slope0 = carrier_temperature_mobility(1360.0, 0.0, 0.35e-12, electrons).d_rise;
    REQUIRE(relative(slope0, -1360.0 * a_rise / 2.0) <= 1e-15);
    const CarrierTemperatureMobility cooled =
        carrier_temperature_mobility(1360.0, -50.0, 0.35e-12, electrons);
    REQUIRE(cooled.mobility > 1360.0);
    REQUIRE(cooled.field == 0.0);
    REQUIRE(std::isfinite(carrier_temperature_mobility(1360.0, -300.0, 0.35e-12, electrons).mobility));
}

TEST_CASE("energy balance: every partial against central differences") {
    // EC2: 1e-7 relative (or 1e-7 of the scale where the partial is near zero).
    const auto check = [](double mu0, double rise, double tau_w, const CanaliParameters& c) {
        const CarrierTemperatureMobility m = carrier_temperature_mobility(mu0, rise, tau_w, c);
        const double hm = 1e-6 * mu0;
        // beta = 2 is analytic through T_L (a step on the scale of 1 K); otherwise x ~ sqrt(rise)
        // and the step must be small against the rise itself.
        const double hr = 1e-4 * (c.beta == 2.0 ? std::max(std::abs(rise), 1.0) : std::abs(rise));
        const auto at = [&](double u, double r) {
            return carrier_temperature_mobility(u, r, tau_w, c);
        };
        const CarrierTemperatureMobility mp = at(mu0 + hm, rise), mm = at(mu0 - hm, rise);
        const CarrierTemperatureMobility rp = at(mu0, rise + hr), rm = at(mu0, rise - hr);
        const auto close = [](double fd, double exact, double scale) {
            return std::abs(fd - exact) <= 1e-7 * std::max(std::abs(exact), scale);
        };
        CAPTURE(mu0, rise, c.beta);
        REQUIRE(close((mp.mobility - mm.mobility) / (2 * hm), m.d_mu0, 1.0));
        REQUIRE(close((rp.mobility - rm.mobility) / (2 * hr), m.d_rise, 1e-6 * mu0));
        const double field_scale = std::max(m.field, 1.0);
        REQUIRE(close((mp.field - mm.field) / (2 * hm), m.field_d_mu0, field_scale / mu0));
        // The field grows like sqrt(rise) for every beta: a step relative to the rise.
        const double he = 1e-4 * std::abs(rise);
        const double fe = (at(mu0, rise + he).field - at(mu0, rise - he).field) / (2 * he);
        REQUIRE(close(fe, m.field_d_rise, 1e-6 * field_scale));
    };
    for (const CanaliParameters c :
         {CanaliParameters{1.07e7, 2.0}, CanaliParameters{8.37e6, 1.0},
          CanaliParameters{1e7, 1.5}}) {
        for (const double rise : {-100.0, 1e-3, 0.5, 30.0, 700.0, 5000.0}) {
            if (rise < 0.0 && c.beta != 2.0) continue;  // constant there; checked above
            check(1360.0, rise, 0.35e-12, c);
            check(200.0, rise, 0.4e-12, c);
        }
    }
}
