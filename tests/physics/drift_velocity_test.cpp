// Unit 68 slice 4c infrastructure (DECISIONS.md amendment 10): the velocity-form Caughey-Thomas and
// piecewise laws against Is-haqzai et al. (2026) Table 11 (an undoped reference fixture, not a
// material), and Rashid's hole mobility on the Caughey-Thomas parameters.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <initializer_list>

#include "NiTCAD/physics/drift_velocity.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD::physics;

namespace {

// Table 11, 300 K: electrons PW, holes CT. 1 V/um = 1e4 V/cm.
constexpr HighFieldLaw electrons{.v_sat_cm_s = 26.59e6, .E_c_V_per_cm = 0.56e4, .beta = 0.41};
constexpr double mu01 = 1880.0;
constexpr HighFieldLaw holes{.v_sat_cm_s = 15.09e6, .E_c_V_per_cm = 0.55e4, .beta = 0.88};
constexpr double mu_h = 2744.0;

double relative(double a, double b) { return std::abs(a / b - 1.0); }

}  // namespace

TEST_CASE("drift velocity: Table 11 at the owner's regression fields", "[drift_velocity][unit68]") {
    struct Row {
        double E_V_per_um, v_n, v_p;  // [V/um], [cm/s]: owner's table, four digits; 0.56 V/um as accepted (A2)
    };
    for (const Row& r : {Row{0.009, 1.692e5, 2.396e5}, Row{0.08, 1.504e6, 1.813e6}, Row{0.56, 4.903e6, 6.927e6},
                         Row{1.0, 6.441e6, 8.903e6}, Row{12.0, 1.443e7, 1.403e7}}) {
        CAPTURE(r.E_V_per_um);
        const double E = r.E_V_per_um * 1e4;
        CHECK(relative(piecewise_velocity(mu01, E, electrons).velocity, r.v_n) < 5e-4);
        CHECK(relative(caughey_thomas_velocity(mu_h, E, holes.v_sat_cm_s, holes.beta).velocity, r.v_p) < 5e-4);
    }
    // mu0,2 = v_s / E_c: Table 11's 4748 and 2744.
    CHECK(relative(electrons.v_sat_cm_s / electrons.E_c_V_per_cm, 4748.0) < 1e-4);
    CHECK(relative(holes.v_sat_cm_s / holes.E_c_V_per_cm, 2744.0) < 1e-3);
}

TEST_CASE("drift velocity: the piecewise switch is where the branches meet", "[drift_velocity][unit68]") {
    const double Et = piecewise_switch_field(mu01, electrons);
    CHECK(relative(Et, 0.0852e4) < 1e-3);  // owner: about 0.0852 V/um, not the paper's 0.08
    const double upper = caughey_thomas_velocity(electrons.v_sat_cm_s / electrons.E_c_V_per_cm, Et,
                                                 electrons.v_sat_cm_s, electrons.beta).velocity;
    CHECK(relative(mu01 * Et, upper) < 1e-13);
    CHECK(relative(piecewise_velocity(mu01, Et * (1 - 1e-12), electrons).velocity,
                   piecewise_velocity(mu01, Et, electrons).velocity) < 1e-9);
    // A doping-dependent mu01 moves the switch and keeps the law continuous.
    for (const double mu : {500.0, 1000.0, 4000.0}) {
        CAPTURE(mu);
        const double e = piecewise_switch_field(mu, electrons);
        CHECK(relative(piecewise_velocity(mu, e * (1 - 1e-12), electrons).velocity,
                       piecewise_velocity(mu, e, electrons).velocity) < 1e-9);
    }
}

TEST_CASE("drift velocity: dv/dE analytic, mu0 at E = 0 for beta < 1", "[drift_velocity][unit68]") {
    for (const double beta : {0.41, 0.88, 1.0, 2.0}) {
        CAPTURE(beta);
        const DriftVelocity zero = caughey_thomas_velocity(mu_h, 0.0, holes.v_sat_cm_s, beta);
        CHECK(zero.velocity == 0.0);
        CHECK(zero.d_dE == mu_h);
        for (const double E : {1.0, 1e2, 1e4, 1e5}) {
            CAPTURE(E);
            const double h = E * 1e-6;
            const double fd = (caughey_thomas_velocity(mu_h, E + h, holes.v_sat_cm_s, beta).velocity -
                               caughey_thomas_velocity(mu_h, E - h, holes.v_sat_cm_s, beta).velocity) / (2 * h);
            CHECK(relative(caughey_thomas_velocity(mu_h, E, holes.v_sat_cm_s, beta).d_dE, fd) < 1e-6);
        }
    }
    CHECK(piecewise_velocity(mu01, 0.0, electrons).d_dE == mu01);
    // Saturation: v -> v_s.
    CHECK(relative(caughey_thomas_velocity(mu_h, 1e9, holes.v_sat_cm_s, holes.beta).velocity, holes.v_sat_cm_s) < 1e-3);
}

TEST_CASE("drift velocity: Rashid's mobility on the Caughey-Thomas parameters", "[drift_velocity][unit68]") {
    // Synthetic coefficients: the mapping itself, over T and N.
    const double mu_min0 = 120.0, n_min = -0.6, mu_max0 = 2000.0, n_max = -2.2, C0 = 3e17, m1 = 1.4, g2 = 0.7;
    SemiconductorParameters p = silicon_parameters;
    p.hole_mobility = rashid_mobility(mu_min0, n_min, mu_max0, n_max, C0, m1, g2);
    const Semiconductor m = *Semiconductor::create(p);
    for (const double T : {200.0, 300.0, 450.0}) {
        for (const double N : {1e14, 1e17, 3e18, 1e20}) {
            CAPTURE(T, N);
            const double t = T / 300.0;
            const double lo = mu_min0 * std::pow(t, n_min), hi = mu_max0 * std::pow(t, n_max);
            const double direct = lo + (hi - lo) / (1.0 + std::pow(t, m1) * std::pow(N / C0, g2));
            CHECK(relative(caughey_thomas_mobility(m, Carrier::hole, N, T), direct) < 1e-13);
        }
    }
}

TEST_CASE("drift velocity: Rashid Table III, T <= 343 K branch, as the diamond hole fixture (owner A4, A5)",
          "[drift_velocity][unit68]") {
    // Rashid et al. (2008) Table III, the T <= 343 K column; the 343-600 K branch (mu_max0 4399,
    // n_max -3.4) is not carried. Valid only to 343 K. Not a registered material.
    constexpr CaugheyThomasParameters rashid = rashid_mobility(55.0, 2.5, 3443.0, -1.55, 1.12e14, -2.0, 0.589);
    CHECK(relative(rashid.N_ref_T_exponent, 2.0 / 0.589) < 1e-15);  // m1 < 0: +3.40, sign kept
    SemiconductorParameters p = silicon_parameters;
    p.hole_mobility = rashid;
    const Semiconductor m = *Semiconductor::create(p);
    // Rashid Fig. 6 at N_A(i') = 1e11: about 3400 at 300 K; 3390 by the formula.
    CHECK(std::abs(caughey_thomas_mobility(m, Carrier::hole, 1e11, 300.0) - 3390.0) < 1.0);
    CHECK(caughey_thomas_mobility(m, Carrier::hole, 0.0, 300.0) == 3443.0);
    // Known discrepancy (A5): undoped 3443 against Is-haqzai Table 11's 2744 (v_s / E_c); not rescaled.
    CHECK(relative(3443.0, holes.v_sat_cm_s / holes.E_c_V_per_cm) > 0.2);
}
