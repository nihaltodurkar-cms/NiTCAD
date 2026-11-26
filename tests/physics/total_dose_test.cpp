// Total-dose oxide charge (Unit 69, DECISIONS.md TDC2): the Co-60 charge yield's limits and the
// trapped-hole sheet density's constant and linearity.
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "NiTCAD/physics/total_dose.hpp"

using namespace NiTCAD::physics;

namespace {

bool close(double a, double b) { return std::abs(a - b) <= 1e-12 * std::abs(b); }

}  // namespace

TEST_CASE("total dose: yield limits", "[total_dose]") {
    CHECK(close(charge_yield(0.0), std::pow(0.1 / 1.35, 0.9)));
    CHECK(close(charge_yield(1e6), std::pow(1.1 / 2.35, 0.9)));
    CHECK(close(charge_yield(1e12), std::pow((1e12 + 0.1e6) / (1e12 + 1.35e6), 0.9)));
    // The approach to 1: 1 - f_y = m (E1 - E0) / E to first order; the second order is 1e-6 relative.
    CHECK(std::abs((1.0 - charge_yield(1e12)) / (0.9 * 1.25e6 / 1e12) - 1.0) < 1e-5);
    CHECK(charge_yield(0.0) < charge_yield(1e5));
    CHECK(charge_yield(1e5) < charge_yield(1e6));
    CHECK(charge_yield(1e6) < charge_yield(1e7));
}

TEST_CASE("total dose: trapped charge", "[total_dose]") {
    CHECK(dose_pair_density_cm3_rad == 8.1e12);
    const double base = oxide_trapped_charge_cm2(1e5, 1e6, 1e-6, 1.0);
    CHECK(close(base, 8.1e12 * 1e5 * 1e-6 * charge_yield(1e6)));
    CHECK(close(oxide_trapped_charge_cm2(2e5, 1e6, 1e-6, 1.0), 2.0 * base));
    CHECK(close(oxide_trapped_charge_cm2(1e5, 1e6, 2e-6, 1.0), 2.0 * base));
    CHECK(close(oxide_trapped_charge_cm2(1e5, 1e6, 1e-6, 0.5), 2.0 * oxide_trapped_charge_cm2(1e5, 1e6, 1e-6, 0.25)));
    CHECK(oxide_trapped_charge_cm2(1e5, 1e6, 1e-6, 0.0) == 0.0);
}
