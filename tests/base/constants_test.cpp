#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "NiTCAD/base/constants.hpp"

using Catch::Approx;
using namespace NiTCAD::base;

TEST_CASE("constants: values are the pinned CODATA 2018 numbers") {
    REQUIRE(q_C == 1.602176634e-19);
    REQUIRE(k_B_J_per_K == 1.380649e-23);
    REQUIRE(eps0_F_per_cm == 8.8541878128e-14);  // 8.8541878128e-12 F/m
    REQUIRE(hbar_J_s == 1.054571817e-34);
    REQUIRE(m0_kg == 9.1093837015e-31);
}

TEST_CASE("constants: k_B in eV/K follows from k_B and q") {
    REQUIRE(k_B_eV_per_K == Approx(8.617333262145e-5).epsilon(1e-12));
    // The legacy project hard-coded the truncated 8.617333262e-5; the difference is ~1.7e-11.
    REQUIRE(k_B_eV_per_K == Approx(8.617333262e-5).epsilon(1e-9));
}

TEST_CASE("constants: thermal voltage at 300 K is 25.852 mV") {
    static_assert(thermal_voltage(300.0) > 0.0258 && thermal_voltage(300.0) < 0.0259);
    REQUIRE(thermal_voltage(300.0) == Approx(0.0258519998).epsilon(1e-9));
    REQUIRE(thermal_voltage(300.0) == Approx(0.025852).epsilon(1e-6));
}

TEST_CASE("constants: thermal voltage scales linearly and satisfies q V_T = k_B T") {
    REQUIRE(thermal_voltage(600.0) == Approx(2.0 * thermal_voltage(300.0)).epsilon(1e-14));
    REQUIRE(q_C * thermal_voltage(350.0) == Approx(k_B_J_per_K * 350.0).epsilon(1e-14));
}
