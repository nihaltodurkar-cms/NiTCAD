// TVB-01, constants and thermal voltage (DECISIONS Unit 55 section 11.8; owner rulings TB-Q1-Q3 in 11.13).
//
// Code under test: include/NiTCAD/base/constants.hpp. Reference: verify/reference/codata2018.hpp, NIST's printed
// CODATA 2018 strings (PA). Criteria, frozen before code:
//   - q, k, eps0 (F/cm, the printed F/m value's exponent shifted by -2) and m0: bit-equal to the correctly rounded
//     double of the printed decimal (std::from_chars);
//   - hbar = h / 2 pi and k in eV/K = k / q, from the exact h, k, q in double-double: relative error <= gamma_3,
//     gamma_n = n u / (1 - n u), u = 2^-53 (the input roundings and the one operation; a correctly rounded value
//     is within u);
//   - V_T(T) = k T / q at T in {1, 4.2, 77, 300, 400, 1000} K (the double T): <= gamma_4 (fl(k), fl(q), the
//     product, the quotient).
// Frozen expected outcome: every comparison passes except hbar, whose base value is NIST's printed 10-digit
// truncation (finding (i) of 11.13): fail. The hidden [.finding] case asserts hbar's criterion and fails.
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/codata2018.hpp"
#include "NiTCAD/verify/reference/contract.hpp"
#include "NiTCAD/verify/reference/double_double.hpp"

#include <catch2/catch_test_macros.hpp>

#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <string>

namespace v = NiTCAD::verify;
namespace ref = NiTCAD::verify::reference;
namespace codata = NiTCAD::verify::reference::codata2018;
namespace base = NiTCAD::base;
using v::Outcome;

namespace {

constexpr double u = 0x1p-53;
constexpr double gamma(int n) { return n * u / (1.0 - n * u); }

ref::DD hbar_reference() { return codata::to_dd(codata::planck) / ref::ldexp(ref::pi(), 1); }

}  // namespace

TEST_CASE("TVB-01: the reference's own conversions") {
    // The printed strings parse to the digits typed, the exponent shift is exact, and 10^n is exact.
    REQUIRE(codata::shifted_decimal(codata::vacuum_permittivity, -2) == "8.8541878128e-14");
    REQUIRE(codata::shifted_decimal(codata::elementary_charge) == "1.602176634e-19");
    REQUIRE(codata::shifted_decimal(codata::reduced_planck) == "1.054571817e-34");
    REQUIRE(codata::is_truncated(codata::reduced_planck));
    REQUIRE(codata::is_truncated(codata::boltzmann_eV));
    REQUIRE_FALSE(codata::is_truncated(codata::planck));
    for (int n = 23; n <= 44; ++n) {
        // 10^n / 10^(n-1) = 10 exactly when both are exact.
        CAPTURE(n);
        REQUIRE(codata::exact_power_of_ten(n) == ref::DD(10.0) * codata::exact_power_of_ten(n - 1));
    }
    // The double-double value rounds to the correctly rounded double of the same decimal.
    for (const auto* c : {&codata::elementary_charge, &codata::boltzmann, &codata::planck, &codata::electron_mass}) {
        CAPTURE(c->quantity);
        REQUIRE(codata::to_dd(*c).to_double() == codata::to_double(codata::shifted_decimal(*c)));
    }
    // A truncated printed value is a lower bound of the exact one it truncates, within one unit of its last digit.
    const ref::DD hbar = hbar_reference();
    const ref::DD printed = codata::to_dd(codata::reduced_planck);
    REQUIRE(printed < hbar);
    REQUIRE((hbar - printed).to_double() < 1e-43);
    const ref::DD kq = codata::to_dd(codata::boltzmann) / codata::to_dd(codata::elementary_charge);
    REQUIRE(codata::to_dd(codata::boltzmann_eV) < kq);
    REQUIRE((kq - codata::to_dd(codata::boltzmann_eV)).to_double() < 1e-14);
}

TEST_CASE("TVB-01: base constants against CODATA 2018 (frozen outcomes)") {
    const auto r = NiTCAD::verify::benchmark::run_tvb01();
    REQUIRE_FALSE(r.run_error.has_value());
    REQUIRE(r.comparisons.size() == 12);
    for (const auto& c : r.comparisons) {
        CAPTURE(c.quantity, c.parameter, c.metric, c.criterion, c.comparison.reason);
        UNSCOPED_INFO(c.quantity << " " << c.parameter << ": " << c.metric << " (" << c.criterion << ")");
        if (c.quantity.starts_with("hbar")) {
            // Finding (i), frozen as a fail; reported as a known finding (RC-Q2).
            CHECK(c.metric > gamma(3));
            CHECK(c.comparison.outcome == Outcome::finding);
            CHECK(c.comparison.reason.find("known finding (i)") != std::string::npos);
        } else {
            CHECK(c.comparison.outcome == Outcome::pass);
        }
    }
    CHECK(r.notes.front().ends_with("yes"));  // hbar is the printed truncation
    REQUIRE(NiTCAD::verify::benchmark::overall(r) == Outcome::finding);
}

TEST_CASE("TVB-01 finding (i): hbar within gamma_3 of h / 2 pi", "[.finding]") {
    CHECK(std::abs((ref::DD(base::hbar_J_s) - hbar_reference()).to_double()) / hbar_reference().to_double() <=
          gamma(3));
}
