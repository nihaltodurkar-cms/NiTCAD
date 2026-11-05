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

v::ReferenceInfo info() {
    v::ReferenceInfo r;
    r.id = "TVB-01";
    r.cls = v::ReferenceClass::E;
    r.independence = v::Independence::I3;
    r.provenance = v::Provenance::PA;
    r.citation =
        "NIST, Fundamental Physical Constants, 2018 CODATA adjustment, complete listing "
        "(physics.nist.gov/cuu/Constants/ArchiveASCII/allascii_2018.txt, read 2026-11-05)";
    r.version = "CODATA 2018";
    r.equations = {1, 4};  // eps0 and q in G1; k T / q in G4
    r.form = "printed decimals; hbar = h / (2 pi), k[eV/K] = k / q, V_T = k T / q from the exact h, k, q";
    r.assumptions = {"SI 2019 exact h, k, q; eps0 and m0 as printed (measured)"};
    r.derivation = "verify/reference/codata2018.hpp: exact decimal-to-double-double conversion";
    r.verification = "TC1 (double-double, pi); std::from_chars correct rounding";
    return r;
}

double relative(double sim, const ref::DD& reference) {
    return std::abs((ref::DD(sim) - reference).to_double()) / std::abs(reference.to_double());
}

std::string bits(double x) {
    char s[32];
    std::snprintf(s, sizeof s, "%016llx", static_cast<unsigned long long>(std::bit_cast<std::uint64_t>(x)));
    return s;
}

v::Comparison bit_equal(double sim, const std::string& decimal) {
    const double r = codata::to_double(decimal);
    UNSCOPED_INFO(decimal << ": base " << bits(sim) << ", from_chars " << bits(r));
    return v::judge(info(), {}, std::bit_cast<std::uint64_t>(sim) == std::bit_cast<std::uint64_t>(r),
                    "bit-equal to " + decimal);
}

v::Comparison within(double sim, const ref::DD& reference, int n, const std::string& what) {
    const double e = relative(sim, reference);
    UNSCOPED_INFO(what << ": relative error " << e << " (gamma_" << n << " = " << gamma(n) << ")");
    return v::judge(info(), {}, e <= gamma(n), what + " relative error <= gamma_" + std::to_string(n));
}

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
    CHECK(bit_equal(base::q_C, codata::shifted_decimal(codata::elementary_charge)).outcome == Outcome::pass);
    CHECK(bit_equal(base::k_B_J_per_K, codata::shifted_decimal(codata::boltzmann)).outcome == Outcome::pass);
    CHECK(bit_equal(base::eps0_F_per_cm, codata::shifted_decimal(codata::vacuum_permittivity, -2)).outcome ==
          Outcome::pass);
    CHECK(bit_equal(base::m0_kg, codata::shifted_decimal(codata::electron_mass)).outcome == Outcome::pass);

    const ref::DD k = codata::to_dd(codata::boltzmann);
    const ref::DD q = codata::to_dd(codata::elementary_charge);
    CHECK(within(base::k_B_eV_per_K, k / q, 3, "k / q [eV/K]").outcome == Outcome::pass);

    // Finding (i): base's hbar is NIST's printed truncation, not h / 2 pi.
    const v::Comparison h = within(base::hbar_J_s, hbar_reference(), 3, "hbar = h / 2 pi");
    CHECK(h.outcome == Outcome::fail);
    // Diagnostic: it is bit-equal to the printed (truncated) digits.
    CHECK(bit_equal(base::hbar_J_s, codata::shifted_decimal(codata::reduced_planck)).outcome == Outcome::pass);

    for (const double T : {1.0, 4.2, 77.0, 300.0, 400.0, 1000.0}) {
        CAPTURE(T);
        CHECK(within(base::thermal_voltage(T), k * ref::DD(T) / q, 4, "V_T(" + std::to_string(T) + " K)").outcome ==
              Outcome::pass);
    }
}

TEST_CASE("TVB-01 finding (i): hbar within gamma_3 of h / 2 pi", "[.finding]") {
    CHECK(within(base::hbar_J_s, hbar_reference(), 3, "hbar = h / 2 pi").outcome == Outcome::pass);
}
