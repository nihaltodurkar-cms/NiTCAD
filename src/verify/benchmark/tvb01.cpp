// TVB-01, constants and thermal voltage against CODATA 2018 (DECISIONS Unit 55 11.8; rulings TB-Q1-Q3 in 11.13).
//
// Code under test: base/constants.hpp. Reference: verify/reference/codata2018.hpp (NIST's printed strings, PA).
// Criteria (frozen): q, k, eps0 (F/cm) and m0 bit-equal to the correctly rounded double of the printed decimal;
// hbar = h / 2 pi and k / q within gamma_3, V_T(T) = k T / q within gamma_4 of the exact double-double values,
// gamma_n = n u / (1 - n u), u = 2^-53. Finding (i) (hbar the printed truncation) was fixed in Unit 56 (FQ1).
#include <bit>
#include <cmath>
#include <cstdint>
#include <string>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/codata2018.hpp"
#include "NiTCAD/verify/reference/double_double.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

namespace ref = verify::reference;
namespace codata = verify::reference::codata2018;

constexpr double u = 0x1p-53;
constexpr double gamma(int n) { return n * u / (1.0 - n * u); }

ReferenceInfo info() {
    ReferenceInfo r;
    r.id = "TVB-01";
    r.cls = ReferenceClass::E;
    r.independence = Independence::I3;
    r.provenance = Provenance::PA;
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

ComparisonRecord bit_equal(const std::string& quantity, const char* unit, double sim, const std::string& decimal) {
    const double r = codata::to_double(decimal);
    const bool equal = std::bit_cast<std::uint64_t>(sim) == std::bit_cast<std::uint64_t>(r);
    ComparisonRecord c{quantity, unit, 0.0, sim, r, 0.0, std::abs(sim - r), "bit-equal to " + decimal, {}};
    c.comparison = judge(info(), {}, equal, c.criterion);
    return c;
}

ComparisonRecord within(const std::string& quantity, const char* unit, double parameter, double sim,
                        const ref::DD& reference, int n) {
    const double e = std::abs((ref::DD(sim) - reference).to_double()) / std::abs(reference.to_double());
    ComparisonRecord c{quantity, unit, parameter, sim, reference.to_double(), 0.0, e,
                       "relative error <= gamma_" + std::to_string(n) + " = " + shortest(gamma(n)), {}};
    c.comparison = judge(info(), {}, e <= gamma(n), c.criterion);
    return c;
}

}  // namespace

BenchmarkResult run_tvb01() {
    BenchmarkResult r;
    r.id = "TVB-01";
    r.title = "Constants and thermal voltage against CODATA 2018";
    r.info = info();
    r.comparisons.push_back(bit_equal("q", "C", base::q_C, codata::shifted_decimal(codata::elementary_charge)));
    r.comparisons.push_back(bit_equal("k_B", "J/K", base::k_B_J_per_K, codata::shifted_decimal(codata::boltzmann)));
    r.comparisons.push_back(
        bit_equal("eps0", "F/cm", base::eps0_F_per_cm, codata::shifted_decimal(codata::vacuum_permittivity, -2)));
    r.comparisons.push_back(bit_equal("m0", "kg", base::m0_kg, codata::shifted_decimal(codata::electron_mass)));
    const ref::DD k = codata::to_dd(codata::boltzmann);
    const ref::DD q = codata::to_dd(codata::elementary_charge);
    r.comparisons.push_back(within("k_B / q", "eV/K", 0.0, base::k_B_eV_per_K, k / q, 3));
    const ref::DD hbar = codata::to_dd(codata::planck) / ref::ldexp(ref::pi(), 1);
    r.comparisons.push_back(within("hbar = h / 2 pi", "J s", 0.0, base::hbar_J_s, hbar, 3));
    for (const double T : {1.0, 4.2, 77.0, 300.0, 400.0, 1000.0}) {
        r.comparisons.push_back(within("V_T", "V", T, base::thermal_voltage(T), k * ref::DD(T) / q, 4));
    }
    r.notes.push_back("diagnostic: base::hbar_J_s is bit-equal to NIST's printed (truncated) " +
                      codata::shifted_decimal(codata::reduced_planck) + ": " +
                      (std::bit_cast<std::uint64_t>(base::hbar_J_s) ==
                               std::bit_cast<std::uint64_t>(
                                   codata::to_double(codata::shifted_decimal(codata::reduced_planck)))
                           ? "yes"
                           : "no"));
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
