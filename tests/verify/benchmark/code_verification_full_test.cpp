// The full tier of TVR-2 (DECISIONS Unit 55 14.7, rulings CV-Q12-CV-Q16 in 14.9): TVB-04 on the tensor ladders and
// the 2D generated, 3D generated and swept families, against its frozen expected outcomes. The 3D levels use MKL
// PARDISO; the test takes the runtime from NITCAD_MKL_RUNTIME_PATH. Local only (owner ruling, option 1): 64^3 takes
// about 12 minutes on 6 threads and 20 GB, more than a CI runner holds, and the 3D generated and swept ladders about
// as long again, so the test exists only with NITCAD_VERIFY_FULL_TIER=ON (Release) and carries the ctest label "full".
#include <catch2/catch_test_macros.hpp>

#include <cstdio>
#include <string>

#include "NiTCAD/linalg/mkl_runtime.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/contract.hpp"

namespace bench = NiTCAD::verify::benchmark;
using NiTCAD::verify::Outcome;

namespace {

Outcome outcome_of(const bench::BenchmarkResult& r, const std::string& quantity) {
    for (const auto& c : r.comparisons) {
        if (c.quantity == quantity) return c.comparison.outcome;
    }
    FAIL("no comparison '" << quantity << "'");
    return Outcome::fail;
}

}  // namespace

TEST_CASE("TVB-04 full tier: tensor ladders and the unstructured families (frozen outcomes)",
          "[tvr2_full][tvb04_full]") {
    const std::string mkl = NITCAD_MKL_RUNTIME_PATH;
    if (mkl.empty()) SKIP("NITCAD_MKL_RUNTIME is not set: the 3D levels need PARDISO");
    REQUIRE(NiTCAD::linalg::load_mkl_runtime(mkl).has_value());
    bench::set_pardiso_threads(6);  // the P-cores of the owner's machine (local runs only)
    const bench::BenchmarkResult r = bench::run_tvb04(bench::Tier::full);
    REQUIRE_FALSE(r.run_error.has_value());
    for (const auto& n : r.notes) UNSCOPED_INFO("note: " << n);
    for (const auto& c : r.comparisons) {
        std::printf("TVB-04 %s: %.4g (%s) %s\n", c.quantity.c_str(), c.metric, c.criterion.c_str(),
                    c.comparison.reason.c_str());
    }
    REQUIRE(r.comparisons.size() == 15);
    for (const char* family : {"1D tensor", "2D tensor", "3D tensor"}) {
        CAPTURE(family);
        CHECK(outcome_of(r, std::string("psi L2 order, ") + family) == Outcome::pass);
    }
    CHECK(outcome_of(r, "psi Linf order, 1D tensor") == Outcome::pass);
    CHECK(outcome_of(r, "psi Linf order, 2D tensor") == Outcome::pass);
    // 3D L-infinity orders 1.60, 1.84, 1.93: spread 0.23 > 0.1, outside the asymptotic range (14.7, CV-Q14).
    CHECK(outcome_of(r, "psi Linf order, 3D tensor") == Outcome::not_evaluated);
    CHECK(outcome_of(r, "psi L2 finest, 1D tensor") == Outcome::pass);
    CHECK(outcome_of(r, "psi L2 finest, 2D tensor") == Outcome::finding);  // known finding (iv), CV-Q12
    CHECK(outcome_of(r, "psi L2 finest, 3D tensor") == Outcome::finding);
    // The unstructured families (CV-Q7): H1 orders p_th = 1 and finest H1 <= 1e-2, all passing. Measured: 2D generated
    // H1 orders 1.21, 1.11, 1.08 (p-hat 1.08, spread 0.099); 3D generated 1.52, 1.29 and swept 1.45, 1.23 (3 levels,
    // judged on the finest pair, 14.7); finest H1 7.9e-4, 3.3e-3, 3.3e-3.
    for (const char* family : {"2D generated", "3D generated", "3D swept"}) {
        CAPTURE(family);
        CHECK(outcome_of(r, std::string("psi H1 order, ") + family) == Outcome::pass);
        CHECK(outcome_of(r, std::string("psi H1 finest, ") + family) == Outcome::pass);
    }
    CHECK(bench::overall(r) == Outcome::not_evaluated);
}

TEST_CASE("TVB-07 full tier: N = 64 ... 512 (frozen outcomes)", "[tvb07_full]") {
    const bench::BenchmarkResult r = bench::run_tvb07(bench::Tier::full);
    REQUIRE_FALSE(r.run_error.has_value());
    for (const auto& n : r.notes) UNSCOPED_INFO("note: " << n);
    // V_bi at 4 levels; psi Linf and E_max orders (measured 1.96 and 2.00); psi Linf finest 6.7e-5, E_max 4.3e-5.
    REQUIRE(r.comparisons.size() == 8);
    for (const auto& c : r.comparisons) {
        CAPTURE(c.quantity, c.metric, c.comparison.reason);
        CHECK(c.comparison.outcome == Outcome::pass);
    }
    CHECK(bench::overall(r) == Outcome::pass);
}

TEST_CASE("TVB-09 full tier: the V_G grid at N = 64 ... 512 (frozen outcomes)", "[tvb09_full]") {
    const bench::BenchmarkResult r = bench::run_tvb09(bench::Tier::full);
    REQUIRE_FALSE(r.run_error.has_value());
    for (const auto& n : r.notes) UNSCOPED_INFO("note: " << n);
    // V_FB at 4 levels, then per V_G (11) the Q_s order and the finest Q_s and C_LF. At V_G = V_FB, Q_s = 0: its
    // three are not_applicable (|Q_s| below 1e-3 of the range). Measured: orders 1.997-2.000, finest Q_s 2e-6 to
    // 6e-6, C_LF 1.5e-6 to 6e-6.
    REQUIRE(r.comparisons.size() == 37);
    for (const auto& c : r.comparisons) {
        CAPTURE(c.quantity, c.metric, c.comparison.reason);
        const bool flat_band = c.quantity.ends_with("V_G - V_FB = 0 V");
        CHECK(c.comparison.outcome == (flat_band ? Outcome::not_applicable : Outcome::pass));
    }
    CHECK(bench::overall(r) == Outcome::pass);
}
