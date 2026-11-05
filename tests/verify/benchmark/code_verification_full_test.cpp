// The full tier of TVR-2 (DECISIONS Unit 55 14.7, rulings CV-Q12-CV-Q13 and the PARDISO ruling in 14.9): TVB-04 on
// the 1D, 2D and 3D tensor ladders, against its frozen expected outcomes. The 3D levels use MKL PARDISO; the test
// takes the runtime from NITCAD_MKL_RUNTIME_PATH. Local only (owner ruling, option 1): 64^3 takes about 12 minutes on
// 6 threads and 20 GB, more than a CI runner holds, so the test exists only with NITCAD_VERIFY_FULL_TIER=ON (Release)
// and carries the ctest label "full".
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

TEST_CASE("TVB-04 full tier: 1D, 2D and 3D tensor ladders (frozen outcomes)", "[tvr2_full]") {
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
    REQUIRE(r.comparisons.size() == 9);
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
    CHECK(bench::overall(r) == Outcome::not_evaluated);
}
