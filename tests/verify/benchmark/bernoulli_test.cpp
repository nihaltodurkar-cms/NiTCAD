// TVB-03, the Bernoulli function B(x) = x / (e^x - 1) (DECISIONS Unit 55 section 11.8; owner ruling TB-Q9 in 11.13).
//
// Code under test: assemble::bernoulli (include/NiTCAD/assemble/bernoulli.hpp, header-only). Reference:
// verify/reference/bernoulli.hpp, double-double, direct and series (PA: the definition and OEIS's B_2k).
// Criterion, frozen before code: relative error <= 4 eps wherever B(x) is a normal double, on a logarithmic grid of
// +-[1e-300, 745] (20 points a decade), x = 0, +-1e-2 and their neighbouring doubles, and every 0.05 on [700, 745].
// The subnormal range is a diagnostic. The reference counts only where its uncertainty is <= 1e-3 of the criterion
// (the TB-Q4 rule): the series and the direct method agree to that on |x| <= 1/2, the identity B(-x) = B(x) + x holds
// to it on (1/2, 670], and the scaled form agrees with the direct one to it on [600, 670].
// Frozen expected outcome: pass at every judged point except where std::expm1(x) overflows (x > ln DBL_MAX =
// 709.78): fail (finding (iii) of 11.13). The hidden [.finding] case asserts pass everywhere and fails. B' is recorded
// as a diagnostic, not judged (TB-Q7).
#include "NiTCAD/assemble/bernoulli.hpp"
#include "NiTCAD/verify/reference/bernoulli.hpp"
#include "NiTCAD/verify/reference/contract.hpp"
#include "NiTCAD/verify/reference/double_double.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>
#include <vector>

namespace v = NiTCAD::verify;
namespace ref = NiTCAD::verify::reference;
using v::Outcome;

namespace {

constexpr double eps = std::numeric_limits<double>::epsilon();
constexpr double criterion = 4.0 * eps;
constexpr double reference_budget = 1e-3 * criterion;

v::ReferenceInfo info() {
    v::ReferenceInfo r;
    r.id = "TVB-03";
    r.cls = v::ReferenceClass::E;
    r.independence = v::Independence::I3;
    r.provenance = v::Provenance::PA;
    r.citation =
        "definition B(x) = x / (e^x - 1); Bernoulli numbers B_2k, OEIS A000367 and A002445 b-files "
        "(oeis.org, read 2026-11-05)";
    r.version = "OEIS b-files n = 1..17";
    r.equations = {3};  // the Scharfetter-Gummel flux of G3
    r.form = "x / expm1(x) in double-double; scaled x e^-r 2^-m above 670; -x below -670";
    r.assumptions = {"none: B is defined by its formula"};
    r.applicability = {"B(x) a normal double"};
    r.derivation = "verify/reference/bernoulli.hpp";
    r.verification = "TC1 (expm1, exp, ln 2); series vs direct; B(-x) = B(x) + x; scaled vs direct";
    r.second_method = "the generating-function series with OEIS's exact B_2k, |x| <= 1/2";
    return r;
}

std::vector<double> grid() {
    std::vector<double> g{0.0};
    for (int i = -300 * 20; i <= static_cast<int>(std::log10(745.0) * 20.0); ++i) {
        const double a = std::pow(10.0, i / 20.0);
        g.push_back(a);
        g.push_back(-a);
    }
    for (const double a : {1e-2, -1e-2}) {
        g.push_back(a);
        g.push_back(std::nextafter(a, 0.0));
        g.push_back(std::nextafter(a, 2.0 * a));
    }
    for (int i = 0; i <= 900; ++i) g.push_back(700.0 + 0.05 * i);
    g.push_back(745.0);
    g.push_back(-745.0);
    return g;
}

bool expm1_overflows(double x) { return x > std::log(DBL_MAX); }

struct Judged {
    double x;
    double error;
    v::Comparison c;
};

std::vector<Judged> judge_all() {
    std::vector<Judged> out;
    for (const double x : grid()) {
        const ref::ScaledDD b = ref::bernoulli_direct(x);
        const v::Predicate normal{"B(x) a normal double", ref::is_normal_double(b) ? 1.0 : 0.0, 1.0,
                                  ref::is_normal_double(b)};
        const double sim = NiTCAD::assemble::bernoulli(x);
        const double e = normal.holds ? ref::relative_error(sim, b) : 0.0;
        const v::Predicate pred[] = {normal};
        out.push_back({x, e, v::judge(info(), pred, e <= criterion, "relative error <= 4 eps")});
    }
    return out;
}

}  // namespace

TEST_CASE("TVB-03: the reference's methods agree to 1e-3 of the criterion") {
    double worst_series = 0.0, worst_identity = 0.0, worst_scaled = 0.0;
    for (const double x : grid()) {
        if (x != 0.0 && std::abs(x) <= 0.5) {
            const ref::DD s = ref::bernoulli_series(ref::DD(x));
            const ref::ScaledDD d = ref::bernoulli_direct(x);
            REQUIRE(d.exponent == 0);
            const double e = std::abs((s - d.mantissa).to_double()) / s.to_double();
            worst_series = std::max(worst_series, e);
            CAPTURE(x);
            REQUIRE(e <= reference_budget);
        }
        if (x > 0.5 && x <= 670.0) {
            const ref::DD bp = ref::bernoulli_direct(x).mantissa;
            const ref::DD bm = ref::bernoulli_direct(-x).mantissa;
            const double e = std::abs((bm - bp - ref::DD(x)).to_double()) / bm.to_double();
            worst_identity = std::max(worst_identity, e);
            CAPTURE(x);
            REQUIRE(e <= reference_budget);
        }
    }
    for (double x = 600.0; x <= 670.0; x += 0.37) {
        // The scaled form, forced, against the direct one.
        const int m = static_cast<int>(std::floor(x / ref::ln2().hi));
        const ref::DD r = ref::DD(x) - ref::DD(static_cast<double>(m)) * ref::ln2();
        const ref::DD scaled = ref::ldexp(ref::DD(x) * ref::exp(-r), -m);
        const ref::DD direct = ref::bernoulli_direct(x).mantissa;
        const double e = std::abs((scaled - direct).to_double()) / direct.to_double();
        worst_scaled = std::max(worst_scaled, e);
        CAPTURE(x);
        REQUIRE(e <= reference_budget);
    }
    UNSCOPED_INFO("worst disagreement in u_dd: series " << worst_series / ref::u_dd << ", identity "
                                                       << worst_identity / ref::u_dd << ", scaled "
                                                       << worst_scaled / ref::u_dd);
    // Anchors: B(0) = 1 exactly; B(1) = 1 / (e - 1) with TC1's e.
    REQUIRE(ref::bernoulli_direct(0.0).mantissa == ref::DD(1.0));
    const ref::DD b1 = ref::DD(1.0) / (ref::e() - ref::DD(1.0));
    REQUIRE(std::abs((ref::bernoulli_direct(1.0).mantissa - b1).to_double()) <= 8.0 * ref::u_dd);
    REQUIRE(std::abs((ref::bernoulli_series(ref::DD(0.5)) - ref::bernoulli_direct(0.5).mantissa).to_double()) <=
            reference_budget);
}

TEST_CASE("TVB-03: assemble::bernoulli against the reference (frozen outcomes)") {
    double worst = 0.0, worst_x = 0.0, worst_subnormal = 0.0;
    int judged = 0, refused = 0, failed = 0;
    double first_fail = 0.0, last_fail = 0.0;
    for (const Judged& j : judge_all()) {
        CAPTURE(j.x, j.error, j.c.reason);
        if (j.c.outcome == Outcome::not_applicable) {
            ++refused;
            const double sim = NiTCAD::assemble::bernoulli(j.x);
            const ref::ScaledDD b = ref::bernoulli_direct(j.x);
            worst_subnormal = std::max(worst_subnormal, ref::relative_error(sim, b));
            continue;
        }
        ++judged;
        if (expm1_overflows(j.x)) {
            CHECK(j.c.outcome == Outcome::fail);
            if (j.c.outcome == Outcome::fail) {
                if (failed++ == 0) first_fail = j.x;
                last_fail = j.x;
            }
        } else {
            CHECK(j.c.outcome == Outcome::pass);
            if (j.error > worst) {
                worst = j.error;
                worst_x = j.x;
            }
        }
    }
    UNSCOPED_INFO("judged " << judged << ", worst " << worst / eps << " eps at x = " << worst_x << "; finding (iii): "
                            << failed << " fails on [" << first_fail << ", " << last_fail << "]; subnormal range: "
                            << refused << " points, worst relative " << worst_subnormal);
    REQUIRE(failed > 0);
}

TEST_CASE("TVB-03: B' diagnostic (recorded, not judged)") {
    double worst = 0.0, worst_x = 0.0;
    for (const double x : grid()) {
        if (std::abs(x) > 670.0) continue;
        ref::DD d;
        if (std::abs(x) <= 0.5) {
            d = ref::bernoulli_series_derivative(ref::DD(x));
        } else {
            // B'(x) = B (1 - B - x) / x, from B(-x) = B(x) + x; no cancellation for |x| > 1/2.
            const ref::DD b = ref::bernoulli_direct(x).mantissa;
            d = b * (ref::DD(1.0) - b - ref::DD(x)) / ref::DD(x);
        }
        if (std::abs(d.hi) < DBL_MIN) continue;
        const double e = std::abs((ref::DD(NiTCAD::assemble::bernoulli_derivative(x)) - d).to_double()) /
                         std::abs(d.to_double());
        if (e > worst) {
            worst = e;
            worst_x = x;
        }
    }
    UNSCOPED_INFO("B' worst relative error " << worst << " (" << worst / eps << " eps) at x = " << worst_x);
    REQUIRE(std::isfinite(worst));
}

TEST_CASE("TVB-03 finding (iii): B within 4 eps where std::expm1 overflows", "[.finding]") {
    for (const Judged& j : judge_all()) {
        CAPTURE(j.x, j.error);
        if (j.c.outcome != Outcome::not_applicable) CHECK(j.c.outcome == Outcome::pass);
    }
}
