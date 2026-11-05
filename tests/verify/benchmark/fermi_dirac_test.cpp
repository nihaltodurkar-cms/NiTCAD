// TVB-02, the Fermi-Dirac integrals F_1/2, F_-1/2, F_3/2 and the inverse of F_1/2 (DECISIONS Unit 55 section 11.8;
// owner rulings TB-Q4-Q8 and the widened finding (ii) in 11.13). The benchmark itself is
// verify::benchmark::run_tvb02 (src/verify/benchmark/tvb02.cpp, RC-Q1); this file holds its outcomes to the frozen
// expectations and checks the reference: its uncertainty within 1e-3 of the criterion (else NOT EVALUATED), the
// polylogarithm series for eta <= -2, the zeta anchors at 0, and the tail bound against the tail it omits.
// Frozen expected outcome: F_1/2, F_3/2 and the inverse pass everywhere. F_-1/2: the frozen "pass for eta <= 15"
// was not met (misses from 13.76, worst 2.0e-11); the owner widened finding (ii) (option A, 2026-11-05), so its
// outcomes are recorded (a miss is reported as known finding (ii), RC-Q2). The hidden [.finding] case asserts
// F_-1/2 everywhere.
#include "NiTCAD/physics/fermi_dirac.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/contract.hpp"
#include "NiTCAD/verify/reference/double_double.hpp"
#include "NiTCAD/verify/reference/fermi_dirac.hpp"
#include "NiTCAD/verify/reference/tanh_sinh.hpp"

#include "../../physics/references.hpp"  // the existing gate's helper: TVD-14 agreement, recorded only
#include "../reference/exact_integer.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <string>

namespace v = NiTCAD::verify;
namespace ref = NiTCAD::verify::reference;
namespace bench = NiTCAD::verify::benchmark;
using ref::DD;
using ref::FermiOrder;
using v::Outcome;

namespace {

constexpr double budget = 1e-3 * 1e-12;

const bench::BenchmarkResult& result() {
    static const bench::BenchmarkResult r = bench::run_tvb02();
    return r;
}

const char* name(FermiOrder j) {
    switch (j) {
        case FermiOrder::minus_half: return "F_-1/2";
        case FermiOrder::half: return "F_1/2";
        case FermiOrder::three_halves: return "F_3/2";
    }
    return "?";
}

}  // namespace

TEST_CASE("TVB-02: the benchmark's outcomes (frozen)") {
    const auto& r = result();
    REQUIRE_FALSE(r.run_error.has_value());
    REQUIRE(r.comparisons.size() == 4 * 324);
    double worst[4] = {0, 0, 0, 0};
    int misses = 0;
    double first_miss = 1e300, worst_eta = 0.0;
    for (const auto& c : r.comparisons) {
        CAPTURE(c.quantity, c.parameter, c.metric, c.comparison.reason);
        REQUIRE(c.comparison.outcome != Outcome::not_evaluated);  // the reference counts everywhere
        const int k = c.quantity == "F_1/2" ? 0 : c.quantity == "F_-1/2" ? 1 : c.quantity == "F_3/2" ? 2 : 3;
        if (k == 1 && c.metric > worst[1]) worst_eta = c.parameter;
        worst[k] = std::max(worst[k], c.metric);
        if (k == 1) {
            CHECK((c.comparison.outcome == Outcome::pass || c.comparison.outcome == Outcome::finding));
            if (c.comparison.outcome == Outcome::finding) {
                ++misses;
                first_miss = std::min(first_miss, c.parameter);
                CHECK(c.comparison.reason.find("known finding (ii)") != std::string::npos);
            }
        } else {
            CHECK(c.comparison.outcome == Outcome::pass);
        }
    }
    UNSCOPED_INFO("worst: F_1/2 " << worst[0] << ", F_-1/2 " << worst[1] << " at eta = " << worst_eta << ", F_3/2 "
                                  << worst[2] << ", inverse " << worst[3] << "; finding (ii): " << misses
                                  << " misses, first at " << first_miss << "; " << r.notes.front());
    REQUIRE(bench::overall(r) == (misses > 0 ? Outcome::finding : Outcome::pass));
}

TEST_CASE("TVB-02: the reference agrees with the polylogarithm series for eta <= -2") {
    double worst = 0.0;
    for (const auto& c : result().comparisons) {
        if (c.parameter > -2.0 || c.quantity == "inverse F_1/2") continue;
        const FermiOrder j = c.quantity == "F_1/2"    ? FermiOrder::half
                             : c.quantity == "F_-1/2" ? FermiOrder::minus_half
                                                      : FermiOrder::three_halves;
        const ref::FermiSeries s = ref::fermi_series(j, DD(c.parameter));
        // q_ref is the reference rounded to a double (half an ulp, 1.1e-16, within the 1e-15 budget).
        const double e = (std::abs((s.value - DD(c.q_ref)).to_double()) + s.truncation_bound) / std::abs(c.q_ref);
        worst = std::max(worst, e);
        CAPTURE(c.parameter, name(j));
        REQUIRE(e <= budget);
    }
    UNSCOPED_INFO("worst series disagreement " << worst);
}

TEST_CASE("TVB-02: the reference at eta = 0 against OEIS's zeta digits") {
    struct Anchor {
        FermiOrder j;
        const char* zeta;
    };
    // zeta(3/2) A078434, zeta(1/2) A059750 (digits of -zeta(1/2)), zeta(5/2) A247041; 40 digits each.
    const Anchor anchors[] = {{FermiOrder::half, "2.612375348685488343348567567924071630570"},
                              {FermiOrder::minus_half, "-1.460354508809586812889499152515298012467"},
                              {FermiOrder::three_halves, "1.341487257250917179756769693348612136623"}};
    for (const Anchor& a : anchors) {
        const ref::FermiQuadrature q = ref::fermi_tanh_sinh(a.j, DD(0.0));
        const DD zeta = q.value / ref::fermi_eta_factor(a.j);
        const exact::LiteralComparison c = exact::compare_with_literal(zeta, a.zeta, 1);
        UNSCOPED_INFO(name(a.j) << "(0) / (1 - 2^-j) against zeta: " << c.measured_udd << " u_dd");
        CAPTURE(name(a.j), c.measured_udd);
        REQUIRE(c.measured_udd * ref::u_dd <= budget);
    }
}

TEST_CASE("TVB-02: the tail bound holds against the tail it omits") {
    // The omitted [U, inf) integrated over [U, sqrt(eta + 200)] (the rest is below e^-200 of the integrand's scale).
    for (const double eta : {-40.0, -2.0, 0.0, 20.0, 40.0}) {
        for (const FermiOrder j : {FermiOrder::minus_half, FermiOrder::half, FermiOrder::three_halves}) {
            const DD x(eta);
            const auto f = [&](const DD& u) { return ref::fermi_integrand(j, x, u); };
            const ref::TanhSinhResult t = ref::tanh_sinh(f, ref::fermi_cut(x), ref::sqrt(x + DD(200.0)), 0.0, 1e-20);
            const double bound = ref::fermi_tanh_sinh(j, x).tail_bound;
            CAPTURE(eta, name(j), t.value.to_double(), bound);
            REQUIRE(t.status == ref::TanhSinhStatus::ok);
            REQUIRE(t.value.to_double() > 0.0);
            REQUIRE(t.value.to_double() <= bound);
        }
    }
}

TEST_CASE("TVB-02: agreement with the existing gate's double-double helper (TVD-14, recorded, not asserted)") {
    // tests/physics/references.hpp, used by tests/physics/fermi_dirac_test.cpp: series below -1, Gauss-Legendre above.
    double w12 = 0.0, wm12 = 0.0;
    for (const auto& c : result().comparisons) {
        if (c.quantity == "F_1/2") {
            w12 = std::max(w12, std::abs(reference::fermi_half(c.parameter).value() - c.q_ref) / c.q_ref);
        } else if (c.quantity == "F_-1/2") {
            wm12 = std::max(wm12, std::abs(reference::fermi_minus_half(c.parameter).value() - c.q_ref) / c.q_ref);
        }
    }
    UNSCOPED_INFO("existing helper vs TVB-02 reference (both rounded to double): F_1/2 " << w12 << ", F_-1/2 "
                                                                                       << wm12);
    REQUIRE(std::isfinite(w12 + wm12));
}

TEST_CASE("TVB-02 finding (ii): F_-1/2 within 1e-12 on all of [-40, 40]", "[.finding]") {
    for (const auto& c : result().comparisons) {
        if (c.quantity != "F_-1/2") continue;
        CAPTURE(c.parameter, c.metric);
        CHECK(c.metric <= 1e-12);
    }
}
