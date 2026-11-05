// TVB-02, the Fermi-Dirac integrals F_1/2, F_-1/2, F_3/2 and the inverse of F_1/2 (DECISIONS Unit 55 section 11.8;
// owner rulings TB-Q4-Q8 in 11.13).
//
// Code under test: physics::fermi_half (value F_1/2, derivative F_-1/2), physics::fermi_three_halves (value F_3/2;
// its derivative, an F_1/2, a diagnostic) and physics::inverse_fermi_half. Reference: verify/reference/fermi_dirac.hpp
// (DLMF 25.12.14, PA): tanh-sinh (TC4) primary, adaptive Gauss-Kronrod (TC3) second method, the polylogarithm series
// for eta <= -2, and the exact anchors F_j(0) = (1 - 2^-j) zeta(j + 1) against OEIS's digits.
// Criteria, frozen before code: <= 1e-12 relative for each value over eta in [-40, 40]; the inverse
// |eta_sim - eta_ref| / max(1, |eta_ref|) <= 1e-12, eta_ref by Brent (TC5) on F_1/2,ref(eta) = nu for each double
// nu = fermi_half(eta_k).value. The reference counts only where its uncertainty (max of the two rules' difference
// and the primary's estimate, plus the tail bound) is <= 1e-3 of the criterion.
// Grid: eta_k = -40 + 0.0123 + 0.25 k inside [-40, 40] (320 points; the ruling's "321" counted the point
// 40.0123, outside the range), plus -2 +- 1e-9, 40 - 1e-9 and 0.
// Frozen expected outcome: F_1/2, F_3/2 and the inverse pass everywhere. F_-1/2: the frozen "pass for eta <= 15"
// was not met (misses from 13.76, worst 2.0e-11); the owner widened finding (ii) (option A, 2026-11-05), so its
// outcomes are recorded, not required. The hidden [.finding] case asserts F_-1/2 everywhere.
#include "NiTCAD/physics/fermi_dirac.hpp"
#include "NiTCAD/verify/reference/brent.hpp"
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
#include <vector>

namespace v = NiTCAD::verify;
namespace ref = NiTCAD::verify::reference;
namespace physics = NiTCAD::physics;
using ref::DD;
using ref::FermiOrder;
using v::Outcome;

namespace {

constexpr double criterion = 1e-12;
constexpr double budget = 1e-3 * criterion;

v::ReferenceInfo info() {
    v::ReferenceInfo r;
    r.id = "TVB-02";
    r.cls = v::ReferenceClass::E;
    r.independence = v::Independence::I3;
    r.provenance = v::Provenance::PA;
    r.citation =
        "NIST DLMF 25.12(iii), eqs. 25.12.14 and 25.12.16, and 25.12.10 (dlmf.nist.gov/25.12, v1.2.8, read "
        "2026-11-05); zeta(3/2), zeta(1/2), zeta(5/2) digits from OEIS A078434, A059750, A247041";
    r.version = "DLMF 1.2.8";
    r.equations = {4};
    r.form = "2 / Gamma(j + 1) int_0^U u^(2j+1) / (1 + e^(u^2 - eta)) du, U^2 = eta + 100, split at sqrt(eta)";
    r.assumptions = {"parabolic band (the integral's definition)"};
    r.derivation = "verify/reference/fermi_dirac.hpp (substitution, split, tail bound)";
    r.verification = "TC3, TC4 quadratures; polylog series for eta <= -2; zeta anchors at 0";
    r.second_method = "adaptive Gauss-Kronrod 7/15 (TC3) of the same form";
    return r;
}

std::vector<double> grid() {
    std::vector<double> g;
    for (int k = 0;; ++k) {
        const double e = -40.0 + 0.0123 + 0.25 * k;
        if (e > 40.0) break;
        g.push_back(e);
    }
    for (const double e : {-2.0 - 1e-9, -2.0 + 1e-9, 40.0 - 1e-9, 0.0}) g.push_back(e);
    return g;
}

struct Reference {
    DD value;
    double u_ref;  // absolute
};

Reference reference(FermiOrder j, double eta) {
    const ref::FermiQuadrature ts = ref::fermi_tanh_sinh(j, DD(eta));
    const ref::FermiQuadrature gk = ref::fermi_gauss_kronrod(j, DD(eta));
    REQUIRE(ts.ok);
    REQUIRE(gk.ok);
    const double diff = std::abs((ts.value - gk.value).to_double());
    return {ts.value, std::max(diff, ts.abs_error) + ts.tail_bound};
}

double relative(double sim, const DD& r) { return std::abs((DD(sim) - r).to_double()) / std::abs(r.to_double()); }

struct Row {
    double eta;
    Reference f12, fm12, f32;
};

const std::vector<Row>& rows() {
    static const std::vector<Row> r = [] {
        std::vector<Row> out;
        for (const double eta : grid()) {
            out.push_back({eta, reference(FermiOrder::half, eta), reference(FermiOrder::minus_half, eta),
                           reference(FermiOrder::three_halves, eta)});
        }
        return out;
    }();
    return r;
}

v::Comparison judge_value(double sim, const Reference& r, const char* what, double& worst) {
    const double e = relative(sim, r.value);
    worst = std::max(worst, e);
    return v::judge(info(), {}, e <= criterion, std::string(what) + " relative error <= 1e-12");
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

TEST_CASE("TVB-02: the reference counts (uncertainty <= 1e-3 of the criterion) and agrees with the series") {
    double worst_u = 0.0, worst_series = 0.0;
    for (const Row& r : rows()) {
        for (const auto* q : {&r.f12, &r.fm12, &r.f32}) {
            const double u = q->u_ref / std::abs(q->value.to_double());
            worst_u = std::max(worst_u, u);
            CAPTURE(r.eta);
            REQUIRE(u <= budget);
        }
        if (r.eta <= -2.0) {
            for (const FermiOrder j : {FermiOrder::minus_half, FermiOrder::half, FermiOrder::three_halves}) {
                const ref::FermiSeries s = ref::fermi_series(j, DD(r.eta));
                const Reference& q = j == FermiOrder::half ? r.f12 : j == FermiOrder::minus_half ? r.fm12 : r.f32;
                const double e = (std::abs((s.value - q.value).to_double()) + s.truncation_bound) /
                                 std::abs(q.value.to_double());
                worst_series = std::max(worst_series, e);
                CAPTURE(r.eta, name(j));
                REQUIRE(e <= budget);
            }
        }
    }
    UNSCOPED_INFO("worst relative reference uncertainty " << worst_u << "; worst series disagreement " << worst_series);
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

TEST_CASE("TVB-02: physics Fermi-Dirac integrals against the reference (frozen outcomes)") {
    double w12 = 0.0, wm12 = 0.0, w32 = 0.0, wm12_low = 0.0, wd32 = 0.0;
    int misses = 0;
    double first_miss = 0.0, worst_miss_eta = 0.0;
    for (const Row& r : rows()) {
        CAPTURE(r.eta);
        const physics::FermiIntegral f = physics::fermi_half(r.eta);
        const physics::FermiIntegral g = physics::fermi_three_halves(r.eta);
        CHECK(judge_value(f.value, r.f12, "F_1/2", w12).outcome == Outcome::pass);
        CHECK(judge_value(g.value, r.f32, "F_3/2", w32).outcome == Outcome::pass);
        const double before = wm12;
        const v::Comparison m = judge_value(f.derivative, r.fm12, "F_-1/2", wm12);
        if (wm12 > before) worst_miss_eta = r.eta;
        // Finding (ii) as widened by the owner (2026-11-05, option A): F_-1/2's outcome is recorded, not required.
        CHECK((m.outcome == Outcome::pass || m.outcome == Outcome::fail));
        if (m.outcome == Outcome::fail && misses++ == 0) first_miss = r.eta;
        if (r.eta <= 15.0) wm12_low = std::max(wm12_low, relative(f.derivative, r.fm12.value));
        wd32 = std::max(wd32, relative(g.derivative, r.f12.value));  // diagnostic
    }
    UNSCOPED_INFO("worst relative error: F_1/2 " << w12 << ", F_3/2 " << w32 << ", F_-1/2 " << wm12 << " at eta = "
                                                 << worst_miss_eta << " (" << wm12_low << " up to 15); finding (ii): "
                                                 << misses << " misses, first at " << first_miss
                                                 << "; diagnostic F_3/2 derivative vs F_1/2 " << wd32);
    REQUIRE(true);
}

TEST_CASE("TVB-02: the inverse of F_1/2 against Brent on the reference (frozen outcome)") {
    double worst = 0.0, worst_eta = 0.0;
    long long evaluations = 0;
    for (const double eta : grid()) {
        const double nu = physics::fermi_half(eta).value;
        const double sim = physics::inverse_fermi_half(nu);
        const double scale = std::max(1.0, std::abs(eta));
        auto f = [&](const DD& x) { return ref::fermi_tanh_sinh(FermiOrder::half, x).value - DD(nu); };
        // A bracket around the code's answer, widened by 100 until it holds the root.
        double delta = 1e-10 * scale;
        DD a = DD(sim) - DD(delta), b = DD(sim) + DD(delta);
        while ((f(a).hi < 0.0) == (f(b).hi < 0.0)) {
            delta *= 100.0;
            a = DD(sim) - DD(delta);
            b = DD(sim) + DD(delta);
            REQUIRE(delta < 1.0);
        }
        const ref::BrentResult root = ref::brent_zero(f, a, b, 1e-18 * scale);
        REQUIRE(root.status == ref::BrentStatus::ok);
        evaluations += root.evaluations;
        const double e = std::abs((DD(sim) - root.root).to_double()) / std::max(1.0, std::abs(root.root.to_double()));
        if (e > worst) {
            worst = e;
            worst_eta = eta;
        }
        CAPTURE(eta, nu, sim, e);
        CHECK(v::judge(info(), {}, e <= criterion, "inverse |d eta| / max(1, |eta|) <= 1e-12").outcome ==
              Outcome::pass);
    }
    UNSCOPED_INFO("inverse: worst " << worst << " at eta = " << worst_eta << "; Brent evaluations " << evaluations);
    REQUIRE(true);
}

TEST_CASE("TVB-02 finding (ii): F_-1/2 within 1e-12 on all of [-40, 40]", "[.finding]") {
    double w = 0.0;
    for (const Row& r : rows()) {
        CAPTURE(r.eta);
        CHECK(judge_value(physics::fermi_half(r.eta).derivative, r.fm12, "F_-1/2", w).outcome == Outcome::pass);
    }
}

TEST_CASE("TVB-02: agreement with the existing gate's double-double helper (TVD-14, recorded, not asserted)") {
    // tests/physics/references.hpp, used by tests/physics/fermi_dirac_test.cpp: series below -1, Gauss-Legendre above.
    double w12 = 0.0, wm12 = 0.0, at = 0.0;
    for (const Row& r : rows()) {
        const double h12 = reference::fermi_half(r.eta).value();
        const double hm12 = reference::fermi_minus_half(r.eta).value();
        w12 = std::max(w12, relative(h12, r.f12.value));
        const double e = relative(hm12, r.fm12.value);
        if (e > wm12) {
            wm12 = e;
            at = r.eta;
        }
    }
    UNSCOPED_INFO("existing helper vs TVB-02 reference (helper rounded to double): F_1/2 " << w12 << ", F_-1/2 " << wm12
                                                                                         << " at eta = " << at);
    REQUIRE(true);
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
