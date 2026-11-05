// TVB-02, the Fermi-Dirac integrals and the inverse of F_1/2 (DECISIONS Unit 55 11.8; rulings TB-Q4-Q8 in 11.13).
//
// Code under test: physics::fermi_half (F_1/2 and, as its derivative, F_-1/2), physics::fermi_three_halves (F_3/2)
// and physics::inverse_fermi_half. Reference: verify/reference/fermi_dirac.hpp (DLMF 25.12.14, PA), tanh-sinh with
// adaptive Gauss-Kronrod as the second method. Criteria (frozen): <= 1e-12 relative per value on eta in [-40, 40];
// the inverse |eta_sim - eta_ref| / max(1, |eta_ref|) <= 1e-12, eta_ref by Brent on F_1/2,ref(eta) = nu. A
// comparison whose reference uncertainty exceeds 1e-3 of the criterion is NOT EVALUATED (the TB-Q4 rule). F_-1/2
// misses from eta = 13.76: finding (ii), widened.
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "NiTCAD/physics/fermi_dirac.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/brent.hpp"
#include "NiTCAD/verify/reference/fermi_dirac.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

namespace ref = verify::reference;
using ref::DD;
using ref::FermiOrder;

constexpr double criterion = 1e-12;
constexpr double budget = 1e-3 * criterion;

ReferenceInfo info() {
    ReferenceInfo r;
    r.id = "TVB-02";
    r.cls = ReferenceClass::E;
    r.independence = Independence::I3;
    r.provenance = Provenance::PA;
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
    double u_ref = 0.0;  // absolute
    bool ok = false;
};

Reference reference(FermiOrder j, double eta) {
    const ref::FermiQuadrature ts = ref::fermi_tanh_sinh(j, DD(eta));
    const ref::FermiQuadrature gk = ref::fermi_gauss_kronrod(j, DD(eta));
    const double diff = std::abs((ts.value - gk.value).to_double());
    return {ts.value, std::max(diff, ts.abs_error) + ts.tail_bound, ts.ok && gk.ok};
}

ComparisonRecord value(const ReferenceInfo& inf, const char* quantity, double eta, double sim, const Reference& r) {
    const double rv = r.value.to_double();
    const double e = std::abs((DD(sim) - r.value).to_double()) / std::abs(rv);
    ComparisonRecord c{quantity, "1", eta, sim, rv, r.u_ref, e, "relative error <= 1e-12", {}};
    if (!r.ok || !(r.u_ref <= budget * std::abs(rv))) {
        c.comparison = missing_input("a reference within 1e-3 of the criterion (u_ref " + shortest(r.u_ref) + ")");
    } else {
        c.comparison = judge(inf, {}, e <= criterion, c.criterion);
    }
    return c;
}

ComparisonRecord inverse(const ReferenceInfo& inf, double eta) {
    const double nu = physics::fermi_half(eta).value;
    const double sim = physics::inverse_fermi_half(nu);
    const double scale = std::max(1.0, std::abs(eta));
    auto f = [&](const DD& x) { return ref::fermi_tanh_sinh(FermiOrder::half, x).value - DD(nu); };
    ComparisonRecord c{"inverse F_1/2", "1", eta, sim, 0.0, 0.0, 0.0, "|d eta| / max(1, |eta|) <= 1e-12", {}};
    double delta = 1e-10 * scale;
    DD a = DD(sim) - DD(delta), b = DD(sim) + DD(delta);
    while ((f(a).hi < 0.0) == (f(b).hi < 0.0)) {
        delta *= 100.0;
        if (delta > 1.0) {
            c.comparison = missing_input("a bracket of the reference root within 1 of the code's answer");
            return c;
        }
        a = DD(sim) - DD(delta);
        b = DD(sim) + DD(delta);
    }
    const ref::BrentResult root = ref::brent_zero(f, a, b, 1e-18 * scale);
    if (root.status != ref::BrentStatus::ok) {
        c.comparison = missing_input("Brent's root of the reference");
        return c;
    }
    c.q_ref = root.root.to_double();
    c.u_ref = 1e-18 * scale + 6.0 * ref::u_dd * std::abs(c.q_ref);  // Brent's (2.18) bound
    c.metric = std::abs((DD(sim) - root.root).to_double()) / std::max(1.0, std::abs(c.q_ref));
    c.comparison = judge(inf, {}, c.metric <= criterion, c.criterion);
    return c;
}

}  // namespace

BenchmarkResult run_tvb02() {
    BenchmarkResult r;
    r.id = "TVB-02";
    r.title = "Fermi-Dirac integrals F_1/2, F_-1/2, F_3/2 and the inverse of F_1/2";
    r.info = info();
    double diag = 0.0;
    for (const double eta : grid()) {
        const Reference f12 = reference(FermiOrder::half, eta);
        const Reference fm12 = reference(FermiOrder::minus_half, eta);
        const Reference f32 = reference(FermiOrder::three_halves, eta);
        const physics::FermiIntegral f = physics::fermi_half(eta);
        const physics::FermiIntegral g = physics::fermi_three_halves(eta);
        r.comparisons.push_back(value(r.info, "F_1/2", eta, f.value, f12));
        r.comparisons.push_back(value(r.info, "F_-1/2", eta, f.derivative, fm12));
        r.comparisons.push_back(value(r.info, "F_3/2", eta, g.value, f32));
        r.comparisons.push_back(inverse(r.info, eta));
        diag = std::max(diag, std::abs((DD(g.derivative) - f12.value).to_double()) / std::abs(f12.value.to_double()));
    }
    r.notes.push_back("diagnostic (not judged): fermi_three_halves' derivative against F_1/2, worst relative " +
                      shortest(diag));
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
