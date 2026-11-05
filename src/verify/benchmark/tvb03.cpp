// TVB-03, the Bernoulli function (DECISIONS Unit 55 11.8; ruling TB-Q9 in 11.13).
//
// Code under test: assemble::bernoulli (header-only). Reference: verify/reference/bernoulli.hpp. Criterion (frozen):
// <= 4 eps relative wherever B(x) is a normal double (the applicability predicate), on a logarithmic grid of
// +-[1e-300, 745] (20 points a decade), x = 0, +-1e-2 and their neighbouring doubles, every 0.05 on [700, 745].
// Finding (iii) (B = 0 where std::expm1 overflows) was fixed in Unit 56 (DECISIONS Unit 55 12, FQ4).
#include <cfloat>
#include <cmath>
#include <limits>
#include <vector>

#include "NiTCAD/assemble/bernoulli.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/bernoulli.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

namespace ref = verify::reference;

ReferenceInfo info() {
    ReferenceInfo r;
    r.id = "TVB-03";
    r.cls = ReferenceClass::E;
    r.independence = Independence::I3;
    r.provenance = Provenance::PA;
    r.citation =
        "definition B(x) = x / (e^x - 1); Bernoulli numbers B_2k, OEIS A000367 and A002445 b-files "
        "(oeis.org, read 2026-11-05)";
    r.version = "OEIS b-files n = 1..17";
    r.equations = {3};
    r.form = "x / expm1(x) in double-double; scaled x e^-r 2^-m above 670; -x below -670";
    r.assumptions = {"none: B is defined by its formula"};
    r.applicability = {"B(x) a normal double"};
    r.derivation = "verify/reference/bernoulli.hpp";
    r.verification = "TC1 (expm1, exp, ln 2); series vs direct; B(-x) = B(x) + x; scaled vs direct";
    r.second_method = "the generating-function series with OEIS's exact B_2k, |x| <= 1/2";
    return r;
}

}  // namespace

static std::vector<double> tvb03_grid() {
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

BenchmarkResult run_tvb03() {
    BenchmarkResult r;
    r.id = "TVB-03";
    r.title = "The Bernoulli function B(x) = x / (e^x - 1)";
    r.info = info();
    const ReferenceInfo inf = r.info;
    constexpr double criterion = 4.0 * std::numeric_limits<double>::epsilon();
    for (const double x : tvb03_grid()) {
        const ref::ScaledDD b = ref::bernoulli_direct(x);
        const bool normal = ref::is_normal_double(b);
        const Predicate pred[] = {{"B(x) a normal double", normal ? 1.0 : 0.0, 1.0, normal}};
        const double sim = assemble::bernoulli(x);
        const double e = normal ? ref::relative_error(sim, b) : 0.0;
        ComparisonRecord c{"B(x)", "1", x, sim, normal ? std::ldexp(b.mantissa.to_double(), b.exponent) : 0.0,
                           0.0, e, "relative error <= 4 eps", {}};
        c.comparison = judge(inf, pred, e <= criterion, c.criterion);
        r.comparisons.push_back(std::move(c));
    }
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
