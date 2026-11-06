// TVB-31.23-31.25, the WKB action, the Tsu-Esaki supply and current (Unit 55, TVR-3 step 8; DECISIONS 15.3, 15.9,
// 15.18).
//
// Code under test: physics::segment_action, wkb_exponent, wkb_transmission, tsu_esaki_supply and tsu_esaki_current
// (gate_tunnelling.hpp). Reference: verify/reference/model_tunnelling.hpp in double-double (model_tunnelling_test).
//   31.23 The segment action L g(a, b) and its partials over a, b in {-1 .. 5} eV (across 0, both ends inside, equal)
//         and near-equal pairs a (1 +- 10^-k); S(E), dS/dE and TC(E) of four barriers (rectangular, trapezoidal,
//         Fowler-Nordheim triangular, two-layer) over E -1 to 3.5 eV in 1/64 eV with every barrier end and its
//         neighbours. 15.9's value rule.
//   31.24 N(E) and its three partials at 300 K, x1 = (F1 - E) / kT from -64 to 64 in 1/8, for F1 - F2 in {0, 1e-9,
//         1e-3, 0.1, 1} eV; |x| > 60 not_applicable.
//   31.25 The current of the four barriers at 300 and 400 K and four biases against a double-double adaptive
//         integral of the same cut range: relative differences reported as diagnostics, never judged (MI-Q8).
#include <algorithm>
#include <array>
#include <cmath>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/physics/gate_tunnelling.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/model_tunnelling.hpp"
#include "tvb31_common.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

namespace mt = verify::reference::model_tunnelling;
using tvb31::DD;
using tvb31::in_range;
using tvb31::value_record;

constexpr std::string_view any_range = "any a, b, E (an exact closed form)";
constexpr std::string_view x_range = "(E - F) / kT within [-60, 60] for both levels";
constexpr std::string_view citation =
    "Gehring, Simulation of Tunneling in Semiconductor Devices, PhD thesis, TU Wien 2003, sec. 3.2, 3.3.1, 3.5.1 "
    "(read in step 1)";
constexpr std::string_view verification =
    "model_tunnelling_test: the action against tanh-sinh (TC4) split at the turning point, partials against central "
    "differences, the current against tanh-sinh; 15.9's rule";

constexpr std::string_view finding_xiv =
    "DECISIONS Unit 55 15.18 MI-Q22 (known finding): tsu_esaki_supply's dN/dE = sigma(x2) - sigma(x1) cancels where "
    "both levels are many kT above E (both sigma -> 1; an absolute error of eps); fix proposed separately (MI-Q10)";

struct TestBarrier {
    const char* name;
    std::vector<physics::BarrierSegment> segments;
};

std::vector<TestBarrier> barriers() {
    return {{"rectangular", {{3.1, 3.1, 1.5e-7, 0.5}}},
            {"trapezoidal", {{3.1, 2.1, 1.5e-7, 0.5}}},
            {"triangular", {{3.1, -0.9, 1.5e-7, 0.5}}},
            {"two-layer", {{3.1, 2.6, 1.0e-7, 0.5}, {1.5, 1.0, 2.0e-7, 0.18}}}};
}

// S, dS/dE and TC of one barrier with NS segments at E, against the reference with the condition sums over E and
// every segment's W0, W1, L and m.
template <std::size_t NS>
void judge_barrier(BenchmarkResult& r, const ReferenceInfo& inf, const TestBarrier& b, double E) {
    constexpr std::size_t N = 1 + 4 * NS;
    using AN = std::array<DD, N>;
    AN x{};
    x[0] = DD(E);
    for (std::size_t k = 0; k < NS; ++k) {
        const physics::BarrierSegment& s = b.segments[k];
        x[1 + 4 * k] = DD(s.W0_eV);
        x[2 + 4 * k] = DD(s.W1_eV);
        x[3 + 4 * k] = DD(s.length_cm);
        x[4 + 4 * k] = DD(s.mass_m0);
    }
    auto barrier_of = [](const AN& y) {
        std::array<mt::Segment, NS> seg{};
        for (std::size_t k = 0; k < NS; ++k) seg[k] = {y[1 + 4 * k], y[2 + 4 * k], y[3 + 4 * k], y[4 + 4 * k]};
        return seg;
    };
    auto fS = [&](const AN& y) {
        const auto seg = barrier_of(y);
        return mt::exponent(std::span<const mt::Segment>(seg), y[0]);
    };
    auto fT = [&](const AN& y) {
        const auto seg = barrier_of(y);
        return mt::transmission(std::span<const mt::Segment>(seg), y[0]);
    };
    double dE = 0.0;
    const double S = physics::wkb_exponent(b.segments, E, &dE);
    const double TC = physics::wkb_transmission(b.segments, E);
    const auto q = fS(x);
    const auto CS = verify::reference::condition_sums<2>(fS, x);
    const std::string ctx = std::string(" (") + b.name + ")";
    const std::vector<Predicate> pr{Predicate{std::string(any_range), 0.0, 0.0, true}};
    r.comparisons.push_back(value_record(inf, "S" + ctx, "1", E, S, q[0], CS[0], pr));
    r.comparisons.push_back(value_record(inf, "dS/dE" + ctx, "1/eV", E, dE, q[1], CS[1], pr));
    r.comparisons.push_back(value_record(inf, "TC" + ctx, "1", E, TC, fT(x), verify::reference::condition_sum(fT, x),
                                         pr));
}

}  // namespace

BenchmarkResult run_tvb31_23(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-31.23";
    r.title = "WKB segment action, exponent and transmission";
    r.info = tvb31::entry_info(r.id, Provenance::PA, std::string(citation),
                               "L g(a, b), g = (2/3)(a+^(3/2) - b+^(3/2)) / (a - b); S = sum sqrt(2 m q) / hbar L g; "
                               "TC = exp(-2 S)",
                               {std::string(any_range)}, {12}, "verify/reference/model_tunnelling.hpp",
                               std::string(verification));
    r.info.cls = ReferenceClass::E;
    r.info.second_method = "tanh-sinh quadrature (TC4) of the root over the segment, split at the turning point (L0)";
    const ReferenceInfo inf = r.info;
    const std::vector<Predicate> pr{Predicate{std::string(any_range), 0.0, 0.0, true}};
    // The segment action.
    const std::array<double, 13> v{-1.0, -1e-3, -1e-12, 0.0, 1e-12, 1e-6, 1e-3, 0.1, 0.5, 1.0, 2.0, 3.1, 5.0};
    std::vector<std::array<double, 2>> pairs;
    for (const double a : v) {
        for (const double b : v) pairs.push_back({a, b});
    }
    for (const double a : {1e-3, 0.5, 3.1}) {
        for (int k = 1; k <= 15; ++k) {
            const double d = std::pow(10.0, -k);
            pairs.push_back({a, a * (1.0 + d)});
            pairs.push_back({a, a * (1.0 - d)});
            pairs.push_back({a * (1.0 + d), a});
        }
    }
    const double L = 1.5e-7;
    using A3 = std::array<DD, 3>;
    auto fa = [](const A3& y) { return mt::segment_action(y[0], y[1], y[2]); };
    for (const auto& p : pairs) {
        const physics::SegmentAction s = physics::segment_action(p[0], p[1], L);
        const A3 x{DD(p[0]), DD(p[1]), DD(L)};
        const auto q = fa(x);
        const auto S = verify::reference::condition_sums<3>(fa, x);
        const std::string ctx = " (a " + shortest(p[0]) + ", b " + shortest(p[1]) + ")";
        r.comparisons.push_back(value_record(inf, "L g" + ctx, "cm eV^1/2", p[0] - p[1], s.value, q[0], S[0], pr));
        r.comparisons.push_back(value_record(inf, "L dg/da" + ctx, "cm eV^-1/2", p[0] - p[1], s.d_a, q[1], S[1], pr));
        r.comparisons.push_back(value_record(inf, "L dg/db" + ctx, "cm eV^-1/2", p[0] - p[1], s.d_b, q[2], S[2], pr));
    }
    // The scale.
    r.comparisons.push_back(value_record(inf, "sqrt(2 m0 q) / hbar", "1/(cm eV^1/2)", 0.0, physics::wkb_action_scale(),
                                         mt::action_scale(), 0.0, pr));
    // Barriers.
    for (const TestBarrier& b : barriers()) {
        std::vector<double> Es;
        for (int i = -64; i <= 224; ++i) {
            if (tvb31::keep(i + 64, 288, tier)) Es.push_back(i / 64.0);
        }
        for (const auto& s : b.segments) {
            for (const double w : {s.W0_eV, s.W1_eV}) {
                for (const double e : {w, std::nextafter(w, -INFINITY), std::nextafter(w, INFINITY)}) Es.push_back(e);
            }
        }
        for (const double E : Es) {
            if (b.segments.size() == 1) {
                judge_barrier<1>(r, inf, b, E);
            } else {
                judge_barrier<2>(r, inf, b, E);
            }
        }
    }
    apply_known_findings(r, known_findings());
    return r;
}

BenchmarkResult run_tvb31_24(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-31.24";
    r.title = "Tsu-Esaki supply function and its partials";
    r.info = tvb31::entry_info(r.id, Provenance::PA, std::string(citation),
                               "N(E) = kT [ln(1 + e^((F1 - E)/kT)) - ln(1 + e^((F2 - E)/kT))]", {std::string(x_range)},
                               {12}, "verify/reference/model_tunnelling.hpp", std::string(verification));
    const ReferenceInfo inf = r.info;
    const double kT = base::k_B_eV_per_K * 300.0, F1 = 0.2;
    using A4 = std::array<DD, 4>;
    auto fs = [](const A4& y) { return mt::supply(y); };
    const std::array<std::string, 4> names{"N", "dN/dF1", "dN/dF2", "dN/dE"};
    const std::array<std::string, 4> units{"eV", "1", "1", "1"};
    tvb31::FloorFinding xiv_floor{"(xiv)", finding_xiv};
    for (const double dF : {0.0, 1e-9, 1e-3, 0.1, 1.0}) {
        const double F2 = F1 - dF;
        for (int i = -512; i <= 512; ++i) {
            if (!tvb31::keep(i + 512, 1024, tier)) continue;
            const double E = F1 - (i / 8.0) * kT;
            const physics::Supply s = physics::tsu_esaki_supply(E, F1, F2, kT);
            const A4 x{DD(E), DD(F1), DD(F2), DD(kT)};
            const auto q = fs(x);
            const auto S = verify::reference::condition_sums<4>(fs, x);
            const double x1 = (F1 - E) / kT, x2 = (F2 - E) / kT;
            Predicate p = in_range(std::string(x_range), std::abs(x1) > std::abs(x2) ? x1 : x2, -60.0, 60.0);
            p.holds = std::abs(x1) <= 60.0 && std::abs(x2) <= 60.0;
            const std::array<double, 4> sim{s.value, s.d_F1, s.d_F2, s.d_E};
            const std::string ctx = " (F1 - F2 " + shortest(dF) + ")";
            for (std::size_t k = 0; k < 4; ++k) {
                r.comparisons.push_back(value_record(inf, names[k] + ctx, units[k], x1, sim[k], q[k], S[k], {p}));
            }
            // Finding (xiv)'s floor: the shipped dN/dE = sigma(x2) - sigma(x1), terms of size sigma each (MI-Q22).
            tvb31::classify_floor<4>(r, {0.0, 0.0, 0.0, std::abs(s.d_F1) + std::abs(s.d_F2)}, xiv_floor);
        }
    }
    apply_known_findings(r, known_findings());
    tvb31::note_if_stale(r, xiv_floor);
    return r;
}

BenchmarkResult run_tvb31_25(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-31.25";
    r.title = "Tsu-Esaki current quadrature (diagnostic, MI-Q8)";
    r.info = tvb31::entry_info(r.id, Provenance::PA, std::string(citation),
                               "J = (A* / k^2) int TC N dE over the cut range; NiTCAD's composite 4-point "
                               "Gauss-Legendre against a double-double adaptive Gauss-Kronrod integral",
                               {"MI-Q8: reported, not judged"}, {12}, "verify/reference/model_tunnelling.hpp",
                               std::string(verification));
    r.info.cls = ReferenceClass::E;
    r.info.second_method = "tanh-sinh quadrature of the same integrand (L0)";
    double worst = 0.0;
    std::string worst_at;
    // The four barriers of 31.23, and one whose ends (0.8 and 0.3 eV) lie inside the integration range, so that the
    // turning point crosses the segment end there: the kinks of MI-Q8 (added after the first run, which had none).
    std::vector<TestBarrier> bs = barriers();
    bs.push_back({"kinked", {{0.8, 0.3, 1.0e-7, 0.4}}});
    for (const TestBarrier& b : bs) {
        std::vector<mt::Segment> seg;
        for (const auto& s : b.segments) seg.push_back({DD(s.W0_eV), DD(s.W1_eV), DD(s.length_cm), DD(s.mass_m0)});
        for (const double T : {300.0, 400.0}) {
            if (tier == Tier::fast && T != 300.0) continue;  // MI-Q9 (15.20): 300 K only in fast
            for (const double V : {0.01, 0.1, 0.5, 1.0}) {
                if (tier == Tier::fast && (V == 0.01 || V == 0.5)) continue;  // and F1 - F2 0.1 and 1 eV
                const double F1 = 0.1, F2 = F1 - V;
                const physics::TunnelCurrent J =
                    physics::tsu_esaki_current(b.segments, {0.0, INFINITY}, F1, F2, T, 120.0);
                const mt::Current q = mt::current(seg, 0.0, INFINITY, F1, F2, T, 120.0);
                if (!q.ok) {
                    r.run_error = std::string("reference current did not converge (") + b.name + ")";
                    return r;
                }
                const double rel = std::abs((DD(J.J) - q.J).hi) / std::abs(q.J.hi);
                const std::string at =
                    std::string(b.name) + ", T " + shortest(T) + " K, F1 - F2 " + shortest(V) + " eV";
                r.notes.push_back("diagnostic (not judged, MI-Q8): J " + shortest(J.J) + " A/cm^2, reference " +
                                  shortest(q.J.to_double()) + ", relative difference " + shortest(rel) + " (" + at +
                                  ")");
                if (rel > worst) {
                    worst = rel;
                    worst_at = at;
                }
            }
        }
    }
    r.notes.push_back("diagnostic (not judged, MI-Q8): worst relative difference " + shortest(worst) + " (" + worst_at +
                      ")");
    return r;
}

}  // namespace NiTCAD::verify::benchmark
