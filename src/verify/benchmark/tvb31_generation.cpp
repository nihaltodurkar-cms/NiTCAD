// TVB-31.20-31.22, impact ionization and band-to-band generation (Unit 55, TVR-3 step 7; DECISIONS 15.3, 15.9, 15.17).
//
// Code under test: physics::impact_ionization_coefficient and impact_ionization_temperature_factor
// (impact_ionization.hpp) with its slope (thermal.hpp); kane_generation and the direct-gap WKB functions
// (band_to_band.hpp). Reference: verify/reference/model_generation.hpp in double-double (model_generation_test).
// Criterion: 15.9's value rule; the documented cuts (gamma B / E > 700, B / F > 700) by 15.9's cut rule.
//   31.20 vOdM, both carriers and every branch, the register's coefficients: E 0 and 1e4-1e7 V/cm, the hole switch and
//         its neighbours, the cut field gamma B / 700 with its neighbours and (1 +- 10^-j) points, at gamma(200, 300,
//         500 K); gamma and d gamma / dT over 200-500 K.
//   31.21 Kane, the shipped (A, B) (grade L, F5): F 0 and 1e4-1e7 V/cm, the cut B / 700 and its neighbours.
//   31.22 The direct-gap WKB functions on two test bands (u = 14.9: E_g 1.42 eV, m_r 0.0335 m0; u = 1.5: 0.5 eV,
//         m0 / 3; no shipped material has direct-gap masses): kappa and both antiderivatives over delta in [0, 1] in
//         1/1024 and toward each edge; segments between ten endpoints in and around the gap and near-coincident
//         pairs D = 10^-k; the path rate over slope, I_k, I_ik and z = k_m^2 I_ik; eq (8) over F 1e7-1e10 V/m.
#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/physics/band_to_band.hpp"
#include "NiTCAD/physics/impact_ionization.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/thermal.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/model_generation.hpp"
#include "tvb31_common.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

namespace mg = verify::reference::model_generation;
using tvb31::classify_floor;
using tvb31::cut_record;
using tvb31::FloorFinding;
using tvb31::note_if_stale;
using tvb31::DD;
using tvb31::in_range;
using tvb31::register_dd;
using tvb31::value_record;

constexpr std::string_view e_range = "E within [1e4, 1e7] V/cm, or 0, or a cut or switch point (15.5)";
constexpr std::string_view t_range = "T within [200, 500] K";
constexpr std::string_view d_range = "delta within [0, 1]";
constexpr std::string_view f_range = "F within [1e7, 1e10] V/m with the exponent at most 700 (a normal rate)";
constexpr std::string_view path_range = "any (the path rate is a formula in its four inputs)";
constexpr std::string_view seg_range = "any endpoints (the integrals count only the part inside the gap, by design)";
constexpr std::string_view verification = "model_generation_test: derivatives against complex step (TC8) and "
                                          "double-double central differences, the antiderivatives against tanh-sinh "
                                          "(TC4), eq (8) against eq (11) on the exact uniform path; 15.9's rule";

constexpr std::string_view finding_xi =
    "DECISIONS Unit 55 15.17 MI-Q18 (known finding): wkb_segment's integrals and partials are differences of two "
    "antiderivatives of size ~u pi / 4 divided by D (and D^2), so a short segment loses eps K / (kappa D) relative; "
    "and below D = 1e-10 its midpoint rule truncates (L f'' D / 12 in the partials); fix proposed separately (MI-Q10)";
constexpr std::string_view finding_xii =
    "DECISIONS Unit 55 15.17 MI-Q18 (known finding): wkb_path_rate's dG/dI_ik = pre t (k_m^2 e / I_ik - (1 - e) / "
    "I_ik^2) cancels as z = k_m^2 I_ik -> 0 (relative eps / z^2); fix proposed separately (MI-Q10)";

// The size of the terms of NiTCAD's two antiderivatives at delta (clipped to [0, 1]), for the rounding floor of
// finding (xi): sqrt(m_r E_g) / (2 u hbar) (u (|theta| + |alpha c|) / 2 + c^3 / 3) and hbar / (2 u sqrt(m_r E_g))
// (u |theta| + c).
struct TermSizes {
    double K, G;
};
TermSizes term_sizes(double delta, const mg::Band& b) {
    const DD u = mg::u_of(b), s = mg::scale_of(b), h = mg::hbar();
    const mg::AlphaCos a = mg::alpha_cos(DD(std::clamp(delta, 0.0, 1.0)), u);
    const double theta = std::abs(std::atan2(a.alpha.hi, a.c.hi)), ac = std::abs(a.alpha.hi * a.c.hi), c = a.c.hi;
    const double kk = (s / (DD(2.0) * u * h)).hi, gg = (h / (DD(2.0) * u * s)).hi;
    return {kk * (0.5 * u.hi * (theta + ac) + c * c * c / 3.0), gg * (u.hi * theta + c)};
}

// The field grid of 31.20 and 31.21: 0, the log grid, `specials` with their neighbours, the cut threshold t with its
// neighbours and t (1 +- 10^-j), and 1e8 beyond.
struct FieldPoint {
    double E;
    bool special;  // a cut or switch point (judged whatever the range)
};
std::vector<FieldPoint> field_grid(int per, std::initializer_list<double> specials, double threshold) {
    std::vector<FieldPoint> g{{0.0, true}};
    for (const double E : tvb31::log_grid(1e4, 1e7, per)) g.push_back({E, false});
    for (const double s : specials) {
        for (const double v : {s, std::nextafter(s, 0.0), std::nextafter(s, INFINITY)}) g.push_back({v, true});
    }
    for (const double v : {threshold, std::nextafter(threshold, 0.0), std::nextafter(threshold, INFINITY)}) {
        g.push_back({v, true});
    }
    for (int j = 1; j <= 12; ++j) {
        g.push_back({threshold * (1.0 + std::pow(10.0, -j)), true});
        g.push_back({threshold * (1.0 - std::pow(10.0, -j)), true});
    }
    g.push_back({1e8, false});
    return g;
}

Predicate field_predicate(double E, bool special) {
    Predicate p = in_range(std::string(e_range), E, 1e4, 1e7);
    if (special) p.holds = true;
    return p;
}

}  // namespace

BenchmarkResult run_tvb31_20(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-31.20";
    r.title = "van Overstraeten-de Man impact ionization and its temperature factor";
    r.info = tvb31::entry_info(r.id, Provenance::PB,
                               "the form alpha = A exp(-B / E) of van Overstraeten and de Man, Solid-State Electron. 13, "
                               "583 (1970), as Allpix Squared 3.3 and Sentaurus print it (MI-Q17); gamma(T) the tanh "
                               "ratio (Sentaurus' form, tag C: the 300 K pin a register row)",
                               "alpha = gamma A exp(-gamma B / E); gamma = tanh(hw / 2kT0) / tanh(hw / 2kT), T0 = 300 K",
                               {std::string(e_range), std::string(t_range)}, {8},
                               "verify/reference/model_generation.hpp", std::string(verification));
    const ReferenceInfo inf = r.info;
    ReferenceInfo ginf = r.info;
    ginf.applicability = {std::string(t_range)};
    ginf.provenance = Provenance::PD;  // a NiTCAD choice (tag C, F4): its documented form, checked as written
    ginf.citation = "the tanh ratio as NiTCAD documents it (Sentaurus' vanOverstraetendeMan form; tag C, F4)";
    auto made = physics::Semiconductor::create(
        tvb31::register_parameters(physics::legacy_silicon_parameters, "legacy_silicon", {"31.20"}));
    if (!made) {
        r.run_error = made.error().message;
        return r;
    }
    const physics::Semiconductor& m = *made;
    const double hw = m.parameters().impact_ionization.phonon_energy_eV;
    const int per = tvb31::per_decade(tier);
    // gamma and its slope.
    std::vector<double> Ts = tvb31::log_grid(200.0, 500.0, per);
    Ts.push_back(300.0);
    Ts.push_back(150.0);
    Ts.push_back(600.0);
    for (const double T : Ts) {
        const std::vector<Predicate> pr{in_range(std::string(t_range), T, 200.0, 500.0)};
        using A2 = std::array<DD, 2>;
        auto f = [](const A2& y) { return mg::vodm_gamma(y); };
        const A2 x{DD(hw), DD(T)};
        const A2 q = f(x);
        const auto S = verify::reference::condition_sums<2>(f, x);
        r.comparisons.push_back(value_record(ginf, "gamma", "1", T, physics::impact_ionization_temperature_factor(hw, T),
                                             q[0], S[0], pr));
        r.comparisons.push_back(value_record(ginf, "dgamma/dT", "1/K", T,
                                             physics::impact_ionization_temperature_factor_slope(hw, T), q[1], S[1],
                                             pr));
    }
    // alpha and its slope.
    for (const bool electron : {true, false}) {
        const auto carrier = electron ? physics::Carrier::electron : physics::Carrier::hole;
        const physics::ImpactIonizationCoefficients& c = physics::impact_ionization(m, carrier);
        for (const double T : {200.0, 300.0, 500.0}) {
            const double gamma = physics::impact_ionization_temperature_factor(hw, T);
            // The cut threshold of the low branch (the one in use at small fields).
            const double threshold = gamma * c.B_low_V_per_cm / 700.0;
            for (const FieldPoint& p : field_grid(per, {c.switch_V_per_cm}, threshold)) {
                const bool high = p.E >= c.switch_V_per_cm;
                const double A = high ? c.A_high_per_cm : c.A_low_per_cm, B = high ? c.B_high_V_per_cm : c.B_low_V_per_cm;
                const physics::ImpactIonizationRate sim = physics::impact_ionization_coefficient(c, gamma, p.E);
                const std::string ctx = std::string(electron ? " electron" : " hole") + " (T " + shortest(T) + ")";
                const std::vector<Predicate> pr{field_predicate(p.E, p.special),
                                                in_range(std::string(t_range), T, 200.0, 500.0)};
                using A4 = std::array<DD, 4>;
                auto f = [](const A4& y) { return mg::vodm(y); };
                if (p.E == 0.0) {
                    r.comparisons.push_back(value_record(inf, "alpha" + ctx, "1/cm", 0.0, sim.alpha, DD(0.0), 0.0, pr));
                    r.comparisons.push_back(value_record(inf, "dalpha/dE" + ctx, "1/V", 0.0, sim.d_dE, DD(0.0), 0.0, pr));
                    continue;
                }
                const A4 x{DD(A), DD(B), DD(gamma), DD(p.E)};
                const auto q = f(x);
                if (gamma * B / p.E > 700.0) {  // the documented cut
                    const auto qt = f(A4{DD(A), DD(B), DD(gamma), DD(gamma * B / 700.0)});
                    r.comparisons.push_back(cut_record(inf, "alpha" + ctx, "1/cm", p.E, sim.alpha, q[0], qt[0], pr));
                    r.comparisons.push_back(cut_record(inf, "dalpha/dE" + ctx, "1/V", p.E, sim.d_dE, q[1], qt[1], pr));
                    continue;
                }
                const auto S = verify::reference::condition_sums<2>(f, x);
                r.comparisons.push_back(value_record(inf, "alpha" + ctx, "1/cm", p.E, sim.alpha, q[0], S[0], pr));
                r.comparisons.push_back(value_record(inf, "dalpha/dE" + ctx, "1/V", p.E, sim.d_dE, q[1], S[1], pr));
            }
        }
    }
    tvb31::table_checks(r, "31.20");
    apply_known_findings(r, known_findings());
    return r;
}

BenchmarkResult run_tvb31_21(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-31.21";
    r.title = "Kane band-to-band generation";
    r.info = tvb31::entry_info(r.id, Provenance::PB,
                               "the Kane form G = A F^2 exp(-B / F) as Hurkx, Klaassen and Knuvers 1992 and Sentaurus' "
                               "public summaries print it; the shipped (A, B) unattributed (F5, grade L)",
                               "G = A F^2 exp(-B / F); dG/dF = G (2 + B / F) / F", {std::string(e_range)}, {8},
                               "verify/reference/model_generation.hpp", std::string(verification));
    r.info.assumptions.push_back("A, B: the shipped values (grade L, F5)");
    const ReferenceInfo inf = r.info;
    const physics::BandToBandParameters& b = physics::legacy_silicon_parameters.band_to_band;
    const double threshold = b.B_V_per_cm / 700.0;
    for (const FieldPoint& p : field_grid(tvb31::per_decade(tier), {}, threshold)) {
        const physics::KaneRate sim = physics::kane_generation(b.A_per_cm3_s, b.B_V_per_cm, p.E);
        const std::vector<Predicate> pr{field_predicate(p.E, p.special)};
        using A3 = std::array<DD, 3>;
        auto f = [](const A3& y) { return mg::kane(y); };
        if (p.E == 0.0) {
            r.comparisons.push_back(value_record(inf, "G", "cm^-3 s^-1", 0.0, sim.rate, DD(0.0), 0.0, pr));
            r.comparisons.push_back(value_record(inf, "dG/dF", "cm^-3 s^-1 cm/V", 0.0, sim.d_dF, DD(0.0), 0.0, pr));
            continue;
        }
        const A3 x{DD(b.A_per_cm3_s), DD(b.B_V_per_cm), DD(p.E)};
        const auto q = f(x);
        if (b.B_V_per_cm / p.E > 700.0) {
            const auto qt = f(A3{x[0], x[1], DD(threshold)});
            r.comparisons.push_back(cut_record(inf, "G", "cm^-3 s^-1", p.E, sim.rate, q[0], qt[0], pr));
            r.comparisons.push_back(cut_record(inf, "dG/dF", "cm^-3 s^-1 cm/V", p.E, sim.d_dF, q[1], qt[1], pr));
            continue;
        }
        const auto S = verify::reference::condition_sums<2>(f, x);
        r.comparisons.push_back(value_record(inf, "G", "cm^-3 s^-1", p.E, sim.rate, q[0], S[0], pr));
        r.comparisons.push_back(value_record(inf, "dG/dF", "cm^-3 s^-1 cm/V", p.E, sim.d_dF, q[1], S[1], pr));
    }
    tvb31::table_checks(r, "31.21");
    apply_known_findings(r, known_findings());
    return r;
}

BenchmarkResult run_tvb31_22(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-31.22";
    r.title = "Direct-gap WKB band-to-band functions";
    r.info = tvb31::entry_info(r.id, Provenance::PA,
                               "Esseni, Pala, Palestri, Alper and Rollo, Semicond. Sci. Technol. 32, 083005 (2017), "
                               "sec. 2.1, eqs (8), (9), (11), (12) (read in step 1)",
                               "kappa(delta), its antiderivatives, segment integrals and partials, the path rate (11) "
                               "and partials, the uniform-field rate (8)",
                               {std::string(d_range), std::string(seg_range), std::string(f_range),
                                std::string(path_range)},
                               {8, 12},
                               "verify/reference/model_generation.hpp", std::string(verification));
    r.info.assumptions.push_back("test bands (no shipped material has direct-gap masses); NiTCAD's base constants");
    ReferenceInfo dinf = r.info;
    dinf.applicability = {std::string(d_range)};
    ReferenceInfo finf = r.info;
    finf.applicability = {std::string(f_range)};
    ReferenceInfo pinf = r.info;
    pinf.applicability = {std::string(path_range)};
    ReferenceInfo sinf = r.info;
    sinf.applicability = {std::string(seg_range)};
    FloorFinding xi_floor{"(xi)", finding_xi};
    FloorFinding xii_floor{"(xii)", finding_xii};
    struct TestBand {
        const char* name;
        double gap_eV, mr_m0;
    };
    for (const TestBand& tb : {TestBand{"u 14.9", 1.42, 0.0335}, TestBand{"u 1.5", 0.5, 1.0 / 3.0}}) {
        const physics::WkbBand band{tb.gap_eV * base::q_C, tb.mr_m0 * base::m0_kg};
        const std::string bn = std::string(" (") + tb.name + ")";
        auto ref_band = [](const DD& g, const DD& mr) { return mg::Band{g, mr}; };
        // kappa and the antiderivatives.
        std::vector<double> ds;
        for (int i = 0; i <= 1024; ++i) ds.push_back(i / 1024.0);
        for (int k = 1; k <= 15; ++k) {
            ds.push_back(std::pow(10.0, -k));
            ds.push_back(1.0 - std::pow(10.0, -k));
        }
        ds.push_back(-0.5);
        ds.push_back(1.5);
        using A3 = std::array<DD, 3>;  // delta, gap, m_r
        const A3 base_x{DD(0.0), DD(band.gap_J), DD(band.reduced_mass_kg)};
        for (const double d : ds) {
            const std::vector<Predicate> pr{in_range(std::string(d_range), d, 0.0, 1.0)};
            A3 x = base_x;
            x[0] = DD(d);
            auto fk = [&](const A3& y) { return mg::kappa(y[0], ref_band(y[1], y[2])); };
            auto fK = [&](const A3& y) { return mg::kappa_antiderivative(y[0], ref_band(y[1], y[2])); };
            auto fG = [&](const A3& y) { return mg::inverse_kappa_antiderivative(y[0], ref_band(y[1], y[2])); };
            r.comparisons.push_back(value_record(dinf, "kappa" + bn, "1/m", d, physics::wkb_kappa(d, band), fk(x),
                                                 verify::reference::condition_sum(fk, x), pr));
            r.comparisons.push_back(value_record(dinf, "int kappa" + bn, "1/m", d,
                                                 physics::wkb_kappa_antiderivative(d, band), fK(x),
                                                 verify::reference::condition_sum(fK, x), pr));
            r.comparisons.push_back(value_record(dinf, "int 1/kappa" + bn, "m", d,
                                                 physics::wkb_inverse_kappa_antiderivative(d, band), fG(x),
                                                 verify::reference::condition_sum(fG, x), pr));
        }
        // Segments.
        struct Pair {
            double a, b;
        };
        std::vector<Pair> pairs;
        const std::array<double, 10> ends{-0.2, 0.0, 1e-9, 0.05, 0.3, 0.5, 0.75, 0.999, 1.0, 1.3};
        for (const double a : ends) {
            for (const double b : ends) {
                if (a != b) pairs.push_back({a, b});
            }
        }
        for (const double base : {1e-6, 0.3, 0.6, 0.999999}) {
            for (int k = 3; k <= 14; ++k) {
                pairs.push_back({base, base + std::pow(10.0, -k)});
                pairs.push_back({base + std::pow(10.0, -k), base});
            }
        }
        const double L = 1e-8;
        const std::array<std::string, 6> sn{"I_k", "I_ik", "dI_k/da", "dI_k/db", "dI_ik/da", "dI_ik/db"};
        const std::array<std::string, 6> su{"1", "m^2", "1", "1", "m^2", "m^2"};
        using A5 = std::array<DD, 5>;  // da, db, L, gap, m_r
        auto fs = [&](const A5& y) { return mg::segment(y[0], y[1], y[2], ref_band(y[3], y[4])).q; };
        for (const Pair& p : pairs) {
            const physics::WkbSegment s = physics::wkb_segment(p.a, p.b, L, band);
            const A5 x{DD(p.a), DD(p.b), DD(L), DD(band.gap_J), DD(band.reduced_mass_kg)};
            const mg::Segment q = mg::segment(x[0], x[1], x[2], ref_band(x[3], x[4]));
            if (!q.ok) {
                r.run_error = "reference segment failed at " + shortest(p.a) + ", " + shortest(p.b);
                return r;
            }
            const auto S = verify::reference::condition_sums<6>(fs, x);
            const std::array<double, 6> sim{s.I_k, s.I_ik, s.dI_k_da, s.dI_k_db, s.dI_ik_da, s.dI_ik_db};
            const std::string ctx = " [" + shortest(p.a) + ", " + shortest(p.b) + "]" + bn;
            for (std::size_t k = 0; k < 6; ++k) {
                Predicate pr{std::string(seg_range), 0.0, 0.0, true};
                r.comparisons.push_back(value_record(sinf, sn[k] + ctx, su[k], p.b - p.a, sim[k], q.q[k], S[k], {pr}));
            }
            // Finding (xi)'s floor, from the shipped form L (K(hi) - K(lo)) / D and L (f(v) -+ dK / D) / D.
            const double D = std::abs(p.b - p.a);
            if (D >= 1e-10) {
                const bool rising = p.b >= p.a;
                const double vmin = rising ? p.a : p.b, vmax = rising ? p.b : p.a;
                const double lo = std::clamp(vmin, 0.0, 1.0), hi = std::clamp(vmax, 0.0, 1.0);
                const bool in_lo = vmin > 0.0 && vmin < 1.0, in_hi = vmax > 0.0 && vmax < 1.0;
                const mg::Band rb = ref_band(x[3], x[4]);
                const TermSizes tl = term_sizes(lo, rb), th = term_sizes(hi, rb);
                const double kl = in_lo ? physics::wkb_kappa(lo, band) : 0.0;
                const double kh = in_hi ? physics::wkb_kappa(hi, band) : 0.0;
                const double il = in_lo ? 1.0 / kl : 0.0, ih = in_hi ? 1.0 / kh : 0.0;
                const double TK = L * (tl.K + th.K) / D, TG = L * (tl.G + th.G) / D;
                const double dk_lo = (L * kl + TK) / D, dk_hi = (L * kh + TK) / D;
                const double di_lo = (L * il + TG) / D, di_hi = (L * ih + TG) / D;
                classify_floor<6>(r,
                                  {TK, TG, rising ? dk_lo : dk_hi, rising ? dk_hi : dk_lo, rising ? di_lo : di_hi,
                                   rising ? di_hi : di_lo},
                                  xi_floor);
            } else {
                // The shipped D < 1e-10 branch, L f(mid) with the partials split equally: its truncation is
                // L |f''| D^2 / 24 in the value and L |f''| D / 12 in each partial (I = L avg f, avg f = f(m) +
                // f''(m) D^2 / 24 + ...). Finding (xi) where the miss is within twice that plus 64 eps |q_ref|.
                const mg::Band rb = ref_band(x[3], x[4]);
                const DD m = (x[0] + x[1]) / DD(2.0);
                const DD h = DD(1e-4 * std::min(m.hi, 1.0 - m.hi));
                auto second = [&](auto&& f) { return std::abs(((f(m + h) - DD(2.0) * f(m) + f(m - h)) / (h * h)).hi); };
                const double k2 = second([&](const DD& v) { return mg::kappa(v, rb); });
                const double i2 = second([&](const DD& v) { return DD(1.0) / mg::kappa(v, rb); });
                const std::array<double, 6> trunc{L * k2 * D * D / 24.0, L * i2 * D * D / 24.0, L * k2 * D / 12.0,
                                                  L * k2 * D / 12.0,     L * i2 * D / 12.0,     L * i2 * D / 12.0};
                for (std::size_t k = 0; k < 6; ++k) {
                    ComparisonRecord& c = r.comparisons[r.comparisons.size() - 6 + k];
                    if (c.comparison.outcome != Outcome::fail) continue;
                    const double bound = 2.0 * trunc[k] + 64.0 * tvb31::eps * std::abs(c.q_ref);
                    if (!(std::abs(c.q_sim - c.q_ref) <= bound)) continue;
                    c.comparison.outcome = Outcome::finding;
                    c.comparison.reason = "known finding (xi) (" + std::string(finding_xi) + "): " +
                                          c.comparison.reason + "; within the midpoint rule's truncation " +
                                          shortest(bound);
                    ++xi_floor.converted;
                }
            }
        }
        // eq (8).
        std::vector<double> Fs = tvb31::log_grid(1e7, 1e10, tvb31::per_decade(tier));
        Fs.push_back(1e6);
        using F3 = std::array<DD, 3>;  // F, gap, m_r
        auto f8 = [&](const F3& y) { return mg::uniform_field_rate(y[0], ref_band(y[1], y[2])); };
        const DD expo_coef = DD(std::acos(-1.0)) * sqrt(DD(band.reduced_mass_kg)) * DD(band.gap_J) *
                             sqrt(DD(band.gap_J)) / (DD(2.0) * DD(base::q_C) * DD(base::hbar_J_s));
        for (const double F : Fs) {
            Predicate pr = in_range(std::string(f_range), F, 1e7, 1e10);
            pr.holds = pr.holds && (expo_coef / DD(F)).hi <= 700.0;
            const F3 x{DD(F), DD(band.gap_J), DD(band.reduced_mass_kg)};
            r.comparisons.push_back(value_record(finf, "eq (8) G" + bn, "m^-3 s^-1", F,
                                                 physics::wkb_uniform_field_rate(F, band), f8(x),
                                                 verify::reference::condition_sum(f8, x), {pr}));
        }
    }
    // The path rate.
    const std::array<std::string, 5> pn{"G", "dG/dslope", "dG/dI_k", "dG/dI_ik", "dG/dk_m^2"};
    const std::array<std::string, 5> pu{"m^-3 s^-1", "m^-2 s^-1 / J", "m^-3 s^-1", "m^-5 s^-1", "m^-1 s^-1"};
    using A4 = std::array<DD, 4>;
    auto fp = [](const A4& y) { return mg::path_rate(y); };
    for (const double F : {-1e8, 1e7, 1e8, 1e9}) {
        for (const double Ik : {0.0, 0.1, 1.0, 10.0, 100.0, 340.0}) {
            for (const double Iik : {0.0, 1e-20, 1e-17, 1e-14}) {
                for (const double z : {1e-10, 1e-4, 0.5, 1.0, 3.0, 50.0, 800.0}) {
                    const double slope = F * base::q_C, km2 = Iik > 0.0 ? z / Iik : z;
                    const physics::WkbPathRate s = physics::wkb_path_rate(slope, Ik, Iik, km2);
                    const A4 x{DD(slope), DD(Ik), DD(Iik), DD(km2)};
                    const auto q = fp(x);
                    const auto S = verify::reference::condition_sums<5>(fp, x);
                    const std::array<double, 5> sim{s.rate, s.d_slope, s.d_I_k, s.d_I_ik, s.d_km2};
                    const std::string ctx = " (F " + shortest(F) + ", I_k " + shortest(Ik) + ", I_ik " +
                                            shortest(Iik) + ", z " + shortest(z) + ")";
                    for (std::size_t k = 0; k < 5; ++k) {
                        Predicate pr{std::string(path_range), 0.0, 0.0, true};
                        r.comparisons.push_back(value_record(pinf, "path " + pn[k] + ctx, pu[k], z, sim[k], q[k], S[k],
                                                             {pr}));
                    }
                    // Finding (xii)'s floor, from the shipped pre t (k_m^2 e / I_ik - (1 - e) / I_ik^2), the
                    // exponential's argument rounded (1 + z).
                    if (Iik > 0.0) {
                        const double pre = std::abs(slope) / (36.0 * base::hbar_J_s), t = std::exp(-2.0 * Ik);
                        const double e = std::exp(-z), br = -std::expm1(-z);
                        const double T = pre * t * (km2 * e / Iik * (1.0 + z) + br / (Iik * Iik));
                        classify_floor<5>(r, {0.0, 0.0, 0.0, T, 0.0}, xii_floor);
                    }
                }
            }
        }
    }
    apply_known_findings(r, known_findings());
    note_if_stale(r, xi_floor);
    note_if_stale(r, xii_floor);
    return r;
}

}  // namespace NiTCAD::verify::benchmark
