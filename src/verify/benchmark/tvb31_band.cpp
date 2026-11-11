// TVB-31.01-31.05, the band set and band-gap narrowing (Unit 55, TVR-3 step 3; DECISIONS 15.3, 15.9).
//
// Code under test: physics::band_gap_eV, conduction_band_dos, valence_band_dos, intrinsic_density,
// intrinsic_level_depth_eV (semiconductor.hpp), band_gap_slope_eV_per_K, dos_mass_log, dos_log_slope,
// intrinsic_density_log_slope (thermal.hpp), bandgap_narrowing_eV, effective_intrinsic_density
// (bandgap_narrowing.hpp), on the shipped sets with the register's printed values written over the fields they bind
// (the formula check), and the shipped fields against the register (the table check; 15.2).
// Reference: verify/reference/model_band.hpp in double-double, its derivatives checked against complex step at L0
// (model_band_test). Criterion: 15.9's value rule for values and partials, frozen before code.
//   31.01 legacy Varshni gap and power-law Nc, Nv, 1-500 K;
//   31.02 Passler gap (Couderc eq (10), Table I), 1-400 K;
//   31.03 Couderc's density-of-states mass fits, Nc, Nv and the log slopes, 1-400 K;
//   31.04 n_i, the intrinsic level's depth and d ln n_i / dT on both sets (legacy 1-500 K, default 1-400 K: its
//         fits' range); n_i judged where it is a normal double (below about 9 K it is not);
//   31.05 Slotboom narrowing, n_ie and d ln n_ie / dT at 300 K, N 1e14-1e21 cm^-3 with N0 and its neighbours: the
//         documented form (bandgap_narrowing.hpp, an implementation check) and the published form (F1, MI-Q4: its
//         miss is the known finding (F1)).
// 31.04's dEg/dT and d ln N/dT are judged in 31.01-31.03, where they belong to their forms.
// Couderc's printed n_i (eq (15)) is a separate fit, not a value of these formulas, so no printed value is judged.
#include <array>
#include <expected>
#include <cfloat>
#include <cmath>
#include <string>
#include <vector>

#include "NiTCAD/physics/bandgap_narrowing.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/thermal.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/model_band.hpp"
#include "tvb31_common.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

namespace mb = verify::reference::model_band;
using tvb31::DD;
using tvb31::in_range;
using tvb31::register_dd;
using tvb31::value_check;
using A4 = std::array<DD, 4>;

constexpr std::string_view legacy_range = "T within the source's range [1, 500] K";
constexpr std::string_view couderc_range = "T within the source's range [1, 400] K";
constexpr std::string_view ni_normal = "n_i a normal double";
constexpr std::string_view doping_range = "N within the source's range [1e14, 1e21] cm^-3";

ReferenceInfo base_info(std::string id, Provenance g, std::string citation, std::string form,
                        std::vector<std::string> applicability) {
    ReferenceInfo r;
    r.id = std::move(id);
    r.cls = ReferenceClass::M;
    r.independence = Independence::I1;
    r.provenance = g;
    r.citation = std::move(citation);
    r.version = "register rows of verify/reference/model_parameters.hpp (DECISIONS 15.10)";
    r.equations = {5};
    r.form = std::move(form);
    r.assumptions = {"parameters: the register's printed decimals, exact in double-double"};
    r.applicability = std::move(applicability);
    r.derivation = "verify/reference/model_band.hpp";
    r.verification = "model_band_test: each analytic derivative against complex step (TC8) to 64 eps S; identities";
    return r;
}

std::expected<physics::Semiconductor, std::string> make(const physics::SemiconductorParameters& p) {
    auto m = physics::Semiconductor::create(p);
    if (!m) return std::unexpected(m.error().message);
    return *m;
}

struct LegacyInputs {
    DD Eg0 = register_dd("legacy_silicon.Eg0_eV");
    DD alpha = register_dd("legacy_silicon.varshni_alpha_eV_per_K");
    DD beta = register_dd("legacy_silicon.varshni_beta_K");
    DD Nc300 = register_dd("legacy_silicon.Nc300");
    DD Nv300 = register_dd("legacy_silicon.Nv300");
};

struct CoudercInputs {
    DD Eg0 = register_dd("silicon.Eg0_eV");
    DD alpha = register_dd("silicon.passler_gap.alpha_eV_per_K");
    DD Theta = register_dd("silicon.passler_gap.Theta_K");
    DD Delta = register_dd("silicon.passler_gap.Delta");
    A4 c{register_dd("silicon.dos_mass.conduction[0]"), register_dd("silicon.dos_mass.conduction[1]"),
         register_dd("silicon.dos_mass.conduction[2]"), register_dd("silicon.dos_mass.conduction[3]")};
    A4 v{register_dd("silicon.dos_mass.valence[0]"), register_dd("silicon.dos_mass.valence[1]"),
         register_dd("silicon.dos_mass.valence[2]"), register_dd("silicon.dos_mass.valence[3]")};
};

physics::SemiconductorParameters legacy_register() {
    return tvb31::register_parameters(physics::legacy_silicon_parameters, "legacy_silicon", {"31.01", "31.05"});
}
physics::SemiconductorParameters default_register() {
    return tvb31::register_parameters(physics::silicon_parameters, "silicon", {"31.02", "31.03"});
}

// The legacy band set as a function of x = {Eg0, alpha, beta, Nc300, Nv300, T, ...}.
DD legacy_ln_ni(const DD* x) {
    return mb::intrinsic_log_density(mb::power_dos(x[3], x[5]), mb::power_dos(x[4], x[5]),
                                     mb::varshni_gap(x[0], x[1], x[2], x[5]), x[5]);
}
DD legacy_ln_ni_slope(const DD* x, const DD& bgn) {
    const DD dlnN = mb::power_dos_log_slope(x[5]);
    return mb::intrinsic_log_slope(dlnN, dlnN, mb::varshni_gap(x[0], x[1], x[2], x[5]),
                                   mb::varshni_gap_slope(x[1], x[2], x[5]), bgn, x[5]);
}

// The default set as a function of x = {Eg0, alpha, Theta, Delta, c0..c3, v0..v3, T}.
A4 slice(const DD* x) { return {x[0], x[1], x[2], x[3]}; }
DD couderc_Eg(const DD* x) { return mb::passler_gap(x[0], x[1], x[2], x[3], x[12]); }
DD couderc_ln_ni(const DD* x) {
    return mb::intrinsic_log_density(mb::couderc_dos(slice(x + 4), x[12]), mb::couderc_dos(slice(x + 8), x[12]),
                                     couderc_Eg(x), x[12]);
}

}  // namespace

BenchmarkResult run_tvb31_01() {
    BenchmarkResult r;
    r.id = "TVB-31.01";
    r.title = "Legacy silicon band set: Varshni gap and power-law densities of states";
    r.info = base_info(r.id, Provenance::PB,
                       "Varshni, Physica 34, 149 (1967), form as reproduced by Palankovski (TU Wien 2000) eq 3.60 and "
                       "Xyce bandgapTbl.tex; N(T) = N300 (T/300)^(3/2), Palankovski eqs 3.91-3.92 (read 2026-11-06)",
                       "Eg = Eg0 - alpha T^2 / (T + beta), its T slope; Nc, Nv = N300 (T / 300)^(3/2), d ln N / dT",
                       {std::string(legacy_range)});
    const ReferenceInfo inf = r.info;
    const auto m = make(legacy_register());
    if (!m) {
        r.run_error = m.error();
        return r;
    }
    const LegacyInputs p;
    for (const double T : tvb31::temperature_grid(1.0, 500.0, {550.0, 600.0})) {
        const std::vector<Predicate> pr{in_range(std::string(legacy_range), T, 1.0, 500.0)};
        const DD t(T);
        r.comparisons.push_back(value_check(
            inf, "Eg", "eV", T, physics::band_gap_eV(*m, T),
            [](const A4& x) { return mb::varshni_gap(x[0], x[1], x[2], x[3]); }, A4{p.Eg0, p.alpha, p.beta, t}, pr));
        r.comparisons.push_back(value_check(
            inf, "dEg/dT", "eV/K", T, physics::band_gap_slope_eV_per_K(*m, T),
            [](const std::array<DD, 3>& x) { return mb::varshni_gap_slope(x[0], x[1], x[2]); },
            std::array<DD, 3>{p.alpha, p.beta, t}, pr));
        for (const bool electron : {true, false}) {
            const auto band = electron ? physics::Carrier::electron : physics::Carrier::hole;
            const double sim = electron ? physics::conduction_band_dos(*m, T) : physics::valence_band_dos(*m, T);
            r.comparisons.push_back(value_check(
                inf, electron ? "Nc" : "Nv", "cm^-3", T, sim,
                [](const std::array<DD, 2>& x) { return mb::power_dos(x[0], x[1]); },
                std::array<DD, 2>{electron ? p.Nc300 : p.Nv300, t}, pr));
            r.comparisons.push_back(value_check(
                inf, electron ? "d ln Nc/dT" : "d ln Nv/dT", "1/K", T, physics::dos_log_slope(*m, band, T),
                [](const std::array<DD, 1>& x) { return mb::power_dos_log_slope(x[0]); }, std::array<DD, 1>{t}, pr));
        }
    }
    tvb31::table_checks(r, "31.01");
    apply_known_findings(r, known_findings());
    return r;
}

BenchmarkResult run_tvb31_02() {
    BenchmarkResult r;
    r.id = "TVB-31.02";
    r.title = "Passler band gap (Couderc et al. 2014, eq (10), Table I)";
    r.info = base_info(r.id, Provenance::PA,
                       "R. Couderc, M. Amara, M. Lemiti, J. Appl. Phys. 115, 093705 (2014), eq (10) and Table I, open "
                       "PDF read in Unit 35",
                       "Eg = Eg0 - alpha Theta [(1 - 3 Delta^2) / (e^(Theta/T) - 1) + (3 Delta^2 / 2)(P^(1/6) - 1)], "
                       "chi = 2T / Theta; its T slope",
                       {std::string(couderc_range)});
    const ReferenceInfo inf = r.info;
    const auto m = make(default_register());
    if (!m) {
        r.run_error = m.error();
        return r;
    }
    const CoudercInputs p;
    for (const double T : tvb31::temperature_grid(1.0, 400.0, {450.0, 500.0})) {
        const std::vector<Predicate> pr{in_range(std::string(couderc_range), T, 1.0, 400.0)};
        const DD t(T);
        r.comparisons.push_back(value_check(
            inf, "Eg", "eV", T, physics::band_gap_eV(*m, T),
            [](const std::array<DD, 5>& x) { return mb::passler_gap(x[0], x[1], x[2], x[3], x[4]); },
            std::array<DD, 5>{p.Eg0, p.alpha, p.Theta, p.Delta, t}, pr));
        r.comparisons.push_back(value_check(
            inf, "dEg/dT", "eV/K", T, physics::band_gap_slope_eV_per_K(*m, T),
            [](const A4& x) { return mb::passler_gap_slope(x[0], x[1], x[2], x[3]); },
            A4{p.alpha, p.Theta, p.Delta, t}, pr));
    }
    tvb31::table_checks(r, "31.02");
    apply_known_findings(r, known_findings());
    return r;
}

BenchmarkResult run_tvb31_03() {
    BenchmarkResult r;
    r.id = "TVB-31.03";
    r.title = "Couderc density-of-states mass fits and Nc(T), Nv(T)";
    r.info = base_info(r.id, Provenance::PA,
                       "Couderc et al. (2014), eqs (4), (5), (18), (19), open PDF read in Unit 35",
                       "N = 4.83e15 M(T) T^(3/2), M = (m_d / m0)^(3/2) = c0 T^3 + c1 T^2 + c2 T + c3; ln M and its slope; "
                       "d ln N / dT = 3 / (2T) + M' / M",
                       {std::string(couderc_range)});
    const ReferenceInfo inf = r.info;
    const auto m = make(default_register());
    if (!m) {
        r.run_error = m.error();
        return r;
    }
    const CoudercInputs p;
    using A5 = std::array<DD, 5>;
    auto with_t = [](const A4& c, const DD& t) { return A5{c[0], c[1], c[2], c[3], t}; };
    for (const double T : tvb31::temperature_grid(1.0, 400.0, {450.0, 500.0})) {
        const std::vector<Predicate> pr{in_range(std::string(couderc_range), T, 1.0, 400.0)};
        const DD t(T);
        for (const bool electron : {true, false}) {
            const auto band = electron ? physics::Carrier::electron : physics::Carrier::hole;
            const A5 x = with_t(electron ? p.c : p.v, t);
            const std::string b = electron ? "conduction" : "valence";
            const double sim = electron ? physics::conduction_band_dos(*m, T) : physics::valence_band_dos(*m, T);
            r.comparisons.push_back(value_check(
                inf, electron ? "Nc" : "Nv", "cm^-3", T, sim,
                [](const A5& y) { return mb::couderc_dos(slice(y.data()), y[4]); }, x, pr));
            const physics::ValueSlope lm = physics::dos_mass_log(*m, band, T);
            r.comparisons.push_back(value_check(
                inf, "ln M " + b, "1", T, lm.value,
                [](const A5& y) { return log(mb::mass_polynomial(slice(y.data()), y[4])); }, x, pr));
            r.comparisons.push_back(value_check(
                inf, "d ln M/dT " + b, "1/K", T, lm.d_dT,
                [](const A5& y) {
                    return mb::mass_polynomial_slope(slice(y.data()), y[4]) / mb::mass_polynomial(slice(y.data()), y[4]);
                },
                x, pr));
            r.comparisons.push_back(value_check(
                inf, electron ? "d ln Nc/dT" : "d ln Nv/dT", "1/K", T, physics::dos_log_slope(*m, band, T),
                [](const A5& y) { return mb::couderc_dos_log_slope(slice(y.data()), y[4]); }, x, pr));
        }
    }
    tvb31::table_checks(r, "31.03");
    apply_known_findings(r, known_findings());
    return r;
}

BenchmarkResult run_tvb31_04() {
    BenchmarkResult r;
    r.id = "TVB-31.04";
    r.title = "Intrinsic density, intrinsic level and d ln n_i / dT on both silicon sets";
    r.info = base_info(r.id, Provenance::PD,
                       "derived here: n_i^2 = Nc Nv e^(-Eg / kT) (Couderc eq (1)) and n = p = n_i under Boltzmann "
                       "statistics; k = CODATA 2018 k / q",
                       "n_i = exp(ln(Nc Nv) / 2 - Eg / 2kT); depth = chi + Eg / 2 + (kT / 2) ln(Nc / Nv); "
                       "d ln n_i / dT = (d ln Nc + d ln Nv) / 2dT + (Eg - T dEg/dT) / (2 k T^2)",
                       {"T within the source's range"});
    r.info.cls = ReferenceClass::E;
    r.info.second_method = "complex step (TC8) of the reference's ln n_i, model_band_test";
    r.info.assumptions.push_back("chi: the shipped value (grade L in the register; the depth's T dependence is judged)");
    const auto legacy = make(legacy_register());
    const auto deflt = make(default_register());
    if (!legacy || !deflt) {
        r.run_error = !legacy ? legacy.error() : deflt.error();
        return r;
    }
    ReferenceInfo inf_ni = r.info;
    inf_ni.applicability.push_back(std::string(ni_normal));
    const ReferenceInfo inf = r.info;
    const double ln_min = std::log(DBL_MIN);
    const LegacyInputs l;
    const CoudercInputs c;
    const DD chi(physics::legacy_silicon_parameters.electron_affinity_eV);
    for (const bool is_legacy : {true, false}) {
        const physics::Semiconductor& m = is_legacy ? *legacy : *deflt;
        const double hi = is_legacy ? 500.0 : 400.0;
        const std::string set = is_legacy ? "legacy: " : "default: ";
        for (const double T : tvb31::temperature_grid(1.0, hi, {hi + 50.0, hi + 100.0})) {
            const Predicate range = in_range("T within the source's range", T, 1.0, hi);
            const DD t(T);
            if (is_legacy) {
                using A6 = std::array<DD, 6>;
                using A7 = std::array<DD, 7>;
                const A6 x{l.Eg0, l.alpha, l.beta, l.Nc300, l.Nv300, t};
                const double ln_ni = legacy_ln_ni(x.data()).to_double();
                const std::vector<Predicate> pn{range, Predicate{std::string(ni_normal), ln_ni, ln_min, ln_ni >= ln_min}};
                r.comparisons.push_back(value_check(
                    inf_ni, set + "n_i", "cm^-3", T, physics::intrinsic_density(m, T),
                    [](const A6& y) { return exp(legacy_ln_ni(y.data())); }, x, pn));
                r.comparisons.push_back(value_check(
                    inf, set + "depth chi + Ec - Ei", "eV", T, physics::intrinsic_level_depth_eV(m, T),
                    [](const A7& y) {
                        return y[6] + mb::conduction_to_intrinsic(mb::power_dos(y[3], y[5]), mb::power_dos(y[4], y[5]),
                                                                  mb::varshni_gap(y[0], y[1], y[2], y[5]), y[5]);
                    },
                    A7{l.Eg0, l.alpha, l.beta, l.Nc300, l.Nv300, t, chi}, {range}));
                r.comparisons.push_back(value_check(
                    inf, set + "d ln n_i/dT", "1/K", T, physics::intrinsic_density_log_slope(m, T),
                    [](const A6& y) { return legacy_ln_ni_slope(y.data(), DD(0.0)); }, x, {range}));
            } else {
                using A13 = std::array<DD, 13>;
                using A14 = std::array<DD, 14>;
                const A13 x{c.Eg0, c.alpha, c.Theta, c.Delta, c.c[0], c.c[1], c.c[2],
                            c.c[3], c.v[0],  c.v[1],  c.v[2],  c.v[3], t};
                const double ln_ni = couderc_ln_ni(x.data()).to_double();
                const std::vector<Predicate> pn{range, Predicate{std::string(ni_normal), ln_ni, ln_min, ln_ni >= ln_min}};
                r.comparisons.push_back(value_check(
                    inf_ni, set + "n_i", "cm^-3", T, physics::intrinsic_density(m, T),
                    [](const A13& y) { return exp(couderc_ln_ni(y.data())); }, x, pn));
                A14 xd{};
                for (std::size_t i = 0; i < 13; ++i) xd[i] = x[i];
                xd[13] = chi;
                r.comparisons.push_back(value_check(
                    inf, set + "depth chi + Ec - Ei", "eV", T, physics::intrinsic_level_depth_eV(m, T),
                    [](const A14& y) {
                        return y[13] + mb::conduction_to_intrinsic(mb::couderc_dos(slice(y.data() + 4), y[12]),
                                                                   mb::couderc_dos(slice(y.data() + 8), y[12]),
                                                                   couderc_Eg(y.data()), y[12]);
                    },
                    xd, {range}));
                r.comparisons.push_back(value_check(
                    inf, set + "d ln n_i/dT", "1/K", T, physics::intrinsic_density_log_slope(m, T),
                    [](const A13& y) {
                        return mb::intrinsic_log_slope(mb::couderc_dos_log_slope(slice(y.data() + 4), y[12]),
                                                       mb::couderc_dos_log_slope(slice(y.data() + 8), y[12]),
                                                       couderc_Eg(y.data()),
                                                       mb::passler_gap_slope(y[1], y[2], y[3], y[12]), DD(0.0), y[12]);
                    },
                    x, {range}));
            }
        }
    }
    apply_known_findings(r, known_findings());
    return r;
}

BenchmarkResult run_tvb31_05() {
    BenchmarkResult r;
    r.id = "TVB-31.05";
    r.title = "Slotboom band-gap narrowing and n_ie (documented and published forms)";
    r.info = base_info(r.id, Provenance::PD,
                       "documented form: bandgap_narrowing.hpp, the published form (Slotboom, de Graaff; Tidy3D and "
                       "Genius, PB) less E0 sqrt(1/2) above N0 and 0 at or below it (F1, DECISIONS 15.4)",
                       "dEg = E0 [ln(N/N0) + sqrt(ln^2(N/N0) + 1/2) - sqrt(1/2)] for N > N0, else 0; "
                       "n_ie = n_i e^(dEg / 2kT); d ln n_ie / dT = d ln n_i / dT - dEg / (2 k T^2); legacy set, 300 K",
                       {std::string(doping_range)});
    ReferenceInfo pub = r.info;
    pub.provenance = Provenance::PB;
    pub.citation = "Slotboom and de Graaff's form as reproduced by Tidy3D (Slotboom model) and Genius (Si BGN), "
                   "values Klaassen, Slotboom, de Graaff 1992 (read 2026-11-06)";
    pub.form = "dEg = E0 [ln(N/N0) + sqrt(ln^2(N/N0) + 1/2)] for every N";
    const ReferenceInfo doc = r.info;
    const auto m = make(legacy_register());
    if (!m) {
        r.run_error = m.error();
        return r;
    }
    const LegacyInputs l;
    const DD E0 = register_dd("legacy_silicon.bandgap_narrowing.E0_eV");
    const DD N0 = register_dd("legacy_silicon.bandgap_narrowing.N0");
    const double n0 = N0.to_double();
    std::vector<double> grid;
    for (int i = 14 * 64; i <= 21 * 64; ++i) grid.push_back(std::pow(10.0, i / 64.0));
    for (const double n : {n0, std::nextafter(n0, 0.0), std::nextafter(n0, INFINITY), 1e13, 3e21, 1e22}) {
        grid.push_back(n);
    }
    // Points approaching the cut from above (a cut moved off N0 by a relative 10^-j is caught there, not only by the
    // nextafter neighbours, where the narrowing is below the bound).
    for (int j = 1; j <= 12; ++j) grid.push_back(n0 * (1.0 + std::pow(10.0, -j)));
    constexpr double T = 300.0;
    const DD t(T);
    using A3 = std::array<DD, 3>;
    using A9 = std::array<DD, 9>;  // Eg0, alpha, beta, Nc300, Nv300, T, E0, N0, N
    for (const double N : grid) {
        const std::vector<Predicate> pr{in_range(std::string(doping_range), N, 1e14, 1e21)};
        const double sim_bgn = physics::bandgap_narrowing_eV(*m, N);
        const A9 x{l.Eg0, l.alpha, l.beta, l.Nc300, l.Nv300, t, E0, N0, DD(N)};
        for (const bool published : {false, true}) {
            const ReferenceInfo& inf = published ? pub : doc;
            const std::string pre = published ? "published form: " : "";
            auto bgn = published ? +[](const DD& e, const DD& n0_, const DD& n) { return mb::slotboom_published(e, n0_, n); }
                                 : +[](const DD& e, const DD& n0_, const DD& n) { return mb::slotboom_documented(e, n0_, n); };
            r.comparisons.push_back(value_check(
                inf, pre + "dEg", "eV", N, sim_bgn, [bgn](const A3& y) { return bgn(y[0], y[1], y[2]); },
                A3{E0, N0, DD(N)}, pr));
            r.comparisons.push_back(value_check(
                inf, pre + "n_ie (300 K)", "cm^-3", N, physics::effective_intrinsic_density(*m, N, T),
                [bgn](const A9& y) { return mb::effective_intrinsic(legacy_ln_ni(y.data()), bgn(y[6], y[7], y[8]), y[5]); },
                x, pr));
            r.comparisons.push_back(value_check(
                inf, pre + "d ln n_ie/dT (300 K)", "1/K", N, physics::intrinsic_density_log_slope(*m, T, sim_bgn),
                [bgn](const A9& y) { return legacy_ln_ni_slope(y.data(), bgn(y[6], y[7], y[8])); }, x, pr));
        }
    }
    tvb31::table_checks(r, "31.05");
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
