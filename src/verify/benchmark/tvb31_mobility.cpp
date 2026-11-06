// TVB-31.06-31.10, mobility (Unit 55, TVR-3 step 4; DECISIONS 15.3, 15.9, 15.14).
//
// Code under test: physics::caughey_thomas_mobility (mobility.hpp) and its T slope (thermal.hpp), masetti_mobility
// and masetti_mobility_slope, canali_mobility (field_mobility.hpp), lombardi_mobility (surface_mobility.hpp),
// carrier_temperature_mobility (energy_balance.hpp). Reference: verify/reference/model_mobility.hpp in double-double,
// its derivatives checked at L0 (model_mobility_test). Criterion: 15.9's value rule for values and partials.
//   31.06 Caughey-Thomas, both carriers: N 0 and 1e13-1e21 at 300 K; T 200-500 K at N 0, 1e15, 1e17, 1e19. The
//         values are grade L (F3): the formula is judged on the shipped numbers, the table rows are findings.
//   31.07 Masetti, boron and phosphorus: N 0 and 1e13-1e21 at 300 K; T 200-500 K at N 0 and 1e17 (the documented
//         mu_max(T), tag C). The published form has no T dependence: its slope is 0, and the miss is (F6).
//   31.08 Canali, both carriers, with the shipped (beta 2, 1) and the law's 300 K parameters: E 0 and 1e-2-1e6 V/cm,
//         mu0 the shipped Caughey-Thomas mu at N = 0, 300 K. The table rows are judged as law_300K rows (F2).
//   31.09 Lombardi, both carriers: E_perp 1-1e7 V/cm at 300 K for N 0 and 1e18; T 200-500 K at 1e5 V/cm, 1e18;
//         mu_b the shipped Caughey-Thomas mu there; the value and its three partials.
//   31.10 The energy-balance mobility, beta 1, 1.5 and 2: a = 0, +-a0's neighbours, a0 / 2, negative a down to -1e-3
//         and 1e-12-1e6 on the log grid; mu0 1360, tau_w 0.35 ps, v_sat 1.07e7; the six returned quantities.
#include <array>
#include <cmath>
#include <expected>
#include <stdexcept>
#include <string>
#include <vector>

#include "NiTCAD/physics/energy_balance.hpp"
#include "NiTCAD/physics/field_mobility.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/surface_mobility.hpp"
#include "NiTCAD/physics/thermal.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/model_mobility.hpp"
#include "tvb31_common.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

namespace mm = verify::reference::model_mobility;
using tvb31::DD;
using tvb31::in_range;
using tvb31::register_dd;
using tvb31::value_check;
using tvb31::zero_or_in_range;

constexpr std::string_view n_range = "N 0 or within the source's range [1e13, 1e21] cm^-3";
constexpr std::string_view t_range = "T within the source's range [200, 500] K";
constexpr std::string_view e_range = "E 0 or within [1e-2, 1e6] V/cm";
constexpr std::string_view ep_range = "E_perp within [1, 1e7] V/cm";
constexpr std::string_view a_range = "a within [-1e-3, 1e6]";
constexpr std::string_view verification = "model_mobility_test: derivatives against complex step (TC8), the balance "
                                          "root against Brent (TC5); 15.9's rule";

std::string carrier_name(bool electron) { return electron ? "electron" : "hole"; }

std::vector<double> doping_grid(int per) {
    std::vector<double> g{0.0};
    for (const double n : tvb31::log_grid(1e13, 1e21, per)) g.push_back(n);
    g.push_back(1e22);  // beyond the range
    return g;
}

}  // namespace

BenchmarkResult run_tvb31_06(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-31.06";
    r.title = "Caughey-Thomas doping and temperature mobility";
    r.info = tvb31::entry_info(r.id, Provenance::PB,
                               "Caughey and Thomas, Proc. IEEE 55, 2192 (1967); the form as reproduced by Xyce and "
                               "Genius (analytic model), read 2026-11-06; mu_max (T/300)^gamma as mobility.hpp "
                               "documents it",
                               "mu = mu_min + (mu_max (T/300)^gamma - mu_min) / (1 + (N/N_ref)^alpha); d mu / dT",
                               {std::string(n_range), std::string(t_range)}, {6}, "verify/reference/model_mobility.hpp",
                               std::string(verification));
    r.info.assumptions.push_back("mu_min, mu_max, N_ref, alpha: the shipped values (grade L, F3); gamma: the register");
    const ReferenceInfo inf = r.info;
    auto made = physics::Semiconductor::create(
        tvb31::register_parameters(physics::legacy_silicon_parameters, "legacy_silicon", {"31.06"}));
    if (!made) {
        r.run_error = made.error().message;
        return r;
    }
    const physics::Semiconductor& m = *made;
    struct Point {
        double N, T;
    };
    std::vector<Point> pts;
    for (const double N : doping_grid(tvb31::per_decade(tier))) pts.push_back({N, 300.0});
    for (const double N : {0.0, 1e15, 1e17, 1e19}) {
        for (const double T : tvb31::log_grid(200.0, 500.0, tvb31::per_decade(tier))) pts.push_back({N, T});
    }
    pts.push_back({1e17, 150.0});
    pts.push_back({1e17, 600.0});
    for (const bool electron : {true, false}) {
        const physics::CaugheyThomasParameters& c =
            electron ? m.parameters().electron_mobility : m.parameters().hole_mobility;
        const auto carrier = electron ? physics::Carrier::electron : physics::Carrier::hole;
        const std::string name = carrier_name(electron);
        const DD gamma = register_dd(std::string("legacy_silicon.") + (electron ? "electron" : "hole") +
                                     "_mobility.T_exponent");
        for (const Point& p : pts) {
            const std::vector<Predicate> pr{zero_or_in_range(std::string(n_range), p.N, 1e13, 1e21),
                                            in_range(std::string(t_range), p.T, 200.0, 500.0)};
            using A7 = std::array<DD, 7>;
            const A7 x{DD(c.mu_min), DD(c.mu_max), DD(c.N_ref), DD(c.alpha), gamma, DD(p.N), DD(p.T)};
            const double param = p.T == 300.0 ? p.N : p.T;
            const std::string at = p.T == 300.0 ? " (300 K)" : " (N " + shortest(p.N) + ")";
            r.comparisons.push_back(value_check(
                inf, name + " mu" + at, "cm^2/(V s)", param, physics::caughey_thomas_mobility(m, carrier, p.N, p.T),
                [](const A7& y) { return mm::ct_mobility(y[0], y[1], y[2], y[3], y[4], y[5], y[6]); }, x, pr));
            r.comparisons.push_back(value_check(
                inf, name + " dmu/dT" + at, "cm^2/(V s K)", param,
                physics::caughey_thomas_mobility_slope(m, carrier, p.N, p.T),
                [](const A7& y) { return mm::ct_mobility_slope(y[1], y[2], y[3], y[4], y[5], y[6]); }, x, pr));
        }
    }
    tvb31::table_checks(r, "31.06");
    apply_known_findings(r, known_findings());
    return r;
}

BenchmarkResult run_tvb31_07(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-31.07";
    r.title = "Masetti doping mobility (boron and phosphorus)";
    r.info = tvb31::entry_info(r.id, Provenance::PB,
                               "Masetti, Severi, Solmi, IEEE TED 30, 764 (1983), eqs (1), (4), Table I, as Allpix "
                               "Squared 3.3 and Palankovski Table 3.24 reproduce them (read 2026-11-06)",
                               "mu = mu_min1 e^(-P_c/N) + (mu_max - mu_min2) / (1 + (N/C_r)^alpha) - mu_1 / (1 + "
                               "(C_s/N)^beta), at 300 K",
                               {std::string(n_range), std::string(t_range)}, {6}, "verify/reference/model_mobility.hpp",
                               std::string(verification));
    ReferenceInfo doc = r.info;
    doc.provenance = Provenance::PD;
    doc.citation = "documented form (mobility.hpp): the published form with mu_max (T/300)^gamma, gamma a NiTCAD "
                   "choice (tag C, F6)";
    doc.form = "the published form with mu_max (T/300)^gamma; d mu / dT";
    const ReferenceInfo pub = r.info;
    const int per = tvb31::per_decade(tier);
    for (const bool boron : {true, false}) {
        const std::string set = boron ? "masetti_boron_silicon." : "masetti_phosphorus_silicon.";
        const physics::MasettiParameters& shipped =
            boron ? physics::masetti_boron_silicon : physics::masetti_phosphorus_silicon;
        // The formula check runs on the register's printed values (P's P_c: 0, eq (1) has no such term).
        const char* names[10] = {"mu_min1", "mu_min2", "mu_max", "mu_1", "P_c", "C_r", "C_s", "alpha", "beta",
                                 "T_exponent"};
        std::array<DD, 10> p{};
        for (int i = 0; i < 10; ++i) {
            p[i] = (!boron && i == 4) ? DD(0.0) : register_dd(set + names[i]);
        }
        physics::MasettiParameters q{};
        double* fields[10] = {&q.mu_min1, &q.mu_min2, &q.mu_max, &q.mu_1, &q.P_c,
                              &q.C_r,     &q.C_s,     &q.alpha,  &q.beta, &q.T_exponent};
        for (int i = 0; i < 10; ++i) *fields[i] = p[i].to_double();
        (void)shipped;
        const std::string name = boron ? "boron" : "phosphorus";
        using A12 = std::array<DD, 12>;
        auto pack = [&p](double N, double T) {
            A12 x{};
            for (int i = 0; i < 10; ++i) x[i] = p[i];
            x[10] = DD(N);
            x[11] = DD(T);
            return x;
        };
        auto params = [](const A12& y) {
            std::array<DD, 10> a{};
            for (int i = 0; i < 10; ++i) a[i] = y[i];
            return a;
        };
        auto mu_f = [params](const A12& y) { return mm::masetti_mobility(params(y), y[10], y[11]); };
        auto slope_f = [params](const A12& y) { return mm::masetti_mobility_slope(params(y), y[10], y[11]); };
        for (const double N : doping_grid(per)) {
            const std::vector<Predicate> pr{zero_or_in_range(std::string(n_range), N, 1e13, 1e21),
                                            in_range(std::string(t_range), 300.0, 200.0, 500.0)};
            const A12 x = pack(N, 300.0);
            r.comparisons.push_back(value_check(pub, name + " mu (300 K)", "cm^2/(V s)", N,
                                                physics::masetti_mobility(q, N, 300.0), mu_f, x, pr));
            const double slope = physics::masetti_mobility_slope(q, N, 300.0);
            r.comparisons.push_back(value_check(doc, name + " dmu/dT (300 K)", "cm^2/(V s K)", N, slope, slope_f, x, pr));
            r.comparisons.push_back(value_check(
                pub, "published form: " + name + " dmu/dT (300 K)", "cm^2/(V s K)", N, slope,
                [](const A12&) { return DD(0.0); }, x, pr));
        }
        for (const double N : {0.0, 1e17}) {
            std::vector<double> ts = tvb31::log_grid(200.0, 500.0, per);
            ts.push_back(150.0);
            ts.push_back(600.0);
            for (const double T : ts) {
                const std::vector<Predicate> pr{zero_or_in_range(std::string(n_range), N, 1e13, 1e21),
                                                in_range(std::string(t_range), T, 200.0, 500.0)};
                const A12 x = pack(N, T);
                const std::string at = " (N " + shortest(N) + ")";
                r.comparisons.push_back(value_check(doc, name + " mu" + at, "cm^2/(V s)", T,
                                                    physics::masetti_mobility(q, N, T), mu_f, x, pr));
                r.comparisons.push_back(value_check(doc, name + " dmu/dT" + at, "cm^2/(V s K)", T,
                                                    physics::masetti_mobility_slope(q, N, T), slope_f, x, pr));
            }
        }
    }
    tvb31::table_checks(r, "31.07");
    apply_known_findings(r, known_findings());
    return r;
}

BenchmarkResult run_tvb31_08(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-31.08";
    r.title = "Canali field-dependent mobility";
    r.info = tvb31::entry_info(r.id, Provenance::PB,
                               "Canali et al., IEEE TED 22, 1045 (1975), the form as reproduced by Allpix Squared 3.3 "
                               "and Ayalew's thesis (TU Wien) sec. 3.3.3, read 2026-11-06",
                               "mu = mu0 / (1 + (mu0 E / v_sat)^beta)^(1/beta); d mu / dE (right-hand at E = 0)",
                               {std::string(e_range)}, {6}, "verify/reference/model_mobility.hpp",
                               std::string(verification));
    const ReferenceInfo inf = r.info;
    std::vector<double> grid{0.0};
    for (const double e : tvb31::log_grid(1e-2, 1e6, tvb31::per_decade(tier))) grid.push_back(e);
    grid.push_back(1e7);
    const auto& legacy = physics::legacy_silicon_parameters;
    for (const bool electron : {true, false}) {
        const physics::CanaliParameters& shipped = electron ? legacy.electron_saturation : legacy.hole_saturation;
        const double mu0 = (electron ? legacy.electron_mobility : legacy.hole_mobility).mu_max;  // CT at N = 0, 300 K
        // The law's 300 K values (the table check's reference), as a second parameter set.
        const auto law = [](std::string_view binding) {
            for (const auto& l : reference::model_parameters::laws) {
                if (l.binding == binding) {
                    return reference::decimal::to_dd(l.coefficient) *
                           pow(DD(300.0), reference::decimal::to_dd(l.exponent));
                }
            }
            throw std::logic_error("no law");
        };
        const std::string b = std::string("legacy_silicon.") + (electron ? "electron" : "hole") + "_saturation.";
        const physics::CanaliParameters law_set{law(b + "v_sat_cm_s").to_double(), law(b + "beta").to_double()};
        for (const bool use_law : {false, true}) {
            const physics::CanaliParameters& c = use_law ? law_set : shipped;
            const std::string name = carrier_name(electron) + (use_law ? " (law at 300 K)" : " (shipped)");
            for (const double E : grid) {
                const std::vector<Predicate> pr{zero_or_in_range(std::string(e_range), E, 1e-2, 1e6)};
                using A4 = std::array<DD, 4>;
                const A4 x{DD(mu0), DD(c.v_sat_cm_s), DD(c.beta), DD(E)};
                const physics::FieldMobility f = physics::canali_mobility(mu0, E, c);
                r.comparisons.push_back(value_check(
                    inf, name + " mu", "cm^2/(V s)", E, f.mobility,
                    [](const A4& y) { return mm::canali_mobility(y[0], y[1], y[2], y[3]); }, x, pr));
                r.comparisons.push_back(value_check(
                    inf, name + " dmu/dE", "cm^3/(V^2 s)", E, f.d_dE,
                    [](const A4& y) { return mm::canali_mobility_slope(y[0], y[1], y[2], y[3]); }, x, pr));
            }
        }
    }
    tvb31::table_checks(r, "31.08");
    apply_known_findings(r, known_findings());
    return r;
}

BenchmarkResult run_tvb31_09(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-31.09";
    r.title = "Lombardi surface mobility and its partials";
    r.info = tvb31::entry_info(r.id, Provenance::PB,
                               "Lombardi, Manzini, Saporito, Vanzi, IEEE TCAD 7, 1164 (1988), as Xyce and Genius "
                               "reproduce it (Unit 46)",
                               "1/mu = 1/mu_b + 1/(B/E + C N^lambda / (T E^(1/3))) + E^2 / delta; d mu / dE, dT, "
                               "d mu_b",
                               {std::string(ep_range), std::string(t_range)}, {6}, "verify/reference/model_mobility.hpp",
                               std::string(verification));
    const ReferenceInfo inf = r.info;
    const auto& legacy = physics::Semiconductor::create(physics::legacy_silicon_parameters);
    if (!legacy) {
        r.run_error = legacy.error().message;
        return r;
    }
    const int per = tvb31::per_decade(tier);
    for (const bool electron : {true, false}) {
        const std::string set = electron ? "lombardi_silicon_electrons." : "lombardi_silicon_holes.";
        const std::array<DD, 4> p{register_dd(set + "B"), register_dd(set + "C"), register_dd(set + "lambda"),
                                  register_dd(set + "delta")};
        const physics::LombardiParameters q{p[0].to_double(), p[1].to_double(), p[2].to_double(), p[3].to_double()};
        const auto carrier = electron ? physics::Carrier::electron : physics::Carrier::hole;
        struct Point {
            double N, T, E;
        };
        std::vector<Point> pts;
        for (const double N : {0.0, 1e18}) {
            std::vector<double> es = tvb31::log_grid(1.0, 1e7, per);
            es.push_back(0.5);
            es.push_back(1e8);
            for (const double E : es) pts.push_back({N, 300.0, E});
        }
        for (const double T : tvb31::log_grid(200.0, 500.0, per)) pts.push_back({1e18, T, 1e5});
        for (const Point& pt : pts) {
            const double mu_b = physics::caughey_thomas_mobility(*legacy, carrier, pt.N, pt.T);
            const std::vector<Predicate> pr{in_range(std::string(ep_range), pt.E, 1.0, 1e7),
                                            in_range(std::string(t_range), pt.T, 200.0, 500.0)};
            using A8 = std::array<DD, 8>;
            const A8 x{p[0], p[1], p[2], p[3], DD(mu_b), DD(pt.N), DD(pt.T), DD(pt.E)};
            const physics::SurfaceMobility s = physics::lombardi_mobility(q, mu_b, pt.N, pt.T, pt.E);
            const bool sweep_T = pt.E == 1e5 && pt.N == 1e18 && pt.T != 300.0;
            const double param = sweep_T ? pt.T : pt.E;
            const std::string at = sweep_T ? " (T sweep)" : " (N " + shortest(pt.N) + ")";
            const std::string name = carrier_name(electron);
            auto q4 = [](const A8& y) { return std::array<DD, 4>{y[0], y[1], y[2], y[3]}; };
            r.comparisons.push_back(value_check(
                inf, name + " mu" + at, "cm^2/(V s)", param, s.mobility,
                [q4](const A8& y) { return mm::lombardi_mobility(q4(y), y[4], y[5], y[6], y[7]); }, x, pr));
            const double sims[3] = {s.d_dE, s.d_dT, s.d_dbulk};
            const char* names[3] = {" dmu/dE", " dmu/dT", " dmu/dmu_b"};
            const char* units[3] = {"cm^3/(V^2 s)", "cm^2/(V s K)", "1"};
            for (int k = 0; k < 3; ++k) {
                r.comparisons.push_back(value_check(
                    inf, name + names[k] + at, units[k], param, sims[k],
                    [q4, k](const A8& y) { return mm::lombardi_partials(q4(y), y[4], y[5], y[6], y[7])[k]; }, x, pr));
            }
        }
    }
    tvb31::table_checks(r, "31.09");
    apply_known_findings(r, known_findings());
    return r;
}

BenchmarkResult run_tvb31_10(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-31.10";
    r.title = "Energy-balance carrier-temperature mobility and effective field";
    r.info = tvb31::entry_info(r.id, Provenance::PD,
                               "derived here (Unit 52 ED4 and the owner's beta = 1 ruling, DECISIONS 2026-10-29): the "
                               "uniform-field balance q tau_w mu_C(E) E^2 = (3/2) k (T - T_L) with Canali's mu_C",
                               "x^2 (1 + x^beta)^(-1/beta) = a, a = mu0 (3/2)(k/q) rise / (tau_w v_sat^2); mu = mu0 m, "
                               "E_eff = x v_sat / mu0; regularized below a0 = 1e-8 (beta != 2)",
                               {std::string(a_range)}, {6, 10}, "verify/reference/model_mobility.hpp",
                               std::string(verification));
    r.info.cls = ReferenceClass::E;
    r.info.second_method = "Brent's zero (TC5) of the balance, model_mobility_test";
    const ReferenceInfo inf = r.info;
    const double mu0 = 1360.0, tau = register_dd("energy_relaxation.electron").to_double(), v_sat = 1.07e7;
    const physics::CanaliParameters base{v_sat, 1.0};
    std::vector<double> as{-1e-3, -1e-6, -1e-9, 0.0, 0.5 * mm::regularization_a0, mm::regularization_a0,
                           std::nextafter(mm::regularization_a0, 0.0), std::nextafter(mm::regularization_a0, 1.0),
                           -2e-3, 2e6, 1e-20, 1e-40, 1e-100};
    for (const double a : tvb31::log_grid(1e-12, 1e6, tvb31::per_decade(tier))) as.push_back(a);
    const DD per_K = mm::balance_a(DD(mu0), DD(1.0), DD(tau), DD(v_sat));
    for (const double beta : {1.0, 1.5, 2.0}) {
        physics::CanaliParameters c = base;
        c.beta = beta;
        const std::string name = "beta " + shortest(beta);
        for (const double a : as) {
            const double rise = (DD(a) / per_K).to_double();
            const std::vector<Predicate> pr{in_range(std::string(a_range), a, -1e-3, 1e6)};
            using A4 = std::array<DD, 4>;
            const A4 x{DD(mu0), DD(rise), DD(tau), DD(v_sat)};
            const physics::CarrierTemperatureMobility s = physics::carrier_temperature_mobility(mu0, rise, tau, c);
            const double sims[6] = {s.mobility, s.d_mu0, s.d_rise, s.field, s.field_d_mu0, s.field_d_rise};
            const char* names[6] = {" mu", " dmu/dmu0", " dmu/drise", " E_eff", " dE_eff/dmu0", " dE_eff/drise"};
            const char* units[6] = {"cm^2/(V s)", "1", "cm^2/(V s K)", "V/cm", "V s/cm^3", "V/(cm K)"};
            auto all = [beta](const A4& y) {
                const mm::BalanceResult b = mm::carrier_temperature_mobility(y, beta);
                return std::array<DD, 6>{b.mobility, b.d_mu0, b.d_rise, b.field, b.field_d_mu0, b.field_d_rise};
            };
            const std::array<DD, 6> q = all(x);
            const std::array<double, 6> S = reference::condition_sums<6>(all, x);
            for (int k = 0; k < 6; ++k) {
                r.comparisons.push_back(
                    tvb31::value_record(inf, name + names[k], units[k], a, sims[k], q[k], S[k], pr));
            }
        }
    }
    tvb31::table_checks(r, "31.10");
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
