// TVB-31.17-31.19, incomplete ionization and the neutral equilibria (Unit 55, TVR-3 step 6; DECISIONS 15.3, 15.9,
// 15.16).
//
// Code under test: physics::ionized_density, dopant_level_eV and bound_fraction, ionized_neutral_equilibrium
// (ionization.hpp); boltzmann_neutral_equilibrium and fermi_dirac_neutral_equilibrium (statistics.hpp). Reference:
// verify/reference/model_statistics.hpp in double-double (model_statistics_test). Criterion: 15.9's value rule.
//   31.17 ionized_density, N = 1e17: eta -40..40 in 1/4 at depths 0, the next double above 0, 0.045, 0.1 and 0.3 eV;
//         depth 0..0.3 eV in 1/256 of the range at eta -40, -5, 0, 5 and 40; g 2 and 4; b 1, 0.5 and 1e-3; kT at
//         300 K. The value and d/d eta. The shipped levels (PC) are table findings.
//   31.18 dopant_level_eV and bound_fraction, phosphorus and boron (the register's Table I, PA): N 0 and 1e14-1e21.
//   31.19 The neutral equilibria (in-model numerics, MI-Q8), legacy silicon at 300 K: Boltzmann and Fermi-Dirac with
//         complete ionization over net doping 0 and +-1-1e21; the ionized equilibrium (both statistics, the legacy
//         fixed levels and the Altermatt levels) over N_D or N_A 1e10-1e21 and compensated pairs about N_A = 1e17.
//         n, p and eta against the Brent root.
#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/physics/ionization.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/statistics.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/model_statistics.hpp"
#include "tvb31_common.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

namespace ms = verify::reference::model_statistics;
using tvb31::DD;
using tvb31::in_range;
using tvb31::register_dd;
using tvb31::value_check;
using tvb31::value_record;
using tvb31::zero_or_in_range;
using tvb31::FloorFinding;
using tvb31::classify_floor;
using tvb31::note_if_stale;

constexpr std::string_view eta_range = "eta within [-40, 40]";
constexpr std::string_view depth_range = "depth within [0, 0.3] eV";
constexpr std::string_view n_range = "N 0 or within the source's range [1e14, 1e21] cm^-3";
constexpr std::string_view c_range = "N_D and N_A within [0, 1e21] cm^-3";
// Finding (x) (DECISIONS Unit 55 15.16, MI-Q16): fermi_dirac_neutral_equilibrium solves ln major = ln c + log1p(minor
// / c), c = |C| / n_ie, whose two right-hand terms (about +-|ln c|) cancel when c << 1; its header documents the root
// as "good to a few 1e-14 in y". Where |eta| is small that absolute error misses 15.9's rule; it is the finding only
// within the residual's rounding floor.
constexpr std::string_view finding_x =
    "DECISIONS Unit 55 15.16 MI-Q16 (known finding): fermi_dirac_neutral_equilibrium's eta has the residual's "
    "absolute floor (its documented few 1e-14) where |eta| is small; fix proposed separately (MI-Q10)";
constexpr std::string_view verification = "model_statistics_test: d/d eta against complex step (TC8); the Boltzmann "
                                          "root against its closed form, the Fermi-Dirac root's balance, the implicit "
                                          "condition sums against re-solved roots; 15.9's rule";

}  // namespace

BenchmarkResult run_tvb31_17(Tier) {
    BenchmarkResult r;
    r.id = "TVB-31.17";
    r.title = "Incomplete ionization";
    r.info = tvb31::entry_info(r.id, Provenance::PD,
                               "derived here from the grand-canonical occupancy of one level in equilibrium with its "
                               "band (the standard form, as Genius reproduces it); the bound fraction b of "
                               "Schenk, Altermatt, Schmithuesen, SISPAD 2006, eq (3)",
                               "N_ion = N ((1 - b) + b / (1 + g e^(eta + E/kT))); complete ionization at E <= 0",
                               {std::string(eta_range), std::string(depth_range)}, {4},
                               "verify/reference/model_statistics.hpp", std::string(verification));
    const ReferenceInfo inf = r.info;
    const double kT = base::thermal_voltage(300.0);
    const double N = 1e17;
    struct Point {
        double eta, depth_eV;
        bool sweep_depth;
    };
    std::vector<Point> pts;
    for (const double depth : {0.0, -1.0, 0.045, 0.1, 0.3}) {  // -1: the next double above 0 in kT units
        for (int i = -160; i <= 160; ++i) pts.push_back({i / 4.0, depth, false});
        pts.push_back({-50.0, depth, false});
        pts.push_back({50.0, depth, false});
    }
    for (const double eta : {-40.0, -5.0, 0.0, 5.0, 40.0}) {
        for (int i = 0; i <= 256; ++i) pts.push_back({eta, 0.3 * i / 256.0, true});
        pts.push_back({eta, 0.35, true});
    }
    const std::array<std::string, 2> units{"cm^-3", "cm^-3"};
    for (const double g : {2.0, 4.0}) {
        for (const double b : {1.0, 0.5, 1e-3}) {
            for (const Point& p : pts) {
                const double depth_kT = p.depth_eV < 0.0 ? std::nextafter(0.0, 1.0) : p.depth_eV / kT;
                const double depth_eV = p.depth_eV < 0.0 ? 0.0 : p.depth_eV;
                const std::vector<Predicate> pr{in_range(std::string(eta_range), p.eta, -40.0, 40.0),
                                                in_range(std::string(depth_range), depth_eV, 0.0, 0.3)};
                const physics::IonizedDensity sim = physics::ionized_density(N, p.eta, depth_kT, g, b);
                const std::string ctx = " (g " + shortest(g) + ", b " + shortest(b) +
                                        (p.sweep_depth ? ", eta " + shortest(p.eta)
                                                       : ", depth " + (p.depth_eV < 0.0 ? std::string("0+")
                                                                                        : shortest(p.depth_eV)) +
                                                             " eV") +
                                        ")";
                using A5 = std::array<DD, 5>;
                const A5 x{DD(N), DD(p.eta), DD(depth_kT), DD(g), DD(b)};
                auto f = [](const A5& y) { return ms::ionized(y); };
                const std::array<DD, 2> q = f(x);
                const std::array<double, 2> S = verify::reference::condition_sums<2>(f, x);
                const double param = p.sweep_depth ? depth_eV : p.eta;
                r.comparisons.push_back(value_record(inf, "N_ion" + ctx, units[0], param, sim.value, q[0], S[0], pr));
                r.comparisons.push_back(value_record(inf, "dN_ion/deta" + ctx, units[1], param, sim.d_eta, q[1], S[1], pr));
            }
        }
    }
    tvb31::table_checks(r, "31.17");
    apply_known_findings(r, known_findings());
    return r;
}

BenchmarkResult run_tvb31_18(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-31.18";
    r.title = "Altermatt-Schenk density-dependent dopant level";
    r.info = tvb31::entry_info(r.id, Provenance::PA,
                               "Schenk, Altermatt, Schmithuesen, SISPAD 2006, pp. 51-54, eqs (1), (5), Table I (read "
                               "in Unit 37)",
                               "E_dop(N) = E_dop0 / (1 + (N/N_ref)^c), b(N) = 1 / (1 + (N/N_b)^d)",
                               {std::string(n_range)}, {4}, "verify/reference/model_statistics.hpp",
                               std::string(verification));
    const ReferenceInfo inf = r.info;
    std::vector<double> grid{0.0};
    for (const double N : tvb31::log_grid(1e14, 1e21, tvb31::per_decade(tier))) grid.push_back(N);
    grid.push_back(1e22);
    for (const char* dopant : {"phosphorus", "boron"}) {
        const std::string b = std::string("altermatt_") + dopant + ".";
        const physics::DensityDependentLevel l{register_dd(b + "E_dop0").to_double(), register_dd(b + "N_ref").to_double(),
                                               register_dd(b + "c").to_double(), register_dd(b + "N_b").to_double(),
                                               register_dd(b + "d").to_double()};
        for (const double N : grid) {
            const std::vector<Predicate> pr{zero_or_in_range(std::string(n_range), N, 1e14, 1e21)};
            using A4 = std::array<DD, 4>;
            using A3 = std::array<DD, 3>;
            r.comparisons.push_back(value_check(
                inf, std::string(dopant) + " E_dop", "eV", N, physics::dopant_level_eV(l, N),
                [](const A4& y) { return ms::altermatt_level(y); }, A4{DD(l.E0_eV), DD(l.N_ref), DD(l.c), DD(N)}, pr));
            r.comparisons.push_back(value_check(inf, std::string(dopant) + " b", "1", N, physics::bound_fraction(l, N),
                                                [](const A3& y) { return ms::altermatt_bound(y); },
                                                A3{DD(l.N_b), DD(l.d), DD(N)}, pr));
        }
    }
    tvb31::table_checks(r, "31.18");
    apply_known_findings(r, known_findings());
    return r;
}

BenchmarkResult run_tvb31_19(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-31.19";
    r.title = "Charge-neutral equilibria (Boltzmann, Fermi-Dirac, ionized)";
    r.info = tvb31::entry_info(r.id, Provenance::PD,
                               "the root in eta of n - p - N_D+ + N_A- = 0 (statistics.hpp, ionization.hpp), by Brent "
                               "(TC5) in double-double; F_1/2 the TVB-02 reference; in-model numerics (MI-Q8)",
                               "n(eta) - p(eta) = N_D+(eta) - N_A-(eta), n = n_ie e^eta or n_ie e^g F_1/2(eta - g)",
                               {std::string(c_range)}, {4}, "verify/reference/model_statistics.hpp",
                               std::string(verification));
    const ReferenceInfo inf = r.info;
    const physics::Semiconductor& m = physics::legacy_silicon();
    const double kT = base::thermal_voltage(300.0);
    const double n_ie = physics::intrinsic_density(m, 300.0);
    const double gn = std::log(m.parameters().Nc300 / n_ie), gp = std::log(m.parameters().Nv300 / n_ie);
    const physics::IonizationParameters& ion = m.parameters().ionization;
    const int per = tvb31::per_decade(tier);
    const std::array<std::string, 3> names{"n", "p", "eta"};
    const std::array<std::string, 3> units{"cm^-3", "cm^-3", "1"};

    FloorFinding x_floor{"(x)", finding_x};
    enum class Model { complete, fixed, altermatt };
    auto run = [&](const std::string& label, double ND, double NA, Model model, bool fd, double parameter) {
        physics::DopantLevels lv{0.0, 0.0, ion.donor_degeneracy, ion.acceptor_degeneracy};
        if (model == Model::fixed) {
            lv.donor_kT = ion.donor_eV / kT;
            lv.acceptor_kT = ion.acceptor_eV / kT;
        } else if (model == Model::altermatt) {
            lv.donor_kT = physics::dopant_level_eV(physics::altermatt_phosphorus, ND) / kT;
            lv.acceptor_kT = physics::dopant_level_eV(physics::altermatt_boron, NA) / kT;
            lv.donor_bound = physics::bound_fraction(physics::altermatt_phosphorus, ND);
            lv.acceptor_bound = physics::bound_fraction(physics::altermatt_boron, NA);
        }
        physics::NeutralEquilibrium sim{};
        if (model == Model::complete) {
            sim = fd ? physics::fermi_dirac_neutral_equilibrium(ND - NA, n_ie, gn, gp)
                     : physics::boltzmann_neutral_equilibrium(ND - NA, n_ie);
        } else {
            sim = physics::ionized_neutral_equilibrium(ND, NA, n_ie, gn, gp, lv, fd);
        }
        const ms::EquilibriumInputs x{DD(ND), DD(NA), DD(n_ie), DD(gn), DD(gp), DD(lv.donor_kT), DD(lv.acceptor_kT),
                                      DD(lv.donor_degeneracy), DD(lv.acceptor_degeneracy), DD(lv.donor_bound),
                                      DD(lv.acceptor_bound)};
        // The bracket starts about NiTCAD's eta (when finite); its ends are verified by sign (model_statistics.hpp).
        const double guess = std::isfinite(sim.eta) ? sim.eta : std::asinh((ND - NA) / (2.0 * n_ie));
        const ms::Equilibrium e = ms::equilibrium(x, fd, guess);
        if (!e.ok) {
            r.run_error = "reference root failed: " + label + " at " + shortest(parameter);
            return false;
        }
        const std::vector<Predicate> pr{[&] {
            Predicate p = in_range(std::string(c_range), ND, 0.0, 1e21);
            p.holds = p.holds && NA <= 1e21;
            return p;
        }()};
        const std::array<double, 3> s{sim.n, sim.p, sim.eta};
        const std::array<DD, 3> q{e.n, e.p, e.eta};
        for (int k = 0; k < 3; ++k) {
            r.comparisons.push_back(
                value_record(inf, names[k] + " (" + label + ")", units[k], parameter, s[k], q[k], e.S[k], pr));
        }
        if (fd && model == Model::complete) {
            // (x): the residual's ln c + log1p(minor / c), c = |C| / n_ie, terms of about 2 |ln c| each.
            const double c = std::abs(ND - NA) / n_ie;
            const double T = std::abs(e.eta.hi) + (c > 0.0 ? 2.0 * std::abs(std::log(c)) : 0.0);
            classify_floor<3>(r, {0.0, 0.0, T}, x_floor);
        }
        return true;
    };

    std::vector<double> C{0.0};
    for (const double c : tvb31::log_grid(1.0, 1e21, per)) {
        C.push_back(c);
        C.push_back(-c);
    }
    C.push_back(1e22);
    C.push_back(-1e22);
    for (const bool fd : {false, true}) {
        const std::string stat = fd ? "Fermi-Dirac" : "Boltzmann";
        for (const double c : C) {
            if (!run(stat + ", complete", c > 0.0 ? c : 0.0, c < 0.0 ? -c : 0.0, Model::complete, fd, c)) return r;
        }
        for (const Model model : {Model::fixed, Model::altermatt}) {
            const std::string lab = stat + (model == Model::fixed ? ", fixed levels" : ", Altermatt levels");
            for (const double N : tvb31::log_grid(1e10, 1e21, per)) {
                if (!run(lab + ", n-type", N, 0.0, model, fd, N)) return r;
                if (!run(lab + ", p-type", 0.0, N, model, fd, -N)) return r;
            }
            // Compensated: N_A = 1e17, N_D = N_A (1 + s 10^-k).
            const double NA = 1e17;
            if (!run(lab + ", compensated", NA, NA, model, fd, 0.0)) return r;
            for (int k = 1; k <= 12; ++k) {
                for (const double s : {1.0, -1.0}) {
                    const double ND = NA * (1.0 + s * std::pow(10.0, -k));
                    if (!run(lab + ", compensated", ND, NA, model, fd, ND - NA)) return r;
                }
            }
        }
    }
    note_if_stale(r, x_floor);
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
