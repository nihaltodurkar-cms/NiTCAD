// TVB-31.11-31.16, lifetimes, recombination and interface traps (Unit 55, TVR-3 step 5; DECISIONS 15.3, 15.9, 15.15).
//
// Code under test: physics::scharfetter_lifetime, srh_recombination, auger_recombination, radiative_recombination
// (recombination.hpp); fermi_occupancy, trap_kinetics and trap_band_levels (interface_traps.hpp). Reference:
// verify/reference/model_recombination.hpp in double-double, its partials checked at L0 against complex step and SRH
// against the Brent root of the capture-emission balance (model_recombination_test). Criterion: 15.9's value rule.
//   31.11 Scharfetter, both carriers: N 0 and 1e13-1e21 (and 1e22 beyond). tau0 L and N_ref PC: the formula is
//         judged on the shipped numbers and the table rows are findings.
//   31.12 SRH, R and both partials: n 1e-10-1e21 on the log grid, p = r n_ie^2 / n for r in {1e-6, 1, 1e6} (r = 1
//         the zero), with the Boltzmann E (n_ie^2) and the Fermi-Dirac E of statistics.hpp (its value and partials
//         as inputs); tau_n, tau_p the shipped Scharfetter values at 1e16; n_ie the shipped n_i at 300 K.
//   31.13 Auger, the same grid, the register's Cn and Cp (PB).
//   31.14 Radiative, the same grid with the Boltzmann E (and the Fermi-Dirac E for GaAs), B of each shipped material
//         with a non-zero B (the rows are
//         read in step 9; here the shipped numbers).
//   31.15 fermi_occupancy over x in [-750, 750] (step 1/2, with 0 and its neighbours); trap_kinetics over n on the
//         log grid, p = r n1 p1 / n, levels tau = (E_t - E_i) / kT in {-20, 0, 20} (and dn1/dn = -0.2 n1 / n,
//         dp1/dp = -0.2 p1 / p at tau 0, r 1), cn = cp / 3 = 1e-8 cm^3/s.
//   31.16 trap_band_levels: the occupied and empty densities, sum over NiTCAD's levels of D_i f(x_i), against the
//         adaptive Gauss-Kronrod integral (TC3, MI-Q8) over four bands up to the gap: E_F -0.5..0.5 eV in 1/256 at
//         300 K, and T 200-500 K at E_F 0 and 0.25 eV.
#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/physics/interface_traps.hpp"
#include "NiTCAD/physics/recombination.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/statistics.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/model_recombination.hpp"
#include "tvb31_common.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

namespace mr = verify::reference::model_recombination;
using tvb31::DD;
using tvb31::in_range;
using tvb31::value_check;
using tvb31::value_record;
using tvb31::zero_or_in_range;
using tvb31::FloorFinding;
using tvb31::classify_floor;
using tvb31::note_if_stale;

constexpr std::string_view n_range = "N 0 or within the source's range [1e13, 1e21] cm^-3";
constexpr std::string_view c_range = "n and p within [1e-10, 1e21] cm^-3";
constexpr std::string_view x_range = "x within [-745, 745] (e^|x| finite)";
constexpr std::string_view ef_range = "E_F within [-0.5, 0.5] eV";
constexpr std::string_view t_range = "T within [200, 500] K";
constexpr std::string_view verification = "model_recombination_test: partials against complex step of cancellation-"
                                          "free value forms (TC8), SRH against the Brent root of the capture-emission "
                                          "balance (TC5), the band closed form against Gauss-Kronrod (TC3); 15.9's "
                                          "rule";

constexpr std::array<double, 3> ratios{1e-6, 1.0, 1e6};

std::string ratio_name(double r) { return "np/E " + shortest(r); }

// n on the log grid of 31.12-31.15, with the range's ends.
std::vector<double> density_grid(int per) { return tvb31::log_grid(1e-10, 1e21, per); }

template <std::size_t N, std::size_t K, class F>
void push_all(BenchmarkResult& r, const ReferenceInfo& inf, const std::string& prefix,
              const std::array<std::string, K>& names, const std::array<std::string, K>& units, double parameter,
              const std::array<double, K>& sim, F&& f, const std::array<DD, N>& x, const std::vector<Predicate>& pr) {
    const std::array<DD, K> q = f(x);
    const std::array<double, K> S = verify::reference::condition_sums<K>(f, x);
    for (std::size_t k = 0; k < K; ++k) {
        const std::string name = prefix.empty() ? names[k] : names[k] + " (" + prefix + ")";
        r.comparisons.push_back(value_record(inf, name, units[k], parameter, sim[k], q[k], S[k], pr));
    }
}

// Findings (viii) and (ix) (DECISIONS Unit 55 15.15, MI-Q15): srh_recombination and trap_kinetics form their partials
// by the bare quotient rule (A D - B C) / D^2, whose two terms cancel when one capture term dominates D or f is near 0
// or 1. The reference expands them so that the cancellation is symbolic. A partial that misses the value rule is the
// finding only where its error is within that form's own rounding floor, 64 eps T with T = (|A| D + |B C|) / D^2 the
// size of its terms; a larger error (a wrong term is of the order of T) stays a fail.
constexpr std::string_view finding_viii =
    "DECISIONS Unit 55 15.15 MI-Q15 (known finding): srh_recombination's partials by the bare quotient rule cancel "
    "when one capture term dominates the denominator; error within the form's rounding floor; fix proposed "
    "separately (MI-Q10)";
constexpr std::string_view finding_ix =
    "DECISIONS Unit 55 15.15 MI-Q15 (known finding): trap_kinetics's partials by the bare quotient rule cancel when "
    "f or 1 - f is near 1 or one capture term dominates; error within the form's rounding floor; fix proposed "
    "separately (MI-Q10)";

}  // namespace

BenchmarkResult run_tvb31_11(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-31.11";
    r.title = "Scharfetter doping-dependent lifetimes";
    r.info = tvb31::entry_info(r.id, Provenance::PB,
                               "the Scharfetter relation tau_min + (tau_max - tau_min) / (1 + (N/N_ref)^gamma) as "
                               "Gritsch's thesis sec. 4.7 prints it (Fossum and Lee 1982), read 2026-11-06; NiTCAD's "
                               "form is its tau_min = 0, gamma = 1 case",
                               "tau = tau0 / (1 + N / N_ref)", {std::string(n_range)}, {7},
                               "verify/reference/model_recombination.hpp", std::string(verification));
    r.info.assumptions.push_back("tau_n0, tau_p0: the shipped values (grade L); N_ref: the register (PC)");
    const ReferenceInfo inf = r.info;
    auto made = physics::Semiconductor::create(
        tvb31::register_parameters(physics::legacy_silicon_parameters, "legacy_silicon", {"31.11"}));
    if (!made) {
        r.run_error = made.error().message;
        return r;
    }
    const physics::Semiconductor& m = *made;
    std::vector<double> grid{0.0};
    for (const double N : tvb31::log_grid(1e13, 1e21, tvb31::per_decade(tier))) grid.push_back(N);
    grid.push_back(1e22);
    const auto& lt = m.parameters().lifetime;
    for (const bool electron : {true, false}) {
        const auto carrier = electron ? physics::Carrier::electron : physics::Carrier::hole;
        const DD tau0(electron ? lt.tau_n0 : lt.tau_p0);
        for (const double N : grid) {
            using A3 = std::array<DD, 3>;
            r.comparisons.push_back(value_check(
                inf, std::string(electron ? "tau_n" : "tau_p"), "s", N, physics::scharfetter_lifetime(m, carrier, N),
                [](const A3& y) { return mr::scharfetter_lifetime(y[0], y[1], y[2]); },
                A3{tau0, DD(lt.N_ref), DD(N)}, {zero_or_in_range(std::string(n_range), N, 1e13, 1e21)}));
        }
    }
    tvb31::table_checks(r, "31.11");
    apply_known_findings(r, known_findings());
    return r;
}

namespace {

// The common grid of 31.12-31.14: n, p = r n_ie^2 / n, and the equilibrium product (Boltzmann or Fermi-Dirac).
struct CarrierPoint {
    double n, p, ratio;
    physics::EquilibriumProduct E;
    bool fd;
};

std::vector<CarrierPoint> carrier_points(int per, double n_ie, double Nc, double Nv, bool with_fd) {
    std::vector<CarrierPoint> pts;
    const double gn = std::log(Nc / n_ie), gp = std::log(Nv / n_ie);
    for (const bool fd : {false, true}) {
        if (fd && !with_fd) continue;
        for (const double ratio : ratios) {
            for (const double n : density_grid(per)) {
                const double p = ratio * n_ie * (n_ie / n);
                physics::EquilibriumProduct E = physics::boltzmann_equilibrium_product(n_ie);
                if (fd) {
                    E = physics::fermi_dirac_equilibrium_product(n_ie, physics::fermi_dirac_degeneracy(n_ie, gn, n),
                                                                 physics::fermi_dirac_degeneracy(n_ie, gp, p));
                }
                pts.push_back({n, p, ratio, E, fd});
            }
        }
    }
    return pts;
}

std::vector<Predicate> carrier_predicates(const CarrierPoint& c) {
    Predicate p = in_range(std::string(c_range), c.n, 1e-10, 1e21);
    p.holds = p.holds && c.p >= 1e-10 && c.p <= 1e21;
    return {p};
}

std::string carrier_prefix(const CarrierPoint& c) {
    return std::string(c.fd ? "Fermi-Dirac E" : "Boltzmann E") + ", " + ratio_name(c.ratio);
}

struct SiliconAt300 {
    double n_ie, Nc, Nv, tau_n, tau_p;
};

SiliconAt300 silicon_at_300() {
    const physics::Semiconductor& m = physics::legacy_silicon();
    return {physics::intrinsic_density(m, 300.0), m.parameters().Nc300, m.parameters().Nv300,
            physics::scharfetter_lifetime(m, physics::Carrier::electron, 1e16),
            physics::scharfetter_lifetime(m, physics::Carrier::hole, 1e16)};
}

const std::array<std::string, 3> rate_names{"R", "dR/dn", "dR/dp"};

}  // namespace

BenchmarkResult run_tvb31_12(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-31.12";
    r.title = "Shockley-Read-Hall recombination and its partials";
    r.info = tvb31::entry_info(r.id, Provenance::PD,
                               "derived here from the capture-emission kinetics of one mid-gap level (Shockley and "
                               "Read 1952, Hall 1952; the standard form, the papers not read)",
                               "R = (n p - E) / (tau_p (n + n_ie) + tau_n (p + n_ie)), the chain through E",
                               {std::string(c_range)}, {7}, "verify/reference/model_recombination.hpp",
                               std::string(verification));
    const ReferenceInfo inf = r.info;
    const SiliconAt300 si = silicon_at_300();
    const std::array<std::string, 3> units{"cm^-3 s^-1", "s^-1", "s^-1"};
    FloorFinding viii{"(viii)", finding_viii};
    for (const CarrierPoint& c : carrier_points(tvb31::per_decade(tier), si.n_ie, si.Nc, si.Nv, true)) {
        const auto R = physics::srh_recombination(c.n, c.p, c.E, si.n_ie, si.tau_n, si.tau_p);
        push_all<8, 3>(r, inf, carrier_prefix(c), rate_names, units, c.n, {R.rate, R.d_dn, R.d_dp},
                       [](const std::array<DD, 8>& y) { return mr::srh(y); },
                       std::array<DD, 8>{DD(c.n), DD(c.p), DD(c.E.value), DD(c.E.d_dn), DD(c.E.d_dp), DD(si.n_ie),
                                         DD(si.tau_n), DD(si.tau_p)},
                       carrier_predicates(c));
        // The shipped form's term sizes: ((p - E_n) D - (n p - E) tau_p) / D^2 and its mirror.
        const double D = si.tau_p * (c.n + si.n_ie) + si.tau_n * (c.p + si.n_ie);
        const double ex = std::abs(c.n * c.p - c.E.value);
        classify_floor<3>(r,
                          {0.0, (std::abs(c.p - c.E.d_dn) * D + ex * si.tau_p) / (D * D),
                           (std::abs(c.n - c.E.d_dp) * D + ex * si.tau_n) / (D * D)},
                          viii);
    }
    note_if_stale(r, viii);
    apply_known_findings(r, known_findings());
    return r;
}

BenchmarkResult run_tvb31_13(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-31.13";
    r.title = "Auger recombination and its partials";
    r.info = tvb31::entry_info(r.id, Provenance::PD,
                               "derived here from detailed balance (a rate proportional to the process's carriers, "
                               "less its reverse at the equilibrium product); Cn, Cp as Palankovski Table 3.36 prints "
                               "them (Dziewior and Schmid 1977, not read)",
                               "R = (Cn n + Cp p)(n p - E), the chain through E", {std::string(c_range)}, {7},
                               "verify/reference/model_recombination.hpp", std::string(verification));
    const ReferenceInfo inf = r.info;
    const SiliconAt300 si = silicon_at_300();
    const double Cn = tvb31::register_dd("legacy_silicon.auger.Cn").to_double();
    const double Cp = tvb31::register_dd("legacy_silicon.auger.Cp").to_double();
    const std::array<std::string, 3> units{"cm^-3 s^-1", "s^-1", "s^-1"};
    for (const CarrierPoint& c : carrier_points(tvb31::per_decade(tier), si.n_ie, si.Nc, si.Nv, true)) {
        const auto R = physics::auger_recombination(c.n, c.p, c.E, Cn, Cp);
        push_all<7, 3>(r, inf, carrier_prefix(c), rate_names, units, c.n, {R.rate, R.d_dn, R.d_dp},
                       [](const std::array<DD, 7>& y) { return mr::auger(y); },
                       std::array<DD, 7>{DD(c.n), DD(c.p), DD(c.E.value), DD(c.E.d_dn), DD(c.E.d_dp), DD(Cn), DD(Cp)},
                       carrier_predicates(c));
    }
    tvb31::table_checks(r, "31.13");
    apply_known_findings(r, known_findings());
    return r;
}

BenchmarkResult run_tvb31_14(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-31.14";
    r.title = "Radiative recombination and its partials";
    r.info = tvb31::entry_info(r.id, Provenance::PD,
                               "derived here from detailed balance; B per material as shipped (the coefficients are "
                               "read in step 9)",
                               "R = B (n p - E)", {std::string(c_range)}, {7},
                               "verify/reference/model_recombination.hpp",
                               std::string(verification));
    r.info.assumptions.push_back("B: the shipped value of each material with a non-zero B");
    const ReferenceInfo inf = r.info;
    const SiliconAt300 si = silicon_at_300();
    const std::array<std::string, 3> units{"cm^-3 s^-1", "s^-1", "s^-1"};
    struct Material {
        const char* name;
        double B;
        bool with_fd;  // the Fermi-Dirac E too, whose partials exercise the chain through E
    };
    const std::array<Material, 4> materials{{{"germanium", physics::germanium_parameters.radiative_cm3_s, false},
                                             {"gallium_arsenide", physics::gallium_arsenide_parameters.radiative_cm3_s, true},
                                             {"indium_gallium_arsenide",
                                              physics::indium_gallium_arsenide_parameters.radiative_cm3_s, false},
                                             {"silicon_carbide_4h", physics::silicon_carbide_4h_parameters.radiative_cm3_s, false}}};
    for (const Material& mat : materials) {
        for (const CarrierPoint& c : carrier_points(tvb31::per_decade(tier), si.n_ie, si.Nc, si.Nv, mat.with_fd)) {
            const auto R = physics::radiative_recombination(c.n, c.p, c.E, mat.B);
            push_all<6, 3>(r, inf, std::string(mat.name) + ", " + carrier_prefix(c), rate_names, units, c.n,
                           {R.rate, R.d_dn, R.d_dp}, [](const std::array<DD, 6>& y) { return mr::radiative(y); },
                           std::array<DD, 6>{DD(c.n), DD(c.p), DD(c.E.value), DD(c.E.d_dn), DD(c.E.d_dp), DD(mat.B)},
                           carrier_predicates(c));
        }
    }
    apply_known_findings(r, known_findings());
    return r;
}

BenchmarkResult run_tvb31_15(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-31.15";
    r.title = "Interface-trap occupancy and capture-emission rate";
    r.info = tvb31::entry_info(r.id, Provenance::PD,
                               "derived here from the capture-emission kinetics of one level (the Shockley-Read "
                               "form, Unit 15b)",
                               "f = (cn n + cp p1) / D, U / N_t = cn cp (n p - n1 p1) / D, D = cn (n + n1) + cp (p + "
                               "p1); f_eq = 1 / (1 + e^x)",
                               {std::string(c_range), std::string(x_range)}, {7},
                               "verify/reference/model_recombination.hpp", std::string(verification));
    // Each half judged under its own applicability.
    ReferenceInfo occ_inf = r.info;
    occ_inf.applicability = {std::string(x_range)};
    ReferenceInfo inf = r.info;
    inf.applicability = {std::string(c_range)};
    // The equilibrium occupancy.
    std::vector<double> xs{0.0, std::nextafter(0.0, 1.0), std::nextafter(0.0, -1.0)};
    for (double x = -750.0; x <= 750.0; x += 0.5) xs.push_back(x);
    const std::array<std::string, 3> occ_names{"f_eq", "1 - f_eq", "df_eq/deta"};
    const std::array<std::string, 3> occ_units{"1", "1", "1"};
    for (const double x : xs) {
        const physics::FermiOccupancy o = physics::fermi_occupancy(x);
        push_all<1, 3>(r, occ_inf, "", occ_names, occ_units, x, {o.occupied, o.empty, o.d_eta},
                       [](const std::array<DD, 1>& y) { return mr::fermi_occupancy(y[0]); }, std::array<DD, 1>{DD(x)},
                       {in_range(std::string(x_range), x, -745.0, 745.0)});
    }
    // The kinetics of one level.
    const double n_ie = 1e10, cn = 1e-8, cp = 3e-8;
    FloorFinding ix{"(ix)", finding_ix};
    const std::array<std::string, 7> names{"f", "1 - f", "df/dn", "df/dp", "U/N_t", "dU/dn", "dU/dp"};
    const std::array<std::string, 7> units{"1", "1", "cm^3", "cm^3", "s^-1", "cm^3 s^-1", "cm^3 s^-1"};
    struct Level {
        double tau, ratio, g;
    };
    std::vector<Level> levels;
    for (const double tau : {-20.0, 0.0, 20.0}) {
        for (const double ratio : ratios) levels.push_back({tau, ratio, 0.0});
    }
    levels.push_back({0.0, 1.0, -0.2});
    for (const Level& l : levels) {
        const double n1 = n_ie * std::exp(l.tau), p1 = n_ie * std::exp(-l.tau);
        const std::string prefix = "level tau " + shortest(l.tau) + ", " + ratio_name(l.ratio) +
                                   (l.g != 0.0 ? ", dn1/dn -0.2 n1/n" : "");
        for (const double n : density_grid(tvb31::per_decade(tier))) {
            const double p = l.ratio * n1 * (p1 / n);
            const double dn1 = l.g * n1 / n, dp1 = l.g * p1 / p;
            const physics::TrapKinetics k = physics::trap_kinetics(n, p, n1, dn1, p1, dp1, cn, cp);
            Predicate pr = in_range(std::string(c_range), n, 1e-10, 1e21);
            pr.holds = pr.holds && p >= 1e-10 && p <= 1e21;
            push_all<8, 7>(r, inf, prefix, names, units, n,
                           {k.occupied, k.empty, k.d_occupied_dn, k.d_occupied_dp, k.rate, k.d_rate_dn, k.d_rate_dp},
                           [](const std::array<DD, 8>& y) { return mr::trap_kinetics(y); },
                           std::array<DD, 8>{DD(n), DD(p), DD(n1), DD(dn1), DD(p1), DD(dp1), DD(cn), DD(cp)}, {pr});
            // The shipped form's term sizes (interface_traps.hpp).
            const double D = cn * (n + n1) + cp * (p + p1), occ = cn * n + cp * p1, D2 = D * D;
            const double dDn = std::abs(cn * (1.0 + dn1)), dDp = std::abs(cp * (1.0 + dp1));
            const double ex = std::abs(n * p - n1 * p1), kk = cn * cp;
            classify_floor<7>(r,
                              {0.0, 0.0, (cn * D + occ * dDn) / D2, (std::abs(cp * dp1) * D + occ * dDp) / D2, 0.0,
                               kk * (std::abs(p - dn1 * p1) * D + ex * dDn) / D2,
                               kk * (std::abs(n - n1 * dp1) * D + ex * dDp) / D2},
                              ix);
        }
    }
    note_if_stale(r, ix);
    apply_known_findings(r, known_findings());
    return r;
}

BenchmarkResult run_tvb31_16(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-31.16";
    r.title = "Trap-band quadrature (trap_band_levels)";
    r.info = tvb31::entry_info(r.id, Provenance::PD,
                               "the occupied density of a uniform band, integral of D f((E - E_F) / kT) dE, by "
                               "adaptive Gauss-Kronrod in double-double (TC3, MI-Q8); in-model numerics",
                               "sum_i D_i f((E_i - E_F) / kT) over trap_band_levels against the integral",
                               {std::string(ef_range), std::string(t_range)}, {7},
                               "verify/reference/model_recombination.hpp", std::string(verification));
    r.info.assumptions.push_back("kT: the reference's CODATA k_B / q times T; NiTCAD's base::thermal_voltage");
    const ReferenceInfo inf = r.info;
    const double D = 1e11;
    const std::array<std::array<double, 2>, 4> bands{{{-0.54, 0.54}, {-0.2, 0.3}, {0.1, 0.101}, {-0.3, -0.05}}};
    struct Point {
        double EF, T;
    };
    std::vector<Point> pts;
    for (int i = -128; i <= 128; ++i) pts.push_back({i / 256.0, 300.0});
    pts.push_back({-0.6, 300.0});
    pts.push_back({0.6, 300.0});
    for (const double EF : {0.0, 0.25}) {
        for (const double T : tvb31::log_grid(200.0, 500.0, tvb31::per_decade(tier))) pts.push_back({EF, T});
        pts.push_back({EF, 150.0});
        pts.push_back({EF, 600.0});
    }
    for (const auto& b : bands) {
        physics::TrapBand band;
        band.density_cm2_eV = D;
        band.energy_low_eV = b[0];
        band.energy_high_eV = b[1];
        const std::string bname = "band [" + shortest(b[0]) + ", " + shortest(b[1]) + "] eV";
        for (const Point& p : pts) {
            const std::vector<physics::TrapLevel> levels = physics::trap_band_levels(band, p.T);
            const double kT = base::thermal_voltage(p.T);
            double occupied = 0.0, empty = 0.0;
            for (const physics::TrapLevel& l : levels) {
                const physics::FermiOccupancy o = physics::fermi_occupancy((l.energy_eV - p.EF) / kT);
                occupied += l.density_cm2 * o.occupied;
                empty += l.density_cm2 * o.empty;
            }
            const std::array<DD, 5> x{DD(D), DD(b[0]), DD(b[1]), DD(p.EF), DD(p.T)};
            const std::vector<Predicate> pr{in_range(std::string(ef_range), p.EF, -0.5, 0.5),
                                            in_range(std::string(t_range), p.T, 200.0, 500.0)};
            const double param = p.T == 300.0 ? p.EF : p.T;
            const std::string at = p.T == 300.0 ? " (300 K)" : " (E_F " + shortest(p.EF) + ")";
            for (const int k : {0, 1}) {
                const mr::BandIntegral q = mr::band_density_quadrature(x, k);
                if (!q.ok) {
                    r.run_error = "Gauss-Kronrod did not converge on " + bname + at;
                    return r;
                }
                auto closed = [k](const std::array<DD, 5>& y) { return mr::band_density_closed(y, k); };
                r.comparisons.push_back(value_record(inf, bname + (k == 0 ? " occupied" : " empty") + at, "cm^-2",
                                                     param, k == 0 ? occupied : empty, q.value,
                                                     verify::reference::condition_sum(closed, x), pr));
            }
        }
    }
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
