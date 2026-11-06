// TVB-31.26-31.29, thermionic emission, MLDA, the lattice conductivity and the thermopower factor (Unit 55, TVR-3
// step 8; DECISIONS 15.3, 15.9, 15.18).
//
// Code under test: physics::emission_velocity_cm_s (all three overloads) and emission_velocity_slope
// (thermionic_emission.hpp, thermal.hpp); mlda_wavelength_cm, mlda_factor and mlda_mean_factor
// (quantum_confinement.hpp); thermal_conductivity, fermi_orders_zero_one and fermi_dirac_thermal_diffusion
// (thermal.hpp). Reference: verify/reference/model_thermal.hpp in double-double (model_thermal_test).
//   31.26 v(N, T) over N 1e15-1e21 at 200, 300 and 500 K; v and dv/dT of both carriers through N(T) on the legacy
//         (power-law) and the default (Couderc) sets, with the density-of-states form and with test Richardson
//         constants (120 A/(cm^2 K^2), no shipped set has them), T 200-500 K; the Schottky overload's v.
//   31.27 lambda(m, T) for both quantization masses; f(z) at z / lambda 0-12 in 1/64; the mean factor at a / lambda
//         1e-8-1e3 with the series switch 1/2 and its neighbours.
//   31.28 kappa and d kappa / dT on the register's printed values, T 200-600 K; the three register rows.
//   31.29 F_-1, F_0, F_1, h and dh/dx at x -750..750 in 1/4, -2..2 in 1/1024, and the switch points +-1 with their
//         neighbours. Where |q_ref| < 2^-969 the reference's low part is subnormal and U_ref is floored at 2^-1074
//         (value_record, MI-Q19): near and below the underflow that exceeds 1e-3 of the bound, not_evaluated (TB-Q4).
#include <array>
#include <cfloat>
#include <cmath>
#include <string>
#include <vector>

#include "NiTCAD/physics/quantum_confinement.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/thermal.hpp"
#include "NiTCAD/physics/thermionic_emission.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/model_band.hpp"
#include "NiTCAD/verify/reference/model_thermal.hpp"
#include "tvb31_common.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

namespace mh = verify::reference::model_thermal;
namespace mb = verify::reference::model_band;
using tvb31::DD;
using tvb31::in_range;
using tvb31::register_dd;
using tvb31::value_record;

constexpr std::string_view t_range = "T within [200, 500] K";
constexpr std::string_view couderc_range = "T within the source's range [1, 400] K";
constexpr std::string_view any_n = "any N > 0 and T > 0 (the formula)";
constexpr std::string_view any_mlda = "any mass and T > 0 (the formula)";
constexpr std::string_view z_range = "z / lambda within [0, 10]";
constexpr std::string_view a_range = "a / lambda within [1e-8, 1e3]";
constexpr std::string_view kappa_range = "T within [200, 600] K";
constexpr std::string_view x_range = "x within [-750, 750]";

constexpr std::string_view finding_xiii =
    "DECISIONS Unit 55 15.18 MI-Q21 (known finding, ruled before the first run): fermi_dirac_thermal_diffusion's "
    "dh/dx = 2 - 2 F_1 F_-1 / F_0^2 cancels as x -> -inf, where dh/dx -> e^x / 2 (an absolute error of eps); fix "
    "proposed separately (MI-Q10)";

Predicate holds(std::string_view name) { return Predicate{std::string(name), 0.0, 0.0, true}; }

}  // namespace

BenchmarkResult run_tvb31_26(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-31.26";
    r.title = "Thermionic emission velocity and its temperature derivative";
    r.info = tvb31::entry_info(r.id, Provenance::PD,
                               "derived here: v = sqrt(kT / (2 pi m)) with m from N = 2 (2 pi m k T / h^2)^(3/2); "
                               "Richardson's A T^2 / (q N)",
                               "v = (kT / h)(2 / N)^(1/3), or A T^2 / (q N); dv/dT = v (1/T - s/3), or v (2/T - s)",
                               {std::string(any_n)}, {11}, "verify/reference/model_thermal.hpp",
                               "model_thermal_test: the closed form against sqrt(kT / (2 pi m)), dv/dT against "
                               "complex step (TC8); 15.9's rule");
    r.info.cls = ReferenceClass::E;
    r.info.second_method = "sqrt(kT / (2 pi m)) with the implied mass; complex step on v(T) (L0)";
    ReferenceInfo dinf = r.info;
    ReferenceInfo tinf = r.info;
    tinf.applicability = {std::string(t_range)};
    ReferenceInfo cinf = r.info;
    cinf.applicability = {std::string(t_range), std::string(couderc_range)};
    // v(N, T).
    using A2 = std::array<DD, 2>;
    auto fv = [](const A2& y) { return mh::emission_dos(y); };
    for (const double T : {200.0, 300.0, 500.0}) {
        for (const double N : tvb31::log_grid(1e15, 1e21, tvb31::per_decade(tier))) {
            const A2 x{DD(N), DD(T)};
            r.comparisons.push_back(value_record(dinf, "v(N, T) (T " + shortest(T) + ")", "cm/s", N,
                                                 physics::emission_velocity_cm_s(N, T), fv(x),
                                                 verify::reference::condition_sum(fv, x), {holds(any_n)}));
        }
    }
    // Through N(T), both sets, both forms.
    auto legacy = tvb31::register_parameters(physics::legacy_silicon_parameters, "legacy_silicon", {"31.01"});
    auto deflt = tvb31::register_parameters(physics::silicon_parameters, "silicon", {"31.02", "31.03"});
    constexpr double A_test = 120.0;
    for (const bool couderc : {false, true}) {
        for (const bool rich : {false, true}) {
            physics::SemiconductorParameters p = couderc ? deflt : legacy;
            if (rich) p.richardson = {.electron = A_test, .hole = A_test};
            auto made = physics::Semiconductor::create(p);
            if (!made) {
                r.run_error = made.error().message;
                return r;
            }
            const physics::Semiconductor& m = *made;
            for (const bool electron : {true, false}) {
                const auto carrier = electron ? physics::Carrier::electron : physics::Carrier::hole;
                const std::string set = couderc ? "silicon" : "legacy_silicon";
                const std::string band = electron ? "conduction" : "valence";
                const std::string ctx = std::string(" (") + (couderc ? "Couderc" : "legacy") + ", " +
                                        (rich ? "Richardson" : "DOS") + ", " + (electron ? "electron" : "hole") + ")";
                for (const double T : tvb31::temperature_grid(200.0, 500.0, {150.0, 600.0})) {
                    std::vector<Predicate> pr{in_range(std::string(t_range), T, 200.0, 500.0)};
                    if (couderc) pr.push_back(in_range(std::string(couderc_range), T, 1.0, 400.0));
                    const ReferenceInfo& inf = couderc ? cinf : tinf;
                    const double v = physics::emission_velocity_cm_s(m, carrier, T);
                    const double dv = physics::emission_velocity_slope(m, carrier, T);
                    if (couderc) {
                        using A6 = std::array<DD, 6>;  // A, T, c0..c3
                        A6 x{DD(rich ? A_test : 0.0), DD(T)};
                        for (int k = 0; k < 4; ++k) {
                            x[2 + k] = register_dd(set + ".dos_mass." + band + "[" + std::to_string(k) + "]");
                        }
                        auto f = [](const A6& y) {
                            const std::array<DD, 4> c{y[2], y[3], y[4], y[5]};
                            return mh::emission(y[0], y[1], mb::couderc_dos(c, y[1]),
                                                mb::couderc_dos_log_slope(c, y[1]));
                        };
                        const auto q = f(x);
                        const auto S = verify::reference::condition_sums<2>(f, x);
                        r.comparisons.push_back(value_record(inf, "v" + ctx, "cm/s", T, v, q[0], S[0], pr));
                        r.comparisons.push_back(value_record(inf, "dv/dT" + ctx, "cm/(s K)", T, dv, q[1], S[1], pr));
                    } else {
                        using A3 = std::array<DD, 3>;  // A, T, N300
                        const A3 x{DD(rich ? A_test : 0.0), DD(T),
                                   register_dd(set + (electron ? ".Nc300" : ".Nv300"))};
                        auto f = [](const A3& y) {
                            return mh::emission(y[0], y[1], mb::power_dos(y[2], y[1]), mb::power_dos_log_slope(y[1]));
                        };
                        const auto q = f(x);
                        const auto S = verify::reference::condition_sums<2>(f, x);
                        r.comparisons.push_back(value_record(inf, "v" + ctx, "cm/s", T, v, q[0], S[0], pr));
                        r.comparisons.push_back(value_record(inf, "dv/dT" + ctx, "cm/(s K)", T, dv, q[1], S[1], pr));
                        if (!rich) {
                            // The Schottky overload with the contact's Richardson constant.
                            const A3 xr{DD(A_test), DD(T), x[2]};
                            const auto qr = f(xr);
                            const auto Sr = verify::reference::condition_sums<2>(f, xr);
                            r.comparisons.push_back(value_record(
                                inf, "v Schottky A*" + ctx, "cm/s", T,
                                physics::emission_velocity_cm_s(m, carrier, T, A_test), qr[0], Sr[0], pr));
                        }
                    }
                }
            }
        }
    }
    apply_known_findings(r, known_findings());
    return r;
}

BenchmarkResult run_tvb31_27(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-31.27";
    r.title = "MLDA wavelength, factor and mean factor";
    r.info = tvb31::entry_info(r.id, Provenance::PB,
                               "Paasch and Ubensee 1982 in Hansch et al. 1989's nondegenerate form, as Ieong et al. "
                               "(SISPAD 1998) and Entner's thesis (eqs 2.106-2.108) print it",
                               "lambda = hbar / sqrt(2 m k T); f = 1 - exp(-z^2 / lambda^2)",
                               {std::string(any_mlda), std::string(z_range), std::string(a_range)}, {13},
                               "verify/reference/model_thermal.hpp",
                               "model_thermal_test: the mean against tanh-sinh (TC4) of f, its series against the "
                               "closed form; 15.9's rule");
    ReferenceInfo linf = r.info;
    linf.applicability = {std::string(any_mlda)};
    ReferenceInfo zinf = r.info;
    zinf.applicability = {std::string(z_range)};
    ReferenceInfo minf = r.info;  // the mean factor: an exact consequence (PD), second method quadrature
    minf.applicability = {std::string(a_range)};
    minf.cls = ReferenceClass::E;
    minf.provenance = Provenance::PD;
    minf.citation = "derived here: (1/a) int_0^a f dz = 1 - sqrt(pi) lambda / (2a) erf(a / lambda)";
    minf.second_method = "tanh-sinh quadrature of f over [0, a] (L0)";
    using A2 = std::array<DD, 2>;
    auto fl = [](const A2& y) { return mh::mlda_wavelength(y); };
    auto ff = [](const A2& y) { return mh::mlda_factor(y); };
    auto fm = [](const A2& y) { return mh::mlda_mean(y); };
    for (const double mass : {physics::mlda_silicon_electron_mass, physics::mlda_silicon_hole_mass}) {
        for (const double T : tvb31::temperature_grid(200.0, 500.0, {})) {
            const A2 x{DD(mass), DD(T)};
            r.comparisons.push_back(value_record(linf, "lambda (m " + shortest(mass) + ")", "cm", T,
                                                 physics::mlda_wavelength_cm(mass, T), fl(x),
                                                 verify::reference::condition_sum(fl, x), {holds(any_mlda)}));
        }
    }
    const double lam = physics::mlda_wavelength_cm(physics::mlda_silicon_electron_mass, 300.0);
    for (int i = 0; i <= 768; ++i) {
        const double z = (i / 64.0) * lam;
        const A2 x{DD(z), DD(lam)};
        r.comparisons.push_back(value_record(zinf, "f(z)", "1", z / lam, physics::mlda_factor(z, lam), ff(x),
                                             verify::reference::condition_sum(ff, x),
                                             {in_range(std::string(z_range), z / lam, 0.0, 10.0)}));
    }
    std::vector<double> xs = tvb31::log_grid(1e-8, 1e3, tvb31::per_decade(tier));
    for (const double s : {0.5, std::nextafter(0.5, 0.0), std::nextafter(0.5, 1.0), 1e-9, 1e4}) xs.push_back(s);
    for (const double s : xs) {
        const double a = s * lam;
        const A2 x{DD(a), DD(lam)};
        r.comparisons.push_back(value_record(minf, "mean factor", "1", s, physics::mlda_mean_factor(a, lam), fm(x),
                                             verify::reference::condition_sum(fm, x),
                                             {in_range(std::string(a_range), s, 1e-8, 1e3)}));
    }
    tvb31::table_checks(r, "31.27");
    apply_known_findings(r, known_findings());
    return r;
}

BenchmarkResult run_tvb31_28(Tier) {
    BenchmarkResult r;
    r.id = "TVB-31.28";
    r.title = "Lattice thermal conductivity power law";
    r.info = tvb31::entry_info(r.id, Provenance::PD,
                               "the power law as NiTCAD documents it (legacy materials.kappa_th; MI-Q20), checked as "
                               "written on the register's printed values (Palankovski Table 3.5, PC)",
                               "kappa = kappa_300 (T / 300)^-a; d kappa / dT = -a kappa / T",
                               {std::string(kappa_range)}, {9}, "verify/reference/model_thermal.hpp",
                               "model_thermal_test: d kappa / dT against central differences; 15.9's rule");
    r.info.assumptions.push_back("kappa_300 and a: the register's printed 1.48 W/(cm K) and 1.3 (shipped 1.33)");
    const ReferenceInfo inf = r.info;
    const physics::SemiconductorParameters p =
        tvb31::register_parameters(physics::legacy_silicon_parameters, "legacy_silicon", {"31.28"});
    using A3 = std::array<DD, 3>;
    auto f = [](const A3& y) { return mh::conductivity(y); };
    for (const double T : tvb31::temperature_grid(200.0, 600.0, {150.0, 700.0})) {
        const physics::ValueSlope k = physics::thermal_conductivity(p.thermal, T);
        const A3 x{register_dd("legacy_silicon.thermal.conductivity_W_cmK"),
                   register_dd("legacy_silicon.thermal.conductivity_exponent"), DD(T)};
        const auto q = f(x);
        const auto S = verify::reference::condition_sums<2>(f, x);
        const std::vector<Predicate> pr{in_range(std::string(kappa_range), T, 200.0, 600.0)};
        r.comparisons.push_back(value_record(inf, "kappa", "W/(cm K)", T, k.value, q[0], S[0], pr));
        r.comparisons.push_back(value_record(inf, "dkappa/dT", "W/(cm K^2)", T, k.d_dT, q[1], S[1], pr));
    }
    tvb31::table_checks(r, "31.28");
    apply_known_findings(r, known_findings());
    return r;
}

BenchmarkResult run_tvb31_29(Tier) {
    BenchmarkResult r;
    r.id = "TVB-31.29";
    r.title = "Fermi orders -1, 0, 1 and the thermopower factor";
    r.info = tvb31::entry_info(r.id, Provenance::PD,
                               "derived here from DLMF 25.12 (F_0 = ln(1 + e^x), F_1 = -Li2(-e^x)); h from Wachutka "
                               "1990 and Lundstrom sec. 3.4 with r = -1/2 as NiTCAD documents it",
                               "F_-1 = 1 / (1 + e^-x), F_0, F_1; h = 2 F_1 / F_0 - 3/2, dh/dx = 2 - 2 F_1 F_-1 / F_0^2",
                               {std::string(x_range)}, {9}, "verify/reference/model_thermal.hpp",
                               "model_thermal_test: F_0, F_1 differences against tanh-sinh (TC4) of F_-1, F_0; the "
                               "series against the closed forms; dh/dx against central differences; 15.9's rule");
    r.info.cls = ReferenceClass::E;
    r.info.second_method = "tanh-sinh quadrature of the next lower order; central differences of h (L0)";
    const ReferenceInfo inf = r.info;
    std::vector<double> xs;
    for (int i = -3000; i <= 3000; ++i) xs.push_back(i / 4.0);
    for (int i = -2048; i <= 2048; ++i) xs.push_back(i / 1024.0);
    for (const double s : {-1.0, 1.0}) {
        xs.push_back(std::nextafter(s, -INFINITY));
        xs.push_back(std::nextafter(s, INFINITY));
    }
    xs.push_back(-760.0);
    xs.push_back(760.0);
    tvb31::FloorFinding xiii_floor{"(xiii)", finding_xiii};
    using A1 = std::array<DD, 1>;
    auto fm1 = [](const A1& y) { return mh::fermi_orders(y[0]).minus_one; };
    auto f0 = [](const A1& y) { return mh::fermi_orders(y[0]).zero; };
    auto f1 = [](const A1& y) { return mh::fermi_orders(y[0]).one; };
    auto fh = [](const A1& y) { return mh::thermal_diffusion(y); };
    for (const double xv : xs) {
        const std::vector<Predicate> pr{in_range(std::string(x_range), xv, -750.0, 750.0)};
        const A1 x{DD(xv)};
        const physics::FermiOrdersZeroOne s = physics::fermi_orders_zero_one(xv);
        const physics::ThermalDiffusionFactor h = physics::fermi_dirac_thermal_diffusion(xv);
        const mh::FermiOrders q = mh::fermi_orders(x[0]);
        r.comparisons.push_back(value_record(inf, "F_-1", "1", xv, s.minus_one, q.minus_one,
                                             verify::reference::condition_sum(fm1, x), pr));
        r.comparisons.push_back(
            value_record(inf, "F_0", "1", xv, s.zero, q.zero, verify::reference::condition_sum(f0, x), pr));
        r.comparisons.push_back(
            value_record(inf, "F_1", "1", xv, s.one, q.one, verify::reference::condition_sum(f1, x), pr));
        const auto qh = fh(x);
        const auto S = verify::reference::condition_sums<2>(fh, x);
        r.comparisons.push_back(value_record(inf, "h", "1", xv, h.value, qh[0], S[0], pr));
        r.comparisons.push_back(value_record(inf, "dh/dx", "1", xv, h.d_x, qh[1], S[1], pr));
        // Finding (xiii)'s floor: the shipped 2 - 2 F_1 F_-1 / F_0^2, terms of size 2 each (MI-Q21).
        tvb31::classify_floor<1>(r, {4.0}, xiii_floor);
    }
    apply_known_findings(r, known_findings());
    tvb31::note_if_stale(r, xiii_floor);
    return r;
}

}  // namespace NiTCAD::verify::benchmark
