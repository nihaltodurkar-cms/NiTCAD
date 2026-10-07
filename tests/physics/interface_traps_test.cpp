// Insulators and interface traps (ARCHITECTURE.md section 11, Unit 15b): the insulator set, the
// Fermi occupancy and the trap-band quadrature against 40-digit values, the steady-state occupancy
// reducing to the Fermi function at equilibrium under either statistics, the exact partials against
// finite differences, and the validation.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/interface_traps.hpp"
#include "NiTCAD/physics/recombination.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/statistics.hpp"

using namespace NiTCAD::physics;

namespace {

constexpr double eps = std::numeric_limits<double>::epsilon();

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

}  // namespace

TEST_CASE("insulator: SiO2 and validation") {
    REQUIRE(silicon_dioxide().parameters().eps_r == 3.9);
    REQUIRE(Insulator::create({.eps_r = 7.5}).has_value());
    for (const double bad : {0.0, -1.0, std::numeric_limits<double>::quiet_NaN(),
                             std::numeric_limits<double>::infinity()}) {
        const auto m = Insulator::create({.eps_r = bad});
        REQUIRE_FALSE(m.has_value());
        REQUIRE(m.error().code == NiTCAD::base::ErrorCode::invalid_input);
    }
}

TEST_CASE("interface traps: Fermi occupancy against 40-digit values") {
    // x, f = 1 / (1 + e^x), 1 - f, f (1 - f) (mpmath, 40 digits).
    struct Ref {
        double x, f, g, d;
    };
    constexpr Ref refs[] = {
        {-40.0, 0.9999999999999999957516, 4.248354255291588977281e-18, 4.248354255291588959232e-18},
        {-3.5, 0.9706877692486436811347, 0.02931223075135631886529, 0.02845302387973555983969},
        {-0.25, 0.5621765008857981040273, 0.4378234991142018959727, 0.2461340827375983475853},
        {0.0, 0.5, 0.5, 0.25},
        {0.7, 0.3318122278318338934692, 0.6681877721681661065308, 0.2217128732931090504051},
        {12.0, 6.144174602214717825574e-6, 0.9999938558253977852822, 6.144136851333175325188e-6},
        {45.0, 2.862518580549393644388e-20, 0.9999999999999999999714, 2.862518580549393644306e-20},
    };
    for (const Ref& r : refs) {
        const FermiOccupancy f = fermi_occupancy(r.x);
        CAPTURE(r.x, f.occupied, f.empty, f.d_eta);
        // Each of f and 1 - f is computed without cancellation, so both keep full relative
        // precision in their tails.
        REQUIRE(close(f.occupied, r.f, 4 * eps));
        REQUIRE(close(f.empty, r.g, 4 * eps));
        REQUIRE(close(f.d_eta, r.d, 8 * eps));
    }
}

TEST_CASE("interface traps: the steady-state occupancy is the Fermi function at equilibrium") {
    // Scaled densities as the assemblers use them. Equilibrium: n = gamma_n n_ie e^eta,
    // p = gamma_p n_ie e^-eta with the same gamma in n1 and p1; Boltzmann is gamma = 1. The
    // cross-sections do not enter the equilibrium occupancy.
    const double n_ie = 1.0674e10 / 1e17;
    for (const double gamma_n : {1.0, 0.71, 0.18}) {
        for (const double gamma_p : {1.0, 0.93, 0.42}) {
            for (const double eta : {-20.0, -6.3, -0.4, 0.0, 1.7, 9.0, 25.0}) {
                for (const double tau : {-18.0, -2.5, 0.0, 0.3, 4.0, 19.0}) {
                    for (const auto [cn, cp] : {std::pair{1e-8, 1e-8}, std::pair{3e-6, 2e-10},
                                                std::pair{1e-12, 5e-7}}) {
                        const double n = gamma_n * n_ie * std::exp(eta);
                        const double p = gamma_p * n_ie * std::exp(-eta);
                        const double n1 = gamma_n * n_ie * std::exp(tau);
                        const double p1 = gamma_p * n_ie * std::exp(-tau);
                        const TrapKinetics t = trap_kinetics(n, p, n1, 0.0, p1, 0.0, cn, cp);
                        const FermiOccupancy f = fermi_occupancy(tau - eta);
                        CAPTURE(gamma_n, gamma_p, eta, tau, cn, cp, t.occupied, f.occupied);
                        // Exact in exact arithmetic; the exponentials round separately, so a few
                        // ulps (measured at most 1.5e-15 for f and 1.9e-15 for 1 - f, relative).
                        REQUIRE(close(t.occupied, f.occupied, 4e-15));
                        REQUIRE(close(t.empty, f.empty, 4e-15));
                        // n p - n1 p1 vanishes up to the rounding of the two products.
                        const double D = cn * (n + n1) + cp * (p + p1);
                        REQUIRE(std::abs(t.rate) <= 8 * eps * cn * cp * (n * p + n1 * p1) / D);
                    }
                }
            }
        }
    }
}

TEST_CASE("interface traps: the partials of the steady-state kinetics") {
    // Away from equilibrium, with density-dependent n1 and p1 (Fermi-Dirac: n1 = gamma_n(n) ...).
    // Central differences in n and p with n1(n), p1(p) moved consistently.
    const double n_ie = 1e-7;
    const double tau = 1.3, cn = 2e-8, cp = 7e-9;
    const std::pair<double, double> states[] = {{1.0, 1e-14},  {1e-14, 1.0}, {0.3, 0.2},
                                                {1e-9, 1e-12}, {2.5, 3.1},   {1e-7, 1e-7}};
    for (const auto [n, p] : states) {
        // ln gamma = w density with w chosen so that it changes by a resolvable amount over the
        // difference step (ln gamma -0.21 and -0.08 at the state).
        const double wn = -0.21 / n, wp = -0.08 / p;
        const auto at = [&](double nn, double pp) {
            const double n1 = n_ie * std::exp(wn * nn + tau), p1 = n_ie * std::exp(wp * pp - tau);
            return trap_kinetics(nn, pp, n1, n1 * wn, p1, p1 * wp, cn, cp);
        };
        const TrapKinetics t = at(n, p);
        const double hn = 1e-6 * n, hp = 1e-6 * p;
        const TrapKinetics np = at(n + hn, p), nm = at(n - hn, p);
        const TrapKinetics pp = at(n, p + hp), pm = at(n, p - hp);
        // f differenced through the smaller of f and 1 - f, which keeps its relative precision.
        const auto f = [&](const TrapKinetics& k) {
            return t.occupied > 0.5 ? -k.empty : k.occupied;
        };
        const double dfn = (f(np) - f(nm)) / ((n + hn) - (n - hn));
        const double dfp = (f(pp) - f(pm)) / ((p + hp) - (p - hp));
        const double drn = (np.rate - nm.rate) / ((n + hn) - (n - hn));
        const double drp = (pp.rate - pm.rate) / ((p + hp) - (p - hp));
        CAPTURE(n, p, t.d_occupied_dn, dfn, t.d_occupied_dp, dfp, t.d_rate_dn, drn, t.d_rate_dp,
                drp);
        const double fscale =
            std::max(std::abs(t.d_occupied_dn) * n, std::abs(t.d_occupied_dp) * p);
        REQUIRE(std::abs(dfn - t.d_occupied_dn) * n <= 1e-8 * fscale);
        REQUIRE(std::abs(dfp - t.d_occupied_dp) * p <= 1e-8 * fscale);
        const double rscale = std::max(std::abs(t.d_rate_dn) * n, std::abs(t.d_rate_dp) * p);
        REQUIRE(std::abs(drn - t.d_rate_dn) * n <= 1e-8 * rscale);
        REQUIRE(std::abs(drp - t.d_rate_dp) * p <= 1e-8 * rscale);
        REQUIRE(close(t.occupied + t.empty, 1.0, 4 * eps));
    }
}

TEST_CASE("interface traps: a mid-gap level is the surface recombination velocity form") {
    // N_t c_n = s_n and N_t c_p = s_p: the trap's U equals SRH with tau = 1 / s per unit area.
    const double n_ie = 1e-7, N = 3e10, cn = 1e-8, cp = 4e-9;
    for (const auto [n, p] :
         {std::pair{1.0, 1e-10}, std::pair{1e-3, 0.4}, std::pair{1e-9, 1e-9}}) {
        const TrapKinetics t = trap_kinetics(n, p, n_ie, 0.0, n_ie, 0.0, cn, cp);
        const RecombinationRate s = srh_recombination(
            n, p, boltzmann_equilibrium_product(n_ie), n_ie, 1.0 / (N * cn), 1.0 / (N * cp));
        CAPTURE(n, p, N * t.rate, s.rate);
        REQUIRE(close(N * t.rate, s.rate, 8 * eps));
        // The partials are differences of nearly equal terms where one density is tiny: compared
        // on the scale of the larger partial.
        const double scale = std::max(std::abs(s.d_dn), std::abs(s.d_dp));
        REQUIRE(std::abs(N * t.d_rate_dn - s.d_dn) <= 16 * eps * scale);
        REQUIRE(std::abs(N * t.d_rate_dp - s.d_dp) <= 16 * eps * scale);
    }
}

TEST_CASE("interface traps: band quadrature against the closed-form equilibrium integral") {
    // A uniform band from -0.55 to 0.53 eV about E_i at 300 K: 42 panels of 6 levels. The
    // occupied density at Fermi level E_F, D_it int f dE, has the closed form
    // D_it kT ln((1 + e^((E_F - lo)/kT)) / (1 + e^((E_F - hi)/kT))) (mpmath, 40 digits).
    const TrapBand band{.type = TrapType::acceptor,
                        .density_cm2_eV = 2e11,
                        .energy_low_eV = -0.55,
                        .energy_high_eV = 0.53};
    const double T = 300.0;
    const double kT = NiTCAD::base::thermal_voltage(T);
    REQUIRE(close(kT, 0.02585199978643553230099098, 2 * eps));
    const std::vector<TrapLevel> levels = trap_band_levels(band, T);
    REQUIRE(levels.size() == 6 * static_cast<std::size_t>(std::ceil(1.08 / kT)));
    double total = 0.0;
    for (const TrapLevel& l : levels) {
        REQUIRE(l.type == TrapType::acceptor);
        REQUIRE(l.energy_eV > band.energy_low_eV);
        REQUIRE(l.energy_eV < band.energy_high_eV);
        total += l.density_cm2;
    }
    REQUIRE(close(total, 2e11 * 1.08, 1e-14));
    struct Ref {
        double EF, integral;
    };
    constexpr Ref refs[] = {
        {-0.6, 0.00349046036301059988093},  {-0.3, 0.2500016317938146783758},
        {0.0, 0.5499999999826137404589},    {0.123, 0.672999996240065518931},
        {0.5, 1.042953339583036665382},     {0.7, 1.07996399876439887113},
    };
    double worst = 0.0;
    for (const Ref& r : refs) {
        double occupied = 0.0;
        for (const TrapLevel& l : levels) {
            occupied += l.density_cm2 * fermi_occupancy((l.energy_eV - r.EF) / kT).occupied;
        }
        const double exact = 2e11 * r.integral;
        worst = std::max(worst, std::abs(occupied - exact) / (2e11 * 1.08));
        CAPTURE(r.EF, occupied, exact);
    }
    CAPTURE(worst);
    // Relative to the band's whole density; measured 7.1e-16.
    REQUIRE(worst < 1e-14);
}

TEST_CASE("interface traps: validation") {
    const Semiconductor si = silicon();
    const double T = 300.0;
    REQUIRE(check_interface_traps({}, si, T).has_value());
    const auto rejected = [&](const InterfaceTraps& t) {
        const auto ok = check_interface_traps(t, si, T);
        if (ok) return false;
        return ok.error().code == NiTCAD::base::ErrorCode::invalid_input;
    };
    const double nan = std::numeric_limits<double>::quiet_NaN();
    REQUIRE(rejected({.levels = {{.density_cm2 = -1.0}}}));
    REQUIRE(rejected({.levels = {{.density_cm2 = nan}}}));
    REQUIRE(rejected({.levels = {{.density_cm2 = 1e10, .energy_eV = 0.6}}}));   // above E_c
    REQUIRE(rejected({.levels = {{.density_cm2 = 1e10, .energy_eV = -0.6}}}));  // below E_v
    REQUIRE(rejected({.levels = {{.density_cm2 = 1e10, .sigma_n_cm2 = 0.0}}}));
    REQUIRE(rejected({.levels = {{.density_cm2 = 1e10, .sigma_p_cm2 = -1e-15}}}));
    REQUIRE(rejected({.levels = {{.type = static_cast<TrapType>(7), .density_cm2 = 1e10}}}));
    REQUIRE(rejected({.bands = {{.density_cm2_eV = 1e11, .energy_low_eV = 0.1,
                                 .energy_high_eV = 0.1}}}));
    REQUIRE(rejected({.bands = {{.density_cm2_eV = 1e11, .energy_low_eV = -0.7,
                                 .energy_high_eV = 0.1}}}));
    REQUIRE(rejected({.bands = {{.density_cm2_eV = -1e11, .energy_low_eV = -0.1,
                                 .energy_high_eV = 0.1}}}));
    REQUIRE(rejected({.thermal_velocity_n_cm_s = 0.0}));
    REQUIRE(rejected({.thermal_velocity_p_cm_s = nan}));
    // The gap edges themselves are inside: E_c - E_i and E_v - E_i of silicon at 300 K.
    const double ec = intrinsic_level_depth_eV(si, T) - si.parameters().electron_affinity_eV;
    const double ev = ec - band_gap_eV(si, T);
    REQUIRE(check_interface_traps({.levels = {{.density_cm2 = 1e10, .energy_eV = ec},
                                              {.density_cm2 = 1e10, .energy_eV = ev}},
                                   .bands = {{.density_cm2_eV = 1e11, .energy_low_eV = ev,
                                              .energy_high_eV = ec}}},
                                  si, T)
                .has_value());
}
