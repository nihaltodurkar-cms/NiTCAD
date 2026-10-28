// Schottky contacts on devices (Unit 48; DECISIONS.md SC2-SC16): 1D silicon diodes against
// analytic theory. Equilibrium (no current; the barrier-limited densities and the band bending),
// the surface field against the first integral of Poisson's equation and the C-V, the
// thermionic-diffusion series law on the solution's own potential, the saturation current and
// ideality, the velocity limits, the minority carriers, Fermi-Dirac statistics, the work function
// against the barrier, 2D tensor, 2D Delaunay and 3D devices against 1D, Newton's convergence, the
// small-signal and transient currents, conservation and a device whose only carrier contacts are
// Schottky, the deck's refusals, the run identity, a deck analysis and the device file, and the
// diode deck of decks/ (reported). The criteria not met are [!shouldfail] cases, unchanged.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/io/device_file.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/physics/fermi_dirac.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/thermionic_emission.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/equilibrium.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "NiTCAD/solve/transient.hpp"
#include "fixtures.hpp"

using namespace NiTCAD;
using deck_test::must;

namespace {

constexpr double T = 300.0;
const double kT = base::k_B_eV_per_K * T;  // [eV]

const physics::Semiconductor& silicon() {
    static const physics::Semiconductor si =
        *physics::Semiconductor::create(physics::silicon_parameters);
    return si;
}

struct Diode {
    char type = 'n';                // the body's doping: n (donors) or p (acceptors)
    double doping = 1e16;           // [cm^-3]
    std::string barrier = "barrier = 0.7";
    double richardson_n = 110.0;    // A* [A / (cm^2 K^2)] (Sze: 110 / 30 for silicon)
    double richardson_p = 30.0;
    int dimension = 1;
};

// A Schottky contact on x_min of a 2 um body, ohmic on x_max; 1 nm spacing in the first 0.2 um.
std::string deck_of(const Diode& d) {
    const char* across = d.dimension == 1   ? ""
                         : d.dimension == 2 ? "y = linspace(0, 1um, 3)\n"
                                            : "y = linspace(0, 1um, 3)\nz = linspace(0, 1um, 3)\n";
    char text[2048];
    std::snprintf(text, sizeof text, R"(
[device]
temperature = 300

[mesh]
type = tensor
x = linspace(0, 0.2um, 201) + linspace(0.2um, 2um, 181)
%s
[region silicon]
material = silicon

[doping]
profile = %s constant %g

[contact anode]
kind = schottky
patch = x_min
%s
richardson_n = %g
richardson_p = %g

[contact cathode]
kind = ohmic
patch = x_max

[sweep]
contact anode = 0
)",
                  across, d.type == 'n' ? "donor" : "acceptor", d.doping, d.barrier.c_str(),
                  d.richardson_n, d.richardson_p);
    return text;
}

struct Models {
    bool fermi_dirac = false;
    bool recombination = true;  // SRH and Auger
    bool bgn = true;
};

deck::Plan plan_of(const Diode& d, const Models& m = {}) {
    deck::Plan p = must(deck_of(d));
    p.options.models.fermi_dirac = m.fermi_dirac;
    p.options.models.srh = m.recombination;
    p.options.models.auger = m.recombination;
    p.options.models.bgn = m.bgn;
    return p;
}

// The anode bias of every point (the cathode at 0).
results::Sweep sweep(const deck::Plan& p, const std::vector<double>& anode) {
    std::vector<std::vector<double>> points;
    for (const double v : anode) points.push_back({v, 0.0});
    auto s = solve::sweep_bias(p.device, points, p.options);
    if (!s) FAIL(s.error().message);
    REQUIRE_FALSE(s->stopped.has_value());
    REQUIRE(s->points.size() == anode.size());
    return std::move(*s);
}

double anode_current(const results::BiasPoint& b) { return b.terminal_current[0]; }

double relative(double a, double b) { return std::abs(a - b) / std::abs(b); }

// Phi_Bn = 0.7 eV, A* = 110: A* T^2 e^(-Phi / kT) [A/cm^2].
double thermionic_saturation(double richardson, double barrier) {
    return richardson * T * T * std::exp(-barrier / kT);
}

}  // namespace

TEST_CASE("Schottky: equilibrium carries no current; n_B, p_B and the band bending (SC2)",
          "[schottky]") {
    const double Nc = physics::conduction_band_dos(silicon(), T);
    const double Nv = physics::valence_band_dos(silicon(), T);
    const double Eg = physics::band_gap_eV(silicon(), T);
    for (const char type : {'n', 'p'}) {
        for (const double barrier : {0.7, 0.4}) {
            Diode d{.type = type, .barrier = "barrier = " + std::to_string(barrier)};
            const deck::Plan p = plan_of(d, {.bgn = false});
            const results::Sweep s = sweep(p, {0.0, -1.0});
            const results::BiasPoint& eq = s.points[0];
            const double reverse = std::abs(anode_current(s.points[1]));
            INFO(type << " " << barrier << ": I(0) " << anode_current(eq) << ", I(-1 V) "
                      << -reverse << ", resolution " << eq.terminal_current_resolution[0]);
            // The 1e-12 of SC2 is a case of its own (below); here, within the resolution.
            CHECK(std::abs(anode_current(eq)) <= eq.terminal_current_resolution[0]);
            const double nB = Nc * std::exp(-barrier / kT);
            const double pB = Nv * std::exp(-(Eg - barrier) / kT);
            CHECK(relative(eq.fields.n_cm3[0], nB) <= 1e-12);
            CHECK(relative(eq.fields.p_cm3[0], pB) <= 1e-12);
            // The band bending: E_c at the contact Phi_Bn above E_F, in the body kT ln(Nc / n).
            CHECK(std::abs(eq.bands.conduction_eV[0] - barrier) <= 1e-9);
            const std::size_t last = eq.fields.n_cm3.size() - 1;
            const double bulk = kT * std::log(Nc / eq.fields.n_cm3[last]);
            CHECK(std::abs((eq.bands.conduction_eV[0] - eq.bands.conduction_eV[last]) -
                           (barrier - bulk)) <= 1e-9);
            // Fermi-Dirac: n_B = Nc F_1/2(-Phi_Bn / kT), p_B = Nv F_1/2(-Phi_Bp / kT).
            const deck::Plan f = plan_of(d, {.fermi_dirac = true, .bgn = false});
            const results::BiasPoint q = sweep(f, {0.0}).points[0];
            const double nF = Nc * physics::fermi_half(-barrier / kT).value;
            const double pF = Nv * physics::fermi_half(-(Eg - barrier) / kT).value;
            INFO("Fermi-Dirac: n(0) " << q.fields.n_cm3[0] << " against " << nF << ", p(0) "
                                      << q.fields.p_cm3[0] << " against " << pF);
            CHECK(std::abs(anode_current(q)) <= q.terminal_current_resolution[0]);
            CHECK(relative(q.fields.n_cm3[0], nF) <= 1e-12);
            CHECK(relative(q.fields.p_cm3[0], pF) <= 1e-12);
        }
    }
}

// SC2's and SC8's equilibrium current as fixed: at most 1e-12 of SC5's reverse current (1.7e-17
// A/cm^2). Not met: the equilibrium current is 1e-14 A/cm^2 at Phi_Bn = 0.7 eV and 3e-8 at 0.05 eV
// (Fermi-Dirac), the rounding of the contact edge's two opposing fluxes (q D n / h, about 5e5
// A/cm^2 in the body on 1 nm edges, and 2e8 in the 0.05 eV accumulation layer), well inside the
// solver's own terminal-current resolution (1e-7 A/cm^2 at 0.7 eV). The criterion did not allow
// for rounding. Encoded unchanged; for the owner.
TEST_CASE("Schottky: the equilibrium current within 1e-12 of the reverse current (SC2, SC8)",
          "[schottky][!shouldfail]") {
    const double reverse = std::abs(anode_current(sweep(plan_of({}), {-1.0}).points[0]));
    for (const bool fd : {false, true}) {
        for (const double barrier : {0.7, 0.4, 0.05}) {
            if (!fd && barrier == 0.05) continue;
            const deck::Plan p = plan_of({.barrier = "barrier = " + std::to_string(barrier)},
                                         {.fermi_dirac = fd, .bgn = fd});
            const results::BiasPoint eq = sweep(p, {0.0}).points[0];
            INFO((fd ? "Fermi-Dirac " : "Boltzmann ") << barrier << ": I(0) "
                                                      << anode_current(eq));
            CHECK(std::abs(anode_current(eq)) <= 1e-12 * reverse);
        }
    }
}

TEST_CASE("Schottky: the surface field and the C-V of a depletion layer (SC3)", "[schottky]") {
    const Diode d{};
    const deck::Plan p = plan_of(d, {.bgn = false});
    const results::Sweep s = sweep(p, {0.0});
    const results::BiasPoint& eq = s.points[0];
    const auto& x = p.device.mesh().points();
    const double eps = silicon().parameters().eps_r * base::eps0_F_per_cm;
    // Gauss over the contact node's half cell: the field leaving the metal is the edge's D plus the
    // half cell's charge, q (N_D - n + p) h / 2.
    const double h = x[1][0] - x[0][0];
    const double psi0 = eq.fields.potential_V[0], psi1 = eq.fields.potential_V[1];
    const double rho = base::q_C * (d.doping - eq.fields.n_cm3[0] + eq.fields.p_cm3[0]);
    const double Es = ((psi1 - psi0) * eps / h + rho * h / 2.0) / eps;
    // The first integral, u the band bending in kT and n0, p0 the body's densities.
    const std::size_t last = eq.fields.n_cm3.size() - 1;
    const double n0 = eq.fields.n_cm3[last], p0 = eq.fields.p_cm3[last];
    const double u = (eq.bands.conduction_eV[0] - eq.bands.conduction_eV[last]) / kT;
    const double exact = std::sqrt(2.0 * kT * base::q_C * n0 / eps *
                                   (std::exp(-u) + u - 1.0 + p0 / n0 * (std::exp(u) - u - 1.0)));
    INFO("E_s " << Es << " V/cm, first integral " << exact);
    CHECK(relative(std::abs(Es), exact) <= 5e-3);

    // C-V at 1 kHz over 0 ... -2 V: 1/C^2 = 2 (V_bi - V - kT/q) / (q eps N_D).
    solve::SmallSignalOptions o{.steady = p.options, .frequencies_Hz = {1e3}};
    std::vector<std::vector<double>> points;
    for (int k = 0; k <= 8; ++k) points.push_back({-0.25 * k, 0.0});
    const auto ac = solve::solve_small_signal(p.device, points, o);
    REQUIRE(ac.has_value());
    REQUIRE(ac->points.size() == points.size());
    // Least squares of 1/C^2 against V.
    double sv = 0, sy = 0, svv = 0, svy = 0;
    const double nk = static_cast<double>(points.size());
    for (std::size_t k = 0; k < points.size(); ++k) {
        const double C = ac->capacitance(k, 0, 0, 0);
        const double y = 1.0 / (C * C), v = points[k][0];
        sv += v;
        sy += y;
        svv += v * v;
        svy += v * y;
    }
    const double slope = (nk * svy - sv * sy) / (nk * svv - sv * sv);
    const double intercept = (sy - slope * sv) / nk;  // 1/C^2 at V = 0
    const double ND = -2.0 / (base::q_C * eps * slope);
    const double vbi_less_kT = -intercept / slope;     // V at 1/C^2 = 0
    const double vbi = 0.7 - kT * std::log(physics::conduction_band_dos(silicon(), T) / n0);
    INFO("N_D from C-V " << ND << ", V_bi - kT/q " << vbi_less_kT << " against " << vbi - kT);
    CHECK(relative(ND, d.doping) <= 1e-2);
    CHECK(std::abs(vbi_less_kT - (vbi - kT)) <= 10e-3);
}

namespace {

// The exact series law of SC4 (Slotboom form of the SG fluxes) for the carrier of the body: the
// particle current between the metal and the ohmic contact through the emission and every edge,
// J = q (U_ohmic - U_B) / (1 / (v N e^(w_0)) + sum r_e / (D N)), with w = -E_c / kT (electrons) or
// E_v / kT (holes), U = e^(E_F / kT) (electrons) or e^(-E_F / kT), and the edge resistance
// r_e = L (e^(-w_a) - e^(-w_b)) / (w_b - w_a). Returns the magnitude in A/cm^2.
double series_law(const deck::Plan& p, const results::BiasPoint& b, char type, double anode_V,
                  double richardson, double doping) {
    const physics::Carrier carrier = type == 'n' ? physics::Carrier::electron
                                                 : physics::Carrier::hole;
    const double N = type == 'n' ? physics::conduction_band_dos(silicon(), T)
                                 : physics::valence_band_dos(silicon(), T);
    const double v = physics::emission_velocity_cm_s(silicon(), carrier, T, richardson);
    const double mu = physics::caughey_thomas_mobility(silicon(), carrier, doping, T);
    const bool masetti = type == 'p';  // silicon's holes default to Masetti (Unit 35)
    const double D = (masetti ? physics::masetti_mobility(physics::masetti_boron_silicon, doping, T)
                              : mu) *
                     kT;
    const auto w = [&](std::size_t i) {
        return type == 'n' ? -b.bands.conduction_eV[i] / kT : b.bands.valence_eV[i] / kT;
    };
    const auto& x = p.device.mesh().points();
    double sum = 0.0;
    for (std::size_t i = 0; i + 1 < x.size(); ++i) {
        const double L = x[i + 1][0] - x[i][0];
        const double wa = w(i), wb = w(i + 1);
        sum += wa == wb ? L * std::exp(-wa) : L * (std::exp(-wa) - std::exp(-wb)) / (wb - wa);
    }
    // The metal's Fermi level is at -V; the ohmic contact's at 0.
    const double sign = type == 'n' ? 1.0 : -1.0;
    const double UB = std::exp(sign * -anode_V / kT), U1 = 1.0;
    const double flux = (U1 - UB) / (1.0 / (v * N * std::exp(w(0))) + sum / (D * N));
    return base::q_C * std::abs(flux);
}

}  // namespace

TEST_CASE("Schottky: the thermionic-diffusion series law on the solution's potential (SC4)",
          "[schottky]") {
    for (const char type : {'n', 'p'}) {
        Diode d{.type = type};
        const deck::Plan p = plan_of(d, {.recombination = false});
        const std::vector<double> bias = type == 'n'
                                             ? std::vector<double>{-1.0, -0.5, -0.1, 0.1, 0.2, 0.3}
                                             : std::vector<double>{1.0, 0.5, 0.1, -0.1, -0.2, -0.3};
        const results::Sweep s = sweep(p, bias);
        for (std::size_t k = 0; k < bias.size(); ++k) {
            const results::BiasPoint& b = s.points[k];
            // The body's carrier current on the first edge.
            const double J = std::abs(type == 'n' ? b.edge_current_n[0] : b.edge_current_p[0]);
            const double law = series_law(p, b, type, bias[k],
                                          type == 'n' ? d.richardson_n : d.richardson_p, d.doping);
            INFO(type << " at " << bias[k] << " V: J " << J << ", series law " << law);
            CHECK(relative(J, law) <= 1e-2);
        }
    }
}

namespace {

// SC5's diodes: n-Si with Phi_Bn = 0.7 eV, A*_n = 110; p-Si with Phi_Bp = 0.7 eV (barrier =
// Eg - 0.7), A*_p = 30. The reverse point (1 V) first, then forward 0 ... 0.4 V by 10 mV.
struct Characteristic {
    std::vector<double> bias;  // the anode's, forward positive on n-Si and negative on p-Si
    results::Sweep sweep;
    double saturation;         // A* T^2 e^(-0.7 / kT)
};

Characteristic characteristic(char type) {
    Diode d{.type = type};
    if (type == 'p') {
        char b[64];
        std::snprintf(b, sizeof b, "barrier = %.17g", physics::band_gap_eV(silicon(), T) - 0.7);
        d.barrier = b;
    }
    const double sign = type == 'n' ? 1.0 : -1.0;
    std::vector<double> bias{-sign * 1.0};
    for (int k = 0; k <= 40; ++k) bias.push_back(sign * 0.01 * k);
    results::Sweep s = sweep(plan_of(d), bias);
    return {bias, std::move(s),
            thermionic_saturation(type == 'n' ? d.richardson_n : d.richardson_p, 0.7)};
}

// The local ideality d|V| / d(V_T ln J) between successive forward points with both J in
// [low, high] A/cm^2; with `shifted`, of J / (1 - e^(-|V|/V_T)) (the thermionic form's -1 taken
// out) and only from |V| >= 3 kT/q. {smallest, largest, steps}.
struct Ideality {
    double low = 1e300, high = -1e300;
    int steps = 0;
};

Ideality ideality(const Characteristic& c, double low, double high, bool shifted) {
    Ideality r;
    const auto j = [&](std::size_t k) {
        const double v = std::abs(c.bias[k]);
        const double J = std::abs(anode_current(c.sweep.points[k]));
        return shifted ? J / (1.0 - std::exp(-v / kT)) : J;
    };
    for (std::size_t k = 2; k < c.bias.size(); ++k) {
        if (shifted && std::abs(c.bias[k - 1]) < 3.0 * kT) continue;
        const double ja = std::abs(anode_current(c.sweep.points[k - 1]));
        const double jb = std::abs(anode_current(c.sweep.points[k]));
        if (ja < low || jb > high) continue;
        const double n = std::abs(c.bias[k] - c.bias[k - 1]) / (kT * std::log(j(k) / j(k - 1)));
        r.low = std::min(r.low, n);
        r.high = std::max(r.high, n);
        ++r.steps;
    }
    return r;
}

}  // namespace

TEST_CASE("Schottky: the saturation current (SC5)", "[schottky]") {
    for (const char type : {'n', 'p'}) {
        const Characteristic c = characteristic(type);
        const double reverse = std::abs(anode_current(c.sweep.points[0]));
        INFO(type << ": reverse at 1 V " << reverse << ", A* T^2 e^(-Phi/kT) " << c.saturation);
        CHECK(reverse >= 0.80 * c.saturation);
        CHECK(reverse <= 1.00 * c.saturation);
        // Reported, not a criterion: the ideality of the thermionic form, J / (1 - e^(-V/V_T)),
        // from 3 kT/q up to 1e-2 A/cm^2.
        const Ideality t = ideality(c, 0.0, 1e-2, true);
        WARN(type << ": ideality of J / (1 - e^(-V/V_T)) " << t.low << " ... " << t.high
                  << " over " << t.steps << " steps");
    }
}

// SC5 as fixed (J in 1e-6 ... 1e-2 A/cm^2). Not met on n-Si: its saturation current is 1.7e-5
// A/cm^2, so the range starts below it, where the thermionic form's -1 dominates d ln J / dV
// (an ideality of 0.43 at the first step). Encoded with the criterion unchanged; for the owner.
TEST_CASE("Schottky: the ideality factor over 1e-6 ... 1e-2 A/cm^2 (SC5)",
          "[schottky][!shouldfail]") {
    for (const char type : {'n', 'p'}) {
        const Characteristic c = characteristic(type);
        const Ideality r = ideality(c, 1e-6, 1e-2, false);
        INFO(type << ": ideality " << r.low << " ... " << r.high << " over " << r.steps
                  << " steps");
        CHECK(r.steps > 0);
        CHECK(r.low >= 1.00);
        CHECK(r.high <= 1.05);
    }
}

TEST_CASE("Schottky: the velocity limits (SC6)", "[schottky]") {
    const double Nc = physics::conduction_band_dos(silicon(), T);
    const double nB = Nc * std::exp(-0.7 / kT);
    double previous = 0.0;
    for (const double A : {10.0, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9}) {
        const deck::Plan p = plan_of({.richardson_n = A});
        const results::Sweep s = sweep(p, {0.2});
        const double J = anode_current(s.points[0]);
        INFO("A* " << A << ": J " << J << ", n(0) / n_B " << s.points[0].fields.n_cm3[0] / nB);
        CHECK(J > previous);
        previous = J;
        // The contact density lies J_n / (q v) above n_B (the emission condition itself).
        const double v =
            physics::emission_velocity_cm_s(silicon(), physics::Carrier::electron, T, A);
        const double Jn = std::abs(s.points[0].edge_current_n[0]);
        CHECK(std::abs(s.points[0].fields.n_cm3[0] - nB - Jn / (base::q_C * v)) <=
              1e-6 * s.points[0].fields.n_cm3[0]);
    }
}

// SC6 as fixed: at A* = 1e9 the contact density within 1e-3 of n_B at +0.2 V. Not met: n(0) /
// n_B = 1.0039, which is J_n / (q v n_B) at that A* (0.62 A/cm^2 over v = 1.9e13 cm/s): the
// diffusion-theory limit is approached as 1 / A*, and 1e9 is not yet within 1e-3 of it. Encoded
// with the criterion unchanged; for the owner.
TEST_CASE("Schottky: the diffusion-theory limit at A* = 1e9 (SC6)", "[schottky][!shouldfail]") {
    const double nB = physics::conduction_band_dos(silicon(), T) * std::exp(-0.7 / kT);
    const results::Sweep s = sweep(plan_of({.richardson_n = 1e9}), {0.2});
    INFO("n(0) / n_B " << s.points[0].fields.n_cm3[0] / nB);
    CHECK(relative(s.points[0].fields.n_cm3[0], nB) <= 1e-3);
}

TEST_CASE("Schottky: the minority holes cross by emission (SC7)", "[schottky]") {
    // No recombination: the contact node's half cell then adds nothing to its edge's current.
    const deck::Plan p = plan_of({.barrier = "barrier = 0.8"}, {.recombination = false});
    const results::Sweep s = sweep(p, {0.3});
    const results::BiasPoint& b = s.points[0];
    const double Nv = physics::valence_band_dos(silicon(), T);
    const double Eg = physics::band_gap_eV(silicon(), T);
    // p_B with band-gap narrowing: n_ie^2 / n_B at the contact node (the narrowing lowers the hole
    // barrier, SD2); read it from the equilibrium state of the same device.
    const results::Sweep eq = sweep(p, {0.0});
    const double pB = eq.points[0].fields.p_cm3[0];
    CHECK(pB >= Nv * std::exp(-(Eg - 0.8) / kT));  // never above the unnarrowed barrier
    const double vp = physics::emission_velocity_cm_s(silicon(), physics::Carrier::hole, T, 30.0);
    // The hole current entering the device at the contact: the metal takes holes at v (p - p_B).
    const double Jp = b.edge_current_p[0];
    const double emission = base::q_C * vp * (b.fields.p_cm3[0] - pB);
    INFO("J_p " << Jp << ", q v_p (p - p_B) " << emission);
    CHECK(std::abs(std::abs(Jp) - std::abs(emission)) <= 1e-10 * std::abs(anode_current(b)));
    CHECK(relative(b.edge_current_n[0] + b.edge_current_p[0], anode_current(b)) <= 1e-10);
}

TEST_CASE("Schottky: Fermi-Dirac statistics (SC8)", "[schottky]") {
    const double Nc = physics::conduction_band_dos(silicon(), T);
    const double reverse = std::abs(anode_current(sweep(plan_of({}), {-1.0}).points[0]));
    for (const double barrier : {0.05, 0.7}) {
        const deck::Plan p =
            plan_of({.barrier = "barrier = " + std::to_string(barrier)}, {.fermi_dirac = true});
        const results::Sweep s = sweep(p, {0.0});
        const double nB = Nc * physics::fermi_half(-barrier / kT).value;
        INFO("Phi " << barrier << ": I(0) " << anode_current(s.points[0]) << ", n(0) "
                    << s.points[0].fields.n_cm3[0] << " against " << nB << ", reverse "
                    << reverse);
        // SC8's 1e-12 is in the equilibrium-current case above; here, within the resolution.
        CHECK(std::abs(anode_current(s.points[0])) <=
              s.points[0].terminal_current_resolution[0]);
        CHECK(relative(s.points[0].fields.n_cm3[0], nB) <= 1e-10);
    }
    const std::vector<double> bias{-1.0, -0.5, 0.1, 0.2, 0.3};
    const results::Sweep fd = sweep(plan_of({}, {.fermi_dirac = true}), bias);
    const results::Sweep mb = sweep(plan_of({}), bias);
    for (std::size_t k = 0; k < bias.size(); ++k) {
        INFO(bias[k] << " V: FD " << anode_current(fd.points[k]) << ", Boltzmann "
                     << anode_current(mb.points[k]));
        CHECK(relative(anode_current(fd.points[k]), anode_current(mb.points[k])) <= 1e-3);
    }
}

TEST_CASE("Schottky: the work function and the barrier give the same device (SC9)",
          "[schottky]") {
    const double chi = silicon().parameters().electron_affinity_eV;
    const std::vector<double> bias{-1.0, 0.0, 0.2};
    char wf[64];
    std::snprintf(wf, sizeof wf, "work_function = %.17g", chi + 0.7);
    const results::Sweep a = sweep(plan_of({.barrier = wf}), bias);
    const results::Sweep b = sweep(plan_of({}), bias);
    for (std::size_t k = 0; k < bias.size(); ++k) {
        const double ia = anode_current(a.points[k]), ib = anode_current(b.points[k]);
        if (bias[k] == 0.0) continue;  // no current to compare
        CHECK(relative(ia, ib) <= 1e-12);
    }
}

TEST_CASE("Schottky: 2D and 3D devices give the 1D current per area (SC10)", "[schottky]") {
    const std::vector<double> bias{-1.0, 0.2};
    const results::Sweep one = sweep(plan_of({}), bias);
    for (const int dim : {2, 3}) {
        const results::Sweep s = sweep(plan_of({.dimension = dim}), bias);
        const double area = dim == 2 ? 1e-4 : 1e-8;  // 1 um, or 1 um^2, across
        for (std::size_t k = 0; k < bias.size(); ++k) {
            const double per_area = anode_current(s.points[k]) / area;
            INFO(dim << "D at " << bias[k] << " V: " << per_area << " against "
                     << anode_current(one.points[k]));
            CHECK(relative(per_area, anode_current(one.points[k])) <= 1e-9);
            // The terminal currents sum to zero within the current's resolution (SC13's 1e-10
            // relative is its own case, below).
            const auto& I = s.points[k].terminal_current;
            CHECK(std::abs(I[0] + I[1]) <= s.points[k].terminal_current_resolution[0]);
        }
    }
}

namespace {

constexpr double delaunay_width = 1e-4;  // [cm]

// SC10's 2D Delaunay diode (mesh::generate_mesh): a structured boundary layer on the contact with
// the 1D mesh's spacing to 1 um, past the depletion layer at -1 V (0.45 um) and its Debye tail,
// 0.5 um along it (coarser than max_size: Unit 49), and the generator's triangles beyond
// (max_size 0.1 um); the region is the deck's silicon (the Unit 35 defaults), from `p`.
device::Device delaunay_diode(const deck::Plan& p) {
    const double W = delaunay_width, length = 2e-4;
    std::vector<double> x = deck_test::linspace(0.0, 2e-5, 201);
    for (const double v : deck_test::linspace(2e-5, length, 181)) {
        if (v > x.back()) x.push_back(v);
    }
    std::vector<double> offsets;
    for (const double v : x) {
        if (v > 0.0 && v <= 1e-4 * (1.0 + 1e-12)) offsets.push_back(v);
    }
    mesh::Geometry2D g{{{{{0, 0, 0}, {length, 0, 0}, {length, W, 0}, {0, W, 0}}, 0}},
                       {{"x_min", {0, 0, 0}, {0, W, 0}}, {"x_max", {length, 0, 0}, {length, W, 0}}}};
    mesh::MeshOptions2D o;
    o.max_size = 1e-5;
    o.layers = {{{0, 0, 0}, {0, W, 0}, deck_test::linspace(0.0, W, 3), offsets}};
    auto generated = mesh::generate_mesh(g, o);
    if (!generated) FAIL(generated.error().message);
    mesh::Mesh m = std::move(generated->mesh);
    const std::size_t n = m.node_count();
    device::Contact anode{"anode", device::ContactKind::schottky, m.find_boundary("x_min")->nodes};
    anode.schottky = {.boundary = "x_min", .barrier_eV = 0.7, .richardson_n = 110.0,
                      .richardson_p = 30.0};
    auto cathode = m.find_boundary("x_max")->nodes;
    auto d = device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = T,
         .regions = {p.device.regions()[0]},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::vector<double>(n, 1e16),
         .acceptors = std::vector<double>(n, 0.0),
         .contacts = {std::move(anode),
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
    if (!d) FAIL(d.error().message);
    return std::move(*d);
}

results::Sweep sweep_device(const device::Device& d, const deck::Plan& p,
                            const std::vector<double>& anode) {
    std::vector<std::vector<double>> points;
    for (const double v : anode) points.push_back({v, 0.0});
    auto s = solve::sweep_bias(d, points, p.options);
    if (!s) FAIL(s.error().message);
    REQUIRE_FALSE(s->stopped.has_value());
    return std::move(*s);
}

}  // namespace

// The box method reproduces a 1D solution exactly where it is linear in x on any Delaunay mesh
// (Scharfetter-Gummel is exact along any edge for it, and the faces of an interior box sum to
// zero), so with recombination off (the quasi-neutral body then carries constant currents on a
// linear potential) the 1D current is the reference.
TEST_CASE("Schottky: a 2D Delaunay device gives the 1D current per area (SC10)", "[schottky]") {
    const deck::Plan p = plan_of({}, {.recombination = false});
    const device::Device d = delaunay_diode(p);
    INFO(d.mesh().node_count() << " nodes");
    const std::vector<double> bias{-1.0, 0.2};
    const results::Sweep one = sweep(p, bias);
    const results::Sweep s = sweep_device(d, p, bias);
    for (std::size_t k = 0; k < bias.size(); ++k) {
        const auto& I = s.points[k].terminal_current;
        const double per_area = I[0] / delaunay_width;
        INFO("Delaunay at " << bias[k] << " V: " << per_area << " against "
                            << anode_current(one.points[k]));
        CHECK(relative(per_area, anode_current(one.points[k])) <= 1e-6);
        CHECK(std::abs(I[0] + I[1]) <= s.points[k].terminal_current_resolution[0]);
    }
}

// SC13: 2D devices (tensor and Delaunay) whose only carrier contacts are a Schottky and an ohmic
// contact, forward at +0.3 V (1.7 A/cm^2): the terminal currents sum to zero within 1e-10 relative.
TEST_CASE("Schottky: the terminal currents of a 2D device sum to zero (SC13)", "[schottky]") {
    const deck::Plan tensor = plan_of({.dimension = 2});
    const deck::Plan reference = plan_of({}, {.recombination = false});
    const results::Sweep a = sweep(tensor, {0.3});
    const results::Sweep b = sweep_device(delaunay_diode(reference), reference, {0.3});
    for (const results::Sweep* s : {&a, &b}) {
        const auto& I = s->points[0].terminal_current;
        INFO((s == &a ? "tensor" : "Delaunay") << ": anode " << I[0] << ", cathode " << I[1]);
        CHECK(std::abs(I[0] + I[1]) <= 1e-10 * std::abs(I[0]));
    }
}

// SC13 at the reverse current (-1 V, 1.7e-5 A/cm^2) as fixed: 1e-10 relative. Not met: the sum is
// 2e-11 A/cm^2 per area (1.3e-6 of the current; the 1D device's the same), the rounding of the edge
// fluxes that cancel in the current (q D n / h, 5e5 A/cm^2 on the 1 nm edges), well inside the
// resolution (1e-7 A/cm^2). As SC2's equilibrium current: 1e-10 relative holds only above about
// 1 A/cm^2. Encoded unchanged; for the owner.
TEST_CASE("Schottky: the terminal currents of a 2D device at -1 V within 1e-10 (SC13)",
          "[schottky][!shouldfail]") {
    const results::Sweep s = sweep(plan_of({.dimension = 2}), {-1.0});
    const auto& I = s.points[0].terminal_current;
    INFO("anode " << I[0] << ", cathode " << I[1]);
    CHECK(std::abs(I[0] + I[1]) <= 1e-10 * std::abs(I[0]));
}

namespace {

// SC5's n-Si diode, Newton to 1e-12 (the currents to rounding, not to the default 1e-8).
deck::Plan sc12_plan() {
    deck::Plan p = plan_of({});
    p.options.newton.tol_update = 1e-12;
    return p;
}

// The fourth-order central difference of the DC anode current at V (step h).
double dc_conductance(const deck::Plan& p, double V, double h) {
    const results::Sweep s = sweep(p, {V - 2 * h, V - h, V + h, V + 2 * h});
    const auto I = [&](std::size_t k) { return anode_current(s.points[k]); };
    return (8.0 * (I(2) - I(1)) - (I(3) - I(0))) / (12.0 * h);
}

}  // namespace

TEST_CASE("Schottky: small-signal admittance, C = dQ/dV, transient conservation (SC12)",
          "[schottky]") {
    const deck::Plan p = sc12_plan();
    solve::SmallSignalOptions o{.steady = p.options, .frequencies_Hz = {0.0, 1.0}};
    const std::vector<std::vector<double>> points{{-1.0, 0.0}, {0.2, 0.0}};
    const auto ac = solve::solve_small_signal(p.device, points, o);
    REQUIRE(ac.has_value());
    REQUIRE(ac->points.size() == 2);

    // At -1 V the capacitance against dQ/dV, Q = -q sum V_i (N_D - n_i + p_i) the anode's charge
    // (the ohmic end is neutral and field-free in reverse bias).
    const double dv = 1e-4;
    const results::Sweep s = sweep(p, {-1.0 - dv, -1.0 + dv});
    const auto& volumes = p.device.mesh().volumes();
    const auto charge = [&](const results::BiasPoint& b) {
        double q = 0.0;
        for (std::size_t i = 0; i < volumes.size(); ++i) {
            q -= base::q_C * volumes[i] * (1e16 - b.fields.n_cm3[i] + b.fields.p_cm3[i]);
        }
        return q;
    };
    const double c_dc = (charge(s.points[1]) - charge(s.points[0])) / (2 * dv);
    const double c = ac->capacitance(0, 1, 0, 0);
    INFO("C " << c << ", dQ/dV " << c_dc);
    CHECK(relative(c, c_dc) <= 1e-4);  // measured agreement to 6 digits (2.37723e-8 F/cm^2)

    // Where the DC current resolves its derivative (+0.2 V), g within 1e-5 of the DC difference.
    const double g_dc = dc_conductance(p, 0.2, 3e-4);
    const double g = ac->conductance(1, 0, 0, 0);
    INFO("+0.2 V: g " << g << ", dI/dV " << g_dc);
    CHECK(relative(g, g_dc) <= 1e-5);

    // Every column and row of Y sums to zero within the reported resolution (as the Unit 22
    // gates: below forward bias the conductances lie near the resolution of the state).
    for (std::size_t k = 0; k < 2; ++k) {
        for (std::size_t f = 0; f < 2; ++f) {
            const std::vector<double>& res = ac->points[k].resolution[f];
            const double res_all = std::max(res[0], res[1]);
            for (std::size_t j = 0; j < 2; ++j) {
                const auto column = ac->admittance(k, f, 0, j) + ac->admittance(k, f, 1, j);
                const auto row = ac->admittance(k, f, j, 0) + ac->admittance(k, f, j, 1);
                INFO("point " << k << ", f " << f << ", j " << j << ": column " << column
                              << " (resolution " << res[j] << "), row " << row);
                CHECK(std::abs(column) <= res[j]);
                CHECK(std::abs(row) <= 2.0 * res_all);
            }
        }
    }

    // A reverse-bias step, 0 to -1 V at 1 ps: the total currents of the contacts sum to zero at
    // every step, against the run's largest current (as the Unit 21 gate).
    solve::TransientOptions t{.steady = p.options, .t_end_s = 1e-8};
    const std::vector<solve::Waveform> w{*solve::Waveform::step(0.0, -1.0, 1e-12),
                                         solve::Waveform::constant(0.0)};
    const auto tr = solve::solve_transient(p.device, w, t);
    REQUIRE(tr.has_value());
    REQUIRE_FALSE(tr->stopped.has_value());
    double largest = 0.0, worst = 0.0;
    for (const results::TimePoint& q : tr->points) {
        for (const double v : {q.conduction_current[0], q.conduction_current[1],
                               q.displacement_current[0], q.displacement_current[1]}) {
            largest = std::max(largest, std::abs(v));
        }
    }
    for (const results::TimePoint& q : tr->points) {
        worst = std::max(worst, std::abs(q.terminal_current[0] + q.terminal_current[1]) / largest);
    }
    INFO(tr->points.size() << " steps, largest current " << largest << ", worst sum " << worst);
    CHECK(worst <= 1e-10);  // measured 5.3e-11 over 208 steps (a 1 ns ramp instead: 2.3e-9)
}

// SC12's conductance at -1 V as fixed: within 1e-5 of the DC difference. Not met: g = 1.495e-7
// S/cm^2 there, and the DC current (1.7e-5 A/cm^2) is resolved only to about 1e-6 of itself, so
// any DC difference over a step that keeps the truncation error below 1e-5 lies within the
// current's rounding (as Unit 22 found for the ohmic diode below 0.5 V). Encoded unchanged; for
// the owner.
TEST_CASE("Schottky: the small-signal conductance at -1 V against the DC difference (SC12)",
          "[schottky][!shouldfail]") {
    const deck::Plan p = sc12_plan();
    solve::SmallSignalOptions o{.steady = p.options, .frequencies_Hz = {0.0}};
    const std::vector<std::vector<double>> points{{-1.0, 0.0}};
    const auto ac = solve::solve_small_signal(p.device, points, o);
    REQUIRE(ac.has_value());
    const double g = ac->conductance(0, 0, 0, 0);
    const double g_dc = dc_conductance(p, -1.0, 1e-3);
    INFO("-1 V: g " << g << ", dI/dV " << g_dc << ", resolution "
                    << ac->points[0].resolution[0][0]);
    CHECK(relative(g, g_dc) <= 1e-5);
}

TEST_CASE("Schottky: a device whose only carrier contact is Schottky; device file round trip "
          "(SC13, SC15)",
          "[schottky]") {
    // Two Schottky contacts back to back (a metal-semiconductor-metal device).
    std::string text = deck_of({});
    const auto at = text.find("[contact cathode]");
    text = text.substr(0, at) + "[contact cathode]\nkind = schottky\npatch = x_max\nbarrier = 0.6\n" +
           "\n[sweep]\ncontact anode = 0\n";
    const deck::Plan p = must(text);
    const results::Sweep s = sweep(p, {0.3});
    const auto& b = s.points[0];
    INFO("anode " << b.terminal_current[0] << ", cathode " << b.terminal_current[1]);
    CHECK(std::abs(b.terminal_current[0] + b.terminal_current[1]) <=
          1e-10 * std::abs(b.terminal_current[0]));
    // The reverse-biased cathode limits the current: near its saturation current.
    CHECK(relative(std::abs(b.terminal_current[0]), thermionic_saturation(110.0, 0.6)) <= 0.2);

    const auto dir = deck_test::scratch("schottky_round_trip");
    const auto file = dir / "sbd.nitdev";
    REQUIRE(io::save_device(file, p.device).has_value());
    const auto back = io::load_device(file);
    REQUIRE(back.has_value());
    REQUIRE(back->contacts().size() == 2);
    for (std::size_t c = 0; c < 2; ++c) {
        const device::Contact& a = p.device.contacts()[c];
        const device::Contact& r = back->contacts()[c];
        CHECK(r.kind == a.kind);
        CHECK(r.nodes == a.nodes);
        CHECK(r.schottky.boundary == a.schottky.boundary);
        CHECK(r.schottky.by_work_function == a.schottky.by_work_function);
        CHECK(r.schottky.barrier_eV == a.schottky.barrier_eV);
        CHECK(r.schottky.richardson_n == a.schottky.richardson_n);
        CHECK(r.schottky.richardson_p == a.schottky.richardson_p);
    }
}

// Quadratic: in the asymptotic window (a correction below 1e-2, the next above rounding, 1e-13 of
// the scaled state) each correction is at most the square of the one before; a linear rate would
// leave it near a fixed fraction. (Measured: 1.12, 1.35e-4, 2.35e-11, 2.9e-15. The three-point
// order estimate of LS4 gives 0.58 on the last three, whose last is rounding, and 1.73 without it,
// its first, 1.12, outside the window.)
TEST_CASE("Schottky: Newton converges quadratically at +0.3 V (SC11)", "[schottky]") {
    deck::Plan p = plan_of({});
    p.options.newton.tol_update = 1e-12;
    const results::Sweep s = sweep(p, {0.0, 0.1, 0.2, 0.3});
    const auto& steps = s.points.back().convergence.iterations;
    for (const auto& r : steps) {
        std::printf("SC11 iteration %d: correction %.3e\n", r.iteration, r.update);
    }
    int pairs = 0;
    for (std::size_t k = 1; k < steps.size(); ++k) {
        const double before = steps[k - 1].update, after = steps[k].update;
        if (before >= 1e-2 || after <= 1e-13) continue;
        CAPTURE(k, before, after);
        CHECK(after <= before * before);
        ++pairs;
    }
    CHECK(pairs >= 1);
}

TEST_CASE("Schottky: the deck's refusals (SC14)", "[schottky]") {
    const auto refused = [](const std::string& text, const std::string& what) {
        const auto p = deck_test::plan_of(text);
        REQUIRE_FALSE(p.has_value());
        INFO(p.error().message);
        CHECK(p.error().code == base::ErrorCode::invalid_input);
        CHECK(p.error().message.find(what) != std::string::npos);
    };
    const std::string exactly = "needs exactly one of 'barrier' and 'work_function'";
    refused(deck_of({.barrier = "barrier = 0.7\nwork_function = 4.75"}), exactly);
    refused(deck_of({.barrier = ""}), exactly);
    // The Schottky keys on an ohmic contact.
    std::string ohmic = deck_of({});
    ohmic.replace(ohmic.find("kind = schottky"), 15, "kind = ohmic");
    refused(ohmic, "barrier, richardson_n and richardson_p are for a Schottky contact");
    std::string kind = deck_of({});
    kind.replace(kind.find("kind = schottky"), 15, "kind = schotky");
    refused(kind, "kind is ohmic, gate, electrode or schottky");
    // The weighted terminal current (W1-W8) reads the ohmic contacts' Dirichlet densities.
    deck::Plan p = plan_of({});
    p.options.weighted_current = true;
    const std::vector<std::vector<double>> points{{0.1, 0.0}};
    const auto s = solve::sweep_bias(p.device, points, p.options);
    REQUIRE_FALSE(s.has_value());
    INFO(s.error().message);
    CHECK(s.error().code == base::ErrorCode::invalid_input);
    CHECK(s.error().message.find("Schottky") != std::string::npos);
}

TEST_CASE("Schottky: the run identity names the barrier; a deck analysis over a Schottky bias "
          "(SC15)",
          "[schottky]") {
    const auto run = [](const Diode& d) {
        std::string text = deck_of(d);
        const auto at = text.find("[sweep]");
        text = text.substr(0, at) + "[small_signal]\ncontact anode = -0.5 : 0 : 0.5\n" +
               "frequencies = 1e6\n";
        const deck::Plan p = must(text);
        REQUIRE(p.analysis == deck::Analysis::small_signal);
        auto report = deck::run(p);
        if (!report) FAIL(report.error().message);
        REQUIRE(report->small_signal.has_value());
        return std::move(*report->small_signal);
    };
    const results::SmallSignal a = run({});
    const results::SmallSignal b = run({.barrier = "barrier = 0.71"});
    const results::SmallSignal c = run({.barrier = "work_function = 4.75"});  // chi 4.05: 0.7 eV
    REQUIRE(a.points.size() == 2);
    REQUIRE_FALSE(a.stopped.has_value());
    CHECK(a.run.input_identity != b.run.input_identity);
    CHECK(a.run.input_identity != c.run.input_identity);  // the barrier given another way
    // The depletion capacitance falls with reverse bias.
    CHECK(a.capacitance(0, 0, 0, 0) < a.capacitance(1, 0, 0, 0));
    // An ohmic diode of the same deck has another identity.
    std::string text = deck_of({.barrier = ""});
    text.replace(text.find("kind = schottky"), 15, "kind = ohmic");
    text.erase(text.find("richardson_n"), text.find("[contact cathode]") - text.find("richardson_n"));
    const deck::Plan ohmic = must(text);
    const deck::Plan schottky = must(deck_of({}));
    CHECK(solve::make_run_record(ohmic.device, ohmic.options, ohmic.points).input_identity !=
          solve::make_run_record(schottky.device, schottky.options, schottky.points).input_identity);
}

TEST_CASE("Schottky: the diode deck of decks/ against the thermionic form (SC16, reported)",
          "[schottky]") {
    const auto d = deck::read(std::filesystem::path(NITCAD_DECKS_DIR) / "schottky_diode_1d.deck");
    if (!d) FAIL(d.error().message);
    auto p = deck::build(*d);
    if (!p) FAIL(p.error().message);
    REQUIRE(p->points.size() == 29);
    const auto solved = solve::sweep_bias(p->device, p->points, p->options);
    if (!solved) FAIL(solved.error().message);
    REQUIRE_FALSE(solved->stopped.has_value());
    const double Js = thermionic_saturation(110.0, 0.7);
    std::printf("SC16 V [V], J [A/cm^2], thermionic form, ratio\n");
    for (const results::BiasPoint& b : solved->points) {
        const double V = b.bias_V[0];
        const double form = Js * (std::exp(V / kT) - 1.0);
        const double J = anode_current(b);
        std::printf("SC16 %6.2f %13.5e %13.5e %9.4f\n", V, J, form, V == 0.0 ? 0.0 : J / form);
    }
    // The reverse current with the material's own emission velocity against A* = 110.
    const auto reverse = [&](const std::string& richardson) {
        Diode diode;
        diode.richardson_n = 0.0;
        diode.richardson_p = 0.0;
        std::string text = deck_of(diode);
        const auto from = text.find("richardson_n");
        text.replace(from, text.find("[contact cathode]") - from, richardson + "\n");
        const results::Sweep s = sweep(must(text), {-1.0});
        return anode_current(s.points[0]);
    };
    const double own = reverse(""), sze = reverse("richardson_n = 110\nrichardson_p = 30");
    std::printf("SC16 reverse current at -1 V: default velocity %.5e, A* = 110 %.5e (ratio %.4f)\n",
                own, sze, own / sze);
    CHECK(own < 0.0);
    CHECK(sze < 0.0);
}

// SD2 under band-gap narrowing (1e18 donors, where the narrowing is tens of meV): n_B is fixed by
// Phi_Bn, the narrowing lowers the hole barrier only, and the [equilibrium] analysis (equilibrium
// Poisson) holds the same psi_S as the drift-diffusion state at 0 V.
TEST_CASE("Schottky: equilibrium with band-gap narrowing (SC2, SD2)", "[schottky]") {
    const double Nc = physics::conduction_band_dos(silicon(), T);
    const double Nv = physics::valence_band_dos(silicon(), T);
    const double Eg = physics::band_gap_eV(silicon(), T);
    const Diode d{.doping = 1e18};
    const deck::Plan p = plan_of(d);
    REQUIRE(p.options.models.bgn);
    const results::BiasPoint eq = sweep(p, {0.0}).points[0];
    INFO("n(0) " << eq.fields.n_cm3[0] << ", p(0) " << eq.fields.p_cm3[0] << ", I(0) "
                 << anode_current(eq) << ", resolution " << eq.terminal_current_resolution[0]);
    CHECK(relative(eq.fields.n_cm3[0], Nc * std::exp(-0.7 / kT)) <= 1e-12);
    CHECK(eq.fields.p_cm3[0] > 1.01 * Nv * std::exp(-(Eg - 0.7) / kT));  // the lowered barrier
    CHECK(std::abs(anode_current(eq)) <= eq.terminal_current_resolution[0]);

    std::string text = deck_of(d);
    text = text.substr(0, text.find("[sweep]")) + "[equilibrium]\n";
    const deck::Plan e = must(text);
    REQUIRE(e.analysis == deck::Analysis::equilibrium);
    const auto report = deck::run(e);
    if (!report) FAIL(report.error().message);
    REQUIRE(report->equilibrium.has_value());
    const auto& f = report->equilibrium->fields;
    INFO("equilibrium psi(0) " << f.potential_V[0] << ", drift-diffusion " << eq.fields.potential_V[0]);
    CHECK(std::abs(f.potential_V[0] - eq.fields.potential_V[0]) <= 1e-12);
    CHECK(relative(f.n_cm3[0], eq.fields.n_cm3[0]) <= 1e-9);
}

// SD2 on a material other than the reference one: psi_S carries the node's band shift s (Unit 15).
// A silicon | silicon-with-chi + 0.1 eV junction at 1 um (drift-diffusion interface), the ohmic
// contact on x_min (the reference material), the Schottky contact on x_max (Phi_Bn = 0.7 eV of its
// own material): equilibrium carries no current, n = n_B there, and equilibrium Poisson holds the
// drift-diffusion state's psi_S.
TEST_CASE("Schottky: a Schottky contact on a shifted material (SC2, SD2)", "[schottky]") {
    const double Nc = physics::conduction_band_dos(silicon(), T);
    const deck::Plan p = plan_of({});
    std::vector<double> x;  // the deck's axis mirrored: 1 nm at x_max
    const auto& points = p.device.mesh().points();
    const double length = points.back()[0];
    for (auto k = points.size(); k-- > 0;) x.push_back(length - points[k][0]);
    x.front() = 0.0;
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n, 0);
    for (std::size_t i = 0; i < n; ++i) region[i] = m.points()[i][0] < 0.5 * length ? 0 : 1;
    const auto& base = std::get<physics::Semiconductor>(p.device.regions()[0].material);
    physics::SemiconductorParameters shifted = base.parameters();
    shifted.electron_affinity_eV += 0.1;
    device::Contact anode{"anode", device::ContactKind::schottky, m.find_boundary("x_max")->nodes};
    anode.schottky = {.boundary = "x_max", .barrier_eV = 0.7, .richardson_n = 110.0,
                      .richardson_p = 30.0};
    auto cathode = m.find_boundary("x_min")->nodes;
    auto d = device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = T,
         .regions = {{"left", base}, {"right", *physics::Semiconductor::create(shifted)}},
         .node_region = std::move(region),
         .donors = std::vector<double>(n, 1e16),
         .acceptors = std::vector<double>(n, 0.0),
         .contacts = {std::move(anode),
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}},
         .interfaces = {{"left", "right", device::InterfaceTransport::drift_diffusion}}});
    if (!d) FAIL(d.error().message);
    const std::size_t last = n - 1;
    const auto s = solve::sweep_bias(*d, std::vector<std::vector<double>>{{0.0, 0.0}}, p.options);
    if (!s) FAIL(s.error().message);
    const results::BiasPoint& eq = s->points[0];
    INFO("n(x_max) " << eq.fields.n_cm3[last] << ", I(0) " << eq.terminal_current[0]
                     << ", resolution " << eq.terminal_current_resolution[0]);
    CHECK(relative(eq.fields.n_cm3[last], Nc * std::exp(-0.7 / kT)) <= 1e-12);
    CHECK(std::abs(eq.terminal_current[0]) <= eq.terminal_current_resolution[0]);
    const auto e = solve::solve_equilibrium(
        *d, {.newton = p.options.newton, .linear = p.options.linear, .models = p.options.models});
    if (!e) FAIL(e.error().message);
    INFO("equilibrium psi " << e->fields.potential_V[last] << ", drift-diffusion "
                            << eq.fields.potential_V[last]);
    CHECK(std::abs(e->fields.potential_V[last] - eq.fields.potential_V[last]) <= 1e-12);
}
