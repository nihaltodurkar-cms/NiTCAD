// L0 for the ideal-diode reference (DECISIONS Unit 55 16.2 (d), 16.3, TVB-16):
//   - without SRH: the short-base law with the QE edges, recomputed here from qe_junction's charge (another code
//     path for the edges) and pn_junction's bulks;
//   - with SRH: the coth law's long- and short-base limits;
//   - SNS: Gauss-Kronrod against tanh-sinh in u, the change of variable against the x-space integral (pn_potential),
//     its sign, its exact 1 / tau scaling, and at deep reverse bias the generation bound
//     q n_i (x_p + x_n) / (tau_n + tau_p) and the depletion-approximation estimate below it;
//   - zero current at zero bias; the refusals.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <string>
#include <vector>

#include "NiTCAD/verify/reference/diode.hpp"
#include "NiTCAD/verify/reference/pn_junction.hpp"
#include "NiTCAD/verify/reference/qe_junction.hpp"
#include "NiTCAD/verify/reference/tanh_sinh.hpp"

using namespace NiTCAD::verify::reference;

namespace {

// Silicon-like inputs, typed here (the reference is base-free).
DiodeInput input(bool srh) {
    DiodeInput in;
    in.junction = {.q_C = 1.602176634e-19, .eps_F_per_cm = 11.7 * 8.8541878128e-14,
                   .thermal_voltage_V = 0.025851999786435, .ni_cm3 = 9.65e9, .acceptors_cm3 = 1e16,
                   .donors_cm3 = 2e16};
    in.p_side = {.length_cm = 5e-4, .minority_diffusivity_cm2_s = 32.0, .tau_n_s = 8e-6, .tau_p_s = 2.5e-6};
    in.n_side = {.length_cm = 5e-4, .minority_diffusivity_cm2_s = 11.0, .tau_n_s = 7e-6, .tau_p_s = 2e-6};
    in.srh = srh;
    return in;
}

double rel(double a, double b) { return std::abs(a - b) / std::abs(b); }

}  // namespace

TEST_CASE("diode: without SRH, the short-base law at the QE edges", "[diode]") {
    const DiodeInput in = input(false);
    const auto eq = pn_junction(in.junction);
    REQUIRE(eq.has_value());
    for (const double V : {0.2, 0.4, 0.6, -1.0}) {
        CAPTURE(V);
        const auto d = diode_current(in, V);
        REQUIRE(d.has_value());
        const auto qe = qe_junction(in.junction, V);
        REQUIRE(qe.has_value());
        const double Q = qe->charge_C_per_cm2.to_double(), q = in.junction.q_C;
        const double xp = Q / (q * in.junction.acceptors_cm3), xn = Q / (q * in.junction.donors_cm3);
        CHECK(rel(d->edge_p_cm.to_double(), xp) <= 1e-14);
        CHECK(rel(d->edge_n_cm.to_double(), xn) <= 1e-14);
        const double J = q * std::expm1(V / in.junction.thermal_voltage_V) *
                         (in.p_side.minority_diffusivity_cm2_s * eq->n_p.to_double() / (in.p_side.length_cm - xp) +
                          in.n_side.minority_diffusivity_cm2_s * eq->p_n.to_double() / (in.n_side.length_cm - xn));
        CHECK(d->n_p0_cm3 == eq->n_p);  // the reported bulk minorities
        CHECK(d->p_n0_cm3 == eq->p_n);
        CHECK(rel(d->diffusion_A_per_cm2.to_double(), J) <= 1e-13);
        CHECK(d->sns_A_per_cm2.hi == 0.0);
        CHECK(d->total_A_per_cm2 == d->diffusion_A_per_cm2);
    }
    const auto zero = diode_current(in, 0.0);
    REQUIRE(zero.has_value());
    CHECK(zero->total_A_per_cm2.hi == 0.0);
}

TEST_CASE("diode: with SRH, the long- and short-base limits of the coth law", "[diode]") {
    const double V = 0.4;
    DiodeInput in = input(true);
    const auto eq = pn_junction(in.junction);
    REQUIRE(eq.has_value());
    const double q = in.junction.q_C, f = std::expm1(V / in.junction.thermal_voltage_V);
    // Long base: both sides 50 diffusion lengths, coth = 1 to e^-100.
    {
        const auto probe = diode_current(in, V);
        REQUIRE(probe.has_value());
        const double Ln = probe->diffusion_length_p_cm.to_double(), Lp = probe->diffusion_length_n_cm.to_double();
        in.p_side.length_cm = 50.0 * Ln;
        in.n_side.length_cm = 50.0 * Lp;
        const auto d = diode_current(in, V);
        REQUIRE(d.has_value());
        const double J = q * f *
                         (in.p_side.minority_diffusivity_cm2_s * eq->n_p.to_double() / Ln +
                          in.n_side.minority_diffusivity_cm2_s * eq->p_n.to_double() / Lp);
        CHECK(rel(d->diffusion_A_per_cm2.to_double(), J) <= 1e-13);
        // L = sqrt(D tau) with tau the exact linearization 1/tau = (n0 + p0) / (tau_p (n0 + n_i) + tau_n (p0 + n_i)):
        // close to 1 / tau_n in the p bulk, and to 1 / tau_p in the n bulk, where the minority is holes.
        CHECK(rel(Ln, std::sqrt(in.p_side.minority_diffusivity_cm2_s * in.p_side.tau_n_s)) <= 1e-5);
        CHECK(rel(Lp, std::sqrt(in.n_side.minority_diffusivity_cm2_s * in.n_side.tau_p_s)) <= 1e-5);
    }
    // Short base: lifetimes so long that W' / L ~ 1e-6; L coth(W'/L) = W' (1 + (W'/L)^2 / 3 + ...).
    {
        DiodeInput s = input(true);
        for (DiodeSide* side : {&s.p_side, &s.n_side}) {
            side->tau_n_s *= 1e8;
            side->tau_p_s *= 1e8;
        }
        const auto d = diode_current(s, V);
        const auto none = diode_current(input(false), V);
        REQUIRE(d.has_value());
        REQUIRE(none.has_value());
        CHECK(rel(d->diffusion_A_per_cm2.to_double(), none->diffusion_A_per_cm2.to_double()) <= 1e-11);
    }
}

TEST_CASE("diode: the SNS current, two methods, sign, 1 / tau scaling, reverse generation bound", "[diode]") {
    const DiodeInput in = input(true);
    for (const double V : {-5.0, 0.5}) {  // two biases: each SRH call is about 0.1 s in Debug
        CAPTURE(V);
        const auto d = diode_current(in, V);
        REQUIRE(d.has_value());
        const double sns = d->sns_A_per_cm2.to_double();
        CAPTURE(sns, d->sns_second_method_difference, d->sns_quadrature_error);
        CHECK((V > 0.0 ? sns > 0.0 : sns < 0.0));
        CHECK(d->sns_second_method_difference <= 1e-12 * std::abs(sns));
        CHECK(d->sns_quadrature_error <= 1e-20 * std::abs(sns));
        CHECK(d->total_A_per_cm2 == d->diffusion_A_per_cm2 + d->sns_A_per_cm2);
        // Every lifetime doubled: U halves exactly in the depletion region (the diffusion term does not).
        DiodeInput slow = in;
        for (DiodeSide* side : {&slow.p_side, &slow.n_side}) {
            side->tau_n_s *= 2.0;
            side->tau_p_s *= 2.0;
        }
        const auto h = diode_current(slow, V);
        REQUIRE(h.has_value());
        CHECK(rel(2.0 * h->sns_A_per_cm2.to_double(), sns) <= 1e-20 + 1e-15);
    }
    // Deep reverse: |U| <= n_i / (tau_n + tau_p) (the denominator is at least n_i (tau_p + tau_n) and the excess
    // magnitude below n_i^2), and nearly equal where n, p << n_i. That excludes, near each edge, the zone where the
    // majority density is above n_i: within phi = V_T ln(N / n_i) of the bulk potential, which the depletion
    // approximation puts at sqrt(2 eps phi / (q N)) from the edge. So J_gen ~ q n_i sum W_eff / (tau_n + tau_p), to
    // about L_D / W_eff (the edges' own smearing).
    const auto d = diode_current(in, -5.0);
    REQUIRE(d.has_value());
    const double q = in.junction.q_C, ni = in.junction.ni_cm3;
    const double bound = q * ni *
                         (d->edge_p_cm.to_double() / (in.p_side.tau_n_s + in.p_side.tau_p_s) +
                          d->edge_n_cm.to_double() / (in.n_side.tau_n_s + in.n_side.tau_p_s));
    const double g = -d->sns_A_per_cm2.to_double();
    CAPTURE(g, bound);
    const double eps = in.junction.eps_F_per_cm, VT = in.junction.thermal_voltage_V;
    const auto excluded = [&](double N) { return std::sqrt(2.0 * eps * VT * std::log(N / ni) / (q * N)); };
    const double estimate =
        q * ni *
        ((d->edge_p_cm.to_double() - excluded(in.junction.acceptors_cm3)) / (in.p_side.tau_n_s + in.p_side.tau_p_s) +
         (d->edge_n_cm.to_double() - excluded(in.junction.donors_cm3)) / (in.n_side.tau_n_s + in.n_side.tau_p_s));
    CAPTURE(estimate, g / estimate);
    CHECK(g <= bound);
    CHECK(std::abs(g / estimate - 1.0) <= 0.1);
}

TEST_CASE("diode: refusals", "[diode]") {
    DiodeInput in = input(true);
    in.p_side.length_cm = 0.0;
    CHECK_FALSE(diode_current(in, 0.3).has_value());
    in = input(true);
    in.n_side.tau_p_s = -1.0;
    CHECK_FALSE(diode_current(in, 0.3).has_value());
    in.srh = false;  // lifetimes are not read without SRH
    CHECK(diode_current(in, 0.3).has_value());
    in = input(false);
    in.p_side.length_cm = 1e-5;  // 0.1 um: the p-side edge at -5 V (about 0.6 um) is past the contact
    const auto past = diode_current(in, -5.0);
    REQUIRE_FALSE(past.has_value());
    CHECK(past.error().find("W' <= 0") != std::string::npos);
    CHECK_FALSE(diode_current(input(false), -50.0).has_value());  // qe_junction's range
    CHECK_FALSE(diode_current(input(false), std::nan("")).has_value());
}

// The change of variable dx = du / |u'| (the first integral) against U integrated in x, with u(x) from pn_potential at
// each tanh-sinh abscissa (a separate profile code path). Each point is a Newton solve: 2 s in Release, 15 s in
// Debug, so it is hidden and run by the Release-only ctest verify.diode_profile (the Debug budget, CV-Q19).
TEST_CASE("diode: the SNS integral in u equals the integral of U(x) over the QE profile", "[.diode_profile]") {
    const DiodeInput in = input(true);
    const double V = -2.0;
    const auto d = diode_current(in, V);
    REQUIRE(d.has_value());
    PnJunctionInput shifted = in.junction;
    shifted.ni_cm3 = d->m_cm3.to_double();
    const auto jm = pn_junction(shifted);
    REQUIRE(jm.has_value());
    const double VT = in.junction.thermal_voltage_V, ni = in.junction.ni_cm3, m = shifted.ni_cm3;
    const double excess = ni * ni * std::expm1(V / VT);
    bool ok = true;
    const auto U_of = [&](const DiodeSide& s) {
        return [&, side = &s](const DD& x) {
            const std::vector<double> at{x.to_double()};
            const auto p = pn_potential(*jm, at);
            if (!p) {
                ok = false;
                return DD(0.0);
            }
            const double u = (*p)[0].psi_V.to_double() / VT;
            const double n = m * std::exp(u), pp = m * std::exp(-u);
            return DD(excess / (side->tau_p_s * (n + ni) + side->tau_n_s * (pp + ni)));
        };
    };
    // tanh-sinh wants an integral of order one: x = t x_edge on t in [0, 1], f divided by U(0).
    const DD xn = d->edge_n_cm, xp = d->edge_p_cm;
    const auto un = U_of(in.n_side), up = U_of(in.p_side);
    const DD u0n = abs(un(DD(0.0))), u0p = abs(up(DD(0.0)));
    const auto rn = tanh_sinh([&](const DD& t) { return un(t * xn) / u0n; }, DD(0.0), DD(1.0), 0.0, 1e-10, 6);
    const auto rp = tanh_sinh([&](const DD& t) { return up(-t * xp) / u0p; }, DD(0.0), DD(1.0), 0.0, 1e-10, 6);
    REQUIRE(ok);
    REQUIRE(rn.status == TanhSinhStatus::ok);
    REQUIRE(rp.status == TanhSinhStatus::ok);
    const double in_x = in.junction.q_C * (rn.value * u0n * xn + rp.value * u0p * xp).to_double();
    CAPTURE(in_x, d->sns_A_per_cm2.to_double(), rn.abs_error, rp.abs_error, static_cast<int>(rn.status), rn.levels,
            rn.evaluations, static_cast<int>(rp.status));
    CHECK(rel(in_x, d->sns_A_per_cm2.to_double()) <= 1e-8);
}
