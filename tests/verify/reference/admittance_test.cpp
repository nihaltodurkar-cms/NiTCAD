// L0 for the diode admittance reference (DECISIONS Unit 55 16.2 (l), 16.3, 16.18 (g), TVB-18):
//   - w = 0: sum (Y_fixed + Y_edge) is d J_diff / dV of diode.hpp (its W(V) edges), against a Richardson-extrapolated
//     double-double central difference, with SRH and without;
//   - the long-base limit G_0 sqrt(1 + i w tau) (W' = 40 L) and the short-base series 1/W' + kappa^2 W' / 3;
//   - the Riccati second method against the closed form over about 1e-3 <= w tau <= 1e4;
//   - the capacitive term i w C_j with qe_junction's C_j; the edge motion's size; the refusals;
//   - the depletion region's minority storage C_SCR (DT-Q34 A): positive, growing with forward bias, against its
//     second method (the fixed-region charge differenced in the bias), negligible at reverse bias.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <complex>
#include <vector>

#include "NiTCAD/verify/reference/admittance.hpp"
#include "NiTCAD/verify/reference/diode.hpp"
#include "NiTCAD/verify/reference/qe_junction.hpp"

using namespace NiTCAD::verify::reference;

namespace {

// Silicon-like inputs, typed here (the reference is base-free); diode_test's set.
DiodeInput input(bool srh, double length_cm) {
    DiodeInput in;
    in.junction = {.q_C = 1.602176634e-19, .eps_F_per_cm = 11.7 * 8.8541878128e-14,
                   .thermal_voltage_V = 0.025851999786435, .ni_cm3 = 9.65e9, .acceptors_cm3 = 1e16,
                   .donors_cm3 = 2e16};
    in.p_side = {.length_cm = length_cm, .minority_diffusivity_cm2_s = 32.0, .tau_n_s = 8e-6, .tau_p_s = 2.5e-6};
    in.n_side = {.length_cm = length_cm, .minority_diffusivity_cm2_s = 11.0, .tau_n_s = 7e-6, .tau_p_s = 2e-6};
    in.srh = srh;
    return in;
}

double rel(std::complex<double> a, std::complex<double> b) { return std::abs(a - b) / std::abs(b); }

// d J_diff / dV by central differences at h and h/2 in double-double, Richardson-extrapolated (O(h^4)).
double conductance(const DiodeInput& in, double V) {
    const auto d = [&](double h) {
        const auto a = diode_current(in, V + h), b = diode_current(in, V - h);
        REQUIRE(a.has_value());
        REQUIRE(b.has_value());
        return (a->diffusion_A_per_cm2 - b->diffusion_A_per_cm2) / DD(2.0 * h);
    };
    const DD d1 = d(1e-4), d2 = d(5e-5);
    return ((DD(4.0) * d2 - d1) / DD(3.0)).to_double();
}

}  // namespace

TEST_CASE("admittance: at w = 0 the edges' motion makes it d J_diff / dV", "[admittance]") {
    for (const bool srh : {true, false}) {
        for (const double L : {5e-4, 5e-2}) {
            const DiodeInput in = input(srh, L);
            // With SRH each diode_current runs the SNS quadrature (the Debug budget, CV-Q19): one bias there.
            const std::vector<double> biases = srh ? std::vector<double>{0.45} : std::vector<double>{0.3, 0.45, 0.6};
            for (const double V : biases) {
                const std::vector<double> w{0.0};
                const auto y = diode_admittance(in, V, w);
                REQUIRE(y.has_value());
                const auto& pt = y->points[0];
                const double G = conductance(in, V);
                CAPTURE(srh, L, V, G, pt.total.real(), pt.shockley.real());
                CHECK(std::abs(pt.total.imag()) == 0.0);
                CHECK(std::abs(pt.total.real() - G) <= 1e-9 * G);
                // Shockley's fixed-edge form misses it by the edge term, which is small but resolved.
                const double edge = (pt.p.edge_motion + pt.n.edge_motion).real();
                CHECK(std::abs(edge) > 1e-9 * G);
                CHECK(std::abs(edge) < 1e-2 * G);
            }
        }
    }
}

TEST_CASE("admittance: the long-base limit G_0 sqrt(1 + i w tau)", "[admittance]") {
    // W' = 40 L on both sides: coth(kappa W') = 1 to e^-80 at every frequency.
    DiodeInput in = input(true, 1.0);
    const std::vector<double> w0{0.0};
    const auto base = diode_admittance(in, 0.5, w0);
    REQUIRE(base.has_value());
    const double Lp = std::sqrt(32.0 * base->lifetime_p_s), Ln = std::sqrt(11.0 * base->lifetime_n_s);
    in.p_side.length_cm = 40.0 * Lp + base->dc.edge_p_cm.to_double();
    in.n_side.length_cm = 40.0 * Ln + base->dc.edge_n_cm.to_double();
    for (const double wt : {1e-3, 1.0, 1e3}) {
        const std::vector<double> w{wt / base->lifetime_p_s};
        const auto y = diode_admittance(in, 0.5, w);
        REQUIRE(y.has_value());
        const auto& pt = y->points[0];
        const auto& dc = y->dc;
        const double G0p = in.junction.q_C * 32.0 * dc.n_p0_cm3.to_double() * std::exp(0.5 / in.junction.thermal_voltage_V) /
                           in.junction.thermal_voltage_V / Lp;
        const std::complex<double> expect = G0p * std::sqrt(std::complex<double>(1.0, w[0] * y->lifetime_p_s));
        CAPTURE(wt, pt.p.fixed_edge, expect);
        CHECK(rel(pt.p.fixed_edge, expect) <= 1e-13);
    }
}

TEST_CASE("admittance: the short-base series at low frequency", "[admittance]") {
    const DiodeInput in = input(false, 5e-4);
    const std::vector<double> w{1e3, 1e4};
    const auto y = diode_admittance(in, 0.4, w);
    REQUIRE(y.has_value());
    const auto& dc = y->dc;
    const double VT = in.junction.thermal_voltage_V, q = in.junction.q_C, Wp = dc.base_p_cm.to_double();
    for (const auto& pt : y->points) {
        // kappa^2 = i w / D: kappa coth(kappa W') = 1/W' + kappa^2 W'/3 - kappa^4 W'^3 / 45 + ...
        const std::complex<double> k2(0.0, pt.omega_rad_s / 32.0);
        const std::complex<double> series = 1.0 / Wp + k2 * Wp / 3.0 - k2 * k2 * Wp * Wp * Wp / 45.0;
        const std::complex<double> expect = q * 32.0 * dc.n_p0_cm3.to_double() * std::exp(0.4 / VT) / VT * series;
        const double next = std::abs(k2 * k2 * k2) * std::pow(Wp, 6) * Wp * 2.0 / 945.0 * Wp;  // the next term's size
        CAPTURE(pt.omega_rad_s, pt.p.fixed_edge, expect, next);
        CHECK(rel(pt.p.fixed_edge, expect) <= 1e-14 + 10.0 * next);
    }
}

TEST_CASE("admittance: the Riccati second method agrees with the closed form", "[admittance]") {
    for (const double L : {5e-4, 5e-2}) {
        const DiodeInput in = input(true, L);
        const std::vector<double> w{1e2, 1e4, 1e5, 1e6, 1e8, 1e9};  // w tau_p about 8e-4 ... 8e3
        const auto y = diode_admittance(in, 0.5, w);
        REQUIRE(y.has_value());
        for (const auto& pt : y->points) {
            CAPTURE(L, pt.omega_rad_s, pt.p.fixed_edge, pt.p.fixed_edge_riccati);
            CHECK(rel(pt.p.fixed_edge_riccati, pt.p.fixed_edge) <= 1e-12);
            CHECK(rel(pt.n.fixed_edge_riccati, pt.n.fixed_edge) <= 1e-12);
            CHECK(pt.second_method_difference <= 1e-12 * std::abs(pt.shockley - pt.capacitive));
        }
    }
}

TEST_CASE("admittance: the edge motion at w > 0", "[admittance]") {
    // The moved edge is an effective junction voltage -u0' dW' V_T / (n_0 e^v) on Y_fixed, plus the curvature term:
    //     Y_edge = dW'/dV [q D U0 u0''/U0 - Y_fixed (u0'/U0) V_T (1 - e^-v)],
    // a second form of the header's expression, from Y_fixed (checked above) and the reported base data.
    for (const double L : {5e-4, 5e-2}) {
        const DiodeInput in = input(true, L);
        const double V = 0.45, VT = in.junction.thermal_voltage_V, q = in.junction.q_C;
        const std::vector<double> w{1e3, 1e5, 1e7};  // w tau_p about 8e-3, 0.8 and 80
        const auto y = diode_admittance(in, V, w);
        REQUIRE(y.has_value());
        const auto& dc = y->dc;
        const double Lp = std::sqrt(32.0 * y->lifetime_p_s), Wp = dc.base_p_cm.to_double();
        const double U0 = dc.n_p0_cm3.to_double() * std::expm1(V / VT);
        const double g = 1.0 / (Lp * std::tanh(Wp / Lp)), h = 1.0 / (Lp * Lp);
        for (const auto& pt : y->points) {
            const std::complex<double> expect =
                y->base_rate_p_cm_per_V * (q * 32.0 * U0 * h - pt.p.fixed_edge * g * VT * -std::expm1(-V / VT));
            CAPTURE(L, pt.omega_rad_s, pt.p.edge_motion, expect);
            CHECK(rel(pt.p.edge_motion, expect) <= 1e-12);
            CHECK(pt.p.edge_motion.imag() != 0.0);
        }
    }
}

TEST_CASE("admittance: the capacitive term and the totals", "[admittance]") {
    const DiodeInput in = input(true, 5e-4);
    const std::vector<double> w{1e5, 1e8};
    const auto y = diode_admittance(in, 0.45, w);
    REQUIRE(y.has_value());
    const auto qe = qe_junction(in.junction, 0.45);
    REQUIRE(qe.has_value());
    CHECK(y->capacitance_F_per_cm2.to_double() == qe->capacitance_F_per_cm2.to_double());
    CHECK(y->base_rate_p_cm_per_V == (qe->capacitance_F_per_cm2 / DD(in.junction.q_C * 1e16)).to_double());
    for (const auto& pt : y->points) {
        CHECK(pt.capacitive == std::complex<double>(0.0, pt.omega_rad_s * qe->capacitance_F_per_cm2.to_double()));
        const std::complex<double> s = pt.p.fixed_edge + pt.n.fixed_edge + pt.capacitive;
        CHECK(rel(pt.shockley, s) <= 1e-15);
        CHECK(pt.scr_storage == std::complex<double>(0.0, pt.omega_rad_s * y->scr_storage_F_per_cm2.to_double()));
        CHECK(rel(pt.total, s + pt.p.edge_motion + pt.n.edge_motion + pt.scr_storage) <= 1e-15);
        CHECK(pt.shockley.imag() > 0.0);  // capacitive and diffusive susceptance both positive
        CHECK(std::abs(1.0 / pt.p.effective_length_cm) > 0.0);
    }
}

TEST_CASE("admittance: the depletion region's minority storage C_SCR (DT-Q34 A)", "[admittance]") {
    // Without SRH (no SNS quadrature in diode_current; C_SCR does not depend on the lifetimes).
    const DiodeInput in = input(false, 5e-4);
    const std::vector<double> w{1e6};
    double last = 0.0;
    for (const double V : {0.35, 0.45}) {  // two biases: the Debug budget (CV-Q19)
        CAPTURE(V);
        const auto y = diode_admittance(in, V, w);
        REQUIRE(y.has_value());
        const double c = y->scr_storage_F_per_cm2.to_double(), cj = y->capacitance_F_per_cm2.to_double();
        CHECK(c > 0.0);
        CHECK(c > last);  // the stored minorities grow with forward bias
        last = c;
        CHECK(c < 0.05 * cj);  // a correction of a few percent at most here (16.37: 0.6-2.7 % for TVB-18)
        if (V == 0.45) {  // the second method: the fixed-region charge, differenced in the bias
            const auto s = scr_storage_second_method(in.junction, V);
            REQUIRE(s.has_value());
            CHECK(std::abs(*s - c) <= 1e-9 * c);
        }
    }
    // Deep reverse bias: no minorities to store.
    const auto r = diode_admittance(in, -2.0, w);
    REQUIRE(r.has_value());
    CHECK(r->scr_storage_F_per_cm2.to_double() < 1e-9 * r->capacitance_F_per_cm2.to_double());
    const auto s = scr_storage_second_method(in.junction, -2.0);
    REQUIRE(s.has_value());
    CHECK(std::abs(*s - r->scr_storage_F_per_cm2.to_double()) <= 1e-9 * r->capacitance_F_per_cm2.to_double());
}

TEST_CASE("admittance: refusals", "[admittance]") {
    const DiodeInput in = input(true, 5e-4);
    const std::vector<double> bad{-1.0}, nan{NAN};
    CHECK_FALSE(diode_admittance(in, 0.4, bad).has_value());
    CHECK_FALSE(diode_admittance(in, 0.4, nan).has_value());
    DiodeInput short_contact = in;
    short_contact.p_side.length_cm = 1e-7;
    const std::vector<double> w{1.0};
    CHECK_FALSE(diode_admittance(short_contact, 0.4, w).has_value());
}
