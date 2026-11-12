// L0 tests of the multiplication reference (verify/reference/multiplication.hpp; DECISIONS Unit 55 17.3, 17.13
// (b), IG-Q3 A, TVB-23).
//   - the field: the barrier and E_max against qe_junction.hpp inside its range (0 to -12 V; 1e16 against 2e17 and
//     1e19, and the mirror 1e19 against 1e16, whose n-side minority term carries the spill-over) to 1e-24; past it
//     (-20, -60, -200 V) against the depletion form with the -2 V_T correction, E_max^2 = 2 q N_eff (V_bi - V -
//     2 V_T) / eps, which differs only by terms exponentially small in the barrier (1e-9), and smooth across
//     qe_junction's limit;
//   - one carrier ionizing: with alpha_p = 0, I_n = 1 - e^-phi(W), M_n = e^phi(W) and M_p = 1; with alpha_n = 0,
//     I_n = 0, M_n = 1 and M_p = e^-phi(W) (phi < 0): closed forms the integrator does not know, to 1e-11;
//   - equal coefficients: phi = 0, I_n = I_p and M_n = M_p;
//   - silicon's two different sets: the identity 1 - I_n = e^-phi(W) - I_p to 1e-11, the tolerance error below
//     1e-10, and the nested Kronrod quadrature (the second method) within 5e-9 at 16 panels;
//   - a weak field in a junction whose two sides both bend by many V_T (1e16 against 2e17, where the field is linear on
//     each side): the integral against the leading Laplace term (W_p + W_n) alpha(E_max) E_max / (gamma B), within 30 %;
//   - a strongly one-sided junction (1e16 against 1e19) at zero bias: the heavy side's carriers spill across the
//     junction, the peak field (1.8e5 V/cm) is well above the depletion form's (5e4 V/cm), and I_n is 4.5e-5: pinned,
//     so that the spill-over layer stays in the reference;
//   - a peak field between the cut gamma B / 200 and ten times it (2e16 against 2e16 at 0 V): the integral against
//     the linear-field series to 2 % (the cut keeps such fields);
//   - the branch switch: swapped branch sets with the switch at 0 and above any field give the single-branch result
//     exactly, and a switch at half the peak a result strictly between the two single-branch ones;
//   - the field window (where E falls to a fraction of its peak) against the linear field's (1 - f) W where the
//     majority carriers are negligible, and its refusals;
//   - no ionization: A = 0 everywhere gives M = 1 exactly;
//   - breakdown of a one-sided 1e16 junction: I_n = 1 within 1 % of 55.3 V (the depletion-approximation value of
//     Unit 19's tests, a sanity bound), and the bias at a given 1 / M returns that 1 / M;
//   - bad inputs are refused.
// The cases that integrate are Release-only ([.tvr5a_reference_release]).
#include "NiTCAD/verify/reference/multiplication.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <cstdio>
#include <limits>

#include "NiTCAD/verify/reference/qe_junction.hpp"

namespace ref = NiTCAD::verify::reference;
using ref::DD;

namespace {

constexpr double q = 1.602176634e-19, VT = 0.025851999786435535, eps = 11.7 * 8.8541878128e-14, ni = 9.673e9;

// Silicon's van Overstraeten-de Man sets (the values of NiTCAD's silicon; any positive data serve here).
constexpr ref::IonizationBranches electrons{7.03e5, 1.231e6, 7.03e5, 1.231e6, 5.0e5};
constexpr ref::IonizationBranches holes{1.582e6, 2.036e6, 6.71e5, 1.693e6, 4.0e5};
constexpr ref::IonizationBranches none{};

ref::MultiplicationInput input(double NA, double ND, const ref::IonizationBranches& e = electrons,
                               const ref::IonizationBranches& h = holes) {
    return {q, eps, VT, ni, NA, ND, e, h, 1.0};
}

double rel(const DD& a, const DD& b) { return std::abs(((a - b) / b).to_double()); }

}  // namespace

TEST_CASE("multiplication: the field is qe_junction's inside its range") {
    // Two dopings: for 1e16 against 1e19 the spill-over term is large, and beyond about -8 V it is formed from its
    // logarithm (the branch that carries the reference past qe_junction's limit).
    // The mirror (1e19 against 1e16) carries the spill-over on the n side, through its minority term; at 0 V it is
    // largest.
    for (const auto& [NA, ND] : {std::array{1e16, 2e17}, std::array{1e16, 1e19}, std::array{1e19, 1e16}})
    for (const double V : {0.0, -1.0, -5.0, -9.0, -12.0}) {
        const auto m = ref::multiplication(input(NA, ND, none, none), V);
        const auto j = ref::qe_junction({q, eps, VT, ni, NA, ND}, V);
        REQUIRE(m.has_value());
        REQUIRE(j.has_value());
        INFO("V " << V << ": barrier " << m->barrier_V.hi << ", E_max " << m->field_max_V_per_cm.hi);
        CHECK(rel(m->barrier_V, j->barrier_V) <= 1e-24);
        CHECK(rel(m->field_max_V_per_cm, j->charge_C_per_cm2 / DD(eps)) <= 1e-24);
        CHECK(rel(m->bending_p_V + m->bending_n_V, m->barrier_V) <= 1e-28);
    }
}

TEST_CASE("multiplication: the field continues past qe_junction's limit") {
    const ref::MultiplicationInput in = input(1e16, 2e17, none, none);
    REQUIRE_FALSE(ref::qe_junction({q, eps, VT, ni, 1e16, 2e17}, -20.0).has_value());  // the range this extends
    const double Vbi = VT * std::log(1e16 * 2e17 / (ni * ni));
    const double Neff = 1.0 / (1.0 / 1e16 + 1.0 / 2e17);
    for (const double V : {-20.0, -60.0, -200.0}) {
        const auto m = ref::multiplication(in, V);
        REQUIRE(m.has_value());
        const double E2 = 2.0 * q * Neff * (Vbi - V - 2.0 * VT) / eps;
        const double E = m->field_max_V_per_cm.to_double();
        INFO("V " << V << ": E_max " << E << " against " << std::sqrt(E2));
        CHECK(std::abs(E * E / E2 - 1.0) <= 1e-9);
        CHECK(std::abs(m->barrier_V.to_double() - (Vbi - V)) <= 1e-9);
    }
    // Smooth across -15 V: the second difference of E_max^2 in V vanishes (E_max^2 is linear in V there).
    const auto e2 = [&](double V) {
        const double E = ref::multiplication(in, V)->field_max_V_per_cm.to_double();
        return E * E;
    };
    const double a = e2(-14.0), b = e2(-15.5), c = e2(-17.0);
    CHECK(std::abs(a - 2.0 * b + c) <= 1e-12 * b);
}

TEST_CASE("multiplication: one carrier ionizing has closed forms", "[.tvr5a_reference_release]") {
    {
        const auto m = ref::multiplication(input(1e16, 1e19, electrons, none), -40.0);
        REQUIRE(m.has_value());
        INFO("phi " << m->phi_total.hi << ", I_n " << m->ionization_integral.hi);
        CHECK(m->phi_total.hi > 0.1);
        CHECK(std::abs((m->ionization_integral - (DD(1.0) - ref::exp(-m->phi_total))).to_double()) <= 1e-11);
        CHECK(rel(m->multiplication_n, ref::exp(m->phi_total)) <= 1e-10);
        CHECK(std::abs(m->multiplication_p.to_double() - 1.0) <= 1e-10);
        CHECK(m->hole_integral.to_double() == 0.0);
    }
    {
        const auto m = ref::multiplication(input(1e16, 1e19, none, holes), -40.0);
        REQUIRE(m.has_value());
        CHECK(m->ionization_integral.to_double() == 0.0);
        CHECK(m->multiplication_n.to_double() == 1.0);
        CHECK(m->phi_total.hi < -0.01);
        CHECK(rel(m->multiplication_p, ref::exp(-m->phi_total)) <= 1e-10);
        // The identity gives I_p = e^-phi - 1.
        CHECK(std::abs((m->hole_integral - (ref::exp(-m->phi_total) - DD(1.0))).to_double()) <= 1e-11);
    }
}

TEST_CASE("multiplication: equal coefficients", "[.tvr5a_reference_release]") {
    const auto m = ref::multiplication(input(1e16, 1e19, electrons, electrons), -40.0);
    REQUIRE(m.has_value());
    CHECK(m->phi_total.to_double() == 0.0);
    CHECK((m->ionization_integral - m->hole_integral).to_double() == 0.0);
    CHECK(rel(m->multiplication_n, m->multiplication_p) <= 1e-28);
    CHECK(m->ionization_integral.hi > 0.05);
}

TEST_CASE("multiplication: the identity, the tolerance error and the second method", "[.tvr5a_reference_release]") {
    for (const auto& [NA, ND, V] : {std::array{1e16, 1e19, -30.0}, std::array{1e16, 1e19, -50.0},
                                    std::array{1e17, 2e17, -14.0}, std::array{1e19, 2e16, -30.0}}) {
        const auto m = ref::multiplication(input(NA, ND), V);
        const auto s = ref::multiplication_quadrature(input(NA, ND), V, 16);
        REQUIRE(m.has_value());
        REQUIRE(s.has_value());
        const double second = std::abs((s->ionization_integral - m->ionization_integral).to_double());
        INFO("N_A " << NA << ", N_D " << ND << ", V " << V << ": I_n " << m->ionization_integral.hi << ", phi "
                    << m->phi_total.hi << ", M_n " << m->multiplication_n.hi << ", M_p " << m->multiplication_p.hi
                    << ", tolerance error " << m->integration_error << ", identity " << m->identity_residual
                    << ", quadrature differs by " << second << ", steps " << m->steps);
        CHECK(m->below_breakdown);
        CHECK(m->ionization_integral.hi > 1e-3);
        CHECK(m->identity_residual <= 1e-11);
        CHECK(m->integration_error <= 1e-10);
        CHECK(second <= 5e-9);
        CHECK(std::abs((s->phi_total - m->phi_total).to_double()) <= 5e-9);
        CHECK(std::abs((s->hole_integral - m->hole_integral).to_double()) <= 5e-9);
        // Electrons ionize more than holes in silicon: phi > 0, so a hole injected from the n side multiplies less.
        CHECK(m->phi_total.hi > 0.0);
        CHECK(m->multiplication_p.to_double() < m->multiplication_n.to_double());
        CHECK(m->multiplication_p.to_double() > 1.0);
    }
}

TEST_CASE("multiplication: a weak field against the leading Laplace term", "[.tvr5a_reference_release]") {
    // 1e16 against 2e17, electrons only, at -5 V: each side bends by many V_T, so the field is linear on each side, E =
    // E_max (1 - x / W), and phi = int A e^(-B / E) dx = (W_p + W_n) A e^-b / b (1 - 2 / b + ...), b = B / E_max. The
    // bracket is between 0.7 and 1 for b near 9.
    const auto m = ref::multiplication(input(1e16, 2e17, electrons, none), -5.0);
    REQUIRE(m.has_value());
    const double E = m->field_max_V_per_cm.to_double(), W = (m->depletion_p_cm + m->depletion_n_cm).to_double();
    const double b = 1.231e6 / E;
    const double leading = W * 7.03e5 * std::exp(-b) / b;
    const double phi = m->phi_total.to_double();  // int alpha_n dx with alpha_p = 0
    INFO("E_max " << E << ", W " << W << ", b " << b << ", phi " << phi << ", leading term " << leading);
    CHECK(phi < leading);
    CHECK(phi > 0.70 * leading);
}

TEST_CASE("multiplication: a peak field just above the cut", "[.tvr5a_reference_release]") {
    // 2e16 against 2e16 at 0 V, electrons only: E_max is near 5e4 V/cm, between the cut gamma B / 200 and ten times
    // it, so the integral is that of the field nearest the peak. For a field linear on each side the integral is
    // (W_p + W_n) A e^-b / b (1 - 2 / b + 6 / b^2 - 24 / b^3 + ...), b = B / E_max near 25 (the next term is 5e-4).
    const auto m = ref::multiplication(input(2e16, 2e16, electrons, none), 0.0);
    REQUIRE(m.has_value());
    const double E = m->field_max_V_per_cm.to_double(), W = (m->depletion_p_cm + m->depletion_n_cm).to_double();
    const double b = 1.231e6 / E;
    const double series = 1.0 - 2.0 / b + 6.0 / (b * b) - 24.0 / (b * b * b);
    const double expected = W * 7.03e5 * std::exp(-b) / b * series;
    const double phi = m->phi_total.to_double();
    INFO("E_max " << E << ", W " << W << ", b " << b << ", phi " << phi << ", linear-field value " << expected);
    REQUIRE(E > 1.231e6 / 200.0);
    REQUIRE(E < 1.231e6 / 20.0);
    CHECK(std::abs(phi / expected - 1.0) <= 2e-2);
}

TEST_CASE("multiplication: each carrier's branch follows its switch field", "[.tvr5a_reference_release]") {
    // Electrons only, 1e16 against 1e19 at -40 V. A set whose high branch is silicon's electron branch and whose switch
    // is 0 uses only that branch; a set with that branch low and a switch above any field, likewise: both give the
    // single-branch result exactly (the same path, the same coefficients at each node). A switch at half the peak
    // field gives a result strictly between the two single-branch results.
    constexpr ref::IonizationBranches weak{1.0e6, 2.0e6, 1.0e6, 2.0e6, 5.0e5};  // smaller than silicon's at these fields
    const auto phi = [](const ref::IonizationBranches& e) {
        const auto m = ref::multiplication(input(1e16, 1e19, e, none), -40.0);
        REQUIRE(m.has_value());
        return m->phi_total;
    };
    const DD strong_only = phi(electrons), weak_only = phi(weak);
    REQUIRE(weak_only.hi < 0.5 * strong_only.hi);
    CHECK((phi({1.0e6, 2.0e6, 7.03e5, 1.231e6, 0.0}) - strong_only).to_double() == 0.0);
    CHECK((phi({7.03e5, 1.231e6, 1.0e6, 2.0e6, 1e300}) - strong_only).to_double() == 0.0);
    const double peak = ref::multiplication(input(1e16, 1e19, none, none), -40.0)->field_max_V_per_cm.to_double();
    const DD mixed = phi({7.03e5, 1.231e6, 1.0e6, 2.0e6, 0.5 * peak});  // low below half the peak, weak above
    INFO("strong " << strong_only.hi << ", weak " << weak_only.hi << ", mixed " << mixed.hi);
    CHECK(mixed.hi < 0.9 * strong_only.hi);
    CHECK(mixed.hi > weak_only.hi);
}

TEST_CASE("multiplication: the field window") {
    // Where the majority carriers are negligible the field is linear, E = E_max (1 - |x| / W), W = eps E_max / (q N):
    // the window's edge at E = f E_max is |x| = (1 - f) W. At the edge of a tenth of the peak the majority density is
    // N e^-(1 + z), z = f^2 eps E_max^2 / (2 q N V_T) (17 on both sides of 5e15 at -100 V), so the linear form holds
    // to about e^-z there. In 1e16 against 2e17 at -5 V the n side bends by only about 10 V_T: its electrons spill
    // into the p side near the junction at about N_D e^-10, 1e-3 of N_A, which bends the p side's field by that much
    // (within 3e-3 of linear), and its own majority tail makes the field there fall more slowly than linearly
    // (x_n between (1 - f) W_n and W_n + 3 L_D, L_D the Debye length at N_D).
    for (const double f : {0.1, 0.5, 0.9}) {
        const auto w = ref::multiplication_field_window(input(5e15, 5e15, none, none), -100.0, f);
        const auto m = ref::multiplication(input(5e15, 5e15, none, none), -100.0);
        REQUIRE(w.has_value());
        REQUIRE(m.has_value());
        INFO("f " << f << ": x_p " << w->x_p_cm.hi << ", x_n " << w->x_n_cm.hi << ", W " << m->depletion_p_cm.hi);
        CHECK(rel(w->field_max_V_per_cm, m->field_max_V_per_cm) == 0.0);
        CHECK(rel(-w->x_p_cm, DD(1.0 - f) * m->depletion_p_cm) <= 1e-6);
        CHECK(rel(w->x_n_cm, DD(1.0 - f) * m->depletion_n_cm) <= 1e-6);
        CHECK(w->abs_error <= 1e-20);
    }
    {
        const auto w = ref::multiplication_field_window(input(1e16, 2e17, none, none), -5.0, 0.1);
        const auto m = ref::multiplication(input(1e16, 2e17, none, none), -5.0);
        REQUIRE(w.has_value());
        REQUIRE(m.has_value());
        const double debye = std::sqrt(eps * VT / (q * 2e17));
        INFO("x_p " << w->x_p_cm.hi << ", x_n " << w->x_n_cm.hi << ", W_p " << m->depletion_p_cm.hi << ", W_n "
                    << m->depletion_n_cm.hi << ", L_D " << debye);
        CHECK(rel(-w->x_p_cm, DD(0.9) * m->depletion_p_cm) <= 3e-3);
        CHECK(w->x_n_cm.hi > 0.9 * m->depletion_n_cm.hi);
        CHECK(w->x_n_cm.hi < m->depletion_n_cm.hi + 3.0 * debye);
    }
    for (const double f : {0.0, 1.0, -0.5, std::numeric_limits<double>::quiet_NaN()}) {
        CHECK_FALSE(ref::multiplication_field_window(input(5e15, 5e15, none, none), -100.0, f).has_value());
    }
    CHECK_FALSE(ref::multiplication_field_window(input(5e15, 5e15, none, none), 0.5, 0.1).has_value());
}

TEST_CASE("multiplication: no ionization") {
    const auto m = ref::multiplication(input(1e16, 1e19, none, none), -40.0);
    REQUIRE(m.has_value());
    CHECK(m->multiplication_n.to_double() == 1.0);
    CHECK(m->multiplication_p.to_double() == 1.0);
    CHECK(m->integration_error == 0.0);
}

TEST_CASE("multiplication: breakdown and the bias at a given 1 / M", "[.tvr5a_reference_release]") {
    const ref::MultiplicationInput in = input(1e16, 1e19);
    const auto bd = ref::bias_at_inverse_multiplication(in, true, 0.0, -80.0, -30.0);
    REQUIRE(bd.has_value());
    INFO("V_BR " << -*bd << " V");
    CHECK(std::abs(-*bd / 55.3 - 1.0) <= 1e-2);
    // 1 / M_n and 1 / M_p vanish at the same bias (both carry the factor 1 - I_n).
    const auto bp = ref::bias_at_inverse_multiplication(in, false, 0.0, -80.0, -30.0);
    REQUIRE(bp.has_value());
    CHECK(std::abs(*bp - *bd) <= 1e-8);
    for (const bool electron : {true, false}) {
        for (const double inverse : {0.5, 0.1, 0.01}) {
            const auto V = ref::bias_at_inverse_multiplication(in, electron, inverse, *bd, -1.0);
            REQUIRE(V.has_value());
            const auto m = ref::multiplication(in, *V);
            REQUIRE(m.has_value());
            const double found = (electron ? m->inverse_n : m->inverse_p).to_double();
            INFO((electron ? "electron" : "hole") << ", 1 / M " << inverse << ": V " << *V << ", found " << found);
            CHECK(std::abs(found - inverse) <= 1e-9);
            CHECK(*V > *bd);
        }
    }
    CHECK_FALSE(ref::bias_at_inverse_multiplication(in, true, 0.5, -10.0, -1.0).has_value());  // not in the bracket
    CHECK_FALSE(ref::bias_at_inverse_multiplication(in, true, 0.5, -1.0, -10.0).has_value());
}

TEST_CASE("multiplication: bad inputs are refused") {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const auto bad = [](const ref::MultiplicationInput& in, double V = -5.0) {
        return !ref::multiplication(in, V).has_value();
    };
    CHECK_FALSE(bad(input(1e16, 1e19, none, none)));
    for (double ref::MultiplicationInput::* m :
         {&ref::MultiplicationInput::q_C, &ref::MultiplicationInput::eps_F_per_cm,
          &ref::MultiplicationInput::thermal_voltage_V, &ref::MultiplicationInput::ni_cm3,
          &ref::MultiplicationInput::acceptors_cm3, &ref::MultiplicationInput::donors_cm3,
          &ref::MultiplicationInput::gamma}) {
        ref::MultiplicationInput in = input(1e16, 1e19, none, none);
        in.*m = 0.0;
        CHECK(bad(in));
        in.*m = nan;
        CHECK(bad(in));
    }
    ref::MultiplicationInput in = input(1e16, 1e19);
    in.electron.B_low_V_per_cm = 0.0;
    CHECK(bad(in));
    in = input(1e16, 1e19);
    in.hole.A_high_per_cm = -1.0;
    CHECK(bad(in));
    CHECK(bad(input(1e16, 1e19, none, none), nan));
    CHECK(bad(input(1e16, 1e19, none, none), 0.1));  // forward bias: this reference is for reverse bias
    CHECK_FALSE(ref::multiplication_quadrature(input(1e16, 1e19), -5.0, 0).has_value());
}

TEST_CASE("multiplication: the spill-over layer of a one-sided junction", "[.tvr5a_reference_release]") {
    const auto m = ref::multiplication(input(1e16, 1e19, electrons, none), 0.0);
    const auto s = ref::multiplication_quadrature(input(1e16, 1e19, electrons, none), 0.0, 16);
    REQUIRE(m.has_value());
    REQUIRE(s.has_value());
    const double Vbi = VT * std::log(1e16 * 1e19 / (ni * ni));
    const double depletion = std::sqrt(2.0 * q * 1e16 * (Vbi - 2.0 * VT) / eps);  // the light side alone
    const double E = m->field_max_V_per_cm.to_double(), phi = m->phi_total.to_double();
    std::printf("one-sided, 0 V: E_max %.6g V/cm (depletion form %.6g), n-side bending %.4g V, phi %.6e (quadrature %.6e)\n",
                E, depletion, m->bending_n_V.to_double(), phi, s->phi_total.to_double());
    // The heavy side bends by only a few V_T, so its electrons reach the junction and their charge adds to the field.
    CHECK(m->bending_n_V.to_double() < 6.0 * VT);
    CHECK(E > 3.0 * depletion);
    CHECK(E < 4.0 * depletion);
    CHECK(phi > 4.4e-5);
    CHECK(phi < 4.6e-5);
    CHECK(std::abs(s->phi_total.to_double() / phi - 1.0) <= 1e-8);
}
