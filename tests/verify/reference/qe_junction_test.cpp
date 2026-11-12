// L0 for the QE junction reference (DECISIONS Unit 55 16.2 (c), 16.3, TVB-17):
//   - at zero bias it is TVR-2's equilibrium junction (pn_junction), Q = eps E_max;
//   - at bias it is the equilibrium junction with n_i -> n_i e^(V / 2 V_T) (the shift identity), through
//     pn_junction's independent code path (asinh bulks, first_integral);
//   - C = -dQ/dV (dual numbers) against complex step and against a Richardson-extrapolated central difference of
//     Q(V) in double-double;
//   - the depletion form eps / W (with -2 V_T) approaches C at reverse bias, to rounding by -5 V.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "NiTCAD/verify/reference/pn_junction.hpp"
#include "NiTCAD/verify/reference/qe_junction.hpp"

using namespace NiTCAD::verify::reference;

namespace {

// Silicon-like inputs, typed here (the reference is base-free).
PnJunctionInput input() {
    return {.q_C = 1.602176634e-19, .eps_F_per_cm = 11.7 * 8.8541878128e-14, .thermal_voltage_V = 0.025851999786435,
            .ni_cm3 = 9.65e9, .acceptors_cm3 = 1e16, .donors_cm3 = 2e16};
}

double rel(const DD& a, const DD& b) { return std::abs(((a - b) / b).hi); }

}  // namespace

TEST_CASE("qe junction: zero bias is the equilibrium junction", "[qe_junction]") {
    const PnJunctionInput in = input();
    const auto q = qe_junction(in, 0.0);
    const auto eq = pn_junction(in);
    REQUIRE(q.has_value());
    REQUIRE(eq.has_value());
    CHECK(rel(q->charge_C_per_cm2, DD(in.eps_F_per_cm) * eq->field_max_V_per_cm) <= 1e-28);
    CHECK(rel(q->barrier_V, eq->built_in_V) <= 1e-28);
    CHECK(rel(q->u_j, eq->u_j) <= 1e-28);
}

TEST_CASE("qe junction: at bias, the equilibrium junction with n_i e^(V / 2 V_T)", "[qe_junction]") {
    const PnJunctionInput in = input();
    for (const double V : {-5.0, -1.0, -0.2, 0.3}) {
        CAPTURE(V);
        const auto q = qe_junction(in, V);
        REQUIRE(q.has_value());
        PnJunctionInput shifted = in;
        shifted.ni_cm3 = q->m_cm3.to_double();  // rounded to a double: 1e-16 relative in m
        const auto eq = pn_junction(shifted);
        REQUIRE(eq.has_value());
        CHECK(rel(q->charge_C_per_cm2, DD(in.eps_F_per_cm) * eq->field_max_V_per_cm) <= 1e-15);
        // The barrier V_T (u_n - u_p) at m (V_bi - V up to the bulks' own m^2 / N^2 dependence).
        CHECK(rel(q->barrier_V, eq->built_in_V) <= 1e-15);
    }
}

TEST_CASE("qe junction: C = -dQ/dV against complex step and a central difference", "[qe_junction]") {
    const PnJunctionInput in = input();
    for (const double V : {-5.0, -2.0, -0.5, 0.0, 0.2}) {
        CAPTURE(V);
        const auto q = qe_junction(in, V);
        REQUIRE(q.has_value());
        const double C = q->capacitance_F_per_cm2.to_double();
        REQUIRE(C > 0.0);
        const auto cs = qe_capacitance_complex_step(in, V);
        REQUIRE(cs.has_value());
        CHECK(std::abs(*cs - C) <= 1e-13 * C);
        // Central differences at d and d / 2, Richardson: O((d / V_T)^4) truncation, rounding about 1e-30 / d.
        const auto Qat = [&](double dv) { return qe_junction(in, V + dv)->charge_C_per_cm2; };
        const double d = 0x1p-20;  // V +- d and V +- d / 2 exact in double
        const DD c1 = (Qat(-d) - Qat(d)) / DD(2.0 * d);
        const DD c2 = (Qat(-0.5 * d) - Qat(0.5 * d)) / DD(d);
        const DD rich = (DD(4.0) * c2 - c1) / DD(3.0);
        CHECK(rel(q->capacitance_F_per_cm2, rich) <= 1e-14);
    }
}

TEST_CASE("qe junction: the depletion form approaches C at reverse bias", "[qe_junction]") {
    const PnJunctionInput in = input();
    double previous = 1.0;
    for (const double V : {-0.5, -2.0, -5.0, -10.0}) {
        CAPTURE(V);
        const auto q = qe_junction(in, V);
        REQUIRE(q.has_value());
        const double d = rel(q->depletion_form_F_per_cm2, q->capacitance_F_per_cm2);
        CAPTURE(d);
        CHECK(d <= std::max(previous, 1e-28));  // falls until it reaches rounding
        previous = d;
    }
    CHECK(previous <= 1e-28);  // exact up to e^(-barrier) and m^2 / N^2 terms (qe_junction.hpp)
    CHECK_FALSE(qe_junction(in, std::nan("")).has_value());
    CHECK_FALSE(qe_junction(in, -50.0).has_value());  // m^2 = n_i^2 e^(V / V_T) would underflow
    const auto far = qe_junction(in, -18.0);  // m about 1e-141: refused by the range guard itself
    REQUIRE_FALSE(far.has_value());
    CHECK(far.error().find("too far") != std::string::npos);
    PnJunctionInput bad = in;
    bad.donors_cm3 = 0.0;
    CHECK_FALSE(qe_junction(bad, -1.0).has_value());
}
