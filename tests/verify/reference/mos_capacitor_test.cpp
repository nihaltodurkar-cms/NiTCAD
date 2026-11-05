// L0 tests of the MOS-capacitor reference (verify/reference/mos_capacitor.hpp) and of the shared first integral
// (first_integral.hpp; DECISIONS Unit 55 14.1 (c), TVB-09).
//   - flat band: u_s = u_b, Q_s = 0, and C_LF = 1 / (1/C_ox + L_D / eps) (the extrinsic Debye length; the bulk's
//     n_b + p_b differs from N_A by 2 n_b, about 2e-14 relative here);
//   - the round trip: V_G recomputed from the returned u_s and Q_s, to 1e-25 relative;
//   - depletion: Q_s = -sqrt(2 q eps N_A V_T (t - 1 + e^-t)), t = psi_s / V_T, where the minority term is
//     negligible (to 1e-6);
//   - C_LF analytic against a central difference of the reference's own Q_G(V_G) in double-double (to 1e-9), and
//     against the complex step where that is evaluated (to 1e-10); C_LF above 0.9 C_ox in strong accumulation
//     (V_FB - 1 V) and strong inversion (V_FB + 2 V);
//   - V_FB on data where it is exact by construction (N_c = N_v, Q_f = 0, N = 2 n_i sinh 1, so u_b = -1);
//   - the first integral's profile: u_b everywhere at u_s = u_b; a negative distance refused; bad inputs refused.
#include "NiTCAD/verify/reference/mos_capacitor.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

#include "NiTCAD/verify/reference/first_integral.hpp"

namespace ref = NiTCAD::verify::reference;
using ref::DD;

namespace {

constexpr double q = 1.602176634e-19, eps_si = 11.7 * 8.8541878128e-14, VT = 0.025851999786435535, ni = 9.673e9;
constexpr double NA = 1e17, cox = 3.9 * 8.8541878128e-14 / 1e-6;

ref::MosInput input(double vfb = -0.35) {
    return {.q_C = q,
            .eps_F_per_cm = eps_si,
            .thermal_voltage_V = VT,
            .ni_cm3 = ni,
            .net_doping_cm3 = -NA,
            .oxide_capacitance_F_per_cm2 = cox,
            .flat_band_V = vfb};
}

}  // namespace

TEST_CASE("MOS reference: flat band") {
    const auto p = ref::mos_point(input(), -0.35);
    REQUIRE(p.has_value());
    CHECK(p->charge_C_per_cm2.to_double() == 0.0);
    CHECK(std::abs(p->band_bending_V.to_double()) <= 1e-28);
    const double LD = std::sqrt(eps_si * VT / (q * NA));
    const double cfb = 1.0 / (1.0 / cox + LD / eps_si);
    INFO("C_LF " << p->capacitance.to_double() << ", 1 / (1/C_ox + L_D/eps) " << cfb);
    CHECK(std::abs(p->capacitance.to_double() / cfb - 1.0) <= 1e-12);
}

TEST_CASE("MOS reference: the returned surface state reproduces V_G") {
    const ref::MosInput in = input();
    for (const double dv : {-1.0, -0.2, 0.05, 0.5, 1.5, 2.0}) {
        const auto p = ref::mos_point(in, in.flat_band_V + dv);
        REQUIRE(p.has_value());
        const DD back = DD(in.flat_band_V) + p->band_bending_V - p->charge_C_per_cm2 / DD(cox);
        const double rel = std::abs((back - p->gate_V).to_double()) / std::abs(p->gate_V.to_double());
        INFO("V_G - V_FB = " << dv << ": relative " << rel << ", Brent evaluations " << p->brent_evaluations);
        CHECK(rel <= 1e-25);
        CHECK((dv > 0.0 ? p->charge_C_per_cm2.hi < 0.0 : p->charge_C_per_cm2.hi > 0.0));
    }
}

TEST_CASE("MOS reference: depletion charge where the minority term is negligible") {
    const ref::MosInput in = input();
    for (const double dv : {0.2, 0.5, 0.8}) {
        const auto p = ref::mos_point(in, in.flat_band_V + dv);
        REQUIRE(p.has_value());
        const double psi = p->band_bending_V.to_double();
        const double tt = psi / VT;  // the majority term exactly: N_A V_T (t - 1 + e^-t)
        const double approx = -std::sqrt(2.0 * q * eps_si * NA * VT * (tt - 1.0 + std::exp(-tt)));
        const double rel = std::abs(p->charge_C_per_cm2.to_double() / approx - 1.0);
        INFO("V_G - V_FB = " << dv << ": psi_s " << psi << " V, Q_s " << p->charge_C_per_cm2.to_double()
                             << ", depletion form " << approx << ", relative " << rel);
        CHECK(rel <= 1e-6);
    }
}

TEST_CASE("MOS reference: C_LF against a difference of Q_G and the complex step; the C_ox limits") {
    const ref::MosInput in = input();
    for (const double dv : {-1.0, -0.2, 0.05, 0.5, 1.0, 1.5, 2.0}) {
        const double V = in.flat_band_V + dv, h = 1e-7;
        const auto p = ref::mos_point(in, V), a = ref::mos_point(in, V + h), b = ref::mos_point(in, V - h);
        REQUIRE(p.has_value());
        REQUIRE(a.has_value());
        REQUIRE(b.has_value());
        // Q_G = -Q_s - q Q_f, so C = -(Q_s(V + h) - Q_s(V - h)) / (2h), in double-double (error ~ h^2 C'' / 6).
        const double fd = (-(a->charge_C_per_cm2 - b->charge_C_per_cm2) / DD(2.0 * h)).to_double();
        const double c = p->capacitance.to_double();
        INFO("V_G - V_FB = " << dv << ": C_LF " << c << ", difference " << fd << ", complex step "
                             << p->capacitance_complex_step);
        CHECK(std::abs(fd / c - 1.0) <= 1e-9);
        if (p->complex_step_evaluated) CHECK(std::abs(p->capacitance_complex_step / c - 1.0) <= 1e-10);
        CHECK(c < cox);
        if (dv == -1.0 || dv == 2.0) CHECK(c > 0.9 * cox);  // strong accumulation and strong inversion
    }
}

TEST_CASE("MOS reference: V_FB where it is exact by construction") {
    const double N = 2.0 * ni * std::sinh(1.0);
    const auto v = ref::flat_band_voltage({.work_function_eV = 4.7,
                                           .electron_affinity_eV = 4.05,
                                           .band_gap_eV = 1.12,
                                           .conduction_dos_cm3 = 2e19,
                                           .valence_dos_cm3 = 2e19,
                                           .thermal_voltage_V = VT,
                                           .ni_cm3 = ni,
                                           .net_doping_cm3 = -N,
                                           .fixed_charge_cm2 = 0.0,
                                           .q_C = q,
                                           .oxide_capacitance_F_per_cm2 = cox});
    REQUIRE(v.has_value());
    const double expected = 4.7 - 4.05 - 0.56 - VT;
    INFO("V_FB " << v->to_double() << ", expected " << expected);
    CHECK(std::abs(v->to_double() - expected) <= 1e-15);
}

TEST_CASE("first integral: the flat profile, a negative distance, and bad MOS inputs") {
    const ref::NeutralBulk b = ref::neutral_bulk(ni, -NA);
    const std::vector<double> d{0.0, 1e-6, 1e-4};
    const auto flat = ref::bulk_profile(b, b.u_b, DD(1.3e-6), d);
    REQUIRE(flat.has_value());
    for (const auto& p : *flat) CHECK(p.u == b.u_b);
    const std::vector<double> bad{-1e-6};
    CHECK_FALSE(ref::bulk_profile(b, b.u_b + DD(1.0), DD(1.3e-6), bad).has_value());
    ref::MosInput in = input();
    in.oxide_capacitance_F_per_cm2 = 0.0;
    CHECK_FALSE(ref::mos_point(in, 0.0).has_value());
    in = input();
    in.net_doping_cm3 = 0.0;
    CHECK_FALSE(ref::mos_point(in, 0.0).has_value());
    CHECK_FALSE(ref::mos_point(input(), std::nan("")).has_value());
}
