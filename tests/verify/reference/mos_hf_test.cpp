// L0 tests of the MOS high-frequency and deep-depletion reference (verify/reference/mos_hf.hpp; DECISIONS Unit 55
// 16.18 (b), 16.40, TVB-10).
//   - v = 0: Q_s(t, 0) is mos_point's Kingston-Neustadter charge (to 1e-28 relative);
//   - the first integral's consistency: Q_s(t, v) = the holes' charge + Q_n, both by quadrature of the charge
//     density over the profile, at v = 0 and v != 0, inversion and accumulation (to 1e-24 relative);
//   - depletion: C_HF = C_LF to 1e-6 (no inversion charge to hold); strong inversion: C_HF well below C_LF (which is
//     near C_ox), flat in V_G to 1 % (the HF minimum), and above the depletion-approximation bound eps / W_max;
//   - the second method (Q_s along Q_n held fixed by Brent, Richardson) runs in TVB-10, as a run-time self-check;
//   - deep depletion: the charge against its quadrature (to 1e-24 relative); equal to the equilibrium charge in
//     depletion (to 1e-8) and smaller in magnitude in strong inversion; the round trip V_G(t_s, Q) (1e-25);
//   - refusals: n-type substrate, flat band exactly.
#include "NiTCAD/verify/reference/mos_hf.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "NiTCAD/verify/reference/first_integral.hpp"
#include "NiTCAD/verify/reference/gauss_kronrod.hpp"

namespace ref = NiTCAD::verify::reference;
using ref::DD;

namespace {

constexpr double q = 1.602176634e-19, eps_si = 11.7 * 8.8541878128e-14, VT = 0.025851999786435535, ni = 9.673e9;
constexpr double NA = 1e17, cox = 3.9 * 8.8541878128e-14 / 1e-6, vfb = -0.9;

ref::MosInput input(double net = -NA) {
    return {.q_C = q,
            .eps_F_per_cm = eps_si,
            .thermal_voltage_V = VT,
            .ni_cm3 = ni,
            .net_doping_cm3 = net,
            .oxide_capacitance_F_per_cm2 = cox,
            .flat_band_V = vfb};
}

double rel(const DD& a, const DD& b) { return std::abs((a - b).to_double()) / std::abs(b.to_double()); }

}  // namespace

TEST_CASE("mos_hf: v = 0 is the Kingston-Neustadter charge", "[mos_hf]") {
    for (const double dv : {-1.0, -0.2, 0.3, 1.0, 2.0}) {
        const auto p = ref::mos_point(input(), vfb + dv);
        REQUIRE(p.has_value());
        const ref::NeutralBulk b = ref::mos_bulk(input());
        const auto c = ref::mos_redistribution_charges(input(), p->u_s - b.u_b, DD(0.0));
        REQUIRE(c.has_value());
        CAPTURE(dv);
        CHECK(rel(c->surface_C_per_cm2, p->charge_C_per_cm2) <= 1e-28);
    }
}

TEST_CASE("mos_hf: Q_s(t, v) is the holes' charge plus Q_n", "[mos_hf]") {
    const ref::NeutralBulk b = ref::mos_bulk(input());
    const DD scale = DD(2.0 * q) / (DD(eps_si) * DD(VT));
    for (const auto& [t, v] : {std::pair{30.0, 0.0}, std::pair{30.0, 0.7}, std::pair{40.0, -0.4}, std::pair{-4.0, 0.0},
                               std::pair{-4.0, 1.0}, std::pair{3.0, 0.2}}) {
        const DD tt(t), vv(v), nv = b.n_b * exp(-vv);
        const auto phi = [&](const DD& s) {
            return b.p_b * ref::poisson_boltzmann_G(-s) + nv * ref::poisson_boltzmann_G(s);
        };
        const auto holes = [&](const DD& s) { return DD(q) * b.p_b * expm1(-s) / sqrt(scale * phi(s)); };
        const DD lo = t < 0 ? tt : DD(0.0), hi = t < 0 ? DD(0.0) : tt;
        const auto Qp = ref::qag15(holes, lo, hi, 0.0, 1e-26);
        const auto c = ref::mos_redistribution_charges(input(), tt, vv);
        REQUIRE(c.has_value());
        CAPTURE(t, v, c->surface_C_per_cm2.to_double(), Qp.value.to_double(), c->minority_C_per_cm2.to_double());
        CHECK(rel(Qp.value + c->minority_C_per_cm2, c->surface_C_per_cm2) <= 1e-24);
        // The electrons' excess has the sign of the inversion (negative) or of their depletion (positive).
        CHECK((t > 0 ? c->minority_C_per_cm2.hi < 0.0 : c->minority_C_per_cm2.hi > 0.0));
    }
}

TEST_CASE("mos_hf: C_HF in depletion and strong inversion", "[mos_hf]") {
    // Depletion: no inversion charge to hold.
    const auto dep = ref::mos_hf_point(input(), vfb + 0.5);
    REQUIRE(dep.has_value());
    CHECK(rel(dep->capacitance_hf, dep->dc.capacitance) <= 1e-6);
    // The second method is not run here (CV-Q19's Debug budget, 16.40): TVB-10 runs it at each of its six gate
    // biases and stops with a run error unless it agrees to 1e-15 C (found 1e-17 .. 1e-28).
    for (const double dv : {1.2, 2.0}) {
        const auto p = ref::mos_hf_point(input(), vfb + dv);
        REQUIRE(p.has_value());
        CHECK_FALSE(p->second_method_evaluated);
        CHECK(p->quadrature_error <= 1e-20);
    }
    const auto s2 = ref::mos_hf_point(input(), vfb + 2.0);
    const auto s3 = ref::mos_hf_point(input(), vfb + 3.0);
    REQUIRE((s2 && s3));
    CHECK(s2->dc.capacitance.to_double() >= 0.9 * cox);  // LF: the inversion layer follows
    CHECK(s2->capacitance_hf.to_double() <= 0.5 * cox);
    CHECK(rel(s3->capacitance_hf, s2->capacitance_hf) <= 1e-2);  // the HF plateau
    // Above the depletion approximation's minimum: the band bending is at least 2 psi_B, so W <= W_max(2 psi_B + V_T).
    const ref::NeutralBulk b = ref::mos_bulk(input());
    const double psi = 2.0 * std::abs(b.u_b.to_double()) * VT;
    const double w_max = std::sqrt(2.0 * eps_si * (psi + VT) / (q * NA));
    const double c_dep = eps_si / w_max * cox / (eps_si / w_max + cox);
    CHECK(s2->capacitance_hf.to_double() >= 0.9 * c_dep);
    CHECK(s2->capacitance_hf.to_double() <= 1.1 * c_dep);
}

TEST_CASE("mos_hf: deep depletion", "[mos_hf]") {
    for (const double dv : {0.3, 1.0, 2.0, 4.0}) {
        const auto d = ref::mos_deep_depletion_point(input(), vfb + dv);
        REQUIRE(d.has_value());
        CAPTURE(dv, d->charge_C_per_cm2.to_double(), d->charge_second_method);
        CHECK(d->second_method_difference <= 1e-24 * std::abs(d->charge_C_per_cm2.to_double()));
        const DD back = DD(vfb) + DD(VT) * d->t_s - d->charge_C_per_cm2 / DD(cox);
        CHECK(rel(back, DD(vfb + dv)) <= 1e-25);
        const auto eq = ref::mos_point(input(), vfb + dv);
        REQUIRE(eq.has_value());
        if (dv < 0.5) {
            CHECK(rel(d->charge_C_per_cm2, eq->charge_C_per_cm2) <= 1e-8);
        } else if (dv >= 2.0) {
            CHECK(std::abs(d->charge_C_per_cm2.to_double()) < 0.9 * std::abs(eq->charge_C_per_cm2.to_double()));
            CHECK(d->t_s.to_double() > (eq->u_s - ref::mos_bulk(input()).u_b).to_double());
        }
    }
}

TEST_CASE("mos_hf: refusals", "[mos_hf]") {
    CHECK_FALSE(ref::mos_hf_point(input(NA), vfb + 1.0).has_value());
    CHECK_FALSE(ref::mos_deep_depletion_point(input(NA), vfb + 1.0).has_value());
    CHECK_FALSE(ref::mos_redistribution_charges(input(NA), DD(1.0), DD(0.0)).has_value());
}
