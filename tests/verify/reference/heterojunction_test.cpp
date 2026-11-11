// L0 tests of the heterojunction reference (verify/reference/heterojunction.hpp; DECISIONS Unit 55 17.3, TVB-08).
//   - one material on both sides: the pn_junction reference's V_bi, junction level and peak field, to 1e-26 relative;
//   - the matching condition eps_1 Phi_1 = eps_2 Phi_2, re-evaluated here with G(t) = expm1(t) - t (not the
//     library's series), to 1e-26 of eps Phi, on an isotype and an anisotype junction of two different materials;
//   - a pure affinity step between otherwise equal sides: V_bi = chi_1 - chi_2 and the band edges' step, to rounding;
//   - a small step (V_bi = 1e-5 V): the linearized (Debye-Hueckel) division phi_0 / V_bi = (eps_2 / L_2) / (eps_1 /
//     L_1 + eps_2 / L_2) and D = eps_1 phi_0 / L_1, to 1e-3 (the nonlinearity is of the order V_bi / V_T = 4e-4);
//   - the profile's ends: phi(0) is the interface level and the profile is monotone between the bulks;
//   - the collocation (I3), Release only ([.tvr5a_reference_release]): two degrees agree to 1e-10 V, the flux is
//     continuous to 1e-9, and the interface level, D and the profile agree with the first integral to within the
//     contact tail it reports;
//   - bad inputs are refused.
#include "NiTCAD/verify/reference/heterojunction.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "NiTCAD/verify/reference/first_integral.hpp"
#include "NiTCAD/verify/reference/pn_junction.hpp"

namespace ref = NiTCAD::verify::reference;
using ref::DD;

namespace {

constexpr double q = 1.602176634e-19, VT = 0.025851999786435535, eps0 = 8.8541878128e-14;

// Two III-V-like materials at 300 K (typical magnitudes; any positive data serve here).
ref::HeteroSide narrow(double C) { return {12.9 * eps0, 4.07, 1.424, 4.7e17, 9.0e18, C}; }
ref::HeteroSide wide(double C) { return {12.05 * eps0, 3.74, 1.798, 6.2e17, 1.1e19, C}; }

ref::HeterojunctionInput input(const ref::HeteroSide& l, const ref::HeteroSide& r) { return {q, VT, l, r}; }

DD G_expm1(const DD& t) { return ref::expm1(t) - t; }

double rel(const DD& a, const DD& b) { return std::abs(((a - b) / b).to_double()); }

}  // namespace

TEST_CASE("heterojunction: one material on both sides is the pn junction") {
    const auto j = ref::heterojunction(input(narrow(-1e16), narrow(2e16)));
    REQUIRE(j.has_value());
    REQUIRE(j->ni_left_cm3 == j->ni_right_cm3);
    const auto p = ref::pn_junction({q, 12.9 * eps0, VT, j->ni_left_cm3, 1e16, 2e16});
    REQUIRE(p.has_value());
    INFO("V_bi " << j->built_in_V.hi << " against " << p->built_in_V.hi);
    CHECK(rel(j->built_in_V, p->built_in_V) <= 1e-26);
    CHECK(rel(j->interface_V, p->psi_j - p->psi_p) <= 1e-26);
    CHECK(rel(j->field_left_V_per_cm, p->field_max_V_per_cm) <= 1e-26);
    CHECK(rel(j->field_right_V_per_cm, p->field_max_V_per_cm) <= 1e-26);
}

TEST_CASE("heterojunction: n_i and the anisotype V_bi against their formulas in plain double") {
    const ref::HeteroSide l = narrow(-1e17), r = wide(1e17);
    const auto j = ref::heterojunction(input(l, r));
    REQUIRE(j.has_value());
    const double ni1 = std::sqrt(l.Nc_cm3 * l.Nv_cm3) * std::exp(-l.band_gap_eV / (2.0 * VT));
    const double ni2 = std::sqrt(r.Nc_cm3 * r.Nv_cm3) * std::exp(-r.band_gap_eV / (2.0 * VT));
    CHECK(std::abs(j->ni_left_cm3 / ni1 - 1.0) <= 1e-13);
    CHECK(std::abs(j->ni_right_cm3 / ni2 - 1.0) <= 1e-13);
    // Majority-carrier limit (n_i << N): the p side's electrons are n_i1^2 / N_A, the N side's are N_D.
    const double vbi = (VT * std::log(1e17 / r.Nc_cm3) - r.affinity_eV) -
                       (VT * std::log(ni1 * ni1 / 1e17 / l.Nc_cm3) - l.affinity_eV);
    INFO("V_bi " << j->built_in_V.to_double() << " against " << vbi);
    CHECK(std::abs(j->built_in_V.to_double() / vbi - 1.0) <= 1e-12);
}

TEST_CASE("heterojunction: the interface level satisfies the matching condition") {
    for (const auto& in : {input(narrow(1e17), wide(1e17)), input(narrow(-1e17), wide(1e17)),
                           input(wide(1e17), narrow(3e16)), input(wide(-2e16), narrow(-5e17))}) {
        const auto j = ref::heterojunction(in);
        REQUIRE(j.has_value());
        const ref::NeutralBulk l = ref::neutral_bulk(j->ni_left_cm3, in.left.net_doping_cm3);
        const ref::NeutralBulk r = ref::neutral_bulk(j->ni_right_cm3, in.right.net_doping_cm3);
        const DD t1 = j->interface_V / DD(VT), t2 = (j->interface_V - j->built_in_V) / DD(VT);
        const DD a = DD(in.left.eps_F_per_cm) * (l.n_b * G_expm1(t1) + l.p_b * G_expm1(-t1));
        const DD b = DD(in.right.eps_F_per_cm) * (r.n_b * G_expm1(t2) + r.p_b * G_expm1(-t2));
        INFO("V_bi " << j->built_in_V.hi << ", phi_0 " << j->interface_V.hi << ", eps Phi " << a.hi);
        CHECK(std::abs((a - b).to_double()) <= 1e-26 * a.to_double());
        // The interface lies strictly between the bulks, and D follows from either side.
        CHECK(j->bend_left_V.hi * j->built_in_V.hi > 0.0);
        CHECK(j->bend_right_V.hi * j->built_in_V.hi > 0.0);
        const DD D2 = DD(2.0) * DD(q) * DD(VT) * b;
        CHECK(rel(j->displacement_C_per_cm2 * j->displacement_C_per_cm2, D2) <= 1e-25);
        CHECK(rel(j->field_left_V_per_cm * DD(in.left.eps_F_per_cm),
                  j->field_right_V_per_cm * DD(in.right.eps_F_per_cm)) <= 1e-28);
    }
}

TEST_CASE("heterojunction: a pure affinity step") {
    ref::HeteroSide a = narrow(1e17), b = narrow(1e17);
    b.affinity_eV = 4.07 - 0.25;
    const auto j = ref::heterojunction(input(a, b));
    REQUIRE(j.has_value());
    // V_bi = chi_1 - chi_2 (equal bulks); the conduction band steps up by the same amount at the interface.
    const DD step = DD(a.affinity_eV) - DD(b.affinity_eV);
    CHECK(rel(j->built_in_V, step) <= 1e-28);
    CHECK(rel(j->ec_right_eV - j->ec_left_eV, step) <= 1e-28);
    CHECK(rel(j->bend_left_V + j->bend_right_V, j->built_in_V) <= 1e-28);
    // V_bi > 0: the potential rises to the right, so the left side's bands bend down toward the interface (its
    // electrons accumulate) and E_c - E_F there lies below the bulk's -V_T ln(n_b / N_c) by the bending.
    const ref::NeutralBulk bulk = ref::neutral_bulk(j->ni_left_cm3, 1e17);
    const DD ec_bulk = -(DD(VT) * ref::log(bulk.n_b / DD(a.Nc_cm3)));
    CHECK(j->bend_left_V.to_double() > 0.0);
    CHECK(j->ec_left_eV.to_double() < ec_bulk.to_double());
    CHECK(rel(ec_bulk - j->ec_left_eV, j->bend_left_V) <= 1e-26);
    // No step at all: flat.
    const auto flat = ref::heterojunction(input(a, a));
    REQUIRE(flat.has_value());
    CHECK(flat->built_in_V.to_double() == 0.0);
    CHECK(flat->displacement_C_per_cm2.to_double() == 0.0);
}

TEST_CASE("heterojunction: a small step divides as the linearized problem") {
    ref::HeteroSide a = narrow(1e17), b = narrow(4e17);
    b.eps_F_per_cm = 9.0 * eps0;
    // Choose chi_2 so that V_bi = 1e-5 V: V_bi = V_T ln(n_b2 / n_b1) + chi_1 - chi_2 for equal N_c.
    const ref::NeutralBulk l = ref::neutral_bulk(1.0, 1e17), r = ref::neutral_bulk(1.0, 4e17);  // n_b ~ N here
    b.affinity_eV = a.affinity_eV + VT * std::log((r.n_b / l.n_b).to_double()) - 1e-5;
    const auto j = ref::heterojunction(input(a, b));
    REQUIRE(j.has_value());
    const double vbi = j->built_in_V.to_double();
    INFO("V_bi " << vbi);
    REQUIRE(std::abs(vbi - 1e-5) <= 1e-9);
    const double g1 = a.eps_F_per_cm / j->debye_left_cm.to_double(), g2 = b.eps_F_per_cm / j->debye_right_cm.to_double();
    const double share = j->interface_V.to_double() / vbi;
    CHECK(std::abs(share - g2 / (g1 + g2)) <= 1e-3 * share);
    const double D = j->displacement_C_per_cm2.to_double();
    CHECK(std::abs(D - g1 * j->interface_V.to_double()) <= 1e-3 * D);
}

TEST_CASE("heterojunction: the profile runs monotonically from bulk to bulk through the interface level") {
    const auto j = ref::heterojunction(input(narrow(-1e17), wide(1e17)));
    REQUIRE(j.has_value());
    const std::vector<double> x{-3e-5, -1e-5, -3e-6, -1e-6, -1e-8, 0.0, 1e-8, 1e-6, 3e-6, 1e-5, 3e-5};
    const auto p = ref::hetero_potential(*j, x);
    REQUIRE(p.has_value());
    CHECK(std::abs(((*p)[5].phi_V - j->interface_V).to_double()) <= 1e-24);
    for (std::size_t i = 1; i < x.size(); ++i) CHECK((*p)[i].phi_V.to_double() > (*p)[i - 1].phi_V.to_double());
    // 0.3 um from the junction (about 13 Debye lengths past the depletion edge) the bulks are reached to 1e-7 V,
    // and the reported tail is that distance.
    CHECK(std::abs((*p)[0].phi_V.to_double()) <= 1e-7);
    CHECK(std::abs(((*p)[10].phi_V - j->built_in_V).to_double()) <= 1e-7);
    CHECK(std::abs((*p)[0].phi_V.to_double()) <= 1.0000001 * VT * (*p)[0].tail.to_double());
    CHECK(std::abs((*p)[0].phi_V.to_double()) >= 0.9999999 * VT * (*p)[0].tail.to_double());
    CHECK(std::abs(((*p)[10].phi_V - j->built_in_V).to_double()) <= 1.0000001 * VT * (*p)[10].tail.to_double());
    // The one-sided slopes at the interface are the fields.
    const double left = ((*p)[5].phi_V - (*p)[4].phi_V).to_double() / 1e-8;
    const double right = ((*p)[6].phi_V - (*p)[5].phi_V).to_double() / 1e-8;
    CHECK(std::abs(left / j->field_left_V_per_cm.to_double() - 1.0) <= 1e-2);
    CHECK(std::abs(right / j->field_right_V_per_cm.to_double() - 1.0) <= 1e-2);
}

TEST_CASE("heterojunction: the collocation agrees with the first integral", "[.tvr5a_reference_release]") {
    for (const auto& in : {input(narrow(-1e17), wide(1e17)), input(narrow(1e17), wide(1e17))}) {
        const auto j = ref::heterojunction(in);
        REQUIRE(j.has_value());
        const double half = 4e-5;
        std::vector<double> x;
        for (int k = 0; k <= 80; ++k) x.push_back(-half + 2.0 * half * k / 80.0);
        const auto fi = ref::hetero_potential(*j, x);
        const auto c1 = ref::hetero_collocation(*j, -half, half, 8, 32, x);
        const auto c2 = ref::hetero_collocation(*j, -half, half, 8, 40, x);
        REQUIRE(fi.has_value());
        REQUIRE(c1.has_value());
        REQUIRE(c2.has_value());
        double degrees = 0.0, against = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i) {
            degrees = std::max(degrees, std::abs(c1->phi_V[i] - c2->phi_V[i]));
            against = std::max(against, std::abs(c2->phi_V[i] - (*fi)[i].phi_V.to_double()));
        }
        const double tail = VT * std::max((*fi)[0].tail.to_double(), (*fi)[x.size() - 1].tail.to_double());
        INFO("degrees differ by " << degrees << " V, collocation against first integral " << against
                                  << " V, contact tail " << tail << " V, estimate " << c2->tail_estimate_V
                                  << " V, Newton " << c2->newton_iterations);
        CHECK(degrees <= 1e-10);
        CHECK(c2->flux_mismatch <= 1e-9);
        CHECK(against <= tail + 1e-10);
        CHECK(std::abs(c2->interface_V - j->interface_V.to_double()) <= tail + 1e-10);
        CHECK(std::abs(c2->displacement_C_per_cm2 / j->displacement_C_per_cm2.to_double() - 1.0) <= 1e-8);
    }
}

TEST_CASE("heterojunction: bad inputs are refused") {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const auto bad = [](ref::HeterojunctionInput in) { return !ref::heterojunction(in).has_value(); };
    ref::HeterojunctionInput in = input(narrow(1e17), wide(1e17));
    CHECK_FALSE(bad(in));
    for (double ref::HeteroSide::* m : {&ref::HeteroSide::eps_F_per_cm, &ref::HeteroSide::band_gap_eV,
                                        &ref::HeteroSide::Nc_cm3, &ref::HeteroSide::Nv_cm3,
                                        &ref::HeteroSide::net_doping_cm3}) {
        ref::HeterojunctionInput b = in;
        b.right.*m = 0.0;
        CHECK(bad(b));
        b = in;
        b.left.*m = nan;
        CHECK(bad(b));
    }
    ref::HeterojunctionInput b = in;
    b.left.eps_F_per_cm = -1.0;
    CHECK(bad(b));
    b = in;
    b.thermal_voltage_V = 0.0;
    CHECK(bad(b));
    b = in;
    b.left.affinity_eV = nan;
    CHECK(bad(b));
    const auto j = ref::heterojunction(in);
    REQUIRE(j.has_value());
    const std::vector<double> x{0.0};
    CHECK_FALSE(ref::hetero_collocation(*j, 1e-5, 2e-5, 4, 16, x).has_value());
    CHECK_FALSE(ref::hetero_collocation(*j, -1e-5, 1e-5, 0, 16, x).has_value());
    CHECK_FALSE(ref::hetero_collocation(*j, -1e-5, 1e-5, 4, 1, x).has_value());
}
