// L0 tests of the Pao-Sah reference (verify/reference/pao_sah.hpp; DECISIONS Unit 55 16.18 (d), (e), 16.46; TVB-20,
// TVB-21):
//   - at V = 0 both forms reduce to the equilibrium MOS-C, so psi_s against TVB-09's Kingston-Neustadter reference
//     (mos_capacitor, separate code), 1e-24 relative;
//   - the bulk ambiguity in strong inversion: forms 1 and 2 differ by less than 1e-10 of Q (16.17: about
//     n_i^2 / N_A^2), the cutoff psi_B against psi_B / 2 by less than 1e-8 (the electrons below psi_B, about
//     n_i V_T / E(psi_B) per area, against the channel's 1e12 cm^-2: some 2e-9);
//   - the second method: dQ/dV_D against mu (-Q_I(V_D)), 1e-12;
//   - limits: the linear region Q -> mu (-Q_I(0)) V_D; the subthreshold slope against ln 10 V_T (1 + C_dep / C_ox)
//     (5 %, the depletion approximation); the charge sheet within 5 % in strong inversion, its closed form against
//     Gauss-Kronrod (1e-24);
//   - refusals.
// The ambiguity and the identity cases (1.5 s in Release) are hidden and run in the Release-only ctest
// verify.pao_sah_release (CV-Q19's Debug budget); TVB-20 also runs the identity at run time.
#include "NiTCAD/verify/reference/pao_sah.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <utility>

#include "NiTCAD/verify/reference/mos_capacitor.hpp"

namespace ref = NiTCAD::verify::reference;
using ref::DD;

namespace {

constexpr double q = 1.602176634e-19, eps0 = 8.8541878128e-14, VT = 0.025851997, ni = 1e10, NA = 1e17;
constexpr double eps = 11.7 * eps0, cox = 3.9 * eps0 / 1e-6, VFB = -0.9, mu = 1400.0;

ref::PaoSahInput input(int form = 1, double cutoff = 1.0) {
    return {.q_C = q,
            .eps_F_per_cm = eps,
            .thermal_voltage_V = VT,
            .ni_cm3 = ni,
            .acceptors_cm3 = NA,
            .oxide_capacitance_F_per_cm2 = cox,
            .flat_band_V = VFB,
            .mobility_cm2_Vs = mu,
            .form = form,
            .cutoff_fraction = cutoff};
}

}  // namespace

TEST_CASE("pao_sah: at V = 0 the surface potential is the MOS-C's", "[pao_sah]") {
    const ref::MosInput mi{.q_C = q,
                           .eps_F_per_cm = eps,
                           .thermal_voltage_V = VT,
                           .ni_cm3 = ni,
                           .net_doping_cm3 = -NA,
                           .oxide_capacitance_F_per_cm2 = cox,
                           .flat_band_V = VFB};
    for (const double VG : {-0.5, 0.0, 0.5, 1.5}) {
        for (const int form : {1, 2}) {
            const auto a = ref::pao_sah_inversion_charge(input(form), VG, 0.0);
            const auto b = ref::mos_point(mi, VG);
            REQUIRE(a.has_value());
            REQUIRE(b.has_value());
            CAPTURE(VG, form);
            CHECK(std::abs(((a->surface_V - b->band_bending_V) / b->band_bending_V).to_double()) <= 1e-24);
        }
    }
}

TEST_CASE("pao_sah: the bulk ambiguity is small in strong inversion", "[pao_sah][.pao_sah_release]") {
    for (const double VD : {0.05, 0.5}) {
        const auto a = ref::pao_sah(input(1, 1.0), 1.5, VD);
        const auto b = ref::pao_sah(input(2, 1.0), 1.5, VD);
        const auto c = ref::pao_sah(input(1, 0.5), 1.5, VD);
        REQUIRE(a.has_value());
        REQUIRE(b.has_value());
        REQUIRE(c.has_value());
        CAPTURE(VD, a->quadrature_error);
        CHECK(a->quadrature_error <= 1e-18);
        CHECK(std::abs(((b->current - a->current) / a->current).to_double()) <= 1e-10);
        CHECK(std::abs(((c->current - a->current) / a->current).to_double()) <= 1e-8);
        CHECK(c->current.to_double() != a->current.to_double());  // the cutoff is read (it adds the weak-inversion tail)
        CHECK(a->current.hi > 0.0);
        CHECK(a->surface_drain_V > a->surface_source_V);  // V_D > 0 raises psi_s toward the drain
    }
}

TEST_CASE("pao_sah: the conductance identity", "[pao_sah][.pao_sah_release]") {
    for (const auto [VG, VD] : {std::pair{0.5, 0.05}, std::pair{1.5, 0.5}}) {  // weak and strong inversion
        const auto e = ref::pao_sah_conductance_identity(input(), VG, VD, 0x1p-17);
        REQUIRE(e.has_value());
        CAPTURE(VG, VD, *e);
        CHECK(*e <= 1e-12);
    }
}

TEST_CASE("pao_sah: limits", "[pao_sah][.tvr4_reference_release]") {
    // Linear region: Q / (mu (-Q_I(0)) V_D) -> 1 as V_D -> 0 (first order in V_D / V_T).
    const auto q0 = ref::pao_sah_inversion_charge(input(), 1.5, 0.0);
    const auto lin = ref::pao_sah(input(), 1.5, 1e-5);
    REQUIRE(q0.has_value());
    REQUIRE(lin.has_value());
    CHECK(std::abs((lin->current / (DD(mu) * -q0->charge * DD(1e-5))).to_double() - 1.0) <= 1e-3);
    // Subthreshold: S = dV_G / dlog10 Q against the depletion formula at psi_s(V_G).
    const double h = 1e-3, VG = 0.1;
    const auto lo = ref::pao_sah(input(), VG - h, 0.05), hi = ref::pao_sah(input(), VG + h, 0.05),
               mid = ref::pao_sah(input(), VG, 0.05);
    REQUIRE(lo.has_value());
    REQUIRE(hi.has_value());
    REQUIRE(mid.has_value());
    const double S = 2.0 * h / std::log10((hi->current / lo->current).to_double());
    const double psi = mid->surface_source_V.to_double();
    CAPTURE(psi, mid->bulk_potential_V.to_double());
    REQUIRE(psi > mid->bulk_potential_V.to_double());  // weak inversion
    REQUIRE(psi < 2.0 * mid->bulk_potential_V.to_double());
    const double Cdep = std::sqrt(q * eps * NA / (2.0 * psi));
    const double S_formula = std::log(10.0) * VT * (1.0 + Cdep / cox);
    CAPTURE(S, S_formula);
    CHECK(std::abs(S / S_formula - 1.0) <= 5e-2);
    // Brews' charge sheet in strong inversion.
    const auto cs = ref::charge_sheet(input(), 1.5, 0.05);
    const auto ps = ref::pao_sah(input(), 1.5, 0.05);
    REQUIRE(cs.has_value());
    REQUIRE(ps.has_value());
    CHECK(cs->second_method <= 1e-24);
    CHECK(std::abs((cs->current / ps->current).to_double() - 1.0) <= 5e-2);
}

TEST_CASE("pao_sah: refusals", "[pao_sah]") {
    CHECK_FALSE(ref::pao_sah(input(), VFB, 0.1).has_value());
    CHECK_FALSE(ref::pao_sah(input(), 1.0, -0.1).has_value());
    CHECK_FALSE(ref::pao_sah(input(3), 1.0, 0.1).has_value());
    CHECK_FALSE(ref::pao_sah(input(1, 0.0), 1.0, 0.1).has_value());
    CHECK_FALSE(ref::pao_sah_conductance_identity(input(), 1.0, 0.1, 0.06).has_value());
    CHECK_FALSE(ref::charge_sheet(input(), VFB + 0.01, 0.05).has_value());  // psi_s <= V_T
}
