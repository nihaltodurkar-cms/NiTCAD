// L0 tests of the band-to-band references on the p-i-n potential (verify/reference/tunnel_generation.hpp; DECISIONS
// Unit 55 17.2 (m)-(o), 17.3, TVB-24). The closed forms below hold where the depletion layers' field is linear in x
// (carriers negligible) and the intrinsic region's field uniform, which is where the generation lives; the Debye
// tails, where they fail, carry G below e^-500 of the peak.
//   - local Kane: P = A E_i^2 e^(-B / E_i) d + (eps / (q N_A) + eps / (q N_D)) int_0^E_i A E^2 e^(-B / E) dE (the
//     sides by dx = eps dE / (q N)), to 1e-10; the second method within 1e-12;
//   - nonlocal Kane, one path: a start inside the i region whose crossing is inside it has F = E_i and l = E_g / E_i;
//     a start at depth a = 0.05 W_p inside the p depletion layer (the holes' share of the field there e^-s ~ 1e-10)
//     has l = a + (E_g - E_i a + q N_A a^2 / (2 eps)) / E_i; both to 1e-10; deep in the p bulk a path at B / F = 675
//     has the full rate and one at 710 none (the cut at 700);
//   - nonlocal Kane pairs: G(E_i)(d - E_g / E_i) <= P <= G(E_i)(d + W_p + W_n) (the mean field of a path never exceeds
//     E_i), the second method within 1e-10; an independent quadrature over the starts of each path's own rate within
//     1e-10, on that device and on a one-sided (1e19/1e18) 2 nm i region (every path crossing a region edge);
//   - WKB, one path inside the i region: I_k = l pi sqrt(m_r E_g) / (4 hbar) and I_ik = l pi hbar / (2 sqrt(m_r E_g))
//     (the full-gap integrals in a uniform field), the slope q E_i, k_m^2 from the flat-field levels' windows, and the
//     rate from these by model_generation's path formula (m^-3 to cm^-3) to 1e-9; a crossing beyond the n side's
//     flat-field level has k_m^2 = 0 and no rate; a
//     start below the p side's flat-field level has k_m^2 = 0 and no rate; E at the flat-field levels;
//   - WKB pairs: the second method within 1e-10; the slope coefficient c_1 P against int G E' / (2 E) dx over the
//     starts, E' = (q N_A / eps)(1 - e^-s) on the p side, -(q N_D / eps)(1 - e^-(B - s)) on the n side and 0 in the i
//     region, by an independent quadrature, to 1e-6;
//   - bad inputs are refused.
// The cases that integrate are Release-only ([.tvr5a_reference_release]).
#include "NiTCAD/verify/reference/tunnel_generation.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/verify/reference/gauss_kronrod.hpp"
#include "NiTCAD/verify/reference/model_generation.hpp"

namespace ref = NiTCAD::verify::reference;
using ref::DD;

namespace {

// The reference takes m0 and hbar from base (model_generation.hpp), so the closed forms here do too.
constexpr double q = 1.602176634e-19, VT = 0.025851999786435535, eps0 = 8.8541878128e-14;
constexpr double hbar = NiTCAD::base::hbar_J_s, m0 = NiTCAD::base::m0_kg;
constexpr double si_eps = 11.7 * eps0, si_ni = 9.673e9, si_gap = 1.1242;
constexpr double gaas_eps = 12.9 * eps0, gaas_ni = 2.1e6, gaas_gap = 1.4223;
constexpr ref::KaneInput kane{3.5e21, 1.03e8};
constexpr ref::WkbInput gaas_wkb{gaas_gap, 0.067, 0.5, 1e3};

ref::PinJunction silicon(double d, double V) {
    return *ref::pin_junction({q, si_eps, VT, si_ni, 1e19, 0.0, 1e19, d}, V);
}
ref::PinJunction gaas(double d, double V) {
    return *ref::pin_junction({q, gaas_eps, VT, gaas_ni, 1e19, 0.0, 1e19, d}, V);
}

double rel(const DD& a, const DD& b) { return std::abs(((a - b) / b).to_double()); }

DD G(const DD& E) { return DD(kane.A_per_cm3_s) * E * E * exp(-(DD(kane.B_V_per_cm) / E)); }

}  // namespace

TEST_CASE("tunnel_generation: local Kane against the depletion form", "[.tvr5a_reference_release]") {
    for (const double V : {-3.0, -6.0}) {
        const double d = 1e-6;
        const ref::PinJunction j = silicon(d, V);
        const auto P = ref::local_kane_pairs(j, kane);
        REQUIRE(P.has_value());
        const DD Ei = j.field_p_edge_V_per_cm;
        const auto side = ref::qag15([&](const DD& E) { return G(E); }, DD(kane.B_V_per_cm / 700.0), Ei, 0.0, 1e-26);
        REQUIRE(side.status == ref::QuadratureStatus::ok);
        const DD expected = G(Ei) * DD(d) + DD(2.0 * si_eps / (q * 1e19)) * side.value;
        CAPTURE(V, P->pairs_cm2_s.to_double(), expected.to_double(), P->second_cm2_s.to_double());
        REQUIRE(rel(P->pairs_cm2_s, expected) < 1e-10);
        REQUIRE(rel(P->second_cm2_s, P->pairs_cm2_s) < 1e-12);
    }
}

TEST_CASE("tunnel_generation: nonlocal Kane paths in the closed forms", "[.tvr5a_reference_release]") {
    const double d = 2e-6;
    const ref::PinJunction j = silicon(d, -3.0);
    const DD Ei = j.field_p_edge_V_per_cm;
    const DD l_i = DD(si_gap) / Ei;
    REQUIRE(l_i.hi < 0.5 * d);
    {
        const auto s = ref::pin_level(j, 0.25 * d);
        REQUIRE(s.has_value());
        const auto p = ref::kane_path(j, si_gap, kane, *s);
        REQUIRE(p.has_value());
        REQUIRE(p->exists);
        CAPTURE(p->field_V_per_cm.to_double(), Ei.to_double());
        REQUIRE(rel(p->field_V_per_cm, Ei) < 1e-10);
        REQUIRE(rel(p->length_cm, l_i) < 1e-10);
        REQUIRE(rel(p->second_length_cm, p->length_cm) < 1e-20);
        REQUIRE(rel(p->rate_cm3_s, G(Ei)) < 1e-9);
    }
    {
        // A start at depth a in the p side, near enough to the i region that the holes are negligible on [-a, 0]:
        // E(x) = E_i + q N_A x / eps there, so the potential rises by E_i a - q N_A a^2 / (2 eps) to x = 0.
        // (Here s_0 is only about 25, so the holes' share of the field, e^-s / (2 (s - 1)), is 1e-7 at a = 0.3 W_p;
        // a = 0.05 W_p keeps it near 1e-12.)
        const double a = 0.05 * j.depletion_p_cm.to_double();
        const auto s = ref::pin_level(j, -a);
        REQUIRE(s.has_value());
        REQUIRE(std::exp(-s->to_double()) < 1e-9);
        const auto p = ref::kane_path(j, si_gap, kane, *s);
        REQUIRE(p.has_value());
        REQUIRE(p->exists);
        const DD rise = Ei * DD(a) - DD(q * 1e19 / (2.0 * si_eps)) * DD(a) * DD(a);
        const DD expected = DD(a) + (DD(si_gap) - rise) / Ei;
        CAPTURE(p->length_cm.to_double(), expected.to_double());
        REQUIRE(rel(p->length_cm, expected) < 1e-10);
        REQUIRE(p->field_V_per_cm < Ei);
    }
    {  // no crossing from the n side
        const auto p = ref::kane_path(j, si_gap, kane, j.barrier - DD(1.0));
        REQUIRE(p.has_value());
        REQUIRE(!p->exists);
    }
    {
        // The cut: deep in the p bulk the path lengthens and F falls; a start with B / F = 675 has the full rate, one
        // with B / F = 710 none. The start by bisection in ln s on the path's own F.
        const auto start_at = [&](double ratio) {
            double lo = std::log(1e-30), hi = std::log(j.s0.to_double());
            for (int it = 0; it < 200; ++it) {
                const double mid = 0.5 * (lo + hi);
                const auto p = ref::kane_path(j, si_gap, kane, DD(std::exp(mid)));
                (kane.B_V_per_cm / p->field_V_per_cm.to_double() > ratio ? lo : hi) = mid;
            }
            return *ref::kane_path(j, si_gap, kane, DD(std::exp(hi)));
        };
        const ref::PathReference in_cut = start_at(675.0), past_cut = start_at(710.0);
        CAPTURE(kane.B_V_per_cm / in_cut.field_V_per_cm.to_double(), kane.B_V_per_cm / past_cut.field_V_per_cm.to_double());
        REQUIRE(std::abs(kane.B_V_per_cm / in_cut.field_V_per_cm.to_double() - 675.0) < 1e-6);
        REQUIRE(in_cut.rate_cm3_s.hi > 0.0);
        REQUIRE(rel(in_cut.rate_cm3_s, G(in_cut.field_V_per_cm)) < 1e-28);
        REQUIRE(past_cut.rate_cm3_s.hi == 0.0);
    }
}

TEST_CASE("tunnel_generation: nonlocal Kane pairs within their bounds","[.tvr5a_reference_release]") {
    const double d = 2e-6;
    const ref::PinJunction j = silicon(d, -3.0);
    const auto P = ref::nonlocal_kane_pairs(j, si_gap, kane);
    if (!P) FAIL(P.error());
    const DD Ei = j.field_p_edge_V_per_cm;
    const DD low = G(Ei) * (DD(d) - DD(si_gap) / Ei);
    const DD high = G(Ei) * (DD(d) + j.depletion_p_cm + j.depletion_n_cm);
    CAPTURE(P->pairs_cm2_s.to_double(), low.to_double(), high.to_double(), P->second_cm2_s.to_double());
    REQUIRE(P->pairs_cm2_s > low);
    REQUIRE(P->pairs_cm2_s < high);
    REQUIRE(rel(P->second_cm2_s, P->pairs_cm2_s) < 1e-10);
}

TEST_CASE("tunnel_generation: nonlocal Kane pairs over the starts, independently", "[.tvr5a_reference_release]") {
    // int G(s_a) dx_a from each path's own rate, dx_a = V_T ds_a / E(s_a), against the pairs: on the 20 nm device and
    // on a one-sided 2 nm one, where every path crosses a region edge (so E at the start and at the crossing differ). Below
    // s = 1e-3 (about 7 Debye lengths into the p bulk) every path's rate is below 1e-30 of G(E_i), and the rate falls as
    // the start moves deeper, so the part left out is below that times the distance.
    // (Weighting by dx at the crossing instead of the start gives the same P on any device: G depends on the start
    // only through l = x_b - x_a, so the two differ by int G(l) dl, which vanishes as l runs from large to large.)
    for (const double d : {2e-6, 2e-7}) {
        const ref::PinJunction j =
            d > 1e-6 ? silicon(d, -3.0) : *ref::pin_junction({q, si_eps, VT, si_ni, 1e19, 0.0, 1e18, d}, -3.0);
        const auto P = ref::nonlocal_kane_pairs(j, si_gap, kane);
        if (!P) FAIL(P.error());
        const DD Ei = j.field_p_edge_V_per_cm;
        const DD g = DD(si_gap) / j.VT, lo(1e-3), top = j.barrier - g;
        REQUIRE(ref::kane_path(j, si_gap, kane, lo)->rate_cm3_s < DD(1e-30) * G(Ei));
        std::vector<DD> at{lo, j.s0 - g, j.s0, j.sd - g, j.sd, top};
        std::erase_if(at, [&](const DD& s) { return s < lo || top < s; });
        std::sort(at.begin(), at.end(), [](const DD& a, const DD& b) { return a < b; });
        DD sum(0.0);
        for (std::size_t k = 0; k + 1 < at.size(); ++k) {
            const auto piece = ref::qag15(
                [&](const DD& s) {
                    const auto p = ref::kane_path(j, si_gap, kane, s);
                    return p->exists ? p->rate_cm3_s * j.VT / ref::pin_field(j, s) : DD(0.0);
                },
                at[k], at[k + 1], 0.0, 1e-12);
            REQUIRE(piece.status == ref::QuadratureStatus::ok);
            sum += piece.value;
        }
        CAPTURE(d, sum.to_double(), P->pairs_cm2_s.to_double());
        REQUIRE(rel(sum, P->pairs_cm2_s) < 1e-10);
    }
}

TEST_CASE("tunnel_generation: a WKB path in a uniform field", "[.tvr5a_reference_release]") {
    const double d = 2e-6;
    const ref::PinJunction j = gaas(d, -4.0);
    const DD Ei = j.field_p_edge_V_per_cm;
    const auto lv = ref::wkb_levels(j, gaas_wkb);
    REQUIRE(lv.has_value());
    REQUIRE(rel(ref::pin_field(j, lv->s_lo), DD(1e3)) < 1e-25);
    REQUIRE(rel(ref::pin_field(j, lv->s_hi), DD(1e3)) < 1e-25);
    REQUIRE(lv->s_lo < j.s0);
    REQUIRE(lv->s_hi > j.sd);
    const DD l_i = DD(gaas_gap) / Ei;
    REQUIRE(l_i.hi < 0.5 * d);
    const auto s = ref::pin_level(j, 0.25 * d);
    REQUIRE(s.has_value());
    const auto p = ref::wkb_path(j, gaas_wkb, *s);
    if (!p) FAIL(p.error());
    REQUIRE(p->exists);
    const DD mc = DD(0.067 * m0), mv = DD(0.5 * m0), mr = mc * mv / (mc + mv);
    const DD Eg = DD(gaas_gap) * DD(q), root = sqrt(mr * Eg), h(hbar), L = l_i * DD(1e-2);
    const DD Ik = L * ref::pi() * root / (DD(4.0) * h), Iik = L * ref::pi() * h / (DD(2.0) * root);
    CAPTURE(p->I_k.to_double(), Ik.to_double(), p->I_ik_m2.to_double(), Iik.to_double());
    REQUIRE(rel(p->I_k, Ik) < 1e-10);
    REQUIRE(rel(p->I_ik_m2, Iik) < 1e-10);
    REQUIRE(rel(p->slope_J_per_m, DD(q) * Ei * DD(100.0)) < 1e-12);  // the i region's carriers: 1e-13
    const DD mcd = DD(0.067) * DD(m0), mvd = DD(0.5) * DD(m0);
    const DD kv = DD(2.0) * mvd * DD(q) * DD(VT) / (h * h), kc = DD(2.0) * mcd * DD(q) * DD(VT) / (h * h);
    const DD wv = kv * (*s - lv->s_lo), wc = kc * (lv->s_hi - *s - DD(gaas_gap) / DD(VT));
    const DD km2 = wv < wc ? wv : wc;
    CAPTURE(p->km2_per_m2.to_double(), km2.to_double());
    REQUIRE(rel(p->km2_per_m2, km2) < 1e-24);
    // The rate from the closed-form path quantities by model_generation's formula, m^-3 to cm^-3.
    const DD G = ref::model_generation::path_rate({DD(q) * Ei * DD(100.0), Ik, Iik, km2})[0] * DD(1e-6);
    CAPTURE(p->rate_cm3_s.to_double(), G.to_double());
    REQUIRE(rel(p->rate_cm3_s, G) < 1e-9);
    REQUIRE(rel(p->second_rate_cm3_s, p->rate_cm3_s) < 1e-12);
    // Below the p side's flat-field level the backward window closes: no rate.
    const auto flat = ref::wkb_path(j, gaas_wkb, lv->s_lo * DD(0.5));
    REQUIRE(flat.has_value());
    REQUIRE(flat->exists);
    REQUIRE(flat->km2_per_m2.hi == 0.0);
    REQUIRE(flat->rate_cm3_s.hi == 0.0);
    // A crossing beyond the n side's flat-field level closes the forward window: no rate either.
    const DD g = DD(gaas_gap) / j.VT;
    const auto beyond = ref::wkb_path(j, gaas_wkb, ldexp(lv->s_hi - g + (j.barrier - g), -1));
    REQUIRE(beyond.has_value());
    REQUIRE(beyond->exists);
    REQUIRE(beyond->km2_per_m2.hi == 0.0);
    REQUIRE(beyond->rate_cm3_s.hi == 0.0);
}

TEST_CASE("tunnel_generation: WKB pairs and the slope coefficient", "[.tvr5a_reference_release]") {
    const double d = 1e-6;
    const ref::PinJunction j = gaas(d, -4.0);
    const auto P = ref::wkb_pairs(j, gaas_wkb);
    if (!P) FAIL(P.error());
    CAPTURE(P->pairs_cm2_s.to_double(), P->second_cm2_s.to_double());
    REQUIRE(P->pairs_cm2_s.hi > 0.0);
    REQUIRE(rel(P->second_cm2_s, P->pairs_cm2_s) < 1e-10);
    const auto c1 = ref::wkb_slope_coefficient(j, gaas_wkb);
    if (!c1) FAIL(c1.error());
    const auto lv = ref::wkb_levels(j, gaas_wkb);
    REQUIRE(lv.has_value());
    // int G E' / (2 E) dx over the starts, dx = V_T ds / E, with E' = (q N_A / eps)(1 - e^-s) on the p side and
    // -(q N_D / eps)(1 - e^-(B - s)) on the n side (the majority carriers' share; the minority's and the i region's
    // carriers are below 1e-30).
    bool ok = true;
    const auto side = [&](const DD& from, const DD& to, double sign, bool p_side) {
        const auto r = ref::qag15(
            [&](const DD& s) {
                const auto p = ref::wkb_path(j, gaas_wkb, s);
                if (!p) ok = false;
                const DD E = ref::pin_field(j, s);
                const DD bulk = p_side ? exp(-s) : exp(s - j.barrier);
                return p ? p->rate_cm3_s * DD(VT) / (E * E) * DD(sign) * (DD(1.0) - bulk) : DD(0.0);
            },
            from, to, 0.0, 1e-14);
        if (r.status != ref::QuadratureStatus::ok) ok = false;
        return r.value;
    };
    const DD sum = side(lv->s_lo, j.s0, 1.0, true) + side(j.sd, lv->s_hi - DD(gaas_gap / VT), -1.0, false);
    REQUIRE(ok);
    const DD expected = DD(q * 1e19 / (2.0 * gaas_eps)) * sum / P->pairs_cm2_s;
    CAPTURE(c1->to_double(), expected.to_double());
    REQUIRE(c1->hi > 0.0);
    REQUIRE(rel(*c1, expected) < 1e-6);
}

// Starts within 1e-9 ... 1e-17 (relative) of each cut of the pairs' quadrature, at -0.5, -1 and -2 V: every path
// forms, its two methods' rates agree, and the pairs form. A start just above s_d - g puts the n edge ~1e-14 below
// delta = 1, where the second method's sliver stops on roundoff (found at TVB-24's rerun, 17.27).
TEST_CASE("tunnel_generation: WKB paths near the quadrature's cuts", "[.tvr5a_reference_release]") {
    for (const double V : {-0.5, -1.0, -2.0}) {
        CAPTURE(V);
        const ref::PinJunction j = gaas(1e-6, V);
        const auto lv = ref::wkb_levels(j, gaas_wkb);
        REQUIRE(lv.has_value());
        const DD g = DD(gaas_gap / VT);
        const DD lo = lv->s_lo, hi = lv->s_hi - g;
        const DD mv(gaas_wkb.hole_mass), mc(gaas_wkb.electron_mass);
        const DD meet = (mv * lv->s_lo + mc * (lv->s_hi - g)) / (mv + mc);
        for (const DD& cut : {lo, hi, j.s0, j.sd, j.s0 - g, j.sd - g, meet}) {
            if (cut < lo || cut > hi) continue;
            for (int e = 9; e <= 17; ++e) {
                for (const double sign : {-1.0, 1.0}) {
                    const DD s = cut + DD(sign * std::pow(10.0, -e)) * cut;
                    if (s < lo || s > hi) continue;
                    CAPTURE(cut.to_double(), sign, e);
                    const auto p = ref::wkb_path(j, gaas_wkb, s);
                    if (!p) FAIL(p.error());
                    if (p->rate_cm3_s.hi > 0.0) REQUIRE(rel(p->second_rate_cm3_s, p->rate_cm3_s) < 1e-10);
                }
            }
        }
        const auto P = ref::wkb_pairs(j, gaas_wkb);
        if (!P) FAIL(P.error());
        REQUIRE(P->pairs_cm2_s.hi > 0.0);
        REQUIRE(rel(P->second_cm2_s, P->pairs_cm2_s) < 1e-10);
    }
}

TEST_CASE("tunnel_generation: bad inputs are refused") {
    const ref::PinJunction j = silicon(1e-6, -1.0);
    REQUIRE(!ref::local_kane_pairs(j, {0.0, 1e8}).has_value());
    REQUIRE(!ref::local_kane_pairs(j, {1e21, -1.0}).has_value());
    REQUIRE(!ref::kane_path(j, -1.0, kane, j.s0).has_value());
    REQUIRE(!ref::kane_path(j, si_gap, kane, DD(0.0)).has_value());
    REQUIRE(!ref::kane_path(j, si_gap, kane, j.barrier).has_value());
    REQUIRE(!ref::wkb_path(j, {gaas_gap, 0.0, 0.5, 1e3}, j.s0).has_value());
    REQUIRE(!ref::wkb_path(j, {gaas_gap, 1.5, 1.5, 1e3}, j.s0).has_value());  // m_r = 0.75 m0: u = 2/3 < 1
    REQUIRE(!ref::wkb_levels(j, {gaas_gap, 0.067, 0.5, 1e9}).has_value());   // the flat field above the edges' field
}
