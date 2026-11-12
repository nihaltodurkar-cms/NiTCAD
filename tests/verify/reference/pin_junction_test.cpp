// L0 tests of the p-i-n reference (verify/reference/pin_junction.hpp; DECISIONS Unit 55 17.3, 17.13 (e), TVB-24).
//   - the limit d -> 0: a 1e-15 cm i region gives qe_junction's barrier and peak field (a different reference, two
//     regions) to 1e-8, at 0 and -1 V, for a symmetric and a one-sided junction;
//   - the n side's field from its bulk distance r equals pin_field(B - r) inside the n side, and is linear in r as
//     r -> 0 (r = 1e-20, where B - r has lost r);
//   - an intrinsic region carries a uniform field: E at the two edges agree to 1e-10 at -3 and -20 V (the latter past
//     qe_junction's range), the length and matching residuals at rounding;
//   - Gauss across a doped i region: E(0) - E(d) = q C_i d / eps to 1e-8, for donors and for acceptors, which also
//     checks that the root honours d;
//   - positions: x(s_0) = 0 and x(s_d) = d; x -> s -> x round trips on both sides and in the i region, x rising with s;
//     deep in a bulk (s = 1e-12 to 1e-8) the distance is the log-linear one, ln(s / s') times the Debye length;
//   - an interior field minimum (zero bias, a long intrinsic region of a symmetric device) at the middle, below both
//     edges;
//   - dE/dx against a central difference of E in the level, times E / V_T, in each region;
//   - the depletion approximation as a sanity bound (the depletion widths within 3 %, a one-sided device);
//   - the second method: Chebyshev collocation on the finite domain agrees on s_0, s_d, the edge and mid fields and the
//     levels at query points to 1e-9 (three devices);
//   - bad inputs are refused.
// The cases that integrate are Release-only ([.tvr5a_reference_release]).
#include "NiTCAD/verify/reference/pin_junction.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <vector>

#include "NiTCAD/verify/reference/qe_junction.hpp"

namespace ref = NiTCAD::verify::reference;
using ref::DD;

namespace {

constexpr double q = 1.602176634e-19, VT = 0.025851999786435535, eps = 11.7 * 8.8541878128e-14, ni = 9.673e9;

ref::PinInput input(double NA, double Ci, double ND, double d) { return {q, eps, VT, ni, NA, Ci, ND, d}; }

double rel(const DD& a, const DD& b) { return std::abs(((a - b) / b).to_double()); }

}  // namespace

TEST_CASE("pin_junction: a vanishing i region is the two-region junction", "[.tvr5a_reference_release]") {
    for (const auto& [NA, ND] : {std::array{1e19, 1e19}, std::array{1e16, 1e19}}) {
        for (const double V : {0.0, -1.0}) {
            const auto p = ref::pin_junction(input(NA, 0.0, ND, 1e-15), V);
            const auto t = ref::qe_junction({q, eps, VT, ni, NA, ND}, V);
            REQUIRE(p.has_value());
            REQUIRE(t.has_value());
            const DD E2 = t->charge_C_per_cm2 / DD(eps);
            CAPTURE(NA, ND, V, p->field_p_edge_V_per_cm.to_double(), E2.to_double());
            REQUIRE(rel(p->VT * p->barrier, t->barrier_V) < 1e-24);
            REQUIRE(rel(p->field_p_edge_V_per_cm, E2) < 1e-8);
            REQUIRE(rel(p->field_n_edge_V_per_cm, E2) < 1e-8);
        }
    }
}

TEST_CASE("pin_junction: the n side's field from its bulk", "[.tvr5a_reference_release]") {
    // pin_field(B - r) at r = 0.2 and 0.9 of B - s_d; E -> V_T r sqrt(scale N_D / 2) as r -> 0 (Phi_n = N_D r^2 / 2 to
    // relative r / 3), where B - r would have lost r.
    const auto p = ref::pin_junction(input(1e16, 0.0, 1e19, 1e-6), -2.0);
    REQUIRE(p.has_value());
    const DD side = p->barrier - p->sd;
    for (const double f : {0.2, 0.9}) {
        const DD r = side * DD(f);
        CAPTURE(f, side.to_double());
        REQUIRE(rel(ref::pin_field_n_bulk(*p, r), ref::pin_field(*p, p->barrier - r)) < 1e-28);
    }
    const DD r(1e-20);
    const DD linear = p->VT * r * sqrt(p->scale * DD(1e19) / DD(2.0));
    REQUIRE(rel(ref::pin_field_n_bulk(*p, r), linear) < 1e-19);
}

TEST_CASE("pin_junction: an intrinsic region carries a uniform field", "[.tvr5a_reference_release]") {
    for (const double V : {-3.0, -20.0}) {
        const auto p = ref::pin_junction(input(1e19, 0.0, 1e19, 1e-6), V);
        REQUIRE(p.has_value());
        CAPTURE(V, p->s0.to_double(), p->sd.to_double(), p->barrier.to_double(), p->field_p_edge_V_per_cm.to_double());
        REQUIRE(p->s0.hi > 0.0);
        REQUIRE(p->sd > p->s0);
        REQUIRE(p->sd < p->barrier);
        REQUIRE(rel(p->field_n_edge_V_per_cm, p->field_p_edge_V_per_cm) < 1e-10);
        REQUIRE(rel(p->field_min_i_V_per_cm, p->field_p_edge_V_per_cm) < 1e-10);
        REQUIRE(p->length_residual_cm < 1e-24 * 1e-6);
        REQUIRE(p->matching_residual < 1e-25);
        // The field times the i region's length is the i region's share of the barrier.
        REQUIRE(rel(p->field_p_edge_V_per_cm * DD(1e-6), p->VT * (p->sd - p->s0)) < 1e-10);
    }
}

TEST_CASE("pin_junction: Gauss across a doped i region", "[.tvr5a_reference_release]") {
    const double d = 1e-6;
    for (const double Ci : {2e17, -2e17}) {
        const auto p = ref::pin_junction(input(1e19, Ci, 1e19, d), -3.0);
        REQUIRE(p.has_value());
        const double expected = q * Ci * d / eps;
        const double got = (p->field_p_edge_V_per_cm - p->field_n_edge_V_per_cm).to_double();
        CAPTURE(Ci, got, expected);
        REQUIRE(std::abs(got / expected - 1.0) < 1e-8);
    }
}

TEST_CASE("pin_junction: positions and levels", "[.tvr5a_reference_release]") {
    const double d = 1e-6;
    const auto p = ref::pin_junction(input(5e18, 0.0, 2e19, d), -2.0);
    REQUIRE(p.has_value());
    const auto x0 = ref::pin_position(*p, p->s0);
    const auto xd = ref::pin_position(*p, p->sd);
    REQUIRE(x0.has_value());
    REQUIRE(xd.has_value());
    REQUIRE(std::abs(x0->to_double()) < 1e-30);
    REQUIRE(std::abs(xd->to_double() / d - 1.0) < 1e-24);
    const DD rn = p->barrier - p->sd;
    const std::vector<DD> levels{p->s0 * DD(1e-6),
                                 p->s0 * DD(0.01),
                                 p->s0 * DD(0.5),
                                 p->s0 + (p->sd - p->s0) * DD(0.3),
                                 p->s0 + (p->sd - p->s0) * DD(0.9),
                                 p->barrier - rn * DD(0.5),
                                 p->barrier - rn * DD(0.01),
                                 p->barrier - rn * DD(1e-6)};
    double last = -1.0;
    for (const DD& s : levels) {
        const auto x = ref::pin_position(*p, s);
        REQUIRE(x.has_value());
        const auto back = ref::pin_level(*p, x->to_double());
        if (!back) FAIL(back.error());
        CAPTURE(s.to_double(), x->to_double(), back->to_double());
        REQUIRE(x->to_double() > last);
        last = x->to_double();
        // x is rounded to a double for pin_level, so the level returns to that rounding times ds/dx.
        const double dsdx = (ref::pin_field(*p, s) / p->VT).to_double();
        REQUIRE(std::abs((*back - s).to_double()) <= 4e-16 * std::abs(x->to_double()) * dsdx + 1e-26 * s.to_double());
    }
    REQUIRE(!ref::pin_position(*p, DD(0.0)).has_value());
    REQUIRE(!ref::pin_position(*p, p->barrier).has_value());
    // Deep in the p bulk the field is linear in the level (Phi = (p_p + n_p) s^2 / 2 to relative O(s)), so
    // x(s) - x(s') = ln(s / s') / sqrt(scale (p_p + n_p) / 2) to relative O(s).
    const auto xa = ref::pin_position(*p, DD(1e-12)), xb = ref::pin_position(*p, DD(1e-8));
    REQUIRE(xa.has_value());
    REQUIRE(xb.has_value());
    const DD debye = DD(1.0) / sqrt(p->scale * (p->p_p + exp(p->ln_np)) / DD(2.0));
    CAPTURE((*xb - *xa).to_double(), (log(DD(1e4)) * debye).to_double());
    REQUIRE(rel(*xb - *xa, log(DD(1e4)) * debye) < 1e-7);
}

TEST_CASE("pin_junction: the field's minimum inside an intrinsic region", "[.tvr5a_reference_release]") {
    // At zero bias the carriers in a long intrinsic region make the field dip between the edges; in a symmetric device
    // the minimum lies at the middle, where n = p.
    const double d = 2e-6;
    const auto p = ref::pin_junction(input(1e17, 0.0, 1e17, d), 0.0);
    REQUIRE(p.has_value());
    const auto mid = ref::pin_level(*p, 0.5 * d);
    if (!mid) FAIL(mid.error());
    const DD Emid = ref::pin_field(*p, *mid);
    CAPTURE(p->field_min_i_V_per_cm.to_double(), Emid.to_double(), p->field_p_edge_V_per_cm.to_double());
    REQUIRE(p->field_min_i_V_per_cm < p->field_p_edge_V_per_cm);
    REQUIRE(p->field_min_i_V_per_cm < p->field_n_edge_V_per_cm);
    REQUIRE(rel(p->field_min_i_V_per_cm, Emid) < 1e-20);
}

TEST_CASE("pin_junction: dE/dx against a difference of the field", "[.tvr5a_reference_release]") {
    const auto p = ref::pin_junction(input(1e19, 1e17, 1e19, 1e-6), -2.0);
    REQUIRE(p.has_value());
    const DD rn = p->barrier - p->sd;
    for (const DD& s : {p->s0 * DD(0.2), p->s0 * DD(0.9), p->s0 + (p->sd - p->s0) * DD(0.5),
                        p->barrier - rn * DD(0.9), p->barrier - rn * DD(0.2)}) {
        const DD h = s * DD(1e-9);
        const DD dEds = (ref::pin_field(*p, s + h) - ref::pin_field(*p, s - h)) / ldexp(h, 1);
        const DD expected = dEds * ref::pin_field(*p, s) / p->VT;
        const DD got = ref::pin_field_slope(*p, s);
        CAPTURE(s.to_double(), got.to_double(), expected.to_double());
        REQUIRE(rel(got, expected) < 1e-12);
    }
}

TEST_CASE("pin_junction: the depletion approximation as a sanity bound", "[.tvr5a_reference_release]") {
    // N_A x_p = N_D x_n, and (q / eps)(N_A x_p^2 / 2 + N_A x_p d + N_D x_n^2 / 2) = V_bi - V - 2 V_T (the majority tails'
    // correction, as qe_junction's depletion form).
    const double NA = 5e18, ND = 2e19, d = 1e-6;
    const auto p = ref::pin_junction(input(NA, 0.0, ND, d), -6.0);
    REQUIRE(p.has_value());
    const double U = (p->VT * p->barrier).to_double() - 2.0 * VT;
    const double a = 0.5 * NA + NA * NA / (2.0 * ND), b = NA * d, c = -U * eps / q;
    const double xp = (-b + std::sqrt(b * b - 4.0 * a * c)) / (2.0 * a);
    CAPTURE(p->depletion_p_cm.to_double(), xp);
    REQUIRE(std::abs(p->depletion_p_cm.to_double() / xp - 1.0) < 0.03);
    REQUIRE(std::abs(p->depletion_n_cm.to_double() / (NA * xp / ND) - 1.0) < 0.03);
}

TEST_CASE("pin_junction: Chebyshev collocation agrees", "[.tvr5a_reference_release]") {
    struct Case {
        double NA, Ci, ND, d, V;
    };
    for (const Case& c : {Case{1e19, 0.0, 1e19, 1e-6, -3.0}, Case{1e19, 2e17, 1e19, 1e-6, -1.0},
                          Case{5e18, -1e17, 2e19, 2e-6, 0.0}}) {
        const auto p = ref::pin_junction(input(c.NA, c.Ci, c.ND, c.d), c.V);
        REQUIRE(p.has_value());
        const double debye = std::sqrt(eps * VT / (q * std::min(c.NA, c.ND)));
        const double xl = -(p->depletion_p_cm.to_double() + 40.0 * debye);
        const double xr = c.d + p->depletion_n_cm.to_double() + 40.0 * debye;
        const std::vector<double> query{0.5 * xl, -0.5 * p->depletion_p_cm.to_double(), 0.25 * c.d,
                                        c.d + 0.5 * p->depletion_n_cm.to_double(), 0.5 * (c.d + xr)};
        const auto col = ref::pin_collocation(*p, xl, xr, 24, 4.0, query);
        REQUIRE(col.has_value());
        CAPTURE(c.NA, c.Ci, c.ND, c.V, col->newton_iterations, col->tail_estimate);
        CAPTURE(p->s0.to_double(), col->s0, p->sd.to_double(), col->sd);
        REQUIRE(std::abs(p->s0.to_double() - col->s0) < 1e-9);
        REQUIRE(std::abs(p->sd.to_double() - col->sd) < 1e-9);
        const double Ep = p->field_p_edge_V_per_cm.to_double(), En = p->field_n_edge_V_per_cm.to_double();
        CAPTURE(Ep, col->field_p_edge_V_per_cm, En, col->field_n_edge_V_per_cm);
        REQUIRE(std::abs(col->field_p_edge_V_per_cm / Ep - 1.0) < 1e-9);
        REQUIRE(std::abs(col->field_n_edge_V_per_cm / En - 1.0) < 1e-9);
        const auto smid = ref::pin_level(*p, 0.5 * c.d);
        if (!smid) FAIL(smid.error());
        REQUIRE(std::abs(col->field_mid_V_per_cm / ref::pin_field(*p, *smid).to_double() - 1.0) < 1e-9);
        for (std::size_t k = 0; k < query.size(); ++k) {
            const auto s = ref::pin_level(*p, query[k]);
            if (!s) FAIL(s.error());
            CAPTURE(query[k], s->to_double(), col->levels[k]);
            REQUIRE(std::abs(s->to_double() - col->levels[k]) < 1e-9);
        }
    }
}

TEST_CASE("pin_junction: bad inputs are refused") {
    REQUIRE(!ref::pin_junction(input(1e19, 0.0, 1e19, 1e-6), 0.1).has_value());
    REQUIRE(!ref::pin_junction(input(1e19, 0.0, 1e19, 0.0), -1.0).has_value());
    REQUIRE(!ref::pin_junction(input(0.0, 0.0, 1e19, 1e-6), -1.0).has_value());
    REQUIRE(!ref::pin_junction(input(1e19, 1e19, 1e19, 1e-6), -1.0).has_value());
    REQUIRE(!ref::pin_junction(input(1e19, -1e19, 1e19, 1e-6), -1.0).has_value());
    REQUIRE(!ref::pin_junction(input(1e19, NAN, 1e19, 1e-6), -1.0).has_value());
    const auto p = ref::pin_junction(input(1e19, 0.0, 1e19, 1e-6), -1.0);
    REQUIRE(p.has_value());
    REQUIRE(!ref::pin_collocation(*p, 1e-6, 2e-6, 24, 4.0, {}).has_value());
    REQUIRE(!ref::pin_collocation(*p, -1e-6, 2e-6, 1, 4.0, {}).has_value());
    REQUIRE(!ref::pin_level(*p, -1.0).has_value());
}
