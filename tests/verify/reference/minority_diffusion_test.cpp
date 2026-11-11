// L0 for the minority-diffusion reference (DECISIONS Unit 55 16.3, TVB-15): the linearized lifetime against a
// difference quotient of the full SRH rate; the Green's function against its sinh form, its symmetry and its
// derivative jump; the Green's-function integral against the infinite-domain erfc closed form on a long bar; and
// the collocation (the second method) against the integral.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "NiTCAD/verify/reference/minority_diffusion.hpp"
#include "NiTCAD/verify/reference/special_functions.hpp"

using namespace NiTCAD::verify::reference;

namespace {

double rel(const DD& a, const DD& b) { return std::abs(((a - b) / b).hi); }

DD sinh_dd(const DD& x) { return ldexp(exp(x) - exp(-x), -1); }

}  // namespace

TEST_CASE("minority diffusion: the linearized SRH lifetime", "[minority_diffusion]") {
    const double ni = 9.65e9, p0 = 1e16, n0 = ni * ni / p0, tn = 3e-6, tp = 1e-6;
    const DD inv = srh_linear_inverse_lifetime(n0, p0, ni, tn, tp);
    // Second method: R(n0 + d, p0 + d) / d for the full rate R = (n p - n_ie^2) / (tau_p (n + n_ie) + tau_n (p + n_ie)),
    // with n0 p0 = n_ie^2 exactly in double-double; the remainder is O(d / p0).
    const DD E = DD(n0) * DD(p0);  // the equilibrium product, so R(0) = 0 exactly
    const DD d = DD(p0 * 1e-15);
    const DD n = DD(n0) + d, p = DD(p0) + d;
    const DD R = (n * p - E) / (DD(tp) * (n + DD(ni)) + DD(tn) * (p + DD(ni)));
    CHECK(rel(inv, R / d) <= 1e-14);
    // Strongly p-type: 1 / tau -> 1 / tau_n (1 + n_ie / p0 terms aside).
    CHECK(rel(inv, DD(1.0) / DD(tn)) <= 2.0 * ni / p0 + 1e-15);
    REQUIRE_THROWS_AS(srh_linear_inverse_lifetime(0.0, p0, ni, tn, tp), std::invalid_argument);
    REQUIRE_THROWS_AS(srh_linear_inverse_lifetime(n0, p0, ni, -1.0, tp), std::invalid_argument);
}

TEST_CASE("minority diffusion: the Green's function", "[minority_diffusion]") {
    const MinorityProblem p{.diffusivity_cm2_s = 2.0, .lifetime_s = 1.0 / 18.0, .length_cm = 1.0,
                            .peak_generation = 1.0, .centre_cm = 0.5, .width_cm = 0.1};  // k = 3
    // k from the tau actually passed (1 / 18 rounded to a double), as the problem poses it.
    const DD D(2.0), L(1.0), k = DD(1.0) / sqrt(D * DD(p.lifetime_s));
    for (const double x : {0.1, 0.37, 0.8}) {
        for (const double s : {0.05, 0.37, 0.6, 0.99}) {
            CAPTURE(x, s);
            const double lo = std::min(x, s), hi = std::max(x, s);
            const DD sinh_form = sinh_dd(k * DD(lo)) * sinh_dd(k * (L - DD(hi))) / (D * k * sinh_dd(k * L));
            CHECK(rel(minority_green(p, x, s), sinh_form) <= 64.0 * u_dd);
            CHECK(minority_green(p, x, s) == minority_green(p, s, x));
        }
    }
    // -D d^2 g / dx^2 + g / tau = delta(x - s): the derivative jumps by -1 / D at x = s.
    const double s = 0.37, h = 1e-12;
    const double hr = (s + h) - s, hl = s - (s - h);  // the steps as rounded (exact differences)
    const DD right = (minority_green(p, s + h, s) - minority_green(p, s, s)) / DD(hr);
    const DD left = (minority_green(p, s, s) - minority_green(p, s - h, s)) / DD(hl);
    CHECK(std::abs((right - left).hi + 0.5) <= 1e-10);
    CHECK(minority_green(p, 0.0, s).hi == 0.0);
    CHECK(minority_green(p, 1.0, s).hi == 0.0);
    REQUIRE_THROWS_AS(minority_green(p, 1.5, s), std::invalid_argument);
}

TEST_CASE("minority diffusion: the integral against the erfc closed form on a long bar",
          "[minority_diffusion][.tvr4_reference_release]") {
    // k L = 200 and the Gaussian 50 widths from each end: the ends change u by e^-100, below double-double.
    // Infinite domain: u = G0 / (2 D k) sigma sqrt(pi / 2) e^(k^2 sigma^2 / 2)
    //   [e^(-k y) erfc((k sigma^2 - y) / (sigma sqrt 2)) + e^(k y) erfc((k sigma^2 + y) / (sigma sqrt 2))], y = x - c.
    const double D = 1.5, k = 200.0, sigma = 0.01, c = 0.5, G0 = 3.0;
    const MinorityProblem p{.diffusivity_cm2_s = D, .lifetime_s = 1.0 / (D * k * k), .length_cm = 1.0,
                            .peak_generation = G0, .centre_cm = c, .width_cm = sigma};
    for (const double x : {0.45, 0.5, 0.503, 0.55}) {
        CAPTURE(x);
        const DD y = DD(x) - DD(c), S(sigma);
        const DD K = DD(1.0) / sqrt(DD(D) * DD(p.lifetime_s));  // tau as passed (rounded), so k is about 200
        const DD a = K * S * S, b = S * sqrt2();
        const DD bracket = exp(-K * y) * erfc((a - y) / b) + exp(K * y) * erfc((a + y) / b);
        const DD closed = DD(G0) / (ldexp(DD(D) * K, 1)) * S * sqrt(ldexp(pi(), -1)) * exp(ldexp(K * K * S * S, -1)) *
                          bracket;
        const MinorityValue u = minority_excess(p, x);
        REQUIRE(u.ok);
        CHECK(rel(u.value, closed) <= 1e-25);
    }
}

TEST_CASE("minority diffusion: collocation against the integral (the second method)",
          "[minority_diffusion][.tvr4_reference_release]") {
    // TVB-15's shape: a Gaussian 10 um wide in a 100 um bar, a diffusion length of 30 um.
    const MinorityProblem p{.diffusivity_cm2_s = 30.0, .lifetime_s = 3e-7, .length_cm = 1e-2,
                            .peak_generation = 1e20, .centre_cm = 5e-3, .width_cm = 1e-3};
    const std::vector<double> xs{0.0, 1e-3, 3.3e-3, 5e-3, 6.1e-3, 9.9e-3, 1e-2};
    const MinorityCollocation col = minority_collocation(p, xs, 8, 40);
    REQUIRE(col.ok);
    const double peak = minority_excess(p, 5e-3).value.hi;
    CHECK(col.tail_estimate > 0.0);  // computed, not a default
    CHECK(col.tail_estimate <= 1e-24 * peak);
    for (std::size_t j = 0; j < xs.size(); ++j) {
        CAPTURE(xs[j]);
        const MinorityValue u = minority_excess(p, xs[j]);
        REQUIRE(u.ok);
        CHECK(std::abs((col.values[j] - u.value).hi) <= 1e-22 * peak);
    }
    REQUIRE_THROWS_AS(minority_collocation(p, xs, 0, 40), std::invalid_argument);
    REQUIRE_THROWS_AS(minority_collocation(p, xs, 8, 3), std::invalid_argument);
    REQUIRE_THROWS_AS(minority_collocation(p, std::vector<double>{2e-2}, 8, 40), std::invalid_argument);
    MinorityProblem bad = p;
    bad.width_cm = 0.0;
    REQUIRE_THROWS_AS(minority_excess(bad, 1e-3), std::invalid_argument);
}
