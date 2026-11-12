// L0 tests of the linear-transient reference (verify/reference/linear_transient.hpp; DECISIONS Unit 55 16.18 (h),
// TVB-19).
//   - tau_d = eps / (q mu N) to rounding; Dormand-Prince against e^{-t / tau_d} (to 1e-12 over five tau_d);
//   - the Haynes-Shockley moments: no recombination (tau infinite) against the uniform age distribution's moments
//     (to 1e-28 relative); the second method (Gauss-Kronrod) to 1e-24; after the source ends, d area / dt =
//     -area / tau, d centre / dt = v and d width^2 / dt = 2 D (dd central differences: 1e-12, the area + (dt / tau)^2); a short pulse's
//     limit; refusals.
#include "NiTCAD/verify/reference/linear_transient.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
#include <numbers>

namespace ref = NiTCAD::verify::reference;
using ref::DD;

namespace {

ref::PulseInput pulse(double tau) {
    return {.diffusivity_cm2_s = 12.0,
            .velocity_cm_s = 4.0e5,
            .lifetime_s = tau,
            .centre_cm = 0.01,
            .sigma_cm = 5e-4,
            .rate_cm3_s = 1e20,
            .t_on_s = 1e-8,
            .t_off_s = 3e-8};
}

double rel(const DD& a, const DD& b) { return std::abs(((a - b) / b).to_double()); }

}  // namespace

TEST_CASE("linear_transient: dielectric relaxation", "[linear_transient]") {
    const auto t = ref::dielectric_relaxation_time(1.0359e-12, 1.602176634e-19, 1350.0, 1e14);
    REQUIRE(t.has_value());
    CHECK(rel(*t, DD(1.0359e-12) / (DD(1.602176634e-19) * DD(1350.0) * DD(1e14))) == 0.0);
    const auto d = ref::relaxation_second_method(*t, 5.0 * t->to_double());
    REQUIRE(d.has_value());
    CHECK(*d <= 1e-12);  // 4000 fixed steps of 1.25e-3 tau_d: far below any criterion
    CHECK_FALSE(ref::dielectric_relaxation_time(0.0, 1.0, 1.0, 1.0).has_value());
    CHECK_FALSE(ref::dielectric_relaxation_time(1.0, 1.0, -1.0, 1.0).has_value());
}

TEST_CASE("linear_transient: Haynes-Shockley moments without recombination", "[linear_transient]") {
    const auto in = pulse(std::numeric_limits<double>::infinity());
    const double t = 1e-7;
    const auto m = ref::haynes_shockley_moments(in, t);
    REQUIRE(m.has_value());
    const DD span = DD(in.t_off_s) - DD(in.t_on_s);
    const DD mid = DD(t) - ldexp(DD(in.t_off_s) + DD(in.t_on_s), -1);
    const DD amp = DD(in.rate_cm3_s) * DD(in.sigma_cm) * sqrt(ldexp(ref::pi(), 1));
    CHECK(rel(m->area_cm2, amp * span) <= 1e-28);
    CHECK(rel(m->centre_cm, DD(in.centre_cm) + DD(in.velocity_cm_s) * mid) <= 1e-28);
    const DD w2 = DD(in.sigma_cm) * DD(in.sigma_cm) + DD(2.0 * in.diffusivity_cm2_s) * mid +
                  DD(in.velocity_cm_s) * DD(in.velocity_cm_s) * span * span / DD(12.0);
    CHECK(rel(m->width_cm * m->width_cm, w2) <= 1e-28);
}

TEST_CASE("linear_transient: Haynes-Shockley moments with recombination", "[linear_transient]") {
    for (const double tau : {1e-8, 2e-7, 1e-5}) {
        const auto in = pulse(tau);
        const double t = 8e-8;
        const auto m = ref::haynes_shockley_moments(in, t);
        REQUIRE(m.has_value());
        CAPTURE(tau, m->area_difference, m->centre_difference, m->width_difference);
        CHECK(m->area_difference <= 1e-24);
        CHECK(m->centre_difference <= 1e-24);
        CHECK(m->width_difference <= 1e-24);
        // The moments' equations after the source ends, by central differences (dt = 1e-12 s).
        const double dt = 1e-12;
        const auto p = ref::haynes_shockley_moments(in, t + dt), q = ref::haynes_shockley_moments(in, t - dt);
        REQUIRE((p && q));
        const DD h = DD(t + dt) - DD(t - dt);  // the doubles actually used
        const DD da = (p->area_cm2 - q->area_cm2) / h, dc = (p->centre_cm - q->centre_cm) / h;
        const DD dw = (p->width_cm * p->width_cm - q->width_cm * q->width_cm) / h;
        CHECK(rel(da, -m->area_cm2 / DD(tau)) <= 1e-12 + (dt / tau) * (dt / tau));  // + the O(dt^2) truncation
        CHECK(rel(dc, DD(in.velocity_cm_s)) <= 1e-12);
        CHECK(rel(dw, DD(2.0 * in.diffusivity_cm2_s)) <= 1e-12);
    }
}

TEST_CASE("linear_transient: a short pulse is a released Gaussian", "[linear_transient]") {
    auto in = pulse(2e-7);
    in.t_off_s = in.t_on_s + 1e-15;
    const double t = 5e-8, age = t - in.t_off_s + 0.5e-15;  // the mean age
    const auto m = ref::haynes_shockley_moments(in, t);
    REQUIRE(m.has_value());
    const double area0 = in.rate_cm3_s * in.sigma_cm * std::sqrt(2.0 * std::numbers::pi) * 1e-15 * std::exp(-age / in.lifetime_s);
    CHECK(std::abs(m->area_cm2.to_double() / area0 - 1.0) <= 1e-7);
    CHECK(std::abs(m->centre_cm.to_double() - (in.centre_cm + in.velocity_cm_s * age)) <= 1e-12);
    CHECK(std::abs(m->width_cm.to_double() / std::sqrt(in.sigma_cm * in.sigma_cm + 2.0 * in.diffusivity_cm2_s * age) -
                   1.0) <= 1e-9);
    CHECK_FALSE(ref::haynes_shockley_moments(in, in.t_on_s).has_value());  // before the end
    auto bad = pulse(1e-7);
    bad.sigma_cm = 0.0;
    CHECK_FALSE(ref::haynes_shockley_moments(bad, 1e-7).has_value());
}
