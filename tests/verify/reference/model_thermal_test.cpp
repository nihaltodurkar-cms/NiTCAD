// TVR-3 step 8 (DECISIONS Unit 55 15.9, 15.18): L0 tests of the references of TVB-31.26-31.29
// (verify/reference/model_thermal.hpp).
//   - Emission: (kT / h)(2 / N)^(1/3) against sqrt(kT / (2 pi m)) with the mass N implies; dv/dT against complex step
//     (TC8) of v(T) through the power-law N.
//   - MLDA: the mean's series against its closed form at x = 1/4; the mean against tanh-sinh (TC4) of f over [0, a].
//   - Conductivity: d kappa / dT against a central difference.
//   - Fermi orders: F_0 and F_1 differences against tanh-sinh of F_-1 and F_0; the series branch against the closed
//     forms at x = -1; F_1 continuous at 0 across the reflection; dh/dx against central differences; h -> 1/2.
#include "NiTCAD/verify/reference/model_thermal.hpp"

#include "NiTCAD/verify/reference/complex_step.hpp"
#include "NiTCAD/verify/reference/tanh_sinh.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <complex>

namespace ref = NiTCAD::verify::reference;
namespace mh = NiTCAD::verify::reference::model_thermal;
using ref::complex_d;
using ref::DD;

TEST_CASE("model_thermal: emission velocity", "[verify][model_thermal]") {
    namespace b = NiTCAD::base;
    for (const double N : {1e17, 2.86e19, 1e21}) {
        for (const double T : {200.0, 300.0, 500.0}) {
            const DD v = mh::emission_dos({DD(N), DD(T)});
            // N [m^-3] = 2 (2 pi m k T / h^2)^(3/2): m = (N / 2)^(2/3) h^2 / (2 pi k T).
            const DD h = DD(2.0) * ref::pi() * DD(b::hbar_J_s), kT = DD(b::k_B_J_per_K) * DD(T);
            const DD Nm = DD(N) * DD(1e6);
            const DD m = exp(log(Nm / DD(2.0)) * DD(2.0) / DD(3.0)) * h * h / (DD(2.0) * ref::pi() * kT);
            const DD vm = sqrt(kT / (DD(2.0) * ref::pi() * m)) * DD(100.0);
            CHECK(std::abs(((v - vm) / vm).hi) <= 1e-28);
        }
    }
    // dv/dT through N = N300 (T / 300)^1.5, both forms.
    for (const double A : {0.0, 120.0}) {
        for (const double T : {200.0, 300.0, 450.0}) {
            auto g = [&](complex_d t) {
                const complex_d N = 2.86e19 * std::pow(t / 300.0, 1.5);
                if (A > 0.0) return A * t * t / (NiTCAD::base::q_C * N);
                const double h = 2.0 * 3.14159265358979323846 * NiTCAD::base::hbar_J_s;
                return NiTCAD::base::k_B_J_per_K * t / h * std::pow(2.0 / N, 1.0 / 3.0);
            };
            const auto cs = ref::complex_step(g, T, 1e-30 * T);
            REQUIRE(cs.status == ref::ComplexStepStatus::ok);
            const DD N = DD(2.86e19) * exp(DD(1.5) * log(DD(T) / DD(300.0)));
            const auto e = mh::emission(DD(A), DD(T), N, DD(1.5) / DD(T));
            INFO("A " << A << " T " << T);
            CHECK(std::abs(e[1].hi - cs.derivative) <= 1e-14 * std::abs(cs.derivative));
            CHECK(std::abs(((e[1] - e[0] / (DD(2.0) * DD(T))) / e[1]).hi) <= 1e-28);
        }
    }
}

TEST_CASE("model_thermal: MLDA", "[verify][model_thermal]") {
    const DD lam(1.0);
    const DD below = mh::mlda_mean({DD(std::nextafter(0.25, 0.0)), lam});
    const DD r(0.25);
    const DD closed = DD(1.0) - sqrt(ref::pi()) / (DD(2.0) * r) * ref::erf(r);
    CHECK(std::abs(((below - closed) / closed).hi) <= 1e-15);  // nextafter moves the value by 2 eps relative
    for (const double a : {1e-3, 0.1, 0.2499, 0.25, 0.7, 2.0, 9.0}) {
        const auto t = ref::tanh_sinh([&](const DD& z) { return mh::mlda_factor({z, lam}); }, DD(0.0), DD(a), 0.0,
                                      1e-30);
        REQUIRE(t.status == ref::TanhSinhStatus::ok);
        const DD q = t.value / DD(a);
        const DD m = mh::mlda_mean({DD(a), lam});
        INFO("a " << a << ": " << m.hi << " vs " << q.hi);
        // The closed form loses u_dd / x^2 to cancellation at small x.
        CHECK(std::abs(((m - q) / q).hi) <= 1e-28 + 4.0 * ref::u_dd / (a * a));
    }
}

TEST_CASE("model_thermal: conductivity", "[verify][model_thermal]") {
    for (const double T : {200.0, 300.0, 600.0}) {
        const auto k = mh::conductivity({DD(1.48), DD(1.3), DD(T)});
        const DD h(T * 0x1p-30);
        const DD cd = (mh::conductivity({DD(1.48), DD(1.3), DD(T) + h})[0] -
                       mh::conductivity({DD(1.48), DD(1.3), DD(T) - h})[0]) /
                      (DD(2.0) * h);
        CHECK(std::abs(((k[1] - cd) / k[1]).hi) <= 1e-15);
    }
}

TEST_CASE("model_thermal: Fermi orders and h", "[verify][model_thermal]") {
    // Differences against quadrature of the next lower order.
    for (const std::array<double, 2> ab : {std::array<double, 2>{-5.0, -1.0}, std::array<double, 2>{-1.0, 0.0},
                                           std::array<double, 2>{-0.3, 1.0}, std::array<double, 2>{1.0, 5.0},
                                           std::array<double, 2>{5.0, 30.0}}) {
        const auto lo = mh::fermi_orders(DD(ab[0])), hi = mh::fermi_orders(DD(ab[1]));
        const auto q0 = ref::tanh_sinh([](const DD& x) { return mh::fermi_orders(x).minus_one; }, DD(ab[0]),
                                       DD(ab[1]), 0.0, 1e-30);
        const auto q1 = ref::tanh_sinh([](const DD& x) { return mh::fermi_orders(x).zero; }, DD(ab[0]), DD(ab[1]),
                                       0.0, 1e-30);
        REQUIRE(q0.status == ref::TanhSinhStatus::ok);
        REQUIRE(q1.status == ref::TanhSinhStatus::ok);
        INFO("[" << ab[0] << ", " << ab[1] << "]");
        CHECK(std::abs((hi.zero - lo.zero - q0.value).hi) <= 1e-29 * (std::abs(hi.zero.hi) + std::abs(lo.zero.hi)));
        CHECK(std::abs((hi.one - lo.one - q1.value).hi) <= 1e-29 * (std::abs(hi.one.hi) + std::abs(lo.one.hi)));
    }
    // The series branch at x = -1 against the closed forms.
    const DD z = exp(DD(-1.0));
    const auto s = mh::fermi_orders(DD(-1.0));
    CHECK(std::abs(((s.zero - ref::log1p_exp(DD(-1.0))) / s.zero).hi) <= 1e-30);
    CHECK(std::abs(((s.one + ref::li2(-z)) / s.one).hi) <= 1e-30);
    CHECK(std::abs(((s.minus_one - z / (DD(1.0) + z)) / s.minus_one).hi) <= 1e-30);
    // F_1 continuous at 0 across the reflection.
    const auto a = mh::fermi_orders(DD(-1e-12)), c = mh::fermi_orders(DD(1e-12));
    CHECK(std::abs((c.one - a.one - DD(2e-12) * mh::fermi_orders(DD(0.0)).zero).hi) <= 1e-28);
    // dh/dx against central differences; h -> 1/2 and dh/dx -> 0 as x -> -inf.
    for (const double x : {-700.0, -40.0, -3.0, -1.0, -0.5, 0.0, 0.7, 1.0, 3.0, 40.0, 700.0}) {
        const DD h(std::max(1.0, std::abs(x)) * 0x1p-30);
        const DD cd = (mh::thermal_diffusion({DD(x) + h})[0] - mh::thermal_diffusion({DD(x) - h})[0]) / (DD(2.0) * h);
        const auto t = mh::thermal_diffusion({DD(x)});
        INFO("x " << x << ": " << t[1].hi << " vs " << cd.hi);
        CHECK(std::abs((t[1] - cd).hi) <= 1e-15 * (1.0 + std::abs(x)));
    }
    CHECK(std::abs((mh::thermal_diffusion({DD(-740.0)})[0] - DD(0.5)).hi) <= 1e-30);
    CHECK(std::abs(mh::thermal_diffusion({DD(-740.0)})[1].hi) <= 1e-30);
    // The cancellation-free dh/dx against the direct 2 - 2 F_1 F_-1 / F_0^2 where that is still accurate (x = -3,
    // dh/dx ~ 0.02), and against its leading term e^x / 2 where it is not (x = -40, relative O(e^x)).
    {
        const auto f = mh::fermi_orders(DD(-3.0));
        const DD direct = DD(2.0) - DD(2.0) * f.one * f.minus_one / (f.zero * f.zero);
        CHECK(std::abs((mh::thermal_diffusion({DD(-3.0)})[1] - direct).hi) <= 1e-29);
        const DD z40 = exp(DD(-40.0));
        CHECK(std::abs(((mh::thermal_diffusion({DD(-40.0)})[1] - z40 / DD(2.0)) / (z40 / DD(2.0))).hi) <= 1e-16);
        const DD z7 = exp(DD(-700.0));  // normal, but its low part subnormal (about 1e-19 relative)
        CHECK(std::abs(((mh::thermal_diffusion({DD(-700.0)})[1] - z7 / DD(2.0)) / (z7 / DD(2.0))).hi) <= 1e-17);
    }
}
