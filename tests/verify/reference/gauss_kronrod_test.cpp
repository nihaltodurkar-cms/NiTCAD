// TC3 (DECISIONS Unit 55 section 11.3): L0 tests of the Gauss-Kronrod 7/15 quadrature.
//   - the computed nodes (1 u_dd) and weights (4 u_dd) against the 33 digits netlib's dqk15.f prints
//     (Fullerton, 80-digit arithmetic), compared exactly; the owner ruled (2026-11-04, after measurement)
//     1 u_dd for the nodes, as for the TC1 constants, and 4 u_dd for the weights, which come from an
//     8 x 8 moment solve (measured max 2.41);
//   - the 15-point rule integrates x^k on [-1, 1] to 4 u_dd for k <= 22 (its degree 3n + 1, n = 7) and misses
//     by more than 4 u_dd at k = 24 (the degree is real);
//   - adaptive: on ten integrals with closed forms (smooth, peaked, kinked, endpoint-singular-derivative,
//     oscillatory, boundary layer) the true error is <= the reported estimate in every case and the estimate
//     is <= the requested tolerance.
#include "NiTCAD/verify/reference/gauss_kronrod.hpp"

#include "exact_integer.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <functional>
#include <string>

namespace ref = NiTCAD::verify::reference;
using ref::DD;
using ref::u_dd;

TEST_CASE("TC3: computed nodes (1 u_dd) and weights (4 u_dd) match dqk15.f's printed digits") {
    const ref::KronrodRule15& r = ref::kronrod15();
    // netlib quadpack/dqk15.f data statements (read 2026-11-04), digit groups joined.
    const char* xgk[8] = {"0.991455371120812639206854697526329", "0.949107912342758524526189684047851",
                          "0.864864423359769072789712788640926", "0.741531185599394439863864773280788",
                          "0.586087235467691130294144838258730", "0.405845151377397166906606412076961",
                          "0.207784955007898467600689403773245", "0.000000000000000000000000000000000"};
    const char* wgk[8] = {"0.022935322010529224963732008058970", "0.063092092629978553290700663189204",
                          "0.104790010322250183839876322541518", "0.140653259715525918745189590510238",
                          "0.169004726639267902826583426598550", "0.190350578064785409913256402421014",
                          "0.204432940075298892414161999234649", "0.209482141084727828012999174891714"};
    const char* wg[4] = {"0.129484966168869693270611432679082", "0.279705391489276667901467771423780",
                         "0.381830050505118944950369775488975", "0.417959183673469387755102040816327"};
    for (int i = 0; i < 8; ++i) {
        const auto cx = exact::compare_with_literal(r.xgk[static_cast<std::size_t>(i)], xgk[i], 1);
        const auto cw = exact::compare_with_literal(r.wgk[static_cast<std::size_t>(i)], wgk[i], 4);
        INFO("xgk(" << i + 1 << ") " << cx.measured_udd << " u_dd, wgk(" << i + 1 << ") " << cw.measured_udd << " u_dd");
        CHECK(cx.within);
        CHECK(cw.within);
    }
    for (int i = 0; i < 4; ++i) {
        const auto c = exact::compare_with_literal(r.wg[static_cast<std::size_t>(i)], wg[i], 4);
        INFO("wg(" << i + 1 << ") " << c.measured_udd << " u_dd");
        CHECK(c.within);
    }
}

TEST_CASE("TC3: the 15-point rule is exact to degree 22 and not 24") {
    for (int k = 0; k <= 24; k += 2) {
        const auto res = ref::qk15(
            [k](const DD& x) {
                DD p(1.0);
                for (int i = 0; i < k; ++i) p *= x;
                return p;
            },
            DD(-1.0), DD(1.0));
        const DD exact = DD(2.0) / DD(k + 1.0);
        const double err = std::abs((res.value - exact).to_double()) / exact.hi / u_dd;
        INFO("k = " << k << ": " << err << " u_dd");
        if (k <= 22)
            CHECK(err <= 4.0);
        else
            CHECK(err > 4.0);
    }
    // Odd monomials vanish by the rule's symmetry.
    for (int k = 1; k <= 23; k += 2) {
        const auto res = ref::qk15(
            [k](const DD& x) {
                DD p(1.0);
                for (int i = 0; i < k; ++i) p *= x;
                return p;
            },
            DD(-1.0), DD(1.0));
        CHECK(res.value.hi == 0.0);
    }
}

namespace {

// T_k(x) by T_{j+1} = 2 x T_j - T_{j-1}.
DD chebyshev(int k, const DD& x) {
    DD t0(1.0), t1 = x;
    if (k == 0) return t0;
    for (int j = 1; j < k; ++j) {
        const DD t2 = ldexp(x * t1, 1) - t0;
        t0 = t1;
        t1 = t2;
    }
    return t1;
}

struct Integral {
    std::string name;
    std::function<DD(const DD&)> f;
    DD a, b, exact;
};

}  // namespace

TEST_CASE("TC3: adaptive 7/15 on ten closed-form integrals: true error <= estimate <= tolerance") {
    const double eps_small = 1e-6;   // boundary peak
    const double c = 0.3, eps_kink = 1e-4;
    const DD one(1.0);
    const Integral cases[] = {
        {"e^x on [0,1]", [](const DD& x) { return ref::exp(x); }, DD(0.0), one, ref::e() - one},
        {"1/(1+x) on [0,1]", [one](const DD& x) { return one / (one + x); }, DD(0.0), one, ref::ln2()},
        {"1/(1+x^2) on [0,1]", [one](const DD& x) { return one / (one + x * x); }, DD(0.0), one, ldexp(ref::pi(), -2)},
        {"sqrt x on [0,1]", [](const DD& x) { return ref::sqrt(x); }, DD(0.0), one, DD(2.0) / DD(3.0)},
        {"log(1+x) on [0,1]", [one](const DD& x) { return ref::log(one + x); }, DD(0.0), one,
         ldexp(ref::ln2(), 1) - one},
        {"1/(x+1e-6) on [0,1]", [eps_small](const DD& x) { return DD(1.0) / (x + DD(eps_small)); }, DD(0.0), one,
         ref::log(DD(1.0) + DD(eps_small)) - ref::log(DD(eps_small))},
        {"1/(|x-0.3|+1e-4) on [0,1]",
         [c, eps_kink](const DD& x) { return DD(1.0) / (abs(x - DD(c)) + DD(eps_kink)); }, DD(0.0), one,
         ref::log((DD(c) + DD(eps_kink)) / DD(eps_kink)) + ref::log((one - DD(c) + DD(eps_kink)) / DD(eps_kink))},
        {"T_40 on [-1,1]", [](const DD& x) { return chebyshev(40, x); }, DD(-1.0), one, DD(-2.0) / DD(1599.0)},
        {"T_100 on [-1,1]", [](const DD& x) { return chebyshev(100, x); }, DD(-1.0), one, DD(-2.0) / DD(9999.0)},
        {"e^(-50x) on [0,1]", [](const DD& x) { return ref::exp(DD(-50.0) * x); }, DD(0.0), one,
         -ref::expm1(DD(-50.0)) / DD(50.0)},
    };
    const double epsrel = 1e-24;
    for (const auto& ic : cases) {
        const auto res = ref::qag15(ic.f, ic.a, ic.b, 0.0, epsrel);
        const double true_err = std::abs((res.value - ic.exact).to_double());
        const double tol = epsrel * std::abs(res.value.to_double());
        INFO(ic.name << ": true error " << true_err << ", estimate " << res.abs_error << ", tolerance " << tol
                     << ", intervals " << res.intervals << ", status " << static_cast<int>(res.status));
        CHECK(res.status == ref::QuadratureStatus::ok);
        CHECK(true_err <= res.abs_error);
        CHECK(res.abs_error <= tol);
    }
}

TEST_CASE("TC3: invalid input is refused") {
    const auto res = ref::qag15([](const DD& x) { return x; }, DD(0.0), DD(1.0), 0.0, 1e-40);
    CHECK(res.status == ref::QuadratureStatus::invalid_input);
}
