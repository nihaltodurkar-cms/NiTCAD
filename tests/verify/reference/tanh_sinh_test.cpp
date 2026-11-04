// TC4 (DECISIONS Unit 55 section 11.3): L0 tests of the tanh-sinh quadrature.
//   - the four frozen integrals: int_0^1 x^-1/2 = 2, int_0^1 ln x = -1, int_0^1 x^-1/2 ln x = -4 (endpoint
//     singularities) and int_0^inf e^-x^2 = sqrt(pi)/2 (mapped to (0, 1] by x = 1/s - 1, as BJL map their
//     problems 11-13), each to 16 u_dd relative, with the error estimate bounding the true error (the
//     estimate is the sources' heuristic with the owner's 8 u_dd int |f| rounding floor and eps2 = eps1^3,
//     both ruled 2026-11-04 after measurement; see the header);
//   - two smooth integrals the same way, and the estimator's published edge rules (E = 1 for n <= 2).
#include "NiTCAD/verify/reference/tanh_sinh.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <functional>
#include <string>

namespace ref = NiTCAD::verify::reference;
using ref::DD;
using ref::u_dd;

namespace {

struct Integral {
    std::string name;
    std::function<DD(const DD&)> f;
    DD a, b, exact;
};

void check_integral(const Integral& ic, double tol_udd) {
    const double epsrel = tol_udd * u_dd;
    const auto res = ref::tanh_sinh(ic.f, ic.a, ic.b, 0.0, epsrel);
    const double true_err = std::abs((res.value - ic.exact).to_double());
    const double measured = true_err / std::abs(ic.exact.hi) / u_dd;
    INFO(ic.name << ": true error " << true_err << " (" << measured << " u_dd), estimate " << res.abs_error
                 << ", levels " << res.levels << ", evaluations " << res.evaluations << ", status "
                 << static_cast<int>(res.status));
    CHECK(res.status == ref::TanhSinhStatus::ok);
    CHECK(measured <= tol_udd);
    CHECK(true_err <= res.abs_error);
}

}  // namespace

TEST_CASE("TC4: the four frozen integrals to 16 u_dd, the estimate bounding the true error") {
    const DD one(1.0);
    const Integral cases[] = {
        {"x^-1/2 on [0,1]", [one](const DD& x) { return one / ref::sqrt(x); }, DD(0.0), one, DD(2.0)},
        {"ln x on [0,1]", [](const DD& x) { return ref::log(x); }, DD(0.0), one, DD(-1.0)},
        {"x^-1/2 ln x on [0,1]", [](const DD& x) { return ref::log(x) / ref::sqrt(x); }, DD(0.0), one, DD(-4.0)},
        {"e^-x^2 on [0,inf), x = 1/s - 1",
         [one](const DD& s) {
             const DD x = one / s - one;
             return ref::exp(-(x * x)) / (s * s);
         },
         DD(0.0), one, ldexp(ref::sqrt(ref::pi()), -1)},
    };
    for (const auto& ic : cases) check_integral(ic, 16.0);
}

TEST_CASE("TC4: smooth integrands the same way") {
    const DD one(1.0);
    check_integral({"1 on [-1,1]", [one](const DD&) { return one; }, DD(-1.0), one, DD(2.0)}, 16.0);
    check_integral({"e^x on [0,1]", [](const DD& x) { return ref::exp(x); }, DD(0.0), one, ref::e() - one}, 16.0);
}

TEST_CASE("TC4: a singularity the table cannot resolve is reported by the endpoint term, not hidden") {
    // int_0^1 x^-0.9 = 10: at the last abscissa (r y ~ 1e-96) |w f| ~ 1e-9 and the truncated tail is
    // 10 (r y)^0.1 ~ 1e-9, so no level can reach 1e-20; the estimate must say so and still bound the error.
    const auto res = ref::tanh_sinh([](const DD& x) { return ref::pow(x, DD(-0.9)); }, DD(0.0), DD(1.0), 0.0, 1e-20);
    const double true_err = std::abs((res.value - DD(10.0)).to_double());
    INFO("true error " << true_err << ", estimate " << res.abs_error << ", levels " << res.levels);
    CHECK(res.status == ref::TanhSinhStatus::limit_reached);
    CHECK(true_err <= res.abs_error);
}

TEST_CASE("TC4: the estimator's published rules and invalid input") {
    const auto two = ref::tanh_sinh([](const DD& x) { return ref::exp(x); }, DD(0.0), DD(1.0), 0.0, 1e-3, 2);
    CHECK(two.status == ref::TanhSinhStatus::limit_reached);
    CHECK(two.abs_error == 1.0);  // E_n = 1 for n <= 2
    // The first two levels are never trusted, however loose the tolerance.
    const auto loose = ref::tanh_sinh([](const DD& x) { return ref::exp(x); }, DD(0.0), DD(1.0), 2.0, 0.0);
    CHECK(loose.status == ref::TanhSinhStatus::ok);
    CHECK(loose.levels == 3);
    const auto bad = ref::tanh_sinh([](const DD& x) { return x; }, DD(0.0), DD(1.0), 0.0, 0.0);
    CHECK(bad.status == ref::TanhSinhStatus::invalid_input);
    const auto nan = ref::tanh_sinh([](const DD&) { return DD(std::nan("")); }, DD(0.0), DD(1.0), 0.0, 1e-20);
    CHECK(nan.status == ref::TanhSinhStatus::bad_integrand);
    const ref::TanhSinhTable& tab = ref::tanh_sinh_table();
    CHECK(tab.y[0] == DD(1.0));
    CHECK(tab.w[0] == ldexp(ref::pi(), -1));
}
