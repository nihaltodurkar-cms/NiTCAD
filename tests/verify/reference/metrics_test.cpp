// Metrics (DECISIONS Unit 55 section 11.5, rulings 11.13): L0 tests.
//   - exact recovery of p from synthetic e_k = C h_k^p, and from synthetic solutions f_k = f* + C h_k^p, to 1e-12;
//   - the NPARC worked example (spatconv.html, Table of Pr on three grids) to half a unit of each printed digit:
//     p = 1.786170, Pr_h=0 = 0.97050 + 0.00080 = 0.97130, GCI12 = 0.103083 %, GCI23 = 0.356249 %; the asymptotic
//     ratio 1.002019 from the page's rounded intermediates, and exactly f1 / f2 = 1.0020237 from the values;
//   - the asymptotic-range check needs 4 grids (owner ruling), with its order spread and GCI ratio;
//   - eps_int = 0 for an additive e(h, tol) and the exact cross term for a multiplicative one;
//   - GCI refused (not_evaluated) for non-monotone solutions, p <= 0, non-constant r, a zero solution.
#include "NiTCAD/verify/reference/metrics.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cfloat>
#include <cmath>
#include <vector>

namespace ref = NiTCAD::verify::reference;
using ref::MetricStatus;

namespace {

constexpr double eps = DBL_EPSILON;

bool within_half_unit(double value, double printed, double unit) { return std::abs(value - printed) <= 0.5 * unit; }

}  // namespace

TEST_CASE("observed order recovered exactly from synthetic errors", "[verify][metrics]") {
    for (const double p : {1.0, 1.5, 2.0, 4.0}) {
        for (const double r : {2.0, 1.5, std::sqrt(2.0), 1.1}) {
            std::vector<double> h, e;
            for (int k = 0; k < 5; ++k) {
                h.push_back(1e-3 * std::pow(r, k));
                e.push_back(0.37 * std::pow(h.back(), p));
            }
            const auto orders = ref::observed_orders_from_errors(h, e);
            REQUIRE(orders.size() == 4);
            for (const auto& m : orders) {
                REQUIRE(m.status == MetricStatus::ok);
                CHECK(std::abs(m.value - p) <= 1e-12);
            }
        }
    }
    // Non-constant ratios are allowed when the errors are exact.
    const std::array<double, 4> h{1.0, 1.3, 2.9, 3.1};
    std::array<double, 4> e{};
    for (int k = 0; k < 4; ++k) e[k] = 2.5 * std::pow(h[k], 2.0);
    for (const auto& m : ref::observed_orders_from_errors(h, e)) CHECK(std::abs(m.value - 2.0) <= 1e-12);
}

TEST_CASE("observed order recovered exactly from synthetic solutions", "[verify][metrics]") {
    for (const double p : {1.0, 1.5, 2.0, 3.0}) {
        for (const double r : {2.0, 1.5}) {
            // f_k = f* + C h_k^p with differences of order 0.1 relative to f* = 1 (no cancellation to speak of).
            const double c = 0.1;
            const double f1 = 1.0 + c, f2 = 1.0 + c * std::pow(r, p), f3 = 1.0 + c * std::pow(r, 2.0 * p);
            const auto m = ref::observed_order_from_solutions(f1, f2, f3, r);
            REQUIRE(m.status == MetricStatus::ok);
            CHECK(std::abs(m.value - p) <= 1e-12);
            const auto x = ref::richardson(f1, f2, r, p);
            REQUIRE(x.status == MetricStatus::ok);
            CHECK(std::abs(x.value - 1.0) <= 1e-12);
        }
    }
}

TEST_CASE("NPARC worked example reproduced to its printed digits", "[verify][metrics]") {
    // spatconv.html: grids 1, 2, 4 (normalised spacing), Pr = 0.97050, 0.96854, 0.96178.
    const std::array<double, 3> h{1.0, 2.0, 4.0};
    const std::array<double, 3> f{0.97050, 0.96854, 0.96178};
    const auto s = ref::grid_convergence(h, f);
    REQUIRE(s.status == MetricStatus::ok);
    CHECK(s.r == 2.0);
    REQUIRE(s.orders.size() == 1);
    CHECK(within_half_unit(s.orders[0], 1.786170, 1e-6));
    const double increment = (f[0] - f[1]) / (std::pow(2.0, s.orders[0]) - 1.0);
    CHECK(within_half_unit(increment, 0.00080, 1e-5));
    CHECK(within_half_unit(s.extrapolated, 0.97130, 1e-5));
    REQUIRE(s.gci.size() == 2);
    CHECK(within_half_unit(100.0 * s.gci[0], 0.103083, 1e-6));
    CHECK(within_half_unit(100.0 * s.gci[1], 0.356249, 1e-6));
    // The page's ratio is the quotient of its rounded intermediates; the exact ratio is f1 / f2 identically.
    CHECK(within_half_unit(0.356249 / (std::pow(2.0, 1.786170) * 0.103083), 1.002019, 1e-6));
    CHECK(std::abs(s.three_grid_ratio - f[0] / f[1]) <= 4.0 * eps);
    CHECK(within_half_unit(s.three_grid_ratio, 1.0020237, 1e-7));
    // Three grids: the asymptotic-range check is not evaluated.
    CHECK(s.asymptotic_status == MetricStatus::not_evaluated);
    CHECK(s.asymptotic_reason == "needs 4 grids");
    CHECK_FALSE(s.in_asymptotic_range);
    // The stand-alone forms agree with the study.
    CHECK(ref::observed_order_from_solutions(f[0], f[1], f[2], 2.0).value == s.orders[0]);
    CHECK(ref::gci(f[0], f[1], 2.0, s.orders[0], 3).value == s.gci[0]);
    CHECK(ref::richardson(f[0], f[1], 2.0, s.orders[0]).value == s.extrapolated);
}

TEST_CASE("GCI safety factors", "[verify][metrics]") {
    const double three = ref::gci(1.0, 1.01, 2.0, 2.0, 3).value;
    const double two = ref::gci(1.0, 1.01, 2.0, 2.0, 2).value;
    CHECK(std::abs(three - 1.25 * 0.01 / 3.0) <= 4.0 * eps * three);
    CHECK(std::abs(two - 3.0 * 0.01 / 3.0) <= 4.0 * eps * two);
    CHECK(ref::gci(1.0, 1.01, 2.0, 2.0, 1).status == MetricStatus::invalid_input);
    CHECK(ref::gci(0.0, 1.01, 2.0, 2.0, 3).status == MetricStatus::not_evaluated);
    CHECK(ref::gci(1.0, 1.01, 1.0, 2.0, 3).status == MetricStatus::invalid_input);
    CHECK(ref::gci(1.0, 1.01, 2.0, 0.0, 3).status == MetricStatus::invalid_input);
}

TEST_CASE("asymptotic-range check on four or more grids", "[verify][metrics]") {
    const std::array<double, 5> h{1.0, 2.0, 4.0, 8.0, 16.0};
    SECTION("pure power law: in range, spread 0, ratio f_{k+1}/f_{k+2}") {
        std::array<double, 5> f{};
        // Dyadic: h powers of two and C = 2^-20, so every f_k and difference is exact.
        for (int k = 0; k < 5; ++k) f[k] = 1.0 + 0x1p-20 * h[k] * h[k];
        const auto s = ref::grid_convergence(h, f);
        REQUIRE(s.status == MetricStatus::ok);
        REQUIRE(s.orders.size() == 3);
        for (double p : s.orders) CHECK(std::abs(p - 2.0) <= 1e-12);
        CHECK(s.asymptotic_status == MetricStatus::ok);
        CHECK(s.order_spread <= 1e-12);
        // rho_k = r^{p_{k+1} - p_k} f_{k+1} / f_{k+2}, k = 0, 1; with equal orders the largest deviation is
        // 1 - f_2 / f_3.
        CHECK(std::abs(s.ratio_deviation - (1.0 - f[2] / f[3])) <= 1e-12);
        CHECK(s.in_asymptotic_range);
    }
    SECTION("an order change of 0.3 between triples is out of range") {
        // Differences d_k = f_{k+1} - f_k growing by r^2 then r^2.3: orders 2 then 2.3.
        const double d0 = 0.001;
        const std::array<double, 4> h4{1.0, 2.0, 4.0, 8.0};
        const std::array<double, 4> f{1.0, 1.0 + d0, 1.0 + d0 + 4.0 * d0,
                                      1.0 + d0 + 4.0 * d0 + 4.0 * d0 * std::pow(2.0, 2.3)};
        const auto s = ref::grid_convergence(h4, f);
        REQUIRE(s.status == MetricStatus::ok);
        CHECK(std::abs(s.orders[0] - 2.0) <= 1e-9);
        CHECK(std::abs(s.orders[1] - 2.3) <= 1e-9);
        CHECK(std::abs(s.order_spread - 0.3) <= 1e-9);
        CHECK(s.ratio_deviation > 0.1);
        CHECK_FALSE(s.in_asymptotic_range);
    }
    SECTION("an order change of 0.05 is in range") {
        const double d0 = 0.001;
        const std::array<double, 4> h4{1.0, 2.0, 4.0, 8.0};
        const std::array<double, 4> f{1.0, 1.0 + d0, 1.0 + d0 + 4.0 * d0,
                                      1.0 + d0 + 4.0 * d0 + 4.0 * d0 * std::pow(2.0, 2.05)};
        const auto s = ref::grid_convergence(h4, f);
        REQUIRE(s.status == MetricStatus::ok);
        CHECK(std::abs(s.order_spread - 0.05) <= 1e-9);
        CHECK(s.ratio_deviation <= 0.1);
        CHECK(s.in_asymptotic_range);
    }
    SECTION("each test can fail alone") {
        // Spread 0.12 > 0.1 while rho = 2^0.12 f_1 / f_2 = 1.083 is inside 10 %: out of range on the spread.
        const double d0 = 0.001;
        const std::array<double, 4> h4{1.0, 2.0, 4.0, 8.0};
        const std::array<double, 4> f{1.0, 1.0 + d0, 1.0 + d0 + 4.0 * d0,
                                      1.0 + d0 + 4.0 * d0 + 4.0 * d0 * std::pow(2.0, 2.12)};
        const auto s = ref::grid_convergence(h4, f);
        REQUIRE(s.status == MetricStatus::ok);
        CHECK(s.order_spread > 0.1);
        CHECK(s.ratio_deviation <= 0.1);
        CHECK_FALSE(s.in_asymptotic_range);
        // Equal orders (spread 0) with large relative differences: rho = f_1 / f_2 = 0.75 is outside 10 %.
        std::array<double, 4> g{};
        for (int k = 0; k < 4; ++k) g[k] = 1.0 + 0x1p-5 * h4[k] * h4[k];  // dyadic, exact
        const auto t = ref::grid_convergence(h4, g);
        REQUIRE(t.status == MetricStatus::ok);
        CHECK(t.order_spread <= 1e-12);
        CHECK(t.ratio_deviation > 0.1);
        CHECK_FALSE(t.in_asymptotic_range);
    }
}

TEST_CASE("grid study refusals", "[verify][metrics]") {
    const std::array<double, 3> h{1.0, 2.0, 4.0};
    SECTION("not monotone") {
        const std::array<double, 3> f{1.0, 1.01, 1.005};
        const auto s = ref::grid_convergence(h, f);
        CHECK(s.status == MetricStatus::not_evaluated);
        CHECK(s.reason == "not monotone");
        CHECK(s.gci.empty());
        const auto m = ref::observed_order_from_solutions(f[0], f[1], f[2], 2.0);
        CHECK(m.status == MetricStatus::not_evaluated);
        CHECK(m.reason == "not monotone");
    }
    SECTION("observed order <= 0 (differences shrinking on coarser grids)") {
        const std::array<double, 3> f{1.0, 1.04, 1.05};
        const auto s = ref::grid_convergence(h, f);
        CHECK(s.status == MetricStatus::not_evaluated);
        CHECK(s.reason == "observed order <= 0");
        CHECK(ref::observed_order_from_solutions(f[0], f[1], f[2], 2.0).status == MetricStatus::not_evaluated);
    }
    SECTION("order exactly 0 (equal differences)") {
        const std::array<double, 3> f{1.0, 1.25, 1.5};
        CHECK(ref::grid_convergence(h, f).reason == "observed order <= 0");
    }
    SECTION("zero difference") {
        const std::array<double, 3> f{1.0, 1.0, 1.5};
        CHECK(ref::grid_convergence(h, f).reason == "zero difference");
    }
    SECTION("non-constant r") {
        const std::array<double, 3> hv{1.0, 2.0, 4.0 * (1.0 + 1e-10)};
        const std::array<double, 3> f{1.0, 1.01, 1.05};
        const auto s = ref::grid_convergence(hv, f);
        CHECK(s.status == MetricStatus::not_evaluated);
        CHECK(s.reason == "non-constant r");
        // Within 1e-12 relative is constant.
        const std::array<double, 3> hc{1.0, 2.0, 4.0 * (1.0 + 1e-14)};
        CHECK(ref::grid_convergence(hc, f).status == MetricStatus::ok);
    }
    SECTION("a zero solution") {
        const std::array<double, 3> f{0.0, 0.01, 0.05};
        CHECK(ref::grid_convergence(h, f).status == MetricStatus::not_evaluated);
    }
    SECTION("malformed") {
        const std::array<double, 2> h2{1.0, 2.0};
        const std::array<double, 2> f2{1.0, 1.01};
        CHECK(ref::grid_convergence(h2, f2).status == MetricStatus::invalid_input);
        const std::array<double, 3> hd{2.0, 1.0, 0.5};
        const std::array<double, 3> f{1.0, 1.01, 1.05};
        CHECK(ref::grid_convergence(hd, f).status == MetricStatus::invalid_input);
        const std::array<double, 3> fn{1.0, NAN, 1.05};
        CHECK(ref::grid_convergence(h, fn).status == MetricStatus::invalid_input);
    }
}

TEST_CASE("factorial interaction contrast", "[verify][metrics]") {
    // Additive e(h, tol) = a(h) + b(tol): eps_int = 0. Dyadic values make the inputs exact, so 0 exactly.
    const double a_h = 0x1.8p-10, a_h2 = 0x1.8p-12, b_t = 0x1.4p-20, b_t2 = 0x1.4p-27;
    CHECK(ref::interaction_contrast(a_h + b_t, a_h + b_t2, a_h2 + b_t, a_h2 + b_t2).value == 0.0);
    // Additive with non-dyadic values: zero to the inputs' own rounding.
    for (const double s : {1.0, 3.7, 1e-8}) {
        const double ah = 0.3 * s, ah2 = 0.077 * s, bt = 1.3e-5 * s, bt2 = 1.1e-7 * s;
        const double c = ref::interaction_contrast(ah + bt, ah + bt2, ah2 + bt, ah2 + bt2).value;
        CHECK(std::abs(c) <= 4.0 * eps * (ah + bt));
    }
    // Multiplicative e = A(h) B(tol): eps_int = (A(h/2) - A(h)) (B(tol/100) - B(tol)).
    const double ah = 0x1.8p-3, ah2 = 0x1.8p-5, bt = 0x1.2p+0, bt2 = 0x1.01p+0;
    const double c = ref::interaction_contrast(ah * bt, ah * bt2, ah2 * bt, ah2 * bt2).value;
    CHECK(c == (ah2 - ah) * (bt2 - bt));
    const double ah_n = 0.31, ah2_n = 0.078, bt_n = 1.37, bt2_n = 1.0042;
    const double cn = ref::interaction_contrast(ah_n * bt_n, ah_n * bt2_n, ah2_n * bt_n, ah2_n * bt2_n).value;
    const double exact = (ah2_n - ah_n) * (bt2_n - bt_n);
    CHECK(std::abs(cn - exact) <= 4.0 * eps * ah_n * bt_n);
    // The four-term sum is exact before its one rounding: 1 - (-2^-60) - 0 + (-1) = 2^-60, which a running
    // double sum loses (1 + 2^-60 rounds to 1).
    CHECK(ref::interaction_contrast(-1.0, 0.0, -0x1p-60, 1.0).value == 0x1p-60);
}

TEST_CASE("combined uncertainty and the interacting flag", "[verify][metrics]") {
    const std::array<double, 5> u{1e-3, -2e-4, 5e-5, 0.0, 1e-6};
    const auto c = ref::combine_uncertainty(u, -1e-4);
    REQUIRE(c.status == MetricStatus::ok);
    CHECK(std::abs(c.value - (1e-3 + 2e-4 + 5e-5 + 1e-6 + 1e-4)) <= 4.0 * eps * c.value);
    CHECK_FALSE(c.interacting);  // 1e-4 / 1.351e-3 = 7.4 %
    const auto d = ref::combine_uncertainty(u, 2e-4);
    CHECK(d.interacting);  // 2e-4 / 1.451e-3 = 13.8 %
    // The rule is a sum of magnitudes, not a root-sum-square.
    const std::array<double, 2> v{3.0, 4.0};
    CHECK(ref::combine_uncertainty(v, 0.0).value == 7.0);
    const std::array<double, 1> w{NAN};
    CHECK(ref::combine_uncertainty(w, 0.0).status == MetricStatus::invalid_input);
}

TEST_CASE("relative and log errors with floors", "[verify][metrics]") {
    CHECK(ref::relative_error(1.01, 1.0, 0.0).value == std::abs(1.01 - 1.0));
    // Below the floor the floor scales: |1e-15 - 2e-15| / 1e-12.
    CHECK(std::abs(ref::relative_error(1e-15, 2e-15, 1e-12).value - 1e-3) <= 4.0 * eps * 1e-3);
    CHECK(ref::relative_error(-2.0, -1.0, 0.0).value == 1.0);
    CHECK(ref::relative_error(1.0, 0.0, 0.0).status == MetricStatus::not_evaluated);
    CHECK(ref::relative_error(1.0, 1.0, -1.0).status == MetricStatus::invalid_input);
    CHECK(ref::relative_error(NAN, 1.0, 0.0).status == MetricStatus::invalid_input);

    CHECK(std::abs(ref::log_error(1e-9, 1e-12, 0.0, 0.0).value - 3.0) <= 4.0 * eps * 3.0);
    CHECK(std::abs(ref::log_error(1e-12, 1e-9, 0.0, 0.0).value - 3.0) <= 4.0 * eps * 3.0);
    CHECK(ref::log_error(1e-15, 1e-9, 1e-14, 0.0).status == MetricStatus::not_evaluated);
    CHECK(ref::log_error(1e-9, 1e-15, 0.0, 1e-14).status == MetricStatus::not_evaluated);
    CHECK(ref::log_error(1e-14, 1e-9, 1e-14, 0.0).status == MetricStatus::not_evaluated);  // at the floor
    CHECK(ref::log_error(-1e-9, 1e-9, 0.0, 0.0).status == MetricStatus::not_evaluated);
    CHECK(ref::log_error(1e-9, 1e-9, -1.0, 0.0).status == MetricStatus::invalid_input);
}

TEST_CASE("field norms", "[verify][metrics]") {
    const std::array<double, 4> e{3.0, -4.0, 0.0, 12.0};
    const std::array<double, 4> w{1.0, 1.0, 7.0, 1.0};
    CHECK(ref::norm_linf(e).value == 12.0);
    const std::array<double, 3> neg{1.0, -7.0, 2.0};
    CHECK(ref::norm_linf(neg).value == 7.0);
    CHECK(ref::norm_l2(e, w).value == 13.0);
    const std::array<double, 4> w2{4.0, 0.25, 1.0, 0.0};  // sqrt(36 + 4) = sqrt(40)
    CHECK(std::abs(ref::norm_l2(e, w2).value - std::sqrt(40.0)) <= 2.0 * eps * std::sqrt(40.0));
    CHECK(ref::seminorm_h1(e, w).value == 13.0);
    // Scaled far from 1: no overflow or underflow in the squares (powers of two keep it exact).
    const std::array<double, 2> big{0x1p+600 * 3.0, 0x1p+600 * 4.0};
    const std::array<double, 2> tiny_w{0x1p-700, 0x1p-700};
    CHECK(ref::norm_l2(big, tiny_w).value == 5.0 * 0x1p+250);
    const std::array<double, 2> small{0x1p-600 * 3.0, 0x1p-600 * 4.0};
    const std::array<double, 2> ones{1.0, 1.0};
    CHECK(ref::norm_l2(small, ones).value == 5.0 * 0x1p-600);
    // An odd weight exponent.
    const std::array<double, 2> wodd{0x1p+3, 0x1p+3};
    const std::array<double, 2> e34{3.0, 4.0};
    CHECK(std::abs(ref::norm_l2(e34, wodd).value -
                   5.0 * std::sqrt(8.0)) <= 2.0 * eps * 5.0 * std::sqrt(8.0));
    // Refusals.
    const std::array<double, 3> w3{1.0, 1.0, 1.0};
    CHECK(ref::norm_l2(e, w3).status == MetricStatus::invalid_input);
    const std::array<double, 4> wneg{1.0, -1.0, 1.0, 1.0};
    CHECK(ref::norm_l2(e, wneg).status == MetricStatus::invalid_input);
    CHECK(ref::seminorm_h1(e, wneg).status == MetricStatus::invalid_input);
    CHECK(ref::norm_linf(std::span<const double>{}).status == MetricStatus::invalid_input);
    const std::array<double, 4> zero{};
    CHECK(ref::norm_l2(zero, w).value == 0.0);
}
