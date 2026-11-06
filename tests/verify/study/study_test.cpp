// verify/study (Unit 55 TVR-2, DECISIONS Unit 55 14.5-14.7, CV-Q6-Q8): order studies against exact solutions, the
// iterative-error rule, the asymptotic-range check from 4 levels and the order band, on synthetic ladders whose
// orders are known exactly; and the unstructured mean spacing (V / N)^(1/D).
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <limits>
#include <vector>

#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/verify/study/study.hpp"

namespace st = NiTCAD::verify::study;

namespace {

struct Ladder {
    std::vector<double> h, e, u;
};

// e = C h^p on h = h0 r^k, k = 0 .. levels - 1; U_iter = fraction e.
Ladder power(double p, int levels, double r = 2.0, double fraction = 0.0) {
    Ladder l;
    for (int k = 0; k < levels; ++k) {
        const double h = 1e-3 * std::pow(r, k);
        l.h.push_back(h);
        l.e.push_back(3.0 * std::pow(h, p));
        l.u.push_back(fraction * l.e.back());
    }
    return l;
}

st::OrderStudy study(const Ladder& l) { return st::order_study(l.h, l.e, l.u); }

}  // namespace

TEST_CASE("study: exact orders of a power-law ladder", "[study]") {
    for (const double p : {1.0, 2.0, 1.5}) {
        for (const double r : {2.0, 1.7}) {  // any ratio: the exact solution removes constant r
            const st::OrderStudy s = study(power(p, 4, r));
            REQUIRE(s.status == st::StudyStatus::ok);
            REQUIRE(s.orders.size() == 3);
            for (const double o : s.orders) CHECK(std::abs(o - p) <= 1e-12);
            CHECK(std::abs(s.observed - p) <= 1e-12);
            CHECK(s.asymptotic_evaluated);
            CHECK(s.in_asymptotic_range);
            CHECK(s.order_spread <= 1e-12);
        }
    }
}

TEST_CASE("study: the order band [p_th - 0.2, p_th + 0.3]", "[study]") {
    CHECK(st::judge_order(study(power(2.0, 4)), {.p_th = 2.0}).pass);
    CHECK(st::judge_order(study(power(1.81, 4)), {.p_th = 2.0}).pass);
    CHECK(st::judge_order(study(power(2.29, 4)), {.p_th = 2.0}).pass);
    const st::OrderVerdict low = st::judge_order(study(power(1.79, 4)), {.p_th = 2.0});
    CHECK(low.judged);
    CHECK_FALSE(low.pass);
    CHECK_FALSE(st::judge_order(study(power(2.31, 4)), {.p_th = 2.0}).pass);
    CHECK_FALSE(st::judge_order(study(power(2.0, 4)), {.p_th = 1.0}).pass);  // 2 > 1.3
    CHECK(st::judge_order(study(power(1.0, 4)), {.p_th = 1.0}).pass);
}

TEST_CASE("study: three levels judge the finest pair without the asymptotic check", "[study]") {
    Ladder l = power(2.0, 3);
    l.e[2] *= 10.0;  // the coarsest pair far off: not checked with 3 levels
    const st::OrderStudy s = study(l);
    REQUIRE(s.status == st::StudyStatus::ok);
    CHECK_FALSE(s.asymptotic_evaluated);
    CHECK(std::abs(s.observed - 2.0) <= 1e-12);
    const st::OrderVerdict v = st::judge_order(s, {.p_th = 2.0});
    CHECK(v.judged);
    CHECK(v.pass);
}

TEST_CASE("study: outside the asymptotic range, not judged", "[study]") {
    // e = h + 300 h^2 on h = 1e-3 .. 8e-3: the orders drift from first towards second by more than 0.1.
    Ladder l;
    for (int k = 0; k < 4; ++k) {
        const double h = 1e-3 * std::pow(2.0, k);
        l.h.push_back(h);
        l.e.push_back(h + 300.0 * h * h);
        l.u.push_back(0.0);
    }
    const st::OrderStudy s = study(l);
    REQUIRE(s.status == st::StudyStatus::ok);
    CHECK(s.asymptotic_evaluated);
    CHECK(s.order_spread > 0.1);
    CHECK_FALSE(s.in_asymptotic_range);
    const st::OrderVerdict v = st::judge_order(s, {.p_th = 1.0});
    CHECK_FALSE(v.judged);
    CHECK(v.reason == "not in the asymptotic range");
}

TEST_CASE("study: the iterative error must be below a tenth of the discretisation error (CV-Q8)", "[study]") {
    CHECK(st::judge_order(study(power(2.0, 4, 2.0, 0.1)), {.p_th = 2.0}).judged);  // exactly a tenth
    const st::OrderVerdict v = st::judge_order(study(power(2.0, 4, 2.0, 0.11)), {.p_th = 2.0});
    CHECK_FALSE(v.judged);
    CHECK(v.reason == "iterative error above a tenth of the discretisation error (CV-Q8)");
}

TEST_CASE("study: invalid ladders", "[study]") {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const std::vector<double> h{1e-3, 2e-3, 4e-3}, e{1.0, 4.0, 16.0}, u{0.0, 0.0, 0.0};
    CHECK(st::order_study(std::vector<double>{1e-3}, std::vector<double>{1.0}, std::vector<double>{0.0}).status ==
          st::StudyStatus::invalid_input);
    CHECK(st::order_study(h, std::vector<double>{1.0, 4.0}, u).status == st::StudyStatus::invalid_input);
    CHECK(st::order_study(std::vector<double>{2e-3, 1e-3, 4e-3}, e, u).status == st::StudyStatus::invalid_input);
    CHECK(st::order_study(h, std::vector<double>{1.0, 0.0, 16.0}, u).status == st::StudyStatus::invalid_input);
    CHECK(st::order_study(h, std::vector<double>{1.0, nan, 16.0}, u).status == st::StudyStatus::invalid_input);
    CHECK(st::order_study(h, e, std::vector<double>{0.0, -1.0, 0.0}).status == st::StudyStatus::invalid_input);
    CHECK_FALSE(st::judge_order(st::order_study(h, std::vector<double>{1.0, 0.0, 16.0}, u), {.p_th = 2.0}).judged);
}

TEST_CASE("study: mean spacing (V / N)^(1/D)", "[study]") {
    std::vector<double> x;
    for (int i = 0; i <= 10; ++i) x.push_back(0.1 * i);
    const auto one = *NiTCAD::mesh::make_tensor_grid(x);
    CHECK(std::abs(st::mean_spacing(one) - 1.0 / 11.0) <= 1e-15);
    const auto two = *NiTCAD::mesh::make_tensor_grid(x, x);
    CHECK(std::abs(st::mean_spacing(two) - 1.0 / 11.0) <= 1e-15);
    const auto three = *NiTCAD::mesh::make_tensor_grid(x, x, x);
    CHECK(std::abs(st::mean_spacing(three) - 1.0 / 11.0) <= 1e-15);
}
