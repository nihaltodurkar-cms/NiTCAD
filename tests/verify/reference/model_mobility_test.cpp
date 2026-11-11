// TVR-3 step 4 (DECISIONS Unit 55 15.9, frozen before code): L0 tests of the mobility references of TVB-31.06-31.10
// (verify/reference/model_mobility.hpp).
//   - Every analytic derivative judged in the benchmark against complex step (TC8) of the reference's value formula,
//     |d - d_cs| <= 64 eps S(d) (15.9, "Partials"); the truncation estimate below 1e-3 of that bound.
//   - 31.10's root is a double-double iteration, so complex step cannot pass through it. Its second method (15.3) is
//     Brent's zero (TC5) of x^2 m(x) - a, to 1e-28 relative; its derivatives are checked against double-double
//     central differences (step 2^-27 relative, error about 2^-54 relative) to the same 64 eps S(d).
//   - Special points: the right-hand Canali slope at E = 0; the regularization's continuity of x and dx/da at a0.
#include "NiTCAD/verify/reference/model_mobility.hpp"

#include "NiTCAD/verify/reference/brent.hpp"
#include "NiTCAD/verify/reference/complex_step.hpp"
#include "NiTCAD/verify/reference/condition.hpp"
#include "NiTCAD/verify/reference/double_double.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

namespace ref = NiTCAD::verify::reference;
namespace mm = NiTCAD::verify::reference::model_mobility;
using ref::complex_d;
using ref::DD;

namespace {

constexpr double eps = std::numeric_limits<double>::epsilon();

std::vector<double> decades(double lo, double hi, int per) {
    std::vector<double> v;
    for (double e = std::log10(lo); e <= std::log10(hi) + 1e-9; e += 1.0 / per) v.push_back(std::pow(10.0, e));
    return v;
}

template <std::size_t N, class Fc, class Fd>
void check_cs(const char* what, double at, Fc&& f_complex, Fd&& d_of, const std::array<DD, N>& inputs) {
    const auto cs = ref::complex_step(f_complex, at);
    REQUIRE(cs.status == ref::ComplexStepStatus::ok);
    const DD d = d_of(inputs);
    const double S = ref::condition_sum(d_of, inputs);
    // 15.9's partial rule in full (MI-Q13): max(1e-13 |d|, 64 eps S(d)).
    const double bound = std::fmax(1e-13 * std::abs(d.to_double()), 64.0 * eps * S);
    INFO(what << " at " << at << ": analytic " << d.to_double() << ", complex step " << cs.derivative << ", bound "
              << bound);
    CHECK(std::abs(d.to_double() - cs.derivative) <= bound);
    CHECK(cs.truncation_estimate <= 1e-3 * bound);
}

template <class S>
std::array<S, 10> masetti_boron() {
    return {S(44.9), S(0.0), S(470.5), S(29.0), S(9.23e16), S(2.23e17), S(6.1e20), S(0.719), S(2.0), S(-2.2)};
}
template <class S>
std::array<S, 10> masetti_phosphorus() {
    return {S(68.5), S(68.5), S(1414.0), S(56.1), S(0.0), S(9.20e16), S(3.41e20), S(0.711), S(1.98), S(-2.5)};
}

}  // namespace

TEST_CASE("model_mobility: Caughey-Thomas and Masetti T slopes against complex step", "[verify][model_mobility]") {
    for (const double N : {0.0, 1e13, 1e16, 1.3e17, 1e19, 1e21}) {
        for (const double T : decades(200.0, 500.0, 16)) {
            const std::array<DD, 6> x{DD(1360.0), DD(1.3e17), DD(0.91), DD(-2.33), DD(N), DD(T)};
            check_cs(
                "CT dmu/dT", T,
                [&](complex_d t) {
                    return mm::ct_mobility(complex_d(92.0), complex_d(1360.0), complex_d(1.3e17), complex_d(0.91),
                                           complex_d(-2.33), complex_d(N), t);
                },
                [](const std::array<DD, 6>& y) { return mm::ct_mobility_slope(y[0], y[1], y[2], y[3], y[4], y[5]); },
                x);
            for (const bool boron : {true, false}) {
                const auto pc = boron ? masetti_boron<complex_d>() : masetti_phosphorus<complex_d>();
                const auto pd = boron ? masetti_boron<DD>() : masetti_phosphorus<DD>();
                std::array<DD, 12> xm{};
                for (int i = 0; i < 10; ++i) xm[i] = pd[i];
                xm[10] = DD(N);
                xm[11] = DD(T);
                check_cs(
                    "Masetti dmu/dT", T, [&](complex_d t) { return mm::masetti_mobility(pc, complex_d(N), t); },
                    [](const std::array<DD, 12>& y) {
                        std::array<DD, 10> p{};
                        for (int i = 0; i < 10; ++i) p[i] = y[i];
                        return mm::masetti_mobility_slope(p, y[10], y[11]);
                    },
                    xm);
            }
        }
    }
}

TEST_CASE("model_mobility: Canali dmu/dE against complex step, and at E = 0", "[verify][model_mobility]") {
    for (const double beta : {1.0, 1.109, 1.213, 2.0}) {
        for (const double E : decades(1e-2, 1e6, 8)) {
            const std::array<DD, 4> x{DD(1360.0), DD(1.07e7), DD(beta), DD(E)};
            check_cs(
                "Canali dmu/dE", E,
                [&](complex_d e) { return mm::canali_mobility(complex_d(1360.0), complex_d(1.07e7), complex_d(beta), e); },
                [](const std::array<DD, 4>& y) { return mm::canali_mobility_slope(y[0], y[1], y[2], y[3]); }, x);
        }
        // At E = 0: mu0, and the right-hand slope's limit (the slope at E = 1e-12 tends to it).
        CHECK(mm::canali_mobility(DD(1360.0), DD(1.07e7), DD(beta), DD(0.0)) == DD(1360.0));
        const double at0 = mm::canali_mobility_slope(DD(1360.0), DD(1.07e7), DD(beta), DD(0.0)).to_double();
        const double near0 = mm::canali_mobility_slope(DD(1360.0), DD(1.07e7), DD(beta), DD(1e-12)).to_double();
        if (beta == 1.0) {
            CHECK(std::abs(near0 - at0) <= 1e-9 * std::abs(at0));
        } else {
            // The slope falls like x^(beta - 1): to 0, slowly for beta near 1.
            CHECK(at0 == 0.0);
            const double nearer = mm::canali_mobility_slope(DD(1360.0), DD(1.07e7), DD(beta), DD(1e-30)).to_double();
            CHECK(std::abs(nearer) < std::abs(near0));
            CHECK(near0 < 0.0);
        }
    }
}

TEST_CASE("model_mobility: Lombardi partials against complex step", "[verify][model_mobility]") {
    const std::array<double, 4> p{4.75e7, 1.74e5, 0.125, 5.82e14};
    auto pc = [&] { return std::array<complex_d, 4>{p[0], p[1], p[2], p[3]}; };
    for (const double N : {0.0, 1e17, 1e19}) {
        for (const double E : decades(1.0, 1e7, 8)) {
            for (const double T : {250.0, 300.0, 400.0}) {
                const std::array<DD, 9> x{DD(p[0]), DD(p[1]), DD(p[2]), DD(p[3]), DD(1000.0), DD(N), DD(T), DD(E),
                                          DD(0.0)};
                auto part = [](int k) {
                    return [k](const std::array<DD, 9>& y) {
                        const std::array<DD, 4> q{y[0], y[1], y[2], y[3]};
                        return mm::lombardi_partials(q, y[4], y[5], y[6], y[7])[k];
                    };
                };
                check_cs(
                    "Lombardi dmu/dE", E,
                    [&](complex_d e) { return mm::lombardi_mobility(pc(), complex_d(1000.0), complex_d(N), complex_d(T), e); },
                    part(0), x);
                check_cs(
                    "Lombardi dmu/dT", T,
                    [&](complex_d t) { return mm::lombardi_mobility(pc(), complex_d(1000.0), complex_d(N), t, complex_d(E)); },
                    part(1), x);
                check_cs(
                    "Lombardi dmu/dmu_b", 1000.0,
                    [&](complex_d b) { return mm::lombardi_mobility(pc(), b, complex_d(N), complex_d(T), complex_d(E)); },
                    part(2), x);
            }
        }
    }
}

TEST_CASE("model_mobility: the balance root against Brent (TC5), its slope, the regularization", "[verify][model_mobility]") {
    for (const double beta : {1.0, 1.5, 2.0}) {
        const DD b(beta);
        for (const double a : decades(1e-8, 1e6, 8)) {
            const mm::Root r = mm::balance_root(DD(a), beta);
            auto f = [&](const DD& x) { return x * x * mm::balance_m(x, b) - DD(a); };
            const double hi = a < 1.0 ? 4.0 * std::sqrt(a) : 4.0 * a;
            const auto br = ref::brent_zero(f, DD(0.0), DD(hi), 1e-30 * std::sqrt(a));
            REQUIRE(br.status == ref::BrentStatus::ok);
            INFO("beta " << beta << " a " << a);
            CHECK(std::abs(((r.x - br.root) / r.x).hi) <= 1e-27);
            // dx/da against a double-double central difference.
            const DD h = ref::ldexp(DD(a), -27);
            const DD cd = (mm::balance_root(DD(a) + h, beta).x - mm::balance_root(DD(a) - h, beta).x) / ref::ldexp(h, 1);
            CHECK(std::abs(((r.d_a - cd) / r.d_a).hi) <= 1e-14);
        }
        // The regularization: x and dx/da continuous at a0; 0 at and below 0.
        const DD a0(mm::regularization_a0);
        const mm::Root at = mm::regularized_root(a0, beta);
        const mm::Root below = mm::regularized_root(DD(std::nextafter(mm::regularization_a0, 0.0)), beta);
        CHECK(std::abs(((at.x - below.x) / at.x).hi) <= 1e-15);
        CHECK(std::abs(((at.d_a - below.d_a) / at.d_a).hi) <= 1e-14);
        CHECK(mm::regularized_root(DD(0.0), beta).x == DD(0.0));
        CHECK(mm::regularized_root(DD(-1e-3), beta).d_a == DD(0.0));
    }
}

TEST_CASE("model_mobility: carrier-temperature partials against central differences", "[verify][model_mobility]") {
    using A4 = std::array<DD, 4>;
    for (const double beta : {1.0, 1.5, 2.0}) {
        const A4 base{DD(1360.0), DD(0.0), DD(0.35e-12), DD(1.07e7)};
        const DD per_K = mm::balance_a(base[0], DD(1.0), base[2], base[3]);
        for (const double a : {-1e-3, -1e-9, 1e-9, 2e-8, 1e-4, 1e-2, 1.0, 1e2, 1e5}) {
            A4 x = base;
            x[1] = DD(a) / per_K;
            const mm::BalanceResult r = mm::carrier_temperature_mobility(x, beta);
            auto field = [&](int k) {
                return [k, beta](const A4& y) {
                    const mm::BalanceResult q = mm::carrier_temperature_mobility(y, beta);
                    return k == 0 ? q.mobility : q.field;
                };
            };
            auto cd = [&](int k, std::size_t i) {
                const DD h = ref::ldexp(ref::abs(x[i]), -27);
                A4 up = x, dn = x;
                up[i] = x[i] + h;
                dn[i] = x[i] - h;
                return ((field(k)(up) - field(k)(dn)) / ref::ldexp(h, 1)).to_double();
            };
            auto d_of = [beta](int which) {
                return [which, beta](const A4& y) {
                    const mm::BalanceResult q = mm::carrier_temperature_mobility(y, beta);
                    const DD v[4] = {q.d_mu0, q.d_rise, q.field_d_mu0, q.field_d_rise};
                    return v[which];
                };
            };
            const double got[4] = {r.d_mu0.to_double(), r.d_rise.to_double(), r.field_d_mu0.to_double(),
                                   r.field_d_rise.to_double()};
            const double want[4] = {cd(0, 0), cd(0, 1), cd(1, 0), cd(1, 1)};
            for (int k = 0; k < 4; ++k) {
                // Below a0 for beta != 2 the regularization's cubic is exact; crossing a0 within a step is avoided by
                // the points chosen (a0 / 2^27 away at least).
                const double S = ref::condition_sum(d_of(k), x);
                INFO("beta " << beta << " a " << a << " partial " << k << ": " << got[k] << " vs " << want[k]);
                CHECK(std::abs(got[k] - want[k]) <= std::fmax(1e-13 * std::abs(got[k]), 64.0 * eps * S));
            }
        }
    }
}
