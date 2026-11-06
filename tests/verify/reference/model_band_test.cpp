// TVR-3 step 3 (DECISIONS Unit 55 15.9, frozen before code): L0 tests of the band-set references of TVB-31.01-31.05
// (verify/reference/model_band.hpp) and of the register's decimal reader (decimal.hpp).
//   - Each analytic derivative the benchmark judges NiTCAD's partials against is checked against complex step (TC8)
//     of the reference's own value formula: |d_analytic - d_cs| <= 64 eps S(d), S(d) = sum_i |x_i dd/dx_i| over every
//     input (15.9, "Partials"). The complex-step truncation estimate is required to be below 1e-3 of that bound.
//   - Identities that tie the forms to their sources: Varshni at 0 K is Eg0; Passler near 0 K is Eg0 less its
//     leading chi^2 term; n_i^2 = Nc Nv e^(-Eg/kT); N(300 K) = N300; k = CODATA k / q; the
//     documented Slotboom form is the published one minus E0 sqrt(1/2) above N0 and 0 at N0.
//   - decimal: parse and the exact M 10^e reading against hand-written cases.
#include "NiTCAD/verify/reference/model_band.hpp"

#include "NiTCAD/verify/reference/complex_step.hpp"
#include "NiTCAD/verify/reference/condition.hpp"
#include "NiTCAD/verify/reference/decimal.hpp"
#include "NiTCAD/verify/reference/double_double.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <numbers>
#include <vector>

namespace ref = NiTCAD::verify::reference;
namespace mb = NiTCAD::verify::reference::model_band;
using ref::complex_d;
using ref::DD;

namespace {

constexpr double eps = std::numeric_limits<double>::epsilon();

// T from 1 K to 600 K, 16 points a decade, and 300 K.
std::vector<double> temperatures() {
    std::vector<double> t;
    for (int i = 0; i <= static_cast<int>(16.0 * std::log10(600.0)); ++i) t.push_back(std::pow(10.0, i / 16.0));
    t.push_back(300.0);
    t.push_back(600.0);
    return t;
}

// Inputs as doubles, so the complex-step and double-double evaluations see the same numbers.
struct Legacy {
    double Eg0 = 1.17, alpha = 4.73e-4, beta = 636.0, Nc300 = 2.86e19, Nv300 = 3.10e19;
};
struct Couderc {
    double Eg0 = 1.17, alpha = 3.23e-4, Theta = 446.0, Delta = 0.51;
    std::array<double, 4> c{-4.609e-10, 6.753e-7, -1.312e-5, 1.094};
    std::array<double, 4> v{2.525e-9, -4.689e-6, 3.376e-3, 0.3426};
};

template <class S>
std::array<S, 4> as(const std::array<double, 4>& a) {
    return {S(a[0]), S(a[1]), S(a[2]), S(a[3])};
}

// Checks d_analytic (double-double) against the complex step of f at T, to 64 eps S(d); S from `d_of`, the analytic
// derivative as a function of all inputs.
template <std::size_t N, class Fc, class Fd>
void check_derivative(const char* what, double T, Fc&& f_complex, const DD& d_analytic, Fd&& d_of,
                      const std::array<DD, N>& inputs) {
    const auto cs = ref::complex_step(f_complex, T);
    REQUIRE(cs.status == ref::ComplexStepStatus::ok);
    const double S = ref::condition_sum(d_of, inputs);
    const double bound = std::fmax(1e-13 * std::abs(d_analytic.to_double()), 64.0 * eps * S);  // 15.9 in full (MI-Q13)
    const double diff = std::abs(d_analytic.to_double() - cs.derivative);
    INFO(what << " at T = " << T << ": analytic " << d_analytic.to_double() << ", complex step " << cs.derivative
              << ", |diff| " << diff << ", bound " << bound);
    CHECK(diff <= bound);
    CHECK(cs.truncation_estimate <= 1e-3 * bound);
}

}  // namespace

TEST_CASE("model_band: Varshni and power-law DOS derivatives against complex step", "[verify][model_band]") {
    const Legacy p;
    for (const double T : temperatures()) {
        const std::array<DD, 4> x{DD(p.Eg0), DD(p.alpha), DD(p.beta), DD(T)};
        check_derivative(
            "dEg/dT Varshni", T,
            [&](complex_d t) { return mb::varshni_gap(complex_d(p.Eg0), complex_d(p.alpha), complex_d(p.beta), t); },
            mb::varshni_gap_slope(DD(p.alpha), DD(p.beta), DD(T)),
            [](const std::array<DD, 4>& y) { return mb::varshni_gap_slope(y[1], y[2], y[3]); }, x);
        const std::array<DD, 2> xn{DD(p.Nc300), DD(T)};
        check_derivative(
            "d ln N/dT power law", T, [&](complex_d t) { return std::log(mb::power_dos(complex_d(p.Nc300), t)); },
            mb::power_dos_log_slope(DD(T)), [](const std::array<DD, 2>& y) { return mb::power_dos_log_slope(y[1]); },
            xn);
    }
}

TEST_CASE("model_band: Passler gap and Couderc DOS derivatives against complex step", "[verify][model_band]") {
    const Couderc p;
    for (const double T : temperatures()) {
        const std::array<DD, 4> x{DD(p.alpha), DD(p.Theta), DD(p.Delta), DD(T)};
        check_derivative(
            "dEg/dT Passler", T,
            [&](complex_d t) {
                return mb::passler_gap(complex_d(p.Eg0), complex_d(p.alpha), complex_d(p.Theta), complex_d(p.Delta), t);
            },
            mb::passler_gap_slope(DD(p.alpha), DD(p.Theta), DD(p.Delta), DD(T)),
            [](const std::array<DD, 4>& y) { return mb::passler_gap_slope(y[0], y[1], y[2], y[3]); }, x);
        for (const auto* c : {&p.c, &p.v}) {
            const std::array<DD, 5> xc{DD((*c)[0]), DD((*c)[1]), DD((*c)[2]), DD((*c)[3]), DD(T)};
            auto slope = [](const std::array<DD, 5>& y) {
                return mb::couderc_dos_log_slope(std::array<DD, 4>{y[0], y[1], y[2], y[3]}, y[4]);
            };
            check_derivative(
                "d ln N/dT Couderc", T,
                [&](complex_d t) { return std::log(mb::couderc_dos(as<complex_d>(*c), t)); },
                mb::couderc_dos_log_slope(as<DD>(*c), DD(T)), slope, xc);
            auto mslope = [](const std::array<DD, 5>& y) {
                return mb::mass_polynomial_slope(std::array<DD, 4>{y[0], y[1], y[2], y[3]}, y[4]);
            };
            check_derivative(
                "dM/dT", T, [&](complex_d t) { return mb::mass_polynomial(as<complex_d>(*c), t); },
                mb::mass_polynomial_slope(as<DD>(*c), DD(T)), mslope, xc);
        }
    }
}

TEST_CASE("model_band: d ln n_i/dT and d ln n_ie/dT against complex step", "[verify][model_band]") {
    const Legacy l;
    const Couderc c;
    const double dEg_bgn = 0.05;  // a fixed narrowing [eV]
    for (const double T : temperatures()) {
        // Legacy set: inputs Eg0, alpha, beta, Nc300, Nv300, dEg_bgn, T.
        auto ln_ni_legacy = [&](complex_d t, double bgn) {
            const complex_d Eg = mb::varshni_gap(complex_d(l.Eg0), complex_d(l.alpha), complex_d(l.beta), t);
            return mb::intrinsic_log_density(mb::power_dos(complex_d(l.Nc300), t), mb::power_dos(complex_d(l.Nv300), t),
                                             Eg, t) +
                   complex_d(bgn) / (2.0 * mb::detail::k_eV<complex_d>() * t);
        };
        auto slope_legacy = [](const std::array<DD, 7>& y) {
            const DD Eg = mb::varshni_gap(y[0], y[1], y[2], y[6]);
            const DD dEg = mb::varshni_gap_slope(y[1], y[2], y[6]);
            const DD dlnN = mb::power_dos_log_slope(y[6]);
            return mb::intrinsic_log_slope(dlnN, dlnN, Eg, dEg, y[5], y[6]);
        };
        for (const double bgn : {0.0, dEg_bgn}) {
            const std::array<DD, 7> x{DD(l.Eg0), DD(l.alpha), DD(l.beta), DD(l.Nc300), DD(l.Nv300), DD(bgn), DD(T)};
            check_derivative(
                bgn == 0.0 ? "d ln n_i/dT legacy" : "d ln n_ie/dT legacy", T,
                [&](complex_d t) { return ln_ni_legacy(t, bgn); }, slope_legacy(x), slope_legacy, x);
        }
        // Default set: inputs Eg0, alpha, Theta, Delta, the eight mass coefficients, T.
        auto slope_couderc = [](const std::array<DD, 13>& y) {
            const std::array<DD, 4> cc{y[4], y[5], y[6], y[7]};
            const std::array<DD, 4> vv{y[8], y[9], y[10], y[11]};
            const DD Eg = mb::passler_gap(y[0], y[1], y[2], y[3], y[12]);
            const DD dEg = mb::passler_gap_slope(y[1], y[2], y[3], y[12]);
            return mb::intrinsic_log_slope(mb::couderc_dos_log_slope(cc, y[12]), mb::couderc_dos_log_slope(vv, y[12]),
                                           Eg, dEg, DD(0.0), y[12]);
        };
        const std::array<DD, 13> x{DD(c.Eg0),  DD(c.alpha), DD(c.Theta), DD(c.Delta), DD(c.c[0]),
                                   DD(c.c[1]), DD(c.c[2]),  DD(c.c[3]),  DD(c.v[0]),  DD(c.v[1]),
                                   DD(c.v[2]), DD(c.v[3]),  DD(T)};
        check_derivative(
            "d ln n_i/dT Couderc", T,
            [&](complex_d t) {
                const complex_d Eg = mb::passler_gap(complex_d(c.Eg0), complex_d(c.alpha), complex_d(c.Theta),
                                                     complex_d(c.Delta), t);
                return mb::intrinsic_log_density(mb::couderc_dos(as<complex_d>(c.c), t),
                                                 mb::couderc_dos(as<complex_d>(c.v), t), Eg, t);
            },
            slope_couderc(x), slope_couderc, x);
    }
}

TEST_CASE("model_band: identities of the forms", "[verify][model_band]") {
    const Legacy l;
    const Couderc c;
    // Both gaps are Eg0 as T -> 0 (Varshni exactly at 0; Passler's first term is e^-446 at 1 K, its second
    // (3 Delta^2 / 2)(P^(1/6) - 1) ~ alpha Theta 0.39 pi^2 chi^2 / (18 (1 + Delta^2)) with chi = 2 / 446).
    CHECK(mb::varshni_gap(DD(l.Eg0), DD(l.alpha), DD(l.beta), DD(0.0)) == DD(l.Eg0));
    const DD at1 = mb::passler_gap(DD(c.Eg0), DD(c.alpha), DD(c.Theta), DD(c.Delta), DD(1.0));
    const double chi = 2.0 / 446.0;
    const double second = c.alpha * c.Theta * 1.5 * c.Delta * c.Delta * (std::pow(std::numbers::pi * chi, 2) /
                                                                          (3.0 * (1.0 + c.Delta * c.Delta)) / 6.0);
    CHECK(std::abs((DD(c.Eg0) - at1).to_double() - second) <= 1e-3 * second);
    // n_i^2 = Nc Nv e^(-Eg / kT) (Couderc eq (1)), at 300 K.
    const DD T(300.0);
    const DD Nc = mb::power_dos(DD(l.Nc300), T);
    const DD Nv = mb::power_dos(DD(l.Nv300), T);
    const DD Eg = mb::varshni_gap(DD(l.Eg0), DD(l.alpha), DD(l.beta), T);
    const DD ni = mb::intrinsic_density(Nc, Nv, Eg, T);
    const DD rhs = Nc * Nv * ref::exp(-Eg / (mb::detail::k_eV_dd() * T));
    CHECK(std::abs(((ni * ni - rhs) / rhs).hi) <= 1e-28);
    // The power law is N300 at 300 K exactly.
    CHECK(mb::power_dos(DD(l.Nc300), T) == DD(l.Nc300));
    // Slotboom: the documented form is 0 at N0 and below, the published one minus E0 sqrt(1/2) above.
    const DD E0(6.92e-3), N0(1.3e17);
    CHECK(mb::slotboom_documented(E0, N0, N0) == DD(0.0));
    CHECK(mb::slotboom_documented(E0, N0, DD(1e16)) == DD(0.0));
    for (const double N : {1.3000000000000002e17, 1e18, 1e20, 1e21}) {
        const DD diff = mb::slotboom_published(E0, N0, DD(N)) - mb::slotboom_documented(E0, N0, DD(N));
        CHECK(std::abs((diff - E0 * ref::sqrt(DD(0.5))).hi) <= 1e-30);
    }
    // The published form at N0 is E0 sqrt(1/2) (4.89 meV), the F1 offset of DECISIONS 15.4.
    CHECK(std::abs(mb::slotboom_published(E0, N0, N0).to_double() - 6.92e-3 * std::sqrt(0.5)) <= 1e-18);
    // k: CODATA 2018's exact k / q, printed truncated as 8.617 333 262... e-5 eV/K: within the truncation.
    const double k = mb::detail::k_eV_dd().to_double();
    CHECK(k >= 8.617333262e-5);
    CHECK(k < 8.617333263e-5);
}

TEST_CASE("decimal: parse and exact reading of printed values", "[verify][model_band]") {
    namespace dec = NiTCAD::verify::reference::decimal;
    const auto p = dec::parse("-4.609e-10");
    CHECK(p.negative);
    CHECK(p.digits == "4609");
    CHECK(p.exponent == -13);
    CHECK(dec::parse("636").exponent == 0);
    CHECK(dec::parse("0.3426").digits == "03426");
    CHECK(dec::parse("0.3426").exponent == -4);
    CHECK(dec::parse("1.3e+17").exponent == 16);
    CHECK_THROWS(dec::parse("1.2.3"));
    CHECK_THROWS(dec::parse("e5"));
    CHECK_THROWS(dec::parse("1e"));
    CHECK(dec::to_double("0.3426") == 0.3426);
    CHECK(dec::to_double("-4.609e-10") == -4.609e-10);
    // M 10^e: integers are exact; 0.3426 = 3426 / 10^4 to one double-double rounding.
    CHECK(dec::to_dd("636") == DD(636.0));
    CHECK(dec::to_dd("2.86e19") == DD(2.86e19));  // 286 10^17 = 2.86e19 exactly (< 2^106, 10^17 exact)
    const DD r = dec::to_dd("0.3426") * DD(10000.0);
    CHECK(std::abs((r - DD(3426.0)).hi) <= 3426.0 * 4.0 * ref::u_dd);
    CHECK(dec::to_dd("-4.609e-10").hi == -4.609e-10);
    CHECK(dec::to_dd("0") == DD(0.0));
    CHECK_THROWS(dec::to_dd("1.234567890123456"));
}
