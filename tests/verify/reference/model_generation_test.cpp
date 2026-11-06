// TVR-3 step 7 (DECISIONS Unit 55 15.9, frozen before code): L0 tests of the generation references of TVB-31.20-31.22
// (verify/reference/model_generation.hpp).
//   - vOdM and Kane derivatives against complex step (TC8) under 15.9's partial rule (MI-Q13); gamma(300 K) = 1 and
//     d gamma / dT against a double-double central difference.
//   - The edge-safe alpha and cos(theta) against the direct formulas inside the gap; the antiderivatives against
//     tanh-sinh quadrature (TC4) of kappa and 1 / kappa, and their whole-gap values pi sqrt(m_r E_g) / (4 hbar) and
//     pi hbar / (2 sqrt(m_r E_g)).
//   - Segment partials against double-double central differences of the segment value; the closed-form brackets
//     against quadrature of the difference.
//   - Path-rate partials against central differences; the bracket series against its closed form at z = 0.5.
//   - eq (8) against eq (11) with the exact uniform-field integrals (k_m^2 I_ik >> 1).
#include "NiTCAD/verify/reference/model_generation.hpp"

#include "NiTCAD/verify/reference/complex_step.hpp"
#include "NiTCAD/verify/reference/condition.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <limits>

namespace ref = NiTCAD::verify::reference;
namespace mg = NiTCAD::verify::reference::model_generation;
using ref::complex_d;
using ref::DD;

namespace {
constexpr double eps = std::numeric_limits<double>::epsilon();

mg::Band band(double gap_eV, double mr_m0) {
    return {DD(gap_eV) * DD(NiTCAD::base::q_C), DD(mr_m0) * DD(NiTCAD::base::m0_kg)};
}

template <std::size_t N, class F>
void cs_check(const char* what, F&& f, const std::array<double, N>& x, std::size_t i) {
    auto g = [&](complex_d t) {
        std::array<complex_d, N> y{};
        for (std::size_t k = 0; k < N; ++k) y[k] = x[k];
        y[i] = t;
        return f(y)[0];
    };
    const auto cs = ref::complex_step(g, x[i], 1e-30 * std::abs(x[i]));
    REQUIRE(cs.status == ref::ComplexStepStatus::ok);
    std::array<DD, N> xd{};
    for (std::size_t k = 0; k < N; ++k) xd[k] = DD(x[k]);
    const DD d = f(xd)[1];
    const auto S = ref::condition_sums<2>(f, xd);
    INFO(what << " at " << x[i] << ": " << d.to_double() << " vs " << cs.derivative);
    CHECK(std::abs(d.to_double() - cs.derivative) <= std::fmax(1e-13 * std::abs(d.to_double()), 64.0 * eps * S[1]));
}
}  // namespace

TEST_CASE("model_generation: vOdM and Kane", "[verify][model_generation]") {
    auto fv = [](const auto& y) { return mg::vodm(y); };
    auto fk = [](const auto& y) { return mg::kane(y); };
    for (double e = 4.0; e <= 7.0; e += 0.125) {
        const double E = std::pow(10.0, e);
        cs_check<4>("vOdM", fv, std::array<double, 4>{7.03e5, 1.231e6, 1.07, E}, 3);
        cs_check<3>("Kane", fk, std::array<double, 3>{3.5e21, 1.03e8, E}, 2);
    }
    const auto g300 = mg::vodm_gamma(std::array<DD, 2>{DD(0.063), DD(300.0)});
    CHECK(std::abs((g300[0] - DD(1.0)).hi) <= 4.0 * ref::u_dd);
    for (const double T : {200.0, 300.0, 450.0, 500.0}) {
        const DD h = ref::ldexp(DD(T), -27);
        const DD cd = (mg::vodm_gamma(std::array<DD, 2>{DD(0.063), DD(T) + h})[0] -
                       mg::vodm_gamma(std::array<DD, 2>{DD(0.063), DD(T) - h})[0]) /
                      ref::ldexp(h, 1);
        const DD d = mg::vodm_gamma(std::array<DD, 2>{DD(0.063), DD(T)})[1];
        INFO("T " << T);
        CHECK(std::abs(((d - cd) / d).hi) <= 1e-14);
    }
}

TEST_CASE("model_generation: kappa and its antiderivatives", "[verify][model_generation]") {
    for (const mg::Band b : {band(1.42, 0.0335), band(0.5, 1.0 / 3.0)}) {
        const DD u = mg::u_of(b);
        for (double d = 0.05; d < 1.0; d += 0.05) {
            const DD S = sqrt(u * (DD(d) - DD(0.5)) + u * u / DD(4.0) + DD(0.25));
            const DD alpha = -u + DD(2.0) * S;
            const mg::AlphaCos a = mg::alpha_cos(DD(d), u);
            CHECK(std::abs((a.alpha - alpha).hi) <= 1e-28 * std::fmax(1.0, u.hi));
            CHECK(std::abs((a.c - sqrt(DD(1.0) - alpha * alpha)).hi) <= 1e-27 * std::fmax(1.0, u.hi));
        }
        for (const std::array<double, 2> ab : {std::array<double, 2>{0.0, 1.0}, std::array<double, 2>{0.1, 0.4}, std::array<double, 2>{0.0, 0.5},
                                               std::array<double, 2>{0.0, 1e-6}, std::array<double, 2>{0.7, 1.0}}) {
            // Each half by t^2 substitution from its end (delta = a + t^2, b - t^2): smooth at the band edges, where
            // kappa ~ sqrt and 1 / kappa ~ 1 / sqrt.
            auto split = [&](auto&& f) {
                const DD a(ab[0]), c(ab[1]), m = (a + c) / DD(2.0);
                const auto l = ref::tanh_sinh([&](const DD& t) { return f(a + t * t) * DD(2.0) * t; }, DD(0.0),
                                              sqrt(m - a), 0.0, 1e-30);
                const auto r = ref::tanh_sinh([&](const DD& t) { return f(c - t * t) * DD(2.0) * t; }, DD(0.0),
                                              sqrt(c - m), 0.0, 1e-30);
                REQUIRE(l.status == ref::TanhSinhStatus::ok);
                REQUIRE(r.status == ref::TanhSinhStatus::ok);
                return l.value + r.value;
            };
            const DD qk = split([&](const DD& v) { return mg::kappa(v, b); });
            const DD qi = split([&](const DD& v) { return DD(1.0) / mg::kappa(v, b); });
            const DD K = mg::kappa_antiderivative(DD(ab[1]), b) - mg::kappa_antiderivative(DD(ab[0]), b);
            const DD G = mg::inverse_kappa_antiderivative(DD(ab[1]), b) - mg::inverse_kappa_antiderivative(DD(ab[0]), b);
            INFO("u " << u.hi << " [" << ab[0] << ", " << ab[1] << "]: K rel " << ((K - qk) / K).hi << ", G rel "
                      << ((G - qi) / G).hi);
            // The closed form is a difference of antiderivatives: its own rounding is 4 u_dd of their size.
            const double sk = std::abs(mg::kappa_antiderivative(DD(ab[1]), b).hi) +
                              std::abs(mg::kappa_antiderivative(DD(ab[0]), b).hi);
            const double sg = std::abs(mg::inverse_kappa_antiderivative(DD(ab[1]), b).hi) +
                              std::abs(mg::inverse_kappa_antiderivative(DD(ab[0]), b).hi);
            CHECK(std::abs((K - qk).hi) <= 1e-27 * std::abs(K.hi) + 4.0 * ref::u_dd * sk);
            // 1 - t^2 is not exact for tiny t, so the substitution cannot reach delta = 1 for the singular 1 / kappa: its
            // check there is the whole-gap identity below.
            if (ab[1] < 1.0) CHECK(std::abs((G - qi).hi) <= 1e-27 * std::abs(G.hi) + 4.0 * ref::u_dd * sg);
        }
        const DD s = mg::scale_of(b), h = mg::hbar();
        const DD K01 = mg::kappa_antiderivative(DD(1.0), b) - mg::kappa_antiderivative(DD(0.0), b);
        const DD G01 = mg::inverse_kappa_antiderivative(DD(1.0), b) - mg::inverse_kappa_antiderivative(DD(0.0), b);
        CHECK(std::abs((K01 / (ref::pi() * s / (DD(4.0) * h)) - DD(1.0)).hi) <= 1e-29);
        CHECK(std::abs((G01 / (ref::pi() * h / (DD(2.0) * s)) - DD(1.0)).hi) <= 1e-29);
    }
}

TEST_CASE("model_generation: segment partials and brackets", "[verify][model_generation]") {
    const mg::Band b = band(1.42, 0.0335);
    const DD L(1e-8);
    for (const std::array<double, 2> s : {std::array<double, 2>{0.1, 0.6}, std::array<double, 2>{0.6, 0.1},
                                          std::array<double, 2>{-0.2, 0.3}, std::array<double, 2>{0.8, 1.3},
                                          std::array<double, 2>{0.3, 0.3 + 1e-3}, std::array<double, 2>{0.5, 0.5 + 1e-6},
                                          std::array<double, 2>{0.2 + 1e-9, 0.2}, std::array<double, 2>{1e-6, 2e-6}}) {
        const mg::Segment seg = mg::segment(DD(s[0]), DD(s[1]), L, b);
        REQUIRE(seg.ok);
        const double D = std::abs(s[1] - s[0]);
        for (int end = 0; end < 2; ++end) {
            const DD h = DD(D * 0x1p-20);
            auto at = [&](const DD& v) {
                return end == 0 ? mg::segment(v, DD(s[1]), L, b) : mg::segment(DD(s[0]), v, L, b);
            };
            const DD x = DD(s[end]);
            const mg::Segment up = at(x + h), dn = at(x - h);
            REQUIRE(up.ok);
            REQUIRE(dn.ok);
            for (int k = 0; k < 2; ++k) {
                const DD cd = (up.q[k] - dn.q[k]) / ref::ldexp(h, 1);
                const DD d = seg.q[2 + 2 * k + end];
                INFO("segment [" << s[0] << ", " << s[1] << "] end " << end << " k " << k << ": " << d.hi << " vs "
                                 << cd.hi);
                // Central difference truncation (h / D)^2 = 2^-40 relative to the derivative's scale.
                CHECK(std::abs((d - cd).hi) <= 1e-11 * (std::abs(d.hi) + std::abs(seg.q[k].hi) / D));
            }
        }
    }
    // The closed-form bracket against quadrature of the difference, at D = 0.05 (both are valid there).
    const DD lo(0.3), hi(0.35);
    const DD fk_hi = mg::kappa(hi, b);
    const DD closed = (hi - lo) * fk_hi - (mg::kappa_antiderivative(hi, b) - mg::kappa_antiderivative(lo, b));
    const auto q = ref::tanh_sinh([&](const DD& v) { return fk_hi - mg::kappa(v, b); }, lo, hi, 0.0, 1e-28);
    REQUIRE(q.status == ref::TanhSinhStatus::ok);
    CHECK(std::abs(((closed - q.value) / q.value).hi) <= 1e-25);
    // The Kronrod-15 bracket (added after the first benchmark run, 15.17) against tanh-sinh at D = 1e-5 inside the
    // gap, where both are valid; and the D -> 0 limit at D = 1e-13: I_k -> L kappa(mid), d I_k / d v -> L kappa' / 2.
    const DD lo5(0.3), hi5 = DD(0.3) + DD(1e-5);
    const DD f5 = mg::kappa(hi5, b);
    const auto q5 = ref::tanh_sinh([&](const DD& v) { return f5 - mg::kappa(v, b); }, lo5, hi5, 0.0, 1e-28);
    REQUIRE(q5.status == ref::TanhSinhStatus::ok);
    const DD k5 = ref::qk15([&](const DD& v) { return f5 - mg::kappa(v, b); }, lo5, hi5).value;
    CHECK(std::abs(((k5 - q5.value) / q5.value).hi) <= 1e-25);
    const DD a13(0.3), b13 = DD(0.3) + DD(1e-13), mid = (a13 + b13) / DD(2.0);
    const mg::Segment tiny = mg::segment(a13, b13, L, b);
    REQUIRE(tiny.ok);
    const DD hk = DD(1e-7);
    const DD dk = (mg::kappa(mid + hk, b) - mg::kappa(mid - hk, b)) / ref::ldexp(hk, 1);
    CHECK(std::abs(((tiny.q[0] - L * mg::kappa(mid, b)) / tiny.q[0]).hi) <= 1e-24);
    for (int end = 2; end < 4; ++end) CHECK(std::abs(((tiny.q[end] - L * dk / DD(2.0)) / tiny.q[end]).hi) <= 1e-12);
}

TEST_CASE("model_generation: path rate and eq (8)", "[verify][model_generation]") {
    for (const double z : {1e-8, 1e-3, 0.5, 0.999, 1.001, 5.0, 50.0}) {
        for (const double Ik : {0.1, 10.0, 200.0}) {
            const std::array<DD, 4> x{DD(1.6e-19 * 1e8), DD(Ik), DD(1e-18), DD(z / 1e-18)};
            const auto r = mg::path_rate(x);
            for (std::size_t i = 0; i < 4; ++i) {
                const DD h = ref::ldexp(ref::abs(x[i]), -40);  // truncation (2h)^2 / 6 of e^(-2 I_k) at I_k = 200 needs a small step
                std::array<DD, 4> up = x, dn = x;
                up[i] = x[i] + h;
                dn[i] = x[i] - h;
                const DD cd = (mg::path_rate(up)[0] - mg::path_rate(dn)[0]) / ref::ldexp(h, 1);
                INFO("z " << z << " I_k " << Ik << " input " << i << ": " << r[1 + i].hi << " vs " << cd.hi << " rel " << ((r[1 + i] - cd) / r[1 + i]).hi);
                // The central difference's rounding, 4 u_dd 2^40 G / x_i, dominates where the partial is a small part of G / x_i
                // (dG/dI_ik at z << 1).
                CHECK(std::abs((r[1 + i] - cd).hi) <=
                      1e-14 * std::abs(r[1 + i].hi) + 4.0 * ref::u_dd * 0x1p40 * std::abs(r[0].hi) / std::abs(x[i].hi));
            }
        }
    }
    // The series of g(z) = 1 - e^-z (1 + z) at z = 0.5 against the closed form (no cancellation there).
    const DD z(0.5);
    const DD closed = -ref::expm1(-z) - z * exp(-z);
    const auto r = mg::path_rate(std::array<DD, 4>{DD(1.0), DD(0.0), DD(1.0), z});
    const DD pre = DD(1.0) / (DD(36.0) * mg::hbar());
    CHECK(std::abs((-r[3] / pre - closed).hi) <= 1e-30);
    // eq (8) against eq (11) on the exact uniform-field path.
    for (const mg::Band b : {band(1.42, 0.0335), band(0.5, 1.0 / 3.0)}) {
        for (const double F : {1e8, 1e9, 1e10}) {
            const DD q(NiTCAD::base::q_C);
            const DD Ltot = b.gap_J / (q * DD(F));
            const DD Ik = Ltot * (mg::kappa_antiderivative(DD(1.0), b) - mg::kappa_antiderivative(DD(0.0), b));
            const DD Iik =
                Ltot * (mg::inverse_kappa_antiderivative(DD(1.0), b) - mg::inverse_kappa_antiderivative(DD(0.0), b));
            const DD path = mg::path_rate(std::array<DD, 4>{q * DD(F), Ik, Iik, DD(1e60)})[0];
            const DD eq8 = mg::uniform_field_rate(DD(F), b);
            INFO("u " << mg::u_of(b).hi << " F " << F);
            CHECK(std::abs(((path - eq8) / eq8).hi) <= 1e-27);
        }
    }
}
