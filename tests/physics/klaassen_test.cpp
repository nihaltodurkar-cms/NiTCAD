// Klaassen (Philips unified) mobility (Unit 66, physics/klaassen_mobility.hpp): against an
// independent transcription of the formulas recorded in DECISIONS.md (Unit 66), the limits (KC2),
// Masetti's majority fits (KC3), minority above majority (KC4) and the partials (KC5).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "NiTCAD/physics/klaassen_mobility.hpp"
#include "NiTCAD/physics/mobility.hpp"

using namespace NiTCAD;
using physics::Carrier;

namespace {

// Written out from DECISIONS.md (Unit 66), not from the implementation: plain doubles, G's
// minimiser by a golden-section search on G itself (the implementation solves G' = 0).
struct Carr {
    double mu_max, mu_min, N_ref, alpha, theta, m;
};
constexpr Carr E{1414.0, 68.5, 9.20e16, 0.711, 2.285, 1.0};
constexpr Carr H{470.5, 44.9, 2.23e17, 0.719, 2.247, 1.258};
constexpr double s1 = 0.89233, s2 = 0.41372, s3 = 0.19778, s4 = 0.28227, s5 = 0.005978,
                 s6 = 1.80618, s7 = 0.72169;
constexpr double r1 = 0.7643, r2 = 2.2999, r3 = 6.5502, r4 = 2.3670, r5 = -0.01552, r6 = 0.6478;

double G_raw(double P, double m, double T) {
    return 1.0 - s1 / std::pow(s2 + std::pow(T / 300.0 / m, s4) * P, s3) +
           s5 / std::pow(std::pow(m * 300.0 / T, s7) * P, s6);
}

double G_ref(double P, double m, double T) {
    double a = std::log(1e-8), b = std::log(1e8);
    const double g = (std::sqrt(5.0) - 1.0) / 2.0;
    for (int k = 0; k < 200; ++k) {
        const double c = b - g * (b - a), d = a + g * (b - a);
        if (G_raw(std::exp(c), m, T) < G_raw(std::exp(d), m, T)) b = d;
        else a = c;
    }
    const double P_min = std::exp(0.5 * (a + b));
    return G_raw(std::max(P, P_min), m, T);
}

double cluster(double N, double c, double N_ref) {
    return N == 0.0 ? 0.0 : N * (1.0 + 1.0 / (c + (N_ref / N) * (N_ref / N)));
}

double reference(Carrier carrier, double ND, double NA, double n, double p, double T) {
    const bool el = carrier == Carrier::electron;
    const Carr& c = el ? E : H;
    const double mj = el ? H.m : E.m;
    const double Ds = cluster(ND, 0.21, 4.0e20), As = cluster(NA, 0.50, 7.2e20);
    const double same = el ? Ds : As, other = el ? As : Ds;  // own-type and opposite dopant
    const double screen = el ? p : n;                         // the opposite carrier
    const double Nsc = Ds + As + screen;
    const double inv_P = 2.459 * std::pow(Nsc, 2.0 / 3.0) / (3.97e13 * std::pow(T / 300.0, 2)) +
                         3.828 * (n + p) / (1.36e20 * c.m * std::pow(T / 300.0, 2));
    const double P = 1.0 / inv_P;
    const double ratio = c.m / mj;
    const double F = std::isinf(P) ? r1
                                   : (r1 * std::pow(P, r6) + r2 + r3 * ratio) /
                                         (std::pow(P, r6) + r4 + r5 * ratio);
    const double G = std::isinf(P) ? 1.0 : G_ref(P, c.m, T);
    const double Neff = same + G * other + screen / F;
    const double muL = c.mu_max * std::pow(T / 300.0, -c.theta);
    const double muN = c.mu_max * c.mu_max / (c.mu_max - c.mu_min) *
                       std::pow(T / 300.0, 3.0 * c.alpha - 1.5);
    const double muc = c.mu_max * c.mu_min / (c.mu_max - c.mu_min) * std::pow(300.0 / T, 0.5);
    if (Nsc == 0.0 && n + p == 0.0) return muL;
    if (Nsc == 0.0) return muL;  // N_sc,eff = 0 too: 1 / mu_DAP = 0 (the recorded limit)
    const double mu_dap = muN * (Nsc / Neff) * std::pow(c.N_ref / Nsc, c.alpha) +
                          muc * (n + p) / Neff;
    return 1.0 / (1.0 / muL + 1.0 / mu_dap);
}

constexpr double ni300 = 9.673e9;

}  // namespace

TEST_CASE("klaassen: the formulas as recorded", "[physics][klaassen]") {
    const auto& K = physics::klaassen_silicon;
    const double Ns[] = {0.0, 1e14, 1e16, 1e17, 1e18, 1e19, 1e20, 5e20};
    const double cs[] = {0.0, 1e10, 1e16, 1e18, 1e20};
    for (const Carrier carrier : {Carrier::electron, Carrier::hole}) {
        for (const double T : {250.0, 300.0, 400.0}) {
            for (const double ND : Ns) {
                for (const double NA : Ns) {
                    for (const double n : cs) {
                        for (const double p : cs) {
                            if (ND + NA + n + p == 0.0) continue;
                            const double mu =
                                physics::klaassen_mobility(K, carrier, ND, NA, n, p, T).mobility;
                            const double ref = reference(carrier, ND, NA, n, p, T);
                            CAPTURE(carrier == Carrier::electron, T, ND, NA, n, p, mu, ref);
                            REQUIRE(std::abs(mu / ref - 1.0) < 1e-10);
                        }
                    }
                }
            }
        }
    }
}

TEST_CASE("klaassen: the lattice limit (KC2)", "[physics][klaassen]") {
    const auto& K = physics::klaassen_silicon;
    for (const double T : {200.0, 300.0, 400.0}) {
        const double e = physics::klaassen_mobility(K, Carrier::electron, 0, 0, 0, 0, T).mobility;
        const double h = physics::klaassen_mobility(K, Carrier::hole, 0, 0, 0, 0, T).mobility;
        CHECK(std::abs(e / (1414.0 * std::pow(T / 300.0, -2.285)) - 1.0) < 1e-14);
        CHECK(std::abs(h / (470.5 * std::pow(T / 300.0, -2.247)) - 1.0) < 1e-14);
    }
    CHECK(physics::klaassen_mobility(K, Carrier::electron, 0, 0, 0, 0, 300.0).mobility == 1414.0);
    CHECK(physics::klaassen_mobility(K, Carrier::hole, 0, 0, 0, 0, 300.0).mobility == 470.5);
}

TEST_CASE("klaassen: majority carriers against Masetti (KC3)", "[physics][klaassen]") {
    const auto& K = physics::klaassen_silicon;
    double worst_e = 0.0, worst_h = 0.0;
    for (int d = 14; d <= 20; ++d) {
        const double N = std::pow(10.0, d);
        const double ke = physics::klaassen_mobility(K, Carrier::electron, N, 0, N,
                                                     ni300 * ni300 / N, 300.0).mobility;
        const double kh = physics::klaassen_mobility(K, Carrier::hole, 0, N, ni300 * ni300 / N,
                                                     N, 300.0).mobility;
        const double me = physics::masetti_mobility(physics::masetti_phosphorus_silicon, N, 300.0);
        const double mh = physics::masetti_mobility(physics::masetti_boron_silicon, N, 300.0);
        std::printf("KC3 N=%.0e  e: klaassen %.2f masetti %.2f (%+.2f%%)  h: klaassen %.2f "
                    "masetti %.2f (%+.2f%%)\n",
                    N, ke, me, 100.0 * (ke / me - 1.0), kh, mh, 100.0 * (kh / mh - 1.0));
        worst_e = std::max(worst_e, std::abs(ke / me - 1.0));
        worst_h = std::max(worst_h, std::abs(kh / mh - 1.0));
    }
    std::printf("KC3 largest deviation: electrons %.2f%%, holes %.2f%%\n", 100.0 * worst_e,
                100.0 * worst_h);
    CHECK(worst_e < 0.10);
    CHECK(worst_h < 0.10);
}

TEST_CASE("klaassen: minority above majority (KC4)", "[physics][klaassen]") {
    const auto& K = physics::klaassen_silicon;
    for (const double N : {1e19, 1e20}) {
        const double tiny = ni300 * ni300 / N;
        const double e_min = physics::klaassen_mobility(K, Carrier::electron, 0, N, tiny, N, 300.0)
                                 .mobility;
        const double e_maj = physics::klaassen_mobility(K, Carrier::electron, N, 0, N, tiny, 300.0)
                                 .mobility;
        const double h_min = physics::klaassen_mobility(K, Carrier::hole, N, 0, N, tiny, 300.0)
                                 .mobility;
        const double h_maj = physics::klaassen_mobility(K, Carrier::hole, 0, N, tiny, N, 300.0)
                                 .mobility;
        std::printf("KC4 N=%.0e  electrons minority/majority %.3f (%.1f / %.1f), holes %.3f "
                    "(%.1f / %.1f)\n",
                    N, e_min / e_maj, e_min, e_maj, h_min / h_maj, h_min, h_maj);
        CHECK(e_min > e_maj);
        CHECK(h_min > h_maj);
    }
}

TEST_CASE("klaassen: the partials (KC5)", "[physics][klaassen]") {
    const auto& K = physics::klaassen_silicon;
    const double Ns[] = {0.0, 1e15, 1e18, 1e20};
    const double cs[] = {0.0, 1e10, 1e16, 1e19};
    // Central differences with a relative step h = 1e-5 x: truncation about 1e-10 relative,
    // rounding about 1e-16 mu / 2e-5 = 5e-12 mu per unit relative change of x. Passing: the
    // partial to 1e-6 relative, or its effect x d mu / dx below 1e-9 mu (under the noise).
    const auto check = [](double analytic, double plus, double minus, double h, double mu,
                          double x) {
        const double fd = (plus - minus) / (2.0 * h);
        const double err = std::abs(analytic - fd) * x;
        return err <= 1e-6 * std::abs(fd) * x || err <= 1e-9 * mu;
    };
    for (const Carrier carrier : {Carrier::electron, Carrier::hole}) {
        for (const double T : {250.0, 300.0, 400.0}) {
            for (const double ND : Ns) {
                for (const double NA : Ns) {
                    for (const double n : cs) {
                        for (const double p : cs) {
                            // Not differentiable without dopants and screening carriers.
                            if (ND + NA == 0.0 && (n == 0.0 || p == 0.0)) continue;
                            const auto at = [&](double nn, double pp, double TT) {
                                return physics::klaassen_mobility(K, carrier, ND, NA, nn, pp, TT);
                            };
                            const auto m = at(n, p, T);
                            CAPTURE(carrier == Carrier::electron, T, ND, NA, n, p, m.mobility,
                                    m.d_dn, m.d_dp, m.d_dT);
                            if (n > 0.0) {
                                const double h = 1e-5 * n;
                                CHECK(check(m.d_dn, at(n + h, p, T).mobility,
                                            at(n - h, p, T).mobility, h, m.mobility, n));
                            }
                            if (p > 0.0) {
                                const double h = 1e-5 * p;
                                CHECK(check(m.d_dp, at(n, p + h, T).mobility,
                                            at(n, p - h, T).mobility, h, m.mobility, p));
                            }
                            const double h = 1e-5 * T;
                            CHECK(check(m.d_dT, at(n, p, T + h).mobility,
                                        at(n, p, T - h).mobility, h, m.mobility, T));
                        }
                    }
                }
            }
        }
    }
}

TEST_CASE("klaassen: the neutral equilibrium variant", "[physics][klaassen]") {
    const auto& K = physics::klaassen_silicon;
    for (const auto [ND, NA] : {std::pair{0.0, 0.0}, std::pair{1e17, 0.0}, std::pair{0.0, 1e19},
                                std::pair{1e18, 3e18}, std::pair{5e15, 5e15}}) {
        for (const Carrier carrier : {Carrier::electron, Carrier::hole}) {
            const double T = 320.0, ni = 3.0e10, dlnni = 0.07;
            const auto eq = physics::neutral_equilibrium(ND, NA, ni);
            CHECK(std::abs(eq.n * eq.p / (ni * ni) - 1.0) < 1e-12);
            CHECK(std::abs((eq.n - eq.p - (ND - NA)) / std::max({eq.n, eq.p, 1.0})) < 1e-12);
            const auto at = physics::klaassen_at(K, carrier, T);
            const auto m = physics::klaassen_equilibrium_mobility(at, ND, NA, ni, dlnni);
            CHECK(m.mobility ==
                  physics::klaassen_mobility(at, ND, NA, eq.n, eq.p).mobility);
            // Total d/dT with n_i(T) = ni exp(dlnni (T' - T)).
            const double h = 1e-5 * T;
            const auto val = [&](double TT) {
                const double nn = ni * std::exp(dlnni * (TT - T));
                return physics::klaassen_equilibrium_mobility(physics::klaassen_at(K, carrier, TT),
                                                              ND, NA, nn, dlnni)
                    .mobility;
            };
            const double fd = (val(T + h) - val(T - h)) / (2.0 * h);
            CAPTURE(ND, NA, carrier == Carrier::electron, m.d_dT, fd);
            CHECK(std::abs(m.d_dT - fd) <= 1e-6 * std::abs(fd) + 1e-12 * m.mobility);
        }
    }
}
