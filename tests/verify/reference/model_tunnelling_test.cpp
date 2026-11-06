// TVR-3 step 8 (DECISIONS Unit 55 15.9, 15.18): L0 tests of the tunnelling references of TVB-31.23-31.25
// (verify/reference/model_tunnelling.hpp).
//   - The segment action against tanh-sinh quadrature (TC4) of sqrt((a + (b - a) t)+) over [0, 1], split at the
//     turning point; its a = b limit sqrt(a); its partials against double-double central differences.
//   - The supply's partials against central differences; exactly 0 at F1 = F2; the Boltzmann limit kT e^((F1 - E)/kT)
//     well above both levels.
//   - The current's Gauss-Kronrod integral against tanh-sinh over the same pieces.
//   - The action scale against another arrangement of its constants; a rectangular barrier's S and TC against
//     kappa L sqrt(V - E) and e^(-2 S) (added after mutants R01 and R03 survived).
#include "NiTCAD/verify/reference/model_tunnelling.hpp"

#include "NiTCAD/verify/reference/tanh_sinh.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <vector>

namespace ref = NiTCAD::verify::reference;
namespace mt = NiTCAD::verify::reference::model_tunnelling;
using ref::DD;

TEST_CASE("model_tunnelling: segment action", "[verify][model_tunnelling]") {
    const DD L(1e-7);
    const std::vector<std::array<double, 2>> pairs{
        {1.0, 2.0}, {2.0, 1.0}, {3.0, 3.0}, {0.5, 0.5 + 1e-9}, {1.0, -1.0}, {-0.5, 2.0}, {1e-6, 3.0}, {3.0, 0.0},
        {0.0, 2.0}, {0.2, -0.0001}, {-2.0, -1.0}, {5.0, 1e-3}};
    for (const auto& p : pairs) {
        const DD a(p[0]), b(p[1]);
        const auto g = mt::segment_action(a, b, L);
        auto root = [&](const DD& t) {
            const DD v = a + (b - a) * t;
            return v.hi > 0.0 ? sqrt(v) : DD(0.0);
        };
        DD q(0.0);
        if (p[0] == p[1]) {
            q = p[0] > 0.0 ? sqrt(a) : DD(0.0);
        } else if (p[0] > 0.0 || p[1] > 0.0) {
            // Only the part inside the barrier, ending at the turning point t* = a / (a - b) (in double-double).
            const DD tstar = a / (a - b);
            const bool turning = (p[0] > 0.0) != (p[1] > 0.0);
            const DD lo = !turning || p[0] > 0.0 ? DD(0.0) : tstar, hi = !turning || p[1] > 0.0 ? DD(1.0) : tstar;
            const auto r = ref::tanh_sinh(root, lo, hi, 0.0, 1e-28);
            REQUIRE(r.status == ref::TanhSinhStatus::ok);
            q = r.value;
        }
        INFO("a " << p[0] << ", b " << p[1] << ": " << g[0].hi << " vs " << (L * q).hi);
        CHECK(std::abs((g[0] - L * q).hi) <= 1e-26 * std::abs(g[0].hi) + 1e-40);
        // Partials, away from a turning point.
        for (int k = 0; k < 2; ++k) {
            const double v = p[k];
            if (std::abs(v) < 1e-3 || (p[0] > 0.0) != (p[1] > 0.0) && std::abs(v) < 0.1) continue;
            const DD h = DD(std::abs(v) * 0x1p-30);
            const auto up = k == 0 ? mt::segment_action(a + h, b, L) : mt::segment_action(a, b + h, L);
            const auto dn = k == 0 ? mt::segment_action(a - h, b, L) : mt::segment_action(a, b - h, L);
            const DD cd = (up[0] - dn[0]) / (DD(2.0) * h);
            INFO("partial " << k << ": " << g[1 + k].hi << " vs " << cd.hi);
            CHECK(std::abs((g[1 + k] - cd).hi) <= 1e-15 * (std::abs(g[1 + k].hi) + std::abs(g[0].hi / v)));
        }
    }
}

TEST_CASE("model_tunnelling: supply", "[verify][model_tunnelling]") {
    const DD kT(0.025852);
    for (const double E : {-1.0, 0.0, 0.3, 1.0, 2.5}) {
        for (const double dF : {0.0, 1e-6, 0.1, 1.0}) {
            const std::array<DD, 4> x{DD(E), DD(0.2), DD(0.2 - dF), kT};
            const auto s = mt::supply(x);
            if (dF == 0.0) CHECK(s[0].hi == 0.0);
            for (int i = 0; i < 3; ++i) {
                const DD h(0x1p-30);
                auto xu = x, xd = x;
                xu[i] = xu[i] + h;
                xd[i] = xd[i] - h;
                const DD cd = (mt::supply(xu)[0] - mt::supply(xd)[0]) / (DD(2.0) * h);
                const DD an = i == 0 ? s[3] : s[i];
                INFO("E " << E << " dF " << dF << " input " << i << ": " << an.hi << " vs " << cd.hi);
                CHECK(std::abs((an - cd).hi) <= 1e-15);
            }
        }
    }
    // Boltzmann limit 40 kT above both levels: N = kT (e^x1 - e^x2) to e^(2 x1) relative.
    const std::array<DD, 4> x{DD(0.2 + 40 * 0.025852), DD(0.2), DD(0.1), kT};
    const DD x1 = (x[1] - x[0]) / kT, x2 = (x[2] - x[0]) / kT;
    const DD boltz = kT * (exp(x1) - exp(x2));
    CHECK(std::abs(((mt::supply(x)[0] - boltz) / boltz).hi) <= 1e-17);
}

TEST_CASE("model_tunnelling: current", "[verify][model_tunnelling]") {
    // A 1.5 nm SiO2-like barrier, 3.1 eV high at its left end, 1 V across it.
    const std::vector<mt::Segment> barrier{{DD(3.1), DD(2.1), DD(1.5e-7), DD(0.5)}};
    const auto J = mt::current(barrier, 0.0, INFINITY, 0.1, -0.9, 300.0, 120.0);
    REQUIRE(J.ok);
    const DD kT = DD(NiTCAD::base::k_B_eV_per_K) * DD(300.0);
    auto f = [&](const DD& E) {
        return mt::transmission(barrier, E) * mt::supply(std::array<DD, 4>{E, DD(0.1), DD(-0.9), kT})[0];
    };
    const double Eh = (DD(0.1) + DD(40.0) * kT).hi;
    const auto t = ref::tanh_sinh(f, DD(0.0), DD(Eh), 0.0, 1e-26);
    REQUIRE(t.status == ref::TanhSinhStatus::ok);
    const DD k(NiTCAD::base::k_B_eV_per_K);
    const DD Jt = DD(120.0) / (k * k) * t.value;
    INFO(J.J.hi << " vs " << Jt.hi);
    CHECK(std::abs(((J.J - Jt) / Jt).hi) <= 1e-22);
    CHECK(J.J.hi > 0.0);
}

TEST_CASE("model_tunnelling: scale and rectangular transmission", "[verify][model_tunnelling]") {
    // The scale, as sqrt(2 m0 q / hbar^2) per metre over 100 (another arrangement of the same constants), and its
    // magnitude 5.1232e7 / (cm eV^(1/2)) from CODATA 2018.
    const DD m0(NiTCAD::base::m0_kg), q(NiTCAD::base::q_C), hbar(NiTCAD::base::hbar_J_s);
    const DD scale = mt::action_scale();
    const DD alt = sqrt(DD(2.0) * m0 * q / (hbar * hbar)) / DD(100.0);
    CHECK(std::abs(((scale - alt) / alt).hi) <= 1e-30);
    CHECK(std::abs(scale.hi / 5.1232e7 - 1.0) <= 1e-4);
    // A rectangular barrier: S = kappa L sqrt(V - E) with kappa = scale sqrt(m), and TC = e^(-2 S).
    const DD V(3.1), L(1.5e-7), m(0.5);
    const std::vector<mt::Segment> barrier{{V, V, L, m}};
    for (const double E : {-1.0, 0.0, 1.0, 3.0}) {
        const DD S = scale * sqrt(m) * L * sqrt(V - DD(E));
        const DD TC = exp(-(S + S));
        INFO("E " << E << ": " << mt::transmission(barrier, DD(E)).hi << " vs " << TC.hi);
        CHECK(std::abs(((mt::exponent(barrier, DD(E))[0] - S) / S).hi) <= 1e-30);
        CHECK(std::abs(((mt::transmission(barrier, DD(E)) - TC) / TC).hi) <= 1e-28);
    }
}
