// TC7 (DECISIONS Unit 55 section 11.3, with the owner's readings of 2026-11-04 in 11.13): L0 tests of Chebyshev
// collocation with domain splitting.
//   - the frozen problem: -(eps u')' = f on [-1, 0] u [0, 1], eps = 1 left and 10 right, exact u = e^x (left) and
//     cos 3x + sin(3x) / 30 (right), so u(0) = 1 and eps u'(0) = 1 on both sides; N = 4, 8, 16, 32, 64 on each
//     subdomain. The error is the fine-grid interpolant error (201 points per subdomain, Clenshaw on the Chebyshev
//     coefficients; owner reading). Decay rule (owner reading): each doubling reduces the error by >= 10^3 or
//     brings it to the floor; once there it stays; the floor is reached by N = 64. The floor is 64 u_dd max|u| N^2
//     (owner ruling after measurement: the frozen 64 u_dd max|u| is below the rounding plateau of any correct
//     second-order collocation, whose matrices have condition O(N^4), Boyd section 7.7);
//   - interfaces (owner reading and ruling): at N = 64 the discrete mismatch |u_L(0) - u_R(0)| <= 64 u_dd max|u|
//     and |eps_L u_L'(0) - eps_R u_R'(0)| <= 64 u_dd max|eps u'| N^2; the errors of u(0) (floor as u) and of each
//     side's flux against the exact 1 (floor 64 u_dd max|eps u'| N^2) follow the decay rule;
//   - additions: cheb(1) ... cheb(5) against the values printed in Spectral Methods in MATLAB p. 55 (4 decimals);
//     the diagonal by (6.6) against Theorem 7's closed forms (6.3)-(6.4), within (N + 8) u_dd sum_j |D_ij|
//     (recursive summation of N terms, Higham's gamma_N, plus the entries' own rounding); exact differentiation of
//     x^k, k <= N, within (N + 2k + 8) u_dd sum_j |D_ij x_j^k| (the k products forming x^k added); these two
//     bounds were 8 u_dd before the first run, which ignored both terms (owner-approved correction);
//     the nodes against an independent cos; coefficients of T_k recovered; the dense solve; p13.m's problem
//     u_xx = e^{4x}, u(+-1) = 0 (exact (e^{4x} - x sinh 4 - cosh 4) / 16) with eps = 1 on both sides.
// cos and sin in the tests are test-local Taylor series (independent of the library's sin), |x| <= 3.2: the
// largest term is 3.2^3 / 6 < 5.5, so the cancellation costs at most a few u_dd.
#include "NiTCAD/verify/reference/chebyshev.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace ref = NiTCAD::verify::reference;
using ref::DD;
using ref::DDVector;
using ref::u_dd;

namespace {

DD cos_t(const DD& x) {
    const DD x2 = x * x;
    DD term(1.0), sum(1.0);
    for (int k = 1; k < 80; ++k) {
        term = -(term * x2) / DD((2.0 * k - 1.0) * (2.0 * k));
        sum = sum + term;
        if (std::abs(term.hi) < 0x1p-112) break;
    }
    return sum;
}
DD sin_t(const DD& x) {
    const DD x2 = x * x;
    DD term = x, sum = x;
    for (int k = 1; k < 80; ++k) {
        term = -(term * x2) / DD((2.0 * k) * (2.0 * k + 1.0));
        sum = sum + term;
        if (std::abs(term.hi) < 0x1p-112) break;
    }
    return sum;
}

const std::vector<int> ladder = {4, 8, 16, 32, 64};

// The condition-scaled floor (owner ruling after measurement): 64 u_dd * scale * N^2.
std::vector<double> floors(double scale) {
    std::vector<double> f;
    for (int n : ladder) f.push_back(64.0 * u_dd * scale * n * n);
    return f;
}

// The decay rule (owner reading): returns an empty string when it holds, else the first violation.
std::string decay_violation(const std::vector<double>& e, const std::vector<double>& floor) {
    for (std::size_t k = 0; k + 1 < e.size(); ++k) {
        const bool at_floor = e[k] <= floor[k];
        const bool ok = at_floor ? e[k + 1] <= floor[k + 1] : (e[k + 1] <= e[k] / 1e3 || e[k + 1] <= floor[k + 1]);
        if (!ok) return "pair " + std::to_string(k) + ": " + std::to_string(e[k]) + " -> " + std::to_string(e[k + 1]);
    }
    if (e.back() > floor.back()) return "floor not reached";
    return {};
}

struct Measured {
    std::vector<double> u_err, uc_err, fl_err, fr_err;
    double mismatch_u = 0, mismatch_flux = 0, tail = 0;
};

Measured run(const ref::TwoDomainProblem& base, const std::vector<int>& ns, DD (*exact)(const DD&),
             const DD& exact_uc, const DD& exact_flux) {
    Measured m;
    for (int n : ns) {
        ref::TwoDomainProblem p = base;
        p.n = n;
        const ref::TwoDomainSolution s = ref::solve_two_domain(p);
        REQUIRE(s.ok);
        double e = 0.0;
        for (int side = 0; side < 2; ++side) {
            const DD lo = side == 0 ? p.a : p.c, hi = side == 0 ? p.c : p.b;
            for (int i = 0; i <= 200; ++i) {
                DD x = lo + (hi - lo) * DD(i) / DD(200.0);
                if (i == 200) x = hi;
                // the right piece is evaluated just right of c through its own coefficients
                const DD v = (side == 1 && i == 0) ? ref::chebyshev_evaluate(s.coef_right, DD(-1.0)) : s.eval(x);
                e = std::max(e, std::abs((v - exact(x)).to_double()));
            }
        }
        m.u_err.push_back(e);
        m.uc_err.push_back(std::abs((s.u_left[0] - exact_uc).to_double()));
        m.fl_err.push_back(std::abs((s.flux_left - exact_flux).to_double()));
        m.fr_err.push_back(std::abs((s.flux_right - exact_flux).to_double()));
        m.mismatch_u = std::abs((s.u_left[0] - s.u_right.back()).to_double());
        m.mismatch_flux = std::abs((s.flux_left - s.flux_right).to_double());
        m.tail = s.tail_estimate;
    }
    return m;
}

std::string list(const std::vector<double>& v) {
    std::string s;
    for (double x : v) {
        char b[32];
        std::snprintf(b, sizeof b, "%.3e ", x);
        s += b;
    }
    return s;
}

DD exact_jump(const DD& x) {
    if (!(DD(0.0) < x)) return ref::exp(x);
    const DD t = DD(3.0) * x;
    return cos_t(t) + sin_t(t) / DD(30.0);
}
DD f_jump(const DD& x) {  // -(eps u')': left -e^x; right -10 u'' = 90 cos 3x + 3 sin 3x
    if (!(DD(0.0) < x)) return -ref::exp(x);
    const DD t = DD(3.0) * x;
    return DD(90.0) * cos_t(t) + DD(3.0) * sin_t(t);
}

DD exact_p13(const DD& x) {
    const DD e4 = ref::exp(DD(4.0)), em4 = ref::exp(DD(-4.0));
    const DD sh = ldexp(e4 - em4, -1), ch = ldexp(e4 + em4, -1);
    return (ref::exp(DD(4.0) * x) - x * sh - ch) / DD(16.0);
}

}  // namespace

TEST_CASE("TC7: cheb(1) ... cheb(5) match the values printed in Spectral Methods in MATLAB, p. 55") {
    const std::vector<std::vector<double>> printed = {
        {0.5, -0.5, 0.5, -0.5},
        {1.5, -2.0, 0.5, 0.5, -0.0, -0.5, -0.5, 2.0, -1.5},
        {3.1667, -4.0, 1.3333, -0.5, 1.0, -0.3333, -1.0, 0.3333, -0.3333, 1.0, 0.3333, -1.0, 0.5, -1.3333, 4.0,
         -3.1667},
        {5.5, -6.8284, 2.0, -1.1716, 0.5, 1.7071, -0.7071, -1.4142, 0.7071, -0.2929, -0.5, 1.4142, -0.0, -1.4142,
         0.5, 0.2929, -0.7071, 1.4142, 0.7071, -1.7071, -0.5, 1.1716, -2.0, 6.8284, -5.5},
        {8.5, -10.4721, 2.8944, -1.5279, 1.1056, -0.5, 2.6180, -1.1708, -2.0, 0.8944, -0.6180, 0.2764,
         -0.7236, 2.0, -0.1708, -1.6180, 0.8944, -0.3820, 0.3820, -0.8944, 1.6180, 0.1708, -2.0, 0.7236,
         -0.2764, 0.6180, -0.8944, 2.0, 1.1708, -2.6180, 0.5, -1.1056, 1.5279, -2.8944, 10.4721, -8.5}};
    double worst = 0.0;
    for (int n = 1; n <= 5; ++n) {
        const ref::DDMatrix d = ref::chebyshev_differentiation(n);
        for (std::size_t i = 0; i < d.v.size(); ++i)
            worst = std::max(worst, std::abs(d.v[i].to_double() - printed[static_cast<std::size_t>(n) - 1][i]));
    }
    INFO("worst difference from the printed 4-decimal values: " << worst);
    CHECK(worst <= 5.0001e-5);
}

TEST_CASE("TC7: nodes, the (6.6) diagonal against Theorem 7, exact differentiation of x^k") {
    double node = 0.0, diag = 0.0, poly = 0.0;
    for (int n : {4, 8, 16, 32, 64}) {
        const DDVector x = ref::chebyshev_points(n);
        const ref::DDMatrix d = ref::chebyshev_differentiation(n);
        const DD pin = ref::pi() / DD(static_cast<double>(n));
        for (int j = 0; j <= n; ++j) {
            const DD ang = pin * DD(static_cast<double>(j));
            node = std::max(node, std::abs((x[static_cast<std::size_t>(j)] - cos_t(ang)).to_double()) / u_dd);
            DD closed;
            if (j == 0 || j == n) {
                closed = DD(2.0 * n * n + 1.0) / DD(6.0);  // (6.3)
                if (j == n) closed = -closed;
            } else {
                const DD s = sin_t(ang);  // 1 - x_j^2 = sin^2(j pi / N), (6.4)
                closed = -cos_t(ang) / ldexp(s * s, 1);
            }
            DD scale(0.0);
            for (int k = 0; k <= n; ++k) scale = scale + abs(d(static_cast<std::size_t>(j), static_cast<std::size_t>(k)));
            diag = std::max(diag, std::abs((d(static_cast<std::size_t>(j), static_cast<std::size_t>(j)) - closed).to_double()) /
                                      (u_dd * scale.to_double() * (n + 8.0)));
        }
        for (int k = 0; k <= n; ++k) {
            for (int i = 0; i <= n; ++i) {
                DD s(0.0), scale(0.0);
                for (int j = 0; j <= n; ++j) {
                    DD xk(1.0);
                    for (int q = 0; q < k; ++q) xk = xk * x[static_cast<std::size_t>(j)];
                    const DD t = d(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) * xk;
                    s = s + t;
                    scale = scale + abs(t);
                }
                DD want(0.0);
                if (k > 0) {
                    want = DD(static_cast<double>(k));
                    for (int q = 0; q < k - 1; ++q) want = want * x[static_cast<std::size_t>(i)];
                }
                poly = std::max(poly, std::abs((s - want).to_double()) / (u_dd * std::max(1.0, scale.to_double()) * (n + 2.0 * k + 8.0)));
            }
        }
    }
    INFO("nodes " << node << " u_dd; diagonal " << diag << " of the (N + 8) u_dd sum|D_ij| bound; D x^k " << poly
                  << " of the (N + 2k + 8) u_dd sum|D_ij x_j^k| bound");
    CHECK(node <= 2.0);
    CHECK(diag <= 1.0);
    CHECK(poly <= 1.0);
}

TEST_CASE("TC7: Chebyshev coefficients of T_k recovered; Clenshaw; the dense solve") {
    double coef = 0.0, clen = 0.0;
    const int n = 16;
    const DDVector x = ref::chebyshev_points(n);
    for (int k = 0; k <= n; ++k) {
        DDVector v(x.size());
        // T_k(x_j) = cos(k j pi / N), by the test's own cos after exact reduction of k j mod 2N.
        for (std::size_t j = 0; j < x.size(); ++j) {
            long long r = (static_cast<long long>(k) * static_cast<long long>(j)) % (2LL * n);
            if (r > n) r = 2LL * n - r;
            v[j] = cos_t(ref::pi() * DD(static_cast<double>(r)) / DD(static_cast<double>(n)));
        }
        const DDVector a = ref::chebyshev_coefficients(v);
        for (int q = 0; q <= n; ++q)
            coef = std::max(coef, std::abs((a[static_cast<std::size_t>(q)] - DD(q == k ? 1.0 : 0.0)).to_double()) / u_dd);
    }
    // Clenshaw on T_5 at xi: 16 xi^5 - 20 xi^3 + 5 xi.
    DDVector a5(6, DD(0.0));
    a5[5] = DD(1.0);
    for (double xd : {-0.9, -0.3, 0.2, 0.77}) {
        const DD xi(xd), x2 = xi * xi;
        const DD want = ((DD(16.0) * x2 - DD(20.0)) * x2 + DD(5.0)) * xi;
        clen = std::max(clen, std::abs((ref::chebyshev_evaluate(a5, xi) - want).to_double()) / u_dd);
    }
    INFO("coefficients " << coef << " u_dd; Clenshaw T_5 " << clen << " u_dd");
    CHECK(coef <= 8.0);
    CHECK(clen <= 16.0);

    // A 6 x 6 system with integer entries and the known solution 1 ... 6, the pivot order exercised by a zero.
    ref::DDMatrix a(6, 6);
    const double m[6][6] = {{0, 2, 1, 0, 0, 3}, {4, 1, 0, 2, 0, 0}, {1, 0, 5, 1, 2, 0},
                            {0, 3, 1, 6, 0, 1}, {2, 0, 0, 1, 7, 2}, {1, 1, 1, 0, 2, 8}};
    DDVector b(6, DD(0.0));
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j) {
            a(static_cast<std::size_t>(i), static_cast<std::size_t>(j)) = DD(m[i][j]);
            b[static_cast<std::size_t>(i)] = b[static_cast<std::size_t>(i)] + DD(m[i][j] * (j + 1));
        }
    REQUIRE(ref::lu_solve(a, b));
    double sol = 0.0;
    for (int i = 0; i < 6; ++i) sol = std::max(sol, std::abs((b[static_cast<std::size_t>(i)] - DD(i + 1.0)).to_double()) / u_dd);
    INFO("dense solve: " << sol << " u_dd");
    CHECK(sol <= 64.0);
    ref::DDMatrix z(2, 2);
    DDVector bz(2, DD(1.0));
    CHECK_FALSE(ref::lu_solve(z, bz));
}

TEST_CASE("TC7: p13.m's problem u_xx = e^{4x}, u(+-1) = 0, split at 0 with eps = 1 on both sides") {
    ref::TwoDomainProblem p;
    p.f = [](const DD& x) { return -ref::exp(DD(4.0) * x); };
    p.u_a = DD(0.0);
    p.u_b = DD(0.0);
    const DD e4 = ref::exp(DD(4.0)), em4 = ref::exp(DD(-4.0));
    const DD flux = (DD(4.0) - ldexp(e4 - em4, -1)) / DD(16.0);  // u'(0) = (4 - sinh 4) / 16
    const Measured m = run(p, ladder, exact_p13, exact_p13(DD(0.0)), flux);
    double umax = 0.0;
    for (int i = 0; i <= 400; ++i) umax = std::max(umax, std::abs(exact_p13(DD(-1.0 + i / 200.0)).to_double()));
    const std::vector<double> floor = floors(umax);
    INFO("errors " << list(m.u_err) << "floors " << list(floor) << "; tail estimate at N = 64 " << m.tail);
    CHECK(decay_violation(m.u_err, floor).empty());
}

TEST_CASE("TC7: the eps-jump problem: geometric decay to 64 u_dd max|u| N^2, interface conditions") {
    ref::TwoDomainProblem p;
    p.eps_left = DD(1.0);
    p.eps_right = DD(10.0);
    p.f = f_jump;
    p.u_a = ref::exp(DD(-1.0));
    p.u_b = exact_jump(DD(1.0));
    const Measured m = run(p, ladder, exact_jump, DD(1.0), DD(1.0));
    double umax = 0.0, fmax = 0.0;
    for (int i = 0; i <= 400; ++i) {
        const DD x(-1.0 + i / 200.0);
        umax = std::max(umax, std::abs(exact_jump(x).to_double()));
        const DD flux = !(DD(0.0) < x) ? ref::exp(x)
                                       : DD(10.0) * (DD(-3.0) * sin_t(DD(3.0) * x) + cos_t(DD(3.0) * x) / DD(10.0));
        fmax = std::max(fmax, std::abs(flux.to_double()));
    }
    const std::vector<double> floor_u = floors(umax), floor_f = floors(fmax);
    INFO("u errors " << list(m.u_err) << "(floors " << list(floor_u) << ")");
    INFO("u(0) errors " << list(m.uc_err));
    INFO("left flux errors " << list(m.fl_err) << "; right " << list(m.fr_err) << "(floors " << list(floor_f) << ")");
    INFO("N = 64 mismatches: u " << m.mismatch_u << ", flux " << m.mismatch_flux << "; tail estimate " << m.tail);
    const std::string vu = decay_violation(m.u_err, floor_u), vc = decay_violation(m.uc_err, floor_u),
                      vl = decay_violation(m.fl_err, floor_f), vr = decay_violation(m.fr_err, floor_f);
    INFO("violations: u [" << vu << "] u(0) [" << vc << "] flux L [" << vl << "] flux R [" << vr << "]");
    CHECK(vu.empty());
    CHECK(vc.empty());
    CHECK(vl.empty());
    CHECK(vr.empty());
    CHECK(m.mismatch_u <= 64.0 * u_dd * umax);
    CHECK(m.mismatch_flux <= floor_f.back());
}

TEST_CASE("TC7: invalid problems are refused") {
    ref::TwoDomainProblem p;
    p.f = [](const DD&) { return DD(0.0); };
    p.n = 1;
    CHECK_FALSE(ref::solve_two_domain(p).ok);
    p.n = 8;
    p.c = DD(2.0);
    CHECK_FALSE(ref::solve_two_domain(p).ok);
    p.c = DD(0.0);
    p.f = nullptr;
    CHECK_FALSE(ref::solve_two_domain(p).ok);
}
