// TC6 (DECISIONS Unit 55 section 11.3, with the owner's readings of 2026-11-04 in 11.13): L0 tests of the
// Dormand-Prince 5(4) pair with Shampine's dense output.
//   - the tableau "to printed digits": every entry against the rationals of the second reproduction (Boost.Odeint,
//     typed here; the library's are typed from dopri5.f) to 1 u_dd; the error weights e = b - bhat against Boost's
//     bhat; Boost's dense output, which is a different order-4 extension (see its test);
//   - observed order 5.0 +- 0.2 at fixed steps, N = 128, 256, 512, 1024, on y' = y and the harmonic oscillator (and,
//     an addition approved by the owner, the nonlinear y' = y^2, since linear problems see only the tall-tree
//     conditions); the ladder was N = 16 ... 128 until measured (owner ruling, DECISIONS 11.13);
//   - adaptive global error proportional to the tolerance over four decades: the least-squares slope of
//     log err against log tol over tol = 1e-14 ... 1e-18 in [0.96, 1.04] (the order band divided by 5, h ~ tol^(1/5);
//     owner reading; the range was 1e-10 ... 1e-14 until measured, owner ruling);
//   - dense output order 4 read as the continuous extension's order (HNW): one step from exact data, the error at
//     theta in {1/4, 1/2, 3/4} falls as h^(4.0 +- 0.2 + 1) (owner reading); the global dense error is printed only;
//   - additions: row sums, the 17 order conditions of order <= 5 for b, the 8 of order <= 4 for bhat and for the
//     dense weights b(theta), bhat failing order 5 (so its order is real), input refusal, the step limit, landing
//     on tend, backward integration and the evaluation count.
// The oscillator's exact solution uses test-local cos and sin series on [0, 2 pi] (terms up to 2 pi^6 / 6! ~ 85, so
// a few hundred u_dd of cancellation, far below the errors measured, which are >= 1e-26).
#include "NiTCAD/verify/reference/dopri5.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

namespace ref = NiTCAD::verify::reference;
using ref::DD;
using ref::DDVector;
using ref::u_dd;

namespace {

DD q(double n, double d) { return DD(n) / DD(d); }

DD cos_series(const DD& x) {
    const DD x2 = x * x;
    DD term(1.0), sum(1.0);
    for (int k = 1; k < 80; ++k) {
        term = -(term * x2) / DD((2.0 * k - 1.0) * (2.0 * k));
        sum = sum + term;
        if (std::abs(term.hi) < 0x1p-112) break;
    }
    return sum;
}
DD sin_series(const DD& x) {
    const DD x2 = x * x;
    DD term = x, sum = x;
    for (int k = 1; k < 80; ++k) {
        term = -(term * x2) / DD((2.0 * k) * (2.0 * k + 1.0));
        sum = sum + term;
        if (std::abs(term.hi) < 0x1p-112) break;
    }
    return sum;
}


struct Problem {
    const char* name;
    DD t_end;
    DDVector y0;
    void (*f)(const DD&, const DDVector&, DDVector&);
    DDVector (*exact)(const DD&);
};

void f_exp(const DD&, const DDVector& y, DDVector& dy) { dy[0] = y[0]; }
DDVector x_exp(const DD& t) { return {ref::exp(t)}; }
void f_osc(const DD&, const DDVector& y, DDVector& dy) {
    dy[0] = y[1];
    dy[1] = -y[0];
}
DDVector x_osc(const DD& t) { return {cos_series(t), -sin_series(t)}; }
void f_sq(const DD&, const DDVector& y, DDVector& dy) { dy[0] = y[0] * y[0]; }
DDVector x_sq(const DD& t) { return {DD(1.0) / (DD(2.0) - t)}; }

std::vector<Problem> problems() {
    return {{"y' = y on [0, 1]", DD(1.0), {DD(1.0)}, f_exp, x_exp},
            {"harmonic oscillator on [0, 2 pi]", ref::ldexp(ref::pi(), 1), {DD(1.0), DD(0.0)}, f_osc, x_osc},
            {"y' = y^2 on [0, 1]", DD(1.0), {DD(0.5)}, f_sq, x_sq}};
}

double max_error(const DDVector& a, const DDVector& b) {
    double e = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) e = std::max(e, std::abs((a[i] - b[i]).to_double()));
    return e;
}

// Residuals of the order conditions for weights w (sum w_i Phi_i(tree) - theta^rho / gamma), trees of order
// <= max_order (Butcher; HNW I, section II.2), rho the tree's order.
std::vector<double> order_residuals(const std::array<DD, 7>& w, const DD& theta, int max_order) {
    const ref::Dopri5Tableau& t = ref::dopri5_tableau();
    using V = std::array<DD, 7>;
    auto A = [&](const V& v) {
        V r{};
        for (int i = 0; i < 7; ++i)
            for (int j = 0; j < i; ++j) r[i] = r[i] + t.a[i][j] * v[j];
        return r;
    };
    auto mul = [](const V& u, const V& v) {
        V r{};
        for (int i = 0; i < 7; ++i) r[i] = u[i] * v[i];
        return r;
    };
    auto dot = [&](const V& v) {
        DD s(0.0);
        for (int i = 0; i < 7; ++i) s = s + w[i] * v[i];
        return s;
    };
    V one;
    one.fill(DD(1.0));
    const V c = t.c, c2 = mul(c, c), c3 = mul(c2, c), c4 = mul(c3, c);
    const V Ac = A(c), Ac2 = A(c2), AAc = A(Ac);
    struct Tree {
        int order;
        V phi;
        double gamma;
    };
    std::vector<Tree> trees = {{1, one, 1},           {2, c, 2},           {3, c2, 3},
                               {3, Ac, 6},            {4, c3, 4},          {4, mul(c, Ac), 8},
                               {4, Ac2, 12},          {4, AAc, 24},        {5, c4, 5},
                               {5, mul(c2, Ac), 10},  {5, mul(c, Ac2), 15}, {5, mul(c, AAc), 30},
                               {5, mul(Ac, Ac), 20},  {5, A(c3), 20},      {5, A(mul(c, Ac)), 40},
                               {5, A(Ac2), 60},       {5, A(AAc), 120}};
    std::vector<double> out;
    for (const Tree& tr : trees) {
        if (tr.order > max_order) continue;
        DD th(1.0);
        for (int k = 0; k < tr.order; ++k) th = th * theta;
        out.push_back(std::abs((dot(tr.phi) - th / DD(tr.gamma)).to_double()));
    }
    return out;
}

double max_of(const std::vector<double>& v) {
    double m = 0.0;
    for (double x : v) m = std::max(m, x);
    return m;
}

// Shampine's dense weights b_i(theta), from the CONT form of dopri5.f (CONTD5).
std::array<DD, 7> dense_weights(const DD& theta) {
    const ref::Dopri5Tableau& t = ref::dopri5_tableau();
    const DD th1 = DD(1.0) - theta;
    std::array<DD, 7> w{};
    for (int i = 0; i < 7; ++i) {
        const DD delta1(i == 0 ? 1.0 : 0.0), delta7(i == 6 ? 1.0 : 0.0);
        const DD bspl = delta1 - t.b[i];
        const DD c3 = -delta7 + t.b[i] - bspl;
        w[i] = theta * t.b[i] + theta * th1 * bspl + theta * theta * th1 * c3 + theta * theta * th1 * th1 * t.d[i];
    }
    return w;
}

// One step from exact data at t = 0: the dense error, max over theta in {1/4, 1/2, 3/4}.
double local_dense_error(const Problem& p, const DD& h) {
    ref::Dopri5Step s;
    s.k[0].assign(p.y0.size(), DD(0.0));
    p.f(DD(0.0), p.y0, s.k[0]);
    ref::dopri5_stages(p.f, DD(0.0), p.y0, h, s);
    const ref::Dopri5Dense dn = ref::dopri5_dense(DD(0.0), p.y0, h, s);
    double e = 0.0;
    for (double th : {0.25, 0.5, 0.75}) {
        const DD t = DD(th) * h;
        e = std::max(e, max_error(dn.eval(t), p.exact(t)));
    }
    return e;
}

}  // namespace

TEST_CASE("TC6: the tableau agrees with the second reproduction (Boost.Odeint) to 1 u_dd") {
    const ref::Dopri5Tableau& t = ref::dopri5_tableau();
    auto rel = [](const DD& x, const DD& y) {
        const double s = std::max(std::abs(y.to_double()), 1e-300);
        return std::abs((x - y).to_double()) / s / u_dd;
    };
    // Boost's runge_kutta_dopri5: a2..a5 (the nodes), b21..b65 (the matrix), c1..c6 (the weights), and the 4th-order
    // weights c - dc.
    const DD nodes[7] = {DD(0.0), q(1, 5), q(3, 10), q(4, 5), q(8, 9), DD(1.0), DD(1.0)};
    const std::vector<std::vector<DD>> rows = {
        {},
        {q(1, 5)},
        {q(3, 40), q(9, 40)},
        {q(44, 45), q(-56, 15), q(32, 9)},
        {q(19372, 6561), q(-25360, 2187), q(64448, 6561), q(-212, 729)},
        {q(9017, 3168), q(-355, 33), q(46732, 5247), q(49, 176), q(-5103, 18656)},
        {q(35, 384), DD(0.0), q(500, 1113), q(125, 192), q(-2187, 6784), q(11, 84)}};
    const DD bhat[7] = {q(5179, 57600), DD(0.0), q(7571, 16695), q(393, 640), q(-92097, 339200), q(187, 2100),
                        q(1, 40)};
    double worst = 0.0;
    for (int i = 0; i < 7; ++i) {
        if (nodes[i].hi != 0.0) worst = std::max(worst, rel(t.c[i], nodes[i]));
        for (int j = 0; j < 7; ++j) {
            const DD boost = j < static_cast<int>(rows[i].size()) ? rows[i][j] : DD(0.0);
            if (boost.hi == 0.0)
                CHECK(t.a[i][j] == DD(0.0));
            else
                worst = std::max(worst, rel(t.a[i][j], boost));
        }
        CHECK(t.b[i] == t.a[6][i]);
    }
    INFO("worst c, a: " << worst << " u_dd");
    CHECK(worst <= 1.0);
    double worst_e = 0.0;
    for (int i = 0; i < 7; ++i)
        worst_e = std::max(worst_e, std::abs((t.e[i] - (t.b[i] - bhat[i])).to_double()) / u_dd);
    INFO("e = b - bhat: worst " << worst_e << " u_dd absolute");
    CHECK(worst_e <= 2.0);
    CHECK(t.d[1] == DD(0.0));
}

TEST_CASE("TC6: row sums, and the order conditions of b (5), bhat (4, not 5) and b(theta) (4)") {
    const ref::Dopri5Tableau& t = ref::dopri5_tableau();
    double worst_row = 0.0;
    for (int i = 1; i < 7; ++i) {
        DD s(0.0);
        for (int j = 0; j < i; ++j) s = s + t.a[i][j];
        worst_row = std::max(worst_row, std::abs((s - t.c[i]).to_double()) / u_dd);
    }
    INFO("row sums: " << worst_row << " u_dd");
    CHECK(worst_row <= 2.0);

    const double b5 = max_of(order_residuals(t.b, DD(1.0), 5)) / u_dd;
    INFO("b, 17 conditions of order <= 5: " << b5 << " u_dd");
    CHECK(b5 <= 4.0);

    std::array<DD, 7> bhat{};
    for (int i = 0; i < 7; ++i) bhat[i] = t.b[i] - t.e[i];
    const double bh4 = max_of(order_residuals(bhat, DD(1.0), 4)) / u_dd;
    const double bh5 = max_of(order_residuals(bhat, DD(1.0), 5));
    INFO("bhat: order <= 4 " << bh4 << " u_dd; order 5 " << bh5 << " absolute");
    CHECK(bh4 <= 4.0);
    CHECK(bh5 > 1e-6);

    double dense4 = 0.0;
    for (int k = 1; k <= 10; ++k) dense4 = std::max(dense4, max_of(order_residuals(dense_weights(DD(k / 10.0)), DD(k / 10.0), 4)));
    INFO("b(theta), 8 conditions of order <= 4, theta = 0.1 ... 1: " << dense4 / u_dd << " u_dd");
    CHECK(dense4 / u_dd <= 8.0);
    double b1 = 0.0;
    const std::array<DD, 7> w1 = dense_weights(DD(1.0));
    for (int i = 0; i < 7; ++i) b1 = std::max(b1, std::abs((w1[i] - t.b[i]).to_double()) / u_dd);
    INFO("b(1) - b: " << b1 << " u_dd");
    CHECK(b1 <= 2.0);
}

TEST_CASE("TC6: the dense weights are a Hermite interpolant: b'(0) = e_1, b'(1) = e_7") {
    // d/dtheta of dense_weights: b + (1 - 2 theta) bspl + (2 theta - 3 theta^2) c3 + 2 theta (1 - theta)(1 - 2 theta) d.
    const ref::Dopri5Tableau& t = ref::dopri5_tableau();
    double worst = 0.0;
    for (double th : {0.0, 1.0}) {
        for (int i = 0; i < 7; ++i) {
            const DD delta1(i == 0 ? 1.0 : 0.0), delta7(i == 6 ? 1.0 : 0.0);
            const DD bspl = delta1 - t.b[i];
            const DD c3 = -delta7 + t.b[i] - bspl;
            const DD x(th);
            const DD dw = t.b[i] + (DD(1.0) - DD(2.0) * x) * bspl + (DD(2.0) * x - DD(3.0) * x * x) * c3 +
                          DD(2.0) * x * (DD(1.0) - x) * (DD(1.0) - DD(2.0) * x) * t.d[i];
            const DD want = th == 0.0 ? delta1 : delta7;
            worst = std::max(worst, std::abs((dw - want).to_double()) / u_dd);
        }
    }
    INFO("b'(0) - e_1, b'(1) - e_7: " << worst << " u_dd");
    CHECK(worst <= 2.0);
}

// The two reproductions print different interpolants (owner ruling 2026-11-04, DECISIONS 11.13): dopri5.f's
// (degree 4 in theta, "Shampine (1986)") and Boost's (degree 5, cited there to HNW pp. 191-192). Both are
// continuous extensions of order 4; they coincide at theta = 0, 1/2 and 1 and differ elsewhere by about
// 1e-4 theta^2 (1 - theta)^2 (theta - 1/2) per weight. The test asserts exactly that.
TEST_CASE("TC6: Boost.Odeint's dense output is a different order-4 extension, equal at theta = 0, 1/2, 1") {
    // Boost's calc_state: b_i(theta) = A b_i -+ C X_i(theta) (+ D for i = 1), b_7(theta) = B + C X_7(theta).
    const DD b1 = q(35, 384), b3 = q(500, 1113), b4 = q(125, 192), b5 = q(-2187, 6784), b6 = q(11, 84);
    auto boost_weights = [&](const DD& th) {
        const DD X1 = DD(5.0) * (DD(2558722523.0) - DD(31403016.0) * th) / DD(11282082432.0);
        const DD X3 = DD(100.0) * (DD(882725551.0) - DD(15701508.0) * th) / DD(32700410799.0);
        const DD X4 = DD(25.0) * (DD(443332067.0) - DD(31403016.0) * th) / DD(1880347072.0);
        const DD X5 = DD(32805.0) * (DD(23143187.0) - DD(3489224.0) * th) / DD(199316789632.0);
        const DD X6 = DD(55.0) * (DD(29972135.0) - DD(7076736.0) * th) / DD(822651844.0);
        const DD X7 = DD(10.0) * (DD(7414447.0) - DD(829305.0) * th) / DD(29380423.0);
        const DD tm1 = th - DD(1.0), th2 = th * th;
        const DD A = th2 * (DD(3.0) - DD(2.0) * th), B = th2 * tm1, C = th2 * tm1 * tm1, D = th * tm1 * tm1;
        return std::array<DD, 7>{A * b1 - C * X1 + D, DD(0.0), A * b3 + C * X3, A * b4 - C * X4,
                                 A * b5 + C * X5,     A * b6 - C * X6, B + C * X7};
    };
    double cond4 = 0.0, equal = 0.0, differ = 1.0;
    for (int k = 0; k <= 16; ++k) {
        const DD th(k / 16.0);
        const std::array<DD, 7> wb = boost_weights(th), wf = dense_weights(th);
        if (k > 0) cond4 = std::max(cond4, max_of(order_residuals(wb, th, 4)) / u_dd);
        double d = 0.0;
        for (int i = 0; i < 7; ++i) d = std::max(d, std::abs((wb[i] - wf[i]).to_double()));
        if (k == 0 || k == 8 || k == 16)
            equal = std::max(equal, d / u_dd);
        else
            differ = std::min(differ, d);
    }
    INFO("Boost: order <= 4 conditions " << cond4 << " u_dd; difference at 0, 1/2, 1 " << equal
                                         << " u_dd; smallest difference elsewhere " << differ);
    CHECK(cond4 <= 8.0);
    CHECK(equal <= 2.0);
    CHECK(differ > 1e-8);
}

TEST_CASE("TC6: the dense output reproduces the step's end values") {
    const Problem p = problems()[2];
    const DD h(0.125);
    ref::Dopri5Step s;
    s.k[0].assign(1, DD(0.0));
    p.f(DD(0.0), p.y0, s.k[0]);
    ref::dopri5_stages(p.f, DD(0.0), p.y0, h, s);
    const ref::Dopri5Dense dn = ref::dopri5_dense(DD(0.0), p.y0, h, s);
    CHECK(dn.eval(DD(0.0))[0] == p.y0[0]);
    CHECK(std::abs((dn.eval(h)[0] - s.y1[0]).to_double()) <= 2.0 * u_dd);
}

TEST_CASE("TC6: observed order 5.0 +- 0.2 at fixed steps") {
    for (const Problem& p : problems()) {
        std::vector<double> err;
        for (long long n : {128, 256, 512, 1024})
            err.push_back(max_error(ref::dopri5_fixed(p.f, DD(0.0), p.y0, p.t_end, n), p.exact(p.t_end)));
        for (std::size_t k = 0; k + 1 < err.size(); ++k) {
            const double order = std::log2(err[k] / err[k + 1]);
            INFO(p.name << ": N = " << (128 << k) << " -> " << (256 << k) << ", errors " << err[k] << ", "
                        << err[k + 1] << ", order " << order);
            CHECK(order >= 4.8);
            CHECK(order <= 5.2);
        }
    }
}

TEST_CASE("TC6: dense output, continuous order 4.0 +- 0.2 (local error ~ h^5 at theta = 1/4, 1/2, 3/4)") {
    for (const Problem& p : problems()) {
        std::vector<double> err;
        for (int k = 3; k <= 7; ++k) err.push_back(local_dense_error(p, DD(std::ldexp(1.0, -k))));
        for (std::size_t k = 0; k + 1 < err.size(); ++k) {
            const double order = std::log2(err[k] / err[k + 1]) - 1.0;
            INFO(p.name << ": h = 2^-" << (k + 3) << " -> 2^-" << (k + 4) << ", errors " << err[k] << ", "
                        << err[k + 1] << ", continuous order " << order);
            CHECK(order >= 3.8);
            CHECK(order <= 4.2);
        }
        // The global dense error over a fixed-step run, printed only (no criterion; it reads about 5).
        std::vector<double> g;
        for (long long n : {16, 32, 64}) {
            double e = 0.0;
            (void)ref::dopri5_fixed(p.f, DD(0.0), p.y0, p.t_end, n, [&](const ref::Dopri5Dense& d) {
                const DD t = d.t_old + ref::ldexp(d.h, -1);
                e = std::max(e, max_error(d.eval(t), p.exact(t)));
            });
            g.push_back(e);
        }
        INFO(p.name << ": global dense error at mid-steps, N = 16, 32, 64: " << g[0] << ", " << g[1] << ", "
                    << g[2] << " (orders " << std::log2(g[0] / g[1]) << ", " << std::log2(g[1] / g[2]) << ")");
        CHECK(g[2] < g[0]);
    }
}

TEST_CASE("TC6: adaptive global error proportional to the tolerance over four decades") {
    for (const Problem& p : problems()) {
        std::vector<double> lt, le;
        for (int k = 14; k <= 18; ++k) {
            ref::Dopri5Options o;
            o.rtol = o.atol = std::pow(10.0, -k);
            const ref::Dopri5Result r = ref::dopri5(p.f, DD(0.0), p.y0, p.t_end, o);
            REQUIRE(r.status == ref::Dopri5Status::ok);
            CHECK(r.t == p.t_end);
            CHECK(r.nfcn == 6 * r.nstep + 2);
            const double e = max_error(r.y, p.exact(p.t_end));
            INFO(p.name << ": tol 1e-" << k << " error " << e << " in " << r.naccpt << " accepted, " << r.nrejct
                        << " rejected");
            CHECK(e > 0.0);
            lt.push_back(std::log10(o.rtol));
            le.push_back(std::log10(e));
        }
        const double n = static_cast<double>(lt.size());
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        for (std::size_t i = 0; i < lt.size(); ++i) {
            sx += lt[i];
            sy += le[i];
            sxx += lt[i] * lt[i];
            sxy += lt[i] * le[i];
        }
        const double slope = (n * sxy - sx * sy) / (n * sxx - sx * sx);
        INFO(p.name << ": slope of log err against log tol " << slope);
        CHECK(slope >= 0.96);
        CHECK(slope <= 1.04);
    }
}

TEST_CASE("TC6: input refusal, the step limit, backward integration and the dense output on accepted steps") {
    const Problem p = problems()[0];
    ref::Dopri5Options o;
    o.rtol = o.atol = 0.0;
    CHECK(ref::dopri5(p.f, DD(0.0), p.y0, DD(1.0), o).status == ref::Dopri5Status::invalid_input);
    o = {};
    CHECK(ref::dopri5(p.f, DD(1.0), p.y0, DD(1.0), o).status == ref::Dopri5Status::invalid_input);
    o.nmax = 3;
    o.rtol = o.atol = 1e-14;
    CHECK(ref::dopri5(p.f, DD(0.0), p.y0, DD(1.0), o).status == ref::Dopri5Status::too_many_steps);

    // Backward: y' = y from t = 1, y = e, to t = 0.
    o = {};
    o.rtol = o.atol = 1e-12;
    const ref::Dopri5Result r = ref::dopri5(p.f, DD(1.0), {ref::exp(DD(1.0))}, DD(0.0), o);
    REQUIRE(r.status == ref::Dopri5Status::ok);
    CHECK(r.t == DD(0.0));
    INFO("backward error " << std::abs((r.y[0] - DD(1.0)).to_double()));
    CHECK(std::abs((r.y[0] - DD(1.0)).to_double()) < 1e-10);

    // The dense output covers [0, 1] without gaps and matches the step ends.
    DD covered(0.0);
    long long calls = 0;
    double worst = 0.0;
    const ref::Dopri5Result a = ref::dopri5(p.f, DD(0.0), p.y0, DD(1.0), o, [&](const ref::Dopri5Dense& d) {
        CHECK(d.t_old == covered);
        covered = d.t_old + d.h;
        ++calls;
        const DD mid = d.t_old + ref::ldexp(d.h, -1);
        worst = std::max(worst, std::abs((d.eval(mid)[0] - ref::exp(mid)).to_double()));
    });
    CHECK(calls == a.naccpt);
    CHECK(std::abs((covered - DD(1.0)).to_double()) <= u_dd);
    INFO("dense error at mid-steps, tol 1e-12: " << worst);
    CHECK(worst < 1e-10);
}

// Additions (no frozen criterion): the stage times c_i h enter only through f(t, y), which no autonomous problem
// sees; the quadrature y' = 5 t^4 is integrated exactly by a method with sum b c^k = 1/(k + 1), k <= 4, and
// y' = 6 t^5 is not (so the test can tell).
TEST_CASE("TC6: non-autonomous quadrature: y' = 5 t^4 exact, y' = 6 t^5 not") {
    auto f4 = [](const DD& t, const DDVector&, DDVector& dy) { dy[0] = DD(5.0) * t * t * t * t; };
    auto f5 = [](const DD& t, const DDVector&, DDVector& dy) { dy[0] = DD(6.0) * t * t * t * t * t; };
    const DD t_end(1.5);
    const DD exact4 = t_end * t_end * t_end * t_end * t_end, exact5 = exact4 * t_end;
    const double e4 = std::abs((ref::dopri5_fixed(f4, DD(0.0), {DD(0.0)}, t_end, 7)[0] - exact4).to_double());
    const double e5 = std::abs((ref::dopri5_fixed(f5, DD(0.0), {DD(0.0)}, t_end, 7)[0] - exact5).to_double());
    INFO("y' = 5 t^4: " << e4 / u_dd / exact4.hi << " u_dd relative; y' = 6 t^5: " << e5 << " absolute");
    CHECK(e4 <= 8.0 * u_dd * exact4.hi);
    CHECK(e5 > 1e-8);
}

// A first step far too large is rejected (dopri5.f does not count first-step rejections) and the run still meets
// the tolerance.
TEST_CASE("TC6: rejected steps: the controller recovers from h0 = 1") {
    const Problem p = problems()[2];
    ref::Dopri5Options o;
    o.rtol = o.atol = 1e-14;
    const ref::Dopri5Result base = ref::dopri5(p.f, DD(0.0), p.y0, p.t_end, o);
    o.h0 = 1.0;
    long long first_rejections = 0;
    const ref::Dopri5Result r = ref::dopri5(p.f, DD(0.0), p.y0, p.t_end, o);
    REQUIRE(r.status == ref::Dopri5Status::ok);
    first_rejections = r.nstep - r.naccpt - r.nrejct;
    const double e = max_error(r.y, p.exact(p.t_end)), e0 = max_error(base.y, p.exact(p.t_end));
    INFO("h0 = 1: " << r.nstep << " steps, " << r.naccpt << " accepted, " << r.nrejct << " counted rejections, "
                    << first_rejections << " first-step rejections; error " << e << " (default h0: " << e0 << ")");
    CHECK(first_rejections >= 1);
    CHECK(e <= 2.0 * e0);
}
