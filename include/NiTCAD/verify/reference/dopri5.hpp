// The Dormand-Prince 5(4) Runge-Kutta pair with Shampine's dense output, in double-double (Unit 55, TVR-1, TC6).
//
// Sources, read in full before code (2026-11-04): E. Hairer and G. Wanner, DOPRI5 (Fortran, version of 25 April
// 1996, corrected 8 August 2005), the code described in Hairer, Norsett and Wanner, Solving Ordinary Differential
// Equations I, 2nd ed. (Springer, 1993), section II.5, Table 5.2; read in the copy distributed with SciPy 1.10.0
// (scipy/integrate/dop/dopri5.f; the authors' page refused the download): CDOPRI (the coefficients, "Runge-Kutta
// coefficients of Dormand and Prince (1980)", and the "dense output of Shampine (1986)"), DOPCOR (the step, the
// error norm, the step-size control with Lund stabilization, acceptance and rejection), HINIT (the initial step)
// and CONTD5 (the dense output). Second reproduction: K. Ahnert and M. Mulansky, Boost.Odeint
// (boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp, develop branch): the same tableau, the 4th-order weights
// themselves (dopri5.f stores only their differences E_i). Boost's dense output (degree 5 in theta, cited there to
// HNW pp. 191-192) is a different order-4 extension from dopri5.f's (degree 4, "Shampine (1986)"): they agree at
// theta = 0, 1/2 and 1 only. Not read: Dormand and Prince (J. Comput. Appl. Math. 6, 1980) and Shampine (Math.
// Comp. 46, 1986), both refused to automated download, nor HNW's book. Provenance PB: the tableau from two agreeing
// reproductions (TVR1-Q2); the dense output from dopri5.f alone plus the tests' independent derivation (its eight
// order-4 conditions for all theta, b(1) = b and the Hermite end slopes; owner ruling 2026-11-04).
//
// The translation follows DOPCOR line by line, in double-double, with UROUND = u_dd. Deviations, each deliberate:
// no stiffness detection (it changes no accepted value, only stops an integration); no SOLOUT call at the initial
// point (there is no step to interpolate); the last step lands on tend exactly (DOPCOR sets X = X + H, which in
// double-double can differ from XEND in the last bit); the step size and time are double-double, the error norm
// and the step-size factors double (they choose h, they do not enter the solution). The tableau's rationals are
// typed from dopri5.f and divided in double-double (each entry within 1 u_dd of its rational).
#pragma once

#include "NiTCAD/verify/reference/double_double.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

namespace NiTCAD::verify::reference {

using DDVector = std::vector<DD>;

// Stage i = 0 ... 6; a[i][j] for j < i; row 6 is b (first same as last); e = b - bhat (dopri5.f's E1 ... E7);
// d are Shampine's dense coefficients D1 ... D7 (dopri5.f's d[6] multiplies the FSAL stage k7).
struct Dopri5Tableau {
    std::array<DD, 7> c{};
    std::array<std::array<DD, 7>, 7> a{};
    std::array<DD, 7> b{};
    std::array<DD, 7> e{};
    std::array<DD, 7> d{};
};

[[nodiscard]] inline const Dopri5Tableau& dopri5_tableau() {
    static const Dopri5Tableau t = [] {
        auto q = [](double n, double d) { return DD(n) / DD(d); };
        Dopri5Tableau r;
        r.c = {DD(0.0), q(1, 5), q(3, 10), q(4, 5), q(8, 9), DD(1.0), DD(1.0)};
        r.a[1] = {q(1, 5)};
        r.a[2] = {q(3, 40), q(9, 40)};
        r.a[3] = {q(44, 45), q(-56, 15), q(32, 9)};
        r.a[4] = {q(19372, 6561), q(-25360, 2187), q(64448, 6561), q(-212, 729)};
        r.a[5] = {q(9017, 3168), q(-355, 33), q(46732, 5247), q(49, 176), q(-5103, 18656)};
        r.a[6] = {q(35, 384), DD(0.0), q(500, 1113), q(125, 192), q(-2187, 6784), q(11, 84)};
        r.b = r.a[6];
        r.e = {q(71, 57600), DD(0.0), q(-71, 16695), q(71, 1920), q(-17253, 339200), q(22, 525), q(-1, 40)};
        r.d = {q(-12715105075.0, 11282082432.0), DD(0.0), q(87487479700.0, 32700410799.0),
               q(-10690763975.0, 1880347072.0), q(701980252875.0, 199316789632.0), q(-1453857185.0, 822651844.0),
               q(69997945.0, 29380423.0)};
        return r;
    }();
    return t;
}

// One step from (t, y) with k[0] = f(t, y) given: the stages, y1 (5th order), the error estimate h * sum e_i k_i
// and k[6] = f(t + h, y1).
struct Dopri5Step {
    std::array<DDVector, 7> k;
    DDVector y1;
    DDVector err;
};

// The dense output of one accepted step on [t_old, t_old + h]: CONTD5's five coefficient vectors.
struct Dopri5Dense {
    DD t_old;
    DD h;
    std::array<DDVector, 5> cont;
    [[nodiscard]] DDVector eval(const DD& t) const {
        const DD theta = (t - t_old) / h;
        const DD theta1 = DD(1.0) - theta;
        DDVector out(cont[0].size());
        for (std::size_t i = 0; i < out.size(); ++i)
            out[i] = cont[0][i] +
                     theta * (cont[1][i] + theta1 * (cont[2][i] + theta * (cont[3][i] + theta1 * cont[4][i])));
        return out;
    }
};

template <class F>
void dopri5_stages(F& f, const DD& t, const DDVector& y, const DD& h, Dopri5Step& s) {
    const Dopri5Tableau& tb = dopri5_tableau();
    const std::size_t n = y.size();
    DDVector yi(n);
    for (int i = 1; i < 7; ++i) {
        for (std::size_t m = 0; m < n; ++m) {
            DD acc(0.0);
            for (int j = 0; j < i; ++j) acc = acc + tb.a[i][j] * s.k[j][m];
            yi[m] = y[m] + h * acc;
        }
        s.k[i].assign(n, DD(0.0));
        f(t + tb.c[i] * h, yi, s.k[i]);
        if (i == 6) s.y1 = yi;
    }
    s.err.assign(n, DD(0.0));
    for (std::size_t m = 0; m < n; ++m) {
        DD acc(0.0);
        for (int j = 0; j < 7; ++j) acc = acc + tb.e[j] * s.k[j][m];
        s.err[m] = h * acc;
    }
}

// DOPCOR's CONT, from an accepted step.
[[nodiscard]] inline Dopri5Dense dopri5_dense(const DD& t, const DDVector& y, const DD& h, const Dopri5Step& s) {
    const Dopri5Tableau& tb = dopri5_tableau();
    Dopri5Dense dn{t, h, {}};
    const std::size_t n = y.size();
    for (auto& v : dn.cont) v.assign(n, DD(0.0));
    for (std::size_t m = 0; m < n; ++m) {
        const DD ydiff = s.y1[m] - y[m];
        const DD bspl = h * s.k[0][m] - ydiff;
        DD acc(0.0);
        for (int j = 0; j < 7; ++j) acc = acc + tb.d[j] * s.k[j][m];
        dn.cont[0][m] = y[m];
        dn.cont[1][m] = ydiff;
        dn.cont[2][m] = bspl;
        dn.cont[3][m] = -(h * s.k[6][m]) + ydiff - bspl;
        dn.cont[4][m] = h * acc;
    }
    return dn;
}

struct NoDenseOutput {
    void operator()(const Dopri5Dense&) const {}
};

// n equal steps from t0 to tend (t_k = t0 + k h, h = (tend - t0) / n): the fixed-step runs of the order tests.
template <class F, class Out = NoDenseOutput>
[[nodiscard]] DDVector dopri5_fixed(F&& f, const DD& t0, DDVector y, const DD& tend, long long n, Out&& out = {}) {
    const DD h = (tend - t0) / DD(static_cast<double>(n));
    Dopri5Step s;
    s.k[0].assign(y.size(), DD(0.0));
    f(t0, y, s.k[0]);
    for (long long step = 0; step < n; ++step) {
        const DD t = t0 + DD(static_cast<double>(step)) * h;
        dopri5_stages(f, t, y, h, s);
        out(dopri5_dense(t, y, h, s));
        y = s.y1;
        s.k[0] = s.k[6];
    }
    return y;
}

struct Dopri5Options {
    double rtol = 1e-10;
    double atol = 1e-10;
    double safe = 0.9;     // WORK(2)
    double fac1 = 0.2;     // WORK(3): hnew / hold >= fac1
    double fac2 = 10.0;    // WORK(4): hnew / hold <= fac2
    double beta = 0.04;    // WORK(5), Lund stabilization
    double hmax = 0.0;     // WORK(6); 0: |tend - t0|
    double h0 = 0.0;       // WORK(7); 0: HINIT
    long long nmax = 100000;
    double uround = u_dd;  // WORK(1)
};

enum class Dopri5Status {
    ok,              // IDID = 1
    invalid_input,   // IDID = -1
    too_many_steps,  // IDID = -2
    step_too_small,  // IDID = -3
};

struct Dopri5Result {
    DDVector y;
    DD t;
    DD h;  // the predicted step after the last accepted one
    Dopri5Status status = Dopri5Status::ok;
    long long nfcn = 0, nstep = 0, naccpt = 0, nrejct = 0;
};

namespace detail {
[[nodiscard]] inline double dopri5_scale(const Dopri5Options& o, const DD& a, const DD& b) {
    return o.atol + o.rtol * std::max(std::abs(a.to_double()), std::abs(b.to_double()));
}

// HINIT, IORD = 5.
template <class F>
[[nodiscard]] double dopri5_hinit(F& f, const DD& x, const DDVector& y, double posneg, const DDVector& f0,
                                  double hmax, const Dopri5Options& o) {
    const std::size_t n = y.size();
    double dnf = 0.0, dny = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double sk = o.atol + o.rtol * std::abs(y[i].to_double());
        dnf += std::pow(f0[i].to_double() / sk, 2);
        dny += std::pow(y[i].to_double() / sk, 2);
    }
    double h = (dnf <= 1e-10 || dny <= 1e-10) ? 1e-6 : std::sqrt(dny / dnf) * 0.01;
    h = std::copysign(std::min(h, hmax), posneg);
    DDVector y1(n), f1(n, DD(0.0));
    for (std::size_t i = 0; i < n; ++i) y1[i] = y[i] + DD(h) * f0[i];
    f(x + DD(h), y1, f1);
    double der2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        const double sk = o.atol + o.rtol * std::abs(y[i].to_double());
        der2 += std::pow((f1[i] - f0[i]).to_double() / sk, 2);
    }
    der2 = std::sqrt(der2) / h;
    const double der12 = std::max(std::abs(der2), std::sqrt(dnf));
    const double h1 = der12 <= 1e-15 ? std::max(1e-6, std::abs(h) * 1e-3) : std::pow(0.01 / der12, 1.0 / 5.0);
    return std::copysign(std::min({100.0 * std::abs(h), h1, hmax}), posneg);
}
}  // namespace detail

// DOPCOR: the adaptive integration from (t0, y) to tend; out(dense) after every accepted step.
template <class F, class Out = NoDenseOutput>
[[nodiscard]] Dopri5Result dopri5(F&& f, DD t0, DDVector y, const DD& tend, const Dopri5Options& o,
                                  Out&& out = {}) {
    Dopri5Result r;
    const std::size_t n = y.size();
    if (n == 0 || !(o.rtol >= 0.0) || !(o.atol >= 0.0) || (o.rtol == 0.0 && o.atol == 0.0) ||
        !(o.uround > 1e-35 && o.uround < 1.0) || !(o.safe > 1e-4 && o.safe < 1.0) || !(o.beta <= 0.2) ||
        o.nmax <= 0 || !std::isfinite(t0.hi) || !std::isfinite(tend.hi) || tend == t0) {
        r.status = Dopri5Status::invalid_input;
        r.y = std::move(y);
        r.t = t0;
        return r;
    }
    const double beta = o.beta < 0.0 ? 0.0 : o.beta;
    double facold = 1e-4;
    const double expo1 = 0.2 - beta * 0.75;
    const double facc1 = 1.0 / o.fac1, facc2 = 1.0 / o.fac2;
    const double posneg = (tend - t0).hi > 0.0 ? 1.0 : -1.0;
    const double hmax = o.hmax == 0.0 ? std::abs((tend - t0).to_double()) : std::abs(o.hmax);
    bool last = false, reject = false;
    DD x = t0;
    Dopri5Step s;
    s.k[0].assign(n, DD(0.0));
    f(x, y, s.k[0]);
    DD h(o.h0);
    if (o.h0 == 0.0) h = DD(detail::dopri5_hinit(f, x, y, posneg, s.k[0], hmax, o));
    r.nfcn += 2;
    for (;;) {
        if (r.nstep > o.nmax) {
            r.status = Dopri5Status::too_many_steps;
            break;
        }
        if (0.1 * std::abs(h.to_double()) <= std::abs(x.to_double()) * o.uround) {
            r.status = Dopri5Status::step_too_small;
            break;
        }
        if (((x + DD(1.01) * h - tend).hi) * posneg > 0.0) {
            h = tend - x;
            last = true;
        }
        ++r.nstep;
        dopri5_stages(f, x, y, h, s);
        r.nfcn += 6;
        double err = 0.0;
        for (std::size_t i = 0; i < n; ++i) err += std::pow(s.err[i].to_double() / detail::dopri5_scale(o, y[i], s.y1[i]), 2);
        err = std::sqrt(err / static_cast<double>(n));
        const double fac11 = std::pow(err, expo1);
        double fac = fac11 / std::pow(facold, beta);  // Lund stabilization
        fac = std::max(facc2, std::min(facc1, fac / o.safe));
        DD hnew = h / DD(fac);
        if (err <= 1.0) {
            facold = std::max(err, 1e-4);
            ++r.naccpt;
            out(dopri5_dense(x, y, h, s));
            y = s.y1;
            s.k[0] = s.k[6];
            x = last ? tend : x + h;
            if (last) {
                r.h = hnew;
                r.status = Dopri5Status::ok;
                break;
            }
            if (std::abs(hnew.to_double()) > hmax) hnew = DD(posneg * hmax);
            if (reject) hnew = std::abs(hnew.to_double()) < std::abs(h.to_double()) ? hnew : h;
            reject = false;
        } else {
            hnew = h / DD(std::min(facc1, fac11 / o.safe));
            reject = true;
            if (r.naccpt >= 1) ++r.nrejct;
            last = false;
        }
        h = hnew;
    }
    r.y = std::move(y);
    r.t = x;
    return r;
}

}  // namespace NiTCAD::verify::reference
