#include "NiTCAD/verify/reference/qe_junction.hpp"

#include <cmath>
#include <complex>

#include "NiTCAD/verify/reference/brent.hpp"
#include "NiTCAD/verify/reference/complex_step.hpp"
#include "NiTCAD/verify/reference/first_integral.hpp"

namespace NiTCAD::verify::reference {

namespace {

bool positive(double v) { return std::isfinite(v) && v > 0.0; }

// A dual number over double-double: value and derivative with respect to v.
struct Dual {
    DD v, d;
};
Dual operator+(const Dual& a, const Dual& b) { return {a.v + b.v, a.d + b.d}; }
Dual operator-(const Dual& a, const Dual& b) { return {a.v - b.v, a.d - b.d}; }
Dual operator*(const Dual& a, const Dual& b) { return {a.v * b.v, a.d * b.v + a.v * b.d}; }
Dual operator/(const Dual& a, const Dual& b) { return {a.v / b.v, (a.d * b.v - a.v * b.d) / (b.v * b.v)}; }
Dual constant(double c) { return {DD(c), DD(0.0)}; }
Dual sqrt(const Dual& a) {
    const DD s = reference::sqrt(a.v);
    return {s, a.d / ldexp(s, 1)};
}
Dual log(const Dual& a) { return {reference::log(a.v), a.d / a.v}; }
Dual G(const Dual& t) { return {poisson_boltzmann_G(t.v), expm1(t.v) * t.d}; }  // G' = e^t - 1
Dual Gp(const Dual& t) { return {expm1(t.v), exp(t.v) * t.d}; }                  // G'
DD G(const DD& t) { return poisson_boltzmann_G(t); }

using reference::complex_d;
complex_d G(const complex_d& t) { return std::exp(t) - 1.0 - t; }
complex_d Gp(const complex_d& t) { return std::exp(t) - 1.0; }

// The bulks at m, in any scalar with +, -, *, /, sqrt, log.
template <class S>
struct Bulks {
    S nn, pn, un, np, pp, up;
};

template <class S, class C>
Bulks<S> bulks(const S& m, C NA, C ND) {
    Bulks<S> b;
    const S half_d = ND * 0.5, half_a = NA * 0.5;
    using std::sqrt, std::log;
    b.nn = half_d + sqrt(half_d * half_d + m * m);
    b.pn = m * m / b.nn;
    b.un = log(b.nn / m);
    b.pp = half_a + sqrt(half_a * half_a + m * m);
    b.np = m * m / b.pp;
    b.up = log(b.np / m);
    return b;
}

// Phi_n(u), Phi_n - Phi_p and its u-derivative.
template <class S>
S phi_n(const Bulks<S>& b, const S& u) { return b.nn * G(u - b.un) + b.pn * G(b.un - u); }
template <class S>
S mismatch(const Bulks<S>& b, const S& u) {
    return phi_n(b, u) - (b.np * G(u - b.up) + b.pp * G(b.up - u));
}
template <class S>
S mismatch_u(const Bulks<S>& b, const S& u) {
    return b.nn * Gp(u - b.un) - b.pn * Gp(b.un - u) - (b.np * Gp(u - b.up) - b.pp * Gp(b.up - u));
}

// Bulks<Dual> from double inputs: the doping is a constant.
Bulks<Dual> dual_bulks(const Dual& m, double NA, double ND) {
    Bulks<Dual> b;
    const Dual half_d = constant(0.5 * ND), half_a = constant(0.5 * NA);
    b.nn = half_d + sqrt(half_d * half_d + m * m);
    b.pn = m * m / b.nn;
    b.un = log(b.nn / m);
    b.pp = half_a + sqrt(half_a * half_a + m * m);
    b.np = m * m / b.pp;
    b.up = log(b.np / m);
    return b;
}

std::expected<void, std::string> check(const PnJunctionInput& in, double bias_V) {
    if (!positive(in.q_C) || !positive(in.eps_F_per_cm) || !positive(in.thermal_voltage_V) || !positive(in.ni_cm3) ||
        !positive(in.acceptors_cm3) || !positive(in.donors_cm3)) {
        return std::unexpected("qe_junction: inputs must be finite and positive");
    }
    if (!std::isfinite(bias_V)) return std::unexpected("qe_junction: bias must be finite");
    return {};
}

}  // namespace

std::expected<QeJunction, std::string> qe_junction(const PnJunctionInput& in, double bias_V) {
    if (auto ok = check(in, bias_V); !ok) return std::unexpected(ok.error());
    const DD VT(in.thermal_voltage_V), v = DD(bias_V) / VT;
    QeJunction j;
    j.m_cm3 = DD(in.ni_cm3) * exp(ldexp(v, -1));
    // m^2 (in p_n and n_p) and its double-double tail must stay normal: about 1e-140 cm^-3 or more.
    if (!(j.m_cm3.hi >= 1e-140) || !std::isfinite(j.m_cm3.hi)) {
        return std::unexpected("qe_junction: bias too far from 0 (n_i e^(V / 2 V_T) outside [1e-140, inf))");
    }
    const Bulks<DD> b = bulks(j.m_cm3, DD(in.acceptors_cm3), DD(in.donors_cm3));
    j.u_n = b.un;
    j.u_p = b.up;
    j.barrier_V = VT * (b.un - b.up);
    const double t = u_dd * (std::abs(b.un.hi) + std::abs(b.up.hi));
    const BrentResult r = brent_zero([&](const DD& u) { return mismatch(b, u); }, b.up, b.un, t);
    if (r.status != BrentStatus::ok) return std::unexpected("qe_junction: Brent did not converge for u_j");
    j.u_j = r.root;
    j.brent_evaluations = r.evaluations;
    const DD scale = DD(2.0) * DD(in.q_C) * DD(in.eps_F_per_cm) * VT;
    // One dual Newton step from the root: its dual part is du_j/dv = -F_v / F_u.
    const Dual m{j.m_cm3, ldexp(j.m_cm3, -1)};  // dm/dv = m / 2
    const Bulks<Dual> bd = dual_bulks(m, in.acceptors_cm3, in.donors_cm3);
    const Dual u0{j.u_j, DD(0.0)};
    const Dual u1 = u0 - mismatch(bd, u0) / mismatch_u(bd, u0);
    const Dual Q = sqrt(Dual{scale, DD(0.0)} * phi_n(bd, u1));
    j.charge_C_per_cm2 = reference::sqrt(scale * phi_n(b, j.u_j));
    j.capacitance_F_per_cm2 = -Q.d / VT;
    const DD eps(in.eps_F_per_cm), q(in.q_C);
    const DD W = reference::sqrt(ldexp(eps, 1) * (j.barrier_V - ldexp(VT, 1)) / q *
                                 (DD(1.0) / DD(in.acceptors_cm3) + DD(1.0) / DD(in.donors_cm3)));
    j.depletion_form_F_per_cm2 = eps / W;
    return j;
}

std::expected<double, std::string> qe_capacitance_complex_step(const PnJunctionInput& in, double bias_V) {
    const auto real = qe_junction(in, bias_V);
    if (!real) return std::unexpected(real.error());
    const double u_start = real->u_j.to_double(), VT = in.thermal_voltage_V;
    const double scale = 2.0 * in.q_C * in.eps_F_per_cm * VT;
    bool converged = true;
    const auto Q = [&](const complex_d& v) {
        const complex_d m = in.ni_cm3 * std::exp(0.5 * v);
        const Bulks<complex_d> b = bulks(m, in.acceptors_cm3, in.donors_cm3);
        complex_d u(u_start, 0.0);
        bool done = false;
        for (int k = 0; k < 20 && !done; ++k) {
            const complex_d du = mismatch(b, u) / mismatch_u(b, u);
            u -= du;
            done = std::abs(du.real()) <= 1e-15 * std::abs(u.real()) && std::abs(du.imag()) <= 1e-15 * std::abs(u.imag());
        }
        converged = converged && done;
        return std::sqrt(scale * phi_n(b, u));
    };
    const ComplexStepResult cs = complex_step(Q, bias_V / VT);
    if (cs.status != ComplexStepStatus::ok) return std::unexpected("qe_capacitance_complex_step: complex step failed");
    if (!converged) return std::unexpected("qe_capacitance_complex_step: Newton did not converge");
    return -cs.derivative / VT;
}

}  // namespace NiTCAD::verify::reference
