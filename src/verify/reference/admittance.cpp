#include "NiTCAD/verify/reference/admittance.hpp"

#include <algorithm>
#include <cmath>

#include "NiTCAD/verify/reference/dopri5.hpp"
#include "NiTCAD/verify/reference/minority_diffusion.hpp"
#include "NiTCAD/verify/reference/pn_junction.hpp"
#include "NiTCAD/verify/reference/qe_junction.hpp"
#include "NiTCAD/verify/reference/special_functions.hpp"

namespace NiTCAD::verify::reference {

namespace {

// Complex double-double, the few operations the closed form needs.
struct CD {
    DD re, im;
};
CD operator+(const CD& a, const CD& b) { return {a.re + b.re, a.im + b.im}; }
CD operator*(const CD& a, const CD& b) { return {a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re}; }
CD operator*(const DD& s, const CD& a) { return {s * a.re, s * a.im}; }
CD operator/(const CD& a, const CD& b) {
    const DD d = b.re * b.re + b.im * b.im;
    return {(a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d};
}
DD norm2(const CD& a) { return a.re * a.re + a.im * a.im; }
std::complex<double> to_complex(const CD& a) { return {a.re.to_double(), a.im.to_double()}; }

// cos b and sin b for any finite b: reduced by 2 pi to [-pi, pi], then halved (special_functions' series hold for
// |x| <= 2): sin r = 2 sin(r/2) cos(r/2), cos r = 1 - 2 sin^2(r/2). Also cos b - 1 = -2 sin^2(b/2) without
// cancellation (sin^2 of the half of the reduced angle).
struct Trig {
    DD c, s, cm1;
};
Trig trig(const DD& b) {
    const DD two_pi = ldexp(pi(), 1);
    const double n = std::nearbyint((b / two_pi).to_double());
    const DD r = b - DD(n) * two_pi;
    const DD h = ldexp(r, -1);
    const DD sh = sin(h), ch = cos(h);
    const DD cm1 = -ldexp(sh * sh, 1);
    return {DD(1.0) + cm1, ldexp(sh * ch, 1), cm1};
}

// e^z - 1 without cancellation: (expm1(a) cos b + (cos b - 1)) + i e^a sin b.
CD expm1(const CD& z) {
    const Trig t = trig(z.im);
    return {reference::expm1(z.re) * t.c + t.cm1, exp(z.re) * t.s};
}

// coth z for Re z > 0: (1 + e^-2z) / (1 - e^-2z) = (2 + expm1(-2z)) / -expm1(-2z); 1 when e^-2Re z < 1e-87.
CD coth(const CD& z) {
    if (z.re.to_double() > 100.0) return {DD(1.0), DD(0.0)};
    const CD m = expm1(CD{-ldexp(z.re, 1), -ldexp(z.im, 1)});
    return CD{DD(2.0) + m.re, m.im} / CD{-m.re, -m.im};
}

// The principal square root for Re z >= 0 (kappa^2 = (1/tau + i w) / D always is), without cancellation:
// x = sqrt((|z| + a) / 2), y = b / 2x.
CD sqrt(const CD& z) {
    const DD r = reference::sqrt(norm2(z));
    if (r.hi == 0.0) return {DD(0.0), DD(0.0)};
    const DD x = reference::sqrt(ldexp(r + z.re, -1));
    return {x, z.im / ldexp(x, 1)};
}

struct SideData {
    const DiodeSide* s = nullptr;
    DD minority;       // n_0
    DD inv_tau;        // SRH's linearization (0 without srh)
    DD base;           // W'
    DD rate;           // dW'/dV
};

// kappa coth(kappa W') by the Riccati form (header), integrated in t = y / l.
std::expected<CD, std::string> riccati(const CD& kappa, const DD& base) {
    const DD mod = reference::sqrt(norm2(kappa));
    const bool zero = mod.hi == 0.0;
    const DD l = zero ? base : DD(1.0) / mod;
    const CD unit = zero ? CD{DD(0.0), DD(0.0)} : CD{kappa.re * l, kappa.im * l};  // kappa l, |kappa l| = 1
    const CD a = unit * unit;
    const DD T = base / l;
    constexpr double reach = 20.0;  // the contact's influence at the edge is e^-2 reach
    const DD t0 = T.to_double() > reach ? T - DD(reach) : DD(0.0);
    DDVector y0{DD(0.0), DD(0.0)};
    if (t0.hi > 0.0) {  // w = 1 / kappa: z = 1 / (kappa l) = conj(unit)
        y0 = {unit.re, -unit.im};
    }
    const auto f = [&](const DD&, const DDVector& y, DDVector& dy) {
        const CD z{y[0], y[1]};
        const CD az = a * (z * z);
        dy[0] = DD(1.0) - az.re;
        dy[1] = -az.im;
    };
    Dopri5Options o;
    o.rtol = 1e-15;
    o.atol = 1e-15;
    o.nmax = 2000000;
    const Dopri5Result r = dopri5(f, t0, y0, T, o);
    if (r.status != Dopri5Status::ok) return std::unexpected("diode_admittance: the Riccati integration failed");
    const CD z{r.y[0], r.y[1]};
    return CD{DD(1.0), DD(0.0)} / CD{z.re * l, z.im * l};
}

}  // namespace

std::expected<DiodeAdmittance, std::string> diode_admittance(const DiodeInput& in, double bias_V,
                                                            std::span<const double> omegas_rad_s) {
    for (const double w : omegas_rad_s) {
        if (!std::isfinite(w) || w < 0.0) return std::unexpected("diode_admittance: a frequency not finite or negative");
    }
    const auto dc = diode_current(in, bias_V);
    if (!dc) return std::unexpected(dc.error());
    const PnJunctionInput& j = in.junction;
    const auto qe = qe_junction(j, bias_V);
    if (!qe) return std::unexpected("diode_admittance: " + qe.error());
    const auto eq = pn_junction(j);
    if (!eq) return std::unexpected("diode_admittance: equilibrium: " + eq.error());

    DiodeAdmittance out;
    out.dc = *dc;
    out.capacitance_F_per_cm2 = qe->capacitance_F_per_cm2;
    const DD q(j.q_C), VT(j.thermal_voltage_V), v = DD(bias_V) / VT;
    const DD ev = exp(v), em1 = reference::expm1(v);
    SideData p{&in.p_side, eq->n_p, DD(0.0), dc->base_p_cm, out.capacitance_F_per_cm2 / (q * DD(j.acceptors_cm3))};
    SideData n{&in.n_side, eq->p_n, DD(0.0), dc->base_n_cm, out.capacitance_F_per_cm2 / (q * DD(j.donors_cm3))};
    if (in.srh) {
        p.inv_tau = srh_linear_inverse_lifetime(eq->n_p.to_double(), eq->p_p.to_double(), j.ni_cm3, in.p_side.tau_n_s,
                                                in.p_side.tau_p_s);
        n.inv_tau = srh_linear_inverse_lifetime(eq->n_n.to_double(), eq->p_n.to_double(), j.ni_cm3, in.n_side.tau_n_s,
                                                in.n_side.tau_p_s);
        out.lifetime_p_s = (DD(1.0) / p.inv_tau).to_double();
        out.lifetime_n_s = (DD(1.0) / n.inv_tau).to_double();
    }
    out.base_rate_p_cm_per_V = p.rate.to_double();
    out.base_rate_n_cm_per_V = n.rate.to_double();

    for (const double w : omegas_rad_s) {
        AdmittancePoint pt;
        pt.omega_rad_s = w;
        CD sum_fixed{DD(0.0), DD(0.0)}, sum_riccati{DD(0.0), DD(0.0)}, sum_edge{DD(0.0), DD(0.0)};
        for (SideData* sd : {&p, &n}) {
            const DD D(sd->s->minority_diffusivity_cm2_s);
            const CD kappa = sqrt(CD{sd->inv_tau / D, DD(w) / D});
            const bool zero = kappa.re.hi == 0.0 && kappa.im.hi == 0.0;
            // kappa coth(kappa W'), and 1 / W' in the limit kappa -> 0.
            const CD kc = zero ? CD{DD(1.0) / sd->base, DD(0.0)} : kappa * coth(kappa * CD{sd->base, DD(0.0)});
            const auto kr = riccati(kappa, sd->base);
            if (!kr) return std::unexpected(kr.error());
            const DD G = q * D * sd->minority * ev / VT;  // per kappa coth
            const CD fixed = G * kc, fixed2 = G * *kr;
            // u0'/U0 and u0''/U0 at the edge.
            DD g, h;
            if (in.srh) {
                const DD L = reference::sqrt(D / sd->inv_tau);
                const CD c = coth(CD{sd->base / L, DD(0.0)});
                g = c.re / L;
                h = DD(1.0) / (L * L);
            } else {
                g = DD(1.0) / sd->base;
                h = DD(0.0);
            }
            const DD U0 = sd->minority * em1;
            const DD scale = q * D * U0 * sd->rate;
            const CD edge = scale * CD{h - kc.re * g, -kc.im * g};
            AdmittanceSide& side = sd == &p ? pt.p : pt.n;
            side.fixed_edge = to_complex(fixed);
            side.fixed_edge_riccati = to_complex(fixed2);
            side.edge_motion = to_complex(edge);
            side.effective_length_cm =
                zero ? std::complex<double>(INFINITY, 0.0) : to_complex(CD{DD(1.0), DD(0.0)} / kappa);
            sum_fixed = sum_fixed + fixed;
            sum_riccati = sum_riccati + fixed2;
            sum_edge = sum_edge + edge;
        }
        const CD cap{DD(0.0), DD(w) * out.capacitance_F_per_cm2};
        pt.capacitive = to_complex(cap);
        pt.shockley = to_complex(sum_fixed + cap);
        pt.total = to_complex(sum_fixed + cap + sum_edge);
        pt.second_method_difference = std::abs(to_complex(sum_fixed + DD(-1.0) * sum_riccati));
        out.points.push_back(pt);
    }
    return out;
}

}  // namespace NiTCAD::verify::reference
