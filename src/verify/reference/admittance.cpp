#include "NiTCAD/verify/reference/admittance.hpp"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include "NiTCAD/verify/reference/dopri5.hpp"
#include "NiTCAD/verify/reference/first_integral.hpp"
#include "NiTCAD/verify/reference/gauss_kronrod.hpp"
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

// Gauss-Kronrod in u for C_SCR and the second method's charges: far below any criterion (C_SCR is a few percent of C_j)
// and of the second method's rounding (1e-16 / h), and cheap enough for the Debug budget (CV-Q19).
constexpr double scr_epsrel = 1e-20;

// The QE junction at the reduced bias v (frame of diode.hpp: n = m e^u, p = m e^-u).
struct ScrFrame {
    DD m, u_j;
    NeutralBulk bn, bp;
    PnJunction jm;
};
std::expected<ScrFrame, std::string> scr_frame(const PnJunctionInput& j, const DD& v) {
    PnJunctionInput shifted = j;
    shifted.ni_cm3 = (DD(j.ni_cm3) * exp(ldexp(v, -1))).to_double();
    const auto jm = pn_junction(shifted);
    if (!jm) return std::unexpected("diode_admittance: C_SCR at bias: " + jm.error());
    return ScrFrame{DD(shifted.ni_cm3), jm->u_j, neutral_bulk(shifted.ni_cm3, j.donors_cm3),
                    neutral_bulk(shifted.ni_cm3, -j.acceptors_cm3), *jm};
}

// The levels at the distances -x_p and x_n from the junction (pn_potential, frame of f).
std::expected<std::pair<DD, DD>, std::string> scr_levels(const PnJunctionInput& j, const ScrFrame& f, const DD& xp,
                                                         const DD& xn) {
    const std::vector<double> at{-xp.to_double(), xn.to_double()};
    const auto u = pn_potential(f.jm, at);
    if (!u) return std::unexpected("diode_admittance: C_SCR edge potential: " + u.error());
    const DD VT(j.thermal_voltage_V);
    return std::pair<DD, DD>{(*u)[0].psi_V / VT, (*u)[1].psi_V / VT};
}

// q times the minority carriers between the levels u_ep (p side) and u_en (n side) of frame f [C/cm^2].
std::expected<DD, std::string> scr_charge(const PnJunctionInput& j, const ScrFrame& f, const DD& u_ep,
                                          const DD& u_en) {
    const DD q(j.q_C), VT(j.thermal_voltage_V);
    const DD k = ldexp(q / (DD(j.eps_F_per_cm) * VT), 1);
    const auto on_n = [&](const DD& u) { return f.m * exp(-u) / sqrt(k * first_integral_phi(f.bn, u)); };
    const auto on_p = [&](const DD& u) { return f.m * exp(u) / sqrt(k * first_integral_phi(f.bp, u)); };
    const AdaptiveResult rn = qag15(on_n, f.u_j, u_en, 0.0, scr_epsrel);
    const AdaptiveResult rp = qag15(on_p, u_ep, f.u_j, 0.0, scr_epsrel);
    if (rn.status != QuadratureStatus::ok || rp.status != QuadratureStatus::ok) {
        return std::unexpected("diode_admittance: C_SCR second-method quadrature did not converge");
    }
    return q * (rn.value + rp.value);
}

// C_SCR = q (int_0^x_n dp/dV|_x dx + int_-x_p^0 dn/dV|_x dx) (DECISIONS Unit 55 16.38, DT-Q34 A; derived in
// admittance.hpp).
std::expected<DD, std::string> scr_storage(const PnJunctionInput& j, double bias_V) {
    const DD q(j.q_C), VT(j.thermal_voltage_V), v = DD(bias_V) / VT;
    const auto f = scr_frame(j, v);
    if (!f) return std::unexpected(f.error());
    const DD xp = f->jm.depletion_p_cm, xn = f->jm.depletion_n_cm;
    const auto edges = scr_levels(j, *f, xp, xn);
    if (!edges) return std::unexpected(edges.error());
    const auto [u_ep, u_en] = *edges;
    const DD m = f->m, uj = f->u_j;
    const DD k = ldexp(q / (DD(j.eps_F_per_cm) * VT), 1);  // u'^2 = k Phi
    const DD NA(j.acceptors_cm3), ND(j.donors_cm3), four_m2 = ldexp(m * m, 2);
    const DD duj = (sqrt(NA * NA + four_m2) - sqrt(ND * ND + four_m2)) / ldexp(NA + ND, 1);  // du_j/dv
    // b(s) = d/dv (k Phi(s))^(-1/2) at fixed s, with dPhi/dv = (n + p - n_b - p_b) / 2.
    const auto b = [&](const NeutralBulk& bulk, const DD& s) {
        const DD kphi = k * first_integral_phi(bulk, s);
        const DD dphi = ldexp(m * (exp(s) + exp(-s)) - (bulk.n_b + bulk.p_b), -1);
        return -ldexp(k * dphi, -1) / (kphi * sqrt(kphi));
    };
    const DD sj = sqrt(k * first_integral_phi(f->bn, uj));
    // One integrand per side: (1/2) carrier / sqrt(k Phi) + b(s) times the carriers beyond s (closed form).
    // n side: holes p = m e^-u from u_j to u_en; p side: electrons n = m e^u from u_ep to u_j.
    const auto on_n = [&](const DD& u) {
        const DD p = m * exp(-u);
        return ldexp(p, -1) / sqrt(k * first_integral_phi(f->bn, u)) + b(f->bn, u) * (p - m * exp(-u_en));
    };
    const auto on_p = [&](const DD& u) {
        const DD n = m * exp(u);
        return ldexp(n, -1) / sqrt(k * first_integral_phi(f->bp, u)) + b(f->bp, u) * (n - m * exp(u_ep));
    };
    const AdaptiveResult rn = qag15(on_n, uj, u_en, 0.0, scr_epsrel);
    const AdaptiveResult rp = qag15(on_p, u_ep, uj, 0.0, scr_epsrel);
    if (rn.status != QuadratureStatus::ok || rp.status != QuadratureStatus::ok) {
        return std::unexpected("diode_admittance: C_SCR quadrature did not converge");
    }
    const DD Pe = m * (exp(-uj) - exp(-u_en));  // int p du over the n side
    const DD Ne = m * (exp(uj) - exp(u_ep));    // int n du over the p side
    const DD per_v = rn.value + rp.value + (-duj / sj) * Pe + (duj / sj) * Ne;
    return q * per_v / VT;
}

}  // namespace

// Second method (admittance.hpp): the charge between the fixed positions -x_p, x_n (at V) as a function of the bias,
// differenced by Richardson's 4-point rule in v (h = 1e-3: truncation ~h^4, rounding of m to double ~1e-16 / h).
std::expected<double, std::string> scr_storage_second_method(const PnJunctionInput& j, double bias_V) {
    const DD VT(j.thermal_voltage_V), v = DD(bias_V) / VT;
    const auto f = scr_frame(j, v);
    if (!f) return std::unexpected(f.error());
    const DD xp = f->jm.depletion_p_cm, xn = f->jm.depletion_n_cm;
    const DD h(1e-3);
    DD g[4];
    const double steps[4] = {-2.0, -1.0, 1.0, 2.0};
    for (int i = 0; i < 4; ++i) {
        const auto fi = scr_frame(j, v + DD(steps[i]) * h);
        if (!fi) return std::unexpected(fi.error());
        const auto ui = scr_levels(j, *fi, xp, xn);
        if (!ui) return std::unexpected(ui.error());
        const auto gi = scr_charge(j, *fi, ui->first, ui->second);
        if (!gi) return std::unexpected(gi.error());
        g[i] = *gi;
    }
    return ((DD(8.0) * (g[2] - g[1]) - (g[3] - g[0])) / (DD(12.0) * h) / VT).to_double();
}

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
    const auto scr = scr_storage(j, bias_V);
    if (!scr) return std::unexpected(scr.error());
    out.scr_storage_F_per_cm2 = *scr;
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
        const CD scr_cap{DD(0.0), DD(w) * out.scr_storage_F_per_cm2};
        pt.scr_storage = to_complex(scr_cap);
        pt.total = to_complex(sum_fixed + cap + sum_edge + scr_cap);
        pt.second_method_difference = std::abs(to_complex(sum_fixed + DD(-1.0) * sum_riccati));
        out.points.push_back(pt);
    }
    return out;
}

}  // namespace NiTCAD::verify::reference
