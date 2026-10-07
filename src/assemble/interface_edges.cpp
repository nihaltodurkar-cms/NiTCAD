#include "NiTCAD/assemble/interface_edges.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

#include "NiTCAD/physics/interface_traps.hpp"
#include "NiTCAD/physics/recombination.hpp"
#include "NiTCAD/physics/statistics.hpp"
#include "scaled_device.hpp"

namespace NiTCAD::assemble {

namespace {

// The root of F(psi) = G psi - b - Q(psi), Q non-increasing and within [q_low, q_high]: Newton
// safeguarded by the bracket [(b + q_low) / G, (b + q_high) / G], to rounding. `eval(psi)` returns
// (Q, dQ/dpsi).
template <class Eval>
double interface_root(double G, double b, double q_low, double q_high, Eval eval) {
    double lo = (b + q_low) / G, hi = (b + q_high) / G;
    if (!(hi > lo)) return lo;  // Q constant (a fixed charge alone): exact
    double psi = std::clamp((b + eval(b / G).first) / G, lo, hi);
    for (int it = 0; it < 100; ++it) {
        const auto [Q, dQ] = eval(psi);
        const double F = G * psi - b - Q;
        if (F == 0.0) break;
        (F > 0.0 ? hi : lo) = psi;
        double next = psi - F / (G - dQ);
        if (!(next > lo && next < hi)) next = 0.5 * (lo + hi);
        const bool done = std::abs(next - psi) <= 2.0 * std::numeric_limits<double>::epsilon() *
                                                       std::max(1.0, std::abs(psi));
        psi = next;
        if (done || !(hi > lo)) break;
    }
    return psi;
}

}  // namespace

InterfaceEdges::InterfaceEdges(std::vector<InterfaceEdge> edges,
                               std::vector<InterfaceLevel> levels, std::size_t interfaces)
    : edges_(std::move(edges)), levels_(std::move(levels)), interfaces_(interfaces) {}

InterfaceEquilibrium InterfaceEdges::equilibrium(const InterfaceEdge& e, double psi_i,
                                                 double psi_s,
                                                 const InterfaceStatistics& s) const {
    // Q(psi_I) by the Fermi function of tau - (psi_I + s); both kinds fall with eta at N f (1 - f).
    const auto Q = [&](double psi) {
        double q = e.fixed_charge_cm2, dq = 0.0;
        for (std::size_t k = e.first_level; k < e.last_level; ++k) {
            const InterfaceLevel& l = levels_[k];
            const physics::FermiOccupancy f =
                physics::fermi_occupancy(l.tau - (psi + s.band_shift));
            q += l.density_cm2 * (l.donor ? f.empty : -f.occupied);
            dq -= l.density_cm2 * f.d_eta;
        }
        return std::pair{e.charge_weight * q, e.charge_weight * dq};
    };
    const double gi = e.g_insulator, gs = e.g_semiconductor, G = gi + gs;
    const double psi = interface_root(G, gi * psi_i + gs * psi_s, e.charge_low, e.charge_high, Q);
    const auto [q, dq] = Q(psi);
    // d psi_I / d psi_i = g_i / D, d psi_I / d psi_s = g_s / D, D = G - dQ/dpsi_I.
    const double D = G - dq;
    const double a = gi / D, b = gs / D;
    InterfaceEquilibrium r{};
    r.psi = psi;
    r.flux_insulator = gi * (psi - psi_i);
    r.flux_semiconductor = gs * (psi - psi_s);
    r.d_flux_insulator[0] = gi * (a - 1.0);
    r.d_flux_insulator[1] = gi * b;
    r.d_flux_semiconductor[0] = gs * a;
    r.d_flux_semiconductor[1] = gs * (b - 1.0);
    r.trapped = q - e.charge_weight * e.fixed_charge_cm2;
    return r;
}

InterfaceEdges::Charge InterfaceEdges::charge_at(const InterfaceEdge& e, double n, double p,
                                                 const InterfaceStatistics& s) const {
    const physics::Degeneracy gn = s.fermi_dirac
                                       ? physics::fermi_dirac_degeneracy(s.n_ie, s.log_dos_n, n)
                                       : physics::Degeneracy{0.0, 0.0};
    const physics::Degeneracy gp = s.fermi_dirac
                                       ? physics::fermi_dirac_degeneracy(s.n_ie, s.log_dos_p, p)
                                       : physics::Degeneracy{0.0, 0.0};
    double q = e.fixed_charge_cm2, q_n = 0.0, q_p = 0.0, r = 0.0, r_n = 0.0, r_p = 0.0;
    for (std::size_t k = e.first_level; k < e.last_level; ++k) {
        const InterfaceLevel& l = levels_[k];
        const double n1 = s.n_ie * std::exp(gn.log_gamma + l.tau);
        const double p1 = s.n_ie * std::exp(gp.log_gamma - l.tau);
        const physics::TrapKinetics t = physics::trap_kinetics(
            n, p, n1, n1 * gn.d_density, p1, p1 * gp.d_density, l.cn, l.cp);
        // Donor N (1 - f), acceptor -N f; d(1 - f) = -df.
        q += l.density_cm2 * (l.donor ? t.empty : -t.occupied);
        q_n -= l.density_cm2 * t.d_occupied_dn;
        q_p -= l.density_cm2 * t.d_occupied_dp;
        r += l.density_cm2 * t.rate;
        r_n += l.density_cm2 * t.d_rate_dn;
        r_p += l.density_cm2 * t.d_rate_dp;
    }
    if (e.velocity_n > 0.0 && e.velocity_p > 0.0) {
        const physics::EquilibriumProduct E =
            s.fermi_dirac ? physics::fermi_dirac_equilibrium_product(s.n_ie, gn, gp)
                          : physics::boltzmann_equilibrium_product(s.n_ie);
        const physics::RecombinationRate v = physics::srh_recombination(
            n, p, E, s.n_ie, 1.0 / e.velocity_n, 1.0 / e.velocity_p);
        r += v.rate;
        r_n += v.d_dn;
        r_p += v.d_dp;
    }
    const double cw = e.charge_weight, rw = e.rate_weight;
    return {cw * q, cw * q_n, cw * q_p, rw * r, rw * r_n, rw * r_p};
}

InterfaceDrift InterfaceEdges::drift(const InterfaceEdge& e, double psi_i, double psi_s,
                                     double n_s, double p_s,
                                     const InterfaceStatistics& s) const {
    // The node's reduced energies (n = gamma n_ie e^xi) and d xi / d density.
    const physics::Degeneracy gn = s.fermi_dirac
                                       ? physics::fermi_dirac_degeneracy(s.n_ie, s.log_dos_n, n_s)
                                       : physics::Degeneracy{0.0, 0.0};
    const physics::Degeneracy gp = s.fermi_dirac
                                       ? physics::fermi_dirac_degeneracy(s.n_ie, s.log_dos_p, p_s)
                                       : physics::Degeneracy{0.0, 0.0};
    const double xi_n = std::log(n_s / s.n_ie) - gn.log_gamma;
    const double xi_p = std::log(p_s / s.n_ie) - gp.log_gamma;
    const double dxi_n = 1.0 / n_s - gn.d_density;
    const double dxi_p = 1.0 / p_s - gp.d_density;
    // The interface densities at psi_I, with d / d eta.
    struct State {
        physics::DensityResult n, p;
        Charge c;
    };
    const auto at = [&](double psi) {
        State st;
        const double shift = psi - psi_s;
        st.n = detail::density(s.fermi_dirac, s.n_ie, s.log_dos_n, xi_n + shift);
        st.p = detail::density(s.fermi_dirac, s.n_ie, s.log_dos_p, xi_p - shift);
        st.c = charge_at(e, st.n.density, st.p.density, s);
        return st;
    };
    // dQ/dpsi_I = Q_n dn_I/dpsi_I + Q_p dp_I/dpsi_I = Q_n Nd - Q_p Pd.
    const auto Q = [&](double psi) {
        const State st = at(psi);
        return std::pair{st.c.value, st.c.d_n * st.n.d_eta - st.c.d_p * st.p.d_eta};
    };
    const double gi = e.g_insulator, gs = e.g_semiconductor, G = gi + gs;
    const double psi = interface_root(G, gi * psi_i + gs * psi_s, e.charge_low, e.charge_high, Q);
    const State st = at(psi);
    const double Nd = st.n.d_eta, Pd = st.p.d_eta;
    // Partials of n_I and p_I holding psi_I, columns (psi_i, psi_s, n_s, p_s).
    const double pn[4] = {0.0, -Nd, Nd * dxi_n, 0.0};
    const double pp[4] = {0.0, Pd, 0.0, Pd * dxi_p};
    // F = G psi_I - g_i psi_i - g_s psi_s - Q: d psi_I / dx = -F_x / D, D = G - dQ/dpsi_I.
    const double D = G - (st.c.d_n * Nd - st.c.d_p * Pd);
    const double unit_i[4] = {1.0, 0.0, 0.0, 0.0}, unit_s[4] = {0.0, 1.0, 0.0, 0.0};
    InterfaceDrift r{};
    r.psi = psi;
    r.flux_insulator = gi * (psi - psi_i);
    r.flux_semiconductor = gs * (psi - psi_s);
    r.rate = st.c.rate;
    for (std::size_t k = 0; k < 4; ++k) {
        const double Fx = -gi * unit_i[k] - gs * unit_s[k] - (st.c.d_n * pn[k] + st.c.d_p * pp[k]);
        const double dpsi = -Fx / D;
        r.d_flux_insulator[k] = gi * (dpsi - unit_i[k]);
        r.d_flux_semiconductor[k] = gs * (dpsi - unit_s[k]);
        const double dn = pn[k] + Nd * dpsi, dp = pp[k] - Pd * dpsi;
        r.d_rate[k] = st.c.rate_n * dn + st.c.rate_p * dp;
    }
    r.trapped = st.c.value - e.charge_weight * e.fixed_charge_cm2;
    return r;
}

}  // namespace NiTCAD::assemble
