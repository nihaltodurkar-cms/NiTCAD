#include "NiTCAD/assemble/interface_edges.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

#include "NiTCAD/base/contract.hpp"
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
    : edges_(std::move(edges)), levels_(std::move(levels)), interfaces_(interfaces) {
    slot_offset_.reserve(edges_.size());
    for (const InterfaceEdge& e : edges_) {
        slot_offset_.push_back(slots_);
        slots_ += e.last_level - e.first_level;
    }
}

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
                                                 const InterfaceStatistics& s,
                                                 const TrapStep* step,
                                                 std::span<double> f) const {
    const physics::Degeneracy gn = s.fermi_dirac
                                       ? physics::fermi_dirac_degeneracy(s.n_ie, s.log_dos_n, n)
                                       : physics::Degeneracy{0.0, 0.0};
    const physics::Degeneracy gp = s.fermi_dirac
                                       ? physics::fermi_dirac_degeneracy(s.n_ie, s.log_dos_p, p)
                                       : physics::Degeneracy{0.0, 0.0};
    double q = e.fixed_charge_cm2, q_n = 0.0, q_p = 0.0, r = 0.0, r_n = 0.0, r_p = 0.0;
    double h = 0.0, h_n = 0.0, h_p = 0.0;  // the hole row's rate, in a time step
    if (step != nullptr) {
        // f = (c + w occ) / (1 + w D); the electron row takes cn (n - f (n + n1)), the hole row
        // cp (f (p + p1) - p1), and Q falls with f as in steady state.
        const double w = step->weight;
        for (std::size_t k = e.first_level; k < e.last_level; ++k) {
            const InterfaceLevel& l = levels_[k];
            const double n1 = s.n_ie * std::exp(gn.log_gamma + l.tau);
            const double p1 = s.n_ie * std::exp(gp.log_gamma - l.tau);
            const double dn1 = n1 * gn.d_density, dp1 = p1 * gp.d_density;
            const double occ = l.cn * n + l.cp * p1;
            const double D = l.cn * (n + n1) + l.cp * (p + p1);
            const double den = 1.0 + w * D;
            const double occupied = (step->history[k - e.first_level] + w * occ) / den;
            const double f_n = w * (l.cn - occupied * l.cn * (1.0 + dn1)) / den;
            const double f_p = w * (l.cp * dp1 - occupied * l.cp * (1.0 + dp1)) / den;
            if (!f.empty()) f[k - e.first_level] = occupied;
            const double N = l.density_cm2;
            q += N * (l.donor ? 1.0 - occupied : -occupied);
            q_n -= N * f_n;
            q_p -= N * f_p;
            r += N * l.cn * (n - occupied * (n + n1));
            r_n += N * l.cn * (1.0 - f_n * (n + n1) - occupied * (1.0 + dn1));
            r_p -= N * l.cn * f_p * (n + n1);
            h += N * l.cp * (occupied * (p + p1) - p1);
            h_n += N * l.cp * f_n * (p + p1);
            h_p += N * l.cp * (f_p * (p + p1) + occupied * (1.0 + dp1) - dp1);
        }
    }
    for (std::size_t k = e.first_level; step == nullptr && k < e.last_level; ++k) {
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
        if (!f.empty()) f[k - e.first_level] = t.occupied;
    }
    if (step == nullptr) {
        h = r;
        h_n = r_n;
        h_p = r_p;
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
        h += v.rate;
        h_n += v.d_dn;
        h_p += v.d_dp;
    }
    const double cw = e.charge_weight, rw = e.rate_weight;
    return {cw * q, cw * q_n, cw * q_p, rw * r, rw * r_n, rw * r_p, rw * h, rw * h_n, rw * h_p};
}

InterfaceEdges::ChargeSmallSignal InterfaceEdges::charge_small_signal(
    const InterfaceEdge& e, double n, double p, const InterfaceStatistics& s,
    std::complex<double> inv_weight) const {
    using Complex = std::complex<double>;
    const physics::Degeneracy gn = s.fermi_dirac
                                       ? physics::fermi_dirac_degeneracy(s.n_ie, s.log_dos_n, n)
                                       : physics::Degeneracy{0.0, 0.0};
    const physics::Degeneracy gp = s.fermi_dirac
                                       ? physics::fermi_dirac_degeneracy(s.n_ie, s.log_dos_p, p)
                                       : physics::Degeneracy{0.0, 0.0};
    Complex q_n, q_p, r_n, r_p, h_n, h_p;
    for (std::size_t k = e.first_level; k < e.last_level; ++k) {
        // As charge_at's time step, linearized about the steady occupancy f: the occupancy moves
        // by df = (f_n dn + f_p dp), (1 / w + D) f_n = cn - f cn (1 + dn1), likewise f_p.
        const InterfaceLevel& l = levels_[k];
        const double n1 = s.n_ie * std::exp(gn.log_gamma + l.tau);
        const double p1 = s.n_ie * std::exp(gp.log_gamma - l.tau);
        const double dn1 = n1 * gn.d_density, dp1 = p1 * gp.d_density;
        const double D = l.cn * (n + n1) + l.cp * (p + p1);
        const double f = physics::trap_kinetics(n, p, n1, dn1, p1, dp1, l.cn, l.cp).occupied;
        const Complex den = inv_weight + D;
        const Complex f_n = (l.cn - f * l.cn * (1.0 + dn1)) / den;
        const Complex f_p = (l.cp * dp1 - f * l.cp * (1.0 + dp1)) / den;
        const double N = l.density_cm2;
        q_n -= N * f_n;
        q_p -= N * f_p;
        r_n += N * l.cn * (1.0 - f_n * (n + n1) - f * (1.0 + dn1));
        r_p -= N * l.cn * f_p * (n + n1);
        h_n += N * l.cp * f_n * (p + p1);
        h_p += N * l.cp * (f_p * (p + p1) + f * (1.0 + dp1) - dp1);
    }
    if (e.velocity_n > 0.0 && e.velocity_p > 0.0) {
        const physics::EquilibriumProduct E =
            s.fermi_dirac ? physics::fermi_dirac_equilibrium_product(s.n_ie, gn, gp)
                          : physics::boltzmann_equilibrium_product(s.n_ie);
        const physics::RecombinationRate v = physics::srh_recombination(
            n, p, E, s.n_ie, 1.0 / e.velocity_n, 1.0 / e.velocity_p);
        r_n += v.d_dn;
        r_p += v.d_dp;
        h_n += v.d_dn;
        h_p += v.d_dp;
    }
    const double cw = e.charge_weight, rw = e.rate_weight;
    return {cw * q_n, cw * q_p, rw * r_n, rw * r_p, rw * h_n, rw * h_p};
}

InterfaceSmallSignal InterfaceEdges::small_signal(const InterfaceEdge& e, double psi_i,
                                                  double psi_s, double n_s, double p_s,
                                                  const InterfaceStatistics& s,
                                                  std::complex<double> inv_weight) const {
    using Complex = std::complex<double>;
    // The steady interface potential, then drift_at's implicit derivatives with the complex
    // partials of the charge and rates.
    const double psi = drift_at(e, psi_i, psi_s, n_s, p_s, s, nullptr, {}).psi;
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
    const physics::DensityResult nI =
        detail::density(s.fermi_dirac, s.n_ie, s.log_dos_n, xi_n + (psi - psi_s));
    const physics::DensityResult pI =
        detail::density(s.fermi_dirac, s.n_ie, s.log_dos_p, xi_p - (psi - psi_s));
    const ChargeSmallSignal c = charge_small_signal(e, nI.density, pI.density, s, inv_weight);
    const double Nd = nI.d_eta, Pd = pI.d_eta;
    const double pn[4] = {0.0, -Nd, Nd * dxi_n, 0.0};
    const double pp[4] = {0.0, Pd, 0.0, Pd * dxi_p};
    const double gi = e.g_insulator, gs = e.g_semiconductor, G = gi + gs;
    const Complex D = G - (c.d_n * Nd - c.d_p * Pd);
    const double unit_i[4] = {1.0, 0.0, 0.0, 0.0}, unit_s[4] = {0.0, 1.0, 0.0, 0.0};
    InterfaceSmallSignal r{};
    for (std::size_t k = 0; k < 4; ++k) {
        const Complex Fx = -gi * unit_i[k] - gs * unit_s[k] - (c.d_n * pn[k] + c.d_p * pp[k]);
        const Complex dpsi = -Fx / D;
        r.d_flux_insulator[k] = gi * (dpsi - unit_i[k]);
        r.d_flux_semiconductor[k] = gs * (dpsi - unit_s[k]);
        const Complex dn = pn[k] + Nd * dpsi, dp = pp[k] - Pd * dpsi;
        r.d_rate[k] = c.rate_n * dn + c.rate_p * dp;
        r.d_rate_p[k] = c.rate_hn * dn + c.rate_hp * dp;
    }
    return r;
}

InterfaceDrift InterfaceEdges::drift(const InterfaceEdge& e, double psi_i, double psi_s,
                                     double n_s, double p_s, const InterfaceStatistics& s,
                                     const TrapStep* step) const {
    return drift_at(e, psi_i, psi_s, n_s, p_s, s, step, {});
}

void InterfaceEdges::occupancies(const InterfaceEdge& e, double psi_i, double psi_s, double n_s,
                                 double p_s, const InterfaceStatistics& s, const TrapStep* step,
                                 std::span<double> f) const {
    NITCAD_EXPECTS(f.size() == e.last_level - e.first_level);
    (void)drift_at(e, psi_i, psi_s, n_s, p_s, s, step, f);
}

InterfaceDrift InterfaceEdges::drift_at(const InterfaceEdge& e, double psi_i, double psi_s,
                                        double n_s, double p_s, const InterfaceStatistics& s,
                                        const TrapStep* step, std::span<double> f) const {
    NITCAD_EXPECTS(step == nullptr || step->history.size() == e.last_level - e.first_level);
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
        st.c = charge_at(e, st.n.density, st.p.density, s, step);
        return st;
    };
    // dQ/dpsi_I = Q_n dn_I/dpsi_I + Q_p dp_I/dpsi_I = Q_n Nd - Q_p Pd.
    const auto Q = [&](double psi) {
        const State st = at(psi);
        return std::pair{st.c.value, st.c.d_n * st.n.d_eta - st.c.d_p * st.p.d_eta};
    };
    const double gi = e.g_insulator, gs = e.g_semiconductor, G = gi + gs;
    // In a time step f lies within [min(c, 0), max(c, 1)], so Q within the bounds of those.
    double low = e.charge_low, high = e.charge_high;
    if (step != nullptr) {
        double ql = e.fixed_charge_cm2, qh = e.fixed_charge_cm2;
        for (std::size_t k = e.first_level; k < e.last_level; ++k) {
            const InterfaceLevel& l = levels_[k];
            const double c = step->history[k - e.first_level];
            const double f_low = std::min(c, 0.0), f_high = std::max(c, 1.0);
            ql += l.density_cm2 * (l.donor ? 1.0 - f_high : -f_high);
            qh += l.density_cm2 * (l.donor ? 1.0 - f_low : -f_low);
        }
        low = e.charge_weight * ql;
        high = e.charge_weight * qh;
    }
    const double psi = interface_root(G, gi * psi_i + gs * psi_s, low, high, Q);
    State st = at(psi);
    if (!f.empty()) st.c = charge_at(e, st.n.density, st.p.density, s, step, f);
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
    r.rate_p = st.c.rate_h;
    for (std::size_t k = 0; k < 4; ++k) {
        const double Fx = -gi * unit_i[k] - gs * unit_s[k] - (st.c.d_n * pn[k] + st.c.d_p * pp[k]);
        const double dpsi = -Fx / D;
        r.d_flux_insulator[k] = gi * (dpsi - unit_i[k]);
        r.d_flux_semiconductor[k] = gs * (dpsi - unit_s[k]);
        const double dn = pn[k] + Nd * dpsi, dp = pp[k] - Pd * dpsi;
        r.d_rate[k] = st.c.rate_n * dn + st.c.rate_p * dp;
        r.d_rate_p[k] = st.c.rate_hn * dn + st.c.rate_hp * dp;
    }
    r.trapped = st.c.value - e.charge_weight * e.fixed_charge_cm2;
    return r;
}

}  // namespace NiTCAD::assemble
