// The electrothermal rows as equation terms (Unit 40 slice 2b; terms.hpp; Unit 23,
// drift_diffusion.hpp "Electrothermal coupling", DECISIONS.md T1-T14). Every node and edge block is
// 4 x 4 (psi, n, p, theta) by construction, so the terms index it as such.
#include <cmath>
#include <cstdint>
#include <utility>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/physics/impact_ionization.hpp"
#include "NiTCAD/physics/recombination.hpp"
#include "NiTCAD/physics/statistics.hpp"
#include "NiTCAD/physics/thermal.hpp"
#include "scaled_device.hpp"
#include "terms.hpp"

namespace NiTCAD::assemble {

namespace {

constexpr FieldId heat_only[] = {FieldId::temperature};
constexpr FieldId psi_only[] = {FieldId::potential};
constexpr FieldId psi_heat[] = {FieldId::potential, FieldId::temperature};
constexpr FieldId densities[] = {FieldId::electrons, FieldId::holes};
constexpr FieldId density_heat[] = {FieldId::electrons, FieldId::holes, FieldId::temperature};
constexpr FieldId all3[] = {FieldId::potential, FieldId::electrons, FieldId::holes};
constexpr FieldId all4[] = {FieldId::potential, FieldId::electrons, FieldId::holes,
                            FieldId::temperature};

}  // namespace

using Evaluation = DriftDiffusion::Terms::Evaluation;

void DriftDiffusion::Terms::thermal_product(const DriftDiffusion& s, const Evaluation& ev,
                                            std::size_t i, physics::EquilibriumProduct& E,
                                            double& dE_t) {
    const std::span<const ThermalState> t = ev.thermal;
    const double nie = t[i].nie, d_nie = t[i].d_nie;
    E = physics::boltzmann_equilibrium_product(nie);
    dE_t = 2.0 * nie * d_nie;
    if (s.fermi_dirac_) {
        const NodeLevels& l = ev.levels[i];
        E = physics::fermi_dirac_equilibrium_product(nie, {l.n.log_gamma, l.n.dg_dc},
                                                     {l.p.log_gamma, l.p.dg_dc});
        dE_t = E.value * (2.0 * d_nie / nie + l.n.dg_dt + l.p.dg_dt);
    }
}

void DriftDiffusion::Terms::add_thermal_rate(const DriftDiffusion& s, const Evaluation& ev,
                                             std::size_t i, double k, double rate, double d_n,
                                             double d_p, double d_t) {
    const std::span<double> f = ev.f, values = ev.values;
    const auto at = [&](std::size_t node, std::size_t r, std::size_t c) -> double& {
        return values[s.block_[16 * node + 4 * r + c]];
    };
    f[4 * i + 1] -= k * rate;
    f[4 * i + 2] += k * rate;
    if (ev.jacobian()) {
        at(i, 1, 1) -= k * d_n;
        at(i, 1, 2) -= k * d_p;
        at(i, 1, 3) -= k * d_t;
        at(i, 2, 1) += k * d_n;
        at(i, 2, 2) += k * d_p;
        at(i, 2, 3) += k * d_t;
    }
}

// Every node's heat row starts as its thermal contact's: an isothermal sink holds the node's rise
// at the contact's (assigned), a thermal resistance R_th adds h (tau - rise_c), any other node 0.
// With sinks false (thermal_contact_heat) every heat row starts at 0.
class DriftDiffusion::Terms::HeatSinks final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"heat sinks", TermDomain::nodes, heat_only, heat_only, true};
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        const std::span<double> f = ev.f, values = ev.values;
        const bool jacobian = ev.jacobian();
        const auto at = [&](std::size_t node, std::size_t r, std::size_t c) -> double& {
            return values[s.block_[16 * node + 4 * r + c]];
        };
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            const double tau = ev.x[4 * i + 3];
            const ThermalNode& tn = s.thermal_nodes_[i];
            if (held(s, ev, i)) {
                f[4 * i + 3] = tau - s.thermal_rise_[static_cast<std::size_t>(tn.contact)];
                if (jacobian) at(i, 3, 3) = 1.0;
            } else if (ev.sinks && tn.h > 0.0) {
                f[4 * i + 3] =
                    tn.h * (tau - s.thermal_rise_[static_cast<std::size_t>(tn.contact)]);
                if (jacobian) at(i, 3, 3) = tn.h;
            } else {
                f[4 * i + 3] = 0.0;
            }
        }
    }
};

// A time step: the lattice's heat storage rho c dT/dt in the heat row of every node not held.
class DriftDiffusion::Terms::LatticeStorage final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"lattice heat storage", TermDomain::nodes, heat_only, heat_only, false};
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        const TimeStep* step = ev.step;
        if (step == nullptr) return;
        const std::span<double> f = ev.f, values = ev.values;
        const bool jacobian = ev.jacobian();
        const auto at = [&](std::size_t node, std::size_t r, std::size_t c) -> double& {
            return values[s.block_[16 * node + 4 * r + c]];
        };
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (held(s, ev, i)) continue;
            const double tau = ev.x[4 * i + 3];
            const ThermalNode& tn = s.thermal_nodes_[i];
            f[4 * i + 3] += tn.capacity * step->rate * (tau - step->storage[3 * i + 2]);
            if (jacobian) at(i, 3, 3) += tn.capacity * step->rate;
        }
    }
};

// An ohmic contact node: psi, n and p held at the contact's values at the node's T (assigned).
class DriftDiffusion::Terms::ThermalContactRows final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"ohmic contact rows (electrothermal)", TermDomain::nodes, all3, all4, true};
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        const std::span<double> f = ev.f, values = ev.values;
        const bool jacobian = ev.jacobian();
        const auto at = [&](std::size_t node, std::size_t r, std::size_t c) -> double& {
            return values[s.block_[16 * node + 4 * r + c]];
        };
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (s.contact_[i] < 0) continue;
            const double psi = ev.x[4 * i], n = ev.x[4 * i + 1], p = ev.x[4 * i + 2];
            const ThermalOhmic o = s.thermal_ohmic(i, ev.thermal[i]);
            f[4 * i] = psi - o.psi;
            f[4 * i + 1] = n - o.n;
            f[4 * i + 2] = p - o.p;
            if (jacobian) {
                at(i, 0, 0) = at(i, 1, 1) = at(i, 2, 2) = 1.0;
                at(i, 0, 3) = -o.d_psi;
                at(i, 1, 3) = -o.d_n;
                at(i, 2, 3) = -o.d_p;
            }
        }
    }
};

// An insulator node: no charge, no carriers; an electrode holds psi at its T (assigned).
class DriftDiffusion::Terms::ThermalInsulatorRows final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"insulator rows (electrothermal)", TermDomain::nodes, all3, all4, true};
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        const std::span<double> f = ev.f, values = ev.values;
        const bool jacobian = ev.jacobian();
        const auto at = [&](std::size_t node, std::size_t r, std::size_t c) -> double& {
            return values[s.block_[16 * node + 4 * r + c]];
        };
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (s.contact_[i] >= 0 || s.insulator_[i] == 0) continue;
            const double psi = ev.x[4 * i], n = ev.x[4 * i + 1], p = ev.x[4 * i + 2];
            const bool electrode = s.electrode_[i] >= 0;
            const Shift e = electrode ? s.electrode_shift(i, ev.thermal[i]) : Shift{0.0, 0.0};
            f[4 * i] = electrode ? psi - s.psi0_[i] - e.value : 0.0;
            f[4 * i + 1] = n;
            f[4 * i + 2] = p;
            if (jacobian) {
                if (electrode) {
                    at(i, 0, 0) = 1.0;
                    at(i, 0, 3) = -e.d_theta;
                }
                at(i, 1, 1) = at(i, 2, 2) = 1.0;
            }
        }
    }
};

// A gate node's oxide term G (psi_G(T) - psi) + S (gate.hpp) in its Poisson row.
class DriftDiffusion::Terms::ThermalGate final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"gate (electrothermal)", TermDomain::nodes, psi_only, psi_heat, false};
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        const std::span<double> f = ev.f, values = ev.values;
        const bool jacobian = ev.jacobian();
        const auto at = [&](std::size_t node, std::size_t r, std::size_t c) -> double& {
            return values[s.block_[16 * node + 4 * r + c]];
        };
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (!semiconductor(s, i) || !s.gates_.on_gate(i)) continue;
            const double psi = ev.x[4 * i];
            const double G = s.gates_.term(i).coupling;
            const Shift e = s.electrode_shift(i, ev.thermal[i]);
            f[4 * i] += s.gates_.row_term(i, psi) + G * e.value;
            if (jacobian) {
                at(i, 0, 0) -= G;
                at(i, 0, 3) += G * e.d_theta;
            }
        }
    }
};

// SRH recombination at the node's T: R = (n p - E) / den, den = tau_p (n + n_ie) +
// tau_n (p + n_ie), d R / d theta = (-dE / den - (n p - E) (tau_p + tau_n) d_nie / den^2).
class DriftDiffusion::Terms::ThermalSrh final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"SRH recombination (electrothermal)", TermDomain::nodes, densities, density_heat,
                false};
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (!semiconductor(s, i)) continue;
            const double n = ev.x[4 * i + 1], p = ev.x[4 * i + 2];
            const double V = s.volume_[i];
            const double nie = ev.thermal[i].nie, d_nie = ev.thermal[i].d_nie;
            const double k = V * s.rate_scale_;
            physics::EquilibriumProduct E;
            double dE_t;
            thermal_product(s, ev, i, E, dE_t);
            const physics::RecombinationRate r =
                physics::srh_recombination(n, p, E, nie, s.tau_n_[i], s.tau_p_[i]);
            const double den = s.tau_p_[i] * (n + nie) + s.tau_n_[i] * (p + nie);
            const double d_t = -dE_t / den - (n * p - E.value) * (s.tau_p_[i] + s.tau_n_[i]) *
                                                  d_nie / (den * den);
            add_thermal_rate(s, ev, i, k, r.rate, r.d_dn, r.d_dp, d_t);
        }
    }
};

// Auger recombination on physical densities (E Ns^2, its partials times Ns): dR / dE =
// -(Cn n + Cp p).
class DriftDiffusion::Terms::ThermalAuger final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"Auger recombination (electrothermal)", TermDomain::nodes, densities,
                density_heat, false};
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        const double Ns = s.Ns_;
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (!semiconductor(s, i)) continue;
            const double n = ev.x[4 * i + 1], p = ev.x[4 * i + 2];
            const double V = s.volume_[i];
            const double k = V * s.rate_scale_;
            physics::EquilibriumProduct E;
            double dE_t;
            thermal_product(s, ev, i, E, dE_t);
            const double ns = n * Ns, ps = p * Ns;
            const physics::EquilibriumProduct Ep{E.value * Ns * Ns, E.d_dn * Ns,
                                                 E.d_dp * Ns};
            const double dEp = dE_t * Ns * Ns;
            const physics::RecombinationRate r =
                physics::auger_recombination(ns, ps, Ep, s.auger_n_[i], s.auger_p_[i]);
            const double d_t = -(s.auger_n_[i] * ns + s.auger_p_[i] * ps) * dEp;
            add_thermal_rate(s, ev, i, k, r.rate / Ns, r.d_dn, r.d_dp, d_t / Ns);
        }
    }
};

// Radiative recombination on physical densities: dR / dE = -B.
class DriftDiffusion::Terms::ThermalRadiative final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"radiative recombination (electrothermal)", TermDomain::nodes, densities,
                density_heat, false};
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        const double Ns = s.Ns_;
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (!semiconductor(s, i) || !(s.radiative_[i] > 0.0)) continue;
            const double n = ev.x[4 * i + 1], p = ev.x[4 * i + 2];
            const double V = s.volume_[i];
            const double k = V * s.rate_scale_;
            physics::EquilibriumProduct E;
            double dE_t;
            thermal_product(s, ev, i, E, dE_t);
            const double ns = n * Ns, ps = p * Ns;
            const physics::EquilibriumProduct Ep{E.value * Ns * Ns, E.d_dn * Ns,
                                                 E.d_dp * Ns};
            const double dEp = dE_t * Ns * Ns;
            const physics::RecombinationRate r =
                physics::radiative_recombination(ns, ps, Ep, s.radiative_[i]);
            add_thermal_rate(s, ev, i, k, r.rate / Ns, r.d_dn, r.d_dp,
                             -s.radiative_[i] * dEp / Ns);
        }
    }
};

// A time step: the carriers' storage in the continuity rows (S_n = n, S_p = p, complete
// ionization, U6), and in the heat row of a node not held the energy they take with them,
// e_n V dn/dt - e_p V dp/dt over K0 (T4: the edges' energy flow counts the carriers a box stores
// as received; they are not heat), e at the node's own state.
class DriftDiffusion::Terms::ThermalCarrierStorage final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"carrier storage (electrothermal)", TermDomain::nodes, density_heat, all4, false};
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        const TimeStep* step = ev.step;
        if (step == nullptr) return;
        const std::span<double> f = ev.f, values = ev.values;
        const bool jacobian = ev.jacobian();
        const std::span<const ThermalState> t = ev.thermal;
        const std::span<const NodeLevels> g = ev.levels;
        const auto at = [&](std::size_t node, std::size_t r, std::size_t c) -> double& {
            return values[s.block_[16 * node + 4 * r + c]];
        };
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            if (!semiconductor(s, i)) continue;
            const double psi = ev.x[4 * i], n = ev.x[4 * i + 1], p = ev.x[4 * i + 2];
            const ThermalNode& tn = s.thermal_nodes_[i];
            const double V = s.volume_[i];
            const double ks = V * step->rate;
            const double dn = n - step->storage[3 * i], dp = p - step->storage[3 * i + 1];
            f[4 * i + 1] -= ks * dn;
            f[4 * i + 2] += ks * dp;
            if (jacobian) {
                at(i, 1, 1) -= ks;
                at(i, 2, 2) += ks;
            }
            if (held(s, ev, i)) continue;
            const Level none{};
            const Level& ln = g.empty() ? none : g[i].n;
            const Level& lp = g.empty() ? none : g[i].p;
            const double theta = t[i].theta;
            const double kn = tn.r_n + 2.5 + ln.extra, kp = tn.r_p + 2.5 + lp.extra;
            const double en = tn.c_n - psi + kn * theta;
            const double ep = tn.c_p - psi + t[i].gap - kp * theta;
            const double inv_K0 = 1.0 / s.heat_scale_;
            f[4 * i + 3] += inv_K0 * ks * (en * dn - ep * dp);
            if (jacobian) {
                at(i, 3, 0) += inv_K0 * ks * (dp - dn);
                at(i, 3, 1) += inv_K0 * ks * (en + dn * theta * ln.dh_dc);
                at(i, 3, 2) += inv_K0 * ks * (-ep + dp * theta * lp.dh_dc);
                at(i, 3, 3) += inv_K0 * ks *
                               (dn * (kn + theta * ln.dh_dt) -
                                dp * (t[i].d_gap - kp - theta * lp.dh_dt));
            }
        }
    }
};

// The declared interfaces with the electrothermal model (fixed charge and surface recombination;
// the model refuses traps, T8). interface_edges.hpp without traps: psi_I solves
// g_i (psi_I - psi_o) + g_s (psi_I - psi_s) = Q, Q = a Q_f / (Ns L_D) constant, so it does not
// depend on T; the half-edge fluxes go to the Poisson rows. Surface recombination at the
// semiconductor node's T: the carriers keep their quasi-Fermi levels across its half-cell, so the
// band reduced energies at the interface are x_n + (psi_I - psi_s) / theta and
// x_p - (psi_I - psi_s) / theta (at fixed x a density moves with N(theta), 1.5 / theta + d ln M /
// d theta in ln), and the rate is SRH with lifetimes 1 / s and the equilibrium product and
// n1 = n_ie of the node's T. Unlike the isothermal Interfaces, always in: the electrothermal
// small_signal_matrix adds no interface terms of its own.
class DriftDiffusion::Terms::ThermalInterfaces final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"interfaces (electrothermal)", TermDomain::interface_edges, all3, all4, false};
    }
    void pattern(const DriftDiffusion& s, std::vector<linalg::Triplet>& t) const override {
        interface_pattern(s, t);
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        if (s.interfaces_.empty()) return;
        const std::span<const double> x = ev.x;
        const std::span<const ThermalState> t = ev.thermal;
        const std::span<double> f = ev.f, values = ev.values;
        const bool jacobian = ev.jacobian();
        for (const InterfaceEdge& v : s.interfaces_.edges()) {
            const std::size_t o = v.insulator, sn = v.semiconductor;
            const EdgeTerm& e = s.edges_[v.edge];
            const ThermalEdge& te = s.thermal_edges_[v.edge];
            const auto pos = [&](std::size_t row, std::size_t r, std::size_t col,
                                 std::size_t c) {
                if (row == col) return s.block_[16 * row + 4 * r + c];
                return (row == e.a ? te.ab : te.ba)[4 * r + c];
            };
            const double gi = v.g_insulator, gs = v.g_semiconductor, G = gi + gs;
            const double psi_o = x[4 * o], psi_s = x[4 * sn];
            const double psi_I =
                (gi * psi_o + gs * psi_s + v.charge_weight * v.fixed_charge_cm2) / G;
            const double dI_o = gi / G, dI_s = gs / G;  // d psi_I / d psi_o, d psi_s
            if (s.electrode_[o] < 0 && s.contact_[o] < 0) {  // the insulator node's Poisson row
                f[4 * o] += gi * (psi_I - psi_o);
                if (jacobian) {
                    values[pos(o, 0, o, 0)] += gi * (dI_o - 1.0);
                    values[pos(o, 0, sn, 0)] += gi * dI_s;
                }
            }
            if (s.contact_[sn] >= 0) continue;  // a Dirichlet semiconductor node
            f[4 * sn] += gs * (psi_I - psi_s);
            if (jacobian) {
                values[pos(sn, 0, o, 0)] += gs * dI_o;
                values[pos(sn, 0, sn, 0)] += gs * (dI_s - 1.0);
            }
            if (!(v.velocity_n > 0.0 && v.velocity_p > 0.0)) continue;
            const ThermalState& ts = t[sn];
            const ThermalNode& tn = s.thermal_nodes_[sn];
            const double theta = ts.theta;
            const double bend = (psi_I - psi_s) / theta;
            // d bend / d (psi_o, psi_s, theta).
            const double db[3] = {dI_o / theta, (dI_s - 1.0) / theta, -bend / theta};
            const Level ln =
                s.level(x[4 * sn + 1], theta, ts.nie, tn.log_nc, ts.mass_n, ts.d_mass_n);
            const Level lp =
                s.level(x[4 * sn + 2], theta, ts.nie, tn.log_nv, ts.mass_p, ts.d_mass_p);
            const double lt = 1.5 * std::log(theta), ln_nie = std::log(ts.nie);
            const double gn_log = tn.log_nc + lt + ts.mass_n - ln_nie;
            const double gp_log = tn.log_nv + lt + ts.mass_p - ln_nie;
            const physics::DensityResult nI =
                detail::density(s.fermi_dirac_, ts.nie, gn_log, ln.x + bend + gn_log);
            const physics::DensityResult pI =
                detail::density(s.fermi_dirac_, ts.nie, gp_log, lp.x - bend + gp_log);
            // n_I and p_I in (psi_o, psi_s, n_s, p_s, theta).
            const double Nd = nI.d_eta, Pd = pI.d_eta;
            const double dn[5] = {Nd * db[0], Nd * db[1], Nd * ln.dx_dc, 0.0,
                                  1.5 * nI.density / theta + Nd * (ln.dx_dt + db[2]) +
                                      ts.d_mass_n * nI.density};
            const double dp[5] = {-Pd * db[0], -Pd * db[1], 0.0, Pd * lp.dx_dc,
                                  1.5 * pI.density / theta + Pd * (lp.dx_dt - db[2]) +
                                      ts.d_mass_p * pI.density};
            // The product at the interface densities and its direct theta derivative (n_ie and,
            // at fixed density, ln gamma through N(theta)).
            physics::EquilibriumProduct E = physics::boltzmann_equilibrium_product(ts.nie);
            double dE_t = 2.0 * ts.nie * ts.d_nie;
            if (s.fermi_dirac_) {
                const Level In =
                    s.level(nI.density, theta, ts.nie, tn.log_nc, ts.mass_n, ts.d_mass_n);
                const Level Ip =
                    s.level(pI.density, theta, ts.nie, tn.log_nv, ts.mass_p, ts.d_mass_p);
                E = physics::fermi_dirac_equilibrium_product(ts.nie, {In.log_gamma, In.dg_dc},
                                                             {Ip.log_gamma, Ip.dg_dc});
                dE_t = E.value * (2.0 * ts.d_nie / ts.nie + In.dg_dt + Ip.dg_dt);
            }
            const double tp = 1.0 / v.velocity_n, tq = 1.0 / v.velocity_p;  // as srh's tau_n, p
            const physics::RecombinationRate r =
                physics::srh_recombination(nI.density, pI.density, E, ts.nie, tp, tq);
            const double den = tq * (nI.density + ts.nie) + tp * (pI.density + ts.nie);
            const double R_t = -dE_t / den - (nI.density * pI.density - E.value) * (tq + tp) *
                                                 ts.d_nie / (den * den);
            const double rw = v.rate_weight;
            f[4 * sn + 1] -= rw * r.rate;
            f[4 * sn + 2] += rw * r.rate;
            if (!jacobian) continue;
            // Columns: psi_o, psi_s, n_s, p_s, theta_s.
            const std::size_t col_node[5] = {o, sn, sn, sn, sn};
            const std::size_t col[5] = {0, 0, 1, 2, 3};
            for (std::size_t c = 0; c < 5; ++c) {
                const double d = r.d_dn * dn[c] + r.d_dp * dp[c] + (c == 4 ? R_t : 0.0);
                if (d == 0.0) continue;
                values[pos(sn, 1, col_node[c], col[c])] -= rw * d;
                values[pos(sn, 2, col_node[c], col[c])] += rw * d;
            }
        }
    }
};

// Every mesh edge's terms with the electrothermal model: the Poisson flux (as the isothermal
// rows), the thermal Scharfetter-Gummel electron and hole fluxes (thermal_fluxes), and in the heat
// rows the conduction and the carriers' energy flow. One term: the heat rows take the carrier
// fluxes and their partials. Its pattern: the edges' full 4 x 4 blocks.
class DriftDiffusion::Terms::ThermalEdges final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"edge fluxes (electrothermal)", TermDomain::mesh_edges, all4, all4, false};
    }
    void pattern(const DriftDiffusion& s, std::vector<linalg::Triplet>& t) const override {
        for (const EdgeTerm& e : s.edges_) {
            for (std::size_t r = 0; r < 4; ++r) {
                for (std::size_t c = 0; c < 4; ++c) {
                    t.push_back({static_cast<linalg::Index>(4 * e.a + r),
                                 static_cast<linalg::Index>(4 * e.b + c), 0.0});
                    t.push_back({static_cast<linalg::Index>(4 * e.b + r),
                                 static_cast<linalg::Index>(4 * e.a + c), 0.0});
                }
            }
        }
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        const std::span<const double> x = ev.x;
        const std::span<const ThermalState> t = ev.thermal;
        const std::span<const NodeLevels> g = ev.levels;
        const std::span<double> f = ev.f, values = ev.values;
        const bool jacobian = ev.jacobian();
        const auto held = [&](std::size_t i) { return Terms::held(s, ev, i); };
        const double inv_K0 = 1.0 / s.heat_scale_;
        for (std::size_t k = 0; k < s.edges_.size(); ++k) {
            const EdgeTerm& e = s.edges_[k];
            const ThermalEdge& te = s.thermal_edges_[k];
            const std::size_t a = e.a, b = e.b;
            const std::size_t ends[2] = {a, b};
            // values position of (row of node `row`, r) with (column of node `col`, c).
            const auto pos = [&](std::size_t row, std::size_t r, std::size_t col,
                                 std::size_t c) {
                if (row == col) return s.block_[16 * row + 4 * r + c];
                return (row == a ? te.ab : te.ba)[4 * r + c];
            };
            // Row (node, r) gains sign * (value, gradient over the ends' (psi, n, p, theta)).
            const auto add_row = [&](std::size_t node, std::size_t r, double sign, double value,
                                     const double (&grad)[2][4]) {
                f[4 * node + r] += sign * value;
                if (!jacobian) return;
                for (std::size_t end = 0; end < 2; ++end) {
                    for (std::size_t c = 0; c < 4; ++c) {
                        if (grad[end][c] == 0.0) continue;
                        values[pos(node, r, ends[end], c)] += sign * grad[end][c];
                    }
                }
            };
            // Conduction, on every edge.
            const Conduction q = s.conduction(k, t);
            double heat[2][4] = {{0, 0, 0, q.d_a}, {0, 0, 0, q.d_b}};
            double heat_value = q.flux;
            // Poisson, as the isothermal rows.
            const double poisson = e.c * (x[4 * b] - x[4 * a]);
            const double dpoisson[2][4] = {{-e.c, 0, 0, 0}, {e.c, 0, 0, 0}};
            if (!e.carriers) {
                // A charged interface edge takes the half-edge fluxes instead
                // (ThermalInterfaces).
                if (!e.charged && s.contact_[a] < 0 && s.electrode_[a] < 0) {
                    add_row(a, 0, 1.0, poisson, dpoisson);
                }
                if (!e.charged && s.contact_[b] < 0 && s.electrode_[b] < 0) {
                    add_row(b, 0, -1.0, poisson, dpoisson);
                }
                if (!held(a)) add_row(a, 3, 1.0, heat_value, heat);
                if (!held(b)) add_row(b, 3, -1.0, heat_value, heat);
                continue;
            }
            const auto [fn, fp] = s.thermal_fluxes(k, x, t, g);
            const double gn[2][4] = {{fn.d_psi[0], fn.d_c[0], 0, fn.d_theta[0]},
                                     {fn.d_psi[1], fn.d_c[1], 0, fn.d_theta[1]}};
            const double gp[2][4] = {{fp.d_psi[0], 0, fp.d_c[0], fp.d_theta[0]},
                                     {fp.d_psi[1], 0, fp.d_c[1], fp.d_theta[1]}};
            for (const std::size_t node : ends) {
                if (s.contact_[node] >= 0) continue;
                const double sign = node == a ? 1.0 : -1.0;
                add_row(node, 0, sign, poisson, dpoisson);
                add_row(node, 1, sign, fn.flux, gn);
                add_row(node, 2, sign, fp.flux, gp);
            }
            // Energy: the mean carrier energies e_n = c_n - psi + (r_n + 5/2 + dh_n) theta and
            // e_p = c_p - psi + g - (r_p + 5/2 + dh_p) theta of the ends (dh the Fermi-Dirac
            // excess of the thermal diffusion factor, 0 under Boltzmann statistics).
            double en = 0.0, ep = 0.0, gen[2][4] = {}, gep[2][4] = {};
            for (std::size_t end = 0; end < 2; ++end) {
                const std::size_t i = ends[end];
                const ThermalNode& tn = s.thermal_nodes_[i];
                const double psi = x[4 * i], theta = t[i].theta;
                const Level none{};
                const Level& ln = g.empty() ? none : g[i].n;
                const Level& lp = g.empty() ? none : g[i].p;
                const double kn = tn.r_n + 2.5 + ln.extra, kp = tn.r_p + 2.5 + lp.extra;
                en += 0.5 * (tn.c_n - psi + kn * theta);
                ep += 0.5 * (tn.c_p - psi + t[i].gap - kp * theta);
                gen[end][0] = gep[end][0] = -0.5;
                gen[end][1] = 0.5 * theta * ln.dh_dc;
                gep[end][2] = -0.5 * theta * lp.dh_dc;
                gen[end][3] = 0.5 * (kn + theta * ln.dh_dt);
                gep[end][3] = 0.5 * (t[i].d_gap - kp - theta * lp.dh_dt);
            }
            // Heat row of each end: + conduction - ((e_n + m) Jn + (e_p + m) Jp) / K0, with the
            // edge's sign for its end, m the end's metal potential.
            for (const std::size_t node : ends) {
                if (held(node)) continue;
                const double sign = node == a ? 1.0 : -1.0;
                const double m = s.metal_[node];
                double grad[2][4];
                for (std::size_t end = 0; end < 2; ++end) {
                    for (std::size_t c = 0; c < 4; ++c) {
                        grad[end][c] = heat[end][c] -
                                       inv_K0 * ((en + m) * gn[end][c] + fn.flux * gen[end][c] +
                                                 (ep + m) * gp[end][c] + fp.flux * gep[end][c]);
                    }
                }
                add_row(node, 3, sign,
                        heat_value - inv_K0 * ((en + m) * fn.flux + (ep + m) * fp.flux), grad);
            }
        }
    }
};

// Impact ionization with the electrothermal model: the ImpactGeneration term with the thermal
// fluxes, whose temperature partials enter through the currents, and gamma at the node's own T:
// d alpha / d gamma = (alpha - E d alpha / dE) / gamma for alpha = gamma A exp(-gamma B / E).
class DriftDiffusion::Terms::ThermalImpact final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"impact ionization (electrothermal)", TermDomain::global, densities, all4, false};
    }
    void pattern(const DriftDiffusion& s, std::vector<linalg::Triplet>& t) const override {
        impact_pattern(s, t);
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        if (s.impact_nodes_.empty()) return;
        const std::span<const double> x = ev.x;
        const std::span<const ThermalState> t = ev.thermal;
        const std::span<const NodeLevels> g = ev.levels;
        const std::span<double> f = ev.f, values = ev.values;
        const bool jacobian = ev.jacobian();
        const int D = s.dimension_;
        const double T0 = s.T0_;
        std::vector<std::pair<ThermalFlux, ThermalFlux>> flux(s.edges_.size());
        for (std::size_t k = 0; k < s.edges_.size(); ++k) {
            if (s.edges_[k].carriers) flux[k] = s.thermal_fluxes(k, x, t, g);
        }
        constexpr double q = base::q_C;
        constexpr double floor2 = impact_current_floor * impact_current_floor;
        const double rho2 = s.impact_resolution_ * s.impact_resolution_;
        for (const ImpactNode& node : s.impact_nodes_) {
            const std::size_t i = node.node;
            const double count = static_cast<double>(node.last - node.first);
            const double T = t[i].theta * T0;
            const double gamma =
                physics::impact_ionization_temperature_factor(node.phonon_eV, T);
            const double d_gamma =
                T0 * physics::impact_ionization_temperature_factor_slope(node.phonon_eV, T);
            double E[3] = {}, jn[3] = {}, jp[3] = {}, R[2] = {0.0, 0.0};
            for (std::size_t k = node.first; k < node.last; ++k) {
                const ImpactEdge& ie = s.impact_edges_[k];
                const EdgeTerm& e = s.edges_[ie.edge];
                const double delta = x[4 * e.b] - x[4 * e.a];
                const double Ek = -ie.field * delta;
                const double jnk = ie.current * flux[ie.edge].first.flux;
                const double jpk = ie.current * flux[ie.edge].second.flux;
                const double root = std::sqrt(1.0 + delta * delta);
                R[0] += ie.current * e.an * (x[4 * e.a + 1] + x[4 * e.b + 1]) * root / count;
                R[1] += ie.current * e.ap * (x[4 * e.a + 2] + x[4 * e.b + 2]) * root / count;
                for (int r = 0; r < D; ++r) {
                    double w = 0.0;
                    for (int c = 0; c < D; ++c) w += node.inverse[r * D + c] * ie.t[c];
                    E[r] += w * Ek;
                    jn[r] += w * jnk;
                    jp[r] += w * jpk;
                }
            }
            double G = 0.0, dG_gamma = 0.0, dE[3] = {}, dj[2][3] = {}, dR[2] = {0.0, 0.0};
            for (int c = 0; c < 2; ++c) {
                const double* j = c == 0 ? jn : jp;
                double mag2 = 0.0, Ej = 0.0;
                for (int r = 0; r < D; ++r) {
                    mag2 += j[r] * j[r];
                    Ej += E[r] * j[r];
                }
                const double eps2 = floor2 + rho2 * R[c] * R[c];
                const double m = std::sqrt(mag2 + eps2);
                const double Epar = std::abs(Ej) / m;
                const physics::ImpactIonizationRate a = physics::impact_ionization_coefficient(
                    c == 0 ? node.n : node.p, gamma, Epar);
                if (a.alpha == 0.0 && a.d_dE == 0.0) continue;
                const double w = mag2 / m;
                G += a.alpha * w / q;
                dG_gamma += (a.alpha - Epar * a.d_dE) / gamma * w / q;
                const double sign = Ej > 0.0 ? 1.0 : Ej < 0.0 ? -1.0 : 0.0;
                const double m3 = m * m * m;
                for (int r = 0; r < D; ++r) {
                    dE[r] += a.d_dE * w * sign * j[r] / m / q;
                    dj[c][r] = (a.d_dE * w * (sign * E[r] / m - std::abs(Ej) * j[r] / m3) +
                                a.alpha * j[r] * (mag2 + 2.0 * eps2) / m3) /
                               q;
                }
                dR[c] = -(a.d_dE * w * std::abs(Ej) + a.alpha * mag2) / m3 / q * rho2 * R[c];
            }
            const double k_row = s.volume_[i] * s.generation_scale_;
            f[4 * i + 1] += k_row * G;
            f[4 * i + 2] -= k_row * G;
            if (!jacobian) continue;
            // Row n of node i gets +v, row p -v, at column comp of node_j through edge `edge`.
            const auto add_entry = [&](std::size_t node_j, std::size_t comp, std::size_t edge,
                                       double v) {
                const EdgeTerm& e = s.edges_[edge];
                const ThermalEdge& te = s.thermal_edges_[edge];
                for (const std::size_t row : {std::size_t{1}, std::size_t{2}}) {
                    const std::size_t pos =
                        node_j == i ? s.block_[16 * i + 4 * row + comp]
                                    : (e.a == i ? te.ab : te.ba)[4 * row + comp];
                    values[pos] += (row == 1 ? k_row : -k_row) * v;
                }
            };
            // gamma's own temperature dependence, on node i's theta column.
            values[s.block_[16 * i + 4 * 1 + 3]] += k_row * dG_gamma * d_gamma;
            values[s.block_[16 * i + 4 * 2 + 3]] -= k_row * dG_gamma * d_gamma;
            for (std::size_t k = node.first; k < node.last; ++k) {
                const ImpactEdge& ie = s.impact_edges_[k];
                const EdgeTerm& e = s.edges_[ie.edge];
                const ThermalFlux& fn = flux[ie.edge].first;
                const ThermalFlux& fp = flux[ie.edge].second;
                double wE = 0.0, wn = 0.0, wp = 0.0;
                for (int r = 0; r < D; ++r) {
                    double w = 0.0;
                    for (int c = 0; c < D; ++c) w += node.inverse[r * D + c] * ie.t[c];
                    wE += w * dE[r];
                    wn += w * dj[0][r];
                    wp += w * dj[1][r];
                }
                wn *= ie.current;
                wp *= ie.current;
                const double delta = x[4 * e.b] - x[4 * e.a];
                const double root = std::sqrt(1.0 + delta * delta);
                const double rn = ie.current * e.an / count, rp = ie.current * e.ap / count;
                const double sn = x[4 * e.a + 1] + x[4 * e.b + 1];
                const double sp = x[4 * e.a + 2] + x[4 * e.b + 2];
                const double dpsi = (dR[0] * rn * sn + dR[1] * rp * sp) * delta / root;
                add_entry(e.a, 0, ie.edge,
                          wE * ie.field + wn * fn.d_psi[0] + wp * fp.d_psi[0] - dpsi);
                add_entry(e.b, 0, ie.edge,
                          -wE * ie.field + wn * fn.d_psi[1] + wp * fp.d_psi[1] + dpsi);
                add_entry(e.a, 1, ie.edge, wn * fn.d_c[0] + dR[0] * rn * root);
                add_entry(e.b, 1, ie.edge, wn * fn.d_c[1] + dR[0] * rn * root);
                add_entry(e.a, 2, ie.edge, wp * fp.d_c[0] + dR[1] * rp * root);
                add_entry(e.b, 2, ie.edge, wp * fp.d_c[1] + dR[1] * rp * root);
                add_entry(e.a, 3, ie.edge, wn * fn.d_theta[0] + wp * fp.d_theta[0]);
                add_entry(e.b, 3, ie.edge, wn * fn.d_theta[1] + wp * fp.d_theta[1]);
            }
        }
    }
};

std::vector<std::shared_ptr<const DriftDiffusion::Term>> DriftDiffusion::Terms::make_thermal(
    const DriftDiffusion& s, const PhysicsModels& models) {
    std::vector<std::shared_ptr<const Term>> t;
    t.push_back(std::make_shared<HeatSinks>());
    t.push_back(std::make_shared<LatticeStorage>());
    t.push_back(std::make_shared<ThermalContactRows>());
    t.push_back(std::make_shared<ThermalInsulatorRows>());
    t.push_back(poisson_charge());  // complete ionization (U6)
    bool gated = false;
    for (std::size_t i = 0; i < s.node_count(); ++i) gated = gated || s.gates_.on_gate(i);
    if (gated) t.push_back(std::make_shared<ThermalGate>());
    if (models.srh) t.push_back(std::make_shared<ThermalSrh>());
    if (models.auger) t.push_back(std::make_shared<ThermalAuger>());
    if (models.radiative) t.push_back(std::make_shared<ThermalRadiative>());
    t.push_back(std::make_shared<ThermalCarrierStorage>());
    if (!s.interfaces_.edges().empty()) t.push_back(std::make_shared<ThermalInterfaces>());
    t.push_back(std::make_shared<ThermalEdges>());
    if (models.impact_ionization) t.push_back(std::make_shared<ThermalImpact>());
    return t;
}

}  // namespace NiTCAD::assemble
