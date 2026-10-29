// Energy-balance transport (Unit 52; DECISIONS.md ED1-ED13; models.energy_balance).
//
// Two more unknowns per node, the carrier temperature rises tau_n = (T_n - T_L) / T_L and tau_p
// (theta = 1 + tau), in the slots after psi, n, p (stride 5). In units of k T_L, densities and
// fluxes scaled as the continuity rows (scaling.hpp), with A the edge's Scharfetter-Gummel factor
// at the carrier-temperature mobility and theta_e the edge's mean theta:
// - Particle flux (ED2, Stratton with an energy-independent momentum relaxation time):
//   J_n = mu (n grad E_c + k grad(n T_n)) is, with w = theta n, a Scharfetter-Gummel flux in w
//   with the driving term delta / theta_e (delta the isothermal one, band and BGN steps included);
//       Fn = A (w_b B(d) - w_a B(-d)),   Fp = -A (w_b B(-d) - w_a B(d)),   d = delta / theta_e,
//   the isothermal fluxes at theta = 1. Fn is the electrons' flow into a, Fp minus the holes'.
// - Energy flux into a (ED3, c = 5/2, EQ3): S = -kappa grad T - c (k T / q) J with
//   kappa = c (k^2 / q) n mu T is, with v = theta^2 n, c times the same flux in v:
//       S_n = c A (v_b B(d) - v_a B(-d)),   S_p = c A (v_b B(-d) - v_a B(d)),
//   so a uniform theta carries c theta per particle.
// - Joule heating: an edge's carriers gain the band edge's drop times their flow,
//   -delta F (electrons and holes alike in this sign convention), half into each end's box, so
//   the summed sources are the total J.E exactly.
// - Rows (assigned by EnergyRows, then the edges add): on a semiconductor node
//       sum_edges (S + Joule / 2 ...) - V k [ (3/2) n tau / tau_w + (3/2) theta U ] = 0,
//   k = rate_scale, U the net recombination rate of the continuity rows (SRH, Auger, radiative),
//   n in the relaxation term floored at 1e-10 max(n, p, |N_D - N_A|) (ED14; p likewise), so a
//   carrier below that share of its node's majority or doping is held at T_L, not undetermined;
//   on an ohmic contact node and an insulator node tau = 0 (ED1). No other boundary carries
//   energy. While hold_temperatures is on, every energy row is tau - tau_held (ED12's first
//   stage) and the edges add to the continuity rows only.
// - Mobility (ED4): mu(T) of physics/energy_balance.hpp at the edge's mean rise, of the edge's
//   low-field mobility (between two materials the harmonic mean of each end's) with node a's
//   relaxation time and Canali parameters (each end's on a two-material edge).
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/bernoulli.hpp"
#include "NiTCAD/base/contract.hpp"
#include "NiTCAD/physics/energy_balance.hpp"
#include "NiTCAD/physics/recombination.hpp"
#include "composition.hpp"
#include "scaled_device.hpp"
#include "terms.hpp"

namespace NiTCAD::assemble {

namespace {

constexpr double energy_flux_factor = 2.5;  // c = 5/2 + r, r = 0 (EQ3)
constexpr FieldId all5[] = {FieldId::potential, FieldId::electrons, FieldId::holes,
                            FieldId::electron_temperature, FieldId::hole_temperature};
constexpr FieldId energies[] = {FieldId::electron_temperature, FieldId::hole_temperature};
constexpr FieldId transport_rows[] = {FieldId::electrons, FieldId::holes,
                                      FieldId::electron_temperature, FieldId::hole_temperature};

// The relaxation floor (ED14, ruled and amended 2026-10-29): each carrier relaxes with density
// max(c, eps max(n, p, |N_D - N_A|)), so a carrier below eps of its node's majority, or of its
// doping in depletion, is held at T_L; with it above, the physics is unchanged. The densities and
// their partials in n and p (0 where the doping sets the floor).
constexpr double relaxation_floor = 1e-10;

struct Relaxing {
    double n, n_n, n_p, p, p_n, p_p;
};

Relaxing relaxing(double n, double p, double doping) {
    const double d = std::abs(doping);
    Relaxing r{n, 1.0, 0.0, p, 0.0, 1.0};
    if (n < relaxation_floor * std::max(p, d)) {
        r.n = relaxation_floor * std::max(p, d);
        r.n_n = 0.0;
        r.n_p = p >= d ? relaxation_floor : 0.0;
    }
    if (p < relaxation_floor * std::max(n, d)) {
        r.p = relaxation_floor * std::max(n, d);
        r.p_p = 0.0;
        r.p_n = n >= d ? relaxation_floor : 0.0;
    }
    return r;
}

}  // namespace

std::optional<base::Error> DriftDiffusion::make_energy_balance(const device::Device& device,
                                                               const PhysicsModels& models) {
    (void)models;
    const auto refuse = [](std::string message) {
        return base::Error{base::ErrorCode::invalid_input, std::move(message), {}};
    };
    // ED11: the device's combinations.
    if (std::ranges::any_of(schottky_, [](std::int32_t c) { return c >= 0; })) {
        return refuse("energy balance with a Schottky contact is not supported");
    }
    if (!gate_tunnels_.empty()) return refuse("energy balance with gate tunnelling is not supported");
    for (const device::Interface& f : device.interfaces()) {
        if (f.transport == device::InterfaceTransport::thermionic_emission) {
            return refuse("energy balance with thermionic interfaces is not supported");
        }
        // Not in DECISIONS.md ED11 before code (found building slice 2): an interface's charge,
        // traps and surface recombination act on its edges' carrier fluxes, which take no
        // carrier temperature yet.
        if (device::has_interface_charge_or_recombination(f)) {
            return refuse("energy balance with interface charge, traps or surface recombination "
                          "is not supported");
        }
    }
    // ED5: each semiconductor node's relaxation times, its region's own or silicon's.
    const std::size_t n = node_count();
    relax_n_.assign(n, 0.0);
    relax_p_.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (insulator_[i] != 0) continue;
        const physics::SemiconductorParameters& p =
            device.material(static_cast<mesh::NodeId>(i)).parameters();
        physics::EnergyRelaxationParameters r = p.energy_relaxation;
        if (r.tau_n_s == 0.0 || r.tau_p_s == 0.0) {
            const bool silicon =
                p.electron_mobility == physics::silicon_parameters.electron_mobility &&
                p.hole_mobility == physics::silicon_parameters.hole_mobility;
            if (!silicon) {
                return refuse("energy balance: a semiconductor region without energy relaxation "
                              "times (only silicon's are built in)");
            }
            if (r.tau_n_s == 0.0) r.tau_n_s = physics::silicon_energy_relaxation.tau_n_s;
            if (r.tau_p_s == 0.0) r.tau_p_s = physics::silicon_energy_relaxation.tau_p_s;
        }
        relax_n_[i] = r.tau_n_s;
        relax_p_[i] = r.tau_p_s;
    }
    energy_edges_.resize(edges_.size());
    for (std::size_t k = 0; k < edges_.size(); ++k) {
        const EdgeTerm& e = edges_[k];
        for (std::size_t r = 0; r < 5; ++r) {
            for (std::size_t c = 0; c < 5; ++c) {
                energy_edges_[k].ab[5 * r + c] = detail::position(pattern_, 5 * e.a + r, 5 * e.b + c);
                energy_edges_[k].ba[5 * r + c] = detail::position(pattern_, 5 * e.b + r, 5 * e.a + c);
            }
        }
    }
    contact_bias_.assign(contact_count_, 0.0);  // set_bias fills it (electrical_power)
    energy_balance_ = true;
    return std::nullopt;
}

std::pair<DriftDiffusion::EnergyFlux, DriftDiffusion::EnergyFlux> DriftDiffusion::energy_fluxes(
    std::size_t k, std::span<const double> x) const {
    const EdgeTerm& e = edges_[k];
    const std::size_t a = e.a, b = e.b;
    const double dpsi = x[5 * b] - x[5 * a];
    const double na = x[5 * a + 1], nb = x[5 * b + 1];
    const double pa = x[5 * a + 2], pb = x[5 * b + 2];
    const double tna = x[5 * a + 3], tnb = x[5 * b + 3];
    const double tpa = x[5 * a + 4], tpb = x[5 * b + 4];
    const double TL = temperature_K_;
    // The edge factor A = (an / mu_edge) mu(T_e) and dA / dtau of either end (the mean's half).
    struct Factor {
        double A, d_t;
    };
    const auto factor = [&](double a0, double mu_edge, double mu_a, double mu_b, double rise,
                            double relax_a, double relax_b, const physics::CanaliParameters& sat_a,
                            const physics::CanaliParameters& sat_b) {
        const double scale = a0 / mu_edge;
        if (!e.mixed) {
            const auto m = physics::carrier_temperature_mobility(mu_edge, rise, relax_a, sat_a);
            return Factor{scale * m.mobility, scale * m.d_rise * 0.5 * TL};
        }
        const auto ma = physics::carrier_temperature_mobility(mu_a, rise, relax_a, sat_a);
        const auto mb = physics::carrier_temperature_mobility(mu_b, rise, relax_b, sat_b);
        const double sum = ma.mobility + mb.mobility;
        const double h = 2.0 * ma.mobility * mb.mobility / sum;
        const double dh = (2.0 * mb.mobility * mb.mobility * ma.d_rise +
                           2.0 * ma.mobility * ma.mobility * mb.d_rise) /
                          (sum * sum);
        return Factor{scale * h, scale * dh * 0.5 * TL};
    };
    // One carrier: sign +1 electrons, -1 holes (the hole flux is the electron form at -d,
    // negated). c_a, c_b its densities, t_a, t_b its rises, delta its driving term.
    const auto carrier = [&](double sign, double c_a, double c_b, double t_a, double t_b,
                             double delta, const Factor& F) {
        EnergyFlux r{};
        const double ta = 1.0 + t_a, tb = 1.0 + t_b;
        const double th = 1.0 + 0.5 * (t_a + t_b);
        const double d = sign * delta / th;
        const double bp = bernoulli(d), bm = bernoulli(-d);
        const double dp = bernoulli_derivative(d), dm = bernoulli_derivative(-d);
        const double wa = c_a * ta, wb = c_b * tb, va = wa * ta, vb = wb * tb;
        // In the electron form: P = w_b B(d) - w_a B(-d), Q likewise in v; their d slopes.
        const double P = wb * bp - wa * bm, Q = vb * bp - va * bm;
        const double Pd = wb * dp + wa * dm, Qd = vb * dp + va * dm;
        const double A = F.A, c = energy_flux_factor;
        const double dd_psi = sign / th;     // d d / d psi_b (minus at a)
        const double dd_t = -0.5 * d / th;  // d d / d tau of either end
        // Particle flux sign A P; energy flux into a c A Q (both carriers).
        r.flux = sign * A * P;
        r.d_psi[1] = sign * A * Pd * dd_psi;
        r.d_psi[0] = -r.d_psi[1];
        r.d_c[0] = sign * A * (-ta * bm);
        r.d_c[1] = sign * A * (tb * bp);
        r.d_t[0] = sign * (F.d_t * P + A * (-c_a * bm) + A * Pd * dd_t);
        r.d_t[1] = sign * (F.d_t * P + A * (c_b * bp) + A * Pd * dd_t);
        r.energy = c * A * Q;
        r.e_psi[1] = c * A * Qd * dd_psi;
        r.e_psi[0] = -r.e_psi[1];
        r.e_c[0] = c * A * (-ta * ta * bm);
        r.e_c[1] = c * A * (tb * tb * bp);
        r.e_t[0] = c * (F.d_t * Q + A * (-2.0 * c_a * ta * bm) + A * Qd * dd_t);
        r.e_t[1] = c * (F.d_t * Q + A * (2.0 * c_b * tb * bp) + A * Qd * dd_t);
        r.delta = delta;
        return r;
    };
    const Factor Fn = factor(e.an, e.mu_n, e.mu_n_a, e.mu_n_b, 0.5 * (tna + tnb) * TL,
                             relax_n_[a], relax_n_[b], e.sat_n, e.sat_n_b);
    const Factor Fp = factor(e.ap, e.mu_p, e.mu_p_a, e.mu_p_b, 0.5 * (tpa + tpb) * TL,
                             relax_p_[a], relax_p_[b], e.sat_p, e.sat_p_b);
    return {carrier(1.0, na, nb, tna, tnb, dpsi + e.shift_n, Fn),
            carrier(-1.0, pa, pb, tpa, tpb, dpsi + e.shift_p, Fp)};
}

void DriftDiffusion::hold_temperatures(std::span<const double> x) {
    NITCAD_EXPECTS(energy_balance_ && x.size() == unknowns());
    held_.resize(2 * node_count());
    for (std::size_t i = 0; i < node_count(); ++i) {
        held_[2 * i] = x[5 * i + 3];
        held_[2 * i + 1] = x[5 * i + 4];
    }
}

// Summing the energy rows over the semiconductor nodes off the contacts leaves every edge's Joule
// heat but the contact ends' halves, minus the energy flux into the contact nodes; the Joule heat
// -sum_e delta F is, by parts with delta = -(E_b - E_a) per carrier, -sum_i (E_c D_n + E_v D_p)
// over all nodes (D the flux sums of the continuity rows): the gap energy of the recombined pairs
// off the contacts, and at the contacts sum_c V_c I_c plus the Peltier term.
DriftDiffusion::EnergyPowers DriftDiffusion::energy_powers(std::span<const double> x) const {
    NITCAD_EXPECTS(energy_balance_ && held_.empty() && x.size() == unknowns());
    const std::size_t nodes = node_count();
    std::vector<double> f(unknowns());
    residual(x, f);
    const BandEdges bands = band_edges(x);
    EnergyPowers p{electrical_power(x), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    std::vector<double> dn(nodes, 0.0), dp(nodes, 0.0);
    for (std::size_t k = 0; k < edges_.size(); ++k) {
        const EdgeTerm& e = edges_[k];
        if (e.charged || !e.carriers) continue;
        const auto [fn, fp] = energy_fluxes(k, x);
        dn[e.a] += fn.flux;
        dn[e.b] -= fn.flux;
        dp[e.a] += fp.flux;
        dp[e.b] -= fp.flux;
        const double joule = -fn.delta * fn.flux - fp.delta * fp.flux;
        if (contact_[e.a] >= 0) {
            p.contact_flux += fn.energy + fp.energy;
            p.contact_joule += 0.5 * joule;
        }
        if (contact_[e.b] >= 0) {
            p.contact_flux -= fn.energy + fp.energy;
            p.contact_joule += 0.5 * joule;
        }
    }
    for (std::size_t i = 0; i < nodes; ++i) {
        if (insulator_[i] != 0) continue;
        if (contact_[i] >= 0) {
            const double V = contact_bias_[static_cast<std::size_t>(contact_[i])];
            p.peltier -= (bands.conduction[i] + V) * dn[i] + (bands.valence[i] + V) * dp[i];
            continue;
        }
        const double n = x[5 * i + 1], h = x[5 * i + 2];
        const double tn = x[5 * i + 3], tp = x[5 * i + 4];
        const double ku = dn[i] - f[5 * i + 1];  // V k U, the electron row's term
        const double k = volume_[i] * rate_scale_;
        const Relaxing rel = relaxing(n, h, doping_[i]);
        p.relaxation += 1.5 * k * (rel.n * tn / relax_n_[i] + rel.p * tp / relax_p_[i]);
        p.recombination +=
            ku * (bands.conduction[i] - bands.valence[i] + 1.5 * (2.0 + tn + tp));
        p.residual += f[5 * i + 3] + f[5 * i + 4] + bands.conduction[i] * f[5 * i + 1] +
                      bands.valence[i] * f[5 * i + 2];
    }
    return p;
}

// The energy rows of every node, assigned (ED1, ED3, ED7).
class DriftDiffusion::Terms::EnergyRows final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"energy rows", TermDomain::nodes, energies, all5, true};
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        const std::span<const double> x = ev.x;
        const bool jacobian = ev.jacobian();
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            const double tn = x[5 * i + 3], tp = x[5 * i + 4];
            if (!semiconductor(s, i) || !s.held_.empty()) {
                const double hn = s.held_.empty() ? 0.0 : s.held_[2 * i];
                const double hp = s.held_.empty() ? 0.0 : s.held_[2 * i + 1];
                ev.f[5 * i + 3] = tn - hn;
                ev.f[5 * i + 4] = tp - hp;
                if (jacobian) at(s, ev, i, 3, 3) = at(s, ev, i, 4, 4) = 1.0;
                continue;
            }
            const double n = x[5 * i + 1], p = x[5 * i + 2];
            const Relaxing rel = relaxing(n, p, s.doping_[i]);
            // U (the continuity rows' net recombination, in their units / k) and its partials.
            double u = 0.0, u_n = 0.0, u_p = 0.0;
            if (s.srh_) {
                const physics::RecombinationRate r = physics::srh_recombination(
                    n, p, product(s, ev, i), s.n_ie_[i], s.tau_n_[i], s.tau_p_[i]);
                u += r.rate;
                u_n += r.d_dn;
                u_p += r.d_dp;
            }
            const auto physical = [&](const physics::RecombinationRate& r) {
                u += r.rate / s.Ns_;
                u_n += r.d_dn;
                u_p += r.d_dp;
            };
            if (s.auger_) {
                physical(physics::auger_recombination(n * s.Ns_, p * s.Ns_,
                                                      physical_product(s, ev, i), s.auger_n_[i],
                                                      s.auger_p_[i]));
            }
            if (s.radiative_on_ && s.radiative_[i] > 0.0) {
                physical(physics::radiative_recombination(
                    n * s.Ns_, p * s.Ns_, physical_product(s, ev, i), s.radiative_[i]));
            }
            const double k = s.volume_[i] * s.rate_scale_;
            const double rn = 1.5 / s.relax_n_[i], rp = 1.5 / s.relax_p_[i];
            ev.f[5 * i + 3] = -k * (rn * rel.n * tn + 1.5 * (1.0 + tn) * u);
            ev.f[5 * i + 4] = -k * (rp * rel.p * tp + 1.5 * (1.0 + tp) * u);
            if (!jacobian) continue;
            at(s, ev, i, 3, 1) = -k * (rn * rel.n_n * tn + 1.5 * (1.0 + tn) * u_n);
            at(s, ev, i, 3, 2) = -k * (rn * rel.n_p * tn + 1.5 * (1.0 + tn) * u_p);
            at(s, ev, i, 3, 3) = -k * (rn * rel.n + 1.5 * u);
            at(s, ev, i, 4, 1) = -k * (rp * rel.p_n * tp + 1.5 * (1.0 + tp) * u_n);
            at(s, ev, i, 4, 2) = -k * (rp * rel.p_p * tp + 1.5 * (1.0 + tp) * u_p);
            at(s, ev, i, 4, 4) = -k * (rp * rel.p + 1.5 * u);
        }
    }
};

// The edges' particle and energy fluxes and Joule heating (ED2, ED3), adding to the continuity
// and energy rows of both ends off the ohmic contacts.
class DriftDiffusion::Terms::EnergyFluxes final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"energy-balance fluxes", TermDomain::mesh_edges, transport_rows, all5, false};
    }
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        const std::span<double> f = ev.f, values = ev.values;
        const bool jacobian = ev.jacobian();
        const bool held = !s.held_.empty();
        for (std::size_t k = 0; k < s.edges_.size(); ++k) {
            const EdgeTerm& e = s.edges_[k];
            if (e.charged || !e.carriers) continue;
            const EnergyEdge& pos = s.energy_edges_[k];
            const auto [fn, fp] = s.energy_fluxes(k, ev.x);
            // Per carrier: its continuity row, its energy row, its density and rise columns.
            const struct {
                const EnergyFlux& r;
                std::size_t row, erow, col, tcol;
            } carriers[] = {{fn, 1, 3, 1, 3}, {fp, 2, 4, 2, 4}};
            for (const auto& c : carriers) {
                const EnergyFlux& r = c.r;
                // Joule J = -delta flux / 2 into each end; its partials.
                const double joule = -0.5 * r.delta * r.flux;
                double j_psi[2], j_c[2], j_t[2];
                j_psi[0] = 0.5 * r.flux - 0.5 * r.delta * r.d_psi[0];
                j_psi[1] = -0.5 * r.flux - 0.5 * r.delta * r.d_psi[1];
                for (std::size_t q = 0; q < 2; ++q) {
                    j_c[q] = -0.5 * r.delta * r.d_c[q];
                    j_t[q] = -0.5 * r.delta * r.d_t[q];
                }
                // End a: +flux, +energy; end b: -flux, -energy; both +Joule. q = 0 is end a.
                for (std::size_t q = 0; q < 2; ++q) {
                    const std::size_t node = q == 0 ? e.a : e.b;
                    if (s.contact_[node] >= 0) continue;
                    const double sign = q == 0 ? 1.0 : -1.0;
                    f[5 * node + c.row] += sign * r.flux;
                    if (held) {  // the energy rows are tau - tau_held (EnergyRows)
                        if (!jacobian) continue;
                        const auto entry = [&](std::size_t o, std::size_t col) -> double& {
                            if (o == q) return at(s, ev, node, c.row, col);
                            return values[(q == 0 ? pos.ab : pos.ba)[5 * c.row + col]];
                        };
                        for (std::size_t o = 0; o < 2; ++o) {
                            entry(o, 0) += sign * r.d_psi[o];
                            entry(o, c.col) += sign * r.d_c[o];
                            entry(o, c.tcol) += sign * r.d_t[o];
                        }
                        continue;
                    }
                    f[5 * node + c.erow] += sign * r.energy + joule;
                    if (!jacobian) continue;
                    // Column of end o (0 = a, 1 = b) in row `row` of this node.
                    const auto entry = [&](std::size_t row, std::size_t o,
                                           std::size_t col) -> double& {
                        if (o == q) return at(s, ev, node, row, col);
                        return values[(q == 0 ? pos.ab : pos.ba)[5 * row + col]];
                    };
                    for (std::size_t o = 0; o < 2; ++o) {
                        entry(c.row, o, 0) += sign * r.d_psi[o];
                        entry(c.row, o, c.col) += sign * r.d_c[o];
                        entry(c.row, o, c.tcol) += sign * r.d_t[o];
                        entry(c.erow, o, 0) += sign * r.e_psi[o] + j_psi[o];
                        entry(c.erow, o, c.col) += sign * r.e_c[o] + j_c[o];
                        entry(c.erow, o, c.tcol) += sign * r.e_t[o] + j_t[o];
                    }
                }
            }
        }
    }
};

std::shared_ptr<const TermAccess::Term> DriftDiffusion::Terms::energy_rows() {
    return std::make_shared<const EnergyRows>();
}
std::shared_ptr<const TermAccess::Term> DriftDiffusion::Terms::energy_fluxes() {
    return std::make_shared<const EnergyFluxes>();
}

}  // namespace NiTCAD::assemble
