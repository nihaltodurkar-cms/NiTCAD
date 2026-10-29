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
//   Fermi-Dirac (ED8 as ruled 2026-10-29): d gains Delta ln gamma - (5/2)(g_e - 1) Delta ln theta,
//   gamma Unit 14's degeneracy factor and g = F_{3/2} / F_{1/2}, both at the density
//   n theta^(-3/2) (N_c at T_n), g_e the edge's mean; at theta = 1 Unit 14's flux exactly.
// - Energy flux into a (ED3, c = 5/2, EQ3): S = -kappa grad T - c (k T / q) J with
//   kappa = c (k^2 / q) n mu T is, with v = theta^2 n, c times the same flux in v:
//       S_n = c A (v_b B(d) - v_a B(-d)),   S_p = c A (v_b B(-d) - v_a B(d)),
//   so a uniform theta carries c theta per particle. Fermi-Dirac: c g_e times it, the same d.
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
// - Time steps and small signal (ED9): the energy row of a semiconductor node off the contacts
//   stores -V rate (W - c_W), W = (3/2) c theta g the carriers' thermal energy (energy_storage;
//   g = 1 with Boltzmann), as the continuity rows store their carriers; the small-signal matrix
//   takes -V s dW. The step's history holds W_n and W_p after S_n and S_p (storage width 4).
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
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/contract.hpp"
#include "NiTCAD/physics/energy_balance.hpp"
#include "NiTCAD/physics/fermi_dirac.hpp"
#include "NiTCAD/physics/impact_ionization.hpp"
#include "NiTCAD/physics/recombination.hpp"
#include "NiTCAD/physics/statistics.hpp"
#include "composition.hpp"
#include "scaled_device.hpp"
#include "terms.hpp"

namespace NiTCAD::assemble {

namespace {

constexpr double energy_flux_factor = 2.5;  // c = 5/2 + r, r = 0 (EQ3)

double harmonic_mean(double a, double b) { return 2.0 * a * b / (a + b); }
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

// A hot carrier's degeneracy (ED8 as ruled, "Nc at T_n"): n = N theta^(3/2) F_{1/2}(eta), so
// Unit 14's ln gamma and g = F_{3/2} / F_{1/2} at the density nu = c theta^(-3/2), with their
// partials in c and tau. Boltzmann: 0 and 1.
struct HotDegeneracy {
    double lg, lg_c, lg_t;  // ln gamma
    double g, g_c, g_t;     // F_{3/2} / F_{1/2}
};

HotDegeneracy hot_degeneracy(double c, double tau, double n_ie, double log_dos) {
    const double theta = 1.0 + tau;
    const double s = 1.0 / (theta * std::sqrt(theta));
    const double nu = c * s, nu_t = -1.5 * nu / theta;  // d nu / d tau
    const double N = n_ie * std::exp(log_dos);
    const double v = nu / N;
    double lg, lg_nu, g, g_nu;
    if (v < 1e-6) {
        // As physics::fermi_dirac_degeneracy's series; g = 1 + 2^(-5/2) v + (1/8 - 2 3^(-5/2)) v^2.
        const physics::Degeneracy D = physics::fermi_dirac_degeneracy(n_ie, log_dos, nu);
        lg = D.log_gamma;
        lg_nu = D.d_density;
        constexpr double k1 = 0.17677669529663688110;    // 2^(-5/2)
        constexpr double k2 = -0.0033000598199168200;   // 1/8 - 2 3^(-5/2)
        g = 1.0 + v * (k1 + k2 * v);
        g_nu = (k1 + 2.0 * k2 * v) / N;
    } else {
        const double x = physics::inverse_fermi_half(v);
        const physics::LogDegeneracy L1 = physics::log_degeneracy(x);
        const physics::LogDegeneracy L3 = physics::log_degeneracy_three_halves(x);
        const double x_nu = 1.0 / (nu * (1.0 + L1.d_eta));  // d ln nu / d x = F_{-1/2} / F_{1/2}
        lg = L1.value;
        lg_nu = L1.d_eta * x_nu;
        g = std::exp(L3.value - L1.value);
        g_nu = g * (L3.d_eta - L1.d_eta) * x_nu;
    }
    return {lg, lg_nu * s, lg_nu * nu_t, g, g_nu * s, g_nu * nu_t};
}

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
        // ED11 as amended (ruled 2026-10-29): an interface's charge, traps and surface
        // recombination act on its edges' carrier fluxes, which take no carrier temperature.
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
    // The edge factor A = (an / mu_edge) mu(T_e), dA / dtau of either end (the mean's half) and
    // dA / dE_perp. Lombardi (Unit 46): the ends' low-field mobilities first take the surface
    // terms in E_perp, and an edge in one material takes mu(T) of their harmonic mean, as
    // edge_fluxes' mu0.
    struct Factor {
        double A, d_t, d_perp;
    };
    const SurfaceField perp = e.surface ? surface_field(e, x) : SurfaceField{};
    const auto factor = [&](physics::Carrier carrier, double a0, double mu_edge, double mu_a,
                            double mu_b, double rise, double relax_a, double relax_b,
                            const physics::CanaliParameters& sat_a,
                            const physics::CanaliParameters& sat_b) {
        const double scale = a0 / mu_edge;
        double s_a = 0.0, s_b = 0.0;  // d mu_a / d E_perp, d mu_b / d E_perp
        if (e.surface) {
            const physics::SurfaceMobility la =
                surface_end(carrier, mu_a, e.N_a, surface_T_, perp.value);
            const physics::SurfaceMobility lb =
                surface_end(carrier, mu_b, e.N_b, surface_T_, perp.value);
            mu_a = la.mobility;
            mu_b = lb.mobility;
            s_a = la.d_dE;
            s_b = lb.d_dE;
        }
        if (!e.mixed) {
            double mu0 = mu_edge, d_mu0 = 0.0;
            if (e.surface) {
                const double low = mu_a + mu_b;
                mu0 = harmonic_mean(mu_a, mu_b);
                d_mu0 = 2.0 * (mu_b * mu_b * s_a + mu_a * mu_a * s_b) / (low * low);
            }
            const auto m = physics::carrier_temperature_mobility(mu0, rise, relax_a, sat_a);
            return Factor{scale * m.mobility, scale * m.d_rise * 0.5 * TL,
                          scale * m.d_mu0 * d_mu0};
        }
        const auto ma = physics::carrier_temperature_mobility(mu_a, rise, relax_a, sat_a);
        const auto mb = physics::carrier_temperature_mobility(mu_b, rise, relax_b, sat_b);
        const double sum = ma.mobility + mb.mobility;
        const double h = 2.0 * ma.mobility * mb.mobility / sum;
        // d h / d mu_a(T) = 2 mu_b(T)^2 / sum^2, and b likewise.
        const double ha = 2.0 * mb.mobility * mb.mobility / (sum * sum);
        const double hb = 2.0 * ma.mobility * ma.mobility / (sum * sum);
        return Factor{scale * h, scale * (ha * ma.d_rise + hb * mb.d_rise) * 0.5 * TL,
                      scale * (ha * ma.d_mu0 * s_a + hb * mb.d_mu0 * s_b)};
    };
    // One carrier: sign +1 electrons, -1 holes (the hole flux is the electron form at -d,
    // negated). c_a, c_b its densities, t_a, t_b its rises, delta its driving term, h_a, h_b its
    // ends' degeneracy (Fermi-Dirac; Boltzmann leaves d and the energy factor as slice 2's).
    const auto carrier = [&](double sign, double c_a, double c_b, double t_a, double t_b,
                             double delta, const Factor& F, const HotDegeneracy& h_a,
                             const HotDegeneracy& h_b) {
        EnergyFlux r{};
        const double ta = 1.0 + t_a, tb = 1.0 + t_b;
        const double th = 1.0 + 0.5 * (t_a + t_b);
        const double d0 = sign * delta / th;
        // ED8 as ruled: d = d0 + Delta ln gamma - (5/2)(g_e - 1) Delta ln theta, the energy flux
        // c g_e times the flux in v. L and its partials in each end's density and rise.
        double L = 0.0, G = 1.0, L_c[2] = {0.0, 0.0}, L_t[2] = {0.0, 0.0};
        double G_c[2] = {0.0, 0.0}, G_t[2] = {0.0, 0.0};
        if (fermi_dirac_) {
            const double ge = 0.5 * (h_a.g + h_b.g), dlt = std::log(tb / ta);
            G = ge;
            L = h_b.lg - h_a.lg - 2.5 * (ge - 1.0) * dlt;
            L_c[0] = -h_a.lg_c - 1.25 * h_a.g_c * dlt;
            L_c[1] = h_b.lg_c - 1.25 * h_b.g_c * dlt;
            L_t[0] = -h_a.lg_t - 1.25 * h_a.g_t * dlt + 2.5 * (ge - 1.0) / ta;
            L_t[1] = h_b.lg_t - 1.25 * h_b.g_t * dlt - 2.5 * (ge - 1.0) / tb;
            G_c[0] = 0.5 * h_a.g_c;
            G_c[1] = 0.5 * h_b.g_c;
            G_t[0] = 0.5 * h_a.g_t;
            G_t[1] = 0.5 * h_b.g_t;
        }
        const double d = d0 + L;
        const double bp = bernoulli(d), bm = bernoulli(-d);
        const double dp = bernoulli_derivative(d), dm = bernoulli_derivative(-d);
        const double wa = c_a * ta, wb = c_b * tb, va = wa * ta, vb = wb * tb;
        // In the electron form: P = w_b B(d) - w_a B(-d), Q likewise in v; their d slopes.
        const double P = wb * bp - wa * bm, Q = vb * bp - va * bm;
        const double Pd = wb * dp + wa * dm, Qd = vb * dp + va * dm;
        const double A = F.A, c = energy_flux_factor;
        const double dd_psi = sign / th;  // d d / d psi_b (minus at a)
        const double dd_t[2] = {-0.5 * d0 / th + L_t[0], -0.5 * d0 / th + L_t[1]};
        // Particle flux sign A P; energy flux into a c G A Q (both carriers).
        r.flux = sign * A * P;
        r.d_psi[1] = sign * A * Pd * dd_psi;
        r.d_psi[0] = -r.d_psi[1];
        r.d_c[0] = sign * A * (-ta * bm + Pd * L_c[0]);
        r.d_c[1] = sign * A * (tb * bp + Pd * L_c[1]);
        r.d_t[0] = sign * (F.d_t * P + A * (-c_a * bm) + A * Pd * dd_t[0]);
        r.d_t[1] = sign * (F.d_t * P + A * (c_b * bp) + A * Pd * dd_t[1]);
        r.energy = c * G * A * Q;
        r.e_psi[1] = c * G * A * Qd * dd_psi;
        r.e_psi[0] = -r.e_psi[1];
        r.e_c[0] = c * (G_c[0] * A * Q + G * A * (-ta * ta * bm + Qd * L_c[0]));
        r.e_c[1] = c * (G_c[1] * A * Q + G * A * (tb * tb * bp + Qd * L_c[1]));
        r.e_t[0] = c * (G_t[0] * A * Q +
                        G * (F.d_t * Q + A * (-2.0 * c_a * ta * bm) + A * Qd * dd_t[0]));
        r.e_t[1] = c * (G_t[1] * A * Q +
                        G * (F.d_t * Q + A * (2.0 * c_b * tb * bp) + A * Qd * dd_t[1]));
        r.delta = delta;
        r.d_perp = sign * F.d_perp * P;
        r.e_perp = c * G * F.d_perp * Q;
        return r;
    };
    const Factor Fn = factor(physics::Carrier::electron, e.an, e.mu_n, e.mu_n_a, e.mu_n_b,
                             0.5 * (tna + tnb) * TL, relax_n_[a], relax_n_[b], e.sat_n, e.sat_n_b);
    const Factor Fp = factor(physics::Carrier::hole, e.ap, e.mu_p, e.mu_p_a, e.mu_p_b,
                             0.5 * (tpa + tpb) * TL, relax_p_[a], relax_p_[b], e.sat_p, e.sat_p_b);
    constexpr HotDegeneracy boltzmann{0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
    HotDegeneracy hna = boltzmann, hnb = boltzmann, hpa = boltzmann, hpb = boltzmann;
    if (fermi_dirac_) {
        hna = hot_degeneracy(na, tna, n_ie_[a], log_dos_n_[a]);
        hnb = hot_degeneracy(nb, tnb, n_ie_[b], log_dos_n_[b]);
        hpa = hot_degeneracy(pa, tpa, n_ie_[a], log_dos_p_[a]);
        hpb = hot_degeneracy(pb, tpb, n_ie_[b], log_dos_p_[b]);
    }
    return {carrier(1.0, na, nb, tna, tnb, dpsi + e.shift_n, Fn, hna, hnb),
            carrier(-1.0, pa, pb, tpa, tpb, dpsi + e.shift_p, Fp, hpa, hpb)};
}

DriftDiffusion::EnergyStorage DriftDiffusion::energy_storage(std::size_t i,
                                                             std::span<const double> x) const {
    // W = (3/2) c theta g per carrier, the pressure's (3/2) (ED8).
    const auto one = [&](double c, double tau, double log_dos, double& w, double& w_c,
                         double& w_t) {
        const double theta = 1.0 + tau;
        HotDegeneracy h{0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
        if (fermi_dirac_) h = hot_degeneracy(c, tau, n_ie_[i], log_dos);
        w = 1.5 * c * theta * h.g;
        w_c = 1.5 * theta * (h.g + c * h.g_c);
        w_t = 1.5 * c * (h.g + theta * h.g_t);
    };
    EnergyStorage s{};
    one(x[5 * i + 1], x[5 * i + 3], log_dos_n_[i], s.n, s.n_c, s.n_t);
    one(x[5 * i + 2], x[5 * i + 4], log_dos_p_[i], s.p, s.p_c, s.p_t);
    return s;
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
// After a time step the continuity rows' terms are V k U plus the storage V rate (S - c): the
// storage's band-edge energy is carrier_storage, the energy rows' own storage thermal_storage.
DriftDiffusion::EnergyPowers DriftDiffusion::energy_powers(std::span<const double> x,
                                                           const TimeStep* step) const {
    NITCAD_EXPECTS(energy_balance_ && held_.empty() && x.size() == unknowns());
    const std::size_t nodes = node_count();
    std::vector<double> f(unknowns());
    if (step != nullptr) {
        residual(x, *step, f);
    } else {
        residual(x, f);
    }
    const std::vector<NodeDegeneracy> g = step != nullptr ? degeneracies(x)
                                                          : std::vector<NodeDegeneracy>{};
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
        double ku = dn[i] - f[5 * i + 1];  // V k U, the electron row's term
        if (step != nullptr) {
            const double ks = volume_[i] * step->rate;
            const NodeStorage st = node_storage(i, x, g);
            const EnergyStorage w = energy_storage(i, x);
            const double sn = ks * (st.n - step->storage[4 * i]);
            const double sp = ks * (st.p - step->storage[4 * i + 1]);
            ku -= sn;
            p.carrier_storage += bands.conduction[i] * sn - bands.valence[i] * sp;
            p.thermal_storage +=
                ks * (w.n - step->storage[4 * i + 2] + w.p - step->storage[4 * i + 3]);
        }
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
            if (jacobian) {
                at(s, ev, i, 3, 1) = -k * (rn * rel.n_n * tn + 1.5 * (1.0 + tn) * u_n);
                at(s, ev, i, 3, 2) = -k * (rn * rel.n_p * tn + 1.5 * (1.0 + tn) * u_p);
                at(s, ev, i, 3, 3) = -k * (rn * rel.n + 1.5 * u);
                at(s, ev, i, 4, 1) = -k * (rp * rel.p_n * tp + 1.5 * (1.0 + tp) * u_n);
                at(s, ev, i, 4, 2) = -k * (rp * rel.p_p * tp + 1.5 * (1.0 + tp) * u_p);
                at(s, ev, i, 4, 4) = -k * (rp * rel.p + 1.5 * u);
            }
            if (ev.step == nullptr) continue;
            // ED9: a time step stores -V rate (W - c_W), as the continuity rows' storage.
            const EnergyStorage w = s.energy_storage(i, x);
            const double ks = s.volume_[i] * ev.step->rate;
            ev.f[5 * i + 3] -= ks * (w.n - ev.step->storage[4 * i + 2]);
            ev.f[5 * i + 4] -= ks * (w.p - ev.step->storage[4 * i + 3]);
            if (!jacobian) continue;
            at(s, ev, i, 3, 1) -= ks * w.n_c;
            at(s, ev, i, 3, 3) -= ks * w.n_t;
            at(s, ev, i, 4, 2) -= ks * w.p_c;
            at(s, ev, i, 4, 4) -= ks * w.p_t;
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
    // Lombardi (Unit 46): the four transport rows of both ends with the psi columns of the
    // edge's E_perp stencil; none without the model.
    void pattern(const DriftDiffusion& s, std::vector<linalg::Triplet>& t) const override {
        surface_pattern(s, t, 4);
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
            // Lombardi: d / d psi_j = d / d E_perp * d E_perp / d psi_j over the stencil.
            const SurfaceField perp =
                jacobian && e.surface ? s.surface_field(e, ev.x) : SurfaceField{};
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
                    if (jacobian && e.surface) {
                        // The continuity row; the energy row too when released. Joule is
                        // linear in the edge factor as the fluxes are.
                        const double ep = sign * r.e_perp - 0.5 * r.delta * r.d_perp;
                        for (std::size_t j = e.surface_first; j < e.surface_last; ++j) {
                            const SurfaceStencil& st = s.surface_stencil_[j];
                            const double w = s.surface_gradient(perp, st);
                            values[st.pos[q][c.row - 1]] += sign * r.d_perp * w;
                            if (!held) values[st.pos[q][c.erow - 1]] += ep * w;
                        }
                    }
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

// Impact ionization driven by the carrier temperatures (ED6 as built): at each node that ionizes,
// per carrier G_c = alpha_c(E_eff(T_c)) |j_c|^2 / m_c / q, Unit 19's current reconstruction and
// resolution (ImpactGeneration in terms.cpp) with E_eff(T_c) of the node's own low-field mobility,
// relaxation time and Canali parameters in place of the field along the current. +V G in the
// electron row, -V G in the hole row and +(3/2) V theta_c G in each energy row (ED7: the
// generated carriers take the mean energy, II as negative U). In a uniform field E_eff = E, so it
// is the field-driven generation (EC3). While the temperatures are held, the continuity rows only.
class DriftDiffusion::Terms::EnergyImpact final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"impact ionization (energy balance)", TermDomain::global, transport_rows, all5,
                false};
    }
    // The full 5 x 5 blocks of every carrier edge (energy balance) hold its entries.
    void add(const DriftDiffusion& s, const Evaluation& ev) const override {
        if (s.impact_nodes_.empty()) return;
        const std::span<const double> x = ev.x;
        const std::span<double> f = ev.f, values = ev.values;
        const bool jacobian = ev.jacobian();
        const bool held = !s.held_.empty();
        const int D = s.dimension_;
        const double TL = s.temperature_K_;
        std::vector<std::pair<EnergyFlux, EnergyFlux>> flux(s.edges_.size());
        for (std::size_t k = 0; k < s.edges_.size(); ++k) {
            if (s.edges_[k].carriers) flux[k] = s.energy_fluxes(k, x);
        }
        constexpr double q = base::q_C;
        constexpr double floor2 = impact_current_floor * impact_current_floor;
        const double rho2 = s.impact_resolution_ * s.impact_resolution_;
        for (const ImpactNode& node : s.impact_nodes_) {
            const std::size_t i = node.node;
            const double count = static_cast<double>(node.last - node.first);
            // The current densities at the node and their resolutions R, as Unit 19's.
            double jn[3] = {}, jp[3] = {}, R[2] = {0.0, 0.0};
            for (std::size_t k = node.first; k < node.last; ++k) {
                const ImpactEdge& ie = s.impact_edges_[k];
                const EdgeTerm& e = s.edges_[ie.edge];
                const double delta = x[5 * e.b] - x[5 * e.a];
                const double root = std::sqrt(1.0 + delta * delta);
                R[0] += ie.current * e.an * (x[5 * e.a + 1] + x[5 * e.b + 1]) * root / count;
                R[1] += ie.current * e.ap * (x[5 * e.a + 2] + x[5 * e.b + 2]) * root / count;
                for (int r = 0; r < D; ++r) {
                    double w = 0.0;  // (P t_k)_r
                    for (int c = 0; c < D; ++c) w += node.inverse[r * D + c] * ie.t[c];
                    jn[r] += w * ie.current * flux[ie.edge].first.flux;
                    jp[r] += w * ie.current * flux[ie.edge].second.flux;
                }
            }
            // Per carrier, w = |j|^2 / m, m = sqrt(|j|^2 + eps^2):
            //     dG/dj   = alpha j (|j|^2 + 2 eps^2) / m^3 / q,
            //     dG/dR   = -alpha |j|^2 / m^3 / q rho^2 R,
            //     dG/dtau = alpha'(E_eff) dE_eff/dtau w / q.
            double G = 0.0, dj[2][3] = {}, dR[2] = {0.0, 0.0}, dt[2] = {0.0, 0.0};
            for (int c = 0; c < 2; ++c) {
                const double* j = c == 0 ? jn : jp;
                double mag2 = 0.0;
                for (int r = 0; r < D; ++r) mag2 += j[r] * j[r];
                const double eps2 = floor2 + rho2 * R[c] * R[c];
                const double m = std::sqrt(mag2 + eps2);
                const physics::CarrierTemperatureMobility hot =
                    physics::carrier_temperature_mobility(
                        c == 0 ? node.mu_n : node.mu_p, x[5 * i + 3 + c] * TL,
                        c == 0 ? s.relax_n_[i] : s.relax_p_[i],
                        c == 0 ? node.sat_n : node.sat_p);
                const physics::ImpactIonizationRate a = physics::impact_ionization_coefficient(
                    c == 0 ? node.n : node.p, node.gamma, hot.field);
                if (a.alpha == 0.0 && a.d_dE == 0.0) continue;
                const double w = mag2 / m;
                G += a.alpha * w / q;
                dt[c] = a.d_dE * hot.field_d_rise * TL * w / q;
                const double m3 = m * m * m;
                for (int r = 0; r < D; ++r) {
                    dj[c][r] = a.alpha * j[r] * (mag2 + 2.0 * eps2) / m3 / q;
                }
                dR[c] = -a.alpha * mag2 / m3 / q * rho2 * R[c];
            }
            const double k_row = s.volume_[i] * s.generation_scale_;
            const double theta_n = 1.0 + x[5 * i + 3], theta_p = 1.0 + x[5 * i + 4];
            f[5 * i + 1] += k_row * G;
            f[5 * i + 2] -= k_row * G;
            if (!held) {
                f[5 * i + 3] += 1.5 * k_row * theta_n * G;
                f[5 * i + 4] += 1.5 * k_row * theta_p * G;
            }
            if (!jacobian) continue;
            // Each row's factor on dG: electron, hole and (released) the two energy rows.
            const double factor[4] = {k_row, -k_row, 1.5 * k_row * theta_n, 1.5 * k_row * theta_p};
            const std::size_t rows = held ? 2 : 4;
            const auto add_entry = [&](std::size_t node_j, std::size_t col, std::size_t edge,
                                       double v) {
                for (std::size_t r = 0; r < rows; ++r) {
                    if (node_j == i) {
                        at(s, ev, i, 1 + r, col) += factor[r] * v;
                        continue;
                    }
                    const EnergyEdge& pos = s.energy_edges_[edge];
                    values[(s.edges_[edge].a == i ? pos.ab : pos.ba)[5 * (1 + r) + col]] +=
                        factor[r] * v;
                }
            };
            add_entry(i, 3, 0, dt[0]);
            add_entry(i, 4, 0, dt[1]);
            if (!held) {
                at(s, ev, i, 3, 3) += 1.5 * k_row * G;
                at(s, ev, i, 4, 4) += 1.5 * k_row * G;
            }
            for (std::size_t k = node.first; k < node.last; ++k) {
                const ImpactEdge& ie = s.impact_edges_[k];
                const EdgeTerm& e = s.edges_[ie.edge];
                const EnergyFlux& fn = flux[ie.edge].first;
                const EnergyFlux& fp = flux[ie.edge].second;
                double wn = 0.0, wp = 0.0;  // dG / d(edge current component)
                for (int r = 0; r < D; ++r) {
                    double w = 0.0;
                    for (int c = 0; c < D; ++c) w += node.inverse[r * D + c] * ie.t[c];
                    wn += w * dj[0][r];
                    wp += w * dj[1][r];
                }
                wn *= ie.current;
                wp *= ie.current;
                // R's partials through this edge, as Unit 19's.
                const double delta = x[5 * e.b] - x[5 * e.a];
                const double root = std::sqrt(1.0 + delta * delta);
                const double rn = ie.current * e.an / count, rp = ie.current * e.ap / count;
                const double sn = x[5 * e.a + 1] + x[5 * e.b + 1];
                const double sp = x[5 * e.a + 2] + x[5 * e.b + 2];
                const double dpsi = (dR[0] * rn * sn + dR[1] * rp * sp) * delta / root;
                for (std::size_t o = 0; o < 2; ++o) {
                    const std::size_t end = o == 0 ? e.a : e.b;
                    add_entry(end, 0, ie.edge,
                              wn * fn.d_psi[o] + wp * fp.d_psi[o] + (o == 0 ? -dpsi : dpsi));
                    add_entry(end, 1, ie.edge, wn * fn.d_c[o] + dR[0] * rn * root);
                    add_entry(end, 2, ie.edge, wp * fp.d_c[o] + dR[1] * rp * root);
                    add_entry(end, 3, ie.edge, wn * fn.d_t[o]);
                    add_entry(end, 4, ie.edge, wp * fp.d_t[o]);
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
std::shared_ptr<const TermAccess::Term> DriftDiffusion::Terms::energy_impact() {
    return std::make_shared<const EnergyImpact>();
}

}  // namespace NiTCAD::assemble
