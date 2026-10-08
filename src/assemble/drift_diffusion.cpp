#include "NiTCAD/assemble/drift_diffusion.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <string>

#include "NiTCAD/assemble/contact_bias.hpp"
#include "NiTCAD/assemble/ohmic.hpp"
#include "NiTCAD/assemble/sg_flux.hpp"
#include "NiTCAD/assemble/thermionic_flux.hpp"
#include "NiTCAD/base/contract.hpp"
#include "NiTCAD/physics/field_mobility.hpp"
#include "NiTCAD/physics/interface_traps.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/recombination.hpp"
#include "NiTCAD/physics/thermionic_emission.hpp"
#include "scaled_device.hpp"

namespace NiTCAD::assemble {

namespace {

// The (row, column) offsets inside a node's three unknowns that an edge couples, in EdgeTerm order:
// Poisson to psi, electron continuity to psi and n, hole continuity to psi and p.
constexpr std::size_t edge_rows[5] = {0, 1, 1, 2, 2};
constexpr std::size_t edge_cols[5] = {0, 0, 1, 0, 2};

double harmonic_mean(double a, double b) { return 2.0 * a * b / (a + b); }

}  // namespace

std::expected<DriftDiffusion, base::Error> DriftDiffusion::create(
    const device::Device& device, const Scaling& scaling, const PhysicsModels& models) {
    auto scaled = detail::make_scaled_device(device, scaling, models);
    if (!scaled) return std::unexpected(std::move(scaled.error()));
    const std::size_t n = scaled->volume.size();
    const double T = scaling.temperature_K;

    DriftDiffusion s;
    s.volume_ = std::move(scaled->volume);
    s.doping_ = std::move(scaled->doping);
    s.n_ie_ = std::move(scaled->n_ie);
    s.log_dos_n_ = std::move(scaled->log_dos_n);
    s.log_dos_p_ = std::move(scaled->log_dos_p);
    s.band_shift_ = std::move(scaled->band_shift);
    s.donors_ = std::move(scaled->donors);
    s.acceptors_ = std::move(scaled->acceptors);
    s.levels_ = std::move(scaled->levels);
    s.radiative_ = std::move(scaled->radiative);
    s.fermi_dirac_ = models.fermi_dirac;
    s.ionization_ = models.incomplete_ionization;
    s.radiative_on_ = models.radiative;
    s.contact_ = std::move(scaled->contact);
    s.insulator_ = std::move(scaled->insulator);
    s.electrode_ = std::move(scaled->electrode);
    s.electrode_potential_ = std::move(scaled->electrode_potential);
    s.interfaces_ = std::move(scaled->interfaces);
    s.gates_ = std::move(scaled->gates);
    s.kinds_ = std::move(scaled->kinds);
    s.contact_count_ = device.contacts().size();
    s.rate_scale_ = scaling.Ns / scaling.R0;
    s.Ns_ = scaling.Ns;
    s.auger_ = models.auger;
    s.V_T_ = scaling.V_T;
    s.srh_ = models.srh;
    s.field_mobility_ = models.field_mobility;
    s.tau_n_.resize(n);
    s.tau_p_.resize(n);
    s.auger_n_.resize(n);
    s.auger_p_.resize(n);
    std::vector<double> mu_n(n), mu_p(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto node = static_cast<mesh::NodeId>(i);
        if (s.insulator_[i] != 0) {  // placeholders; no carrier equation reads them
            mu_n[i] = mu_p[i] = 1.0;
            s.tau_n_[i] = s.tau_p_[i] = 1.0;
            s.auger_n_[i] = s.auger_p_[i] = 0.0;
            continue;
        }
        const physics::Semiconductor& m = device.material(node);
        const double N = device.total_impurity(node);
        const double N_mobility = models.doping_mobility ? N : 0.0;
        mu_n[i] = physics::caughey_thomas_mobility(m, physics::Carrier::electron, N_mobility, T);
        mu_p[i] = physics::caughey_thomas_mobility(m, physics::Carrier::hole, N_mobility, T);
        s.tau_n_[i] = physics::scharfetter_lifetime(m, physics::Carrier::electron, N);
        s.tau_p_[i] = physics::scharfetter_lifetime(m, physics::Carrier::hole, N);
        s.auger_n_[i] = m.parameters().auger.Cn;
        s.auger_p_[i] = m.parameters().auger.Cp;
    }

    std::vector<linalg::Triplet> triplets;
    triplets.reserve(9 * n + 10 * scaled->edges.size());
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t r = 0; r < 3; ++r) {
            for (std::size_t c = 0; c < 3; ++c) {
                triplets.push_back({static_cast<linalg::Index>(3 * i + r),
                                    static_cast<linalg::Index>(3 * i + c), 0.0});
            }
        }
    }
    for (const detail::ScaledEdge& e : scaled->edges) {
        for (std::size_t k = 0; k < 5; ++k) {
            triplets.push_back({static_cast<linalg::Index>(3 * e.i + edge_rows[k]),
                                static_cast<linalg::Index>(3 * e.j + edge_cols[k]), 0.0});
            triplets.push_back({static_cast<linalg::Index>(3 * e.j + edge_rows[k]),
                                static_cast<linalg::Index>(3 * e.i + edge_cols[k]), 0.0});
        }
    }
    // An interface edge: the insulator's Poisson row reads the semiconductor's n and p.
    for (const InterfaceEdge& f : s.interfaces_.edges()) {
        for (const std::size_t c : {std::size_t{1}, std::size_t{2}}) {
            triplets.push_back({static_cast<linalg::Index>(3 * f.insulator),
                                static_cast<linalg::Index>(3 * f.semiconductor + c), 0.0});
        }
    }
    const auto size = static_cast<linalg::Index>(3 * n);
    auto pattern = linalg::SparseMatrix::from_triplets(size, size, triplets);
    if (!pattern) return std::unexpected(std::move(pattern.error()));
    s.pattern_ = std::move(*pattern);

    s.block_.resize(9 * n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t r = 0; r < 3; ++r) {
            for (std::size_t c = 0; c < 3; ++c) {
                s.block_[9 * i + 3 * r + c] = detail::position(s.pattern_, 3 * i + r, 3 * i + c);
            }
        }
    }
    s.edges_.reserve(scaled->edges.size());
    for (const detail::ScaledEdge& e : scaled->edges) {
        EdgeTerm t{};
        t.a = e.i;
        t.b = e.j;
        t.c = e.et * e.geometry;
        t.carriers = e.carriers;
        t.charged = e.charged;
        for (std::size_t k = 0; k < 5; ++k) {
            t.ab[k] = detail::position(s.pattern_, 3 * e.i + edge_rows[k], 3 * e.j + edge_cols[k]);
            t.ba[k] = detail::position(s.pattern_, 3 * e.j + edge_rows[k], 3 * e.i + edge_cols[k]);
        }
        if (!t.carriers) {  // Poisson only
            s.edges_.push_back(t);
            continue;
        }
        t.mu_n = harmonic_mean(mu_n[e.i], mu_n[e.j]);
        t.mu_p = harmonic_mean(mu_p[e.i], mu_p[e.j]);
        t.an = t.mu_n * scaling.V_T / scaling.D0 * e.geometry;
        t.ap = t.mu_p * scaling.V_T / scaling.D0 * e.geometry;
        t.field = scaling.V_T / e.length_cm;
        const physics::Semiconductor& m = device.material(static_cast<mesh::NodeId>(e.i));
        t.sat_n = physics::saturation(m, physics::Carrier::electron);
        t.sat_p = physics::saturation(m, physics::Carrier::hole);
        t.mixed = e.interface;
        t.mu_n_a = mu_n[e.i];
        t.mu_n_b = mu_n[e.j];
        t.mu_p_a = mu_p[e.i];
        t.mu_p_b = mu_p[e.j];
        const physics::Semiconductor& mb = device.material(static_cast<mesh::NodeId>(e.j));
        t.sat_n_b = physics::saturation(mb, physics::Carrier::electron);
        t.sat_p_b = physics::saturation(mb, physics::Carrier::hole);
        const double dln = std::log(s.n_ie_[e.j] / s.n_ie_[e.i]);
        const double ds = s.band_shift_[e.j] - s.band_shift_[e.i];
        t.shift_n = ds + dln;
        t.shift_p = ds - dln;
        t.thermionic = e.thermionic;
        if (t.thermionic) {
            // K = hmean(v_a, v_b) times the scaled interface area: coupling area / L_D^(D-1),
            // which is geometry length / L_D; in units of D0 / L_D (the scaled flux of v n).
            const physics::Semiconductor& ma = m;
            const double nca = physics::conduction_band_dos(ma, T);
            const double ncb = physics::conduction_band_dos(mb, T);
            const double nva = physics::valence_band_dos(ma, T);
            const double nvb = physics::valence_band_dos(mb, T);
            const double area = e.geometry * e.length_cm / scaling.D0;
            const auto electron = physics::Carrier::electron;
            const auto hole = physics::Carrier::hole;
            t.te_kn = harmonic_mean(physics::emission_velocity_cm_s(ma, electron, T),
                                    physics::emission_velocity_cm_s(mb, electron, T)) *
                      area;
            t.te_kp = harmonic_mean(physics::emission_velocity_cm_s(ma, hole, T),
                                    physics::emission_velocity_cm_s(mb, hole, T)) *
                      area;
            t.te_log_nc = std::log(ncb / nca);
            t.te_ratio_nc = nca / ncb;
            t.te_log_nv = std::log(nvb / nva);
            t.te_ratio_nv = nva / nvb;
        }
        s.edges_.push_back(t);
    }

    for (const InterfaceEdge& f : s.interfaces_.edges()) {
        s.interface_np_.emplace_back(
            detail::position(s.pattern_, 3 * f.insulator, 3 * f.semiconductor + 1),
            detail::position(s.pattern_, 3 * f.insulator, 3 * f.semiconductor + 2));
    }
    s.psi0_.assign(n, 0.0);
    s.n0_.assign(n, 0.0);
    s.p0_.assign(n, 0.0);
    const std::vector<double> zero(s.contact_count_, 0.0);
    if (auto ok = s.set_bias(zero); !ok) return std::unexpected(std::move(ok.error()));
    return s;
}

std::expected<void, base::Error> DriftDiffusion::set_bias(std::span<const double> bias_V) {
    if (auto ok = check_contact_bias(kinds_, bias_V, false); !ok) return ok;
    gates_.set_bias(bias_V);
    for (std::size_t i = 0; i < node_count(); ++i) {
        if (electrode_[i] >= 0) {
            psi0_[i] = bias_V[static_cast<std::size_t>(electrode_[i])] / V_T_ +
                       electrode_potential_[i];
            continue;
        }
        if (contact_[i] < 0) continue;
        const double bias = bias_V[static_cast<std::size_t>(contact_[i])] / V_T_;
        const OhmicValue v = ohmic_contact_value(
            detail::neutral_equilibrium(fermi_dirac_, ionization_, doping_[i], donors_[i],
                                        acceptors_[i], n_ie_[i], log_dos_n_[i], log_dos_p_[i],
                                        levels_[i]),
            bias);
        psi0_[i] = v.psi - band_shift_[i];
        n0_[i] = v.n;
        p0_[i] = v.p;
    }
    return {};
}

std::vector<double> DriftDiffusion::state_from_potential(std::span<const double> psi) const {
    NITCAD_EXPECTS(psi.size() == node_count());
    std::vector<double> x(unknowns());
    for (std::size_t i = 0; i < node_count(); ++i) {
        const double eta = psi[i] + band_shift_[i];
        x[3 * i] = psi[i];
        if (insulator_[i] != 0) continue;  // stamped below
        x[3 * i + 1] = detail::density(fermi_dirac_, n_ie_[i], log_dos_n_[i], eta).density;
        x[3 * i + 2] = detail::density(fermi_dirac_, n_ie_[i], log_dos_p_[i], -eta).density;
    }
    stamp_contacts(x);
    return x;
}

void DriftDiffusion::stamp_contacts(std::span<double> x) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    for (std::size_t i = 0; i < node_count(); ++i) {
        if (insulator_[i] != 0) {
            if (electrode_[i] >= 0) x[3 * i] = psi0_[i];
            x[3 * i + 1] = x[3 * i + 2] = 0.0;
            continue;
        }
        if (contact_[i] < 0) continue;
        x[3 * i] = psi0_[i];
        x[3 * i + 1] = n0_[i];
        x[3 * i + 2] = p0_[i];
    }
}

linalg::SparseMatrix DriftDiffusion::make_jacobian() const { return pattern_; }

std::vector<DriftDiffusion::NodeDegeneracy> DriftDiffusion::degeneracies(
    std::span<const double> x) const {
    std::vector<NodeDegeneracy> g;
    if (!fermi_dirac_) return g;
    g.reserve(node_count());
    for (std::size_t i = 0; i < node_count(); ++i) {
        if (insulator_[i] != 0) {
            g.push_back({{0.0, 0.0}, {0.0, 0.0}});
            continue;
        }
        g.push_back({physics::fermi_dirac_degeneracy(n_ie_[i], log_dos_n_[i], x[3 * i + 1]),
                     physics::fermi_dirac_degeneracy(n_ie_[i], log_dos_p_[i], x[3 * i + 2])});
    }
    return g;
}

std::pair<EdgeFlux, EdgeFlux> DriftDiffusion::edge_fluxes(
    const EdgeTerm& e, std::span<const double> x, std::span<const NodeDegeneracy> g) const {
    // The flux functions take the driving term as psi2 - psi1; its partials are those with
    // respect to psi_a and psi_b, since ln(n_ie) and s do not depend on the state.
    const double dpsi = x[3 * e.b] - x[3 * e.a];
    const double na = x[3 * e.a + 1], nb = x[3 * e.b + 1];
    const double pa = x[3 * e.a + 2], pb = x[3 * e.b + 2];
    double delta_n = dpsi + e.shift_n, delta_p = dpsi + e.shift_p;
    if (!g.empty()) {
        delta_n += g[e.b].n.log_gamma - g[e.a].n.log_gamma;
        delta_p -= g[e.b].p.log_gamma - g[e.a].p.log_gamma;
    }
    // Unit edge factors when the Canali factor follows (never on a thermionic edge).
    EdgeFlux fn = e.thermionic ? thermionic_electron_flux(e.te_kn, delta_n, e.te_log_nc,
                                                          e.te_ratio_nc, na, nb)
                               : sg_electron_flux(field_mobility_ ? 1.0 : e.an, 0.0, delta_n,
                                                  na, nb);
    EdgeFlux fp = e.thermionic ? thermionic_hole_flux(e.te_kp, delta_p, e.te_log_nv,
                                                      e.te_ratio_nv, pa, pb)
                               : sg_hole_flux(field_mobility_ ? 1.0 : e.ap, 0.0, delta_p, pa, pb);
    if (!g.empty()) {
        // d flux / d delta is d_psi2; ln gamma_n,a enters delta_n with -, ln gamma_p,a delta_p
        // with +.
        fn.d_c1 -= fn.d_psi2 * g[e.a].n.d_density;
        fn.d_c2 += fn.d_psi2 * g[e.b].n.d_density;
        fp.d_c1 += fp.d_psi2 * g[e.a].p.d_density;
        fp.d_c2 -= fp.d_psi2 * g[e.b].p.d_density;
    }
    if (!field_mobility_ || e.thermionic) return {fn, fp};
    // Edge factor a(E) = a0 mu(E) / mu0, E = field |dpsi|; the flux is linear in a, so
    // d flux / d dpsi gains (flux / a) da/d dpsi. mu(E) is the Canali mobility of mu0 (Unit 13),
    // or between two materials the harmonic mean of each end's Canali mobility (Unit 15).
    const double E = e.field * std::abs(dpsi);
    const double dE = dpsi > 0.0 ? e.field : dpsi < 0.0 ? -e.field : 0.0;  // dE / d dpsi
    const auto mobility = [&](double mu0, double mu_a, double mu_b,
                              const physics::CanaliParameters& sat_a,
                              const physics::CanaliParameters& sat_b) {
        if (!e.mixed) return physics::canali_mobility(mu0, E, sat_a);
        const physics::FieldMobility ma = physics::canali_mobility(mu_a, E, sat_a);
        const physics::FieldMobility mb = physics::canali_mobility(mu_b, E, sat_b);
        const double sum = ma.mobility + mb.mobility;
        return physics::FieldMobility{
            2.0 * ma.mobility * mb.mobility / sum,
            2.0 * (mb.mobility * mb.mobility * ma.d_dE + ma.mobility * ma.mobility * mb.d_dE) /
                (sum * sum)};
    };
    const auto scaled = [&](EdgeFlux unit, double a0, double mu0, physics::FieldMobility mu) {
        const double a = a0 * (mu.mobility / mu0);
        const double da = a0 * (mu.d_dE / mu0) * dE;
        return EdgeFlux{a * unit.flux, a * unit.d_psi1 - unit.flux * da,
                        a * unit.d_psi2 + unit.flux * da, a * unit.d_c1, a * unit.d_c2};
    };
    return {scaled(fn, e.an, e.mu_n, mobility(e.mu_n, e.mu_n_a, e.mu_n_b, e.sat_n, e.sat_n_b)),
            scaled(fp, e.ap, e.mu_p, mobility(e.mu_p, e.mu_p_a, e.mu_p_b, e.sat_p, e.sat_p_b))};
}

DriftDiffusion::NodeStorage DriftDiffusion::node_storage(
    std::size_t i, std::span<const double> x, std::span<const NodeDegeneracy> g) const {
    const double n = x[3 * i + 1], p = x[3 * i + 2];
    NodeStorage s{n, 1.0, p, 1.0};
    if (!ionization_) return s;
    // S_n = n - N_D+(eta_c), S_p = p - N_A-(eta_v), with eta_c and eta_v as the Poisson row's.
    const physics::DopantLevels& l = levels_[i];
    if (donors_[i] > 0.0) {
        const double ln = fermi_dirac_ ? g[i].n.log_gamma : 0.0;
        const double wn = fermi_dirac_ ? g[i].n.d_density : 0.0;
        const physics::IonizedDensity d = physics::ionized_density(
            donors_[i], std::log(n / n_ie_[i]) - ln - log_dos_n_[i], l.donor_kT,
            l.donor_degeneracy);
        s.n -= d.value;
        s.d_n -= d.d_eta * (1.0 / n - wn);
    }
    if (acceptors_[i] > 0.0) {
        const double lp = fermi_dirac_ ? g[i].p.log_gamma : 0.0;
        const double wp = fermi_dirac_ ? g[i].p.d_density : 0.0;
        const physics::IonizedDensity a = physics::ionized_density(
            acceptors_[i], std::log(p / n_ie_[i]) - lp - log_dos_p_[i], l.acceptor_kT,
            l.acceptor_degeneracy);
        s.p -= a.value;
        s.d_p -= a.d_eta * (1.0 / p - wp);
    }
    return s;
}

std::optional<TrapStep> DriftDiffusion::trap_step(std::size_t k, const TimeStep* step) const {
    if (step == nullptr) return std::nullopt;
    const InterfaceEdge& v = interfaces_.edges()[k];
    return TrapStep{rate_scale_ * Ns_ / step->rate,
                    step->traps.subspan(interfaces_.slot_offset(k), v.last_level - v.first_level)};
}

InterfaceDrift DriftDiffusion::interface_at(std::size_t k, std::span<const double> x,
                                            const TimeStep* step) const {
    const InterfaceEdge& v = interfaces_.edges()[k];
    const std::size_t o = v.insulator, s = v.semiconductor;
    const std::optional<TrapStep> t = trap_step(k, step);
    return interfaces_.drift(v, x[3 * o], x[3 * s], x[3 * s + 1], x[3 * s + 2], statistics(s),
                             t ? &*t : nullptr);
}

void DriftDiffusion::assemble(std::span<const double> x, std::span<double> f,
                              std::span<double> values, const TimeStep* step,
                              bool interfaces) const {
    NITCAD_EXPECTS(x.size() == unknowns() && f.size() == unknowns());
    NITCAD_EXPECTS(step == nullptr || (step->storage.size() == 2 * node_count() &&
                                       step->traps.size() == trap_slots()));
    const bool jacobian = !values.empty();
    const auto at = [&](std::size_t node, std::size_t r, std::size_t c) -> double& {
        return values[block_[9 * node + 3 * r + c]];
    };
    const std::vector<NodeDegeneracy> g = degeneracies(x);
    // The equilibrium product of node i, scaled.
    const auto product = [&](std::size_t i) {
        return fermi_dirac_ ? physics::fermi_dirac_equilibrium_product(n_ie_[i], g[i].n, g[i].p)
                            : physics::boltzmann_equilibrium_product(n_ie_[i]);
    };
    for (std::size_t i = 0; i < node_count(); ++i) {
        const double psi = x[3 * i], n = x[3 * i + 1], p = x[3 * i + 2];
        if (contact_[i] >= 0) {
            f[3 * i] = psi - psi0_[i];
            f[3 * i + 1] = n - n0_[i];
            f[3 * i + 2] = p - p0_[i];
            if (jacobian) at(i, 0, 0) = at(i, 1, 1) = at(i, 2, 2) = 1.0;
            continue;
        }
        if (insulator_[i] != 0) {  // no charge, no carriers; an electrode holds psi
            const bool electrode = electrode_[i] >= 0;
            f[3 * i] = electrode ? psi - psi0_[i] : 0.0;
            f[3 * i + 1] = n;
            f[3 * i + 2] = p;
            if (jacobian) {
                if (electrode) at(i, 0, 0) = 1.0;
                at(i, 1, 1) = at(i, 2, 2) = 1.0;
            }
            continue;
        }
        const double V = volume_[i];
        f[3 * i] = -V * (n - p - doping_[i]);
        f[3 * i + 1] = 0.0;
        f[3 * i + 2] = 0.0;
        if (jacobian) {
            at(i, 0, 1) = -V;
            at(i, 0, 2) = V;
        }
        if (ionization_) {
            // eta_c = ln(n / n_ie) - ln gamma_n - g_n, d eta_c / dn = 1/n - d ln gamma_n / dn;
            // eta_v likewise from p.
            const double ln = fermi_dirac_ ? g[i].n.log_gamma : 0.0;
            const double lp = fermi_dirac_ ? g[i].p.log_gamma : 0.0;
            const double wn = fermi_dirac_ ? g[i].n.d_density : 0.0;
            const double wp = fermi_dirac_ ? g[i].p.d_density : 0.0;
            const double eta_c = donors_[i] > 0.0 ? std::log(n / n_ie_[i]) - ln - log_dos_n_[i]
                                                  : 0.0;
            const double eta_v = acceptors_[i] > 0.0
                                     ? std::log(p / n_ie_[i]) - lp - log_dos_p_[i]
                                     : 0.0;
            const detail::IonizedCharge c =
                detail::ionized_charge(donors_[i], acceptors_[i], eta_c, eta_v, levels_[i]);
            f[3 * i] = -V * (n - p - c.value);
            if (jacobian) {
                at(i, 0, 1) = -V * (1.0 - c.d_eta_c * (1.0 / n - wn));
                at(i, 0, 2) = V * (1.0 + c.d_eta_v * (1.0 / p - wp));
            }
        }
        if (gates_.on_gate(i)) {
            f[3 * i] += gates_.row_term(i, psi);
            if (jacobian) at(i, 0, 0) -= gates_.term(i).coupling;
        }
        if (srh_) {
            const physics::RecombinationRate r =
                physics::srh_recombination(n, p, product(i), n_ie_[i], tau_n_[i], tau_p_[i]);
            const double k = V * rate_scale_;
            f[3 * i + 1] -= k * r.rate;
            f[3 * i + 2] += k * r.rate;
            if (jacobian) {
                at(i, 1, 1) -= k * r.d_dn;
                at(i, 1, 2) -= k * r.d_dp;
                at(i, 2, 1) += k * r.d_dn;
                at(i, 2, 2) += k * r.d_dp;
            }
        }
        // Auger and radiative recombination take physical densities; R / R0 = rate_scale (R / Ns),
        // d(R / R0)/dn' = rate_scale dR/dn. The scaled product E' = E / Ns^2 has
        // dE/dn = Ns dE'/dn'.
        const auto physical_product = [&] {
            const double nie = n_ie_[i] * Ns_;
            physics::EquilibriumProduct E = physics::boltzmann_equilibrium_product(nie);
            if (fermi_dirac_) {
                const physics::EquilibriumProduct s = product(i);
                E = {s.value * Ns_ * Ns_, s.d_dn * Ns_, s.d_dp * Ns_};
            }
            return E;
        };
        const auto add_physical = [&](const physics::RecombinationRate& r) {
            const double k = V * rate_scale_;
            f[3 * i + 1] -= k * (r.rate / Ns_);
            f[3 * i + 2] += k * (r.rate / Ns_);
            if (jacobian) {
                at(i, 1, 1) -= k * r.d_dn;
                at(i, 1, 2) -= k * r.d_dp;
                at(i, 2, 1) += k * r.d_dn;
                at(i, 2, 2) += k * r.d_dp;
            }
        };
        if (auger_) {
            add_physical(physics::auger_recombination(n * Ns_, p * Ns_, physical_product(),
                                                      auger_n_[i], auger_p_[i]));
        }
        if (radiative_on_ && radiative_[i] > 0.0) {
            add_physical(physics::radiative_recombination(n * Ns_, p * Ns_, physical_product(),
                                                          radiative_[i]));
        }
        if (step != nullptr) {
            const NodeStorage st = node_storage(i, x, g);
            const double k = V * step->rate;
            f[3 * i + 1] -= k * (st.n - step->storage[2 * i]);
            f[3 * i + 2] += k * (st.p - step->storage[2 * i + 1]);
            if (jacobian) {
                at(i, 1, 1) -= k * st.d_n;
                at(i, 2, 2) += k * st.d_p;
            }
        }
    }
    for (std::size_t k = 0; interfaces && k < interfaces_.edges().size(); ++k) {
        const InterfaceEdge& v = interfaces_.edges()[k];
        const std::size_t o = v.insulator, s = v.semiconductor;
        const InterfaceDrift t = interface_at(k, x, step);
        // Edge positions: row of o with column psi of s, and the rows of s with column psi of o.
        const EdgeTerm& e = edges_[v.edge];
        const bool o_first = e.a == o;
        const std::size_t os = o_first ? e.ab[0] : e.ba[0];
        const std::size_t so[3] = {o_first ? e.ba[0] : e.ab[0], o_first ? e.ba[1] : e.ab[1],
                                   o_first ? e.ba[3] : e.ab[3]};
        if (electrode_[o] < 0) {  // the insulator node's Poisson row
            f[3 * o] += t.flux_insulator;
            if (jacobian) {
                at(o, 0, 0) += t.d_flux_insulator[0];
                values[os] += t.d_flux_insulator[1];
                values[interface_np_[k].first] += t.d_flux_insulator[2];
                values[interface_np_[k].second] += t.d_flux_insulator[3];
            }
        }
        if (contact_[s] >= 0) continue;  // a Dirichlet semiconductor node
        f[3 * s] += t.flux_semiconductor;
        f[3 * s + 1] -= t.rate;
        f[3 * s + 2] += t.rate_p;
        if (jacobian) {
            values[so[0]] += t.d_flux_semiconductor[0];
            values[so[1]] -= t.d_rate[0];
            values[so[2]] += t.d_rate_p[0];
            for (std::size_t c = 1; c < 4; ++c) {
                at(s, 0, c - 1) += t.d_flux_semiconductor[c];
                at(s, 1, c - 1) -= t.d_rate[c];
                at(s, 2, c - 1) += t.d_rate_p[c];
            }
        }
    }
    for (const EdgeTerm& e : edges_) {
        const std::size_t a = e.a, b = e.b;
        const double psi_a = x[3 * a], psi_b = x[3 * b];
        const double poisson = e.c * (psi_b - psi_a);
        if (e.charged) continue;  // the interface's half-edge fluxes, above
        if (!e.carriers) {  // Poisson only, into the rows that are not Dirichlet
            if (contact_[a] < 0 && electrode_[a] < 0) {
                f[3 * a] += poisson;
                if (jacobian) {
                    at(a, 0, 0) -= e.c;
                    values[e.ab[0]] += e.c;
                }
            }
            if (contact_[b] < 0 && electrode_[b] < 0) {
                f[3 * b] -= poisson;
                if (jacobian) {
                    at(b, 0, 0) -= e.c;
                    values[e.ba[0]] += e.c;
                }
            }
            continue;
        }
        const auto [fn, fp] = edge_fluxes(e, x, g);
        if (contact_[a] < 0) {
            f[3 * a] += poisson;
            f[3 * a + 1] += fn.flux;
            f[3 * a + 2] += fp.flux;
            if (jacobian) {
                at(a, 0, 0) -= e.c;
                values[e.ab[0]] += e.c;
                at(a, 1, 0) += fn.d_psi1;
                values[e.ab[1]] += fn.d_psi2;
                at(a, 1, 1) += fn.d_c1;
                values[e.ab[2]] += fn.d_c2;
                at(a, 2, 0) += fp.d_psi1;
                values[e.ab[3]] += fp.d_psi2;
                at(a, 2, 2) += fp.d_c1;
                values[e.ab[4]] += fp.d_c2;
            }
        }
        if (contact_[b] < 0) {
            f[3 * b] -= poisson;
            f[3 * b + 1] -= fn.flux;
            f[3 * b + 2] -= fp.flux;
            if (jacobian) {
                at(b, 0, 0) -= e.c;
                values[e.ba[0]] += e.c;
                values[e.ba[1]] -= fn.d_psi1;
                at(b, 1, 0) -= fn.d_psi2;
                values[e.ba[2]] -= fn.d_c1;
                at(b, 1, 1) -= fn.d_c2;
                values[e.ba[3]] -= fp.d_psi1;
                at(b, 2, 0) -= fp.d_psi2;
                values[e.ba[4]] -= fp.d_c1;
                at(b, 2, 2) -= fp.d_c2;
            }
        }
    }
}

void DriftDiffusion::evaluate(std::span<const double> x, std::span<double> residual,
                              linalg::SparseMatrix& jacobian) const {
    NITCAD_EXPECTS(jacobian.has_same_pattern(pattern_));
    const std::span<double> values = jacobian.values();
    std::fill(values.begin(), values.end(), 0.0);
    assemble(x, residual, values, nullptr);
}

void DriftDiffusion::residual(std::span<const double> x, std::span<double> residual) const {
    assemble(x, residual, {}, nullptr);
}

void DriftDiffusion::evaluate(std::span<const double> x, const TimeStep& step,
                              std::span<double> residual, linalg::SparseMatrix& jacobian) const {
    NITCAD_EXPECTS(jacobian.has_same_pattern(pattern_));
    const std::span<double> values = jacobian.values();
    std::fill(values.begin(), values.end(), 0.0);
    assemble(x, residual, values, &step);
}

void DriftDiffusion::residual(std::span<const double> x, const TimeStep& step,
                              std::span<double> residual) const {
    assemble(x, residual, {}, &step);
}

std::vector<double> DriftDiffusion::storage(std::span<const double> x) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    const std::vector<NodeDegeneracy> g = degeneracies(x);
    std::vector<double> s(2 * node_count(), 0.0);
    for (std::size_t i = 0; i < node_count(); ++i) {
        if (insulator_[i] != 0) continue;
        const NodeStorage st = node_storage(i, x, g);
        s[2 * i] = st.n;
        s[2 * i + 1] = st.p;
    }
    return s;
}

std::vector<double> DriftDiffusion::trap_occupancies(std::span<const double> x,
                                                     const TimeStep* step) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    NITCAD_EXPECTS(step == nullptr || step->traps.size() == trap_slots());
    std::vector<double> f(trap_slots(), 0.0);
    for (std::size_t k = 0; k < interfaces_.edges().size(); ++k) {
        const InterfaceEdge& v = interfaces_.edges()[k];
        const std::size_t o = v.insulator, s = v.semiconductor;
        const std::optional<TrapStep> t = trap_step(k, step);
        interfaces_.occupancies(
            v, x[3 * o], x[3 * s], x[3 * s + 1], x[3 * s + 2], statistics(s), t ? &*t : nullptr,
            std::span(f).subspan(interfaces_.slot_offset(k), v.last_level - v.first_level));
    }
    return f;
}

std::vector<double> DriftDiffusion::contact_charges(std::span<const double> x,
                                                    const TimeStep* step) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    std::vector<double> charge = gates_.charges(x, 3, contact_count_);
    const auto owner = [&](std::size_t i) {
        return contact_[i] >= 0 ? contact_[i] : electrode_[i];
    };
    for (const EdgeTerm& e : edges_) {
        const std::int32_t ca = owner(e.a), cb = owner(e.b);
        if (ca == cb || e.charged) continue;
        const double flux = e.c * (x[3 * e.a] - x[3 * e.b]);  // from a to b
        if (ca >= 0) charge[static_cast<std::size_t>(ca)] += flux;
        if (cb >= 0) charge[static_cast<std::size_t>(cb)] -= flux;
    }
    // A contact node on an interface edge: the flux leaving it is minus its half-edge flux.
    for (std::size_t k = 0; k < interfaces_.edges().size(); ++k) {
        const InterfaceEdge& v = interfaces_.edges()[k];
        const std::int32_t ci = owner(v.insulator), cs = owner(v.semiconductor);
        if (ci < 0 && cs < 0) continue;
        const InterfaceDrift t = interface_at(k, x, step);
        if (ci >= 0) charge[static_cast<std::size_t>(ci)] -= t.flux_insulator;
        if (cs >= 0) charge[static_cast<std::size_t>(cs)] -= t.flux_semiconductor;
    }
    return charge;
}

std::vector<double> DriftDiffusion::conduction_currents(std::span<const double> x,
                                                        const TimeStep& step) const {
    std::vector<double> terminal = terminal_currents(x);
    // The traps of an interface edge on an ohmic node exchange carriers with the contact:
    // electron capture takes rate from it, hole capture rate_p.
    for (std::size_t k = 0; k < interfaces_.edges().size(); ++k) {
        const std::int32_t c = contact_[interfaces_.edges()[k].semiconductor];
        if (c < 0) continue;
        const InterfaceDrift t = interface_at(k, x, &step);
        terminal[static_cast<std::size_t>(c)] += t.rate_p - t.rate;
    }
    return terminal;
}

double DriftDiffusion::update_size(std::span<const double> x, std::span<const double> dx) const {
    NITCAD_EXPECTS(x.size() == unknowns() && dx.size() == unknowns());
    constexpr double infinity = std::numeric_limits<double>::infinity();
    // A density correction is resolved by the linear solve only to about eps max|dx|, and near
    // convergence max|dx| is itself about eps times the largest density: densities below about
    // 1e-23 of the largest cannot reach a relative 1e-8. They are measured against this floor.
    double largest = 0.0;
    for (std::size_t i = 0; i < node_count(); ++i) {
        largest = std::max({largest, x[3 * i + 1], x[3 * i + 2]});
    }
    const double floor = 1e-20 * largest;
    double size = 0.0;
    for (std::size_t i = 0; i < node_count(); ++i) {
        if (insulator_[i] != 0) {
            size = std::max(size, std::abs(dx[3 * i]));
            continue;
        }
        const double n = x[3 * i + 1], p = x[3 * i + 2];
        size = std::max({size, std::abs(dx[3 * i]),
                         n > 0.0 ? std::abs(dx[3 * i + 1]) / std::max(n, floor) : infinity,
                         p > 0.0 ? std::abs(dx[3 * i + 2]) / std::max(p, floor) : infinity});
    }
    return size;
}

void DriftDiffusion::apply_update(std::span<double> x, std::span<const double> dx,
                                  double max_update) const {
    NITCAD_EXPECTS(x.size() == unknowns() && dx.size() == unknowns());
    for (std::size_t i = 0; i < node_count(); ++i) {
        x[3 * i] += std::clamp(dx[3 * i], -max_update, max_update);
        if (insulator_[i] != 0) continue;  // the densities stay 0
        for (const std::size_t k : {3 * i + 1, 3 * i + 2}) {
            x[k] = std::clamp(x[k] + dx[k], 0.1 * x[k], 10.0 * x[k]);
        }
    }
}

std::vector<std::pair<double, double>> DriftDiffusion::edge_currents(
    std::span<const double> x) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    std::vector<std::pair<double, double>> currents;
    currents.reserve(edges_.size());
    const std::vector<NodeDegeneracy> g = degeneracies(x);
    for (const EdgeTerm& e : edges_) {
        if (!e.carriers) {
            currents.emplace_back(0.0, 0.0);
            continue;
        }
        const auto [fn, fp] = edge_fluxes(e, x, g);
        currents.emplace_back(fn.flux, fp.flux);
    }
    return currents;
}

std::vector<double> DriftDiffusion::terminal_currents(std::span<const double> x) const {
    const auto currents = edge_currents(x);
    std::vector<double> terminal(contact_count_, 0.0);
    for (std::size_t k = 0; k < edges_.size(); ++k) {
        const std::int32_t ca = contact_[edges_[k].a], cb = contact_[edges_[k].b];
        if (ca == cb) continue;  // both inside one contact, or both in the device
        const double total = currents[k].first + currents[k].second;
        if (ca >= 0) terminal[static_cast<std::size_t>(ca)] += total;  // leaves contact a
        if (cb >= 0) terminal[static_cast<std::size_t>(cb)] -= total;  // flows into contact b
    }
    return terminal;
}

std::vector<double> DriftDiffusion::terminal_current_resolution(std::span<const double> x) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    // The current through any cut of the device differs from a terminal current by the total
    // current's imbalance at the nodes between them, which is at most the sum over the nodes of
    // |F_n + F_p| (the residual) plus the rounding of the terms that cancel in it: the one-sided
    // flux terms (density times its flux coefficient) and the recombination terms.
    const std::size_t n = node_count();
    std::vector<double> f(unknowns());
    residual(x, f);
    const std::vector<NodeDegeneracy> g = degeneracies(x);
    std::vector<double> terms(n, 0.0), flux_n(n, 0.0);
    for (const EdgeTerm& e : edges_) {
        if (!e.carriers) continue;
        const auto [fn, fp] = edge_fluxes(e, x, g);
        const double t =
            std::abs(fn.d_c1 * x[3 * e.a + 1]) + std::abs(fn.d_c2 * x[3 * e.b + 1]) +
            std::abs(fp.d_c1 * x[3 * e.a + 2]) + std::abs(fp.d_c2 * x[3 * e.b + 2]);
        terms[e.a] += t;
        terms[e.b] += t;
        flux_n[e.a] += fn.flux;
        flux_n[e.b] -= fn.flux;
    }
    const double eps8 = 8.0 * std::numeric_limits<double>::epsilon();
    double bound = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        if (contact_[i] >= 0) {
            bound += eps8 * terms[i];  // the terminal current's own sum
            continue;
        }
        const double recombination = std::abs(flux_n[i] - f[3 * i + 1]);  // V R
        bound += std::abs(f[3 * i + 1] + f[3 * i + 2]) + eps8 * (terms[i] + 2.0 * recombination);
    }
    std::vector<double> resolution(contact_count_, 0.0);
    for (std::size_t c = 0; c < contact_count_; ++c) {
        if (kinds_[c] == device::ContactKind::ohmic) resolution[c] = bound;
    }
    return resolution;
}

BandEdges DriftDiffusion::band_edges(std::span<const double> x) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    const std::size_t n = node_count();
    BandEdges b{std::vector<double>(n), std::vector<double>(n), std::vector<double>(n),
                std::vector<double>(n)};
    for (std::size_t i = 0; i < n; ++i) {
        if (insulator_[i] != 0) {
            const double nan = std::numeric_limits<double>::quiet_NaN();
            b.conduction[i] = b.valence[i] = b.electron_fermi[i] = b.hole_fermi[i] = nan;
            continue;
        }
        const detail::NodeBands e =
            detail::node_bands(fermi_dirac_, n_ie_[i], log_dos_n_[i], log_dos_p_[i],
                               x[3 * i] + band_shift_[i], x[3 * i + 1], x[3 * i + 2]);
        b.conduction[i] = e.conduction;
        b.valence[i] = e.valence;
        b.electron_fermi[i] = e.electron_fermi;
        b.hole_fermi[i] = e.hole_fermi;
    }
    return b;
}

std::vector<double> DriftDiffusion::gate_charges(std::span<const double> x,
                                                 const TimeStep* step) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    std::vector<double> charge = gates_.charges(x, 3, contact_count_);
    for (const EdgeTerm& e : edges_) {
        const std::int32_t ca = electrode_[e.a], cb = electrode_[e.b];
        if (ca == cb || e.charged) continue;
        const double flux = e.c * (x[3 * e.a] - x[3 * e.b]);  // from a to b
        if (ca >= 0) charge[static_cast<std::size_t>(ca)] += flux;
        if (cb >= 0) charge[static_cast<std::size_t>(cb)] -= flux;
    }
    // An electrode node on an interface edge: the flux leaving it is minus its half-edge flux.
    for (std::size_t k = 0; k < interfaces_.edges().size(); ++k) {
        const std::int32_t c = electrode_[interfaces_.edges()[k].insulator];
        if (c < 0) continue;
        charge[static_cast<std::size_t>(c)] -= interface_at(k, x, step).flux_insulator;
    }
    return charge;
}

std::vector<double> DriftDiffusion::interface_trap_charges(std::span<const double> x,
                                                           const TimeStep* step) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    std::vector<double> q(interfaces_.interface_count(), 0.0);
    for (std::size_t k = 0; k < interfaces_.edges().size(); ++k) {
        q[interfaces_.edges()[k].interface] += interface_at(k, x, step).trapped;
    }
    return q;
}

linalg::ComplexSparseMatrix DriftDiffusion::make_small_signal_matrix() const {
    std::vector<linalg::ComplexTriplet> t;
    t.reserve(pattern_.nonzeros());
    for (linalg::Index r = 0; r < pattern_.rows(); ++r) {
        for (linalg::Index k = pattern_.row_offsets()[r]; k < pattern_.row_offsets()[r + 1]; ++k) {
            t.push_back({r, pattern_.col_indices()[k], {}});
        }
    }
    // The pattern was built from triplets already, so this cannot fail.
    return *linalg::ComplexSparseMatrix::from_triplets(pattern_.rows(), pattern_.cols(), t);
}

void DriftDiffusion::small_signal_matrix(std::span<const double> x, std::complex<double> s,
                                         linalg::ComplexSparseMatrix& a) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    NITCAD_EXPECTS(a.rows() == pattern_.rows() &&
                   std::ranges::equal(a.row_offsets(), pattern_.row_offsets()) &&
                   std::ranges::equal(a.col_indices(), pattern_.col_indices()));
    using Complex = std::complex<double>;
    // The steady rows without the interface edges, then the storage and the interfaces' terms.
    std::vector<double> f(unknowns()), real(pattern_.nonzeros(), 0.0);
    assemble(x, f, real, nullptr, false);
    const std::span<Complex> values = a.values();
    for (std::size_t k = 0; k < real.size(); ++k) values[k] = real[k];
    const auto at = [&](std::size_t node, std::size_t r, std::size_t c) -> Complex& {
        return values[block_[9 * node + 3 * r + c]];
    };
    const std::vector<NodeDegeneracy> g = degeneracies(x);
    for (std::size_t i = 0; i < node_count(); ++i) {
        if (contact_[i] >= 0 || insulator_[i] != 0) continue;
        const NodeStorage st = node_storage(i, x, g);
        at(i, 1, 1) -= s * (volume_[i] * st.d_n);
        at(i, 2, 2) += s * (volume_[i] * st.d_p);
    }
    // The traps' weight Ns / s_physical, s_physical = s / t0.
    const Complex inv_weight = s / (rate_scale_ * Ns_);
    for (std::size_t k = 0; k < interfaces_.edges().size(); ++k) {
        const InterfaceEdge& v = interfaces_.edges()[k];
        const std::size_t o = v.insulator, sn = v.semiconductor;
        const InterfaceSmallSignal t = interfaces_.small_signal(
            v, x[3 * o], x[3 * sn], x[3 * sn + 1], x[3 * sn + 2], statistics(sn), inv_weight);
        // As assemble's interface terms.
        const EdgeTerm& e = edges_[v.edge];
        const bool o_first = e.a == o;
        const std::size_t os = o_first ? e.ab[0] : e.ba[0];
        const std::size_t so[3] = {o_first ? e.ba[0] : e.ab[0], o_first ? e.ba[1] : e.ab[1],
                                   o_first ? e.ba[3] : e.ab[3]};
        if (electrode_[o] < 0) {
            at(o, 0, 0) += t.d_flux_insulator[0];
            values[os] += t.d_flux_insulator[1];
            values[interface_np_[k].first] += t.d_flux_insulator[2];
            values[interface_np_[k].second] += t.d_flux_insulator[3];
        }
        if (contact_[sn] >= 0) continue;
        values[so[0]] += t.d_flux_semiconductor[0];
        values[so[1]] -= t.d_rate[0];
        values[so[2]] += t.d_rate_p[0];
        for (std::size_t c = 1; c < 4; ++c) {
            at(sn, 0, c - 1) += t.d_flux_semiconductor[c];
            at(sn, 1, c - 1) -= t.d_rate[c];
            at(sn, 2, c - 1) += t.d_rate_p[c];
        }
    }
}

std::vector<double> DriftDiffusion::bias_derivative(std::size_t contact) const {
    NITCAD_EXPECTS(contact < contact_count_);
    const auto c = static_cast<std::int32_t>(contact);
    std::vector<double> d(unknowns(), 0.0);
    for (std::size_t i = 0; i < node_count(); ++i) {
        // A Dirichlet row psi - psi0 with psi0 = V / V_T + const; a gate row G (psi_G - psi),
        // psi_G = V / V_T - offset. An ohmic node's densities do not depend on the bias.
        if (contact_[i] == c || electrode_[i] == c) {
            d[3 * i] = -1.0 / V_T_;
        } else if (gates_.on_gate(i) && gates_.contact(i) == contact) {
            d[3 * i] = gates_.term(i).coupling / V_T_;
        }
    }
    return d;
}

std::vector<DriftDiffusion::CurrentRow> DriftDiffusion::small_signal_currents(
    std::span<const double> x, std::complex<double> s) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    using Complex = std::complex<double>;
    std::vector<std::vector<std::pair<std::size_t, Complex>>> terms(contact_count_);
    std::vector<CurrentRow> rows(contact_count_);
    const auto add = [&](std::int32_t c, std::size_t column, Complex value) {
        if (c >= 0) terms[static_cast<std::size_t>(c)].emplace_back(column, value);
    };
    // Conduction (terminal_currents): the total flux on the edges leaving each ohmic contact.
    const std::vector<NodeDegeneracy> g = degeneracies(x);
    for (const EdgeTerm& e : edges_) {
        const std::int32_t ca = contact_[e.a], cb = contact_[e.b];
        if (ca == cb || !e.carriers) continue;
        const auto [fn, fp] = edge_fluxes(e, x, g);
        const std::pair<std::size_t, double> d[6] = {
            {3 * e.a, fn.d_psi1 + fp.d_psi1}, {3 * e.b, fn.d_psi2 + fp.d_psi2},
            {3 * e.a + 1, fn.d_c1},           {3 * e.b + 1, fn.d_c2},
            {3 * e.a + 2, fp.d_c1},           {3 * e.b + 2, fp.d_c2}};
        for (const auto& [column, value] : d) {
            add(ca, column, value);
            add(cb, column, -value);
        }
    }
    // The traps of an interface edge on an ohmic node (conduction_currents): rate_p - rate; and
    // the charge of a contact node on an interface edge: minus its half-edge flux.
    const auto owner = [&](std::size_t i) {
        return contact_[i] >= 0 ? contact_[i] : electrode_[i];
    };
    const Complex inv_weight = s / (rate_scale_ * Ns_);
    for (std::size_t k = 0; k < interfaces_.edges().size(); ++k) {
        const InterfaceEdge& v = interfaces_.edges()[k];
        const std::size_t o = v.insulator, sn = v.semiconductor;
        const std::int32_t ci = owner(o), cs = owner(sn);
        if (ci < 0 && cs < 0) continue;
        const InterfaceSmallSignal t = interfaces_.small_signal(
            v, x[3 * o], x[3 * sn], x[3 * sn + 1], x[3 * sn + 2], statistics(sn), inv_weight);
        const std::size_t columns[4] = {3 * o, 3 * sn, 3 * sn + 1, 3 * sn + 2};
        for (std::size_t c = 0; c < 4; ++c) {
            add(contact_[sn], columns[c], t.d_rate_p[c] - t.d_rate[c]);
            add(ci, columns[c], -s * t.d_flux_insulator[c]);
            add(cs, columns[c], -s * t.d_flux_semiconductor[c]);
        }
    }
    // Displacement: s times the contact charges (contact_charges).
    for (std::size_t i = 0; i < node_count(); ++i) {
        if (!gates_.on_gate(i)) continue;
        const std::size_t c = gates_.contact(i);
        const double G = gates_.term(i).coupling;
        terms[c].emplace_back(3 * i, -s * G);
        rows[c].bias += s * (G / V_T_);
    }
    for (const EdgeTerm& e : edges_) {
        const std::int32_t ca = owner(e.a), cb = owner(e.b);
        if (ca == cb || e.charged) continue;
        add(ca, 3 * e.a, s * e.c);
        add(ca, 3 * e.b, -s * e.c);
        add(cb, 3 * e.a, -s * e.c);
        add(cb, 3 * e.b, s * e.c);
    }
    // Merge each row's terms by column, in increasing column order.
    for (std::size_t c = 0; c < contact_count_; ++c) {
        std::ranges::stable_sort(terms[c], {}, &std::pair<std::size_t, Complex>::first);
        for (const auto& [column, value] : terms[c]) {
            if (!rows[c].columns.empty() && rows[c].columns.back() == column) {
                rows[c].values.back() += value;
            } else {
                rows[c].columns.push_back(column);
                rows[c].values.push_back(value);
            }
        }
    }
    return rows;
}

InterfaceStatistics DriftDiffusion::statistics(std::size_t node) const noexcept {
    return {fermi_dirac_, n_ie_[node], log_dos_n_[node], log_dos_p_[node], band_shift_[node]};
}

}  // namespace NiTCAD::assemble
