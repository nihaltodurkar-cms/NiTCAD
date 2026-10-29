// The electrothermal rows of DriftDiffusion (Unit 23; drift_diffusion.hpp, "Electrothermal
// coupling"; DECISIONS.md T1-T14).
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <variant>

#include "NiTCAD/assemble/bernoulli.hpp"
#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/observer.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/contract.hpp"
#include "NiTCAD/physics/bandgap_narrowing.hpp"
#include "NiTCAD/physics/field_mobility.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/recombination.hpp"
#include "NiTCAD/physics/statistics.hpp"
#include "NiTCAD/physics/thermal.hpp"
#include "composition.hpp"
#include "scaled_device.hpp"
#include "terms.hpp"

namespace NiTCAD::assemble {

namespace {

base::Error refused(std::string message) {
    return {base::ErrorCode::invalid_input, std::move(message), {}};
}

// (theta^(1 - a) - 1) / (1 - a), and ln theta at a = 1, with theta = 1 + tau: the Kirchhoff
// transform of (T / T0)^-a, from T0 to theta T0, in units of T0; accurate to round-off of tau.
double kirchhoff(double tau, double a) {
    const double s = 1.0 - a;
    const double l = std::log1p(tau);
    return std::abs(s) < 1e-12 ? l : std::expm1(s * l) / s;
}

double harmonic_mean(double a, double b) { return 2.0 * a * b / (a + b); }

}  // namespace

std::optional<base::Error> DriftDiffusion::make_thermal(const device::Device& device,
                                                        const Scaling& scaling,
                                                        const PhysicsModels& models,
                                                        std::span<const double> edge_geometry) {
    // DECISIONS.md T8: traps and band-to-band tunnelling are not coupled to the lattice
    // temperature; freezing them at T0 would be a silent approximation.
    if (models.btbt_local || models.btbt_nonlocal != NonlocalTunnelling::off) {
        return refused("electrothermal does not support band-to-band tunnelling (DECISIONS.md T8)");
    }
    for (const device::Interface& f : device.interfaces()) {
        if (!f.traps.levels.empty() || !f.traps.bands.empty()) {
            return refused("electrothermal does not support interface traps (DECISIONS.md T8)");
        }
        // DECISIONS.md U3 (deferred): thermionic emission across a temperature step has no
        // validated non-isothermal form yet.
        if (f.transport == device::InterfaceTransport::thermionic_emission) {
            return refused("electrothermal does not support thermionic emission (DECISIONS.md U3, "
                           "deferred: no validated non-isothermal form)");
        }
    }
    // DECISIONS.md U6: the energy of carriers bound to dopants is not accounted for yet.
    if (models.incomplete_ionization) {
        return refused("electrothermal does not support incomplete ionization (DECISIONS.md U6: "
                       "the bound carriers' energy is not derived)");
    }
    // Every material needs thermal data (T11); none is invented.
    const auto regions = device.regions();
    std::vector<std::int32_t> region_material(regions.size(), -1);
    std::vector<physics::ThermalParameters> region_thermal(regions.size());
    for (std::size_t r = 0; r < regions.size(); ++r) {
        const device::Region& region = regions[r];
        if (device::is_insulator(region)) {
            region_thermal[r] = std::get<physics::Insulator>(region.material).parameters().thermal;
        } else {
            const auto& m = std::get<physics::Semiconductor>(region.material);
            region_thermal[r] = m.parameters().thermal;
            region_material[r] = static_cast<std::int32_t>(thermal_materials_.size());
            thermal_materials_.push_back(m);
            masetti_holes_.push_back(masetti_holes(models.hole_mobility, m.parameters()) ? 1 : 0);
        }
        if (!physics::has_thermal_data(region_thermal[r])) {
            return refused("electrothermal: region '" + region.name +
                           "' has no thermal data (a conductivity and heat capacity)");
        }
        // T3: the Fermi-Dirac thermopower is built for r = -1/2 (F_1 / F_0) alone.
        const physics::ThermalParameters& t = region_thermal[r];
        if (models.fermi_dirac && !device::is_insulator(region) &&
            (t.thermopower_exponent_n != -0.5 || t.thermopower_exponent_p != -0.5)) {
            return refused("electrothermal: region '" + region.name +
                           "': Fermi-Dirac statistics need thermopower exponents of -1/2 "
                           "(DECISIONS.md T3)");
        }
    }

    const mesh::Mesh& mesh = device.mesh();
    const std::size_t n = mesh.node_count();
    const int D = mesh.dimension();
    T0_ = scaling.temperature_K;
    thermal_bgn_ = models.bgn;
    masetti_electrons_ = models.electron_mobility == ElectronMobility::masetti;  // Unit 38
    const auto node_region = device.node_region();
    kappa_ref_ = region_thermal[node_region[device.reference_node()]].conductivity_W_cmK;
    heat_scale_ = kappa_ref_ * T0_ / (scaling.V_T * scaling.J0 * scaling.L_D);

    thermal_nodes_.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto node = static_cast<mesh::NodeId>(i);
        const device::RegionId r = node_region[i];
        ThermalNode& t = thermal_nodes_[i];
        t = ThermalNode{region_material[r], region_thermal[r], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        -1, 0.0, 0.0, 0};
        // rho c V L_D^2 / (t0 kappa_ref): the heat row's storage per unit rate and change of tau.
        t.capacity = t.thermal.heat_capacity_J_cm3K * volume_[i] * scaling.L_D * scaling.L_D /
                     (rate_scale_ * kappa_ref_);
        if (t.material < 0) continue;
        const physics::Semiconductor& m = thermal_materials_[static_cast<std::size_t>(t.material)];
        t.impurity = device.total_impurity(node);
        t.impurity_mobility = models.doping_mobility ? t.impurity : 0.0;
        t.narrowing_eV = models.bgn ? physics::bandgap_narrowing_eV(m, t.impurity) : 0.0;
        t.c_n = log_dos_n_[i] - band_shift_[i];
        t.c_p = -(log_dos_p_[i] + band_shift_[i]);
        t.gap0_eV = physics::band_gap_eV(m, T0_);
        t.r_n = t.thermal.thermopower_exponent_n;
        t.r_p = t.thermal.thermopower_exponent_p;
        t.log_nc = std::log(n_ie_[i]) + log_dos_n_[i];
        t.log_nv = std::log(n_ie_[i]) + log_dos_p_[i];
    }

    // Gates and electrodes with a p-polysilicon work function (chi + Eg(T)).
    for (const device::Contact& c : device.contacts()) {
        const bool gate = c.kind == device::ContactKind::gate &&
                          c.gate.electrode == device::GateElectrode::p_poly;
        const bool electrode = c.kind == device::ContactKind::electrode &&
                               c.electrode.kind == device::GateElectrode::p_poly;
        if (!gate && !electrode) continue;
        for (const mesh::NodeId v : c.nodes) {
            thermal_nodes_[static_cast<std::size_t>(v)].ppoly = gate ? 1 : 2;
        }
    }
    // Thermal contacts: the nodes of each patch; R_th as h = A L_D / (R_th kappa_ref) per node.
    const auto contacts = device.thermal_contacts();
    const double area_scale = std::pow(scaling.L_D, D - 1);
    for (std::size_t c = 0; c < contacts.size(); ++c) {
        const device::ThermalContact& tc = contacts[c];
        const mesh::BoundaryPatch* patch = mesh.find_boundary(tc.boundary);
        NITCAD_EXPECTS(patch != nullptr);  // checked by the device
        const bool isothermal = tc.kind == device::ThermalContactKind::isothermal;
        thermal_isothermal_.push_back(isothermal ? 1 : 0);
        thermal_rise_.push_back((tc.temperature_K - T0_) / T0_);
        for (std::size_t k = 0; k < patch->nodes.size(); ++k) {
            ThermalNode& t = thermal_nodes_[static_cast<std::size_t>(patch->nodes[k])];
            t.contact = static_cast<std::int32_t>(c);
            if (!isothermal) {
                t.h = patch->areas[k] / area_scale * scaling.L_D /
                      (tc.resistance_K_cm2_W * kappa_ref_);
            }
        }
    }
    // Every connected part needs a sink for a steady state (T7).
    std::vector<std::size_t> parent(n);
    std::iota(parent.begin(), parent.end(), std::size_t{0});
    const auto root = [&](std::size_t v) {
        while (parent[v] != v) v = parent[v] = parent[parent[v]];
        return v;
    };
    for (const EdgeTerm& e : edges_) parent[root(e.a)] = root(e.b);
    std::vector<char> sunk(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        if (thermal_nodes_[i].contact >= 0) sunk[root(i)] = 1;
    }
    sinks_complete_ = true;
    for (std::size_t i = 0; i < n; ++i) {
        if (sunk[root(i)] == 0) sinks_complete_ = false;
    }
    // Per node, whether its connected thermal domain has a sink (DECISIONS.md U5).
    reaches_sink_.resize(n);
    for (std::size_t i = 0; i < n; ++i) reaches_sink_[i] = sunk[root(i)];

    thermal_edges_.resize(edges_.size());
    for (std::size_t k = 0; k < edges_.size(); ++k) {
        const EdgeTerm& e = edges_[k];
        const ThermalNode& a = thermal_nodes_[e.a];
        const ThermalNode& b = thermal_nodes_[e.b];
        ThermalEdge& t = thermal_edges_[k];
        t.coupling = scaling.V_T / scaling.D0 * edge_geometry[k];
        t.log_nc = e.carriers ? b.log_nc - a.log_nc : 0.0;
        t.log_nv = e.carriers ? b.log_nv - a.log_nv : 0.0;
        t.r_n = 0.5 * (a.r_n + b.r_n);
        t.r_p = 0.5 * (a.r_p + b.r_p);
        t.kirchhoff = a.thermal == b.thermal;
        t.conduction = edge_geometry[k] / kappa_ref_;
        for (std::size_t r = 0; r < 4; ++r) {
            for (std::size_t c = 0; c < 4; ++c) {
                t.ab[4 * r + c] = detail::position(pattern_, 4 * e.a + r, 4 * e.b + c);
                t.ba[4 * r + c] = detail::position(pattern_, 4 * e.b + r, 4 * e.a + c);
            }
        }
    }
    metal_.assign(n, 0.0);
    contact_bias_.assign(contact_count_, 0.0);
    return std::nullopt;
}

bool DriftDiffusion::reaches_heat_sink(std::size_t node) const {
    NITCAD_EXPECTS(electrothermal() && node < node_count());
    return reaches_sink_[node] != 0;
}

std::expected<void, base::Error> DriftDiffusion::set_thermal_bias(
    std::span<const double> temperature_K) {
    NITCAD_EXPECTS(electrothermal());
    if (temperature_K.size() != thermal_rise_.size()) {
        return std::unexpected(base::Error{
            base::ErrorCode::invalid_input, "one temperature per thermal contact is needed",
            base::ErrorContext{.index = std::nullopt,
                               .value = static_cast<double>(temperature_K.size())}});
    }
    for (std::size_t c = 0; c < temperature_K.size(); ++c) {
        if (!(std::isfinite(temperature_K[c]) && temperature_K[c] > 0.0)) {
            return std::unexpected(base::Error{
                base::ErrorCode::invalid_input,
                "a thermal contact's temperature must be finite and positive",
                base::ErrorContext{.index = c, .value = temperature_K[c]}});
        }
    }
    for (std::size_t c = 0; c < temperature_K.size(); ++c) {
        thermal_rise_[c] = (temperature_K[c] - T0_) / T0_;
    }
    return {};
}

DriftDiffusion::ThermalState DriftDiffusion::thermal_state(std::size_t i, double tau) const {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double theta = 1.0 + tau;
    ThermalState s{theta, tau, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, 0.0, 0.0, 0.0, 0.0};
    if (!(std::isfinite(tau) && theta > 0.0)) return s;
    const ThermalNode& t = thermal_nodes_[i];
    const double T = theta * T0_;
    const physics::ValueSlope kappa = physics::thermal_conductivity(t.thermal, T);
    s.kappa = kappa.value;
    s.d_kappa = kappa.d_dT * T0_;
    // u = kappa(T0) (theta^(1-a) - 1) / (1 - a); kappa(T0) = kappa(T) theta^a.
    const double a = t.thermal.conductivity_exponent;
    s.u = kappa.value * std::pow(theta, a) * kirchhoff(tau, a);
    if (t.material < 0) {
        s.nie = s.d_nie = s.mu_n = s.d_mu_n = s.mu_p = s.d_mu_p = s.gap = s.d_gap = 0.0;
        return s;
    }
    const physics::Semiconductor& m = thermal_materials_[static_cast<std::size_t>(t.material)];
    const double nie = thermal_bgn_ ? physics::effective_intrinsic_density(m, t.impurity, T)
                                    : physics::intrinsic_density(m, T);
    s.nie = nie / Ns_;
    s.d_nie = s.nie * T0_ * physics::intrinsic_density_log_slope(m, T, t.narrowing_eV);
    const auto electron = physics::Carrier::electron, hole = physics::Carrier::hole;
    if (masetti_electrons_) {  // Unit 38
        const auto& masetti = physics::masetti_phosphorus_silicon;
        s.mu_n = physics::masetti_mobility(masetti, t.impurity_mobility, T);
        s.d_mu_n = T0_ * physics::masetti_mobility_slope(masetti, t.impurity_mobility, T);
    } else {
        s.mu_n = physics::caughey_thomas_mobility(m, electron, t.impurity_mobility, T);
        s.d_mu_n = T0_ * physics::caughey_thomas_mobility_slope(m, electron, t.impurity_mobility, T);
    }
    if (masetti_holes_[static_cast<std::size_t>(t.material)] != 0) {  // Units 34, 35
        const auto& masetti = physics::masetti_boron_silicon;
        s.mu_p = physics::masetti_mobility(masetti, t.impurity_mobility, T);
        s.d_mu_p = T0_ * physics::masetti_mobility_slope(masetti, t.impurity_mobility, T);
    } else {
        s.mu_p = physics::caughey_thomas_mobility(m, hole, t.impurity_mobility, T);
        s.d_mu_p = T0_ * physics::caughey_thomas_mobility_slope(m, hole, t.impurity_mobility, T);
    }
    s.gap = -(physics::band_gap_eV(m, T) - t.gap0_eV) / V_T_;
    s.d_gap = -T0_ * physics::band_gap_slope_eV_per_K(m, T) / V_T_;
    // Unit 35 (DECISIONS.md S4): the density-of-states masses' change from T0; 0 without a fit.
    if (physics::has_dos_mass_fit(m.parameters())) {
        const physics::ValueSlope cn = physics::dos_mass_log(m, physics::Carrier::electron, T);
        const physics::ValueSlope cp = physics::dos_mass_log(m, physics::Carrier::hole, T);
        s.mass_n = cn.value - physics::dos_mass_log(m, physics::Carrier::electron, T0_).value;
        s.d_mass_n = T0_ * cn.d_dT;
        s.mass_p = cp.value - physics::dos_mass_log(m, physics::Carrier::hole, T0_).value;
        s.d_mass_p = T0_ * cp.d_dT;
    }
    return s;
}

std::vector<DriftDiffusion::ThermalState> DriftDiffusion::thermal_states(
    std::span<const double> x) const {
    std::vector<ThermalState> s;
    s.reserve(node_count());
    for (std::size_t i = 0; i < node_count(); ++i) s.push_back(thermal_state(i, x[4 * i + 3]));
    return s;
}

DriftDiffusion::Level DriftDiffusion::level(double c, double theta, double nie, double log_N,
                                            double mass, double d_mass) const {
    const double l15 = 1.5 * std::log(theta);
    const double log_N_theta = log_N + l15 + mass;  // N(theta) = N0 theta^1.5 M(theta)
    Level v{std::log(c) - log_N_theta, 1.0 / c, -1.5 / theta - d_mass, 0, 0, 0, 0, 0, 0};
    if (!fermi_dirac_) return v;
    // ln gamma depends on c / N(theta) alone: d / d theta = -(1.5 / theta + d_mass) c d / dc.
    const physics::Degeneracy g =
        physics::fermi_dirac_degeneracy(nie, log_N_theta - std::log(nie), c);
    v.log_gamma = g.log_gamma;
    v.dg_dc = g.d_density;
    v.dg_dt = -1.5 / theta * c * g.d_density - d_mass * c * g.d_density;
    v.x -= g.log_gamma;
    v.dx_dc -= v.dg_dc;
    v.dx_dt -= v.dg_dt;
    const physics::ThermalDiffusionFactor h = physics::fermi_dirac_thermal_diffusion(v.x);
    v.extra = h.value - 0.5;
    v.dh_dc = h.d_x * v.dx_dc;
    v.dh_dt = h.d_x * v.dx_dt;
    return v;
}

std::vector<DriftDiffusion::NodeLevels> DriftDiffusion::thermal_levels(
    std::span<const double> x, std::span<const ThermalState> t) const {
    std::vector<NodeLevels> g;
    if (!fermi_dirac_) return g;
    g.resize(node_count());
    for (std::size_t i = 0; i < node_count(); ++i) {
        if (insulator_[i] != 0) continue;
        const ThermalNode& tn = thermal_nodes_[i];
        g[i] = {level(x[4 * i + 1], t[i].theta, t[i].nie, tn.log_nc, t[i].mass_n, t[i].d_mass_n),
                level(x[4 * i + 2], t[i].theta, t[i].nie, tn.log_nv, t[i].mass_p, t[i].d_mass_p)};
    }
    return g;
}

std::pair<DriftDiffusion::ThermalFlux, DriftDiffusion::ThermalFlux>
DriftDiffusion::thermal_fluxes(std::size_t k, std::span<const double> x,
                               std::span<const ThermalState> t,
                               std::span<const NodeLevels> g, SurfaceSlope* slope) const {
    const EdgeTerm& e = edges_[k];
    const ThermalEdge& te = thermal_edges_[k];
    const std::size_t a = e.a, b = e.b;
    const double dpsi = x[4 * b] - x[4 * a];
    const double ta = t[a].theta, tb = t[b].theta;
    const double th = 0.5 * (ta + tb);
    // Fermi-Dirac: delta_n gains (ln gamma_b - ln gamma_a) - hbar D and delta_p loses it, with
    // hbar the ends' mean excess of the thermal diffusion factor and D = (theta_b - theta_a) / th
    // (the nu-factor scheme of Unit 14 with the thermopower's degeneracy term). Its value and
    // partials in (c_a, c_b) and (theta_a, theta_b), per carrier.
    struct Extra {
        double value = 0.0, d_c[2] = {0.0, 0.0}, d_t[2] = {0.0, 0.0};
    };
    const auto extra = [&](const Level& la, const Level& lb) {
        Extra v{};
        const double D = (t[b].tau - t[a].tau) / th;
        const double dD[2] = {-1.0 / th - 0.5 * D / th, 1.0 / th - 0.5 * D / th};
        const double hbar = 0.5 * (la.extra + lb.extra);
        v.value = lb.log_gamma - la.log_gamma - hbar * D;
        v.d_c[0] = -la.dg_dc - 0.5 * la.dh_dc * D;
        v.d_c[1] = lb.dg_dc - 0.5 * lb.dh_dc * D;
        v.d_t[0] = -la.dg_dt - 0.5 * la.dh_dt * D - hbar * dD[0];
        v.d_t[1] = lb.dg_dt - 0.5 * lb.dh_dt * D - hbar * dD[1];
        return v;
    };
    Extra xn{}, xp{};
    if (!g.empty()) {
        xn = extra(g[a].n, g[b].n);
        xp = extra(g[a].p, g[b].p);
        xp.value = -xp.value;
        for (std::size_t j = 0; j < 2; ++j) {
            xp.d_c[j] = -xp.d_c[j];
            xp.d_t[j] = -xp.d_t[j];
        }
    }
    // Unit 35 (DECISIONS.md S4): -k T grad ln M of a temperature-dependent density-of-states mass
    // is a band-parameter step, delta_n gains ln M_c(theta_b) - ln M_c(theta_a) as it has L_c,
    // delta_p loses ln M_v(theta_b) - ln M_v(theta_a) as it has -L_v. 0 without a mass fit.
    if (t[a].d_mass_n != 0.0 || t[b].d_mass_n != 0.0) {
        xn.value += t[b].mass_n - t[a].mass_n;
        xn.d_t[0] -= t[a].d_mass_n;
        xn.d_t[1] += t[b].d_mass_n;
    }
    if (t[a].d_mass_p != 0.0 || t[b].d_mass_p != 0.0) {
        xp.value -= t[b].mass_p - t[a].mass_p;
        xp.d_t[0] += t[a].d_mass_p;
        xp.d_t[1] -= t[b].d_mass_p;
    }
    // One carrier: flux = sign A (w_b B(sign delta) - w_a B(-sign delta)) with sign +1 for
    // electrons (A_n, delta_n) and -1 for holes, where the hole flux -A (w_b B(-delta) -
    // w_a B(delta)) is the electron form at -delta, negated. base: delta = base / th + L + extra,
    // and dbase_a, dbase_b its theta partials beyond th's (the valence band's shift, holes).
    // The edge mobility from the ends' low-field mobilities at their T: their harmonic mean, or
    // with field mobility the Canali mobility of it at the edge's field E = field |dpsi| (on an
    // edge between two materials the harmonic mean of each end's Canali mobility), as the
    // isothermal edges. Canali's mu = mu0 f(mu0 E), so d mu / d mu0 = (mu + E d mu / dE) / mu0.
    // Value and partials in theta_a, theta_b and dpsi.
    struct Mobility {
        double value, d_a, d_b, d_psi;
    };
    const double E = e.field * std::abs(dpsi);
    const double dE = dpsi > 0.0 ? e.field : dpsi < 0.0 ? -e.field : 0.0;  // dE / d dpsi
    const auto mobility = [&](double mu_a, double d_mu_a, double mu_b, double d_mu_b,
                              const physics::CanaliParameters& sat_a,
                              const physics::CanaliParameters& sat_b) {
        const auto hmean = [](double u, double v) {  // and its partials
            const double sum = u + v;
            return std::array<double, 3>{2.0 * u * v / sum, 2.0 * v * v / (sum * sum),
                                         2.0 * u * u / (sum * sum)};
        };
        if (!field_mobility_) {
            const auto h = hmean(mu_a, mu_b);
            return Mobility{h[0], h[1] * d_mu_a, h[2] * d_mu_b, 0.0};
        }
        if (!e.mixed) {
            const auto h = hmean(mu_a, mu_b);
            const physics::FieldMobility c = physics::canali_mobility(h[0], E, sat_a);
            const double d_mu0 = (c.mobility + E * c.d_dE) / h[0];
            return Mobility{c.mobility, d_mu0 * h[1] * d_mu_a, d_mu0 * h[2] * d_mu_b,
                            c.d_dE * dE};
        }
        const physics::FieldMobility ca = physics::canali_mobility(mu_a, E, sat_a);
        const physics::FieldMobility cb = physics::canali_mobility(mu_b, E, sat_b);
        const auto h = hmean(ca.mobility, cb.mobility);
        return Mobility{h[0], h[1] * (ca.mobility + E * ca.d_dE) / mu_a * d_mu_a,
                        h[2] * (cb.mobility + E * cb.d_dE) / mu_b * d_mu_b,
                        (h[1] * ca.d_dE + h[2] * cb.d_dE) * dE};
    };
    const auto carrier = [&](double sign, double r, const Mobility& mu, double c_a, double c_b,
                             double base, double L, double dbase_a, double dbase_b,
                             const Extra& xd) {
        ThermalFlux f{};
        const double delta = base / th + L + xd.value;
        const double s = sign * delta;  // the electron form's argument
        const double ex = 1.0 + r;
        const double pa = std::pow(ta, ex), pb = std::pow(tb, ex);
        const double wa = pa * c_a, wb = pb * c_b;
        const double q = std::pow(th, -r);
        const double A = mu.value * te.coupling * q;
        const double dA_a = te.coupling * q * mu.d_a - 0.5 * r * A / th;
        const double dA_b = te.coupling * q * mu.d_b - 0.5 * r * A / th;
        const double dA_psi = te.coupling * q * mu.d_psi;
        const double bp = bernoulli(s), bm = bernoulli(-s);
        const double unit = sign * (wb * bp - wa * bm);
        // d flux / d delta = sign A (wb B'(s) + wa B'(-s)) sign.
        const double d_delta =
            A * (wb * bernoulli_derivative(s) + wa * bernoulli_derivative(-s));
        const double dw_a = -sign * A * bm, dw_b = sign * A * bp;
        f.flux = A * unit;
        f.d_psi[0] = -d_delta / th - dA_psi * unit;
        f.d_psi[1] = d_delta / th + dA_psi * unit;
        f.d_c[0] = dw_a * pa + d_delta * xd.d_c[0];
        f.d_c[1] = dw_b * pb + d_delta * xd.d_c[1];
        const double d_th = -0.5 * base / (th * th);  // d delta / d theta_a (and _b), via th
        f.d_theta[0] = dA_a * unit + d_delta * (d_th + dbase_a / th + xd.d_t[0]) +
                       dw_a * ex * pa / ta * c_a;
        f.d_theta[1] = dA_b * unit + d_delta * (d_th + dbase_b / th + xd.d_t[1]) +
                       dw_b * ex * pb / tb * c_b;
        f.delta = delta;
        f.theta_e = th;
        f.w[0] = wa;
        f.w[1] = wb;
        return f;
    };
    const double base_n = dpsi + e.shift_n - te.log_nc;
    const double base_p = dpsi + e.shift_p + te.log_nv - (t[b].gap - t[a].gap);
    if (slope != nullptr) *slope = SurfaceSlope{};
    if (!e.surface) {
        const Mobility mu_n =
            mobility(t[a].mu_n, t[a].d_mu_n, t[b].mu_n, t[b].d_mu_n, e.sat_n, e.sat_n_b);
        const Mobility mu_p =
            mobility(t[a].mu_p, t[a].d_mu_p, t[b].mu_p, t[b].d_mu_p, e.sat_p, e.sat_p_b);
        return {carrier(1.0, te.r_n, mu_n, x[4 * a + 1], x[4 * b + 1], base_n, te.log_nc, 0.0,
                        0.0, xn),
                carrier(-1.0, te.r_p, mu_p, x[4 * a + 2], x[4 * b + 2], base_p, -te.log_nv,
                        t[a].d_gap, -t[b].d_gap, xp)};
    }
    // Lombardi (Unit 46): each end's low-field mobility is its bulk mobility with the surface
    // terms at its own T in the edge's E_perp; d mu_end / d theta = (mu_end / mu_b)^2 d mu_b /
    // d theta + T0 d mu_end / dT. The edge mobility is linear in its ends' changes, so its
    // E_perp slope is mobility() with the ends' E_perp slopes in place of their theta slopes.
    const SurfaceField perp = surface_field(e, x);
    const auto with_surface = [&](physics::Carrier c, double mu_a, double d_mu_a, double mu_b,
                                  double d_mu_b, const physics::CanaliParameters& sat_a,
                                  const physics::CanaliParameters& sat_b) {
        const physics::SurfaceMobility la = surface_end(c, mu_a, e.N_a, ta * T0_, perp.value);
        const physics::SurfaceMobility lb = surface_end(c, mu_b, e.N_b, tb * T0_, perp.value);
        const Mobility m = mobility(la.mobility, la.d_dbulk * d_mu_a + T0_ * la.d_dT,
                                    lb.mobility, lb.d_dbulk * d_mu_b + T0_ * lb.d_dT, sat_a,
                                    sat_b);
        const Mobility s = mobility(la.mobility, la.d_dE, lb.mobility, lb.d_dE, sat_a, sat_b);
        return std::pair{m, s.d_a + s.d_b};
    };
    const auto [mu_n, perp_n] = with_surface(physics::Carrier::electron, t[a].mu_n, t[a].d_mu_n,
                                             t[b].mu_n, t[b].d_mu_n, e.sat_n, e.sat_n_b);
    const auto [mu_p, perp_p] = with_surface(physics::Carrier::hole, t[a].mu_p, t[a].d_mu_p,
                                             t[b].mu_p, t[b].d_mu_p, e.sat_p, e.sat_p_b);
    const ThermalFlux fn =
        carrier(1.0, te.r_n, mu_n, x[4 * a + 1], x[4 * b + 1], base_n, te.log_nc, 0.0, 0.0, xn);
    const ThermalFlux fp = carrier(-1.0, te.r_p, mu_p, x[4 * a + 2], x[4 * b + 2], base_p,
                                   -te.log_nv, t[a].d_gap, -t[b].d_gap, xp);
    // The flux is linear in the edge mobility: d flux / d E_perp = flux / mu d mu / d E_perp.
    if (slope != nullptr) {
        *slope = SurfaceSlope{perp, fn.flux / mu_n.value * perp_n, fp.flux / mu_p.value * perp_p};
    }
    return {fn, fp};
}

DriftDiffusion::Conduction DriftDiffusion::conduction(std::size_t k,
                                                      std::span<const ThermalState> t) const {
    const EdgeTerm& e = edges_[k];
    const ThermalEdge& te = thermal_edges_[k];
    const ThermalState& a = t[e.a];
    const ThermalState& b = t[e.b];
    const double g = te.conduction;
    if (te.kirchhoff) return {g * (a.u - b.u), g * a.kappa, -g * b.kappa};
    const double sum = a.kappa + b.kappa;
    const double hm = 2.0 * a.kappa * b.kappa / sum;
    const double drop = a.tau - b.tau;
    return {g * hm * drop,
            g * (hm + drop * 2.0 * b.kappa * b.kappa / (sum * sum) * a.d_kappa),
            g * (-hm + drop * 2.0 * a.kappa * a.kappa / (sum * sum) * b.d_kappa)};
}

DriftDiffusion::ThermalOhmic DriftDiffusion::thermal_ohmic(std::size_t i,
                                                           const ThermalState& t) const {
    const ThermalNode& tn = thermal_nodes_[i];
    const double theta = t.theta;
    const double lt = 1.5 * std::log(theta), ln_nie = std::log(t.nie);
    // Complete ionization (incomplete ionization is refused, DECISIONS.md U6).
    const physics::NeutralEquilibrium q =
        fermi_dirac_ ? physics::fermi_dirac_neutral_equilibrium(
                           doping_[i], t.nie, tn.log_nc + lt + t.mass_n - ln_nie,
                           tn.log_nv + lt + t.mass_p - ln_nie)
                     : physics::boltzmann_neutral_equilibrium(doping_[i], t.nie);
    // With E_F = -m (m = V / V_T) the band reduced energies are x_n = (psi - m - c_n) / theta and
    // x_p = (c_p + g + m - psi) / theta, so psi0 = m + c_n + theta x_n. Neutrality
    // G = n(x_n, theta) - p(x_p, theta) - C = 0 fixes psi0(theta): dpsi0 / dtheta = -G_t / G_psi,
    // with dn / dx = n / (1 - n d ln gamma / dn) (n under Boltzmann statistics) and
    // d n / d theta = (1.5 / theta + d ln M / d theta) n at fixed x.
    const Level ln = level(q.n, theta, t.nie, tn.log_nc, t.mass_n, t.d_mass_n);
    const Level lp = level(q.p, theta, t.nie, tn.log_nv, t.mass_p, t.d_mass_p);
    const double sn = q.n / (1.0 - q.n * ln.dg_dc), sp = q.p / (1.0 - q.p * lp.dg_dc);
    const double G_psi = (sn + sp) / theta;
    const double G_t = 1.5 * q.n / theta - sn * ln.x / theta - 1.5 * q.p / theta -
                       sp * (t.d_gap - lp.x) / theta + t.d_mass_n * q.n - t.d_mass_p * q.p;
    const double d_psi_t = -G_t / G_psi;  // of theta x_n, which psi0 - m - c_n is
    const double dn = 1.5 * q.n / theta + sn * (d_psi_t - ln.x) / theta + t.d_mass_n * q.n;
    const double dp = 1.5 * q.p / theta + sp * (t.d_gap - d_psi_t - lp.x) / theta +
                      t.d_mass_p * q.p;
    return {metal_[i] + tn.c_n + theta * ln.x, q.n, q.p, d_psi_t, dn, dp};
}

DriftDiffusion::Shift DriftDiffusion::electrode_shift(std::size_t i,
                                                      const ThermalState& t) const {
    // psi_E = (V - phi_m(T) + reference depth at T0) / V_T: a p-polysilicon phi_m = chi + Eg(T)
    // moves it by -(Eg(T) - Eg(T0)) / k T0; for a gate that is the node's own g(theta) (the gate's
    // polysilicon is the node's material, gate.hpp), for an electrode silicon's.
    const std::int32_t kind = thermal_nodes_[i].ppoly;
    if (kind == 1) return {t.gap, t.d_gap};
    if (kind == 2) {
        const physics::Semiconductor poly = physics::silicon();
        const double T = t.theta * T0_;
        return {-(physics::band_gap_eV(poly, T) - physics::band_gap_eV(poly, T0_)) / V_T_,
                -T0_ * physics::band_gap_slope_eV_per_K(poly, T) / V_T_};
    }
    return {0.0, 0.0};
}

void DriftDiffusion::add_gate_shift(std::span<const double> x, std::span<double> charge) const {
    for (std::size_t i = 0; i < node_count(); ++i) {
        if (!gates_.on_gate(i) || thermal_nodes_[i].ppoly != 1) continue;
        charge[gates_.contact(i)] +=
            gates_.term(i).coupling * electrode_shift(i, thermal_state(i, x[4 * i + 3])).value;
    }
}

void DriftDiffusion::assemble_heat(std::span<const double> x, std::span<double> f,
                                   std::span<double> values, const TimeStep* step,
                                   bool sinks) const {
    NITCAD_EXPECTS(x.size() == unknowns() && f.size() == unknowns());
    NITCAD_EXPECTS(step == nullptr || (step->storage.size() == 3 * node_count() &&
                                       step->traps.size() == trap_slots()));
    // The electrothermal terms in assembly order (terms.hpp, thermal_terms.cpp). Every one, the
    // interfaces' among them, also for small_signal_matrix: its electrothermal form adds no
    // interface terms of its own.
    const std::vector<ThermalState> t = thermal_states(x);
    const std::vector<NodeLevels> g = thermal_levels(x, t);
    composition_->assemble(*this, {x, f, values, step, {}, t, g, sinks}, true);
}

void DriftDiffusion::add_thermal_storage(std::span<const double> x, std::complex<double> s,
                                         std::span<std::complex<double>> values) const {
    // assemble_heat's step terms with rate s and the history equal to the state (a steady
    // state), linearized: d(S - c) = dS, so only the terms linear in the change remain.
    const std::vector<ThermalState> t = thermal_states(x);
    const std::vector<NodeLevels> g = thermal_levels(x, t);
    const auto at = [&](std::size_t node, std::size_t r, std::size_t c) -> std::complex<double>& {
        return values[block_[16 * node + 4 * r + c]];
    };
    const double inv_K0 = 1.0 / heat_scale_;
    for (std::size_t i = 0; i < node_count(); ++i) {
        const ThermalNode& tn = thermal_nodes_[i];
        const bool held = tn.contact >= 0 &&
                          thermal_isothermal_[static_cast<std::size_t>(tn.contact)] != 0;
        if (!held) at(i, 3, 3) += s * tn.capacity;
        if (contact_[i] >= 0 || insulator_[i] != 0) continue;
        const double V = volume_[i];
        at(i, 1, 1) -= s * V;  // S_n = n, S_p = p (complete ionization, U6)
        at(i, 2, 2) += s * V;
        if (held) continue;
        const Level none{};
        const Level& ln = g.empty() ? none : g[i].n;
        const Level& lp = g.empty() ? none : g[i].p;
        const double psi = x[4 * i], theta = t[i].theta;
        const double en = tn.c_n - psi + (tn.r_n + 2.5 + ln.extra) * theta;
        const double ep = tn.c_p - psi + t[i].gap - (tn.r_p + 2.5 + lp.extra) * theta;
        at(i, 3, 1) += s * (inv_K0 * V * en);
        at(i, 3, 2) -= s * (inv_K0 * V * ep);
    }
}

std::vector<double> DriftDiffusion::bias_derivative(std::size_t contact,
                                                    std::span<const double> x) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    std::vector<double> d = bias_derivative(contact);
    if (!electrothermal()) {
        // Gate tunnelling through a lumped gate (Unit 51): its node's rows -T_n and +T_p read the
        // gate's bias (per volt).
        if (!gate_tunnels_.empty()) {
            const std::vector<NodeDegeneracy> g = degeneracies(x);
            GateTunnelFlux flux;
            for (const GateTunnel& t : gate_tunnels_) {
                if (!t.lumped || static_cast<std::size_t>(t.contact) != contact) continue;
                gate_tunnel_flux(t, x, g, flux);
                d[m_ * t.node + 1] -= flux.d_n[3];
                d[m_ * t.node + 2] += flux.d_p[3];
            }
        }
        detail::observe("bias_derivative.state", d);
        return d;
    }
    // The heat row of a node of the contact gains -m (Jn + Jp) / K0 per edge it is the first node
    // of and + per edge it is the second of (assemble_heat), m = V / V_T; isothermal nodes hold T.
    const auto c = static_cast<std::int32_t>(contact);
    const std::vector<ThermalState> t = thermal_states(x);
    const std::vector<NodeLevels> g = thermal_levels(x, t);
    const auto loose = [&](std::size_t i) {
        const std::int32_t s = thermal_nodes_[i].contact;
        return contact_[i] == c &&
               !(s >= 0 && thermal_isothermal_[static_cast<std::size_t>(s)] != 0);
    };
    const double k = 1.0 / (V_T_ * heat_scale_);
    for (std::size_t q = 0; q < edges_.size(); ++q) {
        const EdgeTerm& e = edges_[q];
        if (!e.carriers || (!loose(e.a) && !loose(e.b))) continue;
        const auto [fn, fp] = thermal_fluxes(q, x, t, g);
        const double J = fn.flux + fp.flux;
        if (loose(e.a)) d[4 * e.a + 3] -= k * J;
        if (loose(e.b)) d[4 * e.b + 3] += k * J;
    }
    detail::observe("bias_derivative.state", d);
    return d;
}

std::vector<std::size_t> DriftDiffusion::bias_rows(std::size_t contact) const {
    NITCAD_EXPECTS(contact < contact_count_);
    const auto c = static_cast<std::int32_t>(contact);
    std::vector<std::size_t> rows;
    for (std::size_t i = 0; i < node_count(); ++i) {
        const bool own = contact_[i] == c || electrode_[i] == c ||
                         (gates_.on_gate(i) && gates_.contact(i) == contact);
        if (!own) continue;
        rows.push_back(m_ * i);
        if (!electrothermal() || contact_[i] != c) continue;
        const std::int32_t s = thermal_nodes_[i].contact;
        if (!(s >= 0 && thermal_isothermal_[static_cast<std::size_t>(s)] != 0)) {
            rows.push_back(4 * i + 3);
        }
    }
    for (const GateTunnel& t : gate_tunnels_) {  // Unit 51: a lumped gate's tunnelling rows
        if (!t.lumped || t.contact != c) continue;
        rows.push_back(m_ * t.node + 1);
        rows.push_back(m_ * t.node + 2);
    }
    return rows;
}

bool DriftDiffusion::bias_nonlinear() const noexcept {
    return electrothermal() ||
           std::ranges::any_of(gate_tunnels_, [](const GateTunnel& t) { return t.lumped; });
}

std::vector<double> DriftDiffusion::thermal_contact_heat(std::span<const double> x,
                                                         const TimeStep* step) const {
    NITCAD_EXPECTS(electrothermal() && x.size() == unknowns());
    std::vector<double> f(unknowns());
    assemble_heat(x, f, {}, step, false);
    std::vector<double> heat(thermal_rise_.size(), 0.0);
    for (std::size_t i = 0; i < node_count(); ++i) {
        const ThermalNode& tn = thermal_nodes_[i];
        if (tn.contact < 0) continue;
        const auto c = static_cast<std::size_t>(tn.contact);
        // The sink takes what the node's balance leaves: conducted out + sink = received.
        heat[c] += thermal_isothermal_[c] != 0 ? -f[4 * i + 3]
                                                : tn.h * (x[4 * i + 3] - thermal_rise_[c]);
    }
    for (double& h : heat) h *= heat_scale_;
    detail::observe("thermal_contact_heat", heat);
    return heat;
}

double DriftDiffusion::electrical_power(std::span<const double> x) const {
    NITCAD_EXPECTS((electrothermal() || energy_balance_) && x.size() == unknowns());
    const std::vector<double> I = terminal_currents(x);
    double power = 0.0;
    for (std::size_t c = 0; c < I.size(); ++c) power += contact_bias_[c] * I[c];
    detail::observe("electrical_power", power);
    return power;
}

std::vector<std::pair<double, double>> DriftDiffusion::edge_joule_heat(
    std::span<const double> x) const {
    NITCAD_EXPECTS(electrothermal() && x.size() == unknowns());
    const std::vector<ThermalState> t = thermal_states(x);
    const std::vector<NodeLevels> g = thermal_levels(x, t);
    std::vector<std::pair<double, double>> joule(edges_.size(), {0.0, 0.0});
    for (std::size_t k = 0; k < edges_.size(); ++k) {
        if (!edges_[k].carriers) continue;
        const auto [fn, fp] = thermal_fluxes(k, x, t, g);
        joule[k] = {fn.flux * fn.theta_e * (std::log(fn.w[1] / fn.w[0]) - fn.delta),
                    fp.flux * fp.theta_e * (std::log(fp.w[0] / fp.w[1]) - fp.delta)};
    }
    detail::observe("edge_joule_heat", joule);
    return joule;
}

}  // namespace NiTCAD::assemble
