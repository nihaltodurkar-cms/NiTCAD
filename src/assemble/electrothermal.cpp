// The electrothermal rows of DriftDiffusion (Unit 23; drift_diffusion.hpp, "Electrothermal
// coupling"; DECISIONS.md T1-T14).
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <string>
#include <utility>
#include <variant>

#include "NiTCAD/assemble/bernoulli.hpp"
#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/contract.hpp"
#include "NiTCAD/physics/bandgap_narrowing.hpp"
#include "NiTCAD/physics/field_mobility.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/recombination.hpp"
#include "NiTCAD/physics/statistics.hpp"
#include "NiTCAD/physics/thermal.hpp"
#include "scaled_device.hpp"

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
    }
    // Models whose temperature coupling is not built yet.
    const auto pending = [](const char* what) {
        return refused(std::string("electrothermal does not yet support ") + what);
    };
    if (models.incomplete_ionization) return pending("incomplete ionization");
    if (models.impact_ionization) return pending("impact ionization");
    for (const device::Contact& c : device.contacts()) {
        if (c.kind != device::ContactKind::ohmic) return pending("gates and electrodes");
    }
    for (const device::Interface& f : device.interfaces()) {
        if (device::has_interface_charge_or_recombination(f)) {
            return pending("interface charge and surface recombination");
        }
        if (f.transport == device::InterfaceTransport::thermionic_emission) {
            return pending("thermionic emission");
        }
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
    const auto node_region = device.node_region();
    kappa_ref_ = region_thermal[node_region[device.reference_node()]].conductivity_W_cmK;
    heat_scale_ = kappa_ref_ * T0_ / (scaling.V_T * scaling.J0 * scaling.L_D);

    thermal_nodes_.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto node = static_cast<mesh::NodeId>(i);
        const device::RegionId r = node_region[i];
        ThermalNode& t = thermal_nodes_[i];
        t = ThermalNode{region_material[r], region_thermal[r], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        -1, 0.0};
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
    ThermalState s{theta, tau, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan, nan};
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
    s.mu_n = physics::caughey_thomas_mobility(m, electron, t.impurity_mobility, T);
    s.mu_p = physics::caughey_thomas_mobility(m, hole, t.impurity_mobility, T);
    s.d_mu_n = T0_ * physics::caughey_thomas_mobility_slope(m, electron, t.impurity_mobility, T);
    s.d_mu_p = T0_ * physics::caughey_thomas_mobility_slope(m, hole, t.impurity_mobility, T);
    s.gap = -(physics::band_gap_eV(m, T) - t.gap0_eV) / V_T_;
    s.d_gap = -T0_ * physics::band_gap_slope_eV_per_K(m, T) / V_T_;
    return s;
}

std::vector<DriftDiffusion::ThermalState> DriftDiffusion::thermal_states(
    std::span<const double> x) const {
    std::vector<ThermalState> s;
    s.reserve(node_count());
    for (std::size_t i = 0; i < node_count(); ++i) s.push_back(thermal_state(i, x[4 * i + 3]));
    return s;
}

DriftDiffusion::Level DriftDiffusion::level(double c, double theta, double nie,
                                            double log_N) const {
    const double l15 = 1.5 * std::log(theta);
    const double log_N_theta = log_N + l15;  // N(theta) = N0 theta^1.5
    Level v{std::log(c) - log_N_theta, 1.0 / c, -1.5 / theta, 0, 0, 0, 0, 0, 0};
    if (!fermi_dirac_) return v;
    // ln gamma depends on c / N(theta) alone: d / d theta = -(1.5 / theta) c d / dc.
    const physics::Degeneracy g =
        physics::fermi_dirac_degeneracy(nie, log_N_theta - std::log(nie), c);
    v.log_gamma = g.log_gamma;
    v.dg_dc = g.d_density;
    v.dg_dt = -1.5 / theta * c * g.d_density;
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
        g[i] = {level(x[4 * i + 1], t[i].theta, t[i].nie, tn.log_nc),
                level(x[4 * i + 2], t[i].theta, t[i].nie, tn.log_nv)};
    }
    return g;
}

std::pair<DriftDiffusion::ThermalFlux, DriftDiffusion::ThermalFlux>
DriftDiffusion::thermal_fluxes(std::size_t k, std::span<const double> x,
                               std::span<const ThermalState> t,
                               std::span<const NodeLevels> g) const {
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
        f.flux = A * unit;
        // d flux / d delta = sign A (wb B'(s) + wa B'(-s)) sign.
        const double d_delta = A * (wb * bernoulli_derivative(s) + wa * bernoulli_derivative(-s));
        const double dw_a = -sign * A * bm, dw_b = sign * A * bp;
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
    const Mobility mu_n =
        mobility(t[a].mu_n, t[a].d_mu_n, t[b].mu_n, t[b].d_mu_n, e.sat_n, e.sat_n_b);
    const Mobility mu_p =
        mobility(t[a].mu_p, t[a].d_mu_p, t[b].mu_p, t[b].d_mu_p, e.sat_p, e.sat_p_b);
    return {carrier(1.0, te.r_n, mu_n, x[4 * a + 1], x[4 * b + 1], base_n, te.log_nc, 0.0, 0.0,
                    xn),
            carrier(-1.0, te.r_p, mu_p, x[4 * a + 2], x[4 * b + 2], base_p, -te.log_nv,
                    t[a].d_gap, -t[b].d_gap, xp)};
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
    const physics::NeutralEquilibrium q =
        fermi_dirac_ ? physics::fermi_dirac_neutral_equilibrium(
                           doping_[i], t.nie, tn.log_nc + lt - ln_nie, tn.log_nv + lt - ln_nie)
                     : physics::boltzmann_neutral_equilibrium(doping_[i], t.nie);
    // With E_F = -m (m = V / V_T) the band reduced energies are x_n = (psi - m - c_n) / theta and
    // x_p = (c_p + g + m - psi) / theta, so psi0 = m + c_n + theta x_n. Neutrality
    // G = n(x_n, theta) - p(x_p, theta) - C = 0 fixes psi0(theta): dpsi0 / dtheta = -G_t / G_psi,
    // with dn / dx = n / (1 - n d ln gamma / dn) (n under Boltzmann statistics) and
    // d n / d theta = 1.5 n / theta at fixed x.
    const Level ln = level(q.n, theta, t.nie, tn.log_nc);
    const Level lp = level(q.p, theta, t.nie, tn.log_nv);
    const double sn = q.n / (1.0 - q.n * ln.dg_dc), sp = q.p / (1.0 - q.p * lp.dg_dc);
    const double G_psi = (sn + sp) / theta;
    const double G_t = 1.5 * q.n / theta - sn * ln.x / theta - 1.5 * q.p / theta -
                       sp * (t.d_gap - lp.x) / theta;
    const double d_psi_t = -G_t / G_psi;  // of theta x_n, which psi0 - m - c_n is
    const double dn = 1.5 * q.n / theta + sn * (d_psi_t - ln.x) / theta;
    const double dp = 1.5 * q.p / theta + sp * (t.d_gap - d_psi_t - lp.x) / theta;
    return {metal_[i] + tn.c_n + theta * ln.x, q.n, q.p, d_psi_t, dn, dp};
}

void DriftDiffusion::assemble_heat(std::span<const double> x, std::span<double> f,
                                   std::span<double> values, const TimeStep* step,
                                   bool sinks) const {
    NITCAD_EXPECTS(x.size() == unknowns() && f.size() == unknowns());
    NITCAD_EXPECTS(step == nullptr);  // thermal transients follow
    const bool jacobian = !values.empty();
    const std::vector<ThermalState> t = thermal_states(x);
    const std::vector<NodeLevels> g = thermal_levels(x, t);
    const auto at = [&](std::size_t node, std::size_t r, std::size_t c) -> double& {
        return values[block_[16 * node + 4 * r + c]];
    };
    // Whether node i's heat row is a sink's Dirichlet row.
    const auto held = [&](std::size_t i) {
        const std::int32_t c = thermal_nodes_[i].contact;
        return sinks && c >= 0 && thermal_isothermal_[static_cast<std::size_t>(c)] != 0;
    };
    for (std::size_t i = 0; i < node_count(); ++i) {
        const double psi = x[4 * i], n = x[4 * i + 1], p = x[4 * i + 2], tau = x[4 * i + 3];
        const ThermalNode& tn = thermal_nodes_[i];
        if (held(i)) {
            f[4 * i + 3] = tau - thermal_rise_[static_cast<std::size_t>(tn.contact)];
            if (jacobian) at(i, 3, 3) = 1.0;
        } else if (sinks && tn.h > 0.0) {
            f[4 * i + 3] =
                tn.h * (tau - thermal_rise_[static_cast<std::size_t>(tn.contact)]);
            if (jacobian) at(i, 3, 3) = tn.h;
        } else {
            f[4 * i + 3] = 0.0;
        }
        if (contact_[i] >= 0) {
            const ThermalOhmic o = thermal_ohmic(i, t[i]);
            f[4 * i] = psi - o.psi;
            f[4 * i + 1] = n - o.n;
            f[4 * i + 2] = p - o.p;
            if (jacobian) {
                at(i, 0, 0) = at(i, 1, 1) = at(i, 2, 2) = 1.0;
                at(i, 0, 3) = -o.d_psi;
                at(i, 1, 3) = -o.d_n;
                at(i, 2, 3) = -o.d_p;
            }
            continue;
        }
        if (insulator_[i] != 0) {  // no charge, no carriers (electrodes are refused)
            f[4 * i] = 0.0;
            f[4 * i + 1] = n;
            f[4 * i + 2] = p;
            if (jacobian) at(i, 1, 1) = at(i, 2, 2) = 1.0;
            continue;
        }
        const double V = volume_[i];
        f[4 * i] = -V * (n - p - doping_[i]);
        f[4 * i + 1] = 0.0;
        f[4 * i + 2] = 0.0;
        if (jacobian) {
            at(i, 0, 1) = -V;
            at(i, 0, 2) = V;
        }
        const double nie = t[i].nie, d_nie = t[i].d_nie;
        // Recombination R / R0 and its partials; d_t its theta derivative through n_ie.
        const auto add = [&](double k, double rate, double d_n, double d_p, double d_t) {
            f[4 * i + 1] -= k * rate;
            f[4 * i + 2] += k * rate;
            if (jacobian) {
                at(i, 1, 1) -= k * d_n;
                at(i, 1, 2) -= k * d_p;
                at(i, 1, 3) -= k * d_t;
                at(i, 2, 1) += k * d_n;
                at(i, 2, 2) += k * d_p;
                at(i, 2, 3) += k * d_t;
            }
        };
        const double k = V * rate_scale_;
        // The equilibrium product at the node's T (scaled): n_ie^2, or n_ie^2 gamma_n gamma_p
        // with its density partials (Fermi-Dirac), and d / d theta.
        physics::EquilibriumProduct E = physics::boltzmann_equilibrium_product(nie);
        double dE_t = 2.0 * nie * d_nie;
        if (fermi_dirac_) {
            const NodeLevels& l = g[i];
            E = physics::fermi_dirac_equilibrium_product(
                nie, {l.n.log_gamma, l.n.dg_dc}, {l.p.log_gamma, l.p.dg_dc});
            dE_t = E.value * (2.0 * d_nie / nie + l.n.dg_dt + l.p.dg_dt);
        }
        if (srh_) {
            // R = (n p - E) / den, den = tau_p (n + n_ie) + tau_n (p + n_ie):
            // d R / d theta = (-dE / den - (n p - E) (tau_p + tau_n) d_nie / den^2).
            const physics::RecombinationRate r =
                physics::srh_recombination(n, p, E, nie, tau_n_[i], tau_p_[i]);
            const double den = tau_p_[i] * (n + nie) + tau_n_[i] * (p + nie);
            const double d_t = -dE_t / den - (n * p - E.value) * (tau_p_[i] + tau_n_[i]) *
                                                  d_nie / (den * den);
            add(k, r.rate, r.d_dn, r.d_dp, d_t);
        }
        // Auger and radiative on physical densities (E Ns^2, its partials times Ns):
        // dR / dE = -(Cn n + Cp p) and -B.
        const double ns = n * Ns_, ps = p * Ns_;
        const physics::EquilibriumProduct Ep{E.value * Ns_ * Ns_, E.d_dn * Ns_, E.d_dp * Ns_};
        const double dEp = dE_t * Ns_ * Ns_;
        if (auger_) {
            const physics::RecombinationRate r =
                physics::auger_recombination(ns, ps, Ep, auger_n_[i], auger_p_[i]);
            const double d_t = -(auger_n_[i] * ns + auger_p_[i] * ps) * dEp;
            add(k, r.rate / Ns_, r.d_dn, r.d_dp, d_t / Ns_);
        }
        if (radiative_on_ && radiative_[i] > 0.0) {
            const physics::RecombinationRate r =
                physics::radiative_recombination(ns, ps, Ep, radiative_[i]);
            add(k, r.rate / Ns_, r.d_dn, r.d_dp, -radiative_[i] * dEp / Ns_);
        }
    }

    const double inv_K0 = 1.0 / heat_scale_;
    for (std::size_t k = 0; k < edges_.size(); ++k) {
        const EdgeTerm& e = edges_[k];
        const ThermalEdge& te = thermal_edges_[k];
        const std::size_t a = e.a, b = e.b;
        const std::size_t ends[2] = {a, b};
        // values position of (row of node `row`, r) with (column of node `col`, c).
        const auto pos = [&](std::size_t row, std::size_t r, std::size_t col, std::size_t c) {
            if (row == col) return block_[16 * row + 4 * r + c];
            return (row == a ? te.ab : te.ba)[4 * r + c];
        };
        // Row (node, r) gains sign * (value, gradient over the ends' (psi, n, p, theta)).
        const auto add = [&](std::size_t node, std::size_t r, double sign, double value,
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
        const Conduction q = conduction(k, t);
        double heat[2][4] = {{0, 0, 0, q.d_a}, {0, 0, 0, q.d_b}};
        double heat_value = q.flux;
        // Poisson, as the isothermal rows.
        const double poisson = e.c * (x[4 * b] - x[4 * a]);
        const double dpoisson[2][4] = {{-e.c, 0, 0, 0}, {e.c, 0, 0, 0}};
        if (!e.carriers) {
            if (contact_[a] < 0) add(a, 0, 1.0, poisson, dpoisson);
            if (contact_[b] < 0) add(b, 0, -1.0, poisson, dpoisson);
            if (!held(a)) add(a, 3, 1.0, heat_value, heat);
            if (!held(b)) add(b, 3, -1.0, heat_value, heat);
            continue;
        }
        const auto [fn, fp] = thermal_fluxes(k, x, t, g);
        const double gn[2][4] = {{fn.d_psi[0], fn.d_c[0], 0, fn.d_theta[0]},
                                 {fn.d_psi[1], fn.d_c[1], 0, fn.d_theta[1]}};
        const double gp[2][4] = {{fp.d_psi[0], 0, fp.d_c[0], fp.d_theta[0]},
                                 {fp.d_psi[1], 0, fp.d_c[1], fp.d_theta[1]}};
        for (const std::size_t node : ends) {
            if (contact_[node] >= 0) continue;
            const double sign = node == a ? 1.0 : -1.0;
            add(node, 0, sign, poisson, dpoisson);
            add(node, 1, sign, fn.flux, gn);
            add(node, 2, sign, fp.flux, gp);
        }
        // Energy: the mean carrier energies e_n = c_n - psi + (r_n + 5/2 + dh_n) theta and
        // e_p = c_p - psi + g - (r_p + 5/2 + dh_p) theta of the ends (dh the Fermi-Dirac
        // excess of the thermal diffusion factor, 0 under Boltzmann statistics).
        double en = 0.0, ep = 0.0, gen[2][4] = {}, gep[2][4] = {};
        for (std::size_t end = 0; end < 2; ++end) {
            const std::size_t i = ends[end];
            const ThermalNode& tn = thermal_nodes_[i];
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
            const double m = metal_[node];
            double grad[2][4];
            for (std::size_t end = 0; end < 2; ++end) {
                for (std::size_t c = 0; c < 4; ++c) {
                    grad[end][c] = heat[end][c] -
                                   inv_K0 * ((en + m) * gn[end][c] + fn.flux * gen[end][c] +
                                             (ep + m) * gp[end][c] + fp.flux * gep[end][c]);
                }
            }
            add(node, 3, sign,
                heat_value - inv_K0 * ((en + m) * fn.flux + (ep + m) * fp.flux), grad);
        }
    }
}

std::vector<double> DriftDiffusion::thermal_contact_heat(std::span<const double> x) const {
    NITCAD_EXPECTS(electrothermal() && x.size() == unknowns());
    std::vector<double> f(unknowns());
    assemble_heat(x, f, {}, nullptr, false);
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
    return heat;
}

double DriftDiffusion::electrical_power(std::span<const double> x) const {
    NITCAD_EXPECTS(electrothermal() && x.size() == unknowns());
    const std::vector<double> I = terminal_currents(x);
    double power = 0.0;
    for (std::size_t c = 0; c < I.size(); ++c) power += contact_bias_[c] * I[c];
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
    return joule;
}

}  // namespace NiTCAD::assemble
