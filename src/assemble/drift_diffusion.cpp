#include "NiTCAD/assemble/drift_diffusion.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <variant>

#include "NiTCAD/assemble/contact_bias.hpp"
#include "NiTCAD/assemble/ohmic.hpp"
#include "NiTCAD/assemble/sg_flux.hpp"
#include "NiTCAD/assemble/thermionic_flux.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/contract.hpp"
#include "NiTCAD/physics/band_to_band.hpp"
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

// The Moore-Penrose pseudo-inverse of a symmetric positive semi-definite D x D matrix (D <= 3,
// row-major), by Jacobi rotations: eigenvalues below 1e-12 of the largest are taken as zero, so a
// node whose edges do not span every direction reconstructs only the directions they span.
void pseudo_inverse(const double* a, int D, double* inverse) {
    double m[3][3] = {}, v[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    for (int r = 0; r < D; ++r) {
        for (int c = 0; c < D; ++c) m[r][c] = a[r * D + c];
    }
    for (int sweep = 0; sweep < 50; ++sweep) {
        double off = 0.0;
        for (int p = 0; p < D; ++p) {
            for (int q = p + 1; q < D; ++q) off += m[p][q] * m[p][q];
        }
        if (off == 0.0) break;
        for (int p = 0; p < D; ++p) {
            for (int q = p + 1; q < D; ++q) {
                if (m[p][q] == 0.0) continue;
                const double theta = 0.5 * std::atan2(2.0 * m[p][q], m[q][q] - m[p][p]);
                const double c = std::cos(theta), s = std::sin(theta);
                for (int k = 0; k < D; ++k) {  // m <- m J
                    const double mkp = m[k][p], mkq = m[k][q];
                    m[k][p] = c * mkp - s * mkq;
                    m[k][q] = s * mkp + c * mkq;
                }
                for (int k = 0; k < D; ++k) {  // m <- J^T m
                    const double mpk = m[p][k], mqk = m[q][k];
                    m[p][k] = c * mpk - s * mqk;
                    m[q][k] = s * mpk + c * mqk;
                }
                for (int k = 0; k < D; ++k) {  // v <- v J
                    const double vkp = v[k][p], vkq = v[k][q];
                    v[k][p] = c * vkp - s * vkq;
                    v[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }
    double largest = 0.0;
    for (int k = 0; k < D; ++k) largest = std::max(largest, std::abs(m[k][k]));
    for (int r = 0; r < D; ++r) {
        for (int c = 0; c < D; ++c) {
            double sum = 0.0;
            for (int k = 0; k < D; ++k) {
                if (m[k][k] > 1e-12 * largest) sum += v[r][k] * v[c][k] / m[k][k];
            }
            inverse[r * D + c] = sum;
        }
    }
}

}  // namespace

std::expected<DriftDiffusion, base::Error> DriftDiffusion::create(
    const device::Device& device, const Scaling& scaling, const PhysicsModels& models) {
    if (models.impact_ionization && !(std::isfinite(models.impact_current_resolution) &&
                                      models.impact_current_resolution >= 0.0)) {
        return std::unexpected(base::Error{
            base::ErrorCode::invalid_input,
            "impact_current_resolution must be finite and not negative",
            base::ErrorContext{.index = std::nullopt, .value = models.impact_current_resolution}});
    }
    const bool nonlocal = models.btbt_nonlocal != NonlocalTunnelling::off;
    const auto refuse = [](const char* message) {
        return std::unexpected(base::Error{base::ErrorCode::invalid_input, message, {}});
    };
    if (nonlocal && models.btbt_local) {
        return refuse("btbt_local and btbt_nonlocal would count the same tunnelling twice");
    }
    if (nonlocal && device.cells() == nullptr) {
        return refuse("nonlocal tunnelling needs the device's tensor cells");
    }
    auto scaled = detail::make_scaled_device(device, scaling, models);
    if (!scaled) return std::unexpected(std::move(scaled.error()));
    const std::size_t n = scaled->volume.size();
    const double T = scaling.temperature_K;
    if (nonlocal) {
        for (const detail::ScaledEdge& e : scaled->edges) {
            if (e.carriers && e.interface) {
                return refuse("nonlocal tunnelling does not cross a heterointerface");
            }
        }
    }
    if (models.btbt_nonlocal == NonlocalTunnelling::direct_wkb) {
        for (const device::Region& r : device.regions()) {
            if (device::is_insulator(r)) continue;
            const physics::BandToBandParameters& b =
                std::get<physics::Semiconductor>(r.material).parameters().band_to_band;
            if (!b.direct_gap || physics::tunnelling_reduced_mass(b) == 0.0) {
                return refuse("direct_wkb needs a direct gap and cited tunnelling masses in every "
                              "semiconductor region (never silicon)");
            }
        }
    }

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

    // With the electrothermal model every node and edge block is full 4 x 4.
    s.m_ = models.electrothermal ? 4 : 3;
    const std::size_t stride = s.m_;
    std::vector<linalg::Triplet> triplets;
    triplets.reserve(stride * stride * n + 2 * stride * stride * scaled->edges.size());
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t r = 0; r < stride; ++r) {
            for (std::size_t c = 0; c < stride; ++c) {
                triplets.push_back({static_cast<linalg::Index>(stride * i + r),
                                    static_cast<linalg::Index>(stride * i + c), 0.0});
            }
        }
    }
    for (const detail::ScaledEdge& e : scaled->edges) {
        if (models.electrothermal) {
            for (std::size_t r = 0; r < 4; ++r) {
                for (std::size_t c = 0; c < 4; ++c) {
                    triplets.push_back({static_cast<linalg::Index>(4 * e.i + r),
                                        static_cast<linalg::Index>(4 * e.j + c), 0.0});
                    triplets.push_back({static_cast<linalg::Index>(4 * e.j + r),
                                        static_cast<linalg::Index>(4 * e.i + c), 0.0});
                }
            }
            continue;
        }
        for (std::size_t k = 0; k < 5; ++k) {
            triplets.push_back({static_cast<linalg::Index>(stride * e.i + edge_rows[k]),
                                static_cast<linalg::Index>(stride * e.j + edge_cols[k]), 0.0});
            triplets.push_back({static_cast<linalg::Index>(stride * e.j + edge_rows[k]),
                                static_cast<linalg::Index>(stride * e.i + edge_cols[k]), 0.0});
        }
    }
    // An interface edge: the insulator's Poisson row reads the semiconductor's n and p.
    for (const InterfaceEdge& f : s.interfaces_.edges()) {
        for (const std::size_t c : {std::size_t{1}, std::size_t{2}}) {
            triplets.push_back({static_cast<linalg::Index>(stride * f.insulator),
                                static_cast<linalg::Index>(stride * f.semiconductor + c), 0.0});
        }
    }
    // Impact ionization: a node's electron generation reads its neighbours' p (through the hole
    // current) and its hole generation their n.
    if (models.impact_ionization) {
        for (const detail::ScaledEdge& e : scaled->edges) {
            if (!e.carriers) continue;
            for (const auto [a, b] : {std::pair{e.i, e.j}, std::pair{e.j, e.i}}) {
                triplets.push_back({static_cast<linalg::Index>(stride * a + 1),
                                    static_cast<linalg::Index>(stride * b + 2), 0.0});
                triplets.push_back({static_cast<linalg::Index>(stride * a + 2),
                                    static_cast<linalg::Index>(stride * b + 1), 0.0});
            }
        }
    }
    const auto size = static_cast<linalg::Index>(stride * n);
    auto pattern = linalg::SparseMatrix::from_triplets(size, size, triplets);
    if (!pattern) return std::unexpected(std::move(pattern.error()));
    s.pattern_ = std::move(*pattern);
    if (nonlocal) s.base_triplets_ = std::move(triplets);  // set_paths adds the paths' entries

    s.block_.resize(stride * stride * n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t r = 0; r < stride; ++r) {
            for (std::size_t c = 0; c < stride; ++c) {
                s.block_[stride * stride * i + stride * r + c] =
                    detail::position(s.pattern_, stride * i + r, stride * i + c);
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
            t.ab[k] = detail::position(s.pattern_, stride * e.i + edge_rows[k],
                                       stride * e.j + edge_cols[k]);
            t.ba[k] = detail::position(s.pattern_, stride * e.j + edge_rows[k],
                                       stride * e.i + edge_cols[k]);
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

    if (models.impact_ionization) {
        // Per node with coefficients: its carrier edges and the pseudo-inverse of sum t t^T.
        const mesh::Mesh& m = device.mesh();
        const int D = m.dimension();
        s.dimension_ = D;
        s.generation_scale_ = 1.0 / scaling.R0;
        s.impact_resolution_ = models.impact_current_resolution;
        std::vector<std::vector<std::size_t>> incident(n);
        for (std::size_t k = 0; k < scaled->edges.size(); ++k) {
            const detail::ScaledEdge& e = scaled->edges[k];
            if (!e.carriers) continue;
            incident[e.i].push_back(k);
            incident[e.j].push_back(k);
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (s.insulator_[i] != 0 || s.contact_[i] >= 0 || incident[i].empty()) continue;
            const physics::Semiconductor& mat = device.material(static_cast<mesh::NodeId>(i));
            const physics::ImpactIonizationParameters& ii = mat.parameters().impact_ionization;
            if (ii.electron.A_low_per_cm == 0.0 && ii.electron.A_high_per_cm == 0.0 &&
                ii.hole.A_low_per_cm == 0.0 && ii.hole.A_high_per_cm == 0.0) {
                continue;
            }
            ImpactNode node{i, s.impact_edges_.size(), 0, {}, ii.electron, ii.hole,
                            physics::impact_ionization_temperature_factor(ii.phonon_energy_eV, T)};
            double outer[9] = {};
            for (const std::size_t k : incident[i]) {
                const detail::ScaledEdge& e = scaled->edges[k];
                const mesh::Point& pa = m.points()[e.i];
                const mesh::Point& pb = m.points()[e.j];
                ImpactEdge ie{k, {}, scaling.V_T / e.length_cm,
                              scaling.J0 / (e.geometry * e.length_cm / scaling.L_D), 0, 0};
                for (int d = 0; d < D; ++d) ie.t[d] = (pb[d] - pa[d]) / e.length_cm;
                for (int r = 0; r < D; ++r) {
                    for (int c = 0; c < D; ++c) outer[r * D + c] += ie.t[r] * ie.t[c];
                }
                const std::size_t other = e.i == i ? e.j : e.i;
                ie.n_p = detail::position(s.pattern_, stride * i + 1, stride * other + 2);
                ie.p_n = detail::position(s.pattern_, stride * i + 2, stride * other + 1);
                s.impact_edges_.push_back(ie);
            }
            pseudo_inverse(outer, D, node.inverse);
            node.last = s.impact_edges_.size();
            s.impact_nodes_.push_back(node);
        }
    }
    if (models.btbt_local) {
        // Per node with a Kane pair: its carrier edges and the pseudo-inverse of sum t t^T.
        const mesh::Mesh& m = device.mesh();
        const int D = m.dimension();
        s.dimension_ = D;
        s.generation_scale_ = 1.0 / scaling.R0;
        std::vector<std::vector<std::size_t>> incident(n);
        for (std::size_t k = 0; k < scaled->edges.size(); ++k) {
            const detail::ScaledEdge& e = scaled->edges[k];
            if (!e.carriers) continue;
            incident[e.i].push_back(k);
            incident[e.j].push_back(k);
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (s.insulator_[i] != 0 || s.contact_[i] >= 0 || incident[i].empty()) continue;
            const physics::BandToBandParameters& b =
                device.material(static_cast<mesh::NodeId>(i)).parameters().band_to_band;
            if (b.A_per_cm3_s == 0.0) continue;
            TunnelNode node{i, s.field_edges_.size(), 0, {}, b.A_per_cm3_s, b.B_V_per_cm};
            double outer[9] = {};
            for (const std::size_t k : incident[i]) {
                const detail::ScaledEdge& e = scaled->edges[k];
                const mesh::Point& pa = m.points()[e.i];
                const mesh::Point& pb = m.points()[e.j];
                FieldEdge fe{k, {}, scaling.V_T / e.length_cm};
                for (int d = 0; d < D; ++d) fe.t[d] = (pb[d] - pa[d]) / e.length_cm;
                for (int r = 0; r < D; ++r) {
                    for (int c = 0; c < D; ++c) outer[r * D + c] += fe.t[r] * fe.t[c];
                }
                s.field_edges_.push_back(fe);
            }
            pseudo_inverse(outer, D, node.inverse);
            node.last = s.field_edges_.size();
            s.tunnel_nodes_.push_back(node);
        }
    }
    if (nonlocal) {
        // Per node: the longest path that may start there, where the rate falls below exp(-690)
        // of its scale: E_g / F_min with F_min = B / 690 (kane), or the field at which eq. (8)'s
        // exponent reaches 690 (direct_wkb).
        s.tunnel_kind_ = models.btbt_nonlocal;
        s.cells_ = *device.cells();
        s.generation_scale_ = 1.0 / scaling.R0;
        for (auto* v : {&s.tunnel_length_, &s.tunnel_gap_, &s.tunnel_mc_, &s.tunnel_mv_,
                        &s.tunnel_A_, &s.tunnel_B_}) {
            v->assign(n, 0.0);
        }
        constexpr double exponent = 690.0;
        for (std::size_t i = 0; i < n; ++i) {
            if (s.insulator_[i] != 0) continue;
            const physics::Semiconductor& mat = device.material(static_cast<mesh::NodeId>(i));
            const physics::BandToBandParameters& b = mat.parameters().band_to_band;
            const double Eg = physics::band_gap_eV(mat, T);
            s.tunnel_gap_[i] = Eg;
            if (s.tunnel_kind_ == NonlocalTunnelling::kane) {
                if (b.A_per_cm3_s == 0.0) continue;
                s.tunnel_A_[i] = b.A_per_cm3_s;
                s.tunnel_B_[i] = b.B_V_per_cm;
                s.tunnel_length_[i] = Eg / (b.B_V_per_cm / exponent);
                continue;
            }
            s.tunnel_mc_[i] = b.electron_mass * base::m0_kg;
            s.tunnel_mv_[i] = b.hole_mass * base::m0_kg;
            const double mr = physics::tunnelling_reduced_mass(b) * base::m0_kg;
            const double Eg_J = Eg * base::q_C;
            const double F_min = std::numbers::pi * std::sqrt(mr) * std::pow(Eg_J, 1.5) /
                                 (2.0 * base::q_C * base::hbar_J_s * exponent);  // [V/m]
            s.tunnel_length_[i] = Eg_J / (base::q_C * F_min) * 100.0;           // [cm]
        }
    }
    s.generation_ = !s.impact_nodes_.empty() || models.btbt_local || nonlocal;
    for (const InterfaceEdge& f : s.interfaces_.edges()) {
        s.interface_np_.emplace_back(
            detail::position(s.pattern_, stride * f.insulator, stride * f.semiconductor + 1),
            detail::position(s.pattern_, stride * f.insulator, stride * f.semiconductor + 2));
    }
    s.psi0_.assign(n, 0.0);
    s.n0_.assign(n, 0.0);
    s.p0_.assign(n, 0.0);
    if (models.electrothermal) {
        std::vector<double> geometry;
        geometry.reserve(scaled->edges.size());
        for (const detail::ScaledEdge& e : scaled->edges) geometry.push_back(e.geometry);
        if (auto e = s.make_thermal(device, scaling, models, geometry)) {
            return std::unexpected(std::move(*e));
        }
    }
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
        if (electrothermal()) {
            metal_[i] = bias;
            contact_bias_[static_cast<std::size_t>(contact_[i])] = bias;
        }
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
        x[m_ * i] = psi[i];
        if (insulator_[i] != 0) continue;  // stamped below
        x[m_ * i + 1] = detail::density(fermi_dirac_, n_ie_[i], log_dos_n_[i], eta).density;
        x[m_ * i + 2] = detail::density(fermi_dirac_, n_ie_[i], log_dos_p_[i], -eta).density;
        if (electrothermal()) x[4 * i + 3] = 0.0;
    }
    stamp_contacts(x);
    return x;
}

void DriftDiffusion::stamp_contacts(std::span<double> x) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    if (electrothermal()) {
        // The isothermal sinks first: an ohmic node's values depend on its temperature.
        for (std::size_t i = 0; i < node_count(); ++i) {
            const std::int32_t c = thermal_nodes_[i].contact;
            if (c >= 0 && thermal_isothermal_[static_cast<std::size_t>(c)] != 0) {
                x[4 * i + 3] = thermal_rise_[static_cast<std::size_t>(c)];
            }
            if (insulator_[i] != 0) {
                x[4 * i + 1] = x[4 * i + 2] = 0.0;
            } else if (contact_[i] >= 0) {
                const ThermalOhmic o = thermal_ohmic(i, thermal_state(i, x[4 * i + 3]));
                x[4 * i] = o.psi;
                x[4 * i + 1] = o.n;
                x[4 * i + 2] = o.p;
            }
        }
        return;
    }
    for (std::size_t i = 0; i < node_count(); ++i) {
        if (insulator_[i] != 0) {
            if (electrode_[i] >= 0) x[m_ * i] = psi0_[i];
            x[m_ * i + 1] = x[m_ * i + 2] = 0.0;
            continue;
        }
        if (contact_[i] < 0) continue;
        x[m_ * i] = psi0_[i];
        x[m_ * i + 1] = n0_[i];
        x[m_ * i + 2] = p0_[i];
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
        g.push_back({physics::fermi_dirac_degeneracy(n_ie_[i], log_dos_n_[i], x[m_ * i + 1]),
                     physics::fermi_dirac_degeneracy(n_ie_[i], log_dos_p_[i], x[m_ * i + 2])});
    }
    return g;
}

std::pair<EdgeFlux, EdgeFlux> DriftDiffusion::edge_fluxes(
    const EdgeTerm& e, std::span<const double> x, std::span<const NodeDegeneracy> g) const {
    // The flux functions take the driving term as psi2 - psi1; its partials are those with
    // respect to psi_a and psi_b, since ln(n_ie) and s do not depend on the state.
    const double dpsi = x[m_ * e.b] - x[m_ * e.a];
    const double na = x[m_ * e.a + 1], nb = x[m_ * e.b + 1];
    const double pa = x[m_ * e.a + 2], pb = x[m_ * e.b + 2];
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
    const double n = x[m_ * i + 1], p = x[m_ * i + 2];
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
    return interfaces_.drift(v, x[m_ * o], x[m_ * s], x[m_ * s + 1], x[m_ * s + 2], statistics(s),
                             t ? &*t : nullptr);
}

void DriftDiffusion::assemble(std::span<const double> x, std::span<double> f,
                              std::span<double> values, const TimeStep* step,
                              bool interfaces) const {
    if (electrothermal()) {
        assemble_heat(x, f, values, step, true);
        return;
    }
    NITCAD_EXPECTS(x.size() == unknowns() && f.size() == unknowns());
    NITCAD_EXPECTS(step == nullptr || (step->storage.size() == 2 * node_count() &&
                                       step->traps.size() == trap_slots()));
    const bool jacobian = !values.empty();
    const auto at = [&](std::size_t node, std::size_t r, std::size_t c) -> double& {
        return values[block_[m_ * m_ * node + m_ * r + c]];
    };
    const std::vector<NodeDegeneracy> g = degeneracies(x);
    // The equilibrium product of node i, scaled.
    const auto product = [&](std::size_t i) {
        return fermi_dirac_ ? physics::fermi_dirac_equilibrium_product(n_ie_[i], g[i].n, g[i].p)
                            : physics::boltzmann_equilibrium_product(n_ie_[i]);
    };
    for (std::size_t i = 0; i < node_count(); ++i) {
        const double psi = x[m_ * i], n = x[m_ * i + 1], p = x[m_ * i + 2];
        if (contact_[i] >= 0) {
            f[m_ * i] = psi - psi0_[i];
            f[m_ * i + 1] = n - n0_[i];
            f[m_ * i + 2] = p - p0_[i];
            if (jacobian) at(i, 0, 0) = at(i, 1, 1) = at(i, 2, 2) = 1.0;
            continue;
        }
        if (insulator_[i] != 0) {  // no charge, no carriers; an electrode holds psi
            const bool electrode = electrode_[i] >= 0;
            f[m_ * i] = electrode ? psi - psi0_[i] : 0.0;
            f[m_ * i + 1] = n;
            f[m_ * i + 2] = p;
            if (jacobian) {
                if (electrode) at(i, 0, 0) = 1.0;
                at(i, 1, 1) = at(i, 2, 2) = 1.0;
            }
            continue;
        }
        const double V = volume_[i];
        f[m_ * i] = -V * (n - p - doping_[i]);
        f[m_ * i + 1] = 0.0;
        f[m_ * i + 2] = 0.0;
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
            f[m_ * i] = -V * (n - p - c.value);
            if (jacobian) {
                at(i, 0, 1) = -V * (1.0 - c.d_eta_c * (1.0 / n - wn));
                at(i, 0, 2) = V * (1.0 + c.d_eta_v * (1.0 / p - wp));
            }
        }
        if (gates_.on_gate(i)) {
            f[m_ * i] += gates_.row_term(i, psi);
            if (jacobian) at(i, 0, 0) -= gates_.term(i).coupling;
        }
        if (srh_) {
            const physics::RecombinationRate r =
                physics::srh_recombination(n, p, product(i), n_ie_[i], tau_n_[i], tau_p_[i]);
            const double k = V * rate_scale_;
            f[m_ * i + 1] -= k * r.rate;
            f[m_ * i + 2] += k * r.rate;
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
            f[m_ * i + 1] -= k * (r.rate / Ns_);
            f[m_ * i + 2] += k * (r.rate / Ns_);
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
            f[m_ * i + 1] -= k * (st.n - step->storage[2 * i]);
            f[m_ * i + 2] += k * (st.p - step->storage[2 * i + 1]);
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
            f[m_ * o] += t.flux_insulator;
            if (jacobian) {
                at(o, 0, 0) += t.d_flux_insulator[0];
                values[os] += t.d_flux_insulator[1];
                values[interface_np_[k].first] += t.d_flux_insulator[2];
                values[interface_np_[k].second] += t.d_flux_insulator[3];
            }
        }
        if (contact_[s] >= 0) continue;  // a Dirichlet semiconductor node
        f[m_ * s] += t.flux_semiconductor;
        f[m_ * s + 1] -= t.rate;
        f[m_ * s + 2] += t.rate_p;
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
        const double psi_a = x[m_ * a], psi_b = x[m_ * b];
        const double poisson = e.c * (psi_b - psi_a);
        if (e.charged) continue;  // the interface's half-edge fluxes, above
        if (!e.carriers) {  // Poisson only, into the rows that are not Dirichlet
            if (contact_[a] < 0 && electrode_[a] < 0) {
                f[m_ * a] += poisson;
                if (jacobian) {
                    at(a, 0, 0) -= e.c;
                    values[e.ab[0]] += e.c;
                }
            }
            if (contact_[b] < 0 && electrode_[b] < 0) {
                f[m_ * b] -= poisson;
                if (jacobian) {
                    at(b, 0, 0) -= e.c;
                    values[e.ba[0]] += e.c;
                }
            }
            continue;
        }
        const auto [fn, fp] = edge_fluxes(e, x, g);
        if (contact_[a] < 0) {
            f[m_ * a] += poisson;
            f[m_ * a + 1] += fn.flux;
            f[m_ * a + 2] += fp.flux;
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
            f[m_ * b] -= poisson;
            f[m_ * b + 1] -= fn.flux;
            f[m_ * b + 2] -= fp.flux;
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
    if (!impact_nodes_.empty()) add_impact_generation(x, g, f, values);
    if (!tunnel_nodes_.empty()) add_local_tunnelling(x, f, values);
    if (!paths_.paths.empty()) add_path_tunnelling(x, f, values);
}

void DriftDiffusion::add_impact_generation(std::span<const double> x,
                                           std::span<const NodeDegeneracy> g,
                                           std::span<double> f,
                                           std::span<double> values) const {
    const bool jacobian = !values.empty();
    const int D = dimension_;
    // The carrier edges' fluxes, once.
    std::vector<std::pair<EdgeFlux, EdgeFlux>> flux(edges_.size());
    for (std::size_t k = 0; k < edges_.size(); ++k) {
        if (edges_[k].carriers) flux[k] = edge_fluxes(edges_[k], x, g);
    }
    constexpr double q = base::q_C;
    constexpr double floor2 = impact_current_floor * impact_current_floor;
    const double rho2 = impact_resolution_ * impact_resolution_;
    for (const ImpactNode& node : impact_nodes_) {
        const std::size_t i = node.node;
        const double count = static_cast<double>(node.last - node.first);
        // The field and the two current densities at the node, reconstructed from their edge
        // components by least squares: v = P sum_k t_k v_k, P the pseudo-inverse of sum t t^T.
        // Edge components: E_k = -field (psi_b - psi_a), j_k = current * flux.
        // And the size of the opposing flux terms whose difference each current is, averaged over
        // the edges: R_c = mean_k current_k a_k (c_a + c_b) sqrt(1 + Delta_k^2), a the edge's
        // low-field factor, c the carrier's density, Delta = psi_b - psi_a.
        double E[3] = {}, jn[3] = {}, jp[3] = {}, R[2] = {0.0, 0.0};
        for (std::size_t k = node.first; k < node.last; ++k) {
            const ImpactEdge& ie = impact_edges_[k];
            const EdgeTerm& e = edges_[ie.edge];
            const double delta = x[m_ * e.b] - x[m_ * e.a];
            const double Ek = -ie.field * delta;
            const double jnk = ie.current * flux[ie.edge].first.flux;
            const double jpk = ie.current * flux[ie.edge].second.flux;
            const double root = std::sqrt(1.0 + delta * delta);
            R[0] += ie.current * e.an * (x[m_ * e.a + 1] + x[m_ * e.b + 1]) * root / count;
            R[1] += ie.current * e.ap * (x[m_ * e.a + 2] + x[m_ * e.b + 2]) * root / count;
            for (int r = 0; r < D; ++r) {
                double w = 0.0;  // (P t_k)_r
                for (int c = 0; c < D; ++c) w += node.inverse[r * D + c] * ie.t[c];
                E[r] += w * Ek;
                jn[r] += w * jnk;
                jp[r] += w * jpk;
            }
        }
        // Per carrier, with m = sqrt(|j|^2 + eps^2), eps^2 = floor^2 + (rho R)^2:
        //     E_par = |E . j| / m,   G = alpha(E_par) |j|^2 / m / q.
        // E_par is the field along the current (|E . j^| for |j| >> eps), at most |E|; |j|^2 / m is
        // |j| for |j| >> eps, at most |j|, and smooth through j = 0, where it and its gradient
        // vanish. eps is the current's resolution at the node: rho =
        // models.impact_current_resolution of the opposing flux terms, above impact_current_floor;
        // a current at or below it (rounding, in the dense spill-over layer at a junction, where it
        // is a difference of fluxes 1e14 times larger and flips sign between Newton iterates) does
        // not ionize and has no kink.
        // Partials, s = sign(E . j), w = |j|^2 / m:
        //     dG/dE   = alpha' w s j / m / q,
        //     dG/dj   = (alpha' w (s E / m - |E . j| j / m^3)
        //                + alpha j (|j|^2 + 2 eps^2) / m^3) / q,
        //     dG/deps = -(alpha' w |E . j| + alpha |j|^2) eps / m^3 / q,  deps/dR = rho^2 R / eps.
        double G = 0.0, dE[3] = {}, dj[2][3] = {}, dR[2] = {0.0, 0.0};
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
                c == 0 ? node.n : node.p, node.gamma, Epar);
            if (a.alpha == 0.0 && a.d_dE == 0.0) continue;
            const double w = mag2 / m;
            G += a.alpha * w / q;
            const double s = Ej > 0.0 ? 1.0 : Ej < 0.0 ? -1.0 : 0.0;
            const double m3 = m * m * m;
            for (int r = 0; r < D; ++r) {
                dE[r] += a.d_dE * w * s * j[r] / m / q;
                dj[c][r] = (a.d_dE * w * (s * E[r] / m - std::abs(Ej) * j[r] / m3) +
                            a.alpha * j[r] * (mag2 + 2.0 * eps2) / m3) /
                           q;
            }
            // dG/deps deps/dR = -(alpha' w |E . j| + alpha |j|^2) / m^3 / q * rho^2 R.
            dR[c] = -(a.d_dE * w * std::abs(Ej) + a.alpha * mag2) / m3 / q * rho2 * R[c];
        }
        const double k_row = volume_[i] * generation_scale_;
        f[m_ * i + 1] += k_row * G;
        f[m_ * i + 2] -= k_row * G;
        if (!jacobian) continue;
        const auto add = [&](std::size_t node_j, int comp, std::size_t k, double v) {
            // Row n of node i gets +v, row p -v, at column comp of node_j through edge k.
            const ImpactEdge& ie = impact_edges_[k];
            const EdgeTerm& e = edges_[ie.edge];
            const bool own = node_j == i;
            const bool i_is_a = e.a == i;
            for (const int row : {1, 2}) {
                const double sign = row == 1 ? k_row : -k_row;
                std::size_t pos;
                if (own) {
                    pos = block_[m_ * m_ * i + m_ * static_cast<std::size_t>(row) +
                                 static_cast<std::size_t>(comp)];
                } else {
                    const std::size_t* p = i_is_a ? e.ab : e.ba;
                    // Edge positions: (psi,psi), (n,psi), (n,n), (p,psi), (p,p).
                    if (row == 1) {
                        pos = comp == 0 ? p[1] : comp == 1 ? p[2] : ie.n_p;
                    } else {
                        pos = comp == 0 ? p[3] : comp == 1 ? ie.p_n : p[4];
                    }
                }
                values[pos] += sign * v;
            }
        };
        for (std::size_t k = node.first; k < node.last; ++k) {
            const ImpactEdge& ie = impact_edges_[k];
            const EdgeTerm& e = edges_[ie.edge];
            const EdgeFlux& fn = flux[ie.edge].first;
            const EdgeFlux& fp = flux[ie.edge].second;
            double wE = 0.0, wn = 0.0, wp = 0.0;  // dG / d(edge component)
            for (int r = 0; r < D; ++r) {
                double w = 0.0;
                for (int c = 0; c < D; ++c) w += node.inverse[r * D + c] * ie.t[c];
                wE += w * dE[r];
                wn += w * dj[0][r];
                wp += w * dj[1][r];
            }
            wn *= ie.current;
            wp *= ie.current;
            // R's own partials through this edge: d/dc_a = d/dc_b = current a root / count,
            // d/dpsi_b = -d/dpsi_a = current a (c_a + c_b) Delta / root / count.
            const double delta = x[m_ * e.b] - x[m_ * e.a];
            const double root = std::sqrt(1.0 + delta * delta);
            const double rn = ie.current * e.an / count, rp = ie.current * e.ap / count;
            const double sn = x[m_ * e.a + 1] + x[m_ * e.b + 1];
            const double sp = x[m_ * e.a + 2] + x[m_ * e.b + 2];
            const double dpsi = (dR[0] * rn * sn + dR[1] * rp * sp) * delta / root;
            // E_k = -field (psi_b - psi_a); the fluxes' own partials.
            add(e.a, 0, k, wE * ie.field + wn * fn.d_psi1 + wp * fp.d_psi1 - dpsi);
            add(e.b, 0, k, -wE * ie.field + wn * fn.d_psi2 + wp * fp.d_psi2 + dpsi);
            add(e.a, 1, k, wn * fn.d_c1 + dR[0] * rn * root);
            add(e.b, 1, k, wn * fn.d_c2 + dR[0] * rn * root);
            add(e.a, 2, k, wp * fp.d_c1 + dR[1] * rp * root);
            add(e.b, 2, k, wp * fp.d_c2 + dR[1] * rp * root);
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
    const std::size_t w = storage_width();
    std::vector<double> s(w * node_count(), 0.0);
    for (std::size_t i = 0; i < node_count(); ++i) {
        if (electrothermal()) s[3 * i + 2] = x[4 * i + 3];  // tau, every node
        if (insulator_[i] != 0) continue;
        const NodeStorage st = node_storage(i, x, g);
        s[w * i] = st.n;
        s[w * i + 1] = st.p;
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
            v, x[m_ * o], x[m_ * s], x[m_ * s + 1], x[m_ * s + 2], statistics(s), t ? &*t : nullptr,
            std::span(f).subspan(interfaces_.slot_offset(k), v.last_level - v.first_level));
    }
    return f;
}

std::vector<double> DriftDiffusion::contact_charges(std::span<const double> x,
                                                    const TimeStep* step) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    std::vector<double> charge = gates_.charges(x, m_, contact_count_);
    const auto owner = [&](std::size_t i) {
        return contact_[i] >= 0 ? contact_[i] : electrode_[i];
    };
    for (const EdgeTerm& e : edges_) {
        const std::int32_t ca = owner(e.a), cb = owner(e.b);
        if (ca == cb || e.charged) continue;
        const double flux = e.c * (x[m_ * e.a] - x[m_ * e.b]);  // from a to b
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
        largest = std::max({largest, x[m_ * i + 1], x[m_ * i + 2]});
    }
    // With impact ionization (Unit 19) a sub-floor density's row also carries the generation of
    // the currents around it, whose rounding (the opposing fluxes of the dense spill-over layers
    // at a junction) limits it far above eps. Measured on the legacy one-sided junction: with
    // the floor at 1e-20 of the largest density Newton chatters in a period two at 1.5e-5 of a
    // density 1e-21 of the largest (-0.5 V); at 1e-12 the trace's corrector stalls at 1.4e-7 of
    // that floor (-9.8 V). Densities below 1e-8 of the largest, the legacy's floor for its
    // impact-ionization paths, are measured against it then.
    const double floor = (generation_ ? 1e-8 : 1e-20) * largest;
    double size = 0.0;
    for (std::size_t i = 0; i < node_count(); ++i) {
        if (insulator_[i] != 0) {
            size = std::max(size, std::abs(dx[m_ * i]));
            continue;
        }
        const double n = x[m_ * i + 1], p = x[m_ * i + 2];
        size = std::max({size, std::abs(dx[m_ * i]),
                         n > 0.0 ? std::abs(dx[m_ * i + 1]) / std::max(n, floor) : infinity,
                         p > 0.0 ? std::abs(dx[m_ * i + 2]) / std::max(p, floor) : infinity});
    }
    if (electrothermal()) {  // theta, an absolute bound in units of T0 (T10)
        for (std::size_t i = 0; i < node_count(); ++i) {
            size = std::max(size, std::abs(dx[4 * i + 3]));
        }
    }
    return size;
}

void DriftDiffusion::apply_update(std::span<double> x, std::span<const double> dx,
                                  double max_update) const {
    NITCAD_EXPECTS(x.size() == unknowns() && dx.size() == unknowns());
    for (std::size_t i = 0; i < node_count(); ++i) {
        x[m_ * i] += std::clamp(dx[m_ * i], -max_update, max_update);
        if (electrothermal()) {  // 50 K (legacy thermal.py) and a factor of 2 at most
            // theta = 1 + tau within [theta / 2, 2 theta]: tau within [(tau - 1) / 2, 1 + 2 tau].
            const double tau = x[4 * i + 3], cap = 50.0 / T0_;
            x[4 * i + 3] = std::clamp(tau + std::clamp(dx[4 * i + 3], -cap, cap),
                                      0.5 * (tau - 1.0), 1.0 + 2.0 * tau);
        }
        if (insulator_[i] != 0) continue;  // the densities stay 0
        for (const std::size_t k : {m_ * i + 1, m_ * i + 2}) {
            x[k] = std::clamp(x[k] + dx[k], 0.1 * x[k], 10.0 * x[k]);
        }
    }
}

std::vector<std::pair<double, double>> DriftDiffusion::edge_currents(
    std::span<const double> x) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    std::vector<std::pair<double, double>> currents;
    currents.reserve(edges_.size());
    if (electrothermal()) {
        const std::vector<ThermalState> t = thermal_states(x);
        const std::vector<NodeLevels> g = thermal_levels(x, t);
        for (std::size_t k = 0; k < edges_.size(); ++k) {
            if (!edges_[k].carriers) {
                currents.emplace_back(0.0, 0.0);
                continue;
            }
            const auto [fn, fp] = thermal_fluxes(k, x, t, g);
            currents.emplace_back(fn.flux, fp.flux);
        }
        return currents;
    }
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
    std::vector<ThermalState> states;
    std::vector<NodeLevels> levels;
    if (electrothermal()) {
        states = thermal_states(x);
        levels = thermal_levels(x, states);
    }
    std::vector<double> terms(n, 0.0), flux_n(n, 0.0);
    for (std::size_t k = 0; k < edges_.size(); ++k) {
        const EdgeTerm& e = edges_[k];
        if (!e.carriers) continue;
        double t, jn;
        if (electrothermal()) {
            const auto [fn, fp] = thermal_fluxes(k, x, states, levels);
            t = std::abs(fn.d_c[0] * x[4 * e.a + 1]) + std::abs(fn.d_c[1] * x[4 * e.b + 1]) +
                std::abs(fp.d_c[0] * x[4 * e.a + 2]) + std::abs(fp.d_c[1] * x[4 * e.b + 2]);
            jn = fn.flux;
        } else {
            const auto [fn, fp] = edge_fluxes(e, x, g);
            t = std::abs(fn.d_c1 * x[m_ * e.a + 1]) + std::abs(fn.d_c2 * x[m_ * e.b + 1]) +
                std::abs(fp.d_c1 * x[m_ * e.a + 2]) + std::abs(fp.d_c2 * x[m_ * e.b + 2]);
            jn = fn.flux;
        }
        terms[e.a] += t;
        terms[e.b] += t;
        flux_n[e.a] += jn;
        flux_n[e.b] -= jn;
    }
    const double eps8 = 8.0 * std::numeric_limits<double>::epsilon();
    double bound = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        if (contact_[i] >= 0) {
            bound += eps8 * terms[i];  // the terminal current's own sum
            continue;
        }
        const double recombination = std::abs(flux_n[i] - f[m_ * i + 1]);  // V R
        bound += std::abs(f[m_ * i + 1] + f[m_ * i + 2]) + eps8 * (terms[i] + 2.0 * recombination);
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
        if (electrothermal()) {
            // E_c = c_n - psi, E_v = c_p - psi + g; E_Fn = E_c + theta x_n, E_Fp = E_v - theta x_p
            // with x the bands' reduced energies at the node's T (level; units of k T0).
            const ThermalNode& tn = thermal_nodes_[i];
            const double psi = x[4 * i];
            const ThermalState s = thermal_state(i, x[4 * i + 3]);
            const double theta = s.theta;
            b.conduction[i] = tn.c_n - psi;
            b.valence[i] = tn.c_p - psi + s.gap;
            b.electron_fermi[i] =
                b.conduction[i] + theta * level(x[4 * i + 1], theta, s.nie, tn.log_nc).x;
            b.hole_fermi[i] = b.valence[i] - theta * level(x[4 * i + 2], theta, s.nie, tn.log_nv).x;
            continue;
        }
        const detail::NodeBands e =
            detail::node_bands(fermi_dirac_, n_ie_[i], log_dos_n_[i], log_dos_p_[i],
                               x[m_ * i] + band_shift_[i], x[m_ * i + 1], x[m_ * i + 2]);
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
    std::vector<double> charge = gates_.charges(x, m_, contact_count_);
    for (const EdgeTerm& e : edges_) {
        const std::int32_t ca = electrode_[e.a], cb = electrode_[e.b];
        if (ca == cb || e.charged) continue;
        const double flux = e.c * (x[m_ * e.a] - x[m_ * e.b]);  // from a to b
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
        return values[block_[m_ * m_ * node + m_ * r + c]];
    };
    if (electrothermal()) {  // the storage of a time step's rows at s, about a steady state
        add_thermal_storage(x, s, values);
        return;  // no interface terms: the model refuses traps and interface charge
    }
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
            v, x[m_ * o], x[m_ * sn], x[m_ * sn + 1], x[m_ * sn + 2], statistics(sn), inv_weight);
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
            d[m_ * i] = -1.0 / V_T_;
        } else if (gates_.on_gate(i) && gates_.contact(i) == contact) {
            d[m_ * i] = gates_.term(i).coupling / V_T_;
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
    if (electrothermal()) {  // with the temperature columns
        const std::vector<ThermalState> t = thermal_states(x);
        const std::vector<NodeLevels> l = thermal_levels(x, t);
        for (std::size_t k = 0; k < edges_.size(); ++k) {
            const EdgeTerm& e = edges_[k];
            const std::int32_t ca = contact_[e.a], cb = contact_[e.b];
            if (ca == cb || !e.carriers) continue;
            const auto [fn, fp] = thermal_fluxes(k, x, t, l);
            const std::pair<std::size_t, double> d[8] = {
                {4 * e.a, fn.d_psi[0] + fp.d_psi[0]},     {4 * e.b, fn.d_psi[1] + fp.d_psi[1]},
                {4 * e.a + 1, fn.d_c[0]},                 {4 * e.b + 1, fn.d_c[1]},
                {4 * e.a + 2, fp.d_c[0]},                 {4 * e.b + 2, fp.d_c[1]},
                {4 * e.a + 3, fn.d_theta[0] + fp.d_theta[0]},
                {4 * e.b + 3, fn.d_theta[1] + fp.d_theta[1]}};
            for (const auto& [column, value] : d) {
                add(ca, column, value);
                add(cb, column, -value);
            }
        }
    }
    const std::vector<NodeDegeneracy> g = degeneracies(x);
    for (const EdgeTerm& e : edges_) {
        const std::int32_t ca = contact_[e.a], cb = contact_[e.b];
        if (electrothermal() || ca == cb || !e.carriers) continue;
        const auto [fn, fp] = edge_fluxes(e, x, g);
        const std::pair<std::size_t, double> d[6] = {
            {m_ * e.a, fn.d_psi1 + fp.d_psi1}, {m_ * e.b, fn.d_psi2 + fp.d_psi2},
            {m_ * e.a + 1, fn.d_c1},           {m_ * e.b + 1, fn.d_c2},
            {m_ * e.a + 2, fp.d_c1},           {m_ * e.b + 2, fp.d_c2}};
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
            v, x[m_ * o], x[m_ * sn], x[m_ * sn + 1], x[m_ * sn + 2], statistics(sn), inv_weight);
        const std::size_t columns[4] = {m_ * o, m_ * sn, m_ * sn + 1, m_ * sn + 2};
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
        terms[c].emplace_back(m_ * i, -s * G);
        rows[c].bias += s * (G / V_T_);
    }
    for (const EdgeTerm& e : edges_) {
        const std::int32_t ca = owner(e.a), cb = owner(e.b);
        if (ca == cb || e.charged) continue;
        add(ca, m_ * e.a, s * e.c);
        add(ca, m_ * e.b, -s * e.c);
        add(cb, m_ * e.a, -s * e.c);
        add(cb, m_ * e.b, s * e.c);
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

void DriftDiffusion::add_local_tunnelling(std::span<const double> x, std::span<double> f,
                                          std::span<double> values) const {
    const bool jacobian = !values.empty();
    const int D = dimension_;
    for (const TunnelNode& node : tunnel_nodes_) {
        const std::size_t i = node.node;
        // The field at the node by least squares from its edges' components
        // E_k = -field (psi_b - psi_a): E = P sum_k t_k E_k.
        double E[3] = {};
        for (std::size_t k = node.first; k < node.last; ++k) {
            const FieldEdge& fe = field_edges_[k];
            const EdgeTerm& e = edges_[fe.edge];
            const double Ek = -fe.field * (x[m_ * e.b] - x[m_ * e.a]);
            for (int r = 0; r < D; ++r) {
                double w = 0.0;  // (P t_k)_r
                for (int c = 0; c < D; ++c) w += node.inverse[r * D + c] * fe.t[c];
                E[r] += w * Ek;
            }
        }
        double F2 = 0.0;
        for (int r = 0; r < D; ++r) F2 += E[r] * E[r];
        const double F = std::sqrt(F2);
        const physics::KaneRate g = physics::kane_generation(node.A, node.B, F);
        if (g.rate == 0.0 && g.d_dF == 0.0) continue;
        const double k_row = volume_[i] * generation_scale_;
        f[m_ * i + 1] += k_row * g.rate;
        f[m_ * i + 2] -= k_row * g.rate;
        if (!jacobian) continue;
        // dG/dE_r = G' E_r / F (F > 0 here: G vanishes at F = 0); an edge's component enters
        // through (P t_k), and E_k through its ends' psi.
        for (std::size_t k = node.first; k < node.last; ++k) {
            const FieldEdge& fe = field_edges_[k];
            const EdgeTerm& e = edges_[fe.edge];
            double dEk = 0.0;  // dG / dE_k
            for (int r = 0; r < D; ++r) {
                double w = 0.0;
                for (int c = 0; c < D; ++c) w += node.inverse[r * D + c] * fe.t[c];
                dEk += w * g.d_dF * E[r] / F;
            }
            const double d_b = -dEk * fe.field;  // dG / dpsi_b; dG / dpsi_a = -d_b
            const bool i_is_a = e.a == i;
            const std::size_t* p = i_is_a ? e.ab : e.ba;  // row i, the other end's psi column
            const double own = i_is_a ? -d_b : d_b, other = i_is_a ? d_b : -d_b;
            values[block_[m_ * m_ * i + m_]] += k_row * own;   // (n, psi) of i
            values[block_[m_ * m_ * i + 2 * m_]] -= k_row * own;   // (p, psi) of i
            values[p[1]] += k_row * other;              // (n, psi) of the other end
            values[p[3]] -= k_row * other;              // (p, psi) of the other end
        }
    }
}

TunnelPaths DriftDiffusion::trace_paths(std::span<const double> x) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    if (!tunnelling()) return {};
    const std::size_t n = node_count();
    std::vector<double> psi(n), valence(n, 0.0), end(n, 0.0);
    std::vector<char> semiconductor(n, 0), contact(n, 0);
    const bool wkb = tunnel_kind_ == NonlocalTunnelling::direct_wkb;
    for (std::size_t i = 0; i < n; ++i) {
        psi[i] = x[m_ * i];
        if (insulator_[i] != 0) continue;
        const double eta = psi[i] + band_shift_[i];
        valence[i] = -(log_dos_p_[i] + eta);
        end[i] = wkb ? valence[i] + tunnel_gap_[i] / V_T_ : log_dos_n_[i] - eta;
        semiconductor[i] = 1;
        contact[i] = contact_[i] >= 0 ? 1 : 0;
    }
    const TunnelTraceInput in{&*cells_, psi, valence, end, semiconductor, contact,
                              tunnel_length_, V_T_, wkb};
    TunnelPaths traced = trace_tunnel_paths(in);
    // Only the paths whose rate at x is at least `matters` of the largest are kept: the others
    // carry no current the state resolves (a start deep in a neutral region crosses after the
    // whole depletion width, at a mean field far below the peak), but each would couple its start
    // to a distant crossing in the Jacobian and fill the factors.
    std::vector<double> rate;
    double largest = 0.0;
    PathTerms t;
    for (const TunnelPath& p : traced.paths) {
        path_terms(p, layout(p), x, t);
        rate.push_back(t.state.rate_cm3_s);
        largest = std::max(largest, t.state.rate_cm3_s);
    }
    TunnelPaths kept;
    for (std::size_t k = 0; k < traced.paths.size(); ++k) {
        if (rate[k] > 0.0 && rate[k] >= matters * largest) {
            kept.paths.push_back(std::move(traced.paths[k]));
        }
    }
    return kept;
}

DriftDiffusion::PathLayout DriftDiffusion::layout(const TunnelPath& p) const {
    NITCAD_EXPECTS(p.crossing + 1 < p.forward.size());
    const auto nodes_of = [](const TunnelSample& s, std::vector<std::size_t>& into) {
        for (int k = 0; k < s.count; ++k) {
            const auto i = static_cast<std::size_t>(k);
            if (s.weights[i] != 0.0) into.push_back(s.nodes[i]);
        }
    };
    const auto unique = [](std::vector<std::size_t>& v) {
        std::ranges::sort(v);
        v.erase(std::unique(v.begin(), v.end()), v.end());
    };
    PathLayout l;
    nodes_of(p.forward[p.crossing], l.deposit);
    nodes_of(p.forward[p.crossing + 1], l.deposit);
    unique(l.deposit);
    l.columns = l.deposit;
    l.columns.push_back(p.start);
    if (tunnel_kind_ == NonlocalTunnelling::direct_wkb) {
        for (const TunnelSample& s : p.forward) nodes_of(s, l.columns);
        for (const TunnelSample& s : p.backward) nodes_of(s, l.columns);
    }
    unique(l.columns);
    l.before_cm.assign(p.forward.size(), 0.0);
    for (std::size_t s = 0; s + 1 < p.forward.size(); ++s) {
        l.before_cm[s + 1] = l.before_cm[s] + p.length_cm[s];
    }
    return l;
}

bool DriftDiffusion::paths_agree(const TunnelPaths& a, const TunnelPaths& b,
                                 std::span<const double> x) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    const auto states = [&](const TunnelPaths& set) {
        std::vector<TunnelPathState> r;
        PathTerms t;
        for (const TunnelPath& p : set.paths) {
            path_terms(p, layout(p), x, t);
            r.push_back(t.state);
        }
        return r;
    };
    const std::vector<TunnelPathState> sa = states(a), sb = states(b);
    double largest = 0.0;
    for (const auto& s : sa) largest = std::max(largest, s.rate_cm3_s);
    for (const auto& s : sb) largest = std::max(largest, s.rate_cm3_s);
    // Every path of `from` that matters has a twin in `to` (both sorted by start) with the same
    // geometry; if their crossing segments differ (by one), the crossing lies at the sample they
    // share, to 1e-6 of a segment, so the two evaluate the same crossing.
    const auto covered = [&](const TunnelPaths& from, const std::vector<TunnelPathState>& sf,
                             const TunnelPaths& to, const std::vector<TunnelPathState>& st) {
        for (std::size_t k = 0; k < from.paths.size(); ++k) {
            if (!(sf[k].rate_cm3_s >= matters * largest) || sf[k].rate_cm3_s == 0.0) continue;
            const TunnelPath& p = from.paths[k];
            const auto it = std::ranges::lower_bound(to.paths, p.start, {}, &TunnelPath::start);
            if (it == to.paths.end() || it->start != p.start || !same_geometry(p, *it)) {
                return false;
            }
            const auto m = static_cast<std::size_t>(it - to.paths.begin());
            if (p.crossing == it->crossing) continue;
            const bool first = p.crossing < it->crossing;  // p's segment ends where it's starts
            const double t_p = sf[k].fraction, t_q = st[m].fraction;
            const bool shared = first ? std::abs(t_p - 1.0) <= 1e-6 && std::abs(t_q) <= 1e-6
                                      : std::abs(t_p) <= 1e-6 && std::abs(t_q - 1.0) <= 1e-6;
            if (!shared) return false;
        }
        return true;
    };
    return covered(a, sa, b, sb) && covered(b, sb, a, sa);
}

void DriftDiffusion::set_paths(TunnelPaths paths) {
    NITCAD_EXPECTS(tunnelling());
    paths_ = std::move(paths);
    layouts_.clear();
    std::vector<linalg::Triplet> t = base_triplets_;
    for (const TunnelPath& p : paths_.paths) {
        PathLayout l = layout(p);
        for (const std::size_t c : l.columns) {
            t.push_back({static_cast<linalg::Index>(m_ * p.start + 2),
                         static_cast<linalg::Index>(m_ * c), 0.0});
            for (const std::size_t d : l.deposit) {
                t.push_back({static_cast<linalg::Index>(m_ * d + 1),
                             static_cast<linalg::Index>(m_ * c), 0.0});
            }
        }
        layouts_.push_back(std::move(l));
    }
    const auto size = static_cast<linalg::Index>(unknowns());
    pattern_ = *linalg::SparseMatrix::from_triplets(size, size, t);  // in range
    locate_positions();
    for (std::size_t k = 0; k < paths_.paths.size(); ++k) {
        PathLayout& l = layouts_[k];
        const std::size_t i = paths_.paths[k].start;
        for (const std::size_t c : l.columns) {
            l.hole.push_back(detail::position(pattern_, m_ * i + 2, m_ * c));
        }
        for (const std::size_t d : l.deposit) {
            for (const std::size_t c : l.columns) {
                l.electron.push_back(detail::position(pattern_, m_ * d + 1, m_ * c));
            }
        }
    }
}

void DriftDiffusion::locate_positions() {
    for (std::size_t i = 0; i < node_count(); ++i) {
        for (std::size_t r = 0; r < 3; ++r) {
            for (std::size_t c = 0; c < 3; ++c) {
                block_[m_ * m_ * i + m_ * r + c] =
                    detail::position(pattern_, m_ * i + r, m_ * i + c);
            }
        }
    }
    for (EdgeTerm& e : edges_) {
        for (std::size_t k = 0; k < 5; ++k) {
            e.ab[k] = detail::position(pattern_, m_ * e.a + edge_rows[k], m_ * e.b + edge_cols[k]);
            e.ba[k] = detail::position(pattern_, m_ * e.b + edge_rows[k], m_ * e.a + edge_cols[k]);
        }
    }
    for (const ImpactNode& node : impact_nodes_) {
        for (std::size_t k = node.first; k < node.last; ++k) {
            ImpactEdge& ie = impact_edges_[k];
            const EdgeTerm& e = edges_[ie.edge];
            const std::size_t other = e.a == node.node ? e.b : e.a;
            ie.n_p = detail::position(pattern_, m_ * node.node + 1, m_ * other + 2);
            ie.p_n = detail::position(pattern_, m_ * node.node + 2, m_ * other + 1);
        }
    }
    for (std::size_t k = 0; k < interfaces_.edges().size(); ++k) {
        const InterfaceEdge& f = interfaces_.edges()[k];
        interface_np_[k] = {detail::position(pattern_, m_ * f.insulator, m_ * f.semiconductor + 1),
                            detail::position(pattern_, m_ * f.insulator, m_ * f.semiconductor + 2)};
    }
}

void DriftDiffusion::path_terms(const TunnelPath& p, const PathLayout& l,
                                std::span<const double> x, PathTerms& t) const {
    const std::size_t nc = l.columns.size(), nd = l.deposit.size();
    const std::size_t i = p.start;
    t.pairs = 0.0;
    t.d_pairs.assign(nc, 0.0);
    t.share.assign(nd, 0.0);
    t.d_share.assign(nd * nc, 0.0);
    t.state = {i, 0.0, 0.0, 0.0, 0.0, false};
    const auto column = [&](std::size_t node) {
        return static_cast<std::size_t>(std::ranges::lower_bound(l.columns, node) -
                                        l.columns.begin());
    };
    const std::size_t ci = column(i);
    // A sample's psi and E_v (units of V_T), each linear in psi with partials w (E_v: -w).
    const auto value = [&](const TunnelSample& s, auto&& f) {
        double v = 0.0;
        for (int q = 0; q < s.count; ++q) {
            const auto m = static_cast<std::size_t>(q);
            if (s.weights[m] != 0.0) v += s.weights[m] * f(s.nodes[m]);
        }
        return v;
    };
    const auto psi_at = [&](std::size_t node) { return x[m_ * node]; };
    const auto valence_at = [&](std::size_t node) {
        return -(log_dos_p_[node] + x[m_ * node] + band_shift_[node]);
    };
    const auto end_at = [&](std::size_t node) {
        return tunnel_kind_ == NonlocalTunnelling::kane
                   ? log_dos_n_[node] - x[m_ * node] - band_shift_[node]
                   : valence_at(node) + tunnel_gap_[i] / V_T_;
    };
    // d/dpsi of a sample's linear value with weights times `sign`, added into d.
    const auto add_weights = [&](const TunnelSample& s, double sign, std::vector<double>& d) {
        for (int q = 0; q < s.count; ++q) {
            const auto m = static_cast<std::size_t>(q);
            if (s.weights[m] != 0.0) d[column(s.nodes[m])] += sign * s.weights[m];
        }
    };

    // The crossing on the frozen segment (a, b): Delta = E_v(i) - E_end, dDelta = -e_i + w.
    const TunnelSample& sa = p.forward[p.crossing];
    const TunnelSample& sb = p.forward[p.crossing + 1];
    const double Ev_i = valence_at(i);
    const double Da = Ev_i - value(sa, end_at), Db = Ev_i - value(sb, end_at);
    const double den = Da - Db;
    if (!(den < 0.0)) return;  // the band does not fall along the segment: no crossing on it
    const double tt = Da / den;
    std::vector<double> dDa(nc, 0.0), dDb(nc, 0.0), dt(nc);
    dDa[ci] = dDb[ci] = -1.0;
    add_weights(sa, 1.0, dDa);
    add_weights(sb, 1.0, dDb);
    for (std::size_t c = 0; c < nc; ++c) dt[c] = (-Db * dDa[c] + Da * dDb[c]) / (den * den);
    const double La = p.length_cm[p.crossing];
    const double length = l.before_cm[p.crossing] + tt * La;
    const double psi_a = value(sa, psi_at), psi_b = value(sb, psi_at);
    const double drop = (1.0 - tt) * psi_a + tt * psi_b - x[m_ * i];
    t.state.length_cm = length;
    t.state.fraction = tt;
    t.state.reached = tt >= -1e-6 && tt <= 1.0 + 1e-6;
    if (length > 0.0) t.state.field_V_per_cm = V_T_ * drop / length;

    double G = 0.0;
    std::vector<double> dG(nc, 0.0);
    if (tunnel_kind_ == NonlocalTunnelling::kane) {
        if (!(length > 0.0 && drop > 0.0)) return;
        // F = V_T (psi_f - psi_i) / l: dF = V_T ((dpsi_f - e_i) / l - drop dl / l^2), with
        // dpsi_f = (1 - t) w_a + t w_b + (psi_b - psi_a) dt and dl = L_a dt.
        const double F = V_T_ * drop / length;
        const physics::KaneRate g = physics::kane_generation(tunnel_A_[i], tunnel_B_[i], F);
        G = g.rate;
        std::vector<double> dpf(nc, 0.0);
        add_weights(sa, 1.0 - tt, dpf);
        add_weights(sb, tt, dpf);
        for (std::size_t c = 0; c < nc; ++c) {
            const double dpsi_f = dpf[c] + (psi_b - psi_a) * dt[c] - (c == ci ? 1.0 : 0.0);
            const double dF = V_T_ * (dpsi_f / length - drop * La * dt[c] / (length * length));
            dG[c] = g.d_dF * dF;
        }
    } else {
        // eq. (11) over the forward samples, in SI: delta = (E_v(i) - E_v(s)) V_T / E_g.
        constexpr double q = base::q_C, hbar = base::hbar_J_s;
        const double Eg = tunnel_gap_[i];
        const double cd = V_T_ / Eg;
        const double mr = tunnel_mc_[i] * tunnel_mv_[i] / (tunnel_mc_[i] + tunnel_mv_[i]);
        const physics::WkbBand band{Eg * q, mr};
        const std::size_t ns = p.forward.size();
        std::vector<double> delta(ns), Ev(ns);
        for (std::size_t s = 0; s < ns; ++s) {
            Ev[s] = value(p.forward[s], valence_at);
            delta[s] = cd * (Ev_i - Ev[s]);
        }
        // dδ_s = cd (-e_i + w_s).
        const auto add_delta = [&](std::size_t s, double scale, std::vector<double>& d) {
            d[ci] -= scale * cd;
            add_weights(p.forward[s], scale * cd, d);
        };
        double Ik = 0.0, Iik = 0.0;
        std::vector<double> dIk(nc, 0.0), dIik(nc, 0.0);
        for (std::size_t s = 0; s + 1 < ns; ++s) {
            const physics::WkbSegment w =
                physics::wkb_segment(delta[s], delta[s + 1], p.length_cm[s] * 1e-2, band);
            Ik += w.I_k;
            Iik += w.I_ik;
            add_delta(s, w.dI_k_da, dIk);
            add_delta(s + 1, w.dI_k_db, dIk);
            add_delta(s, w.dI_ik_da, dIik);
            add_delta(s + 1, w.dI_ik_db, dIik);
        }
        // |dE_v/dx| at the start over the first segment [J/m]: q V_T (E_v(i) - E_v(1)) / L_0.
        const double L0 = p.length_cm[0] * 1e-2;
        const double slope = q * V_T_ * (Ev_i - Ev[1]) / L0;
        std::vector<double> dslope(nc, 0.0);
        dslope[ci] -= q * V_T_ / L0;
        add_weights(p.forward[1], q * V_T_ / L0, dslope);
        // k_m^2 from E_vmax over the start and the backward samples and E_cmin = min(E_v + E_g)
        // over the forward ones.
        double Evmax = Ev_i;
        const TunnelSample* at_max = nullptr;  // nullptr: the start node
        for (const TunnelSample& s : p.backward) {
            const double v = value(s, valence_at);
            if (v > Evmax) {
                Evmax = v;
                at_max = &s;
            }
        }
        std::size_t at_min = 0;
        for (std::size_t s = 1; s < ns; ++s) {
            if (Ev[s] < Ev[at_min]) at_min = s;
        }
        const double Ecmin = Ev[at_min] + Eg / V_T_;
        const double kv = 2.0 * tunnel_mv_[i] * q * V_T_ / (hbar * hbar);
        const double kc = 2.0 * tunnel_mc_[i] * q * V_T_ / (hbar * hbar);
        const double win_v = kv * (Evmax - Ev_i), win_c = kc * (Ev_i - Ecmin);
        const double km2 = std::max(std::min(win_v, win_c), 0.0);
        std::vector<double> dkm2(nc, 0.0);
        if (km2 > 0.0) {
            if (win_v <= win_c) {  // kv (dE_vmax - dE), dE = -e_i
                if (at_max != nullptr) {
                    add_weights(*at_max, -kv, dkm2);
                    dkm2[ci] += kv;
                }  // at the start dE_vmax = dE: no dependence
            } else {  // kc (dE - dE_cmin), dE_cmin = -w_min
                dkm2[ci] -= kc;
                add_weights(p.forward[at_min], kc, dkm2);
            }
        }
        const physics::WkbPathRate r = physics::wkb_path_rate(slope, Ik, Iik, km2);
        G = r.rate * 1e-6;  // [cm^-3 s^-1]
        for (std::size_t c = 0; c < nc; ++c) {
            dG[c] = 1e-6 * (r.d_slope * dslope[c] + r.d_I_k * dIk[c] + r.d_I_ik * dIik[c] +
                            r.d_km2 * dkm2[c]);
        }
    }
    t.state.rate_cm3_s = G;
    const double k_row = volume_[i] * generation_scale_;
    t.pairs = k_row * G;
    for (std::size_t c = 0; c < nc; ++c) t.d_pairs[c] = k_row * dG[c];
    // Deposit: (1 - t) on sample a's stencil, t on b's.
    const auto deposit = [&](const TunnelSample& s, double share, double sign) {
        for (int q = 0; q < s.count; ++q) {
            const auto m = static_cast<std::size_t>(q);
            if (s.weights[m] == 0.0) continue;
            const auto d = static_cast<std::size_t>(
                std::ranges::lower_bound(l.deposit, s.nodes[m]) - l.deposit.begin());
            const double w = s.weights[m];
            t.share[d] += share * w;
            for (std::size_t c = 0; c < nc; ++c) t.d_share[d * nc + c] += sign * w * dt[c];
        }
    };
    deposit(sa, 1.0 - tt, -1.0);
    deposit(sb, tt, 1.0);
}

void DriftDiffusion::add_path_tunnelling(std::span<const double> x, std::span<double> f,
                                         std::span<double> values) const {
    const bool jacobian = !values.empty();
    PathTerms t;
    for (std::size_t k = 0; k < paths_.paths.size(); ++k) {
        const PathLayout& l = layouts_[k];
        path_terms(paths_.paths[k], l, x, t);
        if (t.pairs == 0.0 && std::ranges::all_of(t.d_pairs, [](double v) { return v == 0.0; })) {
            continue;
        }
        const std::size_t nc = l.columns.size();
        const std::size_t i = paths_.paths[k].start;
        f[m_ * i + 2] -= t.pairs;
        for (std::size_t d = 0; d < l.deposit.size(); ++d) {
            f[m_ * l.deposit[d] + 1] += t.share[d] * t.pairs;
        }
        if (!jacobian) continue;
        for (std::size_t c = 0; c < nc; ++c) {
            values[l.hole[c]] -= t.d_pairs[c];
            for (std::size_t d = 0; d < l.deposit.size(); ++d) {
                values[l.electron[d * nc + c]] +=
                    t.share[d] * t.d_pairs[c] + t.pairs * t.d_share[d * nc + c];
            }
        }
    }
}

std::vector<TunnelPathState> DriftDiffusion::path_states(std::span<const double> x) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    std::vector<TunnelPathState> states;
    PathTerms t;
    for (std::size_t k = 0; k < paths_.paths.size(); ++k) {
        path_terms(paths_.paths[k], layouts_[k], x, t);
        states.push_back(t.state);
    }
    return states;
}

InterfaceStatistics DriftDiffusion::statistics(std::size_t node) const noexcept {
    return {fermi_dirac_, n_ie_[node], log_dos_n_[node], log_dos_p_[node], band_shift_[node]};
}

}  // namespace NiTCAD::assemble
