#include "NiTCAD/assemble/drift_diffusion.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numbers>
#include <optional>
#include <string>
#include <variant>

#include "NiTCAD/assemble/contact_bias.hpp"
#include "NiTCAD/assemble/observer.hpp"
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
#include "composition.hpp"
#include "scaled_device.hpp"
#include "terms.hpp"

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
    // Unit 34: Masetti's coefficients are boron-doped silicon's, so every semiconductor region must
    // carry silicon's hole mobility (its Caughey-Thomas set stands for "silicon" here).
    if (models.hole_mobility == HoleMobility::masetti) {
        for (const device::Region& r : device.regions()) {
            const auto* s = std::get_if<physics::Semiconductor>(&r.material);
            if (s && !(s->parameters().hole_mobility == physics::silicon_parameters.hole_mobility)) {
                return std::unexpected(base::Error{
                    base::ErrorCode::invalid_input,
                    "hole_mobility masetti: region '" + r.name +
                        "' does not have silicon's hole mobility (Masetti's coefficients are "
                        "boron-doped silicon's)",
                    {}});
            }
        }
    } else if (models.hole_mobility != HoleMobility::caughey_thomas &&
               models.hole_mobility != HoleMobility::automatic) {
        return refuse("unknown hole mobility model");
    }
    // Unit 38: likewise for Masetti's phosphorus-doped silicon electrons.
    if (models.electron_mobility == ElectronMobility::masetti) {
        for (const device::Region& r : device.regions()) {
            const auto* s = std::get_if<physics::Semiconductor>(&r.material);
            if (s && !(s->parameters().electron_mobility ==
                       physics::silicon_parameters.electron_mobility)) {
                return std::unexpected(base::Error{
                    base::ErrorCode::invalid_input,
                    "electron_mobility masetti: region '" + r.name +
                        "' does not have silicon's electron mobility (Masetti's coefficients are "
                        "phosphorus-doped silicon's)",
                    {}});
            }
        }
    } else if (models.electron_mobility != ElectronMobility::caughey_thomas) {
        return refuse("unknown electron mobility model");
    }
    // Unit 46 (DECISIONS.md LD6, LD8): Lombardi's coefficients are silicon's; impact ionization
    // reads the edge fluxes' partials, which with E_perp reach past its pattern.
    if (models.surface_mobility == SurfaceMobility::lombardi) {
        for (const device::Region& r : device.regions()) {
            const auto* s = std::get_if<physics::Semiconductor>(&r.material);
            if (s && !(s->parameters().electron_mobility ==
                           physics::silicon_parameters.electron_mobility &&
                       s->parameters().hole_mobility ==
                           physics::silicon_parameters.hole_mobility)) {
                return std::unexpected(base::Error{
                    base::ErrorCode::invalid_input,
                    "surface_mobility lombardi: region '" + r.name +
                        "' does not have silicon's mobility (Lombardi's coefficients are "
                        "silicon's)",
                    {}});
            }
        }
        if (models.impact_ionization) {
            return refuse("surface_mobility lombardi with impact_ionization is not supported");
        }
    } else if (models.surface_mobility != SurfaceMobility::off) {
        return refuse("unknown surface mobility model");
    }
    auto scaled = detail::make_scaled_device(device, scaling, models);
    if (!scaled) return std::unexpected(std::move(scaled.error()));
    const std::size_t n = scaled->volume.size();
    const double T = scaling.temperature_K;
    // Unit 48 (DECISIONS.md SD7): a Schottky contact is isothermal and on a device without
    // thermionic heterointerfaces.
    if (std::ranges::any_of(scaled->schottky, [](std::int32_t c) { return c >= 0; })) {
        if (models.electrothermal) {
            return refuse("a Schottky contact with the electrothermal model is not supported");
        }
        if (std::ranges::any_of(scaled->edges,
                                [](const detail::ScaledEdge& e) { return e.thermionic; })) {
            return refuse(
                "a Schottky contact on a device with a thermionic-emission interface is not "
                "supported");
        }
    }
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
    s.schottky_ = std::move(scaled->schottky);
    s.schottky_eta_ = std::move(scaled->schottky_eta);
    s.schottky_kn_ = std::move(scaled->schottky_kn);
    s.schottky_kp_ = std::move(scaled->schottky_kp);
    s.owner_.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        s.owner_[i] = s.contact_[i] >= 0 ? s.contact_[i] : s.schottky_[i];
    }
    s.insulator_ = std::move(scaled->insulator);
    s.electrode_ = std::move(scaled->electrode);
    s.electrode_potential_ = std::move(scaled->electrode_potential);
    s.interfaces_ = std::move(scaled->interfaces);
    s.gates_ = std::move(scaled->gates);
    s.kinds_ = std::move(scaled->kinds);
    s.contact_count_ = device.contacts().size();
    s.temperature_K_ = T;
    // Unit 51 (DECISIONS.md GD9): gate tunnelling is isothermal.
    if (std::ranges::any_of(device.contacts(),
                            [](const device::Contact& c) { return c.tunnelling != 0; })) {
        if (models.electrothermal) {
            return refuse("gate tunnelling with the electrothermal model is not supported");
        }
        if (auto e = s.make_gate_tunnels(device, scaling)) return std::unexpected(std::move(*e));
    }
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
    std::vector<double> mu_n(n), mu_p(n), N_bulk(n, 0.0);
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
        N_bulk[i] = N_mobility;
        mu_n[i] = models.electron_mobility == ElectronMobility::masetti
                      ? physics::masetti_mobility(physics::masetti_phosphorus_silicon, N_mobility, T)
                      : physics::caughey_thomas_mobility(m, physics::Carrier::electron, N_mobility,
                                                         T);
        mu_p[i] = masetti_holes(models.hole_mobility, m.parameters())
                      ? physics::masetti_mobility(physics::masetti_boron_silicon, N_mobility, T)
                      : physics::caughey_thomas_mobility(m, physics::Carrier::hole, N_mobility, T);
        s.tau_n_[i] = physics::scharfetter_lifetime(m, physics::Carrier::electron, N);
        s.tau_p_[i] = physics::scharfetter_lifetime(m, physics::Carrier::hole, N);
        s.auger_n_[i] = m.parameters().auger.Cn;
        s.auger_p_[i] = m.parameters().auger.Cp;
    }

    // The system's equations and unknowns from the device and the models (Unit 40,
    // composition.hpp): psi, n, p per node, and with the electrothermal model the temperature rise
    // last; every node and edge block is then full 4 x 4.
    s.composition_ = Composition::compose(s, models);
    s.m_ = s.composition_->dofs().block();
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
        if (models.electrothermal) continue;  // full 4 x 4 blocks (ThermalEdges, terms.hpp)
        for (std::size_t k = 0; k < 5; ++k) {
            triplets.push_back({static_cast<linalg::Index>(stride * e.i + edge_rows[k]),
                                static_cast<linalg::Index>(stride * e.j + edge_cols[k]), 0.0});
            triplets.push_back({static_cast<linalg::Index>(stride * e.j + edge_rows[k]),
                                static_cast<linalg::Index>(stride * e.i + edge_cols[k]), 0.0});
        }
    }
    // The edges' data first: the terms' pattern entries (Term::pattern) read it.
    s.edges_.reserve(scaled->edges.size());
    for (const detail::ScaledEdge& e : scaled->edges) {
        EdgeTerm t{};
        t.a = e.i;
        t.b = e.j;
        t.c = e.et * e.geometry;
        t.carriers = e.carriers;
        t.charged = e.charged;
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
        t.surface = models.surface_mobility == SurfaceMobility::lombardi && !e.thermionic;
        t.N_a = N_bulk[e.i];
        t.N_b = N_bulk[e.j];
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

    if (models.surface_mobility == SurfaceMobility::lombardi) {
        // Lombardi (Unit 46): each surface edge's E_perp stencil. Node i's field is
        // E_i = sum_k (P_i t_k) E_k over its carrier edges k, E_k = -field_k (psi_kb - psi_ka),
        // P_i the pseudo-inverse of sum_k t_k t_k^T (as impact ionization's); the edge's field is
        // the ends' mean, and v its part normal to the edge, (I - t t^T) of it.
        s.surface_ = true;
        s.surface_T_ = T;
        const mesh::Mesh& m = device.mesh();
        const int D = m.dimension();
        const std::size_t edge_count = scaled->edges.size();
        std::vector<std::vector<std::size_t>> incident(n);
        std::vector<std::array<double, 3>> unit(edge_count, {0.0, 0.0, 0.0});
        for (std::size_t k = 0; k < edge_count; ++k) {
            const detail::ScaledEdge& e = scaled->edges[k];
            if (!e.carriers) continue;
            incident[e.i].push_back(k);
            incident[e.j].push_back(k);
            const mesh::Point& pa = m.points()[e.i];
            const mesh::Point& pb = m.points()[e.j];
            for (int d = 0; d < D; ++d) unit[k][d] = (pb[d] - pa[d]) / e.length_cm;
        }
        std::vector<std::array<double, 9>> inverse(n);
        for (std::size_t i = 0; i < n; ++i) {
            inverse[i].fill(0.0);
            if (s.insulator_[i] != 0 || incident[i].empty()) continue;
            double outer[9] = {};
            for (const std::size_t k : incident[i]) {
                for (int r = 0; r < D; ++r) {
                    for (int c = 0; c < D; ++c) outer[r * D + c] += unit[k][r] * unit[k][c];
                }
            }
            pseudo_inverse(outer, D, inverse[i].data());
        }
        using Coefficient = std::pair<std::size_t, std::array<double, 3>>;
        std::vector<Coefficient> c;
        const auto accumulate = [&](std::size_t node, const double* w, double scale) {
            auto it = std::ranges::find(c, node, &Coefficient::first);
            if (it == c.end()) {
                c.push_back({node, {0.0, 0.0, 0.0}});
                it = c.end() - 1;
            }
            for (int d = 0; d < D; ++d) it->second[d] += scale * w[d];
        };
        const auto add_field = [&](std::size_t i) {  // half of node i's field
            for (const std::size_t k : incident[i]) {
                const detail::ScaledEdge& e = scaled->edges[k];
                const double field = scaling.V_T / e.length_cm;
                double w[3] = {};
                for (int r = 0; r < D; ++r) {
                    for (int q = 0; q < D; ++q) w[r] += inverse[i][r * D + q] * unit[k][q];
                }
                accumulate(e.j, w, -0.5 * field);
                accumulate(e.i, w, 0.5 * field);
            }
        };
        for (std::size_t k = 0; k < s.edges_.size(); ++k) {
            EdgeTerm& t = s.edges_[k];
            t.surface_first = t.surface_last = s.surface_stencil_.size();
            if (!t.surface) continue;
            c.clear();
            add_field(t.a);
            add_field(t.b);
            std::ranges::sort(c, {}, &Coefficient::first);
            const std::array<double, 3>& u = unit[k];
            for (const auto& [node, v] : c) {
                double along = 0.0;
                for (int d = 0; d < D; ++d) along += v[d] * u[d];
                SurfaceStencil j{node, {0.0, 0.0, 0.0}, {}};
                bool zero = true;
                for (int d = 0; d < D; ++d) {
                    j.g[d] = v[d] - along * u[d];
                    zero = zero && j.g[d] == 0.0;
                }
                if (!zero) s.surface_stencil_.push_back(j);
            }
            t.surface_last = s.surface_stencil_.size();
        }
    }

    // The terms' own entries (terms.hpp); set_paths adds the paths' to the base again.
    if (nonlocal) s.base_triplets_ = triplets;
    Terms::pattern(s, triplets);
    const auto size = static_cast<linalg::Index>(stride * n);
    auto pattern = linalg::SparseMatrix::from_triplets(size, size, triplets);
    if (!pattern) return std::unexpected(std::move(pattern.error()));
    s.pattern_ = std::move(*pattern);

    s.block_.resize(stride * stride * n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t r = 0; r < stride; ++r) {
            for (std::size_t c = 0; c < stride; ++c) {
                s.block_[stride * stride * i + stride * r + c] =
                    detail::position(s.pattern_, stride * i + r, stride * i + c);
            }
        }
    }
    s.locate_gate_tunnels();
    for (EdgeTerm& t : s.edges_) {
        for (std::size_t k = 0; k < 5; ++k) {
            t.ab[k] = detail::position(s.pattern_, stride * t.a + edge_rows[k],
                                       stride * t.b + edge_cols[k]);
            t.ba[k] = detail::position(s.pattern_, stride * t.b + edge_rows[k],
                                       stride * t.a + edge_cols[k]);
        }
    }
    s.locate_surface();

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
                            physics::impact_ionization_temperature_factor(ii.phonon_energy_eV, T),
                            ii.phonon_energy_eV};
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
        if (schottky_[i] >= 0) {  // SD2: psi_S, and n_B, p_B at eta_B (bias-independent)
            const double eta = schottky_eta_[i];
            psi0_[i] =
                bias_V[static_cast<std::size_t>(schottky_[i])] / V_T_ + eta - band_shift_[i];
            n0_[i] = detail::density(fermi_dirac_, n_ie_[i], log_dos_n_[i], eta).density;
            p0_[i] = detail::density(fermi_dirac_, n_ie_[i], log_dos_p_[i], -eta).density;
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
    detail::observe("state_from_potential", x);
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
                if (electrode_[i] >= 0) {
                    x[4 * i] = psi0_[i] + electrode_shift(i, thermal_state(i, x[4 * i + 3])).value;
                }
                x[4 * i + 1] = x[4 * i + 2] = 0.0;
            } else if (contact_[i] >= 0) {
                const ThermalOhmic o = thermal_ohmic(i, thermal_state(i, x[4 * i + 3]));
                x[4 * i] = o.psi;
                x[4 * i + 1] = o.n;
                x[4 * i + 2] = o.p;
            }
        }
        detail::observe("stamp_contacts", x);
        return;
    }
    for (std::size_t i = 0; i < node_count(); ++i) {
        if (insulator_[i] != 0) {
            if (electrode_[i] >= 0) x[m_ * i] = psi0_[i];
            x[m_ * i + 1] = x[m_ * i + 2] = 0.0;
            continue;
        }
        if (schottky_[i] >= 0) x[m_ * i] = psi0_[i];  // its densities are unknowns
        if (contact_[i] < 0) continue;
        x[m_ * i] = psi0_[i];
        x[m_ * i + 1] = n0_[i];
        x[m_ * i + 2] = p0_[i];
    }
    detail::observe("stamp_contacts", x);
}

const DofMap& DriftDiffusion::dofs() const noexcept {
    NITCAD_EXPECTS(composition_ != nullptr);  // not a moved-from system
    return composition_->dofs();
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
    const EdgeTerm& e, std::span<const double> x, std::span<const NodeDegeneracy> g,
    SurfaceSlope* slope) const {
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
    // Unit edge factors when the Canali or Lombardi factor follows (never on a thermionic edge).
    const bool unit_factor = field_mobility_ || e.surface;
    EdgeFlux fn = e.thermionic ? thermionic_electron_flux(e.te_kn, delta_n, e.te_log_nc,
                                                          e.te_ratio_nc, na, nb)
                               : sg_electron_flux(unit_factor ? 1.0 : e.an, 0.0, delta_n,
                                                  na, nb);
    EdgeFlux fp = e.thermionic ? thermionic_hole_flux(e.te_kp, delta_p, e.te_log_nv,
                                                      e.te_ratio_nv, pa, pb)
                               : sg_hole_flux(unit_factor ? 1.0 : e.ap, 0.0, delta_p, pa, pb);
    if (!g.empty()) {
        // d flux / d delta is d_psi2; ln gamma_n,a enters delta_n with -, ln gamma_p,a delta_p
        // with +.
        fn.d_c1 -= fn.d_psi2 * g[e.a].n.d_density;
        fn.d_c2 += fn.d_psi2 * g[e.b].n.d_density;
        fp.d_c1 += fp.d_psi2 * g[e.a].p.d_density;
        fp.d_c2 -= fp.d_psi2 * g[e.b].p.d_density;
    }
    if (slope != nullptr) *slope = SurfaceSlope{};
    if (!unit_factor || e.thermionic) return {fn, fp};
    // Edge factor a(E) = a0 mu(E) / mu0, E = field |dpsi|; the flux is linear in a, so
    // d flux / d dpsi gains (flux / a) da/d dpsi. mu(E) is the Canali mobility of mu0 (Unit 13),
    // or between two materials the harmonic mean of each end's Canali mobility (Unit 15). With
    // Lombardi (Unit 46) the ends' low-field mobilities first take the surface terms in E_perp,
    // and mu0 is their harmonic mean; d_perp is d mu / d E_perp.
    const double E = e.field * std::abs(dpsi);
    const double dE = dpsi > 0.0 ? e.field : dpsi < 0.0 ? -e.field : 0.0;  // dE / d dpsi
    const SurfaceField perp = e.surface ? surface_field(e, x) : SurfaceField{};
    struct Mobility {
        physics::FieldMobility mu;
        double d_perp;
    };
    const auto mobility = [&](physics::Carrier carrier, double mu0, double mu_a, double mu_b,
                              const physics::CanaliParameters& sat_a,
                              const physics::CanaliParameters& sat_b) {
        double s_a = 0.0, s_b = 0.0;  // d mu_a / d E_perp, d mu_b / d E_perp
        if (e.surface) {
            const physics::SurfaceMobility la = surface_end(carrier, mu_a, e.N_a, surface_T_,
                                                            perp.value);
            const physics::SurfaceMobility lb = surface_end(carrier, mu_b, e.N_b, surface_T_,
                                                            perp.value);
            mu_a = la.mobility;
            mu_b = lb.mobility;
            s_a = la.d_dE;
            s_b = lb.d_dE;
            mu0 = harmonic_mean(mu_a, mu_b);
        }
        const double low = mu_a + mu_b;  // d hmean / d mu_a = 2 mu_b^2 / low^2
        const double d_mu0 =
            e.surface ? 2.0 * (mu_b * mu_b * s_a + mu_a * mu_a * s_b) / (low * low) : 0.0;
        if (!field_mobility_) return Mobility{{mu0, 0.0}, d_mu0};
        if (!e.mixed) {
            const physics::FieldMobility c = physics::canali_mobility(mu0, E, sat_a);
            return Mobility{c, (c.mobility + E * c.d_dE) / mu0 * d_mu0};
        }
        const physics::FieldMobility ma = physics::canali_mobility(mu_a, E, sat_a);
        const physics::FieldMobility mb = physics::canali_mobility(mu_b, E, sat_b);
        const double sum = ma.mobility + mb.mobility;
        // Canali's mu = mu0 f(mu0 E): d mu / d mu0 = (mu + E d mu / dE) / mu0.
        const double d_perp =
            e.surface ? 2.0 *
                            (mb.mobility * mb.mobility * (ma.mobility + E * ma.d_dE) / mu_a * s_a +
                             ma.mobility * ma.mobility * (mb.mobility + E * mb.d_dE) / mu_b * s_b) /
                            (sum * sum)
                      : 0.0;
        return Mobility{
            physics::FieldMobility{
                2.0 * ma.mobility * mb.mobility / sum,
                2.0 * (mb.mobility * mb.mobility * ma.d_dE + ma.mobility * ma.mobility * mb.d_dE) /
                    (sum * sum)},
            d_perp};
    };
    const auto scaled = [&](EdgeFlux unit, double a0, double mu0, physics::FieldMobility mu) {
        const double a = a0 * (mu.mobility / mu0);
        const double da = a0 * (mu.d_dE / mu0) * dE;
        return EdgeFlux{a * unit.flux, a * unit.d_psi1 - unit.flux * da,
                        a * unit.d_psi2 + unit.flux * da, a * unit.d_c1, a * unit.d_c2};
    };
    const Mobility mn = mobility(physics::Carrier::electron, e.mu_n, e.mu_n_a, e.mu_n_b, e.sat_n,
                                 e.sat_n_b);
    const Mobility mp =
        mobility(physics::Carrier::hole, e.mu_p, e.mu_p_a, e.mu_p_b, e.sat_p, e.sat_p_b);
    if (slope != nullptr && e.surface) {
        *slope = SurfaceSlope{perp, e.an * (mn.d_perp / e.mu_n) * fn.flux,
                              e.ap * (mp.d_perp / e.mu_p) * fp.flux};
    }
    return {scaled(fn, e.an, e.mu_n, mn.mu), scaled(fp, e.ap, e.mu_p, mp.mu)};
}

DriftDiffusion::SurfaceField DriftDiffusion::surface_field(const EdgeTerm& e,
                                                           std::span<const double> x) const {
    SurfaceField f{0.0, {0.0, 0.0, 0.0}};
    for (std::size_t k = e.surface_first; k < e.surface_last; ++k) {
        const SurfaceStencil& j = surface_stencil_[k];
        const double psi = x[m_ * j.node];
        for (std::size_t d = 0; d < 3; ++d) f.v[d] += j.g[d] * psi;
    }
    constexpr double E0 = 1.0;  // [V/cm] (DECISIONS.md LD2)
    f.value = std::sqrt(f.v[0] * f.v[0] + f.v[1] * f.v[1] + f.v[2] * f.v[2] + E0 * E0);
    return f;
}

double DriftDiffusion::surface_gradient(const SurfaceField& f, const SurfaceStencil& j) const {
    return (f.v[0] * j.g[0] + f.v[1] * j.g[1] + f.v[2] * j.g[2]) / f.value;
}

physics::SurfaceMobility DriftDiffusion::surface_end(physics::Carrier carrier, double mu_b,
                                                     double N, double T, double E_perp) {
    return physics::lombardi_mobility(carrier == physics::Carrier::electron
                                          ? physics::lombardi_silicon_electrons
                                          : physics::lombardi_silicon_holes,
                                      mu_b, N, T, E_perp);
}

void DriftDiffusion::locate_surface() {
    // Rows electron, hole and (electrothermal) heat of each end with psi_j's column.
    const std::size_t rows = electrothermal() ? 3 : 2;
    for (const EdgeTerm& e : edges_) {
        for (std::size_t k = e.surface_first; k < e.surface_last; ++k) {
            SurfaceStencil& j = surface_stencil_[k];
            for (std::size_t end = 0; end < 2; ++end) {
                const std::size_t node = end == 0 ? e.a : e.b;
                for (std::size_t r = 0; r < rows; ++r) {
                    j.pos[end][r] = detail::position(pattern_, m_ * node + 1 + r, m_ * j.node);
                }
            }
        }
    }
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
            l.donor_degeneracy, l.donor_bound);
        s.n -= d.value;
        s.d_n -= d.d_eta * (1.0 / n - wn);
    }
    if (acceptors_[i] > 0.0) {
        const double lp = fermi_dirac_ ? g[i].p.log_gamma : 0.0;
        const double wp = fermi_dirac_ ? g[i].p.d_density : 0.0;
        const physics::IonizedDensity a = physics::ionized_density(
            acceptors_[i], std::log(p / n_ie_[i]) - lp - log_dos_p_[i], l.acceptor_kT,
            l.acceptor_degeneracy, l.acceptor_bound);
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
    // The terms in assembly order (composition.hpp).
    const std::vector<NodeDegeneracy> g = degeneracies(x);
    composition_->assemble(*this, {x, f, values, step, g}, interfaces);
}

void DriftDiffusion::evaluate(std::span<const double> x, std::span<double> residual,
                              linalg::SparseMatrix& jacobian) const {
    NITCAD_EXPECTS(jacobian.has_same_pattern(pattern_));
    const std::span<double> values = jacobian.values();
    std::fill(values.begin(), values.end(), 0.0);
    assemble(x, residual, values, nullptr);
    detail::observe("evaluate.x", x);
    detail::observe("evaluate.residual", residual);
    detail::observe("evaluate.jacobian", values);
}

void DriftDiffusion::residual(std::span<const double> x, std::span<double> residual) const {
    assemble(x, residual, {}, nullptr);
    detail::observe("residual.x", x);
    detail::observe("residual.residual", residual);
}

namespace {

// The step's inputs, for the observer (observer.hpp).
void observe_step(const DriftDiffusion::TimeStep& step) {
    detail::observe("step.rate", step.rate);
    detail::observe("step.storage", step.storage);
    detail::observe("step.traps", step.traps);
}

}  // namespace

void DriftDiffusion::evaluate(std::span<const double> x, const TimeStep& step,
                              std::span<double> residual, linalg::SparseMatrix& jacobian) const {
    NITCAD_EXPECTS(jacobian.has_same_pattern(pattern_));
    const std::span<double> values = jacobian.values();
    std::fill(values.begin(), values.end(), 0.0);
    assemble(x, residual, values, &step);
    observe_step(step);
    detail::observe("evaluate.x", x);
    detail::observe("evaluate.residual", residual);
    detail::observe("evaluate.jacobian", values);
}

void DriftDiffusion::residual(std::span<const double> x, const TimeStep& step,
                              std::span<double> residual) const {
    assemble(x, residual, {}, &step);
    observe_step(step);
    detail::observe("residual.x", x);
    detail::observe("residual.residual", residual);
}

std::vector<double> DriftDiffusion::storage(std::span<const double> x) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    const std::size_t w = storage_width();
    std::vector<double> s(w * node_count(), 0.0);
    if (electrothermal()) {  // (n, p, tau): complete ionization (DECISIONS.md U6)
        for (std::size_t i = 0; i < node_count(); ++i) {
            s[3 * i + 2] = x[4 * i + 3];
            if (insulator_[i] != 0) continue;
            s[3 * i] = x[4 * i + 1];
            s[3 * i + 1] = x[4 * i + 2];
        }
        detail::observe("storage", s);
        return s;
    }
    const std::vector<NodeDegeneracy> g = degeneracies(x);
    for (std::size_t i = 0; i < node_count(); ++i) {
        if (insulator_[i] != 0) continue;
        const NodeStorage st = node_storage(i, x, g);
        s[w * i] = st.n;
        s[w * i + 1] = st.p;
    }
    detail::observe("storage", s);
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
    detail::observe("trap_occupancies", f);
    return f;
}

std::vector<double> DriftDiffusion::contact_charges(std::span<const double> x,
                                                    const TimeStep* step) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    std::vector<double> charge = gates_.charges(x, m_, contact_count_);
    if (electrothermal()) add_gate_shift(x, charge);
    const auto owner = [&](std::size_t i) {
        return owner_[i] >= 0 ? owner_[i] : electrode_[i];
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
    detail::observe("contact_charges", charge);
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
    detail::observe("conduction_currents", terminal);
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
    // A minority density below eps of the majority at its node does not change the node's charge
    // in double, and the solve resolves its correction only to its own rounding there: measured on
    // Unit 46's long channel, holes at 2.7e-19 of the n+ electrons kept corrections of 1-4e-8 of
    // themselves for up to 14 iterations. With the minority floor on (opt-in: SiC's minorities,
    // near 1e-35 cm^-3, converge quadratically without it and the floor stalls its diode at 3 V)
    // it is measured against eps of the majority (Unit 50, DECISIONS.md RD2); the majority, never
    // below it, as before.
    constexpr double eps = std::numeric_limits<double>::epsilon();
    double size = 0.0;
    for (std::size_t i = 0; i < node_count(); ++i) {
        if (insulator_[i] != 0) {
            size = std::max(size, std::abs(dx[m_ * i]));
            continue;
        }
        const double n = x[m_ * i + 1], p = x[m_ * i + 2];
        const double local = minority_floor_ ? std::max(floor, eps * std::max(n, p)) : floor;
        size = std::max({size, std::abs(dx[m_ * i]),
                         n > 0.0 ? std::abs(dx[m_ * i + 1]) / std::max(n, local) : infinity,
                         p > 0.0 ? std::abs(dx[m_ * i + 2]) / std::max(p, local) : infinity});
    }
    if (electrothermal()) {  // theta, an absolute bound in units of T0 (T10)
        for (std::size_t i = 0; i < node_count(); ++i) {
            size = std::max(size, std::abs(dx[4 * i + 3]));
        }
    }
    detail::observe("update_size", size);
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
    detail::observe("apply_update", x);
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
        detail::observe("edge_currents", currents);
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
    detail::observe("edge_currents", currents);
    return currents;
}

std::vector<double> DriftDiffusion::terminal_currents(std::span<const double> x) const {
    const auto currents = edge_currents(x);
    std::vector<double> terminal(contact_count_, 0.0);
    for (std::size_t k = 0; k < edges_.size(); ++k) {
        const std::int32_t ca = owner_[edges_[k].a], cb = owner_[edges_[k].b];
        if (ca == cb) continue;  // both inside one contact, or both in the device
        const double total = currents[k].first + currents[k].second;
        if (ca >= 0) terminal[static_cast<std::size_t>(ca)] += total;  // leaves contact a
        if (cb >= 0) terminal[static_cast<std::size_t>(cb)] -= total;  // flows into contact b
    }
    // Gate tunnelling (Unit 51, GD7): what a tunnel's node sends into the device through its
    // edges, T_n - T_p, enters through its gate.
    if (!gate_tunnels_.empty()) {
        const std::vector<NodeDegeneracy> g = degeneracies(x);
        std::vector<GateTunnelFlux> fluxes;
        gate_tunnel_fluxes(x, g, fluxes);
        for (std::size_t k = 0; k < gate_tunnels_.size(); ++k) {
            const auto c = static_cast<std::size_t>(gate_tunnels_[k].contact);
            terminal[c] += fluxes[k].T_n - fluxes[k].T_p;
        }
    }
    detail::observe("terminal_currents", terminal);
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
        if (owner_[i] >= 0) {
            bound += eps8 * terms[i];  // the terminal current's own sum
            continue;
        }
        const double recombination = std::abs(flux_n[i] - f[m_ * i + 1]);  // V R
        bound += std::abs(f[m_ * i + 1] + f[m_ * i + 2]) + eps8 * (terms[i] + 2.0 * recombination);
    }
    std::vector<double> resolution(contact_count_, 0.0);
    for (std::size_t c = 0; c < contact_count_; ++c) {
        if (kinds_[c] == device::ContactKind::ohmic || kinds_[c] == device::ContactKind::schottky) {
            resolution[c] = bound;
        }
    }
    detail::observe("terminal_current_resolution", resolution);
    return resolution;
}

std::vector<double> DriftDiffusion::edge_flux_terms(std::span<const double> x) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    std::vector<double> terms(edges_.size(), 0.0);
    if (electrothermal()) {
        const std::vector<ThermalState> states = thermal_states(x);
        const std::vector<NodeLevels> levels = thermal_levels(x, states);
        for (std::size_t k = 0; k < edges_.size(); ++k) {
            const EdgeTerm& e = edges_[k];
            if (!e.carriers) continue;
            const auto [fn, fp] = thermal_fluxes(k, x, states, levels);
            terms[k] = std::abs(fn.d_c[0] * x[4 * e.a + 1]) + std::abs(fn.d_c[1] * x[4 * e.b + 1]) +
                       std::abs(fp.d_c[0] * x[4 * e.a + 2]) + std::abs(fp.d_c[1] * x[4 * e.b + 2]);
        }
        detail::observe("edge_flux_terms", terms);
        return terms;
    }
    const std::vector<NodeDegeneracy> g = degeneracies(x);
    for (std::size_t k = 0; k < edges_.size(); ++k) {
        const EdgeTerm& e = edges_[k];
        if (!e.carriers) continue;
        const auto [fn, fp] = edge_fluxes(e, x, g);
        terms[k] = std::abs(fn.d_c1 * x[m_ * e.a + 1]) + std::abs(fn.d_c2 * x[m_ * e.b + 1]) +
                   std::abs(fp.d_c1 * x[m_ * e.a + 2]) + std::abs(fp.d_c2 * x[m_ * e.b + 2]);
    }
    detail::observe("edge_flux_terms", terms);
    return terms;
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
                b.conduction[i] +
                theta * level(x[4 * i + 1], theta, s.nie, tn.log_nc, s.mass_n, s.d_mass_n).x;
            b.hole_fermi[i] =
                b.valence[i] -
                theta * level(x[4 * i + 2], theta, s.nie, tn.log_nv, s.mass_p, s.d_mass_p).x;
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
    detail::observe("band_edges.conduction", b.conduction);
    detail::observe("band_edges.valence", b.valence);
    detail::observe("band_edges.electron_fermi", b.electron_fermi);
    detail::observe("band_edges.hole_fermi", b.hole_fermi);
    return b;
}

std::vector<double> DriftDiffusion::gate_charges(std::span<const double> x,
                                                 const TimeStep* step) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    std::vector<double> charge = gates_.charges(x, m_, contact_count_);
    if (electrothermal()) add_gate_shift(x, charge);
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
    detail::observe("gate_charges", charge);
    return charge;
}

std::vector<double> DriftDiffusion::interface_trap_charges(std::span<const double> x,
                                                           const TimeStep* step) const {
    NITCAD_EXPECTS(x.size() == unknowns());
    std::vector<double> q(interfaces_.interface_count(), 0.0);
    for (std::size_t k = 0; k < interfaces_.edges().size(); ++k) {
        q[interfaces_.edges()[k].interface] += interface_at(k, x, step).trapped;
    }
    detail::observe("interface_trap_charges", q);
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
        detail::observe("small_signal_matrix", values);
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
        const bool poisson = schottky_[sn] < 0;  // a Schottky node's Poisson row is Dirichlet
        if (poisson) values[so[0]] += t.d_flux_semiconductor[0];
        values[so[1]] -= t.d_rate[0];
        values[so[2]] += t.d_rate_p[0];
        for (std::size_t c = 1; c < 4; ++c) {
            if (poisson) at(sn, 0, c - 1) += t.d_flux_semiconductor[c];
            at(sn, 1, c - 1) -= t.d_rate[c];
            at(sn, 2, c - 1) += t.d_rate_p[c];
        }
    }
    detail::observe("small_signal_matrix", values);
}

std::vector<double> DriftDiffusion::bias_derivative(std::size_t contact) const {
    NITCAD_EXPECTS(contact < contact_count_);
    const auto c = static_cast<std::int32_t>(contact);
    std::vector<double> d(unknowns(), 0.0);
    for (std::size_t i = 0; i < node_count(); ++i) {
        // A Dirichlet row psi - psi0 with psi0 = V / V_T + const; a gate row G (psi_G - psi),
        // psi_G = V / V_T - offset. An ohmic node's densities, and a Schottky node's n_B and
        // p_B, do not depend on the bias.
        if (contact_[i] == c || electrode_[i] == c || schottky_[i] == c) {
            d[m_ * i] = -1.0 / V_T_;
        } else if (gates_.on_gate(i) && gates_.contact(i) == contact) {
            d[m_ * i] = gates_.term(i).coupling / V_T_;
        }
    }
    detail::observe("bias_derivative", d);
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
            const std::int32_t ca = owner_[e.a], cb = owner_[e.b];
            if (ca == cb || !e.carriers) continue;
            SurfaceSlope slope{};
            const auto [fn, fp] = thermal_fluxes(k, x, t, l, &slope);
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
            for (std::size_t j = e.surface_first; j < e.surface_last; ++j) {  // Unit 46
                const SurfaceStencil& st = surface_stencil_[j];
                const double value = (slope.n + slope.p) * surface_gradient(slope.field, st);
                add(ca, 4 * st.node, value);
                add(cb, 4 * st.node, -value);
            }
        }
    }
    const std::vector<NodeDegeneracy> g = degeneracies(x);
    for (const EdgeTerm& e : edges_) {
        const std::int32_t ca = owner_[e.a], cb = owner_[e.b];
        if (electrothermal() || ca == cb || !e.carriers) continue;
        SurfaceSlope slope{};
        const auto [fn, fp] = edge_fluxes(e, x, g, &slope);
        const std::pair<std::size_t, double> d[6] = {
            {m_ * e.a, fn.d_psi1 + fp.d_psi1}, {m_ * e.b, fn.d_psi2 + fp.d_psi2},
            {m_ * e.a + 1, fn.d_c1},           {m_ * e.b + 1, fn.d_c2},
            {m_ * e.a + 2, fp.d_c1},           {m_ * e.b + 2, fp.d_c2}};
        for (const auto& [column, value] : d) {
            add(ca, column, value);
            add(cb, column, -value);
        }
        for (std::size_t k = e.surface_first; k < e.surface_last; ++k) {  // Lombardi (Unit 46)
            const SurfaceStencil& j = surface_stencil_[k];
            const double value = (slope.n + slope.p) * surface_gradient(slope.field, j);
            add(ca, m_ * j.node, value);
            add(cb, m_ * j.node, -value);
        }
    }
    // The traps of an interface edge on an ohmic node (conduction_currents): rate_p - rate; and
    // the charge of a contact node on an interface edge: minus its half-edge flux.
    const auto owner = [&](std::size_t i) {
        return owner_[i] >= 0 ? owner_[i] : electrode_[i];
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
    // Gate tunnelling (Unit 51): the linearized T_n - T_p of each tunnel, a lumped gate's own
    // bias among its columns.
    if (!gate_tunnels_.empty()) {
        GateTunnelFlux flux;
        for (const GateTunnel& t : gate_tunnels_) {
            gate_tunnel_flux(t, x, g, flux);
            for (std::size_t c = 0; c < 3; ++c) {
                add(t.contact, m_ * t.node + c, flux.d_n[c] - flux.d_p[c]);
            }
            for (std::size_t k = 0; k < t.nodes.size(); ++k) {
                add(t.contact, m_ * t.nodes[k], flux.d_n[3 + k] - flux.d_p[3 + k]);
            }
            if (t.lumped) {
                rows[static_cast<std::size_t>(t.contact)].bias += flux.d_n[3] - flux.d_p[3];
            }
        }
    }
    // Displacement: s times the contact charges (contact_charges).
    for (std::size_t i = 0; i < node_count(); ++i) {
        if (!gates_.on_gate(i)) continue;
        const std::size_t c = gates_.contact(i);
        const double G = gates_.term(i).coupling;
        terms[c].emplace_back(m_ * i, -s * G);
        if (electrothermal()) {  // a p-polysilicon gate's potential moves with the node's T
            const double d = electrode_shift(i, thermal_state(i, x[4 * i + 3])).d_theta;
            if (d != 0.0) terms[c].emplace_back(4 * i + 3, s * G * d);
        }
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
    for (std::size_t c = 0; detail::observer != nullptr && c < rows.size(); ++c) {
        std::vector<double> columns;
        for (const std::size_t k : rows[c].columns) columns.push_back(static_cast<double>(k));
        detail::observe("small_signal_currents.columns", columns);
        detail::observe("small_signal_currents.values", rows[c].values);
        detail::observe("small_signal_currents.bias",
                        std::span<const std::complex<double>>(&rows[c].bias, 1));
    }
    return rows;
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
    for (const TunnelPath& p : paths_.paths) layouts_.push_back(layout(p));
    // The base and every term's entries, the paths' among them (PathTunnelling, terms.cpp).
    std::vector<linalg::Triplet> t = base_triplets_;
    Terms::pattern(*this, t);
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
    locate_gate_tunnels();
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
    locate_surface();
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
