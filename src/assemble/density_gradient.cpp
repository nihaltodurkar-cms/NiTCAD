// Density-gradient quantum confinement (Unit 62; DECISIONS.md QP2-QP5, QP8): the rows of the
// quantum potentials, shared by the drift-diffusion and the thermal-equilibrium tables, and their
// pattern entries and positions.
//
// With c a carrier's density and L = Lambda / kT its quantum potential, Lambda / kT = -beta
// (laplacian sqrt c) / sqrt c. On node i's box, divided by sqrt c_i (the ratio form of Wettstein,
// Schenk and Fichtner, IEEE TED 48 (2001) 279; scaled by L_D):
//
//     V_i L_i + beta sum_j g_ij (sqrt(c_j / c_i) - 1) - beta B_i = 0,
//
// the sum over the node's carrier edges (both ends semiconductor), g the edge's geometry and B_i
// the barrier term of its insulator interfaces and lumped gate (scaled_device.hpp, QP3, QP5).
// Every other boundary is natural. On an ohmic contact node and an insulator node the row is
// L = 0 (QP4). The row is linear in L, and a density at rounding never forms its square root
// alone, only the ratio.
//
// In drift-diffusion c is the unknown density (n, p). At thermal equilibrium the carriers are slaved
// to psi and L (QP8): n = density(eta - L_n) and p = density(-eta - L_p) with eta = psi + s, under
// the selected statistics, so the edge sums read psi and L of both ends.
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

#include "NiTCAD/base/contract.hpp"
#include "scaled_device.hpp"
#include "terms.hpp"

namespace NiTCAD::assemble {

namespace {

using Evaluation = SystemData::Terms::Evaluation;

constexpr FieldId quantum_rows[] = {FieldId::quantum_potential_n, FieldId::quantum_potential_p};
constexpr FieldId quantum_columns[] = {FieldId::potential, FieldId::electrons, FieldId::holes,
                                       FieldId::quantum_potential_n, FieldId::quantum_potential_p};

// Whether node i's quantum potentials are held at 0 (an ohmic contact or an insulator node).
bool quantum_held(const SystemData& s, std::size_t i) noexcept {
    return s.insulator_[i] != 0 || s.contact_[i] >= 0;
}

// sqrt(cb / ca) - 1 without the cancellation of the difference near 1 (two close densities, a fine
// mesh): through log1p and expm1 of the relative difference while it is below 1/2, so the edge sums
// keep their digits as h -> 0 (QV3).
double ratio_minus_one(double ca, double cb) noexcept {
    const double d = cb - ca;
    if (std::abs(d) < 0.5 * ca) return std::expm1(0.5 * std::log1p(d / ca));
    return std::sqrt(cb / ca) - 1.0;
}

}  // namespace

// The node terms: L on a held node, V L - beta B elsewhere (assigned).
class SystemData::Terms::DensityGradientRows final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"density gradient rows", TermDomain::nodes, quantum_rows, quantum_rows, true};
    }
    void add(const SystemData& s, const SystemState&, const Evaluation& ev) const override {
        const std::size_t m = s.m_, sn = s.lambda_n_slot(), sp = s.lambda_p_slot();
        for (std::size_t i = 0; i < s.node_count(); ++i) {
            const double ln = ev.x[m * i + sn], lp = ev.x[m * i + sp];
            if (quantum_held(s, i)) {
                ev.f[m * i + sn] = ln;
                ev.f[m * i + sp] = lp;
                if (ev.jacobian()) at(s, ev, i, sn, sn) = at(s, ev, i, sp, sp) = 1.0;
                continue;
            }
            const double V = s.volume_[i];
            ev.f[m * i + sn] = V * ln - s.dg_beta_n_ * s.dg_barrier_n_[i];
            ev.f[m * i + sp] = V * lp - s.dg_beta_p_ * s.dg_barrier_p_[i];
            if (ev.jacobian()) at(s, ev, i, sn, sn) = at(s, ev, i, sp, sp) = V;
        }
    }
};

// The edge sums beta g (sqrt(c_j / c_i) - 1) into the rows of the ends that are not held (added).
class SystemData::Terms::DensityGradientEdges final : public Term {
public:
    TermInfo info() const noexcept override {
        return {"density gradient edges", TermDomain::mesh_edges, quantum_rows, quantum_columns,
                false};
    }
    // Every entry density gradient adds beyond the node blocks and the edge stencil.
    void pattern(const SystemData& s, const SystemState&,
                 std::vector<linalg::Triplet>& t) const override {
        s.density_gradient_pattern(t);
    }
    void add(const SystemData& s, const SystemState&, const Evaluation& ev) const override {
        if (s.equilibrium_) {
            add_equilibrium(s, ev);
            return;
        }
        const std::size_t m = s.m_;
        const bool jacobian = ev.jacobian();
        for (std::size_t k = 0; k < s.edges_.size(); ++k) {
            const EdgeTerm& e = s.edges_[k];
            if (!e.carriers) continue;
            const DgEdge& d = s.dg_edges_[k];
            const std::size_t a = e.a, b = e.b;
            // Carrier q: its density's slot, its quantum potential's, beta, and its DgEdge entry
            // (row Lambda with the other end's density).
            for (std::size_t q = 0; q < 2; ++q) {
                const std::size_t c = 1 + q, row = 3 + q;
                const double w = (q == 0 ? s.dg_beta_n_ : s.dg_beta_p_) * s.dg_geometry_[k];
                const double ca = ev.x[m * a + c], cb = ev.x[m * b + c];
                const double rm1 = ratio_minus_one(ca, cb), r = 1.0 + rm1;
                if (!quantum_held(s, a)) {
                    ev.f[m * a + row] += w * rm1;
                    if (jacobian) {
                        ev.values[d.ab[q]] += 0.5 * w * r / cb;
                        at(s, ev, a, row, c) -= 0.5 * w * r / ca;
                    }
                }
                if (!quantum_held(s, b)) {
                    ev.f[m * b + row] -= w * rm1 / r;  // 1 / r - 1
                    if (jacobian) {
                        ev.values[d.ba[q]] += 0.5 * w / (r * ca);
                        at(s, ev, b, row, c) -= 0.5 * w / (r * cb);
                    }
                }
            }
        }
    }

private:
    // Equilibrium: the densities slaved to psi and L. Per node and carrier ln c and
    // d ln c / d eta_c (eta_c = eta - L_n for electrons, -eta - L_p for holes).
    static void add_equilibrium(const SystemData& s, const Evaluation& ev) {
        const std::size_t m = s.m_, n = s.node_count();
        std::vector<double> log_c[2], slope[2];
        for (std::size_t q = 0; q < 2; ++q) {
            log_c[q].assign(n, 0.0);
            slope[q].assign(n, 0.0);
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (s.insulator_[i] != 0) continue;
            const double eta = ev.x[m * i] + s.band_shift_[i];
            const double eta_c[2] = {eta - ev.x[m * i + 1], -eta - ev.x[m * i + 2]};
            const double log_dos[2] = {s.log_dos_n_[i], s.log_dos_p_[i]};
            for (std::size_t q = 0; q < 2; ++q) {
                if (!s.fermi_dirac_) {  // n = n_ie e^eta_c, exactly in logarithms
                    log_c[q][i] = std::log(s.n_ie_[i]) + eta_c[q];
                    slope[q][i] = 1.0;
                    continue;
                }
                const physics::DensityResult c =
                    detail::density(true, s.n_ie_[i], log_dos[q], eta_c[q]);
                log_c[q][i] = std::log(c.density);
                slope[q][i] = c.d_eta / c.density;
            }
        }
        const bool jacobian = ev.jacobian();
        for (std::size_t k = 0; k < s.edges_.size(); ++k) {
            const EdgeTerm& e = s.edges_[k];
            if (!e.carriers) continue;
            const DgEdge& d = s.dg_edges_[k];
            const std::size_t a = e.a, b = e.b;
            for (std::size_t q = 0; q < 2; ++q) {
                const std::size_t row = 1 + q;
                // d ln c / d psi = +slope (electrons), -slope (holes); d ln c / d L = -slope.
                const double sign = q == 0 ? 1.0 : -1.0;
                const double w = (q == 0 ? s.dg_beta_n_ : s.dg_beta_p_) * s.dg_geometry_[k];
                const double rm1 = std::expm1(0.5 * (log_c[q][b] - log_c[q][a])), r = 1.0 + rm1;
                if (!quantum_held(s, a)) {
                    ev.f[m * a + row] += w * rm1;
                    if (jacobian) {  // d/d ln c_b = w r / 2, d/d ln c_a = -w r / 2
                        const double hb = 0.5 * w * r * slope[q][b];
                        const double ha = 0.5 * w * r * slope[q][a];
                        ev.values[d.ab[2 * q]] += sign * hb;
                        ev.values[d.ab[2 * q + 1]] -= hb;
                        at(s, ev, a, row, 0) -= sign * ha;
                        at(s, ev, a, row, row) += ha;
                    }
                }
                if (!quantum_held(s, b)) {
                    ev.f[m * b + row] -= w * rm1 / r;  // 1 / r - 1
                    if (jacobian) {  // d/d ln c_a = w / (2 r), d/d ln c_b = -w / (2 r)
                        const double ha = 0.5 * w / r * slope[q][a];
                        const double hb = 0.5 * w / r * slope[q][b];
                        ev.values[d.ba[2 * q]] += sign * ha;
                        ev.values[d.ba[2 * q + 1]] -= ha;
                        at(s, ev, b, row, 0) -= sign * hb;
                        at(s, ev, b, row, row) += hb;
                    }
                }
            }
        }
    }
};

std::shared_ptr<const SystemData::Term> SystemData::Terms::density_gradient_rows() {
    return std::make_shared<DensityGradientRows>();
}

std::shared_ptr<const SystemData::Term> SystemData::Terms::density_gradient_edges() {
    return std::make_shared<DensityGradientEdges>();
}

void SystemData::density_gradient_pattern(std::vector<linalg::Triplet>& t) const {
    const std::size_t m = m_;
    const auto add = [&](std::size_t node, std::size_t r, std::size_t other, std::size_t c) {
        t.push_back({static_cast<linalg::Index>(m * node + r),
                     static_cast<linalg::Index>(m * other + c), 0.0});
    };
    for (const EdgeTerm& e : edges_) {
        if (!e.carriers) continue;
        for (const auto [a, b] : {std::pair{e.a, e.b}, std::pair{e.b, e.a}}) {
            if (equilibrium_) {  // (Lambda_n, psi), (Lambda_n, Lambda_n), (Lambda_p, psi), (Lambda_p, Lambda_p)
                add(a, 1, b, 0);
                add(a, 1, b, 1);
                add(a, 2, b, 0);
                add(a, 2, b, 2);
                continue;
            }
            add(a, 3, b, 1);  // (Lambda_n, n), (Lambda_p, p), (n, Lambda_n), (p, Lambda_p)
            add(a, 4, b, 2);
            add(a, 1, b, 3);
            add(a, 2, b, 4);
            if (dg_impact_) {  // impact ionization: (n, Lambda_p), (p, Lambda_n)
                add(a, 1, b, 4);
                add(a, 2, b, 3);
            }
        }
    }
    if (equilibrium_) return;
    for (const InterfaceEdge& f : interfaces_.edges()) {  // the insulator's Poisson row
        add(f.insulator, 0, f.semiconductor, 3);
        add(f.insulator, 0, f.semiconductor, 4);
    }
}

void SystemData::locate_density_gradient(const linalg::SparseMatrix& pattern) {
    const std::size_t m = m_;
    const auto at = [&](std::size_t node, std::size_t r, std::size_t other, std::size_t c) {
        return detail::position(pattern, m * node + r, m * other + c);
    };
    dg_edges_.assign(edges_.size(), DgEdge{});
    for (std::size_t k = 0; k < edges_.size(); ++k) {
        const EdgeTerm& e = edges_[k];
        if (!e.carriers) continue;
        DgEdge& d = dg_edges_[k];
        for (std::size_t dir = 0; dir < 2; ++dir) {
            const std::size_t a = dir == 0 ? e.a : e.b, b = dir == 0 ? e.b : e.a;
            std::size_t* p = dir == 0 ? d.ab : d.ba;
            if (equilibrium_) {
                p[0] = at(a, 1, b, 0);
                p[1] = at(a, 1, b, 1);
                p[2] = at(a, 2, b, 0);
                p[3] = at(a, 2, b, 2);
                continue;
            }
            p[0] = at(a, 3, b, 1);
            p[1] = at(a, 4, b, 2);
            p[2] = at(a, 1, b, 3);
            p[3] = at(a, 2, b, 4);
        }
    }
    if (equilibrium_) return;
    NITCAD_EXPECTS(impact_nodes_.empty() || dg_impact_);
    for (const ImpactNode& node : impact_nodes_) {
        for (std::size_t k = node.first; k < node.last; ++k) {
            ImpactEdge& ie = impact_edges_[k];
            const EdgeTerm& e = edges_[ie.edge];
            const std::size_t other = e.a == node.node ? e.b : e.a;
            ie.n_lp = at(node.node, 1, other, 4);
            ie.p_ln = at(node.node, 2, other, 3);
        }
    }
    interface_lambda_.clear();
    for (const InterfaceEdge& f : interfaces_.edges()) {
        interface_lambda_.emplace_back(at(f.insulator, 0, f.semiconductor, 3),
                                       at(f.insulator, 0, f.semiconductor, 4));
    }
}

}  // namespace NiTCAD::assemble
