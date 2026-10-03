// Coupled drift-diffusion system: Poisson, electron and hole continuity, three scaled unknowns per
// node interleaved as x[3i] = psi, x[3i+1] = n, x[3i+2] = p (legacy Device1D::residual_jacobian,
// baseline models, generalised from 1D to the mesh graph).
//
// Rows of node i not on a contact, summing over its edges e = (a, b) (+ at a, - at b):
//     Poisson    sum  et c_e (psi_b - psi_a)                    - V_i (n_i - p_i - C_i)
//     electrons  sum  Jn_e                                        - V_i R_i
//     holes      sum  Jp_e                                        + V_i R_i
// with c_e the scaled coupling (EquilibriumPoisson), Jn_e and Jp_e the Scharfetter-Gummel fluxes of
// sg_flux.hpp with edge factor D_e c_e, D_e = hmean(mu_a, mu_b) V_T / D0 (harmonic mean of the
// nodes' mobilities, legacy dn_edge), and R = (R_SRH + R_Auger) / R0. SRH is evaluated on scaled
// densities (it is homogeneous of degree one, so R_SRH / R0 = R_SRH(n', p', ...) Ns / R0); Auger,
// which is cubic, on physical ones.
// With models.field_mobility the edge factor carries the Canali mobility of the edge's own field:
// D_e = mu_C(mu0_e, E_e) V_T / D0, with mu0_e the harmonic mean above and
// E_e = V_T |psi_b - psi_a| / length the electrostatic field along the edge (the field parallel to
// the edge's current, as in the legacy). The Jacobian includes d D_e / d psi, so Newton stays
// quadratic (OLD / NEW / REASON in ARCHITECTURE.md 6.2, Unit 13: the legacy averaged |E| to the
// nodes, applied Canali per node, took the harmonic mean, and lagged it).
// The driving term of the fluxes is delta_n = psi_b - psi_a + ln(n_ie,b / n_ie,a) for electrons and
// delta_p = psi_b - psi_a - ln(n_ie,b / n_ie,a) for holes (legacy delta, delta_p): with band-gap
// narrowing n_ie varies in space, and only these make the equilibrium carry no current.
// Contact rows (ohmic): psi - psi0, n - n0, p - p0, from ohmic_contact_value at the contact's bias.
// A gate node keeps its three box rows; its Poisson row gains the oxide term of gate.hpp at the
// gate's bias, G_i (psi_G,i - psi_i) + S_i, and its continuity rows no boundary flux (legacy
// Device2D GateBC: Robin on psi only).
//
// Jn and Jp are the electron and hole current densities in units of J0 (conventional current, along
// the edge from a to b); their sum is divergence-free at every node off the ohmic contacts (gate
// nodes included), so the terminal current of an ohmic contact is the total flux on the edges
// leaving it, and that of a gate is zero.
//
// Newton update (legacy Device1D::newton): psi is clipped to +-max_update; n and p are clamped to
// [0.1, 10] times their current value, which keeps them positive. The convergence measure is the
// largest of |dpsi|, |dn| / n and |dp| / p of the FULL correction (6.2; the legacy measured the
// damped one). Only during iteration: the converged state is never clamped (6.7).
#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/gate.hpp"
#include "NiTCAD/assemble/models.hpp"
#include "NiTCAD/assemble/sg_flux.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::assemble {

class DriftDiffusion {
public:
    // Errors: those of EquilibriumPoisson::create (scaling temperature, heterojunctions).
    [[nodiscard]] static std::expected<DriftDiffusion, base::Error> create(
        const device::Device& device, const Scaling& scaling,
        const PhysicsModels& models = {});

    [[nodiscard]] std::size_t unknowns() const noexcept { return 3 * node_count(); }
    [[nodiscard]] std::size_t node_count() const noexcept { return volume_.size(); }
    [[nodiscard]] std::size_t contact_count() const noexcept { return contact_count_; }

    // Applied bias of each contact in V, in the order of device.contacts(); sets the Dirichlet
    // values of the ohmic contacts and the electrode potential of the gates. Errors: those of
    // check_contact_bias (size, a value not finite). Initially every bias is 0.
    [[nodiscard]] std::expected<void, base::Error> set_bias(std::span<const double> bias_V);

    // Scaled state of every node in Boltzmann equilibrium at the scaled potential psi
    // (n = n_ie e^psi, p = n_ie e^-psi), with the ohmic contact nodes set to their Dirichlet
    // values. Precondition (NITCAD_EXPECTS): psi has node_count() entries.
    [[nodiscard]] std::vector<double> state_from_potential(std::span<const double> psi) const;
    // Sets the ohmic contact nodes of a state to their Dirichlet values.
    void stamp_contacts(std::span<double> x) const;

    [[nodiscard]] linalg::SparseMatrix make_jacobian() const;
    // Preconditions (NITCAD_EXPECTS): x and residual have unknowns() entries; jacobian has the
    // pattern of make_jacobian().
    void evaluate(std::span<const double> x, std::span<double> residual,
                  linalg::SparseMatrix& jacobian) const;
    void residual(std::span<const double> x, std::span<double> residual) const;

    // Newton hooks (see the header comment).
    [[nodiscard]] double update_size(std::span<const double> x,
                                     std::span<const double> dx) const;
    void apply_update(std::span<double> x, std::span<const double> dx, double max_update) const;

    // Scaled (Jn, Jp) of every mesh edge, in units of J0 times the edge's scaled coupling area
    // (so a 1D edge gives the current density in units of J0).
    [[nodiscard]] std::vector<std::pair<double, double>> edge_currents(
        std::span<const double> x) const;
    // Scaled total current entering the device through each contact (conventional current from
    // the contact into the device), in device.contacts() order; zero for a gate. Physical value:
    // times J0 L_D^(D-1), in A / cm^(3-D) (A/cm^2 in 1D, A/cm in 2D, A in 3D).
    [[nodiscard]] std::vector<double> terminal_currents(std::span<const double> x) const;
    // Scaled charge on each gate electrode, the sum over its nodes of G_i (psi_G,i - psi_i)
    // (gate.hpp), in device.contacts() order; zero for an ohmic contact. Physical value: times
    // q Ns L_D^D, in C / cm^(3-D).
    [[nodiscard]] std::vector<double> gate_charges(std::span<const double> x) const;

private:
    DriftDiffusion() = default;

    struct EdgeTerm {
        std::size_t a, b;   // end nodes
        double c;           // et * scaled coupling, Poisson
        double an, ap;      // electron and hole SG edge factors (low field)
        double mu_n, mu_p;  // low-field edge mobilities (harmonic means) [cm^2/(V s)]
        double field;       // [V/cm] per unit |psi_b - psi_a|: V_T / length
        physics::CanaliParameters sat_n, sat_p;  // of the edge's material
        double dln;         // ln(n_ie,b / n_ie,a)
        // Positions in the Jacobian values: row a with columns of b, row b with columns of a,
        // in the order (psi, psi), (n, psi), (n, n), (p, psi), (p, p).
        std::size_t ab[5], ba[5];
    };

    void assemble(std::span<const double> x, std::span<double> residual,
                  std::span<double> values) const;

    // The electron and hole fluxes of an edge at state x, with field-dependent edge factors when
    // the model is on.
    [[nodiscard]] std::pair<EdgeFlux, EdgeFlux> edge_fluxes(
        const EdgeTerm& e, std::span<const double> x) const;

    std::vector<double> volume_, doping_, n_ie_, tau_n_, tau_p_, auger_n_, auger_p_;
    std::vector<std::int32_t> contact_;   // ohmic contact index per node, or -1
    GateNodes gates_;
    std::vector<device::ContactKind> kinds_;  // per contact
    std::vector<double> psi0_, n0_, p0_;  // Dirichlet values per node (contact nodes only)
    std::vector<EdgeTerm> edges_;
    std::vector<std::size_t> block_;      // per node, 9 positions of its 3x3 diagonal block
    std::size_t contact_count_ = 0;
    double rate_scale_ = 0.0;             // Ns / R0: scaled SRH call to R / R0
    double Ns_ = 0.0;
    double V_T_ = 0.0;
    bool srh_ = true;
    bool auger_ = true;
    bool field_mobility_ = false;
    linalg::SparseMatrix pattern_;
};

}  // namespace NiTCAD::assemble
