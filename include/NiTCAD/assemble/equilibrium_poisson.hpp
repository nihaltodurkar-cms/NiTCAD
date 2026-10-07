// Equilibrium Poisson system: residual F(psi) and Jacobian J(psi) for the scaled potential, one
// unknown per mesh node, carriers slaved to psi by Boltzmann or, with models.fermi_dirac,
// Fermi-Dirac statistics (legacy Device1D::solve_equilibrium_boltzmann and
// solve_equilibrium_fd, generalised from 1D to the mesh graph).
//
// Scaled box-method row of node i not on a contact:
//     F_i = sum over edges e = (i, j):  c_e (psi_j - psi_i)  -  V_i (n_i - p_i - C_i)
// with c_e = coupling area / length / L_D^(D-2) times the edge's permittivity over the reference
// one (the harmonic mean of its ends'), V_i = control volume / L_D^D, n = n_ie e^(psi + s),
// p = n_ie e^-(psi + s) (physics::boltzmann_density; Fermi-Dirac: physics::fermi_dirac_density with
// eta = psi + s and -(psi + s)), C = (N_D - N_A) / Ns and n_ie / Ns from the node's material (with
// band-gap narrowing when models.bgn is set), and s the node's band shift (Unit 15; 0 for a single
// material): the intrinsic-level depth of its material below the vacuum level, less node 0's, over
// V_T. psi is then the electrostatic potential, continuous across a heterointerface. In 1D this is
// the legacy row et (psi[i+1] - psi[i]) / h - ... - dV (n - p - C) (with the M33 affinity shift).
// With models.incomplete_ionization C is the ionized doping N_D+ - N_A- (physics/ionization.hpp)
// at eta_c = psi + s - g_n and eta_v = -(psi + s) - g_p, with its psi derivative.
// Contact row (ohmic, Dirichlet): F_i = psi_i - psi0_i, psi0 from ohmic_contact_value at zero bias,
// less s.
// A gate node keeps its box row and gains the oxide term of gate.hpp, G_i (psi_G,i - psi_i) + S_i,
// at the gate's bias (0 unless set_bias says otherwise; legacy Device2D.solve_equilibrium). A gate
// carries no current, so a biased gate keeps the device in thermal equilibrium as long as every
// ohmic contact is at 0 V: the legacy MOS-C solve (moscap.MOSCapacitor.solve_psi), whose C-V is
// the quasi-static one.
// Insulators (Unit 15b): an insulator node's row is the box row with no charge (its carrier
// densities are 0), and an electrode node's (an insulator node) is Dirichlet, F_i = psi_i - psi_E
// with psi_E the electrode potential at its bias (scaled_device.hpp). A semiconductor node at a
// semiconductor-insulator interface gains the interface term of interface_nodes.hpp: the fixed
// charge and the traps' charge with Fermi occupancy at eta = psi + s.
//
// The Jacobian pattern is built once (diagonal plus both directions of every edge; contact rows
// keep their off-diagonal entries as explicit zeros), and evaluate() rewrites only the values, so a
// LinearSolver reuses its analysis across Newton iterations (6.10, "Known limits").
#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <vector>

#include "NiTCAD/assemble/band_edges.hpp"
#include "NiTCAD/assemble/gate.hpp"
#include "NiTCAD/assemble/interface_nodes.hpp"
#include "NiTCAD/assemble/models.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"
#include "NiTCAD/physics/ionization.hpp"

namespace NiTCAD::assemble {

class EquilibriumPoisson {
public:
    // Errors (invalid_input): scaling.temperature_K differs from the device's.
    // Of the models only `bgn` (the effective n_ie of every node), `fermi_dirac` and
    // `incomplete_ionization` matter here.
    [[nodiscard]] static std::expected<EquilibriumPoisson, base::Error> create(
        const device::Device& device, const Scaling& scaling, const PhysicsModels& models = {});

    [[nodiscard]] std::size_t unknowns() const noexcept { return volume_.size(); }

    // Applied bias of each contact in V, in the order of device.contacts(). Errors: those of
    // check_contact_bias with thermal_equilibrium (an ohmic contact not at 0 V among them). Sets
    // the electrode potential of the gates; initially every bias is 0.
    [[nodiscard]] std::expected<void, base::Error> set_bias(std::span<const double> bias_V);

    // A zero-valued matrix with the Jacobian's pattern, for evaluate() to fill.
    [[nodiscard]] linalg::SparseMatrix make_jacobian() const;

    // Residual and Jacobian at psi. Preconditions (NITCAD_EXPECTS): psi and residual have
    // unknowns() entries; jacobian has the pattern of make_jacobian().
    void evaluate(std::span<const double> psi, std::span<double> residual,
                  linalg::SparseMatrix& jacobian) const;
    // Residual only.
    void residual(std::span<const double> psi, std::span<double> residual) const;

    // The scaled carrier densities slaved to psi (n = n_ie e^(psi + s) and p = n_ie e^-(psi + s)
    // under Boltzmann statistics; 0 on an insulator node). Precondition (NITCAD_EXPECTS): all three
    // spans have unknowns() entries.
    void carriers(std::span<const double> psi, std::span<double> n, std::span<double> p) const;

    // The charge-neutral potential of every node, asinh(C / 2 n_ie) - s under Boltzmann
    // statistics, the root of physics::fermi_dirac_neutral_equilibrium less s under Fermi-Dirac
    // (ohmic contact nodes: their Dirichlet value, which is the same at zero bias); the legacy
    // initial guess. An electrode node takes its Dirichlet value at the current bias, another
    // insulator node 0.
    [[nodiscard]] std::vector<double> charge_neutral_potential() const;

    // The band diagram at psi (band_edges.hpp); the quasi-Fermi levels are 0 (equilibrium). NaN on
    // insulator nodes, which have no bands in the model. Precondition (NITCAD_EXPECTS): psi has
    // unknowns() entries.
    [[nodiscard]] BandEdges band_edges(std::span<const double> psi) const;

    // Per node: the Dirichlet value psi0 on an ohmic contact or electrode node (the electrode's
    // at its current bias); not meaningful elsewhere.
    [[nodiscard]] std::span<const double> contact_potential() const noexcept { return psi0_; }
    // Per node: 1 on an ohmic contact or electrode node (a Dirichlet row), else 0; gate nodes are
    // 0.
    [[nodiscard]] std::span<const char> is_contact() const noexcept { return contact_; }

    // Scaled charge on each gate and electrode, as DriftDiffusion::gate_charges (zero for an ohmic
    // contact). Precondition (NITCAD_EXPECTS): psi has unknowns() entries.
    [[nodiscard]] std::vector<double> gate_charges(std::span<const double> psi) const;
    // Scaled trapped charge of each declared interface (interface_nodes.hpp, without the fixed
    // charge; zero for an interface without traps), in units of q Ns L_D^D. Precondition
    // (NITCAD_EXPECTS): psi has unknowns() entries.
    [[nodiscard]] std::vector<double> interface_trap_charges(std::span<const double> psi) const;

private:
    EquilibriumPoisson() = default;

    struct EdgeTerm {
        std::size_t i, j;          // end nodes
        double c;                  // scaled coupling
        std::size_t ij, ji;        // positions of (i, j) and (j, i) in the Jacobian values
    };

    // Writes the residual and, if jacobian_values is not empty, the charge derivative onto the
    // diagonal of non-contact rows.
    void residual_into(std::span<const double> psi, std::span<double> residual,
                       std::span<double> jacobian_values) const;

    std::vector<double> volume_;    // scaled control volume
    std::vector<double> doping_;    // scaled net doping C
    std::vector<double> n_ie_;      // scaled n_ie
    std::vector<double> log_dos_n_, log_dos_p_;  // ln(Nc / n_ie), ln(Nv / n_ie)
    std::vector<double> band_shift_;             // s
    std::vector<double> donors_, acceptors_;     // N_D / Ns, N_A / Ns
    std::vector<physics::DopantLevels> levels_;  // dopant levels in units of kT
    bool fermi_dirac_ = false;
    bool ionization_ = false;
    std::vector<double> psi0_;      // Dirichlet value on contact nodes
    std::vector<char> contact_;     // 1 on ohmic contact and electrode nodes
    std::vector<char> insulator_;   // 1 on insulator nodes
    std::vector<std::int32_t> electrode_;      // electrode contact per node, or -1
    std::vector<double> electrode_potential_;  // psi_E at zero bias
    double V_T_ = 1.0;
    GateNodes gates_;
    InterfaceNodes interfaces_;
    std::vector<device::ContactKind> kinds_;  // per contact
    std::vector<EdgeTerm> edges_;
    std::vector<std::size_t> diag_; // position of (i, i) in the Jacobian values
    linalg::SparseMatrix pattern_;
};

}  // namespace NiTCAD::assemble
