// Equilibrium Poisson system: residual F(psi) and Jacobian J(psi) for the scaled potential, one
// unknown per mesh node, carriers slaved to psi by Boltzmann statistics (legacy
// Device1D::solve_equilibrium_boltzmann, generalised from 1D to the mesh graph).
//
// Scaled box-method row of node i not on a contact:
//     F_i = sum over edges e = (i, j):  c_e (psi_j - psi_i)  -  V_i (n_i - p_i - C_i)
// with c_e = coupling area / length / L_D^(D-2), V_i = control volume / L_D^D, n = n_ie e^psi,
// p = n_ie e^-psi (physics::boltzmann_density), C = (N_D - N_A) / Ns and n_ie / Ns from the node's
// material. In 1D this is the legacy row et (psi[i+1] - psi[i]) / h - ... - dV (n - p - C), with
// et = 1 for one material.
// Contact row (ohmic, Dirichlet): F_i = psi_i - psi0_i, psi0 from ohmic_contact_value at zero bias.
//
// The Jacobian pattern is built once (diagonal plus both directions of every edge; contact rows
// keep their off-diagonal entries as explicit zeros), and evaluate() rewrites only the values, so a
// LinearSolver reuses its analysis across Newton iterations (6.10, "Known limits").
#pragma once

#include <cstddef>
#include <expected>
#include <span>
#include <vector>

#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"

namespace NiTCAD::assemble {

class EquilibriumPoisson {
public:
    // Errors (invalid_input):
    // - scaling.temperature_K differs from the device's;
    // - a heterojunction: an edge between regions whose material parameters differ (band offsets
    //   and permittivity steps are deferred); the context index is the edge.
    [[nodiscard]] static std::expected<EquilibriumPoisson, base::Error> create(
        const device::Device& device, const Scaling& scaling);

    [[nodiscard]] std::size_t unknowns() const noexcept { return volume_.size(); }

    // A zero-valued matrix with the Jacobian's pattern, for evaluate() to fill.
    [[nodiscard]] linalg::SparseMatrix make_jacobian() const;

    // Residual and Jacobian at psi. Preconditions (NITCAD_EXPECTS): psi and residual have
    // unknowns() entries; jacobian has the pattern of make_jacobian().
    void evaluate(std::span<const double> psi, std::span<double> residual,
                  linalg::SparseMatrix& jacobian) const;
    // Residual only.
    void residual(std::span<const double> psi, std::span<double> residual) const;

    // The scaled carrier densities slaved to psi, n = n_ie e^psi and p = n_ie e^-psi. Precondition
    // (NITCAD_EXPECTS): all three spans have unknowns() entries.
    void carriers(std::span<const double> psi, std::span<double> n, std::span<double> p) const;

    // The charge-neutral potential of every node, asinh(C / 2 n_ie) (contact nodes: their Dirichlet
    // value, which is the same at zero bias); the legacy initial guess.
    [[nodiscard]] std::vector<double> charge_neutral_potential() const;

    // Per node: the Dirichlet value psi0 on a contact node; not meaningful elsewhere.
    [[nodiscard]] std::span<const double> contact_potential() const noexcept { return psi0_; }
    [[nodiscard]] std::span<const char> is_contact() const noexcept { return contact_; }

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
    std::vector<double> psi0_;      // Dirichlet value on contact nodes
    std::vector<char> contact_;     // 1 on contact nodes
    std::vector<EdgeTerm> edges_;
    std::vector<std::size_t> diag_; // position of (i, i) in the Jacobian values
    linalg::SparseMatrix pattern_;
};

}  // namespace NiTCAD::assemble
