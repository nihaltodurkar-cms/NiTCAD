// CPR-style two-stage preconditioner (Unit 58 SD4, SQ2 a) as the preconditioner of the gmres
// backend with preconditioner = cpr. Private to linalg; real matrices only; Eigen-free (ILUT is held
// through Backend).
//
// The unknowns carry field labels (BasicLinearSolver::set_field_labels): 0 for the potential, a
// positive value for a field coupled to it at the same node, negative for an unknown outside every
// node (a bordering row's). A node is a potential unknown and the positive-labelled unknowns right
// after it (the node-interleaved layout, ARCHITECTURE.md 6.12.1).
//
// The potential system (stage 1) is the matrix A restricted to the potential rows with the carrier
// columns reduced into the potential columns of their node: A_p = R A P, R the injection of the
// potential rows and P the prolongation
//     (P y)_u = e_u y_node(u),  e = 1 on a potential unknown and, on the other unknowns of a node,
//     e_c = -(B_cc^-1 B_cs)_c,
// with B the node's own block of A (s its potential, c the others): the response of the node's
// other unknowns to its potential with their own couplings to the neighbours frozen. For drift-
// diffusion at equilibrium e_n = n and e_p = -p exactly (Scharfetter-Gummel), so A_p is the Jacobian
// of the nonlinear Poisson equation, Laplacian plus the screening term (DECISIONS.md Unit 58,
// slice 3); with A_ψψ alone the screening term is lost. A node whose block B_cc is singular (an
// exactly zero pivot) or gives a non-finite e is left uncoupled (e_c = 0). Unknowns outside every
// node have e = 0.
//
// One application, z = M^-1 v:
//     y = AMG(R v)          one V-cycle on A_p (amg.hpp, the slice-2 hierarchy, its settings)
//     z1 = P y
//     z = z1 + ILUT(v - A z1)   (stage 2: ILUT of the full matrix, Unit 44's settings)
// With an exact stage 1, R (v - A z1) = 0: the potential rows' residual is removed before ILUT.
// The matrix is the equilibrated R A C as for ILUT; A_p is then the equilibrated reduced matrix,
// R_ψ A_p C_ψ, and its AMG is given those scales (set_scales). Every kernel computes each output
// entry in its own fixed order (the product row-parallel on the pool), so the bits are the same for
// every thread count (6.8). The setup is serial.
#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <memory>
#include <span>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"
#include "amg.hpp"
#include "backend.hpp"

namespace NiTCAD::base {
class ThreadPool;
}

namespace NiTCAD::linalg::detail {

class CprBackend final : public Backend<double> {
public:
    CprBackend(const AmgSettings& amg, double drop_tolerance, int fill_factor,
               base::ThreadPool* pool);
    ~CprBackend() override;

    // Copies the labels; the next analyze uses them. True if they differ from the last ones.
    [[nodiscard]] bool set_fields(std::span<const std::int32_t> labels) override;
    // Builds the nodes from the labels, the pattern of A_p and ILUT's ordering. invalid_input if the
    // labels do not have one entry per unknown, a positive label has no potential before it in its
    // node, or there is no potential unknown.
    [[nodiscard]] std::expected<void, base::Error> analyze(const SparseMatrix& a) override;
    // The responses e, A_p, its AMG hierarchy, and ILUT of the full matrix. An AMG or ILUT that
    // cannot be built is inaccurate_solve (the fallback's trigger). fill: the stored nonzeros of
    // A_p, the hierarchy and ILUT's factors.
    [[nodiscard]] std::expected<Factorization, base::Error> factorize(
        std::span<const double> csr_values) override;
    // z = M^-1 v: one two-stage application.
    void solve(std::span<const double> b, std::span<double> x) override;
    void set_scales(std::span<const double> rows, std::span<const double> columns,
                    bool applied) override;

    // For tests: the potential system's matrix (equilibrated values) and size, the responses e of
    // every unknown, and the AMG backend.
    [[nodiscard]] const SparseMatrix& potential_matrix() const noexcept { return reduced_; }
    [[nodiscard]] std::span<const double> responses() const noexcept { return response_; }
    [[nodiscard]] const AmgBackend& amg() const noexcept { return *amg_; }

private:
    base::ThreadPool* pool_;
    std::unique_ptr<AmgBackend> amg_;
    std::unique_ptr<Backend<double>> ilut_;
    std::vector<std::int32_t> labels_;
    // The full matrix (its equilibrated values).
    std::vector<Index> offsets_;
    std::vector<Index> columns_;
    std::vector<double> values_;
    // Nodes: node g's unknowns are potential_[g] .. node_end_[g] - 1; node_of_[u] is u's node or -1.
    std::vector<std::size_t> potential_;
    std::vector<std::size_t> node_end_;
    std::vector<Index> node_of_;
    // A_p, and for every entry of a potential row of A its entry of A_p (-1: none).
    SparseMatrix reduced_;
    std::vector<Index> target_;
    std::vector<double> response_;
    std::vector<double> rows_p_, cols_p_;  // the scales of the potential unknowns
    bool scaled_ = false;
    // Work: the restricted residual, the AMG correction, z1, the stage-2 residual and correction.
    std::vector<double> v_p_, y_p_, z1_, r2_, z2_;
    std::vector<double> block_;  // a node's B_cc and right-hand side (factorize)
    bool analyzed_ = false;
};

}  // namespace NiTCAD::linalg::detail
