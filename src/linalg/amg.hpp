// Smoothed-aggregation algebraic multigrid (Unit 58 SD3) as the preconditioner of the gmres backend
// with preconditioner = amg. Private to linalg; real matrices only (small signal stays on SparseLU,
// LQ3a). Eigen-free: the coarsest level is solved by a SparseLuBackend held through Backend.
//
// The equilibrated matrix R A C (the one ILUT gets) is not symmetric when A is, and Galerkin with
// R = P^T needs a symmetric-like operator: measured, a heterogeneous Poisson problem took 120
// iterations on R A C against 46 unscaled. So the hierarchy is built on the similar matrix
// U (R A C) U^-1 = G A G, U = sqrt(C / R), G = sqrt(R C) (a symmetric scaling of A), and one
// application is M^-1 v = U^-1 V(U v), V the V-cycle on G A G. The near-null vector of G A G is
// the constant of the original unknowns, 1 / G. Without equilibration U = G = 1. Every step below
// is invariant under a symmetric diagonal scaling (DECISIONS.md Unit 58, slice 2).
//
// The hierarchy (Vanek, Mandel, Brezina 1996):
// - strength: j is a strong neighbour of i (i != j, a_ij != 0) if |a_ij| >= theta sqrt(|a_ii a_jj|);
//   the graph is made symmetric (strong one way is strong both ways). A row with no off-diagonal
//   value of its own (a Dirichlet row) is isolated: in no aggregate and no strong neighbour.
// - aggregation, serial in node order (deterministic): (1) a node whose strong neighbours are all
//   free starts an aggregate with them; (2) a free node joins the aggregate (of pass 1) of its
//   strongest aggregated neighbour (by |a_ij| / sqrt(|a_ii a_jj|)), the lowest column on a tie; (3) a node still free starts an
//   aggregate with its free strong neighbours (a node without any is a singleton).
// - tentative prolongator: one column per aggregate, the near-null vector (the constant of the
//   original unknowns, 1 / G) normalized on the aggregate.
// - smoothed prolongator P = (I - omega D^-1 A_F) P0, omega = 4 / (3 lambda), with A_F the filtered
//   matrix (weak entries lumped into the diagonal weighted by the near-null vector b, a_ii +=
//   a_ij b_j / b_i, so A_F b = A b), D its diagonal, and lambda the estimate of the spectral radius
//   of D^-1 A_F: 15 power iterations from a fixed pseudo-random vector, times 1.1.
// - Galerkin coarse operator R A P with R = P^T.
// Coarsening stops at coarse_size unknowns or at max_levels levels, and the coarsest level is then
// factorized by SparseLU; or it stalls (no aggregate, or a level shrinking by less than 10%), and
// that level is then only smoothed (twice), since a level that does not coarsen is one whose
// couplings are weak (measured: the equilibrium Poisson rows, diagonally dominant, which
// otherwise left SparseLU a coarsest level of nearly the whole matrix).
//
// The strength masks and aggregates are built at the first factorization after an analysis and
// kept while the pattern holds, so every later factorization has the same level patterns; the
// operator values are rebuilt from every matrix.
//
// One application of M^-1 is one V-cycle from x = 0: Chebyshev pre- and post-smoothing of degree
// smoother_degree on D^-1 A (D the level matrix's diagonal) over [lambda / 10, lambda], the
// restricted residual, the coarse correction. Every kernel computes each output entry in its own
// fixed order (the products row-parallel on the pool in fixed chunks), so the bits are the same
// for every thread count (6.8, SD2). The setup is serial.
#pragma once

#include <cstddef>
#include <cstdint>
#include <expected>
#include <memory>
#include <span>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"
#include "backend.hpp"

namespace NiTCAD::base {
class ThreadPool;
}

namespace NiTCAD::linalg::detail {

struct AmgSettings {
    double strength_threshold;  // theta, in [0, 1)
    int smoother_degree;        // Chebyshev degree, at least 1
    std::size_t coarse_size;    // a level of at most this many unknowns is the coarsest
    int max_levels;             // at least 2
};

// A real CSR matrix with sorted, unique columns per row.
struct Csr {
    std::size_t rows = 0;
    std::size_t cols = 0;
    std::vector<Index> offsets = std::vector<Index>(1, 0);
    std::vector<Index> columns;
    std::vector<double> values;
    [[nodiscard]] std::size_t nonzeros() const noexcept { return values.size(); }
};

namespace amg {

// Fixed constants of the method (DECISIONS.md Unit 58, slice 2).
inline constexpr std::size_t row_chunk = 2048;       // pool chunk of the row kernels
inline constexpr std::size_t vector_chunk = 8192;    // pool chunk of the vector kernels
inline constexpr double chebyshev_ratio = 10.0;      // lambda_max / lambda_min of the smoother
inline constexpr double prolongator_damping = 4.0 / 3.0;  // over lambda_max(D^-1 A_F)
inline constexpr int power_iterations = 15;          // of the lambda_max estimates
inline constexpr double spectral_safety = 1.1;       // factor on them
inline constexpr double min_coarsening = 0.9;        // a level must shrink below this fraction

// y = A x, row-parallel on pool (nullptr: serial).
void multiply(base::ThreadPool* pool, const Csr& a, std::span<const double> x, std::span<double> y);
// A^T with sorted columns.
[[nodiscard]] Csr transpose(const Csr& a);
// A B (Gustavson); each entry summed in A's row order, then B's. Structural zeros are kept, so the
// pattern depends on the patterns only.
[[nodiscard]] Csr product(const Csr& a, const Csr& b);

// One level's coarsening structure: the strong off-diagonal entries of A (a mask over A's entries)
// and each node's aggregate (-1: none).
struct Coarsening {
    std::vector<std::uint8_t> strong;
    std::vector<Index> aggregate;
    std::size_t aggregates = 0;
};
[[nodiscard]] Coarsening coarsen(const Csr& a, double theta);

// Algebraic symmetrization (Unit 63 MD2): the diagonal D = diag(d) for which D A D^-1 is symmetric
// wherever A allows it. A Scharfetter-Gummel continuity block has a_ij / a_ji = c_i / c_j (c the
// equilibrium density's shape, times the rows' and columns' scales), so d = c^-1/2 makes it exactly
// symmetric, found from the values alone:
//     ln d_j - ln d_i = (ln|a_ij| - ln|a_ji|) / 2
// along a breadth-first spanning forest of the usable off-diagonal pairs (both entries stored,
// non-zero and of one sign), each component started at its lowest unknown and its neighbours taken
// in column order (serial, deterministic). Each component's mean ln d is then shifted to 0 and ln d
// clamped to +-symmetrization_clamp, so a reverse-biased block (ln d spanning hundreds) stays finite
// and its squares do not overflow. An unknown with no usable pair is its own component (d = 1).
inline constexpr double symmetrization_clamp = 230.0;  // |ln d| at most: d within 1e+-100
struct Symmetrization {
    std::vector<double> d;
    std::size_t components = 0;
    std::size_t unpaired = 0;  // off-diagonal entries not usable (partner absent, zero or of the other sign)
    std::size_t clamped = 0;   // unknowns whose ln d was clamped
};
[[nodiscard]] Symmetrization symmetrize(const Csr& a);
// The symmetric part of S = D A D^-1 on A's pattern: (S_ij + S_ji) / 2 on the usable pairs, S_ij
// elsewhere (the diagonal, an absent, zero or opposite-signed partner: a contact row's stored zero
// must not halve its neighbour's coupling, measured on the 3D diode's blocks). Equal to S when exact.
[[nodiscard]] std::vector<double> symmetric_part(const Csr& a, std::span<const double> d);

}  // namespace amg

class AmgBackend final : public Backend<double> {
public:
    AmgBackend(const AmgSettings& settings, base::ThreadPool* pool);
    ~AmgBackend() override;

    // Copies A's pattern as the finest level and forgets the coarsening. Requires a square,
    // non-empty matrix.
    [[nodiscard]] std::expected<void, base::Error> analyze(const SparseMatrix& a) override;
    // Builds the hierarchy from values in the CSR order of the analyzed pattern. A coarsest level
    // SparseLU cannot factorize, or a non-finite operator entry, is inaccurate_solve (the fallback's
    // trigger, as for an ILUT that cannot be built). Reports no pivots; fill is the stored
    // nonzeros of the hierarchy beyond the matrix itself (coarse operators, P, R, the coarse LU).
    [[nodiscard]] std::expected<Factorization, base::Error> factorize(
        std::span<const double> csr_values) override;
    // x = M^-1 b: one V-cycle.
    void solve(std::span<const double> b, std::span<double> x) override;
    void set_scales(std::span<const double> rows, std::span<const double> columns,
                    bool applied) override;

    // The unknowns of every level, finest first (for tests and reports).
    [[nodiscard]] std::vector<std::size_t> level_sizes() const;
    // Level l's matrix, and its prolongator and restrictor (l below the coarsest), for tests.
    [[nodiscard]] const Csr& level_matrix(std::size_t l) const { return levels_.at(l).a; }
    [[nodiscard]] const Csr& prolongator(std::size_t l) const { return levels_.at(l).p; }
    [[nodiscard]] const Csr& restrictor(std::size_t l) const { return levels_.at(l).r; }

private:
    struct Level {
        Csr a;
        std::vector<double> inv_d;  // 1 / diagonal of a (inverse_diagonal in amg.cpp)
        double lambda = 1.0;        // estimated spectral radius of D^-1 a
        Csr p;                      // to this level from the next one
        Csr r;                      // P^T
        std::vector<double> b, x, residual, d;
    };
    void build_level_operators(std::size_t l);
    void smooth(Level& level, bool from_zero);
    void cycle(std::size_t l);

    AmgSettings settings_;
    base::ThreadPool* pool_;
    std::vector<Level> levels_;
    std::vector<amg::Coarsening> coarsening_;  // per level but the coarsest, once per pattern
    bool coarsened_ = false;
    std::vector<double> candidate_;            // near-null vector of the finest level
    std::vector<double> similarity_;           // U = sqrt(C / R)
    std::unique_ptr<Backend<double>> coarse_;
    bool coarse_smoothed_ = false;             // the coarsening stalled: the last level is smoothed
    bool coarse_analyzed_ = false;
    bool built_ = false;
};

// One field's block preconditioner (Unit 63 MD2, MD3): the AMG above on the symmetric part of
// S = D A D^-1 (amg::symmetrize), applied as z = D^-1 V(D r). The AMG gets S with the scales
// R = C = d, which makes its own similarity the identity and its near-null vector 1/d: S's, since
// A's is the equilibrium density c = d^-2 (zero flux) and S (D c) = D A c.
class SymmetrizedAmg {
public:
    SymmetrizedAmg(const AmgSettings& settings, base::ThreadPool* pool);

    // Copies A's pattern (square, non-empty, sorted unique columns).
    [[nodiscard]] std::expected<void, base::Error> analyze(const SparseMatrix& a);
    // d from these values, then the hierarchy on S's symmetric part (AmgBackend::factorize).
    [[nodiscard]] std::expected<Factorization, base::Error> factorize(std::span<const double> csr_values);
    // z = D^-1 V(D r): one V-cycle.
    void solve(std::span<const double> r, std::span<double> z);

    [[nodiscard]] const amg::Symmetrization& symmetrization() const noexcept { return symmetrization_; }
    [[nodiscard]] const AmgBackend& amg() const noexcept { return amg_; }

private:
    base::ThreadPool* pool_;
    AmgBackend amg_;
    Csr a_;
    amg::Symmetrization symmetrization_;
    std::vector<double> b_, x_;
};

}  // namespace NiTCAD::linalg::detail
