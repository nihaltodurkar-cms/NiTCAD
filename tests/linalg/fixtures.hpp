// Matrices shared by the linear-solver tests (Units 3 and 22, and the PARDISO backend of Unit 18):
// every backend is tested on the same systems.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <span>
#include <utility>
#include <vector>

#include "NiTCAD/linalg/sparse_matrix.hpp"

namespace fixtures {

using namespace NiTCAD::linalg;

inline SparseMatrix tridiagonal(Index n, double sub, double diag, double super) {
    std::vector<Triplet> t;
    for (Index i = 0; i < n; ++i) {
        if (i > 0) t.push_back({i, i - 1, sub});
        t.push_back({i, i, diag});
        if (i + 1 < n) t.push_back({i, i + 1, super});
    }
    return SparseMatrix::from_triplets(n, n, t).value();
}

inline SparseMatrix from(Index n, std::vector<Triplet> t) {
    return SparseMatrix::from_triplets(n, n, t).value();
}
inline double max_abs_diff(std::span<const double> a, std::span<const double> b) {
    double m = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) m = std::max(m, std::abs(a[i] - b[i]));
    return m;
}

inline double max_abs(std::span<const double> a) {
    double m = 0.0;
    for (const double v : a) m = std::max(m, std::abs(v));
    return m;
}

// Deterministic value in [0, 1) for integer k (no library RNG, so identical on every build).
inline double hash01(std::uint64_t k) {
    k = (k ^ (k >> 30)) * 0xbf58476d1ce4e5b9ULL;
    k = (k ^ (k >> 27)) * 0x94d049bb133111ebULL;
    k ^= k >> 31;
    return static_cast<double>(k >> 11) * 0x1.0p-53;
}

// A nonsymmetric, diagonally dominant system with a chosen solution; b = A x_exact.
struct System {
    SparseMatrix a;
    std::vector<double> x_exact;
    std::vector<double> b;
};

inline System nonsymmetric(Index n) {
    System s{tridiagonal(n, -1.3, 2.5, -0.7), {}, {}};
    for (Index i = 0; i < n; ++i) s.x_exact.push_back(1.0 + std::sin(0.01 * i));
    s.b.resize(static_cast<std::size_t>(n));
    s.a.multiply(s.x_exact, s.b);
    return s;
}

// 1D box-method Laplacian on n nodes with edge conductances g; Dirichlet rows at both ends add
// 1 to the diagonal. Without them the matrix is singular (pure Neumann, a floating region).
inline SparseMatrix laplace_1d(const std::vector<double>& g, bool dirichlet) {
    const auto n = static_cast<Index>(g.size()) + 1;
    std::vector<Triplet> t;
    for (Index e = 0; e + 1 < n; ++e) {
        const double w = g[static_cast<std::size_t>(e)];
        t.insert(t.end(), {{e, e, w}, {e + 1, e + 1, w}, {e, e + 1, -w}, {e + 1, e, -w}});
    }
    if (dirichlet) t.insert(t.end(), {{0, 0, 1.0}, {n - 1, n - 1, 1.0}});
    return from(n, std::move(t));
}

inline std::vector<double> random_conductances(std::size_t edges) {
    std::vector<double> g(edges);
    for (std::size_t e = 0; e < edges; ++e) g[e] = 1.0 / (0.5 + hash01(e));
    return g;
}

// Three unknowns per node on an m x m grid, neighbours coupled by dense 3x3 blocks (the shape of a
// psi/n/p Jacobian), then rows scaled by 10^[-20, 20] and columns by 10^[-10, 10], as carrier
// densities spanning many decades do. x_exact is chosen and b = A x_exact.
inline System badly_scaled_coupled(int m) {
    const int nodes = m * m;
    const Index n = 3 * nodes;
    std::vector<Triplet> t;
    std::uint64_t k = 0;
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < m; ++i) {
            const Index p = j * m + i;
            const std::array<std::array<int, 2>, 4> nb{{{i - 1, j}, {i + 1, j}, {i, j - 1}, {i, j + 1}}};
            for (const auto& [a, c] : nb) {
                if (a < 0 || a >= m || c < 0 || c >= m) continue;
                const Index q = c * m + a;
                for (Index r = 0; r < 3; ++r) {
                    for (Index s = 0; s < 3; ++s) {
                        const double w = 0.1 + 0.9 * hash01(++k);
                        t.push_back({3 * p + r, 3 * q + s, -w});
                        t.push_back({3 * p + r, 3 * p + s, (r == s ? 1.6 : 0.3) * w});
                    }
                }
            }
            for (Index r = 0; r < 3; ++r) t.push_back({3 * p + r, 3 * p + r, 1.0});
        }
    }
    std::vector<double> row(static_cast<std::size_t>(n)), col(static_cast<std::size_t>(n));
    for (std::size_t i = 0; i < row.size(); ++i) {
        row[i] = std::pow(10.0, 40.0 * hash01(1'000'000 + i) - 20.0);
        col[i] = std::pow(10.0, 20.0 * hash01(2'000'000 + i) - 10.0);
    }
    for (Triplet& e : t) {
        e.value *= row[static_cast<std::size_t>(e.row)] * col[static_cast<std::size_t>(e.col)];
    }
    System s{from(n, std::move(t)), {}, std::vector<double>(static_cast<std::size_t>(n))};
    for (std::size_t i = 0; i < col.size(); ++i) {
        s.x_exact.push_back((1.0 + 0.5 * std::sin(0.1 * static_cast<double>(i))) / col[i]);
    }
    s.a.multiply(s.x_exact, s.b);
    return s;
}

inline double relative_error(std::span<const double> x, std::span<const double> exact) {
    return max_abs_diff(x, exact) / max_abs(exact);
}

// n-node chain with unit conductances, a zeroth-order term delta on every node, and (if w > 0)
// edge (n/2 - 1, n/2) weakened to w and a Dirichlet row at node 0.
inline SparseMatrix anchored_chain(Index n, double delta, double w) {
    std::vector<Triplet> t;
    for (Index i = 0; i + 1 < n; ++i) {
        const double g = (w > 0.0 && i == n / 2 - 1) ? w : 1.0;
        t.insert(t.end(), {{i, i, g}, {i + 1, i + 1, g}, {i, i + 1, -g}, {i + 1, i, -g}});
    }
    for (Index i = 0; i < n; ++i) t.push_back({i, i, delta});
    if (w > 0.0) t.push_back({0, 0, 1.0});
    return from(n, std::move(t));
}

using Complex = std::complex<double>;

inline ComplexSparseMatrix complex_from(Index n, const std::vector<ComplexTriplet>& t) {
    return ComplexSparseMatrix::from_triplets(n, n, t).value();
}

inline double max_abs_diff(std::span<const Complex> a, std::span<const Complex> b) {
    double m = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) m = std::max(m, std::abs(a[i] - b[i]));
    return m;
}

inline double max_abs(std::span<const Complex> a) {
    double m = 0.0;
    for (const Complex v : a) m = std::max(m, std::abs(v));
    return m;
}

inline ComplexSparseMatrix to_complex(const SparseMatrix& a, double imag_scale = 0.0) {
    std::vector<ComplexTriplet> t;
    for (Index r = 0; r < a.rows(); ++r) {
        for (Index k = a.row_offsets()[r]; k < a.row_offsets()[r + 1]; ++k) {
            const double v = a.values()[k];
            t.push_back({r, a.col_indices()[k],
                         imag_scale == 0.0 ? Complex{v, 0.0} : Complex{0.0, imag_scale * v}});
        }
    }
    return complex_from(a.rows(), t);
}

// G + i omega C from a real G and a diagonal C (the AC system's shape).
inline ComplexSparseMatrix ac_matrix(const SparseMatrix& g, const std::vector<double>& c_diag,
                              double omega) {
    std::vector<ComplexTriplet> t;
    for (Index r = 0; r < g.rows(); ++r) {
        for (Index k = g.row_offsets()[r]; k < g.row_offsets()[r + 1]; ++k) {
            const Index col = g.col_indices()[k];
            const double imag = col == r ? omega * c_diag[static_cast<std::size_t>(r)] : 0.0;
            t.push_back({r, col, Complex{g.values()[k], imag}});
        }
    }
    return complex_from(g.rows(), t);
}


}  // namespace fixtures
