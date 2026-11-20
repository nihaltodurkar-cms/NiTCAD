// The cpr preconditioner of the gmres backend (Unit 58 SD4, slice 3): the local carrier responses
// and the reduced potential system against an independent Scharfetter-Gummel derivation, the
// configuration and label checks, GMRES-CPR on drift-diffusion Jacobians through the public
// interface (an independent backward error), the bordering unknown, the reuse of the analysis, the
// fallback, and bit-identical results for every thread count (SG2).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "../../src/linalg/cpr.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"
#include "fixtures.hpp"

using namespace NiTCAD::linalg;
using NiTCAD::base::ErrorCode;
using namespace fixtures;

namespace {

constexpr SolverConfig cpr_strict{.backend = SolverBackend::gmres_ilut,
                                  .fallback = LinearFallback::none,
                                  .preconditioner = Preconditioner::cpr};

// Bernoulli function B(x) = x / (e^x - 1) and its derivative.
double bernoulli(double x) { return std::abs(x) < 1e-8 ? 1.0 - x / 2.0 : x / std::expm1(x); }
double bernoulli_derivative(double x) {
    if (std::abs(x) < 1e-5) return -0.5 + x / 6.0;
    const double e = std::expm1(x);
    return (e - x * (e + 1.0)) / (e * e);
}

// The drift-diffusion Jacobian (psi, n, p per node, node-interleaved, de Mari units, unit boxes) of
// a 2D m x m grid with a pn junction at mid x: Poisson lambda2 sum (psi_i - psi_j) + n - p - C,
// Scharfetter-Gummel continuity rows. Contacts at x = 0 and x = m - 1 are identity rows for every
// unknown. The state: psi the charge-neutral potential, n = exp(psi - phi_n), p = exp(-psi), with
// phi_n rising linearly by bias across x (0: thermal equilibrium, every flux zero).
struct DriftDiffusion {
    SparseMatrix a;
    std::vector<double> n, p;      // per node
    std::vector<int> contact;      // per node
    double lambda2 = 0.0;
    Index m = 0;
};

DriftDiffusion drift_diffusion(Index m, double doping, double lambda2, double bias,
                               bool heterogeneous = false) {
    DriftDiffusion d;
    d.m = m;
    d.lambda2 = lambda2;
    const auto nodes = static_cast<std::size_t>(m * m);
    std::vector<double> psi(nodes);
    d.n.resize(nodes);
    d.p.resize(nodes);
    d.contact.assign(nodes, 0);
    for (Index j = 0; j < m; ++j) {
        for (Index i = 0; i < m; ++i) {
            const auto k = static_cast<std::size_t>(j * m + i);
            const double c = 2 * i < m ? doping : -doping;
            psi[k] = std::asinh(c / 2.0);
            const double phi_n = bias * static_cast<double>(i) / static_cast<double>(m - 1);
            d.n[k] = std::exp(psi[k] - phi_n);
            d.p[k] = std::exp(-psi[k]);
            d.contact[k] = (i == 0 || i + 1 == m) ? 1 : 0;
        }
    }
    std::vector<Triplet> t;
    const auto at = [](std::size_t node, int field) { return static_cast<Index>(3 * node + static_cast<std::size_t>(field)); };
    for (Index j = 0; j < m; ++j) {
        for (Index i = 0; i < m; ++i) {
            const auto k = static_cast<std::size_t>(j * m + i);
            if (d.contact[k]) {
                for (int f = 0; f < 3; ++f) t.push_back({at(k, f), at(k, f), 1.0});
                continue;
            }
            t.push_back({at(k, 1), at(k, 2), 0.0});  // explicit zeros of the node block
            t.push_back({at(k, 2), at(k, 1), 0.0});
            t.push_back({at(k, 0), at(k, 1), 1.0});
            t.push_back({at(k, 0), at(k, 2), -1.0});
            const std::pair<Index, Index> neighbours[] = {{i - 1, j}, {i + 1, j}, {i, j - 1}, {i, j + 1}};
            for (const auto& [ni, nj] : neighbours) {
                if (ni < 0 || nj < 0 || ni >= m || nj >= m) continue;
                const auto l = static_cast<std::size_t>(nj * m + ni);
                const double delta = psi[l] - psi[k];
                // heterogeneous: permittivities per edge from 2^-10 to 2^10 (symmetric in k, l).
                const double g = heterogeneous ? std::exp2(std::floor(20.0 * hash01(std::min(k, l) * 7 + std::max(k, l))) - 10.0) : 1.0;
                t.push_back({at(k, 0), at(k, 0), lambda2 * g});
                t.push_back({at(k, 0), at(l, 0), -lambda2 * g});
                const double bm = bernoulli(-delta), bp = bernoulli(delta);
                const double dm = bernoulli_derivative(-delta), dp = bernoulli_derivative(delta);
                // Electrons: n_i B(-delta) - n_j B(delta).
                t.push_back({at(k, 1), at(k, 1), bm});
                t.push_back({at(k, 1), at(l, 1), -bp});
                const double gn = d.n[k] * dm + d.n[l] * dp;
                t.push_back({at(k, 1), at(k, 0), gn});
                t.push_back({at(k, 1), at(l, 0), -gn});
                // Holes: p_i B(delta) - p_j B(-delta).
                t.push_back({at(k, 2), at(k, 2), bp});
                t.push_back({at(k, 2), at(l, 2), -bm});
                const double gp = d.p[k] * dp + d.p[l] * dm;
                t.push_back({at(k, 2), at(k, 0), -gp});
                t.push_back({at(k, 2), at(l, 0), gp});
            }
        }
    }
    d.a = from(static_cast<Index>(3 * nodes), t);
    return d;
}

std::vector<std::int32_t> labels_for(std::size_t nodes, std::size_t extra = 0) {
    std::vector<std::int32_t> l;
    for (std::size_t k = 0; k < nodes; ++k) {
        for (std::int32_t f = 0; f < 3; ++f) l.push_back(f);
    }
    l.insert(l.end(), extra, -1);
    return l;
}

std::vector<double> rhs(std::size_t n) {
    std::vector<double> b(n);
    for (std::size_t i = 0; i < n; ++i) b[i] = std::sin(0.37 * static_cast<double>(i)) + 1e-3;
    return b;
}

// Componentwise backward error, written independently of the solver's.
double backward_error(const SparseMatrix& a, std::span<const double> b, std::span<const double> x) {
    double worst = 0.0;
    for (Index i = 0; i < a.rows(); ++i) {
        double ax = 0.0;
        double abs_ax = 0.0;
        for (Index k = a.row_offsets()[i]; k < a.row_offsets()[i + 1]; ++k) {
            const double term = a.values()[k] * x[static_cast<std::size_t>(a.col_indices()[k])];
            ax += term;
            abs_ax += std::abs(term);
        }
        const auto u = static_cast<std::size_t>(i);
        const double denominator = abs_ax + std::abs(b[u]);
        worst = std::max(worst, denominator > 0.0 ? std::abs(b[u] - ax) / denominator : 0.0);
    }
    return worst;
}

struct Solved {
    std::vector<double> x;
    SolveReport report;
};

Solved solve_with(const SparseMatrix& a, std::span<const double> b, const SolverConfig& config,
                  std::span<const std::int32_t> labels) {
    auto solver = LinearSolver::create(config).value();
    solver.set_field_labels(labels);
    const auto f = solver.factorize(a);
    INFO("factorize: " << (f ? std::string("ok") : f.error().message));
    REQUIRE(f.has_value());
    Solved s{std::vector<double>(b.size()), {}};
    const auto r = solver.solve(b, s.x);
    INFO("solve: " << (r ? std::string("ok") : r.error().message));
    REQUIRE(r.has_value());
    s.report = *r;
    REQUIRE(backward_error(a, b, s.x) <= config.max_backward_error);
    return s;
}

}  // namespace

TEST_CASE("cpr: the configuration and label checks") {
    SolverConfig sparse_lu_cpr{.preconditioner = Preconditioner::cpr};
    const auto direct = LinearSolver::create(sparse_lu_cpr);
    REQUIRE_FALSE(direct.has_value());
    REQUIRE(direct.error().code == ErrorCode::invalid_input);
    const auto complex = ComplexLinearSolver::create(cpr_strict);
    REQUIRE_FALSE(complex.has_value());
    REQUIRE(complex.error().code == ErrorCode::invalid_input);

    const DriftDiffusion d = drift_diffusion(6, 100.0, 0.1, 0.0);
    const std::size_t nodes = 36;
    {
        auto solver = LinearSolver::create(cpr_strict).value();
        const auto f = solver.factorize(d.a);  // no labels
        REQUIRE_FALSE(f.has_value());
        REQUIRE(f.error().code == ErrorCode::invalid_input);
    }
    {
        auto solver = LinearSolver::create(cpr_strict).value();
        std::vector<std::int32_t> short_labels = labels_for(nodes);
        short_labels.pop_back();
        solver.set_field_labels(short_labels);
        REQUIRE(solver.factorize(d.a).error().code == ErrorCode::invalid_input);
    }
    {
        auto solver = LinearSolver::create(cpr_strict).value();
        std::vector<std::int32_t> labels = labels_for(nodes);
        labels[0] = 1;  // a coupled field before any potential
        solver.set_field_labels(labels);
        REQUIRE(solver.factorize(d.a).error().code == ErrorCode::invalid_input);
        labels = labels_for(nodes);
        labels[3] = -1;  // node 1's potential outside every node: its carriers have none before them
        solver.set_field_labels(labels);
        REQUIRE(solver.factorize(d.a).error().code == ErrorCode::invalid_input);
        solver.set_field_labels(std::vector<std::int32_t>(3 * nodes, -1));  // no potential
        REQUIRE(solver.factorize(d.a).error().code == ErrorCode::invalid_input);
        solver.set_field_labels(labels_for(nodes));
        REQUIRE(solver.factorize(d.a).has_value());
    }
    // Labels are ignored by every other configuration.
    auto ilut = LinearSolver::create({.backend = SolverBackend::gmres_ilut}).value();
    ilut.set_field_labels(std::vector<std::int32_t>(1, 7));
    REQUIRE(ilut.factorize(d.a).has_value());
}

TEST_CASE("cpr: at equilibrium the responses are n and -p, and A_p is the nonlinear Poisson Jacobian") {
    const Index m = 12;
    const double lambda2 = 0.05;
    const DriftDiffusion d = drift_diffusion(m, 1e4, lambda2, 0.0);
    const auto nodes = static_cast<std::size_t>(m * m);
    NiTCAD::linalg::detail::CprBackend cpr({.strength_threshold = 0.08, .smoother_degree = 2,
                                            .coarse_size = 1000, .max_levels = 10},
                                           1e-4, 10, nullptr);
    (void)cpr.set_fields(labels_for(nodes));
    REQUIRE(cpr.analyze(d.a).has_value());
    const std::vector<double> ones(3 * nodes, 1.0);
    cpr.set_scales(ones, ones, false);
    REQUIRE(cpr.factorize(d.a.values()).has_value());
    const std::span<const double> e = cpr.responses();
    const SparseMatrix& ap = cpr.potential_matrix();
    REQUIRE(ap.rows() == static_cast<Index>(nodes));
    for (Index j = 0; j < m; ++j) {
        for (Index i = 0; i < m; ++i) {
            const auto k = static_cast<std::size_t>(j * m + i);
            CAPTURE(i, j);
            REQUIRE(e[3 * k] == 1.0);
            const double en = d.contact[k] ? 0.0 : d.n[k];
            const double ep = d.contact[k] ? 0.0 : -d.p[k];
            REQUIRE(std::abs(e[3 * k + 1] - en) <= 1e-12 * std::max(1.0, d.n[k]));
            REQUIRE(std::abs(e[3 * k + 2] - ep) <= 1e-12 * std::max(1.0, d.p[k]));
            // Row k of A_p: lambda2 times the Laplacian, plus n + p on the diagonal.
            const int degree = (i > 0) + (i + 1 < m) + (j > 0) + (j + 1 < m);
            for (Index q = ap.row_offsets()[static_cast<Index>(k)];
                 q < ap.row_offsets()[static_cast<Index>(k) + 1]; ++q) {
                const auto c = static_cast<std::size_t>(ap.col_indices()[q]);
                double expected = 0.0;
                if (d.contact[k]) {
                    expected = c == k ? 1.0 : 0.0;
                } else if (c == k) {
                    expected = lambda2 * degree + d.n[k] + d.p[k];
                } else {
                    expected = -lambda2;
                }
                REQUIRE(std::abs(ap.values()[q] - expected) <= 1e-12 * std::max(1.0, std::abs(expected)));
            }
            REQUIRE(ap.row_offsets()[static_cast<Index>(k) + 1] - ap.row_offsets()[static_cast<Index>(k)] ==
                    (d.contact[k] ? 1 : degree + 1));
        }
    }
}

TEST_CASE("cpr: GMRES-CPR solves drift-diffusion Jacobians, equilibrium and biased") {
    for (const Index m : {Index{24}, Index{48}, Index{96}}) {
        for (const double bias : {0.0, 20.0}) {
            const DriftDiffusion d = drift_diffusion(m, 1e4, 0.05, bias);
            const auto nodes = static_cast<std::size_t>(m * m);
            const std::vector<double> b = rhs(3 * nodes);
            CAPTURE(m, bias);
            const Solved cpr = solve_with(d.a, b, cpr_strict, labels_for(nodes));
            REQUIRE_FALSE(cpr.report.fell_back);
            // GMRES-ILUT, strict, needs more iterations or does not converge at all.
            SolverConfig ilut = cpr_strict;
            ilut.preconditioner = Preconditioner::ilut;
            auto plain = LinearSolver::create(ilut).value();
            REQUIRE(plain.factorize(d.a).has_value());
            std::vector<double> x(b.size());
            const auto r = plain.solve(b, x);
            if (r) REQUIRE(cpr.report.iterations.value() < r->iterations.value());
            std::printf("cpr: %zu unknowns, bias %g: GMRES-CPR %zu, GMRES-ILUT %s%zu iterations\n",
                        b.size(), bias, cpr.report.iterations.value_or(0), r ? "" : "failed after ",
                        r ? r->iterations.value_or(0) : std::size_t{1000});
        }
    }
}

TEST_CASE("cpr: an unknown in no node (a bordering row) is left to stage 2") {
    const Index m = 16;
    const DriftDiffusion d = drift_diffusion(m, 1e3, 0.05, 5.0);
    const auto nodes = static_cast<std::size_t>(m * m);
    const auto n = static_cast<Index>(3 * nodes);
    // One more unknown: its column on a contact's potential rows, its row a dense-ish sum.
    std::vector<Triplet> t;
    for (Index r = 0; r < n; ++r) {
        for (Index k = d.a.row_offsets()[r]; k < d.a.row_offsets()[r + 1]; ++k) {
            t.push_back({r, d.a.col_indices()[k], d.a.values()[k]});
        }
    }
    for (Index j = 0; j < m; ++j) t.push_back({3 * j * m, n, -1.0});
    for (Index r = 0; r < n; r += 7) t.push_back({n, r, 0.01});
    t.push_back({n, n, 1.0});
    const SparseMatrix bordered = from(n + 1, t);
    const std::vector<double> b = rhs(static_cast<std::size_t>(n) + 1);
    const Solved s = solve_with(bordered, b, cpr_strict, labels_for(nodes, 1));
    REQUIRE_FALSE(s.report.fell_back);
}

TEST_CASE("cpr: the same labels keep the analysis, new labels redo it") {
    const DriftDiffusion d = drift_diffusion(10, 100.0, 0.1, 1.0);
    auto solver = LinearSolver::create(cpr_strict).value();
    const std::vector<std::int32_t> labels = labels_for(100);
    solver.set_field_labels(labels);
    REQUIRE(solver.factorize(d.a).has_value());
    REQUIRE(solver.analyses() == 1);
    solver.set_field_labels(labels);
    REQUIRE(solver.is_factorized());
    REQUIRE(solver.factorize(d.a).has_value());
    REQUIRE(solver.analyses() == 1);
    REQUIRE(solver.backend_counts().analyses == 1);
    std::vector<std::int32_t> other = labels;
    other[1] = -1;  // node 0's electrons left to stage 2
    other[2] = -1;
    solver.set_field_labels(other);
    REQUIRE_FALSE(solver.is_factorized());
    REQUIRE(solver.factorize(d.a).has_value());
    REQUIRE(solver.analyses() == 2);
    REQUIRE(solver.backend_counts().analyses == 2);
}

TEST_CASE("cpr: a system GMRES-CPR does not solve in time falls back, or fails in strict mode") {
    const DriftDiffusion d = drift_diffusion(16, 1e4, 0.05, 10.0);
    const std::vector<double> b = rhs(3 * 256);
    SolverConfig config = cpr_strict;
    config.max_iterations = 1;
    {
        auto solver = LinearSolver::create(config).value();
        solver.set_field_labels(labels_for(256));
        REQUIRE(solver.factorize(d.a).has_value());
        std::vector<double> x(b.size());
        const auto r = solver.solve(b, x);
        REQUIRE_FALSE(r.has_value());
        REQUIRE(r.error().code == ErrorCode::inaccurate_solve);
    }
    config.fallback = LinearFallback::sparse_lu;
    const Solved s = solve_with(d.a, b, config, labels_for(256));
    REQUIRE(s.report.fell_back);
}

TEST_CASE("cpr: a singular node block leaves that node uncoupled") {
    // Node 9's electron row reads only its potential (every other entry an explicit zero): its block
    // B_cc is singular, so its responses are 0; its neighbours keep theirs.
    DriftDiffusion d = drift_diffusion(8, 100.0, 0.1, 0.0);
    const Index row = 3 * 9 + 1;  // an interior node's electron row
    std::vector<Triplet> t;
    for (Index r = 0; r < d.a.rows(); ++r) {
        for (Index k = d.a.row_offsets()[r]; k < d.a.row_offsets()[r + 1]; ++k) {
            const Index c = d.a.col_indices()[k];
            double v = d.a.values()[k];
            if (r == row) v = c == row - 1 ? 1.0 : 0.0;  // the row reads only its potential
            t.push_back({r, c, v});
        }
    }
    const SparseMatrix a = from(d.a.rows(), t);
    NiTCAD::linalg::detail::CprBackend cpr({.strength_threshold = 0.08, .smoother_degree = 2,
                                            .coarse_size = 1000, .max_levels = 10},
                                           1e-4, 10, nullptr);
    (void)cpr.set_fields(labels_for(64));
    REQUIRE(cpr.analyze(a).has_value());
    const std::vector<double> ones(192, 1.0);
    cpr.set_scales(ones, ones, false);
    REQUIRE(cpr.factorize(a.values()).has_value());
    REQUIRE(cpr.responses()[3 * 9 + 1] == 0.0);
    REQUIRE(cpr.responses()[3 * 9 + 2] == 0.0);
    REQUIRE(cpr.responses()[3 * 10 + 1] != 0.0);
}

TEST_CASE("cpr: every thread count gives the one-thread bits (Unit 58 SG2)") {
    // 64^2 nodes, 12,288 unknowns: several row chunks (2048) and vector chunks (8192).
    const Index m = 64;
    const DriftDiffusion d = drift_diffusion(m, 1e4, 0.05, 20.0);
    const auto nodes = static_cast<std::size_t>(m * m);
    const std::vector<double> b = rhs(3 * nodes);
    const std::vector<std::int32_t> labels = labels_for(nodes);
    SolverConfig config = cpr_strict;
    config.restart = 15;  // several cycles
    const Solved one = solve_with(d.a, b, config, labels);
    for (const int threads : {2, 3, 4, 6, 8}) {
        CAPTURE(threads);
        config.threads = threads;
        const Solved s = solve_with(d.a, b, config, labels);
        REQUIRE(s.report.iterations == one.report.iterations);
        REQUIRE(std::bit_cast<std::uint64_t>(s.report.backward_error) ==
                std::bit_cast<std::uint64_t>(one.report.backward_error));
        std::size_t differ = 0;
        for (std::size_t i = 0; i < b.size(); ++i) {
            if (std::bit_cast<std::uint64_t>(s.x[i]) != std::bit_cast<std::uint64_t>(one.x[i])) ++differ;
        }
        REQUIRE(differ == 0);
    }
    std::printf("cpr threads: %zu unknowns, %zu iterations, identical at 2, 3, 4, 6, 8 threads\n",
                b.size(), one.report.iterations.value_or(0));
}

TEST_CASE("cpr: the potential system's AMG is built on the symmetric scaling of A_p") {
    // Given the equilibrated values R A C, the AMG's finest level must be G A_p G (G = sqrt(R C) on
    // the potential unknowns, slice 2's similarity), symmetric where A_p is: between interior nodes
    // at equilibrium (lambda2 times a Laplacian with permittivities 2^-10 to 2^10, plus n + p).
    // Row and column scales are chosen apart by decades, as the carrier columns make them.
    const Index m = 10;
    const DriftDiffusion d = drift_diffusion(m, 1e4, 0.05, 0.0, true);
    const auto nodes = static_cast<std::size_t>(m * m);
    const std::size_t n = 3 * nodes;
    std::vector<double> rows(n), cols(n);
    for (std::size_t u = 0; u < n; ++u) {
        rows[u] = std::exp2(std::floor(16.0 * hash01(3 * u + 1)) - 8.0);
        cols[u] = std::exp2(std::floor(16.0 * hash01(3 * u + 2)) - 8.0);
    }
    std::vector<double> scaled(d.a.nonzeros());
    for (Index r = 0; r < d.a.rows(); ++r) {
        for (Index k = d.a.row_offsets()[r]; k < d.a.row_offsets()[r + 1]; ++k) {
            const auto ku = static_cast<std::size_t>(k);
            scaled[ku] = rows[static_cast<std::size_t>(r)] * d.a.values()[ku] *
                         cols[static_cast<std::size_t>(d.a.col_indices()[k])];
        }
    }
    NiTCAD::linalg::detail::CprBackend cpr({.strength_threshold = 0.08, .smoother_degree = 2,
                                            .coarse_size = 10, .max_levels = 10},
                                           1e-4, 10, nullptr);
    (void)cpr.set_fields(labels_for(nodes));
    REQUIRE(cpr.analyze(d.a).has_value());
    cpr.set_scales(rows, cols, true);
    REQUIRE(cpr.factorize(scaled).has_value());
    const NiTCAD::linalg::detail::Csr& fine = cpr.amg().level_matrix(0);
    REQUIRE(fine.rows == nodes);
    const auto value = [&](std::size_t i, std::size_t j) {
        for (auto k = static_cast<std::size_t>(fine.offsets[i]); k < static_cast<std::size_t>(fine.offsets[i + 1]); ++k) {
            if (static_cast<std::size_t>(fine.columns[k]) == j) return fine.values[k];
        }
        return 0.0;
    };
    std::size_t pairs = 0;
    for (std::size_t i = 0; i < nodes; ++i) {
        if (d.contact[i]) continue;
        for (auto k = static_cast<std::size_t>(fine.offsets[i]); k < static_cast<std::size_t>(fine.offsets[i + 1]); ++k) {
            const auto j = static_cast<std::size_t>(fine.columns[k]);
            if (j == i || d.contact[j]) continue;
            CAPTURE(i, j);
            REQUIRE(std::abs(fine.values[k] - value(j, i)) <= 1e-12 * std::abs(fine.values[k]));
            ++pairs;
        }
    }
    REQUIRE(pairs > 100);
}
