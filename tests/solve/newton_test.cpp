// The Newton driver on one-unknown systems, where every iterate can be checked by hand:
// convergence, damping, the full-correction criterion and each error path.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstddef>
#include <functional>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"
#include "NiTCAD/solve/newton.hpp"

using namespace NiTCAD;
using base::ErrorCode;
using solve::NewtonOptions;

namespace {

// F(x) = f(x), J = f'(x), one unknown.
struct Scalar {
    std::function<double(double)> f;
    std::function<double(double)> df;

    std::size_t unknowns() const { return 1; }
    linalg::SparseMatrix make_jacobian() const {
        const linalg::Triplet t{0, 0, 0.0};
        return *linalg::SparseMatrix::from_triplets(1, 1, std::span(&t, 1));
    }
    void evaluate(std::span<const double> x, std::span<double> r, linalg::SparseMatrix& j) const {
        r[0] = f(x[0]);
        j.values()[0] = df(x[0]);
    }
};

const Scalar cube{[](double x) { return x * x * x - 8.0; }, [](double x) { return 3.0 * x * x; }};

linalg::LinearSolver solver() { return *linalg::LinearSolver::create({}); }

}  // namespace

TEST_CASE("newton: converges quadratically to the root") {
    std::vector<double> x{3.0};
    auto s = solver();
    const auto report = solve::newton_solve(cube, x, {}, s);
    REQUIRE(report.has_value());
    REQUIRE(std::abs(x[0] - 2.0) <= 1e-15);
    REQUIRE(report->updates.back() < 1e-8);
    REQUIRE(report->iterations == static_cast<int>(report->updates.size()));
    // Quadratic: each update is at most a modest multiple of the previous one squared.
    for (std::size_t k = 1; k < report->updates.size(); ++k) {
        CAPTURE(k);
        REQUIRE(report->updates[k] <= 2.0 * report->updates[k - 1] * report->updates[k - 1]);
    }
    REQUIRE(s.analyses() == 1);
    REQUIRE(s.factorizations() == static_cast<std::size_t>(report->iterations));
}

TEST_CASE("newton: the step is clipped, the criterion uses the full correction") {
    // From x = 0.1: F = -7.999, J = 0.03, so the full correction is 266.6; it is clipped to 5.
    std::vector<double> x{0.1};
    auto s = solver();
    const auto one = solve::newton_solve(cube, x, NewtonOptions{.max_iterations = 1}, s);
    REQUIRE_FALSE(one.has_value());
    REQUIRE(one.error().code == ErrorCode::non_convergence);
    REQUIRE(one.error().context->index == 1);
    REQUIRE(std::abs(*one.error().context->value - 7.999 / 0.03) <= 1e-9);
    REQUIRE(x[0] == 0.1 + 5.0);  // the last iterate is kept
    // A tolerance above the clipped step but below the full correction does not stop it.
    std::vector<double> y{0.1};
    const auto wide = solve::newton_solve(
        cube, y, NewtonOptions{.max_iterations = 1, .tol_update = 6.0, .max_update = 5.0}, s);
    REQUIRE_FALSE(wide.has_value());
    // Damped run still converges.
    std::vector<double> z{0.1};
    const auto full = solve::newton_solve(cube, z, NewtonOptions{.max_update = 0.5}, s);
    REQUIRE(full.has_value());
    REQUIRE(std::abs(z[0] - 2.0) <= 1e-15);
}

TEST_CASE("newton: options are validated") {
    std::vector<double> x{3.0};
    auto s = solver();
    for (const NewtonOptions& o : {NewtonOptions{.max_iterations = 0},
                                   NewtonOptions{.tol_update = 0.0},
                                   NewtonOptions{.tol_update = std::nan("")},
                                   NewtonOptions{.max_update = -1.0},
                                   NewtonOptions{.max_update = INFINITY}}) {
        REQUIRE(solve::newton_solve(cube, x, o, s).error().code == ErrorCode::invalid_input);
    }
    REQUIRE(x[0] == 3.0);
}

TEST_CASE("newton: a non-finite residual is non-convergence, not a linear-solver error") {
    const Scalar overflow{[](double x) { return std::exp(x) - 1.0; },
                          [](double x) { return std::exp(x); }};
    std::vector<double> x{800.0};
    auto s = solver();
    const auto r = solve::newton_solve(overflow, x, {}, s);
    REQUIRE_FALSE(r.has_value());
    REQUIRE(r.error().code == ErrorCode::non_convergence);
    REQUIRE(r.error().message.find("not finite") != std::string::npos);
}

TEST_CASE("newton: a singular Jacobian is reported with the iteration") {
    const Scalar flat{[](double x) { return x * x + 1.0; }, [](double x) { return 2.0 * x; }};
    std::vector<double> x{0.0};
    auto s = solver();
    const auto r = solve::newton_solve(flat, x, {}, s);
    REQUIRE_FALSE(r.has_value());
    REQUIRE(r.error().code == ErrorCode::singular_system);
    REQUIRE(r.error().message.starts_with("Newton iteration 1: "));
}
