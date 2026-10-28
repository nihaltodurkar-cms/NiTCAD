// The Newton driver on one-unknown systems, where every iterate can be checked by hand:
// convergence, damping, the full-correction criterion and each error path.
#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"
#include "NiTCAD/results/convergence.hpp"
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
    results::ConvergenceRecord record;
    REQUIRE(solve::newton_solve(cube, x, {}, s, record).has_value());
    REQUIRE(std::abs(x[0] - 2.0) <= 1e-15);
    REQUIRE(record.converged);
    const auto& it = record.iterations;
    REQUIRE(it.back().update < 1e-8);
    // Quadratic: each update is at most a modest multiple of the previous one squared.
    for (std::size_t k = 1; k < it.size(); ++k) {
        CAPTURE(k);
        REQUIRE(it[k].iteration == static_cast<int>(k) + 1);
        REQUIRE(it[k].update <= 2.0 * it[k - 1].update * it[k - 1].update);
    }
    REQUIRE(it.front().residual == 19.0);  // F(3) = 27 - 8
    REQUIRE(s.analyses() == 1);
    REQUIRE(s.factorizations() == it.size());
}

TEST_CASE("newton: the step is clipped, the criterion uses the full correction") {
    // From x = 0.1: F = -7.999, J = 0.03, so the full correction is 266.6; it is clipped to 5.
    std::vector<double> x{0.1};
    auto s = solver();
    results::ConvergenceRecord record;
    const auto one = solve::newton_solve(cube, x, NewtonOptions{.max_iterations = 1}, s, record);
    REQUIRE_FALSE(one.has_value());
    REQUIRE(one.error().code == ErrorCode::non_convergence);
    REQUIRE(one.error().context->index == 1);
    REQUIRE(std::abs(*one.error().context->value - 7.999 / 0.03) <= 1e-9);
    REQUIRE(x[0] == 0.1 + 5.0);  // the last iterate is kept
    REQUIRE(record.iterations.size() == 1);  // and so is the history
    REQUIRE_FALSE(record.converged);
    // A tolerance above the clipped step but below the full correction does not stop it.
    std::vector<double> y{0.1};
    const auto wide = solve::newton_solve(
        cube, y, NewtonOptions{.max_iterations = 1, .tol_update = 6.0, .max_update = 5.0}, s,
        record);
    REQUIRE_FALSE(wide.has_value());
    // Damped run still converges.
    std::vector<double> z{0.1};
    const auto full = solve::newton_solve(cube, z, NewtonOptions{.max_update = 0.5}, s, record);
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
        results::ConvergenceRecord record;
        REQUIRE(solve::newton_solve(cube, x, o, s, record).error().code ==
                ErrorCode::invalid_input);
    }
    REQUIRE(x[0] == 3.0);
}

TEST_CASE("newton: a non-finite residual is non-convergence, not a linear-solver error") {
    const Scalar overflow{[](double x) { return std::exp(x) - 1.0; },
                          [](double x) { return std::exp(x); }};
    std::vector<double> x{800.0};
    auto s = solver();
    results::ConvergenceRecord record;
    const auto r = solve::newton_solve(overflow, x, {}, s, record);
    REQUIRE_FALSE(r.has_value());
    REQUIRE(r.error().code == ErrorCode::non_convergence);
    REQUIRE(r.error().message.find("not finite") != std::string::npos);
}

TEST_CASE("newton: a singular Jacobian is reported with the iteration") {
    const Scalar flat{[](double x) { return x * x + 1.0; }, [](double x) { return 2.0 * x; }};
    std::vector<double> x{0.0};
    auto s = solver();
    results::ConvergenceRecord record;
    const auto r = solve::newton_solve(flat, x, {}, s, record);
    REQUIRE_FALSE(r.has_value());
    REQUIRE(r.error().code == ErrorCode::singular_system);
    REQUIRE(r.error().message.starts_with("Newton iteration 1: "));
}

// --- Line search (Unit 30, solve::LineSearch) ----------------------------------------------------

namespace {

// arctan: Newton from |x| > 1.3917 overshoots further each iteration and diverges.
const Scalar arctan{[](double x) { return std::atan(x); },
                    [](double x) { return 1.0 / (1.0 + x * x); }};

NewtonOptions searched(int iterations = 100,
                       solve::LineSearchTest test = solve::LineSearchTest::natural) {
    NewtonOptions o{.max_iterations = iterations, .tol_update = 1e-12, .max_update = 1e6};
    o.line_search.enabled = true;
    o.line_search.test = test;
    return o;
}

constexpr std::array tests{solve::LineSearchTest::residual, solve::LineSearchTest::natural};

}  // namespace

TEST_CASE("line search: Newton on arctan from x = 2 diverges without it and converges with it") {
    auto s = solver();
    results::ConvergenceRecord record;
    std::vector<double> x{2.0};
    const auto plain = solve::newton_solve(
        arctan, x, NewtonOptions{.max_iterations = 8, .tol_update = 1e-12, .max_update = 1e6}, s,
        record);
    REQUIRE_FALSE(plain.has_value());
    CHECK(std::abs(x[0]) > 1e4);  // the iterates run away (1.2e5 after 8 iterations)

    for (const auto test : tests) {
        CAPTURE(static_cast<int>(test));
        std::vector<double> y{2.0};
        REQUIRE(solve::newton_solve(arctan, y, searched(100, test), s, record).has_value());
        CHECK(std::abs(y[0]) < 1e-12);
        // The first step is shortened; the last ones are full, so convergence stays quadratic.
        CHECK(record.iterations.front().step < 1.0);
        CHECK(record.iterations.back().step == 1.0);
        for (std::size_t k = 1; k < record.iterations.size(); ++k) {
            const auto& before = record.iterations[k - 1];
            const auto& after = record.iterations[k];
            if (test == solve::LineSearchTest::residual) {
                // Each step lowered the merit |F|^2 by the sufficient-decrease factor.
                CHECK(after.residual * after.residual <=
                      (1.0 - 2.0 * 1e-4 * before.step) * before.residual * before.residual);
            } else {
                // Each step shrank the next correction (the natural monotonicity test).
                CHECK(after.update <= (1.0 - 0.25 * before.step) * before.update);
            }
        }
    }
}

TEST_CASE("line search: where every full step is accepted the run is bit-identical to none") {
    for (const double start : {3.0, 1.5, 10.0}) {
        auto s = solver();
        results::ConvergenceRecord off_record, on_record;
        std::vector<double> off{start}, on{start};
        NewtonOptions plain{.max_iterations = 100, .tol_update = 1e-12, .max_update = 1e6};
        REQUIRE(solve::newton_solve(cube, off, plain, s, off_record).has_value());
        REQUIRE(solve::newton_solve(cube, on, searched(), s, on_record).has_value());
        CAPTURE(start);
        bool full = true;
        for (const auto& r : on_record.iterations) full = full && r.step == 1.0;
        if (!full) continue;  // a start that backtracks is not this test's case
        CHECK(on[0] == off[0]);
        REQUIRE(on_record.iterations.size() == off_record.iterations.size());
        for (std::size_t k = 0; k < on_record.iterations.size(); ++k) {
            CHECK(on_record.iterations[k].update == off_record.iterations[k].update);
            CHECK(on_record.iterations[k].residual == off_record.iterations[k].residual);
        }
    }
}

TEST_CASE("line search: with no descent it takes the smallest step, never zero, and reports "
          "non-convergence") {
    // A wrong-signed Jacobian: the correction points uphill, so no trial decreases the merit.
    const Scalar uphill{[](double x) { return x; }, [](double) { return -1.0; }};
    auto s = solver();
    results::ConvergenceRecord record;
    std::vector<double> x{1.0};
    auto o = searched(5);
    o.line_search.max_backtracks = 4;
    const auto r = solve::newton_solve(uphill, x, o, s, record);
    REQUIRE_FALSE(r.has_value());
    CHECK(r.error().code == ErrorCode::non_convergence);
    REQUIRE(record.iterations.size() == 5);
    for (const auto& it : record.iterations) CHECK(it.step == 1.0 / 16.0);
    CHECK(x[0] > 1.0);  // it moved (uphill), it did not stall at lambda = 0
}

TEST_CASE("line search: a trial with a non-finite residual counts as no decrease") {
    // F = x - 1, NaN beyond x = 4, with a Jacobian of 0.25 (a quarter of the slope): from -0.5 the
    // correction is +6, to 5.5, where F is NaN. Without a search the run stops there; with one
    // that trial is rejected, 2.5 (merit 2.25, no decrease) too, and 1 (the root) accepted.
    const Scalar overshoot{[](double x) { return x > 4.0 ? std::nan("") : x - 1.0; },
                           [](double) { return 0.25; }};
    auto s = solver();
    results::ConvergenceRecord record;
    std::vector<double> plain_x{-0.5};
    const auto plain = solve::newton_solve(
        overshoot, plain_x,
        NewtonOptions{.max_iterations = 10, .tol_update = 1e-12, .max_update = 1e6}, s, record);
    REQUIRE_FALSE(plain.has_value());
    CHECK(plain.error().message.find("not finite") != std::string::npos);

    std::vector<double> x{-0.5};
    REQUIRE(solve::newton_solve(overshoot, x, searched(), s, record).has_value());
    CHECK(x[0] == 1.0);
    CHECK(record.iterations.front().step == 0.25);
}

TEST_CASE("line search: invalid options are errors") {
    auto s = solver();
    results::ConvergenceRecord record;
    std::vector<double> x{3.0};
    for (const auto& bad : {solve::LineSearch{true, -1, 1e-4}, solve::LineSearch{true, 61, 1e-4},
                            solve::LineSearch{true, 10, 0.0}, solve::LineSearch{true, 10, 0.5},
                            solve::LineSearch{true, 10, std::nan("")},
                            solve::LineSearch{true, 10, 1e-4,
                                              static_cast<solve::LineSearchTest>(7)}}) {
        NewtonOptions o;
        o.line_search = bad;
        const auto r = solve::newton_solve(cube, x, o, s, record);
        REQUIRE_FALSE(r.has_value());
        CHECK(r.error().code == ErrorCode::invalid_input);
    }
}

namespace {

// Two decoupled rows of very different scales (Unit 50, RD1): F_0 = 1e14 (x_0 - 1e7) + noise(x_1),
// a row of scale |J x| = 1e21 holding a perturbation the Jacobian leaves out, and x_0 = 1e7, where
// its correction (about 1e-10, below half an ulp of 1e7) cannot move it: the row stays at its
// rounding, as a device's large rows do. F_1 = c (d + curvature d^2), d = x_1 - target, a row 1e21
// smaller or more that needs a real correction (a minority density). J = diag(1e14, F_1').
struct TwoScales {
    std::function<double(double)> noise;
    double c = 1.0;
    double target = 1e-3;
    double curvature = 0.0;  // F_1 = c (d + curvature d^2), d = x_1 - target

    std::size_t unknowns() const { return 2; }
    linalg::SparseMatrix make_jacobian() const {
        const std::array<linalg::Triplet, 2> t{linalg::Triplet{0, 0, 0.0}, linalg::Triplet{1, 1, 0.0}};
        return *linalg::SparseMatrix::from_triplets(2, 2, t);
    }
    void evaluate(std::span<const double> x, std::span<double> r, linalg::SparseMatrix& j) const {
        r[0] = 1e14 * (x[0] - 1e7) + noise(x[1]);
        const double d = x[1] - target;
        r[1] = c * (d + curvature * d * d);
        j.values()[0] = 1e14;
        j.values()[1] = c * (1.0 + 2.0 * curvature * d);
    }
};

}  // namespace

TEST_CASE("line search: the residual test sees a small row past the large rows' rounding (Unit 50)") {
    // From x = (1e7, 0): the large row sits at 1e4 and grows by 2e-4 of itself as x_1 rises by
    // 1e-3, so the first full step raises sum F^2 by about 5e4 while the small row's share, 8.1e-7,
    // falls: the plain Armijo test rejects every trial and the iteration would creep by 2^-10 of
    // its step. The small row is curved (curvature 100), so the step leaves it at 1.3e-4, far
    // above its rounding: only the row-scaled merit, (1e4 / 1e21)^2 + (9e-4 / 1)^2 (|J x| of the
    // small row is 0, taken as 1) falling to about 1.6e-8, can pass it.
    TwoScales system{[](double y) { return 1e4 * (1.0 + 0.2 * y); }};
    system.curvature = 100.0;
    auto s = solver();
    results::ConvergenceRecord record;
    std::vector<double> x{1e7, 0.0};
    auto o = searched(100, solve::LineSearchTest::residual);
    o.tol_update = 1e-6;
    REQUIRE(solve::newton_solve(system, x, o, s, record).has_value());
    CHECK(x[0] == 1e7);  // pinned at its rounding
    CHECK(std::abs(x[1] - 1e-3) <= 1e-9);
    CHECK(record.iterations.front().step == 1.0);
    CHECK(record.iterations.size() <= 5);
}

TEST_CASE("line search: a trial with every row at its rounding passes the residual test (Unit 50)") {
    // The small row is c = 1e-30 times (x_1 - 1): at x_1 = 0 its residual, 1e-30, is within 64 eps
    // of its scale (|J x| = 0, taken as 1), yet its correction is 1. The full step zeroes it while
    // the large row's perturbation grows from 1e4 to 1.01e6, still within 64 eps of its scale
    // 1e21 (2.8e7): neither merit falls, and the trial passes as rounding in every row.
    const TwoScales system{[](double y) { return 1e4 * (1.0 + 100.0 * y); }, 1e-30, 1.0};
    auto s = solver();
    results::ConvergenceRecord record;
    std::vector<double> x{1e7, 0.0};
    auto o = searched(100, solve::LineSearchTest::residual);
    o.tol_update = 1e-6;
    REQUIRE(solve::newton_solve(system, x, o, s, record).has_value());
    CHECK(x[1] == 1.0);
    CHECK(record.iterations.front().step == 1.0);
    CHECK(record.iterations.size() == 2);
}
