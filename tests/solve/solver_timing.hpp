// Linear-solver cost on a device's drift-diffusion Jacobian (Unit 18 measurements): the first
// factorization (with the analysis), the best of the numeric refactorizations of the same pattern,
// and the best solve, in seconds, for one solver configuration.
#pragma once

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <span>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"

struct SolverTiming {
    double first_s;     // analysis and first factorization
    double refactor_s;  // best numeric refactorization
    double solve_s;     // best solve
    std::vector<double> x;
};

// The Jacobian at the state of the given potential (V) under the given bias.
inline NiTCAD::linalg::SparseMatrix device_jacobian(const NiTCAD::device::Device& d,
                                                    std::span<const double> bias_V,
                                                    std::span<const double> potential_V,
                                                    std::vector<double>& rhs) {
    using namespace NiTCAD;
    const auto scaling = *assemble::make_scaling(d);
    auto system = *assemble::DriftDiffusion::create(d, scaling);
    REQUIRE(system.set_bias(bias_V).has_value());
    std::vector<double> psi(potential_V.size());
    for (std::size_t i = 0; i < psi.size(); ++i) psi[i] = potential_V[i] / scaling.V_T;
    const auto x = system.state_from_potential(psi);
    linalg::SparseMatrix j = system.make_jacobian();
    rhs.assign(x.size(), 0.0);
    system.evaluate(x, rhs, j);
    return j;
}

inline SolverTiming time_solver(const NiTCAD::linalg::SparseMatrix& a, std::span<const double> b,
                                const NiTCAD::linalg::SolverConfig& config, int repeats = 3) {
    using clock = std::chrono::steady_clock;
    const auto seconds = [](clock::time_point t0) {
        return std::chrono::duration<double>(clock::now() - t0).count();
    };
    auto solver = *NiTCAD::linalg::LinearSolver::create(config);
    SolverTiming t{0.0, 1e300, 1e300, std::vector<double>(b.size())};
    auto t0 = clock::now();
    REQUIRE(solver.factorize(a).has_value());
    t.first_s = seconds(t0);
    for (int k = 0; k < repeats; ++k) {
        t0 = clock::now();
        REQUIRE(solver.factorize(a).has_value());
        t.refactor_s = std::min(t.refactor_s, seconds(t0));
        t0 = clock::now();
        REQUIRE(solver.solve(b, t.x).has_value());
        t.solve_s = std::min(t.solve_s, seconds(t0));
    }
    REQUIRE(solver.backend_counts().analyses == 1);
    REQUIRE(solver.backend_counts().factorizations == static_cast<std::size_t>(repeats) + 1);
    return t;
}

// Unit 44 (LG4): the cost of one Newton step's linear solve with a fresh solver (analysis,
// factorization or preconditioner build, and one solve), its GMRES iterations and whether it fell
// back, and memory: the private bytes the solver holds after the solve, and how far the process's
// peak working set rose during it (a lower bound on the step's own peak; exact when the step sets
// a new peak, so run the cheaper configurations first).
struct StepCost {
    double seconds;
    std::size_t iterations;  // 0 for a direct backend
    bool fell_back;
    double held_MB;
    double peak_rise_MB;
    double backward_error;
    std::vector<double> x;
};

struct ProcessMemory {
    double private_MB;
    double peak_working_set_MB;
};

ProcessMemory process_memory();  // process_memory.cpp: keeps windows.h out of the tests

inline StepCost newton_step_cost(const NiTCAD::linalg::SparseMatrix& a, std::span<const double> b,
                                 const NiTCAD::linalg::SolverConfig& config) {
    using clock = std::chrono::steady_clock;
    const ProcessMemory before = process_memory();
    StepCost c{0.0, 0, false, 0.0, 0.0, 0.0, std::vector<double>(b.size())};
    {
        const auto t0 = clock::now();
        auto solver = *NiTCAD::linalg::LinearSolver::create(config);
        REQUIRE(solver.factorize(a).has_value());
        const auto r = solver.solve(b, c.x);
        c.seconds = std::chrono::duration<double>(clock::now() - t0).count();
        REQUIRE(r.has_value());
        c.iterations = r->iterations.value_or(0);
        c.fell_back = r->fell_back;
        c.backward_error = r->backward_error;
        const ProcessMemory after = process_memory();
        c.held_MB = after.private_MB - before.private_MB;
        c.peak_rise_MB = after.peak_working_set_MB - before.peak_working_set_MB;
    }
    return c;
}
