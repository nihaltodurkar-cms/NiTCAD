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
