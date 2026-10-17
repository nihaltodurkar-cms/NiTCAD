// Assembly cost (Unit 40, DECISIONS.md FG4): the median time of DriftDiffusion::evaluate (residual
// and Jacobian) and of residual alone at a state, for the hidden [.assembly_benchmark] cases.
// Each repetition repeats the call until it lasts at least 20 ms (calibrated first, which also
// warms the caches), so the timer's resolution does not matter; the median of the repetitions is
// reported per call.
#pragma once

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdio>
#include <span>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/results/solution.hpp"

struct AssemblyTiming {
    double evaluate_us;  // median per evaluate call
    double residual_us;  // median per residual call
};

// The scaled state of a solved point's fields (the inverse of the solve layer's
// drift_diffusion_point): psi / V_T, n / Ns, p / Ns and, with the electrothermal model,
// tau = (T - T0) / T0.
inline std::vector<double> scaled_state(const NiTCAD::assemble::DriftDiffusion& system,
                                        const NiTCAD::results::NodeFields& f,
                                        const NiTCAD::assemble::Scaling& scaling) {
    const std::size_t m = system.stride();
    std::vector<double> x(system.unknowns());
    for (std::size_t i = 0; i < system.node_count(); ++i) {
        x[m * i] = f.potential_V[i] / scaling.V_T;
        x[m * i + 1] = f.n_cm3[i] / scaling.Ns;
        x[m * i + 2] = f.p_cm3[i] / scaling.Ns;
        if (system.electrothermal()) {
            const double T0 = scaling.temperature_K;
            x[m * i + 3] = (f.temperature_K[i] - T0) / T0;
        }
    }
    return x;
}

template <class Call>
double median_call_us(Call&& call, int repetitions) {
    using clock = std::chrono::steady_clock;
    const auto seconds = [](clock::duration d) { return std::chrono::duration<double>(d).count(); };
    long calls = 1;
    for (;;) {  // calibration
        const auto t0 = clock::now();
        for (long k = 0; k < calls; ++k) call();
        if (seconds(clock::now() - t0) >= 0.02) break;
        calls *= 2;
    }
    std::vector<double> t;
    for (int r = 0; r < repetitions; ++r) {
        const auto t0 = clock::now();
        for (long k = 0; k < calls; ++k) call();
        t.push_back(1e6 * seconds(clock::now() - t0) / static_cast<double>(calls));
    }
    std::ranges::nth_element(t, t.begin() + t.size() / 2);
    return t[t.size() / 2];
}

inline AssemblyTiming time_assembly(const char* name,
                                    const NiTCAD::assemble::DriftDiffusion& system,
                                    std::span<const double> x, int repetitions = 31) {
    std::vector<double> f(system.unknowns());
    NiTCAD::linalg::SparseMatrix j = system.make_jacobian();
    AssemblyTiming t{};
    t.evaluate_us = median_call_us([&] { system.evaluate(x, f, j); }, repetitions);
    t.residual_us = median_call_us([&] { system.residual(x, f); }, repetitions);
    std::printf("assembly %s: %zu unknowns, %zu nonzeros: evaluate %.3f us, residual %.3f us "
                "(median of %d)\n",
                name, system.unknowns(), j.nonzeros(), t.evaluate_us, t.residual_us, repetitions);
    return t;
}
