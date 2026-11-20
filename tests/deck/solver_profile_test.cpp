// Unit 58 slice 0 (SD1): where a Newton iteration's time goes, before any solver change. Hidden,
// Release, run by hand pinned to the P cores (DECISIONS.md Unit 58, slice 0):
//   start /b /wait /affinity 3C3 nitcad_deck_test.exe "[.scale_profile]"
//   start /b /wait /affinity 3C3 nitcad_deck_test.exe "[.scale_series]"
//   start /b /wait /affinity 3C3 nitcad_deck_test.exe "[.scale_threads]"   (slice 1)
//   start /b /wait /affinity 3C3 nitcad_deck_test.exe "[.scale_amg]"       (slice 2)
//   start /b /wait /affinity 3C3 nitcad_deck_test.exe "[.scale_cpr]"       (slice 3)
//   start /b /wait /affinity 3C3 nitcad_deck_test.exe "[.scale_carrier_blocks]"  (Unit 63 slice 0)
// Measurement only: nothing here is gated beyond the cases solving, and no product code changes.
//
// Two measurements per case and linear-solver configuration, from the same start (the case's bias
// stamped on its equilibrium state, as sweep_bias's first point):
// - The phase split: a plain Newton loop written here (full damped steps, no line search) that
//   times each phase: assembly (System::evaluate), the right-hand side and residual scan (rest),
//   the symbolic analysis (once), the first factorization, later factorizations, the solve, and the
//   update (update_size and apply_update). Its iterations are the production loop's without
//   backtracking.
// - The production loop: solve::newton_solve with the case's own Newton options (line search
//   included) on a forwarding system that times evaluate and apply_update; the remainder is the
//   linear solver's factorizations and solves plus the loop's own O(n) work.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <bit>
#include <cstdint>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <filesystem>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "../linalg/block_probe.hpp"
#include "NiTCAD/assemble/energy.hpp"
#include "NiTCAD/assemble/equilibrium.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/assemble/system.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/equilibrium.hpp"
#include "NiTCAD/solve/newton.hpp"
#include "../../src/solve/fields.hpp"
#include "../linalg/mkl_test_runtime.hpp"
#include "../solve/solver_timing.hpp"

using namespace NiTCAD;
namespace fs = std::filesystem;

namespace {

using clock_type = std::chrono::steady_clock;
double since(clock_type::time_point t0) {
    return std::chrono::duration<double>(clock_type::now() - t0).count();
}

struct Configuration {
    const char* name;
    linalg::SolverConfig config;
};

// The configurations of SD1, cheapest memory first (the peak working set only rises).
std::vector<Configuration> configurations(bool sparse_lu, bool mkl) {
    std::vector<Configuration> c;
    c.push_back({"GMRES-ILUT", {.backend = linalg::SolverBackend::gmres_ilut,
                                .fallback = linalg::LinearFallback::none}});
    if (mkl) {
        c.push_back({"PARDISO 1 thread", {.backend = linalg::SolverBackend::mkl_pardiso}});
        c.push_back({"PARDISO 6 threads",
                     {.backend = linalg::SolverBackend::mkl_pardiso, .threads = 6}});
    }
    if (sparse_lu) c.push_back({"SparseLU", {}});
    return c;
}

struct Split {
    double assembly = 0.0, rest = 0.0, analysis = 0.0, first_factorization = 0.0,
           factorization = 0.0, solve = 0.0, update = 0.0;
    int iterations = 0;
    bool converged = false;
    std::size_t gmres_iterations = 0;
    std::optional<std::size_t> fill;
    double held_MB = 0.0, peak_rise_MB = 0.0, total = 0.0;
    std::string error;
};

// The plain Newton loop with each phase timed (see the header).
Split split_newton(const assemble::System& system, std::vector<double> x,
                   const solve::NewtonOptions& newton, const linalg::SolverConfig& config,
                   int max_iterations) {
    Split s;
    const ProcessMemory before = process_memory();
    const auto start = clock_type::now();
    auto made = linalg::LinearSolver::create(config);
    if (!made) {
        s.error = made.error().message;
        return s;
    }
    linalg::LinearSolver& solver = *made;
    const std::size_t n = system.unknowns();
    linalg::SparseMatrix j = system.make_jacobian();
    std::vector<double> f(n), rhs(n), dx(n);
    for (int it = 1; it <= max_iterations; ++it) {
        auto t0 = clock_type::now();
        system.evaluate(x, f, j);
        s.assembly += since(t0);
        t0 = clock_type::now();
        for (std::size_t i = 0; i < n; ++i) rhs[i] = -f[i];
        s.rest += since(t0);
        if (it == 1) {
            t0 = clock_type::now();
            if (auto a = solver.analyze(j); !a) {
                s.error = "analyze: " + a.error().message;
                break;
            }
            s.analysis = since(t0);
        }
        t0 = clock_type::now();
        auto factored = solver.factorize(j);
        (it == 1 ? s.first_factorization : s.factorization) += since(t0);
        if (!factored) {
            s.error = "iteration " + std::to_string(it) + ": " + factored.error().message;
            break;
        }
        s.fill = factored->fill;
        if (it == 1) {
            const ProcessMemory after = process_memory();
            s.held_MB = after.private_MB - before.private_MB;
        }
        t0 = clock_type::now();
        auto solved = solver.solve(rhs, dx);
        s.solve += since(t0);
        if (!solved) {
            s.error = "iteration " + std::to_string(it) + ": " + solved.error().message;
            break;
        }
        s.gmres_iterations += solved->iterations.value_or(0);
        t0 = clock_type::now();
        const double size = system.update_size(x, dx);
        system.apply_update(x, dx, newton.max_update);
        s.update += since(t0);
        s.iterations = it;
        if (size < newton.tol_update) {
            s.converged = true;
            break;
        }
    }
    s.total = since(start);
    s.peak_rise_MB = process_memory().peak_working_set_MB - before.peak_working_set_MB;
    return s;
}

// A system that forwards to the real one and times its assembly and update calls.
struct TimedSystem {
    const assemble::System& system;
    mutable double assembly = 0.0, update = 0.0;
    mutable std::size_t evaluations = 0;
    std::size_t unknowns() const { return system.unknowns(); }
    linalg::SparseMatrix make_jacobian() const { return system.make_jacobian(); }
    void evaluate(std::span<const double> x, std::span<double> f, linalg::SparseMatrix& j) const {
        const auto t0 = clock_type::now();
        system.evaluate(x, f, j);
        assembly += since(t0);
        ++evaluations;
    }
    double update_size(std::span<const double> x, std::span<const double> dx) const {
        const auto t0 = clock_type::now();
        const double s = system.update_size(x, dx);
        update += since(t0);
        return s;
    }
    void apply_update(std::span<double> x, std::span<const double> dx, double max_update) const {
        const auto t0 = clock_type::now();
        system.apply_update(x, dx, max_update);
        update += since(t0);
    }
    double natural_size(std::span<const double> x_ref, std::span<const double> dx) const {
        return system.natural_size(x_ref, dx);
    }
};

struct Production {
    double total = 0.0, assembly = 0.0, update = 0.0;
    std::size_t iterations = 0, evaluations = 0;
    std::string error;
};

Production production_newton(const assemble::System& system, std::vector<double> x,
                             const solve::NewtonOptions& newton,
                             const linalg::SolverConfig& config) {
    Production p;
    auto solver = linalg::LinearSolver::create(config);
    if (!solver) {
        p.error = solver.error().message;
        return p;
    }
    const TimedSystem timed{system};
    results::ConvergenceRecord record;
    const auto t0 = clock_type::now();
    const auto r = solve::newton_solve(timed, x, newton, *solver, record);
    p.total = since(t0);
    if (!r) p.error = r.error().message;
    p.iterations = record.iterations.size();
    p.evaluations = timed.evaluations;
    p.assembly = timed.assembly;
    p.update = timed.update;
    return p;
}

struct Case {
    std::string name;
    device::Device device;
    solve::BiasOptions options;
    std::vector<double> bias;
};

bool mkl_loaded() { return !mkl_test::configured_path.empty() && mkl_test::runtime().has_value(); }

// The fastest available solver for the untimed parts (equilibrium, energy balance's first stage).
linalg::SolverConfig helper_solver() {
    if (mkl_loaded()) return {.backend = linalg::SolverBackend::mkl_pardiso, .threads = 6};
    return {.backend = linalg::SolverBackend::gmres_ilut};
}

void profile(const Case& c, bool sparse_lu, int sparse_lu_iterations = 40) {
    const bool mkl = mkl_loaded();
    const solve::BiasOptions& o = c.options;
    const auto scaling = assemble::make_scaling(c.device, o.Ns_override);
    REQUIRE(scaling.has_value());
    auto t0 = clock_type::now();
    const auto eq = solve::solve_equilibrium(
        c.device, {.newton = o.newton, .linear = helper_solver(), .Ns_override = o.Ns_override,
                   .models = o.models});
    if (!eq) FAIL(c.name << ": equilibrium: " << eq.error().message);
    const double eq_s = since(t0);
    auto made = assemble::System::create(c.device, *scaling, o.models);
    if (!made) FAIL(c.name << ": " << made.error().message);
    assemble::System& system = *made;
    assemble::set_minority_floor(system, o.minority_floor);
    std::vector<double> x = solve::detail::state_from_equilibrium(system, eq->fields, scaling->V_T);
    REQUIRE(system.set_bias(c.bias).has_value());
    system.stamp_contacts(x);
    if (assemble::energy_balance(system)) {
        // ED12's first stage (untimed): psi, n, p at the carrier temperatures the point starts from.
        assemble::hold_temperatures(system, x);
        auto solver = *linalg::LinearSolver::create(helper_solver());
        results::ConvergenceRecord held;
        const auto first = solve::newton_solve(system, x, o.newton, solver, held);
        assemble::release_temperatures(system);
        if (!first) FAIL(c.name << ": energy balance first stage: " << first.error().message);
    }
    linalg::SparseMatrix j = system.make_jacobian();
    std::printf("\n=== %s: %zu nodes, %zu unknowns (stride %zu), %zu Jacobian nonzeros; equilibrium "
                "%.2f s\n",
                c.name.c_str(), system.node_count(), system.unknowns(), system.stride(),
                j.nonzeros(), eq_s);
    std::fflush(stdout);
    for (const Configuration& k : configurations(sparse_lu, mkl)) {
        const bool lu = k.config.backend == linalg::SolverBackend::eigen_sparse_lu;
        const Split s = split_newton(system, x, o.newton, k.config, lu ? sparse_lu_iterations : 40);
        const int its = std::max(s.iterations, 1);
        const double later = s.iterations > 1 ? s.factorization / (s.iterations - 1) : 0.0;
        const double sum = s.assembly + s.rest + s.analysis + s.first_factorization +
                           s.factorization + s.solve + s.update;
        const auto pct = [&](double v) { return sum > 0.0 ? 100.0 * v / sum : 0.0; };
        std::printf("  %-18s split: %d iterations (%s)%s%s; per iteration assembly %.4f s, "
                    "factorization %.4f s (first %.4f s, analysis %.4f s), solve %.4f s, update "
                    "%.5f s, rest %.5f s; shares of %.2f s: assembly %.1f%%, analysis %.1f%%, "
                    "factorizations %.1f%%, solves %.1f%%, update %.2f%%, rest %.2f%%",
                    k.name, s.iterations, s.converged ? "converged" : "not converged",
                    s.error.empty() ? "" : ", stopped: ", s.error.c_str(), s.assembly / its,
                    later, s.first_factorization, s.analysis, s.solve / its, s.update / its,
                    s.rest / its, sum, pct(s.assembly), pct(s.analysis),
                    pct(s.first_factorization + s.factorization), pct(s.solve), pct(s.update),
                    pct(s.rest));
        if (s.gmres_iterations > 0) {
            std::printf("; GMRES %zu products (%.1f per solve)", s.gmres_iterations,
                        static_cast<double>(s.gmres_iterations) / its);
        }
        std::printf("; fill %zu; held %.0f MB, peak rise %.0f MB\n", s.fill.value_or(0), s.held_MB,
                    s.peak_rise_MB);
        std::fflush(stdout);
        if (lu) continue;  // the production loop only for the threaded and iterative backends
        const Production p = production_newton(system, x, o.newton, k.config);
        const double linear = p.total - p.assembly - p.update;
        std::printf("  %-18s newton_solve: %zu iterations, %zu evaluations%s%s; %.2f s: assembly "
                    "%.1f%%, update %.2f%%, linear and the loop's own work %.1f%%\n",
                    k.name, p.iterations, p.evaluations, p.error.empty() ? "" : ", stopped: ",
                    p.error.c_str(), p.total, 100.0 * p.assembly / p.total,
                    100.0 * p.update / p.total, 100.0 * linear / p.total);
        std::fflush(stdout);
    }
}

std::vector<double> uniform_axis(double length, int nodes) {
    std::vector<double> a(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) a[static_cast<std::size_t>(i)] = length * i / (nodes - 1);
    return a;
}

// Unit 44's LG4 diode (bias_test.cpp): p below x = 1 um, n above, 1e17, legacy silicon.
device::Device diode(mesh::Mesh m, const physics::Semiconductor& si, double NA, double ND,
                     double junction_cm) {
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        (m.points()[i][0] < junction_cm ? acceptors[i] : donors[i]) =
            m.points()[i][0] < junction_cm ? NA : ND;
    }
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    return *device::Device::create({.mesh = std::move(m),
                                    .temperature_K = 300.0,
                                    .regions = {{"silicon", si}},
                                    .node_region = std::vector<device::RegionId>(n, 0),
                                    .donors = std::move(donors),
                                    .acceptors = std::move(acceptors),
                                    .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                                                 {"cathode", device::ContactKind::ohmic,
                                                  std::move(cathode)}}});
}

device::Device diode_3d(int nodes) {
    const auto a = uniform_axis(2e-4, nodes), b = uniform_axis(1e-4, nodes);
    return diode(*mesh::make_tensor_grid(a, b, b), physics::legacy_silicon(), 1e17, 1e17, 1e-4);
}

// A deck's device, options and first biased point.
Case deck_case(const std::string& name, const fs::path& path) {
    const auto d = deck::read(path);
    if (!d) FAIL(name << ": " << d.error().message);
    auto p = deck::build(*d);
    if (!p) FAIL(name << ": " << p.error().message);
    std::vector<double> bias = p->points.front();
    for (const auto& point : p->points) {
        if (std::ranges::any_of(point, [](double v) { return v != 0.0; })) {
            bias = point;
            break;
        }
    }
    return {name, std::move(p->device), p->options, std::move(bias)};
}

}  // namespace

TEST_CASE("solver profile: the phase split of a Newton iteration (Unit 58 SD1)",
          "[.scale_profile]") {
    // SparseLU only where it fits: up to 81k unknowns (30^3 diode: Unit 3 measured 70 s per
    // factorization), and with at most 3 iterations above 50k.
    SECTION("2D MOSFET example") {
        profile(deck_case("2D MOSFET example", fs::path(NITCAD_EXAMPLES_DIR) / "mosfet_2d.deck"),
                true);
    }
    SECTION("3D diode 30^3") {
        profile({"3D diode 30^3", diode_3d(30), {}, {0.5, 0.0}}, false);
        // SparseLU's split alone, capped (each factorization about a minute).
        const Case c{"3D diode 30^3, SparseLU", diode_3d(30), {}, {0.5, 0.0}};
        profile(c, true, 3);
    }
    SECTION("3D diode 40^3") { profile({"3D diode 40^3", diode_3d(40), {}, {0.5, 0.0}}, false); }
    SECTION("3D energy balance") {
        // EC19's pn diode (p 1e17 / n 2e17, 1 um, published silicon) on 41 x 21 x 21, at -0.4 V on
        // the cathode, energy balance (its own Canali law; field_mobility is refused with it).
        const auto x = uniform_axis(1e-4, 41), y = uniform_axis(0.2e-4, 21);
        solve::BiasOptions o;
        o.models.energy_balance = true;
        profile({"3D energy balance diode 41x21x21",
                 diode(*mesh::make_tensor_grid(x, y, y), physics::silicon(), 1e17, 2e17, 0.5e-4),
                 o,
                 {0.0, -0.4}},
                false);
    }
    SECTION("reduced swept FinFET") {
        profile(deck_case("reduced swept FinFET",
                          fs::path(NITCAD_DECKS_DIR) / ".." / "tests" / "decks" /
                              "finfet_3d_swept_reduced.deck"),
                false);
    }
    SECTION("tensor FinFET deck") {
        profile(deck_case("tensor FinFET deck", fs::path(NITCAD_DECKS_DIR) / "finfet_3d.deck"),
                false);
    }
}

TEST_CASE("solver profile: one Newton step's linear solve on a 3D diode size series (Unit 58 SD1)",
          "[.scale_series]") {
    // LG4's measurement (bias_test.cpp) extended towards 1e6 unknowns: the Jacobian at 0.5 V on
    // the equilibrium potential, a fresh solver (analysis, factorization or ILUT, one solve).
    // GMRES-ILUT strict (no fallback) and PARDISO at 6 threads at every size, PARDISO at 1 thread
    // up to 40^3. Cheapest memory first in each size.
    const bool mkl = mkl_loaded();
    for (const int nodes : {20, 30, 40, 50, 60, 70}) {
        const auto d = diode_3d(nodes);
        const auto t0 = clock_type::now();
        const auto eq = solve::solve_equilibrium(d, {.linear = helper_solver()});
        REQUIRE(eq.has_value());
        const double eq_s = since(t0);
        std::vector<double> rhs;
        const linalg::SparseMatrix j =
            device_jacobian(d, std::vector<double>{0.5, 0.0}, eq->fields.potential_V, rhs);
        std::printf("\n=== size series %d^3: %zu unknowns, %zu nonzeros (equilibrium %.1f s)\n",
                    nodes, static_cast<std::size_t>(j.rows()), j.nonzeros(), eq_s);
        std::fflush(stdout);
        std::vector<Configuration> runs{
            {"GMRES-ILUT", {.backend = linalg::SolverBackend::gmres_ilut,
                            .fallback = linalg::LinearFallback::none}}};
        if (mkl) {
            if (nodes <= 40) {
                runs.push_back({"PARDISO 1 thread", {.backend = linalg::SolverBackend::mkl_pardiso}});
            }
            runs.push_back({"PARDISO 6 threads",
                            {.backend = linalg::SolverBackend::mkl_pardiso, .threads = 6}});
        }
        for (const Configuration& k : runs) {
            const ProcessMemory before = process_memory();
            const auto s0 = clock_type::now();
            auto solver = *linalg::LinearSolver::create(k.config);
            std::vector<double> x(rhs.size());
            std::string error;
            std::size_t iterations = 0;
            double backward = 0.0, factor_s = 0.0, held = 0.0;
            std::optional<std::size_t> fill;
            if (auto f = solver.factorize(j); !f) {
                error = f.error().message;
            } else {
                fill = f->fill;
                factor_s = since(s0);
                held = process_memory().private_MB - before.private_MB;
                if (auto r = solver.solve(rhs, x); !r) {
                    error = r.error().message;
                } else {
                    iterations = r->iterations.value_or(0);
                    backward = r->backward_error;
                }
            }
            const double total = since(s0);
            std::printf("  %-18s %.2f s (factorization or ILUT %.2f s), %zu GMRES products, "
                        "backward error %.1e, fill %zu, held %.0f MB, peak rise %.0f MB%s%s\n",
                        k.name, total, factor_s, iterations, backward, fill.value_or(0), held,
                        process_memory().peak_working_set_MB - before.peak_working_set_MB,
                        error.empty() ? "" : "; failed: ", error.c_str());
            std::fflush(stdout);
        }
    }
}

TEST_CASE("solver profile: GMRES-ILUT with the thread pool (Unit 58 slice 1)", "[.scale_threads]") {
    // One Newton step's linear solve as [.scale_series], GMRES-ILUT strict at 1 and 6 threads: the
    // ILUT build (serial), the GMRES solve (pooled products and updates) and whether the bits agree.
    for (const int nodes : {40, 60}) {
        const auto d = diode_3d(nodes);
        const auto eq = solve::solve_equilibrium(d, {.linear = helper_solver()});
        REQUIRE(eq.has_value());
        std::vector<double> rhs;
        const linalg::SparseMatrix j =
            device_jacobian(d, std::vector<double>{0.5, 0.0}, eq->fields.potential_V, rhs);
        std::vector<double> first;
        for (const int threads : {1, 6}) {
            auto solver = *linalg::LinearSolver::create({.backend = linalg::SolverBackend::gmres_ilut,
                                                         .threads = threads,
                                                         .fallback = linalg::LinearFallback::none});
            auto t0 = clock_type::now();
            REQUIRE(solver.factorize(j).has_value());
            const double build = since(t0);
            std::vector<double> x(rhs.size());
            t0 = clock_type::now();
            const auto r = solver.solve(rhs, x);
            const double solve = since(t0);
            REQUIRE(r.has_value());
            if (first.empty()) first = x;
            std::size_t differ = 0;
            for (std::size_t i = 0; i < x.size(); ++i) differ += std::bit_cast<std::uint64_t>(x[i]) !=
                                                                 std::bit_cast<std::uint64_t>(first[i]);
            std::printf("threads %d^3 (%zu unknowns), %d threads: ILUT %.2f s, GMRES solve %.2f s "
                        "(%zu products), entries differing from 1 thread: %zu\n",
                        nodes, x.size(), threads, build, solve, r->iterations.value_or(0), differ);
            std::fflush(stdout);
            REQUIRE(differ == 0);
        }
    }
}

namespace {

// One linear solve with a fresh solver: build (analysis and factorization or preconditioner), solve,
// GMRES products, held memory; a failure is printed, not fatal (strict runs show their failures).
void one_solve(const char* name, const linalg::SparseMatrix& a, std::span<const double> b,
               const linalg::SolverConfig& config,
               std::span<const std::int32_t> labels = {}) {
    const ProcessMemory before = process_memory();
    auto solver = *linalg::LinearSolver::create(config);
    solver.set_field_labels(labels);
    std::vector<double> x(b.size());
    std::string error;
    std::size_t iterations = 0;
    double backward = 0.0, build_s = 0.0, solve_s = 0.0, held = 0.0;
    std::optional<std::size_t> fill;
    auto t0 = clock_type::now();
    if (auto f = solver.factorize(a); !f) {
        error = f.error().message;
    } else {
        build_s = since(t0);
        fill = f->fill;
        held = process_memory().private_MB - before.private_MB;
        t0 = clock_type::now();
        if (auto r = solver.solve(b, x); !r) {
            error = r.error().message;
        } else {
            iterations = r->iterations.value_or(0);
            backward = r->backward_error;
        }
        solve_s = since(t0);
    }
    std::printf("  %-24s build %.2f s, solve %.2f s, total %.2f s, %zu GMRES products, backward "
                "error %.1e, fill %zu, held %.0f MB%s%s\n",
                name, build_s, solve_s, build_s + solve_s, iterations, backward, fill.value_or(0),
                held, error.empty() ? "" : "; failed: ", error.c_str());
    std::fflush(stdout);
}

// The equilibrium Poisson Jacobian (Formulation::equilibrium, one unknown per node) at the given
// potential (V), and its residual as the right-hand side.
linalg::SparseMatrix equilibrium_jacobian(const device::Device& d, std::span<const double> potential_V,
                                          std::vector<double>& rhs) {
    const auto scaling = *assemble::make_scaling(d);
    auto system = *assemble::create_equilibrium(d, scaling);
    std::vector<double> psi(potential_V.size());
    for (std::size_t i = 0; i < psi.size(); ++i) psi[i] = potential_V[i] / scaling.V_T;
    linalg::SparseMatrix j = system.make_jacobian();
    rhs.assign(psi.size(), 0.0);
    system.evaluate(psi, rhs, j);
    return j;
}

}  // namespace

TEST_CASE("solver profile: GMRES-AMG on the 3D diode series (Unit 58 slice 2)", "[.scale_amg]") {
    // SG4's AMG part: one solve of the equilibrium Poisson Jacobian, at the charge-neutral guess (the
    // first Newton step) and at the converged potential, GMRES-AMG against GMRES-ILUT, strict. Then
    // SD3's comparison: GMRES-AMG on the whole drift-diffusion Jacobian at 0.5 V (as
    // [.scale_series]), up to 40^3. Threads: AMG at 6 threads on the converged Poisson system.
    const linalg::SolverConfig ilut{.backend = linalg::SolverBackend::gmres_ilut,
                                    .fallback = linalg::LinearFallback::none};
    linalg::SolverConfig amg = ilut;
    amg.preconditioner = linalg::Preconditioner::amg;
    linalg::SolverConfig amg6 = amg;
    amg6.threads = 6;
    for (const int nodes : {20, 30, 40, 50, 60, 70}) {
        const auto d = diode_3d(nodes);
        const auto eq = solve::solve_equilibrium(d, {.linear = helper_solver()});
        REQUIRE(eq.has_value());
        std::vector<double> guess_V;
        {
            const auto scaling = *assemble::make_scaling(d);
            const auto system = *assemble::create_equilibrium(d, scaling);
            for (const double psi : assemble::charge_neutral_potential(system)) {
                guess_V.push_back(psi * scaling.V_T);
            }
        }
        std::vector<double> rhs;
        const linalg::SparseMatrix first = equilibrium_jacobian(d, guess_V, rhs);
        std::printf("\n=== amg %d^3: Poisson %zu unknowns, %zu nonzeros; first Newton step\n", nodes,
                    static_cast<std::size_t>(first.rows()), first.nonzeros());
        one_solve("GMRES-ILUT", first, rhs, ilut);
        one_solve("GMRES-AMG", first, rhs, amg);
        std::vector<double> rhs_eq;
        const linalg::SparseMatrix converged = equilibrium_jacobian(d, eq->fields.potential_V, rhs_eq);
        // At the converged state the residual is at rounding level: solve for a unit-size one.
        for (std::size_t i = 0; i < rhs_eq.size(); ++i) {
            rhs_eq[i] = std::sin(0.37 * static_cast<double>(i)) + 1e-3;
        }
        std::printf("  converged potential:\n");
        one_solve("GMRES-ILUT", converged, rhs_eq, ilut);
        one_solve("GMRES-AMG", converged, rhs_eq, amg);
        one_solve("GMRES-AMG 6 threads", converged, rhs_eq, amg6);
        if (nodes <= 40) {
            std::vector<double> rhs_dd;
            const linalg::SparseMatrix j =
                device_jacobian(d, std::vector<double>{0.5, 0.0}, eq->fields.potential_V, rhs_dd);
            std::printf("  drift-diffusion at 0.5 V, %zu unknowns:\n", static_cast<std::size_t>(j.rows()));
            one_solve("GMRES-ILUT", j, rhs_dd, ilut);
            one_solve("GMRES-AMG (whole system)", j, rhs_dd, amg);
        }
    }
}

TEST_CASE("solver profile: GMRES-CPR on the 3D diode series (Unit 58 slice 3)", "[.scale_cpr]") {
    // SG4's CPR part, SG5 and SG7: one Newton step's linear solve as [.scale_series] (the Jacobian at
    // 0.5 V on the equilibrium potential, a fresh solver), GMRES-CPR against GMRES-ILUT, strict, at
    // one thread; SG6's linear part: GMRES-CPR at 6 threads. Then [.scale_amg]'s 70^3 six-thread
    // anomaly, repeated.
    const linalg::SolverConfig ilut{.backend = linalg::SolverBackend::gmres_ilut,
                                    .fallback = linalg::LinearFallback::none};
    linalg::SolverConfig cpr = ilut;
    cpr.preconditioner = linalg::Preconditioner::cpr;
    linalg::SolverConfig cpr6 = cpr;
    cpr6.threads = 6;
    for (const int nodes : {20, 30, 40, 50, 60, 70}) {
        const auto d = diode_3d(nodes);
        const auto eq = solve::solve_equilibrium(d, {.linear = helper_solver()});
        REQUIRE(eq.has_value());
        std::vector<double> rhs;
        const linalg::SparseMatrix j =
            device_jacobian(d, std::vector<double>{0.5, 0.0}, eq->fields.potential_V, rhs);
        REQUIRE(j.rows() % 3 == 0);
        std::vector<std::int32_t> labels(static_cast<std::size_t>(j.rows()));
        for (std::size_t u = 0; u < labels.size(); ++u) labels[u] = static_cast<std::int32_t>(u % 3);
        std::printf("\n=== cpr %d^3: %zu unknowns, %zu nonzeros\n", nodes,
                    static_cast<std::size_t>(j.rows()), j.nonzeros());
        one_solve("GMRES-ILUT", j, rhs, ilut);
        one_solve("GMRES-CPR", j, rhs, cpr, labels);
        one_solve("GMRES-CPR 6 threads", j, rhs, cpr6, labels);
    }
    {
        const auto d = diode_3d(70);
        const auto eq = solve::solve_equilibrium(d, {.linear = helper_solver()});
        REQUIRE(eq.has_value());
        std::vector<double> rhs;
        const linalg::SparseMatrix converged = equilibrium_jacobian(d, eq->fields.potential_V, rhs);
        for (std::size_t i = 0; i < rhs.size(); ++i) rhs[i] = std::sin(0.37 * static_cast<double>(i)) + 1e-3;
        linalg::SolverConfig amg = ilut;
        amg.preconditioner = linalg::Preconditioner::amg;
        linalg::SolverConfig amg6 = amg;
        amg6.threads = 6;
        std::printf("\n=== amg 70^3 converged Poisson, repeated\n");
        for (int repeat = 0; repeat < 3; ++repeat) {
            one_solve("GMRES-AMG", converged, rhs, amg);
            one_solve("GMRES-AMG 6 threads", converged, rhs, amg6);
        }
    }
}

TEST_CASE("solver profile: the carrier blocks alone on the 3D diode series (Unit 63 slice 0)",
          "[.scale_carrier_blocks]") {
    // MD1: each field's diagonal block of the 0.5 V Jacobian on the equilibrium potential (as
    // [.scale_cpr]), unequilibrated (the symmetrization absorbs any row and column scaling), solved
    // by GMRES(30) to a backward error of 1e-8 with the slice-2 AMG on the block as it is and with
    // the AMG on its symmetrization (MD2). One thread.
    for (const int nodes : {20, 30, 40, 50, 60, 70}) {
        const auto d = diode_3d(nodes);
        const auto eq = solve::solve_equilibrium(d, {.linear = helper_solver()});
        REQUIRE(eq.has_value());
        std::vector<double> rhs;
        const linalg::SparseMatrix j =
            device_jacobian(d, std::vector<double>{0.5, 0.0}, eq->fields.potential_V, rhs);
        REQUIRE(j.rows() % 3 == 0);
        std::printf("\n=== carrier blocks %d^3: %zu unknowns\n", nodes, static_cast<std::size_t>(j.rows()));
        const char* names[] = {"psi block", "electron block", "hole block"};
        for (const std::size_t field : {std::size_t{0}, std::size_t{1}, std::size_t{2}}) {
            (void)block_probe::compare(names[field], block_probe::block(j, 3, field));
        }
    }
}

namespace {

// AMG on S = D A D^-1's symmetric part, with a chosen near-null vector (1/d, d or 1), as a
// preconditioner of A (Unit 63 slice 0 diagnostic).
struct CandidateAmg {
    int mode;  // 0: 1/d, 1: d, 2: 1
    linalg::detail::AmgBackend amg{block_probe::default_amg(), nullptr};
    linalg::detail::Csr a;
    std::vector<double> d, b, x;
    std::expected<void, base::Error> analyze(const linalg::SparseMatrix& m) {
        a = block_probe::to_csr(m);
        return amg.analyze(m);
    }
    std::expected<linalg::detail::Factorization, base::Error> factorize(std::span<const double>) {
        d = linalg::detail::amg::symmetrize(a).d;
        const std::vector<double> s = linalg::detail::amg::symmetric_part(a, d);
        std::vector<double> inv(d.size());
        for (std::size_t i = 0; i < d.size(); ++i) inv[i] = 1.0 / d[i];
        if (mode == 0) amg.set_scales(d, d, true);
        if (mode == 1) amg.set_scales(inv, inv, true);
        b.assign(d.size(), 0.0);
        x.assign(d.size(), 0.0);
        return amg.factorize(s);
    }
    void solve(std::span<const double> r, std::span<double> z) {
        for (std::size_t i = 0; i < r.size(); ++i) b[i] = d[i] * r[i];
        amg.solve(b, x);
        for (std::size_t i = 0; i < z.size(); ++i) z[i] = x[i] / d[i];
    }
};

}  // namespace

TEST_CASE("solver profile: the near-null vector of the symmetrized blocks (Unit 63 slice 0)",
          "[.carrier_diag]") {
    for (const int nodes : {20, 40}) {
        const auto dev = diode_3d(nodes);
        const auto eq = solve::solve_equilibrium(dev, {.linear = helper_solver()});
        REQUIRE(eq.has_value());
        std::vector<double> rhs;
        const linalg::SparseMatrix j =
            device_jacobian(dev, std::vector<double>{0.5, 0.0}, eq->fields.potential_V, rhs);
        std::printf("\n=== %d^3\n", nodes);
        const char* names[] = {"psi", "electrons", "holes"};
        for (const std::size_t field : {std::size_t{0}, std::size_t{1}, std::size_t{2}}) {
            const linalg::SparseMatrix blk = block_probe::block(j, 3, field);
            const auto csr = block_probe::to_csr(blk);
            const auto s = linalg::detail::amg::symmetrize(csr);
            double lo = 0.0, hi = 0.0, asym = 0.0, rowsum = 0.0, diag_neg = 0.0, offdiag_pos = 0.0;
            for (const double v : s.d) {
                lo = std::min(lo, std::log(v));
                hi = std::max(hi, std::log(v));
            }
            for (std::size_t i = 0; i < csr.rows; ++i) {
                double sum = 0.0, scale = 0.0;
                for (auto k = static_cast<std::size_t>(csr.offsets[i]); k < static_cast<std::size_t>(csr.offsets[i + 1]); ++k) {
                    const auto c = static_cast<std::size_t>(csr.columns[k]);
                    const double v = csr.values[k];
                    sum += v;
                    scale += std::abs(v);
                    if (c == i && v < 0) diag_neg += 1;
                    if (c != i && v > 0) offdiag_pos += 1;
                    if (c != i) {
                        for (auto t = static_cast<std::size_t>(csr.offsets[c]); t < static_cast<std::size_t>(csr.offsets[c + 1]); ++t) {
                            if (static_cast<std::size_t>(csr.columns[t]) == i) {
                                asym = std::max(asym, std::abs(v - csr.values[t]) / std::max(std::abs(v), std::abs(csr.values[t])));
                            }
                        }
                    }
                }
                if (scale > 0) rowsum = std::max(rowsum, std::abs(sum) / scale);
            }
            std::printf("  %s: ln d in [%.2f, %.2f], max rel asymmetry %.3g, max |row sum|/row l1 %.3g, "
                        "negative diagonals %.0f, positive off-diagonals %.0f\n",
                        names[field], lo, hi, asym, rowsum, diag_neg, offdiag_pos);
            const std::vector<double> c = block_probe::rhs(static_cast<std::size_t>(blk.rows()));
            const char* modes[] = {"near-null 1/d", "near-null d", "near-null 1"};
            for (int mode = 0; mode < 3; ++mode) {
                CandidateAmg m{mode};
                block_probe::print(modes[mode], block_probe::run(blk, c, m));
            }
        }
    }
}

TEST_CASE("solver profile: the stalled coarsening of the symmetrized carrier blocks (Unit 63 slice 0b)",
          "[.carrier_0b]") {
    // MB1 (diagnosis, default settings, 40^3-70^3) and MB2's candidates: theta 0.08 (default), 0.04,
    // 0.02, 0, each with the stalled level smoothed (as slice 2) or factorized (factorize_stalled).
    const double thetas[] = {0.08, 0.04, 0.02, 0.0};
    const char* names[] = {"psi", "electrons", "holes"};
    for (const int nodes : {20, 30, 40, 50, 60, 70}) {
        const auto dev = diode_3d(nodes);
        const auto eq = solve::solve_equilibrium(dev, {.linear = helper_solver()});
        REQUIRE(eq.has_value());
        std::vector<double> rhs;
        const linalg::SparseMatrix j =
            device_jacobian(dev, std::vector<double>{0.5, 0.0}, eq->fields.potential_V, rhs);
        std::printf("\n=== %d^3: %zu unknowns\n", nodes, static_cast<std::size_t>(j.rows()));
        for (const std::size_t field : {std::size_t{0}, std::size_t{1}, std::size_t{2}}) {
            const linalg::SparseMatrix blk = block_probe::block(j, 3, field);
            std::printf("  %s\n", names[field]);
            for (const bool stalled : {false, true}) {
                for (const double theta : thetas) {
                    linalg::detail::AmgSettings s = block_probe::default_amg();
                    s.strength_threshold = theta;
                    s.factorize_stalled = stalled;
                    linalg::detail::SymmetrizedAmg* m = nullptr;
                    const auto o = block_probe::symmetrized(blk, s, &m);
                    char label[64];
                    std::snprintf(label, sizeof label, "theta %.2f %s", theta, stalled ? "factorized" : "smoothed");
                    block_probe::print(label, o);
                    if (nodes >= 40 && field > 0 && theta == 0.08 && !stalled && o.built) {
                        block_probe::diagnose(m->amg(), s);
                    }
                }
            }
        }
    }
}
