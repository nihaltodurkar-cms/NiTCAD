#include "NiTCAD/solve/equilibrium.hpp"

#include <cstddef>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/equilibrium_poisson.hpp"
#include "NiTCAD/assemble/scaling.hpp"

namespace NiTCAD::solve {

std::expected<results::EquilibriumResult, base::Error> solve_equilibrium(
    const device::Device& device, const EquilibriumOptions& options, const RunControl& control) {
    auto scaling = assemble::make_scaling(device, options.Ns_override);
    if (!scaling) return std::unexpected(std::move(scaling.error()));
    auto system = assemble::EquilibriumPoisson::create(device, *scaling, options.models);
    if (!system) return std::unexpected(std::move(system.error()));
    auto solver = linalg::LinearSolver::create(options.linear);
    if (!solver) return std::unexpected(std::move(solver.error()));

    std::vector<double> psi = system->charge_neutral_potential();
    results::EquilibriumResult result;
    const IterationObserver observe = [&](const results::IterationRecord& r, bool converged) {
        if (control.progress) {
            control.progress(
                {Phase::equilibrium, 0, 1, r.iteration, r.update, r.residual, converged});
        }
    };
    if (auto ok = newton_solve(*system, psi, options.newton, *solver, result.convergence,
                               control.stop, observe);
        !ok) {
        return std::unexpected(std::move(ok.error()));
    }

    const std::size_t n = system->unknowns();
    result.fields = {std::vector<double>(n), std::vector<double>(n), std::vector<double>(n)};
    system->carriers(psi, result.fields.n_cm3, result.fields.p_cm3);
    for (std::size_t i = 0; i < n; ++i) {
        result.fields.potential_V[i] = psi[i] * scaling->V_T;
        result.fields.n_cm3[i] *= scaling->Ns;
        result.fields.p_cm3[i] *= scaling->Ns;
    }
    return result;
}

}  // namespace NiTCAD::solve
