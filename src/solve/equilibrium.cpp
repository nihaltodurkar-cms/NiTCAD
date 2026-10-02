#include "NiTCAD/solve/equilibrium.hpp"

#include <cstddef>
#include <utility>

#include "NiTCAD/assemble/equilibrium_poisson.hpp"

namespace NiTCAD::solve {

std::expected<EquilibriumSolution, base::Error> solve_equilibrium(
    const device::Device& device, const EquilibriumOptions& options) {
    auto scaling = assemble::make_scaling(device, options.Ns_override);
    if (!scaling) return std::unexpected(std::move(scaling.error()));
    auto system = assemble::EquilibriumPoisson::create(device, *scaling);
    if (!system) return std::unexpected(std::move(system.error()));
    auto solver = linalg::LinearSolver::create(options.linear);
    if (!solver) return std::unexpected(std::move(solver.error()));

    std::vector<double> psi = system->charge_neutral_potential();
    auto report = newton_solve(*system, psi, options.newton, *solver);
    if (!report) return std::unexpected(std::move(report.error()));

    const std::size_t n = system->unknowns();
    EquilibriumSolution solution{.scaling = *scaling,
                                 .potential_V = std::vector<double>(n),
                                 .n_cm3 = std::vector<double>(n),
                                 .p_cm3 = std::vector<double>(n),
                                 .newton = std::move(*report)};
    system->carriers(psi, solution.n_cm3, solution.p_cm3);
    for (std::size_t i = 0; i < n; ++i) {
        solution.potential_V[i] = psi[i] * scaling->V_T;
        solution.n_cm3[i] *= scaling->Ns;
        solution.p_cm3[i] *= scaling->Ns;
    }
    return solution;
}

}  // namespace NiTCAD::solve
