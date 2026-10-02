#include "NiTCAD/solve/equilibrium.hpp"

#include <cstddef>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/equilibrium_poisson.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "fields.hpp"

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

    result.fields = detail::equilibrium_fields(*system, psi, *scaling);
    result.gate_charge = system->gate_charges(psi);
    const double charge = detail::charge_scale(*scaling, device.mesh().dimension());
    for (double& Q : result.gate_charge) Q *= charge;
    return result;
}

}  // namespace NiTCAD::solve
