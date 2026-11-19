#include "NiTCAD/solve/equilibrium.hpp"

#include <cstddef>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/equilibrium.hpp"
#include "NiTCAD/assemble/quantum_confinement.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "density_gradient.hpp"
#include "fields.hpp"
#include "observe.hpp"

namespace NiTCAD::solve {

namespace {

// Unit 62 (QP8, QP9): the density-gradient equilibrium, from the classical one, its results read from
// the drift-diffusion state it slaves (zero bias).
std::expected<results::EquilibriumResult, base::Error> quantum_equilibrium(
    const device::Device& device, const EquilibriumOptions& options, const RunControl& control,
    const assemble::Scaling& scaling) {
    auto solver = linalg::LinearSolver::create(options.linear);
    if (!solver) return std::unexpected(std::move(solver.error()));
    const std::vector<double> zero(device.contacts().size(), 0.0);
    results::EquilibriumResult result;
    const IterationObserver observe = [&](const results::IterationRecord& r, bool converged) {
        if (control.progress) {
            control.progress({Phase::equilibrium, 0, 1, r.iteration, r.update, r.residual, converged});
        }
    };
    auto psi = detail::classical_equilibrium(device, scaling, options.models, zero, options.newton,
                                             *solver, result.convergence, control.stop, observe);
    if (!psi) return std::unexpected(std::move(psi.error()));
    auto q = detail::solve_quantum_equilibrium(device, scaling, options.models, zero, *psi,
                                               options.newton, *solver, result.convergence,
                                               control.stop, observe);
    if (!q) return std::unexpected(std::move(q.error()));
    result.quantum_continuation = q->rungs;
    auto dd = assemble::System::create(device, scaling, options.models);
    if (!dd) return std::unexpected(std::move(dd.error()));
    if (auto ok = dd->set_bias(zero); !ok) return std::unexpected(std::move(ok.error()));
    const auto x = dd->state_from_quantum_potential(q->psi, q->lambda_n, q->lambda_p);
    const std::size_t n = device.mesh().node_count(), m = dd->stride();
    result.fields = {std::vector<double>(n), std::vector<double>(n), std::vector<double>(n)};
    for (std::size_t i = 0; i < n; ++i) {
        result.fields.potential_V[i] = x[m * i] * scaling.V_T;
        result.fields.n_cm3[i] = x[m * i + 1] * scaling.Ns;
        result.fields.p_cm3[i] = x[m * i + 2] * scaling.Ns;
    }
    detail::quantum_potentials(*dd, x, scaling.V_T, result.fields);
    result.bands = detail::band_diagram(*dd, x, scaling.V_T);
    const double charge = detail::charge_scale(scaling, device.mesh().dimension());
    result.gate_charge = gate_charges(*dd, x);
    for (double& Q : result.gate_charge) Q *= charge;
    result.interface_trap_charge = interface_trap_charges(*dd, x);
    for (double& Q : result.interface_trap_charge) Q *= charge;
    return result;
}

std::expected<results::EquilibriumResult, base::Error> solve_equilibrium_unobserved(
    const device::Device& device, const EquilibriumOptions& options, const RunControl& control) {
    auto scaling = assemble::make_scaling(device, options.Ns_override);
    if (!scaling) return std::unexpected(std::move(scaling.error()));
    if (options.models.quantum_confinement == assemble::QuantumConfinement::density_gradient) {
        return quantum_equilibrium(device, options, control, *scaling);
    }
    auto system = assemble::System::create(device, *scaling, options.models,
                                           assemble::Formulation::equilibrium);
    if (!system) return std::unexpected(std::move(system.error()));
    auto solver = linalg::LinearSolver::create(options.linear);
    if (!solver) return std::unexpected(std::move(solver.error()));

    std::vector<double> psi = charge_neutral_potential(*system);
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
    result.bands = detail::band_diagram(*system, psi, scaling->V_T);
    result.gate_charge = gate_charges(*system, psi);
    const double charge = detail::charge_scale(*scaling, device.mesh().dimension());
    for (double& Q : result.gate_charge) Q *= charge;
    result.interface_trap_charge = interface_trap_charges(*system, psi);
    for (double& Q : result.interface_trap_charge) Q *= charge;
    return result;
}

}  // namespace

// The public solve: the body above, its result handed to the test observer (observe.hpp).
std::expected<results::EquilibriumResult, base::Error> solve_equilibrium(
    const device::Device& device, const EquilibriumOptions& options, const RunControl& control) {
    auto r = solve_equilibrium_unobserved(device, options, control);
    if (r) detail::observe_result(*r);
    return r;
}

}  // namespace NiTCAD::solve
