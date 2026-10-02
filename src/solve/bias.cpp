#include "NiTCAD/solve/bias.hpp"

#include <cmath>
#include <cstddef>
#include <utility>

#include "NiTCAD/solve/equilibrium.hpp"

namespace NiTCAD::solve {

namespace {

bool valid_state(const DeviceState& s, std::size_t n) {
    if (s.potential_V.size() != n || s.n_cm3.size() != n || s.p_cm3.size() != n) return false;
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(s.potential_V[i])) return false;
        if (!(std::isfinite(s.n_cm3[i]) && s.n_cm3[i] > 0.0)) return false;
        if (!(std::isfinite(s.p_cm3[i]) && s.p_cm3[i] > 0.0)) return false;
    }
    return true;
}

}  // namespace

std::expected<BiasSolution, base::Error> solve_bias(const device::Device& device,
                                                    std::span<const double> bias_V,
                                                    const BiasOptions& options,
                                                    const DeviceState* initial) {
    auto scaling = assemble::make_scaling(device, options.Ns_override);
    if (!scaling) return std::unexpected(std::move(scaling.error()));
    auto system = assemble::DriftDiffusion::create(device, *scaling, options.models);
    if (!system) return std::unexpected(std::move(system.error()));
    if (auto ok = system->set_bias(bias_V); !ok) return std::unexpected(std::move(ok.error()));
    auto solver = linalg::LinearSolver::create(options.linear);
    if (!solver) return std::unexpected(std::move(solver.error()));

    const std::size_t nodes = system->node_count();
    std::vector<double> x;
    if (initial != nullptr) {
        if (!valid_state(*initial, nodes)) {
            return std::unexpected(base::Error{
                base::ErrorCode::invalid_input,
                "initial state needs a finite potential and positive densities for every node",
                std::nullopt});
        }
        x.resize(system->unknowns());
        for (std::size_t i = 0; i < nodes; ++i) {
            x[3 * i] = initial->potential_V[i] / scaling->V_T;
            x[3 * i + 1] = initial->n_cm3[i] / scaling->Ns;
            x[3 * i + 2] = initial->p_cm3[i] / scaling->Ns;
        }
        system->stamp_contacts(x);
    } else {
        auto equilibrium = solve_equilibrium(
            device, {.newton = options.newton, .linear = options.linear,
                     .Ns_override = options.Ns_override});
        if (!equilibrium) return std::unexpected(std::move(equilibrium.error()));
        std::vector<double> psi(nodes);
        for (std::size_t i = 0; i < nodes; ++i) psi[i] = equilibrium->potential_V[i] / scaling->V_T;
        x = system->state_from_potential(psi);
    }

    auto report = newton_solve(*system, x, options.newton, *solver);
    if (!report) return std::unexpected(std::move(report.error()));

    const double current_scale =
        scaling->J0 * std::pow(scaling->L_D, device.mesh().dimension() - 1);
    BiasSolution solution{.scaling = *scaling,
                          .bias_V = std::vector<double>(bias_V.begin(), bias_V.end()),
                          .state = {std::vector<double>(nodes), std::vector<double>(nodes),
                                    std::vector<double>(nodes)},
                          .terminal_current = system->terminal_currents(x),
                          .edge_current_n = {},
                          .edge_current_p = {},
                          .newton = std::move(*report)};
    for (std::size_t i = 0; i < nodes; ++i) {
        solution.state.potential_V[i] = x[3 * i] * scaling->V_T;
        solution.state.n_cm3[i] = x[3 * i + 1] * scaling->Ns;
        solution.state.p_cm3[i] = x[3 * i + 2] * scaling->Ns;
    }
    for (double& I : solution.terminal_current) I *= current_scale;
    for (const auto& [jn, jp] : system->edge_currents(x)) {
        solution.edge_current_n.push_back(jn * current_scale);
        solution.edge_current_p.push_back(jp * current_scale);
    }
    return solution;
}

}  // namespace NiTCAD::solve
