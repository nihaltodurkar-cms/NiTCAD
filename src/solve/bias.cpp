#include "NiTCAD/solve/bias.hpp"

#include <cmath>
#include <cstddef>
#include <string>
#include <utility>

#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/solve/equilibrium.hpp"

namespace NiTCAD::solve {

namespace {

base::Error invalid(std::string message, std::optional<std::size_t> index = std::nullopt) {
    return {base::ErrorCode::invalid_input, std::move(message),
            base::ErrorContext{.index = index, .value = std::nullopt}};
}

bool valid_fields(const results::NodeFields& s, std::size_t n) {
    if (s.potential_V.size() != n || s.n_cm3.size() != n || s.p_cm3.size() != n) return false;
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(s.potential_V[i])) return false;
        if (!(std::isfinite(s.n_cm3[i]) && s.n_cm3[i] > 0.0)) return false;
        if (!(std::isfinite(s.p_cm3[i]) && s.p_cm3[i] > 0.0)) return false;
    }
    return true;
}

}  // namespace

std::expected<results::Sweep, base::Error> sweep_bias(const device::Device& device,
                                                      std::span<const std::vector<double>> points,
                                                      const BiasOptions& options,
                                                      const results::NodeFields* initial,
                                                      const RunControl& control) {
    if (points.empty()) return std::unexpected(invalid("a sweep needs at least one bias point"));
    const std::size_t contacts = device.contacts().size();
    for (std::size_t k = 0; k < points.size(); ++k) {
        if (points[k].size() != contacts) {
            return std::unexpected(invalid("a bias point needs one bias per contact", k));
        }
        for (const double v : points[k]) {
            if (!std::isfinite(v)) return std::unexpected(invalid("a bias is not finite", k));
        }
    }
    const std::size_t nodes = device.mesh().node_count();
    if (initial != nullptr && !valid_fields(*initial, nodes)) {
        return std::unexpected(invalid(
            "initial state needs a finite potential and positive densities for every node"));
    }
    auto scaling = assemble::make_scaling(device, options.Ns_override);
    if (!scaling) return std::unexpected(std::move(scaling.error()));
    auto system = assemble::DriftDiffusion::create(device, *scaling, options.models);
    if (!system) return std::unexpected(std::move(system.error()));
    auto solver = linalg::LinearSolver::create(options.linear);
    if (!solver) return std::unexpected(std::move(solver.error()));

    results::Sweep sweep;
    sweep.run = make_run_record(device, options, points, initial);

    // Starting state, scaled.
    std::vector<double> x(system->unknowns());
    if (initial != nullptr) {
        for (std::size_t i = 0; i < nodes; ++i) {
            x[3 * i] = initial->potential_V[i] / scaling->V_T;
            x[3 * i + 1] = initial->n_cm3[i] / scaling->Ns;
            x[3 * i + 2] = initial->p_cm3[i] / scaling->Ns;
        }
    } else {
        auto equilibrium = solve_equilibrium(
            device,
            {.newton = options.newton, .linear = options.linear,
             .Ns_override = options.Ns_override, .models = options.models},
            control);
        if (!equilibrium) {
            sweep.stopped = std::move(equilibrium.error());
            return sweep;
        }
        std::vector<double> psi(nodes);
        for (std::size_t i = 0; i < nodes; ++i) {
            psi[i] = equilibrium->fields.potential_V[i] / scaling->V_T;
        }
        x = system->state_from_potential(psi);
    }

    const double current_scale =
        scaling->J0 * std::pow(scaling->L_D, device.mesh().dimension() - 1);
    for (std::size_t k = 0; k < points.size(); ++k) {
        if (control.stop.stop_requested()) {
            sweep.stopped = base::Error{
                base::ErrorCode::cancelled,
                "cancelled before bias point " + std::to_string(k),
                base::ErrorContext{.index = k, .value = std::nullopt}};
            return sweep;
        }
        if (auto ok = system->set_bias(points[k]); !ok) {  // validated above; cannot fail
            sweep.stopped = std::move(ok.error());
            return sweep;
        }
        system->stamp_contacts(x);
        results::BiasPoint point;
        const IterationObserver observe = [&](const results::IterationRecord& r, bool converged) {
            if (control.progress) {
                control.progress({Phase::bias, k, points.size(), r.iteration, r.update,
                                  r.residual, converged});
            }
        };
        if (auto ok = newton_solve(*system, x, options.newton, *solver, point.convergence,
                                   control.stop, observe);
            !ok) {
            sweep.stopped = std::move(ok.error());
            sweep.unfinished = std::move(point.convergence);
            return sweep;
        }
        point.bias_V = points[k];
        point.fields = {std::vector<double>(nodes), std::vector<double>(nodes),
                        std::vector<double>(nodes)};
        for (std::size_t i = 0; i < nodes; ++i) {
            point.fields.potential_V[i] = x[3 * i] * scaling->V_T;
            point.fields.n_cm3[i] = x[3 * i + 1] * scaling->Ns;
            point.fields.p_cm3[i] = x[3 * i + 2] * scaling->Ns;
        }
        point.terminal_current = system->terminal_currents(x);
        for (double& I : point.terminal_current) I *= current_scale;
        for (const auto& [jn, jp] : system->edge_currents(x)) {
            point.edge_current_n.push_back(jn * current_scale);
            point.edge_current_p.push_back(jp * current_scale);
        }
        sweep.points.push_back(std::move(point));
    }
    return sweep;
}

std::expected<results::BiasPoint, base::Error> solve_bias(const device::Device& device,
                                                          std::span<const double> bias_V,
                                                          const BiasOptions& options,
                                                          const results::NodeFields* initial,
                                                          const RunControl& control) {
    const std::vector<double> point(bias_V.begin(), bias_V.end());
    auto sweep = sweep_bias(device, std::span(&point, 1), options, initial, control);
    if (!sweep) return std::unexpected(std::move(sweep.error()));
    if (sweep->stopped) return std::unexpected(std::move(*sweep->stopped));
    return std::move(sweep->points.front());
}

}  // namespace NiTCAD::solve
