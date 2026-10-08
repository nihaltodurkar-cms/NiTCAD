#include "NiTCAD/solve/bias.hpp"

#include <cmath>
#include <optional>
#include <cstddef>
#include <string>
#include <utility>

#include "NiTCAD/assemble/contact_bias.hpp"
#include "NiTCAD/assemble/equilibrium_poisson.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/solve/equilibrium.hpp"
#include "fields.hpp"

namespace NiTCAD::solve {

namespace {

base::Error invalid(std::string message, std::optional<std::size_t> index = std::nullopt) {
    return {base::ErrorCode::invalid_input, std::move(message),
            base::ErrorContext{.index = index, .value = std::nullopt}};
}

// A usable starting state; potential_only for the quasi-static sweep, which reads nothing else.
// The densities of insulator nodes are not read (they are 0).
bool valid_fields(const results::NodeFields& s, const device::Device& device,
                  bool potential_only) {
    const std::size_t n = device.mesh().node_count();
    if (s.potential_V.size() != n) return false;
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(s.potential_V[i])) return false;
    }
    if (potential_only) return true;
    if (s.n_cm3.size() != n || s.p_cm3.size() != n) return false;
    for (std::size_t i = 0; i < n; ++i) {
        if (device.is_insulator(static_cast<mesh::NodeId>(i))) continue;
        if (!(std::isfinite(s.n_cm3[i]) && s.n_cm3[i] > 0.0)) return false;
        if (!(std::isfinite(s.p_cm3[i]) && s.p_cm3[i] > 0.0)) return false;
    }
    return true;
}

// What the two equation sets share: the bias loop with its stop checks, progress events and
// Newton runs. `solve_point(k, point)` solves point k into `point` (convergence filled) and returns
// newton_solve's result; on success it also fills the point's fields and terminal quantities.
template <class SolvePoint>
void run_points(std::span<const std::vector<double>> points, const RunControl& control,
                results::Sweep& sweep, SolvePoint solve_point) {
    for (std::size_t k = 0; k < points.size(); ++k) {
        if (control.stop.stop_requested()) {
            sweep.stopped = base::Error{
                base::ErrorCode::cancelled,
                "cancelled before bias point " + std::to_string(k),
                base::ErrorContext{.index = k, .value = std::nullopt}};
            return;
        }
        results::BiasPoint point;
        point.bias_V = points[k];
        const IterationObserver observe = [&](const results::IterationRecord& r, bool converged) {
            if (control.progress) {
                control.progress({Phase::bias, k, points.size(), r.iteration, r.update,
                                  r.residual, converged});
            }
        };
        if (auto ok = solve_point(k, point, observe); !ok) {
            sweep.stopped = std::move(ok.error());
            sweep.unfinished = std::move(point.convergence);
            return;
        }
        sweep.points.push_back(std::move(point));
    }
}

}  // namespace

std::expected<results::Sweep, base::Error> sweep_bias(const device::Device& device,
                                                      std::span<const std::vector<double>> points,
                                                      const BiasOptions& options,
                                                      const results::NodeFields* initial,
                                                      const RunControl& control) {
    if (points.empty()) return std::unexpected(invalid("a sweep needs at least one bias point"));
    const auto contacts = device.contacts();
    const bool equilibrium = options.equations == Equations::equilibrium_poisson;
    if (!equilibrium && options.equations != Equations::drift_diffusion) {
        return std::unexpected(invalid("unknown equation set"));
    }
    // Every point is checked by the assemblers' own rule before anything is solved; the error
    // names the point (context index) and the contact (message), with the contact's bias as value.
    std::vector<device::ContactKind> kinds;
    for (const device::Contact& c : contacts) kinds.push_back(c.kind);
    for (std::size_t k = 0; k < points.size(); ++k) {
        auto ok = assemble::check_contact_bias(kinds, points[k], equilibrium);
        if (ok) continue;
        std::string where = "bias point " + std::to_string(k);
        std::optional<double> value;
        if (const auto& at = ok.error().context; at && at->index) {
            where += ", contact '" + contacts[*at->index].name + "'";
            value = at->value;
        }
        return std::unexpected(base::Error{base::ErrorCode::invalid_input,
                                           where + ": " + ok.error().message,
                                           base::ErrorContext{.index = k, .value = value}});
    }
    const std::size_t nodes = device.mesh().node_count();
    if (initial != nullptr && !valid_fields(*initial, device, equilibrium)) {
        return std::unexpected(invalid(
            equilibrium ? "initial state needs a finite potential for every node"
                        : "initial state needs a finite potential for every node and positive "
                          "densities for every semiconductor node"));
    }
    auto scaling = assemble::make_scaling(device, options.Ns_override);
    if (!scaling) return std::unexpected(std::move(scaling.error()));
    auto solver = linalg::LinearSolver::create(options.linear);
    if (!solver) return std::unexpected(std::move(solver.error()));

    const int D = device.mesh().dimension();
    const double charge_scale = detail::charge_scale(*scaling, D);
    const std::size_t edges = device.mesh().edges().size();

    if (equilibrium) {
        auto system = assemble::EquilibriumPoisson::create(device, *scaling, options.models);
        if (!system) return std::unexpected(std::move(system.error()));
        results::Sweep sweep;
        sweep.run = make_run_record(device, options, points, initial);
        std::vector<double> psi = system->charge_neutral_potential();
        if (initial != nullptr) {
            for (std::size_t i = 0; i < nodes; ++i) {
                psi[i] = initial->potential_V[i] / scaling->V_T;
            }
        }
        run_points(points, control, sweep,
                   [&](std::size_t k, results::BiasPoint& point,
                       const IterationObserver& observe) -> std::expected<void, base::Error> {
                       if (auto ok = system->set_bias(points[k]); !ok) return ok;  // checked
                       if (auto ok = newton_solve(*system, psi, options.newton, *solver,
                                                  point.convergence, control.stop, observe);
                           !ok) {
                           return ok;
                       }
                       point.fields = detail::equilibrium_fields(*system, psi, *scaling);
                       point.bands =
                           detail::band_diagram(system->band_edges(psi), scaling->V_T);
                       point.terminal_current.assign(contacts.size(), 0.0);
                       point.terminal_current_resolution.assign(contacts.size(), 0.0);
                       point.gate_charge = system->gate_charges(psi);
                       for (double& Q : point.gate_charge) Q *= charge_scale;
                       point.interface_trap_charge = system->interface_trap_charges(psi);
                       for (double& Q : point.interface_trap_charge) Q *= charge_scale;
                       point.edge_current_n.assign(edges, 0.0);
                       point.edge_current_p.assign(edges, 0.0);
                       return {};
                   });
        return sweep;
    }

    auto system = assemble::DriftDiffusion::create(device, *scaling, options.models);
    if (!system) return std::unexpected(std::move(system.error()));
    // The same system without impact ionization, for the retry of a failed point (below).
    std::optional<assemble::DriftDiffusion> plain;
    std::optional<linalg::LinearSolver> plain_solver;
    if (options.models.impact_ionization) {
        assemble::PhysicsModels models = options.models;
        models.impact_ionization = false;
        plain.emplace(*assemble::DriftDiffusion::create(device, *scaling, models));  // as above
        plain_solver.emplace(*linalg::LinearSolver::create(options.linear));        // as above
    }
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
        system->stamp_contacts(x);  // insulator densities 0
    } else {
        auto start = solve_equilibrium(
            device,
            {.newton = options.newton, .linear = options.linear,
             .Ns_override = options.Ns_override, .models = options.models},
            control);
        if (!start) {
            sweep.stopped = std::move(start.error());
            return sweep;
        }
        std::vector<double> psi(nodes);
        for (std::size_t i = 0; i < nodes; ++i) {
            psi[i] = start->fields.potential_V[i] / scaling->V_T;
        }
        x = system->state_from_potential(psi);
    }

    run_points(points, control, sweep,
               [&](std::size_t k, results::BiasPoint& point,
                   const IterationObserver& observe) -> std::expected<void, base::Error> {
                   if (auto ok = system->set_bias(points[k]); !ok) return ok;  // checked
                   system->stamp_contacts(x);
                   const std::vector<double> stamped = x;
                   auto solved = newton_solve(*system, x, options.newton, *solver,
                                              point.convergence, control.stop, observe);
                   // With impact ionization a point that fails is retried from its stamped start,
                   // first converged without the generation (legacy device.py, its 0.0 stage):
                   // stamping a new bias puts the whole step across the cell beside each contact,
                   // and the generation at that transient field can swamp the iterations. Not
                   // first, so a state on a high-current branch (a snapback) is solved where it
                   // is; the retry's iterations are recorded but not reported.
                   if (!solved && plain && solved.error().code != base::ErrorCode::cancelled) {
                       x = stamped;
                       (void)plain->set_bias(points[k]);  // checked
                       results::ConvergenceRecord unused;
                       if (auto ok = newton_solve(*plain, x, options.newton, *plain_solver, unused,
                                                  control.stop);
                           !ok) {
                           return ok;
                       }
                       solved = newton_solve(*system, x, options.newton, *solver,
                                             point.convergence, control.stop);
                   }
                   if (!solved) return solved;
                   detail::drift_diffusion_point(*system, x, *scaling, D, point);
                   return {};
               });
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
