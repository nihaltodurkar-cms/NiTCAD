#include "NiTCAD/solve/bias.hpp"

#include <algorithm>
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
#include "observe.hpp"
#include "wall_clock.hpp"
#include "relocation.hpp"
#include "weighted_current.hpp"
#include "thermal.hpp"

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

// The step control's options (bias.hpp).
std::optional<base::Error> check_step_control(const StepControl& s,
                                              std::span<const device::ContactKind> kinds,
                                              bool equilibrium) {
    if (s.max_cuts < 0 || s.max_cuts > 52) {
        return invalid("step_control.max_cuts must lie in 0..52");
    }
    if (!(std::isfinite(s.growth) && s.growth >= 1.0)) {
        return invalid("step_control.growth must be finite and at least 1");
    }
    if (s.max_iterations < 1) return invalid("step_control.max_iterations must be at least 1");
    if (!s.initial_bias_V.empty() &&
        !assemble::check_contact_bias(kinds, s.initial_bias_V, equilibrium)) {
        return invalid("step_control.initial_bias_V needs one finite bias per contact (and, in the "
                       "quasi-static sweep, every ohmic contact at 0 V)");
    }
    return std::nullopt;
}

// Whether a failed solve may be retried with a shorter bias step.
bool retryable(const base::Error& e) {
    return e.code == base::ErrorCode::non_convergence ||
           e.code == base::ErrorCode::singular_system ||
           e.code == base::ErrorCode::inaccurate_solve;
}

// a + t (b - a), with t = 1 giving b exactly.
std::vector<double> along(std::span<const double> a, std::span<const double> b, double t) {
    std::vector<double> v(b.begin(), b.end());
    if (t >= 1.0) return v;
    for (std::size_t c = 0; c < v.size(); ++c) v[c] = a[c] + t * (b[c] - a[c]);
    return v;
}

// What the two equation sets share: the bias loop with its stop checks, progress events, Newton
// runs and step cutting (Unit 28).
// - `attempt(bias, thermal, point, observe)` solves `state` in place at the given contact biases
//   (and thermal-contact temperatures, empty without the electrothermal model), filling
//   point.convergence, and returns the Newton result;
// - `finish(k, point)` fills the point's fields and terminal quantities from the converged state
//   (and may still fail, e.g. a temperature outside the models' range).
// `previous` is the bias of the starting state when it is known (it is the start of point 0's
// step), `thermal` the thermal-contact temperatures per point (empty without the model).
// With cutting, the progress events count a point's iterations across its attempts, so they stay
// strictly increasing.
template <class Attempt, class Finish>
void run_points(std::span<const std::vector<double>> points,
                std::span<const std::vector<double>> thermal,
                std::optional<std::vector<double>> previous, const StepControl& steps,
                const RunControl& control, std::vector<double>& state, results::Sweep& sweep,
                Attempt attempt, Finish finish) {
    std::vector<double> previous_thermal = thermal.empty() ? std::vector<double>{} : thermal[0];
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
        const std::vector<double> no_thermal;
        const std::vector<double>& target_thermal = thermal.empty() ? no_thermal : thermal[k];
        int reported = 0;  // iterations reported in earlier attempts at this point
        const IterationObserver observe = [&](const results::IterationRecord& r, bool converged) {
            if (control.progress) {
                control.progress({Phase::bias, k, points.size(), reported + r.iteration,
                                  r.update, r.residual, converged});
            }
        };
        const std::vector<double> start = state;
        auto solved = attempt(points[k], target_thermal, point, observe);
        reported += static_cast<int>(point.convergence.iterations.size());
        if (!solved && steps.max_cuts > 0 && previous && retryable(solved.error())) {
            // From the last converged state, halve the step until a sub-step converges, then let
            // the next grow; t is the fraction of the step reached.
            const double smallest = std::ldexp(1.0, -steps.max_cuts);
            std::vector<double> good = start;
            double t = 0.0;
            double h = 0.5;
            state = good;
            while (true) {
                if (control.stop.stop_requested()) {
                    solved = std::unexpected(base::Error{
                        base::ErrorCode::cancelled,
                        "cancelled while cutting the step to bias point " + std::to_string(k),
                        base::ErrorContext{.index = k, .value = std::nullopt}});
                    break;
                }
                const double next = t + h >= 1.0 ? 1.0 : t + h;
                auto sub = attempt(along(*previous, points[k], next),
                                   thermal.empty() ? no_thermal
                                                   : along(previous_thermal, target_thermal, next),
                                   point, observe);
                reported += static_cast<int>(point.convergence.iterations.size());
                if (sub) {
                    t = next;
                    good = state;
                    if (t >= 1.0) {
                        solved = std::move(sub);
                        break;
                    }
                    ++point.bias_substeps;
                    h = std::min(h * steps.growth, 1.0 - t);
                    continue;
                }
                if (!retryable(sub.error()) || h * 0.5 < smallest) {
                    base::Error e = std::move(sub.error());
                    if (retryable(e)) {
                        e.message = "bias point " + std::to_string(k) + ": no convergence with "
                                    "the step cut to 2^-" + std::to_string(steps.max_cuts) +
                                    " of the full step (" + e.message + ")";
                    }
                    solved = std::unexpected(std::move(e));
                    break;
                }
                state = good;
                h *= 0.5;
            }
        }
        if (solved) solved = finish(k, point);
        if (!solved) {
            sweep.stopped = std::move(solved.error());
            sweep.unfinished = std::move(point.convergence);
            return;
        }
        previous = points[k];
        previous_thermal = target_thermal;
        sweep.points.push_back(std::move(point));
    }
}

std::expected<results::Sweep, base::Error> sweep_bias_unobserved(const device::Device& device,
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
    if (auto e = check_step_control(options.step_control, kinds, equilibrium)) {
        return std::unexpected(std::move(*e));
    }
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
    // Electrothermal (Unit 23): drift-diffusion only (T12), and the thermal contacts' temperatures
    // of every point checked before anything is solved.
    const bool thermal = options.models.electrothermal;
    std::vector<std::vector<double>> thermal_bias;
    if (thermal) {
        if (equilibrium) {
            return std::unexpected(invalid(
                "electrothermal needs the drift-diffusion equations (DECISIONS.md T12)"));
        }
        auto t = detail::thermal_points(device, options.thermal_bias_K, points.size());
        if (!t) return std::unexpected(std::move(t.error()));
        thermal_bias = std::move(*t);
        if (initial != nullptr && !detail::valid_temperature(*initial, nodes)) {
            return std::unexpected(invalid(
                "initial temperature needs one finite, positive value per node, or none"));
        }
    }
    if (initial != nullptr && !valid_fields(*initial, device, equilibrium)) {
        return std::unexpected(invalid(
            equilibrium ? "initial state needs a finite potential for every node"
                        : "initial state needs a finite potential for every node and positive "
                          "densities for every semiconductor node"));
    }
    if (auto e = detail::check_weighted_current(options)) return std::unexpected(std::move(*e));
    auto scaling = assemble::make_scaling(device, options.Ns_override);
    if (!scaling) return std::unexpected(std::move(scaling.error()));
    auto solver = linalg::LinearSolver::create(options.linear);
    if (!solver) return std::unexpected(std::move(solver.error()));
    std::optional<detail::WeightedCurrents> weighted;
    if (options.weighted_current) {
        auto w = detail::WeightedCurrents::create(device, options.linear);
        if (!w) return std::unexpected(std::move(w.error()));
        weighted.emplace(std::move(*w));
    }

    // With step control, every attempt at a bias (the full step included) fails fast: a
    // diverging Newton run is cut after step_control.max_iterations rather than run to the end.
    NewtonOptions newton = options.newton;
    if (options.step_control.max_cuts > 0) {
        newton.max_iterations =
            std::min(newton.max_iterations, options.step_control.max_iterations);
    }
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
        std::optional<std::vector<double>> previous;
        if (initial != nullptr && !options.step_control.initial_bias_V.empty()) {
            previous = options.step_control.initial_bias_V;
        }
        run_points(points, {}, std::move(previous), options.step_control, control, psi, sweep,
                   [&](std::span<const double> bias, std::span<const double>,
                       results::BiasPoint& point,
                       const IterationObserver& observe) -> std::expected<void, base::Error> {
                       if (auto ok = system->set_bias(bias); !ok) return ok;  // checked
                       return newton_solve(*system, psi, newton, *solver,
                                           point.convergence, control.stop, observe);
                   },
                   [&](std::size_t, results::BiasPoint& point) -> std::expected<void, base::Error> {
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
    system->set_minority_floor(options.minority_floor);
    if (thermal && !system->heat_sinks_complete()) {
        return std::unexpected(invalid(
            "electrothermal: a connected part of the device has no thermal contact, so a steady "
            "state has no heat sink (DECISIONS.md T7)"));
    }
    // The same system without generation (impact ionization, tunnelling), for the retry of a
    // failed point (below).
    std::optional<assemble::DriftDiffusion> plain;
    std::optional<linalg::LinearSolver> plain_solver;
    if (options.models.impact_ionization || options.models.btbt_local ||
        options.models.btbt_nonlocal != assemble::NonlocalTunnelling::off) {
        assemble::PhysicsModels models = options.models;
        models.impact_ionization = false;
        models.btbt_local = false;
        models.btbt_nonlocal = assemble::NonlocalTunnelling::off;
        plain.emplace(*assemble::DriftDiffusion::create(device, *scaling, models));  // as above
        plain->set_minority_floor(options.minority_floor);
        plain_solver.emplace(*linalg::LinearSolver::create(options.linear));        // as above
    }
    results::Sweep sweep;
    sweep.run = make_run_record(device, options, points, initial);

    // Starting state, scaled.
    std::vector<double> x(system->unknowns());
    if (initial != nullptr) {
        detail::scaled_state(*system, *initial, *scaling, x);
        if (thermal) (void)system->set_thermal_bias(thermal_bias.front());  // checked
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

    // The bias of the starting state: zero on every contact at equilibrium, or the given one.
    std::optional<std::vector<double>> previous;
    if (initial == nullptr) {
        previous.emplace(contacts.size(), 0.0);
    } else if (!options.step_control.initial_bias_V.empty()) {
        previous = options.step_control.initial_bias_V;
    }
    run_points(points, thermal_bias, std::move(previous), options.step_control, control, x, sweep,
               [&](std::span<const double> bias, std::span<const double> temperatures,
                   results::BiasPoint& point,
                   const IterationObserver& observe) -> std::expected<void, base::Error> {
                   if (auto ok = system->set_bias(bias); !ok) return ok;  // checked
                   if (thermal) {
                       if (auto ok = system->set_thermal_bias(temperatures); !ok) return ok;
                   }
                   system->stamp_contacts(x);
                   const std::vector<double> stamped = x;
                   // With nonlocal tunnelling the paths are re-traced until they agree with the
                   // solution (relocation.hpp); they carry over from the previous point.
                   const auto steady = [&](const IterationObserver& report) {
                       return detail::solve_relocating(*system, x, [&] {
                           return newton_solve(*system, x, newton, *solver,
                                               point.convergence, control.stop, report);
                       }, point.path_relocations);
                   };
                   auto solved = steady(observe);
                   // With generation a point that fails is retried from its stamped start, first
                   // converged without the generation (legacy device.py, its 0.0 stage):
                   // stamping a new bias puts the whole step across the cell beside each contact,
                   // and the generation at that transient field can swamp the iterations. Not
                   // first, so a state on a high-current branch (a snapback) is solved where it
                   // is; the retry's iterations are recorded but not reported.
                   if (!solved && plain && solved.error().code != base::ErrorCode::cancelled) {
                       x = stamped;
                       (void)plain->set_bias(bias);  // checked
                       if (thermal) (void)plain->set_thermal_bias(temperatures);  // likewise
                       results::ConvergenceRecord unused;
                       if (auto ok = newton_solve(*plain, x, newton, *plain_solver, unused,
                                                  control.stop);
                           !ok) {
                           return ok;
                       }
                       if (system->tunnelling()) system->set_paths(system->trace_paths(x));
                       solved = steady({});
                   }
                   return solved;
               },
               [&](std::size_t k, results::BiasPoint& point) -> std::expected<void, base::Error> {
                   if (thermal) {
                       if (auto e = detail::converged_temperature(device, system->dofs(), x,
                                                                  scaling->temperature_K)) {
                           return std::unexpected(std::move(*e));
                       }
                   }
                   detail::drift_diffusion_point(*system, x, *scaling, D, point,
                                                 options.store_rates);
                   if (thermal) point.thermal_bias_K = thermal_bias[k];
                   if (weighted) {
                       return weighted->fill(*system, x, newton.max_update, *scaling, D, point);
                   }
                   return {};
               });
    return sweep;
}

}  // namespace

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

// The public sweep: the body above, its result handed to the test observer (observe.hpp).
std::expected<results::Sweep, base::Error> sweep_bias(const device::Device& device,
                                                      std::span<const std::vector<double>> points,
                                                      const BiasOptions& options,
                                                      const results::NodeFields* initial,
                                                      const RunControl& control) {
    const detail::WallClock clock;
    auto r = sweep_bias_unobserved(device, points, options, initial, control);
    if (r) {
        r->run.provenance.wall_s = clock.seconds();
        detail::observe_result(*r);
    }
    return r;
}

}  // namespace NiTCAD::solve
