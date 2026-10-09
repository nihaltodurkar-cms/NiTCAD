#include "NiTCAD/solve/transient.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <string>
#include <utility>

#include "NiTCAD/assemble/contact_bias.hpp"
#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/solve/equilibrium.hpp"
#include "NiTCAD/solve/newton.hpp"
#include "fields.hpp"
#include "relocation.hpp"

namespace NiTCAD::solve {

namespace {

using TimeStep = assemble::DriftDiffusion::TimeStep;

base::Error invalid(std::string message, std::optional<std::size_t> index = std::nullopt,
                    std::optional<double> value = std::nullopt) {
    return {base::ErrorCode::invalid_input, std::move(message),
            base::ErrorContext{.index = index, .value = value}};
}

// The system of one time step, for newton_solve.
struct StepSystem {
    const assemble::DriftDiffusion& dd;
    TimeStep step;

    [[nodiscard]] std::size_t unknowns() const { return dd.unknowns(); }
    [[nodiscard]] linalg::SparseMatrix make_jacobian() const { return dd.make_jacobian(); }
    void evaluate(std::span<const double> x, std::span<double> f,
                  linalg::SparseMatrix& j) const {
        dd.evaluate(x, step, f, j);
    }
    [[nodiscard]] double update_size(std::span<const double> x,
                                     std::span<const double> dx) const {
        return dd.update_size(x, dx);
    }
    void apply_update(std::span<double> x, std::span<const double> dx, double max_update) const {
        dd.apply_update(x, dx, max_update);
    }
};

// An accepted state and what later BDF formulas read of it (scaled).
struct State {
    double t;
    std::vector<double> x, storage, traps, charges;
};

// y_new - (a1 y_n - a2 y_n-1) = beta h F(y_new).
struct Formula {
    int order;
    double beta, a1, a2;
};

constexpr Formula backward_euler{1, 1.0, 1.0, 0.0};

// Variable-step BDF2 for a step h after a step h_prev.
Formula bdf2(double h, double h_prev) {
    const double w = h / h_prev, d = 1.0 + 2.0 * w;
    return {2, (1.0 + w) / d, (1.0 + w) * (1.0 + w) / d, w * w / d};
}

std::vector<double> combine(const Formula& f, const std::vector<double>& y_n,
                            const std::vector<double>* y_prev) {
    std::vector<double> c(y_n.size());
    for (std::size_t i = 0; i < c.size(); ++i) {
        c[i] = f.order == 1 ? y_n[i] : f.a1 * y_n[i] - f.a2 * (*y_prev)[i];
    }
    return c;
}

// Per node: what the error estimate measures.
enum class Measure : unsigned char { carriers, potential, none };

bool valid_fields(const results::NodeFields& s, const device::Device& device) {
    const std::size_t n = device.mesh().node_count();
    if (s.potential_V.size() != n || s.n_cm3.size() != n || s.p_cm3.size() != n) return false;
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(s.potential_V[i])) return false;
        if (device.is_insulator(static_cast<mesh::NodeId>(i))) continue;
        if (!(std::isfinite(s.n_cm3[i]) && s.n_cm3[i] > 0.0)) return false;
        if (!(std::isfinite(s.p_cm3[i]) && s.p_cm3[i] > 0.0)) return false;
    }
    return true;
}

bool positive(double v) { return std::isfinite(v) && v > 0.0; }

}  // namespace

std::expected<results::Transient, base::Error> solve_transient(
    const device::Device& device, std::span<const Waveform> waveforms,
    const TransientOptions& options, const results::NodeFields* initial,
    const RunControl& control) {
    const auto contacts = device.contacts();
    if (waveforms.size() != contacts.size()) {
        return std::unexpected(invalid("a transient run needs one waveform per contact",
                                       std::nullopt, static_cast<double>(waveforms.size())));
    }
    if (options.steady.equations != Equations::drift_diffusion) {
        return std::unexpected(invalid("a transient run solves the drift-diffusion equations"));
    }
    if (options.steady.models.electrothermal) {  // Unit 23, in progress
        return std::unexpected(invalid("electrothermal does not yet support transient runs"));
    }
    const double t_end = options.t_end_s;
    if (!positive(t_end)) {
        return std::unexpected(invalid("t_end must be finite and positive", std::nullopt, t_end));
    }
    const double dt0 = options.dt_initial_s.value_or(t_end * 1e-6);
    const double dt_min = options.dt_min_s.value_or(dt0 * 1e-6);
    const double dt_max = options.dt_max_s.value_or(t_end);
    for (const auto& [name, v] : {std::pair{"dt_initial", dt0}, std::pair{"dt_min", dt_min},
                                  std::pair{"dt_max", dt_max}}) {
        if (!positive(v)) {
            return std::unexpected(invalid(std::string(name) + " must be finite and positive",
                                           std::nullopt, v));
        }
    }
    if (dt_min > dt0 || dt0 > dt_max) {
        return std::unexpected(invalid("time steps need dt_min <= dt_initial <= dt_max"));
    }
    if (!positive(options.rtol)) {
        return std::unexpected(
            invalid("rtol must be finite and positive", std::nullopt, options.rtol));
    }
    if (!(std::isfinite(options.density_ref_cm3) && options.density_ref_cm3 >= 0.0)) {
        return std::unexpected(invalid("density_ref must be finite and not negative",
                                       std::nullopt, options.density_ref_cm3));
    }
    for (std::size_t k = 0; k < options.output_times_s.size(); ++k) {
        const double t = options.output_times_s[k];
        if (!(std::isfinite(t) && t > 0.0 && t <= t_end)) {
            return std::unexpected(invalid("output time outside (0, t_end]", k, t));
        }
    }
    // The waveforms at the start and at every breakpoint (both sides), by the assembler's rule.
    std::vector<device::ContactKind> kinds;
    for (const device::Contact& c : contacts) kinds.push_back(c.kind);
    std::vector<double> restarts;  // waveform breakpoints in (0, t_end]
    for (const Waveform& w : waveforms) {
        for (const double t : w.breakpoints(0.0, t_end)) restarts.push_back(t);
    }
    std::ranges::sort(restarts);
    restarts.erase(std::unique(restarts.begin(), restarts.end()), restarts.end());
    const auto bias_at = [&](double t, bool left) {
        std::vector<double> b;
        for (const Waveform& w : waveforms) b.push_back(left ? w.left_value(t) : w.value(t));
        return b;
    };
    {
        std::vector<double> checked{0.0};
        checked.insert(checked.end(), restarts.begin(), restarts.end());
        for (const double t : checked) {
            for (const bool left : {true, false}) {
                auto ok = assemble::check_contact_bias(kinds, bias_at(t, left), false);
                if (ok) continue;
                std::string where = "waveform at t = " + std::to_string(t) + " s";
                std::optional<std::size_t> index;
                std::optional<double> value;
                if (const auto& at = ok.error().context; at && at->index) {
                    where += ", contact '" + contacts[*at->index].name + "'";
                    index = at->index;
                    value = at->value;
                }
                return std::unexpected(base::Error{base::ErrorCode::invalid_input,
                                                   where + ": " + ok.error().message,
                                                   base::ErrorContext{index, value}});
            }
        }
    }
    if (initial != nullptr && !valid_fields(*initial, device)) {
        return std::unexpected(
            invalid("initial state needs a finite potential for every node and positive "
                    "densities for every semiconductor node"));
    }
    const BiasOptions& steady = options.steady;
    auto scaling = assemble::make_scaling(device, steady.Ns_override);
    if (!scaling) return std::unexpected(std::move(scaling.error()));
    auto solver = linalg::LinearSolver::create(steady.linear);
    if (!solver) return std::unexpected(std::move(solver.error()));
    auto system = assemble::DriftDiffusion::create(device, *scaling, steady.models);
    if (!system) return std::unexpected(std::move(system.error()));
    const assemble::DriftDiffusion& dd = *system;

    const std::size_t nodes = device.mesh().node_count();
    const int D = device.mesh().dimension();
    const double current_scale = scaling->J0 * std::pow(scaling->L_D, D - 1);
    const double charge_scale = detail::charge_scale(*scaling, D);
    const double t0 = dd.time_scale();  // seconds per scaled time
    const double density_ref = options.density_ref_cm3 / scaling->Ns;

    // What the error estimate measures per node: the carriers of a semiconductor node off the
    // contacts, the potential alone of an insulator node, nothing on a Dirichlet node.
    std::vector<Measure> measure(nodes, Measure::carriers);
    for (std::size_t i = 0; i < nodes; ++i) {
        if (device.is_insulator(static_cast<mesh::NodeId>(i))) measure[i] = Measure::potential;
    }
    for (const device::Contact& c : contacts) {
        if (c.kind == device::ContactKind::gate) continue;  // box rows
        for (const mesh::NodeId i : c.nodes) measure[static_cast<std::size_t>(i)] = Measure::none;
    }

    results::Transient run;
    run.run = make_run_record(device, options, waveforms, initial);

    // The steady state at t = 0.
    std::vector<double> x(dd.unknowns());
    if (initial != nullptr) {
        for (std::size_t i = 0; i < nodes; ++i) {
            x[3 * i] = initial->potential_V[i] / scaling->V_T;
            x[3 * i + 1] = initial->n_cm3[i] / scaling->Ns;
            x[3 * i + 2] = initial->p_cm3[i] / scaling->Ns;
        }
    } else {
        auto start = solve_equilibrium(
            device,
            {.newton = steady.newton, .linear = steady.linear,
             .Ns_override = steady.Ns_override, .models = steady.models},
            control);
        if (!start) {
            run.stopped = std::move(start.error());
            return run;
        }
        std::vector<double> psi(nodes);
        for (std::size_t i = 0; i < nodes; ++i) {
            psi[i] = start->fields.potential_V[i] / scaling->V_T;
        }
        x = dd.state_from_potential(psi);
    }
    const std::vector<double> bias0 = bias_at(0.0, true);
    (void)system->set_bias(bias0);  // checked
    dd.stamp_contacts(x);
    results::TimePoint first{0.0, 0.0, 0, std::numeric_limits<double>::quiet_NaN(), bias0};
    {
        const IterationObserver observe = [&](const results::IterationRecord& r, bool converged) {
            if (control.progress) {
                control.progress(
                    {Phase::bias, 0, 1, r.iteration, r.update, r.residual, converged, 0.0});
            }
        };
        results::ConvergenceRecord record;
        int relocations = 0;
        if (auto ok = detail::solve_relocating(
                *system, x,
                [&] {
                    return newton_solve(dd, x, steady.newton, *solver, record, control.stop,
                                        observe);
                },
                relocations);
            !ok) {
            run.stopped = std::move(ok.error());
            run.unfinished = std::move(record);
            return run;
        }
    }

    const auto scale = [](std::vector<double> v, double s) {
        for (double& e : v) e *= s;
        return v;
    };
    const auto snapshot = [&](const State& s) {
        results::TransientSnapshot snap;
        snap.time_s = s.t * t0;
        snap.point = run.points.size() - 1;
        snap.fields = {std::vector<double>(nodes), std::vector<double>(nodes),
                       std::vector<double>(nodes)};
        for (std::size_t i = 0; i < nodes; ++i) {
            snap.fields.potential_V[i] = s.x[3 * i] * scaling->V_T;
            snap.fields.n_cm3[i] = s.x[3 * i + 1] * scaling->Ns;
            snap.fields.p_cm3[i] = s.x[3 * i + 2] * scaling->Ns;
        }
        snap.bands = detail::band_diagram(dd.band_edges(s.x), scaling->V_T);
        for (const auto& [jn, jp] : dd.edge_currents(s.x)) {
            snap.edge_current_n.push_back(jn * current_scale);
            snap.edge_current_p.push_back(jp * current_scale);
        }
        snap.trap_occupancy = s.traps;
        return snap;
    };

    // The starting point.
    std::vector<State> history;  // the accepted states since the last (re)start, newest last
    history.push_back({0.0, x, dd.storage(x), dd.trap_occupancies(x), dd.contact_charges(x)});
    first.conduction_current = scale(dd.terminal_currents(x), current_scale);
    first.displacement_current.assign(contacts.size(), 0.0);
    first.terminal_current = first.conduction_current;
    first.contact_charge = scale(history.back().charges, charge_scale);
    first.interface_trap_charge = scale(dd.interface_trap_charges(x), charge_scale);
    run.points.push_back(std::move(first));
    run.snapshots.push_back(snapshot(history.back()));

    // Landing times: breakpoints (restarts), output times and t_end, in scaled time.
    std::vector<double> landings = restarts;
    landings.insert(landings.end(), options.output_times_s.begin(), options.output_times_s.end());
    landings.push_back(t_end);
    std::ranges::sort(landings);
    landings.erase(std::unique(landings.begin(), landings.end()), landings.end());
    const auto is_restart = [&](double t) { return std::ranges::binary_search(restarts, t); };
    const auto is_output = [&](double t) {
        return t == t_end || std::ranges::find(options.output_times_s, t) !=
                                 options.output_times_s.end();
    };

    const bool use_bdf2 = options.integrator == Integrator::bdf2;
    double t = 0.0;                // seconds
    double h = dt0;                // the next step to try, seconds
    std::size_t landing = 0;       // index of the next landing time
    std::size_t attempt = 0;
    bool repeating = false;        // the step is being repeated with re-traced tunnel paths
    std::vector<double> f(dd.unknowns());
    while (landing < landings.size()) {
        if (control.stop.stop_requested()) {
            run.stopped = base::Error{base::ErrorCode::cancelled,
                                      "cancelled before time step " + std::to_string(attempt),
                                      base::ErrorContext{.index = attempt, .value = t}};
            return run;
        }
        // The step: land on the next landing time if it is within 1.1 steps.
        const double target = landings[landing];
        const double h_try = std::min(h, dt_max);
        const bool lands = target - t <= 1.1 * h_try;
        const double t_new = lands ? target : t + h_try;
        const double step_s = t_new - t;
        // Order: backward Euler for the first two steps after a (re)start (the second has an
        // estimate), then the integrator's.
        const std::size_t known = history.size();
        const State& now = history.back();
        const State* prev = known >= 2 ? &history[known - 2] : nullptr;
        const double h_prev = prev != nullptr ? (now.t - prev->t) * t0 : 0.0;
        const Formula formula =
            use_bdf2 && known >= 3 ? bdf2(step_s, h_prev) : backward_euler;
        const std::vector<double> storage =
            combine(formula, now.storage, prev != nullptr ? &prev->storage : nullptr);
        const std::vector<double> traps =
            combine(formula, now.traps, prev != nullptr ? &prev->traps : nullptr);
        const TimeStep step{1.0 / (formula.beta * step_s / t0), storage, traps};
        const std::vector<double> bias = bias_at(t_new, true);
        (void)system->set_bias(bias);  // checked
        x = now.x;
        dd.stamp_contacts(x);
        const IterationObserver observe = [&](const results::IterationRecord& r, bool converged) {
            if (control.progress) {
                control.progress({Phase::transient, attempt, 0, r.iteration, r.update,
                                  r.residual, converged, t_new});
            }
        };
        results::ConvergenceRecord record;
        auto solved = newton_solve(StepSystem{dd, step}, x, steady.newton, *solver, record,
                                   control.stop, observe);
        ++attempt;
        if (!solved) {
            const base::ErrorCode code = solved.error().code;
            if (code == base::ErrorCode::non_convergence ||
                code == base::ErrorCode::singular_system ||
                code == base::ErrorCode::inaccurate_solve) {
                if (step_s / 4.0 >= dt_min) {
                    h = step_s / 4.0;
                    ++run.rejected_steps;
                    repeating = false;
                    continue;
                }
                base::Error e = std::move(solved.error());
                e.message = "time step below dt_min at t = " + std::to_string(t) + " s: " +
                            e.message;
                run.stopped = std::move(e);
            } else {
                run.stopped = std::move(solved.error());
            }
            run.unfinished = std::move(record);
            return run;
        }
        // Nonlocal tunnelling: the paths re-traced at the new state. If they differ, the step is
        // repeated once with them; if they differ again it is accepted with the paths lagging
        // and the next step takes the new ones.
        if (dd.tunnelling()) {
            assemble::TunnelPaths traced = dd.trace_paths(x);
            if (!dd.paths_agree(traced, dd.paths(), x)) {
                system->set_paths(std::move(traced));
                if (!repeating) {
                    repeating = true;
                    ++run.retraced_steps;
                    h = step_s;
                    continue;
                }
                ++run.lagged_steps;
            }
        }
        repeating = false;
        State next{t_new / t0, x, dd.storage(x), dd.trap_occupancies(x, &step),
                   dd.contact_charges(x, &step)};

        // The local error estimate (Milne's device), when the history allows one.
        double ratio = std::numeric_limits<double>::quiet_NaN();
        if (options.adaptive && known >= 2) {
            const double H = h_prev;
            std::vector<double> c(3, 0.0);  // predictor weights of now, prev, prev2
            double factor;
            if (formula.order == 1) {
                // Linear through prev and now: y^P = now + (now - prev) h / H.
                c = {1.0 + step_s / H, -step_s / H, 0.0};
                factor = -step_s / (2.0 * step_s + H);
            } else {
                const State& prev2 = history[known - 3];
                const double H2 = (prev->t - prev2.t) * t0;
                // Lagrange weights at t_new of (t_n, t_n-1, t_n-2), offsets from t_new.
                const double a = step_s, b = step_s + H, d = step_s + H + H2;
                c = {b * d / ((b - a) * (d - a)), a * d / ((a - b) * (d - b)),
                     a * b / ((a - d) * (b - d))};
                const double w = step_s / h_prev;
                const double C =
                    -(1.0 + w) * (1.0 + w) * step_s * step_s * step_s / (6.0 * w * (1.0 + 2.0 * w));
                const double P = a * b * d / 6.0;
                factor = C / (P - C);
            }
            const auto predicted = [&](const auto& field, std::size_t k) {
                double y = c[0] * (now.*field)[k] + c[1] * (prev->*field)[k];
                if (formula.order == 2) y += c[2] * (history[known - 3].*field)[k];
                return y;
            };
            double err = 0.0;
            for (std::size_t i = 0; i < nodes; ++i) {
                if (measure[i] == Measure::none) continue;
                err = std::max(err, std::abs(factor * (next.x[3 * i] -
                                                       predicted(&State::x, 3 * i))));
                if (measure[i] != Measure::carriers) continue;
                for (const std::size_t k : {3 * i + 1, 3 * i + 2}) {
                    const double e = factor * (next.x[k] - predicted(&State::x, k));
                    err = std::max(err, std::abs(e) / (std::abs(next.x[k]) + density_ref));
                }
            }
            for (std::size_t k = 0; k < next.traps.size(); ++k) {
                err = std::max(err, std::abs(factor * (next.traps[k] -
                                                       predicted(&State::traps, k))));
            }
            ratio = err / options.rtol;
            if (ratio > 1.0) {
                const double shrink =
                    std::max(0.2, 0.9 * std::pow(ratio, -1.0 / (formula.order + 1)));
                if (step_s * shrink < dt_min) {
                    run.stopped = base::Error{
                        base::ErrorCode::non_convergence,
                        "time step below dt_min at t = " + std::to_string(t) +
                            " s: the error estimate is not met",
                        base::ErrorContext{.index = attempt - 1, .value = ratio}};
                    run.unfinished = std::move(record);
                    return run;
                }
                h = step_s * shrink;
                ++run.rejected_steps;
                repeating = false;
                continue;
            }
        }

        // Accepted.
        results::TimePoint point{t_new, step_s, formula.order, ratio, bias};
        point.conduction_current = scale(dd.conduction_currents(x, step), current_scale);
        point.displacement_current.resize(contacts.size());
        const std::vector<double> q_old =
            combine(formula, now.charges, prev != nullptr ? &prev->charges : nullptr);
        for (std::size_t c = 0; c < contacts.size(); ++c) {
            point.displacement_current[c] =
                (next.charges[c] - q_old[c]) * step.rate * current_scale;
        }
        point.terminal_current.resize(contacts.size());
        for (std::size_t c = 0; c < contacts.size(); ++c) {
            point.terminal_current[c] =
                point.conduction_current[c] + point.displacement_current[c];
        }
        point.contact_charge = scale(next.charges, charge_scale);
        point.interface_trap_charge = scale(dd.interface_trap_charges(x, &step), charge_scale);
        point.convergence = std::move(record);
        run.points.push_back(std::move(point));

        // The next step. A step shortened to land is not the one the control asked for: the next
        // is at least the asked one (as far as the estimate allows).
        if (options.adaptive) {
            const double grow = std::isnan(ratio)
                                    ? 1.0
                                    : std::clamp(0.9 * std::pow(std::max(ratio, 1e-10),
                                                                -1.0 / (formula.order + 1)),
                                                 0.2, 2.0);
            h = step_s * grow;
            if (step_s < h_try) h = std::max(h, h_try * std::min(grow, 1.0));
            h = std::min(h, dt_max);
        } else {
            h = dt0;
        }
        t = t_new;
        const bool restart = lands && is_restart(t_new);
        if (lands) ++landing;
        // After a breakpoint the solution changes fast again (a jump): the first step, which has
        // no estimate, is taken at most at dt_initial.
        if (restart) {
            history.clear();
            h = std::min(h, dt0);
        }
        history.push_back(std::move(next));
        if (history.size() > 3) history.erase(history.begin());
        if (options.fields_every_step || (lands && is_output(t_new))) {
            run.snapshots.push_back(snapshot(history.back()));
        }
    }
    return run;
}

}  // namespace NiTCAD::solve
