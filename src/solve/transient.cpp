#include "NiTCAD/solve/transient.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <limits>
#include <optional>
#include <string>
#include <utility>

#include "NiTCAD/assemble/charges.hpp"
#include "NiTCAD/assemble/contact_bias.hpp"
#include "NiTCAD/assemble/currents.hpp"
#include "NiTCAD/assemble/energy.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/assemble/system.hpp"
#include "NiTCAD/assemble/thermal.hpp"
#include "NiTCAD/assemble/tunnelling.hpp"
#include "NiTCAD/solve/equilibrium.hpp"
#include "NiTCAD/solve/newton.hpp"
#include "density_gradient.hpp"
#include "fields.hpp"
#include "observe.hpp"
#include "wall_clock.hpp"
#include "relocation.hpp"
#include "thermal.hpp"

namespace NiTCAD::solve {

namespace {

using TimeStep = assemble::System::TimeStep;

base::Error invalid(std::string message, std::optional<std::size_t> index = std::nullopt,
                    std::optional<double> value = std::nullopt) {
    return {base::ErrorCode::invalid_input, std::move(message),
            base::ErrorContext{.index = index, .value = value}};
}

// The system of one time step, for newton_solve.
struct StepSystem {
    const assemble::System& dd;
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

// The steady rows with the temperature of the nodes `held` marks held at `hold` (the
// electrothermal start from a given temperature, or of a domain without a heat sink, Unit 23,
// DECISIONS.md U5): each such heat row becomes tau - hold.
struct HeldTemperature {
    const assemble::System& dd;
    std::span<const double> hold;  // tau per node
    std::span<const char> held;    // per node: 1 if held

    [[nodiscard]] std::size_t unknowns() const { return dd.unknowns(); }
    [[nodiscard]] linalg::SparseMatrix make_jacobian() const { return dd.make_jacobian(); }
    void evaluate(std::span<const double> x, std::span<double> f,
                  linalg::SparseMatrix& j) const {
        dd.evaluate(x, f, j);
        const auto offsets = j.row_offsets();
        const auto columns = j.col_indices();
        const std::span<double> values = j.values();
        const std::size_t m = dd.dofs().block(),
                          t = dd.dofs().slot(assemble::FieldId::temperature).value_or(0);
        for (std::size_t i = 0; i < dd.node_count(); ++i) {
            if (held[i] == 0) continue;
            const std::size_t r = m * i + t;
            f[r] = x[r] - hold[i];
            for (auto k = static_cast<std::size_t>(offsets[r]);
                 k < static_cast<std::size_t>(offsets[r + 1]); ++k) {
                values[k] = static_cast<std::size_t>(columns[k]) == r ? 1.0 : 0.0;
            }
        }
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

// A resume state that belongs to this run and system (Unit 61, RD3): its identity, one to three
// states of the system's sizes in increasing time, the last of them the point's time.
std::expected<void, base::Error> check_resume(const results::TransientResume& r, const assemble::System& dd,
                                              std::uint64_t identity, std::size_t contacts) {
    if (r.identity != identity) {
        return std::unexpected(invalid("the resume state belongs to another run (its device, options or waveforms "
                                       "differ)"));
    }
    if (r.history.empty() || r.history.size() > 3) return std::unexpected(invalid("a resume state holds 1 to 3 states"));
    const std::size_t n = dd.unknowns();
    for (std::size_t k = 0; k < r.history.size(); ++k) {
        const auto& s = r.history[k];
        if (s.x.size() != n) return std::unexpected(invalid("a resume state of another system size", k));
        if (s.storage.size() != dd.storage(s.x).size() || s.traps.size() != trap_occupancies(dd, s.x).size() ||
            s.charges.size() != contacts) {
            return std::unexpected(invalid("a resume state's stored quantities do not match the system", k));
        }
        if (!std::isfinite(s.t) || (k > 0 && !(s.t > r.history[k - 1].t))) {
            return std::unexpected(invalid("a resume state's times must increase", k));
        }
    }
    if (!positive(r.time_s) || !positive(r.next_step_s) || r.point.time_s != r.time_s) {
        return std::unexpected(invalid("a resume state needs a positive time, the point's, and a positive next step"));
    }
    return {};
}

std::expected<results::Transient, base::Error> solve_transient_unobserved(
    const device::Device& device, std::span<const Waveform> waveforms,
    const TransientOptions& options, const results::NodeFields* initial,
    const RunControl& control, const results::TransientResume* resume) {
    const auto contacts = device.contacts();
    if (resume != nullptr && initial != nullptr) {
        return std::unexpected(invalid("a continued run takes no initial state: it starts from its resume state"));
    }
    if (waveforms.size() != contacts.size()) {
        return std::unexpected(invalid("a transient run needs one waveform per contact",
                                       std::nullopt, static_cast<double>(waveforms.size())));
    }
    if (options.steady.equations != Equations::drift_diffusion) {
        return std::unexpected(invalid("a transient run solves the drift-diffusion equations"));
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
    // The manufactured sources' profile jumps (verification only, DT-Q3 A) are breakpoints too:
    // steps land on them and restart, so each step sees one weight.
    if (const assemble::ManufacturedSources* ms = options.steady.models.manufactured) {
        for (const double t : ms->profile_times_s) {
            if (t > 0.0 && t <= t_end) restarts.push_back(t);
        }
    }
    // Heavy-ion tracks (Unit 65, HD5): steps land on t0 - 4 sigma, t0 and t0 + 4 sigma and restart
    // there (the source turns on and off on the pulse's scale, not the state's).
    for (const physics::HeavyIonTrack& tr : options.heavy_ions) {
        for (const double t : {tr.t0_s - 4.0 * tr.sigma_s, tr.t0_s, tr.t0_s + 4.0 * tr.sigma_s}) {
            if (t > 0.0 && t <= t_end) restarts.push_back(t);
        }
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
    // Electrothermal (Unit 23): the thermal contacts' temperatures (steady.thermal_bias_K, one
    // list, or the device's) for the whole run.
    const bool thermal = steady.models.electrothermal;
    std::vector<double> thermal_bias;
    if (thermal) {
        auto t = detail::thermal_points(device, steady.thermal_bias_K, 1);
        if (!t) return std::unexpected(std::move(t.error()));
        thermal_bias = std::move(t->front());
        if (initial != nullptr &&
            !detail::valid_temperature(*initial, device.mesh().node_count())) {
            return std::unexpected(invalid(
                "initial temperature needs one finite, positive value per node, or none"));
        }
    }
    auto scaling = assemble::make_scaling(device, steady.Ns_override);
    if (!scaling) return std::unexpected(std::move(scaling.error()));
    auto solver = linalg::LinearSolver::create(steady.linear);
    if (!solver) return std::unexpected(std::move(solver.error()));
    assemble::PhysicsModels models = steady.models;  // with the heavy-ion tracks (Unit 65)
    if (!options.heavy_ions.empty()) models.heavy_ions = &options.heavy_ions;
    auto system = assemble::System::create(device, *scaling, models);
    if (!system) return std::unexpected(std::move(system.error()));
    set_minority_floor(*system, steady.minority_floor);
    if (thermal) (void)set_thermal_bias(*system, thermal_bias);  // checked
    const assemble::System& dd = *system;

    const std::size_t nodes = device.mesh().node_count();
    const detail::Slots slot = detail::slots(dd.dofs());
    const std::size_t m = slot.block;
    const double T0 = scaling->temperature_K;
    const int D = device.mesh().dimension();
    const double current_scale = scaling->J0 * std::pow(scaling->L_D, D - 1);
    const double charge_scale = detail::charge_scale(*scaling, D);
    const double t0 = dd.time_scale();  // seconds per scaled time
    const double density_ref = options.density_ref_cm3 / scaling->Ns;
    // Energy balance (Unit 52): the carrier temperature rises' slots.
    const bool energy = energy_balance(dd);
    const std::size_t rise_n = dd.dofs().slot(assemble::FieldId::electron_temperature).value_or(0);
    const std::size_t rise_p = dd.dofs().slot(assemble::FieldId::hole_temperature).value_or(0);
    const auto lambda_n = dd.dofs().slot(assemble::FieldId::quantum_potential_n);  // Unit 62
    const auto lambda_p = dd.dofs().slot(assemble::FieldId::quantum_potential_p);
    const bool quantum = lambda_n.has_value();

    // What the error estimate measures per node: the carriers of a semiconductor node off the
    // contacts, the potential alone of an insulator node, nothing on a Dirichlet node.
    std::vector<Measure> measure(nodes, Measure::carriers);
    for (std::size_t i = 0; i < nodes; ++i) {
        if (device.is_insulator(static_cast<mesh::NodeId>(i))) measure[i] = Measure::potential;
    }
    for (const device::Contact& c : contacts) {
        // Box rows; a Schottky node's continuity rows too (its psi is held, its error 0).
        if (c.kind == device::ContactKind::gate || c.kind == device::ContactKind::schottky) {
            continue;
        }
        for (const mesh::NodeId i : c.nodes) measure[static_cast<std::size_t>(i)] = Measure::none;
    }
    // Electrothermal: the temperature rise of every node off the isothermal sinks, absolute in
    // units of T0 (as the potential in V_T).
    std::vector<char> heated(nodes, thermal ? 1 : 0);
    if (thermal) {
        for (const device::ThermalContact& c : device.thermal_contacts()) {
            if (c.kind != device::ThermalContactKind::isothermal) continue;
            for (const mesh::NodeId i : device.mesh().find_boundary(c.boundary)->nodes) {
                heated[static_cast<std::size_t>(i)] = 0;
            }
        }
    }

    results::Transient run;
    run.run = make_run_record(device, options, waveforms, initial);
    // Heavy-ion resolution report (Unit 65, HD7): the mesh deposit against LET L (over the area
    // in 1D, the depth in 2D); a warning when they differ by more than 5%.
    const std::vector<double> deposits = dd.heavy_ion_deposits();
    // The tracks' charge generated from 0 to t (HD8); absent without tracks.
    const auto ion_charge = [&](double t_s) -> std::optional<double> {
        if (deposits.empty()) return std::nullopt;
        double q = 0.0;
        for (std::size_t k = 0; k < deposits.size(); ++k) {
            q += deposits[k] * physics::heavy_ion_time_fraction(options.heavy_ions[k], 0.0, t_s);
        }
        return q;
    };
    {
        for (std::size_t k = 0; k < deposits.size(); ++k) {
            const physics::HeavyIonTrack& tr = options.heavy_ions[k];
            const double per = D == 1 ? tr.area_cm2 : D == 2 ? tr.depth_cm : 1.0;
            const double analytic = tr.let_C_per_cm * tr.length_cm / per;
            if (std::abs(deposits[k] / analytic - 1.0) > 0.05) {
                char text[200];
                std::snprintf(text, sizeof text,
                              "heavy-ion track %zu: the mesh deposits %.4g C, LET x length %.4g C (the mesh "
                              "under-resolves the track, or the track leaves the semiconductor)",
                              k, deposits[k], analytic);
                run.run.warnings.emplace_back(text);
            }
        }
    }
    // The identity a continuation is checked against (Unit 61, RD3): the run's without its initial
    // state.
    const std::uint64_t identity =
        initial == nullptr ? run.run.input_identity : make_run_record(device, options, waveforms).input_identity;
    if (resume != nullptr) {
        if (auto ok = check_resume(*resume, dd, identity, contacts.size()); !ok) return std::unexpected(ok.error());
        if (tunnelling(dd)) {
            return std::unexpected(invalid("a transient with nonlocal tunnelling paths is not continued: the paths "
                                           "are not part of the resume state"));
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
            snap.fields.potential_V[i] = s.x[m * i + slot.potential] * scaling->V_T;
            snap.fields.n_cm3[i] = s.x[m * i + slot.electrons] * scaling->Ns;
            snap.fields.p_cm3[i] = s.x[m * i + slot.holes] * scaling->Ns;
        }
        // The lattice (electrothermal) and carrier (Unit 52) temperatures, T0 (1 + tau).
        detail::temperatures(dd, s.x, snap.fields);
        detail::quantum_potentials(dd, s.x, scaling->V_T, snap.fields);  // Unit 62
        snap.bands = detail::band_diagram(dd, s.x, scaling->V_T);
        detail::append_edge_currents(dd, s.x, current_scale, snap.edge_current_n,
                                     snap.edge_current_p);
        snap.trap_occupancy = s.traps;
        if (steady.store_rates) snap.rates = detail::rate_fields(dd, s.x);
        if (!deposits.empty()) {  // Unit 65, HD8: the tracks' rate at this time, as a rate (-G)
            std::vector<double> g = dd.heavy_ion_generation(snap.time_s);
            for (double& v : g) v = v == 0.0 ? 0.0 : -v;
            snap.rates.heavy_ion = std::move(g);
        }
        return snap;
    };

    std::vector<double> x(dd.unknowns());
    std::vector<State> history;  // the accepted states since the last (re)start, newest last
    if (resume != nullptr) {
        // A continuation (Unit 61, RD3): the history, its last step as the first point, and the
        // snapshot of that state.
        for (const auto& s : resume->history) history.push_back({s.t, s.x, s.storage, s.traps, s.charges});
        run.points.push_back(resume->point);
        run.snapshots.push_back(snapshot(history.back()));
    } else {
    // The steady state at t = 0.
    if (initial != nullptr) {
        detail::scaled_state(dd, *initial, *scaling, x);
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
        x = detail::state_from_equilibrium(dd, start->fields, scaling->V_T);  // Lambda too (Unit 62)
        if (auto w = detail::continuation_warning(start->quantum_continuation); !w.empty()) {
            run.run.warnings.push_back(std::move(w));
        }
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
        // Electrothermal (DECISIONS.md U5): a given initial temperature is the temperature at
        // t = 0 everywhere, held while the electrical rows are solved. Otherwise each connected
        // thermal domain with a sink (a region without one of its own included, if joined to one)
        // starts from its steady state, and only a domain with no sink at all (no steady
        // temperature, T7) is held, at T0.
        const bool given = initial != nullptr && !initial->temperature_K.empty();
        std::vector<char> mask(nodes, 0);
        bool hold = false;
        for (std::size_t i = 0; thermal && i < nodes; ++i) {
            mask[i] = given || !reaches_heat_sink(dd, i) ? 1 : 0;
            hold = hold || mask[i] != 0;
        }
        std::vector<double> held(nodes);
        for (std::size_t i = 0; hold && i < nodes; ++i) held[i] = x[m * i + slot.temperature];
        // Energy balance (Unit 52, ED12 as ruled): psi, n and p first at the carrier temperatures
        // of the start (the lattice's from equilibrium, or the given ones), as sweep_bias's
        // points; the first stage's iterations are not recorded.
        if (energy_balance(dd)) {
            hold_temperatures(*system, x);
            results::ConvergenceRecord unused;
            auto held_stage = newton_solve(dd, x, steady.newton, *solver, unused, control.stop);
            release_temperatures(*system);
            if (!held_stage) {
                run.stopped = std::move(held_stage.error());
                run.unfinished = std::move(unused);
                return run;
            }
        }
        if (auto ok = detail::solve_relocating(
                *system, x,
                [&] {
                    if (hold) {
                        return newton_solve(HeldTemperature{dd, held, mask}, x, steady.newton,
                                            *solver, record, control.stop, observe);
                    }
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

    // The starting point.
    history.push_back({0.0, x, dd.storage(x), trap_occupancies(dd, x), contact_charges(dd, x)});
    first.conduction_current = scale(terminal_currents(dd, x), current_scale);
    first.displacement_current.assign(contacts.size(), 0.0);
    first.terminal_current = first.conduction_current;
    first.contact_charge = scale(history.back().charges, charge_scale);
    first.interface_trap_charge = scale(interface_trap_charges(dd, x), charge_scale);
    if (thermal) {
        first.thermal_contact_heat =
            scale(thermal_contact_heat(dd, x), scaling->V_T * current_scale);
    }
    first.heavy_ion_charge = ion_charge(0.0);
    run.points.push_back(std::move(first));
    run.snapshots.push_back(snapshot(history.back()));
    }

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
    double t = resume != nullptr ? resume->time_s : 0.0;        // seconds
    double h = resume != nullptr ? resume->next_step_s : dt0;   // the next step to try, seconds
    std::size_t landing = resume != nullptr ? resume->landing : 0;  // index of the next landing time
    std::size_t attempt = resume != nullptr ? resume->attempt : 0;
    if (resume != nullptr && landing >= landings.size()) {
        return std::unexpected(invalid("the resume state is the run's final step: nothing is left to solve"));
    }
    bool repeating = false;        // the step is being repeated with re-traced tunnel paths
    std::vector<double> f(dd.unknowns());
    // The resume state at the top of the loop: after the last accepted step (Unit 61, RD3). Kept
    // whenever the run returns once a step has been accepted.
    const auto keep = [&] {
        if (run.points.size() < 2 && resume == nullptr) return;  // nothing accepted yet
        results::TransientResume r;
        for (const State& s : history) r.history.push_back({s.t, s.x, s.storage, s.traps, s.charges});
        r.point = run.points.back();
        r.time_s = t;
        r.next_step_s = h;
        r.landing = landing;
        r.attempt = attempt;
        r.identity = identity;
        r.fields = snapshot(history.back()).fields;
        run.resume = std::move(r);
    };
    // State at the top of an attempt, restored before keeping it on a return from inside one.
    double h_top = h;
    std::size_t attempt_top = attempt;
    const auto stop_here = [&] {
        h = h_top;
        attempt = attempt_top;
        keep();
        return run;
    };
    while (landing < landings.size()) {
        h_top = h;
        attempt_top = attempt;
        if (control.stop.stop_requested()) {
            run.stopped = base::Error{base::ErrorCode::cancelled,
                                      "cancelled before time step " + std::to_string(attempt),
                                      base::ErrorContext{.index = attempt, .value = t}};
            return stop_here();
        }
        // The step: land on the next landing time if it is within 1.1 steps.
        const double target = landings[landing];
        // Inside a heavy-ion pulse's window [t0 - 4 sigma, t0 + 4 sigma) a step is at most sigma / 4
        // and lands only when it reaches the landing time (no 1.1 slack; Unit 65, HD5).
        double cap = dt_max;
        for (const physics::HeavyIonTrack& tr : options.heavy_ions) {
            if (t >= tr.t0_s - 4.0 * tr.sigma_s && t < tr.t0_s + 4.0 * tr.sigma_s) {
                cap = std::min(cap, 0.25 * tr.sigma_s);
            }
        }
        const double h_try = std::min(h, cap);
        const bool lands = target - t <= (cap < dt_max ? 1.0 : 1.1) * h_try;
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
        const TimeStep step{1.0 / (formula.beta * step_s / t0), storage, traps, t_new, t};
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
            return stop_here();
        }
        // Nonlocal tunnelling: the paths re-traced at the new state. If they differ, the step is
        // repeated once with them; if they differ again it is accepted with the paths lagging
        // and the next step takes the new ones.
        if (tunnelling(dd)) {
            assemble::TunnelPaths traced = trace_paths(dd, x);
            if (!paths_agree(dd, traced, paths(dd), x)) {
                set_paths(*system, std::move(traced));
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
        State next{t_new / t0, x, dd.storage(x), trap_occupancies(dd, x, &step),
                   contact_charges(dd, x, &step)};

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
                const std::size_t psi = m * i + slot.potential, tau = m * i + slot.temperature;
                if (heated[i] != 0) {
                    const double e = factor * (next.x[tau] - predicted(&State::x, tau));
                    err = std::max(err, std::abs(e));
                }
                if (measure[i] == Measure::none) continue;
                err = std::max(err, std::abs(factor * (next.x[psi] - predicted(&State::x, psi))));
                if (measure[i] != Measure::carriers) continue;
                for (const std::size_t k : {m * i + slot.electrons, m * i + slot.holes}) {
                    const double e = factor * (next.x[k] - predicted(&State::x, k));
                    err = std::max(err, std::abs(e) / (std::abs(next.x[k]) + density_ref));
                }
                // Density gradient (Unit 62, QP11): Lambda_n and Lambda_p, absolute in kT (as psi).
                for (std::size_t k = 0; quantum && k < 2; ++k) {
                    const std::size_t r = m * i + (k == 0 ? *lambda_n : *lambda_p);
                    err = std::max(err, std::abs(factor * (next.x[r] - predicted(&State::x, r))));
                }
                // Energy balance: each carrier's rise, absolute in units of T_L (as the lattice's).
                for (std::size_t k = 0; energy && k < 2; ++k) {
                    const std::size_t r = m * i + (k == 0 ? rise_n : rise_p);
                    err = std::max(err, std::abs(factor * (next.x[r] - predicted(&State::x, r))));
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
                    return stop_here();
                }
                h = step_s * shrink;
                ++run.rejected_steps;
                repeating = false;
                continue;
            }
        }

        // Accepted (an electrothermal step whose temperature leaves a material's range stops the
        // run, T10).
        if (thermal) {
            if (auto e = detail::converged_temperature(device, dd.dofs(), x, T0)) {
                run.stopped = std::move(*e);
                return stop_here();
            }
        }
        results::TimePoint point{t_new, step_s, formula.order, ratio, bias};
        if (thermal) {
            point.thermal_contact_heat =
                scale(thermal_contact_heat(dd, x, &step), scaling->V_T * current_scale);
        }
        point.conduction_current = scale(conduction_currents(dd, x, step), current_scale);
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
        point.heavy_ion_charge = ion_charge(t_new);
        point.interface_trap_charge = scale(interface_trap_charges(dd, x, &step), charge_scale);
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
    keep();
    return run;
}

}  // namespace

// The public run: the body above, its result handed to the test observer (observe.hpp).
std::expected<results::Transient, base::Error> solve_transient(
    const device::Device& device, std::span<const Waveform> waveforms,
    const TransientOptions& options, const results::NodeFields* initial,
    const RunControl& control, const results::TransientResume* resume) {
    const detail::WallClock clock;
    auto r = solve_transient_unobserved(device, waveforms, options, initial, control, resume);
    if (r) {
        if (options.steady.models.quantum_confinement == assemble::QuantumConfinement::density_gradient) {
            detail::CoarseScan scan(device);  // Unit 62, QP15
            for (const results::TransientSnapshot& s : r->snapshots) scan.add(s.fields);
            if (auto w = scan.warning(); !w.empty()) r->run.warnings.push_back(std::move(w));
        }
        r->run.provenance.wall_s = clock.seconds();
        detail::observe_result(*r);
    }
    return r;
}

}  // namespace NiTCAD::solve
