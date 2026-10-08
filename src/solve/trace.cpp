#include "NiTCAD/solve/trace.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <span>
#include <stop_token>
#include <string>
#include <utility>

#include "NiTCAD/assemble/contact_bias.hpp"
#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/solve/newton.hpp"
#include "fields.hpp"

namespace NiTCAD::solve {

namespace {

base::Error invalid(std::string message, std::optional<double> value = std::nullopt) {
    return {base::ErrorCode::invalid_input, std::move(message),
            base::ErrorContext{.index = std::nullopt, .value = value}};
}

bool positive(double v) { return std::isfinite(v) && v > 0.0; }

// The bordered system of the trace (trace.hpp) at y = (x, lambda):
//     B = [ J  dF ]        B0 = [ J        dF ]
//         [ c^T   ]             [ c_k e_k^T    ]
// B's last row is dense; factorizing it fills the factors (measured 50 times a factorization of
// J on a 1D junction). B0 keeps only the entry k of that row where the tangent is largest, so it is
// as sparse as J and regular wherever B is (the tangent, J's null direction extended, has a
// nonzero k component); B = B0 + e_n u^T with u = c - c_k e_k, and Sherman-Morrison gives
//     B^-1 r = z - w (u . z) / (1 + u . w),   B0 z = r,   B0 w = e_n.
// The steady rows depend on lambda exactly linearly: F(x, lambda) = F(x, lambda_0) + dF (lambda -
// lambda_0).
class Bordered {
public:
    Bordered(const assemble::DriftDiffusion& dd, std::vector<double> dF, double lambda_ref,
             const linalg::SolverConfig& config)
        : dd_(&dd), dF_(std::move(dF)), lambda_ref_(lambda_ref),
          solver_(*linalg::LinearSolver::create(config)) {  // checked by the caller
        jacobian_ = dd.make_jacobian();
    }

    void set_reference(double lambda_ref) { lambda_ref_ = lambda_ref; }

    // The arc-length row c . (y - y0) = ds, with the sparse row's entry at k.
    void set_row(std::vector<double> c, std::vector<double> y0, double ds, std::size_t k) {
        c_ = std::move(c);
        y0_ = std::move(y0);
        ds_ = ds;
        if (k != k_ || pattern_.rows() == 0) build(k);
    }

    [[nodiscard]] std::size_t unknowns() const { return dd_->unknowns() + 1; }

    // F at y (the steady rows and the arc row) and B0's factorization. Errors: the solver's.
    [[nodiscard]] std::expected<double, base::Error> evaluate(std::span<const double> y,
                                                              std::span<double> f,
                                                              results::ConvergenceRecord& record) {
        const std::size_t n = dd_->unknowns();
        dd_->evaluate(y.first(n), f.first(n), jacobian_);
        const double dl = y[n] - lambda_ref_;
        const std::span<double> v = pattern_.values();
        std::fill(v.begin(), v.end(), 0.0);
        for (std::size_t q = 0; q < to_bordered_.size(); ++q) {
            v[to_bordered_[q]] = jacobian_.values()[q];
        }
        for (std::size_t r = 0; r < n; ++r) {
            if (column_[r] == none) continue;
            f[r] += dF_[r] * dl;
            v[column_[r]] = dF_[r];
        }
        v[v.size() - 1] = c_[k_];  // the last row's single entry
        double arc = -ds_;
        for (std::size_t q = 0; q <= n; ++q) arc += c_[q] * (y[q] - y0_[q]);
        f[n] = arc;
        auto factored = solver_.factorize(pattern_);
        if (!factored) return std::unexpected(std::move(factored.error()));
        if (const auto ratio = factored->pivot_ratio) {
            record.smallest_pivot_ratio =
                std::min(record.smallest_pivot_ratio.value_or(*ratio), *ratio);
        }
        double residual = 0.0;
        for (const double e : f) residual = std::max(residual, std::abs(e));
        return residual;
    }

    // x = B^-1 r with the last factorization. Errors: the solver's.
    [[nodiscard]] std::expected<void, base::Error> solve(std::span<const double> r,
                                                         std::span<double> x) {
        const std::size_t n = dd_->unknowns();
        std::vector<double> e(n + 1, 0.0), w(n + 1);
        e[n] = 1.0;
        if (auto ok = solver_.solve(r, x); !ok) return std::unexpected(std::move(ok.error()));
        if (auto ok = solver_.solve(e, w); !ok) return std::unexpected(std::move(ok.error()));
        double uz = 0.0, uw = 0.0;
        for (std::size_t q = 0; q <= n; ++q) {
            if (q == k_) continue;
            uz += c_[q] * x[q];
            uw += c_[q] * w[q];
        }
        const double s = uz / (1.0 + uw);
        for (std::size_t q = 0; q <= n; ++q) x[q] -= s * w[q];
        return {};
    }

    [[nodiscard]] double update_size(std::span<const double> y,
                                     std::span<const double> dy) const {
        const std::size_t n = dd_->unknowns();
        return std::max(dd_->update_size(y.first(n), dy.first(n)), std::abs(dy[n]));
    }
    void apply_update(std::span<double> y, std::span<const double> dy, double max_update) const {
        const std::size_t n = dd_->unknowns();
        dd_->apply_update(y.first(n), dy.first(n), max_update);
        y[n] += std::clamp(dy[n], -max_update, max_update);
    }

private:
    static constexpr std::size_t none = std::numeric_limits<std::size_t>::max();

    // B0's pattern: J's, the column dF on the contact's rows, the last row's entry at k.
    void build(std::size_t k) {
        k_ = k;
        const linalg::SparseMatrix& j = jacobian_;
        const auto n = static_cast<linalg::Index>(dd_->unknowns());
        std::vector<linalg::Triplet> t;
        for (linalg::Index r = 0; r < n; ++r) {
            for (auto q = j.row_offsets()[r]; q < j.row_offsets()[r + 1]; ++q) {
                t.push_back({r, j.col_indices()[q], 0.0});
            }
            if (dF_[static_cast<std::size_t>(r)] != 0.0) t.push_back({r, n, 0.0});
        }
        t.push_back({n, static_cast<linalg::Index>(k), 0.0});
        pattern_ = *linalg::SparseMatrix::from_triplets(n + 1, n + 1, t);  // in range
        // J's columns come first in every row (dF's column, n, is the last).
        to_bordered_.clear();
        column_.clear();
        for (linalg::Index r = 0; r < n; ++r) {
            const auto start = pattern_.row_offsets()[r];
            for (auto q = j.row_offsets()[r]; q < j.row_offsets()[r + 1]; ++q) {
                to_bordered_.push_back(static_cast<std::size_t>(start + (q - j.row_offsets()[r])));
            }
            column_.push_back(dF_[static_cast<std::size_t>(r)] != 0.0
                                  ? static_cast<std::size_t>(pattern_.row_offsets()[r + 1] - 1)
                                  : none);
        }
    }

    const assemble::DriftDiffusion* dd_;
    std::vector<double> dF_;
    double lambda_ref_;
    linalg::LinearSolver solver_;
    linalg::SparseMatrix jacobian_, pattern_;
    std::vector<std::size_t> to_bordered_, column_;
    std::vector<double> c_, y0_;
    double ds_ = 0.0;
    std::size_t k_ = none;
};

// Newton on the bordered system (newton.hpp's rules: the full correction measured before the
// damping cap, the stop token before every iteration).
std::expected<void, base::Error> correct(Bordered& b, std::span<double> y,
                                         const NewtonOptions& options,
                                         results::ConvergenceRecord& record,
                                         const std::stop_token& stop,
                                         const IterationObserver& observe) {
    record = {};
    const std::size_t n = b.unknowns();
    std::vector<double> f(n), rhs(n), dy(n);
    for (int iteration = 1; iteration <= options.max_iterations; ++iteration) {
        if (stop.stop_requested()) {
            return std::unexpected(base::Error{
                base::ErrorCode::cancelled,
                "cancelled before Newton iteration " + std::to_string(iteration),
                base::ErrorContext{.index = static_cast<std::size_t>(iteration),
                                   .value = std::nullopt}});
        }
        auto residual = b.evaluate(y, f, record);
        if (!residual) return std::unexpected(std::move(residual.error()));
        if (!std::isfinite(*residual)) {
            return std::unexpected(base::Error{
                base::ErrorCode::non_convergence,
                "Newton residual is not finite at iteration " + std::to_string(iteration),
                base::ErrorContext{.index = static_cast<std::size_t>(iteration),
                                   .value = std::nullopt}});
        }
        for (std::size_t q = 0; q < n; ++q) rhs[q] = -f[q];
        if (auto ok = b.solve(rhs, dy); !ok) return std::unexpected(std::move(ok.error()));
        const double update = b.update_size(y, dy);
        b.apply_update(y, dy, options.max_update);
        record.iterations.push_back({iteration, update, *residual});
        const bool converged = update < options.tol_update;
        record.converged = converged;
        if (observe) observe(record.iterations.back(), converged);
        if (converged) return {};
    }
    return std::unexpected(base::Error{
        base::ErrorCode::non_convergence,
        "Newton did not converge in " + std::to_string(options.max_iterations) + " iterations",
        base::ErrorContext{.index = static_cast<std::size_t>(options.max_iterations),
                           .value = record.iterations.back().update}});
}

}  // namespace

std::expected<results::Sweep, base::Error> trace_bias(const device::Device& device,
                                                      const TraceOptions& options,
                                                      const results::NodeFields* initial,
                                                      const RunControl& control) {
    const auto contacts = device.contacts();
    const BiasOptions& steady = options.steady;
    if (steady.equations != Equations::drift_diffusion) {
        return std::unexpected(invalid("a trace solves the drift-diffusion equations"));
    }
    if (options.contact >= contacts.size()) {
        return std::unexpected(invalid("the traced contact does not exist",
                                       static_cast<double>(options.contact)));
    }
    std::vector<device::ContactKind> kinds;
    for (const device::Contact& c : contacts) kinds.push_back(c.kind);
    if (auto ok = assemble::check_contact_bias(kinds, options.start_V, false); !ok) {
        base::Error e = std::move(ok.error());
        e.message = "start bias: " + e.message;
        return std::unexpected(std::move(e));
    }
    const double start = options.start_V[options.contact];
    if (!std::isfinite(options.end_V) || options.end_V == start) {
        return std::unexpected(
            invalid("end_V must be finite and differ from the start", options.end_V));
    }
    for (const auto& [name, v] :
         {std::pair{"step_V", options.step_V}, std::pair{"max_step_V", options.max_step_V},
          std::pair{"min_step_V", options.min_step_V}}) {
        if (!positive(v)) {
            return std::unexpected(invalid(std::string(name) + " must be finite and positive", v));
        }
    }
    if (options.min_step_V > options.step_V || options.step_V > options.max_step_V) {
        return std::unexpected(invalid("steps need min_step_V <= step_V <= max_step_V"));
    }
    if (options.current_limit && !positive(*options.current_limit)) {
        return std::unexpected(
            invalid("current_limit must be finite and positive", *options.current_limit));
    }
    if (options.max_points == 0) return std::unexpected(invalid("max_points must be positive"));
    if (!(std::isfinite(options.density_floor_cm3) && options.density_floor_cm3 >= 0.0)) {
        return std::unexpected(invalid("density_floor must be finite and not negative",
                                       options.density_floor_cm3));
    }
    auto scaling = assemble::make_scaling(device, steady.Ns_override);
    if (!scaling) return std::unexpected(std::move(scaling.error()));
    if (auto solver = linalg::LinearSolver::create(steady.linear); !solver) {
        return std::unexpected(std::move(solver.error()));
    }
    auto system = assemble::DriftDiffusion::create(device, *scaling, steady.models);
    if (!system) return std::unexpected(std::move(system.error()));
    assemble::DriftDiffusion& dd = *system;

    // The starting point (sweep_bias checks the initial state).
    const std::vector<std::vector<double>> first{options.start_V};
    auto sweep = sweep_bias(device, first, steady, initial, control);
    if (!sweep) return std::unexpected(std::move(sweep.error()));
    results::Sweep run;
    run.run = make_run_record(device, options, initial);
    if (sweep->stopped) {
        run.stopped = std::move(sweep->stopped);
        run.unfinished = std::move(sweep->unfinished);
        return run;
    }
    run.points.push_back(std::move(sweep->points.front()));

    const double V_T = scaling->V_T;
    const int D = device.mesh().dimension();
    const std::size_t nodes = device.mesh().node_count();
    const std::size_t n = dd.unknowns();
    const double floor = options.density_floor_cm3 / scaling->Ns;
    const double rms = 1.0 / std::sqrt(static_cast<double>(nodes));
    const double direction = options.end_V > start ? 1.0 : -1.0;
    const double lambda_end = options.end_V / V_T;

    std::vector<double> bias = options.start_V;
    std::vector<double> y(n + 1);
    {
        const results::NodeFields& f = run.points.front().fields;
        for (std::size_t i = 0; i < nodes; ++i) {
            y[3 * i] = f.potential_V[i] / V_T;
            y[3 * i + 1] = f.n_cm3[i] / scaling->Ns;
            y[3 * i + 2] = f.p_cm3[i] / scaling->Ns;
        }
    }
    (void)dd.set_bias(bias);  // checked
    dd.stamp_contacts(std::span(y).first(n));
    y[n] = start / V_T;
    std::vector<double> dF = dd.bias_derivative(options.contact);
    for (double& v : dF) v *= V_T;  // per unit lambda
    Bordered bordered(dd, std::move(dF), y[n], steady.linear);

    // The arc metric, in volts: the swept bias; the swept contact's current at 1 V per decade of
    // |I| + I0, I0 ten times the starting point's current resolution (a current below the
    // resolution does not drive the step); and, weighted by 0.1, the state (the potential in V_T
    // and each density relative to itself above density_floor, as root mean squares over the
    // nodes). Through a snapback the current grows by decades while the bias hardly moves: measured
    // in the state alone (averaged over the nodes) such a step looks like a bias step, and its
    // hyperplane can cut the branch below.
    const double current_scale = scaling->J0 * std::pow(scaling->L_D, D - 1);
    const double I0 =
        std::max(10.0 * run.points.front().terminal_current_resolution[options.contact], 1e-30);
    constexpr double state_weight = 0.1;
    std::vector<double> s(n + 1, 0.0), g(n + 1, 0.0);  // state scales; d log10(|I| + I0) / dy
    const auto metric = [&]() {
        for (std::size_t i = 0; i < nodes; ++i) {
            s[3 * i] = V_T * rms;
            if (device.is_insulator(static_cast<mesh::NodeId>(i))) continue;
            s[3 * i + 1] = V_T * rms / std::max(y[3 * i + 1], floor);
            s[3 * i + 2] = V_T * rms / std::max(y[3 * i + 2], floor);
        }
        s[n] = V_T;
        std::fill(g.begin(), g.end(), 0.0);
        const std::span<const double> x = std::span(y).first(n);
        const double I = dd.terminal_currents(x)[options.contact] * current_scale;
        const auto rows = dd.small_signal_currents(x, 0.0);
        const auto& row = rows[options.contact];
        const double factor =
            (I >= 0.0 ? 1.0 : -1.0) * current_scale / ((std::abs(I) + I0) * std::log(10.0));
        for (std::size_t q = 0; q < row.columns.size(); ++q) {
            g[row.columns[q]] = factor * row.values[q].real();
        }
    };
    const auto dot = [&](const std::vector<double>& a, const std::vector<double>& b) {
        double state = 0.0, ga = 0.0, gb = 0.0;
        for (std::size_t q = 0; q < n; ++q) state += s[q] * a[q] * s[q] * b[q];
        for (std::size_t q = 0; q <= n; ++q) {
            ga += g[q] * a[q];
            gb += g[q] * b[q];
        }
        return state_weight * state_weight * state + s[n] * a[n] * s[n] * b[n] + ga * gb;
    };
    // The arc row of tangent v: c = M v, so that c . v = |v|^2.
    const auto row_of = [&](const std::vector<double>& v) {
        std::vector<double> c(n + 1);
        double gv = 0.0;
        for (std::size_t q = 0; q <= n; ++q) gv += g[q] * v[q];
        for (std::size_t q = 0; q < n; ++q) {
            c[q] = state_weight * state_weight * s[q] * s[q] * v[q] + gv * g[q];
        }
        c[n] = s[n] * s[n] * v[n] + gv * g[n];
        return c;
    };
    // The entry of B0's last row: where the tangent, in the state scales, is largest; the previous
    // one is kept while its component is at least 0.3 of the largest (a new entry changes B0's
    // pattern, and the solver analyzes it afresh).
    std::size_t kept = n;
    const auto pivot = [&](const std::vector<double>& v) {
        std::size_t k = n;
        for (std::size_t q = 0; q < n; ++q) {
            if (std::abs(s[q] * v[q]) > std::abs(s[k] * v[k])) k = q;
        }
        if (std::abs(s[kept] * v[kept]) >= 0.3 * std::abs(s[k] * v[k])) return kept;
        kept = k;
        return k;
    };
    // The unit tangent at y: B v = e_n with the arc row of the previous tangent, normalized in the
    // metric; the first is along the bias.
    std::vector<double> v(n + 1, 0.0);
    metric();
    v[n] = direction / V_T;
    const auto tangent = [&]() -> std::expected<void, base::Error> {
        std::vector<double> f(n + 1), e(n + 1, 0.0), w(n + 1);
        bordered.set_row(row_of(v), y, 0.0, pivot(v));
        results::ConvergenceRecord unused;
        if (auto ok = bordered.evaluate(y, f, unused); !ok) return std::unexpected(ok.error());
        e[n] = 1.0;
        if (auto ok = bordered.solve(e, w); !ok) return std::unexpected(ok.error());
        const double norm = std::sqrt(dot(w, w));
        for (std::size_t q = 0; q <= n; ++q) v[q] = w[q] / norm;
        return {};
    };
    if (auto ok = tangent(); !ok) {
        run.stopped = std::move(ok.error());
        return run;
    }

    double ds = options.step_V;
    std::size_t attempt = 0;
    while (run.points.size() < options.max_points) {
        if (control.stop.stop_requested()) {
            run.stopped = base::Error{base::ErrorCode::cancelled,
                                      "cancelled before trace step " + std::to_string(attempt),
                                      base::ErrorContext{.index = attempt, .value = y[n] * V_T}};
            return run;
        }
        // Predict; land on the end bias if the step would pass it.
        const double dlambda = ds * v[n];
        const bool lands =
            (y[n] + dlambda - lambda_end) * direction >= 0.0 && v[n] * direction > 0.0;
        const double fraction = lands ? (lambda_end - y[n]) / dlambda : 1.0;
        std::vector<double> dy(n + 1), y_next = y;
        for (std::size_t q = 0; q <= n; ++q) dy[q] = fraction * ds * v[q];
        bordered.apply_update(y_next, dy, std::numeric_limits<double>::max());
        if (lands) {
            y_next[n] = lambda_end;
            std::vector<double> c(n + 1, 0.0);
            c[n] = 1.0;
            bordered.set_row(std::move(c), y, lambda_end - y[n], n);
        } else {
            bordered.set_row(row_of(v), y, ds, pivot(v));
        }
        results::BiasPoint point;
        const IterationObserver observe = [&](const results::IterationRecord& r, bool converged) {
            if (control.progress) {
                control.progress({Phase::bias, attempt + 1, 0, r.iteration, r.update, r.residual,
                                  converged});
            }
        };
        auto solved =
            correct(bordered, y_next, steady.newton, point.convergence, control.stop, observe);
        ++attempt;
        if (!solved) {
            if (solved.error().code == base::ErrorCode::cancelled) {
                run.stopped = std::move(solved.error());
                run.unfinished = std::move(point.convergence);
                return run;
            }
            ds *= 0.5;
            if (ds < options.min_step_V) {
                base::Error e = std::move(solved.error());
                e.message = "trace step below min_step_V at " + std::to_string(y[n] * V_T) +
                            " V: " + e.message;
                run.stopped = std::move(e);
                run.unfinished = std::move(point.convergence);
                return run;
            }
            continue;
        }
        // The new tangent may turn by at most acos(0.9) = 26 degrees: a larger turn means the
        // corrector found the branch somewhere else (near a fold two branches lie close, and the
        // hyperplane can cut the other one); the step is then retried at half.
        const std::vector<double> y_old = y, v_old = v, s_old = s, g_old = g;
        y = std::move(y_next);
        bias[options.contact] = y[n] * V_T;
        (void)dd.set_bias(bias);  // finite
        bordered.set_reference(y[n]);
        if (!lands) {
            metric();
            if (auto ok = tangent(); !ok) {
                run.stopped = std::move(ok.error());
                return run;
            }
            const double turn = dot(v, v_old) / std::sqrt(dot(v_old, v_old));
            if (turn < 0.9) {
                y = y_old;
                v = v_old;
                s = s_old;
                g = g_old;
                bias[options.contact] = y[n] * V_T;
                (void)dd.set_bias(bias);
                bordered.set_reference(y[n]);
                ds *= 0.5;
                if (ds < options.min_step_V) {
                    run.stopped = base::Error{
                        base::ErrorCode::non_convergence,
                        "trace step below min_step_V at " + std::to_string(y[n] * V_T) +
                            " V: the tangent turns too fast",
                        base::ErrorContext{.index = attempt, .value = turn}};
                    return run;
                }
                continue;
            }
        }
        // Accepted: the point and the next step.
        point.bias_V = bias;
        detail::drift_diffusion_point(dd, std::span(y).first(n), *scaling, D, point);
        const double current = point.terminal_current[options.contact];
        run.points.push_back(std::move(point));
        if (lands) return run;
        if (options.current_limit && std::abs(current) >= *options.current_limit) return run;
        const std::size_t iterations = run.points.back().convergence.iterations.size();
        if (iterations <= 7) ds = std::min(1.5 * ds, options.max_step_V);
        if (iterations >= 12) ds = std::max(0.7 * ds, options.min_step_V);
    }
    return run;
}

}  // namespace NiTCAD::solve
