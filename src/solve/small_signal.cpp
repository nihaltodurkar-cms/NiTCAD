#include "NiTCAD/solve/small_signal.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <numbers>
#include <string>
#include <utility>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"

namespace NiTCAD::solve {

namespace {

using Complex = std::complex<double>;

base::Error invalid(std::string message, std::optional<std::size_t> index = std::nullopt,
                    std::optional<double> value = std::nullopt) {
    return {base::ErrorCode::invalid_input, std::move(message),
            base::ErrorContext{.index = index, .value = value}};
}

}  // namespace

std::expected<results::SmallSignal, base::Error> solve_small_signal(
    const device::Device& device, std::span<const std::vector<double>> points,
    const SmallSignalOptions& options, const results::NodeFields* initial,
    const RunControl& control) {
    if (options.steady.models.electrothermal) {  // Unit 23, in progress
        return std::unexpected(invalid("electrothermal does not yet support small-signal runs"));
    }
    const std::vector<double>& frequencies = options.frequencies_Hz;
    if (frequencies.empty()) {
        return std::unexpected(invalid("a small-signal run needs at least one frequency"));
    }
    for (std::size_t k = 0; k < frequencies.size(); ++k) {
        if (!(std::isfinite(frequencies[k]) && frequencies[k] >= 0.0)) {
            return std::unexpected(
                invalid("a frequency must be finite and not negative", k, frequencies[k]));
        }
    }
    for (std::size_t k = 0; k < options.field_frequencies_Hz.size(); ++k) {
        const double f = options.field_frequencies_Hz[k];
        if (std::ranges::find(frequencies, f) == frequencies.end()) {
            return std::unexpected(
                invalid("a field frequency is not one of the frequencies", k, f));
        }
    }
    const BiasOptions& steady = options.steady;
    auto scaling = assemble::make_scaling(device, steady.Ns_override);
    if (!scaling) return std::unexpected(std::move(scaling.error()));
    auto solver = linalg::ComplexLinearSolver::create(steady.linear);
    if (!solver) return std::unexpected(std::move(solver.error()));
    auto system = assemble::DriftDiffusion::create(device, *scaling, steady.models);
    if (!system) return std::unexpected(std::move(system.error()));

    // The operating points. sweep_bias checks the points and the initial state.
    auto sweep = sweep_bias(device, points, steady, initial, control);
    if (!sweep) return std::unexpected(std::move(sweep.error()));

    const std::size_t contacts = device.contacts().size();
    const std::size_t nodes = device.mesh().node_count();
    const int D = device.mesh().dimension();
    const double current_scale = scaling->J0 * std::pow(scaling->L_D, D - 1);
    const double t0 = system->time_scale();
    const std::size_t n = system->unknowns();

    results::SmallSignal run;
    run.run = make_run_record(device, options, points, initial);
    run.contacts = contacts;
    run.frequency_Hz = frequencies;

    std::vector<std::vector<Complex>> forcing(contacts, std::vector<Complex>(n));
    for (std::size_t j = 0; j < contacts; ++j) {
        const std::vector<double> d = system->bias_derivative(j);
        for (std::size_t q = 0; q < n; ++q) forcing[j][q] = -d[q];
    }
    linalg::ComplexSparseMatrix a = system->make_small_signal_matrix();
    std::vector<std::vector<Complex>> response(contacts, std::vector<Complex>(n));
    std::vector<Complex> product(n);  // A dx, for the residual

    for (std::size_t p = 0; p < sweep->points.size(); ++p) {
        results::BiasPoint& dc = sweep->points[p];
        (void)system->set_bias(dc.bias_V);  // checked by sweep_bias
        std::vector<double> x(n);
        for (std::size_t i = 0; i < nodes; ++i) {
            x[3 * i] = dc.fields.potential_V[i] / scaling->V_T;
            x[3 * i + 1] = dc.fields.n_cm3[i] / scaling->Ns;
            x[3 * i + 2] = dc.fields.p_cm3[i] / scaling->Ns;
        }
        system->stamp_contacts(x);
        // Nonlocal tunnelling: the paths at the operating point, frozen (their geometry's own
        // dependence on the state is not in the small-signal matrix; ARCHITECTURE.md 6.2).
        if (system->tunnelling()) {
            system->set_paths(system->trace_paths(x));
            a = system->make_small_signal_matrix();
        }

        results::SmallSignalPoint point;
        for (std::size_t k = 0; k < frequencies.size(); ++k) {
            if (control.stop.stop_requested()) {
                run.stopped = base::Error{
                    base::ErrorCode::cancelled,
                    "cancelled at bias point " + std::to_string(p) + ", frequency " +
                        std::to_string(k),
                    base::ErrorContext{.index = p, .value = frequencies[k]}};
                return run;
            }
            const double omega = 2.0 * std::numbers::pi * frequencies[k];
            const Complex s{0.0, omega * t0};
            system->small_signal_matrix(x, s, a);
            const auto factored = solver->factorize(a);
            if (!factored) {
                base::Error e = std::move(factored.error());
                e.message = "small-signal system at bias point " + std::to_string(p) + ", " +
                            std::to_string(frequencies[k]) + " Hz: " + e.message;
                run.stopped = std::move(e);
                return run;
            }
            double backward = 0.0;
            std::vector<double> resolution(contacts, 0.0);
            for (std::size_t j = 0; j < contacts; ++j) {
                const auto solved = solver->solve(forcing[j], response[j]);
                if (!solved) {
                    base::Error e = std::move(solved.error());
                    e.message = "small-signal solve at bias point " + std::to_string(p) + ", " +
                                std::to_string(frequencies[k]) + " Hz, contact '" +
                                device.contacts()[j].name + "': " + e.message;
                    run.stopped = std::move(e);
                    return run;
                }
                backward = std::max(backward, solved->backward_error);
                // The currents of all contacts sum to the residual left in the continuity rows
                // and s times that in the Poisson rows (Gauss): their size bounds the miss.
                a.multiply(response[j], product);
                double left = 0.0;
                for (std::size_t i = 0; i < nodes; ++i) {
                    left += std::abs(forcing[j][3 * i + 1] - product[3 * i + 1]) +
                            std::abs(forcing[j][3 * i + 2] - product[3 * i + 2]) +
                            std::abs(s) * std::abs(forcing[j][3 * i] - product[3 * i]);
                }
                resolution[j] = left;
            }
            const auto rows = system->small_signal_currents(x, s);
            std::vector<Complex> y(contacts * contacts);
            constexpr double eps8 = 8.0 * std::numeric_limits<double>::epsilon();
            for (std::size_t j = 0; j < contacts; ++j) {
                double terms = 0.0;  // the largest sum of term magnitudes over the currents
                for (std::size_t i = 0; i < contacts; ++i) {
                    Complex sum = i == j ? rows[i].bias : Complex{};
                    double magnitude = std::abs(sum);
                    for (std::size_t q = 0; q < rows[i].columns.size(); ++q) {
                        const Complex term = rows[i].values[q] * response[j][rows[i].columns[q]];
                        sum += term;
                        magnitude += std::abs(term);
                    }
                    y[i * contacts + j] = sum * current_scale;
                    terms = std::max(terms, magnitude);
                }
                resolution[j] =
                    (resolution[j] + eps8 * static_cast<double>(contacts) * terms) * current_scale;
            }
            point.resolution.push_back(std::move(resolution));
            point.admittance.push_back(std::move(y));
            point.pivot_ratio.push_back(factored->pivot_ratio.value_or(0.0));
            point.backward_error.push_back(backward);
            if (std::ranges::find(options.field_frequencies_Hz, frequencies[k]) !=
                options.field_frequencies_Hz.end()) {
                for (std::size_t j = 0; j < contacts; ++j) {
                    results::SmallSignalFields fields{frequencies[k], j, {}, {}, {}};
                    for (std::size_t i = 0; i < nodes; ++i) {
                        fields.potential.push_back(response[j][3 * i] * scaling->V_T);
                        fields.n_cm3.push_back(response[j][3 * i + 1] * scaling->Ns);
                        fields.p_cm3.push_back(response[j][3 * i + 2] * scaling->Ns);
                    }
                    point.fields.push_back(std::move(fields));
                }
            }
            if (control.progress) {
                control.progress({Phase::small_signal, p, sweep->points.size(),
                                  static_cast<int>(k + 1), 0.0, backward, true, 0.0,
                                  frequencies[k]});
            }
        }
        point.dc = std::move(dc);
        run.points.push_back(std::move(point));
    }
    run.stopped = std::move(sweep->stopped);
    run.unfinished = std::move(sweep->unfinished);
    return run;
}

}  // namespace NiTCAD::solve
