#include "NiTCAD/deck/run.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <string>
#include <system_error>
#include <utility>

#include "NiTCAD/io/device_file.hpp"
#include "NiTCAD/io/files.hpp"
#include "NiTCAD/linalg/mkl_runtime.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/equilibrium.hpp"

namespace NiTCAD::deck {

namespace {

base::Error invalid(std::string message) {
    return {base::ErrorCode::invalid_input, std::move(message), std::nullopt};
}

// The same bias to 1e-12 relative (1e-12 V near zero): a series written as 0 : 0.3 : 0.1 ends at
// 0.3 exactly, the fourth value of 0 : 0.6 : 0.1 is 3 x 0.1 = 0.30000000000000004 (DK9).
bool same_bias(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (std::abs(a[i] - b[i]) > 1e-12 * std::max({1.0, std::abs(a[i]), std::abs(b[i])})) return false;
    }
    return true;
}

std::string format(double v) {
    char buffer[32];
    const int n = std::snprintf(buffer, sizeof buffer, "%.17g", v);
    return std::string(buffer, static_cast<std::size_t>(n));
}

}  // namespace

std::expected<void, base::Error> write_table(const std::filesystem::path& path,
                                             const device::Device& device,
                                             const results::Sweep& sweep) {
    // With the weighted current (Unit 32) every point has it: two more columns per contact.
    const bool weighted = !sweep.points.empty() && !sweep.points.front().weighted_current.empty();
    std::string text = "point";
    for (const device::Contact& c : device.contacts()) text += "," + c.name + " V";
    for (const device::Contact& c : device.contacts()) text += "," + c.name + " A";
    if (weighted) {
        for (const device::Contact& c : device.contacts()) text += "," + c.name + " weighted A";
        for (const device::Contact& c : device.contacts()) {
            text += "," + c.name + " weighted resolution A";
        }
    }
    text += ",substeps\n";
    for (std::size_t k = 0; k < sweep.points.size(); ++k) {
        const results::BiasPoint& p = sweep.points[k];
        text += std::to_string(k);
        for (const double v : p.bias_V) text += "," + format(v);
        for (const double i : p.terminal_current) text += "," + format(i);
        if (weighted) {
            for (const double i : p.weighted_current) text += "," + format(i);
            for (const double i : p.weighted_current_resolution) text += "," + format(i);
        }
        text += "," + std::to_string(p.bias_substeps) + "\n";
    }
    // Atomic like the binary files: a temporary beside the target, renamed over it.
    std::filesystem::path temporary = path;
    temporary += ".tmp";
    {
        std::ofstream f(temporary, std::ios::binary | std::ios::trunc);
        f.write(text.data(), static_cast<std::streamsize>(text.size()));
        f.flush();
        if (!f) {
            f.close();
            std::error_code ignored;
            std::filesystem::remove(temporary, ignored);
            return std::unexpected(base::Error{base::ErrorCode::resource_exhausted,
                                               "cannot write '" + path.string() + "'", std::nullopt});
        }
    }
    std::error_code ec;
    std::filesystem::rename(temporary, path, ec);
    if (ec) {
        std::error_code ignored;
        std::filesystem::remove(temporary, ignored);
        return std::unexpected(base::Error{base::ErrorCode::resource_exhausted,
                                           "cannot write '" + path.string() + "': " + ec.message(),
                                           std::nullopt});
    }
    return {};
}

std::expected<RunReport, base::Error> run(const Plan& plan, const solve::RunControl& control) {
    RunReport report;
    if (plan.options.linear.backend == linalg::SolverBackend::mkl_pardiso) {
        auto mkl = linalg::load_mkl_runtime(plan.mkl_runtime);
        if (!mkl) return std::unexpected(std::move(mkl.error()));
    }
    const auto written = [&report](const std::filesystem::path& p,
                                   std::expected<void, base::Error> r) -> std::expected<void, base::Error> {
        if (!r) return r;
        report.written.push_back(p);
        return {};
    };
    const Outputs& out = plan.outputs;
    if (!out.device.empty()) {
        if (auto r = written(out.device, io::save_device(out.device, plan.device)); !r) return std::unexpected(r.error());
    }
    if (!out.mesh.empty()) {
        auto r = plan.generated ? io::save_generated_mesh(out.mesh, *plan.generated)
                                : io::save_mesh(out.mesh, plan.device.mesh());
        if (auto w = written(out.mesh, std::move(r)); !w) return std::unexpected(w.error());
    }

    if (plan.analysis == Analysis::equilibrium) {
        const solve::EquilibriumOptions options{.newton = plan.options.newton,
                                                .linear = plan.options.linear,
                                                .Ns_override = plan.options.Ns_override,
                                                .models = plan.options.models};
        auto eq = solve::solve_equilibrium(plan.device, options, control);
        if (!eq) return std::unexpected(std::move(eq.error()));
        if (!out.checkpoint.empty()) {
            const std::vector<std::vector<double>> zero{std::vector<double>(plan.device.contacts().size(), 0.0)};
            results::BiasPoint point;
            point.bias_V = zero.front();
            point.fields = eq->fields;
            const auto record = solve::make_run_record(plan.device, plan.options, zero);
            auto r = io::save_checkpoint(out.checkpoint, plan.device.mesh(), record, point);
            if (auto w = written(out.checkpoint, std::move(r)); !w) return std::unexpected(w.error());
        }
        report.equilibrium = std::move(*eq);
        return report;
    }

    solve::BiasOptions options = plan.options;
    std::vector<std::vector<double>> points = plan.points;
    std::optional<io::Checkpoint> restart;
    if (!plan.restart.empty()) {
        auto c = io::load_checkpoint(plan.restart, plan.device.mesh());
        if (!c) return std::unexpected(std::move(c.error()));
        std::size_t k = 0;
        while (k < points.size() && !same_bias(points[k], c->bias_V)) ++k;
        if (k == points.size()) {
            return std::unexpected(invalid("the restart checkpoint '" + plan.restart.string() +
                                           "' holds a bias that is not a point of the sweep"));
        }
        if (k + 1 == points.size()) {
            return std::unexpected(invalid("the restart checkpoint '" + plan.restart.string() +
                                           "' holds the sweep's last point: nothing is left to solve"));
        }
        report.skipped_points = k + 1;
        points.erase(points.begin(), points.begin() + static_cast<std::ptrdiff_t>(k + 1));
        if (!options.thermal_bias_K.empty()) {
            options.thermal_bias_K.erase(options.thermal_bias_K.begin(),
                                         options.thermal_bias_K.begin() + static_cast<std::ptrdiff_t>(k + 1));
        }
        options.step_control.initial_bias_V = c->bias_V;
        restart = std::move(*c);
    }
    auto sweep = solve::sweep_bias(plan.device, points, options, restart ? &restart->fields : nullptr, control);
    if (!sweep) return std::unexpected(std::move(sweep.error()));
    if (!out.sweep.empty()) {
        if (auto w = written(out.sweep, io::save_sweep(out.sweep, *sweep)); !w) return std::unexpected(w.error());
    }
    if (!out.table.empty()) {
        if (auto w = written(out.table, write_table(out.table, plan.device, *sweep)); !w) return std::unexpected(w.error());
    }
    if (!out.checkpoint.empty() && !sweep->points.empty()) {
        auto r = io::save_checkpoint(out.checkpoint, plan.device.mesh(), sweep->run, sweep->points.back());
        if (auto w = written(out.checkpoint, std::move(r)); !w) return std::unexpected(w.error());
    }
    report.sweep = std::move(*sweep);
    return report;
}

}  // namespace NiTCAD::deck
