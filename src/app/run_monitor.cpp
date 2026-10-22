// The run monitor's content (run_monitor.hpp).
#include "NiTCAD/app/run_monitor.hpp"

#include <algorithm>
#include <array>
#include <format>

#include "NiTCAD/analysis/curve.hpp"

namespace NiTCAD::app {

namespace {

constexpr std::array<std::string_view, 8> columns{"#", "Phase", "Point", "Iterations", "Update", "Residual",
                                                  "Converged", "Seconds"};

}  // namespace

std::span<const std::string_view> monitor_columns() noexcept { return columns; }

std::vector<std::string> monitor_cells(const PointRow& r, std::size_t index) {
    std::string phase, point;
    switch (r.phase) {
    case solve::Phase::equilibrium:
        phase = "equilibrium";
        break;
    case solve::Phase::bias:
        phase = "bias";
        point = r.point_count > 0 ? std::format("{} of {}", r.point + 1, r.point_count) : std::format("{}", r.point + 1);
        break;
    case solve::Phase::transient:
        phase = "transient";
        point = std::format("t = {:.4g} s", r.time_s);
        break;
    case solve::Phase::small_signal:
        phase = "small signal";
        point = std::format("{} of {}, {:.4g} Hz", r.point + 1, r.point_count, r.frequency_Hz);
        break;
    }
    return {std::to_string(index + 1),
            phase,
            point,
            std::to_string(r.iterations),
            std::format("{:.3g}", r.update),
            std::format("{:.3g}", r.residual),
            r.converged ? "yes" : "",
            std::format("{:.2f}", r.seconds)};
}

std::expected<render::PlotScene, base::Error> convergence_plot(std::span<const Iteration> iterations) {
    std::vector<double> x, y;
    for (const Iteration& i : iterations) {
        if (!(i.update > 0.0)) continue;
        x.push_back(static_cast<double>(i.number));
        y.push_back(i.update);
    }
    if (x.size() < 2) {
        return std::unexpected(base::Error{base::ErrorCode::invalid_input, "fewer than two iterations to plot", std::nullopt});
    }
    auto c = analysis::make_curve(std::move(x), std::move(y));
    if (!c) return std::unexpected(std::move(c.error()));
    const std::vector<render::SeriesInput> series{{"largest scaled update", std::move(*c)}};
    return render::curve_plot(series, "Newton iteration", "update", render::AxisScale::linear, render::AxisScale::log10,
                              "Convergence");
}

double progress_fraction(std::span<const PointRow> rows, std::size_t planned) noexcept {
    if (planned == 0) return -1.0;
    std::size_t done = 0;
    for (const PointRow& r : rows) {
        if (r.phase == solve::Phase::bias && r.converged) ++done;  // a small-signal run's operating points too
    }
    return std::min(1.0, static_cast<double>(done) / static_cast<double>(planned));
}

}  // namespace NiTCAD::app
