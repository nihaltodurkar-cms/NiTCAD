// A dimension sweep's results (Unit 26n, DT7): app/sweep_results.hpp.
#include "NiTCAD/app/sweep_results.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <format>
#include <limits>

#include "NiTCAD/analysis/curve.hpp"
#include "NiTCAD/app/dimensions.hpp"
#include "NiTCAD/deck/dimension_sweep.hpp"
#include "NiTCAD/deck/syntax.hpp"
#include "NiTCAD/render/color.hpp"

namespace NiTCAD::app {

namespace {

namespace fs = std::filesystem;

base::Error invalid(std::string what) { return {base::ErrorCode::invalid_input, std::move(what), std::nullopt}; }

// A sweep value's number and unit ("20nm" -> 20, "nm"; "1e18" -> 1e18, ""); the value sweep_values accepted.
std::pair<double, std::string_view> split(std::string_view v) {
    double x = 0.0;
    const char* first = v.data() + (v.starts_with('+') ? 1 : 0);
    const auto [end, ec] = std::from_chars(first, v.data() + v.size(), x);
    return {ec == std::errc{} ? x : 0.0, v.substr(static_cast<std::size_t>(end - v.data()))};
}

double scale_of(std::string_view unit) { return unit == "nm" ? 1e-7 : unit == "um" ? 1e-4 : 1.0; }

std::string current_unit(int dim) { return dim == 1 ? "A/cm\xC2\xB2" : dim == 2 ? "A/cm" : "A"; }

}  // namespace

std::expected<SweepResults, std::string> load_sweep(std::string_view text, const fs::path& path, std::string_view name,
                                                    const std::vector<std::string>& values) {
    const auto d = deck::parse(text, path.filename().string(), path.parent_path());
    if (!d) return std::unexpected(d.error().message);
    const auto params = deck::parameters(*d);
    if (!params) return std::unexpected(params.error().message);
    const auto p = std::ranges::find(*params, name, &deck::Parameter::name);
    if (p == params->end()) return std::unexpected(std::format("no parameter '{}' in [parameters]", name));
    SweepResults s;
    s.deck = path;
    s.name = std::string(name);
    s.length = p->length;
    s.directory = deck::sweep_directory(path, name);
    // The files the deck's [output] names, relative to its folder.
    deck::Outputs named;
    for (const deck::Section& sec : d->sections) {
        if (sec.kind != "output") continue;
        for (const deck::Entry& e : sec.entries) {
            const fs::path f = path.parent_path() / deck::unquote(e.value);
            if (e.key == "sweep") named.sweep = f;
            else if (e.key == "transient") named.transient = f;
            else if (e.key == "small_signal") named.small_signal = f;
            else if (e.key == "device") named.device = f;
            else if (e.key == "mesh") named.mesh = f;
        }
    }
    for (const std::string& v : values) {
        if (const auto ok = deck::sweep_values(v, s.length); !ok) return std::unexpected(std::format("{}: {}", name, ok.error()));
        const auto [number, unit] = split(v);
        SweepRun r{.value = v, .numeric = number * scale_of(unit), .folder = deck::sweep_point_directory(s.directory, name, v),
                   .workspace = std::nullopt, .error = {}};
        const deck::Outputs o = deck::with_outputs_in(named, path.parent_path(), r.folder);
        const fs::path result = !o.sweep.empty() ? o.sweep : !o.small_signal.empty() ? o.small_signal : o.transient;
        std::error_code ec;
        if (result.empty() || o.device.empty()) {
            r.error = "the deck's [output] names no result and device files";
        } else if (!fs::exists(result, ec) || !fs::exists(o.device, ec)) {
            r.error = "not run";
        } else {
            const fs::path mesh = !o.mesh.empty() && fs::exists(o.mesh, ec) ? o.mesh : fs::path();
            auto w = workspace_from_result(result, o.device, mesh);
            if (w) r.workspace = std::move(*w);
            else r.error = w.error().message;
        }
        s.runs.push_back(std::move(r));
    }
    return s;
}

std::string run_label(const SweepResults& s, const SweepRun& r) {
    if (!s.length) return s.name + " = " + r.value;
    const auto [number, unit] = split(r.value);
    const std::string_view shown_unit = unit == "um" ? std::string_view("\xC2\xB5m") : unit;
    return s.name + " = " + std::string(r.value.substr(0, r.value.size() - unit.size())) + " " + std::string(shown_unit);
}

// From the text, so 1.5um shows as 1500 nm exactly (1.5e-4 cm times 1e7 is not).
double shown_value(const SweepResults& s, const SweepRun& r) noexcept {
    if (!s.length) return r.numeric;
    const auto [number, unit] = split(r.value);
    return number * (unit == "nm" ? 1.0 : unit == "um" ? 1e3 : 1e7);
}

std::string value_title(const SweepResults& s) { return s.length ? s.name + " (nm)" : s.name; }

std::expected<render::PlotScene, base::Error> sweep_overlay(const SweepResults& s, bool log_y) {
    std::vector<render::SeriesInput> input;
    std::string x_label, y_label;
    for (const SweepRun& r : s.runs) {
        if (!r.workspace) continue;
        const Workspace& w = *r.workspace;
        const auto fc = figure_curve(w);
        if (!fc) continue;
        auto c = analysis::current_curve(*w.sweep, fc->swept, fc->sense);
        if (!c) continue;
        const auto contacts = w.device->contacts();
        if (x_label.empty()) {
            x_label = "V(" + contacts[fc->swept].name + ") (V)";
            y_label = "I(" + contacts[fc->sense].name + ") (" + current_unit(w.device->mesh().dimension()) + ")";
        }
        input.push_back({run_label(s, r), std::move(*c)});
    }
    if (input.empty()) return std::unexpected(invalid("no run of the sweep of " + s.name + " has a curve to show"));
    auto plot = render::curve_plot(input, std::move(x_label), std::move(y_label), render::AxisScale::linear,
                                   log_y ? render::AxisScale::log10 : render::AxisScale::linear, s.name + " sweep");
    if (!plot) return plot;
    // Viridis in sweep order, the first darkest; its last yellow tenth left out (too pale on the plot's white).
    const std::size_t n = plot->series.size();
    for (std::size_t k = 0; k < n; ++k) {
        const double t = n == 1 ? 0.0 : 0.9 * static_cast<double>(k) / static_cast<double>(n - 1);
        plot->series[k].color = render::sample(render::Colormap::viridis, t);
    }
    plot->legend = true;
    return plot;
}

SweepFigures sweep_figures(const SweepResults& s) {
    SweepFigures f;
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (std::size_t i = 0; i < s.runs.size(); ++i) {
        if (!s.runs[i].workspace) continue;
        const auto figures = extracted_figures(*s.runs[i].workspace);
        std::vector<double> row(f.names.size(), nan);
        for (const Figure& g : figures) {
            auto it = std::ranges::find(f.names, g.name);
            if (it == f.names.end()) {
                f.names.push_back(g.name);
                f.units.push_back(g.unit);
                for (auto& other : f.rows) other.push_back(nan);
                row.push_back(nan);
                it = f.names.end() - 1;
            }
            row[static_cast<std::size_t>(it - f.names.begin())] = g.value;
        }
        f.runs.push_back(i);
        f.rows.push_back(std::move(row));
    }
    return f;
}

std::string figures_csv(const SweepResults& s, const SweepFigures& f) {
    std::string out = value_title(s);
    for (std::size_t j = 0; j < f.names.size(); ++j) {
        out += ',';
        out += f.units[j].empty() ? f.names[j] : f.names[j] + " (" + f.units[j] + ")";
    }
    out += '\n';
    for (std::size_t i = 0; i < f.rows.size(); ++i) {
        out += std::format("{}", shown_value(s, s.runs[f.runs[i]]));
        for (const double v : f.rows[i]) {
            out += ',';
            if (std::isfinite(v)) out += std::format("{}", v);
        }
        out += '\n';
    }
    return out;
}

fs::path figures_csv_path(const SweepResults& s) { return s.directory / "figures.csv"; }

std::expected<render::PlotScene, base::Error> figure_plot(const SweepResults& s, const SweepFigures& f, std::size_t j) {
    if (j >= f.names.size()) return std::unexpected(invalid("no such figure"));
    analysis::Curve c;
    for (std::size_t i = 0; i < f.rows.size(); ++i) {
        if (!std::isfinite(f.rows[i][j])) continue;
        c.x.push_back(shown_value(s, s.runs[f.runs[i]]));
        c.y.push_back(f.rows[i][j]);
    }
    if (c.x.empty()) return std::unexpected(invalid("no run has the figure " + f.names[j]));
    const std::string y = f.units[j].empty() ? f.names[j] : f.names[j] + " (" + f.units[j] + ")";
    const render::SeriesInput in[] = {{f.names[j], std::move(c)}};
    auto plot = render::curve_plot(in, value_title(s), y, render::AxisScale::linear, render::AxisScale::linear,
                                   f.names[j] + " vs " + s.name);
    if (!plot) return plot;
    plot->series[0].markers = true;
    plot->series[0].color = render::sample(render::Colormap::viridis, 0.0);
    return plot;
}

}  // namespace NiTCAD::app

namespace NiTCAD::app {

std::expected<std::vector<std::string>, std::string> typed_sweep_values(std::string_view deck_text, std::string_view name,
                                                                        std::string_view typed) {
    const auto rows = dims::read_dimensions(deck_text);
    if (!rows) return std::unexpected(rows.error());
    const auto row = std::ranges::find(*rows, name, &dims::Dimension::name);
    if (row == rows->end()) return std::unexpected(std::format("no dimension '{}'", name));
    if (row->expr.find('$') != std::string::npos) {
        return std::unexpected(std::format("{} follows other dimensions ({}): sweep those", name, row->expr));
    }
    std::string list(typed);
    if (row->length) {
        const std::string_view expr = row->expr;
        const std::string unit = expr.ends_with("nm") || expr.ends_with("um") || expr.ends_with("cm")
                                     ? std::string(expr.substr(expr.size() - 2))
                                     : std::string(std::abs(row->value) < 1e-4 ? "nm" : "um");
        const char sep = list.find(':') != std::string::npos ? ':' : ',';
        std::string out;
        std::size_t at = 0;
        while (true) {
            const std::size_t end = list.find(sep, at);
            std::string item = list.substr(at, end == std::string::npos ? std::string::npos : end - at);
            const auto b = item.find_first_not_of(" \t"), e = item.find_last_not_of(" \t");
            item = b == std::string::npos ? std::string() : item.substr(b, e - b + 1);
            double v = 0.0;
            const auto [p, ec] = std::from_chars(item.data(), item.data() + item.size(), v);
            if (!item.empty() && ec == std::errc{} && p == item.data() + item.size()) item += unit;
            if (!out.empty() || at > 0) out += sep;
            out += item;
            if (end == std::string::npos) break;
            at = end + 1;
        }
        list = std::move(out);
    }
    auto values = deck::sweep_values(list, row->length);
    if (!values) return std::unexpected(std::format("{}: {}", name, values.error()));
    return values;
}

}  // namespace NiTCAD::app
