// The console's text (console_log.hpp).
#include "NiTCAD/app/console_log.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <format>

namespace NiTCAD::app {

namespace {

constexpr std::array<std::string_view, 4> filters{"All", "No Newton detail", "Warnings and errors", "Errors"};

std::string plural(std::size_t n, std::string_view one) { return std::format("{} {}{}", n, one, n == 1 ? "" : "s"); }

// The current unit of a device of dimension d (BiasPoint::terminal_current).
std::string_view current_unit(int d) { return d == 1 ? "A/cm^2" : d == 2 ? "A/cm" : "A"; }

// UTF-8 to UTF-16 code units (an invalid byte becomes U+FFFD).
std::u16string utf16(std::string_view s) {
    std::u16string out;
    for (std::size_t i = 0; i < s.size();) {
        const auto b = static_cast<unsigned char>(s[i]);
        const int n = b < 0x80 ? 1 : (b >> 5) == 6 ? 2 : (b >> 4) == 14 ? 3 : (b >> 3) == 30 ? 4 : 0;
        if (n == 0 || i + n > s.size()) {
            out += static_cast<char16_t>(0xFFFD);
            ++i;
            continue;
        }
        char32_t c = n == 1 ? b : b & (0x7F >> n);
        for (int k = 1; k < n; ++k) c = (c << 6) | (static_cast<unsigned char>(s[i + k]) & 0x3F);
        if (c >= 0x10000) {
            c -= 0x10000;
            out += static_cast<char16_t>(0xD800 + (c >> 10));
            out += static_cast<char16_t>(0xDC00 + (c & 0x3FF));
        } else {
            out += static_cast<char16_t>(c);
        }
        i += static_cast<std::size_t>(n);
    }
    return out;
}

constexpr std::string_view rule = "  ------------------------------------------------------------";

}  // namespace

bool shows(LogFilter filter, LogLevel level) noexcept {
    switch (filter) {
    case LogFilter::all: return true;
    case LogFilter::summary: return level != LogLevel::detail;
    case LogFilter::warnings: return level == LogLevel::warning || level == LogLevel::error;
    case LogFilter::errors: return level == LogLevel::error;
    }
    return true;
}

std::span<const std::string_view> filter_names() noexcept { return filters; }

std::vector<LogLine> plan_banner(const PlanSummary& plan) {
    std::vector<LogLine> out;
    out.push_back({LogLevel::header, std::string(rule)});
    for (const std::string& d : plan.description) out.push_back({LogLevel::info, "   " + d});
    out.push_back({LogLevel::header, std::string(rule)});
    return out;
}

std::vector<LogLine> report_summary(const deck::Plan& plan, const deck::RunReport& report) {
    std::vector<LogLine> out;
    const auto& contacts = plan.device.contacts();
    if (report.sweep && !report.sweep->points.empty()) {
        const auto& pts = report.sweep->points;
        // The swept contacts' voltages (those that change), and the currents of the contacts that carry one.
        std::vector<std::size_t> swept, carrying;
        for (std::size_t c = 0; c < contacts.size(); ++c) {
            const bool varies = std::ranges::any_of(pts, [&](const auto& p) {
                return c < p.bias_V.size() && c < pts.front().bias_V.size() && p.bias_V[c] != pts.front().bias_V[c];
            });
            if (varies) swept.push_back(c);
            if (contacts[c].kind != device::ContactKind::gate) carrying.push_back(c);
        }
        out.push_back({LogLevel::header, std::format("  Terminal currents ({})", current_unit(plan.device.mesh().dimension()))});
        std::string head = "     #";
        for (std::size_t c : swept) head += std::format(" {:>12}", "V(" + contacts[c].name + ")");
        for (std::size_t c : carrying) head += std::format(" {:>12}", "I(" + contacts[c].name + ")");
        out.push_back({LogLevel::info, head});
        const std::size_t shown = std::min(pts.size(), max_table_rows);
        for (std::size_t k = 0; k < shown; ++k) {
            std::string row = std::format("  {:>4}", k + 1);
            for (std::size_t c : swept) row += std::format(" {:>12.4f}", c < pts[k].bias_V.size() ? pts[k].bias_V[c] : NAN);
            for (std::size_t c : carrying) {
                row += std::format(" {:>12.4e}", c < pts[k].terminal_current.size() ? pts[k].terminal_current[c] : NAN);
            }
            out.push_back({LogLevel::info, row});
        }
        if (pts.size() > shown) out.push_back({LogLevel::info, std::format("  ... {} more", pts.size() - shown)});
    }
    if (report.transient && !report.transient->points.empty()) {
        out.push_back({LogLevel::info, std::format("  Transient: {} accepted, t = {:.4g} s",
                                                   plural(report.transient->points.size(), "time point"),
                                                   report.transient->points.back().time_s)});
    }
    if (report.small_signal) {
        out.push_back({LogLevel::info, std::format("  Small signal: {}",
                                                   plural(report.small_signal->points.size(), "frequency"))});
    }
    for (const auto& f : report.written) out.push_back({LogLevel::info, "  Wrote " + f.filename().string()});
    return out;
}

void ConsoleLog::add(LogLevel level, std::string text) {
    lines_.push_back({level, std::move(text)});
    if (lines_.size() > max_lines) {
        const std::size_t drop = lines_.size() - keep_lines;
        lines_.erase(lines_.begin(), lines_.begin() + static_cast<std::ptrdiff_t>(drop));
        dropped_ += drop;
    }
}

void ConsoleLog::begin_run() {
    banner_ = false;
    next_iteration_ = next_summary_ = 0;
    header_row_ = static_cast<std::size_t>(-1);
    following_ = true;
}

void ConsoleLog::clear() {
    dropped_ += lines_.size();
    lines_.clear();
}

std::string ConsoleLog::text(LogFilter filter) const {
    std::string out;
    for (const LogLine& l : lines_) {
        if (!shows(filter, l.level)) continue;
        out += l.text;
        out += '\n';
    }
    return out;
}

void ConsoleLog::header(const RunView& run, std::size_t index) {
    const PointRow& r = run.rows[index];
    std::string text;
    switch (r.phase) {
    case solve::Phase::equilibrium:
        text = "  Equilibrium";
        break;
    case solve::Phase::bias: {
        text = r.point_count > 0 ? std::format("  Bias point {}/{}", r.point + 1, r.point_count)
                                 : std::format("  Bias point {}", r.point + 1);
        // The point's contact voltages, when the plan lists the points the solver steps through.
        if (run.summary && r.point < run.summary->points.size() && r.point_count == run.summary->points.size()) {
            const auto& v = run.summary->points[r.point];
            const auto& contacts = run.summary->contacts;
            std::string biases;
            for (std::size_t c = 0; c < v.size() && c < contacts.size(); ++c) {
                biases += std::format("{}V({}) = {:.4g} V", biases.empty() ? "" : ", ", contacts[c], v[c]);
            }
            if (!biases.empty()) text += "   " + biases;
        }
        break;
    }
    case solve::Phase::transient:
        text = std::format("  Time step {}   t = {:.4g} s", r.point + 1, r.time_s);
        break;
    case solve::Phase::small_signal:
        text = std::format("  Frequency {}/{}   f = {:.4g} Hz", r.point + 1, r.point_count, r.frequency_Hz);
        break;
    }
    add(LogLevel::info, std::move(text));
    add(LogLevel::detail, "     Iter      |update|    |residual|");
    header_row_ = index;
}

void ConsoleLog::summarise(const RunView& run, std::size_t index) {
    const PointRow& r = run.rows[index];
    if (header_row_ != index) header(run, index);
    if (r.converged) {
        add(LogLevel::success, std::format("     Converged in {}: |update| {:.2e}, |residual| {:.2e}, {:.2f} s",
                                           plural(static_cast<std::size_t>(r.iterations), "iteration"), r.update,
                                           r.residual, r.seconds));
    } else {
        add(LogLevel::warning, std::format("     Not converged after {}: |update| {:.2e}, |residual| {:.2e}, {:.2f} s",
                                           plural(static_cast<std::size_t>(r.iterations), "iteration"), r.update,
                                           r.residual, r.seconds));
    }
}

void ConsoleLog::summarise_before(const RunView& run, std::size_t row) {
    for (; next_summary_ < row && next_summary_ < run.rows.size(); ++next_summary_) summarise(run, next_summary_);
}

void ConsoleLog::follow(const RunView& run) {
    if (!following_) return;
    if (!banner_ && run.summary) {
        for (LogLine& l : plan_banner(*run.summary)) add(l.level, std::move(l.text));
        banner_ = true;
    }
    for (; next_iteration_ < run.iterations.size(); ++next_iteration_) {
        const Iteration& i = run.iterations[next_iteration_];
        if (i.row >= run.rows.size()) continue;
        summarise_before(run, i.row);  // the points before this one are over
        if (header_row_ != i.row) header(run, i.row);
        add(LogLevel::detail, std::format("     {:>4}    {:.3e}    {:.3e}", i.in_row, i.update, i.residual));
    }
    // A converged point is over; so is every point before the last one.
    while (next_summary_ < run.rows.size() &&
           (run.rows[next_summary_].converged || next_summary_ + 1 < run.rows.size())) {
        summarise(run, next_summary_);
        ++next_summary_;
    }
}

void ConsoleLog::finish(const RunView& run) {
    if (!following_) return;
    follow(run);
    summarise_before(run, run.rows.size());
    if (!run.rows.empty() || !run.iterations.empty()) {
        add(LogLevel::header, "  Summary");
        const double per = run.rows.empty() ? 0.0
                                            : static_cast<double>(run.iterations.size()) / static_cast<double>(run.rows.size());
        const auto converged = static_cast<std::size_t>(std::ranges::count_if(run.rows, &PointRow::converged));
        add(LogLevel::info, std::format("  {} ({} converged), {}, {:.1f} per point, {:.2f} s", plural(run.rows.size(), "point"),
                                        converged, plural(run.iterations.size(), "Newton iteration"), per, run.seconds));
        if (run.plan && run.report) {
            for (LogLine& l : report_summary(*run.plan, *run.report)) add(l.level, std::move(l.text));
        }
    }
    if (run.stopped) add(LogLevel::warning, "  Stopped early: " + run.stopped->message);
    if (run.error && run.diagnostics.empty()) add(LogLevel::error, "  " + run.error->message);
    for (const Diagnostic& d : run.diagnostics) {
        add(LogLevel::error, d.line > 0 ? std::format("  line {}: {}", d.line, d.message) : "  " + d.message);
    }
    following_ = false;
}

void ConsoleLog::follow(const RunController& r) {
    follow(RunView{.rows = r.points(),
                   .iterations = r.iterations(),
                   .summary = r.plan_summary() ? &*r.plan_summary() : nullptr,
                   .seconds = r.elapsed_seconds(),
                   .diagnostics = r.diagnostics()});
}

void ConsoleLog::finish(const RunController& r) {
    const auto stopped = r.stopped_early();
    finish(RunView{.rows = r.points(),
                   .iterations = r.iterations(),
                   .summary = r.plan_summary() ? &*r.plan_summary() : nullptr,
                   .plan = r.plan() ? &*r.plan() : nullptr,
                   .report = r.report() ? &*r.report() : nullptr,
                   .seconds = r.elapsed_seconds(),
                   .diagnostics = r.diagnostics(),
                   .stopped = stopped ? &*stopped : nullptr,
                   .error = r.state() == RunState::failed && r.error() ? &*r.error() : nullptr});
}

std::string console_rtf(const std::deque<LogLine>& lines, std::size_t from, LogFilter filter, bool continues,
                        int points) {
    // The colour table: 1 detail, 2 info, 3 header, 4 success, 5 warning, 6 error (a dark theme).
    std::string out = "{\\rtf1\\ansi\\deff0{\\fonttbl{\\f0\\fmodern Consolas;}}"
                      "{\\colortbl;\\red140\\green140\\blue140;\\red212\\green212\\blue212;\\red86\\green156\\blue214;"
                      "\\red106\\green200\\blue106;\\red229\\green192\\blue123;\\red244\\green71\\blue71;}";
    out += std::format("\\f0\\fs{} ", 2 * points);
    bool first = !continues;
    for (std::size_t k = from; k < lines.size(); ++k) {
        const LogLine& l = lines[k];
        if (!shows(filter, l.level)) continue;
        if (!first) out += "\\par\n";
        first = false;
        const bool bold = l.level == LogLevel::header || l.level == LogLevel::error;
        out += std::format("\\cf{}{} ", static_cast<int>(l.level) + 1, bold ? "\\b" : "\\b0");
        for (const char16_t u : utf16(l.text)) {
            if (u == u'\\' || u == u'{' || u == u'}') {
                out += '\\';
                out += static_cast<char>(u);
            } else if (u < 0x80) {
                out += static_cast<char>(u);
            } else {
                out += std::format("\\u{}?", static_cast<int>(static_cast<std::int16_t>(u)));
            }
        }
    }
    out += "}";
    return out;
}

}  // namespace NiTCAD::app
