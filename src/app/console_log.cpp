// The console's text (console_log.hpp).
#include "NiTCAD/app/console_log.hpp"

#include <format>

#include "NiTCAD/app/run_monitor.hpp"

namespace NiTCAD::app {

void ConsoleLog::add(std::string line) {
    lines_.push_back(std::move(line));
    if (lines_.size() > max_lines) {
        const std::size_t drop = lines_.size() - keep_lines;
        lines_.erase(lines_.begin(), lines_.begin() + static_cast<std::ptrdiff_t>(drop));
        dropped_ += drop;
    }
}

void ConsoleLog::begin_run() {
    next_iteration_ = next_summary_ = 0;
    header_row_ = static_cast<std::size_t>(-1);
    row_iterations_ = 0;
    following_ = true;
}

void ConsoleLog::clear() {
    dropped_ += lines_.size();
    lines_.clear();
}

void ConsoleLog::header(const PointRow& row, std::size_t index) {
    const auto cells = monitor_cells(row, index);  // number, phase, point, ...
    const std::string_view phase = row.phase == solve::Phase::bias ? "bias point" : std::string_view(cells[1]);
    add(cells[2].empty() ? std::format("  {}", phase) : std::format("  {} {}", phase, cells[2]));
    header_row_ = index;
    row_iterations_ = 0;
}

void ConsoleLog::summarise(const PointRow& row, std::size_t index) {
    if (header_row_ != index) header(row, index);
    add(std::format("    {} {} iteration{}, update {:.3g}, residual {:.3g}, {:.2f} s",
                    row.converged ? "converged in" : "not converged after", row.iterations,
                    row.iterations == 1 ? "" : "s", row.update, row.residual, row.seconds));
}

void ConsoleLog::summarise_before(std::span<const PointRow> rows, std::size_t row) {
    for (; next_summary_ < row && next_summary_ < rows.size(); ++next_summary_) summarise(rows[next_summary_], next_summary_);
}

void ConsoleLog::follow(std::span<const PointRow> rows, std::span<const Iteration> its) {
    if (!following_) return;
    for (; next_iteration_ < its.size(); ++next_iteration_) {
        const Iteration& i = its[next_iteration_];
        if (i.row >= rows.size()) continue;
        summarise_before(rows, i.row);  // the points before this one are over
        if (header_row_ != i.row) header(rows[i.row], i.row);
        ++row_iterations_;
        add(std::format("    Newton {:>3}  update {:.3e}", row_iterations_, i.update));
    }
    // A converged point is over; so is every point before the last one.
    while (next_summary_ < rows.size() && (rows[next_summary_].converged || next_summary_ + 1 < rows.size())) {
        summarise(rows[next_summary_], next_summary_);
        ++next_summary_;
    }
}

void ConsoleLog::finish(std::span<const PointRow> rows, std::span<const Iteration> its, double seconds,
                        std::span<const Diagnostic> diagnostics) {
    if (!following_) return;
    follow(rows, its);
    summarise_before(rows, rows.size());
    if (!rows.empty() || !its.empty()) {
        add(std::format("  {} point{}, {} Newton iteration{}, {:.2f} s", rows.size(), rows.size() == 1 ? "" : "s",
                        its.size(), its.size() == 1 ? "" : "s", seconds));
    }
    for (const Diagnostic& d : diagnostics) {
        add(d.line > 0 ? std::format("  line {}: {}", d.line, d.message) : "  " + d.message);
    }
    following_ = false;
}

}  // namespace NiTCAD::app
