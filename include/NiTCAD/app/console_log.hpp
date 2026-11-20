// The console under the inspector (Unit 26j): a running, levelled log of what the window and the solver
// did, as a device simulator's log reads —
//   a banner per run (analysis, equations, device, linear solver, Newton tolerances, points planned);
//   per solver point a header (the bias point with its contact voltages, the transient time, the
//   small-signal frequency) and a Newton table, one row per iteration: |update| and |residual|;
//   each point's outcome (converged in n iterations, or not), its last |update|, |residual| and seconds;
//   at the end a summary: totals, iterations per point, the terminal I-V table, the files written, why
//   the run stopped early; a failed check's errors by line; and the window's status messages.
// The text lives here (UI thread); the window shows it, coloured by level and filtered, in a rich edit.
//
// follow() turns the RunController's state into lines incrementally: the banner once the plan is built,
// each iteration once, each point summarised once when it converges, when a later point starts, or when
// the run ends (finish()). begin_run() starts following a new check or run (the controller clears its rows).
#pragma once
#include <cstddef>
#include <deque>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/app/run_controller.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"

namespace NiTCAD::app {

// detail: a Newton iteration row; header: a banner or a section title.
enum class LogLevel : unsigned char { detail, info, header, success, warning, error };

struct LogLine {
    LogLevel level = LogLevel::info;
    std::string text;
    bool operator==(const LogLine&) const = default;
};

// What the console shows: everything, everything but the Newton rows, warnings and errors, errors.
enum class LogFilter : unsigned char { all, summary, warnings, errors };
[[nodiscard]] bool shows(LogFilter filter, LogLevel level) noexcept;
[[nodiscard]] std::span<const std::string_view> filter_names() noexcept;  // in LogFilter order

// The state of a check or run the console reads (a RunController's, or a test's).
struct RunView {
    std::span<const PointRow> rows;
    std::span<const Iteration> iterations;
    const PlanSummary* summary = nullptr;      // once built
    const deck::Plan* plan = nullptr;          // once done
    const deck::RunReport* report = nullptr;   // once done
    double seconds = 0.0;
    std::span<const Diagnostic> diagnostics;
    const base::Error* stopped = nullptr;      // why a done run stopped early
    const base::Error* error = nullptr;        // why it failed
};

// The banner of a built plan, and the summary lines of a finished run's report (the I-V table at most
// max_table_rows rows, then how many more).
inline constexpr std::size_t max_table_rows = 200;
[[nodiscard]] std::vector<LogLine> plan_banner(const PlanSummary& plan);
[[nodiscard]] std::vector<LogLine> report_summary(const deck::Plan& plan, const deck::RunReport& report);

class ConsoleLog {
public:
    // Past max_lines the oldest lines are dropped down to keep_lines, so the shown text is rebuilt
    // rarely rather than on every new line.
    static constexpr std::size_t max_lines = 20000, keep_lines = 15000;

    void add(LogLevel level, std::string text);
    void begin_run();
    void follow(const RunView& run);
    void follow(const RunController& runner);
    // follow(), then the points not yet summarised and the summary; finished, the next follow adds nothing.
    void finish(const RunView& run);
    void finish(const RunController& runner);
    void clear();

    [[nodiscard]] const std::deque<LogLine>& lines() const noexcept { return lines_; }
    // Lines ever added, and lines dropped from the front: lines()[k] is line number dropped() + k.
    [[nodiscard]] std::size_t total() const noexcept { return dropped_ + lines_.size(); }
    [[nodiscard]] std::size_t dropped() const noexcept { return dropped_; }
    // The lines as plain text, one per line ("\n"), those the filter shows.
    [[nodiscard]] std::string text(LogFilter filter = LogFilter::all) const;

private:
    void header(const RunView& run, std::size_t index);
    void summarise(const RunView& run, std::size_t index);
    void summarise_before(const RunView& run, std::size_t row);

    std::deque<LogLine> lines_;
    std::size_t dropped_ = 0;
    // The cursor on the controller's state: whether the banner is out, the next iteration to log, the next
    // row to summarise, the row whose header was written last (npos: none).
    bool banner_ = false;
    std::size_t next_iteration_ = 0, next_summary_ = 0, header_row_ = static_cast<std::size_t>(-1);
    bool following_ = false;
};

// Rich text of lines (those the filter shows), for a rich edit: a colour per level on a dark background,
// Consolas at `points` size; `continues` starts with a paragraph break (appending to shown text).
// Lines from index `from` on.
[[nodiscard]] std::string console_rtf(const std::deque<LogLine>& lines, std::size_t from, LogFilter filter,
                                      bool continues, int points = 9);

}  // namespace NiTCAD::app
