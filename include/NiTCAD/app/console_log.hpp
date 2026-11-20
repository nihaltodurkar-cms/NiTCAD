// The console under the inspector (Unit 26j): a running log of what the window and the solver did —
// the runs and checks started, every Newton iteration with its update, each solver point's outcome
// (iterations, last update and residual, seconds), the run's totals, build errors and the status
// messages. The text lives here (UI thread); the window shows it in a read-only edit control.
//
// follow() turns the RunController's rows and iterations into lines incrementally: each iteration once,
// each point summarised once when it converges, when a later point starts, or when the run ends
// (finish()). begin_run() starts following a new check or run (the controller clears its rows).
#pragma once
#include <cstddef>
#include <deque>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/app/run_controller.hpp"

namespace NiTCAD::app {

class ConsoleLog {
public:
    // Past max_lines the oldest lines are dropped down to keep_lines, so the shown text is rebuilt
    // rarely rather than on every new line.
    static constexpr std::size_t max_lines = 20000, keep_lines = 15000;

    void add(std::string line);
    void begin_run();
    // The controller's points() and iterations() so far (RunController::process applied).
    void follow(std::span<const PointRow> rows, std::span<const Iteration> iterations);
    void follow(const RunController& runner) { follow(runner.points(), runner.iterations()); }
    // follow(), then the points not yet summarised, the totals (points, Newton iterations, seconds) and the
    // build errors.
    void finish(std::span<const PointRow> rows, std::span<const Iteration> iterations, double seconds,
                std::span<const Diagnostic> diagnostics);
    void finish(const RunController& runner) {
        finish(runner.points(), runner.iterations(), runner.elapsed_seconds(), runner.diagnostics());
    }
    void clear();

    [[nodiscard]] const std::deque<std::string>& lines() const noexcept { return lines_; }
    // Lines ever added, and lines dropped from the front: lines()[k] is line number dropped() + k.
    [[nodiscard]] std::size_t total() const noexcept { return dropped_ + lines_.size(); }
    [[nodiscard]] std::size_t dropped() const noexcept { return dropped_; }

private:
    void header(const PointRow& row, std::size_t index);
    void summarise(const PointRow& row, std::size_t index);
    void summarise_before(std::span<const PointRow> rows, std::size_t row);

    std::deque<std::string> lines_;
    std::size_t dropped_ = 0;
    // The cursor on the controller's state: the next iteration to log, the next row to summarise, the row
    // whose header was written last (npos: none) and its iterations so far.
    std::size_t next_iteration_ = 0, next_summary_ = 0, header_row_ = static_cast<std::size_t>(-1);
    int row_iterations_ = 0;
    bool following_ = false;
};

}  // namespace NiTCAD::app
