// Running a deck off the UI thread (Unit 26, DECISIONS.md VA3, WG1, WG2, WG3; ARCHITECTURE.md 6.9).
//
// A RunController owns one worker thread per check or run. The worker builds the deck
// (deck::build) and, for a run, runs it (deck::run, exactly what `nitcad run` does) with a
// RunControl whose stop token the controller holds. The solver's progress events, the build's
// outcome and the run's report are queued on the worker; the worker then calls the wake function,
// which must be cheap and thread-safe (the window posts itself a message). The UI thread calls
// process() to apply the queued events: every observable member is read and changed on the UI
// thread only, so the window never waits on the solve and the solver never calls UI code.
//
// States: idle -> checking -> idle (built) | failed (the deck does not build);
//         idle -> running -> done (finished, or stopped early: cancelled, not converged) | failed;
//         running -> stopping (request_stop) -> done (the completed part kept) | failed;
//         done | failed -> checking | running (a new check or run).
// A run that stops early is still done: its report holds the completed points and why it stopped,
// and its outputs are written, as `nitcad run` writes them.
#pragma once

#include <chrono>
#include <cstddef>
#include <filesystem>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

#include "NiTCAD/app/document.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/deck/build.hpp"
#include "NiTCAD/deck/run.hpp"
#include "NiTCAD/solve/control.hpp"

namespace NiTCAD::app {

enum class RunState : unsigned char { idle, checking, running, stopping, done, failed };

[[nodiscard]] std::string_view to_string(RunState s) noexcept;

// One solver point as the run monitor shows it: a bias point, a transient step attempt or a
// small-signal frequency, with its last Newton iteration so far.
struct PointRow {
    solve::Phase phase = solve::Phase::equilibrium;
    std::size_t point = 0;        // the solver's point index (0-based)
    std::size_t point_count = 0;  // 0 when not known in advance (transient, trace)
    int iterations = 0;
    double update = 0.0;          // the last iteration's largest scaled correction
    double residual = 0.0;
    bool converged = false;
    double time_s = 0.0;          // transient: the step's end time
    double frequency_Hz = 0.0;    // small signal: the frequency
    double seconds = 0.0;         // wall time spent on the point so far
};

// One Newton iteration for the live convergence plot: its number over the whole run (from 1), its
// update, and the row it belongs to.
struct Iteration {
    std::size_t number = 0;
    double update = 0.0;
    std::size_t row = 0;
};

class RunController {
public:
    using Wake = std::function<void()>;

    explicit RunController(Wake wake = {});
    ~RunController();  // stops a run in progress and waits for its worker
    RunController(const RunController&) = delete;
    RunController& operator=(const RunController&) = delete;

    // Starts building (check) or building and running (run) the deck text as the file `path`
    // (check_deck). Refused (false) while checking, running or stopping.
    bool start_check(std::string text, std::filesystem::path path);
    bool start_run(std::string text, std::filesystem::path path);
    // Asks a run to stop at the solver's next safe point (running -> stopping); otherwise nothing.
    void request_stop();

    // UI thread: applies the queued events. Returns whether anything changed.
    bool process();
    // Blocks until the worker has finished (tests and the self-test); process() still applies its
    // events.
    void wait();

    [[nodiscard]] RunState state() const noexcept { return state_; }
    [[nodiscard]] bool busy() const noexcept {
        return state_ == RunState::checking || state_ == RunState::running || state_ == RunState::stopping;
    }
    // The last check's or run's build errors (empty when it built).
    [[nodiscard]] const std::vector<Diagnostic>& diagnostics() const noexcept { return diagnostics_; }
    // The plan the last check or run built (kept after the run).
    [[nodiscard]] const std::optional<deck::Plan>& plan() const noexcept { return plan_; }
    // The last run's report (done) or the error it failed with (failed).
    [[nodiscard]] const std::optional<deck::RunReport>& report() const noexcept { return report_; }
    [[nodiscard]] const std::optional<base::Error>& error() const noexcept { return error_; }
    // Why a done run stopped early (cancelled, not converged, ...), if it did.
    [[nodiscard]] std::optional<base::Error> stopped_early() const;
    [[nodiscard]] const std::vector<PointRow>& points() const noexcept { return points_; }
    [[nodiscard]] const std::vector<Iteration>& iterations() const noexcept { return iterations_; }
    // The points the plan asks for (a sweep's or small-signal run's), 0 when not known in advance.
    [[nodiscard]] std::size_t planned_points() const noexcept { return planned_points_; }
    // Wall time of the last run or check, so far while busy (UI thread).
    [[nodiscard]] double elapsed_seconds() const;
    // Takes the plan and report out (for the results workspace); the state stays done.
    struct Taken {
        deck::Plan plan;
        deck::RunReport report;
    };
    [[nodiscard]] std::optional<Taken> take_results();
    // Takes the plan of the last check out (a deck opened in the window: its device for the views).
    [[nodiscard]] std::optional<deck::Plan> take_plan();

private:
    struct Event;
    struct Shared;
    void start(bool run, std::string text, std::filesystem::path path);
    void apply(Event& e);

    Wake wake_;
    std::shared_ptr<Shared> shared_;
    std::jthread worker_;
    RunState state_ = RunState::idle;
    std::vector<Diagnostic> diagnostics_;
    std::optional<deck::Plan> plan_;
    std::optional<deck::RunReport> report_;
    std::optional<base::Error> error_;
    std::vector<PointRow> points_;
    std::vector<Iteration> iterations_;
    std::size_t planned_points_ = 0;
    std::chrono::steady_clock::time_point started_{};
    std::chrono::steady_clock::time_point finished_{};
};

}  // namespace NiTCAD::app
