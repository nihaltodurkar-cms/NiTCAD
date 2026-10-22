// The run controller (run_controller.hpp).
#include "NiTCAD/app/run_controller.hpp"

#include <mutex>
#include <utility>

namespace NiTCAD::app {

std::string_view to_string(RunState s) noexcept {
    switch (s) {
    case RunState::idle: return "idle";
    case RunState::checking: return "checking";
    case RunState::running: return "running";
    case RunState::stopping: return "stopping";
    case RunState::done: return "done";
    case RunState::failed: return "failed";
    }
    return "?";
}

struct RunController::Event {
    enum class Kind : unsigned char { built, progress, finished } kind = Kind::progress;
    bool run = false;
    std::size_t planned = 0;
    solve::Progress progress{};
    double seconds = 0.0;  // on the point so far
    std::vector<Diagnostic> diagnostics;
    std::optional<deck::Plan> plan;
    std::optional<deck::RunReport> report;
    std::optional<base::Error> error;
};

// The queue between the worker and the UI thread. The worker wakes the UI thread only when it puts an
// event into an empty queue, so a fast solver does not flood the window's message queue.
struct RunController::Shared {
    std::mutex mutex;
    std::vector<Event> queue;
};

namespace {

std::size_t plan_point_count(const deck::Plan& p) {
    return p.analysis == deck::Analysis::sweep || p.analysis == deck::Analysis::small_signal ? p.points.size() : 0;
}

}  // namespace

RunController::RunController(Wake wake) : wake_(std::move(wake)), shared_(std::make_shared<Shared>()) {}

RunController::~RunController() {
    if (worker_.joinable()) {
        worker_.request_stop();
        worker_.join();
    }
}

bool RunController::start_check(std::string text, std::filesystem::path path) {
    if (busy()) return false;
    start(false, std::move(text), std::move(path));
    return true;
}

bool RunController::start_run(std::string text, std::filesystem::path path) {
    if (busy()) return false;
    start(true, std::move(text), std::move(path));
    return true;
}

void RunController::start(bool run, std::string text, std::filesystem::path path) {
    if (worker_.joinable()) worker_.join();  // the last worker has finished (not busy)
    process();                               // and nothing of it is left queued
    state_ = run ? RunState::running : RunState::checking;
    diagnostics_.clear();
    plan_.reset();
    report_.reset();
    error_.reset();
    points_.clear();
    iterations_.clear();
    planned_points_ = 0;
    started_ = finished_ = std::chrono::steady_clock::now();
    worker_ = std::jthread([shared = shared_, wake = wake_, run, text = std::move(text),
                            path = std::move(path)](std::stop_token stop) mutable {
        const auto post = [&](Event e) {
            bool was_empty = false;
            {
                std::lock_guard lock(shared->mutex);
                was_empty = shared->queue.empty();
                shared->queue.push_back(std::move(e));
            }
            if (was_empty && wake) wake();
        };
        CheckResult checked = check_deck(text, path);
        if (!checked.plan) {
            Event e{.kind = Event::Kind::finished, .run = run};
            e.diagnostics = std::move(checked.diagnostics);
            e.error = base::Error{base::ErrorCode::invalid_input,
                                  e.diagnostics.empty() ? std::string("the deck does not build") : e.diagnostics.front().message,
                                  std::nullopt};
            post(std::move(e));
            return;
        }
        if (!run) {
            Event e{.kind = Event::Kind::finished, .run = false};
            e.plan = std::move(checked.plan);
            post(std::move(e));
            return;
        }
        post(Event{.kind = Event::Kind::built, .run = true, .planned = plan_point_count(*checked.plan)});
        // Progress: one event per Newton iteration, with the wall time spent on its point.
        auto point_started = std::chrono::steady_clock::now();
        solve::Phase phase = solve::Phase::equilibrium;
        std::size_t point = static_cast<std::size_t>(-1);
        solve::RunControl control{.stop = stop, .progress = {}};
        control.progress = [&](const solve::Progress& p) {
            const auto now = std::chrono::steady_clock::now();
            if (p.phase != phase || p.point != point) {
                // A new point starts where the last iteration of the previous one ended.
                if (point != static_cast<std::size_t>(-1)) point_started = now;
                phase = p.phase;
                point = p.point;
            }
            post(Event{.kind = Event::Kind::progress,
                       .progress = p,
                       .seconds = std::chrono::duration<double>(now - point_started).count()});
        };
        auto report = deck::run(*checked.plan, control);
        Event e{.kind = Event::Kind::finished, .run = true};
        e.plan = std::move(checked.plan);
        if (report) e.report = std::move(*report);
        else e.error = std::move(report.error());
        post(std::move(e));
    });
}

void RunController::request_stop() {
    if (state_ != RunState::running) return;
    worker_.request_stop();
    state_ = RunState::stopping;
}

void RunController::apply(Event& e) {
    switch (e.kind) {
    case Event::Kind::built:
        planned_points_ = e.planned;
        break;
    case Event::Kind::progress: {
        const solve::Progress& p = e.progress;
        if (points_.empty() || points_.back().phase != p.phase || points_.back().point != p.point) {
            points_.push_back({.phase = p.phase, .point = p.point, .point_count = p.point_count});
        }
        PointRow& row = points_.back();
        row.iterations = p.iteration;
        row.update = p.update;
        row.residual = p.residual;
        row.converged = p.converged;
        row.time_s = p.time_s;
        row.frequency_Hz = p.frequency_Hz;
        row.seconds = e.seconds;
        iterations_.push_back({iterations_.size() + 1, p.update, points_.size() - 1});
        break;
    }
    case Event::Kind::finished:
        finished_ = std::chrono::steady_clock::now();
        diagnostics_ = std::move(e.diagnostics);
        plan_ = std::move(e.plan);
        report_ = std::move(e.report);
        error_ = std::move(e.error);
        if (!e.run) state_ = diagnostics_.empty() ? RunState::idle : RunState::failed;
        else state_ = error_ ? RunState::failed : RunState::done;
        break;
    }
}

bool RunController::process() {
    std::vector<Event> events;
    {
        std::lock_guard lock(shared_->mutex);
        events.swap(shared_->queue);
    }
    for (Event& e : events) apply(e);
    return !events.empty();
}

void RunController::wait() {
    if (worker_.joinable()) worker_.join();
}

std::optional<base::Error> RunController::stopped_early() const {
    if (!report_) return std::nullopt;
    if (report_->sweep && report_->sweep->stopped) return report_->sweep->stopped;
    if (report_->transient && report_->transient->stopped) return report_->transient->stopped;
    if (report_->small_signal && report_->small_signal->stopped) return report_->small_signal->stopped;
    return std::nullopt;
}

double RunController::elapsed_seconds() const {
    const auto end = busy() ? std::chrono::steady_clock::now() : finished_;
    return std::chrono::duration<double>(end - started_).count();
}

std::optional<RunController::Taken> RunController::take_results() {
    if (state_ != RunState::done || !plan_ || !report_) return std::nullopt;
    Taken t{std::move(*plan_), std::move(*report_)};
    plan_.reset();
    report_.reset();
    return t;
}

}  // namespace NiTCAD::app
