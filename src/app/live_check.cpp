// Checking the deck while it is edited (live_check.hpp).
#include "NiTCAD/app/live_check.hpp"

#include <condition_variable>
#include <mutex>
#include <utility>

namespace NiTCAD::app {

struct LiveChecker::Shared {
    std::mutex mutex;
    std::condition_variable_any changed;
    struct Job {
        std::uint64_t generation = 0;
        std::string text;
        std::filesystem::path path;
    };
    std::optional<Job> next;       // the newest text not started
    std::uint64_t running = 0;     // the generation being checked, 0 when idle
    struct Done {
        std::uint64_t generation = 0;
        CheckResult result;
    };
    std::optional<Done> done;      // the last check ended, not yet taken by process()
    std::uint64_t ended = 0;       // the newest generation whose check has ended
    std::uint64_t superseded = 0;  // ended checks overwritten before process() saw them
};

LiveChecker::LiveChecker(Wake wake, Check check) : shared_(std::make_shared<Shared>()) {
    if (!check) check = [](const std::string& text, const std::filesystem::path& path) { return check_deck(text, path); };
    worker_ = std::jthread([s = shared_, wake = std::move(wake), check = std::move(check)](std::stop_token stop) {
        for (;;) {
            Shared::Job job;
            {
                std::unique_lock lock(s->mutex);
                if (!s->changed.wait(lock, stop, [&] { return s->next.has_value(); })) return;
                job = std::move(*s->next);
                s->next.reset();
                s->running = job.generation;
            }
            CheckResult result = check(job.text, job.path);
            {
                std::scoped_lock lock(s->mutex);
                if (s->done) ++s->superseded;
                s->done = Shared::Done{job.generation, std::move(result)};
                s->ended = job.generation;
                s->running = 0;
            }
            s->changed.notify_all();
            if (wake) wake();
        }
    });
}

LiveChecker::~LiveChecker() {
    worker_.request_stop();
    // jthread joins; a check in progress runs to its end first.
}

std::uint64_t LiveChecker::submit(std::string text, std::filesystem::path path) {
    {
        std::scoped_lock lock(shared_->mutex);
        shared_->next = Shared::Job{++submitted_, std::move(text), std::move(path)};
    }
    shared_->changed.notify_all();
    return submitted_;
}

bool LiveChecker::process() {
    std::optional<Shared::Done> done;
    {
        std::scoped_lock lock(shared_->mutex);
        done = std::move(shared_->done);
        shared_->done.reset();
        dropped_ += shared_->superseded;
        shared_->superseded = 0;
    }
    if (!done) return false;
    if (done->generation != submitted_) {
        ++dropped_;  // a newer text was submitted: this result describes a deck no longer shown
        return false;
    }
    result_ = std::move(done->result);
    applied_ = done->generation;
    return true;
}

void LiveChecker::wait() {
    std::unique_lock lock(shared_->mutex);
    shared_->changed.wait(lock, [&] { return shared_->ended >= submitted_; });
}

std::optional<deck::Plan> LiveChecker::take_plan() {
    auto plan = std::move(result_.plan);
    result_.plan.reset();
    return plan;
}

}  // namespace NiTCAD::app
