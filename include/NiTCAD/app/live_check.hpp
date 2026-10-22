// Checking the deck while it is edited (Unit 26b, DECISIONS.md XD1, XG1): its own worker thread, apart
// from the run controller, so the deck is checked during a run too.
//
// The window submits the text after a pause in typing; the worker checks the newest text submitted
// (check_deck: parsed and built, nothing solved). A build cannot be stopped part way, so a check of a
// text since replaced runs to its end and its result is dropped: process() applies a result only if
// it is of the newest submission. The worker calls the wake function (cheap, thread-safe) when a
// check ends; process() and every observer run on the UI thread.
#pragma once

#include <cstdint>
#include <filesystem>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <thread>

#include "NiTCAD/app/document.hpp"

namespace NiTCAD::app {

class LiveChecker {
public:
    using Wake = std::function<void()>;
    // The check a worker makes; check_deck unless a test holds the worker.
    using Check = std::function<CheckResult(const std::string& text, const std::filesystem::path& path)>;

    explicit LiveChecker(Wake wake = {}, Check check = {});
    ~LiveChecker();  // waits for a check in progress to end
    LiveChecker(const LiveChecker&) = delete;
    LiveChecker& operator=(const LiveChecker&) = delete;

    // A new text to check; returns its generation (1, 2, ...). A text not yet started is replaced.
    std::uint64_t submit(std::string text, std::filesystem::path path);
    // UI thread: applies the newest submission's result if its check has ended. Returns whether it did.
    bool process();
    // Blocks until the newest submission's check has ended (tests); process() still applies it.
    void wait();

    // The generation of the newest submission, and of the result applied (0 before any).
    [[nodiscard]] std::uint64_t submitted() const noexcept { return submitted_; }
    [[nodiscard]] std::uint64_t applied() const noexcept { return applied_; }
    // Whether the newest submission's result has not been applied yet.
    [[nodiscard]] bool pending() const noexcept { return applied_ != submitted_; }
    // The applied result: the deck's errors, or its plan when it built.
    [[nodiscard]] const CheckResult& result() const noexcept { return result_; }
    // Takes the applied result's plan out (the structure preview).
    [[nodiscard]] std::optional<deck::Plan> take_plan();
    // Results of older submissions dropped so far.
    [[nodiscard]] std::uint64_t dropped() const noexcept { return dropped_; }

private:
    struct Shared;
    std::shared_ptr<Shared> shared_;
    std::jthread worker_;
    std::uint64_t submitted_ = 0;
    std::uint64_t applied_ = 0;
    std::uint64_t dropped_ = 0;
    CheckResult result_;
};

}  // namespace NiTCAD::app
