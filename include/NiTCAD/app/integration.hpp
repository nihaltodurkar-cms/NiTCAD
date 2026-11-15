// The window's integration with the Windows shell (Unit 26d, DECISIONS.md ZD4-ZD6), its decisions as plain
// functions the window follows: what a dropped file opens, the taskbar button's progress through a run, and
// whether the title bar is dark.
#pragma once

#include <expected>
#include <filesystem>
#include <optional>
#include <span>
#include <string>

#include "NiTCAD/app/run_controller.hpp"

namespace NiTCAD::app {

// What a path opens, as the recent-files list opens it (ZD4): a .deck as a deck; any other file as a result
// with the .nitdev of its name beside it and, when present, the .nitmesh.
struct OpenTarget {
    bool deck = false;
    std::filesystem::path path;  // the deck, or the result file
    std::filesystem::path device, mesh;
};
[[nodiscard]] OpenTarget open_target(const std::filesystem::path& path);

// The file a drop opens (ZD4): exactly one, and not a folder; the reason otherwise, for the status bar.
[[nodiscard]] std::expected<std::filesystem::path, std::string> dropped_file(std::span<const std::filesystem::path> files);

// The taskbar button's progress (ZD5, ZQ3). While a run goes on, its fraction (progress_fraction; indeterminate
// when the run does not know its length); a failed run shows the error state until the next run or check; a run
// that ends (done, stopped early or by Stop) and a check clear it. A check's own failure (the deck does not
// build) is shown in the deck, not on the taskbar.
enum class TaskbarFlag : unsigned char { none, normal, indeterminate, error };
struct TaskbarProgress {
    TaskbarFlag flag = TaskbarFlag::none;
    double fraction = 0.0;  // [0, 1], for normal and error
    friend bool operator==(const TaskbarProgress&, const TaskbarProgress&) = default;
};
// The state after the run controller's state `s`; `run` when the controller's work is a run (not a check or an
// open); `fraction` as progress_fraction gives it (-1: unknown).
[[nodiscard]] TaskbarProgress taskbar_progress(RunState s, bool run, double fraction, const TaskbarProgress& previous) noexcept;

// Whether the title bar is dark (ZD6, ZQ4): when Windows' app theme is dark, the AppsUseLightTheme value 0;
// light when it is any other value or missing.
[[nodiscard]] bool dark_title_bar(std::optional<unsigned long> apps_use_light_theme) noexcept;

}  // namespace NiTCAD::app
