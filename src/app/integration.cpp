// The window's integration with the shell (integration.hpp).
#include "NiTCAD/app/integration.hpp"

#include <algorithm>
#include <system_error>

namespace NiTCAD::app {

namespace fs = std::filesystem;

OpenTarget open_target(const fs::path& path) {
    OpenTarget t;
    t.path = path;
    if (path.extension() == ".deck") {
        t.deck = true;
        return t;
    }
    t.device = fs::path(path).replace_extension(".nitdev");
    fs::path mesh = fs::path(path).replace_extension(".nitmesh");
    std::error_code ec;
    if (fs::exists(mesh, ec)) t.mesh = std::move(mesh);
    return t;
}

std::expected<fs::path, std::string> dropped_file(std::span<const fs::path> files) {
    if (files.size() != 1) return std::unexpected(std::string("Drop one deck or result file at a time."));
    std::error_code ec;
    if (fs::is_directory(files.front(), ec)) {
        return std::unexpected(std::string("A folder does not open: drop a deck or a result file."));
    }
    return files.front();
}

TaskbarProgress taskbar_progress(RunState s, bool run, double fraction, const TaskbarProgress& previous) noexcept {
    switch (s) {
    case RunState::checking:
        return {};  // a check clears an earlier failure
    case RunState::running:
    case RunState::stopping:
        if (!run) return {};
        if (fraction < 0.0) return {TaskbarFlag::indeterminate, 0.0};
        return {TaskbarFlag::normal, std::clamp(fraction, 0.0, 1.0)};
    case RunState::failed:
        if (!run) return previous;
        // The bar kept where the run stopped (full when its length was not known), in the error state.
        return {TaskbarFlag::error, previous.flag == TaskbarFlag::normal && previous.fraction > 0.0 ? previous.fraction : 1.0};
    case RunState::done:
        return run ? TaskbarProgress{} : previous;
    case RunState::idle:
        return previous;
    }
    return previous;
}

bool dark_title_bar(std::optional<unsigned long> apps_use_light_theme) noexcept {
    return apps_use_light_theme.has_value() && *apps_use_light_theme == 0;
}

}  // namespace NiTCAD::app
