// What the application remembers per user (Unit 26, DECISIONS.md VA1): the recently opened files and
// the window's layout, kept as "key = value" lines in a text file (%APPDATA%\NiTCAD\app.ini in the
// window). Unknown and malformed lines are ignored, so an old or damaged file only loses those
// settings.
#pragma once

#include <cstddef>
#include <expected>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/base/error.hpp"

namespace NiTCAD::app {

struct Layout {
    int x = -1, y = -1, width = 0, height = 0;  // the window's normal rectangle (0: the default)
    bool maximized = false;
    int tree_width = 0;       // the scene tree's width, 96-DPI pixels (0: the default)
    int inspector_width = 0;  // the inspector's
    int deck_height = 0;      // the deck pane's share, unused when 0
};

struct AppState {
    static constexpr std::size_t max_recent = 8;
    std::vector<std::filesystem::path> recent;  // most recent first, no duplicates
    Layout layout;

    // Moves (or adds) the path to the front, dropping the oldest beyond max_recent.
    void add_recent(const std::filesystem::path& path);
    void remove_recent(const std::filesystem::path& path);
};

[[nodiscard]] std::string serialize(const AppState& state);
[[nodiscard]] AppState parse_state(std::string_view text);

// The state in a file; a missing file is the default state. Errors (resource_exhausted): save cannot
// write the file.
[[nodiscard]] AppState load_state(const std::filesystem::path& file);
[[nodiscard]] std::expected<void, base::Error> save_state(const AppState& state, const std::filesystem::path& file);

}  // namespace NiTCAD::app
