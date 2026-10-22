// The application's commands (Unit 26, DECISIONS.md VA1, WG5): every command has a menu item, and the
// frequent ones an accelerator, so everything is reachable from the keyboard (Alt opens the menus).
// The window builds its menu bar and accelerator table from this list; the test checks it.
#pragma once

#include <cstdint>
#include <span>
#include <string_view>

namespace NiTCAD::app {

enum class Command : std::uint16_t {
    // File
    new_deck = 100, open_deck, open_result, save, save_as, recent_first, exit = recent_first + 16,
    // Edit
    undo = 200, cut, copy, paste, select_all, find_error,
    // Simulation (run)
    check = 300, run, stop,
    // View
    fit = 400, orbit_left, orbit_right, orbit_up, orbit_down, zoom_in, zoom_out, pan_left, pan_right, pan_up,
    pan_down, perspective, probe_mode, measure_mode, next_tab, previous_tab, next_pane,
    tab_view, tab_curves, tab_deck,
    // The Device, Mesh, Physics and Simulation menus' sections and keys (VQ13) are numbered from here.
    section_first = 1000,
    key_first = 2000,
};

// Keys as the window's accelerator table takes them: a virtual-key code and modifiers.
struct Accelerator {
    std::uint16_t key = 0;  // a Win32 virtual-key code (VK_F5, 'S', ...); 0: none
    bool ctrl = false;
    bool shift = false;
    bool alt = false;
};

struct CommandInfo {
    Command id;
    std::string_view menu;  // the top-level menu it is in ("&File", ...)
    std::string_view text;  // the item's text, with its '&' mnemonic
    Accelerator accelerator;
    std::string_view accelerator_text;  // shown in the menu, e.g. "Ctrl+S"
};

// The view's camera commands have no accelerator: arrows, Shift+arrows, + and - and Home act in the
// focused view (accelerators would take them from the deck editor); their keys are shown as hints.
// Every fixed command, in menu order (the recent files, the deck sections and keys are added by the
// window in their places).
[[nodiscard]] std::span<const CommandInfo> commands() noexcept;
[[nodiscard]] const CommandInfo* find_command(Command id) noexcept;
// The top-level menus in order: File, Edit, Device, Mesh, Physics, Simulation, View (VQ13).
[[nodiscard]] std::span<const std::string_view> menu_bar() noexcept;

}  // namespace NiTCAD::app
