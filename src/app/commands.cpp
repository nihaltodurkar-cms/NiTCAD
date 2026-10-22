// The command list (commands.hpp).
#include "NiTCAD/app/commands.hpp"

#include <array>

namespace NiTCAD::app {

namespace {

// Virtual-key codes (winuser.h), named here so this file needs no Windows header.
constexpr std::uint16_t vk_tab = 0x09, vk_f4 = 0x73, vk_f5 = 0x74, vk_f6 = 0x75, vk_f7 = 0x76, vk_f8 = 0x77;

constexpr std::array<std::string_view, 7> bar{"&File", "&Edit", "&Device", "&Mesh", "&Physics", "&Simulation", "&View"};

constexpr std::array<CommandInfo, 35> list{{
    {Command::new_deck, "&File", "&New deck", {'N', true}, "Ctrl+N"},
    {Command::open_deck, "&File", "&Open deck...", {'O', true}, "Ctrl+O"},
    {Command::open_result, "&File", "Open &result...", {'R', true, true}, "Ctrl+Shift+R"},
    {Command::save, "&File", "&Save", {'S', true}, "Ctrl+S"},
    {Command::save_as, "&File", "Save &as...", {'S', true, true}, "Ctrl+Shift+S"},
    {Command::exit, "&File", "E&xit", {vk_f4, false, false, true}, "Alt+F4"},
    {Command::undo, "&Edit", "&Undo", {'Z', true}, "Ctrl+Z"},
    {Command::cut, "&Edit", "Cu&t", {'X', true}, "Ctrl+X"},
    {Command::copy, "&Edit", "&Copy", {'C', true}, "Ctrl+C"},
    {Command::paste, "&Edit", "&Paste", {'V', true}, "Ctrl+V"},
    {Command::select_all, "&Edit", "Select &all", {'A', true}, "Ctrl+A"},
    {Command::find_error, "&Edit", "Go to the &error", {vk_f8}, "F8"},
    {Command::check, "&Simulation", "&Check deck", {vk_f7}, "F7"},
    {Command::run, "&Simulation", "&Run", {vk_f5}, "F5"},
    {Command::stop, "&Simulation", "S&top", {vk_f5, false, true}, "Shift+F5"},
    {Command::fit, "&View", "&Fit", {}, "Home"},
    {Command::orbit_left, "&View", "Orbit &left", {}, "Left"},
    {Command::orbit_right, "&View", "Orbit &right", {}, "Right"},
    {Command::orbit_up, "&View", "Orbit &up", {}, "Up"},
    {Command::orbit_down, "&View", "Orbit &down", {}, "Down"},
    {Command::zoom_in, "&View", "Zoom &in", {}, "+"},
    {Command::zoom_out, "&View", "Zoom &out", {}, "-"},
    {Command::pan_left, "&View", "Pan left", {}, "Shift+Left"},
    {Command::pan_right, "&View", "Pan right", {}, "Shift+Right"},
    {Command::pan_up, "&View", "Pan up", {}, "Shift+Up"},
    {Command::pan_down, "&View", "Pan down", {}, "Shift+Down"},
    {Command::perspective, "&View", "&Perspective", {'P', true, true}, "Ctrl+Shift+P"},
    {Command::probe_mode, "&View", "Pro&be", {'B', true, true}, "Ctrl+Shift+B"},
    {Command::measure_mode, "&View", "&Measure", {'M', true, true}, "Ctrl+Shift+M"},
    {Command::next_tab, "&View", "Next &tab", {vk_tab, true}, "Ctrl+Tab"},
    {Command::previous_tab, "&View", "Pre&vious tab", {vk_tab, true, true}, "Ctrl+Shift+Tab"},
    {Command::next_pane, "&View", "Ne&xt pane", {vk_f6}, "F6"},
    {Command::tab_view, "&View", "Device vie&w", {'1', true}, "Ctrl+1"},
    {Command::tab_curves, "&View", "&Curves", {'2', true}, "Ctrl+2"},
    {Command::tab_deck, "&View", "D&eck", {'3', true}, "Ctrl+3"},
}};

constexpr std::span<const CommandInfo> fixed{list};

}  // namespace

std::span<const CommandInfo> commands() noexcept { return fixed; }

const CommandInfo* find_command(Command id) noexcept {
    for (const CommandInfo& c : fixed) {
        if (c.id == id) return &c;
    }
    return nullptr;
}

std::span<const std::string_view> menu_bar() noexcept { return bar; }

}  // namespace NiTCAD::app
