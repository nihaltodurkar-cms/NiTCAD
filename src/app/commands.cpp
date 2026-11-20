// The command list (commands.hpp).
#include "NiTCAD/app/commands.hpp"

#include <array>

namespace NiTCAD::app {

namespace {

// Virtual-key codes (winuser.h), named here so this file needs no Windows header.
constexpr std::uint16_t vk_tab = 0x09, vk_space = 0x20, vk_end = 0x23, vk_home = 0x24, vk_left = 0x25, vk_right = 0x27,
                        vk_f4 = 0x73, vk_f5 = 0x74, vk_f6 = 0x75, vk_f7 = 0x76, vk_f8 = 0x77;

constexpr std::array<std::string_view, 7> bar{"&File", "&Edit", "&Device", "&Mesh", "&Physics", "&Simulation", "&View"};

constexpr std::array<CommandInfo, 68> list{{
    {Command::new_deck, "&File", "&New deck", {'N', true}, "Ctrl+N"},
    {Command::open_deck, "&File", "&Open deck...", {'O', true}, "Ctrl+O"},
    {Command::open_result, "&File", "Open &result...", {'R', true, true}, "Ctrl+Shift+R"},
    {Command::compare_result, "&File", "Co&mpare with a result...", {'C', true, true}, "Ctrl+Shift+C"},
    {Command::save, "&File", "&Save", {'S', true}, "Ctrl+S"},
    {Command::save_as, "&File", "Save &as...", {'S', true, true}, "Ctrl+Shift+S"},
    {Command::exit, "&File", "E&xit", {vk_f4, false, false, true}, "Alt+F4"},
    {Command::undo, "&Edit", "&Undo", {'Z', true}, "Ctrl+Z"},
    {Command::cut, "&Edit", "Cu&t", {'X', true}, "Ctrl+X"},
    {Command::copy, "&Edit", "&Copy", {'C', true}, "Ctrl+C"},
    {Command::paste, "&Edit", "&Paste", {'V', true}, "Ctrl+V"},
    {Command::select_all, "&Edit", "Select &all", {'A', true}, "Ctrl+A"},
    {Command::find_error, "&Edit", "Go to the &error", {vk_f8}, "F8"},
    {Command::complete, "&Edit", "Co&mplete", {vk_space, true}, "Ctrl+Space"},
    {Command::redo, "&Edit", "&Redo", {'Y', true}, "Ctrl+Y"},
    {Command::duplicate, "&Edit", "&Duplicate", {'D', true}, "Ctrl+D"},
    {Command::delete_selection, "&Edit", "De&lete", {}, "Del", {}},
    {Command::tool_select, "&Edit", "&Select", {}, "V", "&Geometry tools"},
    {Command::tool_rectangle, "&Edit", "&Rectangle", {}, "R", "&Geometry tools"},
    {Command::tool_polygon, "&Edit", "P&olygon", {}, "P", "&Geometry tools"},
    {Command::tool_contact, "&Edit", "Co&ntact", {}, "C", "&Geometry tools"},
    {Command::tool_pan, "&Edit", "&Pan", {}, "H", "&Geometry tools"},
    {Command::tool_measure, "&Edit", "&Measure", {}, "M", "&Geometry tools"},
    {Command::check, "&Simulation", "&Check deck", {vk_f7}, "F7"},
    {Command::run, "&Simulation", "&Run", {vk_f5}, "F5"},
    {Command::stop, "&Simulation", "S&top", {vk_f5, false, true}, "Shift+F5"},
    {Command::compare_previous, "&Simulation", "Compare with the &previous run", {'V', true, true}, "Ctrl+Shift+V"},
    {Command::stop_comparing, "&Simulation", "Stop co&mparing", {'X', true, true}, "Ctrl+Shift+X"},
    // Transient playback (Unit 26c, YD5).
    {Command::play_pause, "&Simulation", "Pla&y / pause snapshots", {vk_space, true, true}, "Ctrl+Shift+Space"},
    {Command::snapshot_first, "&Simulation", "&First snapshot", {vk_home, false, false, true}, "Alt+Home"},
    {Command::snapshot_previous, "&Simulation", "Pre&vious snapshot", {vk_left, false, false, true}, "Alt+Left"},
    {Command::snapshot_next, "&Simulation", "&Next snapshot", {vk_right, false, false, true}, "Alt+Right"},
    {Command::snapshot_last, "&Simulation", "&Last snapshot", {vk_end, false, false, true}, "Alt+End"},
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
    {Command::cut_line_mode, "&View", "Dr&aw cut line", {'L', true, true}, "Ctrl+Shift+L"},
    {Command::mesh_quality, "&View", "Mesh &quality", {'Q', true, true}, "Ctrl+Shift+Q"},
    {Command::show_regions, "&View", "S&how every region", {'H', true, true}, "Ctrl+Shift+H"},
    {Command::hover_readout, "&View", "Displa&y hover readout", {}, ""},
    // Unit 26h (JD4): the numpad keys are the view's own (viewport.cpp), as the arrows are, so the deck
    // editor still types digits.
    {Command::view_front, "&View", "&Front", {}, "Num 1", "Loo&k from"},
    {Command::view_back, "&View", "&Back", {}, "Ctrl+Num 1", "Loo&k from"},
    {Command::view_left, "&View", "&Left", {}, "Ctrl+Num 3", "Loo&k from"},
    {Command::view_right, "&View", "&Right", {}, "Num 3", "Loo&k from"},
    {Command::view_top, "&View", "&Top", {}, "Num 7", "Loo&k from"},
    {Command::view_bottom, "&View", "B&ottom", {}, "Ctrl+Num 7", "Loo&k from"},
    {Command::view_isometric, "&View", "&Isometric", {}, "Num 0", "Loo&k from"},
    {Command::next_tab, "&View", "Next &tab", {vk_tab, true}, "Ctrl+Tab"},
    {Command::previous_tab, "&View", "Pre&vious tab", {vk_tab, true, true}, "Ctrl+Shift+Tab"},
    {Command::next_pane, "&View", "Ne&xt pane", {vk_f6}, "F6"},
    {Command::tab_start, "&View", "Start pa&ge", {'1', true}, "Ctrl+1"},
    {Command::tab_geometry, "&View", "Geometry editor (&2)", {'2', true}, "Ctrl+2"},
    {Command::tab_view, "&View", "Device vie&w", {'3', true}, "Ctrl+3"},
    {Command::tab_curves, "&View", "&Curves", {'4', true}, "Ctrl+4"},
    {Command::tab_deck, "&View", "D&eck", {'5', true}, "Ctrl+5"},
    {Command::tab_run, "&View", "Ru&n monitor", {'6', true}, "Ctrl+6"},
    {Command::tab_results, "&View", "Re&sults", {'7', true}, "Ctrl+7"},
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
