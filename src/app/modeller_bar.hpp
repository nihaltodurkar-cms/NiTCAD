// The modeller's controls on Win32 (Unit 26m slice 0; DECISIONS.md MT1-MT3), shared by the 3D modeller and the 2D
// editor: the group switcher in the page's strip, the floating bar of the armed group's tools with its line, the
// camera's buttons in the view's corner (Camera free / fixed, Reset view, Fit, Measure, Mesh), and the value chip at
// the cursor. The bar and the chip are children of the page's panel over its canvas (the canvas clips siblings).
// Commands come back to the page's panel procedure as WM_COMMAND with the ids below.
#pragma once

#include <windows.h>

#include <optional>
#include <string>
#include <vector>

#include "NiTCAD/app/modeller.hpp"

namespace NiTCAD::app {

class ModellerBar {
public:
    // Ids, from `first`: the groups, then the bar's tools, then the corner's buttons.
    static constexpr int group0 = 0, tool0 = group0 + modeller::group_count, corner0 = tool0 + 8;
    static constexpr int camera = 0, reset = 1, fit = 2, measure = 3, mesh = 4, corner_count = 5;

    // `strip`: the page's tool strip; `panel`: the window over the canvas; `proc`: the page's panel procedure (the
    // floating windows forward their commands to it), `self` its user data.
    ModellerBar(HWND strip, HWND panel, int first, bool three_d, WNDPROC proc, void* self);
    ModellerBar(const ModellerBar&) = delete;
    ModellerBar& operator=(const ModellerBar&) = delete;

    void set_font(HFONT font, int dpi);
    // The strip's buttons from `x`: the groups at the left, the corner's at the right of `right`.
    void layout_strip(int right, int pad, int height);

    [[nodiscard]] modeller::Group group() const noexcept { return group_; }
    void set_group(modeller::Group g);
    void set_tool(int tool);  // the armed tool: its button pressed, the bar's line its how
    void set_camera_fixed(bool fixed);
    void set_mesh(bool on);
    void set_measure(bool on);
    // The bar placed over the viewport (`viewport`, panel pixels) beside `selection` (panel pixels).
    void place(const modeller::Box& viewport, const std::optional<modeller::Box>& selection);
    [[nodiscard]] modeller::Placement placement() const noexcept { return placed_; }
    [[nodiscard]] HWND bar() const noexcept { return bar_; }
    [[nodiscard]] std::string line() const { return line_; }

    // A command's meaning; false when the id is not the bar's.
    struct Command {
        enum class Kind { group, tool, duplicate, corner } kind = Kind::group;
        int index = 0;  // the group, the tool's page number, or the corner's button
    };
    [[nodiscard]] std::optional<Command> command(int id) const;

    // The chip at a panel point, or hidden (empty text).
    void chip(int x, int y, const std::string& text);
    [[nodiscard]] std::string chip_text() const { return chip_text_; }
    [[nodiscard]] HWND chip_window() const noexcept { return chip_; }

    // The floating windows above the canvas again (after the canvas is moved).
    void raise();

private:
    void rebuild();
    void size_bar();

    HWND strip_ = nullptr, panel_ = nullptr, bar_ = nullptr, hint_ = nullptr, chip_ = nullptr;
    std::vector<HWND> groups_, tools_, corner_;
    int first_ = 0, dpi_ = 96, tool_ = 0;
    bool three_d_ = false;
    HFONT font_ = nullptr;
    modeller::Group group_ = modeller::Group::solid;
    modeller::Placement placed_;
    modeller::Box viewport_{};
    std::optional<modeller::Box> selection_;
    std::string line_, chip_text_;
};

}  // namespace NiTCAD::app
