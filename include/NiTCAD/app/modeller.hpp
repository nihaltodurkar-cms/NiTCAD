// The modeller's way of working, without Win32 (Unit 26m slice 0; DECISIONS.md MT1-MT3), shared by the 3D
// modeller (solid_page) and the 2D editor (geometry_page):
//
// - Tool groups (MT2): Solid, Edges, Faces, Transform, Doping, Mesh. What is picked opens its group (nothing:
//   Solid; a solid: Transform; an edge or a vertex: Edges; a face handle or a contact: Faces); the group's tools
//   float in a small bar beside the selection, placed by place_bar. Groups list only the tools built so far; the
//   others say which slice brings theirs.
// - Mouse navigation (MT1, MQ2, MQ3): with the camera free, a left drag on empty background (or on a region's
//   background, which no drag moves) orbits in 3D and pans in 2D; Shift+left pans in 3D (in 2D it keeps the
//   rubber band); the middle button pans, the right one orbits (3D), the wheel zooms. With the camera fixed no drag
//   and no wheel moves the view: every drag goes to the model.
// - The live deck line (MT3): the lines an edit writes, as deck_line finds them between the deck and the edited
//   deck; shown while dragging, the same function on the applied text (MC4).
// - The value chip (MT3): digits typed during a drag (with nm, um or cm; µm when none) set the dragged value
//   exactly; Enter applies, Esc cancels.
// - Slice 1 (MT4, MT5): the Solid group's Fill / Replace (a region's `overlap`, Action::overlap) and the Transform
//   group's modes (Copy, Mirror, Array, Align: Action::mode); an armed mode is shown as armed_key(tool, mode).
#pragma once

#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/app/geometry.hpp"

namespace NiTCAD::app::modeller {

enum class Group : unsigned char { solid, edges, faces, transform, doping, mesh };
inline constexpr int group_count = 6;
[[nodiscard]] std::string_view group_name(Group g) noexcept;

enum class Picked : unsigned char { nothing, solid, edge, face, contact };
[[nodiscard]] Group group_for(Picked p) noexcept;

// A tool of a group's bar: one of the page's tools (the pages' Tool numbers: select 0, box or rectangle 1, prism or
// polygon 2, contact 3, pan 4, measure 5), a Transform mode of the Select tool (MT5: `tool` the Mode), or an action
// done at once on the selection (duplicate; overlap: the picked regions' rule switched, MT4).
enum class Action : unsigned char { tool, mode, duplicate, overlap };
// The Select tool's Transform modes (Unit 26m slice 1, MT5): how a press on a solid acts. Move drags it (none);
// Copy drags a copy (also Alt+drag in any mode); Mirror clicks a face (its centre plane; Ctrl: the face's own plane;
// Shift: a quarter turn about the view's axis); Array drags the pitch, the wheel setting the count; Align clicks a
// face, the other picked solids' matching faces moving to its plane.
// Unit 26m slice 2 (MT6, MT7): the Faces group's Inset (drag the margin in from a face, release, then the depth, a
// click ending it) and Extrude (drag the depth off a face); the Edges group's Chamfer and Fillet (drag in from an edge
// along a prism's axis: the cut, or the radius), Round vertex (the same from a prism's vertex handle) and Taper (drag
// the edge sideways). Each mode belongs to its group (group_of).
enum class Mode : unsigned char { none, copy, mirror, array, align, inset, extrude, chamfer, fillet, round, taper,
                                  profile, implant, junction, max_size, refine };  // slice 4: Doping, Mesh
[[nodiscard]] Group group_of(Mode m) noexcept;
// The bar's armed key: a tool's number, or 100 + a mode.
[[nodiscard]] constexpr int armed_key(int tool, Mode m) noexcept { return m == Mode::none ? tool : 100 + static_cast<int>(m); }
struct ToolEntry {
    std::string_view label;  // with its key
    Action action = Action::tool;
    int tool = 0;
    std::string_view how_2d;  // the bar's line: how to drag it
    std::string_view how_3d;
};
[[nodiscard]] std::span<const ToolEntry> tools(Group g) noexcept;
// The tools of a group on one page: those with no how for that page (3D only, MQ6) left out.
[[nodiscard]] std::vector<ToolEntry> tools_for(Group g, bool three_d);
// What a group gets in a later slice ("" when it has all its tools).
[[nodiscard]] std::string_view later(Group g) noexcept;

// --- MT1: what a press does to the view -----------------------------------------------------------------------------
enum class Button : unsigned char { left, middle, right };
struct Press {
    Button button = Button::left;
    bool background = false;   // nothing a drag could move lies under it
    bool shift = false;
    bool select_tool = true;   // the Select tool armed (a drawing tool draws wherever it is pressed)
    bool pan_tool = false;     // the Pan tool armed or Space held
};
enum class Nav : unsigned char { model, orbit, pan, none };
[[nodiscard]] Nav navigation(bool fixed, bool three_d, const Press& p) noexcept;
[[nodiscard]] constexpr bool wheel_zooms(bool fixed) noexcept { return !fixed; }

// --- MT2: where the floating bar goes ------------------------------------------------------------------------------
struct Box {
    double left = 0, top = 0, right = 0, bottom = 0;
    [[nodiscard]] double width() const noexcept { return right - left; }
    [[nodiscard]] double height() const noexcept { return bottom - top; }
    [[nodiscard]] bool inside(const Box& outer) const noexcept {
        return left >= outer.left && top >= outer.top && right <= outer.right && bottom <= outer.bottom;
    }
    [[nodiscard]] bool overlaps(const Box& o) const noexcept { return left < o.right && o.left < right && top < o.bottom && o.top < bottom; }
    friend bool operator==(const Box&, const Box&) = default;
};
struct Placement {
    Box bar;
    bool corner = false;  // in the viewport's top-left corner (nothing picked, or no room beside the selection)
};
// A bar of `width` × `height` px beside the selection's screen bounds (right, below, left, above, the first that fits
// `gap` px from it and inside the viewport, slid along that side to stay inside), else in the viewport's top-left
// corner `gap` px in.
[[nodiscard]] Placement place_bar(const Box& viewport, const std::optional<Box>& selection, double width, double height, double gap) noexcept;

// --- MT3: the deck line and the typed value -------------------------------------------------------------------------
// The lines of `after` that are not in `before` (a common head and tail of lines set aside; carriage returns and
// blank lines dropped), joined by "  ·  "; when only lines went, "removed: " and those. "" when the same.
[[nodiscard]] std::string deck_line(std::string_view before, std::string_view after);

class Entry {
public:
    // A key typed during a drag: digits, '.', '-', and the unit letters n, u, m, c; false when not taken.
    bool type(wchar_t c);
    bool backspace();  // false when nothing was typed
    void clear() noexcept { text_.clear(); }
    [[nodiscard]] bool empty() const noexcept { return text_.empty(); }
    [[nodiscard]] const std::string& text() const noexcept { return text_; }
    // The length typed, on the quantum (µm when no unit); none while it does not parse.
    [[nodiscard]] std::optional<geo::Coord> value() const;

private:
    std::string text_;
};

// A length for the chip: "0.12 µm".
[[nodiscard]] std::string chip_length(geo::Coord q);

}  // namespace NiTCAD::app::modeller
