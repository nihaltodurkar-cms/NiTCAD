// The modeller's way of working (modeller.hpp; DECISIONS.md Unit 26m, MT1-MT3).
#include "NiTCAD/app/modeller.hpp"

#include <algorithm>
#include <array>
#include <vector>

namespace NiTCAD::app::modeller {

namespace {

constexpr std::array<std::string_view, group_count> names{"Solid", "Edges", "Faces", "Transform", "Doping", "Mesh"};

constexpr std::array<ToolEntry, 4> solid_tools{{
    {"Select (V)", Action::tool, 0, "Click a shape to pick it; drag it to move it.", "Click a solid to pick it; drag it to move it."},
    {"Box (B)", Action::tool, 1, "Drag the rectangle; type a size for it exactly.",
     "Drag the base on the work plane, release, move for the height, click. Type a size for it exactly."},
    {"Prism (P)", Action::tool, 2, "Click the polygon's corners; Enter or the first corner closes it.",
     "Click the polygon's corners on the work plane; Enter closes it; then its height, a click ending it."},
    {"Fill / Replace", Action::overlap, 0, "The picked shapes' regions fill only empty space, or replace again.",
     "The picked solids' regions fill only empty space, or replace again."},
}};
constexpr std::array<ToolEntry, 2> face_tools{{
    {"Push/Pull", Action::tool, 0, "Drag a side's handle along its normal; type the distance.",
     "Drag a face's handle along its normal; type the distance."},
    {"Contact (C)", Action::tool, 3, "Drag along a side of the domain.", "Drag a rectangle on a face of the domain."},
}};
constexpr std::array<ToolEntry, 5> transform_tools{{
    {"Move", Action::tool, 0, "Drag the shape; type the distance along the larger offset.",
     "Drag the solid, or the gizmo's arrow or square; type the distance."},
    {"Copy", Action::mode, static_cast<int>(Mode::copy), "Drag a copy of the picked shapes to its place (or Alt+drag).",
     "Drag a copy of the picked solids to its place (or Alt+drag)."},
    {"Mirror", Action::mode, static_cast<int>(Mode::mirror),
     "Click a shape: mirrored left-right about its centre; Ctrl: up-down; on a side's handle: across that side; Shift: a quarter turn.",
     "Click a face: mirrored about the solid's centre across it; Ctrl: across the face's own plane; Shift: a quarter turn about the view's axis."},
    {"Array", Action::mode, static_cast<int>(Mode::array), "Drag the pitch from the shape; the wheel sets the count.",
     "Drag the pitch from the solid; the wheel sets the count."},
    {"Align", Action::mode, static_cast<int>(Mode::align), "Click a side's handle: the other picked shapes' matching sides move to its line.",
     "Click a face: the other picked solids' matching faces move to its plane."},
}};

}  // namespace

std::string_view group_name(Group g) noexcept { return names[static_cast<std::size_t>(g)]; }

Group group_for(Picked p) noexcept {
    switch (p) {
    case Picked::solid: return Group::transform;
    case Picked::edge: return Group::edges;
    case Picked::face:
    case Picked::contact: return Group::faces;
    case Picked::nothing: break;
    }
    return Group::solid;
}

std::span<const ToolEntry> tools(Group g) noexcept {
    switch (g) {
    case Group::solid: return solid_tools;
    case Group::faces: return face_tools;
    case Group::transform: return transform_tools;
    default: return {};
    }
}

std::string_view later(Group g) noexcept {
    switch (g) {
    case Group::solid: return "Cylinder comes in slice 2.";
    case Group::edges: return "Chamfer, Fillet, Round vertex and Taper come in slice 2.";
    case Group::faces: return "Inset and Extrude come in slice 2.";
    case Group::transform: return "";
    case Group::doping: return "Region profile, Implant and Junction come in slice 4.";
    case Group::mesh: return "Max size, Interface layers and Refine box come in slice 4.";
    }
    return "";
}

Nav navigation(bool fixed, bool three_d, const Press& p) noexcept {
    switch (p.button) {
    case Button::middle: return fixed ? Nav::none : Nav::pan;
    case Button::right: return fixed || !three_d ? Nav::none : Nav::orbit;
    case Button::left: break;
    }
    if (p.pan_tool) return fixed ? Nav::none : Nav::pan;
    if (fixed || !p.select_tool || !p.background) return Nav::model;
    if (three_d) return p.shift ? Nav::pan : Nav::orbit;
    return p.shift ? Nav::model : Nav::pan;  // 2D: Shift keeps the rubber band (MQ3)
}

Placement place_bar(const Box& viewport, const std::optional<Box>& selection, double width, double height, double gap) noexcept {
    const Box corner{viewport.left + gap, viewport.top + gap, viewport.left + gap + width, viewport.top + gap + height};
    if (!selection) return {corner, true};
    const Box& s = *selection;
    const auto slide = [](double lo, double size, double min, double max) { return std::clamp(lo, min, std::max(min, max - size)); };
    const double y = slide(s.top, height, viewport.top + gap, viewport.bottom - gap);
    const double x = slide(s.left, width, viewport.left + gap, viewport.right - gap);
    const std::array<Box, 4> candidates{{
        {s.right + gap, y, s.right + gap + width, y + height},     // right
        {x, s.bottom + gap, x + width, s.bottom + gap + height},   // below
        {s.left - gap - width, y, s.left - gap, y + height},       // left
        {x, s.top - gap - height, x + width, s.top - gap},         // above
    }};
    for (const Box& b : candidates) {
        if (b.inside(viewport) && !b.overlaps(s)) return {b, false};
    }
    return {corner, true};
}

namespace {

std::vector<std::string_view> lines_of(std::string_view text) {
    std::vector<std::string_view> out;
    std::size_t at = 0;
    while (at < text.size()) {
        std::size_t end = text.find('\n', at);
        if (end == std::string_view::npos) end = text.size();
        std::string_view line = text.substr(at, end - at);
        if (line.ends_with('\r')) line.remove_suffix(1);
        out.push_back(line);
        at = end + 1;
    }
    return out;
}

std::string joined(std::span<const std::string_view> lines) {
    std::string out;
    for (const std::string_view l : lines) {
        if (l.find_first_not_of(" \t") == std::string_view::npos) continue;
        if (!out.empty()) out += "  \xC2\xB7  ";
        out += l;
    }
    return out;
}

}  // namespace

std::string deck_line(std::string_view before, std::string_view after) {
    const auto a = lines_of(before), b = lines_of(after);
    std::size_t head = 0;
    while (head < a.size() && head < b.size() && a[head] == b[head]) ++head;
    std::size_t tail = 0;
    while (tail < a.size() - head && tail < b.size() - head && a[a.size() - 1 - tail] == b[b.size() - 1 - tail]) ++tail;
    const std::span<const std::string_view> added(b.data() + head, b.size() - head - tail);
    if (std::string out = joined(added); !out.empty()) return out;
    const std::span<const std::string_view> removed(a.data() + head, a.size() - head - tail);
    if (std::string out = joined(removed); !out.empty()) return "removed: " + out;
    return "";
}

bool Entry::type(wchar_t c) {
    if ((c >= L'0' && c <= L'9') || c == L'.' || c == L'-') {
        text_ += static_cast<char>(c);
        return true;
    }
    if (c == L'n' || c == L'u' || c == L'm' || c == L'c') {
        if (text_.empty()) return false;
        text_ += static_cast<char>(c);
        return true;
    }
    return false;
}

bool Entry::backspace() {
    if (text_.empty()) return false;
    text_.pop_back();
    return true;
}

std::optional<geo::Coord> Entry::value() const {
    std::string s = text_;
    if (s.empty()) return std::nullopt;
    if (!(s.ends_with("um") || s.ends_with("nm") || s.ends_with("cm"))) {
        if (s.back() < '0' || s.back() > '9') {
            if (s.back() != '.') return std::nullopt;  // a unit half typed
        }
        s += "um";
    }
    return geo::parse_length(s);
}

std::string chip_length(geo::Coord q) {
    std::string s = geo::format_length(q);
    if (q == 0) return "0 \xC2\xB5m";
    if (s.ends_with("um")) s.replace(s.size() - 2, 2, " \xC2\xB5m");
    else if (s.ends_with("nm") || s.ends_with("cm")) s.insert(s.size() - 2, " ");
    return s;
}

}  // namespace NiTCAD::app::modeller
