// The 3D modeller on the Geometry tab (solid_page.hpp; DECISIONS.md Unit 26l, SD4).
#include "solid_page.hpp"

#include <commctrl.h>
#include <windowsx.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <format>
#include <numbers>
#include <sstream>

#include "NiTCAD/app/deck_vocabulary.hpp"
#include "NiTCAD/render/axes3d.hpp"
#include "win32.hpp"

namespace NiTCAD::app {

namespace {

using namespace solid;

constexpr const wchar_t* canvas_class = L"NiTCADSolidCanvas";
constexpr const wchar_t* panel_class = L"NiTCADSolidPanel";
constexpr int id_strip = 1, id_canvas = 2, id_bar = 100;
constexpr std::array<const char*, 3> plane_names{"YZ (normal x)", "ZX (normal y)", "XY (normal z)"};

HWND make(const wchar_t* cls, const wchar_t* text, DWORD style, HWND parent, int id, DWORD ex = 0) {
    return CreateWindowExW(ex, cls, text, style | WS_CHILD | WS_VISIBLE, 0, 0, 10, 10, parent,
                           reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)), GetModuleHandleW(nullptr), nullptr);
}

std::string ui_length(Coord q) {
    std::string s = geo::format_length(q);
    if (s.ends_with("um")) s.resize(s.size() - 2);
    return s;
}

std::optional<Coord> ui_parse(std::string s) {
    while (!s.empty() && (s.back() == ' ' || s.back() == '\t' || s.back() == '\r')) s.pop_back();
    if (s.empty()) return std::nullopt;
    if (!(s.ends_with("um") || s.ends_with("nm") || s.ends_with("cm"))) s += "um";
    return geo::parse_length(s);
}

std::vector<std::string> choices_of(std::string_view section, std::string_view key) {
    std::vector<std::string> out;
    if (const SectionInfo* s = find_section_info(section)) {
        if (const KeyInfo* k = find_key_info(*s, key, "")) {
            for (const auto c : k->choices) out.emplace_back(c);
        }
    }
    return out;
}

const Solid* solid_at(const Model& m, SolidRef s) {
    if (s.region < 0 || s.region >= static_cast<int>(m.regions.size())) return nullptr;
    const auto& r = m.regions[static_cast<std::size_t>(s.region)];
    if (s.solid < 0 || s.solid >= static_cast<int>(r.solids.size())) return nullptr;
    return &r.solids[static_cast<std::size_t>(s.solid)];
}

// The axis a face's normal runs along (move_face's numbering), or -1 for a slanted side.
int face_axis(const Solid& s, int face) {
    if (s.kind == SolidKind::box) return face / 2;
    if (face < 2) return s.axis;
    const std::size_t n = s.polygon.size(), k = static_cast<std::size_t>(face - 2);
    if (k >= n) return -1;
    const Pt a = s.polygon[k], b = s.polygon[(k + 1) % n];
    if (a.y == b.y) return (s.axis + 2) % 3;
    if (a.x == b.x) return (s.axis + 1) % 3;
    return -1;
}

std::string fresh_contact(const Model& m) {
    for (int n = 1;; ++n) {
        std::string name = "contact_" + std::to_string(n);
        if (std::ranges::none_of(m.contacts, [&](const Contact& c) { return c.name == name; })) return name;
    }
}

}  // namespace

SolidPage::SolidPage(HWND parent, Host& host, int id, int inspector_id) : host_(host), inspector_id_(inspector_id) {
    win32::register_class(panel_class, DefWindowProcW, GetSysColorBrush(COLOR_WINDOW));
    win32::register_class(canvas_class, canvas_proc, nullptr, LoadCursorW(nullptr, IDC_ARROW), CS_DBLCLKS);
    panel_ = CreateWindowExW(WS_EX_CONTROLPARENT, panel_class, L"", WS_CHILD | WS_CLIPCHILDREN, 0, 0, 10, 10, parent,
                             reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)), GetModuleHandleW(nullptr), nullptr);
    SetWindowLongPtrW(panel_, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(this));
    SetWindowLongPtrW(panel_, GWLP_WNDPROC, reinterpret_cast<LONG_PTR>(panel_proc));
    strip_ = CreateWindowExW(WS_EX_CONTROLPARENT, panel_class, L"", WS_CHILD | WS_VISIBLE | WS_CLIPCHILDREN, 0, 0, 10, 10, panel_,
                             reinterpret_cast<HMENU>(static_cast<INT_PTR>(id_strip)), GetModuleHandleW(nullptr), nullptr);
    SetWindowLongPtrW(strip_, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(this));
    SetWindowLongPtrW(strip_, GWLP_WNDPROC, reinterpret_cast<LONG_PTR>(panel_proc));
    canvas_ = CreateWindowExW(0, canvas_class, L"", WS_CHILD | WS_VISIBLE | WS_TABSTOP | WS_CLIPSIBLINGS, 0, 0, 10, 10, panel_,
                              reinterpret_cast<HMENU>(static_cast<INT_PTR>(id_canvas)), GetModuleHandleW(nullptr), this);
    bar_ = std::make_unique<ModellerBar>(strip_, panel_, id_bar, true, panel_proc, this);
    start_camera_ = v_.camera;

    inspector_ = CreateWindowExW(WS_EX_CONTROLPARENT, panel_class, L"", WS_CHILD | WS_CLIPCHILDREN | WS_VSCROLL, 0, 0, 10, 10, parent,
                                 reinterpret_cast<HMENU>(static_cast<INT_PTR>(inspector_id - 1)), GetModuleHandleW(nullptr), nullptr);
    SetWindowLongPtrW(inspector_, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(this));
    SetWindowLongPtrW(inspector_, GWLP_WNDPROC, reinterpret_cast<LONG_PTR>(inspector_proc));
    controls_.resize(c_count);
    const auto add = [&](int c, const wchar_t* cls, const wchar_t* text, DWORD style, DWORD ex = 0) {
        controls_[static_cast<std::size_t>(c)] = make(cls, text, style, inspector_, inspector_id + c, ex);
    };
    const DWORD edit = WS_TABSTOP | ES_AUTOHSCROLL, combo = WS_TABSTOP | WS_VSCROLL | CBS_DROPDOWNLIST;
    add(c_heading, L"STATIC", L"GEOMETRY (3D)", SS_LEFT);
    add(c_hint, L"STATIC", L"", SS_LEFT);
    add(c_name_label, L"STATIC", L"Name", SS_LEFT);
    add(c_name, L"EDIT", L"", edit, WS_EX_CLIENTEDGE);
    add(c_material_label, L"STATIC", L"Material", SS_LEFT);
    add(c_material, WC_COMBOBOXW, L"", combo);
    add(c_raise, L"BUTTON", L"Raise layer", WS_TABSTOP | BS_PUSHBUTTON);
    add(c_lower, L"BUTTON", L"Lower layer", WS_TABSTOP | BS_PUSHBUTTON);
    const wchar_t* range_labels[] = {L"x from (µm)", L"x to", L"y from", L"y to", L"z from", L"z to"};
    for (int k = 0; k < 6; ++k) {
        add(c_x0_label + 2 * k, L"STATIC", range_labels[k], SS_LEFT);
        add(c_x0 + 2 * k, L"EDIT", L"", edit, WS_EX_CLIENTEDGE);
    }
    add(c_axis_label, L"STATIC", L"Axis", SS_LEFT);
    add(c_axis, WC_COMBOBOXW, L"", combo);
    add(c_from_label, L"STATIC", L"From (µm)", SS_LEFT);
    add(c_from, L"EDIT", L"", edit, WS_EX_CLIENTEDGE);
    add(c_to_label, L"STATIC", L"To (µm)", SS_LEFT);
    add(c_to, L"EDIT", L"", edit, WS_EX_CLIENTEDGE);
    add(c_vertices_label, L"STATIC", L"Polygon (u v, µm)", SS_LEFT);
    add(c_vertices, L"EDIT", L"", WS_TABSTOP | WS_VSCROLL | ES_MULTILINE | ES_AUTOVSCROLL | ES_WANTRETURN, WS_EX_CLIENTEDGE);
    add(c_apply, L"BUTTON", L"Apply", WS_TABSTOP | BS_PUSHBUTTON);
    add(c_kind_label, L"STATIC", L"Kind", SS_LEFT);
    add(c_kind, WC_COMBOBOXW, L"", combo);
    add(c_draw_label, L"STATIC", L"Draw with", SS_LEFT);
    add(c_draw, WC_COMBOBOXW, L"", combo);
    add(c_plane_label, L"STATIC", L"Work plane", SS_LEFT);
    add(c_plane, WC_COMBOBOXW, L"", combo);
    add(c_offset_label, L"STATIC", L"At (µm)", SS_LEFT);
    add(c_offset, L"EDIT", L"", edit, WS_EX_CLIENTEDGE);
    add(c_section_label, L"STATIC", L"Section", SS_LEFT);
    add(c_section, WC_COMBOBOXW, L"", combo);
    add(c_cut_label, L"STATIC", L"Cut at (µm)", SS_LEFT);
    add(c_cut, L"EDIT", L"", edit, WS_EX_CLIENTEDGE);
    add(c_translucent, L"BUTTON", L"Translucent (T)", WS_TABSTOP | BS_AUTOCHECKBOX);
    add(c_convert, L"BUTTON", L"Convert to layered generate3d", WS_TABSTOP | BS_PUSHBUTTON);
    add(c_overlap, L"BUTTON", L"Fills only empty space (overlap = fill)", WS_TABSTOP | BS_AUTOCHECKBOX);
    add(c_resolution, L"BUTTON", L"Use the largest interface_resolution", WS_TABSTOP | BS_PUSHBUTTON);  // MT12
    add(c_params_label, L"STATIC", L"Parameters (name = value; $name in solids)", SS_LEFT);  // MT8
    add(c_params, L"EDIT", L"", WS_TABSTOP | WS_VSCROLL | ES_MULTILINE | ES_AUTOVSCROLL | ES_WANTRETURN, WS_EX_CLIENTEDGE);
    add(c_params_apply, L"BUTTON", L"Apply parameters", WS_TABSTOP | BS_PUSHBUTTON);
    for (const int c : {c_material, c_draw}) {
        for (const std::string& mat : choices_of("region", "material")) SendMessageW(control(c), CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(win32::widen(mat).c_str()));
    }
    for (const std::string& k : choices_of("contact", "kind")) SendMessageW(control(c_kind), CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(win32::widen(k).c_str()));
    for (const wchar_t* a : {L"x", L"y", L"z"}) SendMessageW(control(c_axis), CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(a));
    for (const char* p : plane_names) SendMessageW(control(c_plane), CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(win32::widen(p).c_str()));
    for (const wchar_t* s : {L"Off", L"Cut above x", L"Cut above y", L"Cut above z"}) SendMessageW(control(c_section), CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(s));
    SendMessageW(control(c_draw), CB_SELECTSTRING, static_cast<WPARAM>(-1), reinterpret_cast<LPARAM>(L"silicon"));
    update_tool_buttons();
    update_inspector();
}

SolidPage::~SolidPage() {
    release();
    for (HWND h : {inspector_, panel_}) {
        if (h) {
            SetWindowLongPtrW(h, GWLP_USERDATA, 0);
            DestroyWindow(h);
        }
    }
}

void SolidPage::release() { target_.reset(); }

HWND SolidPage::control(int c) const { return c >= 0 && c < c_count ? controls_[static_cast<std::size_t>(c)] : nullptr; }

void SolidPage::set_dpi(int dpi, HFONT font, HFONT bold) {
    dpi_ = dpi;
    v_.dpi = dpi;
    bar_->set_font(font, dpi);
    for (HWND h : controls_) SendMessageW(h, WM_SETFONT, reinterpret_cast<WPARAM>(font), TRUE);
    SendMessageW(control(c_heading), WM_SETFONT, reinterpret_cast<WPARAM>(bold), TRUE);
    layout();
    layout_inspector();
}

void SolidPage::layout() {
    RECT r{};
    GetClientRect(panel_, &r);
    const int pad = MulDiv(4, dpi_, 96), bh = MulDiv(24, dpi_, 96);
    MoveWindow(strip_, 0, 0, r.right, bh + 2 * pad, TRUE);
    bar_->layout_strip(r.right, pad, bh);
    canvas_top_ = bh + 2 * pad;
    MoveWindow(canvas_, 0, canvas_top_, r.right, std::max<int>(1, r.bottom - canvas_top_), TRUE);
    bar_->raise();
    place_bar();
    invalidate();
}

void SolidPage::layout_inspector() {
    RECT r{};
    GetClientRect(inspector_, &r);
    const int pad = MulDiv(8, dpi_, 96), row = MulDiv(24, dpi_, 96), label = MulDiv(84, dpi_, 96);
    const int w = std::max<int>(MulDiv(80, dpi_, 96), r.right - 2 * pad);
    int y = pad;
    const auto visible = [&](int c) { return (GetWindowLongPtrW(control(c), GWL_STYLE) & WS_VISIBLE) != 0; };
    const auto line = [&](int c, int height = 0) {
        if (!visible(c)) return;
        MoveWindow(control(c), pad, y, w, height ? height : row - 4, TRUE);
        y += height ? height + 4 : row;
    };
    const auto pair = [&](int l, int c, int extra = 0) {
        if (!visible(c)) return;
        MoveWindow(control(l), pad, y + MulDiv(3, dpi_, 96), label - 4, row - 6, TRUE);
        MoveWindow(control(c), pad + label, y, w - label, row - 4 + extra, TRUE);
        y += row;
    };
    line(c_heading);
    line(c_hint, 5 * (row - 6));
    pair(c_name_label, c_name);
    pair(c_material_label, c_material, MulDiv(200, dpi_, 96));
    if (visible(c_raise)) {
        MoveWindow(control(c_raise), pad, y, w / 2 - 2, row - 4, TRUE);
        MoveWindow(control(c_lower), pad + w / 2 + 2, y, w - w / 2 - 2, row - 4, TRUE);
        y += row;
    }
    line(c_overlap);
    for (int k = 0; k < 6; ++k) pair(c_x0_label + 2 * k, c_x0 + 2 * k);
    pair(c_axis_label, c_axis, MulDiv(80, dpi_, 96));
    pair(c_from_label, c_from);
    pair(c_to_label, c_to);
    line(c_vertices_label);
    line(c_vertices, MulDiv(120, dpi_, 96));
    line(c_apply);
    pair(c_kind_label, c_kind, MulDiv(160, dpi_, 96));
    line(c_convert);
    line(c_resolution);
    line(c_params_label);
    line(c_params, MulDiv(80, dpi_, 96));
    line(c_params_apply);
    y += pad;
    pair(c_draw_label, c_draw, MulDiv(200, dpi_, 96));
    pair(c_plane_label, c_plane, MulDiv(90, dpi_, 96));
    pair(c_offset_label, c_offset);
    pair(c_section_label, c_section, MulDiv(110, dpi_, 96));
    pair(c_cut_label, c_cut);
    line(c_translucent);
}

// --- the deck ------------------------------------------------------------------------------------------------------

void SolidPage::reload(bool fit_view) {
    const std::string text = host_.solid_deck();
    if (text == text_seen_ && !fit_view && m_.mesh != MeshKind::unsupported) {
        update_inspector();
        invalidate();
        return;
    }
    if (text != text_seen_) history_.clear();  // edited elsewhere (or a new deck): the page's undo starts again
    const bool first = text_seen_.empty();
    text_seen_ = text;
    read_mesh_keys();
    auto m = read_model(text);
    if (!m) {
        m_ = {};
        message_ = "The deck does not parse: " + m.error();
    } else {
        m_ = std::move(*m);
        message_ = m_.mesh == MeshKind::unsupported ? m_.why : std::string();
    }
    if (first || fit_view) {
        // The work plane on the domain's floor.
        v_.plane = {2, m_.domain.lo.z};
        fit();
    }
    prune_selection();
    locate_problems();
    update_inspector();
    invalidate();
}

void SolidPage::prune_selection() {
    std::erase_if(sel_.solids, [&](SolidRef s) {
        if (s.region < 0 || s.region >= static_cast<int>(m_.regions.size())) return true;
        const auto& r = m_.regions[static_cast<std::size_t>(s.region)];
        return s.solid == -1 ? !r.solids.empty() : s.solid >= static_cast<int>(r.solids.size());
    });
    if (sel_.contact >= static_cast<int>(m_.contacts.size())) sel_.contact = -1;
}

void SolidPage::apply_text(const std::string& text) {
    host_.solid_edited(text);
    text_seen_ = text;
    read_mesh_keys();
    if (auto m = read_model(text)) {
        m_ = std::move(*m);
        message_ = m_.mesh == MeshKind::unsupported ? m_.why : std::string();
    }
    prune_selection();
    locate_problems();
    update_inspector();
    invalidate();
}

std::string SolidPage::parameters_text() const {
    std::string out;
    for (const solid::Param& p : m_.parameters) out += (out.empty() ? "" : "\r\n") + p.name + " = " + p.expr;
    return out;
}

bool SolidPage::apply_parameters(const std::string& lines) {
    if (!m_.editable()) return false;
    Model next = m_;
    std::istringstream in(lines);
    for (std::string line; std::getline(in, line);) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const auto eq = line.find('=');
        const auto trim = [](std::string s) {
            const auto a = s.find_first_not_of(" \t"), b = s.find_last_not_of(" \t");
            return a == std::string::npos ? std::string{} : s.substr(a, b - a + 1);
        };
        if (trim(line).empty()) continue;
        if (eq == std::string::npos) {
            host_.solid_status("Parameters: '" + line + "' is not name = value");
            return false;
        }
        const std::string name = trim(line.substr(0, eq)), expr = trim(line.substr(eq + 1));
        const auto it = std::ranges::find(next.parameters, name, &solid::Param::name);
        if (it != next.parameters.end() && it->expr == expr) continue;
        if (const auto r = solid::set_parameter(next, name, expr); !r) {
            host_.solid_status("Parameters: " + r.error());
            return false;
        }
    }
    return commit_model(std::move(next));
}

bool SolidPage::commit_model(Model next) {
    const std::string text = host_.solid_deck();
    // MT8: a coordinate the edit moved that a parameter ties edits the parameter (every tied face following).
    if (const auto was = read_model(text); was && was->parameters == next.parameters) tie_edits(*was, next);
    auto out = write_model(text, next);
    if (!out) {
        host_.solid_status(out.error());
        return false;
    }
    if (*out == text) {
        m_ = std::move(next);
        invalidate();
        return false;
    }
    history_.record(text);
    apply_text(*out);
    return true;
}

bool SolidPage::undo() {
    cancel();
    auto t = history_.undo(host_.solid_deck());
    if (!t) return false;
    apply_text(*t);
    return true;
}

bool SolidPage::redo() {
    cancel();
    auto t = history_.redo(host_.solid_deck());
    if (!t) return false;
    apply_text(*t);
    return true;
}

bool SolidPage::duplicate() {
    if (sel_.solids.size() != 1 || sel_.solids[0].solid < 0) return false;
    Model next = m_;
    const Coord step = grid();
    const int r = duplicate_solid(next, sel_.solids[0], {step, step, 0});
    if (r < 0 || !commit_model(std::move(next))) return false;
    sel_ = {{{r, 0}}, -1};
    update_inspector();
    return true;
}

bool SolidPage::remove() {
    Model next = m_;
    bool any = false;
    if (sel_.contact >= 0) {
        any = delete_contact(next, sel_.contact);
    } else {
        auto s = sel_.solids;
        std::ranges::sort(s, [](SolidRef a, SolidRef b) { return std::tie(a.region, a.solid) > std::tie(b.region, b.solid); });
        for (const SolidRef r : s) {
            if (r.solid >= 0) any = delete_solid(next, r) || any;
        }
    }
    sel_ = {};
    if (!any) return false;
    const bool done = commit_model(std::move(next));
    update_inspector();
    return done;
}

bool SolidPage::convert() {
    if (m_.mesh != MeshKind::extruded) return false;
    Model next = m_;
    const auto note = convert_to_layered(next);
    if (!note) {
        host_.solid_status("The extrusion cannot become layered: " + note.error());
        return false;
    }
    const bool done = commit_model(std::move(next));
    if (done) host_.solid_status("Converted: " + *note);
    return done;
}

bool SolidPage::select_region(const std::string& name) {
    for (std::size_t r = 0; r < m_.regions.size(); ++r) {
        if (m_.regions[r].name != name) continue;
        sel_ = {{{static_cast<int>(r), m_.regions[r].solids.empty() ? -1 : 0}}, -1};
        update_inspector();
        invalidate();
        return true;
    }
    return false;
}

// --- the view ------------------------------------------------------------------------------------------------------

void SolidPage::set_camera(const render::Camera& c) {
    v_.camera = c;
    invalidate();
}

void SolidPage::fit() {
    v_.camera.zoom = 1.0f;
    v_.camera.pan = {0.0f, 0.0f};
    invalidate();
}

void SolidPage::named_view(int index) {
    const auto& views = render::named_views();
    if (index < 0 || index >= static_cast<int>(views.size())) return;
    v_.camera.azimuth_deg = views[static_cast<std::size_t>(index)].azimuth_deg;
    v_.camera.elevation_deg = views[static_cast<std::size_t>(index)].elevation_deg;
    invalidate();
}

void SolidPage::set_work_plane(WorkPlane p) {
    v_.plane = p;
    update_inspector();
    invalidate();
}

void SolidPage::set_translucent(bool on) {
    v_.translucent = on;
    SendMessageW(control(c_translucent), BM_SETCHECK, on ? BST_CHECKED : BST_UNCHECKED, 0);
    invalidate();
}

void SolidPage::set_section(Section s) {
    v_.section = s;
    update_inspector();
    invalidate();
}

void SolidPage::set_tool(Tool t) {
    cancel();
    tool_ = t;
    mode_ = modeller::Mode::none;
    update_tool_buttons();
    invalidate();
}

void SolidPage::update_tool_buttons() {
    // A tool armed by its key opens a group that holds it (Pan and Measure are in none: the group stays).
    const auto in = [&](modeller::Group g) {
        return std::ranges::any_of(modeller::tools_for(g, true), [&](const modeller::ToolEntry& e) {
            return e.action == modeller::Action::tool && e.tool == static_cast<int>(tool_);
        });
    };
    if (!in(bar_->group())) {
        for (int g = 0; g < modeller::group_count; ++g) {
            if (in(static_cast<modeller::Group>(g))) {
                bar_->set_group(static_cast<modeller::Group>(g));
                break;
            }
        }
    }
    bar_->set_tool(modeller::armed_key(static_cast<int>(tool_), mode_));
    bar_->set_measure(tool_ == Tool::measure);
}

// --- Unit 26m slice 0: the camera lock, the floating bar, the deck line and the chip (MT1-MT3) ---------------------

void SolidPage::set_camera_fixed(bool fixed) {
    camera_fixed_ = fixed;
    bar_->set_camera_fixed(fixed);
    if (fixed && (drag_ == Drag::orbit || drag_ == Drag::pan)) cancel();
}

void SolidPage::reset_view() {
    v_.camera = start_camera_;
    invalidate();
}

void SolidPage::set_group(modeller::Group g) {
    bar_->set_group(g);
    if (mode_ != modeller::Mode::none && modeller::group_of(mode_) != g) mode_ = modeller::Mode::none;
    const auto list = modeller::tools_for(g, true);
    const bool has = std::ranges::any_of(list, [&](const modeller::ToolEntry& e) { return e.action == modeller::Action::tool && e.tool == static_cast<int>(tool_); });
    if (!has && !list.empty() && list.front().action == modeller::Action::tool) {
        set_tool(static_cast<Tool>(list.front().tool));
    } else {
        bar_->set_tool(modeller::armed_key(static_cast<int>(tool_), mode_));
    }
    place_bar();
}


// --- Unit 26m slice 1: the overlap rule and the Transform modes (MT4, MT5) -----------------------------------------

void SolidPage::set_mode(modeller::Mode m) {
    if (tool_ != Tool::select) set_tool(Tool::select);
    cancel();
    mode_ = m;
    array_count_ = 3;
    if (m != modeller::Mode::none && bar_->group() != modeller::group_of(m)) bar_->set_group(modeller::group_of(m));
    bar_->set_tool(modeller::armed_key(static_cast<int>(tool_), mode_));
    place_bar();
    invalidate();
}

bool SolidPage::toggle_overlap() {
    std::vector<int> regions;
    for (const SolidRef s : sel_.solids) {
        if (s.solid >= 0 && std::ranges::find(regions, s.region) == regions.end()) regions.push_back(s.region);
    }
    if (regions.empty()) {
        host_.solid_status("Pick a solid: its region then fills only empty space, or replaces again.");
        return false;
    }
    Model next = m_;
    // One rule for them all: fill unless every one fills already.
    const bool fill = !std::ranges::all_of(regions, [&](int r) { return m_.regions[static_cast<std::size_t>(r)].fill; });
    bool any = false;
    for (const int r : regions) any = set_overlap(next, r, fill) || any;
    if (!any || !commit_model(std::move(next))) return false;
    update_inspector();
    return true;
}

int SolidPage::view_axis() const {
    // The axis whose unit step moves least on the screen runs most nearly along the view.
    const Box3 b = m_.domain;
    const V3 c{0.5 * (um_of(b.lo)[0] + um_of(b.hi)[0]), 0.5 * (um_of(b.lo)[1] + um_of(b.hi)[1]), 0.5 * (um_of(b.lo)[2] + um_of(b.hi)[2])};
    const auto p0 = screen_of_um(c);
    int best = 2;
    double shortest = std::numeric_limits<double>::infinity();
    for (int a = 0; a < 3; ++a) {
        V3 q = c;
        q[static_cast<std::size_t>(a)] += 1.0;
        const auto p1 = screen_of_um(q);
        const double d = std::hypot(p1[0] - p0[0], p1[1] - p0[1]);
        if (d < shortest) {
            shortest = d;
            best = a;
        }
    }
    return best;
}

bool SolidPage::transform_click(const Pick& p, bool shift, bool ctrl) {
    const bool on_solid = (p.kind == PickKind::face || (p.kind == PickKind::handle && p.handle == HandleKind::face)) && p.solid.solid >= 0;
    const Solid* target = on_solid ? solid_at(m_, p.solid) : nullptr;
    // The plane of the face clicked: its axis, its coordinate and whether it is the solid's far side on that axis.
    int axis = -1;
    Coord at = 0;
    bool far_side = false;
    if (target) {
        axis = face_axis(*target, p.index);
        if (axis >= 0) {
            at = quanta_of(p.at)[axis];
            const Box3 b = target->bounds();
            far_side = 2 * at >= b.lo[axis] + b.hi[axis];
        }
    } else if (p.domain_face >= 0) {
        axis = p.domain_face / 2;
        far_side = p.domain_face % 2 == 1;
        at = far_side ? m_.domain.hi[axis] : m_.domain.lo[axis];
    }
    Model next = m_;
    bool ok = false;
    if (mode_ == modeller::Mode::mirror) {
        if (!target) {
            host_.solid_status("Mirror: click a face of a solid.");
            return false;
        }
        if (shift) {
            ok = rotate_solid(next, p.solid, view_axis());
        } else if (axis < 0) {
            host_.solid_status("Mirror: a slanted face has no plane on the lattice; click another face.");
            return false;
        } else {
            const Box3 b = target->bounds();
            ok = mirror_solid_at(next, p.solid, axis, ctrl ? 2 * at : b.lo[axis] + b.hi[axis]);
        }
        if (!ok) {
            host_.solid_status("That mirror leaves the solid invalid.");
            return false;
        }
        sel_ = {{p.solid}, -1};
        const bool done = commit_model(std::move(next));  // false too when nothing changes (a box about its centre)
        update_inspector();
        return done;
    }
    // Align: the picked solids other than the one clicked.
    if (axis < 0) {
        host_.solid_status("Align: click a face (of a solid or the domain) whose plane the others' faces go to.");
        return false;
    }
    std::vector<SolidRef> others;
    for (const SolidRef s : sel_.solids) {
        if (s.solid >= 0 && !(on_solid && s == p.solid)) others.push_back(s);
    }
    if (others.empty()) {
        host_.solid_status("Align: pick the solids to move first (Shift+click), then click the face to align them to.");
        return false;
    }
    ok = align_solids(next, others, 2 * axis + (far_side ? 1 : 0), at);
    if (ok && commit_model(std::move(next))) {
        update_inspector();
        return true;
    }
    host_.solid_status("Align: nothing moves (the faces are there already), or a solid would leave the limits.");
    return false;
}
// --- Unit 26m slice 2: the Faces and Edges modes, Cylinder, meshable by construction (MT6, MT7, MT12) --------------

void SolidPage::begin_drag(Drag d) {
    drag_ = d;
    button_ = true;
    m_start_ = m_;
    SetCapture(canvas_);
}

void SolidPage::read_mesh_keys() {
    const auto s = mesh_length(text_seen_, "interface_resolution");
    res_q_ = s ? static_cast<Coord>(std::llround(*s * static_cast<double>(per_um) * 1e4)) : 0;
    max_size_ = mesh_length(text_seen_, "max_size").value_or(0.0);
    profiles_ = solid::profiles(text_seen_);  // slice 4: drawn each frame, read once per deck
    size_boxes_ = solid::size_boxes(text_seen_);
}

Coord SolidPage::gap_up(Coord v) const {
    const Coord g = 4 * res_q_;
    return v != 0 && std::abs(v) < g ? (v < 0 ? -g : g) : v;
}

std::expected<double, std::string> SolidPage::automatic_resolution() const {
    if (!m_.editable()) return std::unexpected(std::string("the deck is not a layered 3D one"));
    if (!(max_size_ > 0.0)) return std::unexpected(std::string("the deck has no max_size"));
    return solid::automatic_resolution(m_, max_size_);
}

bool SolidPage::use_automatic_resolution() {
    const auto s = automatic_resolution();
    if (!s) {
        host_.solid_status("No interface_resolution passes the rules: " + s.error());
        return false;
    }
    if (std::isinf(*s)) {
        host_.solid_status("The layout has no interfaces: any interface_resolution goes.");
        return false;
    }
    const Coord q = static_cast<Coord>(std::floor(*s * static_cast<double>(per_um) * 1e4));
    if (q <= 0) return false;
    const std::string text = host_.solid_deck();
    const auto out = set_mesh_length(text, "interface_resolution", q);
    const auto e = extent(m_);
    const std::string about = e ? std::format("; about {:.0f} nodes (MT12 estimate)", node_estimate(*e, max_size_, static_cast<double>(q) * 1e-8)) : std::string();
    host_.solid_status("interface_resolution = " + geo::format_length(q) + ", the largest the layout's rules allow" + about);
    if (!out || *out == text) return false;
    history_.record(text);
    apply_text(*out);
    update_inspector();
    return true;
}

// The edge cut of size r on the model the drag started from (unchanged when refused); chords_ the fillet's chords.
Model SolidPage::cut_model(Coord r) {
    Model next = m_start_;
    chords_ = 0;
    const Solid* s = solid_at(m_start_, ref_);
    if (!s || r <= 0) return next;
    const auto ec = edge_corner(*s, edge_);
    if (!ec) return next;
    if (mode_ == modeller::Mode::chamfer) {
        chamfer_edge(next, ref_, edge_, r);
        return next;
    }
    // MQ4, MT12 (2): ⌈θ / 15°⌉ chords, fewer while one is under the generator's 6 s.
    const double turn = geo::corner_turn(ec->prism.polygon, ec->corner);
    chords_ = geo::arc_segments(turn, static_cast<double>(r), 6.0 * static_cast<double>(res_q_));
    if (chords_ > 0) fillet_edge(next, ref_, edge_, r, chords_);
    return next;
}

bool SolidPage::tool_press(const Pick& p, double x, double y) {
    using modeller::Mode;
    if (mode_ == Mode::none || modeller::group_of(mode_) == modeller::Group::transform) return false;
    if (const auto g = modeller::group_of(mode_); g == modeller::Group::doping || g == modeller::Group::mesh) return doping_press(p, x, y);
    const bool on_face = (p.kind == PickKind::face || (p.kind == PickKind::handle && p.handle == HandleKind::face)) && p.solid.solid >= 0;
    if (mode_ == Mode::inset || mode_ == Mode::extrude) {
        if (!on_face) return false;
        const Solid* s = solid_at(m_, p.solid);
        const auto n = s ? face_normal(*s, p.index) : std::nullopt;
        if (!n) {
            host_.solid_status("A slanted face takes no inset or extrusion; press another face.");
            return true;
        }
        ref_ = p.solid;
        index_ = p.index;
        normal_ = *n;
        origin_ = p.at;  // a face's hit point, or a handle's own place on the face (pick sets it)
        sel_ = {{p.solid}, -1};
        update_inspector();
        if (mode_ == Mode::extrude) {
            start_s_ = along_axis(m_, v_, x, y, origin_, n->axis);
            begin_drag(Drag::extrude);
        } else {
            const auto q = on_plane(m_, v_, x, y, n->axis, origin_[static_cast<std::size_t>(n->axis)]);
            if (!q) return true;
            start_p_ = *q;
            margin_ = 0;
            begin_drag(Drag::inset);
        }
        return true;
    }
    // The Edges modes: an edge along a solid's axis (Round vertex: a prism's vertex handle); a press on a solid picks it,
    // so its edges can be picked, and moves nothing.
    const bool on_edge = p.kind == PickKind::edge && p.solid.solid >= 0 && mode_ != Mode::round;
    const bool on_vertex = p.kind == PickKind::handle && p.handle == HandleKind::vertex && mode_ == Mode::round;
    if (!on_edge && !on_vertex) {
        if (p.kind == PickKind::face && p.solid.solid >= 0) {
            sel_ = {{p.solid}, -1};
            update_inspector();
            invalidate();
            host_.solid_status(mode_ == Mode::round ? "Now drag in from one of its vertex handles." : "Now drag in from one of its edges along its axis.");
            return true;
        }
        return false;
    }
    const Solid* s = solid_at(m_, p.solid);
    if (!s) return true;
    edge_ = on_vertex ? 3 * p.index + 2 : p.index;
    const auto ec = edge_corner(*s, edge_);
    if (!ec) {
        host_.solid_status("That edge is a cap's: chamfers, fillets and tapers take edges along a solid's axis.");
        return true;
    }
    ref_ = p.solid;
    axis_ = ec->prism.axis;
    origin_ = um_of(from_uv(axis_, ec->prism.polygon[static_cast<std::size_t>(ec->corner)], quanta_of(p.at)[axis_]));
    const auto q = on_plane(m_, v_, x, y, axis_, origin_[static_cast<std::size_t>(axis_)]);
    if (!q) return true;
    start_p_ = *q;
    chords_ = 0;
    begin_drag(mode_ == Mode::taper ? Drag::taper : Drag::cut);
    return true;
}

void SolidPage::finish_inset() {
    Model done = m_;
    const bool changed = !(done == m_start_);
    m_ = m_start_;
    drag_ = Drag::none;
    button_ = false;
    end_drag_ui();
    if (GetCapture() == canvas_) ReleaseCapture();
    if (changed) commit_model(std::move(done));
    else host_.solid_status("The inset is not taken (it needs a depth, and a margin that leaves it a size).");
    update_inspector();
    invalidate();
}

std::optional<modeller::Box> SolidPage::selection_bounds() const {
    std::vector<Box3> boxes;
    for (const SolidRef s : sel_.solids) {
        if (const Solid* sol = solid_at(m_, s)) boxes.push_back(sol->bounds());
        else if (s.region >= 0) boxes.push_back(m_.domain);  // a region's background
    }
    if (sel_.contact >= 0 && sel_.contact < static_cast<int>(m_.contacts.size())) boxes.push_back(m_.contacts[static_cast<std::size_t>(sel_.contact)].box);
    if (boxes.empty()) return std::nullopt;
    modeller::Box out{1e300, 1e300, -1e300, -1e300};
    for (const Box3& b : boxes) {
        for (int k = 0; k < 8; ++k) {
            const P3 c{k & 1 ? b.hi.x : b.lo.x, k & 2 ? b.hi.y : b.lo.y, k & 4 ? b.hi.z : b.lo.z};
            const auto p = screen_of(c);
            out.left = std::min(out.left, p[0]);
            out.right = std::max(out.right, p[0]);
            out.top = std::min(out.top, p[1]);
            out.bottom = std::max(out.bottom, p[1]);
        }
    }
    return out;
}

void SolidPage::place_bar() {
    RECT c{};
    GetClientRect(canvas_, &c);
    const double top = canvas_top_;
    const modeller::Box viewport{0.0, top, static_cast<double>(c.right), top + c.bottom};
    std::optional<modeller::Box> s = selection_bounds();
    if (s) *s = {s->left, s->top + top, s->right, s->bottom + top};
    bar_->place(viewport, s);
}

modeller::Placement SolidPage::bar_placement() const {
    modeller::Placement p = bar_->placement();
    p.bar.top -= canvas_top_;
    p.bar.bottom -= canvas_top_;
    return p;
}

// The model the drag would write now (the in-progress model, or it with the solid or contact being drawn).
std::optional<Model> SolidPage::pending() const {
    switch (drag_) {
    case Drag::move:
    case Drag::array:
    case Drag::axis:
    case Drag::plane:
    case Drag::face:
    case Drag::cut:
    case Drag::taper:
    case Drag::extrude:
    case Drag::inset:
    case Drag::inset_depth:
    case Drag::vertex:
    case Drag::contact_move: return m_;
    case Drag::radius: {
        if (preview_.solids.empty()) return std::nullopt;
        Model next = m_;
        const Solid& s = preview_.solids.front();
        if (add_prism(next, s.axis, s.from, s.to, s.polygon, material_) < 0) return std::nullopt;
        return next;
    }
    case Drag::base: {
        // As released now: the height one grid step (on_up).
        Model next = m_;
        Box3 b = normalized(a_, b_);
        b.hi[v_.plane.normal] = b.lo[v_.plane.normal] + grid();
        if (add_box(next, b, material_) < 0) return std::nullopt;
        return next;
    }
    case Drag::height: {
        if (preview_.solids.empty()) return std::nullopt;
        Model next = m_;
        const Solid& s = preview_.solids.front();
        const int r = s.kind == SolidKind::box ? add_box(next, s.box, material_) : add_prism(next, s.axis, s.from, s.to, s.polygon, material_);
        if (r < 0) return std::nullopt;
        return next;
    }
    case Drag::contact: {
        Model next = m_;
        if (add_contact(next, fresh_contact(m_), "ohmic", normalized(a_, b_)) < 0) return std::nullopt;
        return next;
    }
    default: return std::nullopt;
    }
}

std::string SolidPage::chip_value() const {
    using modeller::chip_length;
    const char* names[] = {"x", "y", "z"};
    const auto offset = [&](P3 d) {
        std::string out;
        for (int a = 0; a < 3; ++a) {
            if (d[a] == 0) continue;
            if (!out.empty()) out += ", ";
            out += std::format("\xCE\x94{} {}", names[a], chip_length(d[a]));
        }
        return out.empty() ? std::string("\xCE\x94 0 \xC2\xB5m") : out;
    };
    std::string v;
    switch (drag_) {
    case Drag::axis:
    case Drag::move:
    case Drag::plane:
    case Drag::vertex:
    case Drag::contact_move: v = offset(drag_d_); break;
    case Drag::face:
        v = "push " + chip_length(drag_s_);
        if (naming_ >= 0) v += ", named $" + naming_name_ + (naming_ == 0 ? " (new; N: next)" : " (N: next)");
        break;
    case Drag::array: v = offset(drag_d_) + std::format(" × {} (wheel)", array_count_); break;
    case Drag::base:
    case Drag::contact: {
        const Box3 b = normalized(a_, b_);
        std::string out;
        for (int a = 0; a < 3; ++a) {
            if (b.hi[a] == b.lo[a]) continue;
            if (!out.empty()) out += " \xC3\x97 ";
            out += chip_length(b.hi[a] - b.lo[a]);
        }
        v = out.empty() ? "0 \xC2\xB5m" : out;
        break;
    }
    case Drag::height: v = "height " + chip_length(height_); break;
    case Drag::cut:
        v = (mode_ == modeller::Mode::chamfer ? "cut " : "radius ") + chip_length(drag_s_);
        if (mode_ != modeller::Mode::chamfer) v += chords_ > 0 ? std::format(", {} chords", chords_) : std::string(", chords under 6 s");
        break;
    case Drag::taper: {
        v = offset(drag_d_);
        if (const Solid* s = solid_at(m_start_, ref_)) {
            if (const auto ec = edge_corner(*s, edge_)) {
                // The wall leaning: the neighbour edge across the move, its length against the offset.
                const auto& w = ec->prism.polygon;
                const std::size_t n = w.size(), k = static_cast<std::size_t>(ec->corner);
                const Pt d = uv_of(axis_, drag_d_);
                for (const std::size_t j : {(k + n - 1) % n, (k + 1) % n}) {
                    const double ex = static_cast<double>(w[j].x - w[k].x), ey = static_cast<double>(w[j].y - w[k].y);
                    if (ex * static_cast<double>(d.x) + ey * static_cast<double>(d.y) == 0.0 && (d.x != 0 || d.y != 0)) {
                        v += std::format(", wall {:.1f}°", std::atan2(std::hypot(static_cast<double>(d.x), static_cast<double>(d.y)), std::hypot(ex, ey)) * 180.0 / std::numbers::pi);
                        break;
                    }
                }
            }
        }
        break;
    }
    case Drag::extrude: v = "depth " + chip_length(drag_s_); break;
    case Drag::implant:
    case Drag::junction:
        v = (drag_ == Drag::implant ? "peak depth " : "junction depth ") + chip_length(drag_s_) +
            std::format(", {} {:.0e} (D, wheel)", donor_ ? "donor" : "acceptor", level_);
        break;
    case Drag::max_size: {
        const auto e = extent(m_);
        v = "max_size " + chip_length(drag_s_) +
            (e && res_q_ > 0 ? std::format(", about {:.0f} nodes", node_estimate(*e, static_cast<double>(drag_s_) * 1e-8, static_cast<double>(res_q_) * 1e-8)) : std::string());
        break;
    }
    case Drag::refine: v = "size " + chip_length(refine_size_) + " (wheel)"; break;
    case Drag::inset: v = "margin " + chip_length(margin_); break;
    case Drag::inset_depth: v = "margin " + chip_length(margin_) + ", depth " + chip_length(drag_s_); break;
    case Drag::radius: v = "radius " + chip_length(radius_) + (preview_.solids.empty() ? std::string(", sides under 6 s") : std::format(", {} sides", preview_.solids.front().polygon.size())); break;
    case Drag::measure: {
        const double d = std::sqrt(std::pow(static_cast<double>(b_.x - a_.x), 2) + std::pow(static_cast<double>(b_.y - a_.y), 2) +
                                   std::pow(static_cast<double>(b_.z - a_.z), 2)) / static_cast<double>(per_um);
        v = std::format("{:.6g} \xC2\xB5m", d);
        break;
    }
    default: return {};
    }
    if (!entry_.empty()) v += "   [" + entry_.text() + (entry_.value() ? "]" : "?]");
    return v;
}

void SolidPage::after_drag_step(double x, double y) {
    if (drag_ == Drag::implant || drag_ == Drag::junction || drag_ == Drag::max_size || drag_ == Drag::refine) return;  // doping_step
    if (auto next = pending()) {
        // As commit_model writes it: a tied coordinate moved edits its parameter (MT8, MC4).
        if (m_start_.parameters == next->parameters && !next->parameters.empty()) tie_edits(m_start_, *next);
        auto out = write_model(text_seen_, *next);
        live_line_ = out ? modeller::deck_line(text_seen_, *out) : out.error();
        host_.solid_status(live_line_);
    } else {
        live_line_.clear();
    }
    const std::string v = chip_value();
    const int off = MulDiv(16, dpi_, 96);
    bar_->chip(static_cast<int>(x) + off, static_cast<int>(y) + canvas_top_ + off, v);
}

void SolidPage::end_drag_ui() {
    entry_.clear();
    typed_axis_ = -1;
    live_line_.clear();
    bar_->chip(0, 0, "");
}

bool SolidPage::on_char(wchar_t c) {
    const bool scalar = drag_ == Drag::axis || drag_ == Drag::face || drag_ == Drag::height || drag_ == Drag::cut || drag_ == Drag::taper ||
                        drag_ == Drag::extrude || drag_ == Drag::inset || drag_ == Drag::inset_depth || drag_ == Drag::radius ||
                        drag_ == Drag::implant || drag_ == Drag::junction || drag_ == Drag::max_size;
    const bool two = drag_ == Drag::move || drag_ == Drag::plane || drag_ == Drag::base || drag_ == Drag::array;
    // Slice 4 (MT9): D swaps the doping tools' dopant, during a drag too.
    if ((c == L'd' || c == L'D') && modeller::group_of(mode_) == modeller::Group::doping && mode_ != modeller::Mode::none) {
        donor_ = !donor_;
        host_.solid_status(std::format("{} {:.0e} cm^-3 (D: the dopant; the wheel: the level)", donor_ ? "donor" : "acceptor", level_));
        if (drag_ != Drag::none) on_move(last_x_, last_y_, false, false);
        return true;
    }
    if (!scalar && !two) return false;
    // MT8: N names the dragged extent (a new parameter, then each existing length one in turn), while nothing is
    // typed (once a number is, n is its unit's).
    if (drag_ == Drag::face && entry_.empty() && (c == L'n' || c == L'N')) {
        const Solid* s = solid_at(m_start_, ref_);
        const bool along = s && (s->kind == SolidKind::box || index_ < 2);  // a box's face or a prism's cap
        if (!along) return true;
        std::vector<std::string> lengths;
        for (const solid::Param& p : m_start_.parameters) {
            if (p.length) lengths.push_back(p.name);
        }
        naming_ = naming_ + 1 > static_cast<int>(lengths.size()) ? 0 : naming_ + 1;
        if (naming_ == 0) {
            for (int k = 1;; ++k) {
                naming_name_ = "L" + std::to_string(k);
                if (std::ranges::find(m_start_.parameters, naming_name_, &solid::Param::name) == m_start_.parameters.end()) break;
            }
        } else {
            naming_name_ = lengths[static_cast<std::size_t>(naming_ - 1)];
        }
        on_move(last_x_, last_y_, false, false);
        return true;
    }
    if (c == L'\b') {
        if (!entry_.backspace()) return true;
    } else if (!entry_.type(c)) {
        return false;
    }
    if (two && typed_axis_ < 0) {
        // The axis of the larger offset so far (the first in-plane axis when none).
        const P3 d = drag_ == Drag::base ? P3{b_.x - a_.x, b_.y - a_.y, b_.z - a_.z} : drag_d_;
        const int skip = drag_ == Drag::base ? v_.plane.normal : axis_;
        for (int a = 0; a < 3; ++a) {
            if (a == skip) continue;
            if (typed_axis_ < 0 || std::abs(d[a]) > std::abs(d[typed_axis_])) typed_axis_ = a;
        }
    }
    on_move(last_x_, last_y_, false, false);
    return true;
}

void SolidPage::finish_height() {
    Model next = m_;
    const int n = v_.plane.normal;
    const Coord h = height_ != 0 ? height_ : grid();
    const Coord lo = std::min(v_.plane.offset, v_.plane.offset + h), hi = std::max(v_.plane.offset, v_.plane.offset + h);
    int r = -1;
    if (polygon_.empty()) {
        Box3 b = normalized(a_, b_);
        b.lo[n] = lo;
        b.hi[n] = hi;
        r = add_box(next, b, material_);
    } else {
        std::vector<Pt> uv;
        for (const P3& p : polygon_) uv.push_back(uv_of(n, p));
        r = add_prism(next, n, lo, hi, uv, material_);
    }
    drag_ = Drag::none;
    polygon_.clear();
    preview_ = {};
    end_drag_ui();
    if (GetCapture() == canvas_) ReleaseCapture();
    if (r >= 0 && commit_model(std::move(next))) sel_ = {{{r, 0}}, -1};
    else host_.solid_status("That solid is not taken (it must have a size on each axis).");
    update_inspector();
    invalidate();
}

void SolidPage::apply_drag() {
    switch (drag_) {
    case Drag::height: finish_height(); return;
    case Drag::inset_depth: finish_inset(); return;
    case Drag::cut:
    case Drag::taper:
    case Drag::extrude:
    case Drag::inset:
    case Drag::radius:
    case Drag::move:
    case Drag::array:
    case Drag::axis:
    case Drag::plane:
    case Drag::face:
    case Drag::vertex:
    case Drag::contact_move:
    case Drag::base:
    case Drag::implant:
    case Drag::junction:
    case Drag::max_size:
    case Drag::refine:
    case Drag::contact: on_up(last_x_, last_y_); return;
    default: return;
    }
}

void SolidPage::cancel() {
    pending_text_.clear();
    pending_volume_.reset();
    if (drag_ == Drag::move || drag_ == Drag::axis || drag_ == Drag::plane || drag_ == Drag::face || drag_ == Drag::vertex ||
        drag_ == Drag::contact_move || drag_ == Drag::array || drag_ == Drag::cut || drag_ == Drag::taper || drag_ == Drag::extrude ||
        drag_ == Drag::inset || drag_ == Drag::inset_depth) {
        m_ = copy_from_ ? *copy_from_ : m_start_;
        if (copy_from_) sel_ = copy_sel_;
    }
    copy_from_.reset();
    drag_ = Drag::none;
    button_ = false;
    polygon_.clear();
    preview_ = {};
    end_drag_ui();
    if (GetCapture() == canvas_) ReleaseCapture();
    invalidate();
}

void SolidPage::invalidate() {
    if (canvas_) InvalidateRect(canvas_, nullptr, FALSE);
}

std::array<double, 2> SolidPage::screen_of_um(const V3& p) const { return solid::screen_of(m_, v_, p); }
std::array<double, 2> SolidPage::screen_of(P3 p) const { return screen_of_um(um_of(p)); }

Coord SolidPage::grid() const {
    const V3 c{0.5 * (um_of(m_.domain.lo)[0] + um_of(m_.domain.hi)[0]), 0.5 * (um_of(m_.domain.lo)[1] + um_of(m_.domain.hi)[1]),
               0.5 * (um_of(m_.domain.lo)[2] + um_of(m_.domain.hi)[2])};
    return view_grid(m_, v_, c);
}

Coord SolidPage::quantize(double um) const {
    const Coord step = grid();
    return static_cast<Coord>(std::llround(um * static_cast<double>(per_um) / static_cast<double>(step))) * step;
}

std::optional<P3> SolidPage::plane_point(double x, double y) const {
    const auto p = on_plane(m_, v_, x, y, v_.plane.normal, static_cast<double>(v_.plane.offset) / static_cast<double>(per_um));
    if (!p) return std::nullopt;
    P3 raw = quanta_of(*p);
    raw[v_.plane.normal] = v_.plane.offset;
    P3 out = snap(m_, v_.plane, raw, snap_tolerance(m_, v_, *p), std::nullopt, grid()).p;
    out[v_.plane.normal] = v_.plane.offset;
    return out;
}

double SolidPage::render() {
    RECT c{};
    GetClientRect(canvas_, &c);
    v_.width = std::max<int>(16, c.right);
    v_.height = std::max<int>(16, c.bottom);
    v_.grid = m_.mesh == MeshKind::unsupported ? 0 : grid();
    if (!renderer_) renderer_ = host_.solid_renderer();
    if (!renderer_ || c.right < 16 || c.bottom < 16) return -1.0;
    const auto t0 = std::chrono::steady_clock::now();
    solid::Overlay overlay;
    if (check() == Check::invalid) {
        for (const solid::Problem& p : problems_) overlay.markers.insert(overlay.markers.end(), p.markers.begin(), p.markers.end());
    }
    if (mesh_shown()) overlay.mesh = mesh_;
    // Slice 4: the profiles (Doping open) and refine boxes (Mesh open) as see-through volumes, the drag's too.
    if (bar_->group() == modeller::Group::doping) {
        for (const auto& p : profiles_) {
            overlay.volumes.push_back({solid::profile_extent(m_.domain, p), solid::doping_color(p.dopant == device::Dopant::donor, p.peak_cm3)});
        }
    } else if (bar_->group() == modeller::Group::mesh) {
        for (const auto& b : size_boxes_) overlay.volumes.push_back({b.box, {0.55f, 0.55f, 0.60f, 0.22f}});
    }
    if (pending_volume_) overlay.volumes.push_back(*pending_volume_);
    scene_ = solid_scene(m_, v_, sel_, preview_, overlay);
    scene_ms_ = 1e3 * std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    if (!target_) {
        auto t = renderer_->window_target(canvas_, v_.width, v_.height);
        if (!t) {
            error_ = t.error().message;
            return -1.0;
        }
        target_ = std::move(*t);
    }
    if (auto r = renderer_->present(scene_, *target_, v_.width, v_.height); !r) {
        error_ = r.error().message;
        target_.reset();
        return -1.0;
    }
    error_.reset();
    ValidateRect(canvas_, nullptr);
    place_bar();
    return 1e3 * std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// --- live validation and the mesh preview (slice 3, SD5) -----------------------------------------------------------

void SolidPage::checked(const std::string& text, std::span<const Diagnostic> diagnostics) {
    checked_text_ = text;
    diagnostics_.assign(diagnostics.begin(), diagnostics.end());
    locate_problems();
    update_inspector();
    invalidate();
}

void SolidPage::locate_problems() {
    problems_.clear();
    if (checked_text_.empty() || checked_text_ != text_seen_ || m_.mesh == MeshKind::unsupported) return;
    for (const Diagnostic& d : diagnostics_) problems_.push_back(locate_problem(m_, text_seen_, d.line, d.message));
}

SolidPage::Check SolidPage::check() const noexcept {
    if (checked_text_.empty() || checked_text_ != text_seen_) return Check::pending;
    return diagnostics_.empty() ? Check::valid : Check::invalid;
}

std::string SolidPage::run_blocked() const {
    if (check() != Check::invalid || m_.mesh == MeshKind::unsupported) return {};
    return "Run is off until the layout builds: " + diagnostics_.front().message;
}

void SolidPage::set_mesh(const std::string& text, std::vector<V3> edges) {
    mesh_text_ = text;
    mesh_ = std::move(edges);
    if (mesh_preview_) invalidate();
}

void SolidPage::set_mesh_preview(bool on) {
    mesh_preview_ = on;
    bar_->set_mesh(on);
    update_inspector();
    invalidate();
}

std::string SolidPage::banner() const {
    if (m_.mesh == MeshKind::unsupported) return {};
    if (check() == Check::invalid) return diagnostics_.front().message + "  (Run is off until the layout builds.)";
    if (mesh_preview_ && !mesh_shown()) return check() == Check::pending ? "Building the mesh..." : "The mesh shows when the layout builds.";
    return {};
}

// --- the inspector -------------------------------------------------------------------------------------------------

void SolidPage::update_inspector() {
    const auto show = [&](std::initializer_list<int> cs, bool on) {
        for (const int c : cs) ShowWindow(control(c), on ? SW_SHOW : SW_HIDE);
    };
    const auto set = [&](int c, const std::string& t) { win32::set_text(control(c), t); };
    const auto choose = [&](int c, const std::string& t) {
        SendMessageW(control(c), CB_SELECTSTRING, static_cast<WPARAM>(-1), reinterpret_cast<LPARAM>(win32::widen(t).c_str()));
    };
    const auto ranges = [&](const Box3& b) {
        for (int a = 0; a < 3; ++a) {
            show({c_x0_label + 4 * a, c_x0 + 4 * a, c_x1_label + 4 * a, c_x1 + 4 * a}, true);
            set(c_x0 + 4 * a, ui_length(b.lo[a]));
            set(c_x1 + 4 * a, ui_length(b.hi[a]));
        }
    };
    show({c_name_label, c_name, c_material_label, c_material, c_raise, c_lower, c_axis_label, c_axis, c_from_label, c_from, c_to_label,
          c_to, c_vertices_label, c_vertices, c_apply, c_kind_label, c_kind, c_convert, c_overlap},
         false);
    for (int k = 0; k < 12; ++k) show({c_x0_label + k}, false);
    const bool usable = m_.mesh != MeshKind::unsupported;
    show({c_resolution, c_params_label, c_params, c_params_apply}, m_.editable());
    if (GetFocus() != control(c_params)) set(c_params, parameters_text());
    show({c_draw_label, c_draw, c_plane_label, c_plane, c_offset_label, c_offset, c_section_label, c_section, c_cut_label, c_cut, c_translucent},
         usable);
    SendMessageW(control(c_plane), CB_SETCURSEL, static_cast<WPARAM>(v_.plane.normal), 0);
    set(c_offset, ui_length(v_.plane.offset));
    SendMessageW(control(c_section), CB_SETCURSEL, static_cast<WPARAM>(v_.section.axis + 1), 0);
    set(c_cut, ui_length(v_.section.at));
    std::string hint;
    if (!usable) {
        hint = message_.empty() ? std::string("Open a 3D deck to model its device.") : message_;
    } else if (m_.mesh == MeshKind::extruded) {
        hint = m_.convert_why.empty() ? "An extrude3d deck: Convert makes it layered generate3d to edit it."
                                      : "An extrude3d deck, shown only: " + m_.convert_why;
        show({c_convert}, m_.convert_why.empty());
    } else if (sel_.contact >= 0) {
        const Contact& c = m_.contacts[static_cast<std::size_t>(sel_.contact)];
        hint = "Contact " + c.name + (c.editable ? ": its patch box" : " (not editable here)");
        show({c_name_label, c_name, c_kind_label, c_kind}, true);
        set(c_name, c.name);
        choose(c_kind, c.kind);
        if (c.editable) ranges(c.box);
    } else if (sel_.solids.size() == 1) {
        const SolidRef s = sel_.solids[0];
        const Region& r = m_.regions[static_cast<std::size_t>(s.region)];
        show({c_name_label, c_name, c_material_label, c_material}, true);
        set(c_name, r.name);
        choose(c_material, r.material);
        if (s.solid < 0) {
            hint = "Background " + r.name + ": the domain";
            ranges(m_.domain);
        } else {
            const Solid& sol = r.solids[static_cast<std::size_t>(s.solid)];
            show({c_raise, c_lower, c_overlap}, true);
            SendMessageW(control(c_overlap), BM_SETCHECK, r.fill ? BST_CHECKED : BST_UNCHECKED, 0);
            if (sol.kind == SolidKind::box) {
                hint = std::format("Region {}, layer {} of {}: a box", r.name, s.region + 1, m_.regions.size());
                ranges(sol.box);
            } else {
                hint = std::format("Region {}, layer {} of {}: a prism of {} vertices", r.name, s.region + 1, m_.regions.size(), sol.polygon.size());
                show({c_axis_label, c_axis, c_from_label, c_from, c_to_label, c_to, c_vertices_label, c_vertices, c_apply}, true);
                SendMessageW(control(c_axis), CB_SETCURSEL, static_cast<WPARAM>(sol.axis), 0);
                set(c_from, ui_length(sol.from));
                set(c_to, ui_length(sol.to));
                std::string v;
                for (const Pt& p : sol.polygon) v += ui_length(p.x) + " " + ui_length(p.y) + "\r\n";
                set(c_vertices, v);
            }
        }
    } else if (sel_.solids.size() > 1) {
        hint = std::format("{} solids selected", sel_.solids.size());
    } else {
        hint = "Click a solid to select it; draw on the work plane with the tools above the view.";
    }
    if (const std::string b = banner(); !b.empty()) hint = b + "\r\n" + hint;
    set(c_hint, hint);
    layout_inspector();
}

bool SolidPage::commit(int c) {
    if (c < 0 || c >= c_count) return false;
    const auto text = [&](int k) { return win32::text_of(control(k)); };
    const auto choice = [&](int k) {
        const auto i = SendMessageW(control(k), CB_GETCURSEL, 0, 0);
        if (i < 0) return std::string();
        wchar_t b[128]{};
        SendMessageW(control(k), CB_GETLBTEXT, static_cast<WPARAM>(i), reinterpret_cast<LPARAM>(b));
        return win32::narrow(b);
    };
    const auto index = [&](int k) { return static_cast<int>(SendMessageW(control(k), CB_GETCURSEL, 0, 0)); };
    if (c == c_draw) {
        material_ = choice(c_draw);
        return true;
    }
    if (c == c_plane || c == c_offset) {
        const auto at = ui_parse(text(c_offset));
        if (index(c_plane) >= 0 && at) set_work_plane({index(c_plane), *at});
        return true;
    }
    if (c == c_section || c == c_cut) {
        const auto at = ui_parse(text(c_cut));
        set_section({index(c_section) - 1, at.value_or(0)});
        return true;
    }
    if (c == c_translucent) {
        set_translucent(SendMessageW(control(c_translucent), BM_GETCHECK, 0, 0) == BST_CHECKED);
        return true;
    }
    if (c == c_convert) {
        convert();
        return true;
    }
    if (c == c_params_apply) {
        apply_parameters(text(c_params));
        return true;
    }
    if (c == c_resolution) {
        use_automatic_resolution();
        return true;
    }
    Model next = m_;
    bool ok = false;
    const auto box_of_fields = [&]() -> std::optional<Box3> {
        std::array<std::optional<Coord>, 6> v;
        for (int k = 0; k < 6; ++k) v[static_cast<std::size_t>(k)] = ui_parse(text(c_x0 + 2 * k));
        if (std::ranges::any_of(v, [](const auto& o) { return !o; })) return std::nullopt;
        return Box3{{*v[0], *v[2], *v[4]}, {*v[1], *v[3], *v[5]}};
    };
    const bool range_field = c >= c_x0 && c <= c_z1 && (c - c_x0) % 2 == 0;
    if (sel_.contact >= 0) {
        if (c == c_name) ok = rename_contact(next, sel_.contact, text(c_name));
        else if (c == c_kind) ok = set_contact_kind(next, sel_.contact, choice(c_kind));
        else if (range_field) {
            const auto b = box_of_fields();
            ok = b && set_contact(next, sel_.contact, *b);
        }
    } else if (sel_.solids.size() == 1) {
        const SolidRef s = sel_.solids[0];
        if (c == c_name) ok = rename_region(next, s.region, text(c_name));
        else if (c == c_material) ok = solid::set_material(next, s.region, choice(c_material));
        else if (c == c_overlap) ok = set_overlap(next, s.region, SendMessageW(control(c_overlap), BM_GETCHECK, 0, 0) == BST_CHECKED);
        else if (c == c_raise && raise_region(next, s.region)) {
            ok = true;
            sel_ = {{{s.region + 1, s.solid}}, -1};
        } else if (c == c_lower && lower_region(next, s.region)) {
            ok = true;
            sel_ = {{{s.region - 1, s.solid}}, -1};
        } else if (range_field) {
            const auto b = box_of_fields();
            ok = b && (s.solid < 0 ? set_domain(next, *b) : set_box(next, s, *b));
        } else if ((c == c_apply || c == c_from || c == c_to || c == c_vertices || c == c_axis) && s.solid >= 0) {
            const Solid* old = solid_at(m_, s);
            std::vector<Pt> pts;
            std::istringstream in(text(c_vertices));
            bool good = true;
            for (std::string a, b; in >> a >> b;) {
                const auto u = ui_parse(a), w = ui_parse(b);
                if (!u || !w) good = false;
                else pts.push_back({*u, *w});
            }
            const auto from = ui_parse(text(c_from)), to = ui_parse(text(c_to));
            if (good && from && to && old && old->kind == SolidKind::prism && index(c_axis) >= 0) {
                // Replaced through the editing rules: a new prism in its place.
                Model probe = next;
                const int r = add_prism(probe, index(c_axis), *from, *to, pts, "silicon");
                if (r >= 0) {
                    next.regions[static_cast<std::size_t>(s.region)].solids[static_cast<std::size_t>(s.solid)] = probe.regions.back().solids[0];
                    ok = true;
                }
            }
        }
    }
    if (!ok) {
        host_.solid_status("That value is not taken (a name must be new and start with a letter; a solid must keep its size, a prism its "
                           "u- and v-parallel edges, a contact a face of the domain).");
        update_inspector();
        return true;
    }
    commit_model(std::move(next));
    update_inspector();
    invalidate();
    return true;
}

// --- the mouse and keys --------------------------------------------------------------------------------------------

void SolidPage::on_down(double x, double y, bool shift, bool ctrl, bool double_click) {
    SetFocus(canvas_);
    last_x_ = x;
    last_y_ = y;
    if (m_.mesh == MeshKind::unsupported) return;
    const bool editable = m_.editable();
    const auto refuse_edit = [&] {
        host_.solid_status("An extrude3d deck is edited once converted: Convert, in the Geometry section.");
    };
    const auto begin = [&](Drag d) {
        drag_ = d;
        button_ = true;
        m_start_ = m_;
        naming_ = -1;
        SetCapture(canvas_);
    };
    if (space_ || tool_ == Tool::pan) {
        // MT1: the Pan tool and Space+drag move the view only while the camera is free.
        if (!camera_fixed_) begin(Drag::pan);
        return;
    }
    // A height being given: the click ends the box or prism.
    if (drag_ == Drag::height) {
        finish_height();
        return;
    }
    if (drag_ == Drag::inset_depth) {  // MT6: the click ends the inset
        finish_inset();
        return;
    }
    entry_.clear();
    typed_axis_ = -1;
    drag_d_ = {};
    drag_s_ = 0;
    switch (tool_) {
    case Tool::select: {
        const Pick p = pick(m_, v_, x, y, sel_);
        // MT5: Mirror and Align act on a click (a face of a solid; Align also a face of the domain).
        const bool solid_face = (p.kind == PickKind::face || (p.kind == PickKind::handle && p.handle == HandleKind::face)) && p.solid.solid >= 0;
        if ((mode_ == modeller::Mode::mirror || mode_ == modeller::Mode::align) && editable &&
            (solid_face || (mode_ == modeller::Mode::align && p.kind == PickKind::face && p.domain_face >= 0))) {
            transform_click(p, shift, ctrl);
            invalidate();
            return;
        }
        // Unit 26m slice 2 (MT6, MT7): the Faces and Edges modes.
        if (editable && tool_press(p, x, y)) {
            invalidate();
            return;
        }
        // MT2: what is picked opens its group.
        const auto open = [&](modeller::Picked k) {
            const modeller::Group g = modeller::group_for(k);
            if (g != bar_->group()) {
                bar_->set_group(g);
                bar_->set_tool(modeller::armed_key(static_cast<int>(tool_), mode_));
            }
        };
        if (p.kind == PickKind::handle) open(p.handle == HandleKind::face ? modeller::Picked::face : modeller::Picked::solid);
        else if (p.kind == PickKind::edge) open(modeller::Picked::edge);
        else if (p.kind == PickKind::contact) open(modeller::Picked::contact);
        else if (p.kind == PickKind::face && p.solid.solid >= 0) open(modeller::Picked::solid);
        else open(modeller::Picked::nothing);
        // MT1, MQ3: empty background (or a region's background, which no drag moves) orbits or pans a free camera.
        const bool background = p.kind == PickKind::none || (p.kind == PickKind::face && p.solid.solid < 0);
        const modeller::Nav nav = modeller::navigation(camera_fixed_, true, {modeller::Button::left, background, shift, true, false});
        if (nav == modeller::Nav::orbit || nav == modeller::Nav::pan) {
            if (!shift) {
                sel_ = p.kind == PickKind::face && p.solid.region >= 0 ? Selection{{p.solid}, -1} : Selection{};
                update_inspector();
                invalidate();
            }
            begin(nav == modeller::Nav::orbit ? Drag::orbit : Drag::pan);
            return;
        }
        if (p.kind == PickKind::handle) {
            if (!editable) return refuse_edit();
            ref_ = p.solid;
            index_ = p.index;
            const Solid* s = solid_at(m_, ref_);
            const auto hs = handles(m_, v_, ref_);
            const auto h = std::ranges::find_if(hs, [&](const Handle& k) { return k.kind == p.handle && k.index == p.index; });
            if (!s || h == hs.end()) return;
            origin_ = h->pts[0];
            alone_ = ctrl;
            if (p.handle == HandleKind::axis) {
                axis_ = p.index;
                start_s_ = along_axis(m_, v_, x, y, origin_, axis_);
                begin(Drag::axis);
            } else if (p.handle == HandleKind::plane) {
                axis_ = p.index;
                const auto q = on_plane(m_, v_, x, y, axis_, origin_[static_cast<std::size_t>(axis_)]);
                if (!q) return;
                start_p_ = *q;
                begin(Drag::plane);
            } else if (p.handle == HandleKind::face) {
                axis_ = face_axis(*s, p.index);
                if (axis_ < 0) return;
                start_s_ = along_axis(m_, v_, x, y, origin_, axis_);
                begin(Drag::face);
            } else {
                axis_ = s->axis;
                const auto q = on_plane(m_, v_, x, y, axis_, origin_[static_cast<std::size_t>(axis_)]);
                if (!q) return;
                start_p_ = *q;
                begin(Drag::vertex);
            }
            return;
        }
        if (p.kind == PickKind::contact) {
            sel_ = {{}, p.contact};
            update_inspector();
            invalidate();
            if (!editable) return;
            const Box3 b = m_.contacts[static_cast<std::size_t>(p.contact)].box;
            axis_ = p.domain_face >= 0 ? p.domain_face / 2 : 2;
            for (int a = 0; a < 3; ++a) {
                if (b.lo[a] == b.hi[a]) axis_ = a;
            }
            start_p_ = p.at;
            begin(Drag::contact_move);
            return;
        }
        if (p.kind == PickKind::face || p.kind == PickKind::edge) {
            const SolidRef s = p.solid;
            if (s.region < 0) {
                sel_ = {};
            } else if (shift && s.solid >= 0) {
                const auto at = std::ranges::find(sel_.solids, s);
                if (at == sel_.solids.end()) sel_.solids.push_back(s);
                else sel_.solids.erase(at);
                sel_.contact = -1;
            } else if (std::ranges::find(sel_.solids, s) == sel_.solids.end()) {
                sel_ = {{s}, -1};
            }
            update_inspector();
            invalidate();
            if (s.solid < 0 || !editable || shift) return;
            const Solid* sol = solid_at(m_, s);
            axis_ = p.kind == PickKind::face && sol ? face_axis(*sol, p.index) : v_.plane.normal;
            if (axis_ < 0) axis_ = v_.plane.normal;
            start_p_ = p.kind == PickKind::face ? p.at : um_of(sol->bounds().lo);
            if (p.kind == PickKind::edge) {
                const auto q = on_plane(m_, v_, x, y, axis_, start_p_[static_cast<std::size_t>(axis_)]);
                if (!q) return;
                start_p_ = *q;
            }
            if (mode_ == modeller::Mode::array) {
                // MT5: the drag is the pitch, the wheel the count; one region per copy on release.
                sel_ = {{s}, -1};
                ref_ = s;
                update_inspector();
                begin(Drag::array);
                return;
            }
            begin(Drag::move);
            if (mode_ == modeller::Mode::copy || GetKeyState(VK_MENU) < 0) {
                // MT5: a copy of each picked solid is dragged (its region's rule kept); the originals stay.
                copy_from_ = m_;
                copy_sel_ = sel_;
                Model with = m_;
                std::vector<SolidRef> copies;
                for (const SolidRef c : sel_.solids) {
                    if (c.solid < 0) continue;
                    const int r = duplicate_solid(with, c, {});
                    if (r < 0) continue;
                    with.regions[static_cast<std::size_t>(r)].fill = with.regions[static_cast<std::size_t>(c.region)].fill;
                    copies.push_back({r, 0});
                }
                m_ = with;
                m_start_ = with;
                sel_ = {copies, -1};
            }
            return;
        }
        sel_ = {};
        update_inspector();
        invalidate();
        return;
    }
    case Tool::box: {
        if (!editable) return refuse_edit();
        const auto p = plane_point(x, y);
        if (!p) return;
        a_ = b_ = *p;
        polygon_.clear();
        begin(Drag::base);
        return;
    }
    case Tool::prism: {
        if (!editable) return refuse_edit();
        if (double_click) {
            close_polygon();
            return;
        }
        auto p = plane_point(x, y);
        if (!p) return;
        const int n = v_.plane.normal;
        if (!polygon_.empty()) {
            // Clicks as they are (slanted edges since slice 4); Shift keeps the edge along u or v from the last point.
            if (shift) *p = from_uv(n, geo::constrain(uv_of(n, polygon_.back()), uv_of(n, *p)), v_.plane.offset);
            if (polygon_.size() >= 3 && *p == polygon_.front()) {
                close_polygon();
                return;
            }
            if (*p == polygon_.back()) return;
        }
        polygon_.push_back(*p);
        drag_ = Drag::polygon;
        preview_.polyline = polygon_;
        invalidate();
        return;
    }
    case Tool::contact: {
        if (!editable) return refuse_edit();
        const Pick p = pick(m_, v_, x, y, {});
        if (p.domain_face < 0) {
            host_.solid_status("A contact is drawn on a face of the domain.");
            return;
        }
        axis_ = p.domain_face / 2;
        const Coord c = p.domain_face % 2 ? m_.domain.hi[axis_] : m_.domain.lo[axis_];
        P3 q = quanta_of(p.at);
        q[axis_] = c;
        a_ = b_ = snap(m_, {axis_, c}, q, snap_tolerance(m_, v_, p.at), std::nullopt, grid()).p;
        a_[axis_] = b_[axis_] = c;
        begin(Drag::contact);
        return;
    }
    case Tool::cylinder: {  // MT7: the centre on the work plane, then the radius, then the height
        if (!editable) return refuse_edit();
        const auto p = plane_point(x, y);
        if (!p) return;
        a_ = b_ = *p;
        radius_ = 0;
        polygon_.clear();
        preview_ = {};
        begin(Drag::radius);
        return;
    }
    case Tool::measure: {
        const Pick p = pick(m_, v_, x, y, {});
        std::optional<P3> q;
        if (p.kind == PickKind::face || p.kind == PickKind::contact) {
            q = snap(m_, v_.plane, quanta_of(p.at), snap_tolerance(m_, v_, p.at), std::nullopt, 0).p;
        } else {
            q = plane_point(x, y);
        }
        if (!q) return;
        a_ = b_ = *q;
        preview_.measure = std::pair{a_, b_};
        begin(Drag::measure);
        return;
    }
    case Tool::pan: return;
    }
}

void SolidPage::close_polygon() {
    if (polygon_.size() < 3) return;
    std::vector<Pt> uv;
    for (const P3& p : polygon_) uv.push_back(uv_of(v_.plane.normal, p));
    if (!prism_ok(uv)) {
        host_.solid_status("A prism's polygon is simple, of non-zero area, with three vertices or more.");
        polygon_.clear();
        preview_ = {};
        drag_ = Drag::none;
        invalidate();
        return;
    }
    a_ = b_ = polygon_.back();
    height_ = grid();
    drag_ = Drag::height;
    preview_.polyline.clear();
    preview_.solids = {Solid{SolidKind::prism, {}, v_.plane.normal, v_.plane.offset, v_.plane.offset + height_, uv}};
    invalidate();
}

void SolidPage::on_move(double x, double y, bool shift, bool ctrl) {
    const double dx = x - last_x_, dy = y - last_y_;
    last_x_ = x;
    last_y_ = y;
    const double um = static_cast<double>(per_um);
    const auto with_drag = [&](Model next) {
        m_ = std::move(next);
        update_inspector();
        invalidate();
    };
    // MT3: a typed value replaces the dragged one, in the drag's direction (a typed minus sign: as typed).
    const auto typed_or = [&](Coord dragged) {
        const auto v = entry_.value();
        if (!v) return dragged;
        return *v >= 0 && dragged < 0 ? -*v : *v;
    };
    switch (drag_) {
    case Drag::implant:
    case Drag::junction:
    case Drag::max_size:
    case Drag::refine: doping_step(x, y); return;
    case Drag::orbit:
        v_.camera.azimuth_deg -= static_cast<float>(dx * 0.4);
        v_.camera.elevation_deg = std::clamp(v_.camera.elevation_deg + static_cast<float>(dy * 0.4), -89.0f, 89.0f);
        invalidate();
        return;
    case Drag::pan: {
        const render::ViewFrame f = render::view_frame(frame_scene(m_, v_), v_.width, v_.height);
        v_.camera.pan[0] += static_cast<float>(2.0 * dx / (f.right_px - f.left));
        v_.camera.pan[1] -= static_cast<float>(2.0 * dy / (f.bottom_px - f.top));
        invalidate();
        return;
    }
    case Drag::axis: {
        Model next = m_start_;
        P3 d{};
        d[axis_] = typed_or(quantize(along_axis(m_start_, v_, x, y, origin_, axis_) - start_s_));
        drag_d_ = d;
        move_solid(next, ref_, d);
        with_drag(std::move(next));
        return;
    }
    case Drag::move:
    case Drag::array:
    case Drag::plane:
    case Drag::contact_move: {
        const auto q = on_plane(m_start_, v_, x, y, axis_, start_p_[static_cast<std::size_t>(axis_)]);
        if (!q) return;
        P3 d{};
        for (int a = 0; a < 3; ++a) {
            if (a != axis_) d[a] = quantize((*q)[static_cast<std::size_t>(a)] - start_p_[static_cast<std::size_t>(a)]);
        }
        if (typed_axis_ >= 0 && drag_ != Drag::contact_move) d[typed_axis_] = typed_or(d[typed_axis_]);
        // MT12 (2): copies keep clear of other solids' faces and the domain's: onto them or 4 s off.
        if (res_q_ > 0 && (drag_ == Drag::array || (drag_ == Drag::move && copy_from_))) {
            std::vector<Solid> moving;
            std::vector<SolidRef> skip;
            if (drag_ == Drag::array) {
                if (const Solid* s = solid_at(m_start_, ref_)) moving.push_back(*s);
            } else {
                for (const SolidRef c : sel_.solids) {
                    if (const Solid* s = solid_at(m_start_, c)) {
                        moving.push_back(*s);
                        skip.push_back(c);
                    }
                }
            }
            for (int a = 0; a < 3; ++a) {
                if (a != axis_ && d[a] != 0) d[a] = keep_clear(m_start_, moving, skip, a, d[a], 4 * res_q_);
            }
        }
        drag_d_ = d;
        Model next = m_start_;
        if (drag_ == Drag::contact_move) {
            move_contact(next, sel_.contact, d);
        } else if (drag_ == Drag::array) {
            if (!(d == P3{})) array_solid(next, ref_, d, array_count_);
        } else if (drag_ == Drag::plane) {
            move_solid(next, ref_, d);
        } else {
            for (const SolidRef s : sel_.solids) move_solid(next, s, d);
        }
        with_drag(std::move(next));
        return;
    }
    case Drag::face: {
        Model next = m_start_;
        drag_s_ = typed_or(quantize(along_axis(m_start_, v_, x, y, origin_, axis_) - start_s_));
        move_face(next, ref_, index_, drag_s_);
        if (naming_ >= 0) static_cast<void>(name_extent(next, ref_, index_ ^ 1, index_, naming_name_));  // MT8: N
        with_drag(std::move(next));
        return;
    }
    case Drag::vertex: {
        const auto q = on_plane(m_start_, v_, x, y, axis_, start_p_[static_cast<std::size_t>(axis_)]);
        const Solid* s = solid_at(m_start_, ref_);
        if (!q || !s) return;
        const Pt d = uv_of(axis_, {quantize((*q)[0] - start_p_[0]), quantize((*q)[1] - start_p_[1]), quantize((*q)[2] - start_p_[2])});
        const Pt p0 = s->polygon[static_cast<std::size_t>(index_)];
        Model next = m_start_;
        drag_d_ = from_uv(axis_, d, 0);
        set_vertex(next, ref_, index_, {p0.x + d.x, p0.y + d.y}, !(alone_ || ctrl));
        with_drag(std::move(next));
        return;
    }
    case Drag::base:
    case Drag::contact: {
        std::optional<P3> p;
        if (drag_ == Drag::base) {
            p = plane_point(x, y);
        } else {
            const Coord c = a_[axis_];
            if (const auto q = on_plane(m_, v_, x, y, axis_, static_cast<double>(c) / um)) {
                P3 raw = quanta_of(*q);
                raw[axis_] = c;
                p = snap(m_, {axis_, c}, raw, snap_tolerance(m_, v_, *q), std::nullopt, grid()).p;
                (*p)[axis_] = c;
            }
        }
        if (!p) return;
        b_ = *p;
        if (drag_ == Drag::base && typed_axis_ >= 0) b_[typed_axis_] = a_[typed_axis_] + typed_or(b_[typed_axis_] - a_[typed_axis_]);
        preview_.solids = {Solid{SolidKind::box, normalized(a_, b_), 2, 0, 0, {}}};
        invalidate();
        return;
    }
    case Drag::height: {
        const int n = v_.plane.normal;
        const V3 corner = um_of(b_);
        height_ = typed_or(quantize(along_axis(m_, v_, x, y, corner, n) - corner[static_cast<std::size_t>(n)]));
        if (height_ == 0) height_ = grid();
        const Coord lo = std::min(v_.plane.offset, v_.plane.offset + height_), hi = std::max(v_.plane.offset, v_.plane.offset + height_);
        if (polygon_.empty()) {
            Box3 b = normalized(a_, b_);
            b.lo[n] = lo;
            b.hi[n] = hi;
            preview_.solids = {Solid{SolidKind::box, b, 2, 0, 0, {}}};
        } else {
            std::vector<Pt> uv;
            for (const P3& p : polygon_) uv.push_back(uv_of(n, p));
            preview_.solids = {Solid{SolidKind::prism, {}, n, lo, hi, uv}};
        }
        invalidate();
        return;
    }
    case Drag::cut:
    case Drag::taper: {
        const auto q = on_plane(m_start_, v_, x, y, axis_, start_p_[static_cast<std::size_t>(axis_)]);
        if (!q) return;
        P3 raw{};
        for (int a = 0; a < 3; ++a) raw[a] = quantize((*q)[static_cast<std::size_t>(a)] - start_p_[static_cast<std::size_t>(a)]);
        Pt d = uv_of(axis_, raw);
        if (drag_ == Drag::taper) {
            // Sideways along one axis of the cross-section: the larger offset's.
            if (std::abs(d.x) >= std::abs(d.y)) {
                d.y = 0;
                d.x = typed_or(d.x);
            } else {
                d.x = 0;
                d.y = typed_or(d.y);
            }
            drag_d_ = from_uv(axis_, d, 0);
            Model next = m_start_;
            if (d.x != 0 || d.y != 0) taper_edge(next, ref_, edge_, d);
            with_drag(std::move(next));
            return;
        }
        const double len = std::hypot(static_cast<double>(d.x), static_cast<double>(d.y)) / um;
        drag_s_ = std::max<Coord>(0, typed_or(quantize(len)));
        with_drag(cut_model(drag_s_));
        return;
    }
    case Drag::extrude:
    case Drag::inset_depth: {
        const Coord raw = quantize(along_axis(m_start_, v_, x, y, origin_, normal_.axis) - start_s_) * normal_.sign;
        drag_s_ = gap_up(typed_or(raw));  // MT12 (2): no layer under 4 s
        Model next = m_start_;
        if (drag_ == Drag::extrude) {
            if (drag_s_ > 0) extrude_face(next, ref_, index_, drag_s_, material_);
        } else if (drag_s_ != 0) {
            inset_face(next, ref_, index_, margin_, drag_s_, material_);
        }
        with_drag(std::move(next));
        return;
    }
    case Drag::inset: {
        const auto q = on_plane(m_start_, v_, x, y, normal_.axis, start_p_[static_cast<std::size_t>(normal_.axis)]);
        if (!q) return;
        double reach = 0;
        for (int a = 0; a < 3; ++a) {
            if (a != normal_.axis) reach = std::max(reach, std::abs((*q)[static_cast<std::size_t>(a)] - start_p_[static_cast<std::size_t>(a)]));
        }
        margin_ = gap_up(std::max<Coord>(0, typed_or(quantize(reach))));  // MT12 (2): a margin of 0 or at least 4 s
        // Shown with a depth of one grid step until the depth is given.
        Model next = m_start_;
        inset_face(next, ref_, index_, margin_, gap_up(grid()), material_);
        with_drag(std::move(next));
        return;
    }
    case Drag::radius: {
        const auto p = plane_point(x, y);
        if (!p) return;
        b_ = *p;
        const int n = v_.plane.normal;
        const Pt c = uv_of(n, a_), e = uv_of(n, b_);
        radius_ = std::max<Coord>(0, typed_or(quantize(std::hypot(static_cast<double>(e.x - c.x), static_cast<double>(e.y - c.y)) / um)));
        // MQ4, MT12 (2): 16 sides, fewer while a side is under 6 s.
        const int sides = res_q_ > 0 ? geo::polygon_sides(static_cast<double>(radius_), 6.0 * static_cast<double>(res_q_)) : 16;
        preview_.solids.clear();
        if (radius_ > 0 && sides >= 3) {
            preview_.solids = {Solid{SolidKind::prism, {}, n, v_.plane.offset, v_.plane.offset + grid(), geo::regular_polygon(c, radius_, sides)}};
        }
        invalidate();
        return;
    }

    case Drag::polygon: {
        if (const auto p = plane_point(x, y)) {
            const int n = v_.plane.normal;
            preview_.polyline = polygon_;
            preview_.polyline.push_back(shift ? from_uv(n, geo::constrain(uv_of(n, polygon_.back()), uv_of(n, *p)), v_.plane.offset) : *p);
            invalidate();
        }
        return;
    }
    case Drag::measure: {
        const Pick p = pick(m_, v_, x, y, {});
        std::optional<P3> q = p.kind == PickKind::face || p.kind == PickKind::contact
                                  ? std::optional(snap(m_, v_.plane, quanta_of(p.at), snap_tolerance(m_, v_, p.at), std::nullopt, 0).p)
                                  : plane_point(x, y);
        if (!q) return;
        b_ = *q;
        preview_.measure = std::pair{a_, b_};
        const double d = std::sqrt(std::pow(static_cast<double>(b_.x - a_.x), 2) + std::pow(static_cast<double>(b_.y - a_.y), 2) +
                                   std::pow(static_cast<double>(b_.z - a_.z), 2)) / um;
        host_.solid_status(std::format("distance {:.6g} µm", d));
        invalidate();
        return;
    }
    case Drag::none: {
        if (m_.mesh == MeshKind::unsupported) return;
        if (const auto p = plane_point(x, y)) {
            host_.solid_status(std::format("x = {} µm, y = {} µm, z = {} µm", ui_length(p->x), ui_length(p->y), ui_length(p->z)));
        }
        return;
    }
    }
    (void)shift;
}

// --- Slice 4 (MT9, MT10): the Doping and Mesh tools, each writing one deck line ---------------------------------------

bool SolidPage::commit_text(const std::string& out) {
    const std::string text = host_.solid_deck();
    if (out.empty() || out == text) return false;
    history_.record(text);
    apply_text(out);
    return true;
}

bool SolidPage::doping_press(const Pick& p, double x, double y) {
    using modeller::Mode;
    const auto dopant = [&] { return donor_ ? device::Dopant::donor : device::Dopant::acceptor; };
    pending_text_.clear();
    pending_volume_.reset();
    drag_s_ = 0;
    entry_.clear();
    if (mode_ == Mode::max_size) {
        if (!(max_size_ > 0.0)) {
            host_.solid_status("The deck has no max_size to change.");
            return true;
        }
        size_from_ = max_size_;
        press_y_ = y;
        begin_drag(Drag::max_size);
        doping_step(x, y);
        return true;
    }
    if (mode_ == Mode::refine) {
        const auto q = plane_point(x, y);
        if (!q) return true;
        a_ = b_ = *q;
        if (refine_size_ <= 0) refine_size_ = std::max<Coord>(1, static_cast<Coord>(std::llround(max_size_ * 1e8 / 4.0)));
        begin_drag(Drag::refine);
        return true;
    }
    const bool on_face = (p.kind == PickKind::face || (p.kind == PickKind::handle && p.handle == HandleKind::face)) && p.solid.solid >= 0;
    if (!on_face) return false;
    if (mode_ == Mode::profile) {
        sel_ = {{p.solid}, -1};
        const auto out = solid::add_profile(host_.solid_deck(), solid::region_profile(m_, p.solid, dopant(), level_));
        if (!out) {
            host_.solid_status(out.error());
            return true;
        }
        host_.solid_status(modeller::deck_line(host_.solid_deck(), *out));
        commit_text(*out);
        return true;
    }
    const Solid* s = solid_at(m_, p.solid);
    const auto n = s ? face_normal(*s, p.index) : std::nullopt;
    if (!n) {
        host_.solid_status("A slanted face takes no implant or junction; press another face.");
        return true;
    }
    ref_ = p.solid;
    index_ = p.index;
    normal_ = *n;
    origin_ = p.at;
    sel_ = {{p.solid}, -1};
    start_s_ = along_axis(m_, v_, x, y, origin_, n->axis);
    begin_drag(mode_ == Mode::implant ? Drag::implant : Drag::junction);
    return true;
}

void SolidPage::doping_step(double x, double y) {
    const auto dopant = donor_ ? device::Dopant::donor : device::Dopant::acceptor;
    std::expected<std::string, std::string> out = std::unexpected(std::string());
    pending_volume_.reset();
    switch (drag_) {
    case Drag::implant:
    case Drag::junction: {
        const Coord raw = quantize(along_axis(m_start_, v_, x, y, origin_, normal_.axis) - start_s_) * -normal_.sign;  // into the solid
        const auto typed = entry_.value();
        drag_s_ = typed ? std::abs(*typed) : std::max<Coord>(0, raw);
        const Solid* s = solid_at(m_start_, ref_);
        if (!s || drag_s_ <= 0) break;
        const Coord at = quanta_of(origin_)[normal_.axis];
        Box3 window = s->bounds();
        window.lo[normal_.axis] = window.hi[normal_.axis] = at;
        if (drag_ == Drag::implant) {
            const auto p = solid::implant_profile(normal_, at, window, drag_s_, 0, dopant, level_);
            out = solid::add_profile(text_seen_, p);
            pending_volume_ = solid::Volume{solid::profile_extent(m_.domain, p), solid::doping_color(donor_, level_)};
        } else {
            const auto p = solid::junction_profile(profiles_, normal_, at, window, drag_s_, dopant, level_);
            if (!p) {
                out = std::unexpected(p.error());
                break;
            }
            out = solid::add_profile(text_seen_, *p);
            pending_volume_ = solid::Volume{solid::profile_extent(m_.domain, *p), solid::doping_color(donor_, level_)};
        }
        break;
    }
    case Drag::max_size: {
        // Up doubles it every 100 px (at 96 dpi), down halves it; a typed value is the size.
        const double px = 100.0 * dpi_ / 96.0;
        const auto typed = entry_.value();
        drag_s_ = typed ? std::abs(*typed) : std::max<Coord>(1, static_cast<Coord>(std::llround(size_from_ * std::pow(2.0, (press_y_ - y) / px) * 1e8)));
        if (drag_s_ > 0) out = solid::set_mesh_length(text_seen_, "max_size", drag_s_);
        break;
    }
    case Drag::refine: {
        if (const auto q = plane_point(x, y)) b_ = *q;
        const int n = v_.plane.normal;
        Box3 b = normalized(a_, b_);
        b.lo[n] = m_.domain.lo[n];
        b.hi[n] = m_.domain.hi[n];
        const solid::SizeBox3 sb{b, refine_size_};
        out = solid::add_size_box(text_seen_, sb);
        pending_volume_ = solid::Volume{b, {0.55f, 0.55f, 0.60f, 0.22f}};
        break;
    }
    default: return;
    }
    pending_text_ = out ? *out : std::string();
    live_line_ = out ? modeller::deck_line(text_seen_, *out) : out.error();
    if (!live_line_.empty()) host_.solid_status(live_line_);
    const int off = MulDiv(16, dpi_, 96);
    bar_->chip(static_cast<int>(x) + off, static_cast<int>(y) + canvas_top_ + off, chip_value());
    invalidate();
}

void SolidPage::on_up(double /*x*/, double /*y*/) {  // the release's place: no step uses it
    button_ = false;
    if (drag_ != Drag::base && drag_ != Drag::inset && drag_ != Drag::radius) end_drag_ui();
    if (GetCapture() == canvas_) ReleaseCapture();
    switch (drag_) {
    case Drag::inset: {  // MT6: the margin given; now the depth, a click ending it
        drag_ = Drag::inset_depth;
        entry_.clear();
        // The depth from the face itself (its pressed point along the normal), not from the release pixel's ray: the
        // margin drag ends off the normal line, and its nearest point there is not the face's.
        start_s_ = origin_[static_cast<std::size_t>(normal_.axis)];
        drag_s_ = 0;
        m_ = m_start_;
        invalidate();
        return;
    }
    case Drag::implant:
    case Drag::junction:
    case Drag::max_size:
    case Drag::refine: {  // slice 4: the line the drag shows, written
        const std::string out = pending_text_;
        drag_ = Drag::none;
        pending_text_.clear();
        pending_volume_.reset();
        if (!out.empty()) commit_text(out);
        update_inspector();
        invalidate();
        return;
    }
    case Drag::radius: {  // MT7: the radius given; now the height, as a prism's
        if (preview_.solids.empty()) {
            host_.solid_status(radius_ > 0 ? "A cylinder that small has sides under 6 s (lower interface_resolution)." : "A cylinder needs a radius.");
            cancel();
            return;
        }
        const int n = v_.plane.normal;
        polygon_.clear();
        for (const Pt& p : preview_.solids.front().polygon) polygon_.push_back(from_uv(n, p, v_.plane.offset));
        entry_.clear();
        typed_axis_ = -1;
        height_ = grid();
        b_ = polygon_.front();
        drag_ = Drag::height;
        invalidate();
        return;
    }
    case Drag::cut:
    case Drag::taper:
    case Drag::extrude:
    case Drag::move:
    case Drag::array:
    case Drag::axis:
    case Drag::plane:
    case Drag::face:
    case Drag::vertex:
    case Drag::contact_move: {
        Model done = m_;
        const bool moved = !(done == m_start_);
        m_ = copy_from_ ? *copy_from_ : m_start_;
        if (copy_from_ && !moved) sel_ = copy_sel_;  // a copy not moved is no copy
        copy_from_.reset();
        drag_ = Drag::none;
        if (moved) commit_model(std::move(done));
        update_inspector();
        invalidate();
        return;
    }
    case Drag::base: {
        const int n = v_.plane.normal, u = (n + 1) % 3, w = (n + 2) % 3;
        if (a_[u] == b_[u] || a_[w] == b_[w]) {
            cancel();
            return;
        }
        drag_ = Drag::height;
        entry_.clear();
        typed_axis_ = -1;
        height_ = grid();
        Box3 b = normalized(a_, b_);
        b.hi[n] = b.lo[n] + height_;
        preview_.solids = {Solid{SolidKind::box, b, 2, 0, 0, {}}};
        invalidate();
        return;
    }
    case Drag::contact: {
        Box3 b = normalized(a_, b_);
        drag_ = Drag::none;
        preview_ = {};
        Model next = m_;
        const std::string name = fresh_contact(m_);
        const int c = add_contact(next, name, "ohmic", b);
        if (c < 0) {
            host_.solid_status("A contact is a rectangle with area on a face of the domain.");
        } else if (commit_model(std::move(next))) {
            sel_ = {{}, c};
        }
        update_inspector();
        invalidate();
        return;
    }
    case Drag::measure:
    case Drag::pan:
    case Drag::orbit:
        drag_ = Drag::none;
        invalidate();
        return;
    default: return;
    }
}

bool SolidPage::on_key(UINT key, bool shift, bool ctrl) {
    if (ctrl) return false;
    // MT3: while a drag takes a typed value, letters and digits are its (WM_CHAR), not the tools' keys.
    const bool typing = drag_ == Drag::axis || drag_ == Drag::face || drag_ == Drag::height || drag_ == Drag::move || drag_ == Drag::plane ||
                        drag_ == Drag::base || drag_ == Drag::array || drag_ == Drag::cut || drag_ == Drag::taper || drag_ == Drag::extrude ||
                        drag_ == Drag::inset || drag_ == Drag::inset_depth || drag_ == Drag::radius || drag_ == Drag::implant ||
                        drag_ == Drag::junction || drag_ == Drag::max_size;
    if (typing && ((key >= 'A' && key <= 'Z') || (key >= '0' && key <= '9') || (key >= VK_NUMPAD0 && key <= VK_DIVIDE) ||
                   key == VK_OEM_PERIOD || key == VK_OEM_MINUS || key == VK_BACK)) {
        return true;
    }
    switch (key) {
    case 'L': set_camera_fixed(!camera_fixed_); return true;
    case 'V': set_tool(Tool::select); return true;
    case 'B':
    case 'R': set_tool(Tool::box); return true;
    case 'P': set_tool(Tool::prism); return true;
    case 'O': set_tool(Tool::cylinder); return true;
    case 'C': set_tool(Tool::contact); return true;
    case 'H': set_tool(Tool::pan); return true;
    case 'M': set_tool(Tool::measure); return true;
    case 'F': fit(); return true;
    case 'T': set_translucent(!v_.translucent); return true;
    case 'W': {
        const Pick p = pick(m_, v_, last_x_, last_y_, {});
        if (p.kind != PickKind::face && p.kind != PickKind::contact) return true;
        int axis = p.domain_face >= 0 ? p.domain_face / 2 : -1;
        if (const Solid* s = solid_at(m_, p.solid); s && p.index >= 0) axis = face_axis(*s, p.index);
        if (axis >= 0) set_work_plane({axis, quanta_of(p.at)[axis]});
        return true;
    }
    case VK_ESCAPE: cancel(); return true;
    case VK_RETURN:
        if (drag_ == Drag::polygon) close_polygon();
        else apply_drag();
        return true;
    case VK_DELETE: remove(); return true;
    case VK_LEFT:
    case VK_RIGHT:
    case VK_UP:
    case VK_DOWN:
    case VK_PRIOR:
    case VK_NEXT: {
        if (sel_.solids.empty() || !m_.editable()) return true;
        const Coord step = grid() * (shift ? 10 : 1);
        const int n = v_.plane.normal, u = (n + 1) % 3, w = (n + 2) % 3;
        P3 d{};
        if (key == VK_LEFT) d[u] = -step;
        if (key == VK_RIGHT) d[u] = step;
        if (key == VK_UP) d[w] = step;
        if (key == VK_DOWN) d[w] = -step;
        if (key == VK_PRIOR) d[n] = step;
        if (key == VK_NEXT) d[n] = -step;
        Model next = m_;
        bool ok = true;
        for (const SolidRef s : sel_.solids) ok = s.solid >= 0 && move_solid(next, s, d) && ok;
        if (ok) commit_model(std::move(next));
        return true;
    }
    default: return false;
    }
}

// --- window procedures ---------------------------------------------------------------------------------------------

LRESULT CALLBACK SolidPage::panel_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l) {
    auto* self = reinterpret_cast<SolidPage*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
    if (self) {
        if (message == WM_SIZE && hwnd == self->panel_) {
            self->layout();
            return 0;
        }
        if (message == WM_COMMAND && HIWORD(w) == BN_CLICKED) {
            if (const auto c = self->bar_->command(LOWORD(w))) {
                using Kind = ModellerBar::Command::Kind;
                switch (c->kind) {
                case Kind::group: self->set_group(static_cast<modeller::Group>(c->index)); break;
                case Kind::tool: self->set_tool(static_cast<Tool>(c->index)); break;
                case Kind::duplicate: self->duplicate(); break;
                case Kind::mode: self->set_mode(static_cast<modeller::Mode>(c->index)); break;
                case Kind::overlap: self->toggle_overlap(); break;
                case Kind::corner:
                    switch (c->index) {
                    case ModellerBar::camera: self->set_camera_fixed(!self->camera_fixed_); break;
                    case ModellerBar::reset: self->reset_view(); break;
                    case ModellerBar::fit: self->fit(); break;
                    case ModellerBar::measure: self->set_tool(self->tool_ == Tool::measure ? Tool::select : Tool::measure); break;
                    case ModellerBar::mesh: self->set_mesh_preview(!self->mesh_preview_); break;
                    default: break;
                    }
                    break;
                }
                SetFocus(self->canvas_);
                return 0;
            }
        }
    }
    return DefWindowProcW(hwnd, message, w, l);
}

LRESULT CALLBACK SolidPage::inspector_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l) {
    auto* self = reinterpret_cast<SolidPage*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
    if (self) {
        if (message == WM_SIZE) {
            self->layout_inspector();
            return 0;
        }
        if (message == WM_COMMAND) {
            const int c = LOWORD(w) - self->inspector_id_;
            const int code = HIWORD(w);
            if ((code == BN_CLICKED && (c == c_raise || c == c_lower || c == c_apply || c == c_translucent || c == c_convert || c == c_overlap || c == c_resolution || c == c_params_apply)) ||
                (code == CBN_SELCHANGE && (c == c_material || c == c_kind || c == c_draw || c == c_axis || c == c_plane || c == c_section))) {
                self->commit(c);
                return 0;
            }
        }
    }
    return DefWindowProcW(hwnd, message, w, l);
}

LRESULT CALLBACK SolidPage::canvas_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l) {
    if (message == WM_NCCREATE) {
        SetWindowLongPtrW(hwnd, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(reinterpret_cast<CREATESTRUCTW*>(l)->lpCreateParams));
    }
    auto* self = reinterpret_cast<SolidPage*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
    if (self && self->canvas_ == hwnd) return self->on_canvas(message, w, l);
    if (self && message == WM_CREATE) {
        self->canvas_ = hwnd;
        return 0;
    }
    return DefWindowProcW(hwnd, message, w, l);
}

LRESULT SolidPage::on_canvas(UINT message, WPARAM w, LPARAM l) {
    const double x = GET_X_LPARAM(l), y = GET_Y_LPARAM(l);
    const bool shift = (w & MK_SHIFT) != 0, ctrl = (w & MK_CONTROL) != 0;
    switch (message) {
    case WM_GETDLGCODE: return DLGC_WANTALLKEYS | DLGC_WANTARROWS | DLGC_WANTCHARS;
    case WM_ERASEBKGND: return 1;
    case WM_PAINT: {
        PAINTSTRUCT ps;
        BeginPaint(canvas_, &ps);
        EndPaint(canvas_, &ps);
        render();
        return 0;
    }
    case WM_SIZE: invalidate(); return 0;
    case WM_LBUTTONDOWN: on_down(x, y, shift, ctrl, false); return 0;
    case WM_LBUTTONDBLCLK: on_down(x, y, shift, ctrl, true); return 0;
    case WM_MOUSEMOVE:
        on_move(x, y, shift, ctrl);
        after_drag_step(x, y);
        return 0;
    case WM_LBUTTONUP:
        if (drag_ != Drag::none && button_) on_up(x, y);
        return 0;
    case WM_RBUTTONDOWN:
    case WM_MBUTTONDOWN:
        SetFocus(canvas_);
        if (drag_ != Drag::none && drag_ != Drag::polygon && drag_ != Drag::height) return 0;
        last_x_ = x;
        last_y_ = y;
        if (drag_ == Drag::polygon || drag_ == Drag::height) return 0;  // a drawing keeps going
        // MT1: a fixed camera takes no middle or right drag.
        if (modeller::navigation(camera_fixed_, true, {message == WM_RBUTTONDOWN ? modeller::Button::right : modeller::Button::middle}) == modeller::Nav::none) return 0;
        drag_ = message == WM_RBUTTONDOWN ? Drag::orbit : Drag::pan;
        button_ = true;
        SetCapture(canvas_);
        return 0;
    case WM_RBUTTONUP:
    case WM_MBUTTONUP:
        if (drag_ == Drag::orbit || drag_ == Drag::pan) on_up(x, y);
        return 0;
    case WM_MOUSEWHEEL: {
        // Slice 4: the wheel sets the doping tools' level by decades, a refine box's size by halves.
        if (mode_ != modeller::Mode::none && modeller::group_of(mode_) == modeller::Group::doping) {
            level_ = std::clamp(level_ * (GET_WHEEL_DELTA_WPARAM(w) > 0 ? 10.0 : 0.1), 1e12, 1e22);
            level_ = std::pow(10.0, std::round(std::log10(level_)));
            host_.solid_status(std::format("{} {:.0e} cm^-3 (D: the dopant; the wheel: the level)", donor_ ? "donor" : "acceptor", level_));
            if (drag_ != Drag::none) on_move(last_x_, last_y_, false, false);
            return 0;
        }
        if (drag_ == Drag::refine) {
            refine_size_ = std::max<Coord>(1, GET_WHEEL_DELTA_WPARAM(w) > 0 ? refine_size_ * 2 : refine_size_ / 2);
            on_move(last_x_, last_y_, false, false);
            return 0;
        }
        if (drag_ == Drag::array) {  // MT5: the wheel sets the count, whatever the camera
            array_count_ = std::clamp(array_count_ + (GET_WHEEL_DELTA_WPARAM(w) > 0 ? 1 : -1), 2, 100);
            on_move(last_x_, last_y_, false, false);
            after_drag_step(last_x_, last_y_);
            return 0;
        }
        if (!modeller::wheel_zooms(camera_fixed_)) return 0;
        POINT p{GET_X_LPARAM(l), GET_Y_LPARAM(l)};
        ScreenToClient(canvas_, &p);
        const double k = std::pow(1.2, GET_WHEEL_DELTA_WPARAM(w) / 120.0);
        const render::ViewFrame f = render::view_frame(frame_scene(m_, v_), v_.width, v_.height);
        const double nx = 2.0 * (p.x - f.left) / (f.right_px - f.left) - 1.0, ny = 1.0 - 2.0 * (p.y - f.top) / (f.bottom_px - f.top);
        // The point under the cursor stays: cx' = cx + nx (half_w - half_w'), pan' = (mid - cx') / half_w'.
        const double hw = f.half_w / k, hh = f.half_h / k;
        const double mid_x = f.cx + v_.camera.pan[0] * f.half_w, mid_y = f.cy + v_.camera.pan[1] * f.half_h;
        const double cx = f.cx + nx * (f.half_w - hw), cy = f.cy + ny * (f.half_h - hh);
        v_.camera.zoom = static_cast<float>(std::clamp(v_.camera.zoom * k, 1e-3, 1e7));
        v_.camera.pan = {static_cast<float>((mid_x - cx) / hw), static_cast<float>((mid_y - cy) / hh)};
        invalidate();
        return 0;
    }
    case WM_KEYDOWN:
        if (w == VK_SPACE) {
            space_ = true;
            return 0;
        }
        if (on_key(static_cast<UINT>(w), GetKeyState(VK_SHIFT) < 0, GetKeyState(VK_CONTROL) < 0)) return 0;
        break;
    case WM_KEYUP:
        if (w == VK_SPACE) space_ = false;
        return 0;
    case WM_CHAR:
        if (on_char(static_cast<wchar_t>(w))) after_drag_step(last_x_, last_y_);
        return 0;
    case WM_SETCURSOR:
        if (LOWORD(l) == HTCLIENT) {
            SetCursor(LoadCursorW(nullptr, tool_ == Tool::select ? IDC_ARROW : tool_ == Tool::pan || space_ ? IDC_SIZEALL : IDC_CROSS));
            return TRUE;
        }
        break;
    case WM_CAPTURECHANGED:
        if (reinterpret_cast<HWND>(l) != canvas_ && button_) cancel();
        return 0;
    default: break;
    }
    return DefWindowProcW(canvas_, message, w, l);
}

}  // namespace NiTCAD::app
