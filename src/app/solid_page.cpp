// The 3D modeller on the Geometry tab (solid_page.hpp; DECISIONS.md Unit 26l, SD4).
#include "solid_page.hpp"

#include <commctrl.h>
#include <windowsx.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <format>
#include <sstream>

#include "NiTCAD/app/deck_vocabulary.hpp"
#include "NiTCAD/render/axes3d.hpp"
#include "win32.hpp"

namespace NiTCAD::app {

namespace {

using namespace solid;

constexpr const wchar_t* canvas_class = L"NiTCADSolidCanvas";
constexpr const wchar_t* panel_class = L"NiTCADSolidPanel";
constexpr int id_strip = 1, id_canvas = 2, id_tool_first = 10, id_fit = 20, id_mesh = 21;
constexpr std::array<const wchar_t*, 6> tool_texts{L"Select (V)", L"Box (B)", L"Prism (P)", L"Contact (C)", L"Pan (H)", L"Measure (M)"};
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
    for (std::size_t k = 0; k < tool_texts.size(); ++k) {
        tool_buttons_.push_back(make(L"BUTTON", tool_texts[k], WS_TABSTOP | BS_AUTORADIOBUTTON | BS_PUSHLIKE | (k == 0 ? WS_GROUP : 0), strip_,
                                     id_tool_first + static_cast<int>(k)));
    }
    tool_buttons_.push_back(make(L"BUTTON", L"Fit (F)", WS_TABSTOP | WS_GROUP | BS_PUSHBUTTON, strip_, id_fit));
    mesh_button_ = make(L"BUTTON", L"Mesh", WS_TABSTOP | BS_AUTOCHECKBOX | BS_PUSHLIKE, strip_, id_mesh);
    tool_buttons_.push_back(mesh_button_);
    canvas_ = CreateWindowExW(0, canvas_class, L"", WS_CHILD | WS_VISIBLE | WS_TABSTOP, 0, 0, 10, 10, panel_,
                              reinterpret_cast<HMENU>(static_cast<INT_PTR>(id_canvas)), GetModuleHandleW(nullptr), this);

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
    for (HWND h : tool_buttons_) SendMessageW(h, WM_SETFONT, reinterpret_cast<WPARAM>(font), TRUE);
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
    int x = pad;
    for (std::size_t k = 0; k < tool_buttons_.size(); ++k) {
        const int w = MulDiv(k < tool_texts.size() ? 92 : 64, dpi_, 96);
        MoveWindow(tool_buttons_[k], x, pad, w, bh, TRUE);
        x += w + (k + 1 == tool_texts.size() ? 3 * pad : pad / 2);
    }
    MoveWindow(canvas_, 0, bh + 2 * pad, r.right, std::max<int>(1, r.bottom - bh - 2 * pad), TRUE);
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
    for (int k = 0; k < 6; ++k) pair(c_x0_label + 2 * k, c_x0 + 2 * k);
    pair(c_axis_label, c_axis, MulDiv(80, dpi_, 96));
    pair(c_from_label, c_from);
    pair(c_to_label, c_to);
    line(c_vertices_label);
    line(c_vertices, MulDiv(120, dpi_, 96));
    line(c_apply);
    pair(c_kind_label, c_kind, MulDiv(160, dpi_, 96));
    line(c_convert);
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
    if (auto m = read_model(text)) {
        m_ = std::move(*m);
        message_ = m_.mesh == MeshKind::unsupported ? m_.why : std::string();
    }
    prune_selection();
    locate_problems();
    update_inspector();
    invalidate();
}

bool SolidPage::commit_model(Model next) {
    const std::string text = host_.solid_deck();
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
    update_tool_buttons();
    invalidate();
}

void SolidPage::update_tool_buttons() {
    for (std::size_t k = 0; k < tool_texts.size(); ++k) {
        SendMessageW(tool_buttons_[k], BM_SETCHECK, static_cast<int>(k) == static_cast<int>(tool_) ? BST_CHECKED : BST_UNCHECKED, 0);
    }
}

void SolidPage::cancel() {
    if (drag_ == Drag::move || drag_ == Drag::axis || drag_ == Drag::plane || drag_ == Drag::face || drag_ == Drag::vertex ||
        drag_ == Drag::contact_move) {
        m_ = m_start_;
    }
    drag_ = Drag::none;
    button_ = false;
    polygon_.clear();
    preview_ = {};
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
    if (mesh_button_) SendMessageW(mesh_button_, BM_SETCHECK, on ? BST_CHECKED : BST_UNCHECKED, 0);
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
          c_to, c_vertices_label, c_vertices, c_apply, c_kind_label, c_kind, c_convert},
         false);
    for (int k = 0; k < 12; ++k) show({c_x0_label + k}, false);
    const bool usable = m_.mesh != MeshKind::unsupported;
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
            show({c_raise, c_lower}, true);
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
        SetCapture(canvas_);
    };
    if (space_ || tool_ == Tool::pan) {
        begin(Drag::pan);
        return;
    }
    // A height being given: the click ends the box or prism.
    if (drag_ == Drag::height) {
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
        if (r >= 0 && commit_model(std::move(next))) sel_ = {{{r, 0}}, -1};
        else host_.solid_status("That solid is not taken (it must have a size on each axis).");
        update_inspector();
        invalidate();
        return;
    }
    switch (tool_) {
    case Tool::select: {
        const Pick p = pick(m_, v_, x, y, sel_);
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
            begin(Drag::move);
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
    switch (drag_) {
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
        d[axis_] = quantize(along_axis(m_start_, v_, x, y, origin_, axis_) - start_s_);
        move_solid(next, ref_, d);
        with_drag(std::move(next));
        return;
    }
    case Drag::move:
    case Drag::plane:
    case Drag::contact_move: {
        const auto q = on_plane(m_start_, v_, x, y, axis_, start_p_[static_cast<std::size_t>(axis_)]);
        if (!q) return;
        P3 d{};
        for (int a = 0; a < 3; ++a) {
            if (a != axis_) d[a] = quantize((*q)[static_cast<std::size_t>(a)] - start_p_[static_cast<std::size_t>(a)]);
        }
        Model next = m_start_;
        if (drag_ == Drag::contact_move) {
            move_contact(next, sel_.contact, d);
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
        move_face(next, ref_, index_, quantize(along_axis(m_start_, v_, x, y, origin_, axis_) - start_s_));
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
        preview_.solids = {Solid{SolidKind::box, normalized(a_, b_), 2, 0, 0, {}}};
        invalidate();
        return;
    }
    case Drag::height: {
        const int n = v_.plane.normal;
        const V3 corner = um_of(b_);
        height_ = quantize(along_axis(m_, v_, x, y, corner, n) - corner[static_cast<std::size_t>(n)]);
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
        host_.solid_status(std::format("height {} µm", ui_length(height_)));
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

void SolidPage::on_up(double x, double y) {
    (void)x;
    (void)y;
    button_ = false;
    if (GetCapture() == canvas_) ReleaseCapture();
    switch (drag_) {
    case Drag::move:
    case Drag::axis:
    case Drag::plane:
    case Drag::face:
    case Drag::vertex:
    case Drag::contact_move: {
        Model done = m_;
        m_ = m_start_;
        drag_ = Drag::none;
        if (!(done == m_start_)) commit_model(std::move(done));
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
    switch (key) {
    case 'V': set_tool(Tool::select); return true;
    case 'B':
    case 'R': set_tool(Tool::box); return true;
    case 'P': set_tool(Tool::prism); return true;
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
            const int id = LOWORD(w);
            if (id >= id_tool_first && id < id_tool_first + static_cast<int>(tool_texts.size())) {
                self->set_tool(static_cast<Tool>(id - id_tool_first));
                SetFocus(self->canvas_);
                return 0;
            }
            if (id == id_mesh) {
                self->set_mesh_preview(SendMessageW(self->mesh_button_, BM_GETCHECK, 0, 0) == BST_CHECKED);
                SetFocus(self->canvas_);
                return 0;
            }
            if (id == id_fit) {
                self->fit();
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
            if ((code == BN_CLICKED && (c == c_raise || c == c_lower || c == c_apply || c == c_translucent || c == c_convert)) ||
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
    case WM_MOUSEMOVE: on_move(x, y, shift, ctrl); return 0;
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
        drag_ = message == WM_RBUTTONDOWN ? Drag::orbit : Drag::pan;
        button_ = true;
        SetCapture(canvas_);
        return 0;
    case WM_RBUTTONUP:
    case WM_MBUTTONUP:
        if (drag_ == Drag::orbit || drag_ == Drag::pan) on_up(x, y);
        return 0;
    case WM_MOUSEWHEEL: {
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
    case WM_CHAR: return 0;
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
