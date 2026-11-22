// The Geometry tab (geometry_page.hpp; DECISIONS.md Unit 26k, GD6-GD8).
#include "geometry_page.hpp"

#include <commctrl.h>
#include <windowsx.h>

#include <algorithm>
#include <array>
#include <tuple>
#include <cmath>
#include <format>
#include <sstream>

#include "NiTCAD/app/deck_vocabulary.hpp"
#include "NiTCAD/app/solid.hpp"
#include "win32.hpp"

// Last: they bring objbase.h's `interface` macro.
#include <d2d1.h>
#include <dwrite.h>

namespace NiTCAD::app {

namespace {

using namespace geo;

constexpr const wchar_t* canvas_class = L"NiTCADGeometryCanvas";
constexpr const wchar_t* panel_class = L"NiTCADGeometryPanel";
constexpr int id_strip = 1, id_canvas = 2, id_bar = 100;

template <class T>
void release(T*& p) {
    if (p) p->Release();
    p = nullptr;
}

HWND make(const wchar_t* cls, const wchar_t* text, DWORD style, HWND parent, int id, DWORD ex = 0) {
    return CreateWindowExW(ex, cls, text, style | WS_CHILD | WS_VISIBLE, 0, 0, 10, 10, parent,
                           reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)), GetModuleHandleW(nullptr), nullptr);
}

D2D1_COLOR_F rgb(float r, float g, float b, float a = 1.0f) { return D2D1::ColorF(r, g, b, a); }

// A material's fill: the device view's beige for silicon, pale blue for oxides, a hue from the name otherwise.
D2D1_COLOR_F material_fill(std::string_view m) {
    if (m == "silicon" || m == "silicon_legacy") return rgb(0.93f, 0.86f, 0.70f);
    if (m == "sio2") return rgb(0.78f, 0.87f, 0.97f);
    if (m == "insulator") return rgb(0.86f, 0.86f, 0.93f);
    if (m == "germanium") return rgb(0.80f, 0.86f, 0.64f);
    if (m == "gaas") return rgb(0.74f, 0.89f, 0.76f);
    if (m == "ingaas") return rgb(0.70f, 0.83f, 0.87f);
    if (m == "sic_4h") return rgb(0.86f, 0.80f, 0.92f);
    if (m == "algaas") return rgb(0.90f, 0.78f, 0.70f);
    std::size_t h = 1469598103934665603ull;
    for (const char c : m) h = (h ^ static_cast<unsigned char>(c)) * 1099511628211ull;
    const float hue = static_cast<float>(h % 360) / 60.0f;
    const auto f = [&](float n) {
        const float k = std::fmod(n + hue, 6.0f);
        return 0.88f - 0.18f * std::max(0.0f, std::min({k, 4.0f - k, 1.0f}));
    };
    return rgb(f(5), f(3), f(1));
}

D2D1_COLOR_F darker(D2D1_COLOR_F c, float by = 0.45f) { return rgb(c.r * by, c.g * by, c.b * by, 1.0f); }

D2D1_COLOR_F contact_color(std::string_view kind) {
    if (kind == "gate") return rgb(0.55f, 0.20f, 0.70f);
    if (kind == "schottky") return rgb(0.10f, 0.50f, 0.25f);
    if (kind == "electrode") return rgb(0.20f, 0.30f, 0.65f);
    return rgb(0.85f, 0.50f, 0.05f);
}

// A length for the inspector: µm without the unit ("0.4", "2", "-0.05").
std::string ui_length(Coord q) {
    std::string s = format_length(q);
    if (s.ends_with("um")) s.resize(s.size() - 2);
    return s;
}

// The inspector's text as a length: µm unless it names its unit.
std::optional<Coord> ui_parse(std::string s) {
    while (!s.empty() && (s.back() == ' ' || s.back() == '\t' || s.back() == '\r')) s.pop_back();
    if (s.empty()) return std::nullopt;
    if (!(s.ends_with("um") || s.ends_with("nm") || s.ends_with("cm"))) s += "um";
    return parse_length(s);
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

std::wstring wide(std::string_view s) { return win32::widen(s); }

// The starter decks of a device drawn on a deck with no mesh (can_start): the drawn box is the domain, silicon its
// background; the mesh sizes follow the box. The 3D one has a body contact on its bottom face, so it builds; the
// Geometry tab then shows the 3D modeller for it (Unit 26l).
std::string layout2d(Rect d) {
    const Coord w = d.hi.x - d.lo.x, h = d.hi.y - d.lo.y;
    const Coord size = std::max<Coord>(10, std::max(w, h) / 20), res = std::clamp<Coord>(std::min(w, h) / 10, 10, 100);
    return std::format(R"(# A 2D device drawn on the Geometry tab: the regions are layers, each painted over the ones above it in
# the deck (Unit 26k).
[mesh]
type = generate2d
layered = on
domain = {} {} {} {}
max_size = {}
interface_resolution = {}

[region substrate]
material = silicon

[doping]
profile = acceptor constant 1e16

[equilibrium]
)", format_length(d.lo.x), format_length(d.hi.x), format_length(d.lo.y), format_length(d.hi.y), format_length(size), format_length(res));
}

std::string layout3d(Rect d, Coord height) {
    const Coord w = d.hi.x - d.lo.x, h = d.hi.y - d.lo.y;
    const Coord size = std::max<Coord>(10, std::max({w, h, height}) / 8);
    const Coord res = std::clamp<Coord>(std::min({w, h, height}) / 10, 10, 100);
    const std::string x0 = format_length(d.lo.x), x1 = format_length(d.hi.x), y0 = format_length(d.lo.y), y1 = format_length(d.hi.y);
    return std::format(R"(# A 3D device drawn on the Geometry tab: the regions are layers, each solid replacing what it covers, in
# deck order (Unit 26l).
[mesh]
type = generate3d
layered = on
domain = {0} {1} {2} {3} 0 {4}
patch body = {0} {2} 0  {1} {3} 0
max_size = {5}
interface_resolution = {6}

[region substrate]
material = silicon

[doping]
profile = acceptor constant 1e16

[contact body]
kind = ohmic
patch = body

[equilibrium]
)", x0, x1, y0, y1, format_length(height), format_length(size), format_length(res));
}

}  // namespace

GeometryPage::GeometryPage(HWND parent, Host& host, int id, int inspector_id, bool software)
    : host_(host), inspector_id_(inspector_id), software_(software) {
    win32::register_class(panel_class, DefWindowProcW, GetSysColorBrush(COLOR_WINDOW));
    win32::register_class(canvas_class, canvas_proc, nullptr, LoadCursorW(nullptr, IDC_ARROW),
                          CS_HREDRAW | CS_VREDRAW | CS_DBLCLKS);
    panel_ = CreateWindowExW(WS_EX_CONTROLPARENT, panel_class, L"", WS_CHILD | WS_CLIPCHILDREN, 0, 0, 10, 10, parent,
                             reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)), GetModuleHandleW(nullptr), nullptr);
    SetWindowLongPtrW(panel_, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(this));
    SetWindowLongPtrW(panel_, GWLP_WNDPROC, reinterpret_cast<LONG_PTR>(panel_proc));
    strip_ = CreateWindowExW(WS_EX_CONTROLPARENT, panel_class, L"", WS_CHILD | WS_VISIBLE | WS_CLIPCHILDREN, 0, 0, 10, 10,
                             panel_, reinterpret_cast<HMENU>(static_cast<INT_PTR>(id_strip)), GetModuleHandleW(nullptr),
                             nullptr);
    SetWindowLongPtrW(strip_, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(this));
    SetWindowLongPtrW(strip_, GWLP_WNDPROC, reinterpret_cast<LONG_PTR>(panel_proc));
    canvas_ = CreateWindowExW(0, canvas_class, L"", WS_CHILD | WS_VISIBLE | WS_TABSTOP | WS_CLIPSIBLINGS, 0, 0, 10, 10, panel_,
                              reinterpret_cast<HMENU>(static_cast<INT_PTR>(id_canvas)), GetModuleHandleW(nullptr), this);
    bar_ = std::make_unique<ModellerBar>(strip_, panel_, id_bar, false, panel_proc, this);

    // The Geometry section, in the inspector's place.
    inspector_ = CreateWindowExW(WS_EX_CONTROLPARENT, panel_class, L"", WS_CHILD | WS_CLIPCHILDREN | WS_VSCROLL, 0, 0, 10, 10,
                                 parent, reinterpret_cast<HMENU>(static_cast<INT_PTR>(inspector_id - 1)),
                                 GetModuleHandleW(nullptr), nullptr);
    SetWindowLongPtrW(inspector_, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(this));
    SetWindowLongPtrW(inspector_, GWLP_WNDPROC, reinterpret_cast<LONG_PTR>(inspector_proc));
    controls_.resize(c_count);
    const auto add = [&](int c, const wchar_t* cls, const wchar_t* text, DWORD style, DWORD ex = 0) {
        controls_[static_cast<std::size_t>(c)] = make(cls, text, style, inspector_, inspector_id + c, ex);
    };
    add(c_heading, L"STATIC", L"GEOMETRY", SS_LEFT);
    add(c_hint, L"STATIC", L"", SS_LEFT);
    add(c_name_label, L"STATIC", L"Name", SS_LEFT);
    add(c_name, L"EDIT", L"", WS_TABSTOP | ES_AUTOHSCROLL, WS_EX_CLIENTEDGE);
    add(c_material_label, L"STATIC", L"Material", SS_LEFT);
    add(c_material, WC_COMBOBOXW, L"", WS_TABSTOP | WS_VSCROLL | CBS_DROPDOWNLIST);
    add(c_raise, L"BUTTON", L"Raise layer", WS_TABSTOP | BS_PUSHBUTTON);
    add(c_lower, L"BUTTON", L"Lower layer", WS_TABSTOP | BS_PUSHBUTTON);
    add(c_x_label, L"STATIC", L"x (µm)", SS_LEFT);
    add(c_x, L"EDIT", L"", WS_TABSTOP | ES_AUTOHSCROLL, WS_EX_CLIENTEDGE);
    add(c_y_label, L"STATIC", L"y (µm)", SS_LEFT);
    add(c_y, L"EDIT", L"", WS_TABSTOP | ES_AUTOHSCROLL, WS_EX_CLIENTEDGE);
    add(c_w_label, L"STATIC", L"Width", SS_LEFT);
    add(c_w, L"EDIT", L"", WS_TABSTOP | ES_AUTOHSCROLL, WS_EX_CLIENTEDGE);
    add(c_h_label, L"STATIC", L"Height", SS_LEFT);
    add(c_h, L"EDIT", L"", WS_TABSTOP | ES_AUTOHSCROLL, WS_EX_CLIENTEDGE);
    add(c_vertices_label, L"STATIC", L"Vertices (x y, µm)", SS_LEFT);
    add(c_vertices, L"EDIT", L"", WS_TABSTOP | WS_VSCROLL | ES_MULTILINE | ES_AUTOVSCROLL | ES_WANTRETURN, WS_EX_CLIENTEDGE);
    add(c_apply, L"BUTTON", L"Apply", WS_TABSTOP | BS_PUSHBUTTON);
    add(c_kind_label, L"STATIC", L"Kind", SS_LEFT);
    add(c_kind, WC_COMBOBOXW, L"", WS_TABSTOP | WS_VSCROLL | CBS_DROPDOWNLIST);
    add(c_draw_label, L"STATIC", L"Draw with", SS_LEFT);
    add(c_draw, WC_COMBOBOXW, L"", WS_TABSTOP | WS_VSCROLL | CBS_DROPDOWNLIST);
    add(c_overlap, L"BUTTON", L"Fills only empty space (overlap = fill)", WS_TABSTOP | BS_AUTOCHECKBOX);
    for (const int c : {c_material, c_draw}) {
        for (const std::string& m : choices_of("region", "material")) {
            SendMessageW(control(c), CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(wide(m).c_str()));
        }
    }
    for (const std::string& k : choices_of("contact", "kind")) {
        SendMessageW(control(c_kind), CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(wide(k).c_str()));
    }
    SendMessageW(control(c_draw), CB_SELECTSTRING, static_cast<WPARAM>(-1), reinterpret_cast<LPARAM>(L"silicon"));

    D2D1_FACTORY_OPTIONS fo{};
    D2D1CreateFactory(D2D1_FACTORY_TYPE_SINGLE_THREADED, __uuidof(ID2D1Factory), &fo, reinterpret_cast<void**>(&d2d_));
    DWriteCreateFactory(DWRITE_FACTORY_TYPE_SHARED, __uuidof(IDWriteFactory), reinterpret_cast<IUnknown**>(&dwrite_));
    dims_host_.page = this;  // Unit 26n (DT1)
    dims_ = std::make_unique<DimensionsPanel>(inspector_, dims_host_, inspector_id + 900);
    ShowWindow(dims_->hwnd(), SW_SHOWNA);
    update_tool_buttons();
    update_inspector();
}

GeometryPage::~GeometryPage() {
    release_target();
    release(text_);
    release(dwrite_);
    release(d2d_);
    for (HWND h : {inspector_, panel_}) {
        if (h) {
            SetWindowLongPtrW(h, GWLP_USERDATA, 0);
            DestroyWindow(h);
        }
    }
}

HWND GeometryPage::control(int c) const {
    return c >= 0 && c < c_count ? controls_[static_cast<std::size_t>(c)] : nullptr;
}

void GeometryPage::set_dpi(int dpi, HFONT font, HFONT bold) {
    dpi_ = dpi;
    font_ = font;
    bold_ = bold;
    view_.dpi = dpi;
    bar_->set_font(font, dpi);
    for (HWND h : controls_) SendMessageW(h, WM_SETFONT, reinterpret_cast<WPARAM>(font), TRUE);
    SendMessageW(control(c_heading), WM_SETFONT, reinterpret_cast<WPARAM>(bold), TRUE);
    dims_->set_dpi(dpi, font, bold);
    release(text_);
    layout();
    layout_inspector();
}

double GeometryPage::ruler_px() const noexcept { return 20.0 * dpi_ / 96.0; }

void GeometryPage::layout() {
    RECT r{};
    GetClientRect(panel_, &r);
    const int pad = MulDiv(4, dpi_, 96), bh = MulDiv(24, dpi_, 96);
    MoveWindow(strip_, 0, 0, r.right, bh + 2 * pad, TRUE);
    bar_->layout_strip(r.right, pad, bh);
    canvas_top_ = bh + 2 * pad;
    MoveWindow(canvas_, 0, bh + 2 * pad, r.right, std::max<int>(1, r.bottom - bh - 2 * pad), TRUE);
    bar_->raise();
    place_bar();
    if (target_) {
        RECT c{};
        GetClientRect(canvas_, &c);
        target_->Resize(D2D1::SizeU(static_cast<UINT32>(std::max<LONG>(1, c.right)), static_cast<UINT32>(std::max<LONG>(1, c.bottom))));
    }
    if (!fitted_) fit();
    invalidate();
}

void GeometryPage::layout_inspector() {
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
    line(c_hint, 2 * (row - 6));
    pair(c_name_label, c_name);
    pair(c_material_label, c_material, MulDiv(200, dpi_, 96));
    if (visible(c_raise)) {
        MoveWindow(control(c_raise), pad, y, w / 2 - 2, row - 4, TRUE);
        MoveWindow(control(c_lower), pad + w / 2 + 2, y, w - w / 2 - 2, row - 4, TRUE);
        y += row;
    }
    line(c_overlap);
    pair(c_x_label, c_x);
    pair(c_y_label, c_y);
    pair(c_w_label, c_w);
    pair(c_h_label, c_h);
    line(c_vertices_label);
    line(c_vertices, MulDiv(120, dpi_, 96));
    line(c_apply);
    pair(c_kind_label, c_kind, MulDiv(160, dpi_, 96));
    y += pad;
    pair(c_draw_label, c_draw, MulDiv(200, dpi_, 96));
    // Unit 26n (DT1): the dimensions panel below.
    y += pad;
    MoveWindow(dims_->hwnd(), 0, y, r.right, std::max(1, dims_->preferred_height(r.right)), TRUE);
}

// --- the deck -------------------------------------------------------------------------------------------------

void GeometryPage::reload(bool fit_view, bool force) {
    const std::string text = host_.geometry_deck();
    if (text == text_seen_ && !fit_view && !force && g_.editable()) {
        update_inspector();
        invalidate();
        return;
    }
    if (text != text_seen_ && !history_.can_undo() && !history_.can_redo()) {
        // nothing to keep
    } else if (text != text_seen_) {
        history_.clear();  // edited elsewhere: the page's undo starts again
    }
    const bool first = text_seen_.empty() && !g_.editable();
    text_seen_ = text;
    auto g = read_page(text);
    if (!g) {
        g_ = {};
        message_ = "The deck does not parse: " + g.error();
    } else {
        g_ = std::move(*g);
        message_ = g_.editable() ? std::string() : g_.why;
    }
    prune_selection();
    if (fit_view || first) fit();
    locate_problems();
    reload_dimensions();
    update_inspector();
    invalidate();
}

// The selection kept to what the model has (the deck read again: an edit elsewhere, an undo, a redo).
void GeometryPage::prune_selection() {
    const auto valid = [&](ShapeRef s) {
        return s.region >= 0 && s.region < static_cast<int>(g_.regions.size()) &&
               (s.shape == -1 ? g_.regions[static_cast<std::size_t>(s.region)].shapes.empty()
                              : s.shape >= 0 && s.shape < static_cast<int>(g_.regions[static_cast<std::size_t>(s.region)].shapes.size()));
    };
    std::erase_if(selection_, [&](ShapeRef s) { return !valid(s); });
    if (contact_ >= static_cast<int>(g_.contacts.size())) contact_ = -1;
}

void GeometryPage::apply_text(const std::string& text) {
    host_.geometry_edited(text);
    text_seen_ = text;
    if (auto g = read_page(text)) {
        g_ = std::move(*g);
        message_ = g_.editable() ? std::string() : g_.why;
    }
    prune_selection();  // an undo can take the selected shape or contact away
    locate_problems();
    reload_dimensions();
    update_inspector();
    invalidate();
}

// --- live validation and the mesh preview (slice 3, GD9) ----------------------------------------------------------

void GeometryPage::checked(const std::string& text, std::span<const Diagnostic> diagnostics) {
    checked_text_ = text;
    diagnostics_.assign(diagnostics.begin(), diagnostics.end());
    locate_problems();
    invalidate();
}

void GeometryPage::locate_problems() {
    problems_.clear();
    if (checked_text_.empty() || checked_text_ != text_seen_ || !g_.editable()) return;
    for (const Diagnostic& d : diagnostics_) problems_.push_back(locate_problem(g_, text_seen_, d.line, d.message));
}

GeometryPage::Check GeometryPage::check() const noexcept {
    if (checked_text_.empty() || checked_text_ != text_seen_) return Check::pending;
    return diagnostics_.empty() ? Check::valid : Check::invalid;
}

std::string GeometryPage::run_blocked() const {
    if (check() != Check::invalid || !g_.editable()) return {};
    return "Run is off until the layout builds: " + diagnostics_.front().message;
}

void GeometryPage::set_mesh(const std::string& text, MeshOverlay mesh) {
    mesh_text_ = text;
    mesh_ = std::move(mesh);
    if (mesh_preview_) invalidate();
}

void GeometryPage::set_mesh_preview(bool on) {
    mesh_preview_ = on;
    bar_->set_mesh(on);
    invalidate();
}

std::string GeometryPage::banner() const {
    if (!g_.editable()) return {};
    if (check() == Check::invalid) {
        return diagnostics_.front().message + "  (Run is off until the layout builds.)";
    }
    if (mesh_preview_ && !mesh_shown()) return check() == Check::pending ? "Building the mesh..." : "The mesh shows when the layout builds.";
    return {};
}

bool GeometryPage::commit_model(Geometry next) {
    if (view_only_) return false;  // Unit 26n: dimensions only
    const std::string text = host_.geometry_deck();
    auto out = write_geometry(text, next);
    if (!out) {
        host_.geometry_status(out.error());
        return false;
    }
    if (*out == text) {
        g_ = std::move(next);
        invalidate();
        return false;
    }
    history_.record(text);
    apply_text(*out);
    return true;
}

bool GeometryPage::undo() {
    auto t = history_.undo(host_.geometry_deck());
    if (!t) return false;
    apply_text(*t);
    return true;
}

bool GeometryPage::redo() {
    auto t = history_.redo(host_.geometry_deck());
    if (!t) return false;
    apply_text(*t);
    return true;
}

bool GeometryPage::copy() {
    if (selection_.size() != 1 || selection_[0].shape < 0) return false;
    const ShapeRef s = selection_[0];
    clipboard_ = g_.regions[static_cast<std::size_t>(s.region)].shapes[static_cast<std::size_t>(s.shape)];
    clipboard_material_ = g_.regions[static_cast<std::size_t>(s.region)].material;
    return true;
}

bool GeometryPage::paste() {
    if (!clipboard_ || !g_.editable()) return false;
    Geometry next = g_;
    const Coord step = grid_step(view_);
    Shape s = *clipboard_;
    for (Pt& p : s.points) p = {p.x + step, p.y - step};
    const int r = s.kind == ShapeKind::box ? add_box(next, {s.points[0], s.points[1]}, clipboard_material_)
                                           : add_polygon(next, s.points, clipboard_material_);
    if (r < 0 || !commit_model(std::move(next))) return false;
    selection_ = {{r, 0}};
    clipboard_ = s;
    update_inspector();
    return true;
}

bool GeometryPage::duplicate() {
    if (selection_.size() != 1 || selection_[0].shape < 0) return false;
    Geometry next = g_;
    const Coord step = grid_step(view_);
    const int r = duplicate_shape(next, selection_[0], step, -step);
    if (r < 0 || !commit_model(std::move(next))) return false;
    selection_ = {{r, 0}};
    update_inspector();
    return true;
}

bool GeometryPage::remove() {
    Geometry next = g_;
    bool any = false;
    if (contact_ >= 0) {
        any = delete_contact(next, contact_);
        contact_ = -1;
    } else {
        // Later shapes first, so the earlier references stay right.
        auto sel = selection_;
        std::ranges::sort(sel, [](ShapeRef a, ShapeRef b) { return std::tie(a.region, a.shape) > std::tie(b.region, b.shape); });
        for (const ShapeRef s : sel) {
            if (s.shape >= 0) any = delete_shape(next, s) || any;
        }
    }
    selection_.clear();
    if (!any) return false;
    const bool done = commit_model(std::move(next));
    update_inspector();
    return done;
}

bool GeometryPage::select_region(const std::string& name) {
    for (std::size_t r = 0; r < g_.regions.size(); ++r) {
        if (g_.regions[r].name != name) continue;
        selection_ = {{static_cast<int>(r), g_.regions[r].shapes.empty() ? -1 : 0}};
        contact_ = -1;
        update_inspector();
        invalidate();
        return true;
    }
    return false;
}

void GeometryPage::set_view(const geo::View& v) {
    view_ = v;
    view_.dpi = dpi_;
    fitted_ = true;
    invalidate();
}

void GeometryPage::set_tool(Tool t) {
    cancel();
    tool_ = t;
    mode_ = modeller::Mode::none;
    update_tool_buttons();
    invalidate();
}

void GeometryPage::update_tool_buttons() {
    // A tool armed by its key opens a group that holds it (Pan and Measure are in none: the group stays).
    const auto in = [&](modeller::Group g) {
        return std::ranges::any_of(modeller::tools_for(g, false), [&](const modeller::ToolEntry& e) {
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

void GeometryPage::set_camera_fixed(bool fixed) {
    camera_fixed_ = fixed;
    bar_->set_camera_fixed(fixed);
    if (fixed && drag_ == Drag::pan) cancel();
}

void GeometryPage::set_group(modeller::Group g) {
    bar_->set_group(g);
    if (mode_ != modeller::Mode::none && modeller::group_of(mode_) != g) mode_ = modeller::Mode::none;
    const auto list = modeller::tools_for(g, false);
    const bool has = std::ranges::any_of(list, [&](const modeller::ToolEntry& e) { return e.action == modeller::Action::tool && e.tool == static_cast<int>(tool_); });
    if (!has && !list.empty() && list.front().action == modeller::Action::tool) {
        set_tool(static_cast<Tool>(list.front().tool));
    } else {
        bar_->set_tool(modeller::armed_key(static_cast<int>(tool_), mode_));
    }
    place_bar();
}


// --- Unit 26m slice 1: the overlap rule and the Transform modes (MT4, MT5) -----------------------------------------

void GeometryPage::set_mode(modeller::Mode m) {
    if (tool_ != Tool::select) set_tool(Tool::select);
    cancel();
    mode_ = m;
    array_count_ = 3;
    if (m != modeller::Mode::none && bar_->group() != modeller::group_of(m)) bar_->set_group(modeller::group_of(m));
    bar_->set_tool(modeller::armed_key(static_cast<int>(tool_), mode_));
    place_bar();
    invalidate();
}

bool GeometryPage::toggle_overlap() {
    std::vector<int> regions;
    for (const ShapeRef s : selection_) {
        if (s.shape >= 0 && std::ranges::find(regions, s.region) == regions.end()) regions.push_back(s.region);
    }
    if (regions.empty()) {
        host_.geometry_status("Pick a shape: its region then fills only empty space, or replaces again.");
        return false;
    }
    Geometry next = g_;
    // One rule for them all: fill unless every one fills already.
    const bool fill = !std::ranges::all_of(regions, [&](int r) { return g_.regions[static_cast<std::size_t>(r)].fill; });
    bool any = false;
    for (const int r : regions) any = set_overlap(next, r, fill) || any;
    if (!any) {
        host_.geometry_status("A tensor deck's regions only replace (overlap is for layered meshes).");
        return false;
    }
    if (!commit_model(std::move(next))) return false;
    update_inspector();
    return true;
}

bool GeometryPage::transform_click(const Hit& hit, Pt at, bool shift, bool ctrl) {
    const bool on_shape = (hit.kind == HitKind::shape || hit.kind == HitKind::edge || hit.kind == HitKind::vertex || hit.kind == HitKind::handle) &&
                          hit.shape.shape >= 0;
    Geometry next = g_;
    if (mode_ == modeller::Mode::mirror) {
        if (!on_shape) {
            host_.geometry_status("Mirror: click a shape.");
            return false;
        }
        const Rect b = g_.regions[static_cast<std::size_t>(hit.shape.region)].shapes[static_cast<std::size_t>(hit.shape.shape)].bounds();
        bool ok = false;
        if (shift) {
            ok = rotate_shape(next, hit.shape);
        } else if (hit.kind == HitKind::handle && hit.index % 2 == 1) {
            // A side's handle (1 bottom, 3 right, 5 top, 7 left): across that side.
            const int h = hit.index;
            ok = h == 1 || h == 5 ? mirror_shape_at(next, hit.shape, false, 2 * (h == 1 ? b.lo.y : b.hi.y))
                                  : mirror_shape_at(next, hit.shape, true, 2 * (h == 3 ? b.hi.x : b.lo.x));
        } else {
            ok = ctrl ? mirror_shape_at(next, hit.shape, false, b.lo.y + b.hi.y) : mirror_shape_at(next, hit.shape, true, b.lo.x + b.hi.x);
        }
        if (!ok) {
            host_.geometry_status("That mirror leaves the shape outside the limits.");
            return false;
        }
        selection_ = {hit.shape};
        contact_ = -1;
        const bool done = commit_model(std::move(next));  // false too when nothing changes
        update_inspector();
        return done;
    }
    // Align: the side of the clicked shape's bounds (or the domain's) nearest the click; the other picked shapes'
    // same sides go to its line.
    const Rect b = on_shape ? g_.regions[static_cast<std::size_t>(hit.shape.region)].shapes[static_cast<std::size_t>(hit.shape.shape)].bounds()
                            : g_.domain;
    const std::array<Coord, 4> gap{std::abs(at.x - b.lo.x), std::abs(at.x - b.hi.x), std::abs(at.y - b.lo.y), std::abs(at.y - b.hi.y)};
    const int side = static_cast<int>(std::ranges::min_element(gap) - gap.begin());
    const Coord line = side == 0 ? b.lo.x : side == 1 ? b.hi.x : side == 2 ? b.lo.y : b.hi.y;
    std::vector<ShapeRef> others;
    for (const ShapeRef s : selection_) {
        if (s.shape >= 0 && !(on_shape && s == hit.shape)) others.push_back(s);
    }
    if (others.empty()) {
        host_.geometry_status("Align: pick the shapes to move first (Ctrl+click), then click near the side to align them to.");
        return false;
    }
    if (align_shapes(next, others, side, line) && commit_model(std::move(next))) {
        update_inspector();
        return true;
    }
    host_.geometry_status("Align: nothing moves (the sides are there already).");
    return false;
}
std::optional<modeller::Box> GeometryPage::selection_bounds() const {
    std::optional<Rect> all;
    const auto add = [&](Rect r) {
        if (!all) all = r;
        else all = Rect{{std::min(all->lo.x, r.lo.x), std::min(all->lo.y, r.lo.y)}, {std::max(all->hi.x, r.hi.x), std::max(all->hi.y, r.hi.y)}};
    };
    for (const ShapeRef s : selection_) {
        if (s.region < 0 || s.region >= static_cast<int>(g_.regions.size())) continue;
        const auto& shapes = g_.regions[static_cast<std::size_t>(s.region)].shapes;
        if (s.shape >= 0 && s.shape < static_cast<int>(shapes.size())) add(shapes[static_cast<std::size_t>(s.shape)].bounds());
        else if (s.shape < 0) add(g_.domain);
    }
    if (contact_ >= 0 && contact_ < static_cast<int>(g_.contacts.size())) {
        const Contact& c = g_.contacts[static_cast<std::size_t>(contact_)];
        add(normalized(c.a, c.b));
    }
    if (!all) return std::nullopt;
    return modeller::Box{view_.sx(all->lo.x), view_.sy(all->hi.y), view_.sx(all->hi.x), view_.sy(all->lo.y)};
}

void GeometryPage::place_bar() {
    RECT c{};
    GetClientRect(canvas_, &c);
    const double top = canvas_top_, ruler = ruler_px();
    const modeller::Box viewport{ruler, top + ruler, static_cast<double>(c.right), top + c.bottom};
    std::optional<modeller::Box> s = selection_bounds();
    if (s) *s = {s->left, s->top + top, s->right, s->bottom + top};
    bar_->place(viewport, s);
}

modeller::Placement GeometryPage::bar_placement() const {
    modeller::Placement p = bar_->placement();
    p.bar.top -= canvas_top_;
    p.bar.bottom -= canvas_top_;
    return p;
}

std::string GeometryPage::fresh_contact() const {
    for (int n = 1;; ++n) {
        std::string name = "contact_" + std::to_string(n);
        // Not a contact's name, nor a patch's left in the deck (a renamed contact keeps its patch).
        if (std::ranges::none_of(g_.contacts, [&](const Contact& c) { return c.name == name; }) &&
            text_seen_.find("patch " + name + " ") == std::string::npos) {
            return name;
        }
    }
}

std::optional<Geometry> GeometryPage::pending() const {
    switch (drag_) {
    case Drag::corner:
    case Drag::move:
    case Drag::array:
    case Drag::handle:
    case Drag::vertex:
    case Drag::domain:
    case Drag::contact_end: return g_;
    case Drag::rectangle: {
        Geometry next = g_;
        if (add_box(next, normalized(anchor_, at_), material_) < 0) return std::nullopt;
        return next;
    }
    case Drag::contact: {
        Geometry next = g_;
        if (add_contact(next, fresh_contact(), "ohmic", anchor_, at_) < 0) return std::nullopt;
        return next;
    }
    default: return std::nullopt;
    }
}

std::string GeometryPage::chip_value() const {
    using modeller::chip_length;
    const auto offset = [](Coord dx, Coord dy) {
        std::string out;
        if (dx != 0) out = "\xCE\x94x " + chip_length(dx);
        if (dy != 0) out += (out.empty() ? "" : ", ") + std::string("\xCE\x94y ") + chip_length(dy);
        return out.empty() ? std::string("\xCE\x94 0 \xC2\xB5m") : out;
    };
    std::string v;
    switch (drag_) {
    case Drag::move:
    case Drag::vertex:
    case Drag::contact_end: v = offset(drag_d_.x, drag_d_.y); break;
    case Drag::array: v = offset(drag_d_.x, drag_d_.y) + std::format(" × {} (wheel)", array_count_); break;
    case Drag::handle:
    case Drag::domain: {
        Rect r = g_.domain;
        if (drag_ == Drag::handle && selection_.size() == 1 && selection_[0].shape >= 0) {
            r = g_.regions[static_cast<std::size_t>(selection_[0].region)].shapes[static_cast<std::size_t>(selection_[0].shape)].bounds();
        }
        v = chip_length(r.hi.x - r.lo.x) + " \xC3\x97 " + chip_length(r.hi.y - r.lo.y);
        break;
    }
    case Drag::rectangle: {
        const Rect r = normalized(anchor_, at_);
        v = chip_length(r.hi.x - r.lo.x) + " \xC3\x97 " + chip_length(r.hi.y - r.lo.y);
        break;
    }
    case Drag::contact: v = chip_length(std::abs(at_.x - anchor_.x) + std::abs(at_.y - anchor_.y)); break;
    case Drag::corner:
        v = (mode_ == modeller::Mode::chamfer ? "cut " : "radius ") + chip_length(corner_r_);
        if (mode_ != modeller::Mode::chamfer && corner_r_ > 0) v += corner_chords_ > 0 ? std::format(", {} chords", corner_chords_) : std::string(", chords under 6 s");
        break;
    case Drag::measure: {
        const double dx = static_cast<double>(at_.x - anchor_.x), dy = static_cast<double>(at_.y - anchor_.y);
        v = std::format("{:.6g} \xC2\xB5m", std::sqrt(dx * dx + dy * dy) / static_cast<double>(per_um));
        break;
    }
    default: return {};
    }
    if (!entry_.empty()) v += "   [" + entry_.text() + (entry_.value() ? "]" : "?]");
    return v;
}

void GeometryPage::after_drag_step(double x, double y) {
    if (const auto next = pending()) {
        auto out = write_geometry(text_seen_, *next);
        live_line_ = out ? modeller::deck_line(text_seen_, *out) : out.error();
        host_.geometry_status(live_line_);
    } else {
        live_line_.clear();
    }
    const int off = MulDiv(16, dpi_, 96);
    bar_->chip(static_cast<int>(x) + off, static_cast<int>(y) + canvas_top_ + off, chip_value());
}

void GeometryPage::end_drag_ui() {
    entry_.clear();
    typed_axis_ = -1;
    live_line_.clear();
    bar_->chip(0, 0, "");
}

Coord GeometryPage::typed_or(Coord dragged) const {
    const auto v = entry_.value();
    if (!v) return dragged;
    return *v >= 0 && dragged < 0 ? -*v : *v;
}

bool GeometryPage::on_char(wchar_t c) {
    if (drag_ != Drag::move && drag_ != Drag::rectangle && drag_ != Drag::array && drag_ != Drag::corner) return false;
    if (c == L'\b') {
        if (!entry_.backspace()) return true;
    } else if (!entry_.type(c)) {
        return false;
    }
    if (typed_axis_ < 0) {
        const Pt d = drag_ == Drag::move || drag_ == Drag::array ? drag_d_ : Pt{at_.x - anchor_.x, at_.y - anchor_.y};
        typed_axis_ = std::abs(d.y) > std::abs(d.x) ? 1 : 0;
    }
    on_move(last_x_, last_y_, false);
    return true;
}

void GeometryPage::cancel() {
    if (drag_ != Drag::none && drag_ != Drag::measure) g_ = drag_ == Drag::rubber || drag_ == Drag::pan ? g_ : g_start_;
    if (copy_from_) {
        g_ = *copy_from_;
        selection_ = copy_sel_;
        copy_from_.reset();
    }
    end_drag_ui();
    drag_ = Drag::none;
    start_ = Start::none;
    polygon_.clear();
    hover_.reset();
    if (GetCapture() == canvas_) ReleaseCapture();
    invalidate();
}

// --- a new device drawn on a deck with no mesh ------------------------------------------------------------------

bool GeometryPage::can_start() const { return !g_.editable() && message_.starts_with("the deck has no [mesh]"); }

Pt GeometryPage::start_point(double x, double y) const {
    const Coord step = grid_step(view_);
    const Pt p = view_.model(x, y);
    const auto round = [&](Coord v) { return static_cast<Coord>(std::llround(static_cast<double>(v) / static_cast<double>(step))) * step; };
    return {round(p.x), round(p.y)};
}

std::string GeometryPage::start_readout() const {
    if (start_ == Start::none) return {};
    const auto um = [](Coord q) { return std::format("{:g}", static_cast<double>(q) / per_um); };
    const Coord w = std::abs(start_b_.x - start_a_.x), h = std::abs(start_b_.y - start_a_.y);
    if (start_ == Start::footprint) return std::format("{} × {} µm", um(w), um(h));
    return std::format("{} × {} × {} µm", um(w), um(h), um(start_h_));
}

void GeometryPage::start_finish() {
    const Rect d{{std::min(start_a_.x, start_b_.x), std::min(start_a_.y, start_b_.y)},
                 {std::max(start_a_.x, start_b_.x), std::max(start_a_.y, start_b_.y)}};
    const bool three = start_ == Start::height;
    const Coord height = start_h_;
    start_ = Start::none;
    if (GetCapture() == canvas_) ReleaseCapture();
    if (three) {
        apply_text(layout3d(d, height));  // the host shows the 3D modeller for it
        return;
    }
    apply_text(layout2d(d));
    fit();
}

void GeometryPage::fit() {
    RECT c{};
    GetClientRect(canvas_, &c);
    const double rp = ruler_px();
    const double w = std::max(1.0, c.right - rp), h = std::max(1.0, c.bottom - rp);
    Rect d = g_.domain;
    if (d.lo == d.hi) d = {{0, 0}, {2 * per_um, per_um}};
    const double dw = std::max(1e-6, static_cast<double>(d.hi.x - d.lo.x) / per_um);
    const double dh = std::max(1e-6, static_cast<double>(d.hi.y - d.lo.y) / per_um);
    const double s = std::clamp(0.88 * std::min(w / dw, h / dh), 1e-2, 1e6);
    view_.px_per_um = s;
    view_.dpi = dpi_;
    view_.origin_x = rp + (w - dw * s) / 2.0 - s * static_cast<double>(d.lo.x) / per_um;
    view_.origin_y = rp + (h + dh * s) / 2.0 + s * static_cast<double>(d.lo.y) / per_um;
    fitted_ = c.right > 0;
    invalidate();
}

// --- the inspector --------------------------------------------------------------------------------------------

void GeometryPage::update_inspector() {
    const auto show = [&](std::initializer_list<int> cs, bool on) {
        for (const int c : cs) ShowWindow(control(c), on ? SW_SHOW : SW_HIDE);
    };
    const auto set = [&](int c, const std::string& t) { win32::set_text(control(c), t); };
    const auto pick = [&](int c, const std::string& t) {
        SendMessageW(control(c), CB_SELECTSTRING, static_cast<WPARAM>(-1), reinterpret_cast<LPARAM>(wide(t).c_str()));
    };
    show({c_name_label, c_name, c_material_label, c_material, c_raise, c_lower, c_x_label, c_x, c_y_label, c_y, c_w_label,
          c_w, c_h_label, c_h, c_vertices_label, c_vertices, c_apply, c_kind_label, c_kind, c_overlap},
         false);
    show({c_draw_label, c_draw}, g_.editable());
    const auto xywh = [&](Rect r, bool contact) {
        show({c_x_label, c_x, c_y_label, c_y, c_w_label, c_w, c_h_label, c_h}, true);
        win32::set_text(control(c_x_label), contact ? "x0 (µm)" : "x (µm)");
        win32::set_text(control(c_y_label), contact ? "y0 (µm)" : "y (µm)");
        win32::set_text(control(c_w_label), contact ? "x1 (µm)" : "Width");
        win32::set_text(control(c_h_label), contact ? "y1 (µm)" : "Height");
        set(c_x, ui_length(r.lo.x));
        set(c_y, ui_length(r.lo.y));
        set(c_w, ui_length(contact ? r.hi.x : r.hi.x - r.lo.x));
        set(c_h, ui_length(contact ? r.hi.y : r.hi.y - r.lo.y));
    };
    std::string hint;
    if (can_start()) {
        hint = "The deck has no mesh yet. Drag on the canvas to draw a 2D device; Ctrl+drag for a 3D one.";
    } else if (!g_.editable()) {
        hint = message_.empty() ? std::string("Open a 2D deck to draw its device.") : message_;
    } else if (contact_ >= 0) {
        const Contact& c = g_.contacts[static_cast<std::size_t>(contact_)];
        hint = "Contact " + c.name + (c.editable ? std::string() : " (not editable here)");
        show({c_name_label, c_name, c_kind_label, c_kind}, true);
        set(c_name, c.name);
        pick(c_kind, c.kind);
        if (c.editable) xywh({c.a, c.b}, true);
    } else if (selection_.size() == 1) {
        const ShapeRef s = selection_[0];
        const Region& r = g_.regions[static_cast<std::size_t>(s.region)];
        show({c_name_label, c_name, c_material_label, c_material}, true);
        set(c_name, r.name);
        pick(c_material, r.material);
        if (s.shape < 0) {
            hint = "Background " + r.name + ": the domain";
            if (g_.mesh == MeshKind::layered) xywh(g_.domain, false);
        } else {
            const Shape& sh = r.shapes[static_cast<std::size_t>(s.shape)];
            hint = std::format("Region {}, layer {} of {}: {}", r.name, s.region + 1, g_.regions.size(),
                               sh.kind == ShapeKind::box ? "a rectangle" : std::format("a polygon of {} vertices", sh.points.size()));
            show({c_raise, c_lower}, true);
            show({c_overlap}, g_.mesh != MeshKind::tensor);  // MT4
            SendMessageW(control(c_overlap), BM_SETCHECK, r.fill ? BST_CHECKED : BST_UNCHECKED, 0);
            if (sh.kind == ShapeKind::box) {
                xywh(normalized(sh.points[0], sh.points[1]), false);
            } else {
                show({c_vertices_label, c_vertices, c_apply}, true);
                std::string v;
                for (const Pt& p : sh.points) v += ui_length(p.x) + " " + ui_length(p.y) + "\r\n";
                set(c_vertices, v);
            }
        }
    } else if (selection_.size() > 1) {
        hint = std::format("{} shapes selected", selection_.size());
    } else {
        hint = "Click a shape to select it; draw with the tools above the canvas.";
    }
    set(c_hint, hint);
    layout_inspector();
}

bool GeometryPage::commit(int c) {
    if (c < 0 || c >= c_count) return false;
    const auto text = [&](int k) { return win32::text_of(control(k)); };
    const auto choice = [&](int k) {
        const auto i = SendMessageW(control(k), CB_GETCURSEL, 0, 0);
        if (i < 0) return std::string();
        wchar_t b[128]{};
        SendMessageW(control(k), CB_GETLBTEXT, static_cast<WPARAM>(i), reinterpret_cast<LPARAM>(b));
        return win32::narrow(b);
    };
    if (c == c_draw) {
        material_ = choice(c_draw);
        return true;
    }
    Geometry next = g_;
    bool ok = false;
    if (contact_ >= 0) {
        if (c == c_name) ok = rename_contact(next, contact_, text(c_name));
        else if (c == c_kind) ok = set_contact_kind(next, contact_, choice(c_kind));
        else if (c == c_x || c == c_y || c == c_w || c == c_h) {
            const auto x0 = ui_parse(text(c_x)), y0 = ui_parse(text(c_y)), x1 = ui_parse(text(c_w)), y1 = ui_parse(text(c_h));
            ok = x0 && y0 && x1 && y1 && set_contact(next, contact_, {*x0, *y0}, {*x1, *y1});
        }
    } else if (selection_.size() == 1) {
        const ShapeRef s = selection_[0];
        if (c == c_name) ok = rename_region(next, s.region, text(c_name));
        else if (c == c_material) ok = geo::set_material(next, s.region, choice(c_material));
        else if (c == c_overlap) ok = set_overlap(next, s.region, SendMessageW(control(c_overlap), BM_GETCHECK, 0, 0) == BST_CHECKED);
        else if (c == c_raise && raise_region(next, s.region)) {
            ok = true;
            selection_ = {{s.region + 1, s.shape}};
        } else if (c == c_lower && lower_region(next, s.region)) {
            ok = true;
            selection_ = {{s.region - 1, s.shape}};
        } else if (c == c_x || c == c_y || c == c_w || c == c_h) {
            const auto x = ui_parse(text(c_x)), y = ui_parse(text(c_y)), w = ui_parse(text(c_w)), h = ui_parse(text(c_h));
            if (x && y && w && h) {
                const Rect r{{*x, *y}, {*x + *w, *y + *h}};
                ok = s.shape < 0 ? set_domain(next, r) : set_box(next, s, r);
            }
        } else if ((c == c_apply || c == c_vertices) && s.shape >= 0) {
            std::vector<Pt> pts;
            std::istringstream in(text(c_vertices));
            bool good = true;
            for (std::string a, b; in >> a >> b;) {
                const auto x = ui_parse(a), y = ui_parse(b);
                if (!x || !y) good = false;
                else pts.push_back({*x, *y});
            }
            Shape& sh = next.regions[static_cast<std::size_t>(s.region)].shapes[static_cast<std::size_t>(s.shape)];
            if (good && sh.kind == ShapeKind::polygon && pts.size() >= 4) {
                Shape old = sh;
                sh.points = pts;
                // Kept to the editing rules: a probe edit of its first vertex to itself checks them.
                ok = set_vertex(next, s, 0, pts[0]);
                if (!ok) sh = old;
            }
        }
    }
    if (!ok) {
        host_.geometry_status("That value is not taken (a name must be new and start with a letter; a shape must keep its size, "
                              "a polygon its horizontal and vertical edges, a contact one side of the domain).");
        update_inspector();
        return true;
    }
    commit_model(std::move(next));
    update_inspector();
    invalidate();
    return true;
}

// --- the mouse and keys ---------------------------------------------------------------------------------------

Pt GeometryPage::snapped(double x, double y, bool shift, std::optional<ShapeRef> skip) const {
    Pt p = snap(g_, view_, view_.model(x, y), skip).p;
    if (shift && drag_ != Drag::none) p = constrain(anchor_, p);
    return p;
}

void GeometryPage::on_down(double x, double y, bool shift, bool ctrl, bool double_click) {
    SetFocus(canvas_);
    press_x_ = last_x_ = x;
    press_y_ = last_y_ = y;
    if (can_start()) {
        if (start_ == Start::height) {
            start_finish();  // the 3D device's height given: a click ends it
            return;
        }
        start_ = Start::footprint;
        start3d_ = ctrl;
        start_a_ = start_b_ = start_point(x, y);
        SetCapture(canvas_);
        invalidate();
        return;
    }
    if (!g_.editable()) return;
    if ((tool_ == Tool::select || view_only_) && dimension_down(x, y)) return;  // Unit 26n (DT2): a dimension line first
    const auto begin = [&](Drag d) {
        drag_ = d;
        g_start_ = g_;
        SetCapture(canvas_);
    };
    if (view_only_) {  // Unit 26n: a parametric deck's shapes move by their dimensions only
        if (!camera_fixed_) begin(Drag::pan);
        host_.geometry_status("This device is shaped by named dimensions: change them in the panel, or drag their lines.");
        return;
    }
    if (space_ || tool_ == Tool::pan) {
        // MT1: the Pan tool and Space+drag move the view only while the camera is free.
        if (!camera_fixed_) begin(Drag::pan);
        return;
    }
    entry_.clear();
    typed_axis_ = -1;
    drag_d_ = {};
    switch (tool_) {
    case Tool::rectangle:
        anchor_ = at_ = snapped(x, y, false);
        begin(Drag::rectangle);
        return;
    case Tool::contact: {
        anchor_ = at_ = snapped(x, y, false);
        begin(Drag::contact);
        return;
    }
    case Tool::measure:
        anchor_ = at_ = snapped(x, y, false);
        begin(Drag::measure);
        return;
    case Tool::polygon: {
        if (double_click) {
            on_key(VK_RETURN, false, false);
            return;
        }
        Pt p = snapped(x, y, false);
        if (!polygon_.empty() && shift) p = constrain(polygon_.back(), p);  // Shift: horizontal or vertical (slice 4)
        if (polygon_.size() >= 3 && p == polygon_.front()) {
            on_key(VK_RETURN, false, false);
            return;
        }
        if (polygon_.empty() || p != polygon_.back()) polygon_.push_back(p);
        invalidate();
        return;
    }
    default: break;
    }
    // Select: the domain's handles when the background is selected (a layered deck's).
    const std::optional<ShapeRef> one =
        selection_.size() == 1 && selection_[0].shape >= 0 ? std::optional(selection_[0]) : std::nullopt;
    if (selection_.size() == 1 && selection_[0].shape < 0 && g_.mesh == MeshKind::layered) {
        const Rect& d = g_.domain;
        const Coord mx = (d.lo.x + d.hi.x) / 2, my = (d.lo.y + d.hi.y) / 2;
        const Pt handles[8] = {d.lo, {mx, d.lo.y}, {d.hi.x, d.lo.y}, {d.hi.x, my}, d.hi, {mx, d.hi.y}, {d.lo.x, d.hi.y}, {d.lo.x, my}};
        for (int k = 0; k < 8; ++k) {
            if (std::abs(view_.sx(handles[k].x) - x) <= view_.px(4.5) && std::abs(view_.sy(handles[k].y) - y) <= view_.px(4.5)) {
                drag_index_ = k;
                anchor_ = handles[k];
                begin(Drag::domain);
                return;
            }
        }
    }
    const Hit hit = hit_test(g_, view_, x, y, one);
    // MT5: Mirror and Align act on a click.
    if ((mode_ == modeller::Mode::mirror || mode_ == modeller::Mode::align) && (hit.kind != HitKind::contact && hit.kind != HitKind::contact_end)) {
        transform_click(hit, view_.model(x, y), shift, ctrl);
        invalidate();
        return;
    }
    // Unit 26m slice 2 (MQ6): Chamfer and Round vertex drag in from a picked shape's corner (a polygon's vertex, a box's
    // corner handle); a press on a shape picks it, so its corners show.
    if (mode_ == modeller::Mode::chamfer || mode_ == modeller::Mode::round) {
        int corner = -1;
        if (hit.kind == HitKind::vertex && hit.shape.shape >= 0) corner = hit.index;
        else if (hit.kind == HitKind::handle && hit.shape.shape >= 0 && hit.index % 2 == 0) corner = hit.index / 2;
        if (corner >= 0) {
            corner_shape_ = hit.shape;
            drag_index_ = corner;
            corner_r_ = 0;
            corner_chords_ = 0;
            anchor_ = view_.model(x, y);
            begin(Drag::corner);
            return;
        }
        if ((hit.kind == HitKind::shape || hit.kind == HitKind::edge) && hit.shape.shape >= 0) {
            selection_ = {hit.shape};
            contact_ = -1;
            update_inspector();
            invalidate();
            host_.geometry_status("Now drag in from one of its corners.");
            return;
        }
    }
    // MT2: what is picked opens its group.
    {
        modeller::Picked k = modeller::Picked::nothing;
        if (hit.kind == HitKind::vertex || hit.kind == HitKind::edge) k = hit.shape.shape >= 0 ? modeller::Picked::edge : k;
        else if (hit.kind == HitKind::contact || hit.kind == HitKind::contact_end) k = modeller::Picked::contact;
        else if (hit.kind == HitKind::handle || (hit.kind == HitKind::shape && hit.shape.shape >= 0)) k = modeller::Picked::solid;
        const modeller::Group g = modeller::group_for(k);
        if (g != bar_->group()) {
            bar_->set_group(g);
            bar_->set_tool(modeller::armed_key(static_cast<int>(tool_), mode_));
        }
    }
    // MT1, MQ3: empty background (or a region's background) pans a free camera; Shift keeps the rubber band.
    const bool background = hit.kind == HitKind::none || ((hit.kind == HitKind::shape || hit.kind == HitKind::edge) && hit.shape.shape < 0);
    if (modeller::navigation(camera_fixed_, false, {modeller::Button::left, background, shift, true, false}) == modeller::Nav::pan) {
        contact_ = -1;
        if (!ctrl) selection_ = hit.kind == HitKind::none ? std::vector<ShapeRef>{} : std::vector<ShapeRef>{hit.shape};
        update_inspector();
        begin(Drag::pan);
        invalidate();
        return;
    }
    switch (hit.kind) {
    case HitKind::handle:
        drag_index_ = hit.index;
        anchor_ = view_.model(x, y);
        begin(Drag::handle);
        return;
    case HitKind::vertex:
        if (double_click) {
            Geometry next = g_;
            if (remove_vertex(next, hit.shape, hit.index)) commit_model(std::move(next));
            return;
        }
        drag_index_ = hit.index;
        anchor_ = g_.regions[static_cast<std::size_t>(hit.shape.region)].shapes[static_cast<std::size_t>(hit.shape.shape)]
                      .points[static_cast<std::size_t>(hit.index)];
        alone_ = ctrl;  // Ctrl: the vertex moves alone (its edges may slant); else its neighbours follow
        begin(Drag::vertex);
        return;
    case HitKind::edge:
        if (double_click) {
            Geometry next = g_;
            if (insert_vertex(next, hit.shape, hit.index, snapped(x, y, false))) commit_model(std::move(next));
            return;
        }
        break;
    case HitKind::contact_end:
        selection_.clear();
        contact_ = hit.index;
        drag_index_ = hit.end;
        anchor_ = hit.end == 0 ? g_.contacts[static_cast<std::size_t>(hit.index)].a : g_.contacts[static_cast<std::size_t>(hit.index)].b;
        update_inspector();
        begin(Drag::contact_end);
        return;
    case HitKind::contact:
        selection_.clear();
        contact_ = hit.index;
        update_inspector();
        invalidate();
        return;
    default: break;
    }
    contact_ = -1;
    if (hit.kind == HitKind::shape || hit.kind == HitKind::edge) {
        const ShapeRef s = hit.shape;
        const bool in = std::ranges::find(selection_, s) != selection_.end();
        if (ctrl && s.shape >= 0) {
            if (in) std::erase(selection_, s);
            else selection_.push_back(s);
        } else if (!in) {
            selection_ = {s};
        }
        update_inspector();
        invalidate();
        if (s.shape >= 0 && std::ranges::find(selection_, s) != selection_.end()) {
            if (mode_ == modeller::Mode::array) {
                // MT5: the drag is the pitch, the wheel the count; one region per copy on release.
                selection_ = {s};
                update_inspector();
                anchor_ = snapped(x, y, false, s);
                begin(Drag::array);
                return;
            }
            anchor_ = snapped(x, y, false, selection_.size() == 1 ? std::optional(s) : std::nullopt);
            begin(Drag::move);
            if (mode_ == modeller::Mode::copy || GetKeyState(VK_MENU) < 0) {
                // MT5: a copy of each picked shape is dragged (its region's rule kept); the originals stay.
                copy_from_ = g_;
                copy_sel_ = selection_;
                Geometry with = g_;
                std::vector<ShapeRef> copies;
                for (const ShapeRef c : selection_) {
                    if (c.shape < 0) continue;
                    const int r = duplicate_shape(with, c, 0, 0);
                    if (r < 0) continue;
                    with.regions[static_cast<std::size_t>(r)].fill = with.regions[static_cast<std::size_t>(c.region)].fill;
                    copies.push_back({r, 0});
                }
                g_ = with;
                g_start_ = with;
                selection_ = copies;
            }
        }
        return;
    }
    if (!ctrl) selection_.clear();
    update_inspector();
    anchor_ = at_ = view_.model(x, y);
    begin(Drag::rubber);
}

void GeometryPage::on_move(double x, double y, bool shift) {
    if (dimension_move(x, y)) return;  // Unit 26n (DT2)
    const Pt raw = view_.model(x, y);
    host_.geometry_status(std::format("x {} µm, y {} µm", ui_length(raw.x), ui_length(raw.y)));
    if (start_ == Start::footprint) {
        start_b_ = start_point(x, y);
        last_x_ = x;
        last_y_ = y;
        invalidate();
        return;
    }
    if (start_ == Start::height) {
        // The height follows the cursor up from where the footprint ended, on the grid, at least one step.
        const Coord step = grid_step(view_);
        const double um = (start_y_ - y) / view_.px_per_um;
        start_h_ = std::max(step, static_cast<Coord>(std::llround(um * per_um / static_cast<double>(step))) * step);
        last_x_ = x;
        last_y_ = y;
        invalidate();
        return;
    }
    const double dx = x - last_x_, dy = y - last_y_;
    last_x_ = x;
    last_y_ = y;
    if (tool_ == Tool::polygon && drag_ == Drag::none && g_.editable()) {
        Pt p = snapped(x, y, false);
        if (!polygon_.empty() && shift) p = constrain(polygon_.back(), p);  // Shift: horizontal or vertical (slice 4)
        hover_ = p;
        invalidate();
        return;
    }
    const std::optional<ShapeRef> one = selection_.size() == 1 ? std::optional(selection_[0]) : std::nullopt;
    switch (drag_) {
    case Drag::none: return;
    case Drag::pan:
        view_.origin_x += dx;
        view_.origin_y += dy;
        break;
    case Drag::array:
    case Drag::move: {
        const Pt p = snapped(x, y, shift, one);
        g_ = g_start_;
        drag_d_ = {p.x - anchor_.x, p.y - anchor_.y};
        if (typed_axis_ == 0) drag_d_.x = typed_or(drag_d_.x);
        if (typed_axis_ == 1) drag_d_.y = typed_or(drag_d_.y);
        if (drag_ == Drag::array) {
            if (one && (drag_d_.x != 0 || drag_d_.y != 0)) array_shape(g_, *one, drag_d_.x, drag_d_.y, array_count_);
            break;
        }
        for (const ShapeRef s : selection_) move_shape(g_, s, drag_d_.x, drag_d_.y);
        break;
    }
    case Drag::handle: {
        if (!one) break;
        const Pt p = snapped(x, y, shift, one);
        g_ = g_start_;
        const Shape& sh = g_.regions[static_cast<std::size_t>(one->region)].shapes[static_cast<std::size_t>(one->shape)];
        set_box(g_, *one, drag_handle(normalized(sh.points[0], sh.points[1]), drag_index_, p));
        break;
    }
    case Drag::vertex: {
        if (!one) break;
        const Pt p = snapped(x, y, shift, one);
        g_ = g_start_;
        drag_d_ = {p.x - anchor_.x, p.y - anchor_.y};
        set_vertex(g_, *one, drag_index_, p, !alone_);
        break;
    }
    case Drag::domain: {
        const Pt p = snapped(x, y, shift);
        g_ = g_start_;
        set_domain(g_, drag_handle(g_start_.domain, drag_index_, p));
        break;
    }
    case Drag::contact_end: {
        g_ = g_start_;
        Contact c = g_.contacts[static_cast<std::size_t>(contact_)];
        Pt p = snapped(x, y, false);
        const Rect& d = g_.domain;
        if (c.a.x == c.b.x) p = {c.a.x, std::clamp(p.y, d.lo.y, d.hi.y)};
        else p = {std::clamp(p.x, d.lo.x, d.hi.x), c.a.y};
        drag_d_ = {p.x - anchor_.x, p.y - anchor_.y};
        (drag_index_ == 0 ? c.a : c.b) = p;
        set_contact(g_, contact_, c.a, c.b);
        break;
    }
    case Drag::rubber: at_ = raw; break;
    case Drag::contact: {
        // Along the side the drag began on.
        const Rect& d = g_.domain;
        Pt p = snapped(x, y, false);
        // At a corner of the domain, the side the drag goes along (the larger offset).
        const bool on_x_side = anchor_.x == d.lo.x || anchor_.x == d.hi.x, on_y_side = anchor_.y == d.lo.y || anchor_.y == d.hi.y;
        const bool vertical = on_x_side && (!on_y_side || std::abs(p.y - anchor_.y) > std::abs(p.x - anchor_.x));
        if (vertical) p = {anchor_.x, std::clamp(p.y, d.lo.y, d.hi.y)};
        else p = {std::clamp(p.x, d.lo.x, d.hi.x), anchor_.y};
        at_ = p;
        break;
    }
    case Drag::rectangle:
        at_ = snapped(x, y, shift);
        if (typed_axis_ == 0) at_.x = anchor_.x + typed_or(at_.x - anchor_.x);
        if (typed_axis_ == 1) at_.y = anchor_.y + typed_or(at_.y - anchor_.y);
        break;
    case Drag::corner: {
        const Coord step = grid_step(view_);
        const double len = std::hypot(static_cast<double>(raw.x - anchor_.x), static_cast<double>(raw.y - anchor_.y));
        corner_r_ = std::max<Coord>(0, typed_or(static_cast<Coord>(std::llround(len / static_cast<double>(step))) * step));
        g_ = g_start_;
        corner_chords_ = 0;
        if (corner_r_ > 0) {
            if (mode_ == modeller::Mode::chamfer) {
                chamfer_shape(g_, corner_shape_, drag_index_, corner_r_);
            } else {
                // MQ4, MT12 (2): 15° a chord, fewer while one is under the generator's 6 s.
                const auto& sh = g_.regions[static_cast<std::size_t>(corner_shape_.region)].shapes[static_cast<std::size_t>(corner_shape_.shape)];
                const auto s = solid::mesh_length(text_seen_, "interface_resolution");
                const double six_s = s ? 6.0 * *s * static_cast<double>(per_um) * 1e4 : 0.0;
                corner_chords_ = arc_segments(corner_turn(sh.outline(), drag_index_), static_cast<double>(corner_r_), six_s);
                if (corner_chords_ > 0) fillet_shape(g_, corner_shape_, drag_index_, corner_r_, corner_chords_);
            }
        }
        break;
    }
    case Drag::measure: at_ = snapped(x, y, shift); break;
    }
    invalidate();
}

void GeometryPage::on_up(double x, double y, bool shift) {
    if (start_ == Start::footprint) {
        on_move(x, y, shift);
        if (start_a_.x == start_b_.x || start_a_.y == start_b_.y) {
            cancel();  // a click, or a line: nothing drawn
            return;
        }
        if (!start3d_) {
            start_finish();
            return;
        }
        start_ = Start::height;  // the height follows the cursor now; a click ends it (capture kept)
        start_y_ = y;
        start_h_ = grid_step(view_);
        invalidate();
        return;
    }
    if (start_ == Start::height) return;  // the click's own release
    on_move(x, y, shift);
    const Drag d = drag_;
    end_drag_ui();
    drag_ = Drag::none;
    if (GetCapture() == canvas_) ReleaseCapture();
    switch (d) {
    case Drag::corner:
    case Drag::move:
    case Drag::array:
    case Drag::handle:
    case Drag::vertex:
    case Drag::domain:
    case Drag::contact_end: {
        Geometry next = g_;
        const bool moved = !(next == g_start_);
        g_ = copy_from_ ? *copy_from_ : g_start_;
        if (copy_from_ && !moved) selection_ = copy_sel_;  // a copy not moved is no copy
        copy_from_.reset();
        if (moved) commit_model(std::move(next));
        break;
    }
    case Drag::rubber: {
        const Rect r = normalized(anchor_, at_);
        if (std::abs(x - press_x_) + std::abs(y - press_y_) < 3) break;
        for (std::size_t k = 0; k < g_.regions.size(); ++k) {
            for (std::size_t s = 0; s < g_.regions[k].shapes.size(); ++s) {
                const Rect b = g_.regions[k].shapes[s].bounds();
                const ShapeRef ref{static_cast<int>(k), static_cast<int>(s)};
                if (r.contains(b.lo) && r.contains(b.hi) && std::ranges::find(selection_, ref) == selection_.end()) {
                    selection_.push_back(ref);
                }
            }
        }
        break;
    }
    case Drag::rectangle: {
        Geometry next = g_;
        const int r = add_box(next, normalized(anchor_, at_), material_);
        if (r >= 0 && commit_model(std::move(next))) {
            selection_ = {{r, 0}};
            contact_ = -1;
        }
        break;
    }
    case Drag::contact: {
        Geometry next = g_;
        std::string name;
        for (int n = 1;; ++n) {
            name = "contact_" + std::to_string(n);
            // Not a contact's name, nor a patch's left in the deck (a renamed contact keeps its patch).
            if (std::ranges::none_of(g_.contacts, [&](const Contact& c) { return c.name == name; }) &&
                text_seen_.find("patch " + name + " ") == std::string::npos) {
                break;
            }
        }
        const int c = add_contact(next, name, "ohmic", anchor_, at_);
        if (c >= 0 && commit_model(std::move(next))) {
            contact_ = c;
            selection_.clear();
        } else if (c < 0 && !(anchor_ == at_)) {
            host_.geometry_status("A contact is drawn along one side of the domain: start the drag on its boundary.");
        }
        break;
    }
    default: break;
    }
    update_inspector();
    invalidate();
}

bool GeometryPage::on_key(UINT key, bool shift, bool ctrl) {
    if (ctrl) return false;
    // MT3: while a drag takes a typed value, letters and digits are its (WM_CHAR), not the tools' keys.
    if ((drag_ == Drag::move || drag_ == Drag::rectangle || drag_ == Drag::array || drag_ == Drag::corner) && ((key >= 'A' && key <= 'Z') || (key >= '0' && key <= '9') ||
                                                          (key >= VK_NUMPAD0 && key <= VK_DIVIDE) || key == VK_OEM_PERIOD || key == VK_OEM_MINUS || key == VK_BACK)) {
        return true;
    }
    switch (key) {
    case 'L': set_camera_fixed(!camera_fixed_); return true;
    case 'V': set_tool(Tool::select); return true;
    case 'R': set_tool(Tool::rectangle); return true;
    case 'P': set_tool(Tool::polygon); return true;
    case 'C': set_tool(Tool::contact); return true;
    case 'H': set_tool(Tool::pan); return true;
    case 'M': set_tool(Tool::measure); return true;
    case 'F':
    case VK_HOME: fit(); return true;
    case VK_DELETE: remove(); return true;
    case VK_ESCAPE:
        if (start_ == Start::none && drag_ == Drag::none && polygon_.empty()) {
            selection_.clear();
            contact_ = -1;
            update_inspector();
        }
        cancel();
        return true;
    case VK_RETURN: {
        if (drag_ != Drag::none && drag_ != Drag::pan && drag_ != Drag::rubber) {  // MT3: Enter applies the drag as it stands
            on_up(last_x_, last_y_, false);
            return true;
        }
        if (tool_ != Tool::polygon || polygon_.size() < 3) return true;
        std::vector<Pt> p = polygon_;
        polygon_.clear();
        hover_.reset();
        Geometry next = g_;
        const int r = add_polygon(next, std::move(p), material_);
        if (r >= 0 && commit_model(std::move(next))) {
            selection_ = {{r, 0}};
            contact_ = -1;
            update_inspector();
        } else if (r < 0) {
            host_.geometry_status("A polygon needs three or more vertices, its edges not crossing (a tensor deck takes "
                                  "rectangles only).");
        }
        invalidate();
        return true;
    }
    case VK_LEFT:
    case VK_RIGHT:
    case VK_UP:
    case VK_DOWN: {
        if (selection_.empty()) return true;
        const Coord step = grid_step(view_) * (shift ? 10 : 1);
        const Coord dx = key == VK_LEFT ? -step : key == VK_RIGHT ? step : 0;
        const Coord dy = key == VK_DOWN ? -step : key == VK_UP ? step : 0;
        Geometry next = g_;
        bool any = false;
        for (const ShapeRef s : selection_) any = move_shape(next, s, dx, dy) || any;
        if (any) commit_model(std::move(next));
        update_inspector();
        return true;
    }
    default: return false;
    }
}

// --- drawing --------------------------------------------------------------------------------------------------

void GeometryPage::invalidate() {
    if (canvas_) InvalidateRect(canvas_, nullptr, FALSE);
}

bool GeometryPage::ensure_target() {
    if (!d2d_ || !dwrite_) return false;
    if (!text_) {
        dwrite_->CreateTextFormat(L"Segoe UI", nullptr, DWRITE_FONT_WEIGHT_NORMAL, DWRITE_FONT_STYLE_NORMAL,
                                  DWRITE_FONT_STRETCH_NORMAL, static_cast<float>(12.0 * dpi_ / 96.0), L"", &text_);
        if (!text_) return false;
    }
    if (target_) return true;
    RECT c{};
    GetClientRect(canvas_, &c);
    const auto props = D2D1::RenderTargetProperties(software_ ? D2D1_RENDER_TARGET_TYPE_SOFTWARE : D2D1_RENDER_TARGET_TYPE_DEFAULT,
                                                    D2D1::PixelFormat(), 96.0f, 96.0f);
    const auto hwnd_props = D2D1::HwndRenderTargetProperties(
        canvas_, D2D1::SizeU(static_cast<UINT32>(std::max<LONG>(1, c.right)), static_cast<UINT32>(std::max<LONG>(1, c.bottom))));
    return SUCCEEDED(d2d_->CreateHwndRenderTarget(props, hwnd_props, &target_)) && target_;
}

void GeometryPage::release_target() { release(target_); }

std::string GeometryPage::ruler_text(bool x_axis) const { return x_axis ? ruler_x_ : ruler_y_; }

double GeometryPage::render() {
    if (!ensure_target()) return -1.0;
    const auto t0 = std::chrono::steady_clock::now();
    RECT c{};
    GetClientRect(canvas_, &c);
    const float W = static_cast<float>(c.right), H = static_cast<float>(c.bottom);
    const float R = static_cast<float>(ruler_px());
    ID2D1SolidColorBrush* brush = nullptr;
    target_->BeginDraw();
    target_->SetTransform(D2D1::Matrix3x2F::Identity());
    target_->Clear(rgb(0.965f, 0.965f, 0.965f));
    target_->CreateSolidColorBrush(rgb(0, 0, 0), &brush);
    if (!brush) {
        target_->EndDraw();
        return -1.0;
    }
    const auto color = [&](D2D1_COLOR_F k) { brush->SetColor(k); return brush; };
    const auto text = [&](std::wstring_view t, float x, float y, float w, float h, D2D1_COLOR_F k) {
        target_->DrawText(t.data(), static_cast<UINT32>(t.size()), text_, D2D1::RectF(x, y, x + w, y + h), color(k));
    };
    const auto sx = [&](Coord x) { return static_cast<float>(view_.sx(x)); };
    const auto sy = [&](Coord y) { return static_cast<float>(view_.sy(y)); };
    const float px = static_cast<float>(dpi_) / 96.0f;
    if (can_start()) {
        // A deck with no mesh: a grid to draw a new device on, the hint, and the device being drawn with its size.
        target_->PushAxisAlignedClip(D2D1::RectF(R, R, W, H), D2D1_ANTIALIAS_MODE_ALIASED);
        target_->FillRectangle(D2D1::RectF(R, R, W, H), color(rgb(1, 1, 1)));
        const Coord step = grid_step(view_);
        const Pt a = view_.model(R, H), b = view_.model(W, R);
        const auto first = [&](Coord v) { return static_cast<Coord>(std::floor(static_cast<double>(v) / static_cast<double>(step))) * step; };
        if ((b.x - a.x) / step < 600 && (b.y - a.y) / step < 600) {
            const D2D1_COLOR_F grid = rgb(0, 0, 0, 0.07f);
            for (Coord x = first(a.x); x <= b.x; x += step) target_->DrawLine({sx(x), R}, {sx(x), H}, color(grid), 1.0f);
            for (Coord y = first(a.y); y <= b.y; y += step) target_->DrawLine({R, sy(y)}, {W, sy(y)}, color(grid), 1.0f);
        }
        target_->DrawLine({sx(0), R}, {sx(0), H}, color(rgb(0, 0, 0, 0.3f)), 1.0f);  // the axes through the origin
        target_->DrawLine({R, sy(0)}, {W, sy(0)}, color(rgb(0, 0, 0, 0.3f)), 1.0f);
        const auto fill = [&](std::initializer_list<D2D1_POINT_2F> pts, D2D1_COLOR_F k) {
            ID2D1PathGeometry* path = nullptr;
            d2d_->CreatePathGeometry(&path);
            if (!path) return;
            ID2D1GeometrySink* sink = nullptr;
            path->Open(&sink);
            if (sink) {
                auto p = pts.begin();
                sink->BeginFigure(*p, D2D1_FIGURE_BEGIN_FILLED);
                for (++p; p != pts.end(); ++p) sink->AddLine(*p);
                sink->EndFigure(D2D1_FIGURE_END_CLOSED);
                sink->Close();
                sink->Release();
                target_->FillGeometry(path, color(k));
                target_->DrawGeometry(path, color(rgb(0.13f, 0.33f, 0.66f)), 1.5f * px);
            }
            path->Release();
        };
        if (start_ != Start::none) {
            const float l = std::min(sx(start_a_.x), sx(start_b_.x)), r = std::max(sx(start_a_.x), sx(start_b_.x));
            const float t = std::min(sy(start_a_.y), sy(start_b_.y)), bt = std::max(sy(start_a_.y), sy(start_b_.y));
            const D2D1_COLOR_F si = material_fill("silicon");
            if (start_ == Start::height) {
                // The box seen obliquely: its top face drawn up and to the right of the footprint by the height.
                const float hp = static_cast<float>(static_cast<double>(start_h_) / per_um * view_.px_per_um) * 0.7071f;
                fill({{l, bt}, {r, bt}, {r, t}, {l, t}}, si);
                fill({{l, t}, {r, t}, {r + hp, t - hp}, {l + hp, t - hp}}, rgb(si.r * 1.08f, si.g * 1.08f, si.b * 1.08f));
                fill({{r, t}, {r, bt}, {r + hp, bt - hp}, {r + hp, t - hp}}, rgb(si.r * 0.8f, si.g * 0.8f, si.b * 0.8f));
            } else {
                fill({{l, bt}, {r, bt}, {r, t}, {l, t}}, rgb(si.r, si.g, si.b, 0.85f));
            }
            // The live size at the cursor, on a white label.
            const std::wstring s = wide(start_readout() + (start_ == Start::height ? "   click to finish, Esc to cancel" : ""));
            const float lx = static_cast<float>(last_x_) + 14 * px, ly = static_cast<float>(last_y_) + 14 * px;
            const float lw = (7.2f * static_cast<float>(s.size()) + 12) * px, lh = 22 * px;
            target_->FillRoundedRectangle(D2D1::RoundedRect(D2D1::RectF(lx, ly, lx + lw, ly + lh), 4 * px, 4 * px), color(rgb(1, 1, 1, 0.95f)));
            target_->DrawRoundedRectangle(D2D1::RoundedRect(D2D1::RectF(lx, ly, lx + lw, ly + lh), 4 * px, 4 * px), color(rgb(0, 0, 0, 0.25f)), 1.0f);
            text(s, lx + 6 * px, ly + 3 * px, lw, lh, rgb(0.1f, 0.1f, 0.1f));
        }
        target_->PopAxisAlignedClip();
        text_->SetTextAlignment(DWRITE_TEXT_ALIGNMENT_CENTER);
        text(L"Drag to draw a 2D device.   Ctrl+drag for 3D: release, move up for its height, click.", R, R + 10 * px, W - R, 24 * px,
             rgb(0.35f, 0.35f, 0.35f));
        text_->SetTextAlignment(DWRITE_TEXT_ALIGNMENT_LEADING);
    } else if (!g_.editable()) {
        const std::wstring m = wide(message_.empty() ? std::string("Open a 2D deck to draw its device.") : message_);
        text_->SetTextAlignment(DWRITE_TEXT_ALIGNMENT_CENTER);
        text(m, 0, H / 2 - 20 * px, W, 40 * px, rgb(0.35f, 0.35f, 0.35f));
        text_->SetTextAlignment(DWRITE_TEXT_ALIGNMENT_LEADING);
    } else {
        target_->PushAxisAlignedClip(D2D1::RectF(R, R, W, H), D2D1_ANTIALIAS_MODE_ALIASED);
        // The domain and the grid inside it.
        const Rect& d = g_.domain;
        const D2D1_RECT_F dom = D2D1::RectF(sx(d.lo.x), sy(d.hi.y), sx(d.hi.x), sy(d.lo.y));
        const int bg = g_.background();
        target_->FillRectangle(dom, color(bg >= 0 ? material_fill(g_.regions[static_cast<std::size_t>(bg)].material) : rgb(1, 1, 1)));
        const Coord step = grid_step(view_);
        const Pt a = view_.model(R, H), b = view_.model(W, R);
        const auto first = [&](Coord v) { return static_cast<Coord>(std::floor(static_cast<double>(v) / static_cast<double>(step))) * step; };
        if ((b.x - a.x) / step < 600 && (b.y - a.y) / step < 600) {
            const D2D1_COLOR_F grid = rgb(0, 0, 0, 0.07f);
            for (Coord x = first(a.x); x <= b.x; x += step) target_->DrawLine({sx(x), R}, {sx(x), H}, color(grid), 1.0f);
            for (Coord y = first(a.y); y <= b.y; y += step) target_->DrawLine({R, sy(y)}, {W, sy(y)}, color(grid), 1.0f);
        }
        // The layers in deck order, each shape over the ones before.
        const auto path_of = [&](const std::vector<Pt>& pts) {
            ID2D1PathGeometry* path = nullptr;
            d2d_->CreatePathGeometry(&path);
            if (!path) return path;
            ID2D1GeometrySink* sink = nullptr;
            path->Open(&sink);
            if (sink) {
                sink->BeginFigure({sx(pts[0].x), sy(pts[0].y)}, D2D1_FIGURE_BEGIN_FILLED);
                for (std::size_t k = 1; k < pts.size(); ++k) sink->AddLine({sx(pts[k].x), sy(pts[k].y)});
                sink->EndFigure(D2D1_FIGURE_END_CLOSED);
                sink->Close();
                sink->Release();
            }
            return path;
        };
        for (const int ri : geo::paint_order(g_)) {  // MT4: as the build layers them
            const auto r = static_cast<std::size_t>(ri);
            const D2D1_COLOR_F fill = material_fill(g_.regions[r].material);
            const std::size_t n = g_.regions[r].shapes.size();
            for (std::size_t k = 0; k < n; ++k) {
                const Shape& s = g_.regions[r].shapes[g_.regions[r].fill ? n - 1 - k : k];
                if (s.kind == ShapeKind::box) {
                    const Rect q = normalized(s.points[0], s.points[1]);
                    const D2D1_RECT_F rr = D2D1::RectF(sx(q.lo.x), sy(q.hi.y), sx(q.hi.x), sy(q.lo.y));
                    target_->FillRectangle(rr, color(fill));
                    target_->DrawRectangle(rr, color(darker(fill)), 1.0f);
                } else if (ID2D1PathGeometry* path = path_of(s.points)) {
                    target_->FillGeometry(path, color(fill));
                    target_->DrawGeometry(path, color(darker(fill)), 1.0f);
                    path->Release();
                }
            }
        }
        // The mesh preview: the last check's mesh, while it is of the text shown.
        if (mesh_shown()) {
            const D2D1_COLOR_F wire = rgb(0.15f, 0.15f, 0.2f, 0.55f);
            for (const auto& [i, j] : mesh_.edges) {
                const Pt& p = mesh_.points[static_cast<std::size_t>(i)];
                const Pt& q = mesh_.points[static_cast<std::size_t>(j)];
                target_->DrawLine({sx(p.x), sy(p.y)}, {sx(q.x), sy(q.y)}, color(wire), 0.75f * px);
            }
        }
        target_->DrawRectangle(dom, color(rgb(0.2f, 0.2f, 0.2f)), 1.5f * px);
        // Contacts, with their names.
        for (std::size_t k = 0; k < g_.contacts.size(); ++k) {
            const Contact& ct = g_.contacts[k];
            const bool sel = static_cast<int>(k) == contact_;
            const D2D1_POINT_2F p0{sx(ct.a.x), sy(ct.a.y)}, p1{sx(ct.b.x), sy(ct.b.y)};
            target_->DrawLine(p0, p1, color(sel ? rgb(0.1f, 0.45f, 0.95f) : contact_color(ct.kind)), (sel ? 6.0f : 4.0f) * px);
            if (sel) {
                for (const auto& q : {p0, p1}) {
                    const D2D1_RECT_F h = D2D1::RectF(q.x - 3.5f * px, q.y - 3.5f * px, q.x + 3.5f * px, q.y + 3.5f * px);
                    target_->FillRectangle(h, color(rgb(1, 1, 1)));
                    target_->DrawRectangle(h, color(rgb(0.1f, 0.45f, 0.95f)), 1.0f);
                }
            }
            text(wide(ct.name), (p0.x + p1.x) / 2 + 4 * px, (p0.y + p1.y) / 2 + 2 * px, 160 * px, 18 * px, rgb(0.15f, 0.15f, 0.15f));
        }
        // The build's problems where they lie (GD9): a shape or region outlined, a contact, a point ringed, the
        // domain.
        const D2D1_COLOR_F red = rgb(0.85f, 0.1f, 0.1f);
        for (const Problem& pr : problems_) {
            for (const Marker& mk : pr.markers) {
                switch (mk.kind) {
                case Marker::Kind::shape: {
                    const Region& rg = g_.regions[static_cast<std::size_t>(mk.shape.region)];
                    for (std::size_t s = 0; s < rg.shapes.size(); ++s) {
                        if (mk.shape.shape >= 0 && static_cast<int>(s) != mk.shape.shape) continue;
                        if (ID2D1PathGeometry* path = path_of(rg.shapes[s].outline())) {
                            target_->DrawGeometry(path, color(red), 2.5f * px);
                            path->Release();
                        }
                    }
                    break;
                }
                case Marker::Kind::contact: {
                    const Contact& ct = g_.contacts[static_cast<std::size_t>(mk.contact)];
                    target_->DrawLine({sx(ct.a.x), sy(ct.a.y)}, {sx(ct.b.x), sy(ct.b.y)}, color(rgb(0.85f, 0.1f, 0.1f, 0.6f)), 9.0f * px);
                    break;
                }
                case Marker::Kind::point: {
                    const D2D1_POINT_2F o{sx(mk.at.x), sy(mk.at.y)};
                    target_->DrawEllipse(D2D1::Ellipse(o, 7.0f * px, 7.0f * px), color(red), 2.0f * px);
                    target_->DrawLine({o.x - 3 * px, o.y - 3 * px}, {o.x + 3 * px, o.y + 3 * px}, color(red), 1.5f * px);
                    target_->DrawLine({o.x - 3 * px, o.y + 3 * px}, {o.x + 3 * px, o.y - 3 * px}, color(red), 1.5f * px);
                    break;
                }
                case Marker::Kind::domain: target_->DrawRectangle(dom, color(red), 3.0f * px); break;
                }
            }
        }
        // The selection: outlines; the handles of one box (or the domain), the vertices of one polygon.
        const D2D1_COLOR_F blue = rgb(0.1f, 0.45f, 0.95f);
        const auto handle = [&](Pt p) {
            const D2D1_RECT_F h = D2D1::RectF(sx(p.x) - 3.5f * px, sy(p.y) - 3.5f * px, sx(p.x) + 3.5f * px, sy(p.y) + 3.5f * px);
            target_->FillRectangle(h, color(rgb(1, 1, 1)));
            target_->DrawRectangle(h, color(blue), 1.0f);
        };
        const auto handles_of = [&](Rect q) {
            const Coord mx = (q.lo.x + q.hi.x) / 2, my = (q.lo.y + q.hi.y) / 2;
            for (const Pt p : {q.lo, Pt{mx, q.lo.y}, Pt{q.hi.x, q.lo.y}, Pt{q.hi.x, my}, q.hi, Pt{mx, q.hi.y}, Pt{q.lo.x, q.hi.y}, Pt{q.lo.x, my}}) {
                handle(p);
            }
        };
        for (const ShapeRef s : selection_) {
            if (s.shape < 0) {
                target_->DrawRectangle(dom, color(blue), 2.0f * px);
                if (selection_.size() == 1 && g_.mesh == MeshKind::layered) handles_of(d);
                continue;
            }
            const Shape& sh = g_.regions[static_cast<std::size_t>(s.region)].shapes[static_cast<std::size_t>(s.shape)];
            if (ID2D1PathGeometry* path = path_of(sh.outline())) {
                target_->DrawGeometry(path, color(blue), 2.0f * px);
                path->Release();
            }
            if (selection_.size() == 1) {
                if (sh.kind == ShapeKind::box) handles_of(normalized(sh.points[0], sh.points[1]));
                else for (const Pt& p : sh.points) handle(p);
            }
        }
        // What a tool is drawing.
        const D2D1_COLOR_F ghost = rgb(0.1f, 0.45f, 0.95f, 0.18f);
        if (drag_ == Drag::rubber || drag_ == Drag::rectangle) {
            const Rect q = normalized(anchor_, at_);
            const D2D1_RECT_F rr = D2D1::RectF(sx(q.lo.x), sy(q.hi.y), sx(q.hi.x), sy(q.lo.y));
            target_->FillRectangle(rr, color(drag_ == Drag::rectangle ? material_fill(material_) : ghost));
            target_->DrawRectangle(rr, color(blue), 1.0f);
        }
        if (drag_ == Drag::contact) target_->DrawLine({sx(anchor_.x), sy(anchor_.y)}, {sx(at_.x), sy(at_.y)}, color(contact_color("ohmic")), 4.0f * px);
        if (drag_ == Drag::measure || (tool_ == Tool::measure && !(anchor_ == at_))) {
            target_->DrawLine({sx(anchor_.x), sy(anchor_.y)}, {sx(at_.x), sy(at_.y)}, color(rgb(0.8f, 0.1f, 0.1f)), 1.5f * px);
            const double dx = static_cast<double>(at_.x - anchor_.x) / per_um, dy = static_cast<double>(at_.y - anchor_.y) / per_um;
            text(wide(std::format("dx {:.4g} µm, dy {:.4g} µm, length {:.4g} µm", dx, dy, std::hypot(dx, dy))), sx(at_.x) + 6 * px,
                 sy(at_.y) + 4 * px, 320 * px, 18 * px, rgb(0.6f, 0.05f, 0.05f));
        }
        if (!polygon_.empty()) {
            std::vector<Pt> pts = polygon_;
            if (hover_) pts.push_back(*hover_);
            for (std::size_t k = 0; k + 1 < pts.size(); ++k) {
                target_->DrawLine({sx(pts[k].x), sy(pts[k].y)}, {sx(pts[k + 1].x), sy(pts[k + 1].y)}, color(blue), 1.5f * px);
            }
            for (const Pt& p : polygon_) handle(p);
        }
        // The foot line: the first problem's message (Run is off), or the preview's state.
        if (const std::string foot = banner(); !foot.empty()) {
            const float bh = 40.0f * px;
            const bool bad = check() == Check::invalid;
            target_->FillRectangle(D2D1::RectF(R, H - bh, W, H), color(bad ? rgb(1.0f, 0.92f, 0.92f, 0.95f) : rgb(1, 1, 1, 0.9f)));
            target_->DrawLine({R, H - bh}, {W, H - bh}, color(bad ? red : rgb(0.6f, 0.6f, 0.6f)), 1.0f);
            text(wide(foot), R + 8 * px, H - bh + 4 * px, std::max(10.0f, W - R - 16 * px), bh - 6 * px, bad ? rgb(0.6f, 0.0f, 0.0f) : rgb(0.3f, 0.3f, 0.3f));
        }
        target_->PopAxisAlignedClip();
    }
    // Unit 26n (DT2): the dimension lines, each 12 px off its extent, with end ticks and its label.
    if (!can_start() && !marks_.empty()) {
        target_->PushAxisAlignedClip(D2D1::RectF(R, R, W, H), D2D1_ANTIALIAS_MODE_ALIASED);
        const auto lines = dimension_lines();
        for (std::size_t k = 0; k < lines.size(); ++k) {
            const auto& s = lines[k];
            const bool lit = marks_[k].name == dims_->selected();
            const D2D1_COLOR_F ink = rgb(0.13f, 0.40f, 0.80f);
            target_->DrawLine({s[0], s[1]}, {s[2], s[3]}, color(ink), lit ? 2.0f * px : 1.0f * px);
            const bool across = marks_[k].axis == 0;  // ticks run across the line
            for (const auto& [x, y] : {std::pair{s[0], s[1]}, std::pair{s[2], s[3]}}) {
                if (across) target_->DrawLine({x, y - 4 * px}, {x, y + 4 * px}, color(ink), 1.0f * px);
                else target_->DrawLine({x - 4 * px, y}, {x + 4 * px, y}, color(ink), 1.0f * px);
            }
            const float mx = 0.5f * (s[0] + s[2]), my = 0.5f * (s[1] + s[3]);
            text_->SetTextAlignment(DWRITE_TEXT_ALIGNMENT_CENTER);
            text(wide(marks_[k].label), mx - 120 * px, my - 18 * px, 240 * px, 16 * px, ink);
            text_->SetTextAlignment(DWRITE_TEXT_ALIGNMENT_LEADING);
        }
        target_->PopAxisAlignedClip();
    }
    // The rulers, in µm, labelled at the 1, 2, 5 step at least 60 px apart.
    target_->FillRectangle(D2D1::RectF(0, 0, W, R), color(rgb(0.92f, 0.92f, 0.92f)));
    target_->FillRectangle(D2D1::RectF(0, 0, R, H), color(rgb(0.92f, 0.92f, 0.92f)));
    ruler_x_.clear();
    ruler_y_.clear();
    {
        View coarse = view_;
        coarse.px_per_um = view_.px_per_um / 5.0;
        const Coord step = grid_step(coarse);
        const Pt a = view_.model(R, H), b = view_.model(W, R);
        const auto first = [&](Coord v) { return static_cast<Coord>(std::ceil(static_cast<double>(v) / static_cast<double>(step))) * step; };
        const D2D1_COLOR_F ink = rgb(0.3f, 0.3f, 0.3f);
        int n = 0;
        for (Coord x = first(a.x); x <= b.x && n < 200; x += step, ++n) {
            target_->DrawLine({sx(x), R - 6 * px}, {sx(x), R}, color(ink), 1.0f);
            const std::string t = ui_length(x);
            ruler_x_ += t + " ";
            text(wide(t), sx(x) + 2 * px, 1 * px, 80 * px, R, ink);
        }
        n = 0;
        for (Coord y = first(a.y); y <= b.y && n < 200; y += step, ++n) {
            target_->DrawLine({R - 6 * px, sy(y)}, {R, sy(y)}, color(ink), 1.0f);
            const std::string t = ui_length(y);
            ruler_y_ += t + " ";
            text(wide(t), 1 * px, sy(y) + 1 * px, R + 40 * px, 16 * px, ink);
        }
        text(L"µm", 2 * px, 2 * px, R, R, ink);
    }
    brush->Release();
    const HRESULT hr = target_->EndDraw();
    if (hr == static_cast<HRESULT>(D2DERR_RECREATE_TARGET)) release_target();
    ValidateRect(canvas_, nullptr);
    place_bar();
    return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

// --- window procedures ----------------------------------------------------------------------------------------

LRESULT CALLBACK GeometryPage::panel_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l) {
    auto* self = reinterpret_cast<GeometryPage*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
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
                    case ModellerBar::reset:
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

LRESULT CALLBACK GeometryPage::inspector_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l) {
    auto* self = reinterpret_cast<GeometryPage*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
    if (self) {
        if (message == WM_SIZE) {
            self->layout_inspector();
            return 0;
        }
        if (message == WM_COMMAND) {
            const int c = LOWORD(w) - self->inspector_id_;
            const int code = HIWORD(w);
            if ((code == BN_CLICKED && (c == c_raise || c == c_lower || c == c_apply || c == c_overlap)) ||
                (code == CBN_SELCHANGE && (c == c_material || c == c_kind || c == c_draw))) {
                self->commit(c);
                return 0;
            }
        }
    }
    return DefWindowProcW(hwnd, message, w, l);
}

LRESULT CALLBACK GeometryPage::canvas_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l) {
    if (message == WM_NCCREATE) {
        SetWindowLongPtrW(hwnd, GWLP_USERDATA,
                          reinterpret_cast<LONG_PTR>(reinterpret_cast<CREATESTRUCTW*>(l)->lpCreateParams));
    }
    auto* self = reinterpret_cast<GeometryPage*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
    if (self && self->canvas_ == hwnd) return self->on_canvas(message, w, l);
    if (self && message == WM_CREATE) {
        self->canvas_ = hwnd;
        return 0;
    }
    return DefWindowProcW(hwnd, message, w, l);
}

LRESULT GeometryPage::on_canvas(UINT message, WPARAM w, LPARAM l) {
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
    case WM_SIZE:
        if (target_) target_->Resize(D2D1::SizeU(std::max(1u, static_cast<UINT>(LOWORD(l))), std::max(1u, static_cast<UINT>(HIWORD(l)))));
        if (!fitted_) fit();
        return 0;
    case WM_LBUTTONDOWN: on_down(x, y, shift, ctrl, false); return 0;
    case WM_LBUTTONDBLCLK: on_down(x, y, shift, ctrl, true); return 0;
    case WM_MOUSEMOVE:
        on_move(x, y, shift);
        after_drag_step(x, y);
        return 0;
    case WM_LBUTTONUP:
        if (dimension_up()) return 0;  // Unit 26n (DT2)
        if (start_ != Start::none || (drag_ != Drag::none && drag_ != Drag::pan)) on_up(x, y, shift);
        else if (drag_ == Drag::pan) {
            drag_ = Drag::none;
            ReleaseCapture();
        }
        return 0;
    case WM_MBUTTONDOWN:
        SetFocus(canvas_);
        if (modeller::navigation(camera_fixed_, false, {modeller::Button::middle}) == modeller::Nav::none) return 0;
        last_x_ = x;
        last_y_ = y;
        drag_ = Drag::pan;
        SetCapture(canvas_);
        return 0;
    case WM_MBUTTONUP:
        if (drag_ == Drag::pan) {
            drag_ = Drag::none;
            ReleaseCapture();
        }
        return 0;
    case WM_MOUSEWHEEL: {
        if (drag_ == Drag::array) {  // MT5: the wheel sets the count, whatever the camera
            array_count_ = std::clamp(array_count_ + (GET_WHEEL_DELTA_WPARAM(w) > 0 ? 1 : -1), 2, 100);
            on_move(last_x_, last_y_, false);
            after_drag_step(last_x_, last_y_);
            return 0;
        }
        if (!modeller::wheel_zooms(camera_fixed_)) return 0;
        POINT p{GET_X_LPARAM(l), GET_Y_LPARAM(l)};
        ScreenToClient(canvas_, &p);
        const double f = std::pow(1.2, GET_WHEEL_DELTA_WPARAM(w) / 120.0);
        const double s = std::clamp(view_.px_per_um * f, 1e-2, 1e6);
        const double mx = (p.x - view_.origin_x) / view_.px_per_um, my = (view_.origin_y - p.y) / view_.px_per_um;
        view_.px_per_um = s;
        view_.origin_x = p.x - mx * s;
        view_.origin_y = p.y + my * s;
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
        if (reinterpret_cast<HWND>(l) != canvas_ && drag_ != Drag::none) cancel();
        return 0;
    default: break;
    }
    return DefWindowProcW(canvas_, message, w, l);
}

// --- Unit 26n (DT1, DT2): the dimensions panel and the dimension lines -------------------------------------------------

// A deck the editor reads as written; else, when its shapes hold $references, read from the text with every reference
// written out and shown view only (the dimensions edit it; shape tools refuse).
std::expected<geo::Geometry, std::string> GeometryPage::read_page(const std::string& text) {
    view_only_ = false;
    auto g = read_geometry(text, host_.geometry_bounds());
    if ((!g || !g->editable()) && text.find('$') != std::string::npos) {
        if (const auto sub = dims::substituted_text(text)) {
            if (auto gs = read_geometry(*sub, host_.geometry_bounds()); gs && gs->editable()) {
                view_only_ = true;
                return gs;
            }
        }
    }
    return g;
}

std::string GeometryPage::DimsHost::dimensions_deck() { return page->host_.geometry_deck(); }

void GeometryPage::DimsHost::dimensions_edited(const std::string& text) {
    page->history_.record(page->host_.geometry_deck());  // one undo for the edit
    page->apply_text(text);
}

void GeometryPage::DimsHost::dimensions_selected(const std::string&) { page->invalidate(); }

void GeometryPage::DimsHost::dimensions_sweep(const std::string& name) { page->host_.geometry_sweep(name); }

void GeometryPage::reload_dimensions() {
    dims_->reload();
    marks_ = dims::dimension_marks(text_seen_);
    layout_inspector();
}

std::vector<std::array<float, 4>> GeometryPage::dimension_lines() const {
    std::vector<std::array<float, 4>> out;
    const auto q = [](double cm) { return static_cast<geo::Coord>(std::llround(cm / geo::quantum_cm)); };
    const auto lift = static_cast<float>(view_.px(12.0));
    for (const dims::Mark& m : marks_) {
        const float x0 = static_cast<float>(view_.sx(q(m.a[0]))), y0 = static_cast<float>(view_.sy(q(m.a[1])));
        const float x1 = static_cast<float>(view_.sx(q(m.b[0]))), y1 = static_cast<float>(view_.sy(q(m.b[1])));
        if (m.axis == 0) out.push_back({x0, y0 - lift, x1, y1 - lift});
        else out.push_back({x0 - lift, y0, x1 - lift, y1});
    }
    return out;
}

namespace {

double segment_distance(const std::array<float, 4>& s, double x, double y) {
    const double dx = s[2] - s[0], dy = s[3] - s[1], len2 = dx * dx + dy * dy;
    const double t = len2 > 0.0 ? std::clamp(((x - s[0]) * dx + (y - s[1]) * dy) / len2, 0.0, 1.0) : 0.0;
    return std::hypot(x - (s[0] + t * dx), y - (s[1] + t * dy));
}

}  // namespace

bool GeometryPage::dimension_down(double x, double y) {
    const auto lines = dimension_lines();
    for (std::size_t k = 0; k < lines.size(); ++k) {
        if (segment_distance(lines[k], x, y) > view_.px(6.0)) continue;
        dim_drag_ = static_cast<int>(k);
        dim_x_ = x;
        dim_y_ = y;
        dim_literal_.clear();
        dims_->select(marks_[k].name);
        SetCapture(canvas_);
        invalidate();
        return true;
    }
    return false;
}

bool GeometryPage::dimension_move(double x, double y) {
    if (dim_drag_ < 0) return false;
    const dims::Mark& m = marks_[static_cast<std::size_t>(dim_drag_)];
    const auto& rows = dims_->dimensions();
    const auto row = std::ranges::find(rows, m.name, &dims::Dimension::name);
    if (row == rows.end()) return true;
    const geo::Pt p = view_.model(x, y), p0 = view_.model(dim_x_, dim_y_);
    const geo::Coord d = m.axis == 0 ? p.x - p0.x : p.y - p0.y;
    dim_literal_ = dims::length_literal(*row, row->value + m.grows * static_cast<double>(d) * geo::quantum_cm);
    host_.geometry_status(row->name + " = " + dim_literal_);
    return true;
}

bool GeometryPage::dimension_up() {
    if (dim_drag_ < 0) return false;
    const std::string name = marks_[static_cast<std::size_t>(dim_drag_)].name;
    dim_drag_ = -1;
    if (GetCapture() == canvas_) ReleaseCapture();
    const auto& rows = dims_->dimensions();
    const auto row = std::ranges::find(rows, name, &dims::Dimension::name);
    if (row != rows.end() && !dim_literal_.empty() && dim_literal_ != row->expr) {
        dims_->commit(static_cast<std::size_t>(row - rows.begin()), dim_literal_);  // one edit, on release
    }
    dim_literal_.clear();
    invalidate();
    return true;
}

}  // namespace NiTCAD::app
