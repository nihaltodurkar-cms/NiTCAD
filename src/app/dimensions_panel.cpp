#include "dimensions_panel.hpp"

#include "win32.hpp"

#include <commctrl.h>

#include <algorithm>

namespace NiTCAD::app {

namespace {

constexpr const wchar_t* panel_class = L"NiTCADDimensionsPanel";
constexpr int id_message = 1, id_first = 100;  // a row k's edit is id_first + 4k; label +1, header +2, error +3
constexpr COLORREF error_ink = RGB(190, 30, 30);

bool derived(const dims::Dimension& d) { return d.expr.find('$') != std::string::npos; }

std::string label_text(const dims::Dimension& d) { return derived(d) ? d.label + " = " + d.value_text : d.label; }

std::string text_of(HWND h) {
    std::wstring w(static_cast<std::size_t>(GetWindowTextLengthW(h)) + 1, L'\0');
    w.resize(static_cast<std::size_t>(GetWindowTextW(h, w.data(), static_cast<int>(w.size()))));
    return win32::narrow(w);
}

}  // namespace

DimensionsPanel::DimensionsPanel(HWND parent, Host& host, int id) : host_(host) {
    win32::register_class(panel_class, proc, GetSysColorBrush(COLOR_WINDOW));
    hwnd_ = CreateWindowExW(WS_EX_CONTROLPARENT, panel_class, L"", WS_CHILD | WS_CLIPCHILDREN, 0, 0, 10, 10, parent,
                            reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)), GetModuleHandleW(nullptr), nullptr);
    SetWindowLongPtrW(hwnd_, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(this));
    message_label_ = CreateWindowExW(0, L"STATIC", L"", WS_CHILD | SS_LEFT, 0, 0, 10, 10, hwnd_,
                                     reinterpret_cast<HMENU>(static_cast<INT_PTR>(id_message)), GetModuleHandleW(nullptr), nullptr);
}

DimensionsPanel::~DimensionsPanel() {
    if (hwnd_) {
        SetWindowLongPtrW(hwnd_, GWLP_USERDATA, 0);
        for (const Row& r : controls_) RemoveWindowSubclass(r.edit, edit_proc, 1);
        DestroyWindow(hwnd_);
    }
}

void DimensionsPanel::reload() {
    const auto r = dims::read_dimensions(host_.dimensions_deck());
    if (!r) {
        enabled_ = false;
        message_ = "The deck does not parse; the dimensions wait: " + r.error();
        show_rows();
        return;
    }
    enabled_ = true;
    message_ = r->empty() ? "No named dimensions — name a value with N while dragging, or add [parameters] on the Deck tab." : "";
    std::vector<dims::Dimension> shown = *r;
    if (!group_.empty()) std::erase_if(shown, [&](const dims::Dimension& d) { return d.group != group_; });
    const bool same = shown.size() == rows_.size() &&
                      std::ranges::equal(shown, rows_, [](const auto& a, const auto& b) { return a.name == b.name && a.group == b.group; });
    rows_ = std::move(shown);
    if (!same) rebuild();
    show_rows();
}

void DimensionsPanel::rebuild() {
    for (const Row& r : controls_) {
        RemoveWindowSubclass(r.edit, edit_proc, 1);
        for (HWND h : {r.header, r.label, r.edit, r.error}) {
            if (h) DestroyWindow(h);
        }
    }
    controls_.assign(rows_.size(), Row{});
    const auto make = [&](const wchar_t* cls, DWORD style, DWORD ex, int id) {
        HWND h = CreateWindowExW(ex, cls, L"", WS_CHILD | style, 0, 0, 10, 10, hwnd_, reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)),
                                 GetModuleHandleW(nullptr), nullptr);
        SendMessageW(h, WM_SETFONT, reinterpret_cast<WPARAM>(font_), FALSE);
        return h;
    };
    std::string group;
    for (std::size_t k = 0; k < rows_.size(); ++k) {
        const int base = id_first + 4 * static_cast<int>(k);
        Row& c = controls_[k];
        if (!rows_[k].group.empty() && rows_[k].group != group) {
            c.header = make(L"STATIC", SS_LEFT, 0, base + 2);
            SendMessageW(c.header, WM_SETFONT, reinterpret_cast<WPARAM>(bold_), FALSE);
            SetWindowTextW(c.header, win32::widen(rows_[k].group).c_str());
        }
        group = rows_[k].group;
        c.label = make(L"STATIC", SS_LEFT | SS_ENDELLIPSIS, 0, base + 1);
        c.edit = make(L"EDIT", WS_TABSTOP | ES_AUTOHSCROLL, WS_EX_CLIENTEDGE, base);
        c.error = make(L"STATIC", SS_LEFT, 0, base + 3);
        SetWindowSubclass(c.edit, edit_proc, 1, reinterpret_cast<DWORD_PTR>(this));
    }
}

void DimensionsPanel::show_rows() {
    for (std::size_t k = 0; k < controls_.size() && k < rows_.size(); ++k) {
        const Row& c = controls_[k];
        SetWindowTextW(c.label, win32::widen(label_text(rows_[k])).c_str());
        if (GetFocus() != c.edit || !enabled_) SetWindowTextW(c.edit, win32::widen(rows_[k].expr).c_str());
        EnableWindow(c.edit, enabled_);
        SetWindowTextW(c.error, win32::widen(c.error_text).c_str());
    }
    SetWindowTextW(message_label_, win32::widen(message_).c_str());
    layout();
}

void DimensionsPanel::layout() {
    RECT rc{};
    GetClientRect(hwnd_, &rc);
    const int width = std::max<int>(rc.right, 1);
    const int pad = MulDiv(8, dpi_, 96), row = MulDiv(24, dpi_, 96), gap = MulDiv(4, dpi_, 96);
    const int label_w = (width - 2 * pad) * 45 / 100;
    int y = pad;
    const auto place = [&](HWND h, bool shown, int x, int w, int hgt) {
        if (!h) return;
        ShowWindow(h, shown ? SW_SHOWNA : SW_HIDE);
        if (shown) MoveWindow(h, x, y, std::max(w, 1), hgt, FALSE);
    };
    place(message_label_, !message_.empty(), pad, width - 2 * pad, 2 * row);
    if (!message_.empty()) y += 2 * row + gap;
    for (std::size_t k = 0; k < controls_.size(); ++k) {
        const Row& c = controls_[k];
        if (c.header) {
            place(c.header, true, pad, width - 2 * pad, row - gap);
            y += row;
        }
        SendMessageW(c.label, WM_SETFONT, reinterpret_cast<WPARAM>(k < rows_.size() && rows_[k].name == selected_ ? bold_ : font_), FALSE);
        place(c.label, true, pad, label_w - gap, row - gap);
        place(c.edit, true, pad + label_w, width - 2 * pad - label_w, row - gap);
        y += row;
        place(c.error, !c.error_text.empty(), pad, width - 2 * pad, row - gap);
        if (!c.error_text.empty()) y += row;
    }
    InvalidateRect(hwnd_, nullptr, TRUE);
}

int DimensionsPanel::preferred_height(int) const {
    const int pad = MulDiv(8, dpi_, 96), row = MulDiv(24, dpi_, 96), gap = MulDiv(4, dpi_, 96);
    int y = 2 * pad;
    if (!message_.empty()) y += 2 * row + gap;
    for (const Row& c : controls_) y += row * (1 + (c.header ? 1 : 0) + (c.error_text.empty() ? 0 : 1));
    return y;
}

void DimensionsPanel::set_dpi(int dpi, HFONT font, HFONT bold) {
    dpi_ = dpi;
    font_ = font;
    bold_ = bold;
    SendMessageW(message_label_, WM_SETFONT, reinterpret_cast<WPARAM>(font), FALSE);
    for (const Row& c : controls_) {
        for (HWND h : {c.label, c.edit, c.error}) SendMessageW(h, WM_SETFONT, reinterpret_cast<WPARAM>(font), FALSE);
        if (c.header) SendMessageW(c.header, WM_SETFONT, reinterpret_cast<WPARAM>(bold), FALSE);
    }
    layout();
}

void DimensionsPanel::set_group(std::string group) {
    group_ = std::move(group);
    rows_.clear();  // rebuilt for the group
    rebuild();
    reload();
}

void DimensionsPanel::select(const std::string& name) {
    selected_ = name;
    layout();
}

std::string DimensionsPanel::row_error(std::size_t row) const { return row < controls_.size() ? controls_[row].error_text : std::string(); }

HWND DimensionsPanel::row_edit(std::size_t row) const { return row < controls_.size() ? controls_[row].edit : nullptr; }

std::size_t DimensionsPanel::row_of(HWND edit) const {
    for (std::size_t k = 0; k < controls_.size(); ++k) {
        if (controls_[k].edit == edit) return k;
    }
    return controls_.size();
}

bool DimensionsPanel::commit(std::size_t row, const std::string& typed) {
    if (!enabled_ || row >= rows_.size()) return false;
    const std::string deck = host_.dimensions_deck();
    const auto out = dims::set_dimension(deck, rows_[row].name, typed);
    if (!out) {
        controls_[row].error_text = out.error();
        show_rows();
        return false;
    }
    if (const std::string why = host_.dimensions_refuse(*out); !why.empty()) {
        controls_[row].error_text = why;
        show_rows();
        return false;
    }
    controls_[row].error_text.clear();
    if (*out != deck) host_.dimensions_edited(*out);
    reload();
    return true;
}

bool DimensionsPanel::wheel(std::size_t row, int notches) {
    if (!enabled_ || row >= rows_.size()) return false;
    const std::string next = dims::step_dimension(rows_[row], notches);
    if (next == rows_[row].expr) return false;
    return commit(row, next);
}

LRESULT CALLBACK DimensionsPanel::proc(HWND hwnd, UINT message, WPARAM w, LPARAM l) {
    auto* self = reinterpret_cast<DimensionsPanel*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
    if (self) {
        switch (message) {
        case WM_COMMAND:
            if (HIWORD(w) == EN_SETFOCUS) {
                const std::size_t k = self->row_of(reinterpret_cast<HWND>(l));
                if (k < self->rows_.size()) {
                    self->selected_ = self->rows_[k].name;
                    self->layout();
                    self->host_.dimensions_selected(self->selected_);
                }
                return 0;
            }
            break;
        case WM_CTLCOLORSTATIC: {
            HDC dc = reinterpret_cast<HDC>(w);
            const auto* h = reinterpret_cast<HWND>(l);
            const bool error = std::ranges::any_of(self->controls_, [&](const Row& c) { return c.error == h; }) ||
                               (h == self->message_label_ && !self->enabled_);
            SetTextColor(dc, error ? error_ink : GetSysColor(COLOR_WINDOWTEXT));
            SetBkColor(dc, GetSysColor(COLOR_WINDOW));
            return reinterpret_cast<LRESULT>(GetSysColorBrush(COLOR_WINDOW));
        }
        case WM_SIZE:
            self->layout();
            return 0;
        default: break;
        }
    }
    return DefWindowProcW(hwnd, message, w, l);
}

LRESULT CALLBACK DimensionsPanel::edit_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l, UINT_PTR, DWORD_PTR data) {
    auto* self = reinterpret_cast<DimensionsPanel*>(data);
    switch (message) {
    case WM_GETDLGCODE: return DefSubclassProc(hwnd, message, w, l) | DLGC_WANTALLKEYS;
    case WM_KEYDOWN:
        if (w == VK_RETURN) {
            self->commit(self->row_of(hwnd), text_of(hwnd));
            return 0;
        }
        if (w == VK_ESCAPE) {
            const std::size_t k = self->row_of(hwnd);
            if (k < self->rows_.size()) SetWindowTextW(hwnd, win32::widen(self->rows_[k].expr).c_str());
            return 0;
        }
        break;
    case WM_CHAR:
        if (w == '\r' || w == 27) return 0;  // no beep
        break;
    case WM_MOUSEWHEEL:
        self->wheel(self->row_of(hwnd), GET_WHEEL_DELTA_WPARAM(w) / WHEEL_DELTA);
        return 0;
    default: break;
    }
    return DefSubclassProc(hwnd, message, w, l);
}

}  // namespace NiTCAD::app
