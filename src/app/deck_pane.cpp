// The deck pane (deck_pane.hpp).
#include "deck_pane.hpp"

#include <richedit.h>

#include <algorithm>
#include <string>

#include "NiTCAD/app/layout.hpp"
#include "win32.hpp"

namespace NiTCAD::app {

namespace {

constexpr const wchar_t* panel_class = L"NiTCADDeckPane";
constexpr const wchar_t* gutter_class = L"NiTCADGutter";
constexpr int edit_id = 1, list_id = 2;
constexpr int errors_height_96 = 72;

// Rich Edit keeps '\r' paragraph ends; the document keeps '\n'.
std::string to_lf(std::string s) {
    std::string out;
    out.reserve(s.size());
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '\r') {
            out += '\n';
            if (i + 1 < s.size() && s[i + 1] == '\n') ++i;
        } else {
            out += s[i];
        }
    }
    return out;
}

std::string to_crlf(std::string_view s) {
    std::string out;
    out.reserve(s.size() + s.size() / 16);
    for (const char c : s) {
        if (c == '\n') out += "\r\n";
        else if (c != '\r') out += c;
    }
    return out;
}

}  // namespace

DeckPane::DeckPane(HWND parent, Host& host, int id) : host_(host) {
    LoadLibraryW(L"Msftedit.dll");
    win32::register_class(panel_class, panel_proc, GetSysColorBrush(COLOR_WINDOW));
    win32::register_class(gutter_class, gutter_proc, GetSysColorBrush(COLOR_BTNFACE));
    panel_ = CreateWindowExW(WS_EX_CONTROLPARENT, panel_class, L"", WS_CHILD | WS_CLIPCHILDREN, 0, 0, 16, 16, parent,
                             reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)), GetModuleHandleW(nullptr), this);
    gutter_ = CreateWindowExW(0, gutter_class, L"", WS_CHILD | WS_VISIBLE, 0, 0, 16, 16, panel_, nullptr,
                              GetModuleHandleW(nullptr), this);
    edit_ = CreateWindowExW(0, MSFTEDIT_CLASS, L"",
                            WS_CHILD | WS_VISIBLE | WS_VSCROLL | WS_HSCROLL | WS_TABSTOP | ES_MULTILINE |
                                ES_AUTOVSCROLL | ES_AUTOHSCROLL | ES_NOHIDESEL | ES_WANTRETURN,
                            0, 0, 16, 16, panel_, reinterpret_cast<HMENU>(static_cast<INT_PTR>(edit_id)),
                            GetModuleHandleW(nullptr), nullptr);
    SendMessageW(edit_, EM_EXLIMITTEXT, 0, 64 * 1024 * 1024);
    SendMessageW(edit_, EM_SETTEXTMODE, TM_PLAINTEXT | TM_MULTILEVELUNDO, 0);
    SendMessageW(edit_, EM_SETEVENTMASK, 0, ENM_CHANGE | ENM_UPDATE | ENM_SCROLL);
    list_ = CreateWindowExW(WS_EX_CLIENTEDGE, L"LISTBOX", L"",
                            WS_CHILD | WS_VISIBLE | WS_VSCROLL | WS_TABSTOP | LBS_NOTIFY | LBS_NOINTEGRALHEIGHT, 0, 0,
                            16, 16, panel_, reinterpret_cast<HMENU>(static_cast<INT_PTR>(list_id)),
                            GetModuleHandleW(nullptr), nullptr);
    set_dpi(win32::dpi_of(parent));
}

DeckPane::~DeckPane() {
    if (panel_) {
        SetWindowLongPtrW(panel_, GWLP_USERDATA, 0);
        SetWindowLongPtrW(gutter_, GWLP_USERDATA, 0);
        DestroyWindow(panel_);
    }
    if (code_) DeleteObject(code_);
    if (ui_) DeleteObject(ui_);
}

std::string DeckPane::text() const {
    GETTEXTLENGTHEX len{GTL_NUMCHARS | GTL_PRECISE, 1200};
    const auto n = static_cast<std::size_t>(SendMessageW(edit_, EM_GETTEXTLENGTHEX, reinterpret_cast<WPARAM>(&len), 0));
    std::wstring w(n + 1, L'\0');
    GETTEXTEX get{static_cast<DWORD>((n + 1) * sizeof(wchar_t)), GT_DEFAULT, 1200, nullptr, nullptr};
    const auto got = static_cast<std::size_t>(SendMessageW(edit_, EM_GETTEXTEX, reinterpret_cast<WPARAM>(&get),
                                                           reinterpret_cast<LPARAM>(w.data())));
    w.resize(std::min(got, n));
    return to_lf(win32::narrow(w));
}

void DeckPane::set_text(std::string_view utf8) {
    setting_ = true;
    const std::wstring w = win32::widen(to_crlf(utf8));
    SETTEXTEX set{ST_DEFAULT, 1200};
    SendMessageW(edit_, EM_SETTEXTEX, reinterpret_cast<WPARAM>(&set), reinterpret_cast<LPARAM>(w.c_str()));
    SendMessageW(edit_, EM_EMPTYUNDOBUFFER, 0, 0);
    setting_ = false;
    InvalidateRect(gutter_, nullptr, TRUE);
}

void DeckPane::set_diagnostics(const std::vector<Diagnostic>& d) {
    diagnostics_ = d;
    SendMessageW(list_, LB_RESETCONTENT, 0, 0);
    for (const Diagnostic& x : d) {
        SendMessageW(list_, LB_ADDSTRING, 0, reinterpret_cast<LPARAM>(win32::widen(x.message).c_str()));
    }
    InvalidateRect(gutter_, nullptr, TRUE);
}

void DeckPane::go_to_line(int line) {
    if (line < 1) return;
    const LRESULT start = SendMessageW(edit_, EM_LINEINDEX, static_cast<WPARAM>(line - 1), 0);
    if (start < 0) return;
    const LRESULT length = SendMessageW(edit_, EM_LINELENGTH, static_cast<WPARAM>(start), 0);
    SendMessageW(edit_, EM_SETSEL, static_cast<WPARAM>(start), static_cast<LPARAM>(start + length));
    SendMessageW(edit_, EM_SCROLLCARET, 0, 0);
    SetFocus(edit_);
    InvalidateRect(gutter_, nullptr, TRUE);
}

int DeckPane::caret_line() const {
    DWORD start = 0;
    SendMessageW(edit_, EM_GETSEL, reinterpret_cast<WPARAM>(&start), 0);
    return static_cast<int>(SendMessageW(edit_, EM_EXLINEFROMCHAR, 0, static_cast<LPARAM>(start))) + 1;
}

void DeckPane::set_dpi(int dpi) {
    dpi_ = dpi;
    if (code_) DeleteObject(code_);
    if (ui_) DeleteObject(ui_);
    code_ = win32::code_font(dpi);
    ui_ = win32::ui_font(dpi);
    SendMessageW(edit_, WM_SETFONT, reinterpret_cast<WPARAM>(code_), TRUE);
    SendMessageW(list_, WM_SETFONT, reinterpret_cast<WPARAM>(ui_), TRUE);
    layout();
}

int DeckPane::gutter_width() const {
    const int lines = static_cast<int>(SendMessageW(edit_, EM_GETLINECOUNT, 0, 0));
    int digits = 3;
    for (int n = lines; n >= 1000; n /= 10) ++digits;
    HDC dc = GetDC(gutter_);
    HGDIOBJ old = SelectObject(dc, code_);
    SIZE s{};
    GetTextExtentPoint32W(dc, L"0", 1, &s);
    SelectObject(dc, old);
    ReleaseDC(gutter_, dc);
    return s.cx * digits + scaled(12, dpi_);
}

void DeckPane::layout() {
    RECT r{};
    GetClientRect(panel_, &r);
    const int w = r.right, h = r.bottom;
    const int errors = std::min(scaled(errors_height_96, dpi_), h / 3);
    const int g = gutter_width();
    MoveWindow(gutter_, 0, 0, g, h - errors, TRUE);
    MoveWindow(edit_, g, 0, std::max(w - g, 0), h - errors, TRUE);
    MoveWindow(list_, 0, h - errors, w, errors, TRUE);
}

void DeckPane::paint_gutter(HDC dc) {
    RECT r{};
    GetClientRect(gutter_, &r);
    FillRect(dc, &r, GetSysColorBrush(COLOR_BTNFACE));
    HGDIOBJ old = SelectObject(dc, code_);
    SetBkMode(dc, TRANSPARENT);
    const int first = static_cast<int>(SendMessageW(edit_, EM_GETFIRSTVISIBLELINE, 0, 0));
    const int count = static_cast<int>(SendMessageW(edit_, EM_GETLINECOUNT, 0, 0));
    RECT er{};
    GetClientRect(edit_, &er);
    for (int line = first; line < count; ++line) {
        const LRESULT index = SendMessageW(edit_, EM_LINEINDEX, static_cast<WPARAM>(line), 0);
        POINTL p{};
        SendMessageW(edit_, EM_POSFROMCHAR, reinterpret_cast<WPARAM>(&p), index);
        if (p.y > er.bottom) break;
        const bool error = std::ranges::any_of(diagnostics_, [&](const Diagnostic& d) { return d.line == line + 1; });
        SetTextColor(dc, error ? RGB(200, 52, 44) : GetSysColor(COLOR_GRAYTEXT));
        const std::wstring n = std::to_wstring(line + 1);
        RECT cell{0, static_cast<LONG>(p.y), r.right - scaled(6, dpi_), static_cast<LONG>(p.y) + 1000};
        DrawTextW(dc, n.c_str(), static_cast<int>(n.size()), &cell, DT_RIGHT | DT_TOP | DT_SINGLELINE | DT_NOCLIP);
    }
    SelectObject(dc, old);
}

LRESULT CALLBACK DeckPane::panel_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l) {
    if (message == WM_NCCREATE) {
        SetWindowLongPtrW(hwnd, GWLP_USERDATA,
                          reinterpret_cast<LONG_PTR>(reinterpret_cast<CREATESTRUCTW*>(l)->lpCreateParams));
    }
    auto* self = reinterpret_cast<DeckPane*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
    if (!self) return DefWindowProcW(hwnd, message, w, l);
    switch (message) {
    case WM_SIZE:
        self->layout();
        return 0;
    case WM_COMMAND: {
        const int id = LOWORD(w), code = HIWORD(w);
        if (id == edit_id && code == EN_CHANGE) {
            const int g = self->gutter_width();
            RECT gr{};
            GetClientRect(self->gutter_, &gr);
            if (g != gr.right) self->layout();
            InvalidateRect(self->gutter_, nullptr, TRUE);
            if (!self->setting_) self->host_.deck_changed();
            return 0;
        }
        if (id == edit_id && (code == EN_UPDATE || code == EN_VSCROLL)) {
            InvalidateRect(self->gutter_, nullptr, TRUE);
            return 0;
        }
        if (id == list_id && code == LBN_DBLCLK) {
            const auto k = static_cast<std::size_t>(SendMessageW(self->list_, LB_GETCURSEL, 0, 0));
            if (k < self->diagnostics_.size()) self->go_to_line(self->diagnostics_[k].line);
            return 0;
        }
        break;
    }
    case WM_VKEYTOITEM:
        // Enter on an error goes to its line.
        if (reinterpret_cast<HWND>(l) == self->list_ && LOWORD(w) == VK_RETURN) {
            const auto k = static_cast<std::size_t>(SendMessageW(self->list_, LB_GETCURSEL, 0, 0));
            if (k < self->diagnostics_.size()) self->go_to_line(self->diagnostics_[k].line);
            return -2;
        }
        return -1;
    default:
        break;
    }
    return DefWindowProcW(hwnd, message, w, l);
}

LRESULT CALLBACK DeckPane::gutter_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l) {
    if (message == WM_NCCREATE) {
        SetWindowLongPtrW(hwnd, GWLP_USERDATA,
                          reinterpret_cast<LONG_PTR>(reinterpret_cast<CREATESTRUCTW*>(l)->lpCreateParams));
    }
    auto* self = reinterpret_cast<DeckPane*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
    if (self && message == WM_PAINT) {
        PAINTSTRUCT ps;
        HDC dc = BeginPaint(hwnd, &ps);
        self->paint_gutter(dc);
        EndPaint(hwnd, &ps);
        return 0;
    }
    return DefWindowProcW(hwnd, message, w, l);
}

}  // namespace NiTCAD::app
