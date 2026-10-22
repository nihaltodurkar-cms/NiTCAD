// The deck pane (deck_pane.hpp).
#include "deck_pane.hpp"

#include <commctrl.h>
#include <richedit.h>

#include <algorithm>
#include <string>

#include "NiTCAD/app/layout.hpp"
#include "win32.hpp"

namespace NiTCAD::app {

namespace {

constexpr const wchar_t* panel_class = L"NiTCADDeckPane";
constexpr const wchar_t* gutter_class = L"NiTCADGutter";
constexpr int edit_id = 1, list_id = 2, completion_id = 3;
constexpr int completion_rows = 8;
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
                            WS_CHILD | WS_VISIBLE | WS_CLIPSIBLINGS | WS_VSCROLL | WS_HSCROLL | WS_TABSTOP | ES_MULTILINE |
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
    // The completion list floats over the editor (which clips it out), shown at the caret.
    completion_ = CreateWindowExW(0, L"LISTBOX", L"",
                                  WS_CHILD | WS_BORDER | WS_VSCROLL | LBS_NOTIFY |
                                      LBS_NOINTEGRALHEIGHT,
                                  0, 0, 16, 16, panel_, reinterpret_cast<HMENU>(static_cast<INT_PTR>(completion_id)),
                                  GetModuleHandleW(nullptr), nullptr);
    SetWindowSubclass(edit_, edit_proc, 1, reinterpret_cast<DWORD_PTR>(this));
    set_dpi(win32::dpi_of(parent));
}

DeckPane::~DeckPane() {
    if (edit_) RemoveWindowSubclass(edit_, edit_proc, 1);
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
    SendMessageW(completion_, WM_SETFONT, reinterpret_cast<WPARAM>(code_), TRUE);
    hide_completion();
    layout();
}

std::wstring DeckPane::line_text(int line) const {
    const LRESULT start = SendMessageW(edit_, EM_LINEINDEX, static_cast<WPARAM>(line - 1), 0);
    if (start < 0) return {};
    const auto n = static_cast<std::size_t>(SendMessageW(edit_, EM_LINELENGTH, static_cast<WPARAM>(start), 0));
    std::wstring w(std::max<std::size_t>(n, 1) + 1, L'\0');
    *reinterpret_cast<WORD*>(w.data()) = static_cast<WORD>(w.size());
    const auto got = static_cast<std::size_t>(
        SendMessageW(edit_, EM_GETLINE, static_cast<WPARAM>(line - 1), reinterpret_cast<LPARAM>(w.data())));
    w.resize(std::min(got, n));
    return w;
}

std::pair<int, std::size_t> DeckPane::caret() const {
    CHARRANGE r{};
    SendMessageW(edit_, EM_EXGETSEL, 0, reinterpret_cast<LPARAM>(&r));
    int line = static_cast<int>(SendMessageW(edit_, EM_EXLINEFROMCHAR, 0, r.cpMin)) + 1;
    LRESULT start = SendMessageW(edit_, EM_LINEINDEX, static_cast<WPARAM>(line - 1), 0);
    // After a final line end Rich Edit names the line before; the empty last line starts there.
    if (const LRESULT next = SendMessageW(edit_, EM_LINEINDEX, static_cast<WPARAM>(line), 0); next >= 0 && r.cpMin >= next) {
        ++line;
        start = next;
    }
    const std::wstring w = line_text(line);
    const auto n = std::min<std::size_t>(static_cast<std::size_t>(r.cpMin - start), w.size());
    return {line, win32::narrow(w.substr(0, n)).size()};
}

LONG DeckPane::char_of(int line, std::size_t column) const {
    const LRESULT start = SendMessageW(edit_, EM_LINEINDEX, static_cast<WPARAM>(line - 1), 0);
    const std::string utf8 = win32::narrow(line_text(line));
    return static_cast<LONG>(start) + static_cast<LONG>(win32::widen(utf8.substr(0, std::min(column, utf8.size()))).size());
}

void DeckPane::set_caret(int line, std::size_t column) {
    const LONG c = char_of(line, column);
    SendMessageW(edit_, EM_SETSEL, static_cast<WPARAM>(c), static_cast<LPARAM>(c));
}

bool DeckPane::show_completion(bool forced) {
    const auto [line, column] = caret();
    auto c = complete(text(), line, column);
    if (!c || c->items.empty()) {
        hide_completion();
        return false;
    }
    if (!forced && c->from == column) {
        // Nothing typed of the word: open only right after '[' or '='.
        const std::string before = win32::narrow(line_text(line)).substr(0, column);
        const auto last = before.find_last_not_of(" \t");
        if (last == std::string::npos || (before[last] != '[' && before[last] != '=')) {
            hide_completion();
            return false;
        }
    }
    SendMessageW(completion_, WM_SETREDRAW, FALSE, 0);
    SendMessageW(completion_, LB_RESETCONTENT, 0, 0);
    // The list is in the code font: the help aligned by padding the names.
    int longest = 0;
    HDC dc = GetDC(completion_);
    HGDIOBJ old_font = SelectObject(dc, code_);
    const auto text_width = [&](const std::wstring& s) {
        SIZE size{};
        GetTextExtentPoint32W(dc, s.c_str(), static_cast<int>(s.size()), &size);
        return static_cast<int>(size.cx);
    };
    std::size_t widest = 0;
    for (const Candidate& item : c->items) widest = std::max(widest, item.label.size());
    for (const Candidate& item : c->items) {
        const std::wstring row = win32::widen(item.label + std::string(widest + 2 - item.label.size(), ' ') + item.help);
        longest = std::max(longest, text_width(row));
        SendMessageW(completion_, LB_ADDSTRING, 0, reinterpret_cast<LPARAM>(row.c_str()));
    }
    SelectObject(dc, old_font);
    ReleaseDC(completion_, dc);
    SendMessageW(completion_, LB_SETCURSEL, 0, 0);
    SendMessageW(completion_, WM_SETREDRAW, TRUE, 0);
    // Below the word's start.
    POINTL p{};
    SendMessageW(edit_, EM_POSFROMCHAR, reinterpret_cast<WPARAM>(&p), char_of(line, c->from));
    const int row = static_cast<int>(SendMessageW(completion_, LB_GETITEMHEIGHT, 0, 0));
    const int rows = std::min<int>(static_cast<int>(c->items.size()), completion_rows);
    RECT er{};
    GetWindowRect(edit_, &er);
    MapWindowPoints(HWND_DESKTOP, panel_, reinterpret_cast<POINT*>(&er), 2);
    const int line_height = row;
    const int w = std::min<int>(longest + GetSystemMetricsForDpi(SM_CXVSCROLL, static_cast<UINT>(dpi_)) + scaled(8, dpi_),
                                er.right - er.left),
              h = rows * row + 4;
    int x = er.left + static_cast<int>(p.x), y = er.top + static_cast<int>(p.y) + line_height;
    x = std::max<int>(er.left, std::min<int>(x, er.right - w));
    if (y + h > er.bottom) y = std::max<int>(er.top, er.top + static_cast<int>(p.y) - h);
    SetWindowPos(completion_, HWND_TOP, x, y, w, h, SWP_SHOWWINDOW | SWP_NOACTIVATE);
    completion_shown_ = std::move(c);
    return true;
}

void DeckPane::hide_completion() {
    if (completion_) ShowWindow(completion_, SW_HIDE);
    completion_shown_.reset();
}

int DeckPane::completion_selection() const {
    return completion_shown_ ? static_cast<int>(SendMessageW(completion_, LB_GETCURSEL, 0, 0)) : -1;
}

void DeckPane::accept_completion(int index) {
    if (!completion_shown_ || index < 0 || static_cast<std::size_t>(index) >= completion_shown_->items.size()) return;
    const Completion c = *completion_shown_;
    hide_completion();
    const LONG from = char_of(c.line, c.from), to = char_of(c.line, c.to);
    SendMessageW(edit_, EM_SETSEL, static_cast<WPARAM>(from), static_cast<LPARAM>(to));
    const std::wstring insert = win32::widen(c.items[static_cast<std::size_t>(index)].insert);
    SendMessageW(edit_, EM_REPLACESEL, TRUE, reinterpret_cast<LPARAM>(insert.c_str()));
    SetFocus(edit_);
    // A key's value has choices to offer next.
    if (!insert.empty() && (insert.back() == L' ')) show_completion(false);
}

void DeckPane::after_typing(wchar_t typed) {
    const bool word = (typed >= L'a' && typed <= L'z') || (typed >= L'A' && typed <= L'Z') ||
                      (typed >= L'0' && typed <= L'9') || typed == L'_' || typed == L'.';
    if (word || typed == L'[' || typed == L'=' || typed == L' ' || (typed == L'\b' && completion_shown_)) {
        show_completion(false);
    } else {
        hide_completion();
    }
}

void DeckPane::paint_underlines() {
    if (diagnostics_.empty()) return;
    HDC dc = GetDC(edit_);
    RECT er{};
    GetClientRect(edit_, &er);
    HPEN pen = CreatePen(PS_SOLID, std::max(1, scaled(1, dpi_)), RGB(200, 52, 44));
    HGDIOBJ old = SelectObject(dc, pen);
    TEXTMETRICW tm{};
    HGDIOBJ old_font = SelectObject(dc, code_);
    GetTextMetricsW(dc, &tm);
    SelectObject(dc, old_font);
    const int step = std::max(2, scaled(2, dpi_));
    for (const Diagnostic& d : diagnostics_) {
        if (d.line < 1) continue;
        const std::wstring w = line_text(d.line);
        const auto first = w.find_first_not_of(L" \t");
        if (first == std::wstring::npos) continue;
        const auto last = w.find_last_not_of(L" \t");
        const LRESULT start = SendMessageW(edit_, EM_LINEINDEX, static_cast<WPARAM>(d.line - 1), 0);
        POINTL a{}, b{};
        SendMessageW(edit_, EM_POSFROMCHAR, reinterpret_cast<WPARAM>(&a), start + static_cast<LPARAM>(first));
        SendMessageW(edit_, EM_POSFROMCHAR, reinterpret_cast<WPARAM>(&b), start + static_cast<LPARAM>(last + 1));
        if (a.y < er.top - tm.tmHeight || a.y > er.bottom) continue;
        if (b.x <= a.x) b.x = a.x + tm.tmAveCharWidth * static_cast<LONG>(last + 1 - first);
        // A wavy line under the text.
        const int y = static_cast<int>(a.y) + tm.tmHeight - step;
        std::vector<POINT> zig;
        for (int x = static_cast<int>(a.x), k = 0; x <= static_cast<int>(b.x); x += step, ++k) {
            zig.push_back({x, y + ((k & 1) ? step : 0)});
        }
        if (zig.size() > 1) Polyline(dc, zig.data(), static_cast<int>(zig.size()));
    }
    SelectObject(dc, old);
    DeleteObject(pen);
    ReleaseDC(edit_, dc);
}

LRESULT CALLBACK DeckPane::edit_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l, UINT_PTR, DWORD_PTR data) {
    auto* self = reinterpret_cast<DeckPane*>(data);
    switch (message) {
    case WM_KEYDOWN:
        if (self->completion_shown_) {
            const int n = static_cast<int>(self->completion_shown_->items.size());
            const int at = self->completion_selection();
            switch (w) {
            case VK_DOWN:
                SendMessageW(self->completion_, LB_SETCURSEL, static_cast<WPARAM>(std::min(at + 1, n - 1)), 0);
                return 0;
            case VK_UP:
                SendMessageW(self->completion_, LB_SETCURSEL, static_cast<WPARAM>(std::max(at - 1, 0)), 0);
                return 0;
            case VK_RETURN:
            case VK_TAB:
                self->swallow_char_ = true;
                self->accept_completion(at);
                return 0;
            case VK_ESCAPE:
                self->swallow_char_ = true;
                self->hide_completion();
                return 0;
            case VK_LEFT:
            case VK_RIGHT:
            case VK_HOME:
            case VK_END:
            case VK_PRIOR:
            case VK_NEXT:
                self->hide_completion();
                break;
            default:
                break;
            }
        }
        break;
    case WM_GETDLGCODE:
        // While the list is shown, Tab, Enter and Escape are the list's, not the dialog navigation's.
        if (self->completion_shown_) return DefSubclassProc(hwnd, message, w, l) | DLGC_WANTALLKEYS;
        break;
    case WM_CHAR: {
        // The character of an Enter, Tab or Escape the list took.
        if (self->swallow_char_) {
            self->swallow_char_ = false;
            if (w == VK_RETURN || w == VK_TAB || w == VK_ESCAPE) return 0;
        }
        const LRESULT r = DefSubclassProc(hwnd, message, w, l);
        self->after_typing(static_cast<wchar_t>(w));
        return r;
    }
    case WM_KILLFOCUS:
        if (reinterpret_cast<HWND>(w) != self->completion_) self->hide_completion();
        break;
    case WM_PAINT: {
        const LRESULT r = DefSubclassProc(hwnd, message, w, l);
        self->paint_underlines();
        return r;
    }
    default:
        break;
    }
    return DefSubclassProc(hwnd, message, w, l);
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
        if (id == completion_id && code == LBN_DBLCLK) {
            self->accept_completion(self->completion_selection());
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
