#include "sweep_dialog.hpp"

#include "win32.hpp"

#include "NiTCAD/app/sweep_results.hpp"
#include "NiTCAD/deck/dimension_sweep.hpp"

#include <algorithm>
#include <format>

namespace NiTCAD::app {

namespace {

constexpr const wchar_t* dialog_class = L"NiTCADSweep";
enum : int { id_title = 1, id_values, id_hint, id_message, id_run, id_cancel };

std::string text_of(HWND h) {
    std::wstring w(static_cast<std::size_t>(GetWindowTextLengthW(h)) + 1, L'\0');
    w.resize(static_cast<std::size_t>(GetWindowTextW(h, w.data(), static_cast<int>(w.size()))));
    return win32::narrow(w);
}

HWND child(HWND parent, const wchar_t* cls, const wchar_t* text, DWORD style, DWORD ex, int id) {
    return CreateWindowExW(ex, cls, text, WS_CHILD | WS_VISIBLE | style, 0, 0, 10, 10, parent,
                           reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)), GetModuleHandleW(nullptr), nullptr);
}

}  // namespace

SweepDialog::SweepDialog(HWND owner, Host& host, std::string deck, std::filesystem::path path, std::string name,
                         const std::string& label, const std::string& current, int dpi, HFONT font, HFONT bold)
    : owner_(owner), host_(host), deck_(std::move(deck)), path_(std::move(path)), name_(std::move(name)), dpi_(dpi),
      font_(font), bold_(bold) {
    win32::register_class(dialog_class, proc, GetSysColorBrush(COLOR_BTNFACE));
    const int w = MulDiv(460, dpi_, 96), h = MulDiv(230, dpi_, 96);
    RECT o{};
    GetWindowRect(owner_, &o);
    hwnd_ = CreateWindowExW(WS_EX_DLGMODALFRAME | WS_EX_CONTROLPARENT, dialog_class, win32::widen("Sweep " + name_).c_str(),
                            WS_POPUP | WS_CAPTION | WS_SYSMENU | WS_CLIPCHILDREN,
                            o.left + std::max<int>(0, (o.right - o.left - w) / 2),
                            o.top + std::max<int>(0, (o.bottom - o.top - h) / 2), w, h, owner_, nullptr,
                            GetModuleHandleW(nullptr), nullptr);
    SetWindowLongPtrW(hwnd_, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(this));
    const std::string title = label == name_ ? "Sweep " + name_ + " (now " + current + ")"
                                             : "Sweep " + label + " (" + name_ + ", now " + current + ")";
    title_ = child(hwnd_, L"STATIC", win32::widen(title).c_str(), SS_LEFT | SS_ENDELLIPSIS, 0, id_title);
    edit_ = child(hwnd_, L"EDIT", win32::widen(current).c_str(), WS_TABSTOP | ES_AUTOHSCROLL, WS_EX_CLIENTEDGE, id_values);
    hint_ = child(hwnd_, L"STATIC", L"Values: a list (20nm, 30nm, 45nm) or a range (20nm : 45nm : 5nm). One run each.",
                  SS_LEFT, 0, id_hint);
    message_label_ = child(hwnd_, L"STATIC", L"", SS_LEFT, 0, id_message);
    run_ = child(hwnd_, L"BUTTON", L"&Run", WS_TABSTOP | BS_DEFPUSHBUTTON, 0, id_run);
    cancel_ = child(hwnd_, L"BUTTON", L"Cancel", WS_TABSTOP | BS_PUSHBUTTON, 0, id_cancel);
    for (HWND c : {title_, edit_, hint_, message_label_, run_, cancel_}) SendMessageW(c, WM_SETFONT, reinterpret_cast<WPARAM>(font_), FALSE);
    SendMessageW(title_, WM_SETFONT, reinterpret_cast<WPARAM>(bold_), FALSE);
    EnableWindow(owner_, FALSE);
    read_values();
    layout();
    if (IsWindowVisible(owner_)) {
        ShowWindow(hwnd_, SW_SHOW);
        SetFocus(edit_);
        SendMessageW(edit_, EM_SETSEL, 0, -1);
    }
}

SweepDialog::~SweepDialog() {
    EnableWindow(owner_, TRUE);
    if (hwnd_) {
        SetWindowLongPtrW(hwnd_, GWLP_USERDATA, 0);
        DestroyWindow(hwnd_);
    }
    if (IsWindowVisible(owner_)) SetActiveWindow(owner_);
}

void SweepDialog::set_values(const std::string& typed) {
    SetWindowTextW(edit_, win32::widen(typed).c_str());  // EN_CHANGE reads them
    read_values();
}

void SweepDialog::read_values() {
    const auto v = typed_sweep_values(deck_, name_, text_of(edit_));
    if (v) {
        values_ = *v;
        const auto folder = deck::sweep_directory(path_, name_);
        message_ = std::format("{} run{}, one after another, into {}\\", values_.size(), values_.size() == 1 ? "" : "s",
                               folder.filename().string());
    } else {
        values_.clear();
        message_ = v.error();
    }
    SetWindowTextW(message_label_, win32::widen(message_).c_str());
    EnableWindow(run_, can_run());
}

bool SweepDialog::run() {
    if (!can_run()) return false;
    host_.sweep_run(name_, values_);
    return true;
}

void SweepDialog::cancel() { host_.sweep_closed(); }

void SweepDialog::layout() {
    RECT rc{};
    GetClientRect(hwnd_, &rc);
    const int pad = MulDiv(12, dpi_, 96), row = MulDiv(24, dpi_, 96), gap = MulDiv(6, dpi_, 96);
    const int width = rc.right - 2 * pad;
    int y = pad;
    MoveWindow(title_, pad, y, width, row - gap, FALSE);
    y += row;
    MoveWindow(edit_, pad, y, width, row - 2, FALSE);
    y += row + gap;
    MoveWindow(hint_, pad, y, width, 2 * row - gap, FALSE);
    y += 2 * row;
    MoveWindow(message_label_, pad, y, width, 2 * row - gap, FALSE);
    const int bw = MulDiv(88, dpi_, 96), by = rc.bottom - pad - row;
    MoveWindow(cancel_, rc.right - pad - bw, by, bw, row, FALSE);
    MoveWindow(run_, rc.right - pad - 2 * bw - gap, by, bw, row, FALSE);
    InvalidateRect(hwnd_, nullptr, TRUE);
}

LRESULT CALLBACK SweepDialog::proc(HWND hwnd, UINT message, WPARAM w, LPARAM l) {
    auto* self = reinterpret_cast<SweepDialog*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
    if (self) {
        switch (message) {
        case WM_COMMAND:
            switch (LOWORD(w)) {
            case id_values:
                if (HIWORD(w) == EN_CHANGE) self->read_values();
                return 0;
            case id_run: self->run(); return 0;  // the host closes the dialog from a posted message
            case id_cancel: self->cancel(); return 0;
            default: break;
            }
            break;
        case WM_CLOSE: self->cancel(); return 0;
        case WM_SIZE: self->layout(); return 0;
        default: break;
        }
    }
    return DefWindowProcW(hwnd, message, w, l);
}

}  // namespace NiTCAD::app
