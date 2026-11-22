#include "new_device_dialog.hpp"

#include "win32.hpp"

#include <algorithm>
#include <fstream>

namespace NiTCAD::app {

namespace fs = std::filesystem;

namespace {

constexpr const wchar_t* dialog_class = L"NiTCADNewDevice";
constexpr const wchar_t* thumb_class = L"NiTCADNewDeviceThumb";
enum : int {
    id_list = 1,
    id_thumb,
    id_message,
    id_folder_label,
    id_folder,
    id_browse,
    id_name_label,
    id_name,
    id_create,
    id_cancel,
    id_panel
};

std::string read_file(const fs::path& p) {
    std::ifstream f(p, std::ios::binary);
    return {std::istreambuf_iterator<char>(f), std::istreambuf_iterator<char>()};
}

std::string text_of(HWND h) {
    std::wstring w(static_cast<std::size_t>(GetWindowTextLengthW(h)) + 1, L'\0');
    w.resize(static_cast<std::size_t>(GetWindowTextW(h, w.data(), static_cast<int>(w.size()))));
    return win32::narrow(w);
}

COLORREF material_color(const std::string& m) {
    if (m == "silicon") return RGB(150, 150, 155);
    if (m == "sio2") return RGB(170, 210, 235);
    if (m == "insulator") return RGB(205, 230, 240);
    return RGB(220, 200, 150);
}

HWND child(HWND parent, const wchar_t* cls, const wchar_t* text, DWORD style, DWORD ex, int id) {
    return CreateWindowExW(ex, cls, text, WS_CHILD | WS_VISIBLE | style, 0, 0, 10, 10, parent,
                           reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)), GetModuleHandleW(nullptr), nullptr);
}

}  // namespace

void NewDeviceDialog::PanelHost::dimensions_edited(const std::string& text) {
    d->text_ = text;
    d->thumb_ = devices::thumbnail(d->text_);
    d->message_ = d->thumb_.error;
    d->show_message();
}

// A value the build refuses is not written: the build's message on the row.
std::string NewDeviceDialog::PanelHost::dimensions_refuse(const std::string& text) { return devices::thumbnail(text).error; }

NewDeviceDialog::NewDeviceDialog(HWND owner, Host& host, std::vector<devices::Family> families, fs::path folder, int dpi,
                                 HFONT font, HFONT bold)
    : owner_(owner), host_(host), families_(std::move(families)), dpi_(dpi), font_(font), bold_(bold) {
    win32::register_class(dialog_class, proc, GetSysColorBrush(COLOR_BTNFACE));
    win32::register_class(thumb_class, thumb_proc, GetSysColorBrush(COLOR_WINDOW));
    const int w = MulDiv(760, dpi_, 96), h = MulDiv(520, dpi_, 96);
    RECT o{};
    GetWindowRect(owner_, &o);
    hwnd_ = CreateWindowExW(WS_EX_DLGMODALFRAME | WS_EX_CONTROLPARENT, dialog_class, L"New device",
                            WS_POPUP | WS_CAPTION | WS_SYSMENU | WS_CLIPCHILDREN, o.left + std::max<int>(0, (o.right - o.left - w) / 2),
                            o.top + std::max<int>(0, (o.bottom - o.top - h) / 2), w, h, owner_, nullptr, GetModuleHandleW(nullptr), nullptr);
    SetWindowLongPtrW(hwnd_, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(this));
    list_ = child(hwnd_, L"LISTBOX", L"", WS_TABSTOP | WS_VSCROLL | LBS_NOTIFY | LBS_NOINTEGRALHEIGHT, WS_EX_CLIENTEDGE, id_list);
    thumb_wnd_ = child(hwnd_, thumb_class, L"", 0, WS_EX_CLIENTEDGE, id_thumb);
    SetWindowLongPtrW(thumb_wnd_, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(this));
    message_label_ = child(hwnd_, L"STATIC", L"", SS_LEFT, 0, id_message);
    folder_label_ = child(hwnd_, L"STATIC", L"Folder", SS_LEFT, 0, id_folder_label);
    folder_ = child(hwnd_, L"EDIT", win32::widen(folder.string()).c_str(), WS_TABSTOP | ES_AUTOHSCROLL, WS_EX_CLIENTEDGE, id_folder);
    browse_ = child(hwnd_, L"BUTTON", L"&Browse...", WS_TABSTOP | BS_PUSHBUTTON, 0, id_browse);
    name_label_ = child(hwnd_, L"STATIC", L"Name", SS_LEFT, 0, id_name_label);
    name_ = child(hwnd_, L"EDIT", L"", WS_TABSTOP | ES_AUTOHSCROLL, WS_EX_CLIENTEDGE, id_name);
    create_ = child(hwnd_, L"BUTTON", L"&Create", WS_TABSTOP | BS_DEFPUSHBUTTON, 0, id_create);
    cancel_ = child(hwnd_, L"BUTTON", L"Cancel", WS_TABSTOP | BS_PUSHBUTTON, 0, id_cancel);
    for (HWND c : {list_, message_label_, folder_label_, folder_, browse_, name_label_, name_, create_, cancel_}) {
        SendMessageW(c, WM_SETFONT, reinterpret_cast<WPARAM>(font_), FALSE);
    }
    for (const devices::Family& f : families_) {
        const std::string row = f.example.title + (f.dimension ? "" : "  (does not build)");
        SendMessageW(list_, LB_ADDSTRING, 0, reinterpret_cast<LPARAM>(win32::widen(row).c_str()));
    }
    panel_host_.d = this;
    panel_ = std::make_unique<DimensionsPanel>(hwnd_, panel_host_, id_panel);
    panel_->set_dpi(dpi_, font_, bold_);
    ShowWindow(panel_->hwnd(), SW_SHOWNA);
    EnableWindow(owner_, FALSE);
    if (!families_.empty()) select(0);
    layout();
    if (IsWindowVisible(owner_)) ShowWindow(hwnd_, SW_SHOW);
}

NewDeviceDialog::~NewDeviceDialog() {
    EnableWindow(owner_, TRUE);
    if (hwnd_) {
        SetWindowLongPtrW(hwnd_, GWLP_USERDATA, 0);
        SetWindowLongPtrW(thumb_wnd_, GWLP_USERDATA, 0);
        panel_.reset();
        DestroyWindow(hwnd_);
    }
    if (IsWindowVisible(owner_)) SetActiveWindow(owner_);
}

void NewDeviceDialog::select(std::size_t family) {
    if (family >= families_.size()) return;
    selected_ = family;
    SendMessageW(list_, LB_SETCURSEL, family, 0);
    text_ = read_file(families_[family].example.path);
    thumb_ = devices::thumbnail(text_);
    message_ = thumb_.error;
    set_name(families_[family].example.path.stem().string());
    const auto rows = dims::read_dimensions(text_);
    panel_->set_group(rows && !rows->empty() ? rows->front().group : std::string());
    show_message();
    layout();
}

void NewDeviceDialog::set_folder(const fs::path& folder) { SetWindowTextW(folder_, folder.wstring().c_str()); }

fs::path NewDeviceDialog::folder() const { return fs::path(win32::widen(text_of(folder_))); }

void NewDeviceDialog::set_name(const std::string& name) { SetWindowTextW(name_, win32::widen(name).c_str()); }

std::string NewDeviceDialog::name() const { return text_of(name_); }

bool NewDeviceDialog::create() {
    if (!can_create()) return false;
    const auto made = devices::create_device(text_, folder(), name());
    if (!made) {
        message_ = made.error();
        show_message();
        return false;
    }
    host_.new_device_created(*made);
    return true;
}

void NewDeviceDialog::cancel() { host_.new_device_closed(); }

void NewDeviceDialog::show_message() {
    SetWindowTextW(message_label_, win32::widen(message_).c_str());
    EnableWindow(create_, can_create());
    InvalidateRect(thumb_wnd_, nullptr, TRUE);
}

void NewDeviceDialog::layout() {
    RECT rc{};
    GetClientRect(hwnd_, &rc);
    const int pad = MulDiv(10, dpi_, 96), row = MulDiv(24, dpi_, 96), button = MulDiv(90, dpi_, 96);
    const int list_w = MulDiv(220, dpi_, 96), thumb_h = MulDiv(170, dpi_, 96), label_w = MulDiv(50, dpi_, 96);
    const int x1 = pad + list_w + pad, right = rc.right - pad, bottom = rc.bottom - pad;
    const int foot = bottom - 3 * row - 2 * pad;  // the folder, name and buttons rows below
    MoveWindow(list_, pad, pad, list_w, std::max(row, foot - pad), FALSE);
    MoveWindow(thumb_wnd_, x1, pad, std::max(1, right - x1), thumb_h, FALSE);
    MoveWindow(message_label_, x1, pad + thumb_h + pad / 2, std::max(1, right - x1), row, FALSE);
    const int py = pad + thumb_h + pad / 2 + row;
    MoveWindow(panel_->hwnd(), x1, py, std::max(1, right - x1), std::max(row, foot - py), FALSE);
    panel_->layout();
    int y = foot + pad;
    MoveWindow(folder_label_, pad, y + 3, label_w, row, FALSE);
    MoveWindow(folder_, pad + label_w, y, std::max(1, right - 2 * pad - label_w - button), row, FALSE);
    MoveWindow(browse_, right - button, y, button, row, FALSE);
    y += row + pad / 2;
    MoveWindow(name_label_, pad, y + 3, label_w, row, FALSE);
    MoveWindow(name_, pad + label_w, y, MulDiv(260, dpi_, 96), row, FALSE);
    y += row + pad;
    MoveWindow(cancel_, right - button, y, button, row, FALSE);
    MoveWindow(create_, right - 2 * button - pad, y, button, row, FALSE);
    InvalidateRect(hwnd_, nullptr, TRUE);
}

void NewDeviceDialog::paint_thumbnail(HDC dc, const RECT& r) const {
    FillRect(dc, &r, GetSysColorBrush(COLOR_WINDOW));
    if (!thumb_.error.empty() || thumb_.boxes.empty()) return;
    const double w = thumb_.x1 - thumb_.x0, h = thumb_.y1 - thumb_.y0;
    if (!(w > 0.0) || !(h > 0.0)) return;
    const double m = MulDiv(8, dpi_, 96);
    const double aw = r.right - r.left - 2.0 * m, ah = r.bottom - r.top - 2.0 * m;
    if (aw <= 0.0 || ah <= 0.0) return;
    const bool bar = thumb_.y0 == 0.0 && thumb_.y1 == 1.0;  // a 1D deck: a bar across
    const double sx = aw / w, sy = bar ? std::min(ah, 40.0 * dpi_ / 96.0) / h : std::min(ah / h, aw / w);
    const double s = bar ? sx : sy;  // 2D/3D: one scale, the aspect kept
    const double dw = w * s, dh = bar ? h * sy : h * s;
    const double ox = r.left + m + (aw - dw) / 2, oy = r.top + m + (ah - dh) / 2;
    const auto px = [&](double x) { return static_cast<LONG>(ox + (x - thumb_.x0) * s); };
    const auto py = [&](double y) { return static_cast<LONG>(oy + dh - (y - thumb_.y0) * (bar ? sy : s)); };
    for (const devices::ThumbBox& b : thumb_.boxes) {
        const RECT q{px(b.x0), py(b.y1), px(b.x1), py(b.y0)};
        HBRUSH brush = CreateSolidBrush(material_color(b.material));
        FillRect(dc, &q, brush);
        DeleteObject(brush);
    }
    const RECT frame{px(thumb_.x0), py(thumb_.y1), px(thumb_.x1), py(thumb_.y0)};
    FrameRect(dc, &frame, GetSysColorBrush(COLOR_WINDOWTEXT));
}

LRESULT CALLBACK NewDeviceDialog::thumb_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l) {
    if (message == WM_PAINT) {
        PAINTSTRUCT ps;
        HDC dc = BeginPaint(hwnd, &ps);
        RECT r{};
        GetClientRect(hwnd, &r);
        if (const auto* self = reinterpret_cast<NewDeviceDialog*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA))) self->paint_thumbnail(dc, r);
        EndPaint(hwnd, &ps);
        return 0;
    }
    return DefWindowProcW(hwnd, message, w, l);
}

LRESULT CALLBACK NewDeviceDialog::proc(HWND hwnd, UINT message, WPARAM w, LPARAM l) {
    auto* self = reinterpret_cast<NewDeviceDialog*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
    if (self) {
        switch (message) {
        case WM_COMMAND:
            switch (LOWORD(w)) {
            case id_list:
                if (HIWORD(w) == LBN_SELCHANGE) {
                    const auto k = SendMessageW(self->list_, LB_GETCURSEL, 0, 0);
                    if (k >= 0) self->select(static_cast<std::size_t>(k));
                }
                return 0;
            case id_browse: {
                const fs::path p = win32::pick_folder(hwnd, L"The folder for the new deck", self->folder());
                if (!p.empty()) self->set_folder(p);
                return 0;
            }
            case id_create: self->create(); return 0;  // the host closes the dialog from a posted message
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
