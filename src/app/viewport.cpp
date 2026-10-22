// The view (viewport.hpp).
#include "viewport.hpp"

#include <windowsx.h>

#include <chrono>
#include <cmath>
#include <variant>

#include "win32.hpp"

namespace NiTCAD::app {

namespace {

constexpr const wchar_t* class_name = L"NiTCADViewport";
constexpr float orbit_per_pixel = 0.4f;  // degrees
constexpr float orbit_per_key = 5.0f;
constexpr float pan_per_key = 0.1f;
constexpr float zoom_per_key = 1.25f;

}  // namespace

Viewport::Viewport(HWND parent, Host& host, int id) : host_(host) {
    win32::register_class(class_name, proc, nullptr, nullptr, CS_DBLCLKS);
    hwnd_ = CreateWindowExW(0, class_name, L"", WS_CHILD | WS_VISIBLE | WS_TABSTOP | WS_CLIPSIBLINGS, 0, 0, 16, 16,
                            parent, reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)), GetModuleHandleW(nullptr), this);
}

Viewport::~Viewport() {
    release();
    if (hwnd_) {
        SetWindowLongPtrW(hwnd_, GWLP_USERDATA, 0);
        DestroyWindow(hwnd_);
    }
}

void Viewport::release() { target_.reset(); }

int Viewport::width() const noexcept {
    RECT r{};
    GetClientRect(hwnd_, &r);
    return r.right - r.left;
}

int Viewport::height() const noexcept {
    RECT r{};
    GetClientRect(hwnd_, &r);
    return r.bottom - r.top;
}

void Viewport::show(render::Renderer* renderer, const Scene* scene) {
    renderer_ = renderer;
    scene_ = scene;
    message_.clear();
    present();
}

void Viewport::show_message(std::string message) {
    scene_ = nullptr;
    message_ = std::move(message);
    target_.reset();  // GDI paints the message; a swap chain would cover it
    InvalidateRect(hwnd_, nullptr, TRUE);
}

void Viewport::prepare(render::Renderer* renderer) {
    const int w = width(), h = height();
    if (target_ || !renderer || w < 16 || h < 16) return;
    renderer_ = renderer;
    auto t = renderer->window_target(hwnd_, w, h);
    if (t) target_ = std::move(*t);
    else error_ = t.error().message;
}

void Viewport::present() {
    const int w = width(), h = height();
    if (!renderer_ || !scene_ || w < 16 || h < 16) {
        InvalidateRect(hwnd_, nullptr, TRUE);
        return;
    }
    const auto t0 = std::chrono::steady_clock::now();
    if (!target_) {
        auto t = renderer_->window_target(hwnd_, w, h);
        if (!t) {
            error_ = t.error().message;
            InvalidateRect(hwnd_, nullptr, TRUE);
            return;
        }
        target_ = std::move(*t);
    }
    auto r = std::visit([&](const auto& s) { return renderer_->present(s, *target_, w, h); }, *scene_);
    if (!r) {
        error_ = r.error().message;
        target_.reset();
        InvalidateRect(hwnd_, nullptr, TRUE);
        return;
    }
    error_.reset();
    last_ms_ = 1e3 * std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    ValidateRect(hwnd_, nullptr);
}

void Viewport::paint_message(HDC dc) {
    RECT r{};
    GetClientRect(hwnd_, &r);
    FillRect(dc, &r, GetSysColorBrush(COLOR_WINDOW));
    const std::string text = error_ ? "The view cannot be drawn: " + *error_ : message_;
    if (text.empty()) return;
    HFONT font = win32::ui_font(win32::dpi_of(hwnd_));
    HGDIOBJ old = SelectObject(dc, font);
    SetBkMode(dc, TRANSPARENT);
    SetTextColor(dc, GetSysColor(COLOR_GRAYTEXT));
    std::wstring w = win32::widen(text);
    InflateRect(&r, -24, -24);
    DrawTextW(dc, w.c_str(), static_cast<int>(w.size()), &r, DT_CENTER | DT_VCENTER | DT_WORDBREAK);
    SelectObject(dc, old);
    DeleteObject(font);
}

LRESULT CALLBACK Viewport::proc(HWND hwnd, UINT message, WPARAM w, LPARAM l) {
    if (message == WM_NCCREATE) {
        auto* self = static_cast<Viewport*>(reinterpret_cast<CREATESTRUCTW*>(l)->lpCreateParams);
        SetWindowLongPtrW(hwnd, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(self));
        self->hwnd_ = hwnd;
    }
    auto* self = reinterpret_cast<Viewport*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
    return self ? self->handle(message, w, l) : DefWindowProcW(hwnd, message, w, l);
}

LRESULT Viewport::handle(UINT message, WPARAM w, LPARAM l) {
    switch (message) {
    case WM_ERASEBKGND:
        return 1;  // the swap chain or the message paint covers everything
    case WM_PAINT: {
        PAINTSTRUCT ps;
        HDC dc = BeginPaint(hwnd_, &ps);
        if (!scene_ || error_ || !renderer_) paint_message(dc);
        EndPaint(hwnd_, &ps);
        if (scene_ && renderer_ && !error_) present();
        return 0;
    }
    case WM_SIZE:
        if (scene_ && renderer_) present();
        return 0;
    case WM_GETDLGCODE:
        return DLGC_WANTARROWS | DLGC_WANTCHARS;
    case WM_SETFOCUS:
    case WM_KILLFOCUS:
        return 0;
    case WM_LBUTTONDOWN:
    case WM_RBUTTONDOWN:
    case WM_MBUTTONDOWN:
        SetFocus(hwnd_);
        SetCapture(hwnd_);
        dragging_ = true;
        moved_ = false;
        button_ = message == WM_LBUTTONDOWN ? 1 : message == WM_RBUTTONDOWN ? 2 : 3;
        last_ = down_ = {GET_X_LPARAM(l), GET_Y_LPARAM(l)};
        return 0;
    case WM_MOUSEMOVE: {
        const POINT p{GET_X_LPARAM(l), GET_Y_LPARAM(l)};
        if (!dragging_) {
            host_.viewport_hover(static_cast<float>(p.x), static_cast<float>(p.y));
            return 0;
        }
        const int dx = p.x - last_.x, dy = p.y - last_.y;
        if (std::abs(p.x - down_.x) + std::abs(p.y - down_.y) > 3) moved_ = true;
        last_ = p;
        if (!moved_ || (dx == 0 && dy == 0)) return 0;
        if (button_ == 1 && host_.viewport_is_3d()) {
            host_.viewport_orbit(-orbit_per_pixel * static_cast<float>(dx), orbit_per_pixel * static_cast<float>(dy));
        } else {
            const float hw = std::max(width(), 1) * 0.5f, hh = std::max(height(), 1) * 0.5f;
            host_.viewport_pan(static_cast<float>(dx) / hw, -static_cast<float>(dy) / hh);
        }
        return 0;
    }
    case WM_LBUTTONUP:
    case WM_RBUTTONUP:
    case WM_MBUTTONUP: {
        const bool click = dragging_ && !moved_ && message == WM_LBUTTONUP;
        dragging_ = false;
        ReleaseCapture();
        if (click) host_.viewport_click(static_cast<float>(GET_X_LPARAM(l)), static_cast<float>(GET_Y_LPARAM(l)));
        return 0;
    }
    case WM_CAPTURECHANGED:
        dragging_ = false;
        return 0;
    case WM_LBUTTONDBLCLK:
        host_.viewport_fit();
        return 0;
    case WM_MOUSEWHEEL: {
        const float notches = static_cast<float>(GET_WHEEL_DELTA_WPARAM(w)) / WHEEL_DELTA;
        host_.viewport_zoom(std::pow(1.1f, notches));
        return 0;
    }
    case WM_KEYDOWN: {
        const bool shift = GetKeyState(VK_SHIFT) < 0;
        const bool three = host_.viewport_is_3d();
        switch (w) {
        case VK_LEFT:
        case VK_RIGHT:
        case VK_UP:
        case VK_DOWN: {
            const float sx = w == VK_LEFT ? -1.0f : w == VK_RIGHT ? 1.0f : 0.0f;
            const float sy = w == VK_UP ? 1.0f : w == VK_DOWN ? -1.0f : 0.0f;
            if (shift || !three) host_.viewport_pan(sx * pan_per_key, sy * pan_per_key);
            else host_.viewport_orbit(sx * orbit_per_key, -sy * orbit_per_key);
            return 0;
        }
        case VK_ADD:
        case VK_OEM_PLUS:
            host_.viewport_zoom(zoom_per_key);
            return 0;
        case VK_SUBTRACT:
        case VK_OEM_MINUS:
            host_.viewport_zoom(1.0f / zoom_per_key);
            return 0;
        case VK_HOME:
            host_.viewport_fit();
            return 0;
        default:
            break;
        }
        break;
    }
    default:
        break;
    }
    return DefWindowProcW(hwnd_, message, w, l);
}

}  // namespace NiTCAD::app
