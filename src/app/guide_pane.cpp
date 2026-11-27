#include "guide_pane.hpp"

#include "win32.hpp"

#include <algorithm>
#include <format>

#include "NiTCAD/app/layout.hpp"

namespace NiTCAD::app {

namespace {

constexpr const wchar_t* pane_class = L"NiTCADGuidePane";
constexpr int id_picker = 1, id_heading = 2, id_prose = 3, id_action = 4, id_report_label = 5, id_report = 6,
              id_back = 7, id_next = 8, id_position = 9;
constexpr COLORREF error_ink = RGB(190, 30, 30);

// An edit control shows "\r\n" line breaks.
std::string crlf(const std::string& s) {
    std::string out;
    out.reserve(s.size() + 16);
    for (const char c : s) {
        if (c == '\n') out += '\r';
        out += c;
    }
    return out;
}

HWND child(HWND parent, const wchar_t* cls, const wchar_t* text, DWORD style, int id, DWORD ex = 0) {
    return CreateWindowExW(ex, cls, text, WS_CHILD | WS_VISIBLE | style, 0, 0, 10, 10, parent,
                           reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)), GetModuleHandleW(nullptr), nullptr);
}

}  // namespace

GuidePane::GuidePane(HWND parent, Host& host, int id) : host_(host) {
    win32::register_class(pane_class, proc, GetSysColorBrush(COLOR_WINDOW));
    hwnd_ = CreateWindowExW(WS_EX_CONTROLPARENT, pane_class, L"", WS_CHILD | WS_CLIPCHILDREN, 0, 0, 10, 10, parent,
                            reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)), GetModuleHandleW(nullptr), nullptr);
    SetWindowLongPtrW(hwnd_, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(this));
    picker_ = child(hwnd_, L"COMBOBOX", L"", WS_TABSTOP | WS_VSCROLL | CBS_DROPDOWNLIST, id_picker);
    heading_ = child(hwnd_, L"STATIC", L"", SS_LEFT | SS_NOPREFIX, id_heading);
    prose_ = child(hwnd_, L"EDIT", L"", WS_TABSTOP | WS_VSCROLL | ES_MULTILINE | ES_READONLY | ES_AUTOVSCROLL, id_prose);
    action_ = child(hwnd_, L"BUTTON", L"", WS_TABSTOP | BS_PUSHBUTTON, id_action);
    report_label_ = child(hwnd_, L"STATIC", L"What happened", SS_LEFT, id_report_label);
    report_ = child(hwnd_, L"EDIT", L"", WS_TABSTOP | WS_VSCROLL | ES_MULTILINE | ES_READONLY | ES_AUTOVSCROLL, id_report);
    back_ = child(hwnd_, L"BUTTON", L"< &Back", WS_TABSTOP | BS_PUSHBUTTON, id_back);
    position_ = child(hwnd_, L"STATIC", L"", SS_CENTER | SS_NOPREFIX, id_position);
    next_ = child(hwnd_, L"BUTTON", L"&Next >", WS_TABSTOP | BS_PUSHBUTTON, id_next);
    show_step();
}

GuidePane::~GuidePane() {
    if (hwnd_) {
        SetWindowLongPtrW(hwnd_, GWLP_USERDATA, 0);
        DestroyWindow(hwnd_);
    }
    if (back_brush_) DeleteObject(back_brush_);
}

void GuidePane::set_lessons(std::vector<Lesson> lessons) {
    lessons_ = std::move(lessons);
    SendMessageW(picker_, CB_RESETCONTENT, 0, 0);
    for (const Lesson& l : lessons_) {
        const std::string item = std::format("{} \xE2\x80\x94 {}. {}", to_string(l.track), l.order, l.title);  // "Basics — 1. ..."
        SendMessageW(picker_, CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(win32::widen(item).c_str()));
    }
    // None shown until one is picked (and its deck opened): a step must not act on whatever deck is open.
    lesson_ = lessons_.size();
    step_ = 0;
    SendMessageW(picker_, CB_SETCURSEL, static_cast<WPARAM>(-1), 0);
    show_step();
}

void GuidePane::show_lesson(std::size_t lesson) {
    if (lesson >= lessons_.size()) return;
    lesson_ = lesson;
    step_ = 0;
    SendMessageW(picker_, CB_SETCURSEL, lesson, 0);
    show_step();
}

void GuidePane::go(int delta) {
    if (lesson_ >= lessons_.size() || lessons_[lesson_].steps.empty()) return;
    const auto last = static_cast<long long>(lessons_[lesson_].steps.size()) - 1;
    const auto to = static_cast<std::size_t>(std::clamp(static_cast<long long>(step_) + delta, 0LL, last));
    if (to == step_) return;
    step_ = to;
    show_step();
}

void GuidePane::set_report(const std::string& text) {
    error_ = false;
    win32::set_text(report_, crlf(text));
}

void GuidePane::set_error(const std::string& text) {
    error_ = true;
    win32::set_text(report_, crlf(text));
}

void GuidePane::set_busy(bool busy) {
    busy_ = busy;
    const LessonStep* s = current();
    EnableWindow(action_, s && s->action.kind != LessonAction::Kind::none && !busy_);
}

void GuidePane::set_colors(std::optional<Colors> colors) {
    colors_ = colors;
    if (back_brush_) DeleteObject(back_brush_);
    back_brush_ = colors_ ? CreateSolidBrush(colors_->back) : nullptr;
    RedrawWindow(hwnd_, nullptr, nullptr, RDW_INVALIDATE | RDW_ERASE | RDW_ALLCHILDREN);
}

void GuidePane::set_dpi(int dpi, HFONT font, HFONT bold) {
    dpi_ = dpi;
    font_ = font;
    bold_ = bold;
    for (HWND h : {picker_, prose_, action_, report_label_, report_, back_, next_, position_}) {
        SendMessageW(h, WM_SETFONT, reinterpret_cast<WPARAM>(font), FALSE);
    }
    SendMessageW(heading_, WM_SETFONT, reinterpret_cast<WPARAM>(bold ? bold : font), FALSE);
    layout();
}

void GuidePane::act() {
    const LessonStep* s = current();
    if (!s || s->action.kind == LessonAction::Kind::none || busy_) return;
    host_.guide_act(*s);
}

const LessonStep* GuidePane::current() const noexcept {
    if (lesson_ >= lessons_.size() || step_ >= lessons_[lesson_].steps.size()) return nullptr;
    return &lessons_[lesson_].steps[step_];
}

std::string GuidePane::action_label() const {
    const LessonStep* s = current();
    if (!s) return {};
    const LessonAction& a = s->action;
    switch (a.kind) {
        case LessonAction::Kind::none: return {};
        case LessonAction::Kind::run: return "Run";
        case LessonAction::Kind::sweep: return "Sweep " + a.name + " = " + a.value;
        case LessonAction::Kind::set:
        case LessonAction::Kind::edit: return "Set " + a.name + " = " + a.value;
    }
    return {};
}

std::string GuidePane::heading_text() const { return win32::text_of(heading_); }
std::string GuidePane::prose_text() const { return win32::text_of(prose_); }
std::string GuidePane::report_text() const { return win32::text_of(report_); }
std::string GuidePane::position_text() const { return win32::text_of(position_); }

void GuidePane::show_step() {
    const LessonStep* s = current();
    win32::set_text(heading_, s ? s->heading : (lessons_.empty() ? "No lessons found" : "Pick a lesson"));
    layout();  // the heading's height follows its text
    win32::set_text(prose_, s ? crlf(s->prose) : std::string());
    const std::string label = action_label();
    win32::set_text(action_, label.empty() ? std::string("(nothing to do: read, then Next)") : label);
    set_busy(busy_);
    error_ = false;
    win32::set_text(report_, "");
    const std::size_t n = lesson_ < lessons_.size() ? lessons_[lesson_].steps.size() : 0;
    win32::set_text(position_, s ? std::format("Step {} of {}", step_ + 1, n) : std::string());
    EnableWindow(back_, s && step_ > 0);
    EnableWindow(next_, s && step_ + 1 < n);
}

void GuidePane::layout() {
    RECT rc{};
    GetClientRect(hwnd_, &rc);
    const int w = rc.right, h = rc.bottom;
    const auto px = [&](int at_96) { return scaled(at_96, dpi_); };
    const int m = px(8), gap = px(6), row = px(26);
    const int inner = std::max(w - 2 * m, 0);
    int y = m;
    MoveWindow(picker_, m, y, inner, px(300), TRUE);  // the drop-down's height; the field itself is one row
    y += row + gap;
    // The heading as tall as its wrapped text at this width (one row at least).
    int heading_h = px(22);
    if (HDC dc = GetDC(heading_)) {
        HGDIOBJ old = SelectObject(dc, reinterpret_cast<HFONT>(SendMessageW(heading_, WM_GETFONT, 0, 0)));
        const std::wstring text = win32::widen(heading_text());
        RECT r{0, 0, inner, 0};
        DrawTextW(dc, text.c_str(), static_cast<int>(text.size()), &r, DT_CALCRECT | DT_WORDBREAK | DT_NOPREFIX);
        SelectObject(dc, old);
        ReleaseDC(heading_, dc);
        heading_h = std::max(heading_h, static_cast<int>(r.bottom) + px(4));
    }
    MoveWindow(heading_, m, y, inner, heading_h, TRUE);
    y += heading_h + gap;
    // From the bottom up: Back / position / Next, the report, its label, the action.
    const int bottom_row = std::max(h - m - row, y);
    const int report_h = px(84);
    const int report_y = std::max(bottom_row - gap - report_h, y);
    const int label_y = std::max(report_y - px(20), y);
    const int action_y = std::max(label_y - gap - row, y);
    const int prose_h = std::max(action_y - gap - y, 0);
    MoveWindow(prose_, m, y, inner, prose_h, TRUE);
    MoveWindow(action_, m, action_y, inner, row, TRUE);
    MoveWindow(report_label_, m, label_y, inner, px(20), TRUE);
    MoveWindow(report_, m, report_y, inner, std::max(bottom_row - gap - report_y, 0), TRUE);
    const int button = std::min(px(80), inner / 3);
    MoveWindow(back_, m, bottom_row, button, row, TRUE);
    MoveWindow(position_, m + button, bottom_row + px(5), std::max(inner - 2 * button, 0), row - px(5), TRUE);
    MoveWindow(next_, m + inner - button, bottom_row, button, row, TRUE);
}

LRESULT CALLBACK GuidePane::proc(HWND hwnd, UINT message, WPARAM w, LPARAM l) {
    auto* self = reinterpret_cast<GuidePane*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
    if (self) {
        switch (message) {
        case WM_COMMAND: {
            HWND from = reinterpret_cast<HWND>(l);
            if (HIWORD(w) == BN_CLICKED) {
                if (from == self->action_) self->act();
                else if (from == self->back_) self->go(-1);
                else if (from == self->next_) self->go(1);
                return 0;
            }
            // A lesson opens when the list closes on it, not on each arrow key in the closed field (each open asks
            // to discard edits and copies a deck); arrows there put the field back to the lesson shown.
            if (from == self->picker_ && (HIWORD(w) == CBN_CLOSEUP || HIWORD(w) == CBN_SELCHANGE)) {
                const LRESULT k = SendMessageW(self->picker_, CB_GETCURSEL, 0, 0);
                const bool dropped = SendMessageW(self->picker_, CB_GETDROPPEDSTATE, 0, 0) != 0;
                if (HIWORD(w) == CBN_SELCHANGE && dropped) return 0;  // browsing the open list
                if (HIWORD(w) == CBN_CLOSEUP && k >= 0 && static_cast<std::size_t>(k) != self->lesson_) {
                    self->host_.guide_open(static_cast<std::size_t>(k));
                }
                // The host shows the lesson when it opened it; otherwise the picker goes back to the one shown.
                SendMessageW(self->picker_, CB_SETCURSEL, self->current() ? self->lesson_ : static_cast<WPARAM>(-1), 0);
                return 0;
            }
            break;
        }
        case WM_CTLCOLORSTATIC: {
            HDC dc = reinterpret_cast<HDC>(w);
            HWND from = reinterpret_cast<HWND>(l);
            const auto& c = self->colors_;
            const bool muted = from == self->report_label_ || from == self->position_;
            SetTextColor(dc, from == self->report_ && self->error_ ? (c ? c->error : error_ink)
                             : c                                    ? (muted ? c->muted : c->ink)
                                                                    : GetSysColor(COLOR_WINDOWTEXT));
            SetBkColor(dc, c ? c->back : GetSysColor(COLOR_WINDOW));
            return reinterpret_cast<LRESULT>(c ? self->back_brush_ : GetSysColorBrush(COLOR_WINDOW));
        }
        case WM_ERASEBKGND: {  // the pane's own background in the theme's colour
            RECT r{};
            GetClientRect(hwnd, &r);
            FillRect(reinterpret_cast<HDC>(w), &r, self->colors_ ? self->back_brush_ : GetSysColorBrush(COLOR_WINDOW));
            return 1;
        }
        case WM_SIZE:
            self->layout();
            return 0;
        default: break;
        }
    }
    return DefWindowProcW(hwnd, message, w, l);
}

}  // namespace NiTCAD::app
