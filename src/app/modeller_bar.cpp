// The modeller's controls on Win32 (modeller_bar.hpp; DECISIONS.md Unit 26m, MT1-MT3).
#include "modeller_bar.hpp"

#include <algorithm>
#include <array>

#include "win32.hpp"

namespace NiTCAD::app {

namespace {

using namespace modeller;

constexpr const wchar_t* bar_class = L"NiTCADModellerBar";
constexpr std::array<const wchar_t*, ModellerBar::corner_count> corner_texts{L"Camera free (L)", L"Reset view", L"Fit (F)", L"Measure (M)",
                                                                             L"Mesh"};

HWND make(const wchar_t* cls, const wchar_t* text, DWORD style, HWND parent, int id, DWORD ex = 0) {
    return CreateWindowExW(ex, cls, text, style | WS_CHILD | WS_VISIBLE, 0, 0, 10, 10, parent,
                           reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)), GetModuleHandleW(nullptr), nullptr);
}

int scaled(int at_96, int dpi) { return MulDiv(at_96, dpi, 96); }

}  // namespace

ModellerBar::ModellerBar(HWND strip, HWND panel, int first, bool three_d, WNDPROC proc, void* self)
    : strip_(strip), panel_(panel), first_(first), three_d_(three_d) {
    win32::register_class(bar_class, DefWindowProcW, GetSysColorBrush(COLOR_INFOBK));
    for (int g = 0; g < group_count; ++g) {
        groups_.push_back(make(L"BUTTON", win32::widen(group_name(static_cast<Group>(g))).c_str(),
                               WS_TABSTOP | BS_AUTORADIOBUTTON | BS_PUSHLIKE | (g == 0 ? WS_GROUP : 0), strip_, first_ + group0 + g));
    }
    for (int c = 0; c < corner_count; ++c) {
        const bool toggle = c == camera || c == measure || c == mesh;
        corner_.push_back(make(L"BUTTON", corner_texts[static_cast<std::size_t>(c)],
                               WS_TABSTOP | (c == 0 ? WS_GROUP : 0) | (toggle ? BS_AUTOCHECKBOX | BS_PUSHLIKE : BS_PUSHBUTTON), strip_,
                               first_ + corner0 + c));
    }
    bar_ = CreateWindowExW(WS_EX_CONTROLPARENT, bar_class, L"", WS_CHILD | WS_VISIBLE | WS_CLIPCHILDREN | WS_CLIPSIBLINGS | WS_BORDER, 0, 0,
                           10, 10, panel_, nullptr, GetModuleHandleW(nullptr), nullptr);
    SetWindowLongPtrW(bar_, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(self));
    SetWindowLongPtrW(bar_, GWLP_WNDPROC, reinterpret_cast<LONG_PTR>(proc));
    hint_ = make(L"STATIC", L"", SS_LEFT | SS_NOPREFIX, bar_, 0);
    chip_ = CreateWindowExW(0, L"STATIC", L"", WS_CHILD | WS_CLIPSIBLINGS | WS_BORDER | SS_CENTER | SS_NOPREFIX, 0, 0, 10, 10, panel_, nullptr,
                            GetModuleHandleW(nullptr), nullptr);
    set_group(Group::solid);
    raise();
}

void ModellerBar::set_font(HFONT font, int dpi) {
    font_ = font;
    dpi_ = dpi;
    for (const auto& list : {groups_, tools_, corner_}) {
        for (HWND h : list) SendMessageW(h, WM_SETFONT, reinterpret_cast<WPARAM>(font), TRUE);
    }
    for (HWND h : {hint_, chip_}) SendMessageW(h, WM_SETFONT, reinterpret_cast<WPARAM>(font), TRUE);
    size_bar();
}

void ModellerBar::layout_strip(int right, int pad, int height) {
    int x = pad;
    for (HWND h : groups_) {
        const int w = scaled(76, dpi_);
        MoveWindow(h, x, pad, w, height, TRUE);
        x += w + pad / 2;
    }
    const std::array<int, corner_count> widths{scaled(112, dpi_), scaled(80, dpi_), scaled(60, dpi_), scaled(88, dpi_), scaled(56, dpi_)};
    int total = 0;
    for (const int w : widths) total += w + pad / 2;
    int cx = std::max(x + 3 * pad, right - pad - total);
    for (int c = 0; c < corner_count; ++c) {
        MoveWindow(corner_[static_cast<std::size_t>(c)], cx, pad, widths[static_cast<std::size_t>(c)], height, TRUE);
        cx += widths[static_cast<std::size_t>(c)] + pad / 2;
    }
}

void ModellerBar::set_group(Group g) {
    group_ = g;
    for (int k = 0; k < group_count; ++k) {
        SendMessageW(groups_[static_cast<std::size_t>(k)], BM_SETCHECK, k == static_cast<int>(g) ? BST_CHECKED : BST_UNCHECKED, 0);
    }
    rebuild();
}

void ModellerBar::rebuild() {
    for (HWND h : tools_) DestroyWindow(h);
    tools_.clear();
    const auto list = tools_for(group_, three_d_);
    for (std::size_t k = 0; k < list.size(); ++k) {
        const bool action = list[k].action == Action::duplicate || list[k].action == Action::overlap;  // pressed once
        tools_.push_back(make(L"BUTTON", win32::widen(list[k].label).c_str(),
                              WS_TABSTOP | (k == 0 ? WS_GROUP : 0) | (action ? BS_PUSHBUTTON : BS_AUTOCHECKBOX | BS_PUSHLIKE), bar_,
                              first_ + tool0 + static_cast<int>(k)));
        if (font_) SendMessageW(tools_.back(), WM_SETFONT, reinterpret_cast<WPARAM>(font_), TRUE);
    }
    set_tool(tool_);
}

void ModellerBar::set_tool(int tool) {
    tool_ = tool;
    const auto list = tools_for(group_, three_d_);
    std::string how;
    std::string_view name;
    for (std::size_t k = 0; k < list.size(); ++k) {
        const bool armable = list[k].action == Action::tool || list[k].action == Action::mode;
        const bool armed = armable && (list[k].action == Action::tool ? list[k].tool : 100 + list[k].tool) == tool;
        if (armable) SendMessageW(tools_[k], BM_SETCHECK, armed ? BST_CHECKED : BST_UNCHECKED, 0);
        if (armed && name.empty()) {
            name = list[k].label;
            how = std::string(three_d_ ? list[k].how_3d : list[k].how_2d);
        }
    }
    if (name.empty()) {
        line_ = std::string(later(group_));
        if (list.empty()) line_ = std::string(group_name(group_)) + ": " + line_;
    } else {
        const auto paren = name.find(" (");
        line_ = std::string(name.substr(0, paren)) + ": " + how;
    }
    win32::set_text(hint_, line_);
    size_bar();
}

void ModellerBar::set_camera_fixed(bool fixed) {
    HWND b = corner_[camera];
    SendMessageW(b, BM_SETCHECK, fixed ? BST_CHECKED : BST_UNCHECKED, 0);
    SetWindowTextW(b, fixed ? L"Camera fixed (L)" : L"Camera free (L)");
}

void ModellerBar::set_mesh(bool on) { SendMessageW(corner_[mesh], BM_SETCHECK, on ? BST_CHECKED : BST_UNCHECKED, 0); }
void ModellerBar::set_measure(bool on) { SendMessageW(corner_[measure], BM_SETCHECK, on ? BST_CHECKED : BST_UNCHECKED, 0); }

void ModellerBar::size_bar() {
    const int pad = scaled(4, dpi_), bh = scaled(24, dpi_), line = scaled(15, dpi_);
    // The buttons narrow to keep the bar inside the viewport (a group of five at 200 %), never under 48 px at 96 dpi.
    int bw = scaled(84, dpi_);
    const int n = static_cast<int>(tools_.size());
    const int room = static_cast<int>(viewport_.width()) - 2 * scaled(8, dpi_) - 2 - pad;
    if (n > 0 && room > 0 && n * (bw + pad / 2) > room) bw = std::max(scaled(48, dpi_), room / n - pad / 2);
    // Too narrow even then (a small window at 200 %): the buttons wrap into rows.
    const int per_row = n == 0 ? 1 : std::clamp((room - pad) / (bw + pad / 2), 1, n);
    int right = pad;
    for (int k = 0; k < n; ++k) {
        const int x = pad + (k % per_row) * (bw + pad / 2);
        MoveWindow(tools_[static_cast<std::size_t>(k)], x, pad + (k / per_row) * (bh + pad / 2), bw, bh, TRUE);
        right = std::max(right, x + bw + pad / 2);
    }
    const int rows = n == 0 ? 0 : (n + per_row - 1) / per_row;
    // The hint is 300 px wide at 96 dpi, narrowed to the viewport's room (wrapping onto more lines).
    const int w = std::max(right + pad, room > 0 ? std::min(scaled(300, dpi_), room + pad) : scaled(300, dpi_));
    const int top = pad + rows * (bh + pad / 2);
    int lines = 2;
    if (const HDC dc = GetDC(hint_)) {
        const HGDIOBJ old = font_ ? SelectObject(dc, font_) : nullptr;
        const std::wstring text = win32::widen(line_);
        RECT t{0, 0, w - 2 * pad, 0}, one{0, 0, w - 2 * pad, 0};
        DrawTextW(dc, text.c_str(), static_cast<int>(text.size()), &t, DT_CALCRECT | DT_WORDBREAK | DT_NOPREFIX);
        DrawTextW(dc, L"M", 1, &one, DT_CALCRECT | DT_NOPREFIX);
        if (one.bottom > 0) lines = std::max(lines, static_cast<int>((t.bottom + one.bottom - 1) / one.bottom));
        if (old) SelectObject(dc, old);
        ReleaseDC(hint_, dc);
    }
    MoveWindow(hint_, pad, top, w - 2 * pad, lines * line, TRUE);
    const int h = top + lines * line + pad;
    RECT r{};
    GetWindowRect(bar_, &r);
    if (r.right - r.left != w + 2 || r.bottom - r.top != h + 2) SetWindowPos(bar_, nullptr, 0, 0, w + 2, h + 2, SWP_NOMOVE | SWP_NOZORDER | SWP_NOACTIVATE);
    place(viewport_, selection_);
}

void ModellerBar::place(const Box& viewport, const std::optional<Box>& selection) {
    const bool resize = viewport.width() != viewport_.width();
    viewport_ = viewport;
    selection_ = selection;
    if (resize) {  // the buttons' width follows the viewport's (size_bar places the bar again)
        size_bar();
        return;
    }
    RECT r{};
    GetWindowRect(bar_, &r);
    const double gap = scaled(8, dpi_);
    const Placement p = place_bar(viewport, selection, r.right - r.left, r.bottom - r.top, gap);
    if (!(p.bar == placed_.bar)) {
        SetWindowPos(bar_, HWND_TOP, static_cast<int>(p.bar.left), static_cast<int>(p.bar.top), 0, 0, SWP_NOSIZE | SWP_NOACTIVATE);
    }
    placed_ = p;
}

std::optional<ModellerBar::Command> ModellerBar::command(int id) const {
    const int k = id - first_;
    if (k >= group0 && k < group0 + group_count) return Command{Command::Kind::group, k - group0};
    if (k >= corner0 && k < corner0 + corner_count) return Command{Command::Kind::corner, k - corner0};
    const auto list = tools_for(group_, three_d_);
    if (k >= tool0 && k < tool0 + static_cast<int>(list.size())) {
        const ToolEntry& t = list[static_cast<std::size_t>(k - tool0)];
        if (t.action == Action::duplicate) return Command{Command::Kind::duplicate, 0};
        if (t.action == Action::overlap) return Command{Command::Kind::overlap, 0};
        if (t.action == Action::mode) return Command{Command::Kind::mode, t.tool};
        return Command{Command::Kind::tool, t.tool};
    }
    return std::nullopt;
}

void ModellerBar::chip(int x, int y, const std::string& text) {
    if (text == chip_text_ && (text.empty() || IsWindowVisible(chip_))) {
        if (!text.empty()) SetWindowPos(chip_, HWND_TOP, x, y, 0, 0, SWP_NOSIZE | SWP_NOACTIVATE);
        return;
    }
    chip_text_ = text;
    if (text.empty()) {
        ShowWindow(chip_, SW_HIDE);
        return;
    }
    const std::wstring w = win32::widen(text);
    SIZE s{};
    if (HDC dc = GetDC(chip_)) {
        HGDIOBJ old = font_ ? SelectObject(dc, font_) : nullptr;
        GetTextExtentPoint32W(dc, w.c_str(), static_cast<int>(w.size()), &s);
        if (old) SelectObject(dc, old);
        ReleaseDC(chip_, dc);
    }
    SetWindowTextW(chip_, w.c_str());
    SetWindowPos(chip_, HWND_TOP, x, y, s.cx + scaled(12, dpi_), s.cy + scaled(6, dpi_), SWP_NOACTIVATE | SWP_SHOWWINDOW);
}

void ModellerBar::raise() {
    SetWindowPos(bar_, HWND_TOP, 0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE | SWP_NOACTIVATE);
    SetWindowPos(chip_, HWND_TOP, 0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE | SWP_NOACTIVATE);
}

}  // namespace NiTCAD::app
