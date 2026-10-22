// The start, run monitor and results pages (pages.hpp).
#include "pages.hpp"

#include <commctrl.h>

#include <algorithm>
#include <format>
#include <variant>

#include "NiTCAD/app/commands.hpp"
#include "NiTCAD/app/layout.hpp"
#include "NiTCAD/app/run_monitor.hpp"
#include "win32.hpp"

namespace NiTCAD::app {

namespace {

constexpr const wchar_t* page_class = L"NiTCADPage";

HWND make(const wchar_t* cls, const wchar_t* text, DWORD style, HWND parent, int id, DWORD ex = 0) {
    return CreateWindowExW(ex, cls, text, style | WS_CHILD | WS_VISIBLE, 0, 0, 10, 10, parent,
                           reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)), GetModuleHandleW(nullptr), nullptr);
}

HWND make_panel(HWND parent, int id, void* self, WNDPROC proc) {
    win32::register_class(page_class, DefWindowProcW, GetSysColorBrush(COLOR_WINDOW));
    HWND p = CreateWindowExW(WS_EX_CONTROLPARENT, page_class, L"", WS_CHILD | WS_CLIPCHILDREN, 0, 0, 10, 10, parent,
                             reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)), GetModuleHandleW(nullptr), nullptr);
    SetWindowLongPtrW(p, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(self));
    SetWindowLongPtrW(p, GWLP_WNDPROC, reinterpret_cast<LONG_PTR>(proc));
    return p;
}

void add_column(HWND list, int index, const wchar_t* text, int width) {
    LVCOLUMNW c{};
    c.mask = LVCF_TEXT | LVCF_WIDTH;
    c.pszText = const_cast<wchar_t*>(text);
    c.cx = width;
    SendMessageW(list, LVM_INSERTCOLUMNW, static_cast<WPARAM>(index), reinterpret_cast<LPARAM>(&c));
}

void add_row(HWND list, const std::vector<std::string>& cells) {
    LVITEMW item{};
    item.mask = LVIF_TEXT;
    item.iItem = static_cast<int>(SendMessageW(list, LVM_GETITEMCOUNT, 0, 0));
    std::wstring first = win32::widen(cells.empty() ? std::string() : cells[0]);
    item.pszText = first.data();
    const int row = static_cast<int>(SendMessageW(list, LVM_INSERTITEMW, 0, reinterpret_cast<LPARAM>(&item)));
    for (std::size_t c = 1; c < cells.size(); ++c) {
        std::wstring t = win32::widen(cells[c]);
        LVITEMW s{};
        s.iSubItem = static_cast<int>(c);
        s.pszText = t.data();
        SendMessageW(list, LVM_SETITEMTEXTW, static_cast<WPARAM>(row), reinterpret_cast<LPARAM>(&s));
    }
}

std::wstring list_cell(HWND list, int row, int column) {
    wchar_t buffer[512]{};
    LVITEMW item{};
    item.iSubItem = column;
    item.pszText = buffer;
    item.cchTextMax = 512;
    SendMessageW(list, LVM_GETITEMTEXTW, static_cast<WPARAM>(row), reinterpret_cast<LPARAM>(&item));
    return buffer;
}

void set_font(HWND parent, HFONT font) {
    for (HWND c = GetWindow(parent, GW_CHILD); c; c = GetWindow(c, GW_HWNDNEXT)) {
        SendMessageW(c, WM_SETFONT, reinterpret_cast<WPARAM>(font), TRUE);
    }
}

}  // namespace

// --- the start page --------------------------------------------------------------------------------

namespace {
constexpr int id_recent = 10, id_examples = 11;
}

StartPage::StartPage(HWND parent, Host& host, int id) : host_(host) {
    panel_ = make_panel(parent, id, this, proc);
    for (const wchar_t* h : {L"Start", L"Recent", L"Examples"}) headings_.push_back(make(L"STATIC", h, 0, panel_, 0));
    for (const auto& [command, text] : {std::pair{Command::new_deck, L"&New deck"}, std::pair{Command::open_deck, L"&Open deck..."},
                                        std::pair{Command::open_result, L"Open &result..."}}) {
        buttons_.push_back(make(L"BUTTON", text, WS_TABSTOP | BS_PUSHBUTTON, panel_, static_cast<int>(command)));
    }
    recent_ = make(L"LISTBOX", L"", WS_TABSTOP | WS_VSCROLL | LBS_NOTIFY | LBS_NOINTEGRALHEIGHT, panel_, id_recent,
                   WS_EX_CLIENTEDGE);
    examples_ = make(WC_LISTVIEWW, L"", WS_TABSTOP | LVS_REPORT | LVS_SINGLESEL | LVS_SHOWSELALWAYS, panel_, id_examples,
                     WS_EX_CLIENTEDGE);
    SendMessageW(examples_, LVM_SETEXTENDEDLISTVIEWSTYLE, 0, LVS_EX_FULLROWSELECT);
    add_column(examples_, 0, L"Example", 180);
    add_column(examples_, 1, L"What it shows", 520);
}

StartPage::~StartPage() {
    if (panel_) {
        SetWindowLongPtrW(panel_, GWLP_USERDATA, 0);
        DestroyWindow(panel_);
    }
}

void StartPage::set(const std::vector<std::filesystem::path>& recent, const std::vector<Example>& examples) {
    SendMessageW(recent_, LB_RESETCONTENT, 0, 0);
    for (const auto& p : recent) SendMessageW(recent_, LB_ADDSTRING, 0, reinterpret_cast<LPARAM>(p.wstring().c_str()));
    SendMessageW(examples_, LVM_DELETEALLITEMS, 0, 0);
    for (const Example& e : examples) add_row(examples_, {e.title, e.description});
}

void StartPage::set_dpi(int dpi, HFONT font, HFONT bold) {
    dpi_ = dpi;
    set_font(panel_, font);
    for (HWND h : headings_) SendMessageW(h, WM_SETFONT, reinterpret_cast<WPARAM>(bold), TRUE);
    SendMessageW(examples_, LVM_SETCOLUMNWIDTH, 0, scaled(180, dpi));
    SendMessageW(examples_, LVM_SETCOLUMNWIDTH, 1, scaled(520, dpi));
    layout();
}

int StartPage::example_count() const { return static_cast<int>(SendMessageW(examples_, LVM_GETITEMCOUNT, 0, 0)); }
std::wstring StartPage::example_title(int index) const { return list_cell(examples_, index, 0); }
int StartPage::recent_count() const { return static_cast<int>(SendMessageW(recent_, LB_GETCOUNT, 0, 0)); }

void StartPage::layout() {
    RECT r{};
    GetClientRect(panel_, &r);
    const int pad = scaled(16, dpi_), row = scaled(26, dpi_), gap = scaled(6, dpi_), left_w = scaled(260, dpi_);
    int y = pad;
    MoveWindow(headings_[0], pad, y, left_w, row, TRUE);
    y += row;
    for (HWND b : buttons_) {
        MoveWindow(b, pad, y, scaled(160, dpi_), row, TRUE);
        y += row + gap;
    }
    y += gap;
    MoveWindow(headings_[1], pad, y, left_w, row, TRUE);
    y += row;
    MoveWindow(recent_, pad, y, left_w, std::max<int>(r.bottom - y - pad, 0), TRUE);
    const int x = pad + left_w + 2 * pad;
    MoveWindow(headings_[2], x, pad, std::max<int>(r.right - x - pad, 0), row, TRUE);
    MoveWindow(examples_, x, pad + row, std::max<int>(r.right - x - pad, 0), std::max<int>(r.bottom - 2 * pad - row, 0), TRUE);
    SendMessageW(examples_, LVM_SETCOLUMNWIDTH, 1, LVSCW_AUTOSIZE_USEHEADER);  // the description fills the rest
}

LRESULT CALLBACK StartPage::proc(HWND hwnd, UINT message, WPARAM w, LPARAM l) {
    auto* self = reinterpret_cast<StartPage*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
    if (self) {
        switch (message) {
        case WM_SIZE:
            self->layout();
            return 0;
        case WM_COMMAND: {
            const int id = LOWORD(w), code = HIWORD(w);
            if (id == id_recent && code == LBN_DBLCLK) {
                const auto k = SendMessageW(self->recent_, LB_GETCURSEL, 0, 0);
                if (k >= 0) self->host_.start_recent(static_cast<std::size_t>(k));
            } else if (code == BN_CLICKED && id >= 100) {
                self->host_.start_command(id);
            }
            return 0;
        }
        case WM_NOTIFY: {
            const auto* n = reinterpret_cast<const NMHDR*>(l);
            if (n->hwndFrom == self->examples_ && (n->code == NM_DBLCLK || n->code == LVN_ITEMACTIVATE)) {
                const int k = static_cast<int>(SendMessageW(self->examples_, LVM_GETNEXTITEM, static_cast<WPARAM>(-1), LVNI_SELECTED));
                if (k >= 0) self->host_.start_example(static_cast<std::size_t>(k));
                return 0;
            }
            break;
        }
        case WM_VKEYTOITEM:
            if (reinterpret_cast<HWND>(l) == self->recent_ && LOWORD(w) == VK_RETURN) {
                const auto k = SendMessageW(self->recent_, LB_GETCURSEL, 0, 0);
                if (k >= 0) self->host_.start_recent(static_cast<std::size_t>(k));
                return -2;
            }
            return -1;
        case WM_CTLCOLORSTATIC:
            SetBkMode(reinterpret_cast<HDC>(w), TRANSPARENT);
            return reinterpret_cast<LRESULT>(GetSysColorBrush(COLOR_WINDOW));
        default:
            break;
        }
    }
    return DefWindowProcW(hwnd, message, w, l);
}

// --- the run monitor -------------------------------------------------------------------------------

namespace {
constexpr int id_rows = 10, id_stop = 11;
}

RunPage::RunPage(HWND parent, Host& host, int id) : host_(host) {
    panel_ = make_panel(parent, id, this, proc);
    label_ = make(L"STATIC", L"No run yet.", SS_LEFTNOWORDWRAP | SS_CENTERIMAGE, panel_, 0);
    progress_ = make(PROGRESS_CLASSW, L"", 0, panel_, 0);
    SendMessageW(progress_, PBM_SETRANGE32, 0, 1000);
    stop_ = make(L"BUTTON", L"S&top", WS_TABSTOP | BS_PUSHBUTTON, panel_, id_stop);
    list_ = make(WC_LISTVIEWW, L"", WS_TABSTOP | LVS_REPORT | LVS_OWNERDATA | LVS_SHOWSELALWAYS, panel_, id_rows, WS_EX_CLIENTEDGE);
    SendMessageW(list_, LVM_SETEXTENDEDLISTVIEWSTYLE, 0, LVS_EX_FULLROWSELECT | LVS_EX_DOUBLEBUFFER);
    const auto columns = monitor_columns();
    for (std::size_t c = 0; c < columns.size(); ++c) {
        add_column(list_, static_cast<int>(c), win32::widen(columns[c]).c_str(), c == 2 ? 150 : 80);
    }
    plot_ = std::make_unique<Viewport>(panel_, *this, 12);
    plot_->show_message("The Newton iterations of a run are plotted here.");
}

RunPage::~RunPage() {
    release();
    plot_.reset();
    if (panel_) {
        SetWindowLongPtrW(panel_, GWLP_USERDATA, 0);
        DestroyWindow(panel_);
    }
}

void RunPage::release() {
    if (plot_) plot_->release();
}

void RunPage::update(const RunController& runner, render::Renderer* renderer, bool force) {
    rows_ = runner.points();
    SendMessageW(list_, LVM_SETITEMCOUNT, static_cast<WPARAM>(rows_.size()), LVSICF_NOSCROLL | LVSICF_NOINVALIDATEALL);
    InvalidateRect(list_, nullptr, FALSE);
    if (!rows_.empty()) SendMessageW(list_, LVM_ENSUREVISIBLE, static_cast<WPARAM>(rows_.size() - 1), FALSE);
    const double f = progress_fraction(rows_, runner.planned_points());
    const bool marquee = f < 0.0 && runner.busy();
    LONG_PTR style = GetWindowLongPtrW(progress_, GWL_STYLE);
    if (((style & PBS_MARQUEE) != 0) != marquee) {
        SetWindowLongPtrW(progress_, GWL_STYLE, marquee ? style | PBS_MARQUEE : style & ~static_cast<LONG_PTR>(PBS_MARQUEE));
        SendMessageW(progress_, PBM_SETMARQUEE, marquee, 30);
    }
    if (!marquee) SendMessageW(progress_, PBM_SETPOS, static_cast<WPARAM>(1000.0 * std::max(f, runner.state() == RunState::done ? 1.0 : 0.0)), 0);
    EnableWindow(stop_, runner.state() == RunState::running);
    std::string text = std::format("{}: {} rows, {} Newton iterations, {:.1f} s", to_string(runner.state()), rows_.size(),
                                   runner.iterations().size(), runner.elapsed_seconds());
    if (runner.planned_points() > 0) text += std::format(", {} points planned", runner.planned_points());
    if (const auto stopped = runner.stopped_early()) text += " - stopped early: " + stopped->message;
    if (runner.error()) text += " - " + runner.error()->message;
    win32::set_text(label_, text);
    // The plot is drawn on the UI thread: at most every quarter second, and never more than a twentieth of
    // the time (a software adapter takes tens of milliseconds a draw), so the window keeps running (WG1).
    const auto now = std::chrono::steady_clock::now();
    const auto interval = std::chrono::duration<double, std::milli>(std::max(250.0, 20.0 * plot_->last_present_ms()));
    if (!force && now - plotted_ < interval) return;
    plotted_ = now;
    auto scene = convergence_plot(runner.iterations());
    if (!scene || !renderer) return;
    scene_ = std::move(*scene);
    plot_->show(renderer, &*scene_);
}

void RunPage::set_dpi(int dpi, HFONT font) {
    dpi_ = dpi;
    set_font(panel_, font);
    const auto columns = monitor_columns();
    for (std::size_t c = 0; c < columns.size(); ++c) {
        SendMessageW(list_, LVM_SETCOLUMNWIDTH, c, scaled(c == 2 ? 150 : 80, dpi));
    }
    layout();
}

int RunPage::row_count() const { return static_cast<int>(rows_.size()); }

std::wstring RunPage::cell(int row, int column) const {
    if (row < 0 || row >= static_cast<int>(rows_.size())) return {};
    const auto cells = monitor_cells(rows_[static_cast<std::size_t>(row)], static_cast<std::size_t>(row));
    return column >= 0 && column < static_cast<int>(cells.size()) ? win32::widen(cells[static_cast<std::size_t>(column)])
                                                                    : std::wstring();
}

void RunPage::layout() {
    RECT r{};
    GetClientRect(panel_, &r);
    const int pad = scaled(8, dpi_), row = scaled(26, dpi_), stop_w = scaled(90, dpi_);
    const int list_w = std::max<int>((r.right - 3 * pad) * 55 / 100, 0);
    MoveWindow(stop_, pad, pad, stop_w, row, TRUE);
    MoveWindow(progress_, 2 * pad + stop_w, pad + row / 4, std::max<int>(scaled(260, dpi_), 0), row / 2, TRUE);
    MoveWindow(label_, 3 * pad + stop_w + scaled(260, dpi_), pad, std::max<int>(r.right - 4 * pad - stop_w - scaled(260, dpi_), 0),
               row, TRUE);
    const int top = 2 * pad + row;
    MoveWindow(list_, pad, top, list_w, std::max<int>(r.bottom - top - pad, 0), TRUE);
    MoveWindow(plot_->hwnd(), 2 * pad + list_w, top, std::max<int>(r.right - 3 * pad - list_w, 0),
               std::max<int>(r.bottom - top - pad, 0), TRUE);
}

LRESULT CALLBACK RunPage::proc(HWND hwnd, UINT message, WPARAM w, LPARAM l) {
    auto* self = reinterpret_cast<RunPage*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
    if (self) {
        switch (message) {
        case WM_SIZE:
            self->layout();
            return 0;
        case WM_COMMAND:
            if (LOWORD(w) == id_stop && HIWORD(w) == BN_CLICKED) self->host_.run_page_stop();
            return 0;
        case WM_NOTIFY: {
            auto* n = reinterpret_cast<NMHDR*>(l);
            if (n->hwndFrom == self->list_ && n->code == LVN_GETDISPINFOW) {
                auto* d = reinterpret_cast<NMLVDISPINFOW*>(l);
                if ((d->item.mask & LVIF_TEXT) && d->item.iItem < static_cast<int>(self->rows_.size())) {
                    const std::wstring t = self->cell(d->item.iItem, d->item.iSubItem);
                    wcsncpy_s(d->item.pszText, static_cast<std::size_t>(d->item.cchTextMax), t.c_str(), _TRUNCATE);
                }
                return 0;
            }
            break;
        }
        case WM_CTLCOLORSTATIC:
            SetBkMode(reinterpret_cast<HDC>(w), TRANSPARENT);
            return reinterpret_cast<LRESULT>(GetSysColorBrush(COLOR_WINDOW));
        default:
            break;
        }
    }
    return DefWindowProcW(hwnd, message, w, l);
}

// --- results ---------------------------------------------------------------------------------------

ResultsPage::ResultsPage(HWND parent, int id) {
    panel_ = make_panel(parent, id, this, proc);
    figures_label_ = make(L"STATIC", L"Extracted figures (Unit 24)", 0, panel_, 0);
    record_label_ = make(L"STATIC", L"Run record and provenance", 0, panel_, 0);
    figures_ = make(WC_LISTVIEWW, L"", WS_TABSTOP | LVS_REPORT | LVS_SHOWSELALWAYS, panel_, 10, WS_EX_CLIENTEDGE);
    record_ = make(WC_LISTVIEWW, L"", WS_TABSTOP | LVS_REPORT | LVS_SHOWSELALWAYS, panel_, 11, WS_EX_CLIENTEDGE);
    for (HWND h : {figures_, record_}) SendMessageW(h, LVM_SETEXTENDEDLISTVIEWSTYLE, 0, LVS_EX_FULLROWSELECT);
    add_column(figures_, 0, L"Figure", 180);
    add_column(figures_, 1, L"Value", 110);
    add_column(figures_, 2, L"Unit", 80);
    add_column(figures_, 3, L"Curve", 160);
    add_column(figures_, 4, L"Notes", 160);
    add_column(record_, 0, L"Item", 200);
    add_column(record_, 1, L"Value", 420);
}

ResultsPage::~ResultsPage() {
    if (panel_) {
        SetWindowLongPtrW(panel_, GWLP_USERDATA, 0);
        DestroyWindow(panel_);
    }
}

void ResultsPage::show(const Workspace* w) {
    SendMessageW(figures_, LVM_DELETEALLITEMS, 0, 0);
    SendMessageW(record_, LVM_DELETEALLITEMS, 0, 0);
    if (!w) return;
    for (const Figure& f : extracted_figures(*w)) {
        std::string notes;
        if (f.extrapolated) notes += "extrapolated";
        if (f.below_resolution) notes += std::string(notes.empty() ? "" : ", ") + "below resolution";
        add_row(figures_, {f.name, std::format("{:.5g}", f.value), f.unit, f.curve, notes});
    }
    for (const Line& l : run_record(*w)) add_row(record_, {l.name, l.value});
}

void ResultsPage::set_dpi(int dpi, HFONT font, HFONT bold) {
    dpi_ = dpi;
    set_font(panel_, font);
    for (HWND h : {figures_label_, record_label_}) SendMessageW(h, WM_SETFONT, reinterpret_cast<WPARAM>(bold), TRUE);
    layout();
}

int ResultsPage::figure_count() const { return static_cast<int>(SendMessageW(figures_, LVM_GETITEMCOUNT, 0, 0)); }
int ResultsPage::record_count() const { return static_cast<int>(SendMessageW(record_, LVM_GETITEMCOUNT, 0, 0)); }
std::wstring ResultsPage::figure_cell(int row, int column) const { return list_cell(figures_, row, column); }

void ResultsPage::layout() {
    RECT r{};
    GetClientRect(panel_, &r);
    const int pad = scaled(8, dpi_), row = scaled(24, dpi_);
    const int h = std::max<int>((r.bottom - 3 * pad - 2 * row) / 2, 0);
    MoveWindow(figures_label_, pad, pad, std::max<int>(r.right - 2 * pad, 0), row, TRUE);
    MoveWindow(figures_, pad, pad + row, std::max<int>(r.right - 2 * pad, 0), std::max(h * 2 / 3, 0), TRUE);
    const int y = pad + row + h * 2 / 3 + pad;
    MoveWindow(record_label_, pad, y, std::max<int>(r.right - 2 * pad, 0), row, TRUE);
    MoveWindow(record_, pad, y + row, std::max<int>(r.right - 2 * pad, 0), std::max<int>(r.bottom - y - row - pad, 0), TRUE);
}

LRESULT CALLBACK ResultsPage::proc(HWND hwnd, UINT message, WPARAM w, LPARAM l) {
    auto* self = reinterpret_cast<ResultsPage*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
    if (self && message == WM_SIZE) {
        self->layout();
        return 0;
    }
    if (self && message == WM_CTLCOLORSTATIC) {
        SetBkMode(reinterpret_cast<HDC>(w), TRANSPARENT);
        return reinterpret_cast<LRESULT>(GetSysColorBrush(COLOR_WINDOW));
    }
    return DefWindowProcW(hwnd, message, w, l);
}

}  // namespace NiTCAD::app
