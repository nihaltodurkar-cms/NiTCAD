// The main window (main_window.hpp).
#include "main_window.hpp"

#include <commctrl.h>
#include <span>
#include <richedit.h>  // Unit 26j: the console
#include <dwmapi.h>     // ZD6: the title bar's theme
#include <shobjidl.h>   // ZD5: ITaskbarList3
#include <windowsx.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <charconv>
#include <format>
#include <fstream>
#include <tuple>
#include <utility>
#include <variant>

#include "NiTCAD/app/compare.hpp"
#include "NiTCAD/app/run_monitor.hpp"
#include "NiTCAD/render/compare.hpp"
#include "NiTCAD/render/quality.hpp"
#include "NiTCAD/render/scene.hpp"
#include "win32.hpp"

namespace NiTCAD::app {

namespace {

namespace fs = std::filesystem;

constexpr const wchar_t* window_class = L"NiTCADMainWindow";
constexpr const wchar_t* panel_class = L"NiTCADInspector";

// Status bar parts.
constexpr int status_state = 0, status_file = 1, status_adapter = 2, status_cursor = 3;

// The tree's item data: a kind in the high word, an index in the low.
enum TreeKind : int { tree_none = 0, tree_field = 1, tree_view = 2, tree_region = 3 };
LPARAM tree_data(TreeKind k, int index) { return static_cast<LPARAM>((k << 16) | (index & 0xFFFF)); }

constexpr std::array<std::string_view, 4> scale_names{"(default)", "linear", "log", "signed-log"};
// ZD8: the app's default (viridis) first.
constexpr std::array<std::string_view, 5> colormap_names{"viridis", "turbo", "cividis", "coolwarm", "soft"};
// Unit 26c (YD1): the vector fields arrows and streamlines follow ("none" first).
constexpr std::array<std::string_view, 5> vector_names{"none", "E", "J", "Jn", "Jp"};
// The labels of the inspector's groups (Unit 26c), after the nine of Unit 26: their ids follow
// id_first_label in the order made.
constexpr int label_arrows = MainWindow::id_first_label + 9, label_lines = label_arrows + 1, label_from = label_arrows + 2,
              label_axis = label_arrows + 3, label_swept = label_arrows + 4, label_sensed = label_arrows + 5,
              label_frequency = label_arrows + 6;
// Unit 26c (YD5, YD6): the snapshot slider's label and the clipping box's per axis.
constexpr int label_time = label_arrows + 7, label_clip_x = label_arrows + 8;
// Unit 26d (ZD3): the colour range's label, after the three of the clipping box's sliders.
constexpr int label_range = label_clip_x + 3;

std::string format_um(double v) { return std::format("{:.3g}", v); }

// The numbers of a comma-separated list; empty when any is not a number.
std::vector<double> numbers(std::string_view text) {
    std::vector<double> out;
    while (true) {
        const auto comma = text.find(',');
        std::string_view item = text.substr(0, comma);
        while (!item.empty() && item.front() == ' ') item.remove_prefix(1);
        while (!item.empty() && item.back() == ' ') item.remove_suffix(1);
        double v = 0.0;
        const auto [p, ec] = std::from_chars(item.data(), item.data() + item.size(), v);
        if (item.empty() || ec != std::errc{} || p != item.data() + item.size()) return {};
        out.push_back(v);
        if (comma == std::string_view::npos) return out;
        text = text.substr(comma + 1);
    }
}

std::string number(double v) {
    if (!std::isfinite(v)) return "-";
    return std::format("{:.4g}", v);
}

// The views a tab offers: the device tab the field views (and a 1D profile), the curves tab the plots (in 2D and
// 3D the profile along a cut line too).
bool belongs(ViewKind k, int tab, int dimension) {
    return (tab == MainWindow::tab_curves) == (is_plot(k) && (k != ViewKind::profile || dimension > 1));
}

HWND make(const wchar_t* cls, const wchar_t* text, DWORD style, HWND parent, int id, DWORD ex = 0) {
    return CreateWindowExW(ex, cls, text, style | WS_CHILD, 0, 0, 10, 10, parent,
                           reinterpret_cast<HMENU>(static_cast<INT_PTR>(id)), GetModuleHandleW(nullptr), nullptr);
}

void combo_fill(HWND combo, const std::vector<std::string>& items, int selected) {
    SendMessageW(combo, CB_RESETCONTENT, 0, 0);
    for (const auto& s : items) SendMessageW(combo, CB_ADDSTRING, 0, reinterpret_cast<LPARAM>(win32::widen(s).c_str()));
    SendMessageW(combo, CB_SETCURSEL, static_cast<WPARAM>(selected), 0);
}

int combo_index(HWND combo) { return static_cast<int>(SendMessageW(combo, CB_GETCURSEL, 0, 0)); }

std::wstring menu_text(const CommandInfo& c) {
    std::wstring t = win32::widen(c.text);
    if (!c.accelerator_text.empty()) t += L"\t" + win32::widen(c.accelerator_text);
    return t;
}

}  // namespace

// --- creation -------------------------------------------------------------------------------------

MainWindow::MainWindow(const WindowOptions& options)
    : options_(options),
      runner_([this] { PostMessageW(hwnd_, run_message, 0, 0); }),
      live_([this] { PostMessageW(hwnd_, live_message, 0, 0); }) {
    animate_ = options.animate;
    // The taskbar list (ZD5) and the shell dialogs are COM objects of this thread's apartment. The application
    // initializes COM first (this call is then a counted no-op); a host that did not, such as a test, would
    // otherwise get them from another thread's implicit apartment and crash releasing them after it ends.
    com_ = SUCCEEDED(CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED));
}

std::unique_ptr<MainWindow> MainWindow::create(const WindowOptions& options) {
    INITCOMMONCONTROLSEX icc{sizeof(icc), ICC_BAR_CLASSES | ICC_TREEVIEW_CLASSES | ICC_TAB_CLASSES | ICC_STANDARD_CLASSES |
                                              ICC_LISTVIEW_CLASSES | ICC_PROGRESS_CLASS};
    InitCommonControlsEx(&icc);
    win32::register_class(window_class, proc, GetSysColorBrush(COLOR_BTNFACE));
    win32::register_class(panel_class, panel_proc, GetSysColorBrush(COLOR_WINDOW));
    std::unique_ptr<MainWindow> w(new MainWindow(options));
    if (!options.state_file.empty()) w->state_ = load_state(options.state_file);
    const Layout& l = w->state_.layout;
    const bool placed = l.width > 0 && l.height > 0;
    HWND h = CreateWindowExW(WS_EX_CONTROLPARENT, window_class, L"NiTCAD", WS_OVERLAPPEDWINDOW | WS_CLIPCHILDREN,
                             placed ? l.x : CW_USEDEFAULT, placed ? l.y : CW_USEDEFAULT, placed ? l.width : 1280,
                             placed ? l.height : 800, nullptr, nullptr, GetModuleHandleW(nullptr), w.get());
    if (!h || !w->viewport_) return nullptr;
    if (options.show) {
        ShowWindow(h, !options.input ? SW_SHOWNOACTIVATE : placed && l.maximized ? SW_SHOWMAXIMIZED : SW_SHOWNORMAL);
        UpdateWindow(h);
    }
    if (!options.input) EnableWindow(h, FALSE);  // a disabled window still paints and runs its timers
    return w;
}

MainWindow::~MainWindow() {
    runner_.request_stop();
    runner_.wait();
    if (viewport_) viewport_->release();  // the swap chains before the renderer
    if (run_page_) run_page_->release();
    if (hwnd_) {
        SetWindowLongPtrW(hwnd_, GWLP_USERDATA, 0);
        DestroyWindow(hwnd_);
    }
    viewport_.reset();
    deck_.reset();
    renderer_.reset();
    if (taskbar_list_) taskbar_list_->Release();
    if (accelerators_) DestroyAcceleratorTable(accelerators_);
    for (HIMAGELIST l : {toolbar_images_, toolbar_disabled_, tree_images_}) {
        if (l) ImageList_Destroy(l);
    }
    if (font_) DeleteObject(font_);
    if (bold_) DeleteObject(bold_);
    if (com_) CoUninitialize();  // after the taskbar list's release
}

bool MainWindow::create_children() {
    toolbar_ = make(TOOLBARCLASSNAME, L"", WS_VISIBLE | TBSTYLE_FLAT | TBSTYLE_LIST | CCS_NODIVIDER | CCS_NORESIZE |
                                               TBSTYLE_TOOLTIPS,
                    hwnd_, id_toolbar);
    SendMessageW(toolbar_, TB_BUTTONSTRUCTSIZE, sizeof(TBBUTTON), 0);
    SendMessageW(toolbar_, TB_SETEXTENDEDSTYLE, 0, TBSTYLE_EX_MIXEDBUTTONS);
    struct Button {
        Command id;
        const wchar_t* text;
        bool check;
    };
    const Button buttons[] = {{Command::new_deck, L"New", false},      {Command::open_deck, L"Open", false},
                              {Command::save, L"Save", false},         {Command{}, nullptr, false},
                              {Command::check, L"Check", false},       {Command::run, L"Run", false},
                              {Command::stop, L"Stop", false},         {Command{}, nullptr, false},
                              {Command::fit, L"Fit", false},           {Command::probe_mode, L"Probe", true},
                              {Command::measure_mode, L"Measure", true}, {Command::cut_line_mode, L"Cut line", true},
                              {Command::perspective, L"Perspective", true}};
    int image = 0;  // look A: each button's icon, in order (make_icons)
    for (const Button& b : buttons) {
        TBBUTTON t{};
        t.iBitmap = b.text ? image++ : I_IMAGENONE;
        if (!b.text) {
            t.fsStyle = BTNS_SEP;
        } else {
            t.idCommand = static_cast<int>(b.id);
            t.fsState = TBSTATE_ENABLED;
            t.fsStyle = static_cast<BYTE>(BTNS_AUTOSIZE | BTNS_SHOWTEXT | (b.check ? BTNS_CHECK : 0));
            t.iString = reinterpret_cast<INT_PTR>(b.text);
        }
        SendMessageW(toolbar_, TB_ADDBUTTONSW, 1, reinterpret_cast<LPARAM>(&t));
    }
    status_ = make(STATUSCLASSNAMEW, L"", WS_VISIBLE | SBARS_SIZEGRIP, hwnd_, id_status);
    // Look A: a run's progress at the right of the status bar's first part, shown while it runs.
    status_progress_ = CreateWindowExW(0, PROGRESS_CLASSW, L"", WS_CHILD | PBS_SMOOTH, 0, 0, 10, 10, status_, nullptr,
                                       GetModuleHandleW(nullptr), nullptr);
    SendMessageW(status_progress_, PBM_SETRANGE32, 0, 1000);
    tree_ = make(WC_TREEVIEWW, L"", WS_VISIBLE | WS_TABSTOP | TVS_HASLINES | TVS_HASBUTTONS | TVS_LINESATROOT |
                                        TVS_SHOWSELALWAYS,
                 hwnd_, id_tree, WS_EX_CLIENTEDGE);
    // Check boxes for the regions (Unit 26c, YD4): set after creation, as the tree view asks.
    SetWindowLongPtrW(tree_, GWL_STYLE, GetWindowLongPtrW(tree_, GWL_STYLE) | TVS_CHECKBOXES);
    tabs_ = make(WC_TABCONTROLW, L"", WS_VISIBLE | WS_TABSTOP | WS_CLIPSIBLINGS, hwnd_, id_tabs);
    for (const wchar_t* t : {L"Start", L"Device view", L"Curves", L"Deck", L"Run monitor", L"Results"}) {
        TCITEMW item{};
        item.mask = TCIF_TEXT;
        item.pszText = const_cast<wchar_t*>(t);
        SendMessageW(tabs_, TCM_INSERTITEMW, static_cast<WPARAM>(SendMessageW(tabs_, TCM_GETITEMCOUNT, 0, 0)),
                     reinterpret_cast<LPARAM>(&item));
    }
    viewport_ = std::make_unique<Viewport>(hwnd_, *this, id_viewport);
    deck_ = std::make_unique<DeckPane>(hwnd_, *this, id_deck);
    start_ = std::make_unique<StartPage>(hwnd_, *this, id_start);
    run_page_ = std::make_unique<RunPage>(hwnd_, *this, id_run_page);
    results_ = std::make_unique<ResultsPage>(hwnd_, id_results);
    results_->on_pass = [this](std::size_t k) { results_show_pass(k); };
    metrics_ = make(L"STATIC", L"", WS_VISIBLE | SS_LEFTNOWORDWRAP | SS_CENTERIMAGE, hwnd_, id_metrics);
    // Unit 26j: the console along the window's bottom, a classic output window: the tab's text (plain, black on
    // white) or the points table over a tab strip (Output, Points, Messages); right-click for Copy, Save, Clear
    // and Hide.
    console_panel_ = CreateWindowExW(WS_EX_CONTROLPARENT, panel_class, L"", WS_CHILD | WS_VISIBLE | WS_CLIPCHILDREN, 0, 0,
                                     10, 10, hwnd_, reinterpret_cast<HMENU>(static_cast<INT_PTR>(id_console)),
                                     GetModuleHandleW(nullptr), this);
    {
        HWND tabs = make(WC_TABCONTROLW, L"", WS_VISIBLE | WS_TABSTOP | WS_CLIPSIBLINGS | TCS_BOTTOM, console_panel_,
                         id_console_tabs);
        for (const wchar_t* t : {L"Output", L"Points", L"Messages"}) {
            TCITEMW item{};
            item.mask = TCIF_TEXT;
            item.pszText = const_cast<wchar_t*>(t);
            SendMessageW(tabs, TCM_INSERTITEMW, static_cast<WPARAM>(SendMessageW(tabs, TCM_GETITEMCOUNT, 0, 0)),
                         reinterpret_cast<LPARAM>(&item));
        }
    }
    console_points_ = make(WC_LISTVIEWW, L"", WS_TABSTOP | LVS_REPORT | LVS_OWNERDATA | LVS_SINGLESEL | LVS_SHOWSELALWAYS,
                           console_panel_, id_console_points, WS_EX_CLIENTEDGE);
    SendMessageW(console_points_, LVM_SETEXTENDEDLISTVIEWSTYLE, 0, LVS_EX_FULLROWSELECT | LVS_EX_DOUBLEBUFFER);
    LoadLibraryW(L"Msftedit.dll");
    console_ = make(MSFTEDIT_CLASS, L"",
                    WS_VISIBLE | WS_VSCROLL | WS_HSCROLL | ES_MULTILINE | ES_READONLY | ES_AUTOVSCROLL | ES_AUTOHSCROLL |
                        ES_NOHIDESEL,
                    console_panel_, id_console_text, WS_EX_CLIENTEDGE);
    SendMessageW(console_, EM_EXLIMITTEXT, 0, 0x7FFFFFFE);
    SendMessageW(console_, EM_SETTARGETDEVICE, 0, 1);  // no wrapping: the tables keep their columns
    inspector_ = CreateWindowExW(WS_EX_CONTROLPARENT, panel_class, L"", WS_CHILD | WS_VISIBLE | WS_CLIPCHILDREN, 0, 0,
                                 10, 10, hwnd_, reinterpret_cast<HMENU>(static_cast<INT_PTR>(id_inspector)),
                                 GetModuleHandleW(nullptr), this);
    // The inspector's controls, each after its label.
    int label = id_first_label;
    const auto labelled = [&](const wchar_t* text, const wchar_t* cls, DWORD style, int id) {
        make(L"STATIC", text, WS_VISIBLE | SS_CENTERIMAGE | SS_ENDELLIPSIS, inspector_, label++);
        make(cls, L"", WS_VISIBLE | WS_TABSTOP | style, inspector_, id, cls == std::wstring_view(L"EDIT") ? WS_EX_CLIENTEDGE : 0);
    };
    const DWORD drop = CBS_DROPDOWNLIST | WS_VSCROLL;
    labelled(L"View", L"COMBOBOX", drop, id_view);
    labelled(L"Field", L"COMBOBOX", drop, id_field);
    labelled(L"Point", L"COMBOBOX", drop, id_point);
    labelled(L"Scale", L"COMBOBOX", drop, id_scale);
    labelled(L"Colormap", L"COMBOBOX", drop, id_colormap);
    for (const auto& [id, text] : {std::pair{id_contours, L"&Contours"}, std::pair{id_edges, L"Mesh &edges"},
                                   std::pair{id_contacts, L"Co&ntacts"}, std::pair{id_opaque, L"&Opaque insulators"},
                                   std::pair{id_perspective, L"Pe&rspective"}}) {
        make(L"BUTTON", text, WS_VISIBLE | WS_TABSTOP | BS_AUTOCHECKBOX, inspector_, id);
    }
    make(L"BUTTON", L"Enhanced &3D look", WS_VISIBLE | WS_TABSTOP | BS_AUTOCHECKBOX, inspector_, id_enhanced);
    make(L"BUTTON", L"A&xis box", WS_VISIBLE | WS_TABSTOP | BS_AUTOCHECKBOX, inspector_, id_axis_box);  // Unit 26h (JD1)
    labelled(L"Slice plane", L"COMBOBOX", drop, id_slice_axis);
    labelled(L"Slice at (µm)", L"EDIT", ES_AUTOHSCROLL, id_slice_at);
    make(L"BUTTON", L"C&lipping box (µm)", WS_VISIBLE | WS_TABSTOP | BS_AUTOCHECKBOX, inspector_, id_clip);
    make(L"EDIT", L"", WS_VISIBLE | WS_TABSTOP | ES_AUTOHSCROLL, inspector_, id_clip_box, WS_EX_CLIENTEDGE);
    labelled(L"Iso levels", L"EDIT", ES_AUTOHSCROLL, id_levels);
    labelled(L"Cut line (µm)", L"EDIT", ES_AUTOHSCROLL, id_cut);
    // Unit 26c (YD1-YD3): a drop-down list and its box side by side, shown with the views they apply to.
    const auto paired = [&](const wchar_t* text, int list, int box) {
        labelled(text, L"COMBOBOX", drop, list);
        make(L"EDIT", L"", WS_VISIBLE | WS_TABSTOP | ES_AUTOHSCROLL, inspector_, box, WS_EX_CLIENTEDGE);
    };
    paired(L"Arrows", id_arrows, id_arrow_count);
    paired(L"Streamlines", id_lines, id_line_count);
    paired(L"Seeds (µm)", id_from, id_seeds);
    paired(L"Explode (µm)", id_axis, id_gap);
    make(L"BUTTON", L"Sho&w the field", WS_VISIBLE | WS_TABSTOP | BS_AUTOCHECKBOX, inspector_, id_exploded_field);
    labelled(L"Swept contact", L"COMBOBOX", drop, id_swept);
    labelled(L"Sensed contact", L"COMBOBOX", drop, id_sensed);
    labelled(L"Frequency", L"COMBOBOX", drop, id_frequency);
    make(L"BUTTON", L"Lo&g axis", WS_VISIBLE | WS_TABSTOP | BS_AUTOCHECKBOX, inspector_, id_log);
    // Unit 26c (YD5): the snapshot slider, its time and Play; (YD6) a low and a high slider per axis of the box.
    const DWORD bar = WS_VISIBLE | WS_TABSTOP | TBS_HORZ | TBS_NOTICKS;
    make(L"STATIC", L"Snapshot", WS_VISIBLE | SS_CENTERIMAGE, inspector_, label++);
    make(TRACKBAR_CLASSW, L"", bar, inspector_, id_time);
    make(L"STATIC", L"", WS_VISIBLE | SS_LEFT | SS_CENTERIMAGE, inspector_, id_time_text);
    make(L"BUTTON", L"Pla&y", WS_VISIBLE | WS_TABSTOP | BS_PUSHBUTTON, inspector_, id_play);
    for (const auto& [text, low, high] : {std::tuple{L"Clip x", id_clip_x0, id_clip_x1},
                                          std::tuple{L"Clip y", id_clip_y0, id_clip_y1},
                                          std::tuple{L"Clip z", id_clip_z0, id_clip_z1}}) {
        make(L"STATIC", text, WS_VISIBLE | SS_CENTERIMAGE, inspector_, label++);
        for (const int id : {low, high}) {
            make(TRACKBAR_CLASSW, L"", bar, inspector_, id);
            SendMessageW(child(id), TBM_SETRANGE, FALSE, MAKELPARAM(0, clip_steps));
            SendMessageW(child(id), TBM_SETPAGESIZE, 0, clip_steps / 20);
        }
    }
    // Unit 26d (ZD3): the colour range, its two boxes and Auto side by side, after the colormap in tab order.
    make(L"STATIC", L"Range", WS_VISIBLE | SS_CENTERIMAGE, inspector_, label++);
    make(L"EDIT", L"", WS_VISIBLE | WS_TABSTOP | ES_AUTOHSCROLL, inspector_, id_range_min, WS_EX_CLIENTEDGE);
    make(L"EDIT", L"", WS_VISIBLE | WS_TABSTOP | ES_AUTOHSCROLL, inspector_, id_range_max, WS_EX_CLIENTEDGE);
    make(L"BUTTON", L"&Auto", WS_VISIBLE | WS_TABSTOP | BS_AUTOCHECKBOX, inspector_, id_range_auto);
    HWND after = child(id_colormap);
    for (const int id : {label_range, static_cast<int>(id_range_min), static_cast<int>(id_range_max),
                         static_cast<int>(id_range_auto)}) {
        SetWindowPos(child(id), after, 0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE | SWP_NOACTIVATE);
        after = child(id);
    }
    make(L"STATIC", L"", WS_VISIBLE | SS_LEFT, inspector_, id_probe_text);
    // Look A: the sections' headers (drawn by draw_section_header; a click folds the section).
    for (int k = 0; k < section_count; ++k) make(L"BUTTON", L"", WS_VISIBLE | BS_OWNERDRAW, inspector_, id_section_first + k);
    combo_fill(child(id_scale), {scale_names.begin(), scale_names.end()}, 0);
    combo_fill(child(id_colormap), {colormap_names.begin(), colormap_names.end()}, 0);
    combo_fill(child(id_slice_axis), {"none", "x", "y", "z"}, 0);
    combo_fill(child(id_arrows), {vector_names.begin(), vector_names.end()}, 0);
    combo_fill(child(id_lines), {vector_names.begin(), vector_names.end()}, 0);
    combo_fill(child(id_axis), {"x", "y", "z"}, 1);
    return toolbar_ && status_ && tree_ && tabs_ && metrics_ && inspector_ && viewport_->hwnd() && deck_->hwnd() &&
           start_->hwnd() && run_page_->hwnd() && results_->hwnd();
}

HWND MainWindow::child(int id) const noexcept {
    if (id >= id_view) return GetDlgItem(inspector_, id);
    return GetDlgItem(hwnd_, id);
}

// --- menus and accelerators ------------------------------------------------------------------------

void MainWindow::build_menu() {
    HMENU bar = CreateMenu();
    const auto names = menu_bar();
    for (std::size_t i = 0; i < names.size(); ++i) {
        HMENU popup = CreatePopupMenu();
        AppendMenuW(bar, MF_POPUP, reinterpret_cast<UINT_PTR>(popup), win32::widen(names[i]).c_str());
        fill_menu(popup, static_cast<int>(i));
    }
    SetMenu(hwnd_, bar);
}

// A top-level menu's items; the File menu's recent files and the deck menus' sections and keys depend on
// the state, so they are made again each time the menu opens (WM_INITMENUPOPUP).
void MainWindow::fill_menu(HMENU menu, int index) {
    while (GetMenuItemCount(menu) > 0) DeleteMenu(menu, 0, MF_BYPOSITION);
    const std::string_view name = menu_bar()[static_cast<std::size_t>(index)];
    std::vector<std::string_view> submenus;
    const auto add_commands = [&](std::string_view top) {
        for (const CommandInfo& c : commands()) {
            if (c.menu != top) continue;
            if (c.id == Command::exit) {
                HMENU recent = CreatePopupMenu();
                for (std::size_t k = 0; k < state_.recent.size(); ++k) {
                    const std::wstring t = std::format(L"&{} {}", k + 1, state_.recent[k].wstring());
                    AppendMenuW(recent, MF_STRING, static_cast<UINT_PTR>(Command::recent_first) + k, t.c_str());
                }
                if (state_.recent.empty()) AppendMenuW(recent, MF_STRING | MF_GRAYED, 0, L"(none)");
                AppendMenuW(menu, MF_SEPARATOR, 0, nullptr);
                AppendMenuW(menu, MF_POPUP, reinterpret_cast<UINT_PTR>(recent), L"Recent &files");
                AppendMenuW(menu, MF_SEPARATOR, 0, nullptr);
            }
            if (c.id == Command::fit || c.id == Command::perspective || c.id == Command::view_front ||
                c.id == Command::next_tab ||
                c.id == Command::tab_start || c.id == Command::find_error) {
                AppendMenuW(menu, MF_SEPARATOR, 0, nullptr);
            }
            if (!c.submenu.empty()) {
                // Unit 26h: a submenu, at its first item's place, holding every item that names it.
                if (std::ranges::find(submenus, c.submenu) != submenus.end()) continue;
                submenus.push_back(c.submenu);
                HMENU sub = CreatePopupMenu();
                for (const CommandInfo& d : commands()) {
                    if (d.menu == top && d.submenu == c.submenu) {
                        AppendMenuW(sub, MF_STRING, static_cast<UINT_PTR>(d.id), menu_text(d).c_str());
                    }
                }
                AppendMenuW(menu, MF_POPUP, reinterpret_cast<UINT_PTR>(sub), win32::widen(c.submenu).c_str());
                continue;
            }
            AppendMenuW(menu, MF_STRING, static_cast<UINT_PTR>(c.id), menu_text(c).c_str());
        }
    };
    // The deck menus (VQ13): "Go to [kind]" per section, then the keys to insert per section.
    const auto add_sections = [&](DeckMenu m) {
        const std::string text = deck_ ? deck_->text() : document_.text();
        for (const SectionSpec& spec : menu_sections(m)) {
            const auto where = find_section(text, spec.kind);
            const std::size_t id = section_items_.size();
            section_items_.push_back({std::string(spec.kind)});
            const std::wstring label = L"Go to " + win32::widen(spec.title) + (where ? L"" : L" (not in the deck)");
            AppendMenuW(menu, MF_STRING | (where ? 0 : MF_GRAYED),
                        static_cast<UINT_PTR>(Command::section_first) + id, label.c_str());
        }
        bool first = true;
        for (const SectionSpec& spec : menu_sections(m)) {
            if (spec.keys.empty()) continue;
            const auto where = find_section(text, spec.kind);
            SectionLocation absent;
            absent.kind = std::string(spec.kind);
            const auto keys = insertable_keys(where ? *where : absent);
            const bool analysis = spec.kind == "transient";
            HMENU sub = CreatePopupMenu();
            for (const KeySpec& k : keys) {
                const std::size_t id = key_items_.size();
                key_items_.push_back({std::string(spec.kind), k});
                const std::wstring t = win32::widen(std::string(k.key) + " = " + std::string(k.value) + "\t" + std::string(k.help));
                AppendMenuW(sub, MF_STRING | (analysis && !where ? MF_GRAYED : 0),
                            static_cast<UINT_PTR>(Command::key_first) + id, t.c_str());
            }
            if (keys.empty()) AppendMenuW(sub, MF_STRING | MF_GRAYED, 0, L"(every key is set)");
            if (first) AppendMenuW(menu, MF_SEPARATOR, 0, nullptr);
            first = false;
            AppendMenuW(menu, MF_POPUP, reinterpret_cast<UINT_PTR>(sub),
                        (L"&Insert into " + win32::widen(spec.title)).c_str());
        }
    };
    if (name == "&Device") add_sections(DeckMenu::device);
    else if (name == "&Mesh") add_sections(DeckMenu::mesh);
    else if (name == "&Physics") add_sections(DeckMenu::physics);
    else if (name == "&Simulation") {
        add_commands(name);
        AppendMenuW(menu, MF_SEPARATOR, 0, nullptr);
        add_sections(DeckMenu::simulation);
    } else {
        add_commands(name);
    }
}

void MainWindow::build_accelerators() {
    std::vector<ACCEL> table;
    for (const CommandInfo& c : commands()) {
        if (c.accelerator.key == 0 || c.id == Command::exit) continue;  // Alt+F4 is the system's
        ACCEL a{};
        a.fVirt = static_cast<BYTE>(FVIRTKEY | (c.accelerator.ctrl ? FCONTROL : 0) | (c.accelerator.shift ? FSHIFT : 0) |
                                    (c.accelerator.alt ? FALT : 0));
        a.key = c.accelerator.key;
        a.cmd = static_cast<WORD>(c.id);
        table.push_back(a);
    }
    accelerators_ = CreateAcceleratorTableW(table.data(), static_cast<int>(table.size()));
}

bool MainWindow::pre_translate(MSG& m) {
    if (!hwnd_) return false;
    if (accelerators_ && TranslateAcceleratorW(hwnd_, accelerators_, &m)) return true;
    // Tab and Shift+Tab between the panes' controls; the deck editor keeps its own keys.
    if (m.message == WM_KEYDOWN && m.hwnd == deck_->editor()) return false;
    // Unit 26d (ZD1): Enter in an inspector box makes it take effect (dialog navigation would take it as OK).
    if (m.message == WM_KEYDOWN && m.wParam == VK_RETURN && m.hwnd && GetParent(m.hwnd) == inspector_ &&
        is_box(GetDlgCtrlID(m.hwnd))) {
        commit_box(GetDlgCtrlID(m.hwnd));
        return true;
    }
    return IsDialogMessageW(hwnd_, &m) != FALSE;
}

void MainWindow::update_command_states() {
    HMENU bar = GetMenu(hwnd_);
    const bool busy = runner_.busy();
    const bool has_scene = scene_.has_value();
    const auto enable = [&](Command c, bool on) {
        EnableMenuItem(bar, static_cast<UINT>(c), MF_BYCOMMAND | (on ? MF_ENABLED : MF_GRAYED));
        SendMessageW(toolbar_, TB_ENABLEBUTTON, static_cast<WPARAM>(c), on);
    };
    enable(Command::check, !busy);
    enable(Command::run, !busy);
    enable(Command::stop, runner_.state() == RunState::running);
    enable(Command::open_deck, !busy);
    enable(Command::new_deck, !busy);
    enable(Command::open_result, !busy);
    enable(Command::find_error, !deck_->diagnostics().empty());
    enable(Command::compare_result, !busy && workspace_.has_value());
    enable(Command::compare_previous, !busy && workspace_.has_value() && previous_.has_value());
    enable(Command::stop_comparing, run_b_.has_value());
    enable(Command::mesh_quality, workspace_ && workspace_->device->mesh().dimension() >= 2 &&
                                      workspace_->geometry->element_count() > 0);
    enable(Command::show_regions, workspace_ && !device_view_.hidden.empty());
    for (const Command c : {Command::play_pause, Command::snapshot_first, Command::snapshot_previous,
                            Command::snapshot_next, Command::snapshot_last}) {
        enable(c, can_play());
    }
    for (const Command c : {Command::fit, Command::orbit_left, Command::orbit_right, Command::orbit_up,
                            Command::orbit_down, Command::zoom_in, Command::zoom_out, Command::pan_left,
                            Command::pan_right, Command::pan_up, Command::pan_down}) {
        enable(c, has_scene);
    }
    const auto check = [&](Command c, bool on) {
        CheckMenuItem(bar, static_cast<UINT>(c), MF_BYCOMMAND | (on ? MF_CHECKED : MF_UNCHECKED));
        SendMessageW(toolbar_, TB_CHECKBUTTON, static_cast<WPARAM>(c), on);
    };
    check(Command::probe_mode, probe_mode_);
    check(Command::measure_mode, measure_mode_);
    check(Command::cut_line_mode, cut_mode_);
    check(Command::perspective, settings().perspective);
    check(Command::hover_readout, state_.hover_readout);
}

// --- layout and DPI --------------------------------------------------------------------------------

void MainWindow::apply_dpi(int dpi) {
    dpi_ = dpi;
    if (font_) DeleteObject(font_);
    if (bold_) DeleteObject(bold_);
    font_ = win32::ui_font(dpi);
    bold_ = win32::ui_font(dpi, true);
    const auto set = [&](HWND h) { SendMessageW(h, WM_SETFONT, reinterpret_cast<WPARAM>(font_), TRUE); };
    for (HWND h : {toolbar_, status_, tree_, tabs_, metrics_}) set(h);
    for (HWND c = GetWindow(inspector_, GW_CHILD); c; c = GetWindow(c, GW_HWNDNEXT)) set(c);
    for (HWND c = GetWindow(console_panel_, GW_CHILD); c; c = GetWindow(c, GW_HWNDNEXT)) {
        if (c != console_) set(c);
    }
    make_icons();
    SendMessageW(toolbar_, TB_AUTOSIZE, 0, 0);
    deck_->set_dpi(dpi);
    start_->set_dpi(dpi, font_, bold_);
    run_page_->set_dpi(dpi, font_);
    results_->set_dpi(dpi, font_, bold_);
    layout();
}

void MainWindow::layout() {
    RECT r{};
    GetClientRect(hwnd_, &r);
    SendMessageW(toolbar_, TB_AUTOSIZE, 0, 0);
    SendMessageW(status_, WM_SIZE, 0, 0);
    RECT tr{}, sr{};
    GetWindowRect(toolbar_, &tr);
    GetWindowRect(status_, &sr);
    const int toolbar_h = std::max<int>(tr.bottom - tr.top, scaled(28, dpi_));
    panes_ = layout_panes(r.right, r.bottom, dpi_, state_.layout, toolbar_h, sr.bottom - sr.top);
    const auto place = [](HWND h, const Rect& x) { MoveWindow(h, x.x, x.y, x.width, x.height, TRUE); };
    place(toolbar_, panes_.toolbar);
    // The status bar's parts: the state, the file, the adapter, the cursor.
    const int w = r.right;
    const int parts[4] = {w * 40 / 100, w * 65 / 100, w * 82 / 100, -1};
    SendMessageW(status_, SB_SETPARTS, 4, reinterpret_cast<LPARAM>(parts));
    place(tree_, panes_.tree);
    place(tabs_, panes_.tabs);
    place(inspector_, panes_.inspector);
    place(console_panel_, panes_.console);
    const bool view = tab_ == tab_view || tab_ == tab_curves;
    Rect page = panes_.page;
    if (tab_ != tab_view) page.height += panes_.metrics.height;  // the metrics strip is the device view's
    const Rect whole{panes_.page.x, panes_.page.y, panes_.page.width, panes_.page.height + panes_.metrics.height};
    place(viewport_->hwnd(), page);
    for (HWND h : {deck_->hwnd(), start_->hwnd(), run_page_->hwnd(), results_->hwnd()}) place(h, whole);
    place(metrics_, panes_.metrics);
    ShowWindow(viewport_->hwnd(), view ? SW_SHOW : SW_HIDE);
    ShowWindow(metrics_, tab_ == tab_view ? SW_SHOW : SW_HIDE);
    ShowWindow(deck_->hwnd(), tab_ == tab_deck ? SW_SHOW : SW_HIDE);
    ShowWindow(start_->hwnd(), tab_ == tab_start ? SW_SHOW : SW_HIDE);
    ShowWindow(run_page_->hwnd(), tab_ == tab_run ? SW_SHOW : SW_HIDE);
    ShowWindow(results_->hwnd(), tab_ == tab_results ? SW_SHOW : SW_HIDE);
    // The tab strip is drawn over the page's top: keep the page above it in Z order.
    SetWindowPos(tabs_, HWND_BOTTOM, 0, 0, 0, 0, SWP_NOMOVE | SWP_NOSIZE | SWP_NOACTIVATE);
    layout_inspector();
}

namespace {

// Look A: the inspector's sections, each the controls it lays out (and hides when folded), labels too.
constexpr int first_label = MainWindow::id_first_label;
constexpr int sec_display = 0, sec_colour = 1, sec_show = 2, sec_slice = 3, sec_vectors = 4, sec_curves = 5,
              sec_playback = 6, sec_probe = 7;
constexpr std::array<const wchar_t*, MainWindow::section_count> section_titles{
    L"DISPLAY", L"COLOUR", L"SHOW", L"SLICE AND CLIP", L"VECTORS", L"CURVES", L"PLAYBACK", L"PROBE"};
using M = MainWindow;
const std::array<std::vector<int>, MainWindow::section_count>& section_members() {
    static const std::array<std::vector<int>, MainWindow::section_count> m{
        std::vector<int>{first_label + 0, M::id_view, first_label + 1, M::id_field, first_label + 2, M::id_point},
        std::vector<int>{first_label + 3, M::id_scale, first_label + 4, M::id_colormap, label_range, M::id_range_min,
                         M::id_range_max, M::id_range_auto, first_label + 7, M::id_levels},
        std::vector<int>{M::id_contours, M::id_edges, M::id_contacts, M::id_opaque, M::id_perspective, M::id_enhanced,
                         M::id_axis_box},
        std::vector<int>{first_label + 5, M::id_slice_axis, first_label + 6, M::id_slice_at, M::id_clip, M::id_clip_box,
                         label_clip_x, M::id_clip_x0, M::id_clip_x1, label_clip_x + 1, M::id_clip_y0, M::id_clip_y1,
                         label_clip_x + 2, M::id_clip_z0, M::id_clip_z1},
        std::vector<int>{label_arrows, M::id_arrows, M::id_arrow_count, label_lines, M::id_lines, M::id_line_count,
                         label_from, M::id_from, M::id_seeds, label_axis, M::id_axis, M::id_gap, M::id_exploded_field},
        std::vector<int>{label_swept, M::id_swept, label_sensed, M::id_sensed, label_frequency, M::id_frequency,
                         M::id_log},
        std::vector<int>{label_time, M::id_time, M::id_time_text, M::id_play},
        std::vector<int>{first_label + 8, M::id_cut, M::id_probe_text}};
    return m;
}

// A glyph of the system's icon font (Segoe Fluent Icons, else Segoe MDL2 Assets) in a colour, or a filled
// square (a region's swatch) when the glyph is 0.
struct Icon {
    wchar_t glyph;
    COLORREF colour;
};

// The icons drawn at `size` pixels into a 32-bit image list: the glyph's grey coverage as the alpha, the colour
// premultiplied. Null when no icon font is installed.
HIMAGELIST icon_list(int size, std::span<const Icon> icons, bool grey, int* missing = nullptr) {
    HDC dc = CreateCompatibleDC(nullptr);
    HFONT font = nullptr;
    for (const wchar_t* face : {L"Segoe Fluent Icons", L"Segoe MDL2 Assets"}) {
        font = CreateFontW(-size, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE, DEFAULT_CHARSET, OUT_DEFAULT_PRECIS,
                           CLIP_DEFAULT_PRECIS, ANTIALIASED_QUALITY, DEFAULT_PITCH, face);
        HGDIOBJ old = SelectObject(dc, font);
        wchar_t got[LF_FACESIZE]{};
        GetTextFaceW(dc, LF_FACESIZE, got);
        SelectObject(dc, old);
        if (std::wstring_view(got) == face) break;
        DeleteObject(font);
        font = nullptr;
    }
    if (!font) {
        DeleteDC(dc);
        return nullptr;
    }
    HIMAGELIST list = ImageList_Create(size, size, ILC_COLOR32, static_cast<int>(icons.size()), 0);
    BITMAPINFO bi{};
    bi.bmiHeader = {sizeof(BITMAPINFOHEADER), size, -size, 1, 32, BI_RGB};
    for (const Icon& icon : icons) {
        void* bits = nullptr;
        HBITMAP bmp = CreateDIBSection(dc, &bi, DIB_RGB_COLORS, &bits, nullptr, 0);
        auto* px = static_cast<std::uint32_t*>(bits);
        const std::size_t n = static_cast<std::size_t>(size) * static_cast<std::size_t>(size);
        std::fill(px, px + n, 0u);
        HGDIOBJ old_bmp = SelectObject(dc, bmp);
        const COLORREF c = grey ? RGB(160, 164, 170) : icon.colour;
        if (icon.glyph != 0) {
            WORD index = 0;
            HGDIOBJ old = SelectObject(dc, font);
            GetGlyphIndicesW(dc, &icon.glyph, 1, &index, GGI_MARK_NONEXISTING_GLYPHS);
            if (index == 0xFFFF && missing) ++*missing;
            if (index != 0xFFFF) {
                SetBkMode(dc, TRANSPARENT);
                SetTextColor(dc, RGB(255, 255, 255));
                RECT r{0, 0, size, size};
                DrawTextW(dc, &icon.glyph, 1, &r, DT_CENTER | DT_VCENTER | DT_SINGLELINE | DT_NOPREFIX);
            }
            SelectObject(dc, old);
            GdiFlush();
            for (std::size_t k = 0; k < n; ++k) {
                const std::uint32_t a = (px[k] >> 8) & 0xFF;  // the green channel: the glyph's coverage
                px[k] = a << 24 | (GetRValue(c) * a / 255) << 16 | (GetGValue(c) * a / 255) << 8 | (GetBValue(c) * a / 255);
            }
        } else {  // a swatch: a rounded-off square with a darker edge
            const int m = std::max(size / 8, 1);
            for (int y = m + 1; y < size - m; ++y) {
                for (int x = m; x < size - m; ++x) {
                    const bool edge = x == m || x == size - m - 1 || y == m + 1 || y == size - m - 1;
                    const int f = edge ? 70 : 100;
                    px[static_cast<std::size_t>(y * size + x)] =
                        0xFFu << 24 | static_cast<std::uint32_t>(GetRValue(c) * f / 100) << 16 |
                        static_cast<std::uint32_t>(GetGValue(c) * f / 100) << 8 | static_cast<std::uint32_t>(GetBValue(c) * f / 100);
                }
            }
        }
        SelectObject(dc, old_bmp);
        ImageList_Add(list, bmp, nullptr);
        DeleteObject(bmp);
    }
    DeleteObject(font);
    DeleteDC(dc);
    return list;
}

constexpr COLORREF icon_ink = RGB(45, 50, 58);
// The toolbar's buttons in order: New, Open, Save, Check, Run, Stop, Fit, Probe, Measure, Cut line, Perspective.
constexpr std::array<Icon, 11> toolbar_icons{{{0xE8A5, icon_ink}, {0xE838, icon_ink}, {0xE74E, icon_ink},
                                              {0xE73E, icon_ink}, {0xE768, RGB(46, 157, 79)}, {0xE71A, RGB(196, 43, 28)},
                                              {0xE740, icon_ink}, {0xE81D, icon_ink}, {0xED5E, icon_ink},
                                              {0xE9D9, icon_ink}, {0xE7B8, icon_ink}}};
// The tree's: the device, a semiconductor's and an insulator's swatch, a contact, the mesh, the views, a view,
// the fields, a field, the results.
constexpr int tree_icon_device = 0, tree_icon_semiconductor = 1, tree_icon_insulator = 2, tree_icon_contact = 3,
              tree_icon_mesh = 4, tree_icon_views = 5, tree_icon_view = 6, tree_icon_fields = 7, tree_icon_field = 8,
              tree_icon_results = 9;
constexpr COLORREF icon_soft = RGB(105, 112, 124);
constexpr std::array<Icon, 10> tree_icons{{{0xE950, RGB(47, 111, 214)}, {0, RGB(155, 183, 217)}, {0, RGB(233, 215, 168)},
                                           {0xE945, RGB(176, 120, 20)}, {0xE80A, icon_ink}, {0xE890, icon_ink},
                                           {0xE7B3, icon_soft}, {0xE9D9, icon_ink}, {0xE76C, icon_soft},
                                           {0xE9F9, icon_ink}}};

}  // namespace

void MainWindow::make_icons() {
    for (HIMAGELIST* l : {&toolbar_images_, &toolbar_disabled_, &tree_images_}) {
        if (*l) ImageList_Destroy(*l);
        *l = nullptr;
    }
    const int size = scaled(16, dpi_);
    missing_icons_ = 0;
    toolbar_images_ = icon_list(size, toolbar_icons, false, &missing_icons_);
    toolbar_disabled_ = icon_list(size, toolbar_icons, true);
    tree_images_ = icon_list(size, tree_icons, false, &missing_icons_);
    SendMessageW(toolbar_, TB_SETBITMAPSIZE, 0, toolbar_images_ ? MAKELPARAM(size, size) : MAKELPARAM(0, 0));
    SendMessageW(toolbar_, TB_SETIMAGELIST, 0, reinterpret_cast<LPARAM>(toolbar_images_));
    SendMessageW(toolbar_, TB_SETDISABLEDIMAGELIST, 0, reinterpret_cast<LPARAM>(toolbar_disabled_));
    TreeView_SetImageList(tree_, tree_images_, TVSIL_NORMAL);
}

void MainWindow::set_section_folded(int section, bool folded) {
    if (section < 0 || section >= section_count) return;
    const int bit = 1 << section;
    state_.layout.inspector_folded = folded ? state_.layout.inspector_folded | bit : state_.layout.inspector_folded & ~bit;
    layout_inspector();
    InvalidateRect(GetDlgItem(inspector_, id_section_first + section), nullptr, TRUE);
}

// A section's header: its title in small capitals on a light band, a chevron (down open, right folded).
void MainWindow::draw_section_header(const DRAWITEMSTRUCT& d) {
    const int section = static_cast<int>(d.CtlID) - id_section_first;
    HDC dc = d.hDC;
    RECT r = d.rcItem;
    HBRUSH band = CreateSolidBrush(RGB(244, 245, 248));
    FillRect(dc, &r, band);
    DeleteObject(band);
    RECT line{r.left, r.bottom - 1, r.right, r.bottom};
    HBRUSH edge = CreateSolidBrush(RGB(217, 221, 227));
    FillRect(dc, &line, edge);
    line = {r.left, r.top, r.right, r.top + 1};
    FillRect(dc, &line, edge);
    DeleteObject(edge);
    SetBkMode(dc, TRANSPARENT);
    SetTextColor(dc, RGB(80, 87, 98));
    HGDIOBJ old = SelectObject(dc, bold_);
    RECT text{r.left + scaled(10, dpi_), r.top, r.right - scaled(24, dpi_), r.bottom};
    DrawTextW(dc, section_titles[static_cast<std::size_t>(section)], -1, &text, DT_LEFT | DT_VCENTER | DT_SINGLELINE | DT_NOPREFIX);
    SelectObject(dc, font_);
    RECT chevron{r.right - scaled(24, dpi_), r.top, r.right - scaled(8, dpi_), r.bottom};
    DrawTextW(dc, section_folded(section) ? L"▸" : L"▾", 1, &chevron, DT_RIGHT | DT_VCENTER | DT_SINGLELINE);
    SelectObject(dc, old);
    if ((d.itemState & ODS_FOCUS) != 0) {
        RECT f = r;
        InflateRect(&f, -2, -2);
        DrawFocusRect(dc, &f);
    }
}

// Look A: the sections one under the other, each under its header; in a section a label left of its control
// (a drop-down list, a box, a pair side by side, a slider), the check boxes two to a row where they fit. A section
// folded shows its header alone; one whose controls the view does not use is hidden. The probe's text fills the
// rest.
void MainWindow::layout_inspector() {
    RECT r{};
    GetClientRect(inspector_, &r);
    const int pad = scaled(8, dpi_), row = scaled(23, dpi_), gap = scaled(4, dpi_), head = scaled(22, dpi_);
    const int w = std::max<int>(r.right - 2 * pad, 0);
    const int lw = std::min(scaled(84, dpi_), w / 2);  // the label column
    const int cx = pad + lw + gap, cw = std::max(w - lw - gap, 0);
    const int half = std::max((cw - gap) / 2, 0);
    const auto used = [&](int id) { return !unused_controls_.contains(id); };
    const auto show = [&](int id, bool on) {
        HWND h = child(id);
        const bool visible = (GetWindowLongPtrW(h, GWL_STYLE) & WS_VISIBLE) != 0;
        if (visible != on) ShowWindow(h, on ? SW_SHOWNA : SW_HIDE);
    };
    const auto move = [&](int id, int x, int y, int width, int height) {
        MoveWindow(child(id), x, y, std::max(width, 0), height, TRUE);
        show(id, true);
    };
    int y = 0;
    // A label and its control; a combo box's height includes its list.
    const auto field = [&](int label, int id, bool combo) {
        move(label, pad, y, lw, row);
        move(id, cx, y, cw, combo ? row * 8 : row);
        y += row + gap;
    };
    const auto pair = [&](int label, int a, int b, bool combo) {
        move(label, pad, y, lw, row);
        move(a, cx, y, half, combo ? row * 8 : row);
        move(b, cx + half + gap, y, cw - half - gap, row);
        y += row + gap;
    };
    const auto full = [&](int id) {
        move(id, pad, y, w, row);
        y += row + gap;
    };
    const auto& members = section_members();
    for (int k = 0; k < section_count; ++k) {
        const auto& ids = members[static_cast<std::size_t>(k)];
        const bool any = std::ranges::any_of(ids, used);
        const int header = id_section_first + k;
        if (!any || section_folded(k)) {
            for (const int id : ids) show(id, false);
        }
        if (!any) {
            show(header, false);
            continue;
        }
        move(header, 0, y, r.right, head);
        y += head;
        if (section_folded(k)) continue;
        y += gap + gap / 2;
        // Every member hidden first: each placement below shows the ones in use.
        for (const int id : ids) {
            if (!used(id)) show(id, false);
        }
        switch (k) {
        case sec_display:
            field(first_label + 0, id_view, true);
            field(first_label + 1, id_field, true);
            field(first_label + 2, id_point, true);
            break;
        case sec_colour:
            field(first_label + 3, id_scale, true);
            field(first_label + 4, id_colormap, true);
            if (used(id_range_min)) {  // ZD3: the minimum, the maximum and Auto in a row
                const int third = std::max((cw - 2 * gap) / 3, 0);
                move(label_range, pad, y, lw, row);
                move(id_range_min, cx, y, third, row);
                move(id_range_max, cx + third + gap, y, third, row);
                move(id_range_auto, cx + 2 * (third + gap), y, cw - 2 * (third + gap), row);
                y += row + gap;
            }
            field(first_label + 7, id_levels, false);
            break;
        case sec_show: {
            std::vector<int> boxes;
            for (const int id : ids) {
                if (used(id)) boxes.push_back(id);
            }
            const int col = (w - gap) / 2;
            const bool two = col >= scaled(108, dpi_);
            for (std::size_t j = 0; j < boxes.size(); ++j) {
                if (two) {
                    move(boxes[j], pad + (j % 2 == 0 ? 0 : col + gap), y, j % 2 == 0 ? col : w - col - gap, row);
                    if (j % 2 == 1 || j + 1 == boxes.size()) y += row + gap / 2;
                } else {
                    move(boxes[j], pad, y, w, row);
                    y += row + gap / 2;
                }
            }
            break;
        }
        case sec_slice:
            field(first_label + 5, id_slice_axis, true);
            field(first_label + 6, id_slice_at, false);
            full(id_clip);
            move(id_clip_box, cx, y, cw, row);
            y += row + gap;
            if (used(id_clip_x0)) {  // YD6: per axis the low and the high slider side by side
                for (int axis = 0; axis < 3; ++axis) pair(label_clip_x + axis, id_clip_x0 + 2 * axis, id_clip_x1 + 2 * axis, false);
            }
            break;
        case sec_vectors:
            for (const auto& [label, list, box] : {std::tuple{label_arrows, id_arrows, id_arrow_count},
                                                   std::tuple{label_lines, id_lines, id_line_count},
                                                   std::tuple{label_from, id_from, id_seeds},
                                                   std::tuple{label_axis, id_axis, id_gap}}) {
                if (used(list)) pair(label, list, box, true);
            }
            if (used(id_exploded_field)) move(id_exploded_field, cx, y, cw, row), y += row + gap;
            break;
        case sec_curves:
            for (const auto& [label, list] : {std::pair{label_swept, id_swept}, std::pair{label_sensed, id_sensed},
                                              std::pair{label_frequency, id_frequency}}) {
                if (used(list)) field(label, list, true);
            }
            if (used(id_log)) move(id_log, cx, y, cw, row), y += row + gap;
            break;
        case sec_playback: {  // YD5: the snapshot slider, then its time beside Play
            field(label_time, id_time, false);
            const int button = std::min(scaled(72, dpi_), cw);
            move(id_time_text, cx, y, cw - button - gap, row);
            move(id_play, cx + cw - button, y, button, row);
            y += row + gap;
            break;
        }
        case sec_probe:
            field(first_label + 8, id_cut, false);
            move(id_probe_text, pad, y, w, std::max<int>(r.bottom - y - pad, row * 3));
            y += std::max<int>(r.bottom - y - pad, row * 3);
            break;
        default:
            break;
        }
        y += gap;
    }
}

// Look A: the run's progress at the right of the status bar's first part while it runs (a marquee when the
// points are not counted), hidden otherwise.
void MainWindow::update_status_progress() {
    if (!status_progress_) return;
    const RunState s = runner_.state();
    const bool running = s == RunState::running || s == RunState::stopping;
    if (!running) {
        if (IsWindowVisible(status_progress_) || (GetWindowLongPtrW(status_progress_, GWL_STYLE) & WS_VISIBLE) != 0) {
            ShowWindow(status_progress_, SW_HIDE);
        }
        return;
    }
    RECT part{};
    SendMessageW(status_, SB_GETRECT, status_state, reinterpret_cast<LPARAM>(&part));
    const int width = std::min<int>(scaled(110, dpi_), (part.right - part.left) / 3);
    const int height = std::min<int>(scaled(10, dpi_), part.bottom - part.top);
    MoveWindow(status_progress_, part.right - width - scaled(6, dpi_), part.top + (part.bottom - part.top - height) / 2,
               width, height, TRUE);
    const double f = progress_fraction(runner_.points(), runner_.planned_points());
    const bool marquee = f < 0.0;
    const LONG_PTR style = GetWindowLongPtrW(status_progress_, GWL_STYLE);
    if (((style & PBS_MARQUEE) != 0) != marquee) {
        SetWindowLongPtrW(status_progress_, GWL_STYLE, marquee ? style | PBS_MARQUEE : style & ~static_cast<LONG_PTR>(PBS_MARQUEE));
        SendMessageW(status_progress_, PBM_SETMARQUEE, marquee, 30);
    }
    if (!marquee) SendMessageW(status_progress_, PBM_SETPOS, static_cast<WPARAM>(std::lround(1000.0 * f)), 0);
    if ((style & WS_VISIBLE) == 0) ShowWindow(status_progress_, SW_SHOWNA);
}

void MainWindow::update_title() {
    const std::string t = "NiTCAD - " + document_.title();
    win32::set_text(hwnd_, t);
    set_status(status_file, document_.is_untitled() ? std::string("untitled deck") : document_.path().string());
}

void MainWindow::set_status(int part, const std::string& text, LogLevel level) {
    show_status(part, text);
    if (part == status_state && !text.empty()) log_line(level, text);
}

void MainWindow::show_status(int part, const std::string& text) {
    SendMessageW(status_, SB_SETTEXTW, static_cast<WPARAM>(part), reinterpret_cast<LPARAM>(win32::widen(text).c_str()));
}

void MainWindow::log_line(LogLevel level, const std::string& text) {
    console_log_.add(level, text);
    flush_console();
}

namespace {

// EM_STREAMIN's reader over a string.
struct RtfSource {
    const std::string* text;
    std::size_t at = 0;
};

DWORD CALLBACK read_rtf(DWORD_PTR cookie, LPBYTE buffer, LONG size, LONG* read) {
    auto* s = reinterpret_cast<RtfSource*>(cookie);
    const std::size_t n = std::min(static_cast<std::size_t>(size), s->text->size() - s->at);
    std::memcpy(buffer, s->text->data() + s->at, n);
    s->at += n;
    *read = static_cast<LONG>(n);
    return 0;
}

constexpr int tab_output = 0, tab_points = 1, tab_messages = 2;  // the console's tabs
constexpr int menu_copy = 1, menu_save = 2, menu_clear = 3, menu_fold = 4;

}  // namespace

LogFilter MainWindow::console_text_filter() const noexcept {
    return console_tab_ == tab_messages ? LogFilter::warnings : LogFilter::all;
}

// The log's lines since the last flush, appended to the rich edit (the whole text again when lines it shows
// were dropped, or rebuild); it follows the end while the view is at the bottom, and stays put scrolled up.
// Output shows every line, Messages the warnings and errors, plain black text. The Messages count too.
void MainWindow::flush_console(bool rebuild) {
    if (!console_) return;
    update_console_footer();
    // Hidden (folded, or the Points tab shown), the text waits: it is appended when the tab is shown.
    if (rebuild) console_rebuild_ = true;
    if (state_.layout.console_collapsed || console_tab_ == tab_points) return;
    rebuild = std::exchange(console_rebuild_, false);
    const LogFilter filter = console_text_filter();
    if (filter != console_filter_shown_) rebuild = true;  // the other text tab's lines are shown
    console_filter_shown_ = filter;
    if (console_shown_ < console_log_.dropped()) rebuild = true;
    if (!rebuild && console_shown_ == console_log_.total()) return;
    const std::size_t from = rebuild ? 0 : console_shown_ - console_log_.dropped();
    const bool any_before = !rebuild && console_lines_shown_ > 0;
    std::size_t adding = 0;
    for (std::size_t k = from; k < console_log_.lines().size(); ++k) {
        if (shows(filter, console_log_.lines()[k].level)) ++adding;
    }
    console_shown_ = console_log_.total();
    if (!rebuild && adding == 0) return;
    const std::string rtf = console_rtf(console_log_.lines(), from, filter, any_before, 9, true);
    console_lines_shown_ = (rebuild ? 0 : console_lines_shown_) + adding;

    SCROLLINFO si{sizeof(si), SIF_ALL};
    GetScrollInfo(console_, SB_VERT, &si);
    const bool at_bottom = rebuild || si.nMax <= 0 || si.nPos + static_cast<int>(si.nPage) >= si.nMax;
    CHARRANGE selection{};
    SendMessageW(console_, EM_EXGETSEL, 0, reinterpret_cast<LPARAM>(&selection));
    const auto first = SendMessageW(console_, EM_GETFIRSTVISIBLELINE, 0, 0);
    SendMessageW(console_, WM_SETREDRAW, FALSE, 0);
    RtfSource source{&rtf};
    EDITSTREAM stream{reinterpret_cast<DWORD_PTR>(&source), 0, read_rtf};
    if (rebuild) {
        SendMessageW(console_, EM_STREAMIN, SF_RTF, reinterpret_cast<LPARAM>(&stream));
    } else {
        CHARRANGE end{-1, -1};
        SendMessageW(console_, EM_EXSETSEL, 0, reinterpret_cast<LPARAM>(&end));
        SendMessageW(console_, EM_STREAMIN, SF_RTF | SFF_SELECTION, reinterpret_cast<LPARAM>(&stream));
    }
    if (at_bottom) {
        const auto last = SendMessageW(console_, EM_GETLINECOUNT, 0, 0) - 1;
        const auto start = static_cast<LONG>(SendMessageW(console_, EM_LINEINDEX, static_cast<WPARAM>(last), 0));
        CHARRANGE caret{start, start};
        SendMessageW(console_, EM_EXSETSEL, 0, reinterpret_cast<LPARAM>(&caret));
        SendMessageW(console_, WM_VSCROLL, SB_BOTTOM, 0);
    } else {
        SendMessageW(console_, EM_EXSETSEL, 0, reinterpret_cast<LPARAM>(&selection));
        const auto now = SendMessageW(console_, EM_GETFIRSTVISIBLELINE, 0, 0);
        SendMessageW(console_, EM_LINESCROLL, 0, static_cast<LPARAM>(first - now));
    }
    SendMessageW(console_, WM_SETREDRAW, TRUE, 0);
    InvalidateRect(console_, nullptr, TRUE);
}

// The Points tab from the controller's state (before a run's results are taken: its currents). Hidden, it
// waits until shown, unless forced (a run's end: the results are taken next).
void MainWindow::update_points(bool force) {
    if (!force && (state_.layout.console_collapsed || console_tab_ != tab_points)) {
        points_stale_ = true;
        return;
    }
    points_stale_ = false;
    const RunState s = runner_.state();
    const auto stopped = runner_.stopped_early();
    PointsTable t = points_table(RunView{.rows = runner_.points(),
                                         .iterations = runner_.iterations(),
                                         .summary = runner_.plan_summary() ? &*runner_.plan_summary() : nullptr,
                                         .plan = runner_.plan() ? &*runner_.plan() : nullptr,
                                         .report = runner_.report() ? &*runner_.report() : nullptr,
                                         .seconds = runner_.elapsed_seconds(),
                                         .diagnostics = runner_.diagnostics(),
                                         .stopped = stopped ? &*stopped : nullptr,
                                         .running = s == RunState::running || s == RunState::stopping});
    if (t == points_) return;
    if (t.columns != points_.columns) {
        SendMessageW(console_points_, WM_SETREDRAW, FALSE, 0);
        while (SendMessageW(console_points_, LVM_DELETECOLUMN, 0, 0)) {
        }
        for (std::size_t c = 0; c < t.columns.size(); ++c) {
            std::wstring name = win32::widen(t.columns[c]);
            LVCOLUMNW col{};
            col.mask = LVCF_TEXT | LVCF_WIDTH | LVCF_FMT;
            col.fmt = t.numeric[c] ? LVCFMT_RIGHT : LVCFMT_LEFT;
            const int chars = c == 0 ? 6 : c + 1 == t.columns.size() ? 13 : std::max<int>(static_cast<int>(name.size()) + 2, 10);
            col.cx = scaled(7 * chars + 8, dpi_);
            col.pszText = name.data();
            SendMessageW(console_points_, LVM_INSERTCOLUMNW, c, reinterpret_cast<LPARAM>(&col));
        }
        SendMessageW(console_points_, WM_SETREDRAW, TRUE, 0);
    }
    points_ = std::move(t);
    SendMessageW(console_points_, LVM_SETITEMCOUNT, points_.rows.size(), LVSICF_NOINVALIDATEALL | LVSICF_NOSCROLL);
    InvalidateRect(console_points_, nullptr, FALSE);
    // While a run goes, its point in view.
    const auto live = std::ranges::find(points_.states, PointState::live);
    if (live != points_.states.end()) {
        SendMessageW(console_points_, LVM_ENSUREVISIBLE, static_cast<WPARAM>(live - points_.states.begin()), FALSE);
    }
    update_console_footer();
}

void MainWindow::update_console_footer() {
    // The Messages tab carries the count of warnings and errors.
    const std::size_t n = console_log_.count(LogLevel::warning) + console_log_.count(LogLevel::error);
    if (n == messages_shown_) return;  // the label only when the count changes
    messages_shown_ = n;
    const std::wstring label = n == 0 ? std::wstring(L"Messages") : std::format(L"Messages ({})", n);
    TCITEMW item{};
    item.mask = TCIF_TEXT;
    item.pszText = const_cast<wchar_t*>(label.c_str());
    SendMessageW(GetDlgItem(console_panel_, id_console_tabs), TCM_SETITEMW, tab_messages, reinterpret_cast<LPARAM>(&item));
}
// The tab's text or table over the tab strip along the bottom (Output, Points, Messages); folded, the strip alone.
void MainWindow::layout_console() {
    RECT r{};
    GetClientRect(console_panel_, &r);
    const int strip = std::min<int>(scaled(22, dpi_), r.bottom);
    const int body = std::max<int>(r.bottom - strip, 0);
    const bool open = !state_.layout.console_collapsed;
    MoveWindow(console_points_, 0, 0, r.right, body, TRUE);
    MoveWindow(console_, 0, 0, r.right, body, TRUE);
    MoveWindow(GetDlgItem(console_panel_, id_console_tabs), 0, body, r.right, strip, TRUE);
    ShowWindow(console_points_, open && console_tab_ == tab_points ? SW_SHOW : SW_HIDE);
    ShowWindow(console_, open && console_tab_ != tab_points ? SW_SHOW : SW_HIDE);
}

void MainWindow::console_menu(int x, int y) {
    HMENU m = CreatePopupMenu();
    AppendMenuW(m, MF_STRING, menu_copy, console_tab_ == tab_points ? L"Copy table" : L"Copy all");
    AppendMenuW(m, MF_STRING, menu_save, console_tab_ == tab_points ? L"Save table..." : L"Save output...");
    AppendMenuW(m, MF_STRING | (console_tab_ == tab_points ? MF_GRAYED : 0), menu_clear, L"Clear");
    AppendMenuW(m, MF_SEPARATOR, 0, nullptr);
    AppendMenuW(m, MF_STRING, menu_fold, state_.layout.console_collapsed ? L"Show" : L"Hide");
    if (x == -1 && y == -1) {  // the keyboard's menu key: at the pane's corner
        RECT r{};
        GetWindowRect(console_panel_, &r);
        x = r.left + scaled(8, dpi_);
        y = r.top + scaled(8, dpi_);
    }
    const int chosen = TrackPopupMenu(m, TPM_RETURNCMD | TPM_RIGHTBUTTON, x, y, 0, hwnd_, nullptr);
    DestroyMenu(m);
    on_console(chosen, BN_CLICKED);
}

void MainWindow::on_console(int id, int code) {
    if (code != BN_CLICKED) return;
    if (id == menu_clear) {
        console_log_.clear();
        console_lines_shown_ = 0;
        flush_console(true);
    } else if (id == menu_copy) {
        copy_console();
    } else if (id == menu_save) {
        const bool table = console_tab_ == tab_points;
        const fs::path p = table ? win32::save_dialog(hwnd_, L"CSV (*.csv)\0*.csv\0All files\0*.*\0", L"Save the points table", L"csv")
                                 : win32::save_dialog(hwnd_, L"Logs (*.log)\0*.log\0Text (*.txt)\0*.txt\0All files\0*.*\0",
                                                      L"Save the output", L"log");
        if (!p.empty()) save_console(p);
    } else if (id == menu_fold) {
        set_console_collapsed(!state_.layout.console_collapsed);
    }
}
LRESULT MainWindow::on_console_notify(const NMHDR& n, LPARAM l) {
    if (n.idFrom == static_cast<UINT_PTR>(id_console_tabs) && n.code == TCN_SELCHANGE) {
        set_console_tab(static_cast<int>(SendMessageW(n.hwndFrom, TCM_GETCURSEL, 0, 0)));
        set_console_collapsed(false);  // a tab clicked on the folded console opens it
        return 0;
    }
    if (n.hwndFrom != console_points_) return 0;
    if (n.code == LVN_GETDISPINFOW) {
        auto& d = reinterpret_cast<NMLVDISPINFOW*>(l)->item;
        if ((d.mask & LVIF_TEXT) && d.iItem >= 0 && static_cast<std::size_t>(d.iItem) < points_.rows.size()) {
            const auto& cells = points_.rows[static_cast<std::size_t>(d.iItem)];
            const std::wstring w = d.iSubItem >= 0 && static_cast<std::size_t>(d.iSubItem) < cells.size()
                                       ? win32::widen(cells[static_cast<std::size_t>(d.iSubItem)])
                                       : std::wstring();
            if (d.cchTextMax > 0) {
                const auto n2 = std::min<std::size_t>(w.size(), static_cast<std::size_t>(d.cchTextMax - 1));
                std::wmemcpy(d.pszText, w.data(), n2);
                d.pszText[n2] = L'\0';
            }
        }
        return 0;
    }
    if (n.code == NM_CUSTOMDRAW) {
        // The rows coloured by state: the point being solved on light blue, an unconverged point amber, the
        // points not yet solved grey; "converged" in green.
        auto* cd = reinterpret_cast<NMLVCUSTOMDRAW*>(l);
        switch (cd->nmcd.dwDrawStage) {
        case CDDS_PREPAINT: return CDRF_NOTIFYITEMDRAW;
        case CDDS_ITEMPREPAINT: return CDRF_NOTIFYSUBITEMDRAW;
        case CDDS_ITEMPREPAINT | CDDS_SUBITEM: {
            const auto k = static_cast<std::size_t>(cd->nmcd.dwItemSpec);
            if (k >= points_.states.size()) return CDRF_DODEFAULT;
            const PointState st = points_.states[k];
            const bool status = static_cast<std::size_t>(cd->iSubItem) + 1 == points_.columns.size();
            cd->clrTextBk = st == PointState::live ? RGB(232, 241, 254) : GetSysColor(COLOR_WINDOW);
            cd->clrText = st == PointState::pending ? RGB(160, 166, 174)
                          : st == PointState::warning ? RGB(168, 115, 12)
                          : st == PointState::live    ? RGB(31, 95, 191)
                          : status                    ? RGB(46, 138, 70)
                                                      : GetSysColor(COLOR_WINDOWTEXT);
            return CDRF_NEWFONT;
        }
        default: return CDRF_DODEFAULT;
        }
    }
    return 0;
}

void MainWindow::set_console_tab(int tab) {
    if (tab < tab_output || tab > tab_messages) return;
    console_tab_ = tab;
    SendMessageW(GetDlgItem(console_panel_, id_console_tabs), TCM_SETCURSEL, static_cast<WPARAM>(tab), 0);
    layout_console();
    // What waited while hidden.
    if (tab != tab_points) flush_console();
    else if (points_stale_) update_points();
}

void MainWindow::set_console_collapsed(bool collapsed) {
    if (state_.layout.console_collapsed == collapsed) return;
    state_.layout.console_collapsed = collapsed;
    layout();
    if (!collapsed) {
        if (console_tab_ == tab_points && points_stale_) update_points();
        else flush_console();
    }
}

bool MainWindow::copy_console() {
    const std::wstring text = win32::widen(console_tab_ == tab_points ? points_text(points_, '\t')
                                                                          : console_log_.text(console_text_filter()));
    if (!OpenClipboard(hwnd_)) return false;
    EmptyClipboard();
    bool ok = false;
    if (HGLOBAL g = GlobalAlloc(GMEM_MOVEABLE, (text.size() + 1) * sizeof(wchar_t))) {
        if (void* p = GlobalLock(g)) {
            std::memcpy(p, text.c_str(), (text.size() + 1) * sizeof(wchar_t));
            GlobalUnlock(g);
            ok = SetClipboardData(CF_UNICODETEXT, g) != nullptr;
        }
        if (!ok) GlobalFree(g);
    }
    CloseClipboard();
    return ok;
}

bool MainWindow::save_console(const fs::path& path) {
    const bool table = console_tab_ == tab_points;
    std::ofstream out(path, std::ios::binary);
    out << (table ? points_text(points_, ',') : console_log_.text(console_text_filter()));
    out.close();
    if (!out) {
        set_status(status_state, "Could not write " + path.filename().string(), LogLevel::error);
        return false;
    }
    set_status(status_state, (table ? "Saved the points table to " : "Saved the output to ") + path.filename().string());
    return true;
}

int MainWindow::console_line_count() const {
    return console_lines_shown_ == 0 ? 0 : static_cast<int>(SendMessageW(console_, EM_GETLINECOUNT, 0, 0));
}

int MainWindow::console_row_count() const {
    return static_cast<int>(SendMessageW(console_points_, LVM_GETITEMCOUNT, 0, 0));
}

bool MainWindow::console_body_shown() const {
    return (GetWindowLongPtrW(console_tab_ == tab_points ? console_points_ : console_, GWL_STYLE) & WS_VISIBLE) != 0;
}

std::string MainWindow::status(int part) const {
    const auto n = LOWORD(SendMessageW(status_, SB_GETTEXTLENGTHW, static_cast<WPARAM>(part), 0));
    std::wstring w(static_cast<std::size_t>(n) + 1, L'\0');
    SendMessageW(status_, SB_GETTEXTW, static_cast<WPARAM>(part), reinterpret_cast<LPARAM>(w.data()));
    w.resize(n);
    return win32::narrow(w);
}

// --- tabs, settings, scenes ------------------------------------------------------------------------

ViewSettings& MainWindow::settings() noexcept { return tab_ == tab_curves ? curve_view_ : device_view_; }

void MainWindow::select_tab(int tab) {
    hide_readout();
    tab_ = std::clamp(tab, 0, tab_count - 1);
    SendMessageW(tabs_, TCM_SETCURSEL, static_cast<WPARAM>(tab_), 0);
    layout();
    switch (tab_) {
    case tab_deck:
        SetFocus(deck_->editor());
        break;
    case tab_start:
        start_->set(state_.recent, examples_);
        break;
    case tab_run:
        run_page_->plot().prepare(renderer());  // its swap chain before a run, not during it
        run_page_->update(runner_, renderer(), true);
        break;
    case tab_results:
        results_->show(workspace(), run_b());
        break;
    default:
        sync_inspector();
        refresh();
        break;
    }
}

bool MainWindow::open_example(std::size_t index) {
    if (index >= examples_.size() || runner_.busy()) return false;
    fs::path to = options_.examples_copy_dir;
    if (to.empty()) {
        to = win32::documents_dir() / L"NiTCAD" / L"examples";
    }
    auto copy = copy_example(examples_[index], to);
    if (!copy) {
        set_status(status_state, copy.error().message);
        return false;
    }
    return open_deck(*copy);
}

void MainWindow::start_command(int command) { execute(static_cast<Command>(command)); }

void MainWindow::start_recent(std::size_t index) { execute(static_cast<Command>(static_cast<int>(Command::recent_first) + static_cast<int>(index))); }

void MainWindow::start_example(std::size_t index) {
    if (confirm_discard()) open_example(index);
}

void MainWindow::run_page_stop() { execute(Command::stop); }

void MainWindow::default_settings() {
    device_view_ = {};
    curve_view_ = {};
    device_view_.enhanced = curve_view_.enhanced = state_.enhanced_look;  // Unit 26f (KD5)
    device_view_.axis_box = curve_view_.axis_box = state_.axis_box;  // Unit 26h (JQ3)
    if (!workspace_) return;
    const auto views = available_views(*workspace_, run_b());
    const bool results = workspace_->has_results();
    const int dim = workspace_->device->mesh().dimension();
    // The device tab: the potential (or, before a run, the regions); 1D: the potential's profile.
    device_view_.kind = dim == 1 ? ViewKind::profile : results ? ViewKind::field : ViewKind::solid;
    if (!results) {
        device_view_.field = "doping";
        if (dim == 1) device_view_.kind = ViewKind::profile;
    }
    if (dim == 1 && !results) device_view_.kind = ViewKind::profile;
    // The curves tab: the run's curve, or a band diagram.
    curve_view_.kind = ViewKind::iv;
    for (const ViewKind k : views) {
        if (k == ViewKind::iv || k == ViewKind::cv || k == ViewKind::current) {
            curve_view_.kind = k;
            break;
        }
        if (k == ViewKind::bands && dim == 1) curve_view_.kind = k;  // 2D and 3D: only along a cut line
    }
    if (curve_view_.kind == ViewKind::iv && workspace_->sweep) {
        // The current of the contact carrying the most (a gate sweep's drain), on a log axis.
        const auto& last = workspace_->sweep->points.back().terminal_current;
        const auto contacts = workspace_->device->contacts();
        std::size_t sense = 0;
        for (std::size_t k = 1; k < last.size() && k < contacts.size(); ++k) {
            if (std::abs(last[k]) > std::abs(last[sense])) sense = k;
        }
        if (sense < contacts.size()) curve_view_.sense = contacts[sense].name;
        curve_view_.log = true;
    }
}

bool MainWindow::results_show_pass(std::size_t k) {
    if (!workspace_ || !workspace_->adapt || k >= workspace_->adapt->passes.size()) return false;
    const std::size_t last = workspace_->adapt->passes.size() - 1;
    const auto shown = [&] {
        if (workspace_->device->mesh().dimension() >= 2) device_view_.kind = ViewKind::mesh;
        select_tab(tab_view);
        refresh();
        return true;
    };
    if (workspace_->pass == k || (k == last && !workspace_->pass)) return shown();
    if (k == last && adapt_result_) {
        // The final mesh: the results, back.
        Workspace r = std::move(*adapt_result_);
        adapt_result_.reset();
        set_workspace(std::move(r), true);
        return shown();
    }
    auto p = workspace_of_pass(adapt_result_ ? *adapt_result_ : *workspace_, k);
    if (!p) {
        set_status(status_state, p.error().message);
        return false;
    }
    if (!workspace_->pass) {
        adapt_result_ = std::move(*workspace_);
        workspace_.reset();
        ++workspace_generation_;  // Unit 26i (ND5): the metrics again
    }
    set_workspace(std::move(*p), true);
    return shown();
}

void MainWindow::set_workspace(Workspace w, bool keep_settings) {
    play(false);  // YD5: a new workspace stops playback; its film is fitted again
    film_key_.clear();
    if (!w.pass) adapt_result_.reset();  // a pass's mesh is shown only beside its own results
    // The results a run (or the structure preview) of the same deck replaces: the previous run (XD6).
    if (workspace_ && workspace_->has_results() && workspace_->source == w.source) previous_ = std::move(*workspace_);
    workspace_ = std::move(w);
    ++workspace_generation_;  // Unit 26i (ND5)
    preview_ = false;
    prober_.emplace(*workspace_);
    probe_.reset();
    measure_from_.reset();
    cut_from_.reset();
    // A cut line stays while the device keeps its dimension (a run of the deck it was drawn on).
    if (cut_ && cut_->from_um.size() != static_cast<std::size_t>(workspace_->device->mesh().dimension())) cut_.reset();
    if (!keep_settings) default_settings();
    // Hidden regions stay while the device has them (YD4).
    const auto regions = workspace_->device->regions();
    std::erase_if(device_view_.hidden, [&](const std::string& name) {
        return workspace_->device->mesh().dimension() < 2 ||
               std::ranges::find(regions, name, &device::Region::name) == regions.end();
    });
    fill_tree();
    sync_inspector();
    results_->show(&*workspace_, run_b());
    refresh();
}

void MainWindow::set_run_b(std::optional<Workspace> b) {
    run_b_ = std::move(b);
    // A difference map without a run B to take it from: the field itself.
    if (device_view_.kind == ViewKind::difference && !(workspace_ && run_b_ && comparable(*workspace_, *run_b_))) {
        device_view_.kind = workspace_ && workspace_->device->mesh().dimension() == 1 ? ViewKind::profile : ViewKind::field;
        device_view_.scale.clear();
        device_view_.range.reset();
    }
    fill_tree();
    sync_inspector();
    results_->show(workspace(), run_b());
    refresh();
}

void MainWindow::refresh() {
    hide_readout();
    if (tab_ != tab_view && tab_ != tab_curves) return;
    scene_.reset();
    values_.clear();
    ++values_generation_;  // Unit 26i (ND5)
    values_title_.clear();
    if (!workspace_) {
        viewport_->show_message("Open a deck (Ctrl+O) or a result file (Ctrl+Shift+R).");
        update_metrics();
        update_command_states();
        return;
    }
    if (!renderer_) {
        viewport_->show_message("No Direct3D 12 device: " + renderer_error_.value_or("unknown"));
        return;
    }
    ViewSettings& s = settings();
    // Playback (YD5): a transient run's field view on one scale over every snapshot, fitted again only when the
    // field, its scale, range or hidden regions change.
    if (tab_ == tab_view && takes_film(*workspace_, s)) {
        std::string key = s.field + "|" + s.scale + "|" +
                          (s.range ? std::format("{},{}", s.range->first, s.range->second) : std::string());
        for (const std::string& r : s.hidden) key += "|" + r;
        if (!s.film || key != film_key_) {
            auto film = film_scale(*workspace_, s);
            s.film = film ? std::optional<Film>(std::move(*film)) : std::nullopt;
            film_key_ = std::move(key);
        }
    } else {
        s.film.reset();
    }
    // The cut line where the view uses it: plotted along, or drawn on a 2D or 3D device view.
    ViewSettings shown = s;
    const bool along = s.kind == ViewKind::profile || s.kind == ViewKind::bands;
    if (along || (!is_plot(s.kind) && workspace_->device->mesh().dimension() > 1)) shown.cut = cut_;
    auto scene = build_scene(*workspace_, shown, run_b());
    if (!scene) {
        viewport_->show_message(scene.error().message);
        if (tab_ == tab_view) sync_range(true);
        update_metrics();
        update_command_states();
        return;
    }
    scene_ = std::move(*scene);
    if (tab_ == tab_view) sync_range(true);  // ZD3: Auto's boxes show the scale drawn
    // Hidden regions are not probed (YD4).
    std::vector<std::int32_t> hidden;
    if (takes_hidden(s.kind)) {
        const auto regions = workspace_->device->regions();
        for (const std::string& name : s.hidden) {
            const auto r = std::ranges::find(regions, name, &device::Region::name);
            if (r != regions.end()) hidden.push_back(static_cast<std::int32_t>(r - regions.begin()));
        }
    }
    prober_->hide(hidden);
    // The shown field's values, for probes.
    if (s.kind == ViewKind::difference) {
        // B - A, as the map shows it.
        const auto f = render::find_field(s.field);
        const bool doping = f && (*f == render::Field::net_doping || *f == render::Field::donors || *f == render::Field::acceptors);
        const results::BiasPoint none;
        const results::BiasPoint* pa = &none;
        const results::BiasPoint* pb = &none;
        if (f && !doping && !workspace_->points.empty()) {
            const std::size_t k = std::min(s.point.value_or(workspace_->points.size() - 1), workspace_->points.size() - 1);
            pa = &workspace_->points[k];
            const auto j = render::matching_point(*pa, run_b()->points);
            pb = j ? &run_b()->points[*j] : nullptr;
        }
        if (f && pb) {
            auto d = render::difference_field(*workspace_->device, *run_b()->device, *workspace_->geometry, *pa, *pb, *f);
            if (d) {
                values_ = std::move(d->values);
                ++values_generation_;  // Unit 26i (ND5)
                values_title_ = std::move(d->title);
            }
        }
    } else if (const auto f = render::find_field(s.field);
               f && !is_plot(s.kind) && s.kind != ViewKind::quality && s.kind != ViewKind::exploded) {
        const bool doping = *f == render::Field::net_doping || *f == render::Field::donors || *f == render::Field::acceptors;
        if (doping || workspace_->has_results()) {
            const std::size_t k = s.point.value_or(workspace_->points.empty() ? 0 : workspace_->points.size() - 1);
            const results::BiasPoint none;
            const results::BiasPoint& p = doping || k >= workspace_->points.size() ? none : workspace_->points[k];
            auto v = doping ? render::field_values(*workspace_->device, p, *f)
                            : render::field_values(*workspace_->device, *workspace_->geometry, p, *f);
            if (v) {
                values_ = std::move(*v);
                ++values_generation_;  // Unit 26i (ND5)
                values_title_ = std::string(render::field_info(*f).title);
            }
        }
    }
    viewport_->show(&*renderer_, &*scene_);
    update_metrics();
    update_command_states();
}

void MainWindow::present_camera() {
    if (!scene_) return;
    if (auto* f = std::get_if<render::FieldScene>(&*scene_)) {
        apply_camera(*f, settings());
        viewport_->present_camera();  // Unit 26i (ND2): the camera only
        update_metrics();
    }
}

void MainWindow::update_metrics() {
    // Unit 26i (ND5): what depends on the workspace or the field is computed once for it, not every frame.
    std::string text;
    if (workspace_) {
        const auto& m = workspace_->device->mesh();
        if (metrics_cache_.workspace != workspace_generation_) {
            metrics_cache_ = {};
            metrics_cache_.workspace = workspace_generation_;
            double lo[3]{1e300, 1e300, 1e300}, hi[3]{-1e300, -1e300, -1e300};
            for (const auto& p : m.points()) {
                for (int k = 0; k < m.dimension(); ++k) {
                    lo[k] = std::min(lo[k], p[static_cast<std::size_t>(k)]);
                    hi[k] = std::max(hi[k], p[static_cast<std::size_t>(k)]);
                }
            }
            std::string& d = metrics_cache_.device;
            d = "  Device ";
            for (int k = 0; k < m.dimension(); ++k) d += (k ? " × " : "") + format_um((hi[k] - lo[k]) * 1e4);
            d += " µm     Mesh " + std::to_string(m.node_count()) + " nodes";
            if (workspace_->geometry->element_count() > 0) {
                d += ", " + std::to_string(workspace_->geometry->element_count()) +
                     (m.dimension() == 3 ? " tetrahedra" : m.dimension() == 2 ? " triangles" : " segments");
            }
        }
        text = metrics_cache_.device;
        const ViewKind kind = settings().kind;
        if (kind == ViewKind::quality || kind == ViewKind::quality_histogram) {
            // XD7: the minimum, the median and the count below the threshold.
            if (!metrics_cache_.quality) {
                metrics_cache_.quality.emplace();
                if (const auto angles = render::element_quality(*workspace_->geometry); angles && !angles->empty()) {
                    const auto q = render::quality_summary(*angles, m.dimension());
                    *metrics_cache_.quality = std::format("     Min angle {:.1f}°, median {:.1f}°, {} below {:g}°", q.min_deg,
                                                    q.median_deg, q.below, q.threshold_deg);
                }
            }
            text += *metrics_cache_.quality;
        } else if (!values_.empty()) {
            // ZD7: the field's symbol and unit, as the tree and the colour bar show them.
            const bool difference = kind == ViewKind::difference;
            if (metrics_cache_.values != values_generation_ || metrics_cache_.difference != difference) {
                metrics_cache_.values = values_generation_;
                metrics_cache_.difference = difference;
                double top = -1e300;
                for (const double v : values_) if (std::isfinite(v)) top = std::max(top, v);
                metrics_cache_.field = "     " + max_metric(values_title_, difference ? render::difference_scale(values_).hi() : top, difference);
            }
            text += metrics_cache_.field;
        }
        if (const Workspace* b = run_b()) {
            // XD6: run B, and the threshold shift when both runs have one.
            text += "     Run B " + b->source.filename().string();
            for (const FigureShift& f : figure_shifts(*workspace_, *b)) {
                if (!f.name.starts_with("Threshold")) continue;
                text += std::format(", ΔVth {:+.4g} {}", f.delta(), f.unit);
                break;
            }
        }
        if (scene_ && viewport_->last_present_ms() > 0.0) text += "     Drawn in " + std::format("{:.1f}", viewport_->last_present_ms()) + " ms";
    }
    win32::set_text(metrics_, text);
}

// The inspector shows the active tab's settings.
void MainWindow::sync_inspector() {
    ViewSettings& s = settings();
    std::vector<std::string> views;
    int selected = 0;
    if (workspace_) {
        for (const ViewKind k : available_views(*workspace_, run_b())) {
            if (!belongs(k, tab_, workspace_->device->mesh().dimension())) continue;
            if (k == s.kind) selected = static_cast<int>(views.size());
            views.emplace_back(to_string(k));
        }
    }
    combo_fill(child(id_view), views, selected);
    std::vector<std::string> fields = workspace_ ? available_fields(*workspace_) : std::vector<std::string>{};
    const auto f = std::ranges::find(fields, s.field);
    combo_fill(child(id_field), fields, f == fields.end() ? 0 : static_cast<int>(f - fields.begin()));
    // ZD2: each point by its varying biases or its time.
    const std::vector<std::string> points = workspace_ ? point_labels(*workspace_) : std::vector<std::string>{};
    combo_fill(child(id_point), points,
               points.empty() ? 0 : static_cast<int>(s.point.value_or(points.size() - 1)));
    const auto sc = std::ranges::find(scale_names, s.scale.empty() ? std::string_view("(default)") : std::string_view(s.scale));
    SendMessageW(child(id_scale), CB_SETCURSEL, static_cast<WPARAM>(sc - scale_names.begin()), 0);
    const auto cm = std::ranges::find(colormap_names, std::string_view(s.colormap));
    SendMessageW(child(id_colormap), CB_SETCURSEL, static_cast<WPARAM>(cm == colormap_names.end() ? 0 : cm - colormap_names.begin()), 0);
    const auto check = [&](int id, bool on) { Button_SetCheck(child(id), on ? BST_CHECKED : BST_UNCHECKED); };
    check(id_contours, s.contours);
    check(id_edges, s.edges);
    check(id_contacts, s.contacts);
    check(id_opaque, s.opaque);
    check(id_perspective, s.perspective);
    check(id_enhanced, s.enhanced);
    check(id_axis_box, s.axis_box);
    check(id_clip, s.clip.has_value());
    SendMessageW(child(id_slice_axis), CB_SETCURSEL, s.slices.empty() ? 0 : static_cast<WPARAM>(s.slices.front().axis + 1), 0);
    // The boxes (ZD1); a clipping box's text stays while there is none, to be checked again.
    for (const int id : {id_slice_at, id_levels, id_cut, id_arrow_count, id_line_count, id_seeds, id_gap}) show_box(id);
    if (s.clip) show_box(id_clip_box);
    sync_range();
    // Unit 26c (YD1-YD3).
    const auto select = [&](int id, std::span<const std::string_view> names, const std::string& name) {
        const auto it = std::ranges::find(names, name.empty() ? std::string_view("none") : std::string_view(name));
        SendMessageW(child(id), CB_SETCURSEL, static_cast<WPARAM>(it == names.end() ? 0 : it - names.begin()), 0);
    };
    select(id_arrows, vector_names, s.vectors);
    select(id_lines, vector_names, s.streamlines);
    std::vector<std::string> contacts;
    if (workspace_) {
        for (const auto& c : workspace_->device->contacts()) contacts.push_back(c.name);
    }
    // A contact list with its first entry for none chosen; the setting's contact selected.
    const auto contact_list = [&](int id, const char* first, const std::string& name) {
        std::vector<std::string> items{first};
        items.insert(items.end(), contacts.begin(), contacts.end());
        const auto it = std::ranges::find(contacts, name);
        combo_fill(child(id), items, it == contacts.end() ? 0 : static_cast<int>(it - contacts.begin()) + 1);
    };
    contact_list(id_from, "(automatic)", s.from);
    SendMessageW(child(id_axis), CB_SETCURSEL, static_cast<WPARAM>(std::clamp(s.axis, 0, 2)), 0);
    check(id_exploded_field, s.exploded_field);
    contact_list(id_swept, s.kind == ViewKind::current ? "(all)" : "(automatic)", s.contact);
    contact_list(id_sensed, "(automatic)", s.sense);
    std::vector<std::string> frequencies;
    if (workspace_ && workspace_->small_signal) {
        for (const double hz : workspace_->small_signal->frequency_Hz) frequencies.push_back(std::format("{:g} Hz", hz));
    }
    combo_fill(child(id_frequency), frequencies, static_cast<int>(s.frequency.value_or(0)));
    check(id_log, s.log);
    sync_clip_sliders();
    sync_playback();
    show_inspector_groups();
}

// YD6: the box's sides on the sliders (the whole extent without a box).
void MainWindow::sync_clip_sliders() {
    if (!workspace_) return;
    const ViewSettings& s = settings();
    const auto extent = extent_um();
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const auto [lo, hi] = extent[axis];
        const int id = id_clip_x0 + 2 * static_cast<int>(axis);
        SendMessageW(child(id), TBM_SETPOS, TRUE, s.clip ? clip_step(lo, hi, s.clip->lo_um[axis]) : 0);
        SendMessageW(child(id + 1), TBM_SETPOS, TRUE, s.clip ? clip_step(lo, hi, s.clip->hi_um[axis]) : clip_steps);
    }
}

void MainWindow::sync_playback() {
    const ViewSettings& s = settings();
    const std::size_t n = workspace_ ? workspace_->points.size() : 0;
    const std::size_t k = n == 0 ? 0 : std::min(s.point.value_or(n - 1), n - 1);
    SendMessageW(child(id_time), TBM_SETRANGE, FALSE, MAKELPARAM(0, n == 0 ? 0 : n - 1));
    SendMessageW(child(id_time), TBM_SETPOS, TRUE, static_cast<LPARAM>(k));
    win32::set_text(child(id_time_text), can_play() ? render::time_label(workspace_->transient->snapshots[k].time_s) +
                                                          std::format(" ({} of {})", k + 1, n)
                                                    : std::string());
    win32::set_text(child(id_play), playing_ ? "Pa&use" : "Pla&y");
    SendMessageW(child(id_point), CB_SETCURSEL, static_cast<WPARAM>(k), 0);
}

bool MainWindow::can_play() const {
    return workspace_ && takes_film(*workspace_, device_view_) && workspace_->points.size() >= 2;
}

void MainWindow::play(bool on) {
    if (on && !can_play()) return;
    if (on) {
        if (tab_ != tab_view) select_tab(tab_view);
        const std::size_t n = workspace_->points.size();
        // From the snapshot shown, or from the first when the last is shown.
        const std::size_t k = std::min(device_view_.point.value_or(n - 1), n - 1);
        play_from_ = k + 1 == n ? 0 : k;
        play_start_ = std::chrono::steady_clock::now();
        playing_ = true;
        show_snapshot(play_from_);
        SetTimer(hwnd_, playback_timer, playback_tick_ms, nullptr);
        set_status(status_state, "Playing the snapshots.");
    } else if (playing_) {
        KillTimer(hwnd_, playback_timer);
        playing_ = false;
        set_status(status_state, "Paused.");
    }
    sync_playback();
}

void MainWindow::show_snapshot(std::size_t k) {
    if (!workspace_ || workspace_->points.empty()) return;
    device_view_.point = std::min(k, workspace_->points.size() - 1);
    if (tab_ == tab_view) refresh();
    sync_playback();
}

void MainWindow::on_playback_timer() {
    if (!playing_ || !can_play() || tab_ != tab_view) {
        play(false);
        return;
    }
    const std::size_t n = workspace_->points.size();
    const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - play_start_).count();
    const std::size_t k = playback_frame(play_from_, n, elapsed);
    if (k != device_view_.point.value_or(n - 1)) show_snapshot(k);
    if (k + 1 == n) play(false);
}

void MainWindow::on_slider(HWND bar) {
    if (!workspace_) return;
    const int position = static_cast<int>(SendMessageW(bar, TBM_GETPOS, 0, 0));
    if (bar == child(id_time)) {
        // A snapshot chosen: playback pauses on it.
        play(false);
        show_snapshot(static_cast<std::size_t>(std::max(position, 0)));
        return;
    }
    const int id = GetDlgCtrlID(bar);
    if (id < id_clip_x0 || id > id_clip_z1 || workspace_->device->mesh().dimension() != 3) return;
    // A side of the clipping box (YD6): the box made over the whole device if there was none, the other side
    // of the axis moved with it rather than crossed.
    ViewSettings& s = device_view_;
    const auto extent = extent_um();
    if (!s.clip) {
        ViewSettings::Clip c;
        for (std::size_t k = 0; k < 3; ++k) {
            c.lo_um[k] = clip_value_um(extent[k].first, extent[k].second, 0);
            c.hi_um[k] = clip_value_um(extent[k].first, extent[k].second, clip_steps);
        }
        s.clip = c;
    }
    const auto axis = static_cast<std::size_t>((id - id_clip_x0) / 2);
    const bool high = (id - id_clip_x0) % 2 == 1;
    const double v = clip_value_um(extent[axis].first, extent[axis].second, position);
    (high ? s.clip->hi_um : s.clip->lo_um)[axis] = v;
    if (high && s.clip->lo_um[axis] > v) s.clip->lo_um[axis] = v;
    if (!high && s.clip->hi_um[axis] < v) s.clip->hi_um[axis] = v;
    sync_inspector();
    refresh();
}

std::array<std::pair<double, double>, 3> MainWindow::extent_um() const {
    std::array<std::pair<double, double>, 3> out{};
    if (!workspace_) return out;
    const auto& m = workspace_->device->mesh();
    for (std::size_t k = 0; k < 3; ++k) out[k] = {1e300, -1e300};
    for (const auto& p : m.points()) {
        for (std::size_t k = 0; k < 3; ++k) {
            out[k].first = std::min(out[k].first, p[k] * 1e4);
            out[k].second = std::max(out[k].second, p[k] * 1e4);
        }
    }
    return out;
}

void MainWindow::show_inspector_groups() {
    const ViewSettings& s = settings();
    const int dim = workspace_ ? workspace_->device->mesh().dimension() : 0;
    // Arrows and streamlines over the field views of a 2D or 3D run's results (YD1).
    const bool vectors = tab_ == tab_view && dim >= 2 && workspace_->has_results() && takes_vectors(s.kind);
    const bool exploded = tab_ == tab_view && s.kind == ViewKind::exploded;
    const bool curves = tab_ == tab_curves && (s.kind == ViewKind::iv || s.kind == ViewKind::cv || s.kind == ViewKind::current);
    // Look A: recorded as used or not; layout_inspector shows the used ones whose section is open.
    const auto show = [&](std::initializer_list<int> ids, bool on) {
        for (const int id : ids) {
            if (on) unused_controls_.erase(id);
            else unused_controls_.insert(id);
        }
    };
    show({label_arrows, id_arrows, id_arrow_count, label_lines, id_lines, id_line_count, label_from, id_from, id_seeds}, vectors);
    show({label_axis, id_axis, id_gap, id_exploded_field}, exploded);
    // ZD3: the colour range on the views drawn on a colour scale of the field.
    show({label_range, id_range_min, id_range_max, id_range_auto},
         tab_ == tab_view && workspace_ &&
             (s.kind == ViewKind::field || s.kind == ViewKind::slice || s.kind == ViewKind::iso || exploded));
    win32::set_text(GetDlgItem(inspector_, label_swept), s.kind == ViewKind::current ? "Contact" : "Swept contact");
    show({label_swept, id_swept, id_log}, curves);
    show({label_sensed, id_sensed}, curves && s.kind != ViewKind::current);
    show({label_frequency, id_frequency}, curves && s.kind == ViewKind::cv);
    // Playback on a transient run's field view (YD5); the box's sliders where a 3D view clips (YD6).
    show({label_time, id_time, id_time_text, id_play}, tab_ == tab_view && can_play());
    show({id_enhanced}, tab_ == tab_view && dim == 3 && !is_plot(s.kind));  // Unit 26f (KD5)
    show({id_axis_box}, tab_ == tab_view && dim == 3 && !is_plot(s.kind));  // Unit 26h (JQ3): as the look's box
    const bool clips = tab_ == tab_view && dim == 3 && (s.kind == ViewKind::solid || s.kind == ViewKind::field);
    show({label_clip_x, label_clip_x + 1, label_clip_x + 2, id_clip_x0, id_clip_x1, id_clip_y0, id_clip_y1, id_clip_z0,
          id_clip_z1},
         clips);
    layout_inspector();
}

bool MainWindow::pickable() const {
    return scene_ && std::holds_alternative<render::FieldScene>(*scene_) && prober_ &&
           (tab_ == tab_curves ? curve_view_ : device_view_).kind != ViewKind::exploded;
}

// A slice view starts with a plane through the middle of the device along z, an iso view with a level.
void MainWindow::view_defaults(ViewSettings& s) {
    if (!workspace_) return;
    if (s.kind == ViewKind::slice && s.slices.empty()) {
        const auto& m = workspace_->device->mesh();
        double lo = 1e300, hi = -1e300;
        for (const auto& p : m.points()) {
            lo = std::min(lo, p[2]);
            hi = std::max(hi, p[2]);
        }
        s.slices = {{2, 0.5 * (lo + hi) * 1e4}};
    }
    if (s.kind == ViewKind::iso && s.levels.empty()) s.levels = {1e16};
}

void MainWindow::on_inspector(int id, int code) {
    ViewSettings& s = settings();
    const bool selection = code == CBN_SELCHANGE;
    const bool clicked = code == BN_CLICKED;
    if (id == id_view && selection) {
        std::vector<ViewKind> kinds;
        const int dim = workspace_->device->mesh().dimension();
        for (const ViewKind k : available_views(*workspace_, run_b())) if (belongs(k, tab_, dim)) kinds.push_back(k);
        const int k = combo_index(child(id_view));
        if (k >= 0 && k < static_cast<int>(kinds.size())) s.kind = kinds[static_cast<std::size_t>(k)];
        view_defaults(s);
        sync_inspector();
    } else if (id == id_field && selection) {
        const auto fields = available_fields(*workspace_);
        const int k = combo_index(child(id_field));
        if (k >= 0 && k < static_cast<int>(fields.size())) s.field = fields[static_cast<std::size_t>(k)];
        if (s.kind == ViewKind::solid || s.kind == ViewKind::mesh) s.kind = ViewKind::field;
        s.scale.clear();
        s.range.reset();
        sync_inspector();
    } else if (id == id_point && selection) {
        if (&s == &device_view_) play(false);
        const int k = combo_index(child(id_point));
        s.point = k >= 0 ? std::optional<std::size_t>(static_cast<std::size_t>(k)) : std::nullopt;
    } else if (id == id_scale && selection) {
        const int k = combo_index(child(id_scale));
        s.scale = k <= 0 ? std::string() : std::string(scale_names[static_cast<std::size_t>(k)]);
    } else if (id == id_colormap && selection) {
        s.colormap = std::string(colormap_names[static_cast<std::size_t>(std::max(combo_index(child(id_colormap)), 0))]);
    } else if (clicked && (id == id_contours || id == id_edges || id == id_contacts || id == id_opaque || id == id_perspective)) {
        const bool on = Button_GetCheck(child(id)) == BST_CHECKED;
        (id == id_contours ? s.contours : id == id_edges ? s.edges : id == id_contacts ? s.contacts
                                                                      : id == id_opaque ? s.opaque : s.perspective) = on;
    } else if (clicked && id == id_enhanced) {
        // Unit 26f (KD5): this view, and the default of new views, remembered.
        s.enhanced = state_.enhanced_look = Button_GetCheck(child(id)) == BST_CHECKED;
    } else if (clicked && id == id_axis_box) {
        // Unit 26h (JQ3): this view, and the default of new views, remembered.
        s.axis_box = state_.axis_box = Button_GetCheck(child(id)) == BST_CHECKED;
    } else if ((id == id_arrows || id == id_lines) && selection) {
        // Unit 26c (YD1): arrows or streamlines of a vector field; on the field view the map shows its own
        // magnitude, as `nitcad plot` draws them without --field.
        const int k = std::max(combo_index(child(id)), 0);
        const std::string name = k == 0 ? std::string() : std::string(vector_names[static_cast<std::size_t>(k)]);
        (id == id_arrows ? s.vectors : s.streamlines) = name;
        if (!name.empty() && s.kind == ViewKind::field) {
            s.field = name;
            s.scale.clear();
            s.range.reset();
        }
        sync_inspector();
    } else if ((id == id_from || id == id_swept || id == id_sensed) && selection) {
        // A contact, the first entry none (nitcad plot's choice, or every contact's current).
        const auto contacts = workspace_->device->contacts();
        const int k = combo_index(child(id));
        const std::string name = k > 0 && k <= static_cast<int>(contacts.size()) ? contacts[static_cast<std::size_t>(k - 1)].name
                                                                                 : std::string();
        (id == id_from ? s.from : id == id_swept ? s.contact : s.sense) = name;
    } else if (id == id_frequency && selection) {
        const int k = combo_index(child(id_frequency));
        s.frequency = k >= 0 ? std::optional<std::size_t>(static_cast<std::size_t>(k)) : std::nullopt;
    } else if (id == id_axis && selection) {
        s.axis = std::clamp(combo_index(child(id_axis)), 0, 2);
    } else if (clicked && (id == id_exploded_field || id == id_log)) {
        (id == id_log ? s.log : s.exploded_field) = Button_GetCheck(child(id)) == BST_CHECKED;
    } else if (id == id_clip && clicked) {
        commit_box(id_clip_box);  // ZD1: the box's text takes effect when checked; unchecked, none
        return;
    } else if (id == id_slice_axis && selection) {
        commit_box(id_slice_at);  // the plane along the axis chosen, at the box's position
        return;
    } else if (id == id_range_auto && clicked) {
        // ZD3: Auto clears the range; unchecked, the scale drawn becomes the range, to be typed over.
        if (Button_GetCheck(child(id)) == BST_CHECKED) {
            s.range.reset();
        } else if (const auto* f = scene_ ? std::get_if<render::FieldScene>(&*scene_) : nullptr; f && f->colorbar) {
            s.range = std::pair{f->colorbar->scale.lo(), f->colorbar->scale.hi()};
        }
    } else if (code == EN_KILLFOCUS && is_box(id)) {
        // ZD1: focus left a box; its text takes effect when it changed.
        if (SendMessageW(child(id), EM_GETMODIFY, 0, 0) != 0) commit_box(id);
        return;
    } else {
        return;
    }
    refresh();
}

// --- the inspector's boxes (Unit 26d, ZD1, ZD3) ------------------------------------------------------

bool MainWindow::is_box(int id) noexcept {
    switch (id) {
    case id_slice_at:
    case id_clip_box:
    case id_levels:
    case id_cut:
    case id_arrow_count:
    case id_line_count:
    case id_seeds:
    case id_gap:
    case id_range_min:
    case id_range_max:
        return true;
    default:
        return false;
    }
}

std::string MainWindow::box_text(int id) const {
    const ViewSettings& s = tab_ == tab_curves ? curve_view_ : device_view_;
    const auto list = [](const std::vector<double>& v) {
        std::string t;
        for (const double x : v) t += (t.empty() ? "" : ",") + std::format("{}", x);
        return t;
    };
    switch (id) {
    case id_slice_at:
        return s.slices.empty() ? std::string() : std::format("{}", s.slices.front().at_um);
    case id_clip_box: {
        if (!s.clip) return {};
        const auto& c = *s.clip;
        return std::format("{},{},{}:{},{},{}", c.lo_um[0], c.lo_um[1], c.lo_um[2], c.hi_um[0], c.hi_um[1], c.hi_um[2]);
    }
    case id_levels:
        return list(s.levels);
    case id_cut:
        return cut_ ? list(cut_->from_um) + ":" + list(cut_->to_um) : std::string();
    case id_arrow_count:
        return s.arrows > 0 ? std::to_string(s.arrows) : std::string();
    case id_line_count:
        return s.lines > 0 ? std::to_string(s.lines) : std::string();
    case id_seeds:
        return s.seeds;
    case id_gap:
        return s.gap_um ? std::format("{}", *s.gap_um) : std::string();
    case id_range_min:
    case id_range_max: {
        // The range set, or (Auto) the scale the device view was drawn on.
        const bool top = id == id_range_max;
        if (s.range) return std::format("{}", top ? s.range->second : s.range->first);
        const auto* f = tab_ == tab_view && scene_ ? std::get_if<render::FieldScene>(&*scene_) : nullptr;
        if (f && f->colorbar) return std::format("{:.6g}", top ? f->colorbar->scale.hi() : f->colorbar->scale.lo());
        return {};
    }
    default:
        return {};
    }
}

void MainWindow::show_box(int id) {
    win32::set_text(child(id), box_text(id));
    SendMessageW(child(id), EM_SETMODIFY, FALSE, 0);
}

void MainWindow::sync_range(bool keep_typed) {
    for (const int id : {id_range_min, id_range_max}) {
        // A redraw (playback) leaves a box being typed in as it is.
        if (keep_typed && GetFocus() == child(id) && SendMessageW(child(id), EM_GETMODIFY, 0, 0) != 0) continue;
        show_box(id);
    }
    Button_SetCheck(child(id_range_auto), settings().range ? BST_UNCHECKED : BST_CHECKED);
}

void MainWindow::commit_box(int id) {
    if (!hwnd_ || !workspace_ || !is_box(id)) return;
    ViewSettings& s = settings();
    const auto blank = [](const std::string& t) { return t.find_first_not_of(" \t") == std::string::npos; };
    const std::string t = win32::text_of(child(id));
    const std::vector<double> v = numbers(t);
    std::string error;  // why the text is not taken; empty: taken
    switch (id) {
    case id_arrow_count:
    case id_line_count: {
        // The arrow and line counts; empty: nitcad plot's defaults.
        const bool arrows = id == id_arrow_count;
        const std::size_t most = arrows ? 100000 : 10000;
        std::size_t& to = arrows ? s.arrows : s.lines;
        if (blank(t)) {
            to = 0;
        } else if (v.size() == 1 && v[0] >= 1.0 && v[0] <= static_cast<double>(most) && v[0] == std::floor(v[0])) {
            to = static_cast<std::size_t>(v[0]);
        } else {
            error = std::format("{} are a count from 1 to {}.", arrows ? "Arrows" : "Streamlines", most);
        }
        break;
    }
    case id_seeds: {
        std::string seeds = t;
        while (!seeds.empty() && seeds.back() == ' ') seeds.pop_back();
        while (!seeds.empty() && seeds.front() == ' ') seeds.erase(seeds.begin());
        s.seeds = seeds;
        break;
    }
    case id_gap:
        // The exploded view's gap; empty: 12% of the extent.
        if (blank(t)) {
            s.gap_um.reset();
        } else if (v.size() == 1 && v[0] >= 0.0) {
            s.gap_um = v[0];
        } else {
            error = "The gap is a length in µm, at least 0.";
        }
        break;
    case id_slice_at: {
        // The plane along the axis chosen; no axis or no position: none.
        const int axis = combo_index(child(id_slice_axis)) - 1;
        if (axis < 0 || blank(t)) {
            s.slices.clear();
        } else if (v.size() == 1) {
            s.slices = {{axis, v[0]}};
        } else {
            error = "A slice plane is at one position in µm.";
        }
        SendMessageW(child(id_slice_axis), CB_SETCURSEL, s.slices.empty() ? 0 : static_cast<WPARAM>(s.slices.front().axis + 1), 0);
        break;
    }
    case id_clip_box: {
        std::string box = t;
        std::ranges::replace(box, ':', ',');
        if (const auto b = numbers(box); Button_GetCheck(child(id_clip)) != BST_CHECKED) {
            s.clip.reset();
        } else if (b.size() == 6) {
            ViewSettings::Clip c;
            for (std::size_t k = 0; k < 3; ++k) {
                c.lo_um[k] = b[k];
                c.hi_um[k] = b[k + 3];
            }
            s.clip = c;
        } else {
            error = "A clipping box is <x0>,<y0>,<z0>:<x1>,<y1>,<z1> in µm.";
        }
        Button_SetCheck(child(id_clip), s.clip ? BST_CHECKED : BST_UNCHECKED);
        sync_clip_sliders();  // YD6: the typed box on the sliders
        break;
    }
    case id_levels:
        if (blank(t)) {
            s.levels.clear();
        } else if (!v.empty()) {
            s.levels = v;
        } else {
            error = "Iso levels are numbers separated by commas.";
        }
        break;
    case id_cut: {
        // "<x0>,<y0>[,<z0>]:<x1>,<y1>[,<z1>]" in um; empty removes it. A new line is drawn by set_cut.
        const auto colon = t.find(':');
        const std::size_t dim = static_cast<std::size_t>(workspace_->device->mesh().dimension());
        const std::vector<double> a = colon == std::string::npos ? std::vector<double>{} : numbers(t.substr(0, colon));
        const std::vector<double> b = colon == std::string::npos ? std::vector<double>{} : numbers(t.substr(colon + 1));
        if (blank(t)) {
            if (cut_) {
                set_cut(std::nullopt);
                return;
            }
        } else if (a.size() == dim && b.size() == dim && a != b) {
            if (!cut_ || cut_->from_um != a || cut_->to_um != b) {
                set_cut(ViewSettings::Cut{a, b});
                return;
            }
        } else {
            error = "A cut line is <x0>,<y0>[,<z0>]:<x1>,<y1>[,<z1>] in µm, in the device's dimension.";
        }
        break;
    }
    case id_range_min:
    case id_range_max: {
        // ZD3: both boxes, in the colour bar's units as `nitcad plot --range` takes them; both empty: Auto.
        const std::string lo_t = win32::text_of(child(id_range_min)), hi_t = win32::text_of(child(id_range_max));
        const auto lo = numbers(lo_t), hi = numbers(hi_t);
        const auto field = render::find_field(s.field);
        const bool log = s.scale == "log" ||
                         (s.scale.empty() && field && render::field_info(*field).default_scale == render::ScaleKind::log10);
        if (blank(lo_t) && blank(hi_t)) {
            s.range.reset();
        } else if (lo.size() == 1 && hi.size() == 1 && std::isfinite(lo[0]) && std::isfinite(hi[0]) && lo[0] < hi[0] &&
                   (!log || lo[0] > 0.0)) {
            s.range = std::pair{lo[0], hi[0]};
        } else {
            error = log ? "The range is a minimum and a larger maximum, both above 0 on a log scale."
                        : "The range is a minimum and a larger maximum.";
        }
        break;
    }
    default:
        return;
    }
    if (!error.empty()) set_status(status_state, error);
    if (id == id_range_min || id == id_range_max) {
        sync_range();
    } else if (id != id_clip_box || s.clip || !error.empty()) {
        show_box(id);  // a clipping box's text stays while unchecked, to be checked again
    } else {
        SendMessageW(child(id), EM_SETMODIFY, FALSE, 0);
    }
    if (error.empty()) refresh();
}

// --- the scene tree --------------------------------------------------------------------------------

void MainWindow::fill_tree() {
    TreeView_DeleteAllItems(tree_);
    if (!workspace_) return;
    // Only the regions of a 2D or 3D device have a check box (YD4): checked while shown. The tree gives every
    // new item one, so the others' is taken away.
    // Look A: each item's icon (a region's swatch in its material's colour).
    const auto add = [&](HTREEITEM parent, const std::string& text, LPARAM data, int check = 0, int icon = tree_icon_view) {
        TVINSERTSTRUCTW t{};
        t.hParent = parent;
        t.hInsertAfter = TVI_LAST;
        t.item.mask = TVIF_TEXT | TVIF_PARAM | TVIF_IMAGE | TVIF_SELECTEDIMAGE;
        t.item.iImage = t.item.iSelectedImage = icon;
        const std::wstring w = win32::widen(text);
        t.item.pszText = const_cast<wchar_t*>(w.c_str());
        t.item.lParam = data;
        HTREEITEM item = TreeView_InsertItem(tree_, &t);
        TreeView_SetItemState(tree_, item, INDEXTOSTATEIMAGEMASK(check), TVIS_STATEIMAGEMASK);
        return item;
    };
    const auto& d = *workspace_->device;
    const bool boxes = d.mesh().dimension() >= 2;
    HTREEITEM device = add(TVI_ROOT, preview_ ? "Device (preview, not solved)" : "Device", tree_data(tree_none, 0), 0,
                           tree_icon_device);
    for (std::size_t k = 0; k < d.regions().size(); ++k) {
        const auto& r = d.regions()[k];
        const bool hidden = std::ranges::find(device_view_.hidden, r.name) != device_view_.hidden.end();
        add(device, r.name + (device::is_insulator(r) ? " (insulator)" : " (semiconductor)"),
            boxes ? tree_data(tree_region, static_cast<int>(k)) : tree_data(tree_none, 0), boxes ? (hidden ? 1 : 2) : 0,
            device::is_insulator(r) ? tree_icon_insulator : tree_icon_semiconductor);
    }
    for (const auto& c : d.contacts()) {
        const char* kind = c.kind == device::ContactKind::ohmic      ? "ohmic"
                           : c.kind == device::ContactKind::gate     ? "gate"
                           : c.kind == device::ContactKind::schottky ? "schottky"
                                                                     : "electrode";
        add(device, "contact " + c.name + " (" + kind + ")", tree_data(tree_none, 0), 0, tree_icon_contact);
    }
    TreeView_Expand(tree_, device, TVE_EXPAND);
    add(TVI_ROOT, "Mesh: " + std::to_string(d.mesh().node_count()) + " nodes", tree_data(tree_none, 0), 0, tree_icon_mesh);
    HTREEITEM views = add(TVI_ROOT, "Views", tree_data(tree_none, 0), 0, tree_icon_views);
    const auto kinds = available_views(*workspace_, run_b());
    for (std::size_t k = 0; k < kinds.size(); ++k) add(views, std::string(to_string(kinds[k])), tree_data(tree_view, static_cast<int>(k)));
    TreeView_Expand(tree_, views, TVE_EXPAND);
    HTREEITEM fields = add(TVI_ROOT, "Fields", tree_data(tree_none, 0), 0, tree_icon_fields);
    const auto names = available_fields(*workspace_);
    for (std::size_t k = 0; k < names.size(); ++k) {
        const auto f = render::find_field(names[k]);
        add(fields, std::string(f ? render::field_info(*f).title : names[k]), tree_data(tree_field, static_cast<int>(k)), 0,
            tree_icon_field);
    }
    TreeView_Expand(tree_, fields, TVE_EXPAND);
    if (workspace_->has_results()) {
        add(TVI_ROOT, "Results: " + std::to_string(workspace_->points.size()) + " points", tree_data(tree_none, 0), 0,
            tree_icon_results);
    }
}

void MainWindow::sync_region_boxes() {
    if (!workspace_) return;
    // The boxes as the user left them (a click or the space bar, posted): a region shown while checked; an item
    // the space bar gave a box loses it.
    const auto regions = workspace_->device->regions();
    for (HTREEITEM root = TreeView_GetRoot(tree_); root; root = TreeView_GetNextSibling(tree_, root)) {
        for (HTREEITEM item = root; item; item = item == root ? TreeView_GetChild(tree_, root) : TreeView_GetNextSibling(tree_, item)) {
            TVITEMW t{};
            t.mask = TVIF_PARAM | TVIF_STATE;
            t.stateMask = TVIS_STATEIMAGEMASK;
            t.hItem = item;
            TreeView_GetItem(tree_, &t);
            const int kind = static_cast<int>(t.lParam >> 16), index = static_cast<int>(t.lParam & 0xFFFF);
            const int check = static_cast<int>((t.state & TVIS_STATEIMAGEMASK) >> 12);
            if (kind != tree_region) {
                if (check != 0) TreeView_SetItemState(tree_, item, INDEXTOSTATEIMAGEMASK(0), TVIS_STATEIMAGEMASK);
            } else if (index < static_cast<int>(regions.size())) {
                set_region_shown(regions[static_cast<std::size_t>(index)].name, check != 1);
            }
        }
    }
}

void MainWindow::set_region_shown(const std::string& name, bool shown) {
    auto& hidden = device_view_.hidden;
    const bool was_hidden = std::ranges::find(hidden, name) != hidden.end();
    if (shown == !was_hidden) return;
    if (shown) {
        std::erase(hidden, name);
    } else {
        // In the device's order, as the tree lists them.
        hidden.push_back(name);
        const auto regions = workspace_->device->regions();
        const auto rank = [&](const std::string& n) { return std::ranges::find(regions, n, &device::Region::name) - regions.begin(); };
        std::ranges::sort(hidden, {}, rank);
    }
    set_status(status_state, (shown ? "Showing the region " : "Hiding the region ") + name + ".");
    if (tab_ == tab_view) refresh();
}

void MainWindow::on_tree_selection(LPARAM data) {
    if (!workspace_) return;
    const int kind = static_cast<int>(data >> 16), index = static_cast<int>(data & 0xFFFF);
    if (kind == tree_field) {
        const auto names = available_fields(*workspace_);
        if (index >= static_cast<int>(names.size())) return;
        device_view_.field = names[static_cast<std::size_t>(index)];
        if (device_view_.kind == ViewKind::solid || device_view_.kind == ViewKind::mesh) {
            device_view_.kind = workspace_->device->mesh().dimension() == 1 ? ViewKind::profile : ViewKind::field;
        }
        device_view_.scale.clear();
        device_view_.range.reset();
        select_tab(tab_view);
    } else if (kind == tree_view) {
        const auto kinds = available_views(*workspace_, run_b());
        if (index >= static_cast<int>(kinds.size())) return;
        const ViewKind k = kinds[static_cast<std::size_t>(index)];
        const int tab = belongs(k, tab_curves, workspace_->device->mesh().dimension()) ? tab_curves : tab_view;
        (tab == tab_curves ? curve_view_ : device_view_).kind = k;
        tab_ = tab;
        view_defaults(tab == tab_curves ? curve_view_ : device_view_);
        select_tab(tab);
    }
}

// --- documents, checks and runs --------------------------------------------------------------------

bool MainWindow::confirm_discard() {
    if (!document_.modified()) return true;
    const std::wstring q = L"Save the changes to " + win32::widen(document_.title()) + L"?";
    const int answer = MessageBoxW(hwnd_, q.c_str(), L"NiTCAD", MB_YESNOCANCEL | MB_ICONQUESTION);
    if (answer == IDCANCEL) return false;
    if (answer == IDYES) return save(false);
    return true;
}

bool MainWindow::save(bool as) {
    document_.set_text(deck_->text());
    if (as || document_.is_untitled()) {
        const fs::path p = win32::save_dialog(hwnd_, L"Decks (*.deck)\0*.deck\0All files\0*.*\0", L"Save the deck", L"deck");
        if (p.empty()) return false;
        if (auto r = document_.save_as(p); !r) {
            MessageBoxW(hwnd_, win32::widen(r.error().message).c_str(), L"NiTCAD", MB_ICONERROR);
            return false;
        }
        state_.add_recent(p);
    } else if (auto r = document_.save(); !r) {
        MessageBoxW(hwnd_, win32::widen(r.error().message).c_str(), L"NiTCAD", MB_ICONERROR);
        return false;
    }
    update_title();
    return true;
}

bool MainWindow::open_deck(const fs::path& path) {
    if (runner_.busy()) return false;
    auto d = DeckDocument::open(path);
    if (!d) {
        set_status(status_state, d.error().message);
        return false;
    }
    document_ = std::move(*d);
    deck_->set_text(document_.text());
    deck_->set_diagnostics({});
    state_.add_recent(path);
    update_title();
    workspace_.reset();
    ++workspace_generation_;  // Unit 26i (ND5): the metrics again
    run_b_.reset();
    previous_.reset();
    cut_.reset();
    prober_.reset();
    scene_.reset();
    fill_tree();
    pending_ = Pending::open;
    console_log_.begin_run();
    log_line(LogLevel::header, "Open " + document_.path().filename().string());
    runner_.start_check(document_.text(), document_.path());
    update_taskbar(false);  // ZD5: a check clears a failed run's state
    show_status(status_state, "Building the deck...");
    refresh();
    update_command_states();
    return true;
}

bool MainWindow::open_result(const fs::path& result, const fs::path& device, const fs::path& mesh) {
    auto w = workspace_from_result(result, device, mesh);
    if (!w) {
        set_status(status_state, w.error().message);
        MessageBoxW(hwnd_, win32::widen(w.error().message).c_str(), L"NiTCAD", MB_ICONERROR);
        return false;
    }
    state_.add_recent(result);
    cut_.reset();  // drawn on another device
    run_b_.reset();
    previous_.reset();
    workspace_.reset();
    ++workspace_generation_;  // Unit 26i (ND5): the metrics again
    set_workspace(std::move(*w));
    set_status(status_state, "Opened " + result.filename().string());
    select_tab(tab_view);
    return true;
}

bool MainWindow::compare_result(const fs::path& result, const fs::path& device, const fs::path& mesh) {
    if (!workspace_) {
        set_status(status_state, "Open a deck or a result file to compare with first.");
        return false;
    }
    auto b = workspace_from_result(result, device, mesh);
    if (!b) {
        set_status(status_state, b.error().message);
        MessageBoxW(hwnd_, win32::widen(b.error().message).c_str(), L"NiTCAD", MB_ICONERROR);
        return false;
    }
    const bool one_mesh = comparable(*workspace_, *b);
    set_run_b(std::move(*b));
    set_status(status_state, "Comparing with " + result.filename().string() +
                                 (one_mesh ? "" : " (another mesh: curves and figures only)"));
    return true;
}

// The seconds since a time point (the UI timings).
static double since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

void MainWindow::on_run_events() {
    const auto t0 = std::chrono::steady_clock::now();
    if (!runner_.process()) return;
    console_log_.follow(runner_);  // the text only: shown with the status, at most every 50 ms
    ++timings_.events;
    timings_.events_s += since(t0);
    const RunState s = runner_.state();
    if (s == RunState::running || s == RunState::stopping) {
        // The events of a fast solve arrive as fast as the window takes them; showing them costs more than
        // taking them, so the status and the monitor follow at most every 50 ms (and on the run timer when
        // the events pause), and the message queue stays free for input and paint (WG1).
        if (std::chrono::steady_clock::now() - run_ui_updated_ < std::chrono::milliseconds(50)) return;
        on_run_timer();
        return;
    }
    KillTimer(hwnd_, run_timer);
    const Pending pending = std::exchange(pending_, Pending::none);
    const auto t1 = std::chrono::steady_clock::now();
    run_page_->update(runner_, renderer(), true);  // before the results are taken: how the run ended
    console_log_.finish(runner_);
    flush_console();
    update_points(true);  // before the results are taken: their currents
    deck_->set_diagnostics(runner_.diagnostics());
    timings_.finish_s = since(t1);
    const auto t2 = std::chrono::steady_clock::now();
    on_run_finished(pending, s);
    timings_.results_s = since(t2);
}

// The run's progress in the status bar and the run monitor: on its events (at most every 50 ms) and on
// a timer while it runs.
void MainWindow::on_run_timer() {
    const RunState s = runner_.state();
    if (s != RunState::running && s != RunState::stopping) {
        KillTimer(hwnd_, run_timer);
        return;
    }
    run_ui_updated_ = std::chrono::steady_clock::now();
    {
        const auto& rows = runner_.points();
        std::string text = s == RunState::stopping ? "Stopping..." : "Running";
        if (!rows.empty()) {
            const PointRow& r = rows.back();
            if (r.phase == solve::Phase::equilibrium) text += ": equilibrium";
            else if (r.phase == solve::Phase::transient) text += std::format(": t = {:.3g} s", r.time_s);
            else if (r.point_count > 0) text += std::format(": point {} of {}", r.point + 1, r.point_count);
            else text += std::format(": point {}", r.point + 1);
            text += std::format(", iteration {}, update {:.2g}", r.iterations, r.update);
        }
        const auto t0 = std::chrono::steady_clock::now();
        show_status(status_state, text);
        auto t = std::chrono::steady_clock::now();
        flush_console();
        timings_.console_s += since(t);
        t = std::chrono::steady_clock::now();
        update_points();
        timings_.points_s += since(t);
        t = std::chrono::steady_clock::now();
        if (tab_ == tab_run) run_page_->update(runner_, renderer());  // hidden, it catches up when shown
        timings_.monitor_s += since(t);
        update_command_states();
        update_taskbar(pending_ == Pending::run);
        const double spent = since(t0);
        ++timings_.updates;
        timings_.update_s += spent;
        timings_.update_max_s = std::max(timings_.update_max_s, spent);
    }
}

void MainWindow::on_run_finished(Pending pending, RunState s) {
    if (s == RunState::failed) {
        const std::string why = runner_.error() ? runner_.error()->message : std::string("failed");
        set_status(status_state, (pending == Pending::run ? "The run failed: " : "The deck does not build: ") + why,
                   LogLevel::error);
        if (!runner_.diagnostics().empty()) {
            select_tab(tab_deck);
            deck_->go_to_line(runner_.diagnostics().front().line);
        }
    } else if (pending == Pending::open && s == RunState::idle) {
        if (auto plan = runner_.take_plan()) {
            auto w = workspace_from_plan(std::move(*plan), document_.path());
            if (w) {
                set_workspace(std::move(*w));
                set_status(status_state, workspace_->has_results() ? "Opened the deck and its results"
                                                                   : "Opened the deck (no results yet: Run, F5)");
                select_tab(tab_view);
            } else {
                set_status(status_state, w.error().message);
            }
        }
    } else if (pending == Pending::check && s == RunState::idle) {
        set_status(status_state, "The deck builds.", LogLevel::success);
    } else if (pending == Pending::run && s == RunState::done) {
        const auto stopped = runner_.stopped_early();
        const double seconds = runner_.elapsed_seconds();
        if (auto taken = runner_.take_results()) {
            auto w = workspace_from_run(std::move(taken->plan), std::move(taken->report));
            if (w) {
                w->source = document_.path();
                set_workspace(std::move(*w));
                select_tab(tab_view);
            } else {
                set_status(status_state, w.error().message);
            }
        }
        set_status(status_state, stopped ? std::format("Stopped early after {:.1f} s: {}", seconds, stopped->message)
                                         : std::format("Done in {:.1f} s", seconds),
                   stopped ? LogLevel::warning : LogLevel::success);
    }
    update_command_states();
    update_taskbar(pending == Pending::run);
    // ZD5: a run that ends while the window is in the background flashes its taskbar button until activated.
    if (pending == Pending::run && (s == RunState::done || s == RunState::failed) && IsWindowVisible(hwnd_) &&
        GetForegroundWindow() != hwnd_) {
        FLASHWINFO f{sizeof(f), hwnd_, FLASHW_ALL | FLASHW_TIMERNOFG, 0, 0};
        FlashWindowEx(&f);
    }
}

// --- the shell (Unit 26d, ZD4-ZD6) -------------------------------------------------------------------

bool MainWindow::open_path(const fs::path& path) {
    const OpenTarget t = open_target(path);
    return t.deck ? open_deck(t.path) : open_result(t.path, t.device, t.mesh);
}

void MainWindow::on_drop(HDROP drop) {
    std::vector<fs::path> files;
    const UINT n = DragQueryFileW(drop, 0xFFFFFFFF, nullptr, 0);
    for (UINT k = 0; k < n; ++k) {
        const UINT length = DragQueryFileW(drop, k, nullptr, 0);
        std::wstring name(static_cast<std::size_t>(length) + 1, L'\0');
        DragQueryFileW(drop, k, name.data(), length + 1);
        name.resize(length);
        files.emplace_back(std::move(name));
    }
    DragFinish(drop);
    const auto file = dropped_file(files);
    if (!file) {
        set_status(status_state, file.error());
        return;
    }
    if (runner_.busy()) {
        set_status(status_state, "A run or check is in progress: drop the file when it ends.");
        return;
    }
    if (confirm_discard()) open_path(*file);
}

void MainWindow::update_taskbar(bool run) {
    update_status_progress();
    const RunState s = runner_.state();
    const double fraction = s == RunState::running || s == RunState::stopping
                                ? progress_fraction(runner_.points(), runner_.planned_points())
                                : -1.0;
    const TaskbarProgress next = app::taskbar_progress(s, run, fraction, taskbar_);
    if (next == taskbar_) return;
    taskbar_ = next;
    show_taskbar();
}

void MainWindow::show_taskbar() {
    if (!taskbar_list_) return;  // no taskbar button (a hidden window), or not made yet
    const TBPFLAG flag = taskbar_.flag == TaskbarFlag::normal          ? TBPF_NORMAL
                         : taskbar_.flag == TaskbarFlag::indeterminate ? TBPF_INDETERMINATE
                         : taskbar_.flag == TaskbarFlag::error         ? TBPF_ERROR
                                                                       : TBPF_NOPROGRESS;
    taskbar_list_->SetProgressState(hwnd_, flag);
    if (taskbar_.flag == TaskbarFlag::normal || taskbar_.flag == TaskbarFlag::error) {
        taskbar_list_->SetProgressValue(hwnd_, static_cast<ULONGLONG>(std::lround(taskbar_.fraction * 1000.0)), 1000);
    }
}

void MainWindow::apply_title_theme() {
    DWORD value = 1, size = sizeof(value);
    const bool read = RegGetValueW(HKEY_CURRENT_USER, L"Software\\Microsoft\\Windows\\CurrentVersion\\Themes\\Personalize",
                                   L"AppsUseLightTheme", RRF_RT_REG_DWORD, nullptr, &value, &size) == ERROR_SUCCESS;
    dark_title_ = dark_title_bar(read ? std::optional<unsigned long>(value) : std::nullopt);
    const BOOL dark = dark_title_ ? TRUE : FALSE;
    // Fails where the attribute is not known (before Windows 10 20H1): the title bar stays light.
    (void)DwmSetWindowAttribute(hwnd_, DWMWA_USE_IMMERSIVE_DARK_MODE, &dark, sizeof(dark));
    ++title_theme_updates_;
}

void MainWindow::deck_changed() {
    document_.set_text(deck_->text());
    update_title();
    // Checked a pause after the last edit (XD1).
    ++edits_;
    SetTimer(hwnd_, live_timer, live_delay_ms, nullptr);
}

void MainWindow::on_live_check() {
    if (!live_.process()) return;
    // Edited since: the next check (its timer is running) will tell.
    if (live_edits_ != edits_) return;
    const std::string text = deck_->text();
    std::vector<Diagnostic> shown = live_.result().diagnostics;
    for (Diagnostic& d : shown) {
        if (const auto s = suggestion(text, d)) d.message += " (" + *s + ")";
    }
    deck_->set_diagnostics(shown);
    if (!shown.empty()) {
        set_status(status_state, shown.front().message);
    } else if (auto plan = live_.take_plan()) {
        // The structure preview (XD4): results stay while the deck keeps their mesh.
        if (workspace_ && workspace_->has_results() && same_mesh(*workspace_, plan->device.mesh())) {
            set_status(status_state, "The deck builds (the results shown are of its last run).");
        } else if (auto w = workspace_from_structure(std::move(*plan), document_.path())) {
            show_structure(std::move(*w));
            set_status(status_state, "The deck builds: its structure is shown, not solved (Run, F5).");
        } else {
            set_status(status_state, w.error().message);
        }
    }
    update_command_states();
}

void MainWindow::show_structure(Workspace w) {
    // While editing, the view and camera stay as they were if the previous view was a preview of the
    // same dimension.
    const bool keep = workspace_ && preview_ &&
                      workspace_->device->mesh().dimension() == w.device->mesh().dimension();
    set_workspace(std::move(w), keep);
    preview_ = true;
    fill_tree();  // marked as a preview
}

void MainWindow::go_to_section(std::string_view kind) {
    const auto s = find_section(deck_->text(), kind);
    select_tab(tab_deck);
    if (s) {
        deck_->go_to_line(s->line);
        set_status(status_state, "[" + std::string(kind) + "] at line " + std::to_string(s->line));
    } else {
        set_status(status_state, "The deck has no [" + std::string(kind) + "] section.");
    }
}

void MainWindow::insert_key(std::size_t index) {
    if (index >= key_items_.size()) return;
    const KeyItem& k = key_items_[index];
    const auto in = app::insert_key(deck_->text(), k.kind, k.key);
    if (!in) {
        set_status(status_state, "The deck has no [" + k.kind + "] section to insert into.");
        return;
    }
    deck_->set_text(in->text);
    deck_changed();
    select_tab(tab_deck);
    deck_->go_to_line(in->line);
    set_status(status_state, "Inserted " + std::string(k.key.key) + " = " + std::string(k.key.value) + " at line " +
                                 std::to_string(in->line));
}

void MainWindow::forward_edit(UINT message) {
    HWND focus = GetFocus();
    if (!focus) return;
    if (message == EM_SETSEL) SendMessageW(focus, EM_SETSEL, 0, -1);
    else SendMessageW(focus, message, 0, 0);
}

void MainWindow::focus_next_pane() {
    HWND focus = GetFocus();
    const HWND page = tab_ == tab_deck      ? deck_->editor()
                      : tab_ == tab_start   ? start_->hwnd()
                      : tab_ == tab_run     ? run_page_->hwnd()
                      : tab_ == tab_results ? results_->hwnd()
                                            : viewport_->hwnd();
    const HWND order[] = {tree_, page, child(id_view)};
    std::size_t k = 0;
    for (std::size_t i = 0; i < 3; ++i) {
        if (focus == order[i] || IsChild(order[i], focus) || (i == 2 && IsChild(inspector_, focus))) k = i + 1;
    }
    // A page that is a panel gives the focus to its first control.
    HWND target = order[k % 3];
    if (HWND first = GetNextDlgTabItem(target, nullptr, FALSE); first && IsChild(target, first)) target = first;
    SetFocus(target);
}

void MainWindow::execute(Command c) {
    const auto id = static_cast<int>(c);
    if (id >= static_cast<int>(Command::recent_first) && id < static_cast<int>(Command::recent_first) + 16) {
        const auto k = static_cast<std::size_t>(id - static_cast<int>(Command::recent_first));
        if (k < state_.recent.size() && confirm_discard()) open_path(fs::path(state_.recent[k]));
        return;
    }
    if (id >= static_cast<int>(Command::key_first)) {
        insert_key(static_cast<std::size_t>(id - static_cast<int>(Command::key_first)));
        return;
    }
    if (id >= static_cast<int>(Command::section_first)) {
        const auto k = static_cast<std::size_t>(id - static_cast<int>(Command::section_first));
        if (k < section_items_.size()) go_to_section(section_items_[k].kind);
        return;
    }
    ViewSettings& s = settings();
    switch (c) {
    case Command::new_deck:
        if (!confirm_discard()) return;
        document_ = DeckDocument::untitled();
        deck_->set_text("");
        deck_->set_diagnostics({});
        update_title();
        select_tab(tab_deck);
        return;
    case Command::open_deck: {
        if (!confirm_discard()) return;
        const fs::path p = win32::open_dialog(hwnd_, L"Decks (*.deck)\0*.deck\0All files\0*.*\0", L"Open a deck");
        if (!p.empty()) open_deck(p);
        return;
    }
    case Command::open_result: {
        const fs::path r = win32::open_dialog(
            hwnd_, L"Results (*.nitsweep;*.nittran;*.nitac)\0*.nitsweep;*.nittran;*.nitac\0All files\0*.*\0",
            L"Open a result file");
        if (r.empty()) return;
        // Its device file: the one of the same name beside it, else asked for; a mesh file likewise.
        fs::path d = fs::path(r).replace_extension(".nitdev");
        if (!fs::exists(d)) d = win32::open_dialog(hwnd_, L"Device files (*.nitdev)\0*.nitdev\0", L"The result's device file");
        if (d.empty()) return;
        fs::path m = fs::path(r).replace_extension(".nitmesh");
        if (!fs::exists(m)) m.clear();
        open_result(r, d, m);
        return;
    }
    case Command::compare_result: {
        if (!workspace_) return;
        const fs::path r = win32::open_dialog(
            hwnd_, L"Results (*.nitsweep;*.nittran;*.nitac)\0*.nitsweep;*.nittran;*.nitac\0All files\0*.*\0",
            L"Compare with a result file");
        if (r.empty()) return;
        // Its device and mesh files as Open result finds them.
        fs::path d = fs::path(r).replace_extension(".nitdev");
        if (!fs::exists(d)) d = win32::open_dialog(hwnd_, L"Device files (*.nitdev)\0*.nitdev\0", L"The result's device file");
        if (d.empty()) return;
        fs::path m = fs::path(r).replace_extension(".nitmesh");
        if (!fs::exists(m)) m.clear();
        compare_result(r, d, m);
        return;
    }
    case Command::save: save(false); return;
    case Command::save_as: save(true); return;
    case Command::exit: PostMessageW(hwnd_, WM_CLOSE, 0, 0); return;
    case Command::undo: forward_edit(EM_UNDO); return;
    case Command::cut: forward_edit(WM_CUT); return;
    case Command::copy: forward_edit(WM_COPY); return;
    case Command::paste: forward_edit(WM_PASTE); return;
    case Command::select_all: forward_edit(EM_SETSEL); return;
    case Command::find_error:
        if (!deck_->diagnostics().empty()) {
            select_tab(tab_deck);
            deck_->go_to_line(deck_->diagnostics().front().line);
        }
        return;
    case Command::complete:
        select_tab(tab_deck);
        SetFocus(deck_->editor());
        if (!deck_->show_completion(true)) set_status(status_state, "Nothing to complete here.");
        return;
    case Command::check:
        if (runner_.busy()) return;
        document_.set_text(deck_->text());
        pending_ = Pending::check;
        console_log_.begin_run();
        log_line(LogLevel::header, "Check " + (document_.is_untitled() ? std::string("Untitled") : document_.path().filename().string()));
        runner_.start_check(document_.text(), document_.path());
        update_taskbar(false);
        show_status(status_state, "Checking...");
        update_command_states();
        return;
    case Command::run:
        if (runner_.busy()) return;
        document_.set_text(deck_->text());
        if (document_.is_untitled() && !save(true)) return;  // a run writes its outputs beside the deck
        select_tab(tab_run);  // the run monitor follows it (shown before the run starts)
        pending_ = Pending::run;
        console_log_.begin_run();
        timings_ = {};
        set_console_collapsed(false);  // Unit 26j: the console opens for a run
        log_line(LogLevel::header, "Run " + document_.path().filename().string());
        runner_.start_run(document_.text(), document_.path());
        update_points();  // the last run's rows cleared
        update_taskbar(true);
        show_status(status_state, "Running...");
        SetTimer(hwnd_, run_timer, 100, nullptr);  // the progress while the events pause
        run_page_->update(runner_, renderer(), true);  // the last run's rows cleared
        update_command_states();
        return;
    case Command::stop:
        runner_.request_stop();
        set_status(status_state, "Stopping at the solver's next safe point...");
        update_command_states();
        return;
    case Command::compare_previous: {
        if (!workspace_ || !previous_) return;
        std::optional<Workspace> b = std::move(previous_);
        previous_.reset();
        set_run_b(std::move(b));
        set_status(status_state, comparable(*workspace_, *run_b_) ? "Comparing with the previous run"
                                                                  : "Comparing with the previous run (another mesh: curves and figures only)");
        return;
    }
    case Command::stop_comparing:
        // The previous run stays to compare with again.
        if (run_b_ && !previous_ && workspace_ && run_b_->source == workspace_->source) previous_ = std::move(run_b_);
        set_run_b(std::nullopt);
        set_status(status_state, "Not comparing");
        return;
    case Command::fit: viewport_fit(); return;
    case Command::view_front: show_named_view("Front"); return;
    case Command::view_back: show_named_view("Back"); return;
    case Command::view_left: show_named_view("Left"); return;
    case Command::view_right: show_named_view("Right"); return;
    case Command::view_top: show_named_view("Top"); return;
    case Command::view_bottom: show_named_view("Bottom"); return;
    case Command::view_isometric: show_named_view("Isometric"); return;
    case Command::orbit_left: viewport_orbit(-15.0f, 0.0f); return;
    case Command::orbit_right: viewport_orbit(15.0f, 0.0f); return;
    case Command::orbit_up: viewport_orbit(0.0f, 10.0f); return;
    case Command::orbit_down: viewport_orbit(0.0f, -10.0f); return;
    case Command::zoom_in: viewport_zoom(1.25f); return;
    case Command::zoom_out: viewport_zoom(0.8f); return;
    case Command::pan_left: viewport_pan(-0.1f, 0.0f); return;
    case Command::pan_right: viewport_pan(0.1f, 0.0f); return;
    case Command::pan_up: viewport_pan(0.0f, 0.1f); return;
    case Command::pan_down: viewport_pan(0.0f, -0.1f); return;
    case Command::perspective:
        s.perspective = !s.perspective;
        sync_inspector();
        refresh();
        return;
    case Command::probe_mode:
        probe_mode_ = !probe_mode_;
        if (probe_mode_) measure_mode_ = cut_mode_ = false;
        update_command_states();
        return;
    case Command::measure_mode:
        measure_mode_ = !measure_mode_;
        if (measure_mode_) probe_mode_ = cut_mode_ = false;
        measure_from_.reset();
        update_command_states();
        return;
    case Command::cut_line_mode:
        cut_mode_ = !cut_mode_;
        if (cut_mode_) probe_mode_ = measure_mode_ = false;
        cut_from_.reset();
        if (cut_mode_) {
            if (tab_ != tab_view) select_tab(tab_view);
            set_status(status_state, "Draw cut line: click the line's start, then its end.");
        }
        update_command_states();
        return;
    case Command::mesh_quality:
        // XD7: the minimum angle map in the device view.
        if (!workspace_ || workspace_->device->mesh().dimension() < 2) return;
        device_view_.kind = ViewKind::quality;
        device_view_.scale.clear();
        device_view_.range.reset();
        select_tab(tab_view);
        return;
    case Command::play_pause: play(!playing_); return;
    case Command::snapshot_first:
    case Command::snapshot_previous:
    case Command::snapshot_next:
    case Command::snapshot_last: {
        // YD5: a step pauses and lands on the snapshot asked for, in the device view.
        if (!can_play()) return;
        play(false);
        const std::size_t n = workspace_->points.size();
        const std::size_t k = std::min(device_view_.point.value_or(n - 1), n - 1);
        if (tab_ != tab_view) select_tab(tab_view);
        show_snapshot(c == Command::snapshot_first      ? 0
                      : c == Command::snapshot_previous ? (k == 0 ? 0 : k - 1)
                      : c == Command::snapshot_next     ? std::min(k + 1, n - 1)
                                                        : n - 1);
        return;
    }
    case Command::show_regions:
        // YD4: every region shown again, the tree's boxes checked.
        if (!workspace_ || device_view_.hidden.empty()) return;
        device_view_.hidden.clear();
        fill_tree();
        set_status(status_state, "Showing every region.");
        if (tab_ == tab_view) refresh();
        return;
    case Command::hover_readout:
        // HD3: on by default, saved in the app state.
        state_.hover_readout = !state_.hover_readout;
        if (!state_.hover_readout) hide_readout();
        update_command_states();
        set_status(status_state, state_.hover_readout ? "Hover readout on." : "Hover readout off.");
        return;
    case Command::next_tab: select_tab((tab_ + 1) % tab_count); return;
    case Command::previous_tab: select_tab((tab_ + tab_count - 1) % tab_count); return;
    case Command::next_pane: focus_next_pane(); return;
    case Command::tab_view: select_tab(tab_view); return;
    case Command::tab_curves: select_tab(tab_curves); return;
    case Command::tab_deck: select_tab(tab_deck); return;
    case Command::tab_start: select_tab(tab_start); return;
    case Command::tab_run: select_tab(tab_run); return;
    case Command::tab_results: select_tab(tab_results); return;
    default: return;
    }
}

// --- the view's host -------------------------------------------------------------------------------

bool MainWindow::viewport_is_3d() const {
    if (!scene_) return false;
    const auto* f = std::get_if<render::FieldScene>(&*scene_);
    return f && f->dimension == 3;
}

void MainWindow::viewport_orbit(float azimuth_deg, float elevation_deg) {
    if (!viewport_is_3d()) return;
    stop_move();  // JD5: a drag or a key ends a move where it is
    orbit(settings(), azimuth_deg, elevation_deg);
    present_camera();
}

void MainWindow::viewport_pan(float dx, float dy) {
    if (!scene_ || !std::holds_alternative<render::FieldScene>(*scene_)) return;
    stop_move();
    pan_by(settings(), dx, dy);
    present_camera();
}

void MainWindow::viewport_zoom(float factor) {
    if (!scene_ || !std::holds_alternative<render::FieldScene>(*scene_)) return;
    stop_move();
    zoom_by(settings(), factor);
    present_camera();
}

void MainWindow::viewport_fit() {
    // JD5: smoothly, the angles kept.
    CameraPose to = camera_pose(settings());
    to.zoom = 1.0f;
    to.pan_x = to.pan_y = 0.0f;
    if (!scene_ || !std::holds_alternative<render::FieldScene>(*scene_)) {
        fit(settings());
        return present_camera();
    }
    move_camera(to);
}

void MainWindow::probe_at(float x, float y) {
    if (!pickable()) return;
    const auto* f = std::get_if<render::FieldScene>(&*scene_);
    if (!f) return;
    const auto hit = prober_->pick(*f, viewport_->width(), viewport_->height(), x, y, values_);
    std::string text;
    if (!hit) {
        text = "No node within reach.";
    } else if (measure_mode_) {
        if (!measure_from_) {
            measure_from_ = hit;
            text = std::format("Measure from node {}\n({}, {}, {}) µm\nClick the second point.", hit->node,
                               format_um(hit->position_um[0]), format_um(hit->position_um[1]), format_um(hit->position_um[2]));
        } else {
            text = std::format("Distance {} µm\nfrom node {} to node {}", format_um(distance_um(*measure_from_, *hit)),
                               measure_from_->node, hit->node);
            measure_from_.reset();
        }
    } else {
        probe_ = hit;
        text = std::format("Probe: node {}\nx {} µm\ny {} µm", hit->node, format_um(hit->position_um[0]),
                           format_um(hit->position_um[1]));
        if (workspace_->device->mesh().dimension() == 3) text += std::format("\nz {} µm", format_um(hit->position_um[2]));
        text += "\n" + settings().field + " = " + number(hit->value);
    }
    win32::set_text(child(id_probe_text), text);
    set_status(status_state, text.substr(0, text.find('\n')));
}

bool MainWindow::viewport_over_cube(float x, float y) const {
    if (!scene_ || !viewport_) return false;
    const auto* f = std::get_if<render::FieldScene>(&*scene_);
    return f && render::view_cube_hit(*f, viewport_->width(), viewport_->height(), x, y).has_value();
}

void MainWindow::snap_to(const render::CameraPreset& preset) {
    if (!viewport_is_3d()) return;
    move_camera({preset.azimuth_deg, preset.elevation_deg, 1.0f, 0.0f, 0.0f});
    set_status(status_state, "View: " + preset.name);
}

void MainWindow::show_named_view(std::string_view name) {
    for (const render::CameraPreset& p : render::named_views()) {
        if (p.name == name) return snap_to(p);
    }
}

void MainWindow::viewport_named_view(int index) {
    const auto& views = render::named_views();
    if (index >= 0 && static_cast<std::size_t>(index) < views.size()) snap_to(views[static_cast<std::size_t>(index)]);
}

void MainWindow::move_camera(const CameraPose& to) {
    stop_move();
    ViewSettings& s = settings();
    const CameraPose from = camera_pose(s);
    BOOL shown = TRUE;
    if (!SystemParametersInfoW(SPI_GETCLIENTAREAANIMATION, 0, &shown, 0)) shown = TRUE;
    if (!smooth_moves(animate_, shown != FALSE, force_animate_) || from == to || !hwnd_) {
        set_camera_pose(s, to);
        present_camera();
        return;
    }
    move_ = CameraMove{from, to, &s, std::chrono::steady_clock::now()};
    SetTimer(hwnd_, camera_timer, camera_tick_ms, nullptr);
}

void MainWindow::stop_move() {
    if (!move_) return;
    move_.reset();
    if (hwnd_) KillTimer(hwnd_, camera_timer);
}

void MainWindow::on_camera_timer() {
    if (!move_ || move_->settings != &settings() || !scene_) return stop_move();
    const float ms = std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - move_->start).count();
    const float t = std::min(ms / camera_move_ms, 1.0f);
    set_camera_pose(*move_->settings, camera_between(move_->from, move_->to, t));
    present_camera();
    if (t >= 1.0f) stop_move();
}

void MainWindow::viewport_click(float x, float y) {
    // Unit 26h (JD3): a click on the view cube snaps the camera, whatever the mode.
    if (scene_) {
        if (const auto* f = std::get_if<render::FieldScene>(&*scene_)) {
            if (const auto p = render::view_cube_hit(*f, viewport_->width(), viewport_->height(), x, y)) {
                snap_to(*p);
                return;
            }
        }
    }
    if (cut_mode_) cut_at(x, y);
    else if (probe_mode_ || measure_mode_) probe_at(x, y);
}

void MainWindow::cut_at(float x, float y) {
    if (!pickable() || !workspace_) return;
    const auto* f = std::get_if<render::FieldScene>(&*scene_);
    if (!f) return;
    const auto p = prober_->point_at(*f, viewport_->width(), viewport_->height(), x, y);
    if (!p) {
        set_status(status_state, "No point of the device there.");
        return;
    }
    // To the nanometre, as typed.
    std::vector<double> at;
    for (int k = 0; k < workspace_->device->mesh().dimension(); ++k) {
        at.push_back(std::round((*p)[static_cast<std::size_t>(k)] * 1e3) / 1e3);
    }
    const auto point = [](const std::vector<double>& v) {
        std::string t;
        for (const double c : v) t += (t.empty() ? "" : ", ") + format_um(c);
        return "(" + t + ") µm";
    };
    if (!cut_from_) {
        cut_from_ = at;
        const std::string text = "Cut line from " + point(at) + "\nClick its end.";
        win32::set_text(child(id_probe_text), text);
        set_status(status_state, text.substr(0, text.find('\n')));
        return;
    }
    ViewSettings::Cut c{*cut_from_, at};
    cut_from_.reset();
    if (c.from_um == c.to_um) {
        set_status(status_state, "The cut line's ends are the same point: click its start again.");
        return;
    }
    win32::set_text(child(id_probe_text), "Cut line from " + point(c.from_um) + "\nto " + point(c.to_um));
    set_cut(std::move(c));
}

void MainWindow::set_cut(std::optional<ViewSettings::Cut> cut) {
    cut_ = std::move(cut);
    cut_from_.reset();
    if (!workspace_) return;
    if (cut_ && workspace_->device->mesh().dimension() > 1) {
        // Plotted along in the Curves tab: the band diagram once solved, else the doping.
        if (workspace_->has_results()) {
            curve_view_.kind = ViewKind::bands;
        } else {
            curve_view_.kind = ViewKind::profile;
            curve_view_.field = "doping";
        }
        curve_view_.point = device_view_.point;
        set_status(status_state, std::string("Cut line drawn: its ") +
                                     (workspace_->has_results() ? "band diagram" : "doping") + " is in the Curves tab.");
        select_tab(tab_curves);  // draws it
        return;
    }
    sync_inspector();
    refresh();
}

void MainWindow::viewport_hover(float x, float y) {
    if (!scene_) return hide_readout();
    const int w = viewport_->width(), h = viewport_->height();
    if (const auto* p = std::get_if<render::PlotScene>(&*scene_)) {
        // HD2: the nearest data point of any series.
        const auto hit = pick_curve(*p, w, h, x, y);
        if (!hit) {
            set_status(status_cursor, {});
            return hide_readout();
        }
        std::string text = curve_text(*p, *hit), line = text;
        std::ranges::replace(line, '\n', ' ');
        set_status(status_cursor, line);
        show_readout(text, hit->screen_x, hit->screen_y);
        return;
    }
    if (!pickable()) return hide_readout();
    const auto& f = std::get<render::FieldScene>(*scene_);
    // HD1: the node Probe would pick, with the field's value.
    const auto hit = prober_->pick(f, w, h, x, y, values_, 8.0f);
    if (!hit) {
        set_status(status_cursor, {});
        return hide_readout();
    }
    std::string text = hover_text(values_title_, *hit, workspace_->device->mesh().dimension());
    if (values_.empty()) text.erase(0, text.find('\n') + 1);  // no field shown: the position and node
    const std::string position = std::format("x {}  y {}  z {} µm", format_um(hit->position_um[0]),
                                             format_um(hit->position_um[1]), format_um(hit->position_um[2]));
    set_status(status_cursor, values_.empty() ? position : text.substr(0, text.find('\n')) + "     " + position);
    const render::ScreenPoint s = render::screen_position(
        f, w, h, {static_cast<float>(hit->position_um[0]), static_cast<float>(hit->position_um[1]),
                  static_cast<float>(hit->position_um[2])});
    show_readout(text, s.x, s.y);
}

void MainWindow::viewport_leave() { hide_readout(); }

void MainWindow::show_readout(const std::string& text, float x, float y) {
    if (!state_.hover_readout) return hide_readout();
    if (!readout_) {
        readout_ = CreateWindowExW(WS_EX_TOPMOST | WS_EX_TRANSPARENT, TOOLTIPS_CLASSW, nullptr,
                                   WS_POPUP | TTS_NOPREFIX | TTS_ALWAYSTIP, CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT,
                                   CW_USEDEFAULT, hwnd_, nullptr, GetModuleHandleW(nullptr), nullptr);
        if (readout_) {
            TTTOOLINFOW ti{sizeof(ti)};
            ti.uFlags = TTF_TRACK | TTF_ABSOLUTE;
            ti.hwnd = hwnd_;
            ti.uId = 1;
            ti.lpszText = const_cast<wchar_t*>(L"");
            SendMessageW(readout_, TTM_ADDTOOLW, 0, reinterpret_cast<LPARAM>(&ti));
            SendMessageW(readout_, TTM_SETMAXTIPWIDTH, 0, 600);  // lets "\n" break lines
        }
    }
    readout_text_ = text;
    readout_shown_ = true;
    // A hidden window (the tests) keeps the readout's state without putting a tooltip on the desktop.
    if (!readout_ || !IsWindowVisible(hwnd_)) return;
    const std::wstring wide = win32::widen(text);
    TTTOOLINFOW ti{sizeof(ti)};
    ti.hwnd = hwnd_;
    ti.uId = 1;
    ti.lpszText = const_cast<wchar_t*>(wide.c_str());
    SendMessageW(readout_, TTM_UPDATETIPTEXTW, 0, reinterpret_cast<LPARAM>(&ti));
    // 16 px right of and below the picked point, kept on the view's monitor.
    POINT at{static_cast<LONG>(std::lround(x)) + 16, static_cast<LONG>(std::lround(y)) + 16};
    ClientToScreen(viewport_->hwnd(), &at);
    SIZE size{};
    if (const LRESULT r = SendMessageW(readout_, TTM_GETBUBBLESIZE, 0, reinterpret_cast<LPARAM>(&ti)); r) {
        size = {LOWORD(r), HIWORD(r)};
    }
    MONITORINFO mi{sizeof(mi)};
    if (GetMonitorInfoW(MonitorFromPoint(at, MONITOR_DEFAULTTONEAREST), &mi)) {
        if (at.x + size.cx > mi.rcWork.right) at.x = std::max(mi.rcWork.left, at.x - 32 - size.cx);
        if (at.y + size.cy > mi.rcWork.bottom) at.y = std::max(mi.rcWork.top, at.y - 32 - size.cy);
    }
    SendMessageW(readout_, TTM_TRACKPOSITION, 0, MAKELPARAM(at.x, at.y));
    SendMessageW(readout_, TTM_TRACKACTIVATE, TRUE, reinterpret_cast<LPARAM>(&ti));
}

void MainWindow::hide_readout() {
    readout_shown_ = false;
    readout_text_.clear();
    if (!readout_) return;
    TTTOOLINFOW ti{sizeof(ti)};
    ti.hwnd = hwnd_;
    ti.uId = 1;
    SendMessageW(readout_, TTM_TRACKACTIVATE, FALSE, reinterpret_cast<LPARAM>(&ti));
}

// --- the window procedure --------------------------------------------------------------------------

void MainWindow::save_state() {
    if (options_.state_file.empty()) return;
    WINDOWPLACEMENT wp{sizeof(wp)};
    if (GetWindowPlacement(hwnd_, &wp)) {
        const RECT& r = wp.rcNormalPosition;
        state_.layout.x = r.left;
        state_.layout.y = r.top;
        state_.layout.width = r.right - r.left;
        state_.layout.height = r.bottom - r.top;
        state_.layout.maximized = wp.showCmd == SW_SHOWMAXIMIZED;
    }
    (void)app::save_state(state_, options_.state_file);
}

LRESULT CALLBACK MainWindow::proc(HWND hwnd, UINT message, WPARAM w, LPARAM l) {
    if (message == WM_NCCREATE) {
        auto* self = static_cast<MainWindow*>(reinterpret_cast<CREATESTRUCTW*>(l)->lpCreateParams);
        SetWindowLongPtrW(hwnd, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(self));
        self->hwnd_ = hwnd;
    }
    auto* self = reinterpret_cast<MainWindow*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
    return self ? self->handle(message, w, l) : DefWindowProcW(hwnd, message, w, l);
}

LRESULT CALLBACK MainWindow::panel_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l) {
    if (message == WM_NCCREATE) {
        SetWindowLongPtrW(hwnd, GWLP_USERDATA,
                          reinterpret_cast<LONG_PTR>(reinterpret_cast<CREATESTRUCTW*>(l)->lpCreateParams));
    }
    auto* self = reinterpret_cast<MainWindow*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
    if (self) {
        switch (message) {
        case WM_COMMAND:
            if (hwnd == self->console_panel_) {  // Unit 26j
                self->on_console(LOWORD(w), HIWORD(w));
                return 0;
            }
            if (LOWORD(w) >= id_section_first && LOWORD(w) <= id_section_last && HIWORD(w) == BN_CLICKED) {  // look A
                const int k = LOWORD(w) - id_section_first;
                self->set_section_folded(k, !self->section_folded(k));
                return 0;
            }
            if (LOWORD(w) == id_play && HIWORD(w) == BN_CLICKED) {
                self->execute(Command::play_pause);
                return 0;
            }
            self->on_inspector(LOWORD(w), HIWORD(w));
            return 0;
        case WM_NOTIFY:
            if (hwnd == self->console_panel_) return self->on_console_notify(*reinterpret_cast<const NMHDR*>(l), l);
            break;
        case WM_CONTEXTMENU:
            if (hwnd == self->console_panel_) {
                self->console_menu(GET_X_LPARAM(l), GET_Y_LPARAM(l));
                return 0;
            }
            break;
        case WM_DRAWITEM: {
            const auto& d = *reinterpret_cast<const DRAWITEMSTRUCT*>(l);
            if (hwnd == self->inspector_ && d.CtlID >= id_section_first && d.CtlID <= id_section_last) {
                self->draw_section_header(d);
                return TRUE;
            }
            break;
        }
        case WM_HSCROLL:  // a slider (Unit 26c, YD5, YD6)
            if (l != 0) self->on_slider(reinterpret_cast<HWND>(l));
            return 0;
        case WM_SIZE:
            if (hwnd == self->console_panel_) self->layout_console();
            else if (hwnd == self->inspector_) self->layout_inspector();
            return 0;
        case WM_CTLCOLORSTATIC:
            SetBkMode(reinterpret_cast<HDC>(w), TRANSPARENT);
            return reinterpret_cast<LRESULT>(GetSysColorBrush(COLOR_WINDOW));
        default:
            break;
        }
    }
    return DefWindowProcW(hwnd, message, w, l);
}

LRESULT MainWindow::handle(UINT message, WPARAM w, LPARAM l) {
    // ZD5: the shell made the taskbar button (again after Explorer restarts): its progress shown from now on.
    if (taskbar_created_ != 0 && message == taskbar_created_) {
        if (!taskbar_list_) {
            ITaskbarList3* list = nullptr;
            if (SUCCEEDED(CoCreateInstance(CLSID_TaskbarList, nullptr, CLSCTX_INPROC_SERVER, IID_PPV_ARGS(&list)))) {
                if (SUCCEEDED(list->HrInit())) {
                    taskbar_list_ = list;
                } else {
                    list->Release();
                }
            }
        }
        show_taskbar();  // the state kept, on the new button
        return 0;
    }
    switch (message) {
    case WM_CREATE: {
        if (!create_children()) return -1;
        build_menu();
        build_accelerators();
        // Unit 26d: files dropped on the window (ZD4); the taskbar button's message (ZD5); the title bar's theme (ZD6).
        DragAcceptFiles(hwnd_, TRUE);
        taskbar_created_ = RegisterWindowMessageW(L"TaskbarButtonCreated");
        apply_title_theme();
        auto r = render::Renderer::create({.adapter = options_.adapter, .debug_layer = options_.debug_layer});
        if (r) {
            renderer_ = std::move(*r);
            set_status(status_adapter, renderer_->adapter_name());
        } else {
            renderer_error_ = r.error().message;
            set_status(status_adapter, "no Direct3D 12 device");
        }
        apply_dpi(win32::dpi_of(hwnd_));
        update_title();
        // The bundled examples (VQ11): beside the executable, else in the source tree of this build.
        if (options_.examples_dir.empty()) {
            wchar_t exe[MAX_PATH]{};
            GetModuleFileNameW(nullptr, exe, MAX_PATH);
            examples_ = list_examples(examples_directory(fs::path(exe).parent_path(), fs::path(NITCAD_EXAMPLES_DIR)));
        } else {
            examples_ = list_examples(options_.examples_dir);
        }
        refresh();  // the view's placeholder message
        select_tab(tab_start);
        return 0;
    }
    case WM_SIZE:
        layout();
        return 0;
    case WM_DPICHANGED: {
        const RECT* r = reinterpret_cast<const RECT*>(l);
        apply_dpi(HIWORD(w));
        SetWindowPos(hwnd_, nullptr, r->left, r->top, r->right - r->left, r->bottom - r->top,
                     SWP_NOZORDER | SWP_NOACTIVATE);
        layout();
        return 0;
    }
    case WM_GETMINMAXINFO: {
        auto* mm = reinterpret_cast<MINMAXINFO*>(l);
        mm->ptMinTrackSize = {scaled(800, dpi_), scaled(520, dpi_)};
        return 0;
    }
    case WM_COMMAND: {
        const int id = LOWORD(w);
        if (id == IDOK || id == IDCANCEL) return 0;  // Enter and Escape through dialog navigation
        if (l == 0 || reinterpret_cast<HWND>(l) == toolbar_) {
            execute(static_cast<Command>(id));
            return 0;
        }
        break;
    }
    case WM_NOTIFY: {
        const auto* n = reinterpret_cast<const NMHDR*>(l);
        if (n->hwndFrom == tabs_ && n->code == TCN_SELCHANGE) {
            select_tab(static_cast<int>(SendMessageW(tabs_, TCM_GETCURSEL, 0, 0)));
            return 0;
        }
        if (n->hwndFrom == tree_ && n->code == TVN_SELCHANGEDW) {
            const auto* t = reinterpret_cast<const NMTREEVIEWW*>(l);
            if (t->action != TVC_UNKNOWN) on_tree_selection(t->itemNew.lParam);
            return 0;
        }
        // A box clicked or toggled by the space bar (YD4): the tree changes it after telling, so the boxes are
        // read once it has.
        if (n->hwndFrom == tree_ && (n->code == NM_CLICK || (n->code == TVN_KEYDOWN &&
                                                             reinterpret_cast<const NMTVKEYDOWN*>(l)->wVKey == VK_SPACE))) {
            PostMessageW(hwnd_, tree_check_message, 0, 0);
            return 0;
        }
        break;
    }
    case WM_INITMENUPOPUP: {
        if (HIWORD(l)) break;  // the system menu
        HMENU popup = reinterpret_cast<HMENU>(w);
        HMENU bar = GetMenu(hwnd_);
        for (int i = 0; i < GetMenuItemCount(bar); ++i) {
            if (GetSubMenu(bar, i) == popup) {
                if (i == 0 || (i >= 2 && i <= 5)) {
                    section_items_.clear();
                    key_items_.clear();
                    fill_menu(popup, i);
                }
                break;
            }
        }
        update_command_states();
        return 0;
    }
    case WM_DROPFILES:  // ZD4
        on_drop(reinterpret_cast<HDROP>(w));
        return 0;
    case WM_SETTINGCHANGE:  // ZD6: the app theme changed
        if (l != 0 && std::wstring_view(reinterpret_cast<const wchar_t*>(l)) == L"ImmersiveColorSet") apply_title_theme();
        break;
    case run_message:
        on_run_events();
        return 0;
    case live_message:
        on_live_check();
        return 0;
    case tree_check_message:
        sync_region_boxes();
        return 0;
    case WM_TIMER:
        if (w == run_timer) {
            on_run_timer();
            return 0;
        }
        if (w == live_timer) {
            KillTimer(hwnd_, live_timer);
            live_edits_ = edits_;
            live_.submit(deck_->text(), document_.path());
            return 0;
        }
        if (w == playback_timer) {
            on_playback_timer();
            return 0;
        }
        if (w == camera_timer) {
            on_camera_timer();
            return 0;
        }
        break;
    case WM_SETCURSOR: {
        if (reinterpret_cast<HWND>(w) == hwnd_ && LOWORD(l) == HTCLIENT) {
            POINT p;
            GetCursorPos(&p);
            ScreenToClient(hwnd_, &p);
            const auto in = [&](const Rect& r) { return p.x >= r.x - 2 && p.x < r.right() + 2 && p.y >= r.y && p.y < r.bottom(); };
            if (in(panes_.left_splitter) || in(panes_.right_splitter)) {
                SetCursor(LoadCursorW(nullptr, IDC_SIZEWE));
                return TRUE;
            }
            const Rect& c = panes_.console_splitter;  // Unit 26j: the console's, along y
            if (p.x >= c.x && p.x < c.right() && p.y >= c.y - 2 && p.y < c.bottom() + 2) {
                SetCursor(LoadCursorW(nullptr, IDC_SIZENS));
                return TRUE;
            }
        }
        break;
    }
    case WM_LBUTTONDOWN: {
        const int x = GET_X_LPARAM(l), y = GET_Y_LPARAM(l);
        const auto in = [&](const Rect& r) { return x >= r.x - 2 && x < r.right() + 2 && y >= r.y && y < r.bottom(); };
        const Rect& c = panes_.console_splitter;
        const bool console = x >= c.x && x < c.right() && y >= c.y - 2 && y < c.bottom() + 2;
        dragging_splitter_ = in(panes_.left_splitter) ? 1 : in(panes_.right_splitter) ? 2 : console ? 3 : 0;
        if (dragging_splitter_) {
            SetCapture(hwnd_);
            drag_origin_ = dragging_splitter_ == 3 ? y : x;
            drag_width_ = dragging_splitter_ == 1   ? panes_.tree.width
                          : dragging_splitter_ == 2 ? panes_.inspector.width
                                                    : panes_.console.height;
        }
        return 0;
    }
    case WM_MOUSEMOVE:
        if (dragging_splitter_ == 3) {  // up makes the console taller (and opens it); never under 40 px
            state_.layout.console_collapsed = false;
            const int height = std::max(drag_width_ - (GET_Y_LPARAM(l) - drag_origin_), scaled(40, dpi_));
            state_.layout.console_height = MulDiv(height, 96, dpi_);
            layout();
        } else if (dragging_splitter_) {
            const int dx = GET_X_LPARAM(l) - drag_origin_;
            const int width = std::max(drag_width_ + (dragging_splitter_ == 1 ? dx : -dx), scaled(min_side_96, dpi_));
            (dragging_splitter_ == 1 ? state_.layout.tree_width : state_.layout.inspector_width) = MulDiv(width, 96, dpi_);
            layout();
        }
        return 0;
    case WM_LBUTTONUP:
        if (dragging_splitter_) {
            dragging_splitter_ = 0;
            ReleaseCapture();
        }
        return 0;
    case WM_CLOSE:
        if (runner_.busy()) {
            if (MessageBoxW(hwnd_, L"A run is in progress. Stop it and close?", L"NiTCAD", MB_OKCANCEL | MB_ICONQUESTION) != IDOK) {
                return 0;
            }
            runner_.request_stop();
            runner_.wait();
        }
        if (!confirm_discard()) return 0;
        save_state();
        viewport_->release();
        run_page_->release();
        DestroyWindow(hwnd_);
        return 0;
    case WM_DESTROY:
        if (options_.quit_on_close) PostQuitMessage(0);
        return 0;
    case WM_NCDESTROY: {
        const HWND gone = hwnd_;
        SetWindowLongPtrW(gone, GWLP_USERDATA, 0);
        hwnd_ = nullptr;
        return DefWindowProcW(gone, message, w, l);
    }
    default:
        break;
    }
    return DefWindowProcW(hwnd_, message, w, l);
}

}  // namespace NiTCAD::app
