// The main window (main_window.hpp).
#include "main_window.hpp"

#include <commctrl.h>
#include <windowsx.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <charconv>
#include <format>
#include <tuple>
#include <utility>
#include <variant>

#include "NiTCAD/app/compare.hpp"
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
constexpr std::array<std::string_view, 5> colormap_names{"turbo", "viridis", "cividis", "coolwarm", "soft"};
// Unit 26c (YD1): the vector fields arrows and streamlines follow ("none" first).
constexpr std::array<std::string_view, 5> vector_names{"none", "E", "J", "Jn", "Jp"};
// The labels of the inspector's groups (Unit 26c), after the nine of Unit 26: their ids follow
// id_first_label in the order made.
constexpr int label_arrows = MainWindow::id_first_label + 9, label_lines = label_arrows + 1, label_from = label_arrows + 2,
              label_axis = label_arrows + 3, label_swept = label_arrows + 4, label_sensed = label_arrows + 5,
              label_frequency = label_arrows + 6;

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
      live_([this] { PostMessageW(hwnd_, live_message, 0, 0); }) {}

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
    if (accelerators_) DestroyAcceleratorTable(accelerators_);
    if (font_) DeleteObject(font_);
    if (bold_) DeleteObject(bold_);
}

bool MainWindow::create_children() {
    toolbar_ = make(TOOLBARCLASSNAME, L"", WS_VISIBLE | TBSTYLE_FLAT | TBSTYLE_LIST | CCS_NODIVIDER | CCS_NORESIZE |
                                               TBSTYLE_TOOLTIPS,
                    hwnd_, id_toolbar);
    SendMessageW(toolbar_, TB_BUTTONSTRUCTSIZE, sizeof(TBBUTTON), 0);
    SendMessageW(toolbar_, TB_SETEXTENDEDSTYLE, 0, TBSTYLE_EX_MIXEDBUTTONS);
    SendMessageW(toolbar_, TB_SETBITMAPSIZE, 0, MAKELPARAM(0, 0));
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
    for (const Button& b : buttons) {
        TBBUTTON t{};
        t.iBitmap = I_IMAGENONE;
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
    metrics_ = make(L"STATIC", L"", WS_VISIBLE | SS_LEFTNOWORDWRAP | SS_CENTERIMAGE, hwnd_, id_metrics);
    inspector_ = CreateWindowExW(WS_EX_CONTROLPARENT, panel_class, L"", WS_CHILD | WS_VISIBLE | WS_CLIPCHILDREN, 0, 0,
                                 10, 10, hwnd_, reinterpret_cast<HMENU>(static_cast<INT_PTR>(id_inspector)),
                                 GetModuleHandleW(nullptr), this);
    // The inspector's controls, each after its label.
    int label = id_first_label;
    const auto labelled = [&](const wchar_t* text, const wchar_t* cls, DWORD style, int id) {
        make(L"STATIC", text, WS_VISIBLE, inspector_, label++);
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
    paired(L"Arrows (count)", id_arrows, id_arrow_count);
    paired(L"Streamlines (count)", id_lines, id_line_count);
    paired(L"Seeds from (or typed, µm)", id_from, id_seeds);
    paired(L"Exploded along (gap, µm)", id_axis, id_gap);
    make(L"BUTTON", L"Sho&w the field", WS_VISIBLE | WS_TABSTOP | BS_AUTOCHECKBOX, inspector_, id_exploded_field);
    labelled(L"Swept contact", L"COMBOBOX", drop, id_swept);
    labelled(L"Sensed contact", L"COMBOBOX", drop, id_sensed);
    labelled(L"Frequency", L"COMBOBOX", drop, id_frequency);
    make(L"BUTTON", L"Lo&g axis", WS_VISIBLE | WS_TABSTOP | BS_AUTOCHECKBOX, inspector_, id_log);
    make(L"BUTTON", L"&Apply", WS_VISIBLE | WS_TABSTOP | BS_PUSHBUTTON, inspector_, id_apply);
    make(L"STATIC", L"", WS_VISIBLE | SS_LEFT, inspector_, id_probe_text);
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
            if (c.id == Command::fit || c.id == Command::perspective || c.id == Command::next_tab ||
                c.id == Command::tab_start || c.id == Command::find_error) {
                AppendMenuW(menu, MF_SEPARATOR, 0, nullptr);
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

void MainWindow::layout_inspector() {
    RECT r{};
    GetClientRect(inspector_, &r);
    const int pad = scaled(10, dpi_), row = scaled(23, dpi_), gap = scaled(4, dpi_), label_h = scaled(16, dpi_);
    const int w = std::max<int>(r.right - 2 * pad, 0);
    int y = pad;
    int label = id_first_label;
    const auto put = [&](int id, int height) {
        MoveWindow(child(id), pad, y, w, height, TRUE);
        y += height + gap;
    };
    const auto put_label = [&]() {
        MoveWindow(GetDlgItem(inspector_, label++), pad, y, w, label_h, TRUE);
        y += label_h;
    };
    for (int id : {id_view, id_field, id_point, id_scale, id_colormap}) {
        put_label();
        put(id, row * 8);  // a drop-down list's height includes its list
        y -= row * 8 - row;
    }
    for (int id : {id_contours, id_edges, id_contacts, id_opaque, id_perspective}) put(id, row);
    put_label();
    put(id_slice_axis, row * 6);
    y -= row * 6 - row;
    put_label();
    put(id_slice_at, row);
    put(id_clip, row);
    put(id_clip_box, row);
    put_label();
    put(id_levels, row);
    put_label();
    put(id_cut, row);
    // Unit 26c: the groups shown for the active view (show_inspector_groups), one under the other.
    const auto shown = [&](int id) { return (GetWindowLongPtrW(child(id), GWL_STYLE) & WS_VISIBLE) != 0; };
    const auto put_at = [&](int label_id) {
        MoveWindow(GetDlgItem(inspector_, label_id), pad, y, w, label_h, TRUE);
        y += label_h;
    };
    const int half = std::max((w - gap) / 2, 0);
    for (const auto& [label_id, list, box] : {std::tuple{label_arrows, id_arrows, id_arrow_count},
                                             std::tuple{label_lines, id_lines, id_line_count},
                                             std::tuple{label_from, id_from, id_seeds},
                                             std::tuple{label_axis, id_axis, id_gap}}) {
        if (!shown(list)) continue;
        put_at(label_id);
        MoveWindow(child(list), pad, y, half, row * 8, TRUE);
        MoveWindow(child(box), pad + half + gap, y, w - half - gap, row, TRUE);
        y += row + gap;
    }
    if (shown(id_exploded_field)) put(id_exploded_field, row);
    for (const auto& [label_id, list] : {std::pair{label_swept, id_swept}, std::pair{label_sensed, id_sensed},
                                        std::pair{label_frequency, id_frequency}}) {
        if (!shown(list)) continue;
        put_at(label_id);
        put(list, row * 8);
        y -= row * 8 - row;
    }
    if (shown(id_log)) put(id_log, row);
    MoveWindow(child(id_apply), pad, y, scaled(80, dpi_), row + gap, TRUE);
    y += row + 2 * gap;
    MoveWindow(child(id_probe_text), pad, y, w, std::max<int>(r.bottom - y - pad, 0), TRUE);
}

void MainWindow::update_title() {
    const std::string t = "NiTCAD - " + document_.title();
    win32::set_text(hwnd_, t);
    set_status(status_file, document_.is_untitled() ? std::string("untitled deck") : document_.path().string());
}

void MainWindow::set_status(int part, const std::string& text) {
    SendMessageW(status_, SB_SETTEXTW, static_cast<WPARAM>(part), reinterpret_cast<LPARAM>(win32::widen(text).c_str()));
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

void MainWindow::set_workspace(Workspace w, bool keep_settings) {
    // The results a run (or the structure preview) of the same deck replaces: the previous run (XD6).
    if (workspace_ && workspace_->has_results() && workspace_->source == w.source) previous_ = std::move(*workspace_);
    workspace_ = std::move(w);
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
    if (tab_ != tab_view && tab_ != tab_curves) return;
    scene_.reset();
    values_.clear();
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
    // The cut line where the view uses it: plotted along, or drawn on a 2D or 3D device view.
    ViewSettings shown = s;
    const bool along = s.kind == ViewKind::profile || s.kind == ViewKind::bands;
    if (along || (!is_plot(s.kind) && workspace_->device->mesh().dimension() > 1)) shown.cut = cut_;
    auto scene = build_scene(*workspace_, shown, run_b());
    if (!scene) {
        viewport_->show_message(scene.error().message);
        update_metrics();
        update_command_states();
        return;
    }
    scene_ = std::move(*scene);
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
            if (d) values_ = std::move(d->values);
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
            if (v) values_ = std::move(*v);
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
        viewport_->present();
        update_metrics();
    }
}

void MainWindow::update_metrics() {
    std::string text;
    if (workspace_) {
        const auto& m = workspace_->device->mesh();
        double lo[3]{1e300, 1e300, 1e300}, hi[3]{-1e300, -1e300, -1e300};
        for (const auto& p : m.points()) {
            for (int k = 0; k < m.dimension(); ++k) {
                lo[k] = std::min(lo[k], p[static_cast<std::size_t>(k)]);
                hi[k] = std::max(hi[k], p[static_cast<std::size_t>(k)]);
            }
        }
        text = "  Device ";
        for (int k = 0; k < m.dimension(); ++k) text += (k ? " × " : "") + format_um((hi[k] - lo[k]) * 1e4);
        text += " µm     Mesh " + std::to_string(m.node_count()) + " nodes";
        if (workspace_->geometry->element_count() > 0) {
            text += ", " + std::to_string(workspace_->geometry->element_count()) +
                    (m.dimension() == 3 ? " tetrahedra" : m.dimension() == 2 ? " triangles" : " segments");
        }
        const ViewKind kind = settings().kind;
        if (kind == ViewKind::quality || kind == ViewKind::quality_histogram) {
            // XD7: the minimum, the median and the count below the threshold.
            if (const auto angles = render::element_quality(*workspace_->geometry); angles && !angles->empty()) {
                const auto q = render::quality_summary(*angles, m.dimension());
                text += std::format("     Min angle {:.1f}°, median {:.1f}°, {} below {:g}°", q.min_deg, q.median_deg, q.below,
                                    q.threshold_deg);
            }
        } else if (!values_.empty()) {
            double top = -1e300;
            for (const double v : values_) if (std::isfinite(v)) top = std::max(top, v);
            text += (kind == ViewKind::difference ? "     Max |B − A| " : "     Max " + settings().field + " ") +
                    number(kind == ViewKind::difference ? render::difference_scale(values_).hi() : top);
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
    std::vector<std::string> points;
    if (workspace_) {
        for (std::size_t k = 0; k < workspace_->points.size(); ++k) points.push_back(std::to_string(k));
    }
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
    check(id_clip, s.clip.has_value());
    SendMessageW(child(id_slice_axis), CB_SETCURSEL, s.slices.empty() ? 0 : static_cast<WPARAM>(s.slices.front().axis + 1), 0);
    win32::set_text(child(id_slice_at), s.slices.empty() ? std::string() : std::format("{}", s.slices.front().at_um));
    if (s.clip) {
        const auto& c = *s.clip;
        win32::set_text(child(id_clip_box), std::format("{},{},{}:{},{},{}", c.lo_um[0], c.lo_um[1], c.lo_um[2], c.hi_um[0],
                                                        c.hi_um[1], c.hi_um[2]));
    }
    std::string levels;
    for (const double v : s.levels) levels += (levels.empty() ? "" : ",") + std::format("{}", v);
    win32::set_text(child(id_levels), levels);
    std::string cut;
    if (cut_) {
        for (const double v : cut_->from_um) cut += (cut.empty() ? "" : ",") + std::format("{}", v);
        cut += ":";
        for (std::size_t k = 0; k < cut_->to_um.size(); ++k) cut += (k ? "," : "") + std::format("{}", cut_->to_um[k]);
    }
    win32::set_text(child(id_cut), cut);
    // Unit 26c (YD1-YD3).
    const auto select = [&](int id, std::span<const std::string_view> names, const std::string& name) {
        const auto it = std::ranges::find(names, name.empty() ? std::string_view("none") : std::string_view(name));
        SendMessageW(child(id), CB_SETCURSEL, static_cast<WPARAM>(it == names.end() ? 0 : it - names.begin()), 0);
    };
    select(id_arrows, vector_names, s.vectors);
    select(id_lines, vector_names, s.streamlines);
    win32::set_text(child(id_arrow_count), s.arrows > 0 ? std::to_string(s.arrows) : std::string());
    win32::set_text(child(id_line_count), s.lines > 0 ? std::to_string(s.lines) : std::string());
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
    win32::set_text(child(id_seeds), s.seeds);
    SendMessageW(child(id_axis), CB_SETCURSEL, static_cast<WPARAM>(std::clamp(s.axis, 0, 2)), 0);
    win32::set_text(child(id_gap), s.gap_um ? std::format("{}", *s.gap_um) : std::string());
    check(id_exploded_field, s.exploded_field);
    contact_list(id_swept, s.kind == ViewKind::current ? "(all)" : "(automatic)", s.contact);
    contact_list(id_sensed, "(automatic)", s.sense);
    std::vector<std::string> frequencies;
    if (workspace_ && workspace_->small_signal) {
        for (const double hz : workspace_->small_signal->frequency_Hz) frequencies.push_back(std::format("{:g} Hz", hz));
    }
    combo_fill(child(id_frequency), frequencies, static_cast<int>(s.frequency.value_or(0)));
    check(id_log, s.log);
    show_inspector_groups();
}

void MainWindow::show_inspector_groups() {
    const ViewSettings& s = settings();
    const int dim = workspace_ ? workspace_->device->mesh().dimension() : 0;
    // Arrows and streamlines over the field views of a 2D or 3D run's results (YD1).
    const bool vectors = tab_ == tab_view && dim >= 2 && workspace_->has_results() && takes_vectors(s.kind);
    const bool exploded = tab_ == tab_view && s.kind == ViewKind::exploded;
    const bool curves = tab_ == tab_curves && (s.kind == ViewKind::iv || s.kind == ViewKind::cv || s.kind == ViewKind::current);
    const auto show = [&](std::initializer_list<int> ids, bool on) {
        for (const int id : ids) ShowWindow(child(id), on ? SW_SHOW : SW_HIDE);  // the labels too (in the inspector)
    };
    show({label_arrows, id_arrows, id_arrow_count, label_lines, id_lines, id_line_count, label_from, id_from, id_seeds}, vectors);
    show({label_axis, id_axis, id_gap, id_exploded_field}, exploded);
    win32::set_text(GetDlgItem(inspector_, label_swept), s.kind == ViewKind::current ? "Contact" : "Swept contact");
    show({label_swept, id_swept, id_log}, curves);
    show({label_sensed, id_sensed}, curves && s.kind != ViewKind::current);
    show({label_frequency, id_frequency}, curves && s.kind == ViewKind::cv);
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
    } else if ((id == id_apply && clicked) || (id == id_clip && clicked) || (id == id_slice_axis && selection)) {
        // Unit 26c: the arrow and line counts (empty: nitcad plot's defaults), the typed seeds and the exploded
        // view's gap (empty: 12% of the extent), each kept as it was when not valid.
        const auto count = [&](int box, std::size_t most, std::size_t& to, const char* what) {
            const std::string t = win32::text_of(child(box));
            if (t.find_first_not_of(" \t") == std::string::npos) {
                to = 0;
                return;
            }
            const auto v = numbers(t);
            if (v.size() == 1 && v[0] >= 1.0 && v[0] <= static_cast<double>(most) && v[0] == std::floor(v[0])) {
                to = static_cast<std::size_t>(v[0]);
            } else {
                set_status(status_state, std::format("{} are a count from 1 to {}.", what, most));
            }
        };
        count(id_arrow_count, 100000, s.arrows, "Arrows");
        count(id_line_count, 10000, s.lines, "Streamlines");
        std::string seeds = win32::text_of(child(id_seeds));
        while (!seeds.empty() && seeds.back() == ' ') seeds.pop_back();
        while (!seeds.empty() && seeds.front() == ' ') seeds.erase(seeds.begin());
        s.seeds = seeds;
        const std::string gap = win32::text_of(child(id_gap));
        if (gap.find_first_not_of(" \t") == std::string::npos) {
            s.gap_um.reset();
        } else if (const auto v = numbers(gap); v.size() == 1 && v[0] >= 0.0) {
            s.gap_um = v[0];
        } else {
            set_status(status_state, "The gap is a length in µm, at least 0.");
        }
        // The slice, the clipping box and the iso levels from their fields.
        const int axis = combo_index(child(id_slice_axis)) - 1;
        const std::string at = win32::text_of(child(id_slice_at));
        s.slices.clear();
        if (const auto v = numbers(at); axis >= 0 && v.size() == 1) s.slices.push_back({axis, v[0]});
        s.clip.reset();
        if (Button_GetCheck(child(id_clip)) == BST_CHECKED) {
            ViewSettings::Clip c;
            std::string box = win32::text_of(child(id_clip_box));
            std::ranges::replace(box, ':', ',');
            if (const auto v = numbers(box); v.size() == 6) {
                for (std::size_t k = 0; k < 3; ++k) {
                    c.lo_um[k] = v[k];
                    c.hi_um[k] = v[k + 3];
                }
                s.clip = c;
            }
        }
        s.levels.clear();
        s.levels = numbers(win32::text_of(child(id_levels)));
        // The cut line typed: "<x0>,<y0>[,<z0>]:<x1>,<y1>[,<z1>]" in um; empty removes it.
        const std::string typed = win32::text_of(child(id_cut));
        const auto colon = typed.find(':');
        const std::size_t dim = static_cast<std::size_t>(workspace_->device->mesh().dimension());
        const std::vector<double> a = colon == std::string::npos ? std::vector<double>{} : numbers(typed.substr(0, colon));
        const std::vector<double> b = colon == std::string::npos ? std::vector<double>{} : numbers(typed.substr(colon + 1));
        if (typed.find_first_not_of(" \t") == std::string::npos) {
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
            set_status(status_state, "A cut line is <x0>,<y0>[,<z0>]:<x1>,<y1>[,<z1>] in µm, in the device's dimension.");
        }
    } else {
        return;
    }
    refresh();
}

// --- the scene tree --------------------------------------------------------------------------------

void MainWindow::fill_tree() {
    TreeView_DeleteAllItems(tree_);
    if (!workspace_) return;
    // Only the regions of a 2D or 3D device have a check box (YD4): checked while shown. The tree gives every
    // new item one, so the others' is taken away.
    const auto add = [&](HTREEITEM parent, const std::string& text, LPARAM data, int check = 0) {
        TVINSERTSTRUCTW t{};
        t.hParent = parent;
        t.hInsertAfter = TVI_LAST;
        t.item.mask = TVIF_TEXT | TVIF_PARAM;
        const std::wstring w = win32::widen(text);
        t.item.pszText = const_cast<wchar_t*>(w.c_str());
        t.item.lParam = data;
        HTREEITEM item = TreeView_InsertItem(tree_, &t);
        TreeView_SetItemState(tree_, item, INDEXTOSTATEIMAGEMASK(check), TVIS_STATEIMAGEMASK);
        return item;
    };
    const auto& d = *workspace_->device;
    const bool boxes = d.mesh().dimension() >= 2;
    HTREEITEM device = add(TVI_ROOT, preview_ ? "Device (preview, not solved)" : "Device", tree_data(tree_none, 0));
    for (std::size_t k = 0; k < d.regions().size(); ++k) {
        const auto& r = d.regions()[k];
        const bool hidden = std::ranges::find(device_view_.hidden, r.name) != device_view_.hidden.end();
        add(device, r.name + (device::is_insulator(r) ? " (insulator)" : " (semiconductor)"),
            boxes ? tree_data(tree_region, static_cast<int>(k)) : tree_data(tree_none, 0), boxes ? (hidden ? 1 : 2) : 0);
    }
    for (const auto& c : d.contacts()) {
        const char* kind = c.kind == device::ContactKind::ohmic ? "ohmic" : c.kind == device::ContactKind::gate ? "gate" : "electrode";
        add(device, "contact " + c.name + " (" + kind + ")", tree_data(tree_none, 0));
    }
    TreeView_Expand(tree_, device, TVE_EXPAND);
    add(TVI_ROOT, "Mesh: " + std::to_string(d.mesh().node_count()) + " nodes", tree_data(tree_none, 0));
    HTREEITEM views = add(TVI_ROOT, "Views", tree_data(tree_none, 0));
    const auto kinds = available_views(*workspace_, run_b());
    for (std::size_t k = 0; k < kinds.size(); ++k) add(views, std::string(to_string(kinds[k])), tree_data(tree_view, static_cast<int>(k)));
    TreeView_Expand(tree_, views, TVE_EXPAND);
    HTREEITEM fields = add(TVI_ROOT, "Fields", tree_data(tree_none, 0));
    const auto names = available_fields(*workspace_);
    for (std::size_t k = 0; k < names.size(); ++k) {
        const auto f = render::find_field(names[k]);
        add(fields, std::string(f ? render::field_info(*f).title : names[k]), tree_data(tree_field, static_cast<int>(k)));
    }
    TreeView_Expand(tree_, fields, TVE_EXPAND);
    if (workspace_->has_results()) add(TVI_ROOT, "Results: " + std::to_string(workspace_->points.size()) + " points", tree_data(tree_none, 0));
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
    run_b_.reset();
    previous_.reset();
    cut_.reset();
    prober_.reset();
    scene_.reset();
    fill_tree();
    pending_ = Pending::open;
    runner_.start_check(document_.text(), document_.path());
    set_status(status_state, "Building the deck...");
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

void MainWindow::on_run_events() {
    if (!runner_.process()) return;
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
    run_page_->update(runner_, renderer(), true);  // before the results are taken: how the run ended
    deck_->set_diagnostics(runner_.diagnostics());
    on_run_finished(pending, s);
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
        set_status(status_state, text);
        run_page_->update(runner_, renderer());
        update_command_states();
    }
}

void MainWindow::on_run_finished(Pending pending, RunState s) {
    if (s == RunState::failed) {
        const std::string why = runner_.error() ? runner_.error()->message : std::string("failed");
        set_status(status_state, (pending == Pending::run ? "The run failed: " : "The deck does not build: ") + why);
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
        set_status(status_state, "The deck builds.");
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
                                         : std::format("Done in {:.1f} s", seconds));
    }
    update_command_states();
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
        if (k < state_.recent.size() && confirm_discard()) {
            const fs::path p = state_.recent[k];
            if (p.extension() == ".deck") {
                open_deck(p);
            } else {
                fs::path m = fs::path(p).replace_extension(".nitmesh");
                if (!fs::exists(m)) m.clear();
                open_result(p, fs::path(p).replace_extension(".nitdev"), m);
            }
        }
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
        runner_.start_check(document_.text(), document_.path());
        set_status(status_state, "Checking...");
        update_command_states();
        return;
    case Command::run:
        if (runner_.busy()) return;
        document_.set_text(deck_->text());
        if (document_.is_untitled() && !save(true)) return;  // a run writes its outputs beside the deck
        select_tab(tab_run);  // the run monitor follows it (shown before the run starts)
        pending_ = Pending::run;
        runner_.start_run(document_.text(), document_.path());
        set_status(status_state, "Running...");
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
    case Command::show_regions:
        // YD4: every region shown again, the tree's boxes checked.
        if (!workspace_ || device_view_.hidden.empty()) return;
        device_view_.hidden.clear();
        fill_tree();
        set_status(status_state, "Showing every region.");
        if (tab_ == tab_view) refresh();
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
    orbit(settings(), azimuth_deg, elevation_deg);
    present_camera();
}

void MainWindow::viewport_pan(float dx, float dy) {
    if (!scene_ || !std::holds_alternative<render::FieldScene>(*scene_)) return;
    pan_by(settings(), dx, dy);
    present_camera();
}

void MainWindow::viewport_zoom(float factor) {
    if (!scene_ || !std::holds_alternative<render::FieldScene>(*scene_)) return;
    zoom_by(settings(), factor);
    present_camera();
}

void MainWindow::viewport_fit() {
    fit(settings());
    present_camera();
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

void MainWindow::viewport_click(float x, float y) {
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
    if (!pickable()) return;
    const auto* f = std::get_if<render::FieldScene>(&*scene_);
    if (!f) return;
    const auto hit = prober_->pick(*f, viewport_->width(), viewport_->height(), x, y, {}, 8.0f);
    set_status(status_cursor, hit ? std::format("x {}  y {}  z {} µm", format_um(hit->position_um[0]),
                                                format_um(hit->position_um[1]), format_um(hit->position_um[2]))
                                  : std::string());
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
            self->on_inspector(LOWORD(w), HIWORD(w));
            return 0;
        case WM_SIZE:
            self->layout_inspector();
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
    switch (message) {
    case WM_CREATE: {
        if (!create_children()) return -1;
        build_menu();
        build_accelerators();
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
        }
        break;
    }
    case WM_LBUTTONDOWN: {
        const int x = GET_X_LPARAM(l), y = GET_Y_LPARAM(l);
        const auto in = [&](const Rect& r) { return x >= r.x - 2 && x < r.right() + 2 && y >= r.y && y < r.bottom(); };
        dragging_splitter_ = in(panes_.left_splitter) ? 1 : in(panes_.right_splitter) ? 2 : 0;
        if (dragging_splitter_) {
            SetCapture(hwnd_);
            drag_origin_ = x;
            drag_width_ = dragging_splitter_ == 1 ? panes_.tree.width : panes_.inspector.width;
        }
        return 0;
    }
    case WM_MOUSEMOVE:
        if (dragging_splitter_) {
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
