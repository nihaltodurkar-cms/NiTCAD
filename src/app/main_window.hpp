// The main window of nitcad-app (Unit 26, DECISIONS.md VA1-VA4, VA8, VQ13, WG5): the workspace of the
// owner's reference — a menu bar, a tool bar, the scene tree on the left, the tabbed centre (the device
// view with its metrics strip, the curves, the deck), the inspector on the right and the status bar —
// over the app logic (document, run controller, workspace, view settings, deck sections, state).
//
// It adds no numerics: decks are checked and run by the RunController on its worker thread, whose
// events arrive as a posted message, so the window never waits on a solve (WG1).
#pragma once

#include <windows.h>
#include <shellapi.h>  // HDROP (ZD4)
#include <commctrl.h>  // HIMAGELIST (look A)

#include <array>
#include <chrono>
#include <cstddef>
#include <filesystem>
#include <memory>
#include <cstdint>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/app/app_state.hpp"
#include "NiTCAD/app/commands.hpp"
#include "NiTCAD/app/console_log.hpp"
#include "NiTCAD/app/deck_sections.hpp"
#include "NiTCAD/app/document.hpp"
#include "NiTCAD/app/layout.hpp"
#include "NiTCAD/app/probe.hpp"
#include "NiTCAD/app/run_controller.hpp"
#include "NiTCAD/app/view_settings.hpp"
#include "NiTCAD/app/workspace.hpp"
#include "NiTCAD/render/axes3d.hpp"
#include "NiTCAD/render/image.hpp"
#include "NiTCAD/app/examples.hpp"
#include "NiTCAD/app/integration.hpp"
#include "NiTCAD/app/live_check.hpp"
#include "deck_pane.hpp"
#include "new_device_dialog.hpp"
#include "sweep_dialog.hpp"
#include "NiTCAD/app/sweep_results.hpp"
#include "geometry_page.hpp"
#include "solid_page.hpp"
#include "pages.hpp"
#include "viewport.hpp"

struct ITaskbarList3;  // shobjidl.h's, included where used (it brings objbase.h's `interface` macro)

namespace NiTCAD::app {

struct WindowOptions {
    render::Adapter adapter = render::Adapter::automatic;
    bool debug_layer = false;            // the Direct3D 12 debug layer (needs Graphics Tools)
    bool show = true;                    // shown at creation (tests make hidden windows)
    // Activated when shown and taking mouse and keyboard input; the self-test's window does neither, so a user
    // at the desktop cannot move, click or hover over it while it measures its message loop (WG1).
    bool input = true;
    bool quit_on_close = true;           // post WM_QUIT when destroyed (the application's window)
    // Smooth camera moves (Unit 26h, JD5), when Windows shows animations too. The application turns them
    // on; tests and the self-test leave them off, so a camera command lands at once.
    bool animate = false;
    std::filesystem::path state_file;    // the per-user state; empty: not loaded or saved
    // The bundled examples (VQ11): where they are (empty: beside the executable, else the source tree's),
    // and where a user's copies go (empty: the user's Documents, NiTCAD\examples).
    std::filesystem::path examples_dir;
    std::filesystem::path examples_copy_dir;
    std::filesystem::path devices_copy_dir;  // Unit 26n slice 1: where New device writes by default (empty: Documents/NiTCAD/devices)
};

// The posted message that wakes the window for the run controller's events.
inline constexpr UINT run_message = WM_APP + 1;
// The timer that shows a run's progress while its events pause.
inline constexpr UINT_PTR run_timer = 0x52;
// The deck is checked a pause after the last edit (Unit 26b, XD1): the timer, the pause and the message
// that wakes the window when a check ends.
inline constexpr UINT_PTR live_timer = 0x53;
inline constexpr UINT live_delay_ms = 400;
inline constexpr UINT live_message = WM_APP + 2;
// Posted when a box in the scene tree may have changed (Unit 26c, YD4).
inline constexpr UINT tree_check_message = WM_APP + 3;
// Unit 26n slice 1: the New device dialog closed (a deck created, or cancelled), handled after its own message.
inline constexpr UINT new_device_message = WM_APP + 4;
// Unit 26n slice 3: the sweep dialog closed (Run or Cancel), handled once it has returned from its own message.
inline constexpr UINT sweep_message = WM_APP + 5;
// Transient playback (Unit 26c, YD5): the timer that shows the snapshot due while playing. Each tick shows the
// snapshot the time since Play has reached (playback_frame), so a frame drawn late skips the ones it missed.
inline constexpr UINT_PTR playback_timer = 0x54;
inline constexpr UINT playback_tick_ms = 40;
// Unit 26h (JD5): the smooth camera move's frames, and its length.
inline constexpr UINT_PTR camera_timer = 0x55;
inline constexpr UINT camera_tick_ms = 16;
inline constexpr float camera_move_ms = 250.0f;

class MainWindow final : public Viewport::Host,
                         public DeckPane::Host,
                         public StartPage::Host,
                         public RunPage::Host,
                         public GeometryPage::Host,
                         public SolidPage::Host,
                         public NewDeviceDialog::Host,
                         public SweepDialog::Host {
public:
    enum Tab : int { tab_start = 0, tab_geometry, tab_view, tab_curves, tab_deck, tab_run, tab_results, tab_count };

    // nullptr when the window cannot be made (the reason on stderr is not available to a GUI: the
    // caller shows a message box).
    [[nodiscard]] static std::unique_ptr<MainWindow> create(const WindowOptions& options);
    ~MainWindow();
    MainWindow(const MainWindow&) = delete;
    MainWindow& operator=(const MainWindow&) = delete;

    [[nodiscard]] HWND hwnd() const noexcept { return hwnd_; }
    // The accelerators and dialog navigation; true when the message was handled.
    bool pre_translate(MSG& message);
    [[nodiscard]] HACCEL accelerators() const noexcept { return accelerators_; }

    void execute(Command c);
    // Opens a deck: its text in the deck pane, then a check on the worker; when it builds, its device
    // (and the results its outputs hold) in the views.
    bool open_deck(const std::filesystem::path& path);
    // A result file with its device file (and a generated mesh's file).
    bool open_result(const std::filesystem::path& result, const std::filesystem::path& device,
                     const std::filesystem::path& mesh = {});

    // State, for the self-test and the window tests.
    [[nodiscard]] RunController& runner() noexcept { return runner_; }
    [[nodiscard]] LiveChecker& live_checker() noexcept { return live_; }
    // Whether the views show the structure of the deck being edited, not solved (XD4).
    [[nodiscard]] bool preview() const noexcept { return preview_; }
    [[nodiscard]] const DeckDocument& document() const noexcept { return document_; }
    [[nodiscard]] const Workspace* workspace() const noexcept { return workspace_ ? &*workspace_ : nullptr; }
    [[nodiscard]] ViewSettings& settings() noexcept;  // of the active view tab
    [[nodiscard]] render::Renderer* renderer() noexcept { return renderer_ ? &*renderer_ : nullptr; }
    [[nodiscard]] const std::optional<std::string>& renderer_error() const noexcept { return renderer_error_; }
    [[nodiscard]] Viewport& viewport() noexcept { return *viewport_; }
    [[nodiscard]] DeckPane& deck() noexcept { return *deck_; }
    [[nodiscard]] StartPage& start_page() noexcept { return *start_; }
    [[nodiscard]] NewDeviceDialog* new_device_dialog() noexcept { return new_device_.get(); }  // Unit 26n slice 1; null when closed
    // Unit 26n slice 3 (DT6, DT7): the sweep dialog (null when closed); the last sweep's results (null before one);
    // the Curves tab's sweep views ("sweep: curves", then "sweep: <figure> vs <name>"), and choosing one.
    [[nodiscard]] SweepDialog* sweep_dialog() noexcept { return sweep_dialog_.get(); }
    [[nodiscard]] const SweepResults* sweep_results() const noexcept { return sweep_results_ ? &*sweep_results_ : nullptr; }
    [[nodiscard]] std::vector<std::string> sweep_views() const;
    void select_sweep_view(std::size_t k);
    [[nodiscard]] RunPage& run_page() noexcept { return *run_page_; }
    [[nodiscard]] ResultsPage& results_page() noexcept { return *results_; }
    [[nodiscard]] GeometryPage& geometry_page() noexcept { return *geometry_; }
    // Unit 26l: the 3D modeller, shown on the Geometry tab for a layered generate3d or an extrude3d deck.
    [[nodiscard]] SolidPage& solid_page() noexcept { return *solid_; }
    [[nodiscard]] bool solid_shown() const noexcept { return solid_shown_; }
    [[nodiscard]] const std::vector<Example>& examples() const noexcept { return examples_; }
    // Opens a copy of the example (in the copies' directory), as the start page does.
    bool open_example(std::size_t index);
    [[nodiscard]] const Panes& panes() const noexcept { return panes_; }
    [[nodiscard]] int dpi() const noexcept { return dpi_; }
    [[nodiscard]] int tab() const noexcept { return tab_; }
    void select_tab(int tab);
    // An adaptive run's pass k (Unit 17, AD8; the results page's passes list): its mesh in the device
    // view, the results kept aside; the last pass brings the results back. False (with the reason in the
    // status bar) when it cannot be shown.
    bool results_show_pass(std::size_t pass);
    // Builds the active view's scene from the workspace and settings, and presents it.
    void refresh();
    [[nodiscard]] const Scene* scene() const noexcept { return scene_ ? &*scene_ : nullptr; }
    [[nodiscard]] std::string status(int part) const;
    [[nodiscard]] HWND child(int id) const noexcept;
    [[nodiscard]] const std::optional<ProbeHit>& probe() const noexcept { return probe_; }
    // Unit 26h (JD4, JD5): a named view (render::named_views) by name; whether a camera move is under way,
    // and smooth moves turned on or off (tests).
    void show_named_view(std::string_view name);
    [[nodiscard]] bool camera_moving() const noexcept { return move_.has_value(); }
    void set_animate(bool on) noexcept { animate_ = force_animate_ = on; }  // on: whatever Windows says
    void probe_at(float x, float y);
    // The Draw cut line tool (Unit 26b, XD5): the first click sets the line's start, the second its end.
    void cut_at(float x, float y);
    // The cut line, shown on the device views and plotted along in the Curves tab (profile or bands); none
    // removes it.
    void set_cut(std::optional<ViewSettings::Cut> cut);
    [[nodiscard]] const std::optional<ViewSettings::Cut>& cut() const noexcept { return cut_; }
    // Run comparison (Unit 26b, XD6): a result file as run B beside the workspace (its device file and a
    // generated mesh's file as open_result takes them); File > Compare with a result asks for them.
    bool compare_result(const std::filesystem::path& result, const std::filesystem::path& device,
                        const std::filesystem::path& mesh = {});
    [[nodiscard]] const Workspace* run_b() const noexcept { return run_b_ ? &*run_b_ : nullptr; }
    // The deck's results that a run (or the structure preview) replaced: Simulation > Compare with the
    // previous run takes them as run B.
    [[nodiscard]] bool has_previous_run() const noexcept { return previous_.has_value(); }
    // Transient playback (Unit 26c, YD5): whether the snapshots are playing.
    [[nodiscard]] bool playing() const noexcept { return playing_; }
    // Unit 26d. A path opened as the recent-files list opens it (ZD4: open_target), after confirm_discard; the
    // taskbar button's progress (ZD5), kept whether or not the taskbar is there; whether the title bar is dark
    // and how many times the theme was given to it (ZD6: at creation and per theme change).
    bool open_path(const std::filesystem::path& path);
    [[nodiscard]] const TaskbarProgress& taskbar_progress() const noexcept { return taskbar_; }
    [[nodiscard]] bool dark_title() const noexcept { return dark_title_; }
    [[nodiscard]] int title_theme_updates() const noexcept { return title_theme_updates_; }
    [[nodiscard]] bool aurora_on() const noexcept { return state_.aurora; }  // View > Aurora theme
    // Unit 26e (HD1-HD3): the hover readout, whether it is shown and its text (kept while the window is hidden,
    // when the tooltip itself is not activated); View > Hover readout.
    [[nodiscard]] bool readout_shown() const noexcept { return readout_shown_; }
    [[nodiscard]] const std::string& readout_text() const noexcept { return readout_text_; }
    [[nodiscard]] bool hover_readout() const noexcept { return state_.hover_readout; }
    [[nodiscard]] const AppState& app_state() const noexcept { return state_; }
    // Unit 26j: the console along the bottom, a classic output window: its log; the tab shown (0 Output,
    // 1 Points, 2 Messages, as clicking the tab strip sets it); the points table and the rows its list shows; the
    // lines the log's rich edit shows; folded or not (Hide on its menu); Copy and Save (the table, or the
    // tab's lines); whether the tab's body is shown.
    [[nodiscard]] const ConsoleLog& console_log() const noexcept { return console_log_; }
    // The UI thread's time on a run (the speed pass): taking its events; the progress updates (at most every
    // 50 ms), and within them the console, the Points table and the run monitor; at its end, the console,
    // table and monitor, then the results shown. Reset when a run starts.
    struct UiTimings {
        std::size_t events = 0, updates = 0;
        double events_s = 0.0, update_s = 0.0, update_max_s = 0.0;
        double console_s = 0.0, points_s = 0.0, monitor_s = 0.0;
        double finish_s = 0.0, results_s = 0.0;
    };
    [[nodiscard]] const UiTimings& ui_timings() const noexcept { return timings_; }
    [[nodiscard]] int console_tab() const noexcept { return console_tab_; }
    void set_console_tab(int tab);
    [[nodiscard]] const PointsTable& points() const noexcept { return points_; }
    [[nodiscard]] int console_row_count() const;
    [[nodiscard]] int console_line_count() const;
    [[nodiscard]] bool console_collapsed() const noexcept { return state_.layout.console_collapsed; }
    // Look A: the inspector's sections (Display, Colour, Show, Slice and clip, Vectors, Curves, Playback, Probe),
    // each under a header that folds it; the folded ones saved with the layout. A section whose controls the view
    // does not use is hidden, header and all.
    static constexpr int section_count = 8;
    void set_section_folded(int section, bool folded);
    [[nodiscard]] bool section_folded(int section) const noexcept {
        return (state_.layout.inspector_folded >> section & 1) != 0;
    }
    // The toolbar's and the tree's icons drawn from the system's icon font (false: none found).
    [[nodiscard]] bool has_icons() const noexcept { return toolbar_images_ != nullptr; }
    [[nodiscard]] int missing_icons() const noexcept { return missing_icons_; }  // glyphs the font lacks
    [[nodiscard]] HWND status_progress() const noexcept { return status_progress_; }
    void set_console_collapsed(bool collapsed);
    bool copy_console();
    bool save_console(const std::filesystem::path& path);
    [[nodiscard]] bool console_body_shown() const;

    // The ids of the inspector's controls and the panes (for the tests).
    enum Control : int {
        id_tree = 4000, id_tabs, id_viewport, id_deck, id_metrics, id_inspector, id_toolbar, id_status, id_start,
        id_run_page, id_results, id_console, id_geometry, id_geometry_inspector, id_solid,
        id_geometry_controls = 7001,  // the Geometry section's controls (GeometryPage::Control from here)
        id_solid_controls = 7101,     // the 3D Geometry section's controls (SolidPage::Control from here)
        id_view = 5001, id_field, id_point, id_scale, id_colormap, id_contours, id_edges, id_contacts, id_opaque,
        id_perspective, id_slice_axis, id_slice_at, id_clip, id_clip_box, id_levels, id_cut, id_probe_text,
        // Unit 26c (YD1-YD3): arrows, streamlines and their seeds; the exploded view; the curves' pickers.
        id_arrows, id_arrow_count, id_lines, id_line_count, id_from, id_seeds, id_axis, id_gap, id_exploded_field,
        id_swept, id_sensed, id_frequency, id_log,
        // Unit 26c (YD5, YD6): the snapshot slider, its time and Play; the clipping box's low and high sliders.
        id_time, id_time_text, id_play, id_clip_x0, id_clip_x1, id_clip_y0, id_clip_y1, id_clip_z0, id_clip_z1,
        // Unit 26d (ZD3): the colour range's minimum and maximum, and Auto.
        id_range_min, id_range_max, id_range_auto,
        // Unit 26f (KD5): the enhanced 3D look.
        id_enhanced,
        // Unit 26h (JD1): the axis box.
        id_axis_box,
        // Unit 26j: the console's rich edit and tool strip.
        id_console_text, id_console_tabs, id_console_points,
        // Unit 26o (BSD5): the band surface's levels, colour by level, height scale and opacity.
        id_band_ec, id_band_ev, id_band_efn, id_band_efp, id_by_level, id_band_fit, id_exaggeration, id_band_opacity,
        // Look A: the inspector's section headers, a section each.
        id_section_first, id_section_last = id_section_first + section_count - 1,
        id_first_label = 5100,
    };

private:
    explicit MainWindow(const WindowOptions& options);
    static LRESULT CALLBACK proc(HWND hwnd, UINT message, WPARAM w, LPARAM l);
    static LRESULT CALLBACK panel_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l);
    LRESULT handle(UINT message, WPARAM w, LPARAM l);
    bool create_children();
    void build_menu();
    void build_accelerators();
    void fill_menu(HMENU menu, int index);
    void update_command_states();
    void apply_dpi(int dpi);
    void layout();
    void layout_inspector();
    void draw_section_header(const DRAWITEMSTRUCT& d);
    void make_icons();
    // The Aurora theme (View > Aurora theme) or the classic look, on every part the window paints.
    void apply_theme();
    LRESULT draw_toolbar(NMTBCUSTOMDRAW& cd);
    static LRESULT CALLBACK tabs_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l, UINT_PTR id, DWORD_PTR self);
    void update_status_progress();
    void update_title();
    // The state part's messages also go to the console, at the level given (Unit 26j); show_status only shows
    // them (the run's live progress).
    void set_status(int part, const std::string& text, LogLevel level = LogLevel::info);
    void show_status(int part, const std::string& text);
    // Unit 26j: a timestamped line in the console; the console's new lines into its rich edit (all of them
    // again on rebuild); its tool strip laid out; its controls' notifications.
    void log_line(LogLevel level, const std::string& text);
    void flush_console(bool rebuild = false);
    void layout_console();
    void on_console(int id, int code);
    void console_menu(int x, int y);  // the right-click menu at a screen point (-1, -1: the menu key)
    LRESULT on_console_notify(const NMHDR& n, LPARAM l);
    void update_points(bool force = false);
    void update_console_footer();
    [[nodiscard]] LogFilter console_text_filter() const noexcept;
    void on_run_events();
    void on_run_timer();
    void on_live_check();
    void show_structure(Workspace w);
    enum class Pending : unsigned char { none, open, check, run };
    void on_run_finished(Pending pending, RunState state);
    // Unit 26d: a drop (ZD4); the taskbar's progress for the controller's state (ZD5); the title bar's theme
    // from the user's app theme (ZD6).
    void on_drop(HDROP drop);
    void update_taskbar(bool run);
    void show_taskbar();  // the kept progress on the taskbar button, when there is one
    void apply_title_theme();
    void on_inspector(int id, int code);
    // Unit 26d (ZD1): the inspector's text boxes. Whether a control is one; one box taking effect (Enter, or
    // focus leaving it changed), its text shown again from the setting with the reason in the status bar when
    // not valid; and a box's text for the active view's settings.
    [[nodiscard]] static bool is_box(int id) noexcept;
    void commit_box(int id);
    [[nodiscard]] std::string box_text(int id) const;
    void show_box(int id);
    // ZD3: the range's boxes show the range set, or (Auto) the scale the scene was drawn on; keep_typed leaves
    // a box with focus and changed text as it is (a redraw while typing).
    void sync_range(bool keep_typed = false);
    void view_defaults(ViewSettings& s);
    void sync_inspector();
    // Shows the inspector's controls that apply to the active view (Unit 26c), hiding the others.
    void show_inspector_groups();
    // Whether a click or hover on the view can pick a node (not on the exploded view, whose regions moved).
    [[nodiscard]] bool pickable() const;
    void fill_tree();
    void on_tree_selection(LPARAM item);
    void sync_region_boxes();  // the regions shown as their boxes are (Unit 26c, YD4)
    void set_region_shown(const std::string& name, bool shown);
    // Playback (YD5): Play or Pause, a snapshot shown (pausing), the timer's tick, and the slider, time and
    // button kept as the settings are.
    [[nodiscard]] bool can_play() const;
    void play(bool on);
    void show_snapshot(std::size_t k);
    void on_playback_timer();
    // Unit 26h (JD5): move the camera to a pose, smoothly if animations are on; stop a move at its frame.
    void move_camera(const CameraPose& to);
    void stop_move();
    void on_camera_timer();
    void sync_playback();
    void sync_clip_sliders();  // YD6: the box's sides on its sliders
    // A slider of the inspector moved (WM_HSCROLL): the snapshot, or a side of the clipping box (YD6).
    void on_slider(HWND bar);
    // The device's extent along each axis, in um.
    [[nodiscard]] std::array<std::pair<double, double>, 3> extent_um() const;
    void update_metrics();
    // A new workspace; its view settings reset unless kept (the structure preview while editing).
    void set_workspace(Workspace w, bool keep_settings = false);
    // Run B shown beside the workspace, or none.
    void set_run_b(std::optional<Workspace> b);
    void default_settings();
    void present_camera();
    bool confirm_discard();
    bool save(bool as);
    void go_to_section(std::string_view kind);
    void insert_key(std::size_t index);
    void forward_edit(UINT message);
    void focus_next_pane();
    void save_state();

    // Viewport::Host
    void viewport_orbit(float azimuth_deg, float elevation_deg) override;
    void viewport_pan(float dx, float dy) override;
    void viewport_zoom(float factor) override;
    void viewport_fit() override;
    void viewport_click(float x, float y) override;
    void viewport_hover(float x, float y) override;
    void viewport_leave() override;
    bool viewport_is_3d() const override;
    bool viewport_over_cube(float x, float y) const override;
    void viewport_named_view(int index) override;
    // Unit 26h (JD3): snap the camera to a preset (zoom and pan reset as in Fit).
    void snap_to(const render::CameraPreset& preset);
    // DeckPane::Host
    void deck_changed() override;
    // GeometryPage::Host (Unit 26k).
    [[nodiscard]] std::string geometry_deck() override;
    void geometry_edited(const std::string& text) override;
    [[nodiscard]] std::optional<geo::Rect> geometry_bounds() override;
    void geometry_status(const std::string& cursor) override;
    void geometry_selected(const std::string& region) override;
    // SolidPage::Host (Unit 26l).
    [[nodiscard]] std::string solid_deck() override;
    void solid_edited(const std::string& text) override;
    void solid_status(const std::string& text) override;
    [[nodiscard]] render::Renderer* solid_renderer() override;
    void update_geometry_page();  // which of the two pages the Geometry tab shows
    // Unit 26n (DT2): the Geometry tab shows the device preview (the viewport) for a deck neither editor draws.
    [[nodiscard]] bool geometry_preview() const;
    // StartPage::Host and RunPage::Host
    void start_command(int command) override;
    void start_recent(std::size_t index) override;
    void start_example(std::size_t index) override;
    void new_device_created(const std::filesystem::path& deck) override;  // Unit 26n slice 1 (DT4)
    void new_device_closed() override;
    void open_new_device();
    void finish_new_device();
    // Unit 26n slice 3 (DT6): Sweep... on a dimension row; the dialog's Run and Cancel; the sweep started.
    void geometry_sweep(const std::string& name) override;
    void solid_sweep(const std::string& name) override;
    void open_sweep_dialog(const std::string& name);
    void sweep_run(const std::string& name, const std::vector<std::string>& values) override;
    void sweep_closed() override;
    void finish_sweep_dialog();
    void sweep_finished();
    void start_deck_run(const std::pair<std::string, std::vector<std::string>>* sweep);
    [[nodiscard]] std::expected<render::PlotScene, base::Error> sweep_scene() const;
    void run_page_stop() override;

    WindowOptions options_;
    HWND hwnd_ = nullptr;
    HWND toolbar_ = nullptr, status_ = nullptr, tree_ = nullptr, tabs_ = nullptr, metrics_ = nullptr,
         inspector_ = nullptr, console_panel_ = nullptr, console_ = nullptr, console_points_ = nullptr;
    // Look A: the icons' image lists (made per DPI), the status bar's progress bar, and the inspector's
    // controls the view does not use (show_inspector_groups; layout_inspector shows the others not folded).
    HIMAGELIST toolbar_images_ = nullptr, toolbar_disabled_ = nullptr, tree_images_ = nullptr;
    HWND status_progress_ = nullptr;
    int missing_icons_ = 0;
    std::set<int> unused_controls_;
    ConsoleLog console_log_;
    UiTimings timings_;
    PointsTable points_;
    int console_tab_ = 0;
    bool console_rebuild_ = false, points_stale_ = false;  // what waits while the console's tab is hidden
    LogFilter console_filter_shown_ = LogFilter::all;      // the lines the rich edit holds
    std::size_t messages_shown_ = static_cast<std::size_t>(-1);  // the count on the Messages tab
    std::size_t console_shown_ = 0;        // the log's lines (by number, ConsoleLog::total) flushed so far
    std::size_t console_lines_shown_ = 0;  // of those, the lines the filter shows in the rich edit
    std::unique_ptr<Viewport> viewport_;
    std::unique_ptr<DeckPane> deck_;
    std::unique_ptr<StartPage> start_;
    std::unique_ptr<RunPage> run_page_;
    std::unique_ptr<ResultsPage> results_;
    std::unique_ptr<GeometryPage> geometry_;
    std::unique_ptr<SolidPage> solid_;
    bool solid_shown_ = false;
    std::vector<Example> examples_;
    std::filesystem::path devices_dir_;  // the library (examples/devices), Unit 26n slice 1
    std::unique_ptr<NewDeviceDialog> new_device_;
    std::filesystem::path new_device_path_;  // the deck created, opened by finish_new_device
    // Unit 26n slice 3: the sweep dialog, the sweep it asked for (started once it has closed), the text swept, the
    // results, their figures and the sweep view shown on the Curves tab (-1: the run's own views).
    std::unique_ptr<SweepDialog> sweep_dialog_;
    std::optional<std::pair<std::string, std::vector<std::string>>> sweep_request_;
    std::string sweep_text_;
    std::optional<SweepResults> sweep_results_;
    SweepFigures sweep_figures_;
    int sweep_view_ = -1;
    HACCEL accelerators_ = nullptr;
    HFONT font_ = nullptr, bold_ = nullptr;
    HBRUSH canvas_brush_ = nullptr, field_brush_ = nullptr;  // the Aurora theme's canvas and boxes
    bool aurora_title_ = false;  // the title bar coloured by the Aurora theme (reset when it is turned off)
    int dpi_ = 96;
    int tab_ = tab_start;
    Panes panes_;
    AppState state_;

    DeckDocument document_ = DeckDocument::untitled();
    RunController runner_;
    LiveChecker live_;
    std::uint64_t edits_ = 0, live_edits_ = 0;  // edits made, and the count at the last live submission
    bool preview_ = false;
    Pending pending_ = Pending::none;
    std::chrono::steady_clock::time_point run_ui_updated_{};
    std::optional<Workspace> workspace_;
    // Unit 26i (ND5): bumped when workspace_ or values_ change; the status metrics cached against them.
    std::uint64_t workspace_generation_ = 0, values_generation_ = 0;
    struct Metrics {
        std::uint64_t workspace = ~std::uint64_t{0}, values = ~std::uint64_t{0};
        bool difference = false;
        std::string device, field;
        std::optional<std::string> quality;
    } metrics_cache_;
    std::optional<Workspace> run_b_, previous_;  // Unit 26b, XD6
    std::optional<Workspace> adapt_result_;      // Unit 17: the results while a pass's mesh is shown
    std::optional<Prober> prober_;
    ViewSettings device_view_, curve_view_;
    std::optional<Scene> scene_;
    std::vector<double> values_;  // the shown field's values, for probes
    std::string values_title_;    // their title, "<symbol> (<unit>)" (ZD7)
    std::optional<ProbeHit> probe_, measure_from_;
    bool probe_mode_ = true, measure_mode_ = false, cut_mode_ = false;
    struct CameraMove {
        CameraPose from, to;
        ViewSettings* settings = nullptr;  // the view it moves; another tab's view stops it
        std::chrono::steady_clock::time_point start;
    };
    std::optional<CameraMove> move_;
    bool animate_ = false, force_animate_ = false;
    std::optional<ViewSettings::Cut> cut_;          // Unit 26b, XD5
    std::optional<std::vector<double>> cut_from_;  // the first click of a cut line, um
    std::optional<render::Renderer> renderer_;
    std::optional<std::string> renderer_error_;
    // Playback (YD5): playing from which snapshot since when; the settings the device view's film was fitted for.
    bool playing_ = false;
    std::size_t play_from_ = 0;
    std::chrono::steady_clock::time_point play_start_{};
    std::string film_key_;
    // Unit 26d (ZD5, ZD6): the taskbar button (none until the shell says it is made), its progress shown, the
    // message that says so; the title bar's theme.
    ITaskbarList3* taskbar_list_ = nullptr;
    UINT taskbar_created_ = 0;
    bool com_ = false;  // this window initialized COM on its thread (balanced in the destructor)
    TaskbarProgress taskbar_;
    bool dark_title_ = false;
    int title_theme_updates_ = 0;
    // Unit 26e (HD3): the hover readout's tracking tooltip.
    HWND readout_ = nullptr;
    bool readout_shown_ = false;
    std::string readout_text_;
    void show_readout(const std::string& text, float x, float y);  // at the view's client pixel
    void hide_readout();

    // The menus' deck sections and keys (VQ13): what each dynamic item id stands for.
    struct SectionItem {
        std::string kind;
    };
    struct KeyItem {
        std::string kind;
        KeySpec key;
    };
    std::vector<SectionItem> section_items_;
    std::vector<KeyItem> key_items_;

    // Splitters.
    int dragging_splitter_ = 0;  // 1 left, 2 right, 3 the console's (a height, along y)
    int drag_origin_ = 0, drag_width_ = 0;
};

}  // namespace NiTCAD::app
