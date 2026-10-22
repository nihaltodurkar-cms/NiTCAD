// The main window of nitcad-app (Unit 26, DECISIONS.md VA1-VA4, VA8, VQ13, WG5): the workspace of the
// owner's reference — a menu bar, a tool bar, the scene tree on the left, the tabbed centre (the device
// view with its metrics strip, the curves, the deck), the inspector on the right and the status bar —
// over the app logic (document, run controller, workspace, view settings, deck sections, state).
//
// It adds no numerics: decks are checked and run by the RunController on its worker thread, whose
// events arrive as a posted message, so the window never waits on a solve (WG1).
#pragma once

#include <windows.h>

#include <chrono>
#include <filesystem>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "NiTCAD/app/app_state.hpp"
#include "NiTCAD/app/commands.hpp"
#include "NiTCAD/app/deck_sections.hpp"
#include "NiTCAD/app/document.hpp"
#include "NiTCAD/app/layout.hpp"
#include "NiTCAD/app/probe.hpp"
#include "NiTCAD/app/run_controller.hpp"
#include "NiTCAD/app/view_settings.hpp"
#include "NiTCAD/app/workspace.hpp"
#include "NiTCAD/render/image.hpp"
#include "NiTCAD/app/examples.hpp"
#include "NiTCAD/app/live_check.hpp"
#include "deck_pane.hpp"
#include "pages.hpp"
#include "viewport.hpp"

namespace NiTCAD::app {

struct WindowOptions {
    render::Adapter adapter = render::Adapter::automatic;
    bool debug_layer = false;            // the Direct3D 12 debug layer (needs Graphics Tools)
    bool show = true;                    // shown at creation (tests make hidden windows)
    // Activated when shown and taking mouse and keyboard input; the self-test's window does neither, so a user
    // at the desktop cannot move, click or hover over it while it measures its message loop (WG1).
    bool input = true;
    bool quit_on_close = true;           // post WM_QUIT when destroyed (the application's window)
    std::filesystem::path state_file;    // the per-user state; empty: not loaded or saved
    // The bundled examples (VQ11): where they are (empty: beside the executable, else the source tree's),
    // and where a user's copies go (empty: the user's Documents, NiTCAD\examples).
    std::filesystem::path examples_dir;
    std::filesystem::path examples_copy_dir;
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

class MainWindow final : public Viewport::Host,
                         public DeckPane::Host,
                         public StartPage::Host,
                         public RunPage::Host {
public:
    enum Tab : int { tab_start = 0, tab_view, tab_curves, tab_deck, tab_run, tab_results, tab_count };

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
    [[nodiscard]] RunPage& run_page() noexcept { return *run_page_; }
    [[nodiscard]] ResultsPage& results_page() noexcept { return *results_; }
    [[nodiscard]] const std::vector<Example>& examples() const noexcept { return examples_; }
    // Opens a copy of the example (in the copies' directory), as the start page does.
    bool open_example(std::size_t index);
    [[nodiscard]] const Panes& panes() const noexcept { return panes_; }
    [[nodiscard]] int dpi() const noexcept { return dpi_; }
    [[nodiscard]] int tab() const noexcept { return tab_; }
    void select_tab(int tab);
    // Builds the active view's scene from the workspace and settings, and presents it.
    void refresh();
    [[nodiscard]] const Scene* scene() const noexcept { return scene_ ? &*scene_ : nullptr; }
    [[nodiscard]] std::string status(int part) const;
    [[nodiscard]] HWND child(int id) const noexcept;
    [[nodiscard]] const std::optional<ProbeHit>& probe() const noexcept { return probe_; }
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

    // The ids of the inspector's controls and the panes (for the tests).
    enum Control : int {
        id_tree = 4000, id_tabs, id_viewport, id_deck, id_metrics, id_inspector, id_toolbar, id_status, id_start,
        id_run_page, id_results,
        id_view = 5001, id_field, id_point, id_scale, id_colormap, id_contours, id_edges, id_contacts, id_opaque,
        id_perspective, id_slice_axis, id_slice_at, id_clip, id_clip_box, id_levels, id_cut, id_apply, id_probe_text,
        // Unit 26c (YD1-YD3): arrows, streamlines and their seeds; the exploded view; the curves' pickers.
        id_arrows, id_arrow_count, id_lines, id_line_count, id_from, id_seeds, id_axis, id_gap, id_exploded_field,
        id_swept, id_sensed, id_frequency, id_log,
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
    void update_title();
    void set_status(int part, const std::string& text);
    void on_run_events();
    void on_run_timer();
    void on_live_check();
    void show_structure(Workspace w);
    enum class Pending : unsigned char { none, open, check, run };
    void on_run_finished(Pending pending, RunState state);
    void on_inspector(int id, int code);
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
    bool viewport_is_3d() const override;
    // DeckPane::Host
    void deck_changed() override;
    // StartPage::Host and RunPage::Host
    void start_command(int command) override;
    void start_recent(std::size_t index) override;
    void start_example(std::size_t index) override;
    void run_page_stop() override;

    WindowOptions options_;
    HWND hwnd_ = nullptr;
    HWND toolbar_ = nullptr, status_ = nullptr, tree_ = nullptr, tabs_ = nullptr, metrics_ = nullptr,
         inspector_ = nullptr;
    std::unique_ptr<Viewport> viewport_;
    std::unique_ptr<DeckPane> deck_;
    std::unique_ptr<StartPage> start_;
    std::unique_ptr<RunPage> run_page_;
    std::unique_ptr<ResultsPage> results_;
    std::vector<Example> examples_;
    HACCEL accelerators_ = nullptr;
    HFONT font_ = nullptr, bold_ = nullptr;
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
    std::optional<Workspace> run_b_, previous_;  // Unit 26b, XD6
    std::optional<Prober> prober_;
    ViewSettings device_view_, curve_view_;
    std::optional<Scene> scene_;
    std::vector<double> values_;  // the shown field's values, for probes
    std::optional<ProbeHit> probe_, measure_from_;
    bool probe_mode_ = true, measure_mode_ = false, cut_mode_ = false;
    std::optional<ViewSettings::Cut> cut_;          // Unit 26b, XD5
    std::optional<std::vector<double>> cut_from_;  // the first click of a cut line, um
    std::optional<render::Renderer> renderer_;
    std::optional<std::string> renderer_error_;

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
    int dragging_splitter_ = 0;  // 1 left, 2 right
    int drag_origin_ = 0, drag_width_ = 0;
};

}  // namespace NiTCAD::app
