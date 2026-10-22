// The main window of nitcad-app (Unit 26, DECISIONS.md VA1-VA4, VA8, VQ13, WG5): the workspace of the
// owner's reference — a menu bar, a tool bar, the scene tree on the left, the tabbed centre (the device
// view with its metrics strip, the curves, the deck), the inspector on the right and the status bar —
// over the app logic (document, run controller, workspace, view settings, deck sections, state).
//
// It adds no numerics: decks are checked and run by the RunController on its worker thread, whose
// events arrive as a posted message, so the window never waits on a solve (WG1).
#pragma once

#include <windows.h>

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
#include "deck_pane.hpp"
#include "viewport.hpp"

namespace NiTCAD::app {

struct WindowOptions {
    render::Adapter adapter = render::Adapter::automatic;
    bool debug_layer = false;            // the Direct3D 12 debug layer (needs Graphics Tools)
    bool show = true;                    // shown at creation (tests make hidden windows)
    bool quit_on_close = true;           // post WM_QUIT when destroyed (the application's window)
    std::filesystem::path state_file;    // the per-user state; empty: not loaded or saved
};

// The posted message that wakes the window for the run controller's events.
inline constexpr UINT run_message = WM_APP + 1;

class MainWindow final : public Viewport::Host, public DeckPane::Host {
public:
    enum Tab : int { tab_view = 0, tab_curves = 1, tab_deck = 2 };

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
    [[nodiscard]] const DeckDocument& document() const noexcept { return document_; }
    [[nodiscard]] const Workspace* workspace() const noexcept { return workspace_ ? &*workspace_ : nullptr; }
    [[nodiscard]] ViewSettings& settings() noexcept;  // of the active view tab
    [[nodiscard]] render::Renderer* renderer() noexcept { return renderer_ ? &*renderer_ : nullptr; }
    [[nodiscard]] const std::optional<std::string>& renderer_error() const noexcept { return renderer_error_; }
    [[nodiscard]] Viewport& viewport() noexcept { return *viewport_; }
    [[nodiscard]] DeckPane& deck() noexcept { return *deck_; }
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

    // The ids of the inspector's controls and the panes (for the tests).
    enum Control : int {
        id_tree = 4000, id_tabs, id_viewport, id_deck, id_metrics, id_inspector, id_toolbar, id_status,
        id_view = 5001, id_field, id_point, id_scale, id_colormap, id_contours, id_edges, id_contacts, id_opaque,
        id_perspective, id_slice_axis, id_slice_at, id_clip, id_clip_box, id_levels, id_apply, id_probe_text,
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
    void on_inspector(int id, int code);
    void view_defaults(ViewSettings& s);
    void sync_inspector();
    void fill_tree();
    void on_tree_selection(LPARAM item);
    void update_metrics();
    void set_workspace(Workspace w);
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

    WindowOptions options_;
    HWND hwnd_ = nullptr;
    HWND toolbar_ = nullptr, status_ = nullptr, tree_ = nullptr, tabs_ = nullptr, metrics_ = nullptr,
         inspector_ = nullptr;
    std::unique_ptr<Viewport> viewport_;
    std::unique_ptr<DeckPane> deck_;
    HACCEL accelerators_ = nullptr;
    HFONT font_ = nullptr, bold_ = nullptr;
    int dpi_ = 96;
    int tab_ = tab_view;
    Panes panes_;
    AppState state_;

    DeckDocument document_ = DeckDocument::untitled();
    RunController runner_;
    enum class Pending : unsigned char { none, open, check, run } pending_ = Pending::none;
    std::optional<Workspace> workspace_;
    std::optional<Prober> prober_;
    ViewSettings device_view_, curve_view_;
    std::optional<Scene> scene_;
    std::vector<double> values_;  // the shown field's values, for probes
    std::optional<ProbeHit> probe_, measure_from_;
    bool probe_mode_ = true, measure_mode_ = false;
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
