// The 3D modeller on the Geometry tab (Unit 26l slice 2; DECISIONS.md SD4, SQ3): for a layered generate3d or an
// extrude3d deck the tab shows this page in the 2D editor's place. A D3D12 viewport through the render layer draws
// the model (NiTCAD/app/solid_view.hpp); its tools edit it, and the Geometry section in the inspector's place shows
// the selection's numbers.
//
// - The view: right drag orbits, middle drag (or the Pan tool, or Space+drag) pans, the wheel zooms at the cursor,
//   F fits; the named views from the view menu's keys (numpad), as the 3D view. T toggles translucency; the
//   inspector sets the section plane and the work plane (XY, YZ, ZX at an offset; W takes the face under the cursor).
// - Tools, each with a key while the view has the focus: Select (V: pick a solid, a contact or a handle; drag a
//   solid along the face it was picked by, an arrow of the gizmo along its axis, a square in its plane, a face handle
//   along the face's normal, a prism's vertex handle in its cap's plane; Shift adds to the selection; Ctrl drags a
//   vertex alone), Box (B or R: a rectangle on the work plane, then its height, a click ending it), Prism (P: the
//   polygon's vertices clicked on the work plane, each kept along u or v from the last; Enter, a double click or the
//   first vertex closes it; then its height), Contact (C: a rectangle dragged on a face of the domain), Pan (H),
//   Measure (M: the 3D distance between two snapped points). Arrows nudge by a grid step in the work plane,
//   Page Up and Down along its normal; Delete; Escape cancels.
// - Every edit goes to the deck at once by minimal edits; undo and redo step through the deck texts the page wrote.
//   An extrude3d deck is shown and selected but not edited until Convert (SQ5) makes it layered.
//
// - Unit 26m slice 0 (MT1-MT3): the fixed tool row is a group switcher; the armed group's tools float beside the
//   selection; Camera free / fixed (L) decides whether the mouse moves the view (modeller.hpp); Reset view; the status
//   bar shows the deck line a drag will write, the chip at the cursor its value, typed digits setting it exactly.
// - Unit 26m slice 1 (MT4, MT5): Fill / Replace and the inspector's overlap switch; the Transform modes: Copy (or
//   Alt+drag) moves a copy, Mirror clicks (about the centre; Ctrl: the face's plane; Shift: a quarter turn about the
//   view's axis), Array drags (pitch; the wheel sets the count), Align clicks (the picked solids' matching face to
//   the clicked face's plane). Each is one deck edit and one undo step.
//
// The tests drive the view by its messages (SE8; no screen automation).
#pragma once

#include <windows.h>

#include <expected>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/app/document.hpp"
#include "NiTCAD/app/solid.hpp"
#include "NiTCAD/app/solid_view.hpp"
#include "NiTCAD/render/image.hpp"
#include "modeller_bar.hpp"

#include <memory>

namespace NiTCAD::app {

class SolidPage {
public:
    struct Host {
        [[nodiscard]] virtual std::string solid_deck() = 0;     // the deck's text now
        virtual void solid_edited(const std::string& text) = 0;  // the deck after an edit, undo or redo
        virtual void solid_status(const std::string& text) = 0;  // the cursor's position, a measure, a refusal
        [[nodiscard]] virtual render::Renderer* solid_renderer() = 0;

    protected:
        ~Host() = default;
    };

    // The 2D editor's tools' order, so the Edit menu's tool commands work on both.
    enum class Tool : int { select, box, prism, contact, pan, measure, cylinder };  // cylinder: 3D only (Unit 26m MT7)

    SolidPage(HWND parent, Host& host, int id, int inspector_id);
    ~SolidPage();
    SolidPage(const SolidPage&) = delete;
    SolidPage& operator=(const SolidPage&) = delete;

    [[nodiscard]] HWND hwnd() const noexcept { return panel_; }
    [[nodiscard]] HWND canvas() const noexcept { return canvas_; }
    [[nodiscard]] HWND inspector() const noexcept { return inspector_; }
    void set_dpi(int dpi, HFONT font, HFONT bold);
    void layout_inspector();
    // The swap chain released (before the renderer is).
    void release();

    void reload(bool fit = false);
    void fit();
    void set_tool(Tool t);
    [[nodiscard]] Tool tool() const noexcept { return tool_; }
    void set_camera(const render::Camera& c);
    void set_work_plane(solid::WorkPlane p);
    void set_translucent(bool on);
    void set_section(solid::Section s);
    void named_view(int index);  // render::named_views()
    bool undo();
    bool redo();
    bool duplicate();
    bool remove();
    bool convert();  // an extrude3d deck made layered (SQ5)
    bool select_region(const std::string& name);
    bool commit(int control);

    // --- Unit 26m slice 0 (MT1-MT3) ---
    void set_camera_fixed(bool fixed);
    [[nodiscard]] bool camera_fixed() const noexcept { return camera_fixed_; }
    void reset_view();  // the camera the page opened with
    void set_group(modeller::Group g);
    [[nodiscard]] modeller::Group group() const noexcept { return bar_->group(); }
    [[nodiscard]] const ModellerBar& bar() const noexcept { return *bar_; }
    // The floating bar's place over the canvas (canvas pixels), and the selection's screen bounds it was put beside.
    [[nodiscard]] modeller::Placement bar_placement() const;
    [[nodiscard]] std::optional<modeller::Box> selection_bounds() const;
    [[nodiscard]] const std::string& live_line() const noexcept { return live_line_; }  // "" when no drag
    [[nodiscard]] std::string chip() const { return bar_->chip_text(); }
    [[nodiscard]] const modeller::Entry& typed() const noexcept { return entry_; }
    // --- Unit 26m slice 1 (MT4, MT5) ---
    void set_mode(modeller::Mode m);  // a Transform mode of the Select tool (modeller.hpp)
    [[nodiscard]] modeller::Mode mode() const noexcept { return mode_; }
    bool toggle_overlap();  // the picked solids' regions: fill <-> replace
    [[nodiscard]] int array_count() const noexcept { return array_count_; }
    [[nodiscard]] int view_axis() const;  // the model axis most nearly along the view
    // --- Unit 26m slice 2 (MT12) ---
    [[nodiscard]] solid::Coord resolution() const noexcept { return res_q_; }  // the deck's s, quanta
    [[nodiscard]] int chords() const noexcept { return chords_; }
    // The largest s the layout allows (cm; +inf without interfaces), or why none does.
    [[nodiscard]] std::expected<double, std::string> automatic_resolution() const;
    bool use_automatic_resolution();  // written to the deck's interface_resolution (on the quantum, rounded down)
    // Unit 26m MT8: the inspector's parameters, a `name = value` line each; applied: each changed or new one set (every
    // tied face following), one undo step. False (the status says why, nothing written) when one does not resolve.
    [[nodiscard]] std::string parameters_text() const;
    bool apply_parameters(const std::string& lines);

    // State, for the tests.
    [[nodiscard]] const solid::Model& model() const noexcept { return m_; }
    [[nodiscard]] const solid::ViewState& view() const noexcept { return v_; }
    [[nodiscard]] const solid::Selection& selection() const noexcept { return sel_; }
    [[nodiscard]] std::string material() const { return material_; }
    void set_material(std::string m) { material_ = std::move(m); }
    [[nodiscard]] bool can_undo() const noexcept { return history_.can_undo(); }
    [[nodiscard]] bool can_redo() const noexcept { return history_.can_redo(); }
    [[nodiscard]] const std::string& message() const noexcept { return message_; }
    // Where a model point lands in the view's client area.
    [[nodiscard]] std::array<double, 2> screen_of(solid::P3 p) const;
    [[nodiscard]] std::array<double, 2> screen_of_um(const solid::V3& p) const;
    // Draws now (not waiting for WM_PAINT): the present's time in ms; -1 when nothing could be drawn.
    double render();
    [[nodiscard]] const render::FieldScene& scene() const noexcept { return scene_; }
    [[nodiscard]] double scene_ms() const noexcept { return scene_ms_; }  // the last render's scene build
    [[nodiscard]] const std::optional<std::string>& error() const noexcept { return error_; }

    // --- live validation and the mesh preview (slice 3, SD5; as the 2D page's GD9) ---
    // A live check of `text` ended: its diagnostics (empty: the deck builds). Kept while `text` is the text shown.
    void checked(const std::string& text, std::span<const Diagnostic> diagnostics);
    enum class Check { pending, valid, invalid };
    [[nodiscard]] Check check() const noexcept;
    [[nodiscard]] const std::vector<solid::Problem>& problems() const noexcept { return problems_; }
    // Why Run is off on this tab ("" when it is not): the layout's first problem.
    [[nodiscard]] std::string run_blocked() const;
    // The surface mesh (boundary and interface triangles' edges, pairs of points in µm) a live check of `text` built.
    void set_mesh(const std::string& text, std::vector<solid::V3> edges);
    void set_mesh_preview(bool on);
    [[nodiscard]] bool mesh_preview() const noexcept { return mesh_preview_; }
    [[nodiscard]] bool mesh_shown() const noexcept { return mesh_preview_ && !mesh_text_.empty() && mesh_text_ == text_seen_; }
    [[nodiscard]] std::size_t mesh_segments() const noexcept { return mesh_.size() / 2; }
    [[nodiscard]] const std::vector<solid::V3>& mesh() const noexcept { return mesh_; }
    // The problem or mesh line shown at the inspector's top ("" none).
    [[nodiscard]] std::string banner() const;

    enum Control : int {
        c_heading = 0, c_hint, c_name_label, c_name, c_material_label, c_material, c_raise, c_lower,
        c_x0_label, c_x0, c_x1_label, c_x1, c_y0_label, c_y0, c_y1_label, c_y1, c_z0_label, c_z0, c_z1_label, c_z1,
        c_axis_label, c_axis, c_from_label, c_from, c_to_label, c_to, c_vertices_label, c_vertices, c_apply,
        c_kind_label, c_kind, c_draw_label, c_draw, c_plane_label, c_plane, c_offset_label, c_offset,
        c_section_label, c_section, c_cut_label, c_cut, c_translucent, c_convert, c_overlap, c_resolution, c_params_label, c_params, c_params_apply, c_count
    };
    [[nodiscard]] HWND control(int c) const;

private:
    static LRESULT CALLBACK panel_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l);
    static LRESULT CALLBACK canvas_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l);
    static LRESULT CALLBACK inspector_proc(HWND hwnd, UINT message, WPARAM w, LPARAM l);
    LRESULT on_canvas(UINT message, WPARAM w, LPARAM l);
    void layout();
    void on_down(double x, double y, bool shift, bool ctrl, bool double_click);
    void on_move(double x, double y, bool shift, bool ctrl);
    void on_up(double x, double y);
    bool on_key(UINT key, bool shift, bool ctrl);
    void cancel();
    bool commit_model(solid::Model next);
    void apply_text(const std::string& text);
    void update_inspector();
    void update_tool_buttons();
    void invalidate();
    void prune_selection();
    void close_polygon();
    // A point on the work plane under the pixel, snapped (grid, vertices, edges, faces), on the plane.
    [[nodiscard]] std::optional<solid::P3> plane_point(double x, double y) const;
    [[nodiscard]] solid::Coord quantize(double um) const;  // a drag's length, to the grid step
    [[nodiscard]] solid::Coord grid() const;

    Host& host_;
    HWND panel_ = nullptr, strip_ = nullptr, canvas_ = nullptr, inspector_ = nullptr;
    std::vector<HWND> controls_;
    int inspector_id_ = 0;
    int dpi_ = 96;

    render::Renderer* renderer_ = nullptr;
    std::optional<render::WindowTarget> target_;
    render::FieldScene scene_;
    double scene_ms_ = 0.0;
    std::optional<std::string> error_;

    solid::Model m_;
    std::string text_seen_;
    std::string message_;
    solid::EditHistory history_;
    solid::ViewState v_;
    solid::Selection sel_;
    Tool tool_ = Tool::select;
    std::string material_ = "silicon";

    // A drag or a drawing in progress.
    enum class Drag { none, orbit, pan, move, axis, plane, face, vertex, contact_move, base, height, polygon, contact, measure, array,
                      cut, taper, extrude, inset, inset_depth, radius };
    Drag drag_ = Drag::none;
    bool button_ = false;  // the left button is down
    solid::Model m_start_;
    solid::SolidRef ref_;
    int index_ = -1;  // the handle's axis, plane normal, face or vertex
    // MT8: N during a face drag, with nothing typed, names the extent it drags: -1 not named, 0 a new parameter
    // (naming_name_), k > 0 the k-th existing length parameter.
    int naming_ = -1;
    std::string naming_name_;
    int axis_ = 2;    // the axis a drag runs along, or the normal of the plane it runs in
    solid::V3 origin_{};
    double start_s_ = 0.0;
    solid::V3 start_p_{};
    solid::P3 a_{}, b_{};
    solid::Coord height_ = 0;
    std::vector<solid::P3> polygon_;
    solid::Preview preview_;
    double last_x_ = 0, last_y_ = 0;
    bool space_ = false, alone_ = false;

    // The live check.
    std::string checked_text_, mesh_text_;
    std::vector<Diagnostic> diagnostics_;
    std::vector<solid::Problem> problems_;
    std::vector<solid::V3> mesh_;
    bool mesh_preview_ = false;
    void locate_problems();

    // Unit 26m slice 0.
    std::unique_ptr<ModellerBar> bar_;
    bool camera_fixed_ = false;  // MQ2: free on open
    render::Camera start_camera_;
    modeller::Entry entry_;
    int typed_axis_ = -1;         // a two-axis drag's axis the typed value sets
    solid::P3 drag_d_{};          // the drag's offset now (move, plane, axis, vertex)
    solid::Coord drag_s_ = 0;     // a face's offset now
    std::string live_line_;
    int canvas_top_ = 0;
    void after_drag_step(double x, double y);  // the deck line and the chip
    [[nodiscard]] std::optional<solid::Model> pending() const;
    [[nodiscard]] std::string chip_value() const;
    bool on_char(wchar_t c);
    void apply_drag();  // Enter: the drag ended as it stands
    void finish_height();
    void place_bar();
    void end_drag_ui();

    // Unit 26m slice 1.
    modeller::Mode mode_ = modeller::Mode::none;
    int array_count_ = 3;
    std::optional<solid::Model> copy_from_;  // a copy drag: the model before the copies were made
    solid::Selection copy_sel_;
    bool transform_click(const solid::Pick& p, bool shift, bool ctrl);  // Mirror and Align

    // Unit 26m slice 2 (MT6, MT7, MT12).
    bool tool_press(const solid::Pick& p, double x, double y);  // the Faces and Edges modes' presses
    void begin_drag(Drag d);
    void finish_inset();
    void read_mesh_keys();
    [[nodiscard]] solid::Coord gap_up(solid::Coord v) const;  // a non-zero length under 4 s raised to 4 s (MT12)
    [[nodiscard]] solid::Model cut_model(solid::Coord r);     // the edge cut (Chamfer, Fillet, Round vertex) of size r
    int edge_ = -1;                // the edge a cut or taper drags
    int chords_ = 0;               // a fillet's chords (0: under 6 s)
    solid::Normal normal_;         // the face an inset or extrude starts on
    solid::Coord margin_ = 0, radius_ = 0;
    solid::Coord res_q_ = 0;       // the deck's interface_resolution, quanta (0: none)
    double max_size_ = 0.0;        // the deck's max_size, cm
};

}  // namespace NiTCAD::app
