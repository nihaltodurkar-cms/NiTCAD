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
// The tests drive the view by its messages (SE8; no screen automation).
#pragma once

#include <windows.h>

#include <optional>
#include <string>
#include <vector>

#include "NiTCAD/app/solid.hpp"
#include "NiTCAD/app/solid_view.hpp"
#include "NiTCAD/render/image.hpp"

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
    enum class Tool : int { select, box, prism, contact, pan, measure };

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

    enum Control : int {
        c_heading = 0, c_hint, c_name_label, c_name, c_material_label, c_material, c_raise, c_lower,
        c_x0_label, c_x0, c_x1_label, c_x1, c_y0_label, c_y0, c_y1_label, c_y1, c_z0_label, c_z0, c_z1_label, c_z1,
        c_axis_label, c_axis, c_from_label, c_from, c_to_label, c_to, c_vertices_label, c_vertices, c_apply,
        c_kind_label, c_kind, c_draw_label, c_draw, c_plane_label, c_plane, c_offset_label, c_offset,
        c_section_label, c_section, c_cut_label, c_cut, c_translucent, c_convert, c_count
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
    std::vector<HWND> tool_buttons_;
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
    enum class Drag { none, orbit, pan, move, axis, plane, face, vertex, contact_move, base, height, polygon, contact, measure };
    Drag drag_ = Drag::none;
    bool button_ = false;  // the left button is down
    solid::Model m_start_;
    solid::SolidRef ref_;
    int index_ = -1;  // the handle's axis, plane normal, face or vertex
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
};

}  // namespace NiTCAD::app
