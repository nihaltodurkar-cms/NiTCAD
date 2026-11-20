// The 3D modeller's view (Unit 26l, slice 2: SD4), without Win32: the render scene of a solid model (each solid
// shaded by its material in layer order, edges outlined, contacts on the boundary, the selection, its gizmo and
// handles, the work plane's grid, what a tool is drawing), and picking by ray tests in model coordinates against
// the same camera the renderer draws with (render::view_frame).
//
// Drawing: every solid is drawn whole (clipped only by the section plane, whose cut faces are left open), its faces
// pushed out along their normals by a thousandth of a percent of the domain per layer, so where two layers' faces
// meet the later one shows, as the deck paints them. Translucent mode draws the solids see-through (the selection
// opaque). Coordinates on the screen: render's scene in µm, as the 3D results view.
//
// Picking, in order: the selection's handles (gizmo arrows within 5 px of their shafts, plane squares, face and
// vertex handles within 7 px; at 96 dpi, scaled), then the selected solid's edges within 4 px, then the nearest of
// the contacts' patch faces and the solids' faces along the ray (exact ray-plane hits, inside the face by exact
// tests in its plane), a later layer winning a tie; points cut away by the section plane do not count.
#pragma once

#include <array>
#include <optional>
#include <vector>

#include "NiTCAD/app/solid.hpp"
#include "NiTCAD/render/image.hpp"
#include "NiTCAD/render/scene.hpp"

namespace NiTCAD::app::solid {

using V3 = std::array<double, 3>;  // µm

[[nodiscard]] V3 um_of(P3 p) noexcept;
[[nodiscard]] P3 quanta_of(const V3& p) noexcept;  // rounded to the quantum

// A section plane across `axis` (-1: none): what lies past `at` along the axis is cut away.
struct Section {
    int axis = -1;
    Coord at = 0;
    friend bool operator==(const Section&, const Section&) = default;
};

struct ViewState {
    render::Camera camera;  // the default: the 3D view's isometric
    int width = 800, height = 600, dpi = 96;
    bool translucent = false;
    Section section;
    WorkPlane plane;
    Coord grid = 0;  // the work plane's grid step drawn (0: none)
};

struct Selection {
    std::vector<SolidRef> solids;  // a region's background is {region, -1}
    int contact = -1;
};

// What a tool is drawing, shown see-through over the model.
struct Preview {
    std::vector<Solid> solids;
    std::vector<P3> polyline;  // a prism's polygon so far
    std::optional<std::pair<P3, P3>> measure;
};

// A scene with only the camera and bounds (the domain, padded by 5 %): what picking and the handles project with.
[[nodiscard]] render::FieldScene frame_scene(const Model& m, const ViewState& v);
[[nodiscard]] render::FieldScene solid_scene(const Model& m, const ViewState& v, const Selection& s, const Preview& p = {});

// Where a point lands on the screen (pixels), and the ray under a pixel.
[[nodiscard]] std::array<double, 2> screen_of(const Model& m, const ViewState& v, const V3& p);
[[nodiscard]] render::ScreenRay ray_at(const Model& m, const ViewState& v, double x, double y);

// The handles of one selected solid: the gizmo's three arrows (index: the axis; `pts` its tail and tip) and three
// plane squares (index: the plane's normal; `pts` its corners) at the solid's bounds' centre, 70 px long; a handle
// at each face (index: the face, as move_face numbers them; `pts[0]` its centre); a prism's vertices on its `to`
// cap (index: the vertex).
enum class HandleKind : unsigned char { axis, plane, face, vertex };
struct Handle {
    HandleKind kind = HandleKind::axis;
    int index = 0;
    std::array<V3, 4> pts{};
};
[[nodiscard]] std::vector<Handle> handles(const Model& m, const ViewState& v, SolidRef s);

enum class PickKind : unsigned char { none, handle, edge, contact, face };
struct Pick {
    PickKind kind = PickKind::none;
    HandleKind handle = HandleKind::axis;
    int index = -1;  // the handle's, the edge's (Solid::edges), or the face's (move_face's numbering; -1 on a domain face)
    SolidRef solid;  // a face's solid ({background, -1} or {-1, -1} for the domain), the handle's or edge's
    int contact = -1;
    int domain_face = -1;  // the face of the domain the hit lies on (0..5 as a box's), else -1
    V3 at{};               // the hit point, µm
    friend bool operator==(const Pick&, const Pick&) = default;
};
[[nodiscard]] Pick pick(const Model& m, const ViewState& v, double x, double y, const Selection& s);

// Drags: the point of the plane across `normal` at `at` (µm) under a pixel (none when the ray runs along it), and
// the coordinate along `axis` of the point of the line through `origin` along that axis nearest the ray.
[[nodiscard]] std::optional<V3> on_plane(const Model& m, const ViewState& v, double x, double y, int normal, double at);
[[nodiscard]] double along_axis(const Model& m, const ViewState& v, double x, double y, const V3& origin, int axis);
// The snapping tolerance (8 px at 96 dpi, scaled) at a point, in quanta.
[[nodiscard]] Coord snap_tolerance(const Model& m, const ViewState& v, const V3& p);
// The grid step for the view's zoom at a point (the 2D editor's rule).
[[nodiscard]] Coord view_grid(const Model& m, const ViewState& v, const V3& p);

// A simple polygon's triangles (indices into it), exact on integer coordinates; collinear vertices allowed.
[[nodiscard]] std::vector<std::array<int, 3>> triangulate(const std::vector<Pt>& polygon);

}  // namespace NiTCAD::app::solid
