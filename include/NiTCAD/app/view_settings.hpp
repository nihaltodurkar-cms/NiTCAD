// What a view shows (Unit 26, DECISIONS.md VA4, VA8): the inspector's settings, turned into a render
// scene with the choices `nitcad plot` makes for the same options, so the window and the command
// draw the same picture (WG4). plot_arguments gives that command's options for a setting.
#pragma once

#include <array>
#include <cstddef>
#include <expected>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "NiTCAD/app/workspace.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/render/geometry.hpp"
#include "NiTCAD/render/scene.hpp"

namespace NiTCAD::app {

enum class ViewKind : unsigned char {
    solid,    // regions and contacts (2D, 3D)
    mesh,     // the mesh over the regions (2D, 3D)
    field,    // a field: 2D map, 3D surface, 1D profile
    slice,    // a field on slice planes (3D)
    iso,      // isosurfaces (3D) or contour levels (2D)
    bands,    // the band diagram (1D, or along a cut line)
    profile,  // a field along a 1D device or a cut line
    iv,       // a sweep's current against the swept bias
    cv,       // a small-signal run's capacitance against the swept bias
    current,  // a transient run's terminal currents against time
};

[[nodiscard]] std::string_view to_string(ViewKind v) noexcept;  // as `nitcad plot --view` names it
[[nodiscard]] std::optional<ViewKind> find_view(std::string_view name) noexcept;
// Whether the view draws a curve (a PlotScene) rather than the device.
[[nodiscard]] bool is_plot(ViewKind v) noexcept;

struct ViewSettings {
    ViewKind kind = ViewKind::solid;
    std::string field = "potential";             // field_info's key
    std::optional<std::size_t> point;            // the bias point, default the last
    std::string scale;                           // "", linear, log, signed-log ("" the field's own)
    std::string colormap = "turbo";
    std::optional<std::pair<double, double>> range;
    bool contours = false;
    bool edges = false;                          // mesh edges over a map or surface
    bool opaque = false;                         // insulators drawn opaque
    bool contacts = true;                        // electrode slabs and contact pads
    bool perspective = false;
    std::optional<float> azimuth, elevation;     // 3D camera, degrees (render's defaults when empty)
    // Slice planes and the clipping box in um, converted to cm as `nitcad plot` converts its options.
    struct Slice {
        int axis = 0;  // 0 x, 1 y, 2 z
        double at_um = 0.0;
    };
    struct Clip {
        std::array<double, 3> lo_um{}, hi_um{};
    };
    std::vector<Slice> slices;
    std::optional<Clip> clip;
    std::vector<double> levels;                  // iso levels
    // A cut line (Unit 26b, XD5), its ends in um in the device's dimension: the profile and bands views plot
    // along it (needed in 2D and 3D), and the 2D and 3D device views draw it.
    struct Cut {
        std::vector<double> from_um, to_um;
    };
    std::optional<Cut> cut;
    std::string contact, sense;                  // curves: the swept and sensed contacts by name
    std::optional<std::size_t> frequency;        // C-V: the frequency index
    bool log = false;                            // curves and profiles: a log y axis
    // The window's camera (VA8, Unit 26): zoom about the view's centre and pan in fractions of its
    // half size (render::Camera). `nitcad plot` has neither: plot_arguments leaves them out.
    float zoom = 1.0f;
    float pan_x = 0.0f, pan_y = 0.0f;
};

// Camera moves (VA8): orbit by degrees (3D; the elevation kept within +-89), zoom by a factor (kept
// within 1/16 to 1024), pan by fractions of the half view, and fit (zoom 1, no pan, the angles kept).
void orbit(ViewSettings& s, float azimuth_deg, float elevation_deg);
void zoom_by(ViewSettings& s, float factor);
void pan_by(ViewSettings& s, float dx, float dy);
void fit(ViewSettings& s);
// A built scene's camera set from the settings (3D: the angles and projection; any: zoom and pan), so
// a camera move redraws without building the scene again; build_scene applies it too.
void apply_camera(render::FieldScene& scene, const ViewSettings& s);

using Scene = std::variant<render::FieldScene, render::PlotScene>;

// The scene for the settings. Errors (invalid_input) as `nitcad plot` reports them: an unknown
// field, scale or colormap; a view the device or the results do not have (no results, a 3D-only
// view in 2D, ...); a point out of range; render's errors.
[[nodiscard]] std::expected<Scene, base::Error> build_scene(const Workspace& w, const ViewSettings& s);

// The views that apply to the workspace, in menu order: by dimension, and the curve each kind of
// result has.
[[nodiscard]] std::vector<ViewKind> available_views(const Workspace& w);
// The fields a view can choose: every field_info key whose values the results provide (the rates
// only when stored; the doping always).
[[nodiscard]] std::vector<std::string> available_fields(const Workspace& w);

// The options of `nitcad plot` that draw the same view (without the input, --out, --size and
// --adapter).
[[nodiscard]] std::vector<std::string> plot_arguments(const ViewSettings& s);

}  // namespace NiTCAD::app
