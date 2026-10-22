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
    // Unit 26b (XD6, XD7).
    quality,            // each element's angle (2D map, 3D surface or slices)
    quality_histogram,  // the elements per 5-degree bin
    difference,         // run B - run A of a field (2D map, 3D surface or slices, 1D profile)
    // Unit 26c (YD2).
    exploded,  // the regions moved apart along an axis (3D), the field on them or not
};

[[nodiscard]] std::string_view to_string(ViewKind v) noexcept;  // as `nitcad plot --view` names it
[[nodiscard]] std::optional<ViewKind> find_view(std::string_view name) noexcept;
// Whether the view draws a curve (a PlotScene) rather than the device.
[[nodiscard]] bool is_plot(ViewKind v) noexcept;
// Whether the view draws a setting's arrows and streamlines (Unit 26c, YD1): the field, slice and iso views;
// the others leave them out (build_scene) and so does plot_arguments.
[[nodiscard]] bool takes_vectors(ViewKind v) noexcept;
// Whether the view leaves out a setting's hidden regions (Unit 26c, YD4): the solid, mesh, field, slice, iso,
// exploded, quality and difference views; the curves ignore them (build_scene) and so does plot_arguments.
[[nodiscard]] bool takes_hidden(ViewKind v) noexcept;

// Transient playback (Unit 26c, YD5): the scale every snapshot of a transient run's field view is drawn on,
// fitted once over all of them (film_scale) and kept while the snapshot changes.
struct Film {
    render::ColorScale scale;
    double floor = 0.0;  // a vector magnitude on a log scale without a range: values below it drawn at it (0: none)
};

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
    // Arrows and streamlines (Unit 26c, YD1) of a vector field (E, J, Jn or Jp; empty: none) over a 2D or 3D
    // device view, as `nitcad plot --vectors --arrows --streamlines --lines --from --seeds` draw them: the counts
    // 0 for its defaults (300 arrows in 2D, 400 in 3D; 24 lines), the lines from a contact (empty: its choice)
    // or from typed seeds (a segment "<x0>,<y0>[,<z0>]:<x1>,<y1>[,<z1>]" or, in 3D, a plane "x=<um>").
    std::string vectors, streamlines;
    std::size_t arrows = 0, lines = 0;
    std::string from, seeds;
    // The exploded view (YD2): the axis (0 x, 1 y, 2 z), the gap in um (empty: 12% of the extent), and
    // whether the regions show the field (`--field`) or their own colors.
    int axis = 1;
    std::optional<double> gap_um;
    bool exploded_field = false;
    // Hidden regions (Unit 26c, YD4), by name, as `nitcad plot --hide`: left out of the device views
    // (takes_hidden) of a 2D or 3D device; the other views ignore them.
    std::vector<std::string> hidden;
    // Playback (YD5): with a film, a transient run's field view of a 2D or 3D device (takes_film) shows its
    // snapshot (`point`) as `nitcad plot --view filmstrip` draws that snapshot's frame (the time in the title)
    // on the film's scale; plot_arguments gives the snapshot's field view.
    std::optional<Film> film;
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
// Run B (Unit 26b, XD6), when given: the difference view needs it (on the same mesh, at the bias point of
// B equal to A's; refused otherwise), and the I-V, C-V and current views draw both runs' curves.
[[nodiscard]] std::expected<Scene, base::Error> build_scene(const Workspace& w, const ViewSettings& s,
                                                           const Workspace* b = nullptr);

// The views that apply to the workspace, in menu order: by dimension, and the curve each kind of
// result has; the difference with a run B on the same mesh.
[[nodiscard]] std::vector<ViewKind> available_views(const Workspace& w, const Workspace* b = nullptr);
// The fields a view can choose: every field_info key whose values the results provide (the rates
// only when stored; the doping always).
[[nodiscard]] std::vector<std::string> available_fields(const Workspace& w);

// Playback (Unit 26c, YD5). Whether the settings' view plays the workspace's snapshots: the field view of a
// transient run of a 2D or 3D device.
[[nodiscard]] bool takes_film(const Workspace& w, const ViewSettings& s) noexcept;
// One scale over every snapshot's field as `nitcad plot --view filmstrip --every 1` fits it (DC4): the
// range if given, over the shown regions' nodes, a vector magnitude on a log scale floored ten decades
// below its largest value over all. Errors (invalid_input): not takes_film; build_scene's for the field.
[[nodiscard]] std::expected<Film, base::Error> film_scale(const Workspace& w, const ViewSettings& s);
// The snapshot shown `elapsed_s` after playing from `from`, at playback_fps; the last one once reached. Frames
// not drawn in time are skipped, not queued.
inline constexpr double playback_fps = 5.0;
[[nodiscard]] std::size_t playback_frame(std::size_t from, std::size_t count, double elapsed_s) noexcept;

// The clipping box's sliders (Unit 26c, YD6): clip_steps steps over the device's extent along an axis
// [lo_um, hi_um], a step's value rounded to the nanometre, and the step nearest a value (clamped).
inline constexpr int clip_steps = 1000;
[[nodiscard]] double clip_value_um(double lo_um, double hi_um, int step) noexcept;
[[nodiscard]] int clip_step(double lo_um, double hi_um, double value_um) noexcept;

// The options of `nitcad plot` that draw the same view (without the input, --out, --size and
// --adapter); with run B, `--compare <its file>` for the views that use it.
[[nodiscard]] std::vector<std::string> plot_arguments(const ViewSettings& s, const Workspace* b = nullptr);

}  // namespace NiTCAD::app
