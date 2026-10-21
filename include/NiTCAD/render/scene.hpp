// The CPU scenes that drawing consumes (Unit 25a; DECISIONS.md VR1-VR4, VR8): plain data built from
// a device, a geometry and results. Drawing reads only these; nothing here draws.
//
// A FieldScene is triangles and line segments in world coordinates, in micrometres (the mesh's cm
// times 1e4), with a camera, an optional colorbar and an axis triad: the 2D maps and the 3D views.
// A PlotScene is curves on two axes: field profiles of 1D devices, band diagrams and Unit 24's
// curves (I-V, C-V, admittance against frequency, transients).
//
// Field triangles carry, per vertex, the colormap position of the value (ColorScale::position), so
// drawing interpolates positions and looks colors up per pixel; flat triangles carry a color.
// Regions are drawn by the control-volume rule of geometry.hpp, in the reference's palette (VR8):
// insulators pale and translucent, heavily doped n red, p channel blue, the lightly doped
// substrate beige, metal dark grey (doping classes below).
#pragma once

#include <array>
#include <cstdint>
#include <expected>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/analysis/curve.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/render/color.hpp"
#include "NiTCAD/render/geometry.hpp"
#include "NiTCAD/results/solution.hpp"

namespace NiTCAD::render {

using Float3 = std::array<float, 3>;

// --- field scenes -------------------------------------------------------------------------------

struct SceneVertex {
    Float3 position{};      // um
    Float3 normal{};        // unit; (0, 0, 1) in 2D
    float value = 0.0f;     // colormap position in [0, 1], NaN for no data (field batches)
    Rgba color{};           // flat batches
};

enum class Shading : std::uint8_t {
    flat,   // each vertex's color
    field,  // the colormap at the interpolated position
};

struct TriangleBatch {
    std::string name;  // what it shows: a region class, "field", "slice", "electrode"
    Shading shading = Shading::flat;
    bool lit = true;           // directional light and ambient (3D); 2D maps are unlit
    bool transparent = false;  // drawn with order-independent transparency (alpha below 1)
    std::vector<SceneVertex> vertices;  // three per triangle, counterclockwise from the front
};

struct LineBatch {
    std::string name;
    Rgba color{};
    float width_px = 1.0f;
    bool on_top = false;          // drawn over the triangles (2D outlines, mesh edges)
    std::vector<Float3> points;  // two per segment, um
};

struct Colorbar {
    Colormap map = Colormap::turbo;
    ColorScale scale;
    std::string title;  // the field and its unit, e.g. "n (cm⁻³)"
};

enum class Projection : std::uint8_t { orthographic, perspective };

// Looks at `target` from the direction given by azimuth (about +y, from +z toward +x) and elevation
// (above the xz plane), with +y up: the device's vertical axis, as in 2D, where the camera looks
// down -z (azimuth and elevation 0). `extent` is the half height of the bounds' view at aspect 1;
// drawing refits the bounds to the image's aspect.
struct Camera {
    Float3 target{};
    float azimuth_deg = 35.0f;
    float elevation_deg = 25.0f;
    float extent_um = 1.0f;
    Projection projection = Projection::orthographic;
    float fov_deg = 35.0f;
};

struct FieldScene {
    int dimension = 2;
    std::string title;
    std::vector<TriangleBatch> triangles;
    std::vector<LineBatch> lines;
    std::optional<Colorbar> colorbar;
    Camera camera;
    Float3 lo{};  // bounds of everything drawn, um
    Float3 hi{};
    bool axis_triad = true;
};

// --- fields -------------------------------------------------------------------------------------

enum class Field : std::uint8_t {
    potential,         // V
    electrons,         // n, cm^-3
    holes,             // p, cm^-3
    conduction_band,   // E_c, eV
    valence_band,      // E_v, eV
    electron_fermi,    // E_Fn, eV
    hole_fermi,        // E_Fp, eV
    temperature,       // K (the device's temperature without the electrothermal model)
    net_doping,        // N_D - N_A, cm^-3
    donors,            // N_D, cm^-3
    acceptors,         // N_A, cm^-3
    // Unit 45: generation-recombination rates (results::RateFields), cm^-3 s^-1; an error when the
    // result holds none (not stored, or the mechanism off).
    recombination,     // R: net R - G of every mechanism (electron rows), signed
    srh,               // SRH net recombination, signed
    auger,             // Auger net recombination, signed
    radiative,         // radiative net recombination, signed
    surface,           // the interfaces' effective volume rate (electron rows), signed
    impact_generation, // Gii: impact-ionization generation, positive
    btbt_generation,   // Gbtbt: band-to-band generation (electron rows), positive
    last = btbt_generation,
};

struct FieldInfo {
    std::string_view key;    // as nitcad plot names it: "potential", "n", "p", "Ec", ...
    std::string_view title;  // UTF-8, with its unit
    ScaleKind default_scale;
};

[[nodiscard]] FieldInfo field_info(Field f) noexcept;
[[nodiscard]] std::optional<Field> find_field(std::string_view key) noexcept;
// The field's value at every node. Errors (invalid_input): the bias point's fields or bands not
// one per mesh node.
[[nodiscard]] std::expected<std::vector<double>, base::Error> field_values(
    const device::Device& device, const results::BiasPoint& point, Field f);

// --- region colors (VR8) -------------------------------------------------------------------------

enum class Material : std::uint8_t {
    insulator,   // pale, translucent
    n_heavy,     // N_D - N_A >= heavy: red
    n_light,     // 0 < N_D - N_A < heavy: salmon
    p_channel,   // channel <= N_A - N_D < heavy: blue
    p_heavy,     // N_A - N_D >= heavy: dark blue
    substrate,   // |N_D - N_A| below channel (p) or 0: beige
    metal,       // electrodes and contacts: dark grey
};

struct Palette {
    double heavy_cm3 = 1e19;
    double channel_cm3 = 1e17;
    float insulator_alpha = 0.35f;
};

[[nodiscard]] Material material_class(const device::Device& device, NodeId node,
                                      const Palette& palette = {}) noexcept;
[[nodiscard]] Rgba material_color(Material m, const Palette& palette = {}) noexcept;
[[nodiscard]] std::string_view material_name(Material m) noexcept;

// --- views ---------------------------------------------------------------------------------------

struct ViewOptions {
    Palette palette;
    bool contacts = true;           // electrode slabs and contact pads
    bool region_outlines = true;    // 2D: interfaces and the boundary
    bool mesh_edges = false;        // element edges (2D) or boundary-face edges (3D)
    bool transparent_insulators = true;
    Projection projection = Projection::orthographic;
    float azimuth_deg = 35.0f;
    float elevation_deg = 25.0f;
};

// The regions (VR2, VR8): 2D a map of the control-volume pieces; 3D the boundary surface and the
// interfaces between regions, insulators translucent, with contacts. Errors (invalid_input): a
// geometry whose mesh is not the device's, or a 1D or graph-only geometry.
[[nodiscard]] std::expected<FieldScene, base::Error> solid_view(const device::Device& device,
                                                               const Geometry& geometry,
                                                               const ViewOptions& options = {});
// The mesh with its regions: element edges (2D) or boundary-face edges (3D) over the solid view.
[[nodiscard]] std::expected<FieldScene, base::Error> mesh_view(const device::Device& device,
                                                              const Geometry& geometry,
                                                              const ViewOptions& options = {});

struct FieldOptions {
    ViewOptions view;
    Colormap map = Colormap::turbo;
    std::optional<ColorScale> scale;  // fitted to the values when empty
    bool contours = false;            // lines at the colorbar's major ticks
    std::string title;                // the colorbar's title
};

// A nodal field (VR2, VR3): 2D a map on the elements; 3D on the boundary surface (lit), with the
// insulators translucent and the contacts. Errors (invalid_input): values not one per node; a
// geometry as in solid_view; no value to fit a scale to.
[[nodiscard]] std::expected<FieldScene, base::Error> field_view(const device::Device& device,
                                                               const Geometry& geometry,
                                                               std::span<const double> values,
                                                               const FieldOptions& options = {});
// Slice planes through a 3D device (VR2): the field on each section, inside the translucent outline
// of the device's boundary. Errors as field_view, and (invalid_input) a 2D geometry or a plane
// that misses the device.
[[nodiscard]] std::expected<FieldScene, base::Error> slice_view(const device::Device& device,
                                                               const Geometry& geometry,
                                                               std::span<const double> values,
                                                               std::span<const Plane> planes,
                                                               const FieldOptions& options = {});

// --- plots ---------------------------------------------------------------------------------------

enum class AxisScale : std::uint8_t { linear, log10 };

struct Axis {
    std::string label;  // with its unit
    AxisScale scale = AxisScale::linear;
    double lo = 0.0;
    double hi = 1.0;
    std::vector<Tick> ticks;
};

struct Series {
    std::string name;
    std::vector<double> x;
    std::vector<double> y;
    Rgb color{};
    bool markers = false;
    // Points flagged below resolution (Unit 24's flag): drawn with an open marker.
    std::vector<std::uint8_t> flagged;
};

struct PlotScene {
    std::string title;
    Axis x;
    Axis y;
    std::vector<Series> series;
    bool legend = true;
};

struct SeriesInput {
    std::string name;
    analysis::Curve curve;
};

// Curves on shared axes (VR4): ranges fitted to the data (a log axis to the decades holding the
// positive values; |y| when every y is negative, as for a reverse current), nice ticks, points
// with |y| at or below their resolution flagged. Colors cycle through a fixed qualitative set.
// Errors (invalid_input): no series; a log axis with no positive value.
[[nodiscard]] std::expected<PlotScene, base::Error> curve_plot(std::span<const SeriesInput> series,
                                                              std::string x_label,
                                                              std::string y_label,
                                                              AxisScale x_scale = AxisScale::linear,
                                                              AxisScale y_scale = AxisScale::linear,
                                                              std::string title = {});
// A nodal field along a 1D device against x in um. Errors (invalid_input): not a 1D mesh, or values
// not one per node.
[[nodiscard]] std::expected<PlotScene, base::Error> profile_plot(const mesh::Mesh& mesh,
                                                                std::span<const double> values,
                                                                std::string y_label,
                                                                AxisScale y_scale = AxisScale::linear);
// The band diagram of a 1D device (E_c, E_v, E_Fn, E_Fp in eV against x in um; insulator nodes are
// gaps). Errors (invalid_input): not a 1D mesh, or bands not one per node.
[[nodiscard]] std::expected<PlotScene, base::Error> band_plot(const mesh::Mesh& mesh,
                                                             const results::BandDiagram& bands);

}  // namespace NiTCAD::render
