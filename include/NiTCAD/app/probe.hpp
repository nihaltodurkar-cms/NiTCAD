// Probe and measure (Unit 26, DECISIONS.md VA8): the mesh node under a point of the view, its position
// and field value, and the distance between two probed nodes.
//
// The nodes that can be picked are those the view shows: every node of a 2D map; in 3D the nodes of
// the boundary surface (the tetrahedra's outer faces), of which the front-most within the pick radius
// wins. A 1D device has none (its views are plots).
#pragma once

#include <array>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/app/view_settings.hpp"
#include "NiTCAD/app/workspace.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/render/scene.hpp"

namespace NiTCAD::app {

struct ProbeHit {
    mesh::NodeId node = 0;
    std::array<double, 3> position_um{};
    double value = 0.0;          // the field's value at the node (NaN when no values are given)
    float distance_px = 0.0f;    // from the pick point
};

class Prober {
public:
    explicit Prober(const Workspace& w);

    [[nodiscard]] std::span<const mesh::NodeId> candidates() const noexcept { return candidates_; }
    // The node drawn nearest the pixel (x, y) of the scene drawn at width x height, within the radius:
    // the nearest in 2D, the front-most in 3D. `values` has one entry per node, or is empty.
    [[nodiscard]] std::optional<ProbeHit> pick(const render::FieldScene& scene, int width, int height, float x,
                                               float y, std::span<const double> values = {},
                                               float radius_px = 12.0f) const;
    // A point of the device under the pixel, um (Unit 26b, a cut line's ends): in 2D the point itself, anywhere
    // in the view's plane; in 3D the front-most surface node within the radius (as pick).
    [[nodiscard]] std::optional<std::array<double, 3>> point_at(const render::FieldScene& scene, int width, int height,
                                                                float x, float y, float radius_px = 12.0f) const;
    // Region visibility (Unit 26c, YD4): the nodes of these regions (ids) are not picked; in 3D the shown nodes
    // of the faces they shared with shown regions are, as the view draws those faces. Empty: every region.
    void hide(std::span<const std::int32_t> regions);

private:
    const Workspace* workspace_ = nullptr;
    const mesh::Mesh* mesh_ = nullptr;
    int dimension_ = 0;
    std::vector<mesh::NodeId> all_;         // with every region shown
    std::vector<mesh::NodeId> candidates_;  // less the hidden regions'
};

[[nodiscard]] double distance_um(const ProbeHit& a, const ProbeHit& b) noexcept;

// --- hover readout (Unit 26e, DECISIONS.md HD1, HD2) ---------------------------------------------------

// The readout of a field view's hit, from the view's title "<symbol> (<unit>)":
//   "<symbol> = <value> <unit>" ("<symbol>: no value" when the value is not finite)
//   "x <x>  y <y> µm" (with "  z <z>" before " µm" in 3D)
//   "node <n>"
// The value as {:.4g}, positions as {:.3g}.
[[nodiscard]] std::string hover_text(std::string_view title, const ProbeHit& hit, int dimension);

// The data point of a plot drawn nearest the pixel (x, y), over every series, within the radius.
struct CurveHit {
    std::size_t series = 0;
    std::size_t index = 0;
    double x = 0.0;
    double y = 0.0;
    float screen_x = 0.0f;  // where the point is drawn
    float screen_y = 0.0f;
    float distance_px = 0.0f;
};
[[nodiscard]] std::optional<CurveHit> pick_curve(const render::PlotScene& scene, int width, int height, float x,
                                                 float y, float radius_px = 12.0f);
// The readout of a curve hit: the series name when the plot has more than one, then "<x label> = <x>" and
// "<y label> = <y>", each value as {:.4g}.
[[nodiscard]] std::string curve_text(const render::PlotScene& scene, const CurveHit& hit);

// --- band surface readout (Unit 26o, BSD7) ---------------------------------------------------------------

// The band surface drawn under a pixel: the nearest shown surface the pixel's ray meets (a ray-triangle pick on
// its drawn batches), the device point under it (on a 3D view, in the plane), and every shown level there read
// in double as the surface draws them (render::band_energy_at; nullopt where a level has no value or is not
// shown), with the region whose control volume holds the point.
struct BandHit {
    int level = 0;  // Ec, Ev, EFn, EFp (render::BandLevel order)
    std::array<std::optional<double>, 4> energy_eV{};
    std::array<double, 3> position_um{};
    std::string region;
};
// Nothing off every surface, or when the settings are not a band surface view.
[[nodiscard]] std::optional<BandHit> band_pick(const Workspace& w, const ViewSettings& s,
                                               const render::FieldScene& scene, int width, int height, float x,
                                               float y);
// "<level> = <E> eV" for the level under the cursor, then the other shown levels in order ("<level>: no value"
// where none), "x <x>  y <y> µm" (with "  z <z>" in 3D), and the region. E as {:.4g}, positions as {:.3g}.
[[nodiscard]] std::string band_hover_text(const BandHit& hit, int dimension);

}  // namespace NiTCAD::app
