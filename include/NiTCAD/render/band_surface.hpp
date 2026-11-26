// The band surface (Unit 26o; DECISIONS.md BSD1-BSD4): the band edges and quasi-Fermi levels of a solved device as
// energy-height surfaces over its 2D mesh or over one plane of its 3D mesh.
//
// Steps are kept (BSD2, BSQ1 (a)): the base polygons are split into the control-volume pieces every view draws
// regions with (Geometry::owner_pieces), and a piece reads a level only from the nodes of its owner's material
// in the element, so a heterojunction's offset is a wall on the region boundary rather than a slope. Insulator
// nodes carry no bands (NaN): their pieces are a gap. Nothing here computes a band: the levels are the solver's
// (results::BandDiagram), and a node's material is compared only for identity (the solver's own test of a
// heterointerface, equal semiconductor parameters).
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/render/color.hpp"
#include "NiTCAD/render/cut_line.hpp"
#include "NiTCAD/render/geometry.hpp"
#include "NiTCAD/render/scene.hpp"
#include "NiTCAD/results/solution.hpp"

namespace NiTCAD::render {

enum class BandLevel : std::uint8_t { conduction, valence, electron_fermi, hole_fermi };

struct BandSurfaceOptions {
    std::array<bool, 4> levels{true, true, false, false};  // Ec, Ev, EFn, EFp (BandLevel order)
    std::optional<Plane> plane;                             // the section of a 3D device (cm); none in 2D
    double exaggeration = 1.0;                              // [0.1, 10]: the heights' scale on the default
    float opacity = 1.0f;                                   // [0.2, 1]
    std::optional<std::pair<double, double>> range_eV;      // lo < hi; empty: fitted to the drawn levels
    bool by_level = false;                                  // BSQ4: colour by energy unless set
    Colormap map = Colormap::viridis;
    ViewOptions view;                                       // projection and camera angles
    // BSD9: a cut line (render::cut_line) drawn on each shown level at its heights, as a polyline over the surfaces
    // (stepping where they step); left out when it does not lie in the base plane.
    std::optional<CutLine> cut;
};

// The scene, with what reading it back needs: a height is um_per_eV * E in scene y.
struct BandSurfaceScene {
    FieldScene scene;
    double um_per_eV = 0.0;
    std::pair<double, double> range_eV{};
    std::size_t dropped = 0;  // pieces not drawn: a level not finite at one of their points, or of zero area
};

// The levels as energy-height surfaces (BSD1-BSD4). The base is the 2D device's elements, or the 3D device's
// section by options.plane (Geometry::slice); its axes (u, v) are scene x and z in um: in 2D x and y, for a plane
// normal to an axis the other two in order, else u the plane's projection of the first world axis not along the
// normal and v = n x u. Scene y is um_per_eV * E: by default 0.6 times the base's larger extent over the
// range's span (1 eV for a single value, the range then widened by 0.5 eV each side), times exaggeration. One
// batch per level shown, named Ec, Ev, EFn, EFp, with walls between materials in the same batch (vertical,
// normals in the base plane); coloured by energy on one linear scale (a colorbar "E (eV)") or each in the band
// diagram's colour. The range clips the pieces. scene.axes names the axes ("x (µm)", "E (eV)", ...).
// Errors (invalid_input): a 1D device; a plane in 2D or none in 3D; a level not one value per node; no level
// shown; exaggeration or opacity out of range; range lo >= hi or not finite; the plane missing the device;
// nothing finite to draw.
[[nodiscard]] std::expected<BandSurfaceScene, base::Error> band_surface_view(const device::Device& device,
                                                                            const Geometry& geometry,
                                                                            const results::BandDiagram& bands,
                                                                            const BandSurfaceOptions& options = {});

// A level at point x (cm) in element e as the surface draws it (BSD2), in double: nullopt where x's control-volume
// owner is an insulator or the value is not finite. The owner is the element's node whose piece holds x (a
// tensor cell's corner sub-box, a simplex's largest barycentric weight).
[[nodiscard]] std::optional<double> band_energy_at(const device::Device& device, const Geometry& geometry,
                                                   const results::BandDiagram& bands, BandLevel level,
                                                   std::size_t element, const Point& x);

// One energy range over a run's points (BSD8): the lowest and highest finite value of the levels set, over the
// semiconductor nodes of every point (a level of the wrong size skipped), widened by 0.5 eV each side about its
// middle when it is one value (span under 1e-9 eV). Nothing when no value is finite. The nodes of hidden regions
// (ViewOptions::hidden_regions) are left out, as the surface leaves them out.
[[nodiscard]] std::optional<std::pair<double, double>> band_energy_range(const device::Device& device,
                                                                         std::span<const results::BiasPoint> points,
                                                                         const std::array<bool, 4>& levels,
                                                                         std::span<const std::int32_t> hidden_regions = {});

// The plane's outline on a device view (BSD6): the outer edges of the device's section by the plane (each edge of
// Geometry::slice's polygons that no other polygon shares), as one on-top line batch "band plane" in dark orange.
// Nothing when the plane misses the device.
void add_plane_outline(FieldScene& scene, const Geometry& geometry, const Plane& plane);

// The node whose control-volume piece of element `element` holds x (cm), as band_energy_at reads it (BSD7: the
// readout's material).
[[nodiscard]] NodeId band_owner(const Geometry& geometry, std::size_t element, const Point& x);

// One key per mesh node: -1 on insulator nodes; semiconductor nodes share a key exactly when their regions'
// semiconductor parameters compare equal. Keys count from 0 in region order.
[[nodiscard]] std::vector<std::int32_t> band_material_keys(const device::Device& device);

// A level at point x (cm) of element `element`, read for the piece of `owner`'s control volume (BSD2): over the
// element's nodes when they are all of one material key (exactly interpolate()), else over the nodes of the
// owner's key with their weights renormalised, else (no weight on them where x lies) the owner's own value.
[[nodiscard]] double band_value(const Geometry& geometry, std::size_t element, const Point& x, NodeId owner,
                                std::span<const std::int32_t> keys, std::span<const double> level);

}  // namespace NiTCAD::render
