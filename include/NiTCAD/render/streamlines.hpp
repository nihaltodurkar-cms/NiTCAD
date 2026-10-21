// Streamlines of a nodal vector field and where they start (Unit 25b; DECISIONS.md DS5, DS6, BQ2).
//
// The field is the nodal vectors (vectors.hpp) interpolated linearly on each element, the same field
// the maps draw. A line follows its direction at unit speed (dx/ds = v / |v|), so lines through
// fields spanning many decades step by length, not by magnitude: adaptive Dormand-Prince 5(4) with
// the local error per step held below `tolerance` times the device's size, and no step longer than
// half the local element's size (the drawn polyline follows the field). A line stops at the
// domain boundary (within the step control's resolution; contacts lie on it), where the field is
// undefined (a NaN node: J in an insulator) or vanishes (|v| at most 1e-12 of the field's largest
// nodal magnitude) or reverses within a step (a sink: the line has reached the field's zero), at its
// length limit, or after max_steps steps.
//
// Seeds: on a contact, spaced so that each line carries an equal share of the current through it
// (BQ2: the line density then shows the current density), or evenly over it; evenly along a segment;
// evenly over a plane's section of the device.
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/render/geometry.hpp"
#include "NiTCAD/results/solution.hpp"

namespace NiTCAD::render {

// Finds the element holding a point: the hint, then the elements sharing a node with it, then a
// bucket grid over the elements' bounding boxes. A Locator refers to its geometry.
class Locator {
public:
    explicit Locator(const Geometry& geometry);

    struct Hit {
        std::size_t element = 0;
        std::array<double, 4> weights{};  // barycentric, in the order of Geometry::element
    };
    // The element holding x (barycentric weights at least -1e-10), or nothing outside the mesh.
    [[nodiscard]] std::optional<Hit> locate(const Point& x, std::optional<std::size_t> hint = {}) const;
    [[nodiscard]] const Geometry& geometry() const noexcept { return *g_; }
    // The largest extent of the mesh's bounding box (cm).
    [[nodiscard]] double size() const noexcept { return size_; }

private:
    [[nodiscard]] std::optional<Hit> inside(std::size_t e, const Point& x) const noexcept;

    const Geometry* g_;
    Point lo_{};
    std::array<std::size_t, 3> n_{1, 1, 1};
    double cell_ = 1.0;
    double size_ = 0.0;
    std::vector<std::uint32_t> bucket_start_, bucket_items_;
    std::vector<std::uint32_t> node_start_, node_items_;
};

enum class StreamlineStop : std::uint8_t {
    boundary,  // left the mesh
    no_field,  // undefined or vanishing field
    length,    // reached max_length
    steps,     // took max_steps
};

struct Streamline {
    std::vector<Point> points;  // cm, from the seed (both directions joined, backward part first)
    std::vector<double> magnitude;  // |v| at each point (NaN at a seed without a field)
    StreamlineStop forward = StreamlineStop::boundary;
    StreamlineStop backward = StreamlineStop::boundary;
};

struct StreamlineOptions {
    double tolerance = 1e-5;   // local error per step, relative to the device's size
    double max_length = 0.0;   // per direction, cm; 0: 20 times the device's size
    std::size_t max_steps = 20000;
    bool forward = true;       // trace along v
    bool backward = true;      // and against it
};

// One line through `seed`. `field` has one vector per mesh node. A seed outside the mesh gives a
// line of the seed alone.
[[nodiscard]] Streamline trace_streamline(const Locator& locator, std::span<const Point> field,
                                          const Point& seed, const StreamlineOptions& options = {});

// The current entering the device through each node of a contact, in its node order (A/cm^(3-D)):
// the stored electron and hole currents of the edges from the node to nodes outside the contact.
// Empty when the point holds no edge currents.
[[nodiscard]] std::vector<double> contact_node_currents(const device::Device& device,
                                                        const results::BiasPoint& point,
                                                        std::size_t contact);

// `count` seeds on a contact (2D: on its boundary segments; 3D: its boundary faces), placed at the
// midpoints of `count` equal shares of the distribution: with `node_current` (one per contact node,
// as contact_node_currents) each node's |current| spread evenly over its part of the contact; empty,
// the contact's length or area. A contact with no boundary segment or face gives its nodes'
// positions, evenly chosen. Seeds are moved 1e-9 of the device's size into the mesh.
[[nodiscard]] std::vector<Point> contact_seeds(const device::Device& device, const Locator& locator,
                                               std::size_t contact, std::size_t count,
                                               std::span<const double> node_current = {});
// --- arrows (Unit 25c, DC1) ---------------------------------------------------------------------

// Where an arrow stands and what it shows: the unit direction of the field there (zero when the
// field vanishes), |v|, the surface's unit normal (zero in 2D) and the interpolant, so any nodal
// field can be read at the site.
struct Arrow {
    Point position{};
    Point direction{};
    double magnitude = 0.0;
    Point normal{};
    Interpolant at;
};

// 2D: arrows at the points of a square grid over the mesh's bounding box, its spacing chosen for
// about `count` points inside the device (cell centres of the grid; points outside the mesh, or where
// the field is undefined or zero, give no arrow). The grid's spacing is returned in `spacing`.
// `excluded` (empty, or one flag per node) drops a site whose nearest node of its element (the
// largest weight) is flagged: an insulator's for currents, which it does not carry.
[[nodiscard]] std::vector<Arrow> grid_arrows(const Locator& locator, std::span<const Point> field,
                                             std::size_t count, double& spacing,
                                             std::span<const std::uint8_t> excluded = {});
// On surface polygons (3D boundary faces, or slice sections): `count` sites spread evenly by area (a
// stratified sample), each with its polygon's normal; spacing sqrt(area / count). `excluded` as above.
[[nodiscard]] std::vector<Arrow> surface_arrows(const Locator& locator, std::span<const Point> field,
                                                std::span<const Polygon> surface, std::size_t count,
                                                double& spacing, std::span<const std::uint8_t> excluded = {});

// `count` seeds evenly along the segment a-b (at the midpoints of equal parts).
[[nodiscard]] std::vector<Point> segment_seeds(const Point& a, const Point& b, std::size_t count);
// `count` seeds evenly by area over a plane's section of a 3D device (Geometry::slice).
[[nodiscard]] std::vector<Point> plane_seeds(const Geometry& geometry, const Plane& plane,
                                             std::size_t count);

}  // namespace NiTCAD::render
