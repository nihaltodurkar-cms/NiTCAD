// Vector fields reconstructed from a result (Unit 25b; DECISIONS.md DS1-DS3, BQ1): the electric
// field and the current densities, per element and averaged to the nodes. Nothing is solved again:
// they are read from the potential and the stored edge currents.
//
// - E = -grad psi of the linear interpolant on each element (segment, triangle, tetrahedron), in
//   V/cm, insulators included; exact for a linear potential.
// - Jn and Jp per element are the vectors whose projections on the element's mesh edges best match
//   the edge densities I / A (the stored edge current over the edge's coupling area, in A/cm^2 in
//   every dimension), by least squares over the edges with a positive coupling area and both ends
//   in a semiconductor: the solver's own Scharfetter-Gummel currents, exact for a uniform current.
//   An element without D independent such edges has none. J = Jn + Jp (conduction current only).
// - At a node, the average of its elements' vectors weighted by the volume of the node's
//   control-volume piece in each (geometry.hpp's region rule): in a tensor cell the sub-box of the
//   corner, which takes binomial(D, m) / 2^D of a Kuhn simplex's volume for the simplex's vertex m
//   steps along its path; in a simplex the barycentric dual, 1 / (D + 1) of it for each vertex.
//   A node with no element that has the vector gets NaN.
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <span>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/render/geometry.hpp"
#include "NiTCAD/results/solution.hpp"

namespace NiTCAD::render {

enum class VectorField : std::uint8_t {
    electric,          // E, V/cm
    current,           // J = Jn + Jp, A/cm^2
    electron_current,  // Jn, A/cm^2
    hole_current,      // Jp, A/cm^2
};

// One vector per element (components beyond the mesh's dimension 0) and whether it has one.
struct ElementVectors {
    std::vector<Point> value;
    std::vector<std::uint8_t> defined;
};

// Errors (invalid_input): a geometry whose mesh is not the device's or that has no elements (2D or
// 3D nodes and edges only); the point's potential not one value per node (E), or its edge currents
// not one per mesh edge (J, Jn, Jp). A quasi-static result's edge currents are zero (bias.hpp), and
// so are its J, Jn and Jp.
[[nodiscard]] std::expected<ElementVectors, base::Error> element_vectors(
    const device::Device& device, const Geometry& geometry, const results::BiasPoint& point,
    VectorField field);

// Element e's volume (length in 1D, area in 2D) in the mesh's unit, and the fraction of it in each of
// its nodes' control volumes, in the order of Geometry::element (fractions sum to 1).
[[nodiscard]] double element_volume(const Geometry& geometry, std::size_t e) noexcept;
[[nodiscard]] std::array<double, 4> control_volume_shares(const Geometry& geometry,
                                                          std::size_t e) noexcept;

// The control-volume-weighted average at every node; NaN components where no element has one.
[[nodiscard]] std::vector<Point> nodal_vectors(const Geometry& geometry, const ElementVectors& v);

// |v| per entry (NaN stays NaN).
[[nodiscard]] std::vector<double> magnitudes(std::span<const Point> v);

}  // namespace NiTCAD::render
