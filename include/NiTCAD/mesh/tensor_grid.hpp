// Tensor-product (structured) mesh producer for 1D, 2D and 3D.
//
// Each axis is a strictly increasing list of node coordinates (cm), at least two. The box method
// on this grid follows the legacy mesh.py / mesh2d.py / mesh3d.py and device2d/3d.py:
// - control-volume width along an axis: half of each adjacent spacing (one half at an end node);
// - node volume: the product of its widths along every axis;
// - edge along axis a: length is the spacing; coupling area is the product of the end nodes'
//   widths along the other axes (1 in 1D);
// - boundary patches "x_min", "x_max", "y_min", "y_max", "z_min", "z_max" (as far as the
//   dimension goes): the nodes on that face, each with the product of its widths along the
//   face's axes (1 in 1D).
// The result is an ordinary Mesh: nothing above this file knows it came from a grid.
#pragma once

#include <expected>
#include <span>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/mesh/mesh.hpp"

namespace NiTCAD::mesh {

// Errors: degenerate_mesh if an axis has fewer than two nodes, a non-finite coordinate, or is not
// strictly increasing (the context index is the position on that axis); resource_exhausted if the
// node count does not fit NodeId; otherwise those of Mesh::from_parts.
[[nodiscard]] std::expected<Mesh, base::Error> make_tensor_grid(std::span<const double> x);
[[nodiscard]] std::expected<Mesh, base::Error> make_tensor_grid(std::span<const double> x,
                                                                std::span<const double> y);
[[nodiscard]] std::expected<Mesh, base::Error> make_tensor_grid(std::span<const double> x,
                                                                std::span<const double> y,
                                                                std::span<const double> z);

// An axis for a material interface at `position` (Unit 15): the assemblers put an interface at the
// midpoint of the edge joining two regions, so the axis needs two nodes at position -+ spacing / 2
// and none between. Returns the axis with every node closer than `spacing` to the position
// removed and those two inserted (a region is then assigned by x < position). Errors
// (invalid_input): the axis not strictly increasing or not finite; spacing not finite and
// positive; the position not inside the axis by more than spacing.
[[nodiscard]] std::expected<std::vector<double>, base::Error> straddle_interface(
    std::vector<double> axis, double position, double spacing);

}  // namespace NiTCAD::mesh
