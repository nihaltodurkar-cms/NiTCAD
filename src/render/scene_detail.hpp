// Private to render: the scene builders' shared helpers (scene.cpp), for scenes built in other files (the band
// surface, Unit 26o). Each forwards to scene.cpp's own, so every existing view draws as before.
#pragma once

#include <expected>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/render/geometry.hpp"
#include "NiTCAD/render/scene.hpp"

namespace NiTCAD::render::detail {

// Appends a polygon (points in cm, drawn at their um) as a triangle fan facing `normal` (reversed when it faces
// away); `values` are colormap positions per point (empty: none, the batch's flat color).
void fan(TriangleBatch& batch, std::vector<Point> pts, const std::vector<float>& values, const Rgba& color,
         const Point& normal);

// Bounds and the default camera of a scene from everything in it.
void frame(FieldScene& s, const ViewOptions& o);

// The device and geometry share a mesh with elements, in 2D or 3D, and every hidden region is the device's.
[[nodiscard]] std::expected<void, base::Error> check_view(const device::Device& d, const Geometry& g,
                                                          const ViewOptions& o);

}  // namespace NiTCAD::render::detail
