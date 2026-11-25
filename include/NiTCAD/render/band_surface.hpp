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

#include <cstddef>
#include <cstdint>
#include <span>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/render/geometry.hpp"

namespace NiTCAD::render {

// One key per mesh node: -1 on insulator nodes; semiconductor nodes share a key exactly when their regions'
// semiconductor parameters compare equal. Keys count from 0 in region order.
[[nodiscard]] std::vector<std::int32_t> band_material_keys(const device::Device& device);

// A level at point x (cm) of element `element`, read for the piece of `owner`'s control volume (BSD2): over the
// element's nodes when they are all of one material key (exactly interpolate()), else over the nodes of the
// owner's key with their weights renormalised, else (no weight on them where x lies) the owner's own value.
[[nodiscard]] double band_value(const Geometry& geometry, std::size_t element, const Point& x, NodeId owner,
                                std::span<const std::int32_t> keys, std::span<const double> level);

}  // namespace NiTCAD::render
