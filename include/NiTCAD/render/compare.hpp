// Run comparison (Unit 26b; DECISIONS.md XD6): the difference B - A of a field between two runs on one
// mesh, at a bias point both have.
//
// Units: the potential in mV and the band edges and quasi-Fermi levels in meV (the fields are energies in
// eV; a meV is the mV of the potential); the densities n, p, N_D and N_A as the difference of their log10
// (decades: +1 is ten times run A's), NaN where either is not positive; every other field linearly in its
// own unit. A difference is drawn on the cool-warm map, symmetric about 0.
#pragma once

#include <cstddef>
#include <expected>
#include <optional>
#include <span>
#include <string>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/render/color.hpp"
#include "NiTCAD/render/geometry.hpp"
#include "NiTCAD/render/scene.hpp"
#include "NiTCAD/results/solution.hpp"

namespace NiTCAD::render {

// The first of run B's points whose bias equals run A's point: as many contacts, each bias within
// 1e-9 V. None when B has no such point. (A transient's snapshots share their bias: the first is taken.)
[[nodiscard]] std::optional<std::size_t> matching_point(const results::BiasPoint& a,
                                                        std::span<const results::BiasPoint> b) noexcept;

struct DifferenceField {
    std::vector<double> values;  // B - A per node, in the unit of the title
    std::string title;           // e.g. "Δψ (mV)", "Δlog₁₀ n (decades)"
};

// B - A of a field: field_values of run A on its device and of run B on its own (the doping may differ),
// both on run A's geometry (one mesh). Errors as field_values.
[[nodiscard]] std::expected<DifferenceField, base::Error> difference_field(const device::Device& device_a,
                                                                         const device::Device& device_b,
                                                                         const Geometry& geometry,
                                                                         const results::BiasPoint& a,
                                                                         const results::BiasPoint& b, Field f);

// The linear scale symmetric about 0 over the finite values: -m to m for the largest |value| m, -1 to 1
// when every value is 0 or none is finite (a run compared with itself).
[[nodiscard]] ColorScale difference_scale(std::span<const double> values);

}  // namespace NiTCAD::render
