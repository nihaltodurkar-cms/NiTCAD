// Analytic doping profiles sampled at mesh nodes (Unit 16; DECISIONS.md M15).
//
// A profile is a depth shape along one axis times a lateral window along the others:
//   N(x) = peak * depth(x_axis) * prod over the other axes b of window(x_b),
// - depth: constant (1), Gaussian exp(-(d - center)^2 / (2 width^2)), or erfc
//   0.5 erfc((d - center) / width) (a diffusion from a surface at `center`, decreasing in d);
// - window along each other axis: [lo, hi] with an erfc fall-off of length `lateral` on both
//   sides, 0.5 [erfc((lo - x) / lateral) - erfc((hi - x) / lateral)]; lateral = 0 is the box
//   (1 inside, 0 outside, 1/2 on its faces). An axis whose window is unbounded is skipped.
// Profiles of one species add. Sampling is point sampling at nodes, as the devices built by hand
// before Unit 16 do; insulator nodes get zero.
#pragma once

#include <array>
#include <cstdint>
#include <expected>
#include <limits>
#include <span>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/mesh.hpp"

namespace NiTCAD::device {

enum class Dopant : std::uint8_t { donor, acceptor };
enum class DepthShape : std::uint8_t { constant, gaussian, erfc };

struct DopingProfile {
    Dopant dopant = Dopant::donor;
    double peak_cm3 = 0.0;  // finite, >= 0
    DepthShape shape = DepthShape::constant;
    int axis = 0;            // the depth axis, 0..dimension-1
    double center = 0.0;     // [cm] the Gaussian's centre, or the erfc's 50 % point
    double width = 0.0;      // [cm] sigma or the erfc length; > 0 unless constant
    // The lateral window per axis (the depth axis's entry is ignored), [cm].
    std::array<double, 3> lo{-std::numeric_limits<double>::infinity(),
                             -std::numeric_limits<double>::infinity(),
                             -std::numeric_limits<double>::infinity()};
    std::array<double, 3> hi{std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::infinity()};
    double lateral = 0.0;  // [cm] erfc fall-off of the window, >= 0
};

// The profile's value at p [cm^-3]. Precondition: the profile is valid (check_profile).
[[nodiscard]] double evaluate(const DopingProfile& profile, const mesh::Point& p);

struct Doping {
    std::vector<double> donors;
    std::vector<double> acceptors;
};

// The sum of the profiles at every node, per species; zero on nodes whose region is an insulator.
// Errors (invalid_input, index = the profile, or the node for node_region): a peak not finite and
// >= 0; an axis outside the mesh dimension; a width not finite and positive for a Gaussian or erfc;
// a non-finite center or lateral, lateral < 0, or a window with lo > hi or a NaN bound;
// node_region not one per node or a region out of range.
[[nodiscard]] std::expected<Doping, base::Error> sample(std::span<const DopingProfile> profiles,
                                                        const mesh::Mesh& mesh,
                                                        std::span<const RegionId> node_region,
                                                        std::span<const Region> regions);

}  // namespace NiTCAD::device
