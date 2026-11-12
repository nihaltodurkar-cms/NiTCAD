// The intelligent decision layer of adaptive refinement (Unit 17, slice 5; DECISIONS.md Unit 17
// ID1-ID7, ID9, AD12, AR5).
//
// The layer decides only where nodes are spent (ID1, OA6): the estimate, the stop rule, AD11's
// verification and the certificate are computed from solves alone and read nothing here.
//
// A priori prediction (ID2), before any solve, from the device and the analysis's biases:
// - junctions: every edge between semiconductor nodes of opposite net doping gives a junction
//   point (where linear interpolation of N_D - N_A crosses zero) with the edge's direction as its
//   normal. N_1, N_2: the largest |N| of each side's sign within the depletion width of the point
//   (a fixed point from the end nodes' values, three iterations); N_eff = N_1 N_2 / (N_1 + N_2);
//   V_bi = V_T ln(N_1 N_2 / n_i^2) (>= 0); W = sqrt(2 eps (V_bi + V) / (q N_eff)) with V the
//   largest difference between non-gate contact biases over the analysis; the target size is
//   L_D / 2, L_D = sqrt(eps V_T / (q N_eff)), for nodes within W of the point.
// - insulator interfaces under gates: a gate contact's nodes (its oxide thickness and
//   permittivity), and semiconductor nodes next to an insulator when the device has electrode
//   contacts (t_ox the distance to the nearest electrode node, the insulator's permittivity): the
//   surface field E_s = eps_ox V / (eps_s t_ox), V the largest difference between the gate (or an
//   electrode) and the other contacts over the analysis, gives the classical inversion length
//   lambda = V_T / E_s (with MLDA on, the larger of it and MLDA's electron wavelength); the target
//   size is lambda for nodes within 5 lambda; and the depletion region under it (AR5): L_D / 2,
//   L_D at the node's |N|, for nodes within W_d = sqrt(2 eps_s psi / (q N)), psi = min(2 phi_F, V).
//   The normal points into the semiconductor (the mean direction to the node's semiconductor
//   neighbours off the interface).
// - contacts: none. Impact ionization and band-to-band tunnelling add nothing beyond their
//   junction's target (AR5).
// A node's target is the smallest over the features covering it, with that feature's normal. On a
// tensor or swept axis an edge needs refinement when its length projected on the normal exceeds the
// target at either end (a junction or an interface is resolved across, not along, itself); on a
// generated 2D mesh or cross-section the target enters the size field as is (isotropic).
//
// Bias prediction (ID3): from one coarse sweep on the start mesh, per point the united indicator
// excess; the points are the last one and those whose field changed most from the previous point
// (the weighted 2-norm of the difference, weights 1), up to max_points, in sweep order.
#pragma once

#include <array>
#include <cstddef>
#include <expected>
#include <functional>
#include <limits>
#include <span>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/solve/adapt_indicators.hpp"
#include "NiTCAD/solve/bias.hpp"

namespace NiTCAD::solve {

enum class Policy : unsigned char { basic, intelligent };

[[nodiscard]] constexpr const char* to_string(Policy p) noexcept {
    return p == Policy::basic ? "basic" : "intelligent";
}

// Per node: a target size [cm] (infinity where nothing is predicted) and the unit normal of the
// feature that gave it (all zero: isotropic).
struct SizeTargets {
    std::vector<double> size;
    std::vector<std::array<double, 3>> normal;
};

struct Prediction {
    SizeTargets targets;
    std::size_t junction_points = 0;
    std::size_t interface_points = 0;
    std::size_t predicted_nodes = 0;  // nodes with a finite target
};

// ID2. Errors (invalid_input): points not one bias per contact.
[[nodiscard]] std::expected<Prediction, base::Error> predict_sizes(const device::Device& device,
                                                                   std::span<const std::vector<double>> points,
                                                                   const BiasOptions& options);

// Whether edge e of the device's mesh is longer than the targets ask (projected on the normal).
[[nodiscard]] bool edge_exceeds(const device::Device& device, const SizeTargets& targets, std::size_t e);

// ID3: the adaptation points from per-point indicator fields (one per sweep point, same size),
// among the eligible points (empty: all; the loop excludes points where a goal is at or below its
// current resolution, which cannot carry a relative tolerance, AR5): the last eligible point and
// the eligible points of largest change. Errors (invalid_input): no fields, fields of different
// sizes, max_points 0, eligibility not one flag per point, no eligible point.
[[nodiscard]] std::expected<std::vector<std::size_t>, base::Error> choose_points(
    std::span<const IndicatorField> fields, std::size_t max_points, std::span<const char> eligible = {});

// The layer's options (adapt.hpp's AdaptOptions::intelligent).
struct PolicyOptions {
    // ID3: with AdaptOptions::at empty, up to this many points chosen from a coarse sweep; 0 keeps
    // AdaptOptions' default points (the first and last).
    std::size_t max_points = 0;
    // A relative fall the estimate must make for a proposal to count as lowering it (ID1, AR5).
    double miss_margin = 0.1;
    // A hand-written prediction in place of the layer's own (ID2 at the first mesh, and every later
    // proposal): per-node targets on the given device. Used for AC17's wrong prediction; empty: the
    // layer's.
    std::function<SizeTargets(const device::Device&)> prediction;
};

}  // namespace NiTCAD::solve
