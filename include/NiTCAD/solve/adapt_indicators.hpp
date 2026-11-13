// Refinement indicators and marking for adaptive refinement (Unit 17, slice 2; DECISIONS.md Unit 17
// AD2, AD3, AR1).
//
// Indicators are computed per mesh edge whose two ends are semiconductor nodes (an insulator's
// potential is linear between its charges, exact for the box method, so insulator edges carry
// none), each a raw dimensionless value r >= 0 normalized by its threshold t into an excess
// e = max(0, r / t - 1), 0 meaning resolved:
// - potential: |psi_a - psi_b| / V_T, V_T = k T / q at the edge's mean lattice temperature (the
//   device's temperature without the electrothermal model); Scharfetter-Gummel loses accuracy for
//   steps well above V_T. Default t = 1.
// - electron_density, hole_density: |ln n_a - ln n_b| / ln 10, the decades along the edge
//   (inversion layers, depletion edges); an edge with a non-positive density at either end is
//   skipped. Default t = 1.
// - doping: |L(N_a) - L(N_b)| / ln 10 with L(N) = sign(N) ln(max(|N|, n_i) / n_i), N = N_D - N_A and
//   n_i the node's intrinsic density at the device temperature (DECISIONS.md AR1): the decades of
//   net doping, counted down to n_i and back up across a junction, so an abrupt junction between
//   equal and opposite dopings is seen; on an edge of one sign above n_i it is AD2's
//   |Delta ln |N_D - N_A|| per decade. Needs no solve. Default t = 1.
// - generation (impact ionization and band-to-band tunnelling, only when their rates are stored):
//   per node, its share of the device's total generation, G_i V_i / sum_k G_k V_k, times the number
//   of semiconductor nodes (1 for a uniform spread), G the generation part of the impact and
//   band-to-band rates. An edge takes the larger of its ends. Default t = 10.
// A node's excess is the largest over its edges and kinds; its dominant kind is the one that gave
// it (the lowest kind on a tie).
//
// Marking (AD3, Doerfler-style bulk marking): node i is marked when e_i w_i > 0 and
// e_i w_i >= theta max_k(e_k w_k), w the goal weight (1 for feature-only marking). Deterministic:
// the rule does not depend on order.
#pragma once

#include <array>
#include <cstddef>
#include <expected>
#include <span>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/results/solution.hpp"

namespace NiTCAD::solve {

enum class Indicator : unsigned char { potential, electron_density, hole_density, doping, generation };
inline constexpr std::size_t indicator_count = 5;

[[nodiscard]] constexpr const char* to_string(Indicator k) noexcept {
    switch (k) {
        case Indicator::potential: return "potential";
        case Indicator::electron_density: return "electron_density";
        case Indicator::hole_density: return "hole_density";
        case Indicator::doping: return "doping";
        case Indicator::generation: return "generation";
    }
    return "?";
}

struct IndicatorOptions {
    // Which kinds are evaluated (index by Indicator). Kinds needing a state are skipped without one.
    std::array<bool, indicator_count> enabled{true, true, true, true, true};
    // Thresholds t (r / t - 1 is the excess).
    double potential_vt = 1.0;
    double density_decades = 1.0;
    double doping_decades = 1.0;
    double generation_share = 10.0;
};

// Per edge (in the mesh's edge order) and kind, the raw value r (0 where the kind does not apply).
struct EdgeIndicators {
    std::array<std::vector<double>, indicator_count> raw;
};

struct IndicatorField {
    std::vector<double> excess;           // per node, e >= 0
    std::vector<Indicator> dominant;      // per node, the kind that gave its excess
    std::array<std::size_t, indicator_count> nodes_with_excess{};  // per kind, nodes it gives e > 0
};

// Errors (invalid_input): a state whose fields are not one per node; options with a threshold not
// finite and positive. `state` may be null (doping only); `rates` may be null or empty.
[[nodiscard]] std::expected<EdgeIndicators, base::Error> edge_indicators(const device::Device& device,
                                                                         const results::NodeFields* state,
                                                                         const results::RateFields* rates,
                                                                         const IndicatorOptions& options = {});

[[nodiscard]] std::expected<IndicatorField, base::Error> compute_indicators(const device::Device& device,
                                                                          const results::NodeFields* state,
                                                                          const results::RateFields* rates,
                                                                          const IndicatorOptions& options = {});

// Folds several fields (one per adaptation point, AQ3) into one: per node the largest excess (its
// dominant kind from the field that gave it, the first on a tie). Errors (invalid_input): no fields,
// or fields of different sizes.
[[nodiscard]] std::expected<IndicatorField, base::Error> unite(std::span<const IndicatorField> fields);

struct Marking {
    std::vector<char> marked;  // per node
    double threshold = 0.0;    // theta max(e w)
    std::size_t count = 0;
    std::array<std::size_t, indicator_count> by_kind{};  // marked nodes by dominant kind
    std::vector<std::size_t> by_region;                  // marked nodes per device region
};

// Errors (invalid_input): theta not in (0, 1]; weights not empty and not one per node, or a weight
// negative or not finite; a field not one entry per device node.
[[nodiscard]] std::expected<Marking, base::Error> mark_nodes(const device::Device& device,
                                                             const IndicatorField& field,
                                                             std::span<const double> weights, double theta = 0.3);

}  // namespace NiTCAD::solve
