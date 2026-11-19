// Comparing two runs by tolerances (ARCHITECTURE.md 6.8, Unit 60; DECISIONS.md XC1-XC2).
//
// Agreement between machines, builds and implementations is defined by tolerances, never by
// bit-identity digests (6.8). compare_sweeps says how far apart two sweeps of the same input are:
// for each quantity, the largest difference, its ratio to the allowance, and where it occurs.
// OLD / NEW / REASON: OLD, the 11.3 row asked for "cross-machine checksums". NEW, a tolerance
// comparison. REASON, a digest of rounded values flips whenever a value lies near a rounding
// boundary, and says nothing about how far apart two results are.
//
// Structure is compared exactly first (contacts, points, nodes, the bias per point, which optional
// fields are present, whether the run stopped); a mismatch makes the sweeps not comparable, with
// the reason. The input identity and the mesh digest are reported, not gated (XQ6): mesh
// generation may call CRT math that takes other paths on another CPU.
#pragma once

#include <cstddef>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/results/run.hpp"

namespace NiTCAD::analysis {

// The allowances (XC2). For a current, |a - b| <= relative * max(|a|, |b|) + max(resolution_a,
// resolution_b): each run's own resolution is its floor (BiasPoint::terminal_current_resolution).
struct CompareTolerances {
    double relative = 1e-6;           // currents, charges and densities
    double charge_floor = 1e-30;      // absolute, in C cm^(D-3)
    double potential_V = 1e-6;        // absolute
    double density_floor = 1e-8;      // of the point's largest density, as Newton's update_size
    double temperature_K = 1e-6;      // absolute; lattice and carrier temperatures
};

// One quantity over every point (and node or contact): the largest difference and the largest
// ratio of a difference to its allowance, where that ratio occurs.
struct QuantityComparison {
    std::string name;            // "terminal_current", "potential_V", ...
    double max_difference = 0.0;  // in the quantity's unit
    double max_ratio = 0.0;       // difference / allowance; <= 1 passes
    std::size_t point = 0;        // of the largest ratio
    std::size_t index = 0;        // its node or contact
    std::size_t compared = 0;     // how many values were compared
    [[nodiscard]] bool pass() const noexcept { return max_ratio <= 1.0; }
};

enum class CompareOutcome { match, differ, not_comparable };

struct SweepComparison {
    CompareOutcome outcome = CompareOutcome::match;
    std::string reason;  // why the sweeps are not comparable; empty otherwise
    bool identity_equal = false;
    bool mesh_digest_equal = false;
    std::vector<QuantityComparison> quantities;  // only quantities both sweeps hold
};

[[nodiscard]] SweepComparison compare_sweeps(const results::Sweep& a, const results::Sweep& b,
                                             const CompareTolerances& tolerances = {});

[[nodiscard]] std::string_view to_string(CompareOutcome outcome) noexcept;

}  // namespace NiTCAD::analysis
