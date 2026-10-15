// The M2 rule for assembled dual totals (Unit 16; DECISIONS.md M2). Internal to the mesh layer.
//
// Only totals are judged: an edge's coupling area, a node's volume, a node's area in a boundary
// patch, each summed over every simplex that contributes to it. A single simplex's piece is never
// judged (an obtuse simplex's pieces are legitimately negative). Each total is compared with
// tau * scale, the scale having the total's own dimension:
// - coupling area (cm^(D-1)): max(L^(D-1), sum of |pieces|), L the edge's length;
// - volume (cm^D): max(l^D, sum of |pieces|), l the mean length of the node's edges;
// - patch area (cm^(D-1)): max(l^(D-1), sum of |pieces|).
// The sum of |pieces| covers cancellation: rounding in a total scales with its pieces, not with L.
// Scaling the coordinates by lambda scales a total and its scale by the same lambda^power, so a
// decision does not depend on units (up to the rounding of the scaled coordinates).
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>

namespace NiTCAD::mesh::detail {

// About 4.5e3 ulps; the precedent of mesh.cpp's edge-length tolerance.
inline constexpr double dual_tau = 1e-12;

enum class TotalClass : std::uint8_t {
    keep,    // above tau * scale
    zero,    // within tau * scale of zero (both ends included): dropped
    refuse,  // below -tau * scale, or not finite
};

// A coupling area or a patch area.
[[nodiscard]] inline TotalClass classify_total(double total, double scale) noexcept {
    const double t = dual_tau * scale;
    if (!std::isfinite(total) || total < -t) return TotalClass::refuse;
    return total <= t ? TotalClass::zero : TotalClass::keep;
}

// A control volume: accepted only above tau * scale (a zero volume is refused, not dropped).
[[nodiscard]] inline bool accept_volume(double volume, double scale) noexcept {
    return std::isfinite(volume) && volume > dual_tau * scale;
}

[[nodiscard]] inline double total_scale(double length, int power, double magnitude) noexcept {
    return std::max(std::pow(length, static_cast<double>(power)), magnitude);
}

}  // namespace NiTCAD::mesh::detail
