// Helpers shared by the analysis sources (private to src/analysis).
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <string>
#include <utility>

#include "NiTCAD/analysis/curve.hpp"
#include "NiTCAD/base/error.hpp"

namespace NiTCAD::analysis::detail {

inline base::Error invalid(std::string message, std::optional<std::size_t> index = std::nullopt) {
    base::Error e{base::ErrorCode::invalid_input, std::move(message), std::nullopt};
    if (index) e.context = base::ErrorContext{.index = index, .value = std::nullopt};
    return e;
}

// Whether a point from first to last (inclusive) has |y| at or below its resolution.
inline bool below(const Curve& c, std::size_t first, std::size_t last) {
    if (c.resolution.empty()) return false;
    for (std::size_t k = first; k <= last; ++k) {
        if (std::abs(c.y[k]) <= c.resolution[k]) return true;
    }
    return false;
}

// Whether x lies between the bounds (either order), the bounds widened by 1e-9 of the distance
// between them (1e-9 of their size when equal): a sweep's biases are sums of steps, and 0.05 * 14
// is 0.7000000000000001.
inline bool between(double x, double from, double to) {
    const double lo = std::min(from, to), hi = std::max(from, to);
    const double slack = 1e-9 * (hi > lo ? hi - lo : std::abs(hi));
    return lo - slack <= x && x <= hi + slack;
}

}  // namespace NiTCAD::analysis::detail
