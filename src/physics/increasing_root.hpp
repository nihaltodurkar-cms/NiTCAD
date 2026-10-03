// Private to the physics layer: the root of an increasing function, shared by the inverse Fermi
// integral and the Fermi-Dirac charge-neutral equilibrium.
#pragma once

#include <algorithm>
#include <cmath>
#include <limits>

namespace NiTCAD::physics::detail {

struct ValueSlope {
    double value;
    double slope;
};

// The root of an increasing function, f(x) returning its value and slope, from the guess x0:
// Newton's method, kept inside the bracket the iterates have established (bisection when a step
// leaves it, outward doubling while one side is still open). It returns one Newton step after a
// step of at most 1e-10 max(1, |x|), which quadratic convergence leaves at rounding, or when the
// bracket has closed to rounding; the caller guarantees a root exists.
template <class F>
[[nodiscard]] double increasing_root(F f, double x0) {
    constexpr double infinity = std::numeric_limits<double>::infinity();
    constexpr double eps = std::numeric_limits<double>::epsilon();
    double lo = -infinity, hi = infinity;  // f(lo) < 0 < f(hi)
    double x = x0;
    for (int iteration = 0; iteration < 400; ++iteration) {
        const ValueSlope v = f(x);
        if (v.value == 0.0) return x;
        if (v.value < 0.0) {
            lo = x;
        } else {
            hi = x;
        }
        const double width = std::max(1.0, std::abs(x));
        double next = x - v.value / v.slope;
        const bool newton = next > lo && next < hi;
        if (newton && std::abs(next - x) <= 1e-10 * width) return next;
        if (!newton) {
            if (lo == -infinity) {
                next = hi - std::max(1.0, std::abs(hi));
            } else if (hi == infinity) {
                next = lo + std::max(1.0, std::abs(lo));
            } else {
                next = 0.5 * (lo + hi);
            }
        }
        if (hi - lo <= 4.0 * eps * width) return next;
        x = next;
    }
    return x;
}

}  // namespace NiTCAD::physics::detail
