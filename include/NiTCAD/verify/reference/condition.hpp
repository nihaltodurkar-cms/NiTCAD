// The condition sum S(q) = sum_i |x_i dq/dx_i| of DECISIONS Unit 55 15.9 (TVR-3), the scale of the M criterion
// |q_sim - q_ref| <= max(1e-13 |q_ref|, 64 eps S(q)).
//
// Each dq/dx_i is a double-double central difference with the step h_i = 2^-27 |x_i| (exact): its truncation error is
// O(h^2) = 2^-54 relative and its rounding u_dd / 2^-27 = 2^-77, so S is good to about 1e-16 relative, which is
// ample for a tolerance scale (15.5: "a tolerance scale only, not a judged quantity"). Inputs equal to 0 add
// nothing (|x_i| = 0). A difference that straddles a branch (a cut, a switch point) still yields a finite scale.
// Base-free (TVD-19).
#pragma once

#include <array>
#include <cmath>
#include <cstddef>

#include "NiTCAD/verify/reference/double_double.hpp"

namespace NiTCAD::verify::reference {

inline constexpr int condition_step_log2 = -27;

// f: std::array<DD, N> -> DD.
template <std::size_t N, class F>
[[nodiscard]] double condition_sum(F&& f, const std::array<DD, N>& x) {
    double s = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        if (x[i].hi == 0.0) continue;
        const DD h = ldexp(abs(x[i]), condition_step_log2);
        std::array<DD, N> up = x;
        std::array<DD, N> down = x;
        up[i] = x[i] + h;
        down[i] = x[i] - h;
        const DD d = (f(up) - f(down)) / ldexp(h, 1);
        s += std::abs((x[i] * d).hi);
    }
    return s;
}

// The same for a vector-valued f: std::array<DD, N> -> std::array<DD, K>, one evaluation per perturbation.
template <std::size_t K, std::size_t N, class F>
[[nodiscard]] std::array<double, K> condition_sums(F&& f, const std::array<DD, N>& x) {
    std::array<double, K> s{};
    for (std::size_t i = 0; i < N; ++i) {
        if (x[i].hi == 0.0) continue;
        const DD h = ldexp(abs(x[i]), condition_step_log2);
        std::array<DD, N> up = x;
        std::array<DD, N> down = x;
        up[i] = x[i] + h;
        down[i] = x[i] - h;
        const std::array<DD, K> fu = f(up);
        const std::array<DD, K> fd = f(down);
        for (std::size_t k = 0; k < K; ++k) s[k] += std::abs((x[i] * ((fu[k] - fd[k]) / ldexp(h, 1))).hi);
    }
    return s;
}

}  // namespace NiTCAD::verify::reference
