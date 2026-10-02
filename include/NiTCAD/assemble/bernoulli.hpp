// The Bernoulli function of the Scharfetter-Gummel flux, B(x) = x / (e^x - 1), B(0) = 1, and its
// derivative B'(x) (legacy kernels.hpp bernoulli / dbernoulli).
//
// Identities: B(-x) = B(x) + x and B'(x) + B'(-x) = -1. Limits: B(x) -> -x as x -> -inf (with
// B' -> -1), B(x) -> x e^-x -> 0 as x -> +inf.
//
// Evaluation, accurate to a few ulp for B and about 1e-13 relative for B' (tested against 50-digit
// values):
// - |x| < 1e-2: Taylor series to x^8 (B) and x^7 (B'); the next terms are below 1e-20.
// - otherwise B = x / expm1(x) and B' = (1 - x - B) / expm1(x), from B'(x) = B(x)(1 - B(-x)) / x,
//   the identity above and B / x = 1 / expm1(x); the subtraction loses at most about 2 digits near
//   |x| = 1e-2.
// No clipping: expm1 gives B(x) = -x for very negative x and B(x) = 0 once x e^-x underflows.
// OLD: the legacy clipped x to [-700, 700], so B(-1000) was 700 instead of 1000, and used the
// series only below 1e-4, where B' lost about 4 digits to cancellation.
#pragma once

#include <cmath>

namespace NiTCAD::assemble {

[[nodiscard]] inline double bernoulli(double x) noexcept {
    if (std::abs(x) < 1e-2) {
        const double x2 = x * x;
        // 1 - x/2 + x^2/12 - x^4/720 + x^6/30240 - x^8/1209600
        return 1.0 - 0.5 * x +
               x2 * (1.0 / 12.0 + x2 * (-1.0 / 720.0 + x2 * (1.0 / 30240.0 - x2 / 1209600.0)));
    }
    return x / std::expm1(x);
}

[[nodiscard]] inline double bernoulli_derivative(double x) noexcept {
    if (std::abs(x) < 1e-2) {
        const double x2 = x * x;
        // -1/2 + x/6 - x^3/180 + x^5/5040 - x^7/151200
        return -0.5 +
               x * (1.0 / 6.0 + x2 * (-1.0 / 180.0 + x2 * (1.0 / 5040.0 - x2 / 151200.0)));
    }
    const double em1 = std::expm1(x);
    return (1.0 - x - x / em1) / em1;
}

}  // namespace NiTCAD::assemble
