// A dot product as accurate as if computed in twice the working precision (Unit 36): Dot2 of
// T. Ogita, S. M. Rump, S. Oishi, "Accurate sum and dot product", SIAM J. Sci. Comput. 26, 1955
// (2005), with each product's error from an FMA (TwoProduct) and each sum's from TwoSum. Its error
// is about eps |a . b| + n^2 eps^2 sum |a_q b_q|: where the terms cancel over many orders (the
// arc row of solve::trace_bias, 13 measured) the plain sum keeps only rounding. It needs value-safe
// floating point (the build's /fp:precise: no reassociation, no contraction of a b + c).
#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <span>

namespace NiTCAD::linalg {

// sum over q of a[q] b[q], leaving out q == skip. Precondition: a.size() == b.size().
[[nodiscard]] inline double compensated_dot(
    std::span<const double> a, std::span<const double> b,
    std::size_t skip = std::numeric_limits<std::size_t>::max()) noexcept {
    double sum = 0.0, error = 0.0;
    for (std::size_t q = 0; q < a.size(); ++q) {
        if (q == skip) continue;
        const double p = a[q] * b[q];
        const double p_error = std::fma(a[q], b[q], -p);  // a b - p, exactly
        const double t = sum + p;
        const double z = t - sum;
        error += p_error + ((sum - (t - z)) + (p - z));  // + (sum + p) - t, exactly
        sum = t;
    }
    return sum + error;
}

}  // namespace NiTCAD::linalg
