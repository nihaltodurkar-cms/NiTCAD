// The Bernoulli function B(x) = x / (e^x - 1), B(0) = 1, in double-double, for TVB-03 (Unit 55, TVR-1; DECISIONS
// 11.8, 11.13 TB-Q9).
//
// Two methods, written from the definition and the generating function x / (e^x - 1) = sum_n B_n x^n / n!:
//   - direct: x / expm1(x) with the double-double expm1 of TC1, for |x| <= 670; above 670 the scaled form
//     B = x e^-x / (1 - e^-x) = (x e^-r) 2^-m, m = floor(x / ln 2), r = x - m ln 2 (the double-double ln 2 of TC1;
//     its 2^-107 relative error times m <= 1075 moves r by < 2^-96, so e^-r keeps ~29 digits), with 1 - e^-x = 1 to
//     far below u_dd (e^-670 < 1e-290); below -670, B = -x / (1 - e^x) = -x to the same order;
//   - series, |x| <= 1/2: 1 - x/2 + sum_{k=1}^{17} B_2k x^2k / (2k)!, the B_2k typed from OEIS A000367 (numerators)
//     and A002445 (denominators), b-files read 2026-11-05 (provenance PA). At |x| = 1/2 the first omitted term,
//     |B_36| / 36! (1/2)^36 ~ 2 (4 pi)^-36, is below 1e-39.
// The value is returned as m 2^e (a double-double mantissa and a binary exponent) so that B is held to double-double
// accuracy where it is far below the smallest double, and rounded only when compared.
// Depends on the standard library and double_double.hpp only (TVD-19).
#pragma once

#include <array>
#include <cmath>
#include <cstdint>

#include "NiTCAD/verify/reference/double_double.hpp"

namespace NiTCAD::verify::reference {

struct ScaledDD {
    DD mantissa;
    int exponent = 0;  // value = mantissa * 2^exponent
};

// B_2k = numerator / denominator for k = 1 .. 17 (OEIS A000367 and A002445, n = k).
inline constexpr std::array<std::int64_t, 17> bernoulli_numerators{
    1, -1, 1, -1, 5, -691, 7, -3617, 43867, -174611, 854513, -236364091, 8553103, -23749461029, 8615841276005,
    -7709321041217, 2577687858367};
inline constexpr std::array<std::int64_t, 17> bernoulli_denominators{
    6, 30, 42, 30, 66, 2730, 6, 510, 798, 330, 138, 2730, 6, 870, 14322, 510, 6};

// c_k = B_2k / (2k)!, k = 1 .. 17 (index k - 1); each to a few u_dd (the numerators and denominators are exact
// doubles; (2k)! is exact as a double-double up to 30! and rounded once per factor above).
[[nodiscard]] inline const std::array<DD, 17>& bernoulli_series_coefficients() {
    static const std::array<DD, 17> c = [] {
        std::array<DD, 17> out{};
        DD factorial(1.0);
        for (int k = 1; k <= 17; ++k) {
            factorial = factorial * DD(static_cast<double>(2 * k - 1)) * DD(static_cast<double>(2 * k));
            out[static_cast<std::size_t>(k - 1)] =
                DD(static_cast<double>(bernoulli_numerators[static_cast<std::size_t>(k - 1)])) /
                (DD(static_cast<double>(bernoulli_denominators[static_cast<std::size_t>(k - 1)])) * factorial);
        }
        return out;
    }();
    return c;
}

// The series, |x| <= 1/2.
[[nodiscard]] inline DD bernoulli_series(const DD& x) {
    const auto& c = bernoulli_series_coefficients();
    const DD x2 = x * x;
    DD s(0.0);
    for (std::size_t k = c.size(); k-- > 0;) s = (s + c[k]) * x2;
    return DD(1.0) - ldexp(x, -1) + s;
}

// Its derivative, |x| <= 1/2 (TVB-03's diagnostic of B').
[[nodiscard]] inline DD bernoulli_series_derivative(const DD& x) {
    const auto& c = bernoulli_series_coefficients();
    const DD x2 = x * x;
    DD s(0.0);
    for (std::size_t k = c.size(); k-- > 0;) s = s * x2 + c[k] * DD(static_cast<double>(2 * (k + 1)));
    return DD(-0.5) + s * x;
}

// The direct method, every finite x.
[[nodiscard]] inline ScaledDD bernoulli_direct(double x) {
    if (x == 0.0) return {DD(1.0), 0};
    if (x < -670.0) return {DD(-x), 0};
    if (x <= 670.0) return {DD(x) / expm1(DD(x)), 0};
    const int m = static_cast<int>(std::floor(x / ln2().hi));
    const DD r = DD(x) - DD(static_cast<double>(m)) * ln2();
    return {DD(x) * exp(-r), -m};
}

// Whether the value is in the normal double range (TVB-03 judges only there; below it is a diagnostic).
[[nodiscard]] inline bool is_normal_double(const ScaledDD& b) {
    if (b.mantissa.hi == 0.0) return false;
    return std::ilogb(b.mantissa.hi) + b.exponent >= -1022 && std::ilogb(b.mantissa.hi) + b.exponent <= 1023;
}

// |sim - B| / |B|, sim scaled exactly by 2^-exponent (sim normal or 0; B nonzero).
[[nodiscard]] inline double relative_error(double sim, const ScaledDD& b) {
    const DD s(std::ldexp(sim, -b.exponent));
    return std::abs((s - b.mantissa).to_double()) / std::abs(b.mantissa.to_double());
}

}  // namespace NiTCAD::verify::reference
