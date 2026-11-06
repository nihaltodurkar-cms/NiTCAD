// Double-double special functions for the model references of TVB-31 (Unit 55, TVR-3 step 2;
// DECISIONS 15.11). They depend on double_double only (TVD-19).
//
// Sources (NIST DLMF, read 2026-11-06; provenance PA):
//   - log1p: 2 atanh(x / (2 + x)) by its series for |x| < 1/8 (DLMF 4.24 / 4.37 series), log(1 + x) otherwise;
//   - erf: erf x = (2 / sqrt(pi)) e^(-x^2) sum_n 2^n x^(2n+1) / (1 3 ... (2n+1)) (DLMF 7.6.2, positive terms)
//     for |x| < 3; erfc by the even-part continued fraction sqrt(pi) e^(x^2) erfc x =
//     2x / (2x^2 + 1 - 1*2 / (2x^2 + 5 - 3*4 / (2x^2 + 9 - ...))) (DLMF 7.9.2) for x >= 3;
//   - sin, cos: Taylor series, for |x| <= 2 only (atan needs |x| <= pi/4, its tests pi/2);
//   - atan: Newton on sin y - x cos y = 0 from std::atan, after reducing |x| > 1 by atan x = pi/2 - atan(1/x);
//   - asin x = atan2(x, sqrt((1 - x)(1 + x)));
//   - li2 (the dilogarithm, real x <= 1): the series (DLMF 25.12.1) on [-1/2, 1/2], Landen's identity
//     Li2(x) = -Li2(x/(x-1)) - ln^2(1-x)/2 (25.12.3) on [-1, -1/2), inversion
//     Li2(x) = -pi^2/6 - ln^2(-x)/2 - Li2(1/x) (25.12.4) below -1, reflection
//     Li2(x) = pi^2/6 - ln x ln(1-x) - Li2(1-x) (25.12.6) on (1/2, 1).
//
// Domains as double_double's: values of magnitude in [1e-280, 1e280]. Accuracy: a few u_dd relative, except
// where the function's own condition number amplifies the input's (the L0 tests' bounds, DECISIONS 15.11).
#pragma once

#include "NiTCAD/verify/reference/double_double.hpp"

namespace NiTCAD::verify::reference {

// ln(1 + x), x > -1; relative-accurate near 0.
[[nodiscard]] DD log1p(const DD& x);

// ln(1 + e^x) without overflow: x + log1p(e^-x) for x > 0. Full accuracy while e^-|x| stays normal in
// double-double (|x| <= 650); below about -745 it returns 0 (the value is below the smallest double).
[[nodiscard]] DD log1p_exp(const DD& x);

// tanh x = expm1(2x) / (expm1(2x) + 2) with the sign of x; +-1 for |x| > 40 (the remainder 2 e^-80 is below
// u_dd).
[[nodiscard]] DD tanh(const DD& x);

// erf and erfc for real x. erfc returns 0 for x > 26.5 (e^(-x^2) below the smallest double). Below x = 3 erfc is
// 1 - erf: accurate to a few u_dd absolutely, so its relative error grows by erf / erfc (4.5e4 at x = 3).
[[nodiscard]] DD erf(const DD& x);
[[nodiscard]] DD erfc(const DD& x);

// sin and cos for |x| <= 2 (precondition; the series is used as is, so cos near pi/2 is accurate absolutely).
[[nodiscard]] DD sin(const DD& x);
[[nodiscard]] DD cos(const DD& x);

// atan for any finite x; atan2(y, x) for (y, x) not both 0, in (-pi, pi]; asin for |x| <= 1.
[[nodiscard]] DD atan(const DD& x);
[[nodiscard]] DD atan2(const DD& y, const DD& x);
[[nodiscard]] DD asin(const DD& x);

// The dilogarithm Li2(x) for real x <= 1.
[[nodiscard]] DD li2(const DD& x);

}  // namespace NiTCAD::verify::reference
