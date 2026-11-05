// Complete Fermi-Dirac integrals F_j, j = -1/2, 1/2, 3/2, in double-double, for TVB-02 (Unit 55, TVR-1; DECISIONS
// 11.8, 11.13 TB-Q4, TB-Q5, TB-Q7).
//
// Definition (provenance PA): NIST DLMF section 25.12(iii), eq. 25.12.14,
//     F_j(x) = 1 / Gamma(j + 1) Integral_0^inf t^j / (e^(t - x) + 1) dt,
// normalised so that F_j(x) -> e^x as x -> -inf and dF_j/dx = F_{j-1}; and eq. 25.12.16, F_j(x) = -Li_{j+1}(-e^x),
// whose power series (25.12.10, Li_s(z) = sum z^n / n^s, |z| < 1) gives
// F_j(x) = sum_{k>=1} (-1)^(k+1) e^(kx) / k^(j+1) for x < 0.
//
// Quadrature form (derived here): t = u^2 gives
//     F_j(x) = 2 / Gamma(j + 1) Integral_0^inf u^(2j+1) / (1 + e^(u^2 - x)) du,   2j + 1 = 0, 2 or 4,
// which is smooth (the t^-1/2 of F_-1/2 is gone). The range is cut at U = sqrt(x + L), L = 100, and split at
// sqrt(x) when x > 0 (the Fermi edge on a piece's endpoint, where both rules cluster their nodes). The tail is
// bounded, not estimated: on u >= U the integrand is below 2 u^p e^(x - u^2), p = 2j + 1. With I(p) =
// Integral_U^inf u^p e^(-u^2) du, integration by parts gives I(p) = U^(p-1) e^(-U^2) / 2 + (p - 1) / 2 I(p - 2). For
// p >= 1, I(p - 2) <= I(p) / U^2 gives I(p) <= U^(p-1) e^(-U^2) / (2 - (p - 1) / U^2) (needs U^2 > (p - 1) / 2; here
// U^2 >= 60); for p = 0 the coefficient is negative and I(0) <= e^(-U^2) / (2U) (an error in the first draft, caught
// by the test that integrates the omitted tail). The tail of F_j is at most 2 / Gamma(j + 1) e^x I(p), e^(x - U^2) =
// e^-L. The logistic factor is evaluated as e^-y / (1 + e^-y) for y = u^2 - x > 0 and 1 / (1 + e^y) otherwise.
// Domain: -100 < x <= 600 (the tail's e^(x - u^2) and the cut need x + L > 0; TVB-02 uses [-40, 40]).
#pragma once

#include "NiTCAD/verify/reference/double_double.hpp"

namespace NiTCAD::verify::reference {

enum class FermiOrder { minus_half, half, three_halves };

struct FermiQuadrature {
    DD value;
    double abs_error = 0.0;   // the rule's own estimate, summed over the pieces, times the normalisation
    double tail_bound = 0.0;  // the cut's bound, normalised
    bool ok = false;          // every piece's rule reported success
    int evaluations = 0;
};

// The tanh-sinh (TC4) and the adaptive Gauss-Kronrod (TC3) quadratures of the form above, epsrel 1e-27 per piece.
[[nodiscard]] FermiQuadrature fermi_tanh_sinh(FermiOrder j, const DD& x);
[[nodiscard]] FermiQuadrature fermi_gauss_kronrod(FermiOrder j, const DD& x);

// The polylogarithm series for x <= -2 (terms until e^(kx) < 2^-120), with its truncation bound.
struct FermiSeries {
    DD value;
    double truncation_bound = 0.0;
    int terms = 0;
};
[[nodiscard]] FermiSeries fermi_series(FermiOrder j, const DD& x);

// 1 / Gamma(j + 1): 1 / sqrt(pi), 2 / sqrt(pi), 4 / (3 sqrt(pi)).
[[nodiscard]] DD fermi_normalisation(FermiOrder j);

// (1 - 2^-j), the factor of F_j(0) = (1 - 2^-j) zeta(j + 1) (the Dirichlet eta function at j + 1).
[[nodiscard]] DD fermi_eta_factor(FermiOrder j);

// The normalised integrand 2 / Gamma(j + 1) u^(2j+1) / (1 + e^(u^2 - x)) and the cut U = sqrt(x + 100), so that a
// test can integrate the omitted tail [U, inf) and check the bound against it.
[[nodiscard]] DD fermi_integrand(FermiOrder j, const DD& x, const DD& u);
[[nodiscard]] DD fermi_cut(const DD& x);

}  // namespace NiTCAD::verify::reference
