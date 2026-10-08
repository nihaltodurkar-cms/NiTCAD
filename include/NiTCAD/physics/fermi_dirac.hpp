// The complete Fermi-Dirac integral of order 1/2 (ARCHITECTURE.md section 5, Unit 14; legacy
// pytcad/fermi.py and core/include/tcad/physics/fermi.hpp):
//     F_{1/2}(eta) = (2 / sqrt(pi)) Integral_0^inf t^(1/2) / (1 + e^(t - eta)) dt,
// normalized so that F_{1/2}(eta) -> e^eta as eta -> -inf, with dF_{1/2}/deta = F_{-1/2}.
//
// Everything is expressed through the degeneracy factor gamma(eta) = F_{1/2}(eta) / e^eta, which is
// 1 in the Boltzmann limit and falls below it as the carriers degenerate (gamma(0) = 0.765,
// gamma(20) = 1.4e-7). ln gamma is evaluated, never F_{1/2} itself, so nothing overflows or
// underflows where the result is representable:
// - eta <= -2: the exact series F_{1/2} = sum_k (-1)^(k+1) e^(k eta) / k^(3/2), summed for ln gamma
//   and its derivative directly (no cancellation; at most 22 terms);
// - -2 < eta <= 40: quintic Hermite interpolation of ln gamma on a table of step 0.025, built once
//   per process from the values and first two derivatives of ln gamma, which come from
//   Gauss-Legendre quadratures of F_{1/2}, F_{-1/2} and F_{-3/2};
// - eta > 40: the Sommerfeld expansion F_{1/2} = 4 / (3 sqrt(pi)) eta^(3/2) (1 + sum_k c_k eta^-2k)
//   with ten terms, whose residual falls like e^-eta (3.7e-18 at eta = 40).
// Measured against double-double references (tests/physics/fermi_dirac_test.cpp): F_{1/2} within
// 1.5e-15 relative everywhere; F_{-1/2} within 1.3e-13 up to eta = 15 and 3e-12 up to 40, where
// the slope F_{-1/2} / F_{1/2} has fallen to 0.04, and 1e-14 above. The derivatives returned are
// those of the function evaluated (the interpolant's own derivative inside the table), so a
// Jacobian built from them matches finite differences of the residual. ln gamma is continuously
// differentiable; at the joins its value jumps by rounding (1e-17 at -2, an ulp of 40 at 40) and
// its slope by about 1e-16.
//
// OLD / NEW / REASON (ARCHITECTURE.md 5, Unit 14):
// - OLD: tabulated F_{1/2} and F_{-1/2} on [-10, 40] (step 0.005, cubic Hermite on ln F, with
//   F_{-1/2} from its own table, so it was not exactly the derivative of F_{1/2}); outside
//   [-40, 40] the functions refused. NEW: one table of ln gamma (quintic, step 0.025), the series
//   below and Sommerfeld above, so every eta is valid. REASON: the derivative used by Newton is
//   then exactly that of the residual; and a refusal at eta = 40 (n near 5e21 cm^-3) gave a solver
//   a hard wall where a parabolic band is merely an approximation.
#pragma once

namespace NiTCAD::physics {

struct FermiIntegral {
    double value;       // F_{1/2}(eta)
    double derivative;  // dF_{1/2}/deta = F_{-1/2}(eta)
};

// F_{1/2}(eta) and F_{-1/2}(eta). Underflows to 0 below eta of about -745 (as e^eta does) and
// overflows only where eta^(3/2) does. A NaN eta gives NaN.
[[nodiscard]] FermiIntegral fermi_half(double eta) noexcept;

struct LogDegeneracy {
    double value;  // ln gamma(eta) = ln F_{1/2}(eta) - eta, <= 0
    double d_eta;  // d ln gamma / d eta = F_{-1/2} / F_{1/2} - 1, in (-1, 0]
};

// ln gamma(eta) and its derivative, finite for every finite eta. A NaN eta gives NaN.
[[nodiscard]] LogDegeneracy log_degeneracy(double eta) noexcept;

// The eta with F_{1/2}(eta) = nu (of the function above, to rounding). Safeguarded Newton on
// ln F_{1/2}(eta) = ln nu, typically 2-5 evaluations. Precondition (NITCAD_EXPECTS): nu is finite
// and positive.
[[nodiscard]] double inverse_fermi_half(double nu);

}  // namespace NiTCAD::physics
