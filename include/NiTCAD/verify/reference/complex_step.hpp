// Complex-step differentiation in double (Unit 55, TVR-1, TC8).
//
// Sources: J. R. R. A. Martins, I. M. Kroo and J. J. Alonso, "An automated method for sensitivity analysis using
// complex variables", AIAA 2000-0689 (eqs. (11)-(14), (19)-(20), Table 5), and J. R. R. A. Martins, P. Sturdza and
// J. J. Alonso, "The connection between the complex-step derivative approximation and algorithmic
// differentiation", AIAA 2001-0921 (eqs. (1)-(9)), both read in full from aero-comlab.stanford.edu/Papers
// (2026-11-04). The cited primaries, Squire and Trapp (SIAM Rev. 40, 1998) and Martins, Sturdza and Alonso (ACM
// TOMS 29, 2003), refused automated download and were not read; grade PB.
//
// For f real on the real axis and analytic, (14): f'(x) = Im f(x + ih) / h + h^2 f'''(x) / 3! + ..., with no
// subtraction, so h can be made tiny: the papers use 1e-20 to 1e-200 and report the estimate accurate down to
// h ~ 1e-307, below which underflow somewhere in f's evaluation ruins it. The default here is 1e-200.
//
// The error estimate is Richardson's on (14): d(h) = f' - h^2 f'''/6 + O(h^4) and d(2h) = f' - 4 h^2 f'''/6 + ...,
// so the truncation of d(h) is (d(h) - d(2h)) / 3 to leading order. The next term makes that fall short of the
// true error by O(h^2) relative when f''' and f^(5) share a sign (measured 3.9e-6 at h = 1e-3), so the reported
// estimate carries the safety factor F_s = 1.25 that DECISIONS 11.5 fixes for GCI. It estimates truncation only:
// the rounding of f's own complex evaluation is not estimated (owner rulings, DECISIONS 11.13).
//
// f must be written so that its complex extension is analytic. Most standard complex functions are (Table 5,
// "mathematical standard"); the ones below are the Table's new definitions, which take their branch from real
// parts so that the imaginary part "tags along" (2000, Implementation): abs (20) with x >= 0 taking +z, so the
// right-hand derivative is returned at 0; max, min, sign and dim. Comparisons in f must likewise use real parts.
// Known pitfalls of complex formulas that mix real and imaginary parts (2001, (2)-(5)): arcsin as
// -i log(iz + sqrt(1 - z^2)) and sin as (e^{iz} - e^{-iz}) / 2i lose the h term by cancellation.
#pragma once

#include <cfloat>
#include <cmath>
#include <complex>

namespace NiTCAD::verify::reference {

using complex_d = std::complex<double>;

inline constexpr double complex_step_default_h = 1e-200;
inline constexpr double complex_step_safety_factor = 1.25;  // F_s on the Richardson estimate (DECISIONS 11.5)

enum class ComplexStepStatus {
    ok,
    invalid_input,  // x not finite, or h not a positive normal double (or 2h not finite)
    non_finite,     // f returned a non-finite value (for example after underflow inside f)
};

struct ComplexStepResult {
    double derivative = 0.0;           // Im f(x + ih) / h
    double value = 0.0;                // Re f(x + ih), f(x) to O(h^2)
    double truncation_estimate = 0.0;  // 1.25 |d(2h) - d(h)| / 3: truncation only, f's rounding not estimated
    ComplexStepStatus status = ComplexStepStatus::ok;
};

template <class F>
[[nodiscard]] ComplexStepResult complex_step(F&& f, double x, double h = complex_step_default_h) {
    ComplexStepResult out;
    if (!std::isfinite(x) || !std::isfinite(h) || !(h >= DBL_MIN) || !std::isfinite(2.0 * h)) {
        out.status = ComplexStepStatus::invalid_input;
        return out;
    }
    const complex_d f1 = f(complex_d(x, h));
    const complex_d f2 = f(complex_d(x, 2.0 * h));
    if (!std::isfinite(f1.real()) || !std::isfinite(f1.imag()) || !std::isfinite(f2.imag())) {
        out.status = ComplexStepStatus::non_finite;
        return out;
    }
    out.value = f1.real();
    out.derivative = f1.imag() / h;
    const double d2 = f2.imag() / (2.0 * h);
    out.truncation_estimate = complex_step_safety_factor * std::abs(d2 - out.derivative) / 3.0;
    return out;
}

// Martins, Kroo and Alonso 2000, Table 5: the functions with new (non-standard) complex definitions.
namespace cs {

// (20), with x >= 0 taking +z.
[[nodiscard]] inline complex_d abs(const complex_d& z) { return z.real() < 0.0 ? -z : z; }

[[nodiscard]] inline complex_d max(const complex_d& z1, const complex_d& z2) {
    return z1.real() >= z2.real() ? z1 : z2;
}

[[nodiscard]] inline complex_d min(const complex_d& z1, const complex_d& z2) {
    return z1.real() <= z2.real() ? z1 : z2;
}

// Fortran SIGN: the magnitude of z1 with the sign of z2. Table 5 writes +|x1| / -|x1|; the magnitude is taken
// with (20) so that z1's imaginary part follows it.
[[nodiscard]] inline complex_d sign(const complex_d& z1, const complex_d& z2) {
    return z2.real() >= 0.0 ? abs(z1) : -abs(z1);
}

// Fortran DIM: the positive difference.
[[nodiscard]] inline complex_d dim(const complex_d& z1, const complex_d& z2) {
    return z1.real() > z2.real() ? z1 - z2 : complex_d(0.0, 0.0);
}

}  // namespace cs

}  // namespace NiTCAD::verify::reference
