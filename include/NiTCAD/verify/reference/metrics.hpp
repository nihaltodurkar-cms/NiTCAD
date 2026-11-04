// Error metrics and numerical-uncertainty formulas (Unit 55, TVR-1, DECISIONS 4.1-4.3, 11.5; TVD-7).
//
// Source for the grid-convergence formulas: NASA NPARC Alliance, "Examining Spatial (Grid) Convergence"
// (grc.nasa.gov/www/wind/valid/tutorial/spatconv.html, updated 10 Feb 2021), read in full with its VERIFY run
// (prD.do, prD.out); PA for the formulas the page prints in text and its worked example, which the L0 test
// reproduces (DECISIONS 11.13, metrics rulings). The page's formulas, with grid 1 the finest, h_k increasing and
// r = h_{k+1} / h_k:
//   p     = ln[(f3 - f2) / (f2 - f1)] / ln r                         (three solutions, constant r)
//   f_h=0 = f1 + (f1 - f2) / (r^p - 1)                               (Richardson, order p)
//   GCI12 = F_s |(f2 - f1) / f1| / (r^p - 1),  F_s = 1.25 for >= 3 grids, 3.0 for 2 grids
// The asymptotic-range check (owner ruling): the page's GCI23 / (r^p GCI12) ~ 1 is, with p from the same three
// grids, identically f1 / f2, so it is computed and reported but not used; the check needs >= 4 grids:
// successive orders within 0.1 of each other, and (second test) rho_k = GCI_{k+2,k+3} / (r^{p_k} GCI_{k+1,k+2})
// within 10 % of 1, both GCIs formed with p_k, the order of grids k, k+1, k+2 (0-based). Algebraically
// rho_k = r^{p_{k+1} - p_k} f_{k+1} / f_{k+2}, so it is not a restatement of the three-grid identity.
//
// The factorial interaction contrast and the combination rule are DECISIONS 4.1 (TVD-17): eps_int = e(h/2, tol/100)
// - e(h/2, tol) - e(h, tol/100) + e(h, tol); U = sum |U_k| + |eps_int|, "interacting errors" when |eps_int| > 0.1 U.
//
// Field norms take vectors and weights only (DECISIONS 11.5: the mesh geometry is supplied by the consistency and
// benchmark layers). The H1 seminorm returned is sqrt(sum_e c_e (de_e)^2) with c_e = A_e / d_e; DECISIONS 4.2
// writes the sum, which is the square of the value returned. Sums are accumulated in double-double and rounded
// once, so they do not depend on summation order beyond that final rounding.
//
// Inputs and results are double: the metrics judge simulation output, which is double.
#pragma once

#include <span>
#include <string_view>
#include <vector>

namespace NiTCAD::verify::reference {

inline constexpr double gci_safety_factor_three_grids = 1.25;
inline constexpr double gci_safety_factor_two_grids = 3.0;
inline constexpr double asymptotic_order_spread = 0.1;  // successive orders within 0.1 (DECISIONS 4.2)
inline constexpr double asymptotic_ratio_band = 0.1;    // GCI ratio within 10 % (DECISIONS 4.2)
inline constexpr double constant_ratio_tolerance = 1e-12;  // relative; constant r (owner ruling)
inline constexpr double interaction_fraction = 0.1;      // |eps_int| > 10 % of U: interacting (DECISIONS 4.1)

enum class MetricStatus {
    ok,
    invalid_input,  // a malformed argument: non-finite, negative floor or weight, size mismatch, r <= 1, ...
    not_evaluated,  // a well-formed argument on which the metric is not defined; the reason says why
};

struct Metric {
    double value = 0.0;
    MetricStatus status = MetricStatus::ok;
    std::string_view reason;  // empty when ok
};

// |q_sim - q_ref| / max(|q_ref|, q_floor). q_floor >= 0 (the current resolution for currents).
[[nodiscard]] Metric relative_error(double q_sim, double q_ref, double q_floor);

// |log10(q_sim / q_ref)|, defined only when q_sim > floor_sim and q_ref > floor_ref (floors >= 0; magnitudes are
// the caller's to take, so a negative current is passed as |I|).
[[nodiscard]] Metric log_error(double q_sim, double q_ref, double floor_sim, double floor_ref);

// max_i |e_i|.
[[nodiscard]] Metric norm_linf(std::span<const double> e);
// sqrt(sum_i w_i e_i^2), w_i >= 0 (dual volumes).
[[nodiscard]] Metric norm_l2(std::span<const double> e, std::span<const double> w);
// sqrt(sum_e c_e de_e^2), c_e >= 0 (A_e / d_e of the dual faces), de_e the difference of e along edge e.
[[nodiscard]] Metric seminorm_h1(std::span<const double> de, std::span<const double> c);

// Observed orders from errors against an exact solution: p_k = ln(e_{k+1} / e_k) / ln(h_{k+1} / h_k), k = 0 ...
// N - 2, with h strictly increasing (finest first) and e > 0. Any ratio is allowed: the exact solution removes the
// constant-r condition of the three-solution form.
[[nodiscard]] std::vector<Metric> observed_orders_from_errors(std::span<const double> h, std::span<const double> e);

// Order from three solutions at constant ratio r > 1, f1 the finest. not_evaluated when f2 - f1 or f3 - f2 is
// zero, when they differ in sign (not monotone), or when p <= 0.
[[nodiscard]] Metric observed_order_from_solutions(double f1, double f2, double f3, double r);

// f1 + (f1 - f2) / (r^p - 1); r > 1, p > 0.
[[nodiscard]] Metric richardson(double f1, double f2, double r, double p);

// F_s |(f_coarse - f_fine) / f_fine| / (r^p - 1) as a fraction (not per cent), F_s = 1.25 when grids >= 3 and
// 3.0 when grids == 2 (then p is the caller's theoretical order). not_evaluated when f_fine = 0 (NPARC: the
// relative error is meaningless there).
[[nodiscard]] Metric gci(double f_fine, double f_coarse, double r, double p, int grids);

struct GridStudy {
    MetricStatus status = MetricStatus::ok;
    std::string_view reason;
    double r = 0.0;
    std::vector<double> orders;  // p_k from grids k, k+1, k+2 (N - 2 values)
    double extrapolated = 0.0;   // Richardson from grids 0, 1 with p_0
    std::vector<double> gci;     // GCI_{k,k+1} with p_0, F_s = 1.25 (N - 1 values), fractions
    double three_grid_ratio = 0.0;  // GCI_{1,2} / (r^{p_0} GCI_{0,1}): identically f_0 / f_1, reported only
    // Asymptotic-range check: ok with >= 4 grids, else not_evaluated ("needs 4 grids").
    MetricStatus asymptotic_status = MetricStatus::not_evaluated;
    std::string_view asymptotic_reason;
    bool in_asymptotic_range = false;
    double order_spread = 0.0;     // max_k |p_{k+1} - p_k|
    double ratio_deviation = 0.0;  // max_k |rho_k - 1|
};

// A grid study over N >= 3 solutions f (finest first) on spacings h (strictly increasing, constant ratio).
// not_evaluated (with every field after r left empty or zero) when r is not constant, the solutions are not
// monotone (every successive difference nonzero and of one sign), any p_k <= 0, or any f_k = 0 (each GCI is
// relative to its finer solution).
[[nodiscard]] GridStudy grid_convergence(std::span<const double> h, std::span<const double> f);

// e(h/2, tol/100) - e(h/2, tol) - e(h, tol/100) + e(h, tol), accumulated in double-double.
[[nodiscard]] Metric interaction_contrast(double e_h_tol, double e_h_tolfine, double e_hfine_tol,
                                          double e_hfine_tolfine);

struct CombinedUncertainty {
    MetricStatus status = MetricStatus::ok;
    std::string_view reason;
    double value = 0.0;        // sum |U_k| + |eps_int|
    bool interacting = false;  // |eps_int| > 0.1 value
};

// U = sum |U_k| + |eps_int| (DECISIONS 4.1(c), a sum of magnitudes, not a root-sum-square).
[[nodiscard]] CombinedUncertainty combine_uncertainty(std::span<const double> components, double eps_int);

}  // namespace NiTCAD::verify::reference
