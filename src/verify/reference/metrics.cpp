// Error metrics and numerical-uncertainty formulas (Unit 55, TVR-1, DECISIONS 11.5). See metrics.hpp.
#include "NiTCAD/verify/reference/metrics.hpp"

#include "NiTCAD/verify/reference/double_double.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace NiTCAD::verify::reference {

namespace {

Metric invalid(std::string_view why) { return {0.0, MetricStatus::invalid_input, why}; }
Metric not_evaluated(std::string_view why) { return {0.0, MetricStatus::not_evaluated, why}; }
Metric ok(double v) { return {v, MetricStatus::ok, {}}; }

bool finite(double x) { return std::isfinite(x); }

bool all_finite(std::span<const double> v) {
    return std::all_of(v.begin(), v.end(), [](double x) { return std::isfinite(x); });
}

// sqrt(sum w_i x_i^2) accumulated in double-double. x is scaled by a power of two (exact) so that the squares
// neither overflow nor underflow, and the weights likewise.
double weighted_root_sum_squares(std::span<const double> x, std::span<const double> w) {
    double xmax = 0.0, wmax = 0.0;
    for (double v : x) xmax = std::max(xmax, std::abs(v));
    for (double v : w) wmax = std::max(wmax, std::abs(v));
    if (xmax == 0.0 || wmax == 0.0) return 0.0;
    int ex = 0, ew = 0;
    (void)std::frexp(xmax, &ex);
    (void)std::frexp(wmax, &ew);
    if (ew % 2 != 0) ++ew;  // an even exponent so that its square root is a power of two
    DD sum(0.0);
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double xs = std::ldexp(x[i], -ex);
        sum += two_prod(xs, xs) * DD(std::ldexp(w[i], -ew));
    }
    return std::ldexp(sqrt(sum).to_double(), ex + ew / 2);
}

}  // namespace

Metric relative_error(double q_sim, double q_ref, double q_floor) {
    if (!finite(q_sim) || !finite(q_ref) || !finite(q_floor)) return invalid("non-finite argument");
    if (q_floor < 0.0) return invalid("negative floor");
    const double scale = std::max(std::abs(q_ref), q_floor);
    if (scale == 0.0) return not_evaluated("reference and floor both zero");
    // The difference of two doubles is exact in double-double, so only the final rounding and division remain.
    return ok(std::abs((DD(q_sim) - DD(q_ref)).to_double()) / scale);
}

Metric log_error(double q_sim, double q_ref, double floor_sim, double floor_ref) {
    if (!finite(q_sim) || !finite(q_ref) || !finite(floor_sim) || !finite(floor_ref)) {
        return invalid("non-finite argument");
    }
    if (floor_sim < 0.0 || floor_ref < 0.0) return invalid("negative floor");
    if (!(q_sim > floor_sim) || !(q_ref > floor_ref)) return not_evaluated("at or below floor");
    return ok(std::abs(std::log10(q_sim / q_ref)));
}

Metric norm_linf(std::span<const double> e) {
    if (e.empty()) return invalid("empty vector");
    if (!all_finite(e)) return invalid("non-finite entry");
    double m = 0.0;
    for (double x : e) m = std::max(m, std::abs(x));
    return ok(m);
}

Metric norm_l2(std::span<const double> e, std::span<const double> w) {
    if (e.empty()) return invalid("empty vector");
    if (e.size() != w.size()) return invalid("size mismatch");
    if (!all_finite(e) || !all_finite(w)) return invalid("non-finite entry");
    if (std::any_of(w.begin(), w.end(), [](double x) { return x < 0.0; })) return invalid("negative weight");
    return ok(weighted_root_sum_squares(e, w));
}

Metric seminorm_h1(std::span<const double> de, std::span<const double> c) {
    if (de.empty()) return invalid("empty vector");
    if (de.size() != c.size()) return invalid("size mismatch");
    if (!all_finite(de) || !all_finite(c)) return invalid("non-finite entry");
    if (std::any_of(c.begin(), c.end(), [](double x) { return x < 0.0; })) return invalid("negative coupling");
    return ok(weighted_root_sum_squares(de, c));
}

std::vector<Metric> observed_orders_from_errors(std::span<const double> h, std::span<const double> e) {
    if (h.size() != e.size() || h.size() < 2) return {invalid("need >= 2 matching spacings and errors")};
    if (!all_finite(h) || !all_finite(e)) return {invalid("non-finite entry")};
    for (std::size_t k = 0; k < h.size(); ++k) {
        if (!(h[k] > 0.0)) return {invalid("spacing not positive")};
        if (k > 0 && !(h[k] > h[k - 1])) return {invalid("spacings not strictly increasing")};
    }
    std::vector<Metric> out;
    out.reserve(h.size() - 1);
    for (std::size_t k = 0; k + 1 < h.size(); ++k) {
        if (!(e[k] > 0.0) || !(e[k + 1] > 0.0)) {
            out.push_back(not_evaluated("error not positive"));
            continue;
        }
        out.push_back(ok(std::log(e[k + 1] / e[k]) / std::log(h[k + 1] / h[k])));
    }
    return out;
}

Metric observed_order_from_solutions(double f1, double f2, double f3, double r) {
    if (!finite(f1) || !finite(f2) || !finite(f3) || !finite(r)) return invalid("non-finite argument");
    if (!(r > 1.0)) return invalid("ratio not > 1");
    const double d12 = f2 - f1, d23 = f3 - f2;
    if (d12 == 0.0 || d23 == 0.0) return not_evaluated("zero difference");
    if ((d12 > 0.0) != (d23 > 0.0)) return not_evaluated("not monotone");
    const double p = std::log(d23 / d12) / std::log(r);
    if (!(p > 0.0)) return not_evaluated("observed order <= 0");
    return ok(p);
}

Metric richardson(double f1, double f2, double r, double p) {
    if (!finite(f1) || !finite(f2) || !finite(r) || !finite(p)) return invalid("non-finite argument");
    if (!(r > 1.0)) return invalid("ratio not > 1");
    if (!(p > 0.0)) return invalid("order not > 0");
    return ok(f1 + (f1 - f2) / (std::pow(r, p) - 1.0));
}

Metric gci(double f_fine, double f_coarse, double r, double p, int grids) {
    if (!finite(f_fine) || !finite(f_coarse) || !finite(r) || !finite(p)) return invalid("non-finite argument");
    if (!(r > 1.0)) return invalid("ratio not > 1");
    if (!(p > 0.0)) return invalid("order not > 0");
    if (grids < 2) return invalid("fewer than 2 grids");
    if (f_fine == 0.0) return not_evaluated("fine solution zero");
    const double fs = grids >= 3 ? gci_safety_factor_three_grids : gci_safety_factor_two_grids;
    return ok(fs * std::abs((f_coarse - f_fine) / f_fine) / (std::pow(r, p) - 1.0));
}

GridStudy grid_convergence(std::span<const double> h, std::span<const double> f) {
    GridStudy s;
    auto refuse = [&s](MetricStatus st, std::string_view why) {
        s.status = st;
        s.reason = why;
        s.asymptotic_status = MetricStatus::not_evaluated;
        s.asymptotic_reason = why;
        return s;
    };
    const std::size_t n = f.size();
    if (h.size() != n || n < 3) return refuse(MetricStatus::invalid_input, "need >= 3 matching spacings and solutions");
    if (!all_finite(h) || !all_finite(f)) return refuse(MetricStatus::invalid_input, "non-finite entry");
    for (std::size_t k = 0; k < n; ++k) {
        if (!(h[k] > 0.0)) return refuse(MetricStatus::invalid_input, "spacing not positive");
        if (k > 0 && !(h[k] > h[k - 1])) {
            return refuse(MetricStatus::invalid_input, "spacings not strictly increasing");
        }
    }
    const double r = h[1] / h[0];
    s.r = r;
    for (std::size_t k = 1; k + 1 < n; ++k) {
        if (std::abs(h[k + 1] / h[k] - r) > constant_ratio_tolerance * r) {
            return refuse(MetricStatus::not_evaluated, "non-constant r");
        }
    }
    std::vector<double> d(n - 1);
    for (std::size_t k = 0; k + 1 < n; ++k) {
        d[k] = f[k + 1] - f[k];
        if (d[k] == 0.0) return refuse(MetricStatus::not_evaluated, "zero difference");
        if ((d[k] > 0.0) != (d[0] > 0.0)) return refuse(MetricStatus::not_evaluated, "not monotone");
    }
    for (double x : f) {
        if (x == 0.0) return refuse(MetricStatus::not_evaluated, "a solution is zero (relative GCI undefined)");
    }
    const double ln_r = std::log(r);
    std::vector<double> p(n - 2);
    for (std::size_t k = 0; k + 2 < n; ++k) {
        p[k] = std::log(d[k + 1] / d[k]) / ln_r;
        if (!(p[k] > 0.0)) return refuse(MetricStatus::not_evaluated, "observed order <= 0");
    }
    s.orders = p;
    const double rp0 = std::pow(r, p[0]);
    s.extrapolated = f[0] + (f[0] - f[1]) / (rp0 - 1.0);
    s.gci.resize(n - 1);
    for (std::size_t k = 0; k + 1 < n; ++k) {
        s.gci[k] = gci_safety_factor_three_grids * std::abs(d[k] / f[k]) / (rp0 - 1.0);
    }
    s.three_grid_ratio = s.gci[1] / (rp0 * s.gci[0]);
    if (n < 4) {
        s.asymptotic_status = MetricStatus::not_evaluated;
        s.asymptotic_reason = "needs 4 grids";
        return s;
    }
    // rho_k = GCI_{k+2,k+3}(p_k) / (r^{p_k} GCI_{k+1,k+2}(p_k)); F_s and (r^{p_k} - 1) cancel.
    for (std::size_t k = 0; k + 3 < n; ++k) {
        s.order_spread = std::max(s.order_spread, std::abs(p[k + 1] - p[k]));
        const double rho = std::abs(d[k + 2] / f[k + 2]) / (std::pow(r, p[k]) * std::abs(d[k + 1] / f[k + 1]));
        s.ratio_deviation = std::max(s.ratio_deviation, std::abs(rho - 1.0));
    }
    s.asymptotic_status = MetricStatus::ok;
    s.in_asymptotic_range = s.order_spread <= asymptotic_order_spread && s.ratio_deviation <= asymptotic_ratio_band;
    return s;
}

Metric interaction_contrast(double e_h_tol, double e_h_tolfine, double e_hfine_tol, double e_hfine_tolfine) {
    if (!finite(e_h_tol) || !finite(e_h_tolfine) || !finite(e_hfine_tol) || !finite(e_hfine_tolfine)) {
        return invalid("non-finite argument");
    }
    DD c(e_hfine_tolfine);
    c -= DD(e_hfine_tol);
    c -= DD(e_h_tolfine);
    c += DD(e_h_tol);
    return ok(c.to_double());
}

CombinedUncertainty combine_uncertainty(std::span<const double> components, double eps_int) {
    CombinedUncertainty u;
    if (!all_finite(components) || !finite(eps_int)) {
        u.status = MetricStatus::invalid_input;
        u.reason = "non-finite argument";
        return u;
    }
    DD sum(std::abs(eps_int));
    for (double x : components) sum += DD(std::abs(x));
    u.value = sum.to_double();
    u.interacting = std::abs(eps_int) > interaction_fraction * u.value;
    return u;
}

}  // namespace NiTCAD::verify::reference
