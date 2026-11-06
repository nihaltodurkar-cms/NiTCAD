// Convergence studies against exact solutions (study.hpp).
#include "NiTCAD/verify/study/study.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

#include "NiTCAD/verify/reference/metrics.hpp"

namespace NiTCAD::verify::study {

namespace ref = verify::reference;

double mean_spacing(const mesh::Mesh& m) {
    const auto v = m.volumes();
    const double total = std::accumulate(v.begin(), v.end(), 0.0);
    return std::pow(total / static_cast<double>(m.node_count()), 1.0 / static_cast<double>(m.dimension()));
}

OrderStudy order_study(std::span<const double> h, std::span<const double> e, std::span<const double> u_iter) {
    OrderStudy s;
    const auto invalid = [&](const char* why) {
        s.status = StudyStatus::invalid_input;
        s.reason = why;
        return s;
    };
    if (h.size() < 2 || e.size() != h.size() || u_iter.size() != h.size()) return invalid("sizes: h, e, U_iter");
    for (std::size_t k = 0; k < h.size(); ++k) {
        if (!(std::isfinite(h[k]) && h[k] > 0.0) || (k > 0 && !(h[k] > h[k - 1]))) {
            return invalid("h must be finite, positive and strictly increasing (finest first)");
        }
        if (!(std::isfinite(e[k]) && e[k] > 0.0)) return invalid("errors must be finite and positive");
        if (!(std::isfinite(u_iter[k]) && u_iter[k] >= 0.0)) return invalid("U_iter must be finite, >= 0");
    }
    s.h.assign(h.begin(), h.end());
    s.error.assign(e.begin(), e.end());
    for (std::size_t k = 0; k < h.size(); ++k) s.iterative_ratio = std::max(s.iterative_ratio, u_iter[k] / e[k]);
    for (const ref::Metric& p : ref::observed_orders_from_errors(h, e)) {
        if (p.status != ref::MetricStatus::ok || !std::isfinite(p.value)) {
            s.status = StudyStatus::not_evaluated;
            s.reason = "an observed order is not defined";
            return s;
        }
        s.orders.push_back(p.value);
    }
    s.observed = s.orders.front();
    for (std::size_t k = 0; k + 1 < s.orders.size(); ++k) {
        s.order_spread = std::max(s.order_spread, std::abs(s.orders[k + 1] - s.orders[k]));
    }
    if (static_cast<int>(h.size()) >= asymptotic_levels) {
        s.asymptotic_evaluated = true;
        // Against an exact solution the errors converge to 0, so TVR-1's GCI-ratio test (relative to a nonzero
        // solution value) does not apply: rho_k would be r^-p, not 1. The check is the order spread alone.
        s.in_asymptotic_range = s.order_spread <= ref::asymptotic_order_spread;
    }
    return s;
}

OrderVerdict judge_order(const OrderStudy& s, const OrderBand& band) {
    OrderVerdict v;
    if (s.status != StudyStatus::ok) {
        v.reason = "study not evaluated: " + s.reason;
        return v;
    }
    if (s.iterative_ratio > iterative_fraction) {
        v.reason = "iterative error above a tenth of the discretisation error (CV-Q8)";
        return v;
    }
    if (s.asymptotic_evaluated && !s.in_asymptotic_range) {
        v.reason = "not in the asymptotic range";
        return v;
    }
    v.judged = true;
    v.pass = s.observed >= band.p_th - band.below && s.observed <= band.p_th + band.above;
    v.reason = v.pass ? "" : "observed order outside [p_th - 0.2, p_th + 0.3]";
    return v;
}

}  // namespace NiTCAD::verify::study
