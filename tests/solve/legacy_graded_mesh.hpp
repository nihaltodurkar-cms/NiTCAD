// Test helper: a port of the legacy pytcad/mesh.py graded_mesh (one focus point), so the Unit 9
// gates run on the legacy fixtures' own meshes. Checked against the legacy output in bias_test.cpp.
//
// Spacing target s(x) = min(h_max, h_min + (ratio - 1) |x - x_focus|); nodes at equal increments of
// the arc length t(x) = integral dx / s(x) (trapezoid on a dense uniform sampling), then the cell
// sizes are gradient-limited to `ratio` by repeated forward/backward cumulative minima in log space
// and rescaled to span [0, L].
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

inline std::vector<double> legacy_graded_mesh(double L, double x_focus, double h_min, double h_max,
                                              double ratio = 1.15) {
    h_max = std::max(h_max, h_min);
    const double g = ratio - 1.0;
    const double m_uncapped = 50.0 * L / std::max(h_min, 1e-30);
    const auto m =
        static_cast<std::size_t>(std::min(2000001.0, std::max(2001.0, m_uncapped) + 1.0));

    // numpy.linspace(0, L, m): i * step, with the last point exactly L.
    const double step = L / static_cast<double>(m - 1);
    std::vector<double> xs(m), t(m, 0.0);
    for (std::size_t i = 0; i < m; ++i) xs[i] = static_cast<double>(i) * step;
    xs[m - 1] = L;
    double previous = 1.0 / std::min(h_max, h_min + g * std::abs(xs[0] - x_focus));
    for (std::size_t i = 1; i < m; ++i) {
        const double inv = 1.0 / std::min(h_max, h_min + g * std::abs(xs[i] - x_focus));
        t[i] = t[i - 1] + 0.5 * (inv + previous) * (xs[i] - xs[i - 1]);
        previous = inv;
    }
    const auto cells = static_cast<std::size_t>(std::max(1.0, std::ceil(t[m - 1])));
    const double factor = static_cast<double>(cells) / t[m - 1];
    for (double& v : t) v *= factor;

    // numpy.interp(arange(cells + 1), t, xs).
    std::vector<double> nodes(cells + 1);
    std::size_t j = 0;
    for (std::size_t k = 0; k <= cells; ++k) {
        const double target = static_cast<double>(k);
        while (j + 2 < m && t[j + 1] <= target) ++j;
        if (target >= t[m - 1]) {
            nodes[k] = xs[m - 1];
        } else {
            const double slope = (xs[j + 1] - xs[j]) / (t[j + 1] - t[j]);
            nodes[k] = slope * (target - t[j]) + xs[j];
        }
    }
    nodes[0] = 0.0;
    nodes[cells] = L;
    if (cells < 3 || g <= 0.0) return nodes;

    std::vector<double> h(cells), lh(cells);
    for (std::size_t k = 0; k < cells; ++k) h[k] = nodes[k + 1] - nodes[k];
    const double lr = std::log(ratio);
    for (int sweep = 0; sweep < 50; ++sweep) {
        for (std::size_t k = 0; k < cells; ++k) lh[k] = std::log(h[k]);
        double low = lh[0];
        for (std::size_t k = 0; k < cells; ++k) {  // forward cumulative minimum
            low = std::min(low, lh[k] - static_cast<double>(k) * lr);
            lh[k] = low + static_cast<double>(k) * lr;
        }
        low = lh[cells - 1] + static_cast<double>(cells - 1) * lr;
        for (std::size_t k = cells; k-- > 0;) {  // backward cumulative minimum
            low = std::min(low, lh[k] + static_cast<double>(k) * lr);
            lh[k] = low - static_cast<double>(k) * lr;
        }
        double sum = 0.0;
        for (std::size_t k = 0; k < cells; ++k) {
            h[k] = std::exp(lh[k]);
            sum += h[k];
        }
        double worst = 0.0;
        for (std::size_t k = 0; k < cells; ++k) {
            h[k] *= L / sum;
            if (k > 0) worst = std::max({worst, h[k] / h[k - 1], h[k - 1] / h[k]});
        }
        if (worst <= ratio * (1.0 + 1e-12)) break;
    }
    for (std::size_t k = 0; k < cells; ++k) nodes[k + 1] = nodes[k] + h[k];
    nodes[cells] = L;
    return nodes;
}
