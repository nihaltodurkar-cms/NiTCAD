#include "NiTCAD/verify/reference/pn_junction.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>

#include "NiTCAD/verify/reference/chebyshev.hpp"
#include "NiTCAD/verify/reference/first_integral.hpp"

namespace NiTCAD::verify::reference {

namespace {

bool positive_finite(double v) { return std::isfinite(v) && v > 0.0; }

}  // namespace

std::expected<PnJunction, std::string> pn_junction(const PnJunctionInput& in) {
    if (!positive_finite(in.q_C) || !positive_finite(in.eps_F_per_cm) || !positive_finite(in.thermal_voltage_V) ||
        !positive_finite(in.ni_cm3) || !positive_finite(in.acceptors_cm3) || !positive_finite(in.donors_cm3)) {
        return std::unexpected("pn_junction: every input must be finite and positive");
    }
    PnJunction j;
    j.in = in;
    const DD VT(in.thermal_voltage_V);
    const NeutralBulk n = neutral_bulk(in.ni_cm3, in.donors_cm3);
    const NeutralBulk p = neutral_bulk(in.ni_cm3, -in.acceptors_cm3);
    j.u_n = n.u_b;
    j.u_p = p.u_b;
    j.n_n = n.n_b;
    j.p_n = n.p_b;
    j.n_p = p.n_b;
    j.p_p = p.p_b;
    j.psi_n = VT * j.u_n;
    j.psi_p = VT * j.u_p;
    j.built_in_V = j.psi_n - j.psi_p;
    // Phi_n(u) - Phi_p(u) in cm^-3: u and eps u' continuous at x = 0.
    const auto mismatch = [&](const DD& u) { return first_integral_phi(n, u) - first_integral_phi(p, u); };
    const double t = u_dd * (std::abs(j.u_n.hi) + std::abs(j.u_p.hi));
    const BrentResult b = brent_zero(mismatch, j.u_p, j.u_n, t);
    if (b.status != BrentStatus::ok) return std::unexpected("pn_junction: Brent did not converge for u_j");
    j.u_j = b.root;
    j.brent_evaluations = b.evaluations;
    j.psi_j = VT * j.u_j;
    const DD eps(in.eps_F_per_cm), q(in.q_C);
    j.field_max_V_per_cm = VT * sqrt(DD(2.0) * q * first_integral_phi(n, j.u_j) / (eps * VT));
    j.depletion_n_cm = eps * j.field_max_V_per_cm / (q * n.doping);
    j.depletion_p_cm = eps * j.field_max_V_per_cm / (q * p.doping);
    j.debye_n_cm = sqrt(eps * VT / (q * n.doping));
    j.debye_p_cm = sqrt(eps * VT / (q * p.doping));
    return j;
}

std::expected<std::vector<PnPoint>, std::string> pn_potential(const PnJunction& j, std::span<const double> x_cm) {
    std::vector<PnPoint> out(x_cm.size());
    const DD VT(j.in.thermal_voltage_V);
    for (const bool n_side : {true, false}) {
        const NeutralBulk bulk = n_side ? neutral_bulk(j.in.ni_cm3, j.in.donors_cm3)
                                        : neutral_bulk(j.in.ni_cm3, -j.in.acceptors_cm3);
        std::vector<std::size_t> index;
        std::vector<double> distance;
        for (std::size_t i = 0; i < x_cm.size(); ++i) {
            if (n_side ? x_cm[i] >= 0.0 : x_cm[i] < 0.0) {
                index.push_back(i);
                distance.push_back(std::abs(x_cm[i]));
            }
        }
        const auto prof = bulk_profile(bulk, j.u_j, n_side ? j.debye_n_cm : j.debye_p_cm, distance);
        if (!prof) return std::unexpected("pn_potential: " + prof.error());
        for (std::size_t k = 0; k < index.size(); ++k) out[index[k]] = {VT * (*prof)[k].u, (*prof)[k].tail};
    }
    return out;
}


namespace {

// Dense LU with partial pivoting in double; false on a zero pivot. Solves a x = b in place.
bool dense_solve(std::vector<double>& a, std::vector<double>& b, std::size_t n) {
    for (std::size_t c = 0; c < n; ++c) {
        std::size_t p = c;
        for (std::size_t r = c + 1; r < n; ++r) {
            if (std::abs(a[r * n + c]) > std::abs(a[p * n + c])) p = r;
        }
        if (a[p * n + c] == 0.0) return false;
        if (p != c) {
            for (std::size_t k = 0; k < n; ++k) std::swap(a[p * n + k], a[c * n + k]);
            std::swap(b[p], b[c]);
        }
        for (std::size_t r = c + 1; r < n; ++r) {
            const double f = a[r * n + c] / a[c * n + c];
            if (f == 0.0) continue;
            for (std::size_t k = c; k < n; ++k) a[r * n + k] -= f * a[c * n + k];
            b[r] -= f * b[c];
        }
    }
    for (std::size_t c = n; c-- > 0;) {
        double s = b[c];
        for (std::size_t k = c + 1; k < n; ++k) s -= a[c * n + k] * b[k];
        b[c] = s / a[c * n + c];
    }
    return true;
}

}  // namespace

std::expected<PnCollocation, std::string> pn_collocation(const PnJunction& j, double x_left_cm, double x_right_cm,
                                                         int subdomains, int degree,
                                                         std::span<const double> query_cm) {
    if (!(x_left_cm < 0.0 && x_right_cm > 0.0) || subdomains < 1 || degree < 2) {
        return std::unexpected("pn_collocation: needs x_left < 0 < x_right, subdomains >= 1, degree >= 2");
    }
    // Scaled: x_hat = x / lambda, lambda the p side's Debye length; u'' = (n_i / N_A) (e^u - e^-u) - C / N_A.
    const double lambda = j.debye_p_cm.to_double();
    const double ni = j.in.ni_cm3, NA = j.in.acceptors_cm3, ND = j.in.donors_cm3;
    const double up = j.u_p.to_double(), un = j.u_n.to_double();
    const int P = 2 * subdomains;  // pieces, left to right
    const int N = degree;
    const std::size_t per = static_cast<std::size_t>(N + 1);
    const std::size_t n = static_cast<std::size_t>(P) * per;
    // Piece ends.
    std::vector<double> ends;
    for (int k = 0; k <= subdomains; ++k) ends.push_back(x_left_cm * (1.0 - static_cast<double>(k) / subdomains));
    for (int k = 1; k <= subdomains; ++k) ends.push_back(x_right_cm * static_cast<double>(k) / subdomains);
    ends[static_cast<std::size_t>(subdomains)] = 0.0;
    // cheb.m points xi_m = cos(m pi / N) (xi_0 = 1: the piece's right end) and D, D^2, in double.
    const DDVector xi_dd = chebyshev_points(N);
    const DDMatrix D_dd = chebyshev_differentiation(N);
    std::vector<double> xi(per), D(per * per), D2(per * per, 0.0);
    for (std::size_t a = 0; a < per; ++a) {
        xi[a] = xi_dd[a].to_double();
        for (std::size_t b = 0; b < per; ++b) D[a * per + b] = D_dd(a, b).to_double();
    }
    for (std::size_t a = 0; a < per; ++a) {
        for (std::size_t b = 0; b < per; ++b) {
            double s = 0.0;
            for (std::size_t c = 0; c < per; ++c) s += D[a * per + c] * D[c * per + b];
            D2[a * per + b] = s;
        }
    }
    std::vector<double> half(static_cast<std::size_t>(P)), mid(static_cast<std::size_t>(P));  // scaled
    for (int p = 0; p < P; ++p) {
        half[static_cast<std::size_t>(p)] = 0.5 * (ends[static_cast<std::size_t>(p) + 1] - ends[static_cast<std::size_t>(p)]) / lambda;
        mid[static_cast<std::size_t>(p)] = 0.5 * (ends[static_cast<std::size_t>(p) + 1] + ends[static_cast<std::size_t>(p)]) / lambda;
    }
    // Initial guess: a tanh step between the bulk levels, width two p-side Debye lengths.
    std::vector<double> u(n);
    for (std::size_t p = 0; p < static_cast<std::size_t>(P); ++p) {
        for (std::size_t m = 0; m < per; ++m) {
            const double x = mid[p] + half[p] * xi[m];
            u[p * per + m] = up + (un - up) * 0.5 * (1.0 + std::tanh(0.5 * x));
        }
    }
    PnCollocation out;
    const auto idx = [&](std::size_t p, std::size_t m) { return p * per + m; };
    for (int it = 1; it <= 50; ++it) {
        std::vector<double> A(n * n, 0.0), r(n, 0.0);
        for (std::size_t p = 0; p < static_cast<std::size_t>(P); ++p) {
            const double s1 = 1.0 / half[p], s2 = s1 * s1;
            const double C = p < static_cast<std::size_t>(subdomains) ? -NA : ND;  // net doping
            // Interior collocation rows m = 1 ... N-1: u'' - f(u) = 0.
            for (std::size_t m = 1; m + 1 < per; ++m) {
                const std::size_t row = idx(p, m);
                double d2 = 0.0;
                for (std::size_t b = 0; b < per; ++b) {
                    d2 += D2[m * per + b] * u[idx(p, b)];
                    A[row * n + idx(p, b)] = s2 * D2[m * per + b];
                }
                const double v = u[row];
                const double f = (ni / NA) * (std::exp(v) - std::exp(-v)) - C / NA;
                r[row] = s2 * d2 - f;
                A[row * n + row] -= (ni / NA) * (std::exp(v) + std::exp(-v));
            }
            // Row m = 0 (right end) and m = N (left end): Dirichlet at the outer ends, patching inside.
            if (p == 0) {  // left end of the first piece: u = u_p
                const std::size_t row = idx(p, per - 1);
                A[row * n + row] = 1.0;
                r[row] = u[row] - up;
            }
            if (p + 1 == static_cast<std::size_t>(P)) {  // right end of the last piece: u = u_n
                const std::size_t row = idx(p, 0);
                A[row * n + row] = 1.0;
                r[row] = u[row] - un;
            } else {
                // Between piece p (its m = 0) and p + 1 (its m = N): continuity in row idx(p, 0), flux in row
                // idx(p + 1, N).
                const std::size_t q = p + 1;
                const std::size_t rc = idx(p, 0), rf = idx(q, per - 1);
                A[rc * n + idx(p, 0)] = 1.0;
                A[rc * n + idx(q, per - 1)] = -1.0;
                r[rc] = u[idx(p, 0)] - u[idx(q, per - 1)];
                const double sp = 1.0 / half[p], sq = 1.0 / half[q];
                double dl = 0.0, dr = 0.0;
                for (std::size_t b = 0; b < per; ++b) {
                    dl += D[0 * per + b] * u[idx(p, b)];
                    dr += D[(per - 1) * per + b] * u[idx(q, b)];
                    A[rf * n + idx(p, b)] += sp * D[0 * per + b];
                    A[rf * n + idx(q, b)] -= sq * D[(per - 1) * per + b];
                }
                r[rf] = sp * dl - sq * dr;
            }
        }
        if (!dense_solve(A, r, n)) return std::unexpected("pn_collocation: singular Newton matrix");
        double worst = 0.0;
        for (const double d : r) worst = std::max(worst, std::abs(d));
        const double damp = worst > 1.0 ? 1.0 / worst : 1.0;  // at most one V_T per step
        for (std::size_t k = 0; k < n; ++k) u[k] -= damp * r[k];
        out.newton_iterations = it;
        out.last_update = worst * damp;
        if (worst <= 1e-13 * std::max(std::abs(up), std::abs(un))) break;
        if (it == 50) return std::unexpected("pn_collocation: Newton did not converge in 50 iterations");
    }
    const double VT = j.in.thermal_voltage_V;
    std::vector<DDVector> coef(static_cast<std::size_t>(P));
    double tail = 0.0;
    for (std::size_t p = 0; p < static_cast<std::size_t>(P); ++p) {
        DDVector vals(per);
        for (std::size_t m = 0; m < per; ++m) vals[m] = DD(u[idx(p, m)]);
        coef[p] = chebyshev_coefficients(vals);
        tail = std::max(tail, std::abs(coef[p][per - 1].hi) + std::abs(coef[p][per - 2].hi));
    }
    out.tail_estimate_V = tail * VT;
    for (const double x : query_cm) {
        std::size_t p = 0;
        while (p + 1 < static_cast<std::size_t>(P) && x > ends[p + 1]) ++p;
        const double z = std::clamp((x / lambda - mid[p]) / half[p], -1.0, 1.0);
        out.psi_V.push_back(VT * chebyshev_evaluate(coef[p], DD(z)).to_double());
    }
    return out;
}

}  // namespace NiTCAD::verify::reference
