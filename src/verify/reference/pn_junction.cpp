#include "NiTCAD/verify/reference/pn_junction.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>

#include "NiTCAD/verify/reference/chebyshev.hpp"
#include "NiTCAD/verify/reference/tanh_sinh.hpp"

namespace NiTCAD::verify::reference {

namespace {

// G(t) = e^t - 1 - t: the Taylor series sum_{k >= 2} t^k / k! for |t| < 1/4 (positive terms for t > 0, alternating
// and decreasing for t < 0, so no cancellation beyond a factor 1 + O(t)); expm1(t) - t elsewhere, where the
// difference loses at most a factor of about 8.
DD G(const DD& t) {
    if (std::abs(t.hi) < 0.25) {
        DD term = t * t * DD(0.5), sum = term;
        for (int k = 3; k < 60; ++k) {
            term = term * t / DD(static_cast<double>(k));
            sum += term;
            if (std::abs(term.hi) <= 0x1p-110 * std::abs(sum.hi)) break;
        }
        return sum;
    }
    return expm1(t) - t;
}

// asinh(a) for a > 0 as log(a + sqrt(a^2 + 1)); returns e^asinh(a) as well.
std::pair<DD, DD> asinh_and_exp(const DD& a) {
    const DD e = a + sqrt(a * a + DD(1.0));
    return {log(e), e};
}

struct Side {
    DD u_b;     // bulk level
    DD n_b, p_b;
    DD N;       // |net doping|
    double sign;  // +1: n side (u increases with x toward u_b), -1: p side
    // Phi_b / N at tail t = |u - u_b|: (1/2) (du/dx_hat)^2 with x_hat = x / debye. Taken from t itself, never from
    // u - u_b (which cancels to about u_dd |u_b| / t relative when t is small, and made x(s) noisy at 1e-18).
    [[nodiscard]] DD phi_tail(const DD& t) const {
        const DD d = sign > 0 ? -t : t;  // u - u_b
        return (n_b * G(d) + p_b * G(-d)) / N;
    }
    // u at tail t: the bulk side of the junction.
    [[nodiscard]] DD u_at(const DD& t) const { return sign > 0 ? u_b - t : u_b + t; }
    // k(s) = e^s / |du/dx_hat| at t = e^s.
    [[nodiscard]] DD k(const DD& s) const {
        const DD t = exp(s);
        return t / sqrt(ldexp(phi_tail(t), 1));
    }
};

bool positive_finite(double v) { return std::isfinite(v) && v > 0.0; }

}  // namespace

std::expected<PnJunction, std::string> pn_junction(const PnJunctionInput& in) {
    if (!positive_finite(in.q_C) || !positive_finite(in.eps_F_per_cm) || !positive_finite(in.thermal_voltage_V) ||
        !positive_finite(in.ni_cm3) || !positive_finite(in.acceptors_cm3) || !positive_finite(in.donors_cm3)) {
        return std::unexpected("pn_junction: every input must be finite and positive");
    }
    PnJunction j;
    j.in = in;
    const DD ni(in.ni_cm3), VT(in.thermal_voltage_V);
    const auto [an, en] = asinh_and_exp(DD(in.donors_cm3) / (ni * DD(2.0)));
    const auto [ap, ep] = asinh_and_exp(DD(in.acceptors_cm3) / (ni * DD(2.0)));
    j.u_n = an;
    j.u_p = -ap;
    j.n_n = ni * en;  // n_i e^{u_n}
    j.p_n = ni / en;
    j.p_p = ni * ep;  // n_i e^{-u_p}
    j.n_p = ni / ep;
    j.psi_n = VT * j.u_n;
    j.psi_p = VT * j.u_p;
    j.built_in_V = j.psi_n - j.psi_p;
    const DD Nn(in.donors_cm3), Np(in.acceptors_cm3);
    // Phi_n(u) - Phi_p(u) in cm^-3 (both sides' first integrals, unscaled).
    const auto mismatch = [&](const DD& u) {
        return (j.n_n * G(u - j.u_n) + j.p_n * G(j.u_n - u)) - (j.n_p * G(u - j.u_p) + j.p_p * G(j.u_p - u));
    };
    const double t = u_dd * (std::abs(j.u_n.hi) + std::abs(j.u_p.hi));
    const BrentResult b = brent_zero(mismatch, j.u_p, j.u_n, t);
    if (b.status != BrentStatus::ok) return std::unexpected("pn_junction: Brent did not converge for u_j");
    j.u_j = b.root;
    j.brent_evaluations = b.evaluations;
    j.psi_j = VT * j.u_j;
    const DD eps(in.eps_F_per_cm), q(in.q_C);
    const DD phi_j = j.n_n * G(j.u_j - j.u_n) + j.p_n * G(j.u_n - j.u_j);
    j.field_max_V_per_cm = VT * sqrt(DD(2.0) * q * phi_j / (eps * VT));
    j.depletion_n_cm = eps * j.field_max_V_per_cm / (q * Nn);
    j.depletion_p_cm = eps * j.field_max_V_per_cm / (q * Np);
    j.debye_n_cm = sqrt(eps * VT / (q * Nn));
    j.debye_p_cm = sqrt(eps * VT / (q * Np));
    return j;
}

std::expected<std::vector<PnPoint>, std::string> pn_potential(const PnJunction& j, std::span<const double> x_cm) {
    std::vector<PnPoint> out(x_cm.size());
    const DD VT(j.in.thermal_voltage_V);
    for (const double sign : {1.0, -1.0}) {
        const Side side = sign > 0 ? Side{j.u_n, j.n_n, j.p_n, DD(j.in.donors_cm3), 1.0}
                                   : Side{j.u_p, j.n_p, j.p_p, DD(j.in.acceptors_cm3), -1.0};
        const DD debye = sign > 0 ? j.debye_n_cm : j.debye_p_cm;
        // This side's points by increasing distance from the junction.
        std::vector<std::size_t> order;
        for (std::size_t i = 0; i < x_cm.size(); ++i) {
            if (sign > 0 ? x_cm[i] >= 0.0 : x_cm[i] < 0.0) order.push_back(i);
        }
        std::sort(order.begin(), order.end(),
                  [&](std::size_t a, std::size_t b) { return std::abs(x_cm[a]) < std::abs(x_cm[b]); });
        // X(s) = int_s^{s_j} k: the scaled distance from the junction. Start at the junction.
        const DD s_j = log(sign > 0 ? j.u_n - j.u_j : j.u_j - j.u_p);
        DD s_prev = s_j, X_prev(0.0);
        for (const std::size_t i : order) {
            const DD target = DD(std::abs(x_cm[i])) / debye;
            // Safeguarded Newton on F(s) = X(s) - target, decreasing in s: bracket [lo, hi] with F(lo) >= 0 >= F(hi).
            DD hi = s_prev, lo;
            bool have_lo = false;
            DD s = s_prev, X = X_prev;
            int it = 0;
            for (; it < 200; ++it) {
                const DD F = X - target;
                if (!(F > DD(0.0))) hi = s;
                if (!(F < DD(0.0))) {
                    lo = s;
                    have_lo = true;
                }
                const DD k = side.k(s);
                if (std::abs(F.hi) <= 1e-22 * std::max(1.0, std::abs(target.hi))) break;  // far below any criterion (1e-7)
                if (have_lo && !(lo < hi)) {  // the bracket has closed above the tolerance: X(s) is noisier than it
                    return std::unexpected("pn_potential: the Newton bracket closed with |X - x| = " +
                                           std::to_string(std::abs(F.hi)) + " at x = " + std::to_string(x_cm[i]));
                }
                DD next = s + F / k;  // dX/ds = -k
                const bool inside = next < hi && (!have_lo || lo < next);
                if (!inside) next = have_lo ? ldexp(lo + hi, -1) : hi - DD(1.0);
                // X at next, integrated from the last accepted point s_prev.
                const auto integrand = [&](const DD& sig) { return side.k(sig); };
                const bool down = next < s_prev;
                const TanhSinhResult q = down ? tanh_sinh(integrand, next, s_prev, 1e-28, 1e-26)
                                              : tanh_sinh(integrand, s_prev, next, 1e-28, 1e-26);
                if (q.status != TanhSinhStatus::ok) {
                    return std::unexpected("pn_potential: tanh-sinh did not converge at x = " +
                                           std::to_string(x_cm[i]));
                }
                X = down ? X_prev + q.value : X_prev - q.value;
                s = next;
            }
            if (it == 200) return std::unexpected("pn_potential: Newton did not converge at x = " + std::to_string(x_cm[i]));
            s_prev = s;
            X_prev = X;
            const DD tail = exp(s);
            out[i] = {VT * side.u_at(tail), tail};
        }
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
