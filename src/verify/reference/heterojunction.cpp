#include "NiTCAD/verify/reference/heterojunction.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>

#include "NiTCAD/verify/reference/brent.hpp"
#include "NiTCAD/verify/reference/chebyshev.hpp"
#include "NiTCAD/verify/reference/first_integral.hpp"

namespace NiTCAD::verify::reference {

namespace {

bool positive_finite(double v) { return std::isfinite(v) && v > 0.0; }

bool valid(const HeteroSide& s) {
    return positive_finite(s.eps_F_per_cm) && std::isfinite(s.affinity_eV) && positive_finite(s.band_gap_eV) &&
           positive_finite(s.Nc_cm3) && positive_finite(s.Nv_cm3) && std::isfinite(s.net_doping_cm3) &&
           s.net_doping_cm3 != 0.0;
}

// n_i = sqrt(N_c N_v) e^(-E_g / (2 V_T)), formed in double-double and rounded once.
double intrinsic(const HeteroSide& s, double VT) {
    return (sqrt(DD(s.Nc_cm3) * DD(s.Nv_cm3)) * exp(-DD(s.band_gap_eV) / (DD(2.0) * DD(VT)))).to_double();
}

struct Sides {
    NeutralBulk left, right;
    DD v_bi;  // V_bi / V_T
};

Sides sides_of(const Heterojunction& j) {
    const DD VT(j.in.thermal_voltage_V);
    Sides s{neutral_bulk(j.ni_left_cm3, j.in.left.net_doping_cm3),
            neutral_bulk(j.ni_right_cm3, j.in.right.net_doping_cm3), DD(0.0)};
    s.v_bi = (log(s.right.n_b / DD(j.in.right.Nc_cm3)) - DD(j.in.right.affinity_eV) / VT) -
             (log(s.left.n_b / DD(j.in.left.Nc_cm3)) - DD(j.in.left.affinity_eV) / VT);
    return s;
}

// Phi_b at the deviation t from the bulk's own level, from t itself.
DD phi_at(const NeutralBulk& b, const DD& t) {
    return b.n_b * poisson_boltzmann_G(t) + b.p_b * poisson_boltzmann_G(-t);
}

// Dense LU with partial pivoting in double; false on a zero pivot. Solves a x = b in place. (pn_junction.cpp keeps
// its own copy in its unnamed namespace; the two references share no code by design, I3.)
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

std::expected<Heterojunction, std::string> heterojunction(const HeterojunctionInput& in) {
    if (!positive_finite(in.q_C) || !positive_finite(in.thermal_voltage_V) || !valid(in.left) || !valid(in.right)) {
        return std::unexpected("heterojunction: an input is not finite, not positive where it must be, or a doping "
                               "is zero");
    }
    Heterojunction j;
    j.in = in;
    j.ni_left_cm3 = intrinsic(in.left, in.thermal_voltage_V);
    j.ni_right_cm3 = intrinsic(in.right, in.thermal_voltage_V);
    if (!positive_finite(j.ni_left_cm3) || !positive_finite(j.ni_right_cm3)) {
        return std::unexpected("heterojunction: an intrinsic density is outside the normal range");
    }
    const Sides s = sides_of(j);
    const DD VT(in.thermal_voltage_V), q(in.q_C);
    const DD e1(in.left.eps_F_per_cm), e2(in.right.eps_F_per_cm);
    DD t1(0.0);
    if (s.v_bi.hi != 0.0) {
        const auto mismatch = [&](const DD& t) { return e1 * phi_at(s.left, t) - e2 * phi_at(s.right, t - s.v_bi); };
        const DD lo = s.v_bi.hi > 0.0 ? DD(0.0) : s.v_bi, hi = s.v_bi.hi > 0.0 ? s.v_bi : DD(0.0);
        const BrentResult b = brent_zero(mismatch, lo, hi, u_dd * std::abs(s.v_bi.hi));
        if (b.status != BrentStatus::ok) return std::unexpected("heterojunction: Brent did not converge");
        t1 = b.root;
        j.brent_evaluations = b.evaluations;
    }
    j.built_in_V = VT * s.v_bi;
    j.interface_V = VT * t1;
    j.bend_left_V = j.interface_V;
    j.bend_right_V = j.built_in_V - j.interface_V;
    j.displacement_C_per_cm2 = sqrt(DD(2.0) * q * VT * e1 * phi_at(s.left, t1));
    j.field_left_V_per_cm = j.displacement_C_per_cm2 / e1;
    j.field_right_V_per_cm = j.displacement_C_per_cm2 / e2;
    // E_c - E_F = -q phi - chi with phi = phi_b1 + bend_left, phi_b1 = V_T ln(n_b1 / N_c1) - chi_1.
    j.ec_left_eV = -(VT * log(s.left.n_b / DD(in.left.Nc_cm3))) - j.bend_left_V;
    j.ec_right_eV = j.ec_left_eV + DD(in.left.affinity_eV) - DD(in.right.affinity_eV);
    j.debye_left_cm = sqrt(e1 * VT / (q * s.left.doping));
    j.debye_right_cm = sqrt(e2 * VT / (q * s.right.doping));
    return j;
}

std::expected<std::vector<HeteroPoint>, std::string> hetero_potential(const Heterojunction& j,
                                                                      std::span<const double> x_cm) {
    std::vector<HeteroPoint> out(x_cm.size());
    const Sides s = sides_of(j);
    const DD VT(j.in.thermal_voltage_V);
    const DD t1 = j.interface_V / VT;
    for (const bool right : {false, true}) {
        const NeutralBulk& bulk = right ? s.right : s.left;
        const DD t_surface = right ? t1 - s.v_bi : t1;
        std::vector<std::size_t> index;
        std::vector<double> distance;
        for (std::size_t i = 0; i < x_cm.size(); ++i) {
            if (right ? x_cm[i] >= 0.0 : x_cm[i] < 0.0) {
                index.push_back(i);
                distance.push_back(std::abs(x_cm[i]));
            }
        }
        const auto prof =
            bulk_profile(bulk, bulk.u_b + t_surface, right ? j.debye_right_cm : j.debye_left_cm, distance);
        if (!prof) return std::unexpected("hetero_potential: " + prof.error());
        const bool up = t_surface.hi > 0.0 || (t_surface.hi == 0.0 && t_surface.lo > 0.0);
        for (std::size_t k = 0; k < index.size(); ++k) {
            const DD t = up ? (*prof)[k].tail : -(*prof)[k].tail;  // the deviation, from the tail itself
            out[index[k]] = {right ? j.built_in_V + VT * t : VT * t, (*prof)[k].tail};
        }
    }
    return out;
}

std::expected<HeteroCollocation, std::string> hetero_collocation(const Heterojunction& j, double x_left_cm,
                                                                 double x_right_cm, int subdomains, int degree,
                                                                 std::span<const double> query_cm) {
    if (!(x_left_cm < 0.0 && x_right_cm > 0.0) || subdomains < 1 || degree < 2) {
        return std::unexpected("hetero_collocation: needs x_left < 0 < x_right, subdomains >= 1, degree >= 2");
    }
    // Scaled: x_hat = x / lambda, lambda the left side's Debye length, densities in units of N_1 = |C_1|:
    // y'' = (eps_1 / eps_k) (n_bk e^(y - y_bk) - p_bk e^-(y - y_bk) - C_k) / N_1, y_b1 = 0, y_b2 = v_bi.
    const Sides s = sides_of(j);
    const double lambda = j.debye_left_cm.to_double();
    const double N1 = s.left.doping.to_double();
    const double vbi = s.v_bi.to_double();
    struct Piece {
        double a;        // eps_1 / eps_k
        double nb, pb;   // / N_1
        double C;        // / N_1
        double yb;
        double eps;      // eps_k / eps_1
    };
    const double er = j.in.right.eps_F_per_cm / j.in.left.eps_F_per_cm;
    const Piece left{1.0, s.left.n_b.to_double() / N1, s.left.p_b.to_double() / N1, j.in.left.net_doping_cm3 / N1,
                     0.0, 1.0};
    const Piece right{1.0 / er, s.right.n_b.to_double() / N1, s.right.p_b.to_double() / N1,
                      j.in.right.net_doping_cm3 / N1, vbi, er};
    const std::size_t P = 2 * static_cast<std::size_t>(subdomains);
    const std::size_t S = static_cast<std::size_t>(subdomains);
    const std::size_t per = static_cast<std::size_t>(degree) + 1;
    const std::size_t n = P * per;
    std::vector<double> ends;
    for (int k = 0; k <= subdomains; ++k) ends.push_back(x_left_cm * (1.0 - static_cast<double>(k) / subdomains));
    for (int k = 1; k <= subdomains; ++k) ends.push_back(x_right_cm * static_cast<double>(k) / subdomains);
    ends[S] = 0.0;
    const DDVector xi_dd = chebyshev_points(degree);
    const DDMatrix D_dd = chebyshev_differentiation(degree);
    std::vector<double> xi(per), D(per * per), D2(per * per, 0.0);
    for (std::size_t a = 0; a < per; ++a) {
        xi[a] = xi_dd[a].to_double();
        for (std::size_t b = 0; b < per; ++b) D[a * per + b] = D_dd(a, b).to_double();
    }
    for (std::size_t a = 0; a < per; ++a) {
        for (std::size_t b = 0; b < per; ++b) {
            double sum = 0.0;
            for (std::size_t c = 0; c < per; ++c) sum += D[a * per + c] * D[c * per + b];
            D2[a * per + b] = sum;
        }
    }
    std::vector<double> half(P), mid(P);
    for (std::size_t p = 0; p < P; ++p) {
        half[p] = 0.5 * (ends[p + 1] - ends[p]) / lambda;
        mid[p] = 0.5 * (ends[p + 1] + ends[p]) / lambda;
    }
    const auto piece = [&](std::size_t p) -> const Piece& { return p < S ? left : right; };
    const auto idx = [&](std::size_t p, std::size_t m) { return p * per + m; };
    // Initial guess: a tanh step between the bulk levels, two left Debye lengths wide.
    std::vector<double> y(n);
    for (std::size_t p = 0; p < P; ++p) {
        for (std::size_t m = 0; m < per; ++m) y[idx(p, m)] = vbi * 0.5 * (1.0 + std::tanh(0.5 * (mid[p] + half[p] * xi[m])));
    }
    HeteroCollocation out;
    const double scale = std::max(std::abs(vbi), 1.0);
    for (int it = 1; it <= 80; ++it) {
        std::vector<double> A(n * n, 0.0), r(n, 0.0);
        for (std::size_t p = 0; p < P; ++p) {
            const Piece& k = piece(p);
            const double s2 = 1.0 / (half[p] * half[p]);
            for (std::size_t m = 1; m + 1 < per; ++m) {  // interior rows: y'' - f(y) = 0
                const std::size_t row = idx(p, m);
                double d2 = 0.0;
                for (std::size_t b = 0; b < per; ++b) {
                    d2 += D2[m * per + b] * y[idx(p, b)];
                    A[row * n + idx(p, b)] = s2 * D2[m * per + b];
                }
                const double t = y[row] - k.yb;
                const double en = k.nb * std::exp(t), ep = k.pb * std::exp(-t);
                r[row] = s2 * d2 - k.a * (en - ep - k.C);
                A[row * n + row] -= k.a * (en + ep);
            }
            if (p == 0) {  // y = 0 at the left end
                const std::size_t row = idx(p, per - 1);
                A[row * n + row] = 1.0;
                r[row] = y[row];
            }
            if (p + 1 == P) {  // y = v_bi at the right end
                const std::size_t row = idx(p, 0);
                A[row * n + row] = 1.0;
                r[row] = y[row] - vbi;
            } else {  // continuity in row (p, 0), flux in row (p + 1, N)
                const std::size_t q = p + 1;
                const std::size_t rc = idx(p, 0), rf = idx(q, per - 1);
                A[rc * n + idx(p, 0)] = 1.0;
                A[rc * n + idx(q, per - 1)] = -1.0;
                r[rc] = y[idx(p, 0)] - y[idx(q, per - 1)];
                const double sp = piece(p).eps / half[p], sq = piece(q).eps / half[q];
                double dl = 0.0, dr = 0.0;
                for (std::size_t b = 0; b < per; ++b) {
                    dl += D[b] * y[idx(p, b)];
                    dr += D[(per - 1) * per + b] * y[idx(q, b)];
                    A[rf * n + idx(p, b)] += sp * D[b];
                    A[rf * n + idx(q, b)] -= sq * D[(per - 1) * per + b];
                }
                r[rf] = sp * dl - sq * dr;
            }
        }
        if (!dense_solve(A, r, n)) return std::unexpected("hetero_collocation: singular Newton matrix");
        double worst = 0.0;
        for (const double d : r) worst = std::max(worst, std::abs(d));
        const double damp = worst > 1.0 ? 1.0 / worst : 1.0;  // at most one V_T per step
        for (std::size_t k = 0; k < n; ++k) y[k] -= damp * r[k];
        out.newton_iterations = it;
        out.last_update = worst * damp;
        if (worst <= 1e-13 * scale) break;
        if (it == 80) return std::unexpected("hetero_collocation: Newton did not converge in 80 iterations");
    }
    const double VT = j.in.thermal_voltage_V;
    std::vector<DDVector> coef(P);
    double tail = 0.0;
    for (std::size_t p = 0; p < P; ++p) {
        DDVector vals(per);
        for (std::size_t m = 0; m < per; ++m) vals[m] = DD(y[idx(p, m)]);
        coef[p] = chebyshev_coefficients(vals);
        tail = std::max(tail, std::abs(coef[p][per - 1].hi) + std::abs(coef[p][per - 2].hi));
    }
    out.tail_estimate_V = tail * VT;
    for (const double x : query_cm) {
        std::size_t p = 0;
        while (p + 1 < P && x > ends[p + 1]) ++p;
        if (x >= 0.0 && p < S) p = S;  // x = 0 with the right side
        const double z = std::clamp((x / lambda - mid[p]) / half[p], -1.0, 1.0);
        out.phi_V.push_back(VT * chebyshev_evaluate(coef[p], DD(z)).to_double());
    }
    // The interface: y at the left piece's right end, and both one-sided derivatives.
    out.interface_V = VT * y[idx(S - 1, 0)];
    double dl = 0.0, dr = 0.0;
    for (std::size_t b = 0; b < per; ++b) {
        dl += D[b] * y[idx(S - 1, b)];
        dr += D[(per - 1) * per + b] * y[idx(S, b)];
    }
    const double fl = dl / half[S - 1], fr = er * dr / half[S];  // eps y' / eps_1, per lambda
    out.displacement_C_per_cm2 = j.in.left.eps_F_per_cm * VT * std::abs(fl) / lambda;
    out.flux_mismatch = fl != 0.0 ? std::abs(fl - fr) / std::abs(fl) : std::abs(fr);
    return out;
}

}  // namespace NiTCAD::verify::reference
