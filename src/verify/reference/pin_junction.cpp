#include "NiTCAD/verify/reference/pin_junction.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <vector>

#include "NiTCAD/verify/reference/brent.hpp"
#include "NiTCAD/verify/reference/chebyshev.hpp"
#include "NiTCAD/verify/reference/first_integral.hpp"
#include "NiTCAD/verify/reference/gauss_kronrod.hpp"

namespace NiTCAD::verify::reference {

namespace {

constexpr double quadrature_rtol = 1e-26;
constexpr double deepest_t = -700.0;  // ln of the smallest distance to a bulk level pin_level reaches

bool positive(double v) { return std::isfinite(v) && v > 0.0; }

// p G(a) for a >= 0 with the density p = e^lp carried as a logarithm (as multiplication.cpp).
DD minority(const DD& lp, const DD& a) {
    if (lp.hi > -300.0) return exp(lp) * poisson_boltzmann_G(a);
    return exp(lp + a) - exp(lp) * (DD(1.0) + a);
}

// p expm1(a) with p = e^lp.
DD scaled_expm1(const DD& lp, const DD& a) {
    if (lp.hi > -300.0) return exp(lp) * expm1(a);
    return exp(lp + a) - exp(lp);
}

DD phi_p(const PinJunction& j, const DD& s) { return j.p_p * poisson_boltzmann_G(-s) + minority(j.ln_np, s); }
DD phi_n(const PinJunction& j, const DD& r) { return j.n_n * poisson_boltzmann_G(-r) + minority(j.ln_pn, r); }
DD phi_i(const PinJunction& j, const DD& s0, const DD& phi0, const DD& s) {
    // The electrons' term without cancellation near s_0 (expm1), and as a plain difference once s - s_0 > 1, where
    // nothing cancels and expm1 would overflow at a large barrier.
    const DD electrons = (s - s0).hi < 1.0 ? j.n_n * exp(s0 - j.barrier) * expm1(s - s0)
                                           : j.n_n * (exp(s - j.barrier) - exp(s0 - j.barrier));
    return phi0 + electrons + j.p_p * exp(-s0) * expm1(s0 - s) - DD(j.in.intrinsic_net_cm3) * (s - s0);
}
DD field_of(const PinJunction& j, const DD& phi) { return j.VT * sqrt(j.scale * phi); }

// The level where dPhi_i/ds = n - p - C_i vanishes: n_p w - p_p / w = C_i with w = e^s (n_p p_p = m^2).
DD stationary(const PinJunction& j) {
    const double C = j.in.intrinsic_net_cm3;
    if (C == 0.0) return ldexp(log(j.p_p) - j.ln_np, -1);
    const DD m2 = exp(ldexp(j.ln_m, 1));
    const DD R = sqrt(DD(C) * DD(C) + DD(4.0) * m2);
    if (C > 0.0) return log(DD(C) + R) - ln2() - j.ln_np;
    return log(ldexp(j.p_p, 1)) - log(R - DD(C));
}

// The n edge for a trial s_0 (header comment), and the i region's length.
struct Trial {
    DD sd, phi0, L;
    double abs_error = 0.0;
    bool ok = false;  // false: Phi_i not positive on [s_0, s_d], s_d outside [s_0, B), or the quadrature failing
};

Trial trial(const PinJunction& j, const DD& s0) {
    Trial t;
    t.phi0 = phi_p(j, s0);
    const DD slope = DD(j.in.donors_cm3) - DD(j.in.intrinsic_net_cm3);
    DD sd = s0 + (phi_n(j, j.barrier - s0) - t.phi0) / slope;
    sd = sd - (phi_i(j, s0, t.phi0, sd) - phi_n(j, j.barrier - sd)) / slope;
    t.sd = sd;
    if (!(t.phi0.hi > 0.0) || sd.hi < s0.hi || !(sd.hi < j.barrier.hi)) return t;
    if (!(phi_i(j, s0, t.phi0, sd).hi > 0.0)) return t;
    const DD st = stationary(j);
    if (st.hi > s0.hi && st.hi < sd.hi && !(phi_i(j, s0, t.phi0, st).hi > 0.0)) return t;
    if (sd.hi == s0.hi && sd.lo == s0.lo) {
        t.L = DD(0.0);
        t.ok = true;
        return t;
    }
    const AdaptiveResult a =
        qag15([&](const DD& s) { return j.VT / field_of(j, phi_i(j, s0, t.phi0, s)); }, s0, sd, 0.0, quadrature_rtol);
    if (a.status != QuadratureStatus::ok) return t;
    t.L = a.value;
    t.abs_error = a.abs_error;
    t.ok = true;
    return t;
}

// G(t) / t^2, by its series where t^2 would underflow or G would cancel.
DD g_over_t2(const DD& t) {
    if (std::abs(t.hi) < 1e-8) return DD(0.5) + t * (DD(1.0) / DD(6.0) + t * (DD(1.0) / DD(24.0) + t / DD(120.0)));
    return poisson_boltzmann_G(t) / (t * t);
}
// Phi / c^2 of a side with majority density maj and minority e^lp, c the distance to its bulk level.
DD phi_over_c2(const DD& maj, const DD& lp, const DD& c) {
    const DD minor = c.hi < 1.0 ? exp(lp) * g_over_t2(c) : minority(lp, c) / (c * c);
    return maj * g_over_t2(-c) + minor;
}
// V_T e^t / E = 1 / sqrt(scale Phi / c^2) on one side, t = ln c the logarithm of the distance to that side's bulk
// level: bounded as c -> 0 (where it tends to the bulk's Debye length), and free of underflow.
DD p_integrand(const PinJunction& j, const DD& t) {
    return DD(1.0) / sqrt(j.scale * phi_over_c2(j.p_p, j.ln_np, exp(t)));
}
DD n_integrand(const PinJunction& j, const DD& t) {
    return DD(1.0) / sqrt(j.scale * phi_over_c2(j.n_n, j.ln_pn, exp(t)));
}

std::expected<DD, std::string> integral(const AdaptiveResult& a) {
    if (a.status != QuadratureStatus::ok) return std::unexpected("pin_junction: a position quadrature did not converge");
    return a.value;
}

// x at t = ln s on the p side (x <= 0), and at t = ln (B - s) on the n side (x >= d).
std::expected<DD, std::string> p_position(const PinJunction& j, const DD& t) {
    const DD top = log(j.s0);
    if (!(t < top)) return DD(0.0);
    const auto v = integral(qag15([&](const DD& z) { return p_integrand(j, z); }, t, top, 0.0, quadrature_rtol));
    if (!v) return v;
    return -*v;
}
std::expected<DD, std::string> n_position(const PinJunction& j, const DD& t) {
    const DD top = log(j.barrier - j.sd);
    const DD d(j.in.i_length_cm);
    if (!(t < top)) return d;
    const auto v = integral(qag15([&](const DD& z) { return n_integrand(j, z); }, t, top, 0.0, quadrature_rtol));
    if (!v) return v;
    return d + *v;
}

// Gaussian elimination with partial pivoting on a dense n x n system in double; b becomes the solution.
bool dense_solve(std::vector<double>& A, std::vector<double>& b, std::size_t n) {
    for (std::size_t c = 0; c < n; ++c) {
        std::size_t piv = c;
        for (std::size_t r = c + 1; r < n; ++r) {
            if (std::abs(A[r * n + c]) > std::abs(A[piv * n + c])) piv = r;
        }
        if (A[piv * n + c] == 0.0) return false;
        if (piv != c) {
            for (std::size_t k = 0; k < n; ++k) std::swap(A[c * n + k], A[piv * n + k]);
            std::swap(b[c], b[piv]);
        }
        const double inv = 1.0 / A[c * n + c];
        for (std::size_t r = c + 1; r < n; ++r) {
            const double f = A[r * n + c] * inv;
            if (f == 0.0) continue;
            for (std::size_t k = c; k < n; ++k) A[r * n + k] -= f * A[c * n + k];
            b[r] -= f * b[c];
        }
    }
    for (std::size_t c = n; c-- > 0;) {
        double sum = b[c];
        for (std::size_t k = c + 1; k < n; ++k) sum -= A[c * n + k] * b[k];
        b[c] = sum / A[c * n + c];
    }
    return true;
}

}  // namespace

std::expected<PinJunction, std::string> pin_junction(const PinInput& in, double bias_V) {
    if (!positive(in.q_C) || !positive(in.eps_F_per_cm) || !positive(in.thermal_voltage_V) || !positive(in.ni_cm3) ||
        !positive(in.acceptors_cm3) || !positive(in.donors_cm3) || !positive(in.i_length_cm) ||
        !std::isfinite(in.intrinsic_net_cm3) || !(in.intrinsic_net_cm3 > -in.acceptors_cm3) ||
        !(in.intrinsic_net_cm3 < in.donors_cm3)) {
        return std::unexpected("pin_junction: an input is not finite, or not positive where it must be, or C_i is "
                               "outside (-N_A, N_D)");
    }
    if (!std::isfinite(bias_V) || bias_V > 0.0) {
        return std::unexpected("pin_junction: the bias must be finite and not forward");
    }
    PinJunction j;
    j.in = in;
    j.VT = DD(in.thermal_voltage_V);
    j.scale = DD(2.0) * DD(in.q_C) / (DD(in.eps_F_per_cm) * j.VT);
    j.ln_m = log(DD(in.ni_cm3)) + DD(bias_V) / ldexp(j.VT, 1);
    const DD m = exp(j.ln_m);  // 0 once it underflows: the majority densities are then the dopings exactly
    const DD half_a = ldexp(DD(in.acceptors_cm3), -1), half_d = ldexp(DD(in.donors_cm3), -1);
    j.p_p = half_a + sqrt(half_a * half_a + m * m);
    j.n_n = half_d + sqrt(half_d * half_d + m * m);
    j.ln_np = ldexp(j.ln_m, 1) - log(j.p_p);
    j.ln_pn = ldexp(j.ln_m, 1) - log(j.n_n);
    j.barrier = log(j.n_n) + log(j.p_p) - ldexp(j.ln_m, 1);
    if (!(j.barrier.hi > 0.0)) return std::unexpected("pin_junction: no barrier at this bias");

    // The two-region level (d = 0): Phi_p(s) = Phi_n(B - s). L vanishes there.
    const BrentResult two = brent_zero([&](const DD& s) { return phi_p(j, s) - phi_n(j, j.barrier - s); }, DD(0.0),
                                       j.barrier, u_dd * j.barrier.hi);
    if (two.status != BrentStatus::ok) return std::unexpected("pin_junction: Brent did not converge for the two-region level");
    const DD d(in.i_length_cm);
    const auto f = [&](const DD& s0) -> std::optional<DD> {
        ++j.evaluations;
        const Trial t = trial(j, s0);
        if (!t.ok) return std::nullopt;
        return t.L - d;
    };
    // The low end of the bracket: halve s_0 until L exceeds d, bisecting towards the valid side once a trial leaves
    // it (L grows without bound as s_0 approaches that side's edge).
    DD good = two.root, bad, lo;
    bool have_bad = false, found = false;
    DD s = ldexp(good, -1);
    for (int k = 0; k < 600; ++k) {
        const auto v = f(s);
        if (v && v->hi > 0.0) {
            lo = s;
            found = true;
            break;
        }
        if (v) good = s;
        else {
            bad = s;
            have_bad = true;
        }
        s = have_bad ? ldexp(bad + good, -1) : ldexp(s, -1);
        if (!have_bad && s.hi < 1e-300) break;
    }
    if (!found) return std::unexpected("pin_junction: the i region's length is not bracketed");
    bool lost = false;
    const BrentResult root = brent_zero(
        [&](const DD& s0) {
            const auto v = f(s0);
            if (!v) lost = true;
            return v ? *v : DD(0.0);
        },
        lo, good, 4.0 * u_dd * good.hi);
    if (lost || root.status != BrentStatus::ok) return std::unexpected("pin_junction: Brent did not converge for s_0");
    j.s0 = root.root;
    const Trial t = trial(j, j.s0);
    if (!t.ok) return std::unexpected("pin_junction: the root's length did not evaluate");
    j.sd = t.sd;
    j.phi0 = t.phi0;
    j.length_residual_cm = std::abs((t.L - d).to_double());
    const DD phi_end = phi_n(j, j.barrier - j.sd);
    j.matching_residual = std::abs(((phi_i(j, j.s0, j.phi0, j.sd) - phi_end) / phi_end).to_double());
    j.field_p_edge_V_per_cm = field_of(j, j.phi0);
    j.field_n_edge_V_per_cm = field_of(j, phi_end);
    DD low = j.field_p_edge_V_per_cm < j.field_n_edge_V_per_cm ? j.field_p_edge_V_per_cm : j.field_n_edge_V_per_cm;
    const DD st = stationary(j);
    if (st > j.s0 && st < j.sd) {
        const DD e = field_of(j, phi_i(j, j.s0, j.phi0, st));
        if (e < low) low = e;
    }
    j.field_min_i_V_per_cm = low;
    const DD eps(in.eps_F_per_cm), q(in.q_C);
    j.depletion_p_cm = eps * j.field_p_edge_V_per_cm / (q * DD(in.acceptors_cm3));
    j.depletion_n_cm = eps * j.field_n_edge_V_per_cm / (q * DD(in.donors_cm3));
    return j;
}

DD pin_field(const PinJunction& j, const DD& s) {
    if (!(s > j.s0)) return field_of(j, phi_p(j, s));
    if (!(s > j.sd)) return field_of(j, phi_i(j, j.s0, j.phi0, s));
    return field_of(j, phi_n(j, j.barrier - s));
}

DD pin_field_slope(const PinJunction& j, const DD& s) {
    const DD k = DD(j.in.q_C) / DD(j.in.eps_F_per_cm);
    // Each side from its bulk, so that neutrality cancels exactly: p side n - p + N_A = n_p expm1(s) - p_p expm1(-s);
    // n side n - p - N_D = n_n expm1(-r) - p_n expm1(r).
    if (!(s > j.s0)) return k * (scaled_expm1(j.ln_np, s) - j.p_p * expm1(-s));
    if (!(s > j.sd)) {
        return k * (exp(log(j.n_n) + s - j.barrier) - exp(log(j.p_p) - s) - DD(j.in.intrinsic_net_cm3));
    }
    const DD r = j.barrier - s;
    return k * (j.n_n * expm1(-r) - scaled_expm1(j.ln_pn, r));
}

std::expected<DD, std::string> pin_position(const PinJunction& j, const DD& s) {
    if (!(s.hi > 0.0 && s.hi < j.barrier.hi)) return std::unexpected("pin_position: the level is outside (0, B)");
    if (!(s > j.s0)) return p_position(j, log(s));
    if (!(s > j.sd)) {
        if (s == j.s0) return DD(0.0);
        return integral(qag15([&](const DD& z) { return j.VT / field_of(j, phi_i(j, j.s0, j.phi0, z)); }, j.s0, s, 0.0,
                              quadrature_rtol));
    }
    return n_position(j, log(j.barrier - s));
}

std::expected<DD, std::string> pin_level(const PinJunction& j, double x_cm) {
    const double d = j.in.i_length_cm;
    if (!std::isfinite(x_cm)) return std::unexpected("pin_level: the position is not finite");
    std::string error;
    const auto check = [&](const std::expected<DD, std::string>& v) {
        if (!v) {
            error = v.error();
            return DD(0.0);
        }
        return *v;
    };
    if (x_cm == 0.0) return j.s0;
    if (x_cm == d) return j.sd;
    if (x_cm < 0.0) {
        const auto g = [&](const DD& t) { return check(p_position(j, t)) - DD(x_cm); };
        const DD lo(deepest_t), hi = log(j.s0);
        if (!(g(lo).hi < 0.0)) return std::unexpected("pin_level: the position lies deeper in the p bulk than e^-700");
        const BrentResult r = brent_zero(g, lo, hi, 1e-28);
        if (!error.empty()) return std::unexpected(error);
        if (r.status != BrentStatus::ok) return std::unexpected("pin_level: Brent did not converge on the p side");
        return exp(r.root);
    }
    if (x_cm < d) {
        const auto g = [&](const DD& s) {
            return check(integral(qag15([&](const DD& z) { return j.VT / field_of(j, phi_i(j, j.s0, j.phi0, z)); },
                                        j.s0, s, 0.0, quadrature_rtol))) -
                   DD(x_cm);
        };
        const BrentResult r = brent_zero(g, j.s0, j.sd, 4.0 * u_dd * j.sd.hi);
        if (!error.empty()) return std::unexpected(error);
        if (r.status != BrentStatus::ok) return std::unexpected("pin_level: Brent did not converge in the i region");
        return r.root;
    }
    const auto g = [&](const DD& t) { return check(n_position(j, t)) - DD(x_cm); };
    const DD lo(deepest_t), hi = log(j.barrier - j.sd);
    if (!(g(lo).hi > 0.0)) return std::unexpected("pin_level: the position lies deeper in the n bulk than e^-700");
    const BrentResult r = brent_zero(g, lo, hi, 1e-28);
    if (!error.empty()) return std::unexpected(error);
    if (r.status != BrentStatus::ok) return std::unexpected("pin_level: Brent did not converge on the n side");
    return j.barrier - exp(r.root);
}

std::expected<PinCollocation, std::string> pin_collocation(const PinJunction& j, double x_left_cm, double x_right_cm,
                                                           int degree, double piece_debye,
                                                           std::span<const double> query_cm) {
    const double d = j.in.i_length_cm;
    if (!(x_left_cm < 0.0 && x_right_cm > d) || degree < 2 || !(piece_debye > 0.0)) {
        return std::unexpected("pin_collocation: needs x_left < 0 < d < x_right, degree >= 2, piece_debye > 0");
    }
    const double NA = j.in.acceptors_cm3, ND = j.in.donors_cm3, Ci = j.in.intrinsic_net_cm3;
    const double Nref = std::max(NA, ND);
    const double VT = j.in.thermal_voltage_V, q = j.in.q_C, eps = j.in.eps_F_per_cm;
    const double lambda = std::sqrt(eps * VT / (q * Nref));
    const double B = j.barrier.to_double();
    const double nn = j.n_n.to_double() / Nref, pp = j.p_p.to_double() / Nref;
    // Pieces in increasing x (physical ends, cm), each with its net doping / N_ref.
    struct Piece {
        double a, b, C;
    };
    std::vector<Piece> pieces;
    std::size_t edge_0 = 0, edge_d = 0;  // the pieces ending at x = 0 and at x = d
    const auto add = [&](double a, double b, double C) {
        const int k = std::max(1, static_cast<int>(std::ceil((b - a) / (piece_debye * lambda))));
        for (int i = 0; i < k; ++i) {
            const double pa = i == 0 ? a : a + (b - a) * i / k;
            const double pb = i + 1 == k ? b : a + (b - a) * (i + 1) / k;
            pieces.push_back({pa, pb, C});
        }
    };
    add(x_left_cm, 0.0, -NA / Nref);
    edge_0 = pieces.size() - 1;
    add(0.0, d, Ci / Nref);
    edge_d = pieces.size() - 1;
    add(d, x_right_cm, ND / Nref);
    const std::size_t P = pieces.size(), per = static_cast<std::size_t>(degree) + 1, n = P * per;
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
    std::vector<double> half(P), mid(P);  // in units of lambda
    for (std::size_t p = 0; p < P; ++p) {
        half[p] = 0.5 * (pieces[p].b - pieces[p].a) / lambda;
        mid[p] = 0.5 * (pieces[p].b + pieces[p].a) / lambda;
    }
    const auto idx = [&](std::size_t p, std::size_t m) { return p * per + m; };
    // Initial guess: the depletion approximation. N_A x_p = N_D x_n + C_i d, and the drops add to V_T B:
    // (q / eps)(N_A x_p^2 / 2 + N_A x_p d - C_i d^2 / 2 + (N_A x_p - C_i d)^2 / (2 N_D)) = V_T B.
    const double k = q / (eps * VT);
    const double qa = 0.5 * NA + NA * NA / (2.0 * ND);
    const double qb = NA * d - NA * Ci * d / ND;
    const double qc = -0.5 * Ci * d * d + Ci * Ci * d * d / (2.0 * ND) - B / k;
    const double xp = std::max((-qb + std::sqrt(std::max(qb * qb - 4.0 * qa * qc, 0.0))) / (2.0 * qa), 0.0);
    const double xn = std::max((NA * xp - Ci * d) / ND, 0.0);
    const double y0 = 0.5 * k * NA * xp * xp;
    const auto guess = [&](double x) {
        double y;
        if (x <= -xp) y = 0.0;
        else if (x <= 0.0) y = 0.5 * k * NA * (x + xp) * (x + xp);
        else if (x <= d) y = y0 + k * (NA * xp * x - 0.5 * Ci * x * x);
        else if (x < d + xn) y = B - 0.5 * k * ND * (d + xn - x) * (d + xn - x);
        else y = B;
        return std::clamp(y, 0.0, B);
    };
    std::vector<double> y(n);
    for (std::size_t p = 0; p < P; ++p) {
        for (std::size_t m = 0; m < per; ++m) y[idx(p, m)] = guess(lambda * (mid[p] + half[p] * xi[m]));
    }
    PinCollocation out;
    const double scale = std::max(B, 1.0);
    for (int it = 1; it <= 200; ++it) {
        std::vector<double> A(n * n, 0.0), r(n, 0.0);
        for (std::size_t p = 0; p < P; ++p) {
            const double s2 = 1.0 / (half[p] * half[p]);
            for (std::size_t m = 1; m + 1 < per; ++m) {  // interior rows: y'' - (n - p - C) = 0
                const std::size_t row = idx(p, m);
                double d2 = 0.0;
                for (std::size_t b = 0; b < per; ++b) {
                    d2 += D2[m * per + b] * y[idx(p, b)];
                    A[row * n + idx(p, b)] = s2 * D2[m * per + b];
                }
                const double en = nn * std::exp(y[row] - B), ep = pp * std::exp(-y[row]);
                r[row] = s2 * d2 - (en - ep - pieces[p].C);
                A[row * n + row] -= en + ep;
            }
            if (p == 0) {  // y = 0 at the left end (node N is xi = -1)
                const std::size_t row = idx(p, per - 1);
                A[row * n + row] = 1.0;
                r[row] = y[row];
            }
            if (p + 1 == P) {  // y = B at the right end
                const std::size_t row = idx(p, 0);
                A[row * n + row] = 1.0;
                r[row] = y[row] - B;
            } else {  // continuity in row (p, 0), the derivative in row (p + 1, N)
                const std::size_t nq = p + 1;
                const std::size_t rc = idx(p, 0), rf = idx(nq, per - 1);
                A[rc * n + idx(p, 0)] = 1.0;
                A[rc * n + idx(nq, per - 1)] = -1.0;
                r[rc] = y[idx(p, 0)] - y[idx(nq, per - 1)];
                const double sp = 1.0 / half[p], sq = 1.0 / half[nq];
                double dl = 0.0, dr = 0.0;
                for (std::size_t b = 0; b < per; ++b) {
                    dl += D[b] * y[idx(p, b)];
                    dr += D[(per - 1) * per + b] * y[idx(nq, b)];
                    A[rf * n + idx(p, b)] += sp * D[b];
                    A[rf * n + idx(nq, b)] -= sq * D[(per - 1) * per + b];
                }
                r[rf] = sp * dl - sq * dr;
            }
        }
        if (!dense_solve(A, r, n)) return std::unexpected("pin_collocation: singular Newton matrix");
        double worst = 0.0;
        for (const double v : r) worst = std::max(worst, std::abs(v));
        const double damp = worst > 1.0 ? 1.0 / worst : 1.0;  // at most one V_T per step
        for (std::size_t c = 0; c < n; ++c) y[c] -= damp * r[c];
        out.newton_iterations = it;
        out.last_update = worst * damp;
        if (worst <= 1e-13 * scale) break;
        if (it == 200) return std::unexpected("pin_collocation: Newton did not converge in 200 iterations");
    }
    std::vector<DDVector> coef(P), dcoef(P);
    for (std::size_t p = 0; p < P; ++p) {
        DDVector vals(per), dvals(per);
        for (std::size_t m = 0; m < per; ++m) {
            vals[m] = DD(y[idx(p, m)]);
            double dv = 0.0;
            for (std::size_t b = 0; b < per; ++b) dv += D[m * per + b] * y[idx(p, b)];
            dvals[m] = DD(dv / half[p]);  // dy / d(x / lambda)
        }
        coef[p] = chebyshev_coefficients(vals);
        dcoef[p] = chebyshev_coefficients(dvals);
        out.tail_estimate = std::max(out.tail_estimate, std::abs(coef[p][per - 1].hi) + std::abs(coef[p][per - 2].hi));
    }
    const double to_field = VT / lambda;
    out.s0 = y[idx(edge_0, 0)];
    out.sd = y[idx(edge_d, 0)];
    out.field_p_edge_V_per_cm = to_field * chebyshev_evaluate(dcoef[edge_0], DD(1.0)).to_double();
    out.field_n_edge_V_per_cm = to_field * chebyshev_evaluate(dcoef[edge_d + 1], DD(-1.0)).to_double();
    const auto locate = [&](double x) {
        std::size_t p = 0;
        while (p + 1 < P && x > pieces[p].b) ++p;
        return p;
    };
    {
        const double xm = 0.5 * d;
        const std::size_t p = locate(xm);
        const double z = std::clamp((xm / lambda - mid[p]) / half[p], -1.0, 1.0);
        out.field_mid_V_per_cm = to_field * chebyshev_evaluate(dcoef[p], DD(z)).to_double();
    }
    for (const double x : query_cm) {
        const std::size_t p = locate(x);
        const double z = std::clamp((x / lambda - mid[p]) / half[p], -1.0, 1.0);
        out.levels.push_back(chebyshev_evaluate(coef[p], DD(z)).to_double());
    }
    return out;
}

}  // namespace NiTCAD::verify::reference
