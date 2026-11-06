#include "NiTCAD/verify/reference/thermal_slab.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "NiTCAD/verify/reference/brent.hpp"
#include "NiTCAD/verify/reference/dopri5.hpp"
#include "NiTCAD/verify/reference/gauss_kronrod.hpp"

namespace NiTCAD::verify::reference {

namespace {

bool positive_finite(double v) { return std::isfinite(v) && v > 0.0; }

bool valid(const Conductivity& c) {
    return positive_finite(c.kappa300_W_cmK) && std::isfinite(c.exponent) && c.exponent != 1.0;
}

DD kappa(const Conductivity& c, const DD& T) { return DD(c.kappa300_W_cmK) * pow(T / DD(300.0), DD(-c.exponent)); }

}  // namespace

DD kirchhoff(const Conductivity& c, const DD& T) {
    const DD s = DD(1.0) - DD(c.exponent);
    return DD(c.kappa300_W_cmK) * DD(300.0) / s * (pow(T / DD(300.0), s) - DD(1.0));
}

DD kirchhoff_inverse(const Conductivity& c, const DD& u) {
    const DD s = DD(1.0) - DD(c.exponent);
    return DD(300.0) * pow(DD(1.0) + s * u / (DD(300.0) * DD(c.kappa300_W_cmK)), DD(1.0) / s);
}

std::expected<std::vector<DD>, std::string> slab_temperature(const Conductivity& c, double T0, double T1,
                                                             const std::vector<double>& s) {
    if (!valid(c) || !positive_finite(T0) || !positive_finite(T1)) {
        return std::unexpected("slab_temperature: invalid conductivity or end temperature");
    }
    const DD u0 = kirchhoff(c, DD(T0)), u1 = kirchhoff(c, DD(T1));
    std::vector<DD> out;
    for (const double f : s) {
        if (!(f >= 0.0 && f <= 1.0)) return std::unexpected("slab_temperature: a fraction outside [0, 1]");
        out.push_back(kirchhoff_inverse(c, u0 + DD(f) * (u1 - u0)));
    }
    return out;
}

std::expected<ResistorProfile, std::string> self_heated_resistor(const ResistorInput& in,
                                                                 const std::vector<double>& x) {
    if (!positive_finite(in.q_C) || !positive_finite(in.k_eV_per_K) || !positive_finite(in.length_cm) ||
        !positive_finite(in.density_cm3) || !std::isfinite(in.voltage_V) || in.voltage_V == 0.0 ||
        !positive_finite(in.ambient_K) || !valid(in.kappa) || !(in.mu_min >= 0.0) ||
        !positive_finite(in.mu_max) || !positive_finite(in.N_ref) || !std::isfinite(in.alpha) ||
        !std::isfinite(in.gamma) || in.steps < 4) {
        return std::unexpected("self_heated_resistor: an input is not finite, or out of range");
    }
    for (const double v : x) {
        if (!(v >= 0.0 && v <= in.length_cm)) return std::unexpected("self_heated_resistor: x outside [0, L]");
    }
    const DD q(in.q_C), k(in.k_eV_per_K), L(in.length_cm), N(in.density_cm3), V(in.voltage_V);
    const DD ratio = pow(N / DD(in.N_ref), DD(in.alpha));
    const auto mu = [&](const DD& T) {
        return DD(in.mu_min) + (DD(in.mu_max) * pow(T / DD(300.0), DD(in.gamma)) - DD(in.mu_min)) / (DD(1.0) + ratio);
    };
    const DD u_amb = kirchhoff(in.kappa, DD(in.ambient_K));

    // One shot at current J and u'(0) = s; the dense segments when wanted.
    struct Shot {
        DD uL, phiL;
    };
    const auto shoot = [&](const DD& J, const DD& s, long long M, std::vector<Dopri5Dense>* dense) {
        const auto f = [&](const DD&, const DDVector& y, DDVector& dy) {
            const DD T = kirchhoff_inverse(in.kappa, y[0]);
            const DD m = mu(T);
            const DD Tp = y[1] / kappa(in.kappa, T);
            dy[0] = y[1];
            dy[1] = -(J * J) / (q * m * N) - DD(1.5) * k * J * Tp;
            dy[2] = DD(1.0) / m;
        };
        const auto out = [&](const Dopri5Dense& d) {
            if (dense) dense->push_back(d);
        };
        const DDVector y = dopri5_fixed(f, DD(0.0), DDVector{u_amb, s, DD(0.0)}, L, M, out);
        return Shot{y[0], y[2]};
    };
    // u'(0) for current J: Brent on u(L) = u_amb, bracketed by 0 and 8 Q L (Q the ambient Joule density).
    bool failed = false;
    const auto slope = [&](const DD& J, long long M) {
        const DD Q = J * J / (q * mu(DD(in.ambient_K)) * N);
        const auto miss = [&](const DD& s) { return (shoot(J, s, M, nullptr).uL - u_amb) / (Q * L * L); };
        const BrentResult b = brent_zero(miss, DD(0.0), DD(8.0) * Q * L, 1e-30 * (Q * L).hi);
        if (b.status != BrentStatus::ok) failed = true;
        return b.root;
    };
    // g(J) = J + q N V / int dx / mu(T_J) = 0 by the secant method from the isothermal J and its fixed-point step.
    const auto solve = [&](long long M, int& iterations) {
        const auto g = [&](const DD& J, DD& s) {
            s = slope(J, M);
            return J + q * N * V / shoot(J, s, M, nullptr).phiL;
        };
        DD s0(0.0), s1(0.0);
        DD J0 = -q * N * V * mu(DD(in.ambient_K)) / L;
        DD g0 = g(J0, s0);
        DD J1 = J0 - g0;
        DD g1 = g(J1, s1);
        for (iterations = 1; iterations <= 100; ++iterations) {
            if (g1.hi == 0.0 || std::abs(((J1 - J0) / J1).to_double()) <= 1e-29) break;
            const DD J2 = J1 - g1 * (J1 - J0) / (g1 - g0);
            J0 = J1;
            g0 = g1;
            s0 = s1;
            J1 = J2;
            g1 = g(J1, s1);
        }
        return std::pair{J1, s1};
    };
    ResistorProfile out;
    const auto [J, s] = solve(in.steps, out.iterations);
    int it2 = 0;
    const auto [J2, s2] = solve(in.steps / 2, it2);
    if (failed || out.iterations > 100 || it2 > 100) {
        return std::unexpected("self_heated_resistor: the shooting or the current iteration did not converge");
    }
    out.current_A_cm2 = J;
    std::vector<Dopri5Dense> dense, dense2;
    (void)shoot(J, s, in.steps, &dense);
    (void)shoot(J2, s2, in.steps / 2, &dense2);
    const auto at = [&](const std::vector<Dopri5Dense>& d, double xv) {
        const DD h = L / DD(static_cast<double>(d.size()));
        const auto i = std::min(d.size() - 1, static_cast<std::size_t>(std::max(0.0, std::floor((DD(xv) / h).hi))));
        return kirchhoff_inverse(in.kappa, d[i].eval(DD(xv))[0]);
    };
    for (const double xv : x) {
        out.temperature.push_back(at(dense, xv));
        out.steps_difference = std::max(out.steps_difference, std::abs((at(dense, xv) - at(dense2, xv)).to_double()));
    }

    if (!in.run_second_method) return out;
    // The second method: Kirchhoff-form finite differences in double, 2^14 and 2^15 cells, Richardson.
    const double qd = in.q_C, kd = in.k_eV_per_K, Ld = in.length_cm, Nd = in.density_cm3, Vd = in.voltage_V;
    const double rd = std::pow(Nd / in.N_ref, in.alpha);
    const auto mud = [&](double T) { return in.mu_min + (in.mu_max * std::pow(T / 300.0, in.gamma) - in.mu_min) / (1.0 + rd); };
    const auto mud_T = [&](double T) { return in.gamma * in.mu_max * std::pow(T / 300.0, in.gamma) / T / (1.0 + rd); };
    const double sk = 1.0 - in.kappa.exponent;
    const auto ud = [&](double T) { return in.kappa.kappa300_W_cmK * 300.0 / sk * (std::pow(T / 300.0, sk) - 1.0); };
    const auto kd_T = [&](double T) { return in.kappa.kappa300_W_cmK * std::pow(T / 300.0, -in.kappa.exponent); };
    const auto fd = [&](int M) {
        const double h = Ld / M;
        std::vector<double> T(static_cast<std::size_t>(M) + 1, in.ambient_K);
        double Jd = -qd * Nd * Vd * mud(in.ambient_K) / Ld;
        for (int outer = 0; outer < 60; ++outer) {
            for (int it = 0; it < 50; ++it) {
                std::vector<double> a(M + 1), b(M + 1), c(M + 1), r(M + 1), cp(M + 1), dp(M + 1);
                const double cc = 1.5 * kd * Jd / (2.0 * h);
                for (int i = 1; i < M; ++i) {
                    const double Q = Jd * Jd / (qd * mud(T[i]) * Nd);
                    r[i] = (ud(T[i + 1]) - 2.0 * ud(T[i]) + ud(T[i - 1])) / (h * h) + Q + cc * (T[i + 1] - T[i - 1]);
                    a[i] = kd_T(T[i - 1]) / (h * h) - cc;
                    c[i] = kd_T(T[i + 1]) / (h * h) + cc;
                    b[i] = -2.0 * kd_T(T[i]) / (h * h) - Q * mud_T(T[i]) / mud(T[i]);
                }
                for (int i = 1; i < M; ++i) {
                    const double den = b[i] - (i > 1 ? a[i] * cp[i - 1] : 0.0);
                    cp[i] = c[i] / den;
                    dp[i] = (-r[i] - (i > 1 ? a[i] * dp[i - 1] : 0.0)) / den;
                }
                double change = 0.0;
                for (int i = M - 1; i >= 1; --i) {
                    const double d = dp[i] - (i < M - 1 ? cp[i] * dp[i + 1] : 0.0);
                    dp[i] = d;
                    T[i] += d;
                    change = std::max(change, std::abs(d));
                }
                if (change < 1e-12) break;
            }
            double phi = 0.0;
            for (int i = 0; i < M; ++i) phi += 0.5 * h * (1.0 / mud(T[i]) + 1.0 / mud(T[i + 1]));
            const double next = -qd * Nd * Vd / phi;
            const bool done = std::abs(next / Jd - 1.0) <= 1e-15;
            Jd = next;
            if (done) break;
        }
        return T;
    };
    const int M1 = 1 << 14;
    const std::vector<double> T1 = fd(M1), T2 = fd(2 * M1);
    for (int i = 0; i <= M1; ++i) {
        const double xr = Ld * i / M1;
        const double rich = (4.0 * T2[static_cast<std::size_t>(2 * i)] - T1[static_cast<std::size_t>(i)]) / 3.0;
        out.second_method = std::max(out.second_method, std::abs(rich - at(dense, std::min(xr, Ld)).to_double()));
    }
    return out;
}

std::expected<SeebeckResult, std::string> seebeck_rise(const SeebeckInput& in) {
    if (!positive_finite(in.k_eV_per_K) || !positive_finite(in.density_cm3) || !positive_finite(in.Nc300_cm3) ||
        !std::isfinite(in.thermopower_exponent) || !positive_finite(in.T0) || !positive_finite(in.T1)) {
        return std::unexpected("seebeck_rise: an input is not finite and positive");
    }
    const DD k(in.k_eV_per_K), T0(in.T0), T1(in.T1);
    const DD A = DD(2.5) + DD(in.thermopower_exponent) - log(DD(in.density_cm3) / DD(in.Nc300_cm3));
    const auto F = [&](const DD& T) { return T * log(T / DD(300.0)) - T; };
    SeebeckResult r;
    r.rise_eV = -k * (A * (T1 - T0) + DD(1.5) * (F(T1) - F(T0)));
    const auto integrand = [&](const DD& T) {
        return DD(2.5) + DD(in.thermopower_exponent) -
               log(DD(in.density_cm3) / (DD(in.Nc300_cm3) * pow(T / DD(300.0), DD(1.5))));
    };
    const AdaptiveResult I = qag15(integrand, T0, T1, 0.0, 1e-28);
    if (I.status != QuadratureStatus::ok && I.status != QuadratureStatus::roundoff) {
        return std::unexpected("seebeck_rise: the quadrature failed");
    }
    r.second_method = std::abs((-k * I.value - r.rise_eV).to_double());
    return r;
}

DD peltier_coefficient(double k_eV_per_K, double density_cm3, double Nc300_cm3, double thermopower_exponent,
                       double T) {
    const DD TT(T);
    return -DD(k_eV_per_K) * TT *
           (DD(2.5) + DD(thermopower_exponent) -
            log(DD(density_cm3) / (DD(Nc300_cm3) * pow(TT / DD(300.0), DD(1.5)))));
}

}  // namespace NiTCAD::verify::reference
