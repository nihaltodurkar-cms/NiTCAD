// Test helper: reference values computed here, in double-double arithmetic
// (reference_arithmetic.hpp), from the defining formulas, independently of the code under test.
// They replace literals that earlier units took from an outside tool.
//
// - The Bernoulli function x / expm1(x) and its derivative.
// - The normalized Fermi-Dirac integrals F_{1/2} and F_{-1/2} (F_j = int t^j / (1 + e^(t - eta))
//   dt / Gamma(j + 1)): for eta < -1 the alternating series sum (-1)^(k+1) e^(k eta) / k^(j+1)
//   (the polylogarithm -Li_{j+1}(-e^eta)); otherwise the integrals after t = u^2 (smooth
//   integrands), F_{1/2} = 4 / sqrt(pi) int u^2 f du and F_{-1/2} = 2 / sqrt(pi) int f du,
//   f = 1 / (1 + e^(u^2 - eta)), by composite 20-point Gauss-Legendre on [0, sqrt(eta+ + 120)]
//   with panels no wider than pi / (4 sqrt(eta)) near a large eta (the integrand's poles lie
//   pi / (2 sqrt(eta)) off the axis there).
// - Roots by bisection in DD arithmetic, 200 halvings.
// - The material formulas of the legacy (Varshni gap, (T/300)^1.5 densities of states, n_i, the
//   intrinsic-level depth, Slotboom narrowing, Caughey-Thomas mobility, the emission velocity),
//   evaluated in DD from this code's parameter sets.
#pragma once

#include <algorithm>
#include <cmath>
#include <functional>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "reference_arithmetic.hpp"

namespace reference {

// B(x) = x / (e^x - 1) and B'(x) = B (1/x - 1 - 1/(e^x - 1)); below 0.1 from the series of B in
// the Bernoulli numbers, B = sum B_n x^n / n!, which avoids the cancellation.
struct Bernoulli {
    DD value, derivative;
};
inline Bernoulli bernoulli(double xd) {
    const DD x = xd;
    if (std::abs(xd) < 0.1) {  // the tenth term is below 1e-35 there
        // B = 1 - x/2 + sum_{k>=1} B_2k x^2k / (2k)!.
        DD b = DD(1.0) - x / DD(2.0), db = DD(-0.5);
        DD power = x * x, factorial = 2.0;  // x^2k, (2k)!
        for (int k = 1; k <= 10; ++k) {
            const DD c = bernoulli_number(2 * k) / factorial;
            b += c * power;
            db += c * DD(2.0 * k) * power / x;
            power = power * x * x;
            factorial = factorial * DD(2.0 * k + 1.0) * DD(2.0 * k + 2.0);
        }
        return {b, db};
    }
    const DD em = expm1(x);
    if (em.hi == std::numeric_limits<double>::infinity()) return {0.0, 0.0};
    const DD B = x / em;
    return {B, B * (DD(1.0) / x - DD(1.0) - DD(1.0) / em)};
}

namespace detail {

inline DD fermi_series(double eta, double order) {  // sum (-1)^(k+1) e^(k eta) / k^order
    const DD z = exp(DD(eta));
    DD sum = 0.0, zk = z;
    for (int k = 1; k < 5000; ++k) {
        const DD term = zk / exp(DD(order) * log(DD(k)));
        sum = (k % 2 == 1) ? sum + term : sum - term;
        if (std::abs(term.hi) < 1e-34 * std::abs(sum.hi)) break;
        zk = zk * z;
        if (zk.hi == 0.0) break;
    }
    return sum;
}

inline const GaussLegendre& gl20() {
    static const GaussLegendre g = gauss_legendre(20);
    return g;
}

// int_0^U g(u) / (1 + e^(u^2 - eta)) du, g(u) = u^2 or 1.
inline DD fermi_quadrature(double eta, bool squared) {
    const double U = std::sqrt(std::max(eta, 0.0) + 120.0);
    const double width = std::min(0.25, pi().hi / (4.0 * std::sqrt(std::max(eta, 1.0))));
    const int panels = static_cast<int>(std::ceil(U / width));
    const GaussLegendre& g = gl20();
    DD sum = 0.0;
    for (int p = 0; p < panels; ++p) {
        const DD a = DD(U) * DD(p) / DD(panels), b = DD(U) * DD(p + 1.0) / DD(panels);
        const DD half = (b - a) / DD(2.0), mid = (a + b) / DD(2.0);
        for (std::size_t i = 0; i < g.x.size(); ++i) {
            for (const double sign : {-1.0, 1.0}) {
                const DD u = mid + DD(sign) * half * g.x[i];
                const DD y = u * u - DD(eta);
                if (y.hi > 300.0) continue;  // e^-300 of the integrand's scale
                const DD f = DD(1.0) / (DD(1.0) + exp(y));
                sum += g.w[i] * half * (squared ? u * u * f : f);
            }
        }
    }
    return sum;
}

}  // namespace detail

// F_{1/2}(eta) and F_{-1/2}(eta), normalized (F_j -> e^eta as eta -> -infinity).
inline DD fermi_half(double eta) {
    if (eta < -1.0) return detail::fermi_series(eta, 1.5);
    return DD(4.0) / sqrt(pi()) * detail::fermi_quadrature(eta, true);
}
inline DD fermi_minus_half(double eta) {
    if (eta < -1.0) return detail::fermi_series(eta, 0.5);
    return DD(2.0) / sqrt(pi()) * detail::fermi_quadrature(eta, false);
}
// Both methods at one point, for checking them against each other where both converge.
inline DD fermi_half_series(double eta) { return detail::fermi_series(eta, 1.5); }
inline DD fermi_half_quadrature(double eta) {
    return DD(4.0) / sqrt(pi()) * detail::fermi_quadrature(eta, true);
}

// The root of an increasing f on [lo, hi] by 200 halvings.
inline double bisect(const std::function<DD(double)>& f, double lo, double hi) {
    for (int i = 0; i < 200; ++i) {
        const double mid = 0.5 * (lo + hi);
        if (mid == lo || mid == hi) break;
        (f(mid).hi < 0.0 ? lo : hi) = mid;
    }
    return 0.5 * (lo + hi);
}

// Charge-neutral equilibria in the n_ie gauge (statistics.hpp): n = n_ie e^g F(eta - g) under
// Fermi-Dirac statistics, n_ie e^eta under Boltzmann; p likewise with -eta. With dopant levels the
// ionized densities are N / (1 + g_d e^(eta_band + E/kT)) (ionization.hpp). Root by bisection on
// [-200, 200], in DD.
struct Neutral {
    double eta;
    DD n, p;
};
struct Dopants {
    double donors = 0.0, acceptors = 0.0;
    double donor_kT = 0.0, acceptor_kT = 0.0;  // <= 0: complete ionization
    double donor_degeneracy = 2.0, acceptor_degeneracy = 4.0;
};
inline DD carrier_density(double ni, double g, double eta, bool fermi_dirac) {
    if (!fermi_dirac) return DD(ni) * exp(DD(eta));
    return DD(ni) * exp(DD(g)) * fermi_half(eta - g);
}
inline DD ionized(double N, double eta_band, double depth_kT, double degeneracy) {
    if (!(depth_kT > 0.0) || N == 0.0) return DD(N);
    return DD(N) / (DD(1.0) + DD(degeneracy) * exp(DD(eta_band) + DD(depth_kT)));
}
inline Neutral neutral_equilibrium(const Dopants& d, double ni, double gn, double gp,
                                   bool fermi_dirac) {
    const auto charge = [&](double eta) {
        return carrier_density(ni, gn, eta, fermi_dirac) -
               carrier_density(ni, gp, -eta, fermi_dirac) -
               ionized(d.donors, eta - gn, d.donor_kT, d.donor_degeneracy) +
               ionized(d.acceptors, -eta - gp, d.acceptor_kT, d.acceptor_degeneracy);
    };
    const double eta = bisect(charge, -200.0, 200.0);
    return {eta, carrier_density(ni, gn, eta, fermi_dirac),
            carrier_density(ni, gp, -eta, fermi_dirac)};
}

// Material formulas (legacy materials.py), in DD from this code's parameters.
inline DD thermal_voltage(double T) {
    return DD(NiTCAD::base::k_B_J_per_K) * DD(T) / DD(NiTCAD::base::q_C);
}
inline DD band_gap(const NiTCAD::physics::SemiconductorParameters& p, double T) {
    return DD(p.Eg0_eV) -
           DD(p.varshni_alpha_eV_per_K) * DD(T) * DD(T) / (DD(T) + DD(p.varshni_beta_K));
}
inline DD dos_factor(double T) {
    const DD r = DD(T) / DD(300.0);
    return r * sqrt(r);
}
inline DD intrinsic_density(const NiTCAD::physics::SemiconductorParameters& p, double T) {
    return sqrt(DD(p.Nc300) * DD(p.Nv300)) * dos_factor(T) *
           exp(-band_gap(p, T) / (DD(2.0) * thermal_voltage(T)));
}
inline DD intrinsic_level_depth(const NiTCAD::physics::SemiconductorParameters& p, double T) {
    return DD(p.electron_affinity_eV) + DD(0.5) * band_gap(p, T) +
           DD(0.5) * thermal_voltage(T) * log(DD(p.Nc300) / DD(p.Nv300));
}
inline DD slotboom_narrowing(const NiTCAD::physics::SemiconductorParameters& p, double N) {
    const DD x = log(DD(N) / DD(p.bandgap_narrowing.N0));
    return DD(p.bandgap_narrowing.E0_eV) * (x + sqrt(x * x + DD(0.5)) - sqrt(DD(0.5)));
}
inline DD caughey_thomas(const NiTCAD::physics::CaugheyThomasParameters& ct, double N, double T) {
    const DD mu_max = DD(ct.mu_max) * pow(DD(T) / DD(300.0), DD(ct.T_exponent));
    return DD(ct.mu_min) +
           (mu_max - DD(ct.mu_min)) / (DD(1.0) + pow(DD(N) / DD(ct.N_ref), DD(ct.alpha)));
}

}  // namespace reference
