// References for TVB-31.17-31.19: incomplete ionization, the Altermatt-Schenk density-dependent level and the
// charge-neutral equilibria (Unit 55, TVR-3 step 6; DECISIONS 15.3, 15.9, 15.16).
//   - Ionized fraction of a dopant with one level in equilibrium with its band (the standard form, Genius and Xyce
//     reproduce it; PD here from the level's grand-canonical occupancy): s = 1 / (1 + g e^(eta + E/kT)), with the
//     Altermatt-Schenk bound fraction b, N_ion = N ((1 - b) + b s) and d N_ion / d eta = -N b s (1 - s). A depth of 0
//     or less is the documented complete ionization (the dopant not modelled; ionization.hpp), N_ion = N.
//   - Altermatt-Schenk, SISPAD 2006 eqs (1), (5): E_dop(N) = E_dop,0 / (1 + (N/N_ref)^c), b(N) = 1 / (1 + (N/N_b)^d)
//     (PA, read in Unit 37).
//   - The neutral equilibrium: the root in eta of n(eta) - p(eta) - N_D+ + N_A- = 0 with n = n_ie e^eta (Boltzmann) or
//     n = n_ie e^g F_1/2(eta - g) (Fermi-Dirac, the n_ie gauge of statistics.hpp, F_1/2 normalised to e^x), and p
//     the mirror, by Brent (TC5) in double-double. F_1/2 is the TVB-02 reference (fermi_dirac.hpp): its series for
//     x <= -2, tanh-sinh above. The condition sums of n, p and eta at the root are by implicit
//     differentiation, d eta / d x_i = -(dG/dx_i) / (dG/d eta), the partials of G and of the densities at fixed eta
//     by double-double central differences, so the root is solved once per point.
#pragma once

#include <array>
#include <cmath>


#include "NiTCAD/verify/reference/brent.hpp"
#include "NiTCAD/verify/reference/double_double.hpp"
#include "NiTCAD/verify/reference/fermi_dirac.hpp"
#include "NiTCAD/verify/reference/model_band.hpp"
#include "NiTCAD/verify/reference/model_mobility.hpp"
#include "NiTCAD/verify/reference/model_recombination.hpp"

namespace NiTCAD::verify::reference::model_statistics {

// x = {N, eta, E/kT, g, b}: {N_ion, d N_ion / d eta}.
template <class S>
[[nodiscard]] std::array<S, 2> ionized(const std::array<S, 5>& x) {
    if (!(model_band::detail::re(x[2]) > 0.0)) return {x[0], S(0.0)};
    const std::array<S, 3> o = model_recombination::fermi_occupancy(x[1] + x[2] + log(x[3]));
    return {x[0] * ((S(1.0) - x[4]) + x[4] * o[0]), -(x[0] * x[4] * o[2])};
}

// x = {E_dop0, N_ref, c, N}: E_dop(N).
template <class S>
[[nodiscard]] S altermatt_level(const std::array<S, 4>& x) {
    if (model_band::detail::re(x[3]) == 0.0) return x[0];
    return x[0] / (S(1.0) + model_mobility::detail::spow(x[3] / x[1], x[2]));
}

// x = {N_b, d, N}: b(N).
template <class S>
[[nodiscard]] S altermatt_bound(const std::array<S, 3>& x) {
    if (model_band::detail::re(x[2]) == 0.0) return S(1.0);
    return S(1.0) / (S(1.0) + model_mobility::detail::spow(x[2] / x[0], x[1]));
}

// F_1/2(x), normalised to e^x as x -> -inf: the series for x <= -2, tanh-sinh (TC4) above (cross-checked against
// Gauss-Kronrod in TVB-02 and model_statistics_test); NaN if the quadrature reports failure.
[[nodiscard]] inline DD fermi_half(const DD& x) {
    if (x.hi <= -2.0) return fermi_series(FermiOrder::half, x).value;
    const FermiQuadrature q = fermi_tanh_sinh(FermiOrder::half, x);
    return q.ok ? q.value : DD(std::nan(""));
}

// The inputs of an equilibrium: {N_D, N_A, n_ie, g_n, g_p, E_D/kT, E_A/kT, g_D, g_A, b_D, b_A}.
using EquilibriumInputs = std::array<DD, 11>;

// n = n_ie e^eta or n_ie e^g_n F_1/2(eta - g_n); p the mirror.
[[nodiscard]] inline DD electron_density(const EquilibriumInputs& x, const DD& eta, bool fd) {
    return fd ? x[2] * exp(x[3]) * fermi_half(eta - x[3]) : x[2] * exp(eta);
}
[[nodiscard]] inline DD hole_density(const EquilibriumInputs& x, const DD& eta, bool fd) {
    return fd ? x[2] * exp(x[4]) * fermi_half(-eta - x[4]) : x[2] * exp(-eta);
}
[[nodiscard]] inline DD donors_ionized(const EquilibriumInputs& x, const DD& eta) {
    return ionized(std::array<DD, 5>{x[0], eta - x[3], x[5], x[7], x[9]})[0];
}
[[nodiscard]] inline DD acceptors_ionized(const EquilibriumInputs& x, const DD& eta) {
    return ionized(std::array<DD, 5>{x[1], -eta - x[4], x[6], x[8], x[10]})[0];
}

// G = (n - p - N_D+ + N_A-) / (N_D + N_A + n_ie), increasing in eta, from given densities.
[[nodiscard]] inline DD charge_with(const EquilibriumInputs& x, const DD& eta, const DD& n, const DD& p) {
    return (n - p - donors_ionized(x, eta) + acceptors_ionized(x, eta)) / (x[0] + x[1] + x[2]);
}
[[nodiscard]] inline DD charge(const EquilibriumInputs& x, const DD& eta, bool fd) {
    return charge_with(x, eta, electron_density(x, eta, fd), hole_density(x, eta, fd));
}

struct Equilibrium {
    DD n, p, eta;
    std::array<double, 3> S{};  // condition sums of n, p, eta over the 11 inputs
    bool ok = false;
};

// The root by Brent (TC5) to t = 1e-30 max(1, |eta|), from a bracket about `guess` whose ends are verified by the sign
// of G and widened (x4) until they bracket: the root is the reference's own balance whatever the guess.
[[nodiscard]] inline Equilibrium equilibrium(const EquilibriumInputs& x, bool fd, double guess) {
    Equilibrium out;
    auto G = [&](const DD& eta) { return charge(x, eta, fd); };
    double w = 1e-6 * std::fmax(1.0, std::abs(guess));
    DD a(guess - w), b(guess + w);
    for (int i = 0; i < 40 && !(G(a).hi < 0.0); ++i) a = DD(guess) - DD(w *= 4.0);
    w = 1e-6 * std::fmax(1.0, std::abs(guess));
    for (int i = 0; i < 40 && !(G(b).hi > 0.0); ++i) b = DD(guess) + DD(w *= 4.0);
    const BrentResult r = brent_zero(G, a, b, 1e-30 * std::fmax(1.0, std::abs(guess)));
    if (r.status != BrentStatus::ok) return out;
    const DD eta = r.root;
    out.eta = eta;
    out.n = electron_density(x, eta, fd);
    out.p = hole_density(x, eta, fd);
    // Implicit differentiation, d eta / d x_i = -(dG/dx_i) / (dG/d eta), by central differences at fixed eta. Only
    // n_ie, g_n and g_p move the densities: n and p are linear in n_ie, and g_n (g_p) moves n (p) alone.
    auto central = [](auto&& f, const DD& v, const DD& h) { return (f(v + h) - f(v - h)) / ldexp(h, 1); };
    const DD he = ldexp(DD(std::fmax(1.0, std::abs(eta.hi))), -27);
    const DD dn_deta = central([&](const DD& e) { return electron_density(x, e, fd); }, eta, he);
    const DD dp_deta = central([&](const DD& e) { return hole_density(x, e, fd); }, eta, he);
    const DD sum = x[0] + x[1] + x[2];
    const DD dG_deta =
        (dn_deta - dp_deta - central([&](const DD& e) { return donors_ionized(x, e); }, eta, he) +
         central([&](const DD& e) { return acceptors_ionized(x, e); }, eta, he)) /
        sum;
    double Sn = 0.0, Sp = 0.0, Se = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        if (x[i].hi == 0.0) continue;  // x_i d/dx_i vanishes
        if ((i == 7 || i == 9) && !(x[5].hi > 0.0)) continue;   // complete ionization: g_D, b_D do not enter
        if ((i == 8 || i == 10) && !(x[6].hi > 0.0)) continue;
        const DD h = ldexp(abs(x[i]), -27);
        auto at = [&](const DD& v) {
            EquilibriumInputs y = x;
            y[i] = v;
            return y;
        };
        DD dn(0.0), dp(0.0), dG;
        if (i == 2) {
            dn = out.n / x[2];
            dp = out.p / x[2];
            dG = central([&](const DD& v) { return charge_with(at(v), eta, out.n * (v / x[2]), out.p * (v / x[2])); },
                         x[2], h);
        } else if (i == 3) {
            dn = central([&](const DD& v) { return electron_density(at(v), eta, fd); }, x[3], h);
            dG = (dn - central([&](const DD& v) { return donors_ionized(at(v), eta); }, x[3], h)) / sum;
        } else if (i == 4) {
            dp = central([&](const DD& v) { return hole_density(at(v), eta, fd); }, x[4], h);
            dG = (-dp + central([&](const DD& v) { return acceptors_ionized(at(v), eta); }, x[4], h)) / sum;
        } else {
            dG = central([&](const DD& v) { return charge_with(at(v), eta, out.n, out.p); }, x[i], h);
        }
        const DD de = -dG / dG_deta;
        Se += std::abs((x[i] * de).hi);
        Sn += std::abs((x[i] * (dn + dn_deta * de)).hi);
        Sp += std::abs((x[i] * (dp + dp_deta * de)).hi);
    }
    out.S = {Sn, Sp, Se};
    out.ok = std::isfinite(out.n.hi) && std::isfinite(out.p.hi) && std::isfinite(Sn + Sp + Se);
    return out;
}

}  // namespace NiTCAD::verify::reference::model_statistics
