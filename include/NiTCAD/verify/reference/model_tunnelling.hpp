// References for TVB-31.23-31.25: the WKB action of a piecewise-linear barrier, the Tsu-Esaki supply and current
// (Unit 55, TVR-3 step 8; DECISIONS 15.3, 15.9, 15.18). Gehring, "Simulation of Tunneling in Semiconductor Devices",
// TU Wien 2003, sec. 3.2, 3.3.1 and 3.5.1 (PA, read in step 1). Energies in eV, lengths in cm, masses in m0.
//   - One segment, W linear from W0 to W1 over L, a = W0 - E, b = W1 - E: int_0^L sqrt((a + (b - a) x / L)+) dx =
//     L g(a, b), g = (2/3) (A^3 - B^3) / (a - b) with A = sqrt(a+), B = sqrt(b+). Both ends inside: g = (2/3) (a + A B
//     + b) / (A + B) (A^3 - B^3 = (A - B)(a + AB + b), a - b = (A - B)(A + B)), with the limit sqrt a at a = b; one end
//     outside (o <= 0 < c): g = (2/3) C^3 / (c - o). The partials differentiate these through A = sqrt a (dA/da = 1 /
//     (2A)): dg/da = (A + 2B) / (3 (A + B)^2) inside; dg/dc = C (c - 3 o) / (3 (c - o)^2), dg/do = g / (c - o)
//     across a turning point. Second method (L0): tanh-sinh quadrature (TC4) of the root over [0, 1], split at the
//     turning point; the partials against central differences.
//   - The exponent S(E) = sum_k sqrt(2 m0 q) / hbar sqrt(m_k) L_k g_k (the 1e-2 for cm), the transmission e^(-2S).
//   - The supply N(E) = kT [ln(1 + e^((F1 - E)/kT)) - ln(1 + e^((F2 - E)/kT))], its partials the logistic
//     sigma(x) = 1 / (1 + e^-x): dN/dF1 = sigma(x1), dN/dF2 = -sigma(x2), dN/dE = sigma(x2) - sigma(x1).
//   - The current J = (A* / k^2) int TC N dE over [E_low, E_high] cut as NiTCAD documents (the barrier top, the
//     highest level + 40 kT): adaptive Gauss-Kronrod (TC3) between the kinks (every barrier end inside the range).
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <span>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/verify/reference/double_double.hpp"
#include "NiTCAD/verify/reference/gauss_kronrod.hpp"
#include "NiTCAD/verify/reference/special_functions.hpp"

namespace NiTCAD::verify::reference::model_tunnelling {

// sqrt(2 m0 q) / hbar [1 / (cm eV^(1/2))].
[[nodiscard]] inline DD action_scale() {
    return sqrt(DD(2.0) * DD(base::m0_kg) * DD(base::q_C)) / DD(base::hbar_J_s) / DD(100.0);
}

// {L g, L dg/da, L dg/db}; zero when a, b <= 0.
[[nodiscard]] inline std::array<DD, 3> segment_action(const DD& a, const DD& b, const DD& L) {
    const bool ai = a.hi > 0.0, bi = b.hi > 0.0;
    if (!ai && !bi) return {DD(0.0), DD(0.0), DD(0.0)};
    if (ai && bi) {
        const DD A = sqrt(a), B = sqrt(b), s = A + B;
        return {L * DD(2.0) * (a + A * B + b) / (DD(3.0) * s), L * (A + DD(2.0) * B) / (DD(3.0) * s * s),
                L * (DD(2.0) * A + B) / (DD(3.0) * s * s)};
    }
    const DD c = ai ? a : b, o = ai ? b : a;
    const DD C = sqrt(c), w = c - o;
    const DD g = DD(2.0) * c * C / (DD(3.0) * w);
    const DD d_in = L * C * (c - DD(3.0) * o) / (DD(3.0) * w * w), d_out = L * g / w;
    return ai ? std::array<DD, 3>{L * g, d_in, d_out} : std::array<DD, 3>{L * g, d_out, d_in};
}

// One segment of a barrier.
struct Segment {
    DD W0, W1, L, m;
};

// {S, dS/dE}.
[[nodiscard]] inline std::array<DD, 2> exponent(std::span<const Segment> barrier, const DD& E) {
    const DD scale = action_scale();
    DD S(0.0), dS(0.0);
    for (const Segment& s : barrier) {
        const auto g = segment_action(s.W0 - E, s.W1 - E, s.L);
        const DD c = scale * sqrt(s.m);
        S = S + c * g[0];
        dS = dS - c * (g[1] + g[2]);
    }
    return {S, dS};
}

[[nodiscard]] inline DD transmission(std::span<const Segment> barrier, const DD& E) {
    return exp(DD(-2.0) * exponent(barrier, E)[0]);
}

[[nodiscard]] inline DD logistic(const DD& x) {
    return x.hi >= 0.0 ? DD(1.0) / (DD(1.0) + exp(-x)) : exp(x) / (DD(1.0) + exp(x));
}

// x = {E, F1, F2, kT}: {N, dN/dF1, dN/dF2, dN/dE}.
[[nodiscard]] inline std::array<DD, 4> supply(const std::array<DD, 4>& x) {
    const DD x1 = (x[1] - x[0]) / x[3], x2 = (x[2] - x[0]) / x[3];
    const DD s1 = logistic(x1), s2 = logistic(x2);
    return {x[3] * (log1p_exp(x1) - log1p_exp(x2)), s1, -s2, s2 - s1};
}

// J [A/cm^2] over [E_low, E_high] with NiTCAD's documented cuts; richardson [A/(cm^2 K^2)]. ok false if the
// quadrature did not converge.
struct Current {
    DD J;
    bool ok = true;
};
[[nodiscard]] inline Current current(std::span<const Segment> barrier, double E_low, double E_high, double F1,
                                     double F2, double T, double richardson) {
    const DD k(base::k_B_eV_per_K);
    const DD kT = k * DD(T);
    double top = -INFINITY;
    for (const Segment& s : barrier) top = std::max({top, s.W0.hi, s.W1.hi});
    const double level = std::max({F1, F2, E_low});
    const double Eh = std::min({E_high, top, (DD(level) + DD(40.0) * kT).hi});
    Current out{DD(0.0), true};
    if (!(Eh > E_low)) return out;
    std::vector<double> cuts{E_low};
    for (const Segment& s : barrier) {
        for (const double w : {s.W0.hi, s.W1.hi}) {
            if (w > E_low && w < Eh) cuts.push_back(w);
        }
    }
    cuts.push_back(Eh);
    std::sort(cuts.begin(), cuts.end());
    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
    auto f = [&](const DD& E) {
        return transmission(barrier, E) * supply(std::array<DD, 4>{E, DD(F1), DD(F2), kT})[0];
    };
    DD sum(0.0);
    for (std::size_t i = 0; i + 1 < cuts.size(); ++i) {
        const auto q = qag15(f, DD(cuts[i]), DD(cuts[i + 1]), 0.0, 1e-24);
        if (q.status != QuadratureStatus::ok) out.ok = false;
        sum = sum + q.value;
    }
    out.J = DD(richardson) / (k * k) * sum;
    return out;
}

}  // namespace NiTCAD::verify::reference::model_tunnelling
