#include "NiTCAD/verify/reference/first_integral.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "NiTCAD/verify/reference/tanh_sinh.hpp"

namespace NiTCAD::verify::reference {

// The Taylor series sum_{k >= 2} t^k / k! for |t| < 1/4 (positive terms for t > 0, alternating and decreasing for
// t < 0, so no cancellation beyond a factor 1 + O(t)); expm1(t) - t elsewhere, where the difference loses at most a
// factor of about 8.
DD poisson_boltzmann_G(const DD& t) {
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

NeutralBulk neutral_bulk(double ni_cm3, double net_doping_cm3) {
    // e^{|u_b|} = a + sqrt(a^2 + 1), a = |C| / (2 n_i): the majority density is n_i e^{|u_b|}, the minority n_i
    // e^{-|u_b|}, so n_b - p_b = 2 n_i sinh(u_b) = C and n_b p_b = n_i^2 hold to double-double.
    const DD ni(ni_cm3), C(std::abs(net_doping_cm3));
    const DD a = C / (ni * DD(2.0));
    const DD e = a + sqrt(a * a + DD(1.0));
    const DD level = log(e);
    NeutralBulk b;
    b.doping = C;
    if (net_doping_cm3 > 0.0) {
        b.u_b = level;
        b.n_b = ni * e;
        b.p_b = ni / e;
    } else {
        b.u_b = -level;
        b.n_b = ni / e;
        b.p_b = ni * e;
    }
    return b;
}

DD first_integral_phi(const NeutralBulk& b, const DD& u) {
    const DD t = u - b.u_b;
    return b.n_b * poisson_boltzmann_G(t) + b.p_b * poisson_boltzmann_G(-t);
}

std::expected<std::vector<ProfilePoint>, std::string> bulk_profile(const NeutralBulk& b, const DD& u_surface,
                                                                   const DD& debye_cm,
                                                                   std::span<const double> distance_cm) {
    std::vector<ProfilePoint> out(distance_cm.size());
    for (const double d : distance_cm) {
        if (!(d >= 0.0)) return std::unexpected("bulk_profile: a distance is negative or not a number");
    }
    const DD t_s = abs(u_surface - b.u_b);
    const double dir = u_surface > b.u_b ? 1.0 : -1.0;  // u = u_b + dir t
    const auto u_at = [&](const DD& t) { return dir > 0 ? b.u_b + t : b.u_b - t; };
    if (t_s.hi == 0.0) {
        for (auto& p : out) p = {b.u_b, DD(0.0)};
        return out;
    }
    // k(s) at t = e^s: Phi / |C| taken from t directly (u - u_b = dir t).
    const auto k = [&](const DD& s) {
        const DD t = exp(s);
        const DD d = dir > 0 ? t : -t;
        const DD phi = (b.n_b * poisson_boltzmann_G(d) + b.p_b * poisson_boltzmann_G(-d)) / b.doping;
        return t / sqrt(ldexp(phi, 1));
    };
    std::vector<std::size_t> order(distance_cm.size());
    for (std::size_t i = 0; i < order.size(); ++i) order[i] = i;
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t c) { return distance_cm[a] < distance_cm[c]; });
    DD s_prev = log(t_s), X_prev(0.0);
    for (const std::size_t i : order) {
        const DD target = DD(distance_cm[i]) / debye_cm;
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
            if (std::abs(F.hi) <= 1e-22 * std::max(1.0, std::abs(target.hi))) break;  // far below any criterion (1e-7)
            if (have_lo && !(lo < hi)) {  // the bracket has closed above the tolerance: X(s) is noisier than it
                return std::unexpected("bulk_profile: the Newton bracket closed with |X - d| = " +
                                       std::to_string(std::abs(F.hi)) + " at d = " + std::to_string(distance_cm[i]));
            }
            DD next = s + F / k(s);  // dX/ds = -k
            const bool inside = next < hi && (!have_lo || lo < next);
            if (!inside) next = have_lo ? ldexp(lo + hi, -1) : hi - DD(1.0);
            // X at next, integrated from the last accepted point s_prev.
            const bool down = next < s_prev;
            const TanhSinhResult q = down ? tanh_sinh(k, next, s_prev, 1e-28, 1e-26) : tanh_sinh(k, s_prev, next, 1e-28, 1e-26);
            if (q.status != TanhSinhStatus::ok) {
                return std::unexpected("bulk_profile: tanh-sinh did not converge at d = " + std::to_string(distance_cm[i]));
            }
            X = down ? X_prev + q.value : X_prev - q.value;
            s = next;
        }
        if (it == 200) return std::unexpected("bulk_profile: Newton did not converge at d = " + std::to_string(distance_cm[i]));
        s_prev = s;
        X_prev = X;
        const DD tail = exp(s);
        out[i] = {u_at(tail), tail};
    }
    return out;
}

}  // namespace NiTCAD::verify::reference
