#include "NiTCAD/verify/reference/multiplication.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include "NiTCAD/verify/reference/brent.hpp"
#include "NiTCAD/verify/reference/chebyshev.hpp"
#include "NiTCAD/verify/reference/dopri5.hpp"
#include "NiTCAD/verify/reference/first_integral.hpp"
#include "NiTCAD/verify/reference/gauss_kronrod.hpp"
#include "NiTCAD/verify/reference/model_generation.hpp"

namespace NiTCAD::verify::reference {

namespace {

bool positive(double v) { return std::isfinite(v) && v > 0.0; }

bool valid(const IonizationBranches& b) {
    const bool none = b.A_low_per_cm == 0.0 && b.A_high_per_cm == 0.0;
    const bool finite = std::isfinite(b.A_low_per_cm) && std::isfinite(b.A_high_per_cm) && b.A_low_per_cm >= 0.0 &&
                        b.A_high_per_cm >= 0.0;
    return finite && (none || (positive(b.B_low_V_per_cm) && positive(b.B_high_V_per_cm) &&
                               std::isfinite(b.switch_V_per_cm) && b.switch_V_per_cm >= 0.0));
}

// p G(a) for a >= 0 with the density p = e^lp carried as a logarithm.
DD minority(const DD& lp, const DD& a) {
    if (lp.hi > -300.0) return exp(lp) * poisson_boltzmann_G(a);
    return exp(lp + a) - exp(lp) * (DD(1.0) + a);
}

// The junction at one bias: both sides from their own bulks.
struct Junction {
    DD VT, scale;            // scale = 2 q / (eps V_T): E = V_T sqrt(scale Phi)
    DD p_p, n_n, ln_np, ln_pn;
    DD barrier, s_j, r_j;    // reduced
    DD phi_p(const DD& s) const { return p_p * poisson_boltzmann_G(-s) + minority(ln_np, s); }
    DD phi_n(const DD& r) const { return n_n * poisson_boltzmann_G(-r) + minority(ln_pn, r); }
    DD field(bool p_side, const DD& c) const { return VT * sqrt(scale * (p_side ? phi_p(c) : phi_n(c))); }
};

std::expected<Junction, std::string> junction(const MultiplicationInput& in, double bias_V) {
    if (!positive(in.q_C) || !positive(in.eps_F_per_cm) || !positive(in.thermal_voltage_V) || !positive(in.ni_cm3) ||
        !positive(in.acceptors_cm3) || !positive(in.donors_cm3) || !valid(in.electron) || !valid(in.hole) ||
        !positive(in.gamma)) {
        return std::unexpected("multiplication: an input is not finite, or not positive where it must be");
    }
    if (!std::isfinite(bias_V) || bias_V > 0.0) return std::unexpected("multiplication: the bias must be finite and not forward");
    Junction j;
    j.VT = DD(in.thermal_voltage_V);
    j.scale = DD(2.0) * DD(in.q_C) / (DD(in.eps_F_per_cm) * j.VT);
    const DD ln_m = log(DD(in.ni_cm3)) + DD(bias_V) / ldexp(j.VT, 1);
    const DD m = exp(ln_m);  // 0 once it underflows: the majority densities are then the dopings exactly
    const DD half_a = ldexp(DD(in.acceptors_cm3), -1), half_d = ldexp(DD(in.donors_cm3), -1);
    j.p_p = half_a + sqrt(half_a * half_a + m * m);
    j.n_n = half_d + sqrt(half_d * half_d + m * m);
    j.ln_np = ldexp(ln_m, 1) - log(j.p_p);
    j.ln_pn = ldexp(ln_m, 1) - log(j.n_n);
    j.barrier = log(j.n_n) + log(j.p_p) - ldexp(ln_m, 1);
    if (!(j.barrier.hi > 0.0)) return std::unexpected("multiplication: no barrier at this bias");
    const BrentResult root = brent_zero([&](const DD& s) { return j.phi_p(s) - j.phi_n(j.barrier - s); }, DD(0.0),
                                        j.barrier, u_dd * j.barrier.hi);
    if (root.status != BrentStatus::ok) return std::unexpected("multiplication: Brent did not converge for the junction");
    j.s_j = root.root;
    j.r_j = j.barrier - j.s_j;
    return j;
}

DD alpha(const IonizationBranches& b, double gamma, const DD& E) {
    const bool high = E.hi >= b.switch_V_per_cm;
    const double A = high ? b.A_high_per_cm : b.A_low_per_cm;
    if (A == 0.0) return DD(0.0);
    const double B = high ? b.B_high_V_per_cm : b.B_low_V_per_cm;
    return model_generation::vodm<DD>({DD(A), DD(B), DD(gamma), E})[0];
}

// A piece of the integration path in a side's level coordinate, from a to b in the direction of +x (a < b on the
// p side, a > b on the n side); dx = sign V_T dc / E.
struct Piece {
    bool p_side;
    DD a, b;
};

// The pieces in order of increasing x: the p side from its cut level up to the junction, the n side from the
// junction down to its cut level, each split at the switch fields. Empty when the peak field is below the cut.
std::expected<std::vector<Piece>, std::string> pieces(const MultiplicationInput& in, const Junction& j) {
    std::vector<Piece> out;
    double b_min = INFINITY;
    for (const IonizationBranches* c : {&in.electron, &in.hole}) {
        if (c->A_low_per_cm > 0.0) b_min = std::min(b_min, c->B_low_V_per_cm);
        if (c->A_high_per_cm > 0.0) b_min = std::min(b_min, c->B_high_V_per_cm);
    }
    if (!std::isfinite(b_min)) return out;
    const double cut = in.gamma * b_min / 200.0;
    const DD peak = j.field(true, j.s_j);
    if (!(peak.hi > cut)) return out;
    // Marks at the cut field and at each doubling of it below the peak: exp(-gamma B / E) changes by e^100 over
    // the first octave and by a few e-folds over the last, so equal panels in the level resolve each octave alike.
    std::vector<double> marks;
    for (double f = cut; f < peak.hi; f *= 2.0) marks.push_back(f);
    for (const double s : {in.electron.switch_V_per_cm, in.hole.switch_V_per_cm}) {
        if (s > cut && s < peak.hi) marks.push_back(s);
    }
    std::sort(marks.begin(), marks.end());
    marks.erase(std::unique(marks.begin(), marks.end()), marks.end());
    for (const bool p_side : {true, false}) {
        const DD top = p_side ? j.s_j : j.r_j;
        std::vector<DD> levels;
        for (const double f : marks) {
            const BrentResult r = brent_zero([&](const DD& c) { return j.field(p_side, c) - DD(f); }, DD(0.0), top,
                                             u_dd * top.hi);
            if (r.status != BrentStatus::ok) return std::unexpected("multiplication: Brent did not converge for a field level");
            levels.push_back(r.root);
        }
        levels.push_back(top);
        if (p_side) {
            for (std::size_t k = 0; k + 1 < levels.size(); ++k) out.push_back({true, levels[k], levels[k + 1]});
        } else {
            for (std::size_t k = levels.size() - 1; k > 0; --k) out.push_back({false, levels[k], levels[k - 1]});
        }
    }
    return out;
}

struct Integrals {
    DD I_n, I_p, phi;
};

Integrals integrate(const MultiplicationInput& in, const Junction& j, const std::vector<Piece>& path, int panels) {
    Integrals s{DD(0.0), DD(0.0), DD(0.0)};
    for (const Piece& piece : path) {
        const DD sign = piece.p_side ? j.VT : -j.VT;  // dx = sign dc / E
        const auto rate = [&](const DD& c, bool electron, bool difference) {
            const DD E = j.field(piece.p_side, c);
            const DD an = alpha(in.electron, in.gamma, E), ap = alpha(in.hole, in.gamma, E);
            return sign * (difference ? an - ap : electron ? an : ap) / E;
        };
        const DD width = (piece.b - piece.a) / DD(static_cast<double>(panels));
        for (int k = 0; k < panels; ++k) {
            const DD a = piece.a + width * DD(static_cast<double>(k));
            const DD b = k + 1 == panels ? piece.b : a + width;
            const auto phi_from = [&](const DD& c) { return qk15([&](const DD& z) { return rate(z, true, true); }, a, c).value; };
            const DD start = s.phi;
            s.I_n += exp(-start) * qk15([&](const DD& c) { return rate(c, true, false) * exp(-phi_from(c)); }, a, b).value;
            s.I_p += exp(-start) * qk15([&](const DD& c) { return rate(c, false, false) * exp(-phi_from(c)); }, a, b).value;
            s.phi += phi_from(b);
        }
    }
    return s;
}

// Dormand-Prince on d(phi, I_n, I_p) / d level across the pieces, at relative tolerance rtol.
std::expected<Integrals, std::string> integrate_ode(const MultiplicationInput& in, const Junction& j,
                                                    const std::vector<Piece>& path, double rtol, long long* steps) {
    DDVector y{DD(0.0), DD(0.0), DD(0.0)};  // phi, I_n, I_p
    Dopri5Options o;
    o.rtol = rtol;
    o.atol = 1e-24;
    for (const Piece& piece : path) {
        const DD sign = piece.p_side ? j.VT : -j.VT;
        const auto f = [&](const DD& c, const DDVector& v, DDVector& d) {
            const DD E = j.field(piece.p_side, c);
            const DD an = alpha(in.electron, in.gamma, E), ap = alpha(in.hole, in.gamma, E);
            const DD w = sign * exp(-v[0]) / E;
            d[0] = sign * (an - ap) / E;
            d[1] = an * w;
            d[2] = ap * w;
        };
        const Dopri5Result r = dopri5(f, piece.a, y, piece.b, o);
        if (r.status != Dopri5Status::ok) return std::unexpected("multiplication: the integrator did not finish");
        y = r.y;
        if (steps != nullptr) *steps += r.naccpt;
    }
    return Integrals{y[1], y[2], y[0]};
}

}  // namespace

std::expected<Multiplication, std::string> multiplication(const MultiplicationInput& in, double bias_V) {
    const auto j = junction(in, bias_V);
    if (!j) return std::unexpected(j.error());
    const auto path = pieces(in, *j);
    if (!path) return std::unexpected(path.error());
    Multiplication m;
    const auto fine = integrate_ode(in, *j, *path, 1e-13, &m.steps);
    const auto coarse = integrate_ode(in, *j, *path, 1e-11, nullptr);
    if (!fine || !coarse) return std::unexpected(!fine ? fine.error() : coarse.error());
    m.barrier_V = j->VT * j->barrier;
    m.bending_p_V = j->VT * j->s_j;
    m.bending_n_V = j->VT * j->r_j;
    m.field_max_V_per_cm = j->field(true, j->s_j);
    const DD eps(in.eps_F_per_cm), q(in.q_C);
    m.depletion_p_cm = eps * m.field_max_V_per_cm / (q * DD(in.acceptors_cm3));
    m.depletion_n_cm = eps * m.field_max_V_per_cm / (q * DD(in.donors_cm3));
    m.ionization_integral = fine->I_n;
    m.hole_integral = fine->I_p;
    m.phi_total = fine->phi;
    m.inverse_n = DD(1.0) - fine->I_n;
    m.inverse_p = m.inverse_n * exp(fine->phi);
    m.below_breakdown = m.inverse_n.hi > 0.0;
    if (m.below_breakdown) {
        m.multiplication_n = DD(1.0) / m.inverse_n;
        m.multiplication_p = DD(1.0) / m.inverse_p;
    }
    m.integration_error = std::abs((fine->I_n - coarse->I_n).to_double());
    m.identity_residual = std::abs((m.inverse_n - (exp(-fine->phi) - fine->I_p)).to_double());
    return m;
}

std::expected<MultiplicationQuadrature, std::string> multiplication_quadrature(const MultiplicationInput& in,
                                                                               double bias_V, int panels) {
    if (panels < 1) return std::unexpected("multiplication_quadrature: panels must be at least 1");
    const auto j = junction(in, bias_V);
    if (!j) return std::unexpected(j.error());
    const auto path = pieces(in, *j);
    if (!path) return std::unexpected(path.error());
    const Integrals s = integrate(in, *j, *path, panels);
    return MultiplicationQuadrature{s.I_n, s.I_p, s.phi};
}

std::expected<double, std::string> bias_at_inverse_multiplication(const MultiplicationInput& in, bool electron,
                                                                 double inverse, double reverse_V, double forward_V) {
    if (!(reverse_V < forward_V) || !(forward_V <= 0.0) || !std::isfinite(reverse_V) || !std::isfinite(inverse)) {
        return std::unexpected("bias_at_inverse_multiplication: needs reverse_V < forward_V <= 0 and a finite target");
    }
    std::string error;
    const auto g = [&](const DD& V) {
        const auto m = multiplication(in, V.to_double());
        if (!m) {
            error = m.error();
            return DD(0.0);
        }
        return (electron ? m->inverse_n : m->inverse_p) - DD(inverse);
    };
    const DD lo = g(DD(reverse_V)), hi = g(DD(forward_V));
    if (!error.empty()) return std::unexpected(error);
    if (!(lo.hi <= 0.0 && hi.hi >= 0.0)) {
        return std::unexpected("bias_at_inverse_multiplication: the bracket does not hold the value (1 / M is " +
                               std::to_string(lo.to_double() + inverse) + " and " +
                               std::to_string(hi.to_double() + inverse) + " at its ends)");
    }
    const BrentResult r = brent_zero(g, DD(reverse_V), DD(forward_V), 1e-10);
    if (!error.empty()) return std::unexpected(error);
    if (r.status != BrentStatus::ok) return std::unexpected("bias_at_inverse_multiplication: Brent did not converge");
    return r.root.to_double();
}

}  // namespace NiTCAD::verify::reference
