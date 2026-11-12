#include "NiTCAD/verify/reference/mos_hf.hpp"

#include <algorithm>
#include <cmath>

#include "NiTCAD/verify/reference/brent.hpp"
#include "NiTCAD/verify/reference/first_integral.hpp"
#include "NiTCAD/verify/reference/gauss_kronrod.hpp"

namespace NiTCAD::verify::reference {

namespace {

constexpr double quad_epsrel = 1e-24;
constexpr double delta = 1e-3;  // the second method's step in t

bool positive_finite(double v) { return std::isfinite(v) && v > 0.0; }

std::expected<void, std::string> check(const MosInput& in) {
    if (!positive_finite(in.q_C) || !positive_finite(in.eps_F_per_cm) || !positive_finite(in.thermal_voltage_V) ||
        !positive_finite(in.ni_cm3) || !std::isfinite(in.net_doping_cm3) || !positive_finite(in.oxide_capacitance_F_per_cm2) ||
        !std::isfinite(in.flat_band_V)) {
        return std::unexpected("mos_hf: an input is not finite, or not positive where it must be");
    }
    if (!(in.net_doping_cm3 < 0.0)) return std::unexpected("mos_hf: the substrate must be p-type (net doping < 0)");
    return {};
}

// The redistribution model on one substrate.
struct Model {
    NeutralBulk b;
    DD q, eps, VT, k, scale;  // scale = 2 q / (eps V_T): g^2 = scale Phi

    explicit Model(const MosInput& in)
        : b(mos_bulk(in)),
          q(in.q_C),
          eps(in.eps_F_per_cm),
          VT(in.thermal_voltage_V),
          k(sqrt(DD(2.0) * q * eps * VT)),
          scale(DD(2.0) * q / (eps * VT)) {}

    [[nodiscard]] DD nv(const DD& v) const { return b.n_b * exp(-v); }
    [[nodiscard]] DD phi(const DD& t, const DD& v) const {
        return b.p_b * poisson_boltzmann_G(-t) + nv(v) * poisson_boltzmann_G(t);
    }
    [[nodiscard]] DD phi_t(const DD& t, const DD& v) const { return -b.p_b * expm1(-t) + nv(v) * expm1(t); }
    [[nodiscard]] DD phi_v(const DD& t, const DD& v) const { return -nv(v) * poisson_boltzmann_G(t); }
    [[nodiscard]] static DD sgn(const DD& t) { return t > DD(0.0) ? DD(1.0) : DD(-1.0); }
    [[nodiscard]] DD surface(const DD& t, const DD& v) const { return -sgn(t) * k * sqrt(phi(t, v)); }
    // h(s) = n_b e^{-v} expm1(s) / g(s, v)
    [[nodiscard]] DD h(const DD& s, const DD& v) const { return nv(v) * expm1(s) / sqrt(scale * phi(s, v)); }
};

struct Integral {
    DD value;
    double abs_error = 0.0;
    bool ok = true;
};

template <class F>
Integral over(F&& f, const DD& t) {
    const DD lo = t < DD(0.0) ? t : DD(0.0), hi = t < DD(0.0) ? DD(0.0) : t;
    const AdaptiveResult r = qag15(f, lo, hi, 0.0, quad_epsrel);
    // Roundoff at the double-double floor is accepted, as elsewhere in verify/reference; other statuses are not.
    return {r.value, r.abs_error, r.status == QuadratureStatus::ok || r.status == QuadratureStatus::roundoff};
}

}  // namespace

std::expected<MosCharges, std::string> mos_redistribution_charges(const MosInput& in, const DD& t, const DD& v) {
    if (const auto c = check(in); !c) return std::unexpected(c.error());
    const Model m(in);
    const Integral I = over([&](const DD& s) { return m.h(s, v); }, t);
    if (!I.ok) return std::unexpected("mos_redistribution_charges: the Q_n quadrature failed");
    return MosCharges{m.surface(t, v), -m.q * I.value, m.q.to_double() * I.abs_error};
}

std::expected<MosHfPoint, std::string> mos_hf_point(const MosInput& in, double gate_V, bool second_method) {
    if (const auto c = check(in); !c) return std::unexpected(c.error());
    const auto dc = mos_point(in, gate_V);
    if (!dc) return std::unexpected(dc.error());
    const Model m(in);
    const DD t = dc->u_s - m.b.u_b, v0(0.0);
    if (t.hi == 0.0) return std::unexpected("mos_hf_point: flat band exactly (t_s = 0)");
    MosHfPoint p;
    p.dc = *dc;
    const Integral I = over([&](const DD& s) { return m.h(s, v0); }, t);
    const Integral I2 = over(
        [&](const DD& s) {
            const DD ph = m.phi(s, v0);
            return m.h(s, v0) * (DD(1.0) + m.phi_v(s, v0) / ldexp(ph, 1));
        },
        t);
    if (!I.ok || !I2.ok) return std::unexpected("mos_hf_point: a Q_n quadrature failed");
    p.minority_charge_C_per_cm2 = -m.q * I.value;
    const DD sq = sqrt(m.phi(t, v0));
    const DD dQs_dt = -Model::sgn(t) * m.k * m.phi_t(t, v0) / ldexp(sq, 1);
    const DD dQs_dv = -Model::sgn(t) * m.k * m.phi_v(t, v0) / ldexp(sq, 1);
    const DD dQn_dt = -m.q * Model::sgn(t) * m.h(t, v0);
    const DD dQn_dv = m.q * I2.value;
    p.dv_dt = -dQn_dt / dQn_dv;
    p.quadrature_error = I2.value.hi != 0.0 ? I2.abs_error / std::abs(I2.value.hi) : 0.0;
    p.surface_capacitance_hf = -(dQs_dt + dQs_dv * p.dv_dt) / m.VT;
    const DD Cox(in.oxide_capacitance_F_per_cm2);
    p.capacitance_hf = p.surface_capacitance_hf * Cox / (p.surface_capacitance_hf + Cox);
    if (!second_method) return p;

    // Q_s along Q_n = Q_n(t_s, 0), at t_s +- d.
    const DD Qn0 = p.minority_charge_C_per_cm2;
    bool failed = false;
    const auto constrained = [&](const DD& tt) {
        const auto mismatch = [&](const DD& v) {
            const Integral J = over([&](const DD& s) { return m.h(s, v); }, tt);
            if (!J.ok) failed = true;
            return (-m.q * J.value - Qn0) / abs(Qn0);
        };
        const DD guess = p.dv_dt * (tt - t);
        const BrentResult r = brent_zero(mismatch, guess - DD(2.0), guess + DD(2.0), 1e-28);
        if (r.status != BrentStatus::ok) failed = true;
        return m.surface(tt, r.root);
    };
    const auto D = [&](double d) {
        return -(constrained(t + DD(d)) - constrained(t - DD(d))) / (DD(2.0 * d) * m.VT);
    };
    const DD Cs2 = (DD(4.0) * D(0.5 * delta) - D(delta)) / DD(3.0);
    if (failed) return std::unexpected("mos_hf_point: the second method's Brent or quadrature failed");
    const DD C2 = Cs2 * Cox / (Cs2 + Cox);
    p.capacitance_hf_second_method = C2.to_double();
    p.second_method_difference = std::abs((C2 - p.capacitance_hf).to_double());
    p.second_method_evaluated = true;
    return p;
}

std::expected<DeepDepletionPoint, std::string> mos_deep_depletion_point(const MosInput& in, double gate_V) {
    if (const auto c = check(in); !c) return std::unexpected(c.error());
    if (!std::isfinite(gate_V)) return std::unexpected("mos_deep_depletion_point: gate bias not finite");
    const Model m(in);
    const DD NA(-in.net_doping_cm3), VT = m.VT, Cox(in.oxide_capacitance_F_per_cm2), VFB(in.flat_band_V);
    const auto charge = [&](const DD& t) {
        if (t.hi == 0.0) return DD(0.0);
        return -Model::sgn(t) * m.k * sqrt(NA * poisson_boltzmann_G(-t));
    };
    const auto mismatch = [&](const DD& t) { return VFB + VT * t - charge(t) / Cox - DD(gate_V); };
    // |V_T t| <= |V_G - V_FB| (the oxide drop has the band bending's sign).
    const DD a = DD(std::abs(gate_V - in.flat_band_V) / in.thermal_voltage_V + 1.0);
    const BrentResult r = brent_zero(mismatch, -a, a, u_dd * a.hi);
    if (r.status != BrentStatus::ok) return std::unexpected("mos_deep_depletion_point: Brent did not converge");
    DeepDepletionPoint p;
    p.gate_V = DD(gate_V);
    p.t_s = r.root;
    p.brent_evaluations = r.evaluations;
    p.charge_C_per_cm2 = charge(p.t_s);
    // The second method: the charge density over the profile.
    const Integral I = over(
        [&](const DD& s) { return m.q * NA * expm1(-s) / sqrt(m.scale * NA * poisson_boltzmann_G(-s)); }, p.t_s);
    if (!I.ok) return std::unexpected("mos_deep_depletion_point: the charge quadrature failed");
    p.charge_second_method = I.value.to_double();
    p.second_method_difference = std::abs((I.value - p.charge_C_per_cm2).to_double());
    return p;
}

}  // namespace NiTCAD::verify::reference
