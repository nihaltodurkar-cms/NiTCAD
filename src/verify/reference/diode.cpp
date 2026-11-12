#include "NiTCAD/verify/reference/diode.hpp"

#include <cmath>
#include <vector>

#include "NiTCAD/verify/reference/first_integral.hpp"
#include "NiTCAD/verify/reference/gauss_kronrod.hpp"
#include "NiTCAD/verify/reference/minority_diffusion.hpp"
#include "NiTCAD/verify/reference/qe_junction.hpp"
#include "NiTCAD/verify/reference/tanh_sinh.hpp"

namespace NiTCAD::verify::reference {

namespace {

bool positive(double v) { return std::isfinite(v) && v > 0.0; }

// coth(x), x > 0, without cancellation: (1 + e^-2x) / (1 - e^-2x).
DD coth(const DD& x) {
    const DD e2 = exp(-ldexp(x, 1));
    return (DD(1.0) + e2) / -expm1(-ldexp(x, 1));
}

// The minority diffusion current of one base per unit (e^v - 1), divided by q: D m0 / L coth(W' / L), or D m0 / W'
// without recombination, m0 the minority density. n0, p0 are the bulk's electron and hole densities (for SRH's
// linearization). Sets the diffusion length (0 without srh).
DD base_term(const DiodeSide& s, const DD& minority, double n0, double p0, double ni, const DD& base, bool srh,
             DD& length) {
    const DD D(s.minority_diffusivity_cm2_s);
    if (!srh) {
        length = DD(0.0);
        return D * minority / base;
    }
    const DD inv_tau = srh_linear_inverse_lifetime(n0, p0, ni, s.tau_n_s, s.tau_p_s);
    length = sqrt(D / inv_tau);
    return D * minority / length * coth(base / length);
}

constexpr double sns_epsrel = 1e-24;    // Gauss-Kronrod in u
constexpr double second_epsrel = 1e-26;  // tanh-sinh in u

}  // namespace

std::expected<DiodeCurrent, std::string> diode_current(const DiodeInput& in, double bias_V) {
    const PnJunctionInput& j = in.junction;
    if (!positive(in.p_side.length_cm) || !positive(in.n_side.length_cm) ||
        !positive(in.p_side.minority_diffusivity_cm2_s) || !positive(in.n_side.minority_diffusivity_cm2_s)) {
        return std::unexpected("diode_current: lengths and diffusivities must be finite and positive");
    }
    if (in.srh && (!positive(in.p_side.tau_n_s) || !positive(in.p_side.tau_p_s) || !positive(in.n_side.tau_n_s) ||
                   !positive(in.n_side.tau_p_s))) {
        return std::unexpected("diode_current: lifetimes must be finite and positive with srh");
    }
    const auto qe = qe_junction(j, bias_V);  // checks the junction's inputs and the bias range
    if (!qe) return std::unexpected("diode_current: " + qe.error());
    const auto eq = pn_junction(j);
    if (!eq) return std::unexpected("diode_current: equilibrium: " + eq.error());
    PnJunctionInput shifted = j;
    shifted.ni_cm3 = qe->m_cm3.to_double();
    const auto jm = pn_junction(shifted);
    if (!jm) return std::unexpected("diode_current: at bias: " + jm.error());

    DiodeCurrent out;
    out.m_cm3 = qe->m_cm3;
    out.edge_p_cm = jm->depletion_p_cm;
    out.edge_n_cm = jm->depletion_n_cm;
    out.base_p_cm = DD(in.p_side.length_cm) - out.edge_p_cm;
    out.base_n_cm = DD(in.n_side.length_cm) - out.edge_n_cm;
    if (!(out.base_p_cm.hi > 0.0) || !(out.base_n_cm.hi > 0.0)) {
        return std::unexpected("diode_current: a depletion edge reaches a contact (W' <= 0)");
    }
    out.debye_p_cm = jm->debye_p_cm;
    out.debye_n_cm = jm->debye_n_cm;
    out.n_p0_cm3 = eq->n_p;
    out.p_n0_cm3 = eq->p_n;
    const DD VT(j.thermal_voltage_V), v = DD(bias_V) / VT, q(j.q_C);
    const DD factor = expm1(v);  // e^v - 1
    const DD per_q = base_term(in.p_side, eq->n_p, eq->n_p.to_double(), eq->p_p.to_double(), j.ni_cm3, out.base_p_cm,
                               in.srh, out.diffusion_length_p_cm) +
                     base_term(in.n_side, eq->p_n, eq->n_n.to_double(), eq->p_n.to_double(), j.ni_cm3, out.base_n_cm,
                               in.srh, out.diffusion_length_n_cm);
    out.diffusion_A_per_cm2 = q * factor * per_q;
    out.total_A_per_cm2 = out.diffusion_A_per_cm2;
    if (!in.srh) return out;

    // SNS: U on the reference's QE profile (frame of jm: n = m e^u, p = m e^-u), np - n_i^2 = n_i^2 (e^v - 1).
    const DD m(shifted.ni_cm3), ni(j.ni_cm3);
    const DD excess = ni * ni * factor;
    const auto U = [&](const DD& u, const DiodeSide& s) {
        const DD n = m * exp(u), p = m * exp(-u);
        return excess / (DD(s.tau_p_s) * (n + ni) + DD(s.tau_n_s) * (p + ni));
    };
    const std::vector<double> edges{-out.edge_p_cm.to_double(), out.edge_n_cm.to_double()};
    const auto at_edges = pn_potential(*jm, edges);
    if (!at_edges) return std::unexpected("diode_current: edge potential: " + at_edges.error());
    const DD u_ep = (*at_edges)[0].psi_V / VT, u_en = (*at_edges)[1].psi_V / VT;
    const NeutralBulk bn = neutral_bulk(shifted.ni_cm3, j.donors_cm3);
    const NeutralBulk bp = neutral_bulk(shifted.ni_cm3, -j.acceptors_cm3);
    const DD k = ldexp(q / (DD(j.eps_F_per_cm) * VT), 1);  // u'^2 = k Phi
    const auto on_n = [&](const DD& u) { return U(u, in.n_side) / sqrt(k * first_integral_phi(bn, u)); };
    const auto on_p = [&](const DD& u) { return U(u, in.p_side) / sqrt(k * first_integral_phi(bp, u)); };
    const AdaptiveResult rn = qag15(on_n, jm->u_j, u_en, 0.0, sns_epsrel);
    const AdaptiveResult rp = qag15(on_p, u_ep, jm->u_j, 0.0, sns_epsrel);
    if (rn.status != QuadratureStatus::ok || rp.status != QuadratureStatus::ok) {
        return std::unexpected("diode_current: SNS quadrature did not converge");
    }
    out.sns_A_per_cm2 = q * (rn.value + rp.value);
    out.sns_quadrature_error = j.q_C * (rn.abs_error + rp.abs_error);

    // Second method (16.3): the same integrals by tanh-sinh. The change of variable dx = du / |u'| is checked
    // against the x-space integral at L0 (diode_test).
    // tanh-sinh's error estimate is absolute and meant for integrals of order one (tanh_sinh.hpp), so each integrand
    // is divided by its magnitude at the junction (the interval in u is O(1) to O(100)).
    const DD sn = abs(on_n(jm->u_j)), sp = abs(on_p(jm->u_j));
    const TanhSinhResult tn = tanh_sinh([&](const DD& u) { return on_n(u) / sn; }, jm->u_j, u_en, 0.0, second_epsrel);
    const TanhSinhResult tp = tanh_sinh([&](const DD& u) { return on_p(u) / sp; }, u_ep, jm->u_j, 0.0, second_epsrel);
    if (tn.status != TanhSinhStatus::ok || tp.status != TanhSinhStatus::ok) {
        return std::unexpected("diode_current: SNS second method (tanh-sinh) did not converge");
    }
    const DD second = q * (tn.value * sn + tp.value * sp);
    out.sns_second_method_difference = std::abs((second - out.sns_A_per_cm2).hi);
    out.total_A_per_cm2 = out.diffusion_A_per_cm2 + out.sns_A_per_cm2;
    return out;
}

}  // namespace NiTCAD::verify::reference
