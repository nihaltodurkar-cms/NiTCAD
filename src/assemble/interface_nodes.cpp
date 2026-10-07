#include "NiTCAD/assemble/interface_nodes.hpp"

#include <cmath>
#include <utility>

#include "NiTCAD/physics/interface_traps.hpp"
#include "NiTCAD/physics/recombination.hpp"

namespace NiTCAD::assemble {

InterfaceNodes::InterfaceNodes(std::vector<InterfaceNode> nodes,
                               std::vector<InterfaceLevel> levels, std::size_t interfaces)
    : nodes_(std::move(nodes)), levels_(std::move(levels)), interfaces_(interfaces) {}

InterfaceNodes::EquilibriumCharge InterfaceNodes::equilibrium_charge(const InterfaceNode& v,
                                                                     double eta) const noexcept {
    double value = v.fixed_charge_cm2, d_eta = 0.0;
    for (std::size_t k = v.first_level; k < v.last_level; ++k) {
        const InterfaceLevel& l = levels_[k];
        const physics::FermiOccupancy f = physics::fermi_occupancy(l.tau - eta);
        // Donor N (1 - f), acceptor -N f: both fall with eta at the rate N f (1 - f).
        value += l.density_cm2 * (l.donor ? f.empty : -f.occupied);
        d_eta -= l.density_cm2 * f.d_eta;
    }
    return {v.charge_weight * value, v.charge_weight * d_eta};
}

InterfaceTerms InterfaceNodes::terms(const InterfaceNode& v, double n, double p, double n_ie,
                                     double log_gamma_n, double d_log_gamma_n,
                                     double log_gamma_p, double d_log_gamma_p, double product,
                                     double product_dn, double product_dp) const noexcept {
    double q = v.fixed_charge_cm2, q_dn = 0.0, q_dp = 0.0;
    double r = 0.0, r_dn = 0.0, r_dp = 0.0;
    for (std::size_t k = v.first_level; k < v.last_level; ++k) {
        const InterfaceLevel& l = levels_[k];
        const double n1 = n_ie * std::exp(log_gamma_n + l.tau);
        const double p1 = n_ie * std::exp(log_gamma_p - l.tau);
        const physics::TrapKinetics t = physics::trap_kinetics(
            n, p, n1, n1 * d_log_gamma_n, p1, p1 * d_log_gamma_p, l.cn, l.cp);
        // Donor N (1 - f), acceptor -N f; d(1 - f) = -df.
        q += l.density_cm2 * (l.donor ? t.empty : -t.occupied);
        q_dn -= l.density_cm2 * t.d_occupied_dn;
        q_dp -= l.density_cm2 * t.d_occupied_dp;
        r += l.density_cm2 * t.rate;
        r_dn += l.density_cm2 * t.d_rate_dn;
        r_dp += l.density_cm2 * t.d_rate_dp;
    }
    if (v.velocity_n > 0.0 && v.velocity_p > 0.0) {
        const physics::RecombinationRate s = physics::srh_recombination(
            n, p, {product, product_dn, product_dp}, n_ie, 1.0 / v.velocity_n,
            1.0 / v.velocity_p);
        r += s.rate;
        r_dn += s.d_dn;
        r_dp += s.d_dp;
    }
    const double cw = v.charge_weight, rw = v.rate_weight;
    return {cw * q, cw * q_dn, cw * q_dp, rw * r, rw * r_dn, rw * r_dp};
}

}  // namespace NiTCAD::assemble
