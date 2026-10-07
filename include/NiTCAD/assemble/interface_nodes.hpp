// Interface charge, traps and surface recombination in scaled variables (Unit 15b), shared by both
// assemblers as GateNodes is (gate.hpp).
//
// A declared semiconductor-insulator interface (device::Interface) acts on each semiconductor node
// i that has an edge to an insulator node of it, with the node's share of the interface area
// A_i = the summed coupling area of those edges (the interface lies at their midpoints). In the
// node's rows (sheet densities in cm^-2, a_i = A_i / L_D^(D-1)):
//     Poisson    + a_i / (Ns L_D) (Q_f + sum_k N_k q_k)
//     electrons  - a_i Ns / (R0 L_D) (sum_k N_k r_k + r_s)
//     holes      + the same
// with, per trap level k (physics/interface_traps.hpp), q_k = 1 - f_k for a donor and -f_k for an
// acceptor, r_k = U_k / N_k on scaled densities, and r_s the surface recombination velocity term
// (both velocities positive; otherwise none). Both are homogeneous of degree one, so U [cm^-2 s^-1]
// = Ns r(n', p') on scaled densities n', p'. Trap bands enter as their quadrature levels.
// - In thermal equilibrium (EquilibriumPoisson) the occupancy is the Fermi function of
//   tau_k - eta, eta = psi + s, and nothing recombines.
// - In drift-diffusion it is the steady-state SRH occupancy at the node's n and p, with
//   n1 = gamma_n n_ie e^tau and p1 = gamma_p n_ie e^-tau (gamma 1 under Boltzmann statistics).
// An ohmic contact node is a Dirichlet node and takes no interface terms.
#pragma once

#include <cstddef>
#include <vector>

namespace NiTCAD::assemble {

struct InterfaceLevel {
    bool donor;          // donor-like (else acceptor-like)
    double density_cm2;  // N_k
    double tau;          // (E_t - E_i) / kT
    double cn, cp;       // sigma v_th [cm^3/s]
};

struct InterfaceNode {
    std::size_t node;
    std::size_t interface;          // index in device.interfaces()
    double charge_weight;           // a_i / (Ns L_D)
    double rate_weight;             // a_i Ns / (R0 L_D)
    double fixed_charge_cm2;        // Q_f
    double velocity_n, velocity_p;  // s_n, s_p [cm/s]
    std::size_t first_level, last_level;  // the interface's levels in InterfaceNodes::levels()
};

// The interface terms of one node, scaled: the Poisson row term (charge) and the recombination
// term r (the continuity rows take -/+ r), with partials in n and p.
struct InterfaceTerms {
    double charge, charge_dn, charge_dp;
    double rate, rate_dn, rate_dp;
};

class InterfaceNodes {
public:
    InterfaceNodes() = default;
    InterfaceNodes(std::vector<InterfaceNode> nodes, std::vector<InterfaceLevel> levels,
                   std::size_t interfaces);

    [[nodiscard]] const std::vector<InterfaceNode>& nodes() const noexcept { return nodes_; }
    [[nodiscard]] const std::vector<InterfaceLevel>& levels() const noexcept { return levels_; }
    [[nodiscard]] bool empty() const noexcept { return nodes_.empty(); }

    // Equilibrium: the Poisson row term and its eta derivative (trapped charge by the Fermi
    // function of tau - eta).
    struct EquilibriumCharge {
        double value, d_eta;
    };
    [[nodiscard]] EquilibriumCharge equilibrium_charge(const InterfaceNode& node,
                                                       double eta) const noexcept;

    // Drift-diffusion: the terms at scaled densities n and p of a node with scaled n_ie. ln gamma
    // and d ln gamma / d density of each density (0 under Boltzmann statistics), and the node's
    // equilibrium product (statistics.hpp, scaled) for the velocity term.
    [[nodiscard]] InterfaceTerms terms(const InterfaceNode& node, double n, double p, double n_ie,
                                       double log_gamma_n, double d_log_gamma_n,
                                       double log_gamma_p, double d_log_gamma_p,
                                       double product, double product_dn,
                                       double product_dp) const noexcept;

    // The scaled trapped charge (without Q_f) of every interface, in units of q Ns L_D^D: the sum
    // over its nodes of charge_weight sum_k N_k q_k. `occupancy(node, level)` returns f.
    template <class Occupancy>
    [[nodiscard]] std::vector<double> trapped_charges(Occupancy occupancy) const {
        std::vector<double> q(interfaces_, 0.0);
        for (const InterfaceNode& v : nodes_) {
            double sum = 0.0;
            for (std::size_t k = v.first_level; k < v.last_level; ++k) {
                const InterfaceLevel& l = levels_[k];
                const auto [f, g] = occupancy(v, l);
                sum += l.density_cm2 * (l.donor ? g : -f);
            }
            q[v.interface] += v.charge_weight * sum;
        }
        return q;
    }

private:
    std::vector<InterfaceNode> nodes_;
    std::vector<InterfaceLevel> levels_;
    std::size_t interfaces_ = 0;
};

}  // namespace NiTCAD::assemble
