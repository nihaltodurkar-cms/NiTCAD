// Interface charge, traps and surface recombination at the interface potential (Unit 15b), shared
// by both assemblers as GateNodes is (gate.hpp).
//
// A declared semiconductor-insulator interface (device::Interface) carrying a fixed charge, traps
// or surface recombination acts on each mesh edge joining its insulator and semiconductor nodes.
// The interface lies at the edge's midpoint (mesh::straddle_interface), so the edge is two
// half-edges in series, of conductance g_i = 2 et_i c and g_s = 2 et_s c (c the edge's scaled
// coupling, et the
// end's permittivity over the reference one), and the interface potential psi_I between them obeys
// Gauss's law on a thin box around the interface patch of the edge (area A, a = A / L_D^(D-1)):
//
//     g_i (psi_I - psi_i) + g_s (psi_I - psi_s) = Q(psi_I),
//     Q = a / (Ns L_D) (Q_f + sum_k N_k q_k),     q_k = 1 - f_k (donor), -f_k (acceptor).
//
// The insulator node's Poisson row takes the half-edge flux g_i (psi_I - psi_i) and the
// semiconductor node's g_s (psi_I - psi_s) in place of the edge's flux; with Q = 0 the two are the
// edge's harmonic-mean flux. The semiconductor node's continuity rows take -/+ a Ns / (R0 L_D)
// (sum_k N_k r_k + r_s) with r the trap and surface-velocity rates (physics/interface_traps.hpp),
// homogeneous of degree one so U [cm^-2 s^-1] = Ns r(n', p') on scaled densities.
// - In thermal equilibrium (EquilibriumPoisson) f_k is the Fermi function of tau_k - eta_I,
//   eta_I = psi_I + s (the semiconductor's band shift), and nothing recombines.
// - In drift-diffusion the carriers keep their quasi-Fermi levels across the semiconductor's
//   half-cell: n_I is the density (of the selected statistics) at the reduced energy
//   ln(n_s / n_ie) - ln gamma_n(n_s) + psi_I - psi_s (Boltzmann: n_I = n_s e^(psi_I - psi_s)), p_I
//   likewise with the opposite sign; f_k is the steady-state SRH occupancy at (n_I, p_I), with
//   n1 = gamma_n(n_I) n_ie e^tau and p1 = gamma_p(p_I) n_ie e^-tau. At equilibrium this is the
//   Fermi function of tau_k - eta_I.
// Q falls as psi_I rises (more electrons, fewer holes at the interface), and it is bounded by the
// fixed charge and the trap densities, so the local equation has one root in a known bracket; it is
// solved by safeguarded Newton to rounding. The global system keeps its unknowns: psi_I is
// eliminated, and its derivatives with respect to the edge's end unknowns follow from the local
// equation (implicit function), so the Jacobian is exact.
// A Dirichlet end (an electrode or ohmic contact node) takes none of the edge's terms; the edge's
// psi_I still depends on it.
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

struct InterfaceEdge {
    std::size_t edge;                     // mesh edge
    std::size_t insulator, semiconductor; // its end nodes
    std::size_t interface;                // index in device.interfaces()
    double g_insulator, g_semiconductor;  // half-edge conductances (scaled)
    double charge_weight;                 // a / (Ns L_D)
    double rate_weight;                   // a Ns / (R0 L_D)
    double fixed_charge_cm2;              // Q_f
    double velocity_n, velocity_p;        // s_n, s_p [cm/s]
    std::size_t first_level, last_level;  // the interface's levels in InterfaceEdges::levels()
    double charge_low, charge_high;       // bounds of Q (fixed charge, all traps charged one way)
};

// The equilibrium interface of an edge at its end potentials: psi_I, the half-edge fluxes into the
// insulator and semiconductor rows and their partials in (psi_insulator, psi_semiconductor), and
// the trapped charge (Q without the fixed charge).
struct InterfaceEquilibrium {
    double psi;
    double flux_insulator, flux_semiconductor;
    double d_flux_insulator[2], d_flux_semiconductor[2];
    double trapped;
};

// The drift-diffusion interface of an edge: as InterfaceEquilibrium with the recombination term r
// (the semiconductor's continuity rows take -/+ r), all partials in (psi_insulator,
// psi_semiconductor, n_semiconductor, p_semiconductor).
struct InterfaceDrift {
    double psi;
    double flux_insulator, flux_semiconductor, rate;
    double d_flux_insulator[4], d_flux_semiconductor[4], d_rate[4];
    double trapped;
};

// The semiconductor node's statistics: Fermi-Dirac or Boltzmann, scaled n_ie, ln(Nc / n_ie) and
// ln(Nv / n_ie), and its band shift s.
struct InterfaceStatistics {
    bool fermi_dirac;
    double n_ie, log_dos_n, log_dos_p, band_shift;
};

class InterfaceEdges {
public:
    InterfaceEdges() = default;
    InterfaceEdges(std::vector<InterfaceEdge> edges, std::vector<InterfaceLevel> levels,
                   std::size_t interfaces);

    [[nodiscard]] const std::vector<InterfaceEdge>& edges() const noexcept { return edges_; }
    [[nodiscard]] const std::vector<InterfaceLevel>& levels() const noexcept { return levels_; }
    [[nodiscard]] bool empty() const noexcept { return edges_.empty(); }
    [[nodiscard]] std::size_t interface_count() const noexcept { return interfaces_; }

    [[nodiscard]] InterfaceEquilibrium equilibrium(const InterfaceEdge& e, double psi_insulator,
                                                   double psi_semiconductor,
                                                   const InterfaceStatistics& s) const;
    [[nodiscard]] InterfaceDrift drift(const InterfaceEdge& e, double psi_insulator,
                                       double psi_semiconductor, double n, double p,
                                       const InterfaceStatistics& s) const;

private:
    struct Charge {  // Q and its partials, at one interface state
        double value, d_n, d_p, rate, rate_n, rate_p;
    };
    // Q (with the fixed charge) and r at interface densities n_I, p_I.
    [[nodiscard]] Charge charge_at(const InterfaceEdge& e, double n, double p,
                                   const InterfaceStatistics& s) const;

    std::vector<InterfaceEdge> edges_;
    std::vector<InterfaceLevel> levels_;
    std::size_t interfaces_ = 0;
};

}  // namespace NiTCAD::assemble
