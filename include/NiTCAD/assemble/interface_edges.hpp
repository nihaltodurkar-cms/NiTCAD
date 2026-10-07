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
//
// Trap dynamics (Unit 21): in a time step the occupancy of each level obeys
//     df/dt = cn n_I (1 - f) - cn n1 f - cp p_I f + cp p1 (1 - f) = occ - f D
// (occ = cn n_I + cp p1, D = cn (n_I + n1) + cp (p_I + p1), physical time and densities), and a
// BDF step f - c = beta h (occ - f D) (c the combination of earlier occupancies) is linear in f:
//     f = (c + k occ) / (1 + k D),   k = beta h Ns (scaled densities in occ and D),
// a weighted mean of c and the steady occupancy occ / D, so f stays within [min(c, 0), max(c, 1)]
// and Q within the bounds that follow (the bracket of the local solve). It is eliminated inside the
// local solve, so the global system keeps its unknowns; as k grows it tends to the steady-state
// occupancy. Electron and hole capture then differ by N df/dt: the electron row takes the net
// electron capture cn (n_I (1 - f) - n1 f) and the hole row the net hole capture
// cp (p_I f - p1 (1 - f)), each summed over N_k (InterfaceDrift::rate and rate_p).
#pragma once

#include <cstddef>
#include <span>
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

// The drift-diffusion interface of an edge: as InterfaceEquilibrium with the recombination terms
// (the semiconductor's electron row takes -rate, its hole row +rate_p; in steady state the two are
// equal), all partials in (psi_insulator, psi_semiconductor, n_semiconductor, p_semiconductor).
struct InterfaceDrift {
    double psi;
    double flux_insulator, flux_semiconductor, rate, rate_p;
    double d_flux_insulator[4], d_flux_semiconductor[4], d_rate[4], d_rate_p[4];
    double trapped;
};

// A time step's trap history for one edge (Unit 21): k = beta h Ns in the units of the levels'
// cn and cp (s / cm^3 times cm^3/s), and c per level of the edge's interface, in level order.
struct TrapStep {
    double weight;
    std::span<const double> history;
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
    // With `step` the traps follow its history (see the header comment), else the steady state.
    // Precondition (NITCAD_EXPECTS): the history has one entry per level of the edge.
    [[nodiscard]] InterfaceDrift drift(const InterfaceEdge& e, double psi_insulator,
                                       double psi_semiconductor, double n, double p,
                                       const InterfaceStatistics& s,
                                       const TrapStep* step = nullptr) const;
    // The electron occupancy f of each level of the edge at the drift() state, into `f` (one entry
    // per level, in level order).
    void occupancies(const InterfaceEdge& e, double psi_insulator, double psi_semiconductor,
                     double n, double p, const InterfaceStatistics& s, const TrapStep* step,
                     std::span<double> f) const;

    // Trap slots: one per (interface edge, level of its interface); edge k's slots start at
    // slot_offset(k).
    [[nodiscard]] std::size_t slot_count() const noexcept { return slots_; }
    [[nodiscard]] std::size_t slot_offset(std::size_t edge) const noexcept {
        return slot_offset_[edge];
    }

private:
    struct Charge {  // Q and its partials, at one interface state; rate_e / rate_h the rows' terms
        double value, d_n, d_p, rate, rate_n, rate_p, rate_h, rate_hn, rate_hp;
    };
    // Q (with the fixed charge) and the rates at interface densities n_I, p_I; with `step` the
    // occupancies of the step, written to `f` when it is not empty.
    [[nodiscard]] Charge charge_at(const InterfaceEdge& e, double n, double p,
                                   const InterfaceStatistics& s, const TrapStep* step = nullptr,
                                   std::span<double> f = {}) const;
    // drift() and occupancies(): the occupancies go to `f` when it is not empty.
    [[nodiscard]] InterfaceDrift drift_at(const InterfaceEdge& e, double psi_i, double psi_s,
                                          double n_s, double p_s, const InterfaceStatistics& s,
                                          const TrapStep* step, std::span<double> f) const;

    std::vector<InterfaceEdge> edges_;
    std::vector<InterfaceLevel> levels_;
    std::vector<std::size_t> slot_offset_;
    std::size_t interfaces_ = 0;
    std::size_t slots_ = 0;
};

}  // namespace NiTCAD::assemble
