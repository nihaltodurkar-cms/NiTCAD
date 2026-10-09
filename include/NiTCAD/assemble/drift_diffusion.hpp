// Coupled drift-diffusion system: Poisson, electron and hole continuity, three scaled unknowns per
// node interleaved as x[3i] = psi, x[3i+1] = n, x[3i+2] = p (legacy Device1D::residual_jacobian,
// baseline models, generalised from 1D to the mesh graph).
//
// Rows of node i not on a contact, summing over its edges e = (a, b) (+ at a, - at b):
//     Poisson    sum  et c_e (psi_b - psi_a)                    - V_i (n_i - p_i - C_i)
//     electrons  sum  Jn_e                                        - V_i R_i
//     holes      sum  Jp_e                                        + V_i R_i
// with c_e the scaled coupling (EquilibriumPoisson), Jn_e and Jp_e the Scharfetter-Gummel fluxes of
// sg_flux.hpp with edge factor D_e c_e, D_e = hmean(mu_a, mu_b) V_T / D0 (harmonic mean of the
// nodes' mobilities, legacy dn_edge), and R = (R_SRH + R_Auger) / R0. SRH is evaluated on scaled
// densities (it is homogeneous of degree one, so R_SRH / R0 = R_SRH(n', p', ...) Ns / R0); Auger,
// which is cubic, on physical ones. Both drive n p towards the equilibrium product of the selected
// statistics (n_ie^2; Fermi-Dirac: n_ie^2 gamma_n gamma_p at the node's densities, with partials).
// With models.field_mobility the edge factor carries the Canali mobility of the edge's own field:
// D_e = mu_C(mu0_e, E_e) V_T / D0, with mu0_e the harmonic mean above and
// E_e = V_T |psi_b - psi_a| / length the electrostatic field along the edge (the field parallel to
// the edge's current, as in the legacy). The Jacobian includes d D_e / d psi, so Newton stays
// quadratic (OLD / NEW / REASON in ARCHITECTURE.md 6.2, Unit 13: the legacy averaged |E| to the
// nodes, applied Canali per node, took the harmonic mean, and lagged it).
// The driving term of the fluxes is delta_n = psi_b - psi_a + s_b - s_a + ln(n_ie,b / n_ie,a) for
// electrons and delta_p = psi_b - psi_a + s_b - s_a - ln(n_ie,b / n_ie,a) for holes (legacy delta,
// delta_p): with band-gap narrowing n_ie varies in space, and across a heterointerface the band
// shift s (EquilibriumPoisson; Unit 15) steps, and only these make the equilibrium carry no
// current.
// The edge mobility is the harmonic mean of the two nodes' also across an interface (legacy); with
// field mobility, on an edge between two materials it is the harmonic mean of each end's Canali
// mobility at the edge's field (each with its own material's parameters), on other edges the
// Canali mobility of the harmonic mean, as Unit 13.
// On the edges joining the regions of an interface the device declares thermionic_emission
// (device::Interface) the fluxes are the thermionic-emission ones of thermionic_flux.hpp instead,
// with the same driving terms, K the harmonic mean of the ends' emission velocities
// (physics::emission_velocity_cm_s of each material) times the scaled interface area, and N the
// ends' Nc (electrons) or Nv (holes); no mobility enters them. Undeclared material steps (a graded
// composition) stay drift-diffusion.
// With models.incomplete_ionization the Poisson row's doping is N_D+ - N_A- at the electrons' and
// holes' own band reduced energies (from n and p), with its n and p derivatives; with
// models.radiative, R gains B (n p - E) (physical densities, as Auger).
// With models.fermi_dirac (Unit 14, the legacy nu-factor scheme) delta_n gains
// ln gamma_n,b - ln gamma_n,a and delta_p loses ln gamma_p,b - ln gamma_p,a, each node's degeneracy
// factor taken from its own density (physics::fermi_dirac_degeneracy). Then n_b / n_a = e^delta_n
// at equilibrium, so the flux again vanishes there, and in the continuum limit the flux is
// mu n grad(phi_n): the generalized Einstein relation D / mu = (V_T n / (dn / d eta)) is implicit.
// ln gamma depends on the density, so the Jacobian's density columns gain
// d flux / d delta times d ln gamma / d density.
// Contact rows (ohmic): psi - psi0, n - n0, p - p0, from ohmic_contact_value at the contact's bias
// (psi0 less s).
// A gate node keeps its three box rows; its Poisson row gains the oxide term of gate.hpp at the
// gate's bias, G_i (psi_G,i - psi_i) + S_i, and its continuity rows no boundary flux (legacy
// Device2D GateBC: Robin on psi only).
//
// Insulators (Unit 15b): an insulator node's continuity rows are n = 0 and p = 0 (F = n, F = p),
// its Poisson row the box row with no charge, or on an electrode node the Dirichlet row
// psi - psi_E (scaled_device.hpp); no carrier flux is assembled on an edge with an insulator end,
// so none crosses a semiconductor-insulator edge. An edge of a semiconductor-insulator interface
// with charge, traps or recombination carries the half-edge fluxes of interface_edges.hpp instead
// of its own, the fixed and trapped charge sitting at the interface potential psi_I, and its
// semiconductor node's continuity rows take the traps' and surface recombination, all at the
// interface densities (the node's quasi-Fermi levels carried to psi_I) with the steady-state SRH
// occupancy. The insulator node's Poisson row then depends on the semiconductor node's n and p
// too (the pattern gains those two entries on such edges).
//
// Jn and Jp are the electron and hole current densities in units of J0 (conventional current, along
// the edge from a to b); their sum is divergence-free at every node off the ohmic contacts (gate
// nodes included), so the terminal current of an ohmic contact is the total flux on the edges
// leaving it, and that of a gate is zero.
//
// Newton update (legacy Device1D::newton): psi is clipped to +-max_update; n and p are clamped to
// [0.1, 10] times their current value, which keeps them positive. The convergence measure is the
// largest of |dpsi|, |dn| / n and |dp| / p of the FULL correction (6.2; the legacy measured the
// damped one); an insulator node's densities stay 0 and are not measured. Only during iteration:
// the converged state is never clamped (6.7).
//
// Time steps (Unit 21): a BDF step from earlier states to time t_new solves the rows above with
// the storage term of each semiconductor node off the ohmic contacts,
//     electrons  ... - V_i r (S_n,i - c_n,i),     holes  ... + V_i r (S_p,i - c_p,i),
// with r = 1 / (beta h) (h in units of t0 = Ns / R0, the time scale of the scaled rates), c the
// step's combination of the earlier storages and S the carriers a node stores: S_n = n and
// S_p = p, or with incomplete ionization S_n = n - N_D+ and S_p = p - N_A- (the electrons and
// holes bound to the dopants are stored too, so dopant ionization, instantaneous in the model,
// keeps the charge balance). Poisson stays algebraic. The interface traps follow their own history
// (interface_edges.hpp, TrapStep with weight t0 Ns / r), their electron and hole capture differing
// by the trapped charge's change. Then the conduction current entering through the contacts is
// the rate of change of the charge stored inside the device (carriers, bound carriers and traps)
// for the same BDF difference, and Gauss's law makes the charges of the contacts
// (contact_charges) balance that charge, so the total currents, conduction plus the BDF difference
// of the contact charges, sum to zero.
//
// Small signal (Unit 22): about a steady state x0, a perturbation dx e^(s t) of the contact biases
// dV e^(s t) obeys (J + s C + T(s)) dx = -dF/dV dV, with J the steady Jacobian, C the storage
// term's Jacobian (the electron row -V_i dS_n/dn, the hole row +V_i dS_p/dp, s in units of
// 1 / t0) and T(s) the interface traps' terms with their occupancy's dynamics
// (InterfaceSmallSignal), in place of the steady ones. A time step's Jacobian with rate r, from
// the history x0 and its steady trap occupancies, is this matrix at s = r: backward Euler is exact
// for the exponential e^(s t) when h = 1 / s. The total current of a contact is its conduction
// current plus s times its charge.
//
// Impact ionization (Unit 19, models.impact_ionization; legacy device.py M15 and ii_grid.py):
// the electron row of a semiconductor node off the ohmic contacts gains + V_i G / R0 and the hole
// row - V_i G / R0, with G = (alpha_n(E_n) |j_n| + alpha_p(E_p) |j_p|) / q. The field E and the
// current densities j_n, j_p at the node are vectors reconstructed by least squares from their
// components along the node's carrier edges (v = P sum_k t_k v_k, t_k the edge's unit vector, P
// the pseudo-inverse of sum_k t_k t_k^T; on a tensor grid the per-axis mean of the two edges, the
// legacy's; on any mesh exact for a uniform field). With m = sqrt(|j_c|^2 + eps_c^2), carrier c
// ionizes at the field along its current, E_c = |E . j_c| / m, and its term is
// alpha_c(E_c) |j_c|^2 / m: a field across the current does not ionize, G is the same for a
// current at any angle to the mesh, and G never exceeds alpha(|E|) |j|. For |j| >> eps these are
// |E . j^| and |j|; through j = 0 they are smooth and vanish with their gradients.
// eps_c^2 = floor^2 + (rho R_c)^2: floor = impact_current_floor (1e-12 A/cm^2), below which a
// current fades out of the generation (no generation without current, no kink as a current at
// rounding level flips its sign between Newton iterates: |j| itself made Newton cycle in a period
// two), and, optionally, rho = models.impact_current_resolution (0 by default) times R_c, the
// opposing flux terms the current is the difference of (the mean over the node's edges of
// a_k (c_a + c_b) sqrt(1 + Delta_k^2), a the edge's low-field factor, c the carrier's density,
// Delta the potential difference, as a current density). A current is a small difference of large
// fluxes where the carriers are dense: in the spill-over layer at a junction, where the field and
// alpha peak, the electron current is rounding, some 1e-14 of the fluxes. Without rho it
// generates there, a negligible amount (1e-11 against 2e-9 A/cm^2 of leakage, measured), but
// with the field at breakdown its sign flips can make a trace's corrector cycle; rho = 1e-12
// removes that and leaves the depleted region (where R is small) alone, but stalls the trace of
// an open-base transistor at its start. Reverse leakage (1e-9 A/cm^2 and more) is unaffected.
// OLD / NEW / REASON: the legacy smoothed |J| as sqrt(J^2 + eps^2) with eps 1e-6 of the largest
// edge current, which generates at zero current and ties every node's generation to that one
// edge (a dense Jacobian column), and in 2D/3D averaged the per-axis magnitudes; NEW, a vector
// reconstruction and a local resolution in a form that vanishes at zero current; REASON, no
// generation from currents the state does not resolve, a local and exact Jacobian (eps's own
// partials included), and a uniform field reconstructed exactly at any angle.
// The Jacobian is exact. A node's electron generation reads its neighbours' p (through the hole
// current) and its hole generation their n: with the model on the pattern gains those two
// entries per edge direction; with it off nothing changes.
//
// Band-to-band tunnelling (Unit 20; legacy btbt.py, btbt_grid.py, nonlocal_path.py). Pure
// generation in both models (no Hurkx D factor, no occupancy factor), so it does not vanish at
// equilibrium (ARCHITECTURE.md 6.2, Unit 20, known limits).
// - Local (models.btbt_local): the electron row of a semiconductor node off the ohmic contacts
//   gains + V_i G / R0 and the hole row - V_i G / R0, G = A F^2 exp(-B / F) with F the size of the
//   node's field reconstructed as for impact ionization (on a tensor grid the per-axis mean of the
//   node's two edges, the legacy's). Psi columns only, already in the pattern.
// - Nonlocal (models.btbt_nonlocal, tunnel_paths.hpp): each path of the frozen set (set_paths) from
//   start node i generates P = V_i G_p / R0 pairs: holes in i's hole row (- P), electrons at the
//   crossing in the electron rows of the crossing segment's two samples' stencils, (1 - t) on the
//   first and t on the second, times their stencil weights, so electrons and holes balance
//   exactly. With Delta(s) = E_v(i) - E_end(s) along the frozen crossing segment (a, b),
//   t = Delta_a / (Delta_a - Delta_b), extrapolated along the segment while the crossing moves
//   off it (smooth in psi; a relocation moves the segment). Rates:
//   - kane: E_end = E_c; the path's length l = (segments before a) + t L_ab and mean field
//     F = (psi_f - psi_i) V_T / l [V/cm], psi_f = (1 - t) psi_a + t psi_b the potential at the
//     crossing (with band-gap narrowing off in one material, psi_f - psi_i = E_g / q); G_p =
//     A F^2 exp(-B / F) with the start's (A, B). In a uniform field F is the field.
//   - direct_wkb: E_end = E_v + E_g (delta = 1, physics/band_to_band.hpp), G_p the path rate of
//     eq. (11) from the exact segment integrals over the forward samples, |dE_v/dx| over the first
//     segment, and k_m^2 from the band extrema over the path's own span (the backward polyline
//     and the start for E_vmax, the forward samples for E_cmin = min (E_v + E_g)). The span must
//     reach past the start: from the start to the crossing E_v is largest at the start itself, so
//     k_m^2 would vanish.
//   The rows of a path read the psi of its columns: the start and the crossing segment's stencils
//   (kane), every sample's stencil (direct_wkb). set_paths rebuilds the pattern with those
//   entries; without paths the pattern is unchanged. The geometry's own dependence on psi is not
//   in the Jacobian (frozen; ARCHITECTURE.md 6.2).
// With any generation model on, update_size measures densities against 1e-8 of the largest (as
// impact ionization, Unit 19).
//
// Electrothermal coupling (Unit 23, models.electrothermal; DECISIONS.md T1-T14; legacy thermal.py,
// thermal_grid.py). Each node has a fourth unknown, the scaled temperature rise
// tau = (T - T0) / T0 (T0 the device's temperature, the scaling's), interleaved as x[4i + 3];
// theta = T / T0 = 1 + tau below. The rise, not theta, is the unknown so that a small heating
// keeps its digits: theta near 1 resolves T only to about 7e-14 K at 300 K, which limits the
// heat conducted to a sink to about 4e-8 W/cm^2 across a micron (measured: the energy balance
// of a diode at 3.7e-6 W/cm^2 closed to 0.2%). Without the model the stride stays 3 and nothing
// changes. Energies below are in units of k T0, densities and fluxes as above.
// - Bands: the electron affinity is temperature independent, so the conduction-band edge
//   E_c = c_n - psi does not move with T (c_n = ln(Nc / n_ie) - s at T0, band-gap narrowing
//   included), and the valence-band edge E_v = c_p - psi + g(theta) carries the gap's change,
//   c_p = -(ln(Nv / n_ie) + s) at T0, g = -(Eg(T) - Eg(T0)) / k T0.
// - Carrier fluxes (T2, Wachutka's thermodynamic model with Boltzmann statistics, Nc and Nv
//   ~ T^(3/2), the momentum relaxation time ~ E^r): J_n = mu_n (n grad E_c' + k T grad n
//   + k (1 + r) n grad T), with E_c' = E_c - k T ln Nc0 (Nc0 the material's Nc at T0: a step of it
//   acts as a band step), and the mirror for holes. With w = theta^(1+r) n this is
//   mu_n k T^(-r) (grad w + w grad E_c' / k T); freezing T in the coefficients at the edge's mean
//   theta_e = (theta_a + theta_b) / 2 gives the Scharfetter-Gummel flux in w,
//       Jn = A_n (w_b B(delta_n) - w_a B(-delta_n)),   A_n = hmean(mu_n,a, mu_n,b) V_T / D0 g_e
//                                                              theta_e^(-r),
//       delta_n = (psi_b - psi_a + shift_n - L_c) / theta_e + L_c,
//       Jp = -A_p (w_b B(-delta_p) - w_a B(delta_p)),
//       delta_p = (psi_b - psi_a + shift_p + L_v - (g_b - g_a)) / theta_e - L_v,
//   with g_e the edge's coupling, mu at each node's own T, L_c = ln(Nc0_b / Nc0_a) and
//   L_v = ln(Nv0_b / Nv0_a) (0 in one material), r the mean of the ends' thermopower exponents.
//   At theta = 1 these are the isothermal fluxes; at any uniform theta they vanish at
//   equilibrium; between two temperatures the open-circuit voltage is the Seebeck voltage.
// - Heat (T4): each carrier carries the energy e_n = E_c + (r_n + 5/2) theta (electrons) or
//   e_p = E_v - (r_p + 5/2) theta (holes) per particle, so the energy flux along an edge, from a to
//   b, is -(ebar_n Jn + ebar_p Jp) (ebar the mean of the ends'). The heat a node's box receives
//   is the energy flowing into it: + (ebar_n Jn + ebar_p Jp) from each edge it is the first node
//   of, - from each it is the second of, plus on an ohmic contact node m (Jn + Jp) with the same
//   signs, m = V / V_T (the metal delivers carriers at its Fermi
//   level -V; that term is the contact's Peltier heat, T5). This is the discrete divergence of
//   the energy flux, so it contains the Joule, recombination ((Eg + (r_n + r_p + 5) kT) R),
//   Peltier and Thomson heat, and in a steady state the heat of all nodes is sum_c V_c I_c to
//   round-off.
// - Heat row of node i (scaled by K0 = kappa_ref T0 / (V_T J0 L_D), kappa_ref the reference
//   material's conductivity at 300 K, so conduction is O(1)): sum over its edges of the
//   conducted heat leaving it minus the heat above / K0, plus h_i (theta_i - theta_amb) on a node
//   of a thermal resistance contact (h_i = A_i L_D / (R_th kappa_ref), A_i its scaled face
//   area). Conduction along an edge of one material is exact for steady 1D conduction (the
//   Kirchhoff transform u(theta) = integral of kappa from 1 to theta):
//   (g_e / kappa_ref) (u(theta_a) - u(theta_b)); between two materials (or to an insulator)
//   (g_e / kappa_ref) hmean(kappa_a, kappa_b) (theta_a - theta_b), each kappa at its own node's T.
//   A node of an isothermal thermal contact has the row theta - theta_c. Insulator nodes conduct
//   (their material's kappa) and carry no carriers.
// - Temperature in the other rows: n_ie(T) (and its product) in SRH, Auger and radiative
//   recombination; the Caughey-Thomas mobility; an ohmic contact's neutral equilibrium and
//   potential psi0 = V / V_T + c_n + theta ln(n0 / Nc(theta)) at its own node's T.
// - Newton: a temperature step is clipped to 50 K (scaled by T0) and to [0.5, 2] times theta;
//   update_size includes |d theta|.
// The Jacobian is exact (with the model on, the node and edge blocks are full 4 x 4).
#pragma once

#include <complex>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/band_edges.hpp"
#include "NiTCAD/assemble/gate.hpp"
#include "NiTCAD/assemble/interface_edges.hpp"
#include "NiTCAD/assemble/models.hpp"
#include "NiTCAD/assemble/sg_flux.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/assemble/tunnel_paths.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"
#include "NiTCAD/physics/impact_ionization.hpp"
#include "NiTCAD/physics/ionization.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/statistics.hpp"

namespace NiTCAD::assemble {

class DriftDiffusion {
public:
    // Errors: those of EquilibriumPoisson::create (scaling temperature); invalid_input for
    // impact_current_resolution not finite and >= 0 with impact ionization on; for btbt_local
    // with btbt_nonlocal; for btbt_nonlocal without the device's cells, on a device with an edge
    // between two different semiconductors (no tunnel path crosses a heterointerface), or, for
    // direct_wkb, with a semiconductor region whose material lacks a direct gap or tunnelling
    // masses (silicon never has them).
    [[nodiscard]] static std::expected<DriftDiffusion, base::Error> create(
        const device::Device& device, const Scaling& scaling,
        const PhysicsModels& models = {});

    // The current density below which impact ionization fades out (see the header comment).
    static constexpr double impact_current_floor = 1e-12;  // [A/cm^2]

    [[nodiscard]] std::size_t unknowns() const noexcept { return m_ * node_count(); }
    [[nodiscard]] std::size_t node_count() const noexcept { return volume_.size(); }
    [[nodiscard]] std::size_t contact_count() const noexcept { return contact_count_; }
    // Unknowns per node: 3 (psi, n, p), or 4 with the electrothermal model (theta last).
    [[nodiscard]] std::size_t stride() const noexcept { return m_; }

    struct TimeStep;  // a BDF time step (below)

    // Electrothermal (Unit 23; see the header comment). Whether the model is on.
    [[nodiscard]] bool electrothermal() const noexcept { return m_ == 4; }
    // The temperature of each thermal contact in K (isothermal: its own; resistance: the ambient
    // behind it), in the order of device.thermal_contacts(); initially the device's values.
    // Errors (invalid_input, the index the thermal contact, the value it): the size is not the
    // thermal contact count; a value not finite and positive. (The solve layer checks a
    // temperature against the materials' range.) Precondition (NITCAD_EXPECTS): the model is on.
    [[nodiscard]] std::expected<void, base::Error> set_thermal_bias(
        std::span<const double> temperature_K);
    // Whether every connected part of the mesh (every edge conducts heat) has a node of a thermal
    // contact. A steady state needs one (DECISIONS.md T7); a transient does not. False without
    // the model.
    [[nodiscard]] bool heat_sinks_complete() const noexcept { return sinks_complete_; }
    // Scaled heat leaving through each thermal contact at state x (the heat its isothermal nodes'
    // rows would balance without their Dirichlet terms, or its R_th flux), in units of
    // V_T J0 L_D^(D-1) (a power in W / cm^(3-D) after multiplying by it), in
    // device.thermal_contacts() order. In a steady state their sum is electrical_power(x).
    // Preconditions (NITCAD_EXPECTS): the model is on; x has unknowns() entries.
    // With a time step, the heat rows' storage terms are in the balance.
    [[nodiscard]] std::vector<double> thermal_contact_heat(std::span<const double> x,
                                                           const TimeStep* step = nullptr) const;
    // The power the ohmic contacts deliver, sum_c V_c I_c, in the same units.
    [[nodiscard]] double electrical_power(std::span<const double> x) const;
    // The Joule part of each carrier edge's heat (Jn, Jp), in the same units: each flux times the
    // drop of its quasi-Fermi level along it, theta_e (ln(w_b / w_a) - delta_n) for electrons and
    // theta_e (ln(w_a / w_b) - delta_p) for holes; each >= 0, since a flux runs down that drop.
    // Zero on edges without carriers. Preconditions as thermal_contact_heat.
    [[nodiscard]] std::vector<std::pair<double, double>> edge_joule_heat(
        std::span<const double> x) const;

    // Applied bias of each contact in V, in the order of device.contacts(); sets the Dirichlet
    // values of the ohmic contacts and the electrode potential of the gates. Errors: those of
    // check_contact_bias (size, a value not finite). Initially every bias is 0.
    [[nodiscard]] std::expected<void, base::Error> set_bias(std::span<const double> bias_V);

    // Scaled state of every node in equilibrium at the scaled potential psi under the selected
    // statistics (Boltzmann: n = n_ie e^(psi + s), p = n_ie e^-(psi + s); 0 on insulator nodes),
    // with the ohmic contact and electrode nodes set to their Dirichlet values. Precondition
    // (NITCAD_EXPECTS): psi has node_count() entries.
    [[nodiscard]] std::vector<double> state_from_potential(std::span<const double> psi) const;
    // Sets the ohmic contact and electrode nodes of a state to their Dirichlet values and the
    // densities of insulator nodes to 0.
    void stamp_contacts(std::span<double> x) const;

    [[nodiscard]] linalg::SparseMatrix make_jacobian() const;
    // Preconditions (NITCAD_EXPECTS): x and residual have unknowns() entries; jacobian has the
    // pattern of make_jacobian().
    void evaluate(std::span<const double> x, std::span<double> residual,
                  linalg::SparseMatrix& jacobian) const;
    void residual(std::span<const double> x, std::span<double> residual) const;

    // A BDF time step (see the header comment): rate = 1 / (beta h), h in units of time_scale();
    // storage holds c_n and c_p of every node (storage_width() node_count() entries, node i from
    // storage_width() i; read only off the contacts and insulators), and with the electrothermal
    // model c of the rise tau third (read off the isothermal sinks); traps c of every trap slot
    // (trap_slots() entries).
    struct TimeStep {
        double rate;
        std::span<const double> storage;
        std::span<const double> traps;
    };
    // t0 = Ns / R0 [s].
    [[nodiscard]] double time_scale() const noexcept { return rate_scale_; }
    // Interface trap occupancies, one per (interface edge, level) (InterfaceEdges slots).
    [[nodiscard]] std::size_t trap_slots() const noexcept { return interfaces_.slot_count(); }
    // Preconditions (NITCAD_EXPECTS) as evaluate(), and step's spans of the sizes above.
    void evaluate(std::span<const double> x, const TimeStep& step, std::span<double> residual,
                  linalg::SparseMatrix& jacobian) const;
    void residual(std::span<const double> x, const TimeStep& step,
                  std::span<double> residual) const;
    // (S_n, S_p) of every node at state x, scaled (0 on insulator nodes), with the electrothermal
    // model (S_n, S_p, tau).
    [[nodiscard]] std::vector<double> storage(std::span<const double> x) const;
    [[nodiscard]] std::size_t storage_width() const noexcept { return m_ - 1; }
    // The occupancy f of every trap slot at state x: after `step`, or without one the steady
    // state.
    [[nodiscard]] std::vector<double> trap_occupancies(std::span<const double> x,
                                                       const TimeStep* step = nullptr) const;
    // Scaled charge of every contact at state x, in device.contacts() order: the displacement flux
    // leaving its nodes along the edges to nodes outside it (a gate: as gate_charges). Unlike
    // gate_charges it is not zero for an ohmic contact. Units as gate_charges.
    [[nodiscard]] std::vector<double> contact_charges(std::span<const double> x,
                                                      const TimeStep* step = nullptr) const;
    // Scaled conduction current entering through each contact after `step`: terminal_currents
    // plus, on an interface edge whose semiconductor node is on an ohmic contact, the current
    // that charges its traps from the contact.
    [[nodiscard]] std::vector<double> conduction_currents(std::span<const double> x,
                                                          const TimeStep& step) const;

    // Small-signal analysis (Unit 22; see the header comment). s is the Laplace variable in units
    // of 1 / time_scale() (s = i omega t0 for a frequency omega).
    // The pattern of small_signal_matrix: make_jacobian()'s, complex.
    [[nodiscard]] linalg::ComplexSparseMatrix make_small_signal_matrix() const;
    // J + s C + T(s) at a steady state x. Preconditions (NITCAD_EXPECTS): x has unknowns() entries;
    // a has the pattern of make_small_signal_matrix().
    void small_signal_matrix(std::span<const double> x, std::complex<double> s,
                             linalg::ComplexSparseMatrix& a) const;
    // dF / dV of contact c, per volt (the residual's dependence on its bias). Precondition
    // (NITCAD_EXPECTS): c < contact_count().
    [[nodiscard]] std::vector<double> bias_derivative(std::size_t contact) const;
    // With the electrothermal model the heat row of an ohmic node also depends on its contact's
    // bias, through the metal's Peltier term m (Jn + Jp) / K0, so dF / dV depends on the state:
    // this overload evaluates it at x (without the model it is bias_derivative(contact)). The
    // rows where it can be nonzero at any state: bias_rows(contact). F stays linear in V at fixed
    // x. Preconditions (NITCAD_EXPECTS): contact < contact_count(); x has unknowns() entries.
    [[nodiscard]] std::vector<double> bias_derivative(std::size_t contact,
                                                      std::span<const double> x) const;
    [[nodiscard]] std::vector<std::size_t> bias_rows(std::size_t contact) const;
    // The small-signal total current entering through a contact, conduction plus s times its
    // charge (conduction_currents and contact_charges, linearized about the steady state x):
    //     dI = sum_k values[k] dx[columns[k]] + bias dV,
    // dV the contact's own bias in V (a gate's charge depends on it directly). Scaled as
    // terminal_currents; one row per contact, in device.contacts() order.
    struct CurrentRow {
        std::vector<std::size_t> columns;
        std::vector<std::complex<double>> values;
        std::complex<double> bias;
    };
    [[nodiscard]] std::vector<CurrentRow> small_signal_currents(std::span<const double> x,
                                                                std::complex<double> s) const;

    // Nonlocal tunnelling (Unit 20; see the header comment). Whether the model is on; the paths
    // traced at state x (tunnel_paths.hpp) whose rate at x is at least `matters` of the largest
    // (the others would only fill the Jacobian); the frozen set the rows use (empty at creation),
    // whose setting rebuilds the Jacobian pattern (make_jacobian and make_small_signal_matrix again);
    // and each path's state at x. Preconditions (NITCAD_EXPECTS): x has unknowns() entries; the
    // paths of set_paths are trace_paths's of this system.
    [[nodiscard]] bool tunnelling() const noexcept {
        return tunnel_kind_ != NonlocalTunnelling::off;
    }
    [[nodiscard]] TunnelPaths trace_paths(std::span<const double> x) const;
    void set_paths(TunnelPaths paths);
    [[nodiscard]] const TunnelPaths& paths() const noexcept { return paths_; }
    [[nodiscard]] std::vector<TunnelPathState> path_states(std::span<const double> x) const;
    // Whether two path sets of this system act alike at state x: every path whose rate at x is at
    // least 1e-10 of the largest rate in either set has a path from the same start in the other
    // with the same geometry (tunnel_paths.hpp same_geometry). Paths below that are left out: a
    // start at the edge of the start criterion, or a crossing near the longest length, flips in
    // and out of the set with the rounding of the state, at a rate near exp(-690) of the scale.
    [[nodiscard]] bool paths_agree(const TunnelPaths& a, const TunnelPaths& b,
                                   std::span<const double> x) const;
    // The share of the largest rate below which a path does not matter (paths_agree), and below
    // which trace_paths drops it.
    static constexpr double matters = 1e-10;

    // Newton hooks (see the header comment).
    [[nodiscard]] double update_size(std::span<const double> x,
                                     std::span<const double> dx) const;
    void apply_update(std::span<double> x, std::span<const double> dx, double max_update) const;

    // Scaled (Jn, Jp) of every mesh edge, in units of J0 times the edge's scaled coupling area
    // (so a 1D edge gives the current density in units of J0).
    [[nodiscard]] std::vector<std::pair<double, double>> edge_currents(
        std::span<const double> x) const;
    // Scaled total current entering the device through each contact (conventional current from
    // the contact into the device), in device.contacts() order; zero for a gate. Physical value:
    // times J0 L_D^(D-1), in A / cm^(3-D) (A/cm^2 in 1D, A/cm in 2D, A in 3D).
    [[nodiscard]] std::vector<double> terminal_currents(std::span<const double> x) const;
    // Per ohmic contact, a bound on how far the current through any cut of the device can differ
    // from terminal_currents (same scale): the sum over the nodes of the total-current residual
    // |F_n + F_p| plus 8 eps times the magnitudes of the terms that cancel in it (the one-sided
    // flux terms, density times its flux coefficient, and the recombination terms). The same
    // value for every ohmic contact; zero for a gate. A current below it is not resolved by the
    // (psi, n, p) state.
    [[nodiscard]] std::vector<double> terminal_current_resolution(std::span<const double> x) const;
    // The band diagram at state x (band_edges.hpp); NaN on insulator nodes. Precondition
    // (NITCAD_EXPECTS): x has unknowns() entries and positive densities off the insulators.
    [[nodiscard]] BandEdges band_edges(std::span<const double> x) const;
    // Scaled charge on each gate, the sum over its nodes of G_i (psi_G,i - psi_i) (gate.hpp), and
    // on each electrode, the displacement flux leaving its nodes along the edges to other nodes,
    // in device.contacts() order; zero for an ohmic contact. Physical value: times q Ns L_D^D, in
    // C / cm^(3-D).
    [[nodiscard]] std::vector<double> gate_charges(std::span<const double> x,
                                                   const TimeStep* step = nullptr) const;
    // Scaled trapped charge of each declared interface at state x (interface_edges.hpp, without the
    // fixed charge), in units of q Ns L_D^D. Precondition (NITCAD_EXPECTS): x has unknowns()
    // entries.
    [[nodiscard]] std::vector<double> interface_trap_charges(
        std::span<const double> x, const TimeStep* step = nullptr) const;

private:
    DriftDiffusion() = default;

    struct EdgeTerm {
        std::size_t a, b;   // end nodes
        double c;           // et * scaled coupling, Poisson
        double an, ap;      // electron and hole SG edge factors (low field)
        double mu_n, mu_p;  // low-field edge mobilities (harmonic means) [cm^2/(V s)]
        double field;       // [V/cm] per unit |psi_b - psi_a|: V_T / length
        physics::CanaliParameters sat_n, sat_p;  // of node a's material
        // An edge between two materials: the ends' low-field mobilities and node b's Canali
        // parameters (field mobility takes the harmonic mean of the ends' Canali mobilities).
        bool mixed;
        double mu_n_a, mu_n_b, mu_p_a, mu_p_b;
        physics::CanaliParameters sat_n_b, sat_p_b;
        double shift_n;     // s_b - s_a + ln(n_ie,b / n_ie,a): delta_n = psi_b - psi_a + shift_n
        double shift_p;     // s_b - s_a - ln(n_ie,b / n_ie,a)
        // Thermionic emission (an interface edge with the model on): factors K and, per band,
        // ln(N_b / N_a) and N_a / N_b.
        bool thermionic;
        bool carriers;      // both ends are semiconductor nodes (else no carrier flux)
        bool charged;       // an interface edge: its Poisson flux is the interface's
        double te_kn, te_kp, te_log_nc, te_ratio_nc, te_log_nv, te_ratio_nv;
        // Positions in the Jacobian values: row a with columns of b, row b with columns of a,
        // in the order (psi, psi), (n, psi), (n, n), (p, psi), (p, p).
        std::size_t ab[5], ba[5];
    };

    // The degeneracy factors of a node's electron and hole densities (Fermi-Dirac).
    struct NodeDegeneracy {
        physics::Degeneracy n, p;
    };

    // Without `interfaces` the interface edges' terms are left out (small_signal_matrix adds its
    // own).
    void assemble(std::span<const double> x, std::span<double> residual,
                  std::span<double> values, const TimeStep* step, bool interfaces = true) const;

    // The storage of a semiconductor node and its derivatives d S_n / dn, d S_p / dp.
    struct NodeStorage {
        double n, d_n, p, d_p;
    };
    [[nodiscard]] NodeStorage node_storage(std::size_t i, std::span<const double> x,
                                           std::span<const NodeDegeneracy> g) const;
    // The trap step of interface edge k, or nullptr without a step.
    [[nodiscard]] std::optional<TrapStep> trap_step(std::size_t k, const TimeStep* step) const;
    // The interface of edge k at state x.
    [[nodiscard]] InterfaceDrift interface_at(std::size_t k, std::span<const double> x,
                                              const TimeStep* step) const;

    // Impact ionization (Unit 19; see the header comment). Per node that ionizes: its carrier
    // edges (in impact_edges_) and the pseudo-inverse of sum t t^T over them (row-major, D x D);
    // per such edge: the unit vector t from its first to its second node, the field per unit
    // psi_b - psi_a and the current density per unit scaled flux, and the Jacobian positions of
    // the node's electron row with the other end's p column and of its hole row with the other
    // end's n column. Empty when the model is off.
    struct ImpactNode {
        std::size_t node;
        std::size_t first, last;
        double inverse[9];
        physics::ImpactIonizationCoefficients n, p;
        double gamma;  // the temperature factor of the node's material
    };
    struct ImpactEdge {
        std::size_t edge;
        double t[3];
        double field;    // [V/cm] per unit (psi_b - psi_a): V_T / length
        double current;  // [A/cm^2] per unit scaled flux: J0 / scaled coupling area
        std::size_t n_p, p_n;
    };
    void add_impact_generation(std::span<const double> x, std::span<const NodeDegeneracy> g,
                               std::span<double> f, std::span<double> values) const;

    // Local band-to-band tunnelling (Unit 20): per node, its carrier edges (in field_edges_) and
    // the pseudo-inverse of sum t t^T, as ImpactNode, and its material's Kane pair.
    struct TunnelNode {
        std::size_t node;
        std::size_t first, last;
        double inverse[9];
        double A, B;
    };
    struct FieldEdge {
        std::size_t edge;
        double t[3];
        double field;  // [V/cm] per unit (psi_b - psi_a)
    };
    void add_local_tunnelling(std::span<const double> x, std::span<double> f,
                              std::span<double> values) const;

    // Nonlocal tunnelling (Unit 20): per frozen path, the nodes whose psi its rows read (columns),
    // the nodes its electrons go to (deposit), the start's hole-row positions per column and the
    // deposit rows' positions (deposit-major), and the forward length before each segment [cm].
    struct PathLayout {
        std::vector<std::size_t> columns, deposit;
        std::vector<std::size_t> hole, electron;
        std::vector<double> before_cm;
    };
    // One path's live terms at x: its pairs P (scaled) and dP / dpsi per column, the share of each
    // deposit node and its partials (deposit-major), and the diagnostics.
    struct PathTerms {
        double pairs = 0.0;
        std::vector<double> d_pairs;
        std::vector<double> share, d_share;
        TunnelPathState state{};
    };
    // The columns, deposit nodes and lengths of a path (no Jacobian positions).
    [[nodiscard]] PathLayout layout(const TunnelPath& p) const;
    void path_terms(const TunnelPath& p, const PathLayout& l, std::span<const double> x,
                    PathTerms& t) const;
    void add_path_tunnelling(std::span<const double> x, std::span<double> f,
                             std::span<double> values) const;
    // The Jacobian positions every term uses, after the pattern changed (set_paths).
    void locate_positions();

    // Per node at state x under Fermi-Dirac statistics; empty under Boltzmann.
    [[nodiscard]] std::vector<NodeDegeneracy> degeneracies(std::span<const double> x) const;

    // The electron and hole fluxes of an edge at state x, with field-dependent edge factors when
    // the model is on and the degeneracy terms g (degeneracies(x)) when it is not empty.
    [[nodiscard]] std::pair<EdgeFlux, EdgeFlux> edge_fluxes(
        const EdgeTerm& e, std::span<const double> x, std::span<const NodeDegeneracy> g) const;

    // Electrothermal (Unit 23; electrothermal.cpp). Per node, constant: the material (index into
    // thermal_materials_, or -1 on an insulator), its thermal data, the impurity densities the
    // mobility and band-gap narrowing read, the band constants c_n, c_p and Eg(T0), the
    // thermopower exponents, ln Nc0 and ln Nv0 (scaled), the thermal contact (-1: none), the
    // R_th coefficient h (0: none) and the heat capacity C = rho c V L_D^2 / (t0 kappa_ref) of the
    // heat row (times the rate and the change of tau).
    struct ThermalNode {
        std::int32_t material;
        physics::ThermalParameters thermal;
        double impurity, impurity_mobility, narrowing_eV;
        double c_n, c_p, gap0_eV, r_n, r_p, log_nc, log_nv;
        std::int32_t contact;
        double h;
        double capacity;
    };
    // Per node at a state: theta and tau = theta - 1 (the unknown), n_ie / Ns, the mobilities
    // [cm^2/(V s)], g (the valence band's shift), kappa [W/(cm K)] and the Kirchhoff transform u,
    // each with d / d theta.
    struct ThermalState {
        double theta, tau, nie, d_nie, mu_n, d_mu_n, mu_p, d_mu_p, gap, d_gap, kappa, d_kappa, u;
    };
    // Per edge, constant: V_T / D0 times the coupling (g_e), L_c, L_v, the mean thermopower
    // exponents, whether conduction takes the Kirchhoff form (both ends in one material), g_e /
    // kappa_ref, and the Jacobian positions of the full 4 x 4 blocks (row a with columns of b,
    // row b with columns of a; row-major).
    struct ThermalEdge {
        double coupling, log_nc, log_nv, r_n, r_p;
        bool kirchhoff;
        double conduction;
        std::size_t ab[16], ba[16];
    };
    // A flux and its partials with respect to psi, the density and theta of each end (a, b).
    struct ThermalFlux {
        double flux;
        double d_psi[2], d_c[2], d_theta[2];
        double delta, theta_e, w[2];  // for the Joule part
    };
    [[nodiscard]] std::vector<ThermalState> thermal_states(std::span<const double> x) const;
    // A carrier's band reduced energy x = ln(c / N(theta)) - ln gamma (Fermi-Dirac; ln gamma 0
    // under Boltzmann statistics) at its node's T, with ln gamma and the Fermi-Dirac thermal
    // diffusion factor's excess over the Boltzmann 1 + r (r = -1/2): h - 1/2, each with its
    // partials in the density c and theta. log_N: ln N0 (scaled) of the band.
    struct Level {
        double x, dx_dc, dx_dt;
        double log_gamma, dg_dc, dg_dt;
        double extra, dh_dc, dh_dt;
    };
    [[nodiscard]] Level level(double density, double theta, double nie, double log_N) const;
    struct NodeLevels {
        Level n, p;
    };
    // Per node at x under Fermi-Dirac statistics; empty under Boltzmann.
    [[nodiscard]] std::vector<NodeLevels> thermal_levels(std::span<const double> x,
                                                         std::span<const ThermalState> t) const;
    [[nodiscard]] std::pair<ThermalFlux, ThermalFlux> thermal_fluxes(
        std::size_t edge, std::span<const double> x, std::span<const ThermalState> t,
        std::span<const NodeLevels> g) const;
    // The heat row's conduction term along an edge, from a to b, and its theta partials.
    struct Conduction {
        double flux, d_a, d_b;
    };
    [[nodiscard]] Conduction conduction(std::size_t edge, std::span<const ThermalState> t) const;
    // An ohmic contact node's Dirichlet values at its T and their theta derivatives.
    struct ThermalOhmic {
        double psi, n, p, d_psi, d_n, d_p;
    };
    [[nodiscard]] ThermalOhmic thermal_ohmic(std::size_t node, const ThermalState& t) const;
    // Creation of the thermal data; errors as create's (see electrothermal.cpp).
    // edge_geometry: each mesh edge's scaled coupling (coupling area / length / L_D^(D-2)).
    [[nodiscard]] std::optional<base::Error> make_thermal(const device::Device& device,
                                                          const Scaling& scaling,
                                                          const PhysicsModels& models,
                                                          std::span<const double> edge_geometry);

    std::size_t m_ = 3;  // unknowns per node
    std::vector<physics::Semiconductor> thermal_materials_;
    std::vector<ThermalNode> thermal_nodes_;
    std::vector<ThermalEdge> thermal_edges_;
    std::vector<double> thermal_rise_;  // per thermal contact: (T_c - T0) / T0
    std::vector<double> metal_;          // per node: V / V_T on an ohmic contact node, else 0
    std::vector<double> contact_bias_;   // per contact: V / V_T
    // A thermal contact's kind (true: isothermal) and its nodes' indices into its list.
    std::vector<char> thermal_isothermal_;
    double T0_ = 0.0, heat_scale_ = 0.0, kappa_ref_ = 0.0;
    bool thermal_bgn_ = false;
    bool sinks_complete_ = false;
    // One node's thermal state at tau (NaN entries unless tau is finite and theta positive).
    [[nodiscard]] ThermalState thermal_state(std::size_t node, double tau) const;
    // s times the Jacobian of the electrothermal storage terms (assemble_heat with a step) about a
    // steady state, added to the small-signal matrix's values.
    void add_thermal_storage(std::span<const double> x, std::complex<double> s,
                             std::span<std::complex<double>> values) const;
    // The electrothermal rows (assemble with the model on); with sinks false, without the thermal
    // contacts' terms, so every heat row is the node's heat balance alone.
    void assemble_heat(std::span<const double> x, std::span<double> f, std::span<double> values,
                       const TimeStep* step, bool sinks) const;

    std::vector<double> volume_, doping_, n_ie_, tau_n_, tau_p_, auger_n_, auger_p_;
    std::vector<double> log_dos_n_, log_dos_p_;  // ln(Nc / n_ie), ln(Nv / n_ie)
    std::vector<double> band_shift_;             // s
    std::vector<double> donors_, acceptors_;     // N_D / Ns, N_A / Ns
    std::vector<physics::DopantLevels> levels_;  // dopant levels in units of kT
    std::vector<double> radiative_;              // B [cm^3/s]
    std::vector<std::int32_t> contact_;   // ohmic contact index per node, or -1
    std::vector<char> insulator_;         // 1 on insulator nodes
    std::vector<std::int32_t> electrode_; // electrode contact index per node, or -1
    std::vector<double> electrode_potential_;  // psi_E at zero bias
    GateNodes gates_;
    InterfaceEdges interfaces_;
    // Per interface edge, the Jacobian positions of the insulator node's Poisson row with the
    // semiconductor node's n and p columns.
    std::vector<std::pair<std::size_t, std::size_t>> interface_np_;

    [[nodiscard]] InterfaceStatistics statistics(std::size_t node) const noexcept;
    std::vector<device::ContactKind> kinds_;  // per contact
    std::vector<double> psi0_, n0_, p0_;  // Dirichlet values per node (contact nodes only)
    std::vector<EdgeTerm> edges_;
    std::vector<std::size_t> block_;      // per node, 9 positions of its 3x3 diagonal block
    std::vector<ImpactNode> impact_nodes_;
    std::vector<ImpactEdge> impact_edges_;
    std::vector<TunnelNode> tunnel_nodes_;
    std::vector<FieldEdge> field_edges_;
    // Nonlocal tunnelling: the model, the device's cells, per node the longest path that may start
    // there (0: none) [cm], the material's gap at T [eV] and tunnelling masses [kg] (direct_wkb),
    // and the Kane pair (kane); the frozen paths, their layouts, and the pattern's triplets
    // without them.
    NonlocalTunnelling tunnel_kind_ = NonlocalTunnelling::off;
    std::optional<mesh::TensorCells> cells_;
    std::vector<double> tunnel_length_, tunnel_gap_, tunnel_mc_, tunnel_mv_, tunnel_A_, tunnel_B_;
    TunnelPaths paths_;
    std::vector<PathLayout> layouts_;
    std::vector<linalg::Triplet> base_triplets_;
    bool generation_ = false;  // any generation model on
    int dimension_ = 1;
    double generation_scale_ = 0.0;        // 1 / R0: [cm^-3 s^-1] to the rows' scaled rate
    double impact_resolution_ = 0.0;       // models.impact_current_resolution
    std::size_t contact_count_ = 0;
    double rate_scale_ = 0.0;             // Ns / R0: scaled SRH call to R / R0
    double Ns_ = 0.0;
    double V_T_ = 0.0;
    bool srh_ = true;
    bool auger_ = true;
    bool field_mobility_ = false;
    bool fermi_dirac_ = false;
    bool ionization_ = false;
    bool radiative_on_ = true;
    linalg::SparseMatrix pattern_;
};

}  // namespace NiTCAD::assemble
