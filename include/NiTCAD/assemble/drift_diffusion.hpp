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
// Surface mobility (Unit 46, models.surface_mobility = lombardi; DECISIONS.md LD1-LD8): each
// carrier edge's ends take Lombardi's mobility (physics/surface_mobility.hpp), the Matthiessen
// sum of the end's bulk mobility and the surface terms at its N and T in the edge's E_perp, before
// the harmonic mean and Canali. E_perp = sqrt(|v|^2 + E_0^2), E_0 = 1 V/cm, with v the part of the
// edge's field normal to the edge: the mean of the ends' fields, each reconstructed by least
// squares from its carrier edges as for impact ionization, less its component along the edge.
// v is linear in psi over both ends' stencils (SurfaceStencil), so the Jacobian is exact: the
// carrier rows of a and b gain the psi columns of those nodes (with the electrothermal model the
// heat rows too); with the model off nothing changes. In 1D v vanishes and E_perp = E_0.
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
// impact ionization, Unit 19). With set_minority_floor on (Unit 50, DECISIONS.md RD2; off by
// default) every density is also measured against at least eps of the larger density at its node:
// a deep minority's correction is rounding there.
//
// Electrothermal coupling (Unit 23, models.electrothermal; DECISIONS.md T1-T14; legacy thermal.py,
// thermal_grid.py). Each node has a fourth unknown, the scaled temperature rise
// tau = (T - T0) / T0 (T0 the device's temperature, the scaling's), interleaved as x[4i + 3];
// theta = T / T0 = 1 + tau below. The rise, not theta, is the unknown (DECISIONS.md U1) so that
// a small heating keeps its digits: a double near 1 is spaced 2.2e-16, so theta resolves T only
// to 6.7e-14 K at 300 K, while tau resolves the rise to 2.2e-16 of itself. The heat reaching a
// sink is the flux kappa dT / h through the edge next to it, so with theta it was resolved only
// to kappa 6.7e-14 K / h, 4e-8 W/cm^2 for h = 25 nm. Measured on a diode at 3.7e-6 W/cm^2
// (0.3 V): with theta the energy balance missed by 9.0e-9 W/cm^2 (0.24%); with tau by 2.4e-12
// (6.4e-7), the rounding floor of the energy flows (U8). Without the model the stride stays 3 and
// nothing changes. Energies below are in units of k T0, densities and fluxes as above.
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
//   equilibrium; between two temperatures the open-circuit voltage is the Seebeck voltage, to
//   second order in the spacing (DECISIONS.md U2). Thermionic-emission interfaces and incomplete
//   ionization are refused with the model (U3 deferred, U6). A material with density-of-states
//   mass fits (Unit 35, DECISIONS.md S4) has N(T) = N(T0) theta^(3/2) M(theta): the extra
//   -k T n grad ln M is a band-parameter step, so delta_n gains ln M_c(theta_b) - ln M_c(theta_a)
//   beside L_c and delta_p loses ln M_v(theta_b) - ln M_v(theta_a) beside -L_v, and a density at
//   fixed reduced energy moves with 3 / (2 theta) + d ln M / d theta; the energy a carrier carries
//   is unchanged (parabolic bands at the local mass, an assumption). M = 1 without a fit.
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
#include <memory>
#include <optional>
#include <span>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/band_edges.hpp"
#include "NiTCAD/assemble/fields.hpp"
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
#include "NiTCAD/physics/surface_mobility.hpp"

namespace NiTCAD::assemble {

struct SystemData;   // the model data and the state (src/assemble/system_data.hpp, Unit 57 GE2)
struct SystemState;

class DriftDiffusion {
public:
    // Errors: those of EquilibriumPoisson::create (scaling temperature); invalid_input for
    // impact_current_resolution not finite and >= 0 with impact ionization on; for btbt_local
    // with btbt_nonlocal; for btbt_nonlocal without the device's cells, on a device with an edge
    // between two different semiconductors (no tunnel path crosses a heterointerface), or, for
    // direct_wkb, with a semiconductor region whose material lacks a direct gap or tunnelling
    // masses (silicon never has them). With models.electrothermal (invalid_input): band-to-band
    // tunnelling or interface traps (DECISIONS.md T8), a thermionic-emission interface (U3,
    // deferred), incomplete ionization (U6), Fermi-Dirac statistics with a thermopower exponent
    // other than -1/2 (T3), a region whose material has no thermal data (T11).
    [[nodiscard]] static std::expected<DriftDiffusion, base::Error> create(
        const device::Device& device, const Scaling& scaling,
        const PhysicsModels& models = {});

    // The current density below which impact ionization fades out (see the header comment).
    static constexpr double impact_current_floor = 1e-12;  // [A/cm^2]

    [[nodiscard]] std::size_t unknowns() const noexcept;
    [[nodiscard]] std::size_t node_count() const noexcept;
    [[nodiscard]] std::size_t contact_count() const noexcept;
    // Unknowns per node: 3 (psi, n, p), or 4 with the electrothermal model (theta last).
    [[nodiscard]] std::size_t stride() const noexcept;
    // Where each field's unknowns are (fields.hpp; Unit 40): node-interleaved, stride() per node.
    // Precondition (NITCAD_EXPECTS): not a moved-from system.
    [[nodiscard]] const DofMap& dofs() const noexcept;

    struct TimeStep;  // a BDF time step (below)

    // Electrothermal (Unit 23; see the header comment). Whether the model is on.
    [[nodiscard]] bool electrothermal() const noexcept;
    // Energy balance (Unit 52, models.energy_balance): tau_n and tau_p are the fourth and fifth
    // unknowns of every node.
    [[nodiscard]] bool energy_balance() const noexcept;
    // The carriers' energy budget of state x (Unit 52, EC7 as ruled 2026-10-29), in the units of
    // electrical_power (k T_L per scaled particle flux): in a steady state
    //     input = relaxation + recombination + contact_flux + contact_joule - peltier + residual
    // to rounding, summed over the semiconductor nodes off the ohmic contacts:
    // - input: sum_c V_c I_c (electrical_power);
    // - relaxation: V k (3/2) (n tau_n / tau_wn + p tau_p / tau_wp), the power to the lattice;
    // - recombination: V k U (E_c - E_v + (3/2)(theta_n + theta_p)), each recombined pair's gap
    //   and thermal energy (U the continuity rows' net rate);
    // - contact_flux: the energy flux the edges carry into the ohmic contact nodes (where tau = 0);
    // - contact_joule: the half of the Joule heat of the contact edges no energy row keeps;
    // - peltier: sum over the contact nodes of the carriers' band-edge energy above their Fermi
    //   level times their flow into the device, (E_c + V) D_n + (E_v + V) D_p;
    // - residual: the Newton residual's power, sum of the energy rows' residuals and of E_c and
    //   E_v times the continuity rows' (zero for an exact solution).
    // After a time step (Unit 52, EC13) the balance is the step's: input is the conduction
    // power, recombination counts the rates alone, and two terms store energy:
    //     input = ... + residual + thermal_storage + carrier_storage
    // - thermal_storage: V rate (W_n + W_p - their history), the carriers' thermal energy (ED9);
    // - carrier_storage: V rate (E_c (S_n - c_n) - E_v (S_p - c_p)), the band-edge energy of the
    //   carriers stored (the continuity rows' storage terms, at the step's band edges).
    // Both are 0 without a step. Preconditions (NITCAD_EXPECTS): the model is on and the
    // temperatures are not held; x has unknowns() entries; step's spans of evaluate's sizes.
    struct EnergyPowers {
        double input, relaxation, recombination, contact_flux, contact_joule, peltier, residual;
        double thermal_storage = 0.0, carrier_storage = 0.0;
    };
    [[nodiscard]] EnergyPowers energy_powers(std::span<const double> x,
                                             const TimeStep* step = nullptr) const;
    // The two-stage start of a bias point (Unit 52, ED12 as ruled 2026-10-29): hold_temperatures
    // makes every energy row tau - tau_held = 0, tau_held the rises of state x (so psi, n and p are
    // solved at the given carrier temperatures), until release_temperatures. Precondition
    // (NITCAD_EXPECTS): the model is on; x has unknowns() entries.
    void hold_temperatures(std::span<const double> x);
    void release_temperatures() noexcept;
    [[nodiscard]] bool temperatures_held() const noexcept;
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
    [[nodiscard]] bool heat_sinks_complete() const noexcept;
    // Whether the connected thermal domain of a node has a node of a thermal contact (DECISIONS.md
    // U5): a region without a sink of its own but joined to one has. Preconditions
    // (NITCAD_EXPECTS): the model is on; node < node_count().
    [[nodiscard]] bool reaches_heat_sink(std::size_t node) const;
    // K0 = kappa_ref T0 / (V_T J0 L_D) (see the header comment): the heat rows are the nodes' heat
    // balances in units of V_T J0 L_D^(D-1) divided by it. 0 without the model.
    [[nodiscard]] double heat_scale() const noexcept;
    // Scaled heat leaving through each thermal contact at state x (the heat its isothermal nodes'
    // rows would balance without their Dirichlet terms, or its R_th flux), in units of
    // V_T J0 L_D^(D-1) (a power in W / cm^(3-D) after multiplying by it), in
    // device.thermal_contacts() order. In a steady state their sum is electrical_power(x).
    // Preconditions (NITCAD_EXPECTS): the model is on; x has unknowns() entries.
    // With a time step, the heat rows' storage terms are in the balance.
    [[nodiscard]] std::vector<double> thermal_contact_heat(std::span<const double> x,
                                                           const TimeStep* step = nullptr) const;
    // The power the ohmic contacts deliver, sum_c V_c I_c, in the same units (also with energy
    // balance, Unit 52).
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
    // model c of the rise tau third (read off the isothermal sinks), with energy balance c of the
    // carriers' thermal energies W_n and W_p third and fourth (Unit 52, ED9); traps c of every
    // trap slot (trap_slots() entries).
    // time_s is the step's end time [s], read only by the manufactured sources' time profile
    // (verification only, models.hpp; DECISIONS.md Unit 55 16.3, DT-Q3 A).
    struct TimeStep {
        double rate;
        std::span<const double> storage;
        std::span<const double> traps;
        double time_s = 0.0;
    };
    // t0 = Ns / R0 [s].
    [[nodiscard]] double time_scale() const noexcept;
    // Interface trap occupancies, one per (interface edge, level) (InterfaceEdges slots).
    [[nodiscard]] std::size_t trap_slots() const noexcept;
    // Preconditions (NITCAD_EXPECTS) as evaluate(), and step's spans of the sizes above.
    void evaluate(std::span<const double> x, const TimeStep& step, std::span<double> residual,
                  linalg::SparseMatrix& jacobian) const;
    void residual(std::span<const double> x, const TimeStep& step,
                  std::span<double> residual) const;
    // (S_n, S_p) of every node at state x, scaled (0 on insulator nodes), with the electrothermal
    // model (S_n, S_p, tau), with energy balance (S_n, S_p, W_n, W_p).
    [[nodiscard]] std::vector<double> storage(std::span<const double> x) const;
    [[nodiscard]] std::size_t storage_width() const noexcept;
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
    // Whether dF / dV depends on the state (bias_derivative's x overload must be used): with the
    // electrothermal model, and with gate tunnelling through a lumped gate (Unit 51), whose
    // tunnelling current depends on the gate's bias exponentially. The rows bias_rows names then
    // include the tunnelling nodes' continuity rows.
    [[nodiscard]] bool bias_nonlinear() const noexcept;
    // Gate tunnelling (Unit 51, GD10): whether any contact tunnels, and per contact the part of
    // its terminal current each process carries at state x (scaled as terminal_currents; zero off
    // a tunnelling contact). They sum to the tunnelling part of terminal_currents. Empty vectors
    // without gate tunnelling. Precondition (NITCAD_EXPECTS): x has unknowns() entries.
    [[nodiscard]] bool gate_tunnelling() const noexcept;
    struct GateTunnellingCurrents {
        std::vector<double> ecb, hvb, evb;
    };
    [[nodiscard]] GateTunnellingCurrents gate_tunnelling_currents(std::span<const double> x) const;
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
    [[nodiscard]] bool tunnelling() const noexcept;
    [[nodiscard]] TunnelPaths trace_paths(std::span<const double> x) const;
    void set_paths(TunnelPaths paths);
    [[nodiscard]] const TunnelPaths& paths() const noexcept;
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
    // The minority floor (Unit 50, DECISIONS.md RD2; off by default): with it on, update_size
    // measures every density against at least eps of the larger density at its node.
    void set_minority_floor(bool on) noexcept;
    [[nodiscard]] bool minority_floor() const noexcept;
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
    // Per mesh edge, the magnitudes of the one-sided terms whose difference is the edge's total
    // current: |d_c1 c_a| + |d_c2 c_b| of each carrier (density times its flux coefficient), in
    // the units of edge_currents; 0 for an edge without carriers. It scales the edge current's
    // rounding (8 eps times it), and half of it is the edge's conductance to a quasi-Fermi
    // potential difference near equilibrium (the weights of the weighted current, Unit 32).
    [[nodiscard]] std::vector<double> edge_flux_terms(std::span<const double> x) const;
    // The band diagram at state x (band_edges.hpp); NaN on insulator nodes. Precondition
    // (NITCAD_EXPECTS): x has unknowns() entries and positive densities off the insulators.
    [[nodiscard]] BandEdges band_edges(std::span<const double> x) const;

    // Unit 45 (DECISIONS.md RS1): each generation-recombination mechanism at state x, per node, in
    // cm^-3 s^-1, read from its own equation term: the term's contribution to the node's electron
    // and hole continuity rows divided by the node's control volume, so the rates are the
    // residual's, operand for operand. Sign: net recombination R - G (a generation is negative).
    // A mechanism is listed only when the system has its term; the local and nonlocal
    // band-to-band terms add into one entry. Interfaces: their edges' carrier exchange (surface
    // recombination and trap capture) over the node's control volume, an effective volume rate at
    // interface nodes; read with steady traps (no time step). Precondition (NITCAD_EXPECTS): x has
    // unknowns() entries.
    enum class RateMechanism : std::uint8_t { srh, auger, radiative, interfaces, impact, band_to_band };
    struct MechanismRate {
        RateMechanism mechanism;
        std::vector<double> electron;  // per node, as the electron row sees it
        std::vector<double> hole;
    };
    [[nodiscard]] std::vector<MechanismRate> rates(std::span<const double> x) const;
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


    // A copy shares the model data and owns its state (Unit 57 GE2).
    DriftDiffusion(const DriftDiffusion& other);
    DriftDiffusion& operator=(const DriftDiffusion& other);
    DriftDiffusion(DriftDiffusion&& other) noexcept;
    DriftDiffusion& operator=(DriftDiffusion&& other) noexcept;
    ~DriftDiffusion();

private:
    DriftDiffusion();

    // The terms, their composition and helpers are internal to the assemble library
    // (src/assemble/system_data.hpp, terms.hpp, composition.hpp); SystemData reads these two.
    friend struct SystemData;

    // Without `interfaces` the interface edges' terms are left out (small_signal_matrix adds its
    // own).
    void assemble(std::span<const double> x, std::span<double> residual,
                  std::span<double> values, const TimeStep* step, bool interfaces = true) const;

    std::shared_ptr<const SystemData> data_;  // fixed at creation, shared by copies
    std::unique_ptr<SystemState> state_;      // owned: biases, paths, pattern
};

}  // namespace NiTCAD::assemble
