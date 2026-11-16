// The generic equation system (Unit 57, DECISIONS.md GE3): a composition of equation terms over a
// map of the unknowns (fields.hpp DofMap), with only the members every equation set has. The
// model-specific queries and mutators are free functions in the topic headers: currents.hpp,
// thermal.hpp, tunnelling.hpp, energy.hpp, charges.hpp and rates.hpp, and band_edges.hpp's
// band_edges. The equation sets of a system created from a device and models are those of the
// coupled drift-diffusion system below.
//
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
#include <expected>
#include <memory>
#include <span>
#include <vector>

#include "NiTCAD/assemble/fields.hpp"
#include "NiTCAD/assemble/models.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/linalg/sparse_matrix.hpp"

namespace NiTCAD::assemble {

struct SystemData;  // the model data and the state (src/assemble/system_data.hpp, Unit 57 GE2)
struct SystemState;
class SystemKey;  // constructible only inside src/assemble (the passkey of data() and state())

class System {
public:
    // Errors: invalid_input for impact_current_resolution not finite and >= 0 with impact
    // ionization on; for btbt_local with btbt_nonlocal; for btbt_nonlocal without the device's
    // cells, on a device with an edge between two different semiconductors (no tunnel path
    // crosses a heterointerface), or, for direct_wkb, with a semiconductor region whose material
    // lacks a direct gap or tunnelling masses (silicon never has them). With models.electrothermal
    // (invalid_input): band-to-band tunnelling or interface traps (DECISIONS.md T8), a
    // thermionic-emission interface (U3, deferred), incomplete ionization (U6), Fermi-Dirac
    // statistics with a thermopower exponent other than -1/2 (T3), a region whose material has
    // no thermal data (T11). Those of the scaled device (scaling temperature).
    [[nodiscard]] static std::expected<System, base::Error> create(
        const device::Device& device, const Scaling& scaling,
        const PhysicsModels& models = {});

    [[nodiscard]] std::size_t unknowns() const noexcept;
    [[nodiscard]] std::size_t node_count() const noexcept;
    [[nodiscard]] std::size_t contact_count() const noexcept;
    // Unknowns per node: 3 (psi, n, p), 4 with the electrothermal model (theta last), 5 with
    // energy balance (tau_n, tau_p last).
    [[nodiscard]] std::size_t stride() const noexcept;
    // Where each field's unknowns are (fields.hpp; Unit 40): node-interleaved, stride() per node.
    // Precondition (NITCAD_EXPECTS): not a moved-from system.
    [[nodiscard]] const DofMap& dofs() const noexcept;

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

    // Small-signal analysis (Unit 22; see the header comment). s is the Laplace variable in units
    // of 1 / time_scale() (s = i omega t0 for a frequency omega).
    // The pattern of small_signal_matrix: make_jacobian()'s, complex.
    [[nodiscard]] linalg::ComplexSparseMatrix make_small_signal_matrix() const;
    // J + s C + T(s) at a steady state x: the steady terms but those a frequency-domain term
    // replaces, the storage at s, then the terms active only in the frequency domain (the
    // interface edges' traps with their occupancy's dynamics). Preconditions (NITCAD_EXPECTS):
    // x has unknowns() entries; a has the pattern of make_small_signal_matrix().
    void small_signal_matrix(std::span<const double> x, std::complex<double> s,
                             linalg::ComplexSparseMatrix& a) const;
    // dF / dV of contact c, per volt (the residual's dependence on its bias). Precondition
    // (NITCAD_EXPECTS): c < contact_count().
    [[nodiscard]] std::vector<double> bias_derivative(std::size_t contact) const;
    // Where dF / dV depends on the state (bias_nonlinear: with the electrothermal model the heat
    // row of an ohmic node depends on its contact's bias through the metal's Peltier term
    // m (Jn + Jp) / K0; with gate tunnelling through a lumped gate, whose tunnelling current
    // depends on the gate's bias exponentially), this overload evaluates it at x (otherwise it is
    // bias_derivative(contact)). The rows where it can be nonzero at any state:
    // bias_rows(contact), which then include the tunnelling nodes' continuity rows. With the
    // electrothermal model F stays linear in V at fixed x. Preconditions (NITCAD_EXPECTS):
    // contact < contact_count(); x has unknowns() entries.
    [[nodiscard]] std::vector<double> bias_derivative(std::size_t contact,
                                                      std::span<const double> x) const;
    [[nodiscard]] std::vector<std::size_t> bias_rows(std::size_t contact) const;
    // Whether dF / dV depends on the state (bias_derivative's x overload must be used).
    [[nodiscard]] bool bias_nonlinear() const noexcept;

    // Newton hooks (see the header comment; the minority floor: set_minority_floor below).
    [[nodiscard]] double update_size(std::span<const double> x,
                                     std::span<const double> dx) const;
    void apply_update(std::span<double> x, std::span<const double> dx, double max_update) const;

    // The internals, for src/assemble alone (the key cannot be made outside it): the model data
    // the copies share and the state each owns.
    [[nodiscard]] const SystemData& data(const SystemKey& key) const noexcept;
    [[nodiscard]] const SystemState& state(const SystemKey& key) const noexcept;
    [[nodiscard]] SystemState& state(const SystemKey& key) noexcept;
    // Replaces the shared data by a changed copy (copy on write: set_paths moves the positions).
    void replace_data(const SystemKey& key, std::shared_ptr<const SystemData> data) noexcept;

    // A copy shares the model data and owns its state (Unit 57 GE2).
    System(const System& other);
    System& operator=(const System& other);
    System(System&& other) noexcept;
    System& operator=(System&& other) noexcept;
    ~System();

private:
    System();

    void assemble(std::span<const double> x, std::span<double> residual,
                  std::span<double> values, const TimeStep* step, bool frequency = false) const;

    std::shared_ptr<const SystemData> data_;  // fixed at creation, shared by copies
    std::unique_ptr<SystemState> state_;      // owned: biases, paths, pattern
};

// The minority floor (Unit 50, DECISIONS.md RD2; off by default): with it on, update_size
// measures every density against at least eps of the larger density at its node.
void set_minority_floor(System& system, bool on) noexcept;
[[nodiscard]] bool minority_floor(const System& system) noexcept;

}  // namespace NiTCAD::assemble
