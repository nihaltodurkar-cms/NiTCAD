<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 51 — gate tunnelling (started 2026-10-28; on `physics/gate-tunnelling`)

**Owner request:** "start GATE tunneling" (2026-10-28). ARCHITECTURE.md 11.3 lists gate tunnelling as Planned (the
nanoscale gap found by the FinFET deck). Numbered 51 because Units 47 (PR #66) and 50 (PR #71) are not merged; this
branch starts from `main` (`7a11d26`).
**Owner rulings (2026-10-28):** (1) both gate kinds tunnel: a lumped `gate` locally through its t_ox, and an
`electrode` on a meshed insulator along a straight path through the meshed stack (tensor grids); (2) self-consistent:
a boundary flux in the continuity rows with an exact Jacobian, counted in the gate's terminal current; (3) the
Tsu–Esaki energy integral with a WKB transmission through the trapezoidal or piecewise multi-layer barrier; (4) three
processes: conduction-band electrons (ECB), valence-band holes (HVB) and valence-band electrons (EVB); (5) SiO₂ gets
built-in, cited band data; any other insulator must give its band data in the deck, else the tunnelling gate is
refused; (6) validation against closed forms, plus a measured J–V comparison only if an accessible source turns up;
(7) excluded and refused: electrothermal, the weighted terminal current (Unit 32), electrode tunnelling on an
unstructured mesh; not built: image-force lowering. The standing conditions hold: no change to any device without
tunnelling switched on, criteria fixed before code, Release and Debug suites at every commit, `/W4 /WX`, hand
mutants, nothing pushed or merged without authorization.

**Sources read (before coding).**
- A. Gehring, "Simulation of Tunneling in Semiconductor Devices", dissertation, TU Wien, 2003 (open access,
  iue.tuwien.ac.at/phd/gehring): §3.2 (the Tsu–Esaki formula J = (4π m q / h³) ∫ TC(E_x) N(E_x) dE_x, parabolic
  bands, conserved parallel momentum, and the ECB, HVB and EVB energy limits), §3.3.1 (the Fermi–Dirac supply
  N = kT ln[(1 + e^((E_F1 − E)/kT)) / (1 + e^((E_F2 − E)/kT))]), §3.5.1 (the WKB transmission
  exp(−(2/ħ) ∫ √(2 m_diel (W − E)) dx) between classical turning points, and its direct and Fowler–Nordheim closed
  forms for a linear barrier), §4.2.2 (a stack is evaluated once per path, its current returned to the end nodes'
  continuity rows or the contact's current), §5.1.2 (SiO₂ barriers of 3.2 eV for electrons and 4.65 eV for holes
  in its simulations; ECB, HVB and EVB in n- and pMOS) and Table 5.1 (permittivity, gap and band offsets to silicon
  of nine dielectrics; SiO₂ 3.9, 9.0, 3.0, 4.9 eV; HfO₂ 25, 5.7, 1.5, 3.1 eV; from Leroy 2003, Wilk 2001). Its
  equations are images; their structure was read from the text and is re-derived below.
- Electron tunnelling mass in SiO₂: Brar, Wilk and Seabaugh, Appl. Phys. Lett. 69, 2728 (1996), 0.30 m₀ with a
  parabolic dispersion and 0.41 m₀ with a non-parabolic one, from 1.6–3.5 nm oxides (abstract only, via a search
  summary; the paper is paywalled). Hole mass: Yang et al., IEEE Trans. Electron Devices 47, 2161 (2000),
  0.51 m₀ parabolic, fitted to p⁺-poly pMOSFETs with 1.23–2.16 nm oxides (abstract only). Neither could be read.
  **Owner ruling (2026-10-28):** take the more common values instead: 0.50 m₀ for electrons and holes, the round
  value usually assumed for SiO₂ (Ng et al., J. Appl. Phys. 2004, quoted in a search summary: "usually assumed to
  be 0.5 m_e or 0.4 m_e"); the hole value is also Yang's 0.51 rounded (GD5).
- Legacy: `workbench/physics/tunneling.py` (M12-S1): Fowler–Nordheim and a rectangular WKB transmission as
  diagnostics at the solved oxide field, not self-consistent. OLD / NEW / REASON:
  - OLD: J_FN = (A/φ) E² e^(−B φ^(3/2) / E) with the free-electron mass, evaluated after the solve.
  - NEW: the Tsu–Esaki integral with the insulator's tunnelling mass, in the residual and Jacobian.
  - REASON: rulings 2 and 3; FN misses direct tunnelling, the regime below about 3 nm, and a diagnostic does not
    conserve current.

**Design:**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| GD1 | Energies on one scale per node: the semiconductor's E_c, E_v (from ψ, χ, the band shift, as `band_edges`), its quasi-Fermi levels E_Fn, E_Fp (from n and p through the node's statistics, Boltzmann or the inverse F_1/2); the electrode's Fermi level E_FG = −q V_G. The barrier is the vacuum level minus χ_ox: at the semiconductor end E_c + χ_s − χ_ox, at the electrode end E_FG + φ_m − χ_ox, the insulator's valence edge E_g,ox below it. | Endpoint vacuum levels make the lumped barrier exact for a fixed interface charge (it changes the field, not the ends) and need no V_ox. | — |
| GD2 | Lumped gate (ruling 1): one linear segment of length t_ox between the two ends of GD1, the gate stack's insulator band data (GD5). Electrode (ruling 1): from each semiconductor node with a face on an insulator, a straight path along that face's grid axis, node to node through insulator nodes, ending at a node of an `electrode` contact; no path if it meets a semiconductor node, the mesh boundary or no electrode. Each mesh edge on the path is one linear segment between its nodes' vacuum levels with the band data of the edge's region. A node with insulator faces on several axes has one path per face, each weighted by that face's area. | Tensor grids make the normal a grid line (no locator); a piecewise-linear barrier is exact for the box method's linear potential along an edge, and per-edge band data gives multi-layer stacks (HfO₂ on SiO₂). | A path is frozen geometry (as Unit 20's); a curved or oblique interface is out (unstructured refused, ruling 7). |
| GD3 | WKB: TC(E) = exp(−(2/ħ) Σ_segments ∫ √(2 m (W(x) − E)) dx) over the segments where W > E, each integral in closed form, (2/(3s)) [(W − E)^(3/2)] between the segment's ends or turning point, s the segment's slope (the limit √(W − E) L when s = 0). Electron barriers W = E_c,ox(x) with the insulator's electron mass; hole barriers W = E − E_v,ox(x) mirrored, with its hole mass; EVB electrons see the conduction-band barrier with the insulator's EVB mass (default its electron mass). Above the barrier's highest point no energy is integrated. | Closed form per segment: no spatial quadrature, exact derivatives. | No over-barrier (thermionic) current: negligible for 3 eV barriers, not for a low-offset high-k at high field; recorded. No image force (ruling 7), no quasi-bound states (Tsu–Esaki's continuum overestimates inversion current, Gehring §3.2); recorded. |
| GD4 | Tsu–Esaki: J = (A* / k²) ∫ TC(E) N(E) dE (N has units of energy; with TC = 1 above E_c and a Boltzmann supply it is A* T² e^(−(E_c − E_F)/kT), thermionic emission), N(E) = kT ln[(1 + e^((E_F,s − E)/kT)) / (1 + e^((E_FG − E)/kT))], A* the semiconductor carrier's Richardson constant (`physics::emission_velocity_cm_s`'s: the material's, or its density-of-states mass). Limits per process (Gehring §3.2): ECB from max(E_c,s, E_c,G) to the barrier top, E_F,s = E_Fn; HVB from the barrier bottom to min(E_v,s, E_v,G), mirrored, E_F,s = E_Fp; EVB between E_c,G and E_v,s (electrons from the semiconductor's valence band to the gate's conduction band) and between E_c,s and E_v,G (the reverse pair), E_F of a valence band its E_Fp. A metal electrode has no band edges: no lower limit for ECB, no upper for HVB, and no EVB (its exchange with the valence band is HVB). Polysilicon: degenerate silicon, E_c,G = E_FG (n⁺) or E_v,G = E_FG (p⁺), the other edge E_g,Si(T) away (as `electrode_work_function_eV`). The integral is a fixed composite Gauss–Legendre rule (4 points, panels of at most kT/2), no adaptivity. | (4π m q k² / h³) = A*: one prefactor shared with Schottky and thermionic emission. A fixed rule keeps the residual smooth and its Jacobian the exact derivative of the discrete sum (moving limits through the panel map). | Hundreds of WKB evaluations per node per assembly; GC14 bounds the cost. |
| GD5 | Insulator band data: `InsulatorParameters` gains χ_ox, E_g,ox and the tunnelling masses (electron, hole, EVB), 0 meaning absent. SiO₂'s: χ_ox = 1.05 eV (Δ E_c = 3.0 eV to silicon's χ = 4.05, Gehring Table 5.1), E_g,ox = 9.0 eV, m_n = m_p = 0.50 (owner ruling above; EVB defaults to m_n). A lumped gate stack takes SiO₂'s unless it gives its own; a stack whose permittivity is not 3.9 must give its own. An electrode path crossing an insulator without band data is refused. | Ruling 5. | The SiO₂ masses are fitted values tied to their barriers; any change is an owner ruling. |
| GD6 | Rows: the semiconductor node's electron row gains the ECB and EVB fluxes into the gate (carriers leave the node when the semiconductor side's supply is the larger) times the face area; the hole row gains the HVB flux and the EVB flux as a hole flux (an electron leaving the valence band leaves a hole). The Poisson row is unchanged. The partials reach ψ, n and p of the node and ψ of every insulator node on its paths, and the electrode's bias (`bias_derivative`). The pattern gains the path entries at creation. | Ruling 2; Gehring §4.2.2: the current returns to the end nodes' continuity rows and the contact. | Path entries widen the electron and hole rows of gate-interface nodes. |
| GD7 | Currents: the gate's (or electrode's) conduction current is the sum of its nodes' (or paths') tunnelling fluxes, in `terminal_currents`' sign; its charge is unchanged. Transient and small-signal take the term's partials (no storage of its own). | Conservation: the semiconductor contacts carry the same current through the edges. | — |
| GD8 | Switch: per contact, `tunnelling = off` (the default) or a list of processes, `ecb`, `hvb`, `evb` (`on` = all three). With it off the term is not created: residual, Jacobian pattern and values, run identities and device files are those of `main`. | Ruling 4 with opt-out per process for comparison; GC1. | — |
| GD9 | Refused (invalid_input): with `electrothermal`; with the weighted terminal current; electrode tunnelling on a device without `DeviceDescription::cells` (unstructured); a crossed insulator without band data; masses or E_g,ox not finite and positive. | Ruling 7. | Follow-ups recorded in ARCHITECTURE 11.3. |
| GD10 | The deck: `[contact name] tunnelling = …`; on a gate contact optional `oxide_affinity`, `oxide_band_gap`, `oxide_mass_n`, `oxide_mass_p`, `oxide_mass_evb`; on an insulator region `affinity`, `band_gap`, `tunnel_mass_n`, `tunnel_mass_p`, `tunnel_mass_evb`. The device file stores them; the run record names the processes and band data; `BiasPoint` gains the per-contact tunnelling current by process; the app vocabulary and CLI name the keys. | Deck parity; the per-process currents are what a user compares to measurements. | The device file's format version rises if the stored contact or insulator needs it. |

**Acceptance criteria, fixed before any implementation:**

| # | Criterion |
|---|---|
| GC1 | Tunnelling off (and every existing deck): every existing test and pin passes unchanged; residuals, Jacobian patterns and values, run identities and the Unit 40 traces unchanged; a device file without tunnelling is byte-identical to `main`'s. |
| GC2 | WKB closed forms: for linear, rectangular and three-layer piecewise-linear barriers, with and without a turning point inside a segment, TC against a 2000-panel Gauss–Legendre quadrature of the same integrand within 1e-10 relative; the rectangular barrier equals exp(−2 √(2 m (W − E)) d / ħ) within 1e-13. A two-layer stack of identical band data equals one layer of the summed thickness within 1e-13. |
| GC3 | Fowler–Nordheim limit: ECB from a degenerate source through a single SiO₂ layer at 10–14 MV/cm, slope of ln(J/E_ox²) against 1/E_ox equal to −4 √(2 m_n) (qΦ_B)^(3/2) / (3 q ħ) within 2%. Direct regime: at a fixed oxide voltage of 0.5 V, the slope of ln J against t_ox over 1.2–2.0 nm equals −2 κ̄ of the trapezoid (its mean-barrier WKB exponent per length) within 3%. |
| GC4 | Supply and equilibrium: the supply equals its Boltzmann closed form within 1e-12 where both levels are 10 kT below E; at thermal equilibrium (every process, n- and p-substrate, metal and both polysilicon electrodes, both statistics) the tunnelling current is zero within 1e-14 of its value at 1 V. |
| GC5 | Quadrature: halving the panel width changes J by at most 1e-6 relative at every GC3, GC4 and GC10 point. |
| GC6 | Energy limits: EVB is exactly zero while E_c,G > E_v,s (and the reverse pair's condition) and positive past it; ECB to a metal and HVB to a metal have no gate-side limit; a polysilicon gate's limits follow its edges (checked on a 1D n⁺-poly nMOS at −1 … +2 V). |
| GC7 | Jacobian: against central differences of the residual with tunnelling on (1D lumped gate; 2D and 3D tensor electrode stacks of one and two layers; Boltzmann and Fermi–Dirac; each process alone and all three), the largest relative column error at most 1e-6, the off-node ψ columns of the path included. Newton converges quadratically on GC8's device at +1.5 V. |
| GC8 | Self-consistent MOS capacitor (1D, p-Si 5e17, n⁺ poly, 1.5 nm SiO₂, lumped, V_G = −2 … +2 V by 0.1 V): every point converges; the gate current equals the substrate contact's current within 1e-10 relative; the ECB current at each point equals the GD4 integral evaluated on the solution's own band edges and quasi-Fermi levels within 1e-12. |
| GC9 | Lumped against meshed: the same capacitor with the oxide meshed (0.05 nm spacing, an `electrode` n⁺ poly) gives the lumped gate current within 2% at +1 V and +2 V. |
| GC10 | Polarity and dominance: on GC8's nMOS in inversion ECB is the largest process; on its pMOS mirror (n-Si 5e17, p⁺ poly) in inversion HVB is; the per-process currents carry the signs of GD6. |
| GC11 | Geometry: a 2D tensor-grid device uniform across the gate gives the 1D current per unit area within 1e-9 relative; a 3D tri-gate fin node at a corner takes one path per face, each with its face's area (checked against the sum of two separately built 2D devices within 1e-9). |
| GC12 | Dynamics: the small-signal conductance at 1 Hz against DC differences of the gate current within 1e-5 (GC8, +1.5 V); after a gate step the total currents of all contacts sum to zero within 1e-10 relative at every time step. |
| GC13 | Refusals of GD9 and of a negative or absent required band datum, a gate stack with permittivity ≠ 3.9 and no band data, an unknown process name: each invalid_input with the offending value in its context. |
| GC14 | Cost: on the FinFET deck, assembly time with all three processes at most 1.5 times without (median of 31, Release, P-core pinned as Unit 40). The FinFET sweep with tunnelling converges at every point that converges without it. |
| GC15 | The deck, device file (round trip), run record, `BiasPoint` currents, CLI and app vocabulary carry GD10's keys and results. |
| GC16 | Reported, not gated: GC8's I_G–V_G against a measured ultrathin-oxide J–V if an accessible source is found (then its tolerance is ruled before the comparison is run); the FinFET deck's gate current against its drain current. |
| GC17 | Full Release and Debug suites; zero warnings; hand mutants on each GD row caught or recorded as equivalent. |

**Owner rulings on the criteria (AskUserQuestion, 2026-10-28), after the first deck-level runs.** Measured on GC8's
capacitor (lumped, all processes, the deck defaults): its gate current is SRH-generation-limited (2.98e-9 A/cm² at
+0.5 V, 8.59e-9 at +2 V); the inversion electrons' quasi-Fermi level lies 1.2e-6 to 9.6e-5 eV above the gate's; the
substrate's `terminal_current_resolution` is 1.5e-5 to 8.7e-7 A/cm², so its current (−3.04e-9, −9.11e-9) is not
resolved. With SRH and Auger off the gate current is 5.6e-12 A/cm² and the two currents agree to about 1e-7.

| # | Ruling |
|---|---|
| GR1 | GC8's conservation check (gate current against the substrate's) runs with SRH and Auger recombination off, so the current is tunnelling-limited and resolved at the substrate; tolerance as GC8 (1e-10 relative), else its measured value goes back to the owner. The rest of GC8 keeps the deck defaults. |
| GR2 | GC8's ECB recomputation: \|ΔJ\| ≤ 1e-12 × the one-way flux, the GD4 integral with the electrode's supply removed (only the semiconductor side supplying), which is the scale the rounding of the inputs resolves; the net current near quasi-equilibrium is a small difference of the two one-way fluxes. |
| GR3 | GC12's small-signal conductance is compared at 0 Hz (not 1 Hz) with the DC difference, same device and 1e-5: the capacitor's inversion layer responds in about 10³ s, so 1 Hz is its high-frequency response (G = 6.5e-8 against the DC 3.2e-9 S/cm², physical). |
| GR4 | GC9's meshed oxide (31 oxide nodes at 0.05 nm) trips `min_pivot_ratio` (1e-11; measured 6.7e-13 at the next silicon node's electrons, scaling with the tunnel entries; 2.7e-4 for the lumped gate and for paths of at most 6 nodes). With the guard off Newton converges in 5–9 iterations and the meshed current is the lumped one within 6e-5. Investigate the elimination first (ordering, structural entries) and fix it in the assembler or solver if possible; ask again if not. |

**GR4 investigation (2026-10-28).** A dense partial-pivoting LU of the same equilibrated Jacobian at the equilibrium
state gives pivot ratios of 2.7e-4 to 3.5e-4 in natural order, and the dense factors estimate ‖A⁻¹‖ (equilibrated,
max norm over sign vectors): meshed without tunnelling 5.7e18, meshed with HVB only 5.7e18, meshed with all processes
1.2e14; lumped without 8.4e17, lumped with all 1.2e14. The capacitor's Jacobian is near-singular with or without
tunnelling: the inversion electrons are almost floating (their link to the substrate contact through the depletion
region is exponentially weak; only generation and the tunnel current anchor them), and tunnelling improves the
conditioning by about four orders. Eigen's SparseLU pivots partially over each column's structural rows
(`SparseLU_pivotL.h`, threshold 1.0), so its smallest pivot is genuine for its COLAMD order; the path's row changes
that order, which exposes the near-singularity to the `min_pivot_ratio` guard, a floating-region detector whose
verdict depends on the order (it passes the equally conditioned lumped device). Without tunnelling the meshed device
passes the guard at 0 V and then does not converge at +0.5 V. No assembler or solver defect: the path entries change
the ordering, not the conditioning. Returned to the owner.

| # | Ruling |
|---|---|
| GR5 | GC9's meshed device runs with `[solver] min_pivot_ratio = 0`; recorded as a limitation: a tunnelling MOS capacitor without a minority-carrier contact is near-singular and may trip the guard (lower it, or give the inversion layer a source). |

**Measured after GR1–GR5, and ruled (AskUserQuestion, 2026-10-28).** GR1's setup (SRH and Auger off): gate and
substrate agree to about 2e-10 relative in inversion; in accumulation an absolute floor of 7.6e-10 A/cm², the same with
Newton's `tol_update` at 1e-8 and 1e-12 (the substrate terminal current's rounding of cancelling hole fluxes at
p = 5e17; 37% of the 2.6e-10 A/cm² at −0.1 V). GR3's 0 Hz admittance: Y_GG against the DC difference 4e-6 relative at
−1.0 V; at +1.5 V unresolved: +2.7e-9 to −4.8e-8 S/cm² against the DC 3.19e-9 as the linear solve's backward error goes
from 1.4e-11 to 1.7e-16, the solver's own reported resolution 6e-7 to 8e-7 S/cm² (the gate current there is a 1e-5
cancellation between its bias term and the inversion electrons, which follow the gate). GC12's transient (1.4 → 1.5 V,
GR1's setup), the worst relative sum over steps with a current above 1e-3 A/cm²: 2.5e-8 without tunnelling (the
transient's own floor), 6.1e-8 with it. GR2 as worded scales by the semiconductor side's one-way flux, about 1e-16
A/cm² in accumulation, where the current comes from the gate (the measured ECB agreement there is 6e-17 A/cm² on
0.126).

| # | Ruling |
|---|---|
| GR6 | GC8's conservation (GR1's setup): \|I_G + I_sub\| ≤ 1e-10 \|I_G\| + 1e-9 A/cm². |
| GR7 | GC12's small-signal comparison (GR3, 0 Hz) at −1.0 V (accumulation), within 1e-5. The 0 Hz gate admittance of a generation-limited capacitor in inversion is recorded as unresolved; its reported resolution covers the miss. |
| GR8 | GC12's transient: \|Σ I\| ≤ 1e-7 × max \|I\| at each step with max \|I\| > 1e-3 A/cm² (GR1's setup applied to the transient). |
| GR9 | GR2's scale is the larger of the two one-way fluxes (the semiconductor's side supplying alone, and the electrode's alone). |
| GR10 | The FinFET deck's HfO₂ band data for GC14 and GC16 (AskUserQuestion, 2026-10-28): gap 5.7 eV and affinity 2.55 eV (conduction offset 1.5 eV to silicon's 4.05; Gehring Table 5.1), tunnelling masses 0.18 m₀ for electrons and holes (owner ruling: no accessible source found for 0.18; the open-access measurement found, Monaghan et al., Solid-State Electronics 2009, https://research.ucc.ie/en/publications/determination-of-electron-effective-mass-and-electron-affinity-in/, gives 0.11 ± 0.03 m₀ and an affinity of 2.0 ± 0.25 eV). The interfacial SiO₂ takes the built-in set. |
| GR11 | GC14's cost (AskUserQuestion, 2026-10-28): the tunnels' fluxes are computed on threads, the first threading inside assembly; each flux is computed whole by one thread and added in tunnel order, so the result does not depend on the thread count. The baseline (without tunnelling) stays single-threaded. Measured before the ruling, single-threaded: 4.7× (about 1,300 tunnels × 2 processes × 320 quadrature points × about 10 segments, near 55 ns a point). |
| GR12 | Gate tunnelling with quantum confinement (owner, 2026-10-29, after the rebase onto `main` brought Unit 47): refused (invalid_input) as GD9's exclusions, tested in GC13. MLDA moves the surface quasi-Fermi level by −kT ln f_n, which the supply reads, so the same density tunnels more; neither unit validated the pair. A later unit may rule on the supply under confinement and lift the refusal. |

**Slice 2 (device, assembly, deck, files, record, app, CLI), as built.**
- `device::Contact::tunnelling` (`TunnelProcess` bits; 0 off) and `GateStack`'s oxide band data (`oxide_affinity_eV`,
  `oxide_band_gap_eV`, `oxide_mass_n`, `oxide_mass_p`, `oxide_mass_evb`; 0 absent). The device refuses tunnelling on
  anything but a gate or an electrode, unknown process bits, a partial or non-finite band-data set, and a tunnelling
  gate whose permittivity is not 3.9 without band data (GD9).
- `physics::silicon_dioxide_tunnelling_parameters` holds SiO₂'s band data; `silicon_dioxide_parameters` is unchanged
  and the assembler substitutes the tunnelling set for it, so device files and run identities without tunnelling are
  unchanged (GC1).
- `assemble/gate_tunnels.cpp`: a lumped gate gives one tunnel per gate node (one piece of t_ox, its patch area); an
  electrode gives one per semiconductor–insulator grid edge whose straight grid line through the insulator reaches one
  of its nodes (the edge's coupling area). The barrier starts at the first edge's midpoint (the box method's interface;
  the semiconductor's half is not barrier), with the permittivity-weighted potential there; an edge between two
  insulator materials is split at its midpoint the same way. Electrode paths need `DeviceDescription::cells` (else
  invalid_input). Each process's value and gradient are carried over the local columns (ψ, n, p of the node, ψ of each
  path node, a lumped gate's bias), so the Jacobian is the exact derivative of the discrete flux.
- The term `GateTunnelling` (continuity rows, global domain) adds −T_n to the electron row and +T_p to the hole row;
  the pattern gains the path's ψ columns in both. `terminal_currents` adds T_n − T_p to the tunnel's contact;
  `small_signal_currents` its partials and a lumped gate's bias term; `bias_derivative(contact, x)` the lumped bias
  partials. A lumped tunnelling gate makes dF/dV depend on the state: `bias_nonlinear()`, read by the bordered trace and
  the small-signal forcing.
- Refusals (GD9): electrothermal with tunnelling (`DriftDiffusion::create`), the weighted terminal current
  (`weighted_current.cpp`).
- Deck: `[contact] tunnelling = off | on | ecb, hvb, evb` and the gate keys `oxide_affinity`, `oxide_band_gap`,
  `oxide_mass_n`, `oxide_mass_p`, `oxide_mass_evb`. An insulator region's band data uses the parameter names, as every
  material parameter key does (`electron_affinity_eV`, `band_gap_eV`, `tunnel_mass_n`, `tunnel_mass_p`,
  `tunnel_mass_evb`), not GD10's shorter names.
- Device file: contact key 7 (`tunnelling`) and gate-record keys 7–11 (band data), insulator-parameter keys 2–6, each
  written only when nonzero, so files without tunnelling are byte-identical (format version unchanged).
- Run record: insulator band data and each tunnelling contact's processes and gate band data hashed only when present;
  each tunnelling contact named as the setting `tunnelling.<contact>` (its process bits).
- `BiasPoint::tunnelling_ecb`, `_hvb`, `_evb` (per contact, terminal_current's unit and sign; empty without
  tunnelling), sweep-file keys 19–21; the CLI's sweep table gains `<contact> ECB A`, `HVB A`, `EVB A` per tunnelling
  contact; the app vocabulary names the six contact keys (the gate's `o…` completion test now lists the five oxide
  keys).
- Tests: `tests/assemble/gate_tunnelling_test.cpp` (GC6, GC7 Jacobians, equilibrium, bias and small-signal partials,
  GC13); `tests/deck/gate_tunnelling_test.cpp` (GC7 Newton, GC8, GC9, GC10, GC11 in 2D and 3D, GC12, GC15; hidden: GC14 cost and sweep, GC16). Measured: GC8
  ECB against the integral on the solution's bands 6e-17 A/cm² on 0.126 (accumulation); GC9 0.000% and 0.006%;
  GC12 4e-6 (−1 V) and the transient 6.1e-8; GC7 corrections 1.0e-2, 3.5e-5, 2.5e-11. With Newton's `tol_update` at
  1e-12 (not the deck's 1e-8) the same transient's sum reaches 1.4e-4 (an absolute 2e-7 A/cm² while the currents decay
  through 1 A/cm²); not investigated.

**GC14's cost, as built (GR11).**
- `physics::tsu_esaki_current` takes the WKB sum inline: each segment's factor scale·√m once per call, each end's root
  once per point (adjacent segments share an end), the supply's softplus and logistic from one exponential, scratch on
  the stack. Bit-identical to slice 2 (the same operations on the same values): GC7's Newton corrections are unchanged.
  A variant with one division per segment instead of two (about 7% faster) changed the rounding and GC7 then failed
  (corrections 3.1e-3, 6.8e-5 at +1.5 V): the capacitor's near-singular Jacobian (GR4) makes GC7's sequence depend
  on the last bit of the barrier action. Not adopted; recorded as a sensitivity of GC7.
- `DriftDiffusion::GateTunnelJob` (`terms.hpp`, `gate_tunnels.cpp`): `Composition::assemble` starts it before the
  terms, worker threads (the process's affinity count less one, at least two 16-tunnel chunks each) claim chunks of
  tunnels in turn, the other terms add theirs meanwhile, and the `GateTunnelling` term's `finish()` has the assembling
  thread take the remaining chunks, joins the workers and scatters the fluxes in tunnel order. `gate_tunnel_fluxes`
  (terminal and per-process currents) runs the same job at once. Each flux's local gradients live in a per-thread
  arena (no allocation once grown).
- Measured (FinFET, Release, P-core pinned 0x3C3, median of 31): evaluate 13.2 ms without, 15.7 ms with (1.195×);
  residual 9.6 and 13.6 ms (1.42×). Steps: 5.6× (slice 2), 4.7× (bit-identical arithmetic), 1.8× (threads after the
  other terms), 1.2× (threads alongside them).

**GC14's sweep, GC16 and GC11 in 3D, as measured (2026-10-29).**
- GC14 (`[.gate_tunnelling_finfet_sweep]`, PARDISO on 6 threads, the deck's residual line search): without tunnelling
  17 of 17 points, with all three gates tunnelling (HfO₂ band data of GR10) 17 of 17. Both sweeps together took
  69.8 min.
- GC16, reported: the FinFET's gate current against its drain current at V_D = 50 mV. I_G is −3.5 fA at 0 V
  (−1.6e-5 of I_D), 2.1e-6 to 5.3e-6 of I_D from 0.05 to 0.4 V, and rises to 2.29 nA against 18.8 µA at 0.8 V
  (1.22e-4). No measured ultrathin-oxide J–V source was found, so that comparison is not run.
- GC11 in 3D (`tests/deck/gate_tunnelling_test.cpp`): a silicon corner with 1.5 nm of SiO₂ on two faces and an
  n+ poly gate on each, against the same section in 2D with one gate tunnelling per device, at one fixed state
  (the section's equilibrium potential under 1.5 V on both gates, copied along z). Each 3D gate carries 1 µm times
  its 2D device's current; equal to the printed seven digits (8.949119e-10 and 1.002478e-09 A), within 1e-9.
  Mutant: a node keeping only its first tunnel (the corner losing its second face's path) fails by 8.7e-4.

**GC1's traces, as measured (2026-10-29, after the rebase and GR12).** Unit 42's PG1 method: the Unit 40 assembly
traces of `main` (`329f203`) and of this branch, Release, each run pinned to the P cores from its start, in
declaration order, with the MKL runtime loaded, over the default cases of the assemble, solve, deck and validation
executables and the hidden `[mosfet]`, `[transient]`, `[breakdown]`, `[btbt]`, `[unstructured]`, `[mlda]` (solve
and deck) and `[lombardi]` entries. Every one of `main`'s 469 test-case blocks is byte-identical in the branch's
traces (assemble 127, solve 220, deck 73, MOSFET 11, breakdown 7, transient 6, BTBT 3, unstructured 3, solve MLDA 4,
deck MLDA 2, Lombardi 5, validation 8). The only blocks without a counterpart are this unit's 19 new test cases.
With tunnelling off nothing the assembler computes changed, including the threaded assembly (GR11).

**GC17's hand mutants, as run (2026-10-29).** 23 mutants (the last for GR12's refusal, after the rebase), each built in Release and run against the gate
tunnelling tests (physics, assemble, deck), then reverted. All are caught; none is recorded as equivalent. Seven
first survived and each gap was closed with a test before the mutant was run again:

| Row | Mutant | Caught by |
|---|---|---|
| GD1 | E_c without the band shift | survived (all-silicon tests: the shift is 0); new: GC8's ECB on a germanium surface with the gate at x_max, off the reference material |
| GD1 | the lumped gate's vacuum level +0.1 eV | GC8, GC15 |
| GD2 | the corner node keeps one tunnel | GC11 in 3D |
| GD2 | the lumped oxide 1.1 t_ox | GC8 |
| GD2 | an electrode path continuing through semiconductor nodes | survived; new: no path through a semiconductor (fixed states, electrode / SiO₂ / Si / SiO₂ / Si) |
| GD3 | the transmission e^(−S) for e^(−2S) | GC3 slopes |
| GD3 | the turning-point segment's action halved | GC3 FN slope, the current's partials |
| GD3 | the two-sided segment's end partials swapped | the segment and current partials |
| GD4 | the supply cut at 10 kT | GC5 |
| GD4 | no derivative through the barrier top at the top cut | survived (no case reached that cut); new: a fifth partials case with the levels near the top |
| GD4 | an end root shared across a discontinuous step | GC5, the current's partials |
| GD4 | the supply logistic's branches swapped | the current's partials |
| GD5 | EVB's mass falling back to the hole mass | survived (every insulator had equal masses); new: unequal masses in the band data test, and SiO₂'s ruled values pinned |
| GD6 | EVB's semiconductor-valence direction with the wrong sign | survived (only != 0 was checked); new: GC6 at +3 V, the rows −T_n and +T_p against the reported EVB current |
| GD6 | EVB's gate-valence direction with the wrong sign | GC6 at −3 V (added with the one above) |
| GD6 | HVB with the electron mass | new: each process takes its own mass (fixed state, unequal masses through the gate's keys) |
| GD7 | HVB's terminal current with the wrong sign | GC10 |
| GD8 | `on` without EVB | the assemble tests (GC7 stacks) |
| GD9 | no electrothermal refusal | GC13 |
| GD9 | no weighted-current refusal | survived; new: the weighted current refused with gate tunnelling |
| GD9 (GR12) | no quantum-confinement refusal | GC13 (added with the refusal; checks the message) |
| GD10 | the gate's `oxide_mass_n` and `oxide_mass_p` swapped | each process takes its own mass |
| GD10 | a region's `electron_affinity_eV` into its band gap | survived (only the hidden FinFET case uses region keys); new: an insulator region's band data from the deck |

The germanium case needs the floating-region guard off, as GC9 (GR5): with the gate at x_max the capacitor's
near-singularity meets the guard. A lumped `n_poly` gate takes its work function from the channel's material (the
gate electrostatics' convention, `assemble::gate_work_function_eV`), so on germanium it is an n+ germanium gate.
