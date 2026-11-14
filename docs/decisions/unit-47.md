<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 47 — quantum confinement by the modified local density approximation, opt-in (started 2026-10-24; on `physics/quantum-confinement`)

**Owner request:** "do #1" (2026-10-24), MLDA quantum confinement, the second of the nanoscale models of
ARCHITECTURE.md 11.3. **Owner rulings (2026-10-24):** (1) the density-of-states factor is applied under both
statistics, exact MLDA under Boltzmann and Hänsch's approximation under Fermi–Dirac; (2) refused with
electrothermal; (3) the quantization masses default to Ioffe's m_l and heavy-hole values and can be overridden in the
deck, silicon only; (4) at a lumped-gate node (the interface at the node) the factor is its mean over the node's
half-cell. The standing conditions hold: opt-in and bit-identical when off, criteria fixed before code, Release and
Debug suites at every commit, `/W4 /WX`, hand mutants, nothing pushed or merged without authorization.

**Sources read (before coding).** G. Paasch, H. Übensee, phys. stat. sol. (b) 113, 165 (1982), and W. Hänsch et al.,
Solid-State Electron. 32, 839 (1989), are paywalled and were not read. The model is read from:
(a) M. Ieong, R. Logan, J. Slinkman, "Efficient quantum correction model for multi-dimensional CMOS simulations",
SISPAD 1998 (in4.iue.tuwien.ac.at/pdfs/sispad1998/pdfs/Ieong_34.pdf; text extracted with `pdftotext`). Its MLDA
density, eq. (1), is an image and was not read; its text gives the ingredients (the distance z to the interface, N_c,
the quasi-Fermi level, the zeroth-order spherical Bessel function j0, the thermal wavelength λ = ħ / sqrt(2 m kT)),
says the expression reduces to Fermi–Dirac statistics for large z, and describes the implementation used here: the
corrected density is turned into a quantum potential that is lumped into the band-structure parameters of the
heterostructure current equation, so equilibrium carries no current;
(b) S. Entner, PhD thesis, TU Wien, section 2.4.1.1 (iue.tuwien.ac.at/phd/entner/node12.html), citing Paasch–Übensee
and Hänsch: eq. (2.106) h = 1 − exp(−(z + z_0)^2 / (ζ^2 λ^2)), (2.107) λ = ħ / sqrt(2 m* k T), (2.108)
N_c = N_c,0 h; the correction does not depend on bias and does not move the band edges; z_0 = 0 makes the density
of states vanish at the interface, "which causes numerical problems" (ζ = 1, z_0 = 1 nm in its figure);
(c) Ioffe's silicon band-structure page (ioffe.ru/SVA/NSM/Semicond/Si/bandstr.html, read over HTTP): m_l = 0.98 m0,
m_t = 0.19 m0; heavy holes 0.49 m0, light 0.16 m0.

The Boltzmann reduction is derived here, not read: with the MLDA kernel 1 − j0(2 z sqrt(ε) / λ),

    ∫0^∞ sqrt(ε) e^(−ε) [1 − j0(2 a sqrt(ε))] dε / (sqrt(π) / 2) = 1 − exp(−a^2),   a = z / λ,

(u = sqrt(ε) and ∫0^∞ u e^(−u^2) sin(b u) du = sqrt(π) b e^(−b^2/4) / 4), which is (b) with z_0 = 0 and ζ = 1. QC2
checks it by quadrature. With the masses of (c), λ_n = 1.2263 nm and λ_p = 1.7343 nm at 300 K (λ for m0 is
1.2140 nm).

**Design:**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| QD1 | `PhysicsModels::quantum_confinement`: `off` (the default) or `mlda`; the masses `mlda_mass_n` (default 0.98) and `mlda_mass_p` (0.49), in m0; the deck's `[models] quantum_confinement`, `mlda_mass_n`, `mlda_mass_p`; the run record names the model and both masses only when `mlda`; the app's deck vocabulary lists the keys. | Opt-in, as Units 34, 37, 38, 46; ruling (3). | — |
| QD2 | Per semiconductor node i and carrier c, f_c,i = 1 − exp(−d_i^2 / λ_c^2), λ_c = ħ / sqrt(2 m_c m0 k T) at the device temperature; d_i the distance from the node to the nearest interface point. Interface points: the midpoint of every edge between a semiconductor node and an insulator node (where the assembler puts that interface; mesh/interfaces.hpp), and every lumped-gate node itself (d = 0). One nearest distance: at a corner the nearer interface decides. | The Boltzmann MLDA (derived above); the interfaces as the assembler already represents them. | Corners and thin bodies get a single-interface correction (no overlap of two walls' factors). Heterojunctions and contacts are not confining interfaces. |
| QD3 | A lumped-gate node's factor is the mean of f over z in [0, ℓ/2], ℓ the node's shortest carrier edge: 1 − (sqrt(π) λ / (2 a)) erf(a / λ), a = ℓ/2. Every other node takes f at its own d. | Ruling (4); positive, no fitted parameter, the node's share of the charge converges as the mesh is refined. | On a tensor grid ℓ/2 is the box depth; on an unstructured mesh it is an estimate of it. |
| QD4 | The factors enter only through the scaled device (`make_scaled_device`): N_c → N_c f_n and N_v → N_v f_p, band edges unchanged. With the existing fields that is n_ie → n_ie sqrt(f_n f_p), s → s + ½ ln(f_n / f_p), ln(N_c / n_ie) → + ½ ln(f_n / f_p), ln(N_v / n_ie) → + ½ ln(f_p / f_n). Every assembler (equilibrium Poisson, drift-diffusion, transient, small-signal, the quasi-static sweep), the contact values, SRH and the trap and surface-recombination terms read these, so the model needs no new unknowns, its terms are constant in the solution and the Jacobian is unchanged in form. | Ieong's lumping into band parameters; ψ stays the electrostatic potential and the gate and electrode offsets keep their meaning; equilibrium is current-free by the existing heterostructure flux (Unit 15). | Insulator nodes untouched. Band-diagram output shows unchanged band edges; the densities carry the correction. |
| QD5 | Ruling (1): the same factors under Fermi–Dirac (the DOS form of Hänsch; Entner (2.108)). | One model for both statistics. | Under Fermi–Dirac this is Hänsch's approximation, not the j0 integral at the quasi-Fermi level; documented. |
| QD6 | Refused (invalid_input): with electrothermal (ruling (2)); on a device with a semiconductor region whose Caughey–Thomas parameters are not silicon's (as LD6); masses not finite and positive; with `mlda` and no confining interface (no insulator node and no lumped gate). | λ is evaluated at the device temperature; the masses are silicon's. | — |
| QD7 | The distance search: squared distances to all interface points, brute force, once per assembler construction. | Simple; measured in QC12. | O(nodes × interface points); a bucket grid only if QC12 shows it matters. |

**Acceptance criteria, fixed before any implementation:**

| # | Criterion |
|---|---|
| QC1 | Off (the default): every existing test and pin passes unchanged; run identities and Jacobian patterns unchanged. |
| QC2 | λ against ħ / sqrt(2 m m0 k T) from the CODATA constants within 1e-14 relative (300 K: 1.2263 and 1.7343 nm within 1e-4 nm); the factor f(d) against an independent quadrature of the Boltzmann MLDA integral above, within 1e-10, at d / λ = 0.01 … 10; the half-cell mean of QD3 against quadrature within 1e-12, and it tends to f's mean as a → 0 and to 1 as a → ∞. |
| QC3 | Distances: on 2D and 3D tensor meshes with a planar meshed oxide, each semiconductor node's d equals its distance to the interface plane within 1e-12 relative; on a 2D Delaunay mesh with a mirrored-pair interface, within 1e-12 for the nodes adjacent to the interface; at a 2D corner (two oxide walls) d is the smaller wall distance; a lumped-gate node has d = 0. |
| QC4 | Scaled device: N_c f_n and N_v f_p reconstructed from the modified fields within 1e-14 relative at every node; band edges (node_bands) unchanged within 1e-14 V_T; insulator nodes and nodes with f = 1 bit-identical to the model off. |
| QC5 | Equilibrium with the model on (a meshed-oxide 1D MOS-C in inversion and accumulation, the 2D MOSFET example at Vg = 1 V, both statistics): every edge current within its rounding bound (as Unit 11's BGN gate); n p = n_ie^2 f_n f_p within 1e-12 relative at every semiconductor node; SRH and surface recombination zero at equilibrium. |
| QC6 | The Jacobian against central differences with the model on (2D, Delaunay, 3D; Boltzmann and Fermi–Dirac; with a lumped gate and with a meshed oxide): largest relative column error at most 1e-6; Newton on the 2D MOSFET in strong inversion converges quadratically. |
| QC7 | 1D MOS-C, meshed 5 nm oxide, 1e17 p-type, surface spacing 0.05 nm (Release only): in strong inversion the electron density peaks inside the silicon at a distance between 0.3 and 2 nm (classically at the first node); the inversion charge at fixed V_G is lower than without the model; the inversion C–V gate capacitance (quasi-static) is lower than without the model; the threshold (Q_inv = 1e-7 C/cm^2) shifts up by 10–150 mV. In accumulation the hole charge is reduced likewise. |
| QC8 | Mesh convergence (QC7's MOS-C, Release only): the inversion charge at V_G − V_FB = 2 V on surface spacings 0.2, 0.1 and 0.05 nm; the last two differ by at most 1%, and less than the first two. |
| QC9 | A lumped-gate MOS-C against the meshed-oxide one of the same oxide capacitance (QC7's device, spacing 0.05 nm): the inversion charge at V_G − V_FB = 2 V within 5%. |
| QC10 | 2D = 1D: a y-uniform 2D MOS-C reproduces the 1D one within 1e-12; the 2D MOSFET example rotated 90° gives the same Id at Vg = 0.5 and 1.5 V within 1e-9; small-signal conductance against DC differences within 1e-5; a transient long step equals the steady state. |
| QC11 | Refusals: electrothermal, a germanium region, a zero or negative mass, no confining interface; accepted on silicon with another parameter changed, and with Lombardi together. |
| QC12 | Reported, not gated: the distance search time and the sweep times of the 2D and 3D MOSFET examples with and without the model; their Id at Vg = 1.5 V (2D) and 1.2 V (3D) with and without; the FinFET deck's V_th and on-current with and without (its mesh convergence is still open); Unit 46's LS6 gap with MLDA and Lombardi together. |
| QC13 | Full Release and Debug suites; zero warnings; hand mutants caught. |

**As built (2026-10-24).** `physics/quantum_confinement` (λ, the factor and the half-cell mean),
`assemble/quantum_confinement` (`mlda_factors`: the refusals, the interface points, the distances, the factors),
the MLDA block of `make_scaled_device` (QD4), the interface-trap level shift, the deck keys, the run record and the
app's vocabulary. Tests: `tests/physics/quantum_confinement_test.cpp` (QC2), `tests/assemble/quantum_confinement_test.cpp`
(QC3, QC4, QC6, QC11), `tests/solve/quantum_confinement_test.cpp` (QC5, QC10; and, Release only under `[mlda]`, ctest
`solve_mlda`: QC7–QC9 and the QC5 n p report), `tests/deck/quantum_confinement_test.cpp` (QC1; and, Release only
under `[mlda]`, ctest `deck_mlda`: the 2D example's equilibrium, Newton's order, QC10 and QC12), and the QC12 LS6
case in `tests/deck/surface_mobility_test.cpp` (also under `deck_mlda`).

| # | Added in the build | Why |
|---|---|---|
| QD8 | A device with a thermionic-emission interface is refused with `mlda`. | That interface's flux reads the band offsets and the densities of states on both sides directly; QD4 moves the latter, and how the factor should enter that flux was not designed. Silicon devices with a meshed oxide have none. |
| QD9 | `InterfaceEdge::level_shift`: an interface trap's level relative to the intrinsic level is moved by the node's ½ ln(f_n / f_p), so a declared trap energy keeps its meaning against the (unmoved) band edges. Zero with the model off (bit-identical: `tau + 0.0`). | Without it, QD4's shift of the intrinsic level would move every trap relative to the bands; QC4's trap-charge check (equal within 1e-12) caught it. |
| QD10 | QC6's finite-difference gate also runs with Lombardi on (2D and 3D). | QC11 accepts the two together; the 1e17 LS6 sweep with both stalled (QC12), and the gate shows the Jacobian is not the cause. |

**Results:**

| # | Result |
|---|---|
| QC1 | Met. With the model off every existing test and pin passes unchanged (one vocabulary pin grows with the keys, under QC13); the run record has no MLDA entries (tested). |
| QC2 | Met: 65 assertions. λ against long-double CODATA within 1e-14; at 300 K λ_n = 1.2263 nm, λ_p = 1.7343 nm. f against Simpson quadrature of the MLDA integral within 1e-10 at d / λ = 0.01 … 10; the half-cell mean within 1e-12, with both limits. |
| QC3 | Met on 2D and 3D tensor meshes, at a 2D corner (the smaller wall distance), on a 2D Delaunay mesh (nodes adjacent to the interface), and at lumped-gate nodes (d = 0, the half-cell mean). |
| QC4 | Met: N_c f_n and N_v f_p within 1e-14 under both statistics; node_bands unchanged; insulator nodes and nodes with f = 1 bit-identical; trap charges equal within 1e-12; the Jacobian pattern unchanged. |
| QC5 | **Met except the n p relation on the drift-diffusion state.** 1D MOS-C (accumulation −2 V, inversion 1.5 V; both statistics) and the 2D MOSFET example at Vg = 1 V: every terminal current within its resolution; SRH and interface rates within 1e-9 of their scale. n p = n_ie² f_n f_p holds within 6.7e-16 on the equilibrium (quasi-static) state at the same bias, and is gated there. On the drift-diffusion state the worst node misses 1e-12: 3.26e-12 in inversion, at p = 8.8e-5 cm⁻³; without the model the same state gives 7.45e-13. That is the rounding of the minority density in the solver's last update (a tighter `tol_update` does not change it), not the model. Reported by a hidden case; the criterion is unchanged. **For the owner.** |
| QC6 | Met: the finite-difference gate within 1e-6 on 1D, 2D, 3D, Delaunay and lumped-gate devices under both statistics, and (added in the build, below) with Lombardi as well in 2D and 3D. Newton on the 2D MOSFET example at Vg = 1.5 V: corrections 1.6e2, 9.0, 1.1e-1, 4.1e-6, 6.1e-15, an order of 1.99. |
| QC7 | Met (h = 0.05 nm, V_FB = −0.9805 V, V_G − V_FB = 2 V): the electron density peaks 0.800 nm inside the silicon (classically 0.025 nm, the first node); Q_inv 4.8610e-7 against 5.1824e-7 C/cm²; quasi-static C/C_ox 0.9153 against 0.9477; threshold shift +25.7 mV. Accumulation (V_G − V_FB = −2 V): Q_G −1.1904e-6 against −1.2466e-6 C/cm². |
| QC8 | Met: Q_inv 4.860640e-7, 4.860840e-7 and 4.860958e-7 C/cm² at 0.2, 0.1 and 0.05 nm; the differences are 4.1e-5 and 2.4e-5 relative. |
| QC9 | Met: lumped 4.86099e-7 against meshed 4.86096e-7 C/cm² (ratio 1.0000). |
| QC10 | Met. The y-uniform 2D MOS-C equals the 1D one within 1e-12 at every bias; the reflected 2D MOSFET example (x and y exchanged, as LS9) gives the same Id at Vg = 0.5 and 1.5 V within 1e-9; small-signal g_d and g_m match DC central differences to 8 digits (6.07386027e-1 and 7.14614488e-2); a 1D MOS-C transient step to −1.5 V ends at the DC state (ψ within 1e-9 V, p within 1e-7). |
| QC11 | Met: electrothermal, germanium, zero and negative masses, no confining interface and (QD8) a thermionic interface are refused with a `quantum_confinement mlda:` message, the same from the equilibrium and drift-diffusion assemblers; accepted on silicon with another parameter changed, and with Lombardi. |
| QC12 | 2D MOSFET example: the nearest silicon node is 5.00 nm from the oxide (its surface spacing), where f_n = 1 − 6e-8 and f_p = 1 − 2.5e-4, so the electron current Id(1.5 V) is unchanged to 6 digits (ratio 1.000000); sweep 0.69 s off, 0.68 s on. 3D example: 4.16 s off and on, Id(1.2 V) ratio 1.000000. The examples' meshes do not resolve the 1–2 nm confinement layer, so the model has no effect on them. FinFET deck (45,375 nodes, nearest silicon node 0.100 nm from the oxide): distance search 0.064 s. FinFET sweep with MLDA (Vd = 50 mV, the deck's Fermi–Dirac statistics so QD5's form, PARDISO on 6 threads, `line_search = natural`; before Unit 50 the deck's `residual` test stopped at 0.35 V, see limitations): 6.7 min, no step cuts. Id against the classical run: ratio 0.83 at 0 V, 0.79 at 0.2 V, 0.77 at 0.35 V (the least), 0.92 at 0.5 V, 0.974 at 0.8 V (on-current 1.835e-5 A against 1.883e-5 A). The threshold at a constant 1 µA rises by about 10 mV (0.257 V against 0.247 V, interpolated in ln Id). LS6 with MLDA and Lombardi: the 1e16 gap is +5.2% at 3.4e5 V/cm, +13.5% at 7.4e5, +18.7% at 9.0e5 and +21.6% at 9.9e5 (without MLDA: +5.4% at 3.5e5, +14.4% at 7.5e5, +20.0% at 9.2e5, and 21.6% at 9.6e5 with 0.25 V steps). The 1e17 gap is +0.8% at 3.2e5 and +10.5% at 9.7e5 (without MLDA: +0.8% at 3.3e5, +11.7% at 9.8e5). The 1e17 sweep stops at the default tolerance (limitations below); its gaps are from a rerun at `tol_update = 1e-7`. MLDA lowers the gap by about 1 point, so it does not explain LS6's 1e16 miss: the Jensen bias remains. |
| QC13 | Met. Release (115 tests, with `solve_mlda` and `deck_mlda`) and Debug (104) suites pass; zero warnings under `/W4 /WX` in both. One existing pin changed with the vocabulary, not loosened: the app's `[models]` completion on an empty line lists 19 keys (20 less `srh`), not 16. Hand mutants, 16 of 16 caught: the trap level shift dropped (QC4); the band shift, n_ie scaling, ln(N_c / n_ie) sign and ln(N_v / n_ie) shift each removed (QC4); the lumped-gate mean over the full cell instead of the half (QC3); the interface point at the insulator node instead of the edge midpoint (QC3); λ without the 2 (QC2); the sign of the mean's series (QC2); the electrothermal, silicon and thermionic refusals each removed (QC11; the electrothermal one survived a weak message check at first, and the test was tightened); the masses not recorded (QC1); the deck's `mlda_mass_p` written into `mlda_mass_n` (QC1); lumped-gate nodes not marked as half cells (QC9); the model not applied in the scaled device (QC7). |

**Limitations found (for the owner):**
- **Two Newton stalls, diagnosed (2026-10-24, owner: "yes" to finding the node and variable).** Neither is in
  the model. A scratch probe restarted Newton from a saved converged point at the next bias and listed the
  largest scaled corrections by node and unknown (`DriftDiffusion::update_size`'s measure).
  - **The FinFET: the residual line search at a rounding residual.** From point 6 (0.30 V) with the gate 0.2 mV
    higher (the run's last cut, 2^-8 of the step), full Newton steps converge in 6 iterations: corrections 3.0e2,
    8.1, 0.15, 1.4e-3, 3.4e-7, 2.8e-9. From iteration 5 the residual is at rounding (4.5e-13), while a correction
    of 3.4e-7 is still left on holes of about 1e-3 cm⁻³ in the n+ source and drain next to the oxide, below the
    density floor. The deck's `line_search = residual` Armijo test sees only noise there, rejects the step,
    backtracks to 1/256 and never finishes the iteration. The run with the residual test was the one reported
    first: it stopped at Vg = 0.35 V after 104 min. A rerun at `tol_update = 1e-7` was worse (304 iterations by
    point 1) and was stopped. **The same deck with `line_search = natural` runs all 17 points in 6.7 min, 8–15
    iterations per point, with no cuts.** Its currents equal those of the residual-test run at the points both
    solved, to 5 digits.
  - **The 1e17 long channel with Lombardi: the minority-hole floor in uncontacted n+.** At the step from 1.0 to
    1.5 V the corrections settle at 1.7–3.8e-8 with the residual at rounding (3.3e-11), at full steps. They sit on
    minority holes (p 10–20 cm⁻³) in the n+ pockets under the oxide, above the contacts' boxes (y 1.99–2.0 µm),
    uniform down whole columns: those holes are nearly decoupled, so rounding sets their correction. **With MLDA
    off the same step stalls at 1.4–2.2e-8**, so the floor exists with Lombardi alone. Lombardi-only sweeps pass
    the 1e-8 test or not by margin, and MLDA moves the floor a little.
  - **Resolved by Unit 50** (`main` 0c14841, PR #71, 2026-10-28), the separate unit proposed here: (a) the residual
    line search takes the full step when the residual is already at its rounding level; (b) the opt-in
    `minority_floor`. The FinFET deck's advice to use `line_search = natural` with MLDA, first given here, is
    withdrawn. Rebased on `main`, the FinFET with MLDA and the deck's `line_search = residual`: all 17 points in
    11.4 min (PARDISO, 6 threads), 8–12 iterations per point, no step cuts; Id(0.8 V) 1.8349e-5 A against the
    natural run's 1.835e-5 A (2026-10-28).
- **QC5's n p relation on the drift-diffusion state** (above): 3.26e-12 against 1e-12, rounding of a minority
  density of 1e-4 cm⁻³.
- **The examples' meshes do not resolve the layer.** The 2D and 3D MOSFET examples' surface spacing (5 nm) is
  far coarser than λ, so the model changes their currents by less than 1e-6. The FinFET's 0.1 nm surface spacing
  does resolve it.
- **Corners, thin bodies, Fermi–Dirac and heterojunctions** are as designed (QD2, QD5): one nearest interface,
  Hänsch's DOS form under Fermi–Dirac, and no confinement at heterojunctions or contacts.

