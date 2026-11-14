<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 46 — Lombardi surface mobility, opt-in (started 2026-10-24; on `physics/surface-mobility`)

**Owner request:** "start Lombardi surface-roughness mobility" (2026-10-24), the first of the nanoscale models of
ARCHITECTURE.md 11.3 (found by the FinFET deck). **Owner rulings (2026-10-24):** (1) E⊥ is the field normal to the
current, (2) the Jacobian is exact, with a wider pattern, (3) the electrothermal assembly is included, (4) the surface
terms combine with the run's selected bulk mobility. The standing conditions hold: opt-in and bit-identical when off,
criteria fixed before code, Release and Debug suites at every commit, `/W4 /WX`, hand mutants, nothing pushed or
merged without authorization.

**Sources read (before coding).** C. Lombardi, S. Manzini, A. Saporito, M. Vanzi, IEEE Trans. CAD 7, 1164 (1988),
doi:10.1109/43.9186, is paywalled and was not read. The model is read from two text reproductions that cite it:
(a) the Xyce Reference Guide (Sandia; `doc/Reference_Guide/Xyce_RG_PDE.tex` and `lombardiTbl.tex` on
github.com/Xyce/Xyce), (b) Genius TCAD Open (`src/material/Si/Si_mob_Lombardi.cc` on
github.com/cogenda/Genius-TCAD-Open, Medici's `.LSM` parameter names). Both give

    1/mu = 1/mu_b + 1/mu_ac + 1/mu_sr,
    mu_ac = B / E⊥ + C N^lambda / (T E⊥^(1/3)),      mu_sr = delta / E⊥^2,

N the total impurity density [cm^-3], T [K], E⊥ [V/cm]. COMSOL's Semiconductor Module guide (doc.comsol.com,
"Lombardi Surface Mobility Model") cites the paper for its silicon values; its equations are images, not read.

| Coefficient | Xyce | Genius | COMSOL (as extracted) | Used |
|---|---|---|---|---|
| B_n [cm/s] | 4.75e7 | 4.75e7 | 4.75e7 | 4.75e7 |
| C_n [cm^(5/3+3λ) K V^(-2/3) s^-1] | 1.74e5 | 1.74e5 | 1.74e5 | 1.74e5 |
| λ_n | 0.125 | 0.125 | 0.125 | 0.125 |
| δ_n [V/s] | 5.82e14 | 5.82e14 | 5.82e14 | 5.82e14 |
| B_p [cm/s] | 9.93e6 | 9.93e6 | 9.93e7 (an extraction artefact or a misprint; two sources disagree) | 9.93e6 |
| C_p | 8.84e5 | 8.84e5 | 8.84e5 | 8.84e5 |
| λ_p | 0.0317 | 0.0317 | 0.0317 | 0.0317 |
| δ_p [V/s] | 2.05e14 | 2.05e14 | 2.05e14 | 2.05e14 |

The E⊥ exponent of μ_sr is 2 in both (Xyce `exn8`, `exp8`). Genius applies Canali's saturation to the Matthiessen
sum (the low-field mobility includes the surface terms), as NiTCAD will. Lombardi's own bulk term is Masetti's
arsenic electron set and boron hole set; by ruling (4) it is not used.

**Design:**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| LD1 | `PhysicsModels::surface_mobility`: `off` (the default) or `lombardi`; the deck's `[models] surface_mobility`; the run record names it only when `lombardi`; the app's deck vocabulary lists it. | Opt-in, as Units 34, 37, 38. | — |
| LD2 | Per edge, the field vector is the mean of the two ends' least-squares fields (the reconstruction of impact ionization: E_i = P_i Σ_k t_k E_k over the node's carrier edges, so only semiconductor fields, never an oxide's), and E⊥ = sqrt(\|E − (E·t)t\|^2 + E_0^2) with t the edge's unit vector and E_0 = 1 V/cm. | Ruling (1); exact for a uniform field at any angle; E_0 keeps B / E⊥ and δ / E⊥^2 finite and smooth where E⊥ vanishes (Genius adds 1 V/cm). | Bulk edges at an angle to a junction field are reduced too. In 1D E⊥ = E_0 (the change is below 1e-4). |
| LD3 | Each end's low-field mobility is the Matthiessen sum of its bulk mobility (the run's selection, at its N and T) and the surface terms at its own N and T with the edge's E⊥; the edge mobility is their harmonic mean, then Canali when on (on a mixed edge each end's Canali, then the mean), as today. N is the bulk model's N (0 with `doping_mobility` off). | Ruling (4); the existing harmonic mean and Canali order unchanged. | — |
| LD4 | Exact Jacobian: E⊥ is linear in ψ before the root, so each edge stores the coefficients g_j (D-vectors, t-projected out) over the nodes j of both ends' stencils; the carrier rows of a and b gain ψ columns of those nodes (distance 2), added to the pattern only with the model on. | Ruling (2); quadratic Newton and exact small-signal (Units 13, 19). | More nonzeros; measured and reported (LS11). |
| LD5 | Electrothermal: the per-end mobility and its θ partial include the surface terms (∂μ_ac/∂T = −C N^λ / (T^2 E⊥^(1/3))); E⊥ is from ψ, which does not depend on θ. | Ruling (3). | — |
| LD6 | `lombardi` is refused (invalid_input) on a device with a semiconductor region whose electron or hole Caughey–Thomas parameters are not silicon's. | The coefficients are silicon's; as Units 34, 38. | SiGe and others not offered. |
| LD7 | Not on thermionic edges (no mobility enters them). | As Canali. | — |

**Acceptance criteria, fixed before any implementation:**

| # | Criterion |
|---|---|
| LS1 | Off (the default): every existing test and pin passes unchanged; run identities and Jacobian patterns unchanged. |
| LS2 | μ_ac, μ_sr and the Matthiessen sum equal the formula and coefficients above, evaluated independently in the test, at E⊥ = 1e3 … 1e7 V/cm (each decade), N = 0, 1e14 … 1e20 cm⁻³ and T = 250, 300, 400 K, within 1e-12; the E⊥ and T slopes within 1e-6 of central differences. |
| LS3 | E⊥ of a uniform field on 2D meshes (tensor and Delaunay) at 0, 30, 45 and 60° and on a 3D tensor mesh: each interior edge's E⊥ equals sqrt(\|E × t\|^2 + E_0^2) within 1e-10 relative. |
| LS4 | The Jacobian against central differences of the residual with the model on (isothermal and electrothermal; 2D and 3D; with and without Canali): the largest relative column error at most 1e-6. Newton on the 2D MOSFET in strong inversion converges quadratically (the last three correction norms). |
| LS5 | 1D diode: the current with `lombardi` within 1e-4 relative of the current without (E⊥ = E_0). |
| LS6 | A long-channel n-MOSFET (Release only), Vd = 10 mV, substrate 1e16 and 1e17 cm⁻³: the effective mobility μ_eff = (L / W) g_d / Q_inv at mid-channel, against μ of LD3 evaluated at E_eff = q (Q_b + Q_inv / 2) / ε_Si and the substrate's N, within 15% for E_eff in 3e5 … 1e6 V/cm. |
| LS7 | Universality: in LS6's devices, μ_eff of the two substrates at equal E_eff (interpolated) within 15% for E_eff in 5e5 … 1e6 V/cm. |
| LS8 | The 2D MOSFET example's Id–Vg (Vd = 0.1 V): at Vg = 1.5 V the current with `lombardi` lies below the current without; below threshold (Vg = 0.2 V) within 10%. The ratio at 1.5 V is pinned. |
| LS9 | The 2D MOSFET example rotated 90° (x and y exchanged): the same Id at Vg = 0.5 and 1.5 V within 1e-9 relative. |
| LS10 | Refused on a germanium region and accepted on silicon with another parameter changed; with electrothermal and 300 K sinks within 1e-3 of the isothermal current (LS8's device, Vg = 1.5 V); small-signal conductance against DC differences within 1e-5. |
| LS11 | Reported, not gated: the pattern's nonzeros and the 2D and 3D MOSFET examples' sweep times with and without the model; the FinFET deck's on-current with and without (the deck's mesh convergence is still open). |
| LS12 | Full Release and Debug suites; zero warnings; hand mutants caught. |

**As built (2026-10-24).** `physics/surface_mobility` (the formula and its partials),
`DriftDiffusion` (the per-edge stencils of LD4, the surface path of `edge_fluxes` and `thermal_fluxes`, the
small-signal current rows), the `CarrierFluxes` and `ThermalEdges` terms (the pattern and the E⊥ partials, in the
heat rows too with electrothermal), the deck key, the run record and the app's vocabulary. Tests:
`tests/physics/surface_mobility_test.cpp` (LS2), `tests/assemble/surface_mobility_test.cpp` (LS1, LS3, LS4, LS10),
`tests/deck/surface_mobility_test.cpp` (LS5, and Release only under `[lombardi]`, ctest `deck_lombardi`: LS4's
Newton, LS6–LS11).

| # | Added in the build | Why |
|---|---|---|
| LD8 | `lombardi` with `impact_ionization` is refused (not supported). | Impact ionization's generation rows use the fluxes' partials over its own pattern; the E⊥ partials reach distance 2, outside it. Combining them needs that pattern widened; not done in this unit. **For the owner.** |

**Results:**

| # | Result |
|---|---|
| LS1 | Met. With the model off the pattern and residual are identical (tested); every existing test and pin passes. |
| LS2 | Met: 1,821 assertions (formula within 1e-12, slopes within 1e-6). |
| LS3 | Met on 2D tensor, 2D Delaunay and 3D tensor meshes at every angle. |
| LS4 | Met. The FD gate holds within 1e-6 in all nine cases (2D, Delaunay, 3D; ±Canali; electrothermal 2D, 2D Canali, 3D Canali). On the 2D MOSFET at Vg = 1.5 V the last three corrections are 1.9e-1, 2.5e-5 and 4.3e-15, an order of 2.5. |
| LS5 | Met: the diode within 1e-4, and the run record names the model only when on. |
| LS6 | **1e17 met; 1e16 not met.** Measured at Vg 0 … 6.5 V by 0.5 V on a mesh with 0.05 nm spacing in the top 5 nm. 1e17: gap +0.8% at 3.3e5 V/cm rising to +11.7% at 9.8e5. 1e16: +5.4% at 3.5e5, +14.4% at 7.5e5, **+20.0% at 9.2e5** (21.6% at 9.6e5 with 0.25 V steps). Without the model the same extraction gives the bulk mobility within 1% (checked), so the extraction is sound. Cause: μ_eff is a current-weighted average of μ(E⊥) over the classical inversion layer, whose field falls from E_s at the interface into the bulk. μ is convex in E⊥ (Jensen), so the average exceeds μ(E_eff). The gap grows with field and is larger at the lighter substrate, where Q_inv/2 dominates E_eff. Encoded as a `[!shouldfail]` case with the criterion unchanged. **For the owner.** |
| LS7 | Met: worst 13.2% at 5e5 V/cm, 11.5% at 1e6. |
| LS8 | **Partly not met.** At 1.5 V the current is lower: 4.5348e-2 A/cm against 6.5274e-2, ratio 0.694738, pinned within 1e-4. **Below threshold (0.2 V) the ratio is 0.7249, not within 10%.** On the same device with 0.05 nm spacing at the surface it is 0.718, so the mesh does not cause it. Cause: in depletion the surface field is about 4e4 V/cm at 0.2 V, and there μ_ac (B/E⊥ ≈ 1,100; the C term ≈ 1,660 at 1e16) already lowers μ from 1,250 to about 430 at the interface. The subthreshold current flows in that layer, and the model has no distance damping. The criterion assumed a small effect; the model says otherwise. Encoded as a `[!shouldfail]` case. **For the owner.** |
| LS9 | Met: the reflected device's Id within 1e-9 at 0.5 and 1.5 V. |
| LS10 | Met. Refused on germanium, accepted on silicon with another parameter changed; electrothermal with 300 K sinks 0.999936 of the isothermal current; small-signal g_d and g_m against DC differences within 1e-8. Added in the build, as the mutants found: the small-signal current rows against central differences of the terminal currents (2D, 3D Canali, 2D electrothermal Canali), within 1e-6 of the largest entry. On the example the drain contact lies in the n+ pocket, where E⊥ is small, so the device-level check alone could not catch a missing E⊥ partial in those rows. |
| LS11 | 2D MOSFET example: Jacobian nonzeros 30,244 off, 37,284 with the model (+23%); sweep 0.94 s off, 1.5 s with it. 3D MOSFET example: 5.6 s off, 12.0 s on; Id(1.2 V) ratio 0.718. FinFET deck (PARDISO, 6 threads): on-current at 0.8 V 1.6586e-5 A with the model, 1.8835e-5 without (ratio 0.881); at 0.2 V 0.971. 30 min with the model against 18 min without, but the with-model time was measured while other builds ran, so it is not a clean comparison. |
| LS12 | Release (110 tests) and Debug (101) suites pass, zero warnings. One existing pin changed with the vocabulary, not loosened: the app's `[models]` completion on an empty line lists 16 keys (17 less `srh`), not 15. Hand mutants, all caught: the E⊥ projection removed (LS3); the electrothermal T term of the θ partial (LS4); the heat-row E⊥ partials (LS4); the thermal continuity E⊥ partials (LS4); the isothermal and thermal small-signal E⊥ terms (LS10's current rows; both survived the device-level check before those rows were added); the germanium refusal (LS10); the run-record entry (LS5); E_0 = 0 (a contract violation at E⊥ = 0 in the 1D diode). |

**Limitations found (for the owner):**
- **Mesh dependence.** In classical drift-diffusion the inversion layer is about kT/(qE_s) thick, 0.26 nm at
  1e6 V/cm. On a mesh with 1 nm surface spacing the LS6 gap reached +58%; it converged near 0.05 nm. The examples'
  and the FinFET's surface spacings (about 10 nm) under-resolve the layer, so their Lombardi currents are
  indicative. LS8's subthreshold ratio, which depends on a depletion layer many nm thick, is not affected. A
  follow-up could take E⊥ at the interface node from the interface's D-field (Gauss) instead of the
  least-squares nodal fields.
- **No distance damping.** The model applies at every carrier edge with E⊥ from the field normal to the current,
  not only near an interface. Bulk edges at an angle to a junction field are reduced too (LD2). Commercial tools
  often damp the surface terms with distance from the interface. Not in the owner's rulings; not done.

**Ruling on LS6 (1e16) and LS8 (below threshold) (2026-10-24).** The owner delegated it ("Decide yourself").
Both are accepted as limitations of the model as specified, not defects. The criteria are unchanged, and the two
`[!shouldfail]` cases stay, so a change that makes either pass is flagged.
- LS8's below-threshold limit assumed a small effect. Lombardi without distance damping lowers the mobility in the
  depletion field at the surface, and the result has converged on the mesh.
- LS6's 1e16 gap is the Jensen bias of comparing μ_eff with μ(E_eff) across the classical inversion layer. The
  extraction is verified, and 1e17 and LS7 pass. Quantum confinement (MLDA, ARCHITECTURE 11.3) narrows the layer's
  field spread, so it is the place to revisit this.
- LD8 (refused with impact ionization) and the mesh dependence stand as documented limitations.

