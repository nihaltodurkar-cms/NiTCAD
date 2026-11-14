<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 35 — published silicon defaults: band set and hole mobility (merged to `main` as `2d1fdc7`, PR #44, 2026-10-16)

**Owner authorization and rulings (2026-10-15)**, asked before any work: change the defaults. (1) n_i from a
published, physically consistent silicon band set (Eg(T), Nc(T), Nv(T)), verified from the source before coding — not
a rescale of Nc/Nv or an Eg0 shift; (2) the default hole mobility becomes "automatic": Masetti in regions with
silicon's hole mobility parameters, the material's Caughey–Thomas elsewhere (heterostructures keep working); (3) the
legacy stays selectable — the legacy silicon set (n_i 1.0674e10 cm⁻³) and `hole_mobility = caughey_thomas` as named
options reproducing old results bit for bit; (4) legacy-port gates keep running on the explicitly selected legacy
options; new-default gates are added with their before/after values recorded and justified against Unit 33's
benchmarks. Incomplete ionization is not touched.

**Source read (before coding).** R. Couderc, M. Amara, M. Lemiti, "Reassessment of the intrinsic carrier density
temperature dependence in crystalline silicon", J. Appl. Phys. 115, 093705 (2014), doi:10.1063/1.4867776, open
access (hal.science/hal-00956479), read from the PDF's own text streams. The extraction drops minus signs (plus signs
survive as a separate glyph); each sign below was fixed by the paper's own cross-checks, recorded in the last column.

| Quantity | As used | Paper | Check |
|---|---|---|---|
| n_i | n_i² = Nc Nv exp(−Eg(T) / kT) | eq (1) | — |
| Nc, Nv | 4.83e15 (m_dc / m0)^(3/2) T^(3/2), likewise Nv [cm⁻³] | eqs (4), (5) (the constant as printed; CODATA gives 4.829e15) | — |
| (m_dc / m0)^(3/2) | A T³ + B T² + C T + D, A = −4.609e-10, B = 6.753e-7, C = −1.312e-5, D = 1.094 | eq (18), 0–400 K | It fits 6 m_t √m_l / m0^(3/2) with m_l = 0.9163, m_t = 0.1905 Eg(0) / Eg(T) (eqs 16, 17) within 0.5% (the paper); these signs give 0.03% from 20 to 400 K, and D = 1.094 = 6 (0.1905) √0.9163 |
| (m_dv / m0)^(3/2) | A T³ + B T² + C T + D, A = 2.525e-9, B = −4.689e-6, C = 3.376e-3, D = 0.3426 | eq (19), from eq (15) | With eq (18) and Eg(T) the resulting n_i is within 1.5% of eq (15) from 50 to 400 K (the paper); these signs give 1.3% (largest, 100 K), 0.12% at 300 K |
| Eg(T) | Eg(0) − αΘ [ (1 − 3Δ²) / (exp(Θ / T) − 1) + (3Δ² / 2) ( ⁶√(1 + π² χ² / (3 (1 + Δ²)) + (3Δ² − 1) χ³ / 4 + 8 χ⁴ / 3 + χ⁶) − 1 ) ], χ = 2T / Θ | eq (10), R. Pässler, Phys. Rev. B 66, 085201 (2002), the paper's chosen gap model | π² / (3(1 + Δ²)) (not π² · 3(1 + Δ²)): the mass fit above holds to 0.03% only so (2.4% otherwise); prefactor αΘ: gives the renormalization αΘ / 2 the paper states for Pässler's model and Eg(300 K) = 1.1242 eV |
| Pässler's parameters | Eg(0) = 1.17 eV, α = 3.23e-4 eV/K, Θ = 446 K, Δ = 0.51 | Table I | — |
| Reference n_i(T) | 1.541e15 T^1.712 exp(−Eg(T) / 2kT) | eq (15); "9.68 × 10⁹ cm⁻³" at 300 K | Reproduced: 9.685e9 |

Resulting silicon values at 300 K: Eg = 1.12416 eV (legacy 1.12452), Nc = 2.857e19 (legacy 2.86e19), Nv = 2.514e19 cm⁻³
(legacy 3.10e19), n_i = 9.673e9 cm⁻³ (legacy 1.0674e10; Misiakos 9.7 ± 0.1e9, Unit 33 V4).

**Design (decided here, before coding):**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| S1 | Two optional parameter groups on `SemiconductorParameters`: a Pässler gap (α, Θ, Δ, with `Eg0_eV` as Eg(0); the Varshni α and β then 0) and density-of-states mass polynomials (conduction and valence; Nc300 and Nv300 then 0). All zero (every existing set): the legacy Varshni and (T / 300)^(3/2) formulas, evaluated exactly as before. | A published set in the material's own data; legacy sets bit for bit; the forms are general, not a silicon switch. | Two ways to state a band; `create` refuses a mixture (both forms given for one quantity). |
| S2 | `physics::silicon_parameters` (deck `silicon`) becomes Couderc's set; every other field is unchanged. The old set is `physics::legacy_silicon_parameters` (`legacy_silicon()`, deck `silicon_legacy`). | Ruling (1), (3). | Code and decks that named silicon change results (the point of the unit). |
| S3 | `HoleMobility::automatic`, the new default: per semiconductor region, Masetti where the region's hole mobility parameters equal silicon's Caughey–Thomas set, the material's Caughey–Thomas elsewhere. `caughey_thomas` (legacy) and `masetti` (refuses non-silicon regions, Unit 34) stay. The run record names the setting unless it is `caughey_thomas`, so legacy-option runs keep their identity. | Ruling (2), (3). | "Silicon" is recognized by its hole mobility parameters (as Unit 34's refusal). |
| S4 | Electrothermal: with temperature-dependent DOS masses, N(T) = N(T0) θ^(3/2) M(θ). The extra −kT ∇ln M in Wachutka's flux is a band-parameter step across an edge: ln M(θ_b) − ln M(θ_a) joins the edge's ln Nc0 step (electrons; −ln Nv for holes), with its θ partials; the densities at fixed reduced energy move with 3 / (2θ) + d ln M / dθ in place of 3 / (2θ). The energy a carrier carries, (r + 5/2) kT above its band edge, is unchanged (an assumption: the thermopower is that of parabolic bands at the local mass). With no mass polynomial, M ≡ 1 and every term is the legacy's, bit for bit. | The default material must keep working electrothermally; refusing it or freezing M at T0 would be a silent approximation (T8's rule). | The thermopower assumption above is recorded, not validated. |
| S5 | The p-polysilicon work function of an `electrode` contact (no semiconductor under it) uses `silicon_parameters`, now Couderc's gap: 0.36 meV below the legacy at 300 K. | One silicon default. | A legacy run with a p-poly electrode contact differs by that shift (gates take the gap of the semiconductor under them, so a legacy-silicon device is bit for bit). |

**Acceptance criteria, fixed before any implementation:**

| # | Criterion |
|---|---|
| DC1 | Legacy selectable, bit for bit: with `silicon_legacy` / `legacy_silicon_parameters` and `hole_mobility = caughey_thomas`, a diode sweep and the boron 1e18 benchmark written by `nitcad run` are byte-identical to the same decks run on `main` (`0768cae`, the deck's `silicon` then), and every legacy-port gate passes on those options with its unchanged tolerance. |
| DC2 | The implemented Eg(T), Nc(T), Nv(T) equal the formulas above evaluated independently in the test (from this table, not the parameter struct) at 50, 77, 100, 200, 300, 400 K within 1e-12 relative; legacy sets unchanged (Varshni and power law, bit for bit). |
| DC3 | n_i(300 K) of the default silicon in [9.6e9, 9.8e9] cm⁻³ (Misiakos's 9.7 ± 0.1e9) and within 1.5% of eq (15) from 50 to 400 K (the paper's stated precision). |
| DC4 | Unit 33's solved n_i benchmark (√(np) in the 1e14 bar at 0 V) with the defaults: \|gap\| ≤ 1.1% against 9.7e9 (within the measurement's uncertainty); the legacy pin (+10.0%) unchanged on legacy options. |
| DC5 | Unit 33/34 resistivity benchmarks with the defaults: boron gaps equal Unit 34's Masetti pins and phosphorus gaps equal the Caughey–Thomas pins, each within 0.2 points (a majority-carrier resistivity does not depend on n_i); the largest boron \|gap\| ≤ 11.2%. |
| DC6 | dEg/dT, d ln Nc / dT, d ln Nv / dT and d ln n_i / dT used by the electrothermal model match central differences within 1e-6 relative at 100, 300 and 400 K (both band forms). |
| DC7 | Electrothermal with the defaults: the existing Jacobian (finite-difference) and energy-balance checks pass on the new silicon; with isothermal sinks at 300 K and a small current, the electrothermal current matches the isothermal one within 1e-3. |
| DC8 | `automatic`: Masetti in a silicon region and Caughey–Thomas in a germanium region of one device (no refusal), each region's hole current matching a run with the explicit model on a single-material device; `masetti` still refuses the germanium region. |
| DC9 | New-default gates: each legacy-port gate's device is also run on the defaults, the value recorded here before/after with its explanation; a diffusion-limited diode current changes by the n_i² ratio (9.673e9 / 1.0674e10)² = 0.821 within 2% (hole mobility held at Caughey–Thomas to isolate it). |
| DC10 | Full Release and Debug suites after each commit; zero new warnings. |

**Implementation.** `physics::PasslerGapParameters` and `DensityOfStatesMassParameters` on `SemiconductorParameters`
(S1; `create` refuses a mixed or incomplete form; `check_temperature` refuses a mass fit that is not positive, silicon's
conduction cubic near 2030 K); `silicon_parameters` is Couderc's set and `legacy_silicon_parameters` /
`legacy_silicon()` the old one (S2); `physics::dos_mass_log` and `dos_log_slope`, and the gap, n_i and
emission-velocity slopes of both forms (`thermal.hpp`); `HoleMobility::automatic` (value 2; the default) with
`assemble::masetti_holes` deciding per material (S3); the electrothermal mass term (S4: `ThermalState::mass_n/p`, the
edges' ln M step, the levels, the ohmic contacts and the interfaces); device-file keys 56–58 and 60–67 (absent: the
legacy forms, so a file written before Unit 35 loads as legacy silicon); the run identity digests the new forms only
when present; the deck's `material = silicon_legacy`, `passler_gap.*` overrides and `hole_mobility = automatic`.
Incomplete ionization is not changed.

**Results against the criteria (measured after the criteria commit, `f250d9c`):**

| # | Result | Met? |
|---|---|---|
| DC1 | `nitcad run` on `silicon_legacy` with `hole_mobility = caughey_thomas`: the diode sweep and the boron 1e18 bar's `.sweep` and `.csv` files are byte-identical to `main`'s (`0768cae`) on the same decks. The legacy-port gates run on the legacy options and pass with their tolerances: J(0.5 V) = 1.2801e-2 (bias, sweep, field mobility, heavy doping), the M33 heterojunction set, the MOS capacitor port and the MOSFET Id–Vg and mesh gates, the legacy diode turn-off, Unit 33's pins and Unit 34's Masetti pins. | yes |
| DC2 | Eg, Nc, Nv and n_i within 1e-12 of the independent evaluation at 50, 77, 100, 200, 300, 400 K; the legacy sets' formulas unchanged (their double-double reference tests pass). | yes |
| DC3 | n_i(300 K) = 9.673e9 cm⁻³; within 1.3% of eq (15) from 50 to 400 K (largest near 100 K). | yes |
| DC4 | Solved n_i 9.6731e9: −0.3% (0.3 times the uncertainty); the legacy pin +10.0% unchanged. | yes |
| DC5 | Boron at 300 K: −1.1%, −1.1%, −5.3%, −11.0%, −3.9%; 23 °C: −1.1%, −1.1%, −5.8%, −12.4%, −4.9% (Unit 34's Masetti pins); phosphorus unchanged; largest boron gap 11.0%. | yes |
| DC6 | Within 1e-6 of central differences at 100, 300, 400 K for both forms (and the emission velocity with and without a Richardson constant). | yes |
| DC7 | On the default silicon: the finite-difference Jacobian gates (steady and time step), the energy identity and the uniform-T equilibrium pass; the Seebeck voltage (300–310 K, both silicon sets) within 1e-5 and the 300–400 K refinement case within 1e-5 of ∫(2 − ln(N_D / Nc(T))) dT with the material's own Nc(T); the Peltier coefficient within 1e-7 at Couderc's Nc(300 K); the isothermal limit within 1e-3. | yes |
| DC8 | A silicon/germanium diode solves with `automatic`, `masetti` refuses it; on silicon `automatic` equals `masetti` bit for bit, on germanium `caughey_thomas`. | yes |
| DC9 | See the before/after table; the n_i² ratio of a diffusion current (Caughey–Thomas holes): 0.8219 at 0.4 V and 0.8217 at 0.5 V against 0.8213 (+0.08%, +0.05%). | yes |
| DC10 | Release and Debug: 18/18 suites; no warnings. | yes |

**Before (legacy options) and after (the defaults):**

| Gate | Before | After | Why |
|---|---|---|---|
| Unit 9 diode, J(0.5 V) | 1.2801e-2 A/cm² | 1.1399e-2 (−11.0%) | n_i² × 0.821; the hole part (23%) up by Masetti's 319 against 232 cm²/(V s) in the 1e17 n side |
| n_i(300 K) | 1.0674e10 cm⁻³ | 9.673e9 | Couderc's set (Nv 3.10e19 → 2.514e19, Nc 2.86e19 → 2.857e19, Eg 1.12452 → 1.12416 eV) |
| MOSFET V_th (max g_m) | 0.1520 V | 0.1527 V | — |
| MOSFET swing | 68.79 mV/decade | 68.61 | — |
| MOSFET Id(Vg = 1 V) | 0.31471 A/cm | 0.31448 (−0.1%) | an electron channel current |
| MOSFET on/off | 1.459e11 | 1.303e11 | the off current, junction leakage, 12% higher (not decomposed) |
| Boron resistivity, largest gap | +48.0% | −11.0% | Masetti (Unit 34) |
| Incomplete ionization, boron 1e18 / 1e19 | +120.9% / +301.7% | +40.3% / +241.3% | Masetti and Nv (the acceptor occupancy); still the fixed 45 meV level (Unit 33 finding) |

Tests that compare against analytic forms written with the legacy Caughey–Thomas hole mobility or the T^(3/2) law now
take the default's (the MOS RC and dielectric-relaxation references, the surface-recombination diode, the field-
mobility resistor, the gate's intrinsic offset) or the material's Nc(T) (Seebeck, Peltier); legacy reproductions run on
the legacy options. Two comparisons were relaxed to the quantity's resolution, each explained where it is: a sweep or a
restart against a chain of one-point solves compares contact currents within their Unit 32 resolution (with the lower
n_i the cancelling majority fluxes round the current by 6e-8 of it at 0.4 V), and a transient snapshot's time is
matched to 1e-12 (it is one ulp off 3 τ through the scaled time). The mid-gap-trap test allows the rows' rounding
(the rate is 2e-8 of rows rounding at 1e-16).

**Unresolved issues (found by this unit, not fixed):**
- **Nonlocal band-to-band trace on the defaults** (resolved by Unit 36: the arc row's rounding). The 5e19 zener's `trace_bias` to −2 V stalls at −1.31 V: the
  corrector's updates stay between 3e-8 and 2.5e-6 against the 1e-8 tolerance with the residual flat, a rounding floor
  of the lagged-path system; the swept solves converge, and with Caughey–Thomas holes the trace completes. The
  mechanics test runs on Caughey–Thomas; a `[!shouldfail]` case records the default.
- **Open-base snap-back trace on the defaults** (resolved by Unit 36). The fold (17.46 V, 1.61e-3 A/cm²) is crossed in about 5000 steps,
  the current moving by its resolution (2.9e-7 A/cm²); the curve and BV_CEO are reached, but the highest traced bias
  falls past the turning-point window. The case runs on the legacy options (not repeated as `[!shouldfail]`: four
  minutes).
- The p-polysilicon work function of an `electrode` contact takes the default silicon's gap (S5): 0.36 meV below the
  legacy's at 300 K, so a legacy run with such an electrode is not byte for byte (none of the legacy gates has one).
- The thermopower of a temperature-dependent mass (S4) is an assumption, checked for consistency (Seebeck, Peltier,
  energy identity), not against a measurement.
- Couderc's fits are stated for 0–400 K (n_i to 1.5% over 50–400 K); above 400 K they are evaluated as written.

**Performance.** The band forms are evaluated once per node at assembly (isothermal) or once per node per assembly
(electrothermal: two polynomials and logarithms per band); not separately timed.

Mutation checks (10), all caught: Pässler's prefactor halved, the mass polynomial's slope missing its constant, the
electron and the hole ln M step's sign, the level's d ln M / dθ dropped, `automatic` inverted, the run identity
skipping the valence fit, the device file dropping it, the deck's `silicon_legacy` mapping, the n_i slope's mass term.

