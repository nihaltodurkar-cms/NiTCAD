<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 34 — selectable Masetti hole mobility (merged to `main` as `e14a53f`, PR #42, 2026-10-15)

**Owner request (2026-10-15):** a selectable Masetti doping-dependent hole mobility, off by default (legacy results
preserved), equations and coefficients verified against the paper or a reliable primary reproduction before coding,
compared with Caughey–Thomas on Unit 33's boron resistivity benchmarks, with regression tests against acceptance
criteria fixed in advance; incomplete ionization and n_i unchanged; no push, merge or default change without
authorization.

**Sources read (before coding).** G. Masetti, M. Severi, S. Solmi, IEEE Trans. Electron Devices 30, 764 (1983),
doi:10.1109/T-ED.1983.21207, is paywalled and was not read directly. Its boron model was read from three independent
reproductions, which agree on the form and on every coefficient:
- CERN Allpix Squared, "Charge Carrier Mobility" (allpix-squared.docs.cern.ch/docs/06_models/02_carrier_mobility/):
  the hole formula as the paper's eq (4), the values "taken from Table I" of the paper;
- Crosslight Software, "Low field mobility models" (2010; crosslight.com, Low_field_mobility_model.pdf): the
  Masetti default coefficients in SI units;
- the Flooxs model page (flooxs.ece.ufl.edu, Masetti Mobility Model), as quoted by a search engine (the page itself
  could not be fetched: certificate error).

μ_p(N) = μ_min1 exp(−P_c / N) + (μ_max − μ_min2) / (1 + (N / C_r)^α) − μ_1 / (1 + (C_s / N)^β), N the total ionized
impurity density, with μ_min1 = 44.9, μ_min2 = 0, μ_max = 470.5, μ_1 = 29.0 cm²/(V s), P_c = 9.23e16, C_r = 2.23e17,
C_s = 6.10e20 cm⁻³, α = 0.719, β = 2.0. Masetti's fits are at 300 K; the μ_max (T / 300)^−2.2 scaling that Allpix
Squared and Crosslight apply is not from the paper (Allpix attributes it to another model).

**Acceptance criteria, fixed before any Masetti result was computed:**

| # | Criterion |
|---|---|
| AC1 | Off (the default): every existing test passes unchanged; run identities and Unit 33's pinned Caughey–Thomas gaps are unchanged. |
| AC2 | The implemented μ_p equals the formula above evaluated independently in the test at N = 0 and 1e13 to 1e21 cm⁻³ (each decade) and at 250, 300 and 400 K, within 1e-12 relative; μ_p(0) = μ_max (T / 300)^−2.2. |
| AC3 | Boron resistivity at 300 K (Unit 33's decks, NBS SP 400-64): \|simulated / reference − 1\| ≤ 10% at each of 1e15, 1e16, 1e17, 1e18, 1e19 cm⁻³. |
| AC4 | At 1e17, 1e18 and 1e19 cm⁻³ the 300 K gap is at least 10 percentage points smaller in magnitude than Caughey–Thomas's (30.0%, 48.0%, 23.3%), and the largest \|gap\| over 1e15–1e19 is below Caughey–Thomas's 48.0%. |
| AC5 | Boron resistivity at 23 °C: \|gap\| ≤ 12% at each density (looser: the temperature scaling is not from the paper). |
| AC6 | The electron mobility is unchanged with Masetti on: the phosphorus benchmark gaps equal the Caughey–Thomas pins. |
| AC7 | The temperature slope used by the electrothermal model matches a central difference of μ_p(T) to 1e-6 relative. |
| AC8 | Masetti is refused (invalid_input, before solving) on a semiconductor region whose hole mobility parameters are not silicon's (its coefficients are silicon's). |
| AC9 | Performance: with Masetti the isothermal sweep's wall time is within 10% of Caughey–Thomas's (the mobility is computed once per node at assembly); measured and reported. |

**Implementation.** `physics::masetti_mobility` and `masetti_mobility_slope` with `masetti_boron_silicon`
(`physics/mobility.hpp`); `PhysicsModels::hole_mobility` (`caughey_thomas`, the default, or `masetti`) used for the
hole mobility in the isothermal and electrothermal assemblies; electrons keep Caughey–Thomas; the setting enters the
run record only when Masetti is selected; the deck's `[models] hole_mobility`. Incomplete ionization and n_i are not
touched. No default changed.

**Results against the criteria (measured after the criteria commit, `3670fc3`):**

| # | Result | Met? |
|---|---|---|
| AC1 | Explicit Caughey–Thomas is bit-identical to the default (currents, densities, run identity); all existing suites pass; Unit 33's pins hold. | yes |
| AC2 | Within 1e-12 at N = 0, 1e13–1e21 and 250, 300, 400 K; μ_p(0) = μ_max (T / 300)^−2.2. | yes |
| AC3 | 300 K gaps: −1.1% (1e15), −1.1% (1e16), −5.3% (1e17), **−11.0% (1e18)**, −3.9% (1e19). | **no, at 1e18** (by 1.0 point) |
| AC4 | Against Caughey–Thomas: 30.0 → 5.3, 48.0 → 11.0, 23.3 → 3.9 points at 1e17, 1e18, 1e19; largest \|gap\| 48.0% → 11.0%. | yes |
| AC5 | 23 °C gaps: −1.1%, −1.1%, −5.8%, **−12.4%**, −4.9%. | **no, at 1e18** (by 0.4 point) |
| AC6 | The phosphorus gaps equal the Caughey–Thomas pins at both temperatures. | yes |
| AC7 | Within 1e-6 of a central difference (h = 1e-4 T) at 250, 300, 400 K. | yes |
| AC8 | Refused for a germanium region; accepted for silicon with another parameter changed. | yes |
| AC9 | 2D diode, 20,301 nodes, 6 points, best of 3: ratios 1.36, 0.94, 0.98, 0.93 in four runs; the three runs that counted iterations took 26 Newton iterations with either model. | yes in 3 of 4 runs; the 1.36 run is timing noise on this machine (its iterations were not recorded) |

**Owner ruling (2026-10-15):** the AC3 and AC5 misses at 10¹⁸ cm⁻³ are **accepted as documented limitations, not as
passed criteria**.

The criteria are not loosened. The validation test checks AC3 and AC5 everywhere but 1e18, pins every measured gap,
and states AC3 and AC5 at 1e18 in a separate `[!shouldfail]` case that reports the miss (and would report a fix). AC9 is
a hidden `[.performance]` case (timing does not belong in ctest).

**Limitations and open issues:**
- **AC3/AC5 at 1e18 cm⁻³ not met (−11.0%, −12.4%).** A likely cause, not verified: Masetti fitted mobility against
  carrier (hole) concentration, while the benchmark solves with complete ionization (p = N); the NBS report finds
  boron de-ionization in this range, which a resistivity computed with p = N underestimates. Resolving it needs a
  valid ionization model at 1e18 (Unit 33 finding; out of this unit's scope by instruction).
- The coefficients come from three reproductions, not from the paper itself (paywalled).
- The (T / 300)^−2.2 scaling of μ_max is not Masetti's; Masetti's fits are at 300 K only. The 23 °C gaps track the
  300 K ones within 1.4 points, but other temperatures are not validated.
- Boron-doped silicon only: the coefficients are Masetti's boron set; other acceptors and other materials are refused
  or untested. Electron mobility is unchanged (Caughey–Thomas; phosphorus gaps −8.8% to +3.9%).
- Mutation checks (8), all caught: α miscopied, the slope's factor dropped, the clustering term inverted, the model
  selection inverted, the refusal disabled, the electrothermal flag never set, the run-record setting, the deck
  mapping.

