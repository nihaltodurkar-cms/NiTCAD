<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 33 — validation against published data (merged to `main` as `a14ac96`, PR #41, 2026-10-15)

**Owner rulings (2026-10-15)**, asked before any work: references are values stated in the text or tables of
published sources (no digitized curves); the benchmarks are decks in a new top-level `benchmarks/` directory, run by a
Release-only ctest against their cited values; a disagreement is reported and pinned, not tuned (no model parameter
changes in this unit); web search and fetch may be used to find and read the sources.

| # | Choice | Why | Trade-off |
|---|---|---|---|
| V1 | Resistivity against dopant density from W. R. Thurber, R. L. Mattis, Y. M. Liu, J. J. Filliben, NBS Special Publication 400-64 (1981), doi:10.6028/NBS.SP.400-64, the basis of ASTM F723. Phosphorus: eq (1), log₁₀(qρN / 1 V s cm⁻²) = (A₀ + A₁X + A₂X² + A₃X³) / (1 + B₁X + B₂X² + B₃X³), X = log₁₀(N / 10¹⁶ cm⁻³), table 9 (300 K: A = −3.0652, 2.1853, −0.61080, 0.056189; B = −0.67642, 0.19542, −0.018100; 23 °C: A = −3.0769, 2.2108, −0.62272, 0.057501; B = −0.68157, 0.19833, −0.018376). Boron: eq (4), qρN = min + (max − min) / (1 + (N / N_ref)^α), table 12 (300 K: 0.00215, 0.02053 V s cm⁻², 4.09e18 cm⁻³, −0.727; 23 °C: 0.00209, 0.02024, 3.88e18, −0.737). q = 1.602e-19 C, as the report says to use with its fits. | A measured, standardized, device-level relation: the simulator solves a resistor and the resistivity follows from its current. | It tests mobility (and ionization) at equilibrium-like fields only. |
| V2 | The source was read from the Internet Archive's OCR text of the scanned report (the NIST PDF is images; the ASTM preview stops before the equations). The OCR dropped the + signs of eq (1) and the exponents of N₀ and of table 12's N_ref: N₀ = 10¹⁶ gives the report's 0.52 Ω cm at 10¹⁶ and makes its two phosphorus fits agree; N_ref = 4.09e18 (3.88e18 at 23 °C) is the only decade that keeps eq (4) within the report's stated 3% of eq (3) (table 11, legible: 0.00220, 0.01973, ρ_ref 0.01782 Ω cm, 1.086 at 300 K). The test checks every transcribed fit against the report's other fit within its stated self-consistency (boron 3%, phosphorus 4.5%). | No number is taken from memory; a miscopied coefficient fails the check (mutation-tested). | Two exponents are inferred, not read; the inference is recorded here. |
| V3 | Eleven decks (`benchmarks/resistivity/`): phosphorus 1e14–1e19, boron 1e15–1e19 cm⁻³ (inside the report's measured range), a 10 µm bar at 300 K, 0 and 0.01 V; ρ = 0.01 V / (J L). Each also at 23 °C (the deck's temperature overridden by the test) and with incomplete ionization. | One deck per device, readable and runnable with `nitcad run`. | Decks have no variables, so each density is its own file. |
| V4 | The intrinsic density from K. Misiakos, D. Tsamakis, J. Appl. Phys. 74, 3293 (1993), doi:10.1063/1.354551: "(9.7 +/- 0.1) X 10(9) cm-3" at 300 K (abstract, read from its open-archives record). The model's is √(np) in the solved 1e14 bar at 0 V. | A measured value with a stated uncertainty, tested through a solved device. | — |
| V5 | Each gap (simulated / reference − 1) is pinned at its measured value within 0.2 percentage points; no parameter is changed (owner ruling). | A model change must show up as a changed benchmark, deliberately. | The tests pass with known disagreements; the findings below say where. |
| V6 | Not included: abrupt-junction breakdown voltages (S. M. Sze, G. Gibbons, Appl. Phys. Lett. 8, 111 (1966)) and MOS C-V landmarks. The Sze–Gibbons values could not be read from an accessible copy; textbook C-V landmarks are analytic, not measured. | Only values actually read are used. | Breakdown stays validated against the analytic one-sided value only (Unit 19). |

Findings (pinned gaps, simulated against the reference):
- **Phosphorus, 300 K:** +3.9% (1e14), +0.5% (1e15), −6.9% (1e16), −8.8% (1e17), +0.6% (1e18), −7.1% (1e19). The 23 °C
  gaps are within 1 point of these (the temperature dependence holds).
- **Boron, 300 K:** −4.4% (1e15), +3.4% (1e16), **+30.0% (1e17), +48.0% (1e18), +23.3% (1e19)**. The legacy hole
  mobility is the Caughey–Thomas fit to Irvin's curve; the report itself finds that "the Caughey and Thomas
  expression is a poor fit to the mobility data" and differs from the older curves by up to 45% near 5e17 for boron.
  A boron mobility from this report's own fits would close the gap; that is a model change (an owner decision).
- **Incomplete ionization** (a fixed 45 meV level, Unit 15) makes the heavily doped bars much worse: phosphorus +31%
  (1e18) and +138% (1e19), boron +121% (1e18) and +302% (1e19). The level does not lower with density (the report
  finds de-ionization only at intermediate densities, and full activity at the highest), so the model should not be
  used for ionization above about 1e17 cm⁻³ until it has a density-dependent level.
- **Intrinsic density:** the legacy silicon n_i, 1.0674e10 cm⁻³, is 10.0% above the measured 9.7e9 (9.7 times its
  uncertainty); currents set by n_i² (diode saturation, generation) are about 21% high by that alone.
- Mutation checks (4): a phosphorus coefficient, the boron N_ref decade, a 23 °C ρ-fit exponent (caught by the
  transcription check), the temperature override; all caught.

**Scope (owner instruction, 2026-10-15): Unit 33 is a validation baseline only.** It changes no physics, model,
parameter or solver code; the findings above are recorded for later units, each an owner decision.

**Reconstructed coefficients (read from the OCR text where legible; reconstructed where not):**

| Quantity | Printed (OCR) | Used | How it was fixed |
|---|---|---|---|
| Eq (1) operators | "Aq - A^X . A^X . A3X" | A₀ + A₁X + A₂X² + A₃X³ over 1 + B₁X + B₂X² + B₃X³ | With + signs the ρ- and N-based phosphorus fits agree within 1.3% at 10¹⁶ cm⁻³ (the report's stated self-consistency) |
| N₀ (eq 1) | "Nq =10^^ cm~^" | 10¹⁶ cm⁻³ | Gives 0.52 Ω cm at 10¹⁶ cm⁻³ phosphorus, consistent with both fits |
| Boron N_ref, 300 K (table 12) | "4.09 x 10 + 0.20 x 10" | 4.09 × 10¹⁸ cm⁻³ | At ρ = ρ_ref of eq (3) (N = 3.84e18), eq (4) agrees within 1.5% only for 10¹⁸; 10¹⁷ is 60% off. The report states the two fits agree within 3% from 10¹⁴ to 10¹⁹ |
| Boron N_ref, 23 °C (table 12) | "3.88 x 10 + 0.21 x 10" | 3.88 × 10¹⁸ cm⁻³ | The same test at 23 °C: 1.8% for 10¹⁸ |
| All other coefficients (tables 9, 11, 12) | legible | as printed | — |

The test checks every fit used against the report's other fit within the stated self-consistency (boron 3%,
phosphorus 4.5%), so an error in a reconstruction would fail it.

**Unsupported regimes (do not rely on the model there; none is fixed by this unit):**
- **Boron (p-type) mobility at 10¹⁷–10¹⁹ cm⁻³:** resistivity 23–48% high (the legacy Caughey–Thomas hole mobility).
- **Incomplete ionization above about 10¹⁷ cm⁻³:** resistivity 31% to 302% high at 10¹⁸–10¹⁹ (a fixed 45 meV level).
- **Absolute n_i-driven currents:** n_i(300 K) 10% high against the measurement, so n_i² currents (diode saturation,
  generation) about 21% high.
- **Outside the benchmarked range:** below 10¹⁴ (phosphorus) or 10¹⁵ (boron) and above 10¹⁹ cm⁻³; temperatures
  other than 300 K and 23 °C; compensated material; dopants other than phosphorus and boron (the report's fits
  apply to them only approximately). The report notes its phosphorus conversion understates density above about
  5 × 10¹⁹ cm⁻³.
- **Not validated against published data at all:** breakdown (only the analytic one-sided value, Unit 19), MOS C-V,
  high-field transport (the benchmarks use 10 V/cm), transient and small-signal behaviour.

