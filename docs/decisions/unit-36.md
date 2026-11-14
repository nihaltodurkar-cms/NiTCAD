<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 36 — trace robustness at the current-resolution floor (merged to `main` as `07cd53a`, PR #46, 2026-10-16)

**Owner ruling (2026-10-16):** Units 36–39 in order (trace robustness, density-dependent ionization, Masetti electron
mobility, breakdown validation), each its own branch, criteria before code and PR; stop for the go-ahead to merge
after each; Units 37 and 38 opt-in.

**The two failures (Unit 35, unresolved):** on the defaults, the 5e19 zener's nonlocal band-to-band `trace_bias` stalls
at −1.31 V, and the open-base transistor's snap-back trace crosses its fold in about 5000 steps.

**Diagnosis (measured with temporary instrumentation, removed before this commit).** In both, the corrector's largest
update is the swept bias λ itself (every potential on the swept side moves with it), alternating in sign at 1e-7 to
1e-5 with the residual flat: λ is not resolved. The arc-length row c is M τ, and its current term is the gradient of
log₁₀(|I| + I0) built from the swept contact's current row (`small_signal_currents`), whose entries are the contact's
majority fluxes' derivatives: they cancel as the contact current does (its Unit 32 resolution). Measured on the zener:
Σ|c_q x_q| up to 7e8 against a row value of 1e-4, a cancellation over 13 orders, so both the arc residual
c·(y − y₀) − ds and the Sherman–Morrison products u·z and u·w of the bordered solve are rounding. Checks: with the current
term removed the zener corrector converges quadratically (4–5 iterations, updates 1e-14) and reaches −2 V; with the
products computed in compensated arithmetic it reaches −2 V with the current term kept. The snap-back with compensated
products crosses the fold in a few steps but then accepts a corrected point on the other branch (16.87 V, 7.9e-4
A/cm² → 13.33 V, 9.2e-2 A/cm²; the step's arc length was 1 V) and runs the curve backwards; the tangent-turn test does
not catch it.

| # | Choice | Why | Trade-off |
|---|---|---|---|
| R1 | Every product with the arc row (the arc residual, u·z and u·w) is a compensated dot product (Ogita, Rump, Oishi, SIAM J. Sci. Comput. 26, 1955 (2005), Dot2: TwoProduct by FMA, TwoSum), as accurate as in twice the working precision. | The residual and the solve then agree with c to the rounding of the result, not of Σ|c_q x_q|; c is ours (any fixed row defines a valid hyperplane), so only consistency matters. | About four flops per entry more, once per iteration (n + 1 entries). |
| R2 | A corrected point farther from its predicted point than the step itself, measured in the arc metric, is rejected and the step retried at half (a distance test of the kind in Allgower and Georg, *Numerical Continuation Methods*, 1990). | The turn test missed a 3.5 V jump to the other branch; a corrector that converges so far from its prediction is not following the branch. | A legitimate large correction is retried smaller. |
| R3 | No iterative refinement of the bordered solve. | Measured with R1 and R2: Newton iterations over the zener trace 250 → 169, over the snap-back 901 → 879; not worth the extra solves. | — |

**Acceptance criteria, fixed before the implementation (the diagnosis above was measured with temporary code):**

| # | Criterion |
|---|---|
| TR1 | The zener nonlocal trace to −2 V on the defaults completes, its end point the swept solve's within the resolution and 1e-6 (the existing assertions); Unit 35's `[!shouldfail]` case becomes an ordinary case. |
| TR2 | The open-base snap-back on the defaults passes every assertion of the existing case (fold, turning-point window, the high-current branch to 20 A/cm², the bias falling past the fold, α_T M_n); the case returns to the defaults. |
| TR3 | Both cases also pass on the legacy options (legacy silicon, Caughey–Thomas holes), as before. |
| TR4 | Every other test passes with its tolerance unchanged (the traces of the breakdown, band-to-band, electrothermal and weighted-current suites among them). |
| TR5 | Cost: the snap-back on the defaults takes at most 200 points (Unit 35: about 5000); points and wall times of both cases are reported. |
| TR6 | The compensated dot product is exact to the rounding of its result on a cancelling pair (terms of 1e16 summing to 1); the distance test and R1 are mutation-checked (each removed, a TR1/TR2 case fails). |
| TR7 | Full Release and Debug suites; zero warnings. |

**Implementation.** `linalg::compensated_dot` (`include/NiTCAD/linalg/compensated_dot.hpp`, Dot2 with FMA; needs the
build's `/fp:precise`); `trace_bias` computes the arc residual and the Sherman–Morrison products u·z, u·w with it (R1)
and rejects a corrected point farther than ds from its prediction in the arc metric (R2; not on the step that lands on
end_V, whose row is the bias itself). `trace.hpp` documents both.

**Results against the criteria (measured after the criteria commit, `7ba51ef`):**

| # | Result | Met? |
|---|---|---|
| TR1 | The zener trace reaches −2 V on the defaults and on Caughey–Thomas holes, its end point the swept solve's within the resolution and 1e-6; the `[!shouldfail]` case is removed, the ordinary case runs both. | yes |
| TR2 | The snap-back on the defaults passes every assertion: fold 17.4625 V, the turning point at it, the high-current branch to 20 A/cm², α_T M_n within 1e-3, every point a steady state. | yes |
| TR3 | On the legacy options: fold 17.2666 V (Unit 19's 17.27), every assertion passes. | yes |
| TR4 | Every other test passes with its tolerance unchanged. One test changed: the electrothermal runaway fold read the fold from the parabola through three traced points, an O(h³) estimate that depended on where the top point fell (Unit 36's corrector takes fewer, differently placed points: 1.4e-7 V against the 1e-7 tolerance, 0.13 K from the fold; before, 4.6e-9 at 0.05 K). It now maximizes the quartic through five points (O(h⁵)): 7.3e-10 V, and −2.7e-10 V on the code before Unit 36; the tolerance is unchanged. | yes |
| TR5 | Snap-back: 59 points in 0.7 s on the defaults (Unit 35: about 5000), 55 points in 0.8 s on the legacy options. Zener: about 20 step attempts, under a second, on each option set. | yes |
| TR6 | `compensated_dot` returns 1 for 1e16 + 1 − 1e16 (the plain sum 0) and 2⁻⁶⁰ for (1 + 2⁻³⁰)² − 1 − 2⁻²⁹. Mutation checks (5), all caught: the arc residual plain, the Sherman–Morrison products plain, the distance test off, Dot2's product error dropped, its sum error dropped. | yes |
| TR7 | Release 18/18, Debug 12/12; no warnings. | yes |

**Limitations.** The current term of the arc metric still comes from the contact current; only its arithmetic is
compensated. A current below its resolution still cannot steer a step (I0, unchanged). The distance test uses the arc
metric at the step's start.

