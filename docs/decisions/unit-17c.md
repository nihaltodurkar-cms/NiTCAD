<!-- Archived from DECISIONS.md on 2026-11-17 (branch docs/archive-17-26d), word for word; DECISIONS.md keeps its summary. -->

### Unit 17c — adaptation for derived goals: V_T, g_m and swing (approved 2026-11-13; on `mesh/derived-goals`)

**Request.** Owner, 2026-11-13: "draft the design for the V_T, g_m and swing refinement". This is the
follow-up that Unit 17's AR8 deferred: "a derived goal needs its own adjoint (a functional of several points) and
resolution rule: a later unit". **Approved** by the owner (2026-11-13: "approved, choose recommended"). Branch `mesh/derived-goals`, stacked on 17b (`fad06b4`; 17b is PR #82).

**What exists.**
- `analysis/dc.hpp` extracts each quantity from a transfer curve:
  - `threshold_max_gm`: tangent at the largest g_m, less V_d / 2;
  - `max_transconductance`: central derivative, with the point and its neighbours as the window;
  - `subthreshold_swing`: the steepest segment below 1e-2 of the largest current.
- Every adaptation pass already solves the whole sweep on one mesh (AQ3).
- AR3's adjoint gives, for each adaptation point k, an estimate D_k of the error in I_D(k) and per-node weights ω_k.
- AD11 verifies on M′, the uniform refinement of the final mesh M.

A derived goal G = F(I_1 … I_N) depends on the mesh only through the currents in its window, so the whole unit is
the chain rule plus a rule for when the window moves.

| ID | Construction | Basis | Gain | Risk / cost |
|---|---|---|---|---|
| DG1 | **Goals.** `goal = vt <contact>`, `goal = gm <contact>` and `goal = swing <contact>` select the extractors above, applied to the goal contact's current against the swept gate bias. They need a `[sweep]` that varies exactly one gate or electrode contact, holds every other contact fixed, and has at least 3 points (swing: 2 in its window). V_T's drain bias is the goal contact's fixed bias. A falling curve (p-channel) is scaled by −1, as dc.hpp documents. Current goals may be mixed with derived goals. | dc.hpp, Unit 24 | The figures a MOSFET user reports become adaptation goals. | — |
| DG2 | **Window and adaptation points.** On each pass, the extractor's window on that pass's sweep is: for g_m and V_T, the max-g_m point and its two neighbours; for swing, the steepest segment's two points. Those window points are the goal's adaptation points, and `at = auto` takes their union with the current goals' points. An explicit `at` combined with a derived goal is refused, because the window decides. | — | The adjoint is solved only where the goal depends on the current. | — |
| DG3 | **Sensitivities.** c_k = ∂F/∂I_k at the pass's currents, with the window held fixed, computed analytically per extractor. g_m: ±1/(V_{k+1} − V_{k−1}) on the neighbours. V_T: through the tangent, V* − I*/g_m − V_d/2. Swing: through ΔV / Δlog10 I. A finite-difference check of each formula is a test (DC2). | Chain rule | No new adjoint system: the existing per-point adjoints are combined. | Valid only while the window does not move (DG6). |
| DG4 | **Estimate and weights.** The estimate is E_G = Σ_k \|c_k\| · E_k, where E_k is AR3's estimate for I_D(k); without weighting, E_k is the last change. Summing absolute values is the triangle-inequality bound on the linearization, so errors at different points cannot cancel in the estimate. Per-node weights are ω = Σ_k \|c_k\| ω_k (each ω_k already normalized, AD3), renormalized to sum to 1, and AD3's marking, ID4's benefit and ID6's allocation use them unchanged. | Dual-weighted residual for a functional of several solutions | The loop is driven by the goal, not by every current in the sweep. | — |
| DG5 | **Allowance (AD13 for derived goals).** V_T has an absolute allowance, `vt_tolerance` in V (default 1e-3). Its value can be near 0, so a relative tolerance is meaningless there. g_m and swing are positive and use `tolerance` relative, tol · max(\|G\|, r_G). The resolution is r_G = Σ_k \|c_k\| r_k, where r_k is the current's resolution (Unit 15). | AD13; Unit 15 | Each figure's tolerance is stated in units that suit it. | One more key. |
| DG6 | **Moving window.** The extractors are piecewise: the argmax point or steepest segment can move to a neighbour, and G jumps when it does. A stop is proposed only when the window is the same on the last two passes, and the certificate requires the same window on M and M′. If the window differs, the result is **not_certified** with the reason "window moved", and that counts as a rejected stop (AR2: continue, marked by the verified difference). | — | A certificate is never issued across a discontinuity. | A device whose g_m peak is flat can end not_certified (reported, never hidden). |
| DG7 | **Verification on M′ (AD11).** M′ solves the window's points plus one neighbour on each side, so a local move of the argmax or steepest segment is seen on M′ too. Richardson extrapolation is applied to G itself, as AD11 does: G_M′ is F applied to M′'s currents, and the order is observed from the loop's passes. | AD11 | An independent check that is still cheaper than solving the whole sweep on M′ (about 5 points against N). | A window move further than one point away on M′ is not seen there; it is seen across the loop's passes (DG6). |
| DG8 | **What the certificate covers.** Only the mesh's discretization error, for the given gate step. The gate step's own error in V_T, g_m and swing (sampling) is not certified. `nitcad run` and the report state this next to each derived certificate. | — | Honest scope. | Users may read "certified V_T" as also covering the step; the text says it does not. |
| DG9 | **Report, CSV and identity.** The goal kind is added to the report (new key, written only for derived goals), `adapt.goal` (the kind and the contact) and the CSV columns `<kind> <c> certificate`, `estimate`, `verified error`, `allowed` and `order`. Units: V, S cm^(D−3) and mV/dec. The window's points are reported per pass. With current goals only, the report, CSV and identity are bit-identical to today's. | AD8 | — | — |
| DG10 | **Learning (17b).** Histories record the goal kind, and schema 1 is unchanged: kind ≠ 0 is a category the model has not seen, so LB8 already sends the device to ID2. Training on derived-goal histories is a later schema. | LB3, LB8 | No new code path for the predictor. | — |
| DG11 | **Unchanged.** The adjoint's refusals (`check_adjoint`), AD11 itself (independent of all predictions, OA6), AD12's fallback, the mesh kinds (AQ2) and one mesh per sweep (AQ3) are all as in Unit 17. | — | — | — |
| DG12 | **Layering.** The layer table forbids `solve` to depend on `analysis` (L6 above L5), so the three figures are restated in `solve/adapt_derived.cpp` from dc.cpp and curve.cpp (derivative, `peak_transconductance`, `threshold_max_gm`, `subthreshold_swing` at its default ceiling); a test holds them equal to `analysis/dc.hpp`'s, bit for bit, on n- and p-channel curves (with DC2). | ARCHITECTURE.md 4 | No layering change. | Two copies of three short formulas; the test catches drift. |

**Acceptance criteria. They are frozen before any run and never loosened after one; a miss gets `[!shouldfail]` with
its measurement and goes to the owner for a ruling.**

| ID | Criterion |
|---|---|
| DC1 | Nothing changes. Every existing test passes, and every deck with current goals only (or without `[adapt]`) has a byte-identical report, CSV and identity. |
| DC2 | **Sensitivities and estimate.** On analytic transfer curves (exponential below threshold, linear above), each extractor's c_k equals a central finite difference of the extractor with the window fixed, to 1e-6 relative. E_G and ω equal values computed by hand on a fixture. |
| DC3 | **Accuracy gate.** Device: the tensor MOSFET (AR6: V_d 0.05 V, V_g from −0.2 to 1.0 V in 0.05 V steps, 25 points), adjoint-weighted, each goal adapted separately: `vt` at 1 mV, and `gm` and `swing` at 1e-2. Each run ends **certified**. Its actual error against the reference is within its allowance. The reference is the extraction on a 4×/8×/16× uniform Richardson ladder, recorded as a constant with a 2× guard solve, as AR4 did. |
| DC4 | **Efficiency (ID9: the goal earns its place).** Baseline: `goal = current drain` on the same window points, with the tolerance tightened through 1e-2, 3e-3, 1e-3 and 3e-4 until its V_T, g_m and swing meet DC3's allowances. The derived-goal run must use ≤ 0.8× the baseline's final DOFs for V_T and for swing, and ≤ 1/3 of the uniform ladder's nodes for the same error. |
| DC5 | **Estimate quality.** At the final pass of each DC3 run, E_G is within a factor 3 of the actual error (AC7's bound). Both E_G and the verified error are reported, with their ratios to the actual error. |
| DC6 | **Moving window.** In a fixture where the max-g_m point moves between M and M′, the result is not_certified with "window moved", never certified. With `certify = required`, the run exits with code 3. |
| DC7 | **Deck and refusals.** `nitcad check` accepts `vt`, `gm` and `swing` and refuses each of the following with its line: a sweep with no single swept gate or electrode contact; fewer than 3 points; an explicit `at` combined with a derived goal; `vt_tolerance` without a `vt` goal; a non-ohmic goal contact with adjoint weighting. |
| DC8 | **Determinism and independence.** Two runs are bit-identical. AC21 is extended: on a forced common mesh, a derived certificate is identical with `policy = basic`, with `policy = intelligent` and with a `predictor`. |
| DC9 | **Hand mutants (AC12).** At least one per row DG1–DG10 and DG12; each survivor is closed by a test or argued equivalent. |
| DC10 | **Suites.** Release and Debug pass at every commit with `/W4 /WX`; Debug grows by ≤ 60 s; the DC3, DC4 and DC5 gates are `[.derived_goals]`. |

**Owner questions, ruled as recommended (owner, 2026-11-13: "approved, choose recommended"): DQ1 (a), DQ2 (a), DQ3 (a), DQ4 (a).**
- **DQ1 (V_T's allowance).**
  - (a) **Recommended:** absolute, `vt_tolerance` in V, default 1 mV.
  - (b) Relative, like the current goals.
  - Reason for (a): V_T can be near 0 V, where a relative allowance is meaningless.
- **DQ2 (M′'s points).**
  - (a) **Recommended:** the window plus one neighbour on each side (DG7).
  - (b) The whole sweep, which also catches far window moves on M′ at about N/5 times the verification cost (25 points against 5 in DC3).
- **DQ3 (p-channel).**
  - (a) **Recommended:** supported by dc.hpp's −1 scaling.
  - (b) Refused in this unit.
- **DQ4 (DC4's baseline).**
  - (a) **Recommended:** current-goal adaptation tightened until it meets the same derived accuracy.
  - (b) Uniform refinement only, as AC5 did.
  - Reason for (a): the question this unit answers is whether a derived goal beats adapting the currents.

**Slices (planned).**
1. DG1–DG5 in `solve::adapt`: sensitivities, estimate, weights, allowance. Tests DC2 and DC1.
2. DG6–DG8 (window rule, verification, certificate scope), plus DC6.
3. DG9 and DG10: deck, report, CSV, identity, app columns. Tests DC7 and DC8.
4. Gates DC3–DC5, then mutants (DC9) and the suites (DC10).

**Owner instruction on tests (2026-11-13).** "Dont start the 17c tests", then "keep writing the code, no tests". The
code is written first. No 17c test is written or run until the owner says so. DC1–DC10 stay as frozen above; none is
measured yet.

**Slices 1–3 as written (code only; compiled in Release, `nitcad` and `nitcad_app`, under `/W4 /WX`; nothing run).**
- **DG1–DG3, DG12.** `solve/adapt_derived.{hpp,cpp}`: `GoalKind`, `TransferCurve`, and `derived_goal`, which
  restates the three extractors and their analytic sensitivities. The p-channel sign rule: a curve whose largest
  |I| is negative is read with currents and biases times −1. dV_T/dI equals the read curve's dG/dy; g_m's and the
  swing's sensitivities are s · dG/dy. `AdaptGoal::kind` and `AdaptOptions::vt_tolerance` (default 1e-3 V) are added.
- **`adapt()`, DG2.** Refused:
  - an explicit `at` with a derived goal;
  - fewer than 3 points;
  - a sweep that moves anything but exactly one contact;
  - a swept contact that is not a gate or an electrode;
  - a derived goal on the swept contact.

  The current goals keep `at`; their defaults apply only when a current goal exists, and ID3's coarse sweep runs
  only then. The entries are the current goals × `at` (goal-major), then one per derived goal. Indicators and
  differences unite `at` with the windows. Without derived goals the points are `at` itself, the same order as
  before, so DC1 holds by construction.
- **DG4.** One `estimate_goals` call per window point covers the derived goals' contacts there. The estimate is
  Σ |c_k| E_k. The weights Σ |c_k| ω_k are renormalized over the semiconductor nodes, and ω is the largest over all
  goals.
- **DG5.** `allowed_of` gives a vt entry its absolute `vt_tolerance`; every other entry keeps
  `allowed_error(tol, G, r)`. It is used by `relative`, the stop rule, AD11 and the certificates.
- **DG6.** A derived goal's stop needs its window unchanged from the previous pass. On M′ a moved window, or a
  failed extraction, rejects the stop. The certificate is then `not_certified` with `window_moved`.
- **DG7.**
  - With derived goals and no current goal, M′ solves the windows plus one point either side, from the first of
    them; `solve_on` and `difference_field` take that offset.
  - The derived values on M′ are extracted from M's curve with M′'s currents at the solved points. The extraction
    runs on that whole curve, so the swing's ceiling (1e-2 of the largest current) sees every point.
  - With a current goal present, M′ solves the whole sweep as before.
- **DG8, DG9.**
  - `nitcad run` prints a derived certificate as "<kind> of <contact> (window at point k)", in its unit, with
    the scope sentence and any moved window.
  - The CSV columns are "<kind> <contact> certificate", "estimate <unit>" and so on; a current goal's header is
    unchanged.
  - The app's certificate lines show the kind and unit.
  - The `.nitadapt` report gains keys only when a derived goal is present: report 29 `goal_kind` (texts) and 30
    `vt_tolerance`, pass 22 `window` (first and last per derived goal), certificate 4 `window_moved` (written only
    when true). The entries' count on load is (goals − derived) × points + derived.
  - `adapt.goal_kind` (per goal) and `adapt.vt_tolerance` enter the settings and identity only with a derived goal.
  - The app's `[adapt]` vocabulary gains vt, gm, swing and `vt_tolerance`. The `predictor` key's stale "refused in
    this unit" text is corrected.
- **DG10.**
  - Feature 12 is the run's goal kind: 0 with current goals only, else its first derived goal's (1 vt, 2 gm,
    3 swing). It is passed to `edge_features` and `HistoryRecorder`.
  - Schema 1 is unchanged. A model trained on current runs has feature 12's range [0, 0] (±5% of 0), so every edge
    of a derived-goal run is out of range and LB8 falls back to ID2. That is DG10, with no new code path.
- **Not yet:** ARCHITECTURE.md (6.12.4, the unit table, 11.3) and the tests (DC1–DC10, held per the owner's instruction).

**Tests allowed (owner, 2026-11-14: "yes" to "May I start the 17c tests?").** Slice 4 findings and rulings follow;
each ruling is the recommended option under the standing instruction (no criterion changed).

| # | Finding | Options | Ruling |
|---|---|---|---|
| DR1 | **DG4's absolute sum defeats DC3 and DC5.** Trial (Release, DC3's runs, `max_passes` 12): no goal certified. On the last pass (5358 nodes) V_T's estimate is 20.7 mV and g_m's 1.9e-2 A/cm·V, against actual errors of a few mV and about 0.5% (the 16× ladder); swing's is 22 mV/dec at 69.7. Σ\|c_k\|E_k loses the cancellation the figures are built on: g_m is a difference of neighbouring currents, so a common relative error ε gives ε·g_m but the bound ε(I_{k+1}+I_{k−1})/ΔV, about (V_g − V_T)/ΔV ≈ 7× more; V_T is invariant to a common scaling, and the swing to a common relative error. | (a) **Signed:** E_G = 2\|Σ c_k D_k\| from the adjoint's signed differences, weights \|Σ c_k ℓ_k(n)\| from the nodes' signed contributions (the linearized dual-weighted estimate of G itself). A cancelled estimate cannot certify alone: AD11's verified error on M′ must also be within the allowance. (b) Keep the bound: DC3 and DC5 fail. (c) The bound to certify, signed to mark. | **(a)**. `GoalEstimate::node_signed` added; DG4's text above is superseded for derived goals. |
| DR2 | **Feature 12 ignored when predicting.** `learned_targets` called `edge_features` without the goal kind, so a derived run's prediction read feature 12 as 0 (current) while the range check and the history used the run's kind. | (a) **Pass the kind** (`learned_targets(..., goal_kind)`, `run.cpp` passes it). (b) Leave it (DG10's fallback hides it for schema-1 models, not for a later schema). | **(a)**; DC8's predictor case is the regression test (its tree splits on feature 12). |
| DR3 | **DC3's uniform ladder is not asymptotic** (as AR9 on the same device). Uniform 1×–16× (275–61 985 nodes): V_T 0.08600, 0.11257, 0.13438, 0.14707, 0.15248 (order 1.23, extrapolated 0.15649, 4 mV beyond 16×); g_m 0.33709, 0.36610, 0.38795, 0.39172, 0.38989 (not monotone); swing 69.393, 69.279, 69.701, 70.126, 70.348 (order 0.94, 70.589). | (a) **Per-axis ladders, as AR9**, recorded with a 2× guard. (b) DC3's 4×/8×/16× as frozen: V_T's reference uncertain by 4× its allowance. (c) 32× uniform (247 k nodes, 25 points: hours). | **(a)**: DC3's reference method, not its allowances. The probe (`[.adapt_derived_ref2]`) solves only the windows' neighbourhoods (points 2..5 and 11..17); at 4×4 it reproduces the full sweep's figures to 12 digits. |
| DR4 | **The adjoint's D_k are not accurate enough for the swing on adapted meshes.** After DR1, g_m certified (3036 nodes; estimate 1.57e-3 against the verified 1.55e-3) but the swing's estimate stayed at 15.5 mV/dec from 3.9 k to 10 k nodes. Probe (`[.adapt_derived_dprobe]`, `_dprobe2`): uniform 4× → 8×, D_k/I at points 2, 3, 4 is −0.243, −0.254, −0.262 against actual −0.247, −0.256, −0.263; on the swing run's pass-6 mesh (55 × 34, x left coarse) +0.209, +0.031, −0.064 against actual −0.261, −0.269, −0.275. The actual errors are nearly common (the swing's error is small); the linearization's remainder is not, and the swing differences it. | (a) **Solved difference on M′:** a derived goal's estimate is AD11's Richardson error of G from the window ±1 points solved on M′ (DG7's set), each pass; the adjoint keeps the weights (DR1). (b) Adjoint only: the swing never certifies. (c) Refuse the swing with adjoint weighting. | **(a)**. Costs one Newton solve of 3–5 points on M′ per pass (the adjoint already factorizes on M′ at each window point). |
| DR5 | **DG7's substituted curve mixes meshes at its ends.** With DR4 the swing's estimate rose to 24.6 mV/dec and g_m's grew with refinement. M′'s currents replace M's only at the window ±1, so the segment (or central slope) that straddles the solved range's end mixes an M′ current with an M current about 26% higher. Such a segment looks steep, and the extractor picks it. The same defect was in AD11's verification (`verify_pass`), where it would also have shown as false window moves. | (a) **Search only inside the solved range:** `derived_goal(kind, curve, from, to)` admits only windows whose points all lie in [from, to]. Local moves (±1, DG7) are still seen; the defaults search the whole curve, so DG12's bit-equality is unchanged. (b) Solve the whole sweep on M′ (DQ2 (b)): about 5× the cost. | **(a)**, used for both the estimate (DR4) and AD11. A test (DR5) checks the range on g_m and the swing. |
| DR6 | **A false swing certificate on a pre-asymptotic mesh.** After DR5 the swing run certifies on pass 1 (338 nodes): 69.347 mV/dec, estimate = verified error 0.126, allowance 0.69. The reference is about 70.4–70.6, so the actual error is about 1.1 mV/dec. Uniform factors kx × ky: 1×1 69.393, 2×1 69.293, 1×2 69.373, 2×2 69.279, 4×2 69.302, 2×4 69.686, 4×4 69.701, 8×4 69.696, 4×8 70.130; y alone ky = 32, 64: 70.42, 70.44. The swing's y error appears only once the y step is below about 10 nm: every one-step difference on the coarse meshes (uniform or per axis) is ≤ 0.1, so no estimate from a single refinement sees it. V_T (12 passes, 6084 nodes: estimate 1.9 mV against 1 mV; the window flips 13..15 ↔ 14..16) has not certified; g_m certified at 3036 nodes (estimate 1.55e-3). | (a) Verify derived goals on two levels, M′ and M″ (2× and 4×, window ±1 only); 1× would then see 69.28 → 69.70. Heuristic: a threshold beyond two levels is still missed. (b) A resolution prerequisite before a derived stop (e.g. AD3's feature indicators without excess near the interface). (c) Keep the method and state the limitation: a derived certificate covers only meshes in the asymptotic range; DC3 fails for the swing. (d) Remove the swing (or its certification) from 17c. | **(a)**, owner 2026-11-15 ("Two-level verify (Recommended)"). AD11 for derived goals also solves the window ±1 on M″ (M′ refined uniformly again); the verified error is the larger of AD11's and the distance from M″'s three-level extrapolation (order from the two differences, clamped to [0.5, 3]); a window moved on M″ also rejects. The limitation (a threshold beyond two levels is still missed) goes to ARCHITECTURE 11.3. |
| DR7 | **The swing's marks follow the adjoint's unreliable weights.** With DR6 the false stop is rejected (pass 1: verified 0.95 against 0.13 single-level), but the swing then stays at 69.69–69.73 from 1.1 k to 9.5 k nodes (reference about 70.44), with its estimate at about 1.0: the marks come from the adjoint's weights, which DR4 showed are wrong in subthreshold. | (a) **A swing goal marks by its solved difference on M′** (DR4's solve, the window ±1; unweighted, as AR2 marks after a rejected stop); V_T and g_m keep the adjoint's weights (strong inversion: g_m certified with them). (b) All derived goals mark by difference. (c) Unchanged. | **(a)**. |

**DC3 trial after DR6–DR7 (Release, 2026-11-15; `[.adapt_derived_trial]`, `max_passes` 20, adjoint, tol 1e-2, vt 1 mV).**
- g_m: certified on pass 8 (3036 nodes): 0.385008, estimate = verified 1.55e-3, allowance 3.85e-3.
- Swing: certified on pass 4 (10 098 nodes): 70.128 mV/dec, estimate 0.24, verified 0.32, allowance 0.70; the false stop of DR6 (pass 1) rejected.
- V_T: certified on pass 13 (12 719 nodes): 153.549 mV, estimate = verified 0.87 mV, allowance 1 mV; window 14..16 from pass 11.
- Per-axis ladder (DR3) at ky = 32, kx = 4, 8, 16, 32: V_T 150.100, 152.504, 153.327, 153.851 mV; g_m 0.378112, 0.384683, 0.386714, 0.387635; swing 70.4217, 70.4178, 70.4170, 70.4168. At kx = 4, ky = 32, 64, 128: V_T 150.100, 150.158, 150.141; g_m 0.378112, 0.376874, 0.376498; swing 70.4217, 70.4400, 70.4446. The free window flips 13..15 → 14..16 between kx = 8 and 16, so V_T's x ladder mixes a jump with the error: a fixed-window ladder (`[.adapt_derived_ref3]`, window 14..16) is being measured for DC3's V_T reference.
- Fixed-window ladder (`[.adapt_derived_ref3]`, window 14..16 held): at ky = 32, kx = 4, 8, 16, 32: V_T 148.894, 152.158, 153.327, 153.851 mV (orders 1.48, 1.16); g_m 0.376850, 0.384335, 0.386714, 0.387635 (1.65, 1.37). At kx = 4, ky = 32, 64, 128: V_T 148.894, 148.936, 148.903 (within 0.04 mV, not monotone); g_m 0.376850, 0.375599, 0.375205. At kx = 8, ky = 64: V_T 152.142, g_m 0.382959 (the y step 32 → 64 is 1.10× that at kx = 4 for g_m).
- **DC3's references (DR3 (a)):** each axis extrapolated at its last observed order, the y correction scaled to the finer x: V_T 154.28 mV (uncertain by about 0.3 mV: the x order still drifts), g_m 0.38612 (about 0.13%), swing 70.441 (about 0.005). The 32 × 32 solve takes 16 minutes, so the values are recorded in the test; the guard re-solves the 2× uniform mesh (0.112573879392, 0.366102000996, 69.2786037546, to 1e-6).

| # | Finding | Options | Ruling |
|---|---|---|---|
| DR8 | **The gates do not fit CI's Release job.** DC10 puts the DC3–DC5 gates under `[.derived_goals]`; their Release times are g_m 113 s, swing 687 s, V_T 1429 s. CI's Release job takes 47–54 min of its 60-minute timeout (ctest `-j 4`). | (a) **Only the g_m gate in Release ctest** (`adapt_derived_gm_gate`, tag `[.derived_goals_ci]` beside `[.derived_goals]`); V_T and the swing run by hand with `[.derived_goals]`, their results recorded here. (b) All three in ctest, the timeout raised to 90 min. (c) None in ctest. | **(a)**. `CMakeLists.txt` registers `adapt_derived_gm_gate` and `adapt_derived_slow` (DC6) in the Release block. |
| DR9 | **Code review of 17c (2026-11-15, owner "Review Unit 17").** Five defects: (1) with a swing goal beside a current, vt or gm goal, every pass marked only by the swing's difference, discarding the others' adjoint weights (against DR7 (a)); (2) under `policy = intelligent` ID6's allocation replaced DR7's marks while the pass still reported `marked_by_difference`; (3) the verification's seconds left out M″ while its nodes were M″'s (the current goals were verified on M′); (4) the ID3 coarse step measured derived goals by the adjoint estimate DR4 replaced, so pass 0's miss test compared unlike estimates; (5) a figure not extractable on M′ fell back to M's value, giving a verified error never measured. Also: the swing's adjoint solves at its window were computed and never used (DR4 replaces the estimate, DR7 the marks). | (a) **Fix all six:** (1) the union of AD3's adjoint-weighted marks and the swing's difference marks; (2) no proposal with a swing goal (ID4's benefit runs through the adjoint, unreliable for the swing, DR4); (3) the whole verification's time, M′'s nodes; (4) no coarse estimate with a derived goal (pass 0 has no miss test); (5) an infinite verified error; the swing takes no part in the adjoint. (b) Fix (1)–(5) only. | **(a)**. Single-goal basic runs (the gates) are unchanged. Tests: a swing beside gm keeps every node gm marks alone; a swing under `intelligent` makes no proposal. Also fixed (owner, 2026-11-15: "yes fix 7 and 10"): M″ is solved only when M′ accepts the stop or on the last pass (whose verification is the certificate's); `solve::swept_contacts` is the one rule for `adapt()` and the deck (whose copy read past a shorter point). Not changed (recorded): verification re-solves M′'s window that DR4 already solved; the window ±1 range and the current overlay are written three to four times. |

**DC3 and DC5 gates (Release, 2026-11-15; `[.derived_goals]`, adjoint, tol 1e-2, V_T 1 mV, `max_passes` 20).** All pass.

| Goal | Certified | Mesh (passes) | Time | Actual error (allowed) | Estimate / actual (DC5: within 3×) |
|---|---|---|---|---|---|
| g_m | 0.385008 | 3036 nodes, 9108 DOFs (9) | 113 s | 1.11e-3 (3.85e-3) | 1.4 |
| Swing | 70.1284 mV/dec | 10 098 nodes, 30 294 DOFs (5) | 687 s | 0.313 (0.70) | 0.77 (verified 0.32) |
| V_T | 153.549 mV | 12 719 nodes, 38 157 DOFs (14) | 1429 s | 0.73 mV (1 mV) | 1.2 (verified 0.87 mV) |

**DC4 (Release, 2026-11-15; `[.adapt_derived_baseline]`, `goal = current drain`, adjoint, `max_passes` 20).**
- Baseline on the V_T/g_m window (points 13, 14, 15): tol 1e-2 certified at 3675 DOFs (7 passes) with V_T 149.90 mV (4.4 mV off: misses DC3's 1 mV) and g_m 0.38128 (1.3%). At 3e-3 it reached the cap: 409 332 DOFs (136 444 nodes, 20 passes, not certified) with V_T 154.18 mV and g_m 0.38800, both within DC3's allowances (6210 s). The 1e-3 and 3e-4 rows repeated it exactly (capped), so the run was stopped before the swing rows.
- Baseline on the swing window (points 3, 4), from the earlier run (`max_passes` 12): tol 1e-2 reached 19 968 DOFs (12 passes, not certified) with 70.1426 mV/dec, 0.30 off, within the 0.70 allowance.
- Against the baseline: V_T 38 157 / 409 332 = **0.093×** (≤ 0.8×: passes). Swing 30 294 / 19 968 = **1.52×: misses**.
- Against the uniform ladder (DR3; nodes for the same error): V_T 12 719 nodes at 0.73 mV, where 16× (61 985 nodes) is still about 1.8 mV off (about 151 k nodes by log-log): about 0.08×, passes. g_m (not required) 3036 nodes; no level up to 16× reaches 1.1e-3. Swing 10 098 nodes at 0.313, where 8× (15 633 nodes) is at 0.315: 0.65× by interpolation (misses 1/3); read level by level, 16× (61 985) is the first at or below it (0.16×).
- **DC4 misses for the swing on both halves.** The current goal on points 3, 4 already weights the subthreshold currents the swing is built from, and the swing's error is in y (DR6), which uniform y refinement reaches cheaply. **Accepted as a limitation (owner, 2026-11-15: "Accept as limitation (Recommended)")**; stated in ARCHITECTURE 11.3. The swing still certifies within its allowance (DC3, DC5).

**DC9 (hand mutants, 2026-11-15; `/c/Users/disha/mut17/mutants17c.txt`, Release, each mutant built and run against
`[adapt_derived],[adapt_learn],[adapt_policy]`, `deck adapt*` and `[.adapt_derived_slow]`, restored after each).** 53
mutants covering DG1–DG10, DG12 and DR2–DR9. **51 caught, 1 equivalent, 1 survivor (DG6, accepted by the owner).**
- Round 1 (the parallel runner): 31 caught, 16 survived, 3 not compiling. Several survivors were false: the runner's
  Catch2 shards ran in random order and missed tests (fixed with `--order decl`). Round 2 ran **serially** (owner,
  2026-11-15: "dont run parallel"): the 16 survivors and compilable variants of the stillborn mutants; 9 caught.
- Tests added for the round-2 survivors (60d4444 and this commit), all caught in round 3:
  - DG2 (indicators united over the window's points, not its point alone): DC6's pass-1 mark count (136; the point
    alone 134).
  - DG4 (adjoint weights summed over the window): DC8's pass-0 marks (3) and pass-1 nodes (312; the point alone 4 and
    364).
  - DG7 / DR5 (M′'s figure from the points solved on M′ alone): DR6's case, verified error 1.33474 (a curve mixing M and
    M′ jumps segment: 43.7).
  - DG9 (no kinds or vt_tolerance in the report of a current-only run): AC9's deck test.
  - DR2 (the deck's learned predictor asked for the run's goal kind): a one-tree model splitting on feature 12 leaves a
    gm run's start mesh unrefined.
  - DR6 (the larger of the single- and two-level errors) and DR9 (M″ only with a derived goal, when M′ accepts or on the
    last pass): `derived_verified_error` and `two_level_verification` are now tested functions in `adapt.hpp`, as
    Unit 17's slice 7b did for the loop's rules; `two_level_error` is declared there too. The DR6 test's case had
    changed since DR6 (now 275 → 480 nodes; M″'s moved window rejects the stop), and its comment is corrected.
  - DR9 (no coarse estimate with a derived goal): an intelligent run with a current and a gm goal, a prediction
    refining one far corner: pass 0 is not a miss.
- **Equivalent: DR9-swing-adjoint** (the swing's window points kept in the adjoint pass). The swing's estimate is
  pushed as infinite (DR4 replaces it) and its weights skip ω (`continue` before ω), and every other goal sums only
  over its own window, so the extra transposed solves change only time (and `estimate_fell_back` when one of those
  discarded solves falls back).
- **Survivor: DG6-moved-on-M1** (M′'s window-move check removed). Not equivalent: a stop could be accepted when M′'s
  window moves but M″'s returns to M's and every error is within its allowance. Not observed in six probed
  configurations (DR6's case, DC8, a V_T run at 1 mV and 0.1 mV, a g_m run without weighting, and DC6's). In those, a
  move on M′ either changes the segment (the single-level error then exceeds the allowance) or is caught again on M″.
  **Accepted as a survivor (owner, 2026-11-15: "Accept as survivor (Recommended)").**

**DC10 (suites, 2026-11-15; on 259c6b2).** Release, `/W4 /WX`, ctest serial (`-j1`, owner ruling): **222/222 pass**
in 2,864 s (log `/c/Users/disha/mut17/dc10_rel.log`); the build is warning-free (`/WX`). The first attempt was stopped
by the owner at about test 125 with all passing, and was rerun from the start. The local build has no MKL runtime
configured, so the PARDISO cases skip here (they are required in CI; PARDISO on oneMKL 2026.2.1 is PR #84).
**Debug is left to CI (owner, 2026-11-15: "skip debug locally, leave it to CI"; confirmed for 17c the same day).**
So DC10's Debug half, and its bound that Debug grows by ≤ 60 s, are judged on this branch's CI run, not locally.
Release was run at the branch head only, not at every commit; the earlier commits ran their own fast tests.

**DC10's Debug bound on CI (2026-11-16): missed, then the two cases moved to Release only (owner: "(a) Run those two
cases in Release only").** On `ca04aa5` the push run's Debug job passed in 2h14m, but `deck.adapt` took 6,231 s
(50 s on `main`); the PR run's Debug job was cancelled at its 150-minute limit with 199/203 passed and `deck.adapt`
still running after 119 min (its runner was slower: `adapt.adapt_derived` 2,003 s against 1,262 s). In Release
`deck.adapt` went from 2.6 s (`main`) to 339 s. The growth is the two deck cases this unit added, "derived goals in
the deck, the CLI, the table and the report" and "a learned predictor sees the run's goal kind (DG10, DR2)": each runs
adaptive loops to completion. They now carry the hidden tag `[.deck_adapt_derived]` and run as the Release-only ctest
entry `deck.adapt_derived_release` (as Unit 55's slow cases). Locally, Release: `deck.adapt` 1.6 s,
`deck.adapt_derived_release` 185 s, both pass. Debug no longer runs the two cases, so it should return to about
`main`'s time plus `adapt.adapt_derived`; to be confirmed on CI.

---

