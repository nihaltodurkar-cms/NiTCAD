<!-- Archived from DECISIONS.md on 2026-11-20 (branch docs/loose-ends), word for word; DECISIONS.md keeps its summary. -->

### Unit 59 — a natural-test norm for reverse bias (proposed 2026-11-18; approved 2026-11-18; started 2026-11-19; on `solve/natural-norm`)

**Owner request (2026-11-18):** "write a proposal for item 1": the 11.3 row "Robust nonlinear solving", next step "A
natural-test norm monotone in reverse bias". The branch starts from `assemble/generic-equations` `3b4024a` (PR #87):
the norm needs the unknowns' layout, which Unit 57 made explicit (Unit 40's FE motivation; Unit 57 named this as a
consumer). It is rebased onto `main` once #87 merges.

**What is measured today (the starting point):**
- **Unit 30 (L1–L6):** `LineSearchTest::natural` passes a trial when the simplified correction J(x)⁻¹F(trial) is at
  most (1 − λ/4) of the full correction, both in `update_size`: the largest of |dψ| (V_T) and |dn|/n, |dp|/p (each
  density against a floor, 1e-20 or 1e-8 of the largest density; Unit 50's opt-in minority floor raises it).
- **Its recorded failure:** on the legacy diode from equilibrium straight to −5 V, plain Newton converges in 45
  iterations and the natural search fails. Along plain Newton's own path the measure jumps between 1e2 and 1e8,
  because it is the relative correction of minority densities near 1e-14. Two variants were measured and dropped:
  scaling the test by the clipped fraction (MOSFET 1/10, −5 V still failing) and the largest scaled component
  (1.0 V then failing too).
- **Unit 50's finding:** on today's defaults the Unit 30 MOSFET benchmark (field mobility, V_G = 1 V, drain
  0.1–1 V) gives natural 2/10 (10/10 recorded) and natural with step cutting 10/10 in 252 iterations. Not bisected.
- **Two details of today's test, read from `newton.hpp`:** the full correction is measured at the iterate x and
  the simplified one at the trial, so the densities' scales differ between the two sides of the inequality
  (Deuflhard's test holds the scaling fixed within a step); and the full correction is unclipped while λ is the
  fraction of the clipped step.

**Standing conditions:** opt-in. The default path (no line search) and the residual test are untouched and
bit-identical (Unit 57's GG1 method). Unit 30's natural test stays available unchanged unless NQ3 rules otherwise.
`solve` gains no dependency on `assemble`: the norm is a system hook found by the concept, like `update_size`.
Criteria are fixed before code. Release and Debug suites pass at every commit with `/W4 /WX`, plus hand mutants.
Built with MSVC 14.51. Nothing is pushed or merged without authorization.

**Left out of this unit:** a Gummel or pseudo-transient fallback (the other half of the 11.3 next step; NQ1); a line
search inside `trace_bias`'s corrector; making any line search the default; changing the convergence test
(`update_size`, which every recorded result depends on).

**Design:**

| # | Slice | Change | Notes |
|---|---|---|---|
| ND1 | 0 | **Measure first.** A `[.benchmark]` diagnostic that replays plain Newton's converging path (diode at −5 V; the MOSFET's first drain point) and, at every iteration, records θ = ‖Δx̄‖ / ‖Δx‖ in today's measure and in each candidate norm (ND3), with both corrections scaled at the iterate. A norm suits the test when θ < 1 on plain Newton's accepted steps. Also the baseline of NG3/NG4 on today's tree, and (NQ5) a bisection of Unit 50's MOSFET finding. | No product code. Decides ND3 by measurement; the result is recorded before slice 1. |
| ND2 | 1 | **The hook.** `newton.hpp`: if the system has `natural_size(x_ref, dx)`, the natural test measures both corrections with it at the same x_ref (the iterate); otherwise it uses `update_size` as now. Selected by `LineSearch::norm` (NQ3). | A system without the hook, and `norm = convergence`, give Unit 30's test exactly. |
| ND3 | 1 | **The norm, from the layout.** `System::natural_size`, built from the DofMap's fields: each component scaled by w = max(\|x_ref\|, a field threshold), combined per field and then over fields. Candidates measured in ND1: **(a)** Deuflhard's weighted RMS (NLEQ's scaled norm), potential threshold 1 V_T, density threshold τ times the field's largest value; **(b)** potential-only (the densities follow ψ); **(c)** log-compressed densities, asinh(dn/w) in V_T units next to dψ. Held and insulator nodes are skipped. Fields added by an extension (GG8) are measured like ψ. | The choice and τ come from ND1 (NQ2). The field classes are kept inside `SystemData`; `Field` is unchanged. |
| ND4 | 1 | **Configuration and deck.** `[solver] natural_norm = convergence \| weighted` (NQ7), `BiasOptions` and the run record (only when not the default, as Unit 30's settings). | Old decks parse and run unchanged. |
| ND5 | 2 | **Results.** NG3–NG6 on the chosen norm; ARCHITECTURE 11.3 "As built" after acceptance. | |

**Criteria (frozen once approved):**

| # | Gate |
|---|---|
| NG1 | Bit-identity: every traced suite entry equal to the base tree's with the line search off, and with the natural test on `norm = convergence` (Unit 57's GG1 method). |
| NG2 | Decks: every deck's output unchanged (GG2's method). |
| NG3 | Diode jumps (`[.diode_jumps]`), natural search with the new norm, from equilibrium: converges at 0.6, 0.8, 1.0 V **and −5 V** (every point plain Newton converges); each state equal to a 0.05 V-step sweep's to 1e-12. 1.2 V and −20 V reported, not gated. |
| NG4 | MOSFET line-search benchmark: natural search with the new norm 10/10, each point equal to the 1 V_T sweep's to 1e-12; not fewer points than the old norm on the same tree. |
| NG5 | Unit tests: the norm's value on a hand-built 3-node system (held nodes skipped, an extension field like ψ); the natural test's λ sequence unchanged when F's rows are rescaled (affine invariance); the scalar `newton_test` cases pass under both norms. |
| NG6 | No regression where the natural test passes today: Unit 30's diode 0–0.6 V sweep and every `newton_test` case still converge, with iteration counts reported. |
| NG7 | Release and Debug suites, `/W4 /WX`; hand mutants of the hook, the scaling point and the thresholds, each caught. |

**Questions for the owner (recommended option first):**

| # | Question | Options |
|---|---|---|
| NQ1 | Scope | **The norm only (recommended)**; Gummel or pseudo-transient fallback as its own later unit. / Both in this unit. |
| NQ2 | How the norm is chosen | **Measure the three candidates in slice 0 and take the one that passes NG3 and NG4, ties to fewer iterations (recommended).** / Fix (a) now. |
| NQ3 | Old or new norm by default for `natural` | **A new option, `natural_norm`, default `convergence` (Unit 30 unchanged); promotion decided after the results (recommended).** / Replace Unit 30's norm (changes opted-in runs). |
| NQ4 | Base branch | **Stack on PR #87 and rebase after it merges (recommended).** / Wait for #87's merge. |
| NQ5 | Unit 50's MOSFET finding (natural 2/10 today, 10/10 recorded) | **Bisect in slice 0, at most about two hours, as NG4's baseline (recommended).** / Record today's numbers only. |
| NQ6 | The line search's default | **Stays off (recommended).** / Turn on with the new norm if NG3–NG6 pass. |
| NQ7 | Deck key | **`[solver] natural_norm = convergence \| weighted` (recommended).** / Another name. |
**Owner ruling (2026-11-18): "choose recommended and dont start unit 59".** NQ1–NQ7 take the recommended option:the norm only; the norm chosen by slice 0's measurement (NG3 and NG4, ties to fewer iterations); a new option`natural_norm`, default `convergence`; stacked on PR #87 and rebased after it merges; Unit 50's MOSFET findingbisected in slice 0 (at most about two hours); the line search stays off by default; deck key `[solver]natural_norm = convergence | weighted`. The criteria NG1–NG7 are frozen. **Not started:** no slice begins until theowner says so.
**Started 2026-11-19 (owner: "now start Unit 59"),** rebased onto `main` `e57918a`.

**Slice 0 results (ND1, 2026-11-19; Release, MSVC 14.51, Eigen SparseLU).** Test code only:
`tests/solve/natural_norm_probe.hpp` and two hidden `[.benchmark][.natural_norm]` cases (`step_control_test.cpp`,
`mosfet_test.cpp`). The probe replays `newton_solve` through `System`'s public hooks; its plain path and its
natural search under today's measure equal `newton_solve`'s iteration for iteration (update and step, checked in
both cases: 1527 + 515 assertions pass). Candidates use the iterate for both corrections, combine per field (RMS)
then over fields (RMS), skip held unknowns (Jacobian row = own diagonal), ψ weighted by max(|ψ|, 1).
- **Natural search per measure (NG3's points; iterations, F = failed in 100):**

| Measure | 0.6 V | 0.8 V | 1.0 V | −5 V | 1.2 V (rep.) | −20 V (rep.) | MOSFET (NG4) |
|---|---|---|---|---|---|---|---|
| today (Unit 30) | 13 | 16 | 20 | F | F | F | 2/10, 185 it |
| (a) RMS τ 1e-4 | 15 | F | F | F | F | F | 4/10 |
| (a) RMS τ 1e-8 | 91 | F | F | 75 | F | F | 10/10, 93 it |
| (a) RMS τ 1e-12 | F | F | F | F | F | F | 10/10, 93 it |
| (b) ψ only | 14 | 16 | 20 | 45 | F | F | 0/10 |
| (c) asinh τ 1e-4 | 14 | 20 | 42 | 100 | 43 | F | 4/10 |
| **(c) asinh τ 1e-8** | **37** | **44** | **50** | **61** | **53** | F | **10/10, 93 it** |
| (c) asinh τ 1e-12 | 25 | 34 | 41 | F | 45 | F | 10/10, 93 it |

  Only (c) asinh at τ = 1e-8 passes NG3 and NG4, so by NQ2 it is the norm for slice 1 (no tie). Every converged
  diode state is within 3.8e-15 (relative) of the 0.05 V sweep's.
- **θ along plain Newton's path (λ = 1, test passes at θ ≤ 3/4):** today's measure jumps (diode −5 V largest θ 25.1
  at the trial, 2.5 at the iterate; MOSFET at 5 V_T up to 7.7e3). Measuring both sides at the iterate alone does not
  fix it (−5 V: 2 of 44 steps pass). The candidates' θ stays near 1 through plain Newton's slow middle phase
  (diode 0.8–1.0 V: θ 0.99–1.07 over many steps), which is why the search shortens steps there; asinh keeps θ within
  1.6 at every gated point (largest 1.53, at 1.0 V). ψ alone contracts on the diode (−5 V: 44/44 steps θ < 1) but not on the MOSFET.
- **MOSFET plain path (first drain point):** at 5 V_T fails in 100 (as recorded); at 1 V_T converges in 11.
- **NG4's comparison:** the largest plain relative difference against the 1 V_T sweep is 0.14 (today's norm) to 0.6
  (the candidates), all in densities below 1e-8 of their field's largest, where both runs (tol_update 1e-8) are at
  rounding. Weighted as the norm weighs them (ψ against max(|ψ|, 1), a density against max(|n|, 1e-8 of the
  largest)) every completed point is within 2.0e-15. Which comparison NG4's "to 1e-12" means is put to the owner.
  **Owner ruling (2026-11-19, "choose recommended and start slice 1"):** NG4's "to 1e-12" is that weighted
  comparison.
- **NQ5 bisection (Unit 50's finding):** natural test on the `[.line_search]` benchmark: `f3e2811` (Unit 30's merge)
  10/10 in 227 iterations, as recorded; `e14a53f` (#42) 10/10, 227; `2d1fdc7` (#44, Unit 35: Couderc band set and
  automatic Masetti holes as the silicon defaults) 2/10 in 185, the same as today. The change is the device's
  physics, not the solver: NG4's baseline on today's tree is 2/10, 185 iterations.

**Slice 1 (ND2–ND4, 2026-11-19).**
- `newton.hpp`: `LineSearch::norm` (`NaturalNorm::convergence`, the default, or `weighted`, validated). With
  `weighted`, the natural test and a system with `natural_size(x_ref, dx)`, both corrections are measured with it at
  the iterate x0 (the full one once per iteration, each trial's simplified one against it); otherwise Unit 30's test
  is unchanged, code path included.
- `System::natural_size` (ND3, candidate (c) at τ = 1e-8): held unknowns as `stamp_contacts` holds them (an ohmic
  contact's ψ, n, p; a Schottky contact's ψ; an electrode's ψ; an insulator's n, p; an isothermal sink's rise; the
  carrier temperatures on contacts and insulators), extension fields measured as ψ, in equilibrium ψ alone.
- `[solver] natural_norm = convergence | weighted` (deck, `deck/build.hpp`, the app's section table and vocabulary);
  the run record names it only when not `convergence`.
- The product reproduces slice 0's measurement of (c) at 1e-8 exactly: diode 37, 44, 50, 61 iterations (1.2 V: 53;
  −20 V fails), each state within 1.7e-15 of the 0.05 V sweep's; MOSFET 10/10 in 93 iterations, currents within
  5.2e-14 of the 1 V_T sweep's.
- NG5 tests: `natural_size` by hand on a four-node bar with the GG8 extension (two free nodes, so the 1e-8 threshold
  is exercised), the MOS's electrode and insulator densities skipped; Newton with the hook (both corrections at the
  iterate, never between trials), affine invariance under power-of-two row scales (steps and updates bit-equal), the
  scalar cases equal under either norm (no hook: Unit 30's test). NG6: the diode's 0–0.6 V sweep under the weighted
  test reaches the plain sweep's states (fields to 1e-9, currents within their resolution) in 39 natural-test
  iterations (31 with the convergence norm); the run identity differs.
- NG7 hand mutants, each caught: the hook ignored; the full correction measured at the damped update; the trial's
  at the trial; τ = 1e-4; ψ's threshold dropped; held unknowns measured; no asinh; an insulator's ψ held.
- Targeted Release runs only (solve `natural norm*`, `line search*`, `newton*`; assemble `[extensibility]`; deck;
  app `app deck*`): all pass. NG1, NG2 and the full Release and Debug suites wait for the owner's word (slice 2).

**Slice 2 (2026-11-19, owner: "start slice 2, choose recommended", then "after NG1 passes just push").** Base:
`main` `e57918a`, built in a separate worktree; both trees Release, MSVC 14.51.
- **NG2:** the 12 decks of `tests/decks/` run with both `nitcad.exe`: all 38 files agree (`[.deck_compare]`).
- **NG1:** the 89 ctest entries of the five trace-linked executables, `NITCAD_TRACE_FILE` set, `--order decl`, pinned
  to the P cores (mask 0x3F): 86 traces byte-identical to `main`'s; 3 identical in every test case `main` has, with
  only this unit's added cases (`assemble.extensibility`: the two NG5 cases; `solve.newton`: the three
  `natural norm` cases; `solve.step_control`: the NG6 case). No test case differs.
- **Not run (owner: push after NG1):** the full Release ctest suite, the Debug build and suite. NG7's suites, ND5's
  ARCHITECTURE 11.3 "As built" and the unit's close remain.

**Close (2026-11-20, owner: "do the loose ends from merged units").** Merged as PR #101, `main` `38163aa`
(2026-11-19, CI still running at the merge).
- **NG7 suites, met after the merge without a new run:** CI on `38163aa` (run 44554600656) passed Release 248/248
  and Debug 227/227 with `/W4 /WX`. Both local full suites later ran on trees containing this unit (Unit 58 on
  `main` `3c36224`/`737cfc3`): Release 250/251 (only `adapt_intelligent` AC11's timing check under `-j 6` load,
  which passes alone) and Debug 230/230.
- **ND5:** ARCHITECTURE 11.3 "As built (Unit 59)" written, after Unit 50's section; the "Robust nonlinear
  solving" row names the norm.
- The unit is closed; this record moved to `docs/decisions/unit-59.md`.

