<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 50 — Newton at the rounding floor: the residual line search and the minority-density floor (started 2026-10-28; on `solve/rounding-floors`)

**Owner request:** "solver unit" (2026-10-28): the unit Unit 47 proposed after diagnosing two Newton stalls
(DECISIONS.md Unit 47 on `physics/quantum-confinement`, PR #66): (a) the residual line search at a rounding
residual (the FinFET with MLDA stalled at 0.35 V); (b) the minority-density floor (the 1e17 long channel with
Lombardi chatters at 1–4e-8).

**Reproduced on `main` (before code), one device for both.** Unit 46's 1e17 long-channel n-MOSFET
(`surface_mobility_test.cpp`, Vg 0–6.5 V in 0.5 V, Vd = 10 mV) with Lombardi, Release:
- Default options: all 14 points, 55 s, but at 0.5, 2.5, 3.5, 4.0, 4.5 and 5.5 V the corrections chatter at
  1.0–3.8e-8 for up to 14 iterations and finish only when one dips below 1e-8. The largest correction is always
  a minority hole in an n+ pocket: p = 2.1–2.8e-19 of the scaled majority (about 27 cm⁻³ against 1e20), its
  correction 1–7e-27 of the majority, so 1–4e-8 of p. The global floor (1e-20 of the largest density) is below p,
  so these count.
- `line_search = residual`: no point solves (100 iterations at Vg = 0). From iteration 4 the largest |F| is
  3e-11, rounding of the large rows, while one minority-hole node (p = 2.75e-20) still needs a factor of 8.7; the
  Armijo test on Σ F² sees only the large rows' noise, takes 2^-10-ish steps, and p creeps by 0.6% an iteration.
- Measured against each row's own scale s_i = Σ_j |J_ij x_j| (the size of the terms it sums), that hole row's
  residual is 1.5e14 times ε s_i: it is not rounding, only invisible in Σ F². Converged states sit at
  0.9–7.5 ε s_i; the chattering states at about 1e3 ε s_i on the minority rows (the solve's rounding, not the rows').
- Prototypes (environment-switched, not committed): Armijo on Σ (F_i / s_i)² alone stalls the equilibrium
  solve at 2^-10 steps (the clipped step is not a descent direction for it, Unit 30's finding); passing a trial
  that satisfies Armijo on **either** merit solves all 14 points. Measuring n and p against
  max(density, global floor, ε max(n, p) at the node) converges every point quadratically to 3e-11–1.4e-10, 47 s,
  drain currents equal to the default run's within 1e-11 (11 of 14 bit-identical).

**Owner rulings (2026-10-28, all as recommended):** (1) the residual test passes a trial on either merit, or
when every row is within 64 ε of its scale; (2) the floor ε × the node's majority density, on by default;
(3) in scope: the FinFET deck as a regression, a recheck of Unit 48's rounding-limited criteria (report only),
and a retest of PR #66's FinFET with MLDA and the residual test once this merges.

**Design.**
- RD1 (`solve::newton_solve`, residual test only): at each searched iteration, s_i = Σ_j |J_ij x_j| from J and x
  at the iterate (1 where it is 0). A trial passes when Σ F² or Σ (F_i / s_i)² falls by the Armijo factor
  (1 − 2 α λ), or every |F_i(trial)| ≤ 64 ε s_i. Every trial the old test passed still passes, so a run changes
  only where the old test rejected the full step. The best trial (none passing) is ranked on Σ F², as before
  (corrected as built: the first text said the scaled merit). The natural test and runs without a search are
  unchanged.
- RD2 (`DriftDiffusion::update_size`, **opt-in**, see the second ruling below): a density correction is measured
  against max(density, global floor, ε max(n, p) at that node), ε = 2^-52. A minority density below ε of the
  majority at its node does not change the node's charge in double, and its correction is resolved only to the
  solve's rounding there. Majority densities (≥ the minority) are measured as before. Insulator nodes and θ
  unchanged. Off, `update_size` is the old one.
- RD3 (as ruled the second time): `solve::BiasOptions::minority_floor` (off), set on the system by
  `DriftDiffusion::set_minority_floor` wherever a bias, trace, transient or small-signal solve builds one (and
  on the sweep's retry system without generation); the deck key `[solver] minority_floor = on | off`; recorded
  in the run's settings only when on, so the identities of runs without it are unchanged. The app's deck
  vocabulary offers the key.

**Found in the suites with the floor on by default (2026-10-28), and the second ruling.** `solve.heterojunction`
failed in Release and Debug: its 4H-SiC p-n diode (incomplete ionization, 0.5 V steps to 3 V) stopped at 3.0 V
after 100 iterations. On `main` that point already needs 69 iterations, swinging to 1e12, but there Newton
resolves SiC's minorities (about 1e-35 cm⁻³, far below ε of the majority) quadratically, so RD2's premise (such
a minority's correction is rounding) does not hold in general; with the floor its corrections no longer count
and the point wanders. The current at 1.5 V moved from 5.9e-13 to −1.19e-12 A/cm², below the test's
resolution. **Owner ruling (2026-10-28, as recommended): make the floor opt-in**, behind a solver option and a
`[solver]` deck key, off by default, so every existing run stays bit-identical, SiC included; the long channel
and similar MOSFET decks turn it on, and RC2 runs with it on. RD1 stays on wherever the residual test is used.
The criteria below are kept as fixed, with RC2's and RC6's conditions as this ruling sets them; RC9 is added.

**Criteria (fixed before code).**
- RC1: the long channel with Lombardi and `line_search = residual` solves all 14 points (0 before).
- RC2: the long channel with Lombardi, default options (with `minority_floor` on, per the second ruling): every
  point converges with no more than 2 iterations whose correction lies in [1e-8, 1e-4) (up to 13 before), and
  its drain currents equal the old run's within 1e-9 relative.
- RC3: `newton_test`: a two-row system with rows scaled 1e20 apart, at a state where the large row is at rounding
  and the small row needs a real correction, converges with the residual test (and fails without RD1); a
  trial at rounding in every row passes; the Unit 30 cases are unchanged.
- RC4: `update_size`: a minority density below ε of the majority at its node is measured against ε × majority;
  above it, and for the majority, against the density as before; insulator nodes unchanged.
- RC5: Release and Debug suites pass. Unit 30's hidden diode-jump and MOSFET line-search benchmarks: every
  outcome (converged or not) as recorded in Unit 30 for the natural test; the residual test's recorded failures
  may only become successes.
- RC6: the FinFET deck (residual test, PARDISO, 6 threads) solves 17/17 before and after, drain currents within
  1e-6 relative (after: the same deck, so `minority_floor` off, per the second ruling).
- RC7 (report only): Unit 48's `[!shouldfail]` SC2/SC8 and SC13 rerun and their values recorded.
- RC8: hand mutants: the scaled merit dropped (RC1 fails), the at-rounding pass dropped (an RC3 case fails), the
  floor removed (RC2, RC4 fail), the floor taken from the global largest density instead of the node's (RC4
  fails); added with the second ruling: the floor applied whatever the switch (RC9 fails).
- RC9 (added with the second ruling): with `minority_floor` off, the default, `update_size` measures a minority
  below ε of its node's majority against itself (as before); `solve.heterojunction` passes; the deck key reads
  `on` and `off` and defaults to off.
- After merge (follow-up, done 2026-10-28): PR #66 rebased on `main` 0c14841; its FinFET with MLDA and `line_search = residual` solves 17/17 in 11.4 min, no cuts (Unit 47, limitations).

**As built (2026-10-28, Release and Debug, this machine).** Every criterion met.
- RC1: the long channel with Lombardi, `line_search = residual`: 14/14 points with the floor off and with it on
  (0 before).
- RC2: with `minority_floor` on, every point has at most 2 iterations in [1e-8, 1e-4); drain currents within 1e-9 of
  the 14 recorded before the change. Hidden test `[.lombardi]`, Release (the three sweeps: 145 s).
- RC3: `newton_test`: the two-scale system (rows 1e14 apart, the large one at rounding) converges with the residual
  test in at most 5 iterations from a full first step; a trial at rounding in every row passes (2 iterations); the
  Unit 30 cases unchanged.
- RC4, RC9: `drift_diffusion_test`: off, the minority p = 1e-19 of n = 1 is measured against itself; on, against
  ε × 1, against ε × 1e-3 at a node whose majority is 1e-3, with holes as the majority likewise; above ε of the
  majority and for the majority, against itself. `build_test`: the key reads on and off, off by default.
- RC5: Release 112/112 and Debug 103/103, `solve.heterojunction` included. Unit 30's `[.diode_jumps]` and
  `[.line_search]` benchmarks: outcomes and iteration counts identical to `main` (`c288391`, built from the same
  tree with the changes stashed). The residual test's MOSFET 0/10 stays 0/10: that case is damping-limited, not a
  rounding floor. **Finding, not Unit 50's:** `main` no longer matches Unit 30's MOSFET record (natural test 2/10
  against 10/10; ARCHITECTURE.md 11.3 "As built (Unit 50)").
- RC6: FinFET deck (residual test, PARDISO, 6 threads, floor off): 17/17 before and after; currents within 9.5e-12
  relative (point 0 within 7e-11); Newton iterations 428 → 182; point 0's substeps 3 → 1; wall time 28 → 10.6 min.
- RC7 (report only): SC2/SC8 9.5e-15, 1.4e-14, 5.6e-14, 2.9e-8 A/cm²; SC13 2.2e-15 A/cm; SC12 1.25e-4; SC6
  1.00395; SC5 0.4276: all as recorded in Unit 48 (the default path is unchanged). These are rounding of the
  current itself, not a Newton floor; the minority floor does not bear on them.
- RC8: M1 (scaled merit dropped): the RC3 two-scale case and RC1 fail; M2 (at-rounding pass dropped): the RC3
  every-row case fails; M3 (floor removed): RC4 fails (and RC2, measured with the floor always on); M4 (floor from
  the largest density in the device): RC4 fails at the 1e-3 node; M5 (floor whatever the switch): RC4's off case
  and `solve.heterojunction` fail. All caught.
- Found on the way, fixed here: `app.live_check`'s stale-result case (XG1) failed once locally under `-j 6` and
  once on CI (`dropped() == 1` read 0). The test waited for its check function to return, but the worker posts
  the result after that, so `process()` could run before it. The test now waits on the checker's wake, called
  after the result is posted; a 50 ms delay injected before the post fails the old test every time and passes the
  new one; the app tests pass 20 times under `-j 6` (Release) and 5 (Debug). The checker is unchanged.
