<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 44 — an iterative linear solver (proposed and authorized 2026-10-17; on `linalg/iterative`)

**Owner request (2026-10-17):** the scope of the "Scalable linear solvers" row of ARCHITECTURE.md 11.3 (planned:
"iterative Krylov with ILU/AMG behind `linalg`; parallel assembly"). Documentation only until the owner approves.

**Why now (audited on `main` at `913dac0`).**
- Every solve is a direct sparse LU: Eigen SparseLU (default) or MKL PARDISO (Unit 18, opt-in). Fill and time grow
  steeply on 3D meshes: one SparseLU refactorization of the 20³-node 3D diode took 5.8 s on this machine (Unit 41's
  cost measurement), and 3D meshing has been on `main` since Unit 16.
- The interface is ready for it: every run creates its solver from `BiasOptions::linear` (`bias`, `equilibrium`,
  `transient`, `small_signal`, `trace`, `weighted_current`); the backend interface (`src/linalg/backend.hpp`) is
  analyze / factorize / solve; `SolveReport::iterations` already exists, documented "for iterative backends".
- The acceptance check is backend-neutral: every solve is accepted on its componentwise backward error
  (`max_backward_error`, 1e-8), with iterative refinement.
- Eigen 5.0.1 (vcpkg, pinned) has `IncompleteLUT` and `BiCGSTAB` (supported) and GMRES, IDR(s), BiCGSTAB(L) and
  ILU(0) in `unsupported/`. No new dependency is needed.

**Design (proposed):**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| LS1 | **A third backend, `SolverBackend::gmres_ilut`, opt-in.** "Factorize" builds the preconditioner (Eigen's `IncompleteLUT`, with its own fill-reducing ordering); "solve" runs restarted GMRES with right preconditioning from x = 0. Real and complex (`ComplexLinearSolver`, the small-signal systems), as the other backends. | One more backend behind the same interface: the runs, the deck and every acceptance check are unchanged; the default backend and its results are untouched. | Eigen's ILUT, not an in-tree one: its ordering and drop rule are Eigen's. |
| LS2 | **GMRES in-tree (LQ1)**, about 200 lines in `linalg`: restarted GMRES(m) with modified Gram–Schmidt and Givens rotations. Its stopping test is the layer's own acceptance measure: each restart cycle ends with the componentwise backward error of the true residual, and the solve stops when it is at most `max_backward_error` or the iteration limit is reached. | The stopping test must be the componentwise backward error the layer accepts on (GMRES's own residual estimate is a 2-norm of the preconditioned system); `unsupported/` code has no stability promise. | Code to own and test; Eigen's GMRES would be less. |
| LS3 | **Configuration** (`SolverConfig`): `restart` (30), `max_iterations` (1000), `ilut_drop_tolerance` (1e-4), `ilut_fill_factor` (10). The existing `equilibrate` applies before the preconditioner, as for the direct backends. `threads` must be 1. | Equilibration matters more for ILU than for LU on drift-diffusion's widely scaled rows. Single-threaded keeps results deterministic on one machine (6.8). | Defaults are guesses until LG4 measures them. |
| LS4 | **Failure:** a solve that does not reach `max_backward_error` within `max_iterations` is `inaccurate_solve`, as a direct solve that fails acceptance; `SolveReport::iterations` reports the iterations. No pivot ratio (ILUT's diagonal is no singularity measure, as PARDISO's), so `min_pivot_ratio` is not applied; a singular system shows as non-convergence. Fallback to a direct solve: LQ2. | One failure path for every backend: Newton, step cutting and the run's `stopped` already handle `inaccurate_solve`. | A floating region is not reported as `singular_system` by this backend. |
| LS5 | **The deck:** `[solver] backend = gmres_ilut`, with `restart`, `max_linear_iterations`, `ilut_drop`, `ilut_fill`. | The batch runner reaches it like PARDISO. | Four keys. |
| LS6 | **Out of scope:** algebraic multigrid; block (field-split) preconditioners (they need Unit 40's unknown map, a later unit); parallel assembly; a GPU solver; changing the default backend. | Keeps the unit to one backend and its evidence. | 3D scaling beyond what ILUT gives waits for AMG. |

**Owner decisions (asked; ruled below):**

| # | Question | Recommendation |
|---|---|---|
| LQ1 | GMRES in-tree (LS2), or Eigen's `unsupported/` GMRES with a componentwise check after it? | **In-tree.** The stopping test is the acceptance measure, and the code stays under this project's tests. |
| LQ2 | When GMRES fails, fall back to Eigen SparseLU for that system (and say so in the report), or fail with `inaccurate_solve`? | **Fail** in this unit (one behaviour, measurable failure rate); a fallback can follow once LG4 shows how often it fails. |
| LQ3 | Complex systems (small signal) in this unit? | **Yes**: the same template; the high-frequency admittance systems are a useful test of the preconditioner. |
| LQ4 | If LG4 shows no device where GMRES-ILUT beats SparseLU in time or memory, merge it anyway as an opt-in (documented as such), or stop? | **Merge as opt-in** with the measurements recorded; it is the base AMG and block preconditioners plug into. |

**Acceptance gates, fixed before any implementation:**

| # | Gate |
|---|---|
| LG1 | **The default is untouched.** With `eigen_sparse_lu` every suite passes unchanged and the assembly traces are byte-identical to `main`'s (`--order decl`, PG1's method). |
| LG2 | **GMRES on its own.** Linear-algebra tests: GMRES with no preconditioner and with ILUT reaches the backward-error target on nonsymmetric real and complex systems (including `tests/linalg`'s shared systems); its iteration count on a matrix with k distinct eigenvalues is at most k; restart and iteration limits are honoured; a singular or unconverged system is `inaccurate_solve` with the iterations reported; `threads != 1` and bad settings are `invalid_input`. |
| LG3 | **It agrees with SparseLU** on the devices, as PARDISO does (Unit 18's `[pardiso]` cases as `[gmres]` cases): the diode bias gates, the MOS quasi-static C-V, the small-signal admittances, the transient switch-on, within the same acceptance tolerances; the Newton iteration counts reported. |
| LG4 | **Measured, not assumed:** a Release benchmark (`[.benchmark]`) of one Newton step's linear solve on the 2D MOSFET and on 3D diodes of growing size (at least 20³, 30³ and 40³ nodes where SparseLU still fits): time, peak working set, GMRES iterations, for SparseLU and GMRES-ILUT (and PARDISO where configured). The results go into DECISIONS.md whatever they show (LQ4). |
| LG5 | **The deck** keys of LS5 reach `SolverConfig`; a 3D deck runs with `backend = gmres_ilut`. |
| LG6 | **Suites:** Release and Debug at every commit; no warnings under `/W4 /WX`. |
| LG7 | **Mutation checks** (Release): GMRES (a skipped reorthogonalization, a wrong Givens sign, the stopping test on the estimate instead of the true residual, a restart that keeps the old basis), the preconditioner applied on the wrong side, a configuration check removed, a deck key ignored. |
| LG8 | **No scope creep:** nothing outside LS1–LS5; ARCHITECTURE.md 6.10 and 11.3 updated to "as built" after acceptance. |

**Commits (proposed):** 1) GMRES and the backend in `linalg`, with LG2; 2) the device agreement cases (LG3) and the
deck keys (LG5); 3) the benchmark and its results (LG4).

**Risks:**
- Drift-diffusion Jacobians are hard for ILU: the continuity rows are strongly nonsymmetric at high fields, and the
  scales of ψ, n and p differ by tens of decades before equilibration. GMRES-ILUT may fail on some devices; LG3 and
  LG4 will show where, and LQ2 decides what a failure does.
- A backward error of 1e-8 componentwise may take many iterations where a direct solve gets 1e-16 at once; the
  benchmark may show a looser `max_backward_error` is needed for speed, which is the user's choice per run, not a
  default change.
- Results differ from SparseLU's at the tolerance, not bit for bit, as PARDISO's do.

**Owner rulings (2026-10-17).** The owner approved the unit and delegated LQ1–LQ4 to the choice best for a
commercial TCAD product ("for decisions choose whats best for Commercial TCAD").

| # | Ruling | Why, for a commercial product |
|---|---|---|
| LQ1 | GMRES in-tree (LS2). | The stopping test is the acceptance measure, and the solver stays under this project's tests and stability promise. |
| LQ2 | **Fall back** (differs from the recommendation). When GMRES-ILUT fails to reach `max_backward_error` within `max_iterations`, or the ILUT build fails, that system is solved by Eigen SparseLU, and the fallback is reported: `SolveReport` gains a field saying the system fell back, with the GMRES iterations spent. The direct solve's own checks then apply as usual (`inaccurate_solve`, `singular_system`). The deck key `linear_fallback = sparse_lu` (default) or `none` (fail with `inaccurate_solve`, the strict behaviour) controls it. | A production run must not abort where a direct solve can finish; the user sees in the report how often it happened. The strict mode keeps the failure rate measurable. |
| LQ3 | Complex systems (small signal) in this unit. | The same template; high-frequency admittances test the preconditioner. |
| LQ4 | Merge as opt-in with the measurements recorded, whatever LG4 shows. | It is the base for AMG and block preconditioners, the path to large 3D. |

**Amendments from LQ2.** LS4: failure falls back as LQ2, with `none` giving the behaviour written there. LS5: a
fifth key, `linear_fallback`. LG2: with `linear_fallback = none` an unconverged system is `inaccurate_solve`; with
`sparse_lu` the same system is solved, the fallback is reported, and its solution meets the backward-error target.
LG3: the agreement cases run with `none`, so a fallback cannot hide a failure; any device that fails is recorded.
LG4: also records each case's fallback count with the default setting. LG7: a further mutant, a fallback that is
not reported.

**Status:** accepted by the owner 2026-10-21 ("go next"); LQ3 amended by the owner (LQ3a); commits 1–3 done;
ARCHITECTURE.md 6.10 and 11.3 updated to "as built" (LG8). Merged to `main` as `b518c99` (PR #54, 2026-10-21).

**Commit 3 (LG4 benchmark), measured 2026-10-21** on the development machine (Release, one thread, 31.6 GB),
`[.benchmark][gmres_benchmark]` in `bias_test.cpp` and `mosfet_test.cpp`; MKL not configured, so no PARDISO
column. One Newton step's linear solve with a fresh solver (analysis, factorization or ILUT build, one solve), at
the default settings (`max_backward_error` 1e-8, restart 30, drop 1e-4, fill 10, fallback `sparse_lu`). "Held" is
the private memory the solver holds after the solve; "peak rise" how far the process's peak working set rose
during the step (GMRES-ILUT ran first and never raised it: its peak stayed below the Jacobian assembly's).

| Case | Unknowns | SparseLU | GMRES-ILUT | GMRES iterations | Fallbacks | Step difference |
|---|---|---|---|---|---|---|
| 2D MOSFET 90 × 50, Vg = 1 V | 28,470 | 0.23 s, held 122 MB | 0.47 s, held 48 MB | 120 | 0 | 3.3e-13 |
| 3D diode 20³, 0.5 V | 24,000 | 7.1 s, held 409 MB, peak +321 MB | 0.77 s, held 50 MB | 40 | 0 | 2.2e-14 |
| 3D diode 30³ | 81,000 | 105 s, held 3.1 GB, peak +2.6 GB | 3.8 s, held 174 MB | 60 | 0 | 7.4e-14 |
| 3D diode 40³ | 192,000 | 1021 s, held 10.9 GB, peak +8.0 GB | 14.1 s, held 422 MB | 90 | 0 | 8.7e-14 |

The MOSFET's Id-Vg sweep (16 points, Vds = 0.05 V): SparseLU 21.1 s and 100 Newton iterations, GMRES-ILUT 29.9 s
and 84, every drain current within max(1e-8 |I|, resolution) of SparseLU's (worst ratio 0.000).
**Reading (LQ4):** on the 2D MOSFET SparseLU is twice as fast and GMRES-ILUT holds 2.5× less memory; on 3D
GMRES-ILUT is 9×, 27× and 72× faster and holds 8×, 18× and 26× less memory, its iterations growing about as the
grid's linear size (40, 60, 90). No fallback was taken in any case. The default backend stays SparseLU (LS6);
these figures are the evidence for recommending `backend = gmres_ilut` for 3D devices.

**LQ3a (owner ruling, 2026-10-21): small-signal systems are solved directly.** Measured in commit 2: with GMRES
the admittances agree with SparseLU only within each run's own reported resolution, and that resolution reached
22% of max |Y| (p+n junction) and 93% (MOS capacitor with an interface trap level), the same at
`max_backward_error` 1e-8 and 1e-9; at 1e-10 GMRES-ILUT did not converge the trap device's operating point in 1000
iterations. An admittance is a near-cancellation in the continuity rows: a componentwise backward error of 1e-9
still leaves an absolute residual comparable to the AC current, and only a rounding-level (direct) solve resolves
it. Asked with three options (direct AC solves; keep GMRES with coarse resolutions; stop the unit), the owner chose
**direct AC solves**: with `backend = gmres_ilut`, `solve_small_signal` takes its operating points with GMRES-ILUT
and solves the complex small-signal systems with Eigen SparseLU (`small_signal.hpp`). The complex GMRES path stays
in `linalg`, tested (LG2), for later preconditioners. Not recorded in the results: the run record's
`linear.backend` says gmres_ilut, and each frequency's `pivot_ratio` is the direct factorization's.

**Commit 2 (device agreement LG3, deck keys LG5), as built.** `[solver]` keys `backend = gmres_ilut`, `restart`,
`max_linear_iterations`, `ilut_drop`, `ilut_fill`, `linear_fallback = sparse_lu|none`, each checked with its line;
the run record digests the five settings only with `gmres_ilut`, so every other run keeps its identity. `[gmres]`
cases (strict mode, in the `solve` ctest) beside Unit 18's `[pardiso]` ones, measured in Release:
- diode: one ordering per Newton run, one ILUT per Jacobian (10 iterations as Eigen, 71 preconditioner
  applications); J(0.5 V) 1.2800773807e-2 as Eigen; the 0–0.7 V sweep within 5.7e-9 relative, 69 Newton
  iterations as Eigen's 69;
- MOS quasi-static C-V: gate charge within 7.9e-17 of the largest, φ_s within 1.1e-16 V;
- transient switch-on (101 BDF2 steps): currents within 4.4e-21 of the largest;
- small signal (LQ3a): junction and trap level within 0.10 and 0.23 of the PARDISO case's bound;
- deck: a 3D diode deck (21 × 6 × 6 nodes) with `backend = gmres_ilut` sweeps to 0.4 V within 4.4e-16 of the
  SparseLU run.
- **Suites:** Release 20/20 and Debug 14/14, no warnings. **Mutation checks (LG7, Release): 12 of 12 caught** —
  each of the six deck keys ignored, each of the five run-record entries dropped, and the small-signal systems
  left on GMRES (LQ3a).

**Commit 1 (GMRES and the backend, LG2), as built.** `SolverBackend::gmres_ilut` and `LinearFallback`;
`SolverConfig` gains `restart`, `max_iterations`, `ilut_drop_tolerance`, `ilut_fill_factor`, `fallback` (validated
for every backend); `SolveReport` gains `fell_back`. The kernel is `src/linalg/gmres.hpp` (no Eigen; MGS with one
full reorthogonalization pass, Givens rotations, right preconditioning, x = 0 start); the preconditioner is
`src/linalg/eigen_ilut.cpp` (Eigen 5.0.1 `IncompleteLUT`, version-asserted; a non-finite factor is a failed
build). GMRES runs on the equilibrated system the ILUT was built from and accepts only on the componentwise backward
error of the true iterate, measured once per cycle; a cycle ends early when the residual estimate falls below
1/2 · β · max_backward_error / (current backward error), so a wrong aim costs a restart, never a wrong answer.
- **Finding, amends LS4.** A singular system does not show as non-convergence. On a floating region with an
  inconsistent right-hand side GMRES reaches an x of norm 1.1e17 with backward error 1.0e-16, exactly as Eigen
  SparseLU does with `min_pivot_ratio = 0` (|x| 3.1e16): the componentwise backward error cannot see singularity.
  `gmres_ilut` therefore does not detect floating regions; the device layer's topology check (`device.cpp`)
  rejects them before any solve, so device runs are unaffected. A fallback still applies the direct pivot check.
- **Finding, LG2's k-eigenvalue test.** Termination in k steps holds in floating point only when the minimal
  polynomial is well conditioned: with 6 eigenvalues over four decades its terms reach 1e19 before cancelling
  (measured) and GMRES needs a 7th step. The test uses 6 eigenvalues in [1, 13] and gets at most 6.
- **Measured (Release, `tests/linalg`, strict mode):** nonsymmetric n = 2000 and the contacted 1D Laplacian in 1
  iteration; the 40 × 40 2D Laplacian in 9 (forward error 3.8e-9 against SparseLU); the badly scaled coupled
  system (rows over 40 decades) in 170 (9.7e-9); complex G + iωC over ω = 1e-6…1e3 within 3.7e-14 of SparseLU.
- **Suites:** Release 20/20 and Debug 14/14, no warnings; `nitcad_linalg_test` +18 test cases (`gmres_kernel_test`,
  `gmres_ilut_test`). **Mutation checks (LG7, Release): 19 of 20 caught** — skipped reorthogonalization, Givens
  sign, stopping on the estimate, a restart keeping the old basis, left instead of right preconditioning, each of
  the five configuration checks, the fallback unreported (two places), the fallback ignoring `none`, the operator
  or the measure not equilibrated, the ILUT finiteness check removed, `direct_` not reset on refactorization. The
  survivor, a Givens rotation without the conjugate of its second entry, is equivalent: that entry is the
  Arnoldi norm, real and positive.
