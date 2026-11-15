<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 43 — transient, small-signal and trace analyses in the deck (proposed, authorized and accepted 2026-10-17; on `deck/analyses`, PR #53)

**Owner request (2026-10-17):** the scope of the next unit named by Unit 42's PD1. Documentation only until the
owner approves.

**Why now (audited on `main` at `145f23c`).**
- The deck (Unit 31, `deck/build.hpp`) has two analyses, `[equilibrium]` and `[sweep]` (`deck::Analysis`). The
  transient run (Unit 21), the small-signal run (Unit 22) and the arc-length trace (Unit 19) are reachable from C++
  only.
- Since Unit 42 every result has a file (`save_transient`, `save_small_signal`; a trace is a `Sweep`), so nothing
  but the deck stands between a user and these runs.
- ARCHITECTURE.md 11.3 names it as the input/output row's remaining step.

**Design (proposed):**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| AN1 | **Three new analysis sections**, `[transient]`, `[small_signal]`, `[trace]`. A deck has exactly one analysis section of the five (DK3 extended). `deck::Analysis` gains the three; `Plan` gains the options of each (`solve::TransientOptions` with its waveforms, `SmallSignalOptions`, `TraceOptions`). Every key maps to one C++ option, with the C++ default when absent. | One run per deck, as now; the deck stays a thin, strict layer over the solver options (DK2, DK3). | Five-way exclusivity to check and test. |
| AN2 | **`[transient]`:** `contact <name> = <waveform>` per contact (a contact not named is held at 0 V); waveforms `constant(v)`, `step(v0, v1, t)`, `ramp(v0, v1, t0, t1)`, `pulse(base, pulse, t_start, width[, rise, fall])`, `pwl(t v, t v, ...)`, `sine(offset, amplitude, f[, delay, phase])`, each one `solve::Waveform` factory. Keys `t_end` (required), `dt_initial`, `dt_min`, `dt_max`, `integrator = bdf2 \| backward_euler`, `adaptive = on \| off`, `rtol`, `density_ref`, `output_times = <series>`, `fields_every_step = on \| off`; `thermal <name> = <K>` (electrothermal: one temperature per thermal contact, as the run takes). | The `Waveform` factories already validate; the deck only parses their arguments. | A small function-call syntax (AQ1). |
| AN3 | **`[small_signal]`:** `contact <name> = <series>` operating points and `thermal <name> = <series>`, exactly as `[sweep]` (DK6); `frequencies = <series>` (required), `field_frequencies = <series>`. `logspace(a, b, n)` joins `linspace` in series (AQ5): frequencies are logarithmic. | Reuses the sweep's point rules unchanged. | — |
| AN4 | **`[trace]`:** `contact = <name>` (the swept contact) and `end = <V>` (required); `start <name> = <V>` per contact (default 0); `step`, `max_step`, `min_step`, `current_limit`, `max_points`, `density_floor`; `thermal <name> = <K>`. | `TraceOptions` field for field. | — |
| AN5 | **Outputs.** `[output]` gains `transient = <file>` and `small_signal = <file>`; a trace writes `sweep =` (it is a `Sweep`). Each analysis accepts only the outputs it can write; another is an error naming the analysis. `table =` writes a CSV per analysis (AQ3). `checkpoint =`: a trace's last point and a small-signal run's last operating point; refused with `[transient]` (AQ4). `device` and `mesh` as now. | The Unit 42 files are written by the runs that make them; a wrong output is a deck error, not a silent no-op. | — |
| AN6 | **`restart` stays `[sweep]`'s** (DK9). The new sections refuse it. | A transient or trace restart needs more than a converged point (history, arc direction). | No restart for long transients (later unit). |
| AN7 | **The runner:** `nitcad check` describes the new analyses (contacts' waveforms, frequencies, trace range); `nitcad run` prints progress per accepted time step, operating point or traced point, writes the outputs, prints the provenance; exit 2 when any run stopped early (its completed part written), as for a sweep. | One behaviour for every analysis. | — |
| AN8 | **Out of scope:** restart of the new analyses; a transient or small-signal run starting from a checkpoint file (`initial`); harmonic balance; new solver behaviour, models or defaults; `ARCHITECTURE.md` beyond 6.5 and 11.3 "as built" after acceptance. | Keeps the unit to the deck and the runner. | — |

**Owner decisions needed:**

| # | Question | Recommendation |
|---|---|---|
| AQ1 | Waveform syntax: function calls (`pulse(0, 1, 1e-9, 5e-9)`) as `linspace` already is, or words (`pulse 0 1 1e-9 5e-9`) as doping profiles are? | **Function calls**: named shapes with positional arguments read like `linspace`, and `pwl` needs pairs. |
| AQ2 | Time suffixes (`ps`, `ns`, `us`, `ms`)? DK5 allows suffixes on lengths only. | **Keep DK5**: times in seconds (`1e-9`), as every other non-length quantity. A later unit can add unit suffixes everywhere at once. |
| AQ3 | CSV tables: transient, one row per accepted step (time, step, order, the bias, total, conduction and displacement current per contact); small-signal, one row per (point, frequency) (the bias per contact, then G_ij and C_ij for every contact pair); trace, the sweep table. | **As stated.** |
| AQ4 | `checkpoint =` with `[transient]`: a transient state at t_end is not a steady state a sweep can restart from. | **Refuse it** (an error naming the analysis). |
| AQ5 | `logspace(a, b, n)` in series: `n` values from a to b evenly spaced in log10, the first a and the last b exactly, both positive; available in every series (axes included, harmlessly). | **Add it.** |

**Acceptance gates, fixed before any implementation:**

| # | Gate |
|---|---|
| AG1 | **Existing decks unchanged.** Every existing deck test passes unchanged; for `[equilibrium]` and `[sweep]` decks the plan, the run identity and the written files' values are as on `main`; the assembly traces of the deck and validation entries are byte-identical to `main`'s (`--order decl`, PG1's method). |
| AG2 | **The deck equals the API.** Each new analysis run from a deck gives, bit for bit, the result of the C++ call with the options the deck states (run identity included): a 1D diode transient with each waveform kind, a MOS-C small-signal run over a `logspace`, a diode trace; with and without the electrothermal model where the run supports it. |
| AG3 | **Every key and every error.** Each key reaches its option (a test per key, changing it from its default); absent keys give the C++ defaults. Each error is `invalid_input` "`<file>:<line>: <what>`" (DK4): unknown key, missing required key, malformed or invalid waveform, unknown contact, two analysis sections, `restart` outside `[sweep]`, an output the analysis cannot write, `checkpoint` with `[transient]`, `logspace` with a non-positive bound. |
| AG4 | **Outputs.** The written files load (`load_transient`, `load_small_signal`, `load_sweep`, `load_checkpoint`) equal to the in-memory results, provenance included; the CSV tables have AQ3's columns, one row per step, (point, frequency) or point, and are written atomically as the sweep table. |
| AG5 | **The runner.** `check` and `run` for each analysis; exit 0, 1 (a deck error) and 2 (a run stopped early: its completed part written) for each; a cancelled transient writes its accepted steps. |
| AG6 | **Suites.** Full Release and Debug suites at every commit; no warnings under `/W4 /WX`. |
| AG7 | **Mutation checks.** Hand mutants in the new parsing, plan building, output selection and tables (a key ignored, a default changed, a waveform argument swapped, an exclusivity check dropped, a table column dropped) are caught; run in Release (Debug mutants open modal assertion dialogs). |
| AG8 | **No scope creep.** Nothing outside AN1–AN7; `ARCHITECTURE.md` 6.5 and 11.3 updated to "as built" only after the owner accepts. |

**Commits (proposed):** 1) the three sections, waveforms and `logspace` in `deck::build`, with the plan tests (AG3);
2) running them and their outputs and tables in `deck::run` (AG2, AG4); 3) the runner (AG5). Each passes AG6 on its
own.

**Risks:**
- `deck/build.cpp` is already large (over 1000 lines); the new sections add a few hundred. A split into a separate
  `analysis.cpp` may be worth it in commit 1.
- A transient table with `fields_every_step` and many steps is large; the table holds terminal quantities only.
- AG2's bit-for-bit comparison needs both sides to see the same options, including defaults the deck fills; any
  default the deck sets differently from C++ would show here, as intended.

**Owner approval (2026-10-17):** "start the scope", taken as approval with every recommendation: AQ1
function-call waveforms; AQ2 times in seconds (DK5 kept); AQ3 the CSV columns as stated; AQ4 `checkpoint` refused
with `[transient]`; AQ5 `logspace`.

**As built, so far:**
- Commit 1: the three sections in `deck::build` (`Plan::transient`, `waveforms`, `small_signal`, `trace`; their
  `steady` left default, the run uses `Plan::options`); waveforms in `src/deck/waveforms.hpp`; `logspace` in the
  series reader. AN3 refined: series never had `linspace` (an axis form), so `logspace` is a series form of its own,
  usable in every series. `[sweep]` and `[small_signal]` share one operating-point reader (the sweep's, moved
  unchanged). Until commit 2, `deck::run` refuses the new analyses.
  - **The one change to existing behaviour (AG1 refined):** the message for a deck without exactly one analysis
    now lists the five sections; the two existing tests that quoted the old message were updated to it.
  - Tests (`tests/deck/analyses_test.cpp`, 6 cases): every key changed from its default reaches its option;
    absent keys equal the C++ defaults; each waveform form equals its factory; `logspace` values and exact bounds
    (also where 10^log10(x) does not round-trip); the outputs each analysis takes and refuses; every mistake
    with its line.
  - Mutants (Release): 14 applied, 14 caught (a key ignored: `rtol`, `max_points`; `t_end` not required; step
    arguments swapped; no analysis accepted; `logspace` bounds not exact, first and last, caught only after the
    inexact-bound test was added; the trace start always on the first contact; a transient checkpoint accepted;
    field frequencies stored as frequencies; a pulse with a rise only; `backward_euler` read as `bdf2`; the sweep
    output refused for a transient only; `[trace]` parsed as `[transient]`). A first form of the sweep-output mutant
    did not compile and was replaced.

- Commit 2 (`b77aacc`): `deck::run` calls `solve_transient`, `solve_small_signal` or `trace_bias` with the
  analysis's options and `Plan::options` as their `steady`; writes the transient, small_signal or sweep file, the
  checkpoint (small signal: the last operating point; trace: the last point) and the AQ3 table
  (`write_transient_table`, `write_small_signal_table`; a trace uses the sweep table). `RunReport` gains
  `transient` and `small_signal`; a trace's result is `sweep`.
  - Tests (`tests/deck/analyses_run_test.cpp`, 3 cases; AG2, AG4): each deck run equals the C++ call on the
    hand-built diode (not the deck's device), bit for bit, run identity included: a transient (a ramp, fixed
    backward-Euler steps, an output time), a small-signal run (two points, a `logspace`, field frequencies), a
    trace (with a solver option); the files reload equal, provenance included; the checkpoints; the tables'
    headers, row counts and a whole small-signal row rebuilt from the run.
  - AG2 narrowed: a diode stands in for the MOS-C the proposal named, and the waveform kinds are checked against
    their factories in the plan tests (commit 1), with one run (a ramp) through the solver; no electrothermal run
    from a deck (the thermal keys are checked in the plan).
- Commit 3 (`892d5fc`): `nitcad check` describes the new analyses; `run` prints progress per accepted step, per
  operating point and frequency, and per traced point; exit 2 when a trace, transient or small-signal run stops
  early. Tests: the CLI for each analysis (check, run, exit 0), a stopped transient and small-signal run (exit 2),
  and a stopped transient's file holding its accepted steps (AG5; the stop is a non-convergence, which takes the
  same path as a cancellation).
- Mutants for commits 2–3 (Release): 11 applied, 11 caught (displacement columns dropped; the small-signal
  checkpoint of the first point; each run without the plan's options, the trace's caught after its test gained a
  solver option; capacitances written as conductances, caught after the whole-row test was added; the trace table
  not written; the step column holding the time; a stopped transient exiting 0; waveform names swapped; trace
  progress printed as a sweep's).
- **AG1 result:** with `--order decl`, the same MKL configuration and P-core pinning, all 18 of `main`'s deck
  test-case blocks and all 8 validation blocks are byte-identical; the only other blocks are this unit's 11 new
  cases.
- **Suites (AG6):** commit 1 Release 20/20, Debug 14/14. Commit 2 changes the `deck` library and its tests only,
  so its gate was the deck suite (Debug, 27 cases), not the full suites, to save about 40 minutes; commit 3
  Release 20/20, Debug 14/14. No warnings under `/W4 /WX`.

**Status:** accepted by the owner (2026-10-17, "accept"); ARCHITECTURE.md 6.5 and 11.3 updated to "as built"
(AG8). Merged to `main` as `913dac0` (PR #53, 2026-10-17).
