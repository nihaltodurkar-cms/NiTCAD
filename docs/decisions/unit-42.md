<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 42 — results and provenance (proposed, authorized and accepted 2026-10-17; on `results/provenance`, PR #52)

**Owner request (2026-10-17):** a scope proposal for the "Results and provenance" row of ARCHITECTURE.md 11.3
(status "partial"; missing "build/version stamp in files; transient and small-signal results on disk"). Documentation
only until the owner approves.

**Why now (audited on `main` at `67e8fe5`).**
- `results::RunRecord` holds `input_identity` (a digest of the inputs) and `settings`. Nothing says which build
  produced a result: `project(NiTCAD VERSION 0.0.0)` is the only version, and no source revision, compiler, linear
  backend, MKL version or start time is recorded anywhere.
- After Unit 41 a result is reproducible on one machine whatever the core. To reproduce it, a reader needs the build
  and the platform it came from, and a file does not carry them.
- `io` saves meshes, sweeps, checkpoints and devices (content kinds 1–4). `results::Transient` (Unit 21) and
  `results::SmallSignal` (Unit 22) have no file: a transient or AC run is lost when the process ends.
- `trace_bias` (Unit 19) returns a `Sweep`, so `save_sweep` already stores it (no test says so).
- The deck (Unit 31) runs `equilibrium` and `sweep` only (`deck::Analysis`).

**Design (proposed):**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| PV1 | **A `results::Provenance` record in `RunRecord`** (plain data; `results` still depends on `base` only). Build: NiTCAD version, source revision (git commit and a `dirty` flag), build type, compiler and version, Eigen version. Run: linear backend in effect, the MKL version string when PARDISO is loaded, Eigen's fixed cache sizes (Unit 41), thread count, OS version, CPU model string, logical CPU count, start time (UTC, ISO 8601) and wall time. | Enough to say which build, backend and platform produced a result and when; the items Unit 41 showed can change last bits are all in it. | A few hundred bytes per run. |
| PV2 | **Provenance is not hashed.** `input_identity` keeps its current definition. | Run identities and checkpoints written before this unit stay valid; identity is "same inputs", provenance is "how and where". | Two runs of one input on different builds share an identity; provenance tells them apart. |
| PV3 | **Build facts from a generated header** (`base`), refreshed at build time by a CMake step that rewrites it only when the revision or dirty state changes; outside a git checkout the revision is `unknown`. Run facts gathered by `solve` when it makes the run record (`make_run_record`), the wall time when the run ends. | A configure-time stamp goes stale after the next commit; rewriting only on change avoids relinking every build. | One custom build step; it needs `git` on the build machine's path (not a runtime dependency). |
| PV4 | **In files as a keyed record inside the run record** (IO2): sweeps, checkpoints and the new kinds carry it; files without it load with an empty provenance. No format version change. | Old files keep loading; new readers see the stamp. | Files written before this unit have no provenance, by necessity. |
| PV5 | **Transient and small-signal files:** content kinds 5 and 6; `save_transient` / `load_transient`, `save_small_signal` / `load_small_signal`; every field bit for bit (complex values as two doubles by bit pattern), reusing the point, field, convergence and error records. A test that a `trace_bias` sweep round-trips through `save_sweep`. | Closes the "on disk" half of the 11.3 row and the R2 deferral for these results. | Snapshot-heavy transients make large files (Unit 21's snapshot policy already bounds them). |
| PV6 | **The runner:** `nitcad --version` prints the build stamp; `nitcad run` prints the provenance with its report. | The first thing to ask of a result is which build made it. | — |
| PV7 | **Out of scope:** deck analyses for transient, small-signal and trace runs (PD1); visualization formats (VTK/HDF5, Unit 25); cross-machine checksums (6.8); any change to numerics, defaults or deck syntax beyond PD1's ruling. | Keeps the unit to records and files. | Transient and AC files are reachable from the API only, unless PD1 brings the deck in. |

**Owner decisions needed:**

| # | Question | Recommendation |
|---|---|---|
| PD1 | Add `[transient]`, `[small_signal]` and `[trace]` deck analyses (with their outputs) in this unit, or as the next one? | **Next unit.** Each brings its own deck keys (time stepping, frequencies, arc length) and checks; this unit stays reviewable. |
| PD2 | Start time and wall time in files: two runs of one deck then give files that differ in those bytes. | **Include them.** Tests compare loaded values with provenance excluded, not file bytes. |
| PD3 | Host name and user name in provenance? | **Exclude.** Files get shared; CPU model and OS are what reproduction needs. |
| PD4 | Provenance in the CSV table (`out.table`)? A `#` header line breaks some readers; a sidecar file is one more output. | **Leave the CSV unchanged**; the binary sweep file and the runner's report carry the provenance. |
| PD5 | The version stays `0.0.0` until a release policy exists? | **Yes**; the revision identifies builds until then. |

**Acceptance gates, fixed before any implementation:**

| # | Gate |
|---|---|
| PG1 | **No numerical change.** Every suite passes unchanged; the Unit 40 assembly traces (P-core pinned) are bit-identical to `main`'s; `input_identity` of fixed fixtures equals the value on `main`. |
| PG2 | **Round trips.** `Transient` and `SmallSignal` round-trip bit for bit, including NaN, −0 and infinities, empty and present `stopped` / `unfinished`, snapshots, traps, electrothermal and tunnelling fields; provenance round-trips in every kind that carries it; a `trace_bias` sweep round-trips. |
| PG3 | **Compatibility.** Sweep, checkpoint and device files written by `main` (committed as fixtures) load unchanged with an empty provenance; a loader given another content kind refuses it (`invalid_input`). |
| PG4 | **Corruption.** Truncation, a checksum mismatch and appended bytes are refused for kinds 5 and 6, as for kinds 1–4. |
| PG5 | **Provenance is true.** The revision equals `git rev-parse HEAD` of the build and `dirty` follows the working tree; the MKL version is present only with PARDISO loaded; the cache sizes equal Unit 41's constants; wall time ≥ 0 and the start time parses. |
| PG6 | **Suites.** Full Release and Debug suites pass at every commit; no warnings under `/W4 /WX`. |
| PG7 | **Mutation checks.** Hand mutants of the new writers, readers and provenance gathering (a dropped field, a swapped key, a wrong content kind, provenance hashed into the identity) are caught. |
| PG8 | **No scope creep.** Nothing outside PV1–PV6 and PD1's ruling; ARCHITECTURE.md 6.6 and 11.3 updated to "as built" only after the owner accepts the result. |

**Commits (proposed):** 1) `Provenance` and the build-info header, filled by `make_run_record` (PG1, PG5); 2) provenance
in sweep and checkpoint files (PG3); 3) transient and small-signal files (PG2, PG4); 4) the runner's `--version` and
report (PV6). Each passes PG6 on its own.

**Risks:**
- Tests that compare whole `Sweep` / `Transient` values for equality would see the new time fields differ between
  runs; they must compare with provenance excluded (PG1 catches any that do not).
- A build from a source archive has no git: the revision is `unknown` and `dirty` false, which PG5 must accept.
- The build-time git step must not force a rebuild when nothing changed (checked by a second, no-op build).

**Owner approval (2026-10-17):** "go ahead", taken as approval of the proposal with every recommendation:
PD1 the deck analyses are the next unit; PD2 start and wall time are in files; PD3 no host or user name; PD4 the
CSV table is unchanged; PD5 the version stays `0.0.0`.

**As built, so far:**
- Commit 1 (`4831f75`): `results::Provenance` in `RunRecord`; `cmake/build_info.cmake` writes the generated
  `NiTCAD/build_info.hpp` at every build, only when its text changes (a second build recompiles nothing);
  `linalg::eigen_build()`; `solve::provenance()`; every `make_run_record` attaches it and the four public runs
  stamp `wall_s`. The backend and its thread count are already settings, so they are not repeated (PV1 refined).
  Release 19/19. Mutants: 8 applied, 8 caught (the identity hashing the start time; a run's wall time dropped,
  for trace and for transient; cache sizes swapped; the start time to whole seconds; a wrong revision; a wrong CPU
  count; the Eigen version dropped). One further mutant, `%T` replaced by `%H:%M:%S`, is equivalent: `%S`
  prints the fractional seconds too.
- Commit 2: provenance as run-record key 3 and checkpoint key 7 (`io::Checkpoint::provenance`), every field
  optional on read. Fixtures `tests/io/fixtures/pre_unit42/` (a 1D bar run by `nitcad` built from `main` at
  `67e8fe5`) load with an empty provenance (PG3). `.gitattributes` keeps the binary fixtures unconverted.
  Mutants: 7 applied, 7 caught (provenance not written in the run record or in the checkpoint; a field not read;
  provenance required on read; the dirty flag inverted; wall time not written; its key colliding with the
  settings').

- `b5940c8`, a test fix: the commit-2 Debug suite failed 4 assertions, all the revision check. It required the
  stamped revision to equal git's HEAD when the test runs, and commit 2 was made after that build. A build
  older than HEAD is normal, so the revision must now be a 40-digit commit of HEAD's history
  (`git merge-base --is-ancestor`), or `unknown` without git (PG5 refined). The same commit adds PV5's test that
  a `trace_bias` sweep round-trips through `save_sweep`.
- Commit 3 (`397b08c`): content kinds 5 and 6, `save_transient` / `load_transient` and `save_small_signal` /
  `load_small_signal`. Each point's admittances are stored as one array (frequencies × contacts² complex values,
  as (real, imaginary) pairs) and its resolutions as one (frequencies × contacts), so a loaded run's
  `admittance(point, frequency, i, j)` stays in bounds. Loading also refuses a snapshot of a missing time point,
  an order outside 0–2, and small-signal fields of a missing contact (PG4). Mutants: 11 applied, 11 caught (error
  ratio, trap occupancy, small-signal fields or resolutions dropped; real and imaginary parts swapped; the
  admittance size check loosened, caught only after a too-long case was added to the test; the snapshot and
  contact bounds off by one; the order unchecked; the wrong content kind; the lagged count written as the
  retraced one). A first form of the snapshot mutant did not compile and was replaced.
- Commit 4 (`d8186bb`): `nitcad --version` prints the build on one line; `nitcad run` prints the build, the
  platform, the start and the wall time after the files it wrote (for an equilibrium run, which has no run
  record, the provenance of now). New ctest `cli_version` checks the line. Mutants: 3 applied, 3 caught.
- **Suites (PG6):** commit 1 Release 19/19; commit 2 Release 19/19, Debug 12 of 13 (the revision check above);
  commit 3 Release 19/19, Debug 13/13; commit 4 Release 20/20, Debug 14/14. No warnings under `/W4 /WX`. After
  the stamp, a second build recompiles nothing.
- **Debug assertion dialogs:** a mutant that makes a test read a failed `std::expected` opens the MSVC Debug
  runtime's modal assertion dialog, which waits for a click. Later mutant runs used the Release build.
- **PG1 trace gate, method:** Catch2 3 orders test cases randomly per run, and fixtures cached in function-local
  statics (for example `fixture_curve()` in `mos_test.cpp`) are traced inside whichever case first uses them. So
  traces are compared with `--order decl`, with both builds configured with the same MKL runtime (PARDISO cases
  run, not skipped), each pinned to P cores, block by block (a block is one test case's trace).
- **PG1 result:** every one of `main`'s 362 blocks over the nine traced entries is byte-identical in this branch's
  traces (assemble 107, deck 18, solve 199, breakdown 7, BTBT 3, MOSFET 11, transient 6, unstructured 3,
  validation 8), `result.identity` included. The only blocks without a counterpart are this unit's three new solve
  cases. No numerical result and no input identity changed.

**Status:** accepted by the owner (2026-10-17, "accepted"); ARCHITECTURE.md 6.6 and 11.3 updated to "as built"
(PG8). Merged to `main` as `145f23c` (PR #52, 2026-10-17).
