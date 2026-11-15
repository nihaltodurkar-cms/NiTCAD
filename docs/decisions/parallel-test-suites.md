<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Parallel test suites (owner, 2026-10-23, "the test is taking too long", then "do the test-speed"; on `test/parallel-suites`)

The Debug suite ran its 19 tests one after another in about 40 minutes (`solve` 769 s, `render` 741 s,
`app_window` 410 s), on a 10-thread workstation: the presets set no job count.
- **One ctest entry per test file** (`nitcad_add_file_tests` in CMakeLists.txt) for the seven long executables
  (`mesh`, `assemble`, `render`, `solve`, `deck`, `app`, `app_window`): `<name>.<file>`, each running
  `"-#" "[#<file>_test]~[.]" --allow-running-no-tests`. Catch2's sharding was not used: with the default random
  order and no fixed seed its shards may overlap or miss tests (Catch2 3.16 warns so), and a fixed seed or
  declaration order would give up the random order. Each file keeps it. Hidden tests stay out; the Release-only
  hidden-tag entries are unchanged; `solve.mosfet` runs nothing (all its tests are hidden). **Coverage:** per
  executable, the files' test-case counts sum to the plain run's: 95, 107, 78, 215, 43, 57, 20.
- **Four jobs** in both test presets (`"execution": {"jobs": 4}`; CI's `ctest --preset` takes it). The Debug
  suite is 2,037 s of work; four jobs already hold the rest beside its longest entry, and fewer jobs load the
  timing gates less. `render.image` (its 2D and 3D MOSFET fixtures solve for about 530 s in Debug before any
  view is drawn; splitting its sections would only repeat them) has `COST 600`, so it starts first without
  ctest's timing record.
- **The WG1 files run alone** (`RUN_SERIAL`): `app.run_controller`, `app_window.self_test`,
  `app_window.analysis_window`, and in Release `app_window_3d`. Run beside the others (six jobs, then four),
  the self-test passed once and failed once with gaps of 494 and 611 ms, its 16-point solve at 95 s against
  about 42 s alone. AC9 (`solve.hole_mobility`, a cost ratio within one process) passed in parallel in every
  run and is left parallel.
- **Suites (working tree on `1be857e`):** Debug 99 of 99 in 947 s (was about 2,400 s); Release 107 of 107 in
  618 s (was 1,013 s at `3d0e78a`). **Clean worktree at `15c1088`:** builds from scratch, Release 94 s and Debug
  56 s (vcpkg's binary cache), no warnings under `/W4 /WX`; Release 107 of 107 in 393 s; Debug 98 of 99 in 850 s:
  `app_window.self_test`, run alone after `render.image`, failed WG1 with one gap of 296.6 ms, its solve at the
  usual speed (16 points in 35.4 s) and the monitor's plot drawn in 1.8 ms: the WG1 spike open above (Unit 26c),
  not load.
- **Hosted CI at `1af3b55` (PR #62; the first CI to finish since Unit 26b — PR #61's runs and `main`'s were
  cancelled):** Release passed (21 min); Debug 97 of 99 (1,929 s), both failures in entries run alone. (1) The
  26c·3 playback test's wait for its transient MOSFET run timed out: `pump_until`'s 120 s, a guard against a hang,
  is shorter than that run on the runner's Debug build. The guard in `analysis_window_test.cpp` is now 600 s; the
  file passes locally (139 s). (2) The GPU self-test (WARP on the runner) failed WG1 with a gap of 114.1 ms, its
  solve at 69.3 s (about twice this workstation's).
- **WG1 on the hosted runner (owner, 2026-10-23, "c": a looser CI-only bound):** 250 ms where the environment
  variable `CI` is set (GitHub Actions sets it), 100 ms everywhere else, as ruled. `nitcad-app --self-test` takes
  `--gap-limit-ms <ms>` (default 100; more than 0) and names the bound in its WG1 failure; the app does not read
  the environment. The tests take the bound from one place, `app_test::wg1_limit_ms()` (tests/app/fixtures.hpp),
  for all three WG1 checks: the self-tests (passed as `--gap-limit-ms`), the run controller and playback. New
  checks: `--gap-limit-ms 0` and a missing value are usage errors (64); a 1D diode run (701 points) under a
  0.001 ms bound fails WG1 naming it. Mutant: the bound ignored (100 ms fixed) — caught. Release: the three WG1
  entries pass with and without `CI` set; no warnings under `/W4 /WX`.
- **Hosted CI at `be4d95a`:** Release passed (push and pull-request runs, 21 and 28 min); the three WG1 entries
  passed under the 250 ms bound, and the playback test under its 600 s guard. Debug 98 of 99 in both runs: in
  `app_window.window` the wait for the example MOSFET run (the run monitor test) timed out at `window_test.cpp`'s
  own `pump_until`, still 120 s (the entry failed at 171.5 s and 178.6 s). That guard is now 600 s, as in
  `analysis_window_test.cpp`.
- Not done: restructuring `render.image` (2D and 3D views as separate test cases) to bring Debug near 10
  minutes; worktrees run at once still share their `%TEMP%` scratch folders.

