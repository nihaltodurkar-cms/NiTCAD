<!-- Archived from DECISIONS.md on 2026-11-20 (branch docs/loose-ends), word for word; DECISIONS.md keeps its summary. -->

### Unit 60 — release version policy (PD5) and cross-machine result comparison (6.8) (proposed 2026-11-19; on `results/release-policy`)

**Owner request (2026-11-19):** "start Release version policy (PD5) and cross-machine checksums", the 11.3 row
"Results and provenance", next step "A release version policy (PD5); cross-machine checksums (6.8)". The branch
starts from `main` `c9e3201`.

**Starting point:**
- **PD5 (Unit 42):** "The version stays `0.0.0` until a release policy exists". The version's only source is
  `project(NiTCAD VERSION 0.0.0)`. `cmake/build_info.cmake` writes it into every run's provenance next to the git
  revision and a `dirty` flag, and `nitcad --version` prints it. No tag, release or changelog exists.
- **6.8:** agreement between machines is defined by tolerances, never by bit-identity digests.
  `RunRecord::input_identity` (FNV-1a over the inputs' bits) "identifies a run on one build; it is not a
  cross-machine checksum". The Eigen blocking is CPU-independent since Unit 41, but no result has ever been compared
  between two machines. MKL PARDISO's CPU dispatch is unverified.
- **CI** runs on every branch push and on PRs (GitHub-hosted `windows-2025`, a different CPU from the development
  machine), with the MKL runtime installed. It does not run on tag pushes.

**OLD / NEW / REASON (the name "checksum"):** OLD: the 11.3 row says "cross-machine checksums". NEW: a
tolerance comparison of two result files, plus committed reference results that CI compares on its own machine.
REASON: 6.8 forbids bit-identity digests. A digest of rounded values flips whenever a value lies near a rounding
boundary: with about 1e5 values at 1e-12 relative noise and 6 digits, the chance that some value straddles is about
0.1 per run. A comparison reports how far apart two results are, which a digest cannot.

**Standing conditions:** no solver, assembly or physics code changes, so every result is bit-identical by
construction. Only the version string changes (VQ3). Criteria are fixed before code. Release and Debug suites with
`/W4 /WX`, plus hand mutants. Built with MSVC 14.51. Nothing is pushed, tagged or merged without authorization.

**Left out of this unit:** comparing transient and small-signal files (XQ2; their loaders exist); edge currents and
the rate fields (they are derived from the compared fields); MKL PARDISO references (the CI threads and its CPU
dispatch would be tested together, so a failure would not say which); publishing a release.

**Design: the release policy (RV):**

| # | Rule |
|---|---|
| RV1 | **Semantic versioning** `MAJOR.MINOR.PATCH`. The only source is `project(NiTCAD VERSION …)`. Before 1.0, MAJOR stays 0 and MINOR takes MAJOR's role. 1.0 is the owner's decision. |
| RV2 | **What bumps what (before 1.0):** **MINOR** for a new capability; a changed default that moves a default result beyond `nitcad compare`'s default tolerances (any regenerated cross-machine reference, XC6, counts); a deck key removed or renamed, or a meaning changed; a file format version change. **PATCH** for everything else: fixes that leave every reference matching, docs, build and CI. After 1.0, anything that breaks an old deck, file or public header is MAJOR. A file that loaded before must still load (IO2), so a refusal is never a PATCH. |
| RV3 | **A release** is the owner's decision. A release PR bumps the version and adds the release's row to ARCHITECTURE 9's release table (version, date, the units it contains). After it merges, the owner pushes an annotated tag `vX.Y.Z` on the merge commit. Claude never creates or pushes a tag. |
| RV4 | **Between releases** `main` carries the last released version. A build is identified by the version, the revision and `dirty`, which every run record already holds. |
| RV5 | **The CI tag check:** CI also runs on `v*` tag pushes. A tag build fails unless the tag equals `v` + the CMake version and the build is not `dirty`. The comparison is a `cmake -P` script, `cmake/check_release_tag.cmake`, so it is tested locally (XG5). |
| RV6 | **Release notes** are the tag's annotation, written from DECISIONS' unit summaries. No changelog file (VQ5). |

**Design: cross-machine comparison (XC):**

| # | Slice | Change | Notes |
|---|---|---|---|
| XC1 | 1 | `analysis/compare.hpp`: `compare_sweeps(a, b, CompareTolerances)` returns a `SweepComparison`. First the structure: the same contact count, point count, node count and bias per point (exactly), and the same `stopped` state. Otherwise the result is `not_comparable` with the reason. Then, for each quantity, the largest difference, its ratio to the allowance, the worst point and node, and pass or fail. Input identity and mesh digest are reported as equal or not; they are not gated (XQ6). | `analysis` already depends on `results` only. |
| XC2 | 1 | **Quantities and default tolerances (XQ3):** terminal currents \|a−b\| ≤ 1e-6·max(\|a\|,\|b\|) + max(resolution_a, resolution_b), each run's own floor (Unit 15); weighted currents the same with their resolution, when present; gate and interface-trap charge 1e-6 relative + 1e-30; ψ \|Δ\| ≤ 1e-6 V; n and p 1e-6 relative against a floor of 1e-8 of the point's largest density (update_size's floor); lattice and carrier temperatures 1e-6 K, when present. | Newton stops at 1e-8 of a scaled correction, and its last step is quadratic, so differences from round-off sit far below these. The margin is printed so that a near miss shows. |
| XC3 | 2 | `nitcad compare <a.nitsweep> <b.nitsweep> [--rtol <x>] [--potential-tol <V>]`: one line per quantity (largest difference, ratio, where), the two runs' CPU and revision lines, then `match` or `differ`. Exit 0 match, **4 differ** (XQ5), 1 unreadable or not comparable, 64 usage. | Any `.nitsweep`, trace sweeps included. |
| XC4 | 2 | **References:** `tests/cross_machine/` holds the `.nitsweep` that `examples/pn_diode_1d.deck` and `examples/mosfet_2d.deck` write on the development machine (Eigen SparseLU, the decks' defaults), regenerated only by the documented command in XC6. | XQ4: small, Eigen only. |
| XC5 | 2 | **The ctests** `cross_machine_<deck>` (both presets): a `cmake -P` script copies the deck into the build tree, runs `nitcad run`, then `nitcad compare` against the reference, and prints both CPU lines. CI therefore compares the development machine's results with its own on every push. | |
| XC6 | 3 | **When a reference stops matching** after an intended change, the PR regenerates it (`nitcad run` on the deck, copy the `.nitsweep`), records the old and new largest differences in its DECISIONS record, and counts as MINOR under RV2. | The comparison's report is the evidence. |

**Criteria (frozen once approved):**

| # | Gate |
|---|---|
| XG1 | `compare_sweeps` unit tests: identical sweeps give zero differences and `match`. For each quantity alone, a perturbation at 0.9× its allowance passes and one at 1.1× fails, the report naming that quantity, point and node. The resolution floor is used: a difference of 0.5× the floor on a zero current passes. Every structural mismatch (contacts, points, nodes, bias, `stopped`) gives `not_comparable` with its reason. |
| XG2 | A sweep compared with its own save-and-load copy: zero differences (bit-for-bit files, Unit 29). |
| XG3 | `nitcad compare` exit statuses 0, 4, 1 and 64, each exercised by a ctest. |
| XG4 | `cross_machine_pn_diode_1d` and `cross_machine_mosfet_2d` pass locally and **on the CI runner**. The CI log's CPU line differs from the reference's (recorded); the largest ratios on CI are recorded in the results. |
| XG5 | `check_release_tag.cmake`: a matching tag on a clean build passes; a mismatched tag, a missing `v` and a dirty build each fail (four ctests). |
| XG6 | Nothing else moves: no file under `src/` outside `src/analysis/` and `src/cli/` changes. `cli_version` passes with the new version. Every existing test passes unchanged. |
| XG7 | `/W4 /WX`. Hand mutants (the tolerance's sign or comparison, the floor dropped, one structural check removed, the exit status for `differ`, the tag check's `dirty` test) are each caught. |

**Questions for the owner (recommended option first):**

| # | Question | Options |
|---|---|---|
| VQ1 | Version scheme | **Semantic versioning, 0.y.z before 1.0 (recommended).** / Calendar versions (2026.11). |
| VQ2 | Version on `main` between releases | **The last released version; the revision tells builds apart (recommended).** / The next version with `-dev`. |
| VQ3 | First version | **0.1.0, set in this unit; the tag waits for the owner (recommended).** / Keep 0.0.0 until the first release PR. |
| VQ4 | Who tags | **The owner, an annotated tag on the release PR's merge commit (recommended).** / Claude, on the owner's word. |
| VQ5 | Release notes | **The tag annotation, from DECISIONS' summaries; no changelog file (recommended).** / A `CHANGELOG.md`. |
| XQ1 | What "cross-machine checksum" means | **A tolerance comparison plus committed references checked by CI (recommended).** / A digest of rounded values (rejected above). |
| XQ2 | File kinds compared | **Sweeps, traces included (recommended)**; transient and small-signal later. / All three now. |
| XQ3 | Default tolerances | **As in XC2 (recommended).** / Tighter (1e-9), with a risk of flaky CI. |
| XQ4 | Reference decks | **`pn_diode_1d` and `mosfet_2d`, Eigen only (recommended).** / Add `mosfet_3d` and a PARDISO run. |
| XQ5 | Exit status for "differ" | **4, a new status (recommended).** / 1, shared with failures. |
| XQ6 | Identity and mesh digest unequal across machines | **Reported, not gated (recommended).** Deck parsing is exact, but mesh generation calls CRT math that may take other paths on another CPU. / Gated. |

**Owner ruling (2026-11-19, standing instruction "always choose recommended", 2026-11-11):** VQ1–VQ5 and XQ1–XQ6
take the recommended option. Pushing, the PR, tags and the full Release/Debug suites still wait for the owner's
word.

**Results (2026-11-19):** slices 1–3 are committed locally: `9c1b269` (`compare_sweeps`), `eec8f76` (`nitcad compare`,
0.1.0, the tag check, the ctests), `4efc99c` (the references, written by `eec8f76`'s Release build; its provenance
says "modified" because the owner's uncommitted `.gitignore` edit and two untracked decks were in the tree), and the
records.
- **XG1, XG2:** `[compare]`, 7 test cases (398 assertions), Release and Debug.
- **XG3:** `compare_cli.*`, 7 ctests (0 self, 4 differ with the diode's body doping doubled, 1 not comparable, 1
  unreadable, 64 for one file, a negative `--rtol` and an unknown option).
- **XG4 local half:** `cross_machine_pn_diode_1d` and `cross_machine_mosfet_2d` pass in Release (bit-identical, 0
  differences) and in **Debug against the Release references**. The diode differs at round-off: largest ratio 2.1e-6
  (terminal current 6.0e-12 A/cm²), ψ 1.1e-10, n 4.1e-10, p 2.9e-10. The MOSFET is bit-identical. **The CI half is
  pending** until the branch is pushed.
- **XG4 CI half (recorded 2026-11-20, met):** from the verbose `cross_machine_` step of CI run 44554600656 on
  `main` `38163aa` (Unit 59's merge; the first green `main` run with that step). The references' CPU line reads
  "Intel(R) Core(TM) Ultra 5 225F" (revision `eec8f76`, Release); the runner's reads "AMD EPYC 7763 64-Core
  Processor", so the comparison crosses vendors. **Release:** both decks bit-identical on the runner, every
  quantity's largest ratio 0 (diode 30 currents and 6015 node values per field; MOSFET 64 and 17 056), input
  identity and mesh digest equal. **Debug against the Release references:** the diode's largest ratios are
  terminal current 2.089e-6 (6.045e-12 A/cm², point 4), ψ 1.110e-10, n 4.147e-10, p 2.905e-10, gate charge 0, the
  same as the local Debug run; the MOSFET is bit-identical. In Debug the input identity is reported different (the
  mesh digest equal); it is reported, not compared. Both ctests pass in both configurations.
- **XG5:** `release_tag.*`, 5 ctests (a match, this tree's version, a mismatch, a missing `v`, a dirty tree).
- **XG6:** the branch touches `src/analysis/compare.cpp` and `src/cli/` only under `src/`. `cli_version` passes and
  prints 0.1.0.
- **XG7:** `/W4 /WX` clean. 8 hand mutants, all caught: the current floor dropped, the bias check removed, the
  allowance floor dropped, `differ` exiting 1, a one-sided NaN passing, the tag check ignoring `dirty`, the
  `stopped` check removed, the density floor from the smallest density.
- **Size:** about 710 lines added in code, tests, scripts, CMake and CI (docs aside), plus 2.9 MB of references
  (`mosfet_2d` 2.2 MB, `pn_diode_1d` 0.7 MB).
- **Full Release and Debug ctest not run** (rule: they wait for the owner's word). Only the analysis, CLI,
  cross-machine and tag tests above were run.

**Close (2026-11-20, owner: "do the loose ends from merged units").** Merged as PR #96, `main` `9093cc8`
(2026-11-19, CI green). The full suites ran on CI with the merge's descendants (run 44554600656: Release 248/248,
Debug 227/227) and locally on later `main` trees (Unit 58: Release 250/251, the AC11 timing check only; Debug
230/230). The version is 0.1.0, untagged (only the owner tags). The unit is closed; this record moved to
`docs/decisions/unit-60.md`.
