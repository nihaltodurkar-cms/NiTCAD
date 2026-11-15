<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 41 — last-bit reproducibility on hybrid CPUs (opened 2026-10-17; fix decided by the owner 2026-10-17; on `linalg/hybrid-cpu-determinism`)

**Owner request (2026-10-17):** track the P/E-core last-bit difference found during Unit 40 (its commit 0
findings) as its own unit, separate from Unit 40. First a minimal reproducer and the root cause. Do not change any
numerical result or implement a fix until the impact, the performance cost and the reproducibility policy have been
reviewed.

**Reproducer (`benchmarks/hybrid_cpu/`, standalone, Eigen only, not in the CMake build).**
- `repro.cpp` prints the executing core's CPUID leaf-4 cache sizes, the sizes Eigen adopted, and bit hashes
  (FNV-1a) of a 300 × 300 dense product and a SparseLU solve of a 22 500-unknown nonsymmetric five-point operator.
- Run pinned to each logical CPU of this machine (Intel Core Ultra 5 225F, 6 P + 4 E cores):

| Cores | CPUID L1d / L2 | Eigen's adopted L1 / L2 | SparseLU hash, `auto` | `fixed` (any one set) |
|---|---|---|---|---|
| P: CPUs 0, 1, 6–9 | 48 KiB / 3 MiB | the same | `317ad696…` | one hash on all ten CPUs |
| E: CPUs 2–5 | 32 KiB / 4 MiB | the same | `50df5dd6…` | the same hash |

- The 300 × 300 dense product agrees on every core. A 1200 × 1200 product does not.
- `diode_2d.deck` (a 121 × 81-node 2D diode, `nitcad run` from `main`): 6 of the table's current values differ
  between a run pinned to CPU 0 and one pinned to CPU 2, by at most 1.1e-15 relative. CPU 0 against CPU 6 (two P
  cores) is identical. Six small decks (1D, and a small 2D oxide structure) are identical on P and E cores. Their
  systems are too small for the blocking to change.

**Root cause.**
- Eigen (5.0.1, vcpkg; the backend's static_assert pins it) detects the cache sizes once per process, with CPUID leaf 4 on the core that runs the first
  matrix product, and keeps them in a static.
- The blocking of its matrix products is computed from them. SparseLU's supernodal updates (`panel_bmod`,
  `kernel_bmod`) are such products.
- On a hybrid CPU the two core types report different L1 data caches, so the products split their sums
  differently and round differently in the last bits. Newton then amplifies the differences from the first linear
  solve on.
- **The L1 size alone decides it.** On a P core, `fixed 32768 3145728 …` (only L1 changed to the E core's) gives
  the E-core hash. On an E core, `fixed 49152 …` gives the P-core hash. L2 at 3 or 4 MiB makes no difference.
- So the result depends on the core type the first product runs on. That is usually the start core, but a thread
  can migrate before then. It does not depend on the machine's load or on the core of later products.
- [unverified] The MKL PARDISO backend (Unit 18) dispatches by CPU too. Whether it differs between core types is
  not measured; MKL's conditional numerical reproducibility settings (`MKL_CBWR`) would be the place to look.

**Measured (2026-10-17, at the owner's request; no code, flag or result changed).** The tooling:
`benchmarks/hybrid_cpu/measurement.patch` is a measurement-only patch, not committed to any build. It applies to
Unit 40's tree, which has the trace listener, and:
- adds full values for the `result.*` tags to the trace;
- adds `NITCAD_EIGEN_CACHE=l1,l2,l3`, which fixes Eigen's cache sizes before any test runs;
- adds a hidden `[.unit41_cost]` case timing SparseLU on the MOSFET and a 20³ 3D diode.

`compare-values.sh` compares two value traces block by block. The suites were traced three ways, each pinned:
P cores (mask 0x3C3), E cores (0x3C), and E cores with the cache sizes fixed to the P-core values.

*Impact on the suites:*
- **Every test passes on the E cores as on the P cores.**
- **Fixing the cache sizes removes the difference completely.** On E cores with the P-core sizes, all nine traced
  entries are bit-identical to the P-core run. Eigen's L1-dependent blocking is the whole cause.
- **Which entries differ.** With Eigen's own detection, 4 of the 9 entries are identical on P and E (assemble,
  deck, transient, validation) and 5 differ. In `solve`, `solve_unstructured` and `solve_breakdown` the
  differences are in values. In `solve_mosfet` and `solve_btbt` they also change the work done. Ten of the MOSFET
  cases and two tunnelling cases ran a different number of evaluations, i.e. Newton iterations or bias cuts. For
  example, the curved-junction tunnelling case made 1825 calls on P and 1558 on E; the mesh-independence MOSFET
  case 224 and 204.
- **Size of the changes.** Fields (potential, densities, band edges) agree to at most 4e-15 relative of each
  array's largest value, even where iteration counts changed.
- **Terminal currents.** Every difference is at most 4.3e-5 of the larger of the current and its own reported
  resolution (`terminal_current_resolution`). Currents more than ten times their resolution change by at most
  3.1e-8 relative: Newton stops at a different iterate within its tolerance. Small-signal admittances: the
  dominant entries agree to 1e-11 relative; entries at the AC resolution (about 1e-11 S against 6e-9 S) change sign.
- **Quantities that are rounding by nature** change by large relative amounts, as expected: backward errors, the
  AC resolution, edge currents and quasi-Fermi levels at equilibrium.

*Cost of fixing the cache sizes* (SparseLU refactorization, fastest of three pinned runs; the noise, between two
configurations with the same blocking, is up to 2%):

| Fixed to | P core, MOSFET / 3D diode | E core, MOSFET / 3D diode | Result bits everywhere |
|---|---|---|---|
| nothing (today) | 155.7 / 5771 ms | 180.9 / 6563 ms | P and E differ |
| P-core sizes (L1 48 KiB, L2 3 MiB) | unchanged | 184.1 / 6809 ms (+1.8% / +3.7%) | today's P-core bits |
| E-core sizes (L1 32 KiB, L2 4 MiB) | 154.8 / 5779 ms | unchanged | today's E-core bits |
| the smaller of each (L1 32 KiB, L2 3 MiB) | 155.2 / 5787 ms | 184.6 / 6526 ms | today's E-core bits |

The cost of any fixed choice is within about 4%, at most at the edge of the noise. The choice of sizes decides
which of today's two results survives: P-core sizes keep the P-core runs' bits, which match Unit 40's recorded
traces; a 32 KiB L1 keeps the E-core runs'.

**Options for the owner (option 1 with the P-core values chosen 2026-10-17; below):**
1. Fix Eigen's cache sizes at start-up to one set, for example in `linalg`'s Eigen backend before the first
   factorization:
   - to the P-core values, so runs that start on an E core change in the last bits;
   - or to a 32 KiB L1, so P-core runs change.

   Either makes one machine reproducible whatever core a run starts on.
2. As 1, but take the smallest L1 over all cores, queried per core with thread affinity. The same on this machine
   as the 32 KiB choice, and it adapts to other hybrid parts. Different machines still differ.
3. A solver option (sizes or "fixed"/"detected"), recorded in the run record.
4. Document only: results are reproducible on one core type, and runs may differ in the last bits between P and E
   cores.

For reproducibility across machines, the blocking must not depend on the hardware at all: one fixed set
everywhere (option 1 with constants). Its cost on other CPUs is not measured.

**Owner decision (2026-10-17):** option 1, with the P-core values. The CPU a run executes on must not change its
result; bit-identity across machines is not promised, and tolerances stay the contract (ARCHITECTURE.md 6.8).

**Fix (branch `linalg/hybrid-cpu-determinism`, from `main` 629bed0 with the two investigation commits).**
- `src/linalg/eigen_sparse_lu.cpp`: `fix_cache_sizes()` calls `Eigen::setCpuCacheSizes` once per process
  (`std::call_once`) with L1 48 KiB, L2 3 MiB, L3 20 MiB, at the start of `SparseLuBackend::analyze`, before any
  product of the backend.
- `tests/linalg/eigen_blocking_probe.cpp`, ctest `linalg_eigen_blocking`: starts itself once per logical CPU of
  its affinity mask, each child suspended, pinned, then resumed; each solves a 22 500-unknown nonsymmetric
  five-point system with the default `LinearSolver` and prints a hash of the solution's bits. It fails unless all
  hashes agree. Its own executable, because Eigen detects the cache sizes once per process.

**Verified (2026-10-17):**
- Probe with the fix: one hash (`130858804b42571a`) on all ten CPUs.
- Mutant (the `setCpuCacheSizes` call removed): the probe fails; CPUs 2–5 (E) give `b1405c79cef29076`, the P cores
  `130858804b42571a`, the fixed build's hash. The P-core results are kept.
- Release suite 19/19 (the probe added).
- Not run, at the owner's direction: a whole-suite trace pinned to P against E cores with the fix in place. The
  investigation's E-core run with the P-core sizes was bit-identical to the P cores on all nine entries, but it
  fixed the sizes before any test ran; the fix does so at the backend's first analysis. An Eigen product large
  enough to be blocked that runs before the first factorization would still use the detected sizes. The probe
  does not cover that.
- Debug suite 13/13, the probe among them (its CPUs agree in Debug too).

**Still open:**
- [unverified] PARDISO (MKL) on P against E cores (`MKL_CBWR`).
- The cost of the fixed sizes on CPUs other than this one is not measured.
- [unverified] Whether any Eigen product large enough to be blocked runs before the first factorization (above).

**Status:** fixed; merged to `main` as `67e8fe5` (PR #51, 2026-10-17).
