<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 40 — generic fields and equations behind the assembler (proposed 2026-10-16; authorized and accepted 2026-10-17; on `assemble/fields-equations`)

**Owner request (2026-10-16):** document a refactor toward ARCHITECTURE.md 6.12.1–6.12.2 (the equation graph, 2.1) in
three slices; documentation only on branch `assemble/fields-equations`; no implementation until explicitly authorized.
Scope: this refactor only. No new physics, models, solver strategies, deck keys or file formats.

**Why now (audited on `main` at `f153d57`).** `DriftDiffusion` is one class with the whole system in it:
`src/assemble/drift_diffusion.cpp` (2036 lines) and `src/assemble/electrothermal.cpp` (1136 lines). The unknown layout
is the hard-coded stride `m_` (3: ψ, n, p; 4 with θ last), read as `m_ * i + k` (and `m_ * m_` blocks) on 114 lines of `drift_diffusion.cpp`
and as literal `4 * i + k` throughout `assemble_heat`; `stride()` is read by `solve/bias.cpp`, `solve/fields.hpp`,
`solve/small_signal.cpp`, `solve/trace.cpp`, `solve/transient.cpp` and `electrothermal_test.cpp`. Each physics unit
added private structs and an `add_*` pass to the class (`ImpactNode`, `TunnelNode`, `PathLayout`, `ThermalNode`, …).
Several targets need the layout and the terms to be explicit: Unit 30's open finding (a natural-test norm monotone in
reverse bias "needs the system's unknown layout"), block preconditioners for iterative solvers (field blocks), state
transfer for mesh refinement (Unit 17; field location and units), and rendering (Unit 25; named fields with units).

**The current accumulation order (to be preserved exactly; `DriftDiffusion::assemble`, isothermal path):**

1. Per node in index order: a Dirichlet node (ohmic contact) sets its three rows and diagonal and skips the rest; an
   insulator node sets its rows; otherwise the Poisson charge `-V (n - p - N)` is **assigned**, then incomplete
   ionization **re-assigns** the Poisson row and its two entries, then the gate term is added, then SRH, Auger,
   radiative, then the BDF storage term.
2. The interface edges, in `InterfaceEdges` order (the insulator node's Poisson row, then the semiconductor node's
   three rows).
3. The mesh edges, in mesh edge order: Poisson only on edges without carriers; otherwise Poisson, electron and hole
   fluxes into row a, then row b, skipping Dirichlet rows; interface (`charged`) edges skipped.
4. Impact ionization, local band-to-band tunnelling, nonlocal path tunnelling, each in its own node or path order.

The electrothermal model takes a separate pass (`assemble_heat`, stride 4, its own order); `small_signal_matrix`
calls `assemble` with `interfaces = false` and adds the interface terms itself; `EquilibriumPoisson` is a separate
assembler. Floating-point addition is not associative, so each value of F and J is a sum whose order is fixed by this
sequence; any slice that changes the order of the additions into one entry changes results in the last bits.

**Design (proposed):**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| FE1 | **Slice 1 — field descriptors and the per-node unknown map.** A `Field` descriptor: name, units, location (node only in this unit), unknown or fixed, region applicability (semiconductor, insulator, all). A `DofMap` built from the active fields: the unknown index of (field, node), the per-node block size, and the field of an unknown index. For every existing configuration it reproduces today's layout exactly: node-interleaved, ψ, n, p, then θ with the electrothermal model, every field at every node (insulator nodes keep their n, p unknowns held at 0). `stride()` stays and returns the block size; the solve layer and tests move to the map only where they read a field's slot. | The layout becomes data instead of `m_ * i + k` and `4 * i + k` literals; the norm (Unit 30), preconditioner blocks and field output can ask for a field's unknowns. Same layout, so checkpoints (Unit 29), run identities and the Jacobian pattern are unchanged. | Insulator nodes still carry unused density unknowns (dropping them changes the layout: deferred, FD5). Edge, face and cell fields are not in this unit. |
| FE2 | **Slice 2 — equation-term interfaces.** A term declares the fields it writes and reads, its domain (nodes, mesh edges, interface edges, contacts, paths), its Jacobian pattern contribution at creation, and an `add(state, residual, values)` that adds its residual and exact partials at positions located once at creation (as today's `ab`/`ba`/`block_`). The existing passes become terms, each keeping its own loop body and order: node terms (Dirichlet, insulator, Poisson charge with the ionization override, gate, SRH, Auger, radiative, storage), interface-edge terms, mesh-edge terms (Poisson, Scharfetter–Gummel electron and hole fluxes), impact ionization, local and nonlocal tunnelling; and the heat equation and thermal variants of `assemble_heat`. | The 2.1 graph's local, edge and interface/boundary terms as objects; a new equation is a new term, not a new pass in the class. Exact Jacobians stay hand-written (6.12.3: no automatic differentiation in this unit). | Dispatch overhead (FG4 measures it): edge and global terms are called once per pass over their domain, node terms once per node inside the fused node loop (FE3, FR2), never per matrix entry. |
| FE3 | **Order rule for slice 2.** Every entry of F and J receives its contributions in today's sequence with today's operands: node terms in the order of step 1 (assign → ionization re-assign → gate → SRH → Auger → radiative → storage), with today's skips (a Dirichlet or insulator node gets no later node term); then interface-edge, mesh-edge and global terms in steps 2, 3, 4, each over its domain in today's index order. A node term only writes its own node's rows, so the node terms may run fused in one node loop or as one loop per term without changing any entry's order; the fused loop is kept unless FG4 says otherwise (it is today's memory access pattern). Assignments stay assignments (the Dirichlet rows, the Poisson charge and its ionization override); no term is rewritten as "zero then add" where today's code assigns. | Bit-identical results (FG1) need every entry's additions in the same order with the same operands. | Sequential assembly only: concurrent assembly, where several edges add into one row, changes the order (deferred, FD6). |
| FE4 | **Slice 3 — composition.** An assembler that builds the active term list and the `DofMap` from the device and `PhysicsModels`, in the FE3 order, and evaluates F, J, the time-step and small-signal forms over it. `DriftDiffusion` stays as a compatibility wrapper over it with every public member unchanged in signature, preconditions and errors (`create`, `evaluate`, `residual`, `storage`, `small_signal_matrix`, `terminal_currents`, the thermal and tunnelling members, the Newton hooks); the solve layer, the deck, `io` and the tests compile and run without change beyond FE1's map reads. | Composition per 6.12.2 (diode = Poisson + continuity; thermal = those + heat) without moving the solve layer; the wrapper keeps the change reviewable and revertible. | Two names for one system until a later unit retires the wrapper (not this unit). |
| FE5 | **Diagnostics stay where they are in this unit.** `edge_currents`, `terminal_current_resolution`, `edge_flux_terms`, `band_edges`, `gate_charges`, `contact_charges`, `trap_occupancies`, `thermal_contact_heat` and the weighted-current inputs keep their current code paths, reading the layout through the map. | They are derived quantities, not residual terms; moving them is not needed for composition and multiplies the bit-identity surface. | They are not yet generic derived fields (6.12.1's "derived" kind is deferred, FD7). |
| FE6 | **`EquilibriumPoisson` is out of scope** unless the owner rules otherwise (FD1). | It is a separate, small (284 lines) assembler with its own layout (ψ only); it composes from the same Poisson terms later. | Two Poisson implementations remain for now. |
| FE7 | **One branch; a preparatory commit, then one commit per slice, each passing every gate on its own.** Commit 0 adds only the measurement tools, with no assembler change: the FG4 benchmark and the FG1/FG2 recorder (hidden Catch2 cases, not in ctest, that write F/J bit patterns and run results to a path given in an environment variable, and compare against such a file when one is given). The baselines are recorded at commit 0, which is numerically `main`. Then slice 1 (FE1), slice 2 (FE2; split into 2a, the isothermal terms, and 2b, the electrothermal terms, so each review is about half the size), slice 3 (FE4). The owner may stop after any slice; slice 1 is useful alone (FG7). | The benchmark and the recorder do not exist on `main` (the existing `tests/solve/solver_timing.hpp` times the linear solver, not assembly), so "before" must be measured with the same tools as "after". Each commit is revertible and reviewable; a failed bit-identity check names the slice that caused it. | Slices 2 and 3 touch most of `drift_diffusion.cpp` and `electrothermal.cpp`; reviews are large even when split. |

**Acceptance gates, fixed before any implementation:**

| # | Gate |
|---|---|
| FG1 | **Bit-identical numerics.** For every configuration the existing suites exercise (Boltzmann and Fermi–Dirac; every recombination, mobility, ionization, band-set and gap option; gates, insulators, interfaces and traps; heterojunctions and thermionic emission; impact ionization, local and nonlocal tunnelling; electrothermal; steady, transient, small-signal and `trace_bias`; tensor-grid and unstructured meshes in 1D, 2D and 3D), F and every Jacobian value equal the pre-refactor values **bit for bit** (compared as `std::bit_cast<std::uint64_t>`, not with a tolerance), as do the time-step and small-signal matrices, `bias_derivative` and `small_signal_currents`. Checked by the commit-0 recorder (FE7): at each configuration it evaluates F and J at fixed states built without the code under change (`EquilibriumPoisson`'s equilibrium, out of scope, FE6, with a set bias, plus a deterministic perturbation of every unknown so no row is at a solution) and writes their bit patterns; after each slice the same cases are compared entry by entry, and the first differing row and column are reported. No second copy of the old assembler is kept in the tree. |
| FG2 | **Identical iteration counts and results.** Every Newton iteration count, step-cut count, line-search λ, arc-length step, transient step sequence and accepted/rejected step, every terminal current, weighted current and field, every small-signal admittance, and every run identity digest equal those recorded at commit 0, compared as bit patterns by the same recorder over the solve-layer cases of the existing suites; every existing pin passes unchanged. |
| FG3 | **FD-Jacobian checks.** Every existing FD-Jacobian gate passes unchanged; in addition each term (FE2) has its own FD check of its partials at a non-equilibrium state, under the same tolerances as the existing gates. |
| FG4 | **No measurable assembly slowdown.** A Release benchmark (`[.benchmark]`, not in ctest) times `evaluate` (residual and Jacobian) and `residual` alone on: the legacy 1D diode, the 2D MOSFET with field mobility, a 3D unstructured device, the electrothermal diode, and a nonlocal-tunnelling diode. Measured at commit 0 and after each slice on this machine, in the same session and power state, the median of at least 31 repetitions after warm-up; the gate is the post-slice median not above the pre-refactor median by more than the threshold of FD2 (proposed: 2%, or the measured run-to-run spread of the medians if that is larger, reported alongside). The end-to-end MOSFET and transient benchmarks' wall times are reported, not gated. |
| FG5 | **Suites.** Full Release and Debug suites pass at every slice commit (18/18 and 12/12 as last recorded on `main` after Unit 38, to be re-run at commit 0 before relying on it); zero warnings under the existing flags (`/W4 /WX`). |
| FG6 | **Mutation checks.** For every slice, mutants applied by hand to the new code are caught by the suites **or by the FG1/FG2 comparison** (a reordering mutant, such as two node terms swapped, changes only the last bits, which the suites' tolerances do not see; the comparison is the gate for those, and the report says which of the two caught each mutant). At least: the `DofMap` swapping two fields, an off-by-one block size, the ionization override turned into an addition, two node terms swapped in order, an edge term skipping row b's Dirichlet check, a term's pattern entry dropped, a thermal term's θ column moved, the small-signal pass including interface terms twice, and the wrapper ignoring a model switch. Each is listed with the test that caught it; a surviving mutant is a missing gate, added before the slice is committed. |
| FG7 | **Slice 1 stands alone.** After slice 1 alone: FG1–FG6 pass; `stride()` and every public signature unchanged; the field of every unknown index and the unknown of every (field, node) are queryable; the assembler's block size comes from the map (one source of truth: `m_` is set from it; the literal `4` and `16` of the electrothermal code, where the block is 4 by construction, move with that code into terms in slice 2b); the solve layer's slot reads (`bias.cpp`, `small_signal.cpp`, `trace.cpp`, `transient.cpp`, `fields.hpp`) go through the map. The per-access `m_ * i + k` arithmetic inside the assembler's loops stays until slice 2 moves those loops into terms (rewriting all of it in slice 1 would make slice 1 as large as slice 2). |
| FG9 | **Layering.** The field, map and term types live in `assemble`; `solve` reads only the public map queries (FD4); no new dependency from `physics`, `mesh` or `linalg` on `assemble`; the existing header-boundary checks pass. |
| FG8 | **No scope creep.** No change to physics, models, defaults, deck syntax, file formats, run-record fields or public `DriftDiffusion` signatures; `ARCHITECTURE.md` 2.2, 6.12.1–2 and 11.3 updated to "as built" only after the owner accepts the result. |

**Risks:**

| # | Risk | Mitigation |
|---|---|---|
| FR1 | Bit identity is fragile: a reordered sum, a hoisted common factor (e.g. computing `V * rate_scale_` once instead of per term), a fused multiply-add the compiler forms differently after inlining across a virtual boundary, or a temporary in a different precision changes the last bits and then the iteration counts. | FE3's order rule; terms keep today's expressions verbatim, operand order included; FG1 compares bit patterns at many states before any iteration-level comparison. The flags stay `/fp:precise` with no `/arch` (x64 SSE2 scalar code, `CMakeLists.txt:44`); MSVC is not expected to form contractions under `/fp:precise` without `/fp:contract` [unverified for this toolset: FG1 is what proves it]. If FG1 fails for a reason other than a reordering in the code, the slice stops and the owner is told, rather than tolerances being loosened. |
| FR2 | Indirection costs: virtual dispatch or `std::function` per node or per edge would slow assembly measurably on large 3D meshes; the fused node loop (FE3) needs a per-node call per node term. | Edge and global terms: dispatch once per term per pass, the loop inside the term. Node terms: a per-node call through a list fixed at creation, or static composition (a fixed tuple or variant of terms) if FG4 shows the calls cost (FD3). |
| FR3 | The ionization override and the Dirichlet assignments are not additive; a generic "terms add" interface would break them silently. | Assignment terms are explicit in the interface (a node term may set its rows); FG6 has the override-as-addition mutant. |
| FR4 | The electrothermal pass (`assemble_heat`) duplicates the isothermal terms with T-dependent variants and its own order; folding both into one term list risks differences in either. | Separate thermal term variants with `assemble_heat`'s own order; FG1 runs on electrothermal configurations explicitly. |
| FR5 | Small-signal assembly calls `assemble(…, interfaces = false)` and adds its own interface terms; transient storage and trap elimination thread through `TimeStep`. A missed special case changes AC or transient results only. | FG1 covers the time-step and small-signal matrices; FG2 covers transient step sequences and AC admittances. |
| FR6 | Nonlocal tunnelling changes the Jacobian pattern after creation (`set_paths` → `locate_positions`); a term list with positions fixed at creation would go stale. | The path term relocates its positions on `set_paths` as today; a mutant that skips relocation (FG6, the dropped pattern entry) must be caught. |
| FR7 | Review size: slices 2 and 3 rewrite most of 3200 lines. | One commit per slice, slice 2 split in two (FE7); the recorded baselines (FG1, FG2) instead of a duplicate old assembler; the owner may stop after any slice. |
| FR9 | The FG1 states are only as good as their coverage: a term that is zero at every recorded state (e.g. impact generation below its current floor, a tunnelling path below `matters`) is not checked by bit comparison. | Each configuration's recorded states include one where each of its optional terms is active (nonzero contribution checked when recording); the recorder fails at commit 0 if a term is inactive at every state of its case. |
| FR8 | Unused density unknowns on insulator nodes are baked into the layout; a later layout change (FD5) will change checkpoints and run identities. | Not done in this unit; recorded so the later unit plans a checkpoint version bump. |

**Unresolved decisions (need an owner ruling before implementation):**

| # | Question | Proposal |
|---|---|---|
| FD1 | Include `EquilibriumPoisson` in slice 3 (built from the same Poisson terms) or leave it separate? | Leave it separate in this unit (FE6). |
| FD2 | The slowdown threshold for FG4. | 2%, or the measured spread of the medians if larger, with both reported. |
| FD3 | Dispatch mechanism: virtual terms, or static composition (variant/tuple) chosen at creation? | Virtual dispatch once per term per pass; static composition only if FG4 fails. |
| FD4 | Public or internal: is the `Field`/`DofMap`/term API public in `include/NiTCAD/assemble/`, or internal to `src/assemble/` with only the map queries public? | Internal terms; public read-only `DofMap` queries (the solve layer and later the norm and preconditioners need them). |
| FD5 | Drop the density unknowns on insulator nodes (a layout change)? | No; a later unit with a checkpoint version bump. |
| FD6 | Parallel assembly (edges or nodes split across threads, which changes the order of additions into shared rows and so is not bit-identical)? | Not in this unit; a later unit with tolerance-based gates the owner sets then. |
| FD7 | Derived fields (currents, band edges, charges) as generic fields? | Not in this unit (FE5). |
| FD8 | Where the recorded baselines (FG1, FG2) live: scratch only, or a committed reference file? | The recorder (hidden test cases) is committed at commit 0; the recorded files stay in scratch, regenerated from commit 0 on demand. Bit patterns are specific to this toolset and flags, so a committed file would fail spuriously on CI or after a toolset change. |
| FD9 | Keep the recorder after the unit (a standing bit-identity check for later refactors), or remove it at the end? | Keep it, hidden and out of ctest. |

**Self-review (2026-10-16, before the owner's review).** Checked against the code on `main` (`f153d57`): the
accumulation order, `small_signal_matrix`'s `assemble(…, false)` call, the five solve-layer `stride()` readers, the
build's `/fp:precise` and the run identity (inputs only, so unchanged by a refactor) are as stated. Corrected from the
first draft: (1) node terms may run as separate loops without reordering any entry (FE3); (2) no assembly benchmark or
recorder exists on `main`, so a commit 0 adds them before the refactor, and baselines come from it (FE7, FG1, FG2,
FG4); (3) a duplicate old assembler kept "test-only" was replaced by recorded bit patterns (FG1, FR7); (4) reordering
mutants are invisible to the suites' tolerances and are judged by the comparison (FG6); (5) FG7 asked slice 1 to
rewrite all 114 lines of index arithmetic, which is slice 2's work, now narrowed; (6) the toolset-specific bit
patterns make a committed reference file unsafe (FD8). Added: slice 2 split into isothermal and electrothermal
commits (FE7), the layering gate (FG9), the coverage risk of the recorded states (FR9), FD9.

**Owner ruling (2026-10-17):** "go ahead" after the self-review, taken as authorization to implement with the
proposals of FD1–FD9 accepted as written.

**Commit 0 as built (the recorder differs from FG1's wording above, and is stronger).** Instead of F and J at fixed
states built for the purpose, the recorder traces **every** value the assembler's public members return during the
existing suites:
- `include/NiTCAD/assemble/observer.hpp`, `src/assemble/observer.cpp`: a null-by-default observer. `DriftDiffusion`'s
  public members hand it their inputs and outputs under a tag: `evaluate` and `residual` (state, residual, Jacobian
  values, and a time step's rate, storage and traps), `small_signal_matrix`, both `bias_derivative`s,
  `small_signal_currents` (columns, values, bias), `state_from_potential`, `stamp_contacts`, `update_size`,
  `apply_update`, `storage`, `trap_occupancies`, `edge_currents`, `terminal_currents`, `terminal_current_resolution`,
  `edge_flux_terms`, `band_edges`, `gate_charges`, `contact_charges`, `conduction_currents`, `interface_trap_charges`,
  `thermal_contact_heat`, `electrical_power`, `edge_joule_heat`. Unset, each call pays one null test; the hooks only
  read, so no value changes.
- `tests/baseline/trace_listener.cpp`, a Catch2 listener linked into the assemble, solve, deck and validation test
  executables: with `NITCAD_TRACE_FILE` set, one line per observed value (index in its section, tag, length, a 64-bit
  FNV-1a hash of its bit patterns) under each test case and section. Child processes (the deck tests' `nitcad.exe`)
  are not traced.
- Each traced ctest entry is run with the variable set (scratch script; the files stay in scratch, FD8); a slice
  passes FG1 and FG2 when its traces equal commit 0's line for line. The solve layer is unchanged by every slice, so
  equal assembler outputs at every call give equal iteration counts and results; the first differing line names the
  test, section and call. It covers every term the suites activate, at the suites' own states.
- **FR9's planned mitigation was not built.** The plan was for the recorder to fail at commit 0 when a term is
  inactive at every recorded state. Because the recorder traces only the states the suites reach, a term the suites
  never activated would go unchecked bit for bit. In its place:
  - The existing per-model FD-gate cases evaluate each model where it is active: impact ionization, local and
    nonlocal tunnelling, gates, interfaces and traps, the electrothermal fixtures, time steps. This was checked by
    reading those tests, not measured.
  - From slice 2a, the terms test requires every term to be active in its fixtures and checks its partials.
- The FG4 benchmark: `tests/solve/assembly_timing.hpp` and five `[.assembly_benchmark]` cases (the legacy 1D diode,
  the 2D MOSFET with field mobility, the 3D generated diode, the electrothermal pn bar, the Zener diode with nonlocal
  paths) in the solve test files that already build those devices.

**Commit 0 findings (2026-10-17, commit `4e8a3af`, this machine: Intel Core Ultra 5 225F, 6 P + 4 E cores):**
- Catch2 v3 runs test cases in random order by default ("Randomness seeded to …"), so traces are recorded with
  `--order decl`; in random order one-time static work (e.g. a cached reference solution) is traced under
  whichever case runs first.
- Determinism: two full traced runs (A, B) of the same commit-0 binaries agreed on every traced entry except
  `solve_unstructured`, whose three cases (an 80 000-node tensor reference, 241 803 unknowns) differed in B from the
  first observed value on (the printed results agreed to every printed digit). Later full runs showed the same
  kind of split in `solve` (after a complex solve), `solve_btbt` (after a Newton linear solve) and `solve_mosfet`:
  always exactly two outcomes per entry, always starting right after a linear solve. **Cause found:** the
  outcome depends on the core type the process runs on. `btbt: a curved junction tunnels before a planar one`,
  pinned to each logical CPU in turn, gives one trace on CPUs 0, 1 and 6–9 (the P cores) and another on CPUs 2–5
  (the E-core cluster); an unpinned process lands on either. (The earlier pinned probes happened to use P cores
  only.) The likely mechanism [unverified]: Eigen sizes its dense kernels' blocking from the cache sizes of the
  core it first queries (the E-core cluster's L2 differs from a P core's), and SparseLU's supernodal updates then
  sum in a different order. **This is a reproducibility defect on `main`** (ARCHITECTURE.md 6.8), outside this
  unit's scope: results depend in the last bits on the core a run starts on, on hybrid CPUs. A fix (e.g. fixing
  Eigen's cache sizes at start-up) changes last bits of existing results, so it needs its own owner-approved unit.
  **For this unit's gates:** every trace run is pinned to the P cores (affinity mask 0x3C3), and a slice's trace
  must equal the pinned baseline's exactly.
- FG4 noise: three back-to-back runs of the benchmark on commit 0 gave medians spread by 20–40% per device (e.g.
  the MOSFET's `evaluate` 3818–4655 µs, the 3D diode's 1706–2274 µs): the machine's run-to-run variation is far
  above FD2's 2%. FG4 is therefore measured by alternating runs of the commit-0 test binary (kept in scratch) and
  the slice's, several times back to back, comparing the minimum of each device's medians; the gate is the slice's
  minimum not above commit 0's by more than 2% plus the spread of commit 0's minima across the alternation, with
  every number reported.
- An existing timing gate, `btbt: the cost of the nonlocal paths` (`t_on < 2 t_off + 0.05`), failed once when nine
  traced test processes ran at once (0.88 s against 0.63 s); it passes under `ctest -j 4`. Not changed in this
  unit; noted as a load-sensitive test.

- Commit 0b (`a7ddc6b`): the public solves also hand their results (fields, bands, currents, charges, admittances,
  convergence records, run identities) to the observer, so the slices' solve-layer conversions are compared bit
  for bit, not only through the tests' tolerances.
- FG4 as measured: commit 0b's and the slice's `[.assembly_benchmark]` runs alternate four times, each pinned to
  CPU 0 (a P core), and the minimum per device is compared. Pinned, the spread of commit 0b's four runs is 1–10%
  per device (unpinned it was 20–40%).

**Slice 1 results (fields and the unknown map):**

| Gate | Result |
|---|---|
| FG1/FG2 | Pinned traces: two commit-0b baselines identical to each other on all nine traced entries; slice 1 identical to them on all nine (the only extra lines are the two new `fields` test cases). |
| FG3 | Every existing FD-Jacobian gate passes unchanged (no term moved in slice 1). |
| FG4 | `evaluate`, slice 1 against commit 0b (minimum of four pinned runs): MOSFET 2056 against 2038 µs (+0.9%), 3D diode 665 against 650 (+2.2%; commit 0b's own spread 9.7%), electrothermal bar +0.4%, 1D diode +0.6%, Zener +0.9%; `residual` −1.8% to −0.1%. Within the gate (no assembly loop changed). |
| FG5 | Release 18/18, Debug 12/12, no warnings. |
| FG6 | 7 mutants, all caught: electron and hole fields swapped in the map (fields test), `DofMap::index` off by one (fields test), the temperature slot pointing at the holes (electrothermal test), `scaled_state` writing electrons into the hole slot (small-signal test; the run also crashed), the transient error estimate reading the potential for the rise (electrothermal test), `converged_temperature` reading the potential (electrothermal test); the small-signal residual bound reading the electron slot for the potential survived the suites and was caught by the trace (only the `result.ac_resolution` lines differ). |
| FG7 | `stride()` and every public signature unchanged; `DriftDiffusion::dofs()` answers the field and node of every unknown and the unknown of every (field, node); the block size is set from the map; the solve layer reads every slot through it (`solve/fields.hpp` `Slots`, `scaled_state`). |
| FG9 | The map is in `assemble` (`fields.hpp`); `solve` reads only its public queries. |

**Slice 2a as built (the isothermal terms; `src/assemble/terms.hpp`, `terms.cpp`, internal per FD4):**
- `DriftDiffusion::Term`: `info()` (name, domain, the fields whose rows it writes and those its partials are taken
  in, whether it assigns its rows), `pattern()` (the entries it needs beyond every node's full block and the edges'
  stencil), `add(system, evaluation)`. Terms hold no data: they read the system's at each call, so a system stays
  copyable and shares its terms (`std::shared_ptr<const Term>`). `DriftDiffusion::Terms` holds the term types, the
  factory, the pattern hook and a one-term evaluation for the tests.
- Fifteen terms, in FE3's order: ohmic contact rows, insulator rows, Poisson charge, incomplete ionization, gate,
  SRH, Auger, radiative, storage (nodes); interfaces (interface edges); Poisson flux, carrier fluxes (mesh edges);
  impact ionization, local and nonlocal band-to-band tunnelling (global). The loop bodies are the old ones, moved;
  `add_impact_generation`, `add_local_tunnelling` and `add_path_tunnelling` became the three global terms.
  `assemble` is the term loop; with `interfaces = false` (small-signal) it skips the interface-edge terms.
- **One loop per node term, not one fused node loop (FE3's default).** FD3's dispatch (one virtual call per term
  per pass) would otherwise have been a call per node per term. FE3 allows separate loops: a node term writes only
  its own node's rows, so every entry receives its additions in the old order. FG4 decided: no slowdown (below).
- **The mesh-edge pass is split** into the Poisson flux and the carrier fluxes. The two write disjoint entries (ψ
  rows; n and p rows), so the order of the additions into each entry is unchanged. The electron and hole fluxes
  stay one term, because `edge_fluxes` computes them together.
- **The pattern.** `create` builds the edge table before the pattern, adds each term's own entries (interfaces'
  insulator-row density columns, impact ionization's other carrier, the paths'), then locates the positions as
  before. `set_paths` rebuilds the pattern from the base plus every term's entries. The pattern is the same set, and
  `from_triplets` sorts, so it is identical. Which terms exist depends on the models alone (the impact and tunnel
  node tables are built after the pattern); a term whose table is empty adds nothing.
- `tests/assemble/terms_test.cpp` (FG3): each term's partials against central differences of its own residual, at
  non-equilibrium states across five fixtures: a diode with every node model under both statistics, steady and with
  a time step; the oxide/electrode/trapped-interface fixture, steady and with a time step; a lumped gate; impact
  ionization at 3e5 V/cm; local and nonlocal tunnelling at 4e6 V/cm. Each runs in 1D and 2D (the gate in 1D only),
  with the tolerance of the Unit 19/20 gates (1e-6, the best of steps 1e-5 to 1e-7). Each fixture requires its terms
  to be active (FR9). A second case checks the order, the domains and which terms assign. The column scale is the
  whole system's Jacobian, the existing gates' scale. Scaled by a term's own column, the vanishing side of a
  Scharfetter–Gummel flux at a drop of hundreds of V_T (an entry near 1e-128 beside a flux near 300) is below what
  any difference step resolves, and reads as error 1.

**Slice 2a results:**

| Gate | Result |
|---|---|
| FG1/FG2 | Pinned traces identical to slice 1's on all nine traced entries. |
| FG3 | Every existing FD-Jacobian gate passes unchanged; the 15 terms pass their own checks (230 assertions). |
| FG4 | Minimum of four alternating pinned runs, slice 2a against slice 1. `evaluate`: MOSFET 2028 against 2033 µs (−0.2%), 3D diode 647 against 664 (−2.7%), electrothermal bar −0.1% (not a term yet), 1D diode 11.9 against 13.2 (−10%), Zener 11.1 against 12.1 (−8.6%). `residual`: MOSFET +1.5% (slice 1's own spread 4.8%), others −3.2% to −0.1%. No slowdown. |
| FG5 | Release 18/18, Debug 12/12, no warnings. |
| FG6 | 9 mutants, all caught by the suites. Ionization override turned into an addition: heterojunction test. SRH and Auger swapped in order: the terms test's order check (the swap itself only changes last bits). Carrier fluxes ignoring row b's Dirichlet check: electrothermal test. The interface pattern without the hole column, and the paths' pattern entries dropped: the contract check in the position lookup. The impact pattern on the wrong edges: impact test. `assemble` ignoring `interfaces = false`, so small-signal adds the interface terms twice: small-signal test. The Poisson flux ignoring an electrode row: insulator test. The gate term left out: gate test. |
| FG9 | The terms are in `src/assemble` (internal); no new dependency; `solve` unchanged. |

**Slice 2b as built (the electrothermal rows; `src/assemble/thermal_terms.cpp`):**
- `assemble_heat` stays the entry point, because `thermal_contact_heat` calls it with the sinks out. It is now the
  loop over the electrothermal term list, which `Terms::make` builds when the model is on. The evaluation carries the
  thermal states, the levels and the sinks flag.
- Twelve terms, in `assemble_heat`'s order:
  - heat sinks (assigns the heat row: the isothermal sink's Dirichlet row, R_th, or 0)
  - lattice heat storage
  - ohmic contact rows at the node's T
  - insulator rows (an electrode at its T)
  - the Poisson charge (the isothermal term, shared; complete ionization, U6)
  - gate with the p-polysilicon shift
  - SRH, Auger and radiative at the node's T
  - carrier storage (the continuity rows and the stored carriers' energy in the heat row)
  - interfaces
  - edge fluxes (Poisson, carriers and heat)
  - impact ionization with gamma(T)
- Per node, the heat row is written before the electrical rows, and the carrier storage adds to it last. With one
  loop per term, every entry still receives its additions in that order.
- **The edges stay one term.** The heat rows take the carrier fluxes and their partials, so splitting would compute
  `thermal_fluxes` twice.
- **The interface term runs even for small-signal.** `assemble_heat` never skipped interfaces, and the electrothermal
  `small_signal_matrix` adds no interface terms of its own.
- The literal 4 and 16 moved with the code into the terms (FG7). Every block is 4 × 4 by construction.
- The edges' full 4 × 4 pattern blocks moved from `create` into the edge term's `pattern()`, so the pattern is the
  same set. `make_thermal`'s positions and the diagnostics (`add_thermal_storage`, `bias_derivative`,
  `thermal_contact_heat`) keep their code (FE5).
- `terms_test` gains:
  - a heated 1D and 2D diode in radiative silicon, with an oxide and a recombining fixed-charge interface in 2D,
    under both statistics, steady and with a time step;
  - a p-polysilicon gate;
  - impact ionization at about 3e5 V/cm;
  - the electrothermal order;
  - a held heat row being its sink's Dirichlet row in a time step too.

  The last was added for two mutants that survived the suites (below).

**Slice 2b results:**

| Gate | Result |
|---|---|
| FG1/FG2 | Pinned traces identical to slice 2a's on all nine traced entries. |
| FG3 | Existing FD-Jacobian gates unchanged; all 27 terms pass their own checks at 1e-6 and are active (478 assertions in 4 cases). |
| FG4 | Minimum of four alternating pinned runs, slice 2b against 2a. `evaluate`: electrothermal bar 194.6 against 192.7 µs (+1.0%, slice 2a's own spread 7.1%); `residual` +0.9%. The isothermal devices −3.0% to −0.2% (unchanged code). Within the gate. |
| FG5 | Release 18/18, Debug 12/12, no warnings. |
| FG6 | 8 mutants, all caught. Heat-sink assignment turned into an addition: electrothermal test. A θ column moved (the contact row's ψ entry into the p column): electrothermal test. Thermal SRH and Auger swapped: the terms test's order check. The edges' pattern without the heat rows: the position lookup's contract check. The interfaces term left out: the terms test's activity check. The gate's p-polysilicon shift dropped: electrothermal test. Lattice storage, and separately the carriers' stored energy, added to a held (isothermal-sink) heat row: **both survived the suites**. At a converged step the held rise equals the sink's before and after, so the term vanishes at the solution, and no fixture held a node that was neither contact nor insulator. The new held-row case (a sink on y_min) catches both. |

**Slice 3 as built (composition; `src/assemble/composition.hpp`, `composition.cpp`):**
- `DriftDiffusion::Composition` owns the unknowns' map and the active terms, with each term's equation sets, and
  runs the evaluation loop over them. `DriftDiffusion` holds one, immutable and shared by copies, in place of
  `dofs_` and the term list, and forwards `dofs()`, `assemble` and `assemble_heat` to it. Every public signature is
  unchanged (`dofs()` moved out of line).
- **Equation sets after ARCHITECTURE.md 6.12.2.** Poisson, continuity, contacts, oxide, gate, interfaces,
  generation and heat. They form a bit set: the electrothermal edge term belongs to Poisson, continuity and heat
  at once.
- **Composition tables.** The two `make` functions became tables in `terms.cpp` and `thermal_terms.cpp`. Each row
  is a term's sets, an activation predicate over the models and the device data `create` has set by then
  (contacts, insulators, gates, interfaces), and a factory. Table order is FE3's assembly order. `compose` picks
  the table by the electrothermal model, keeps the active rows, and builds the map from the fields the equations
  solve for.
  - Diode = Poisson + continuity + contacts.
  - MOS adds oxide, gate and interfaces.
  - Thermal adds heat and a fourth unknown.
  - The generation models add generation.
- The contact-row and insulator-row terms are now active only when the device has such nodes. Before, they were
  always in, with empty loops. No value changes.
- **Not done:** the model data the terms read (per node, edge and contact) stays in `DriftDiffusion`, because its
  diagnostics read it too (FE5). `DriftDiffusion` is therefore both the compatibility wrapper (FE4) and that
  data's store. Moving the data out and retiring the wrapper is a later unit's work, as FE4 already says.
- `terms_test` gains three composition cases:
  - the equation sets and block sizes of a diode, a MOS structure with an oxide and interface, a lumped-gate
    MOS-C, a generation model, and a thermal diode;
  - every model switch (SRH, Auger, radiative, incomplete ionization, impact, local and nonlocal tunnelling; the
    electrothermal SRH, Auger, radiative and impact) adds exactly its term;
  - the electrothermal small-signal matrix at s = 0 equals the steady Jacobian entry for entry, so its interface
    terms stay in.

**Slice 3 results:**

| Gate | Result |
|---|---|
| FG1/FG2 | Pinned traces identical to slice 2b's on all nine traced entries, every test case except the `terms:` and `composition:` cases. Those trace their own calls to `evaluate`, so they change with the term lists of their fixtures (an oxide-free fixture no longer carries an empty insulator-rows term). |
| FG3 | Existing FD-Jacobian gates unchanged; the term checks pass. |
| FG4 | Minimum of four alternating pinned runs, slice 3 against 2b, measured twice. `evaluate` −0.9% to +2.3% (1D diode, 11.4 against 11.7 µs; +0.9% in the first measurement), `residual` −2.1% to +0.3%. The first measurement's 3D-diode `residual` +2.5% was +0.3% when repeated. Within the gate; the only code-path change is one indirect call per evaluation. |
| FG5 | Release 18/18, Debug 12/12, no warnings. |
| FG6 | 7 mutants. Six caught: the Auger switch ignored (drift-diffusion test); the electrothermal table never chosen (crash); the `interfaces` flag ignored, so small-signal adds the interface terms twice (small-signal test); radiative and storage swapped in the table (the order check); the temperature field dropped from the map (electrothermal test); the gate predicate false (gate test). The seventh (the electrothermal interface term skipped without `interfaces`) **survived and was equivalent**: `assemble_heat` always passes `true`, so the check it mutated was dead code. The dead check was removed. The new electrothermal small-signal case guards that behaviour. |
| FG8 | No change to physics, models, defaults, deck syntax, file formats, run records or public signatures. `ARCHITECTURE.md` 2.2, 6.12.1–2 and 11.3 are to be updated to "as built" after the owner accepts the result. |
| FG9 | Fields and the map public in `assemble`; terms and composition internal (`src/assemble`); no new dependency; `solve` reads only the map. |

**Review (2026-10-17, owner request: review all eight commits against `main` before push authorization):**

| Area | How | Result |
|---|---|---|
| Term order and bodies | `main`'s `assemble`, `add_*`, `assemble_heat` and `add_thermal_*` against the 27 terms, statement by statement after normalizing only the mechanical renames (`s.`, `ev.`, `stride`/`m` for `m_`, `at(s, ev, …)`) | Every expression unchanged operand for operand. The only other lines are skip conditions equivalent to the old if/continue chain, helpers with the old lambdas' bodies, and wrapping. The tables' order is `main`'s sequence |
| Activation predicates | Each against the guard it replaced | Equivalent. The impact, tunnel and path tables are filled only under the same model switches. `InterfaceEdges::empty()` is `edges_.empty()`. The pattern entries are the same set |
| `create`, `set_paths` | Read line by line | Positions are located the same way in the same pattern; paths rebuild from the same entry set |
| Public API | Headers diffed against `main` | Nothing public removed or changed. Additions: `dofs()`, `fields.hpp`, `observer.hpp` |
| Bit identity against `main` | `main` built in a separate worktree. Six decks run with `nitcad run` from `main` and from HEAD, pinned to a P core: a diode with Fermi–Dirac, field mobility and incomplete ionization; an oxide with an electrode and a trapped interface; a lumped gate; electrothermal; impact ionization; nonlocal tunnelling | Output tables byte-identical. This adds direct steady-state evidence for the link from commit 0 back to `main`, which until then rested on the hooks only reading. Transient and small-signal (which the deck runner does not run) still rest on that, and on the suites |

Six findings, none a blocker, fixed at the owner's word ("do 1"):
1. **The terms' field declarations were unchecked.** New test case: every term writes only the rows and Jacobian
   columns of the fields it declares, across every fixture. Before it was committed, it caught a deliberately
   wrong declaration (the gate term declaring the density rows).
2. The five `*_unobserved` solves had external linkage and no declaration. They are now in anonymous namespaces.
3. `Term`, `Terms` and `Composition` were public nested names. They are now private to `DriftDiffusion`. The
   library and its tests name them through `TermAccess`, declared only in the internal `src/assemble/terms.hpp`.
4. `dofs()` on a moved-from system dereferenced a null composition. It now has a contract check
   (`NITCAD_EXPECTS`).
5. Electrothermal under Fermi–Dirac computes the equilibrium product in each of SRH, Auger and radiative (once
   before). Measured with a new `[.assembly_benchmark]` case, the 401-node pn bar under Fermi–Dirac, against
   `main` (the case was added there temporarily, then removed): `evaluate` 339.9 against 342.8 µs (−0.8%),
   `residual` −3.0%. No slowdown.
6. The FR9 wording above was corrected: the planned recorder check was not built, and the text now says what
   covers its concern instead.

After the fixes: pinned traces identical to slice 3's on all nine traced entries (the `terms:` and `composition:`
cases aside, as before); Release 18/18, Debug 12/12, no warnings; FG4 against slice 3 with six alternating pinned
runs, `evaluate` −0.6% to +0.7%, `residual` −1.6% to +2.8% (the 3D diode; slice 3's own spread 14.6%). A first
four-run measurement read up to +6.3%, but with both binaries about 15% slower than usual, and did not repeat.

**Status:** commit 0, commit 0b and slices 1, 2a, 2b, 3 done. Unit 40 is complete on `assemble/fields-equations`,
reviewed, and the review's findings fixed. **Accepted by the owner (2026-10-17, "accept and push Unit 40")**;
`ARCHITECTURE.md` 2.2, 6.12.1–2 and 11.3 updated to "as built" (FG8). The hybrid-CPU reproducibility
defect (commit 0 findings) is tracked separately as Unit 41 (`investigate/hybrid-cpu-determinism`).
