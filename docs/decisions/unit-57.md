<!-- Archived from DECISIONS.md on 2026-11-18 (branch docs/unit57-as-built), word for word; DECISIONS.md keeps its summary. -->

### Unit 57 — generic fields and equations, completed (proposed 2026-11-15; approved 2026-11-15; reviewed and amended 2026-11-15; started 2026-11-16; on `assemble/generic-equations`)

**Owner request (2026-11-15):** "What next to implement … as per architecture", then "Plan for finish generic
fields and equations (the Unit 40 follow-up)". This unit does what ARCHITECTURE.md 6.12.1–2 and 11.3 list as "Not
yet" after Unit 40: move the model data out of `DriftDiffusion` and retire it as the wrapper, add the derived field
kind with edge locations, and build `EquilibriumPoisson` from the same terms. The branch starts from `main`
`beb494f`. The number is 57: it was freed when the cuDSS plan was dropped (2026-11-15); 56 is reserved by Unit 55.

**Why now:** 11.1's first priority. Every model since Unit 40 (Units 46–54) was added through terms, but also as
members and data in `DriftDiffusion` (header 1,011 lines, source 2,069, about 53 data members the terms read
through `TermAccess`). Density-gradient confinement (11.3, nanoscale) needs a new unknown per carrier. Unit 30's
open natural-test norm needs the layout of each equation. Both are cheaper on a generic system.

**Standing conditions (as Unit 40):** numerically a pure refactor. F, J and every result are bit-identical
(P-core-pinned traces, Unit 41). There is no change to physics, defaults, the deck, file formats, run identities or
checkpoints. Criteria are fixed before code. Release and Debug suites pass at every slice commit with `/W4 /WX`,
plus hand mutants. Nothing is pushed or merged without authorization.

**Left out of this unit:** dropping the insulator nodes' n and p unknowns (FD5: a layout and checkpoint change);
parallel assembly (FD6: not bit-identical); automatic differentiation (6.12.3); the natural-test norm itself
(Unit 30's follow-up, a consumer of this unit); density gradient and any other new physics; the render layer's
\|E\| and \|J\| (Unit 25b; they read results, not the system).

**Design:**

| # | Slice | Change | Notes |
|---|---|---|---|
| GE1 | 0 | **Baselines.** (i) Extend the FG1/FG2 recorder's traced entries to the configurations added after Unit 40: Lombardi, MLDA, Altermatt levels, Schottky contacts, gate tunnelling, energy balance, swept and generated 3D meshes, and an `[adapt]` deck. (ii) **Add trace hooks (`detail::observe`) for what slices 3 and 4 change and the trace does not yet see** (review R3): `EquilibriumPoisson::evaluate` and `residual` (x, F, J), the per-mechanism rates (`RateFields`), the carrier temperatures in K, the gate-tunnelling currents and the energy powers. The hooks only copy values when an observer is set, so no number changes. (iii) **Write GG2's decks that do not exist yet** (R4): the repository has only `decks/finfet_3d*.deck`, `decks/schottky_diode_1d.deck` and the four `examples/`. The new ones go in `tests/decks/` (committed, small, each under a minute in Release). (iv) Extend the FG4 benchmark with an energy-balance and a gate-tunnelling device. Record (i), (ii) and the benchmark on `main`, pinned to P cores. | No assembler change. The recorded files stay in scratch (FD8). |
| GE2 | 1 | **The model data store.** The data fixed at creation moves into an internal `SystemData`, held as `shared_ptr<const SystemData>`. That is the per-node, per-edge and per-contact data, the scaled tables, the node blocks' positions `block_`, the edges' `EdgeTerm::ab/ba`, and the gate tunnels and their positions. The data that changes after creation moves into a value-type `SystemState`. That is the contact biases and metal values, the electrode potentials, gate biases, thermal rises, held carrier temperatures and the minority floor. **It also holds the tunnel paths, their layouts, the base triplets and the Jacobian pattern, because `set_paths` rebuilds all four** (`src/assemble/drift_diffusion.cpp:1824-1841`; review R1). Terms read `(const SystemData&, const SystemState&, const Evaluation&)` and no longer `DriftDiffusion`. The private helpers they call (`edge_fluxes`, `thermal_fluxes`, `level`, `interface_at`, `node_storage`, `path_terms`, `gate_tunnel_fluxes`, `surface_gradient`, …) become functions over the data, in the source file of the model they belong to. `TermAccess` is removed. | A copy shares the data and owns its state, so a relocation in one copy never reaches another. `DriftDiffusion`'s public signatures are unchanged in this slice. |
| GE3 | 2 | **The generic system.** A public `assemble::System`: a composition (terms + `DofMap`) over `SystemData` and `SystemState`, with only the members every equation set has: `create`, `dofs`, `unknowns`, `set_bias`, `make_jacobian`, `evaluate` / `residual` (steady and with a time step), `storage`, the small-signal matrix, `bias_derivative` / `bias_rows`, `update_size` / `apply_update`, `stamp_contacts`, `state_from_potential`. The small-signal path's own interface pass becomes a term, active only in the frequency domain. The model-specific members become free functions in topic headers: `currents.hpp` (terminal, edge, conduction and small-signal currents, resolution, flux terms), `thermal.hpp` (contact heat, electrical power, Joule heat, thermal bias), `tunnelling.hpp` (paths, path states, gate-tunnelling currents), `energy.hpp` (energy powers, hold and release temperatures), `charges.hpp` (contact, gate and trap charges, occupancies). **Queries take `(const System&, x)`. Mutators take `System&`** (`set_thermal_bias`, `hold_temperatures`, `release_temperatures`, `set_paths`, `set_minority_floor`; R6). They reach `SystemData` and `SystemState` through accessors whose parameter is a passkey type only `src/assemble` can construct, so the public header exposes no internal type. Callers in `solve` and the tests move to these (GQ1). The trace labels (`detail::observe` tags) keep their names, so the slice-1 and slice-2 traces compare line for line. | Only `solve` and the tests use `DriftDiffusion`; `deck`, `io` and `app` do not. |
| GE4 | 3 | **Derived fields (FD7).** `Field` gains a kind (unknown, fixed, derived) and the locations edge and interface edge besides node. Shapes stay scalar only: no listed field is a vector, so a vector shape waits for its first user (R9). A derived-field table in `assemble` gives a descriptor and an evaluator `(const System&, x) → values` for each existing diagnostic: band edges and quasi-Fermi levels (node), per-mechanism rates and net rates (node), the electron and hole edge currents and flux terms (edge), trap occupancies (interface edge), and lattice and carrier temperatures in K (node). The solve layer fills `results::BiasPoint`, `RateFields` and `BandDiagram` through the table. | Results structs and file formats are unchanged (GQ4). |
| GE5 | 4 | **Equilibrium from the same terms (FD1).** Equilibrium becomes a composition with a ψ-only `DofMap` and the Poisson, contact, oxide, gate and interface sets. It shares `SystemData`'s scaled tables (both already use `make_scaled_device`). **Shared terms:** only the Poisson flux and the gate term, which compute the same expressions at stride 1. **Equilibrium's own term variants** (R2; `src/assemble/equilibrium_poisson.cpp:140-235`), each keeping its present expressions and order of additions: (i) the charge with carriers slaved to ψ (Boltzmann, Fermi–Dirac, incomplete ionization), assigning the residual and the diagonal; (ii) the insulator rows, which assign 0 to the residual and the diagonal; (iii) the interface edges at Fermi occupancy (`InterfaceEdges::equilibrium`), not the steady-state occupancy; (iv) the contact rows, ohmic at ψ₀ and Schottky at ψ_S, whose diagonal is set to 1 in a final pass after every addition (today's post-pass in `evaluate`). Per entry, the order stays: charge → gate → interfaces → edge flux, then the contact diagonals. `solve::equilibrium`'s signature is unchanged (GQ3). | Removes the separate assembler and its own edge type. The saving is about half of `EquilibriumPoisson`'s 442 lines, since the variants remain. |
| GE6 | — | **Order rule.** As FE3: every entry of F and J receives its contributions in the order and with the operands it has on `main`. `/fp:precise` stays, so moving code between functions or files cannot change a value. The traces check this anyway. | |
| GE7 | — | **Commits.** One branch. Commit 0 is GE1, then one commit per slice, each passing every gate on its own. The owner may stop after any slice; slice 1 is useful alone. Then the docs: ARCHITECTURE.md 2.2, 6.12.1–2 and 11.3 updated to "as built" only after the owner accepts the result. | |

**Criteria (fixed before code):**

| # | Criterion |
|---|---|
| GG1 | **Bit-identical numerics:** on every traced entry (Unit 40's nine plus GE1's), F, every Jacobian value, the time-step and small-signal matrices, `bias_derivative`, every Newton count, step cut, λ, arc step, transient step sequence, current, admittance, field and run identity equal commit 0's bit for bit (P-core-pinned). |
| GG2 | **Byte-identical decks:** the `nitcad run` output tables and every written file (`.nitsweep`, `.nitckpt`, `.nitadapt`) of 12 decks are byte-identical to `main`'s. The provenance record (revision, start time, wall time, R4) is excluded: the files are compared record by record through the `io` reader with that record skipped, and the CSV tables line by line. The decks (GE1 iii): diode, oxide with traps, lumped gate, electrothermal, impact ionization, nonlocal tunnelling, Lombardi + MLDA MOSFET, Schottky diode, gate tunnelling, energy balance, swept 3D FinFET (reduced), and an `[adapt]` deck. |
| GG3 | **FD-Jacobian:** every existing FD gate and per-term check passes unchanged; GE5's slaved-charge term has its own FD check under Boltzmann, Fermi–Dirac and incomplete ionization. |
| GG4 | **No slowdown:** FG4's benchmark with GE1's two devices added; `evaluate` and `residual` medians within 2% of commit 0's, or of its measured spread if larger (FD2). |
| GG5 | **Suites:** full Release and Debug suites pass at every slice commit, with no warnings. |
| GG6 | **Mutants:** at least 8 per slice, each caught by the suites or the traces (the report says which). At least these:
  - a term reading the state's old bias after `set_bias`;
  - **the bias state shared between copies by mistake, so that `set_bias` or `set_paths` on one copy changes another** (R10);
  - the pattern kept in `SystemData`, so that a path relocation is ignored;
  - two terms swapped;
  - the frequency-domain interface term active in steady state;
  - a derived field reading the wrong slot;
  - the equilibrium charge term's Fermi–Dirac branch using Boltzmann;
  - equilibrium missing the interface set;
  - the equilibrium contact diagonal set before the edge flux, not after it;
  - a topic function reading the composition's map of another system. |
| GG7 | **Retired:** `drift_diffusion.hpp` and `TermAccess` no longer exist (GQ1 a). `System`'s public header names no model. A test parses the header's member declarations, not its comments (R11), and requires that none is named for a model: no `thermal`, `heat`, `tunnel`, `path`, `energy`, `temperature`, `schottky` or `gate` in a member name. Adding a term changes only its own file and the composition table. |
| GG8 | **Extensibility (GQ5), two-way coupled (R5):** a test-only equation adds one node field χ to the drift-diffusion composition without editing `System` or any existing term. χ is at every node (`DofMap`'s layout, FD5), held at 0 on insulator and contact nodes. Its row is ∇²χ = ψ − ψ_ref. A second test-only term feeds χ back by adding a source k·χ to the electron and hole continuity rows. The test checks: (i) the `DofMap` and pattern grow by one field; (ii) both terms' partials pass the FD check, including the cross blocks (n, p rows in χ; χ row in ψ); (iii) Newton converges quadratically in its last iterations; (iv) at k = 0 the ψ, n, p fields equal the system without χ to 1e-12; (v) at k ≠ 0 on a 1D bar with a manufactured source, the error falls at second order under two refinements. |
| GG9 | **Layering and scope:** terms, data and state stay internal (`src/assemble`); `System`, `Field`, `DofMap` and the derived-field table are public. No new dependency from `physics`, `mesh` or `linalg`. No change to physics, defaults, deck syntax, file formats or run records. |

**Questions (each with a recommendation):**

| # | Question | Options | Recommendation |
|---|---|---|---|
| GQ1 | The fate of `DriftDiffusion` | (a) Replace it with `assemble::System`, migrate every caller in slice 2, and remove `DriftDiffusion`. (b) Keep `DriftDiffusion` as a deprecated alias for one unit. (c) Keep the name for the generic class. | **(a)**: the callers are `solve` (14 files) and the tests (30), all in-repo; an alias would leave two names for one thing. |
| GQ2 | Where the model-specific queries live | (a) Free functions in topic headers (GE3). (b) Members of `System` (today's shape). | **(a)**: `System` stays generic; a new model adds a header, not members. |
| GQ3 | Equilibrium's public type | (a) `System` created with the equilibrium equation sets; `EquilibriumPoisson` removed, its two users in `solve` migrated. (b) Keep `EquilibriumPoisson` as a thin type over a composition. | **(a)**: one numerical architecture (11.2). Its `carriers` and `charge_neutral_potential` become free functions. |
| GQ4 | Derived fields in results | (a) The table feeds today's result structs; formats unchanged. (b) Results carry a generic list of named fields too (a file-format version bump). | **(a)**: no format change in a refactor; (b) is a later unit, if the app or exports need it. |
| GQ5 | The extensibility test (GG8) | (a) Yes. (b) No. | **(a)**: it is the evidence that density gradient can be added without a second architecture. |
| GQ6 | Slice order | (a) As GE2–GE5. (b) Equilibrium (GE5) before the system (GE3). | **(a)**: GE5 is simpler on the finished `System`. |

**Rulings (owner, 2026-11-15: "choose recommended dont start"):** GQ1–GQ6 all **(a)**, as recommended. GE1–GE7
and GG1–GG9 stand as written. Implementation is not started and waits for the owner's word to begin slice 0.

**Review (2026-11-15; amendments applied on the owner's "do necessary, dont start yet"):** the plan was checked
against the code, and 11 findings were folded into the rows above. No criterion was loosened.
- **R1, GE2:** `set_paths` rebuilds the pattern, base triplets and path layouts, so these belong in `SystemState`,
  not the fixed data.
- **R2, GE5:** equilibrium's contact, insulator, charge and interface rows differ from drift-diffusion's at stride 1,
  so they are kept as variants. The "about 440 lines" claim is corrected.
- **R3, GE1 ii:** the trace did not see the equilibrium assembly, the rates, carrier temperatures, gate-tunnelling
  currents or energy powers, which slices 3 and 4 change.
- **R4, GE1 iii and GG2:** most of GG2's decks did not exist, and the byte comparison had to skip provenance.
- **R5, GG8:** χ was one-way coupled, so it did not show the coupling density gradient needs, and "per semiconductor
  node" contradicted `DofMap`'s layout.
- **R6, GE3:** mutators cannot be free functions over `const System&`. They take `System&` and reach the internals
  through a passkey.
- **R7, start condition (below):** PR #85 changes 10 `solve` files with 29 `DriftDiffusion` references.
- **R8, GE2 notes and GG4:** "copies get cheaper" had no user, since nothing in `src` copies a `DriftDiffusion`. The
  claim and GG4's copy report are removed.
- **R9, GE4:** the vector shape had no field, so it was dropped.
- **R10, GG6:** the "data mutated through a copy" mutant was impossible on `const` data, so it was replaced.
- **R11, GG7:** the header scan reads member names, not comments.

**Start condition (R7):** slice 0 starts after PR #85 (Unit 17c) is merged to `main`, with this branch rebased onto
it. If #85 is still open when the owner says start, the owner is asked whether to rebase onto `mesh/derived-goals`
instead.

**Start (owner, 2026-11-16):** "start Unit 57" while #85 was open (its Release jobs green, Debug running). Asked,
the owner first chose to wait for #85, then asked "cant you do it in Unit 57 without rebase". Checked against
`origin/mesh/derived-goals`: **R7 overstated the overlap.** #85's diff adds or removes no `DriftDiffusion` or
`EquilibriumPoisson` reference; of the files it changes only `src/solve/adapt_adjoint.cpp` (4 lines) uses
`DriftDiffusion`. The rest of the overlap is `DECISIONS.md`, `ARCHITECTURE.md` and `CMakeLists.txt`. So slice 0
starts on this branch as it stands (from `beb494f`, no rebase). When #85 merges, `main` is merged into this branch
and the baselines of what #85 changes (the `[adapt]` deck and the `adapt` traces) are recorded again on the merge
before the next slice is compared against them.

**Slice 0 as built (commit `55f7669`, 2026-11-16; this machine: Core Ultra 5 225F).** No assembler change.
- GE1 ii hooks: `equilibrium.evaluate.{x,residual,jacobian}`, `equilibrium.residual.{x,residual}`,
  `rates.{mechanism,electron,hole}`, `energy_powers`, `gate_tunnelling.{ecb,hvb,evb}`; on the solve side
  `result.{electron,hole}_temperature`, `result.ac_{electron,hole}_temperature`, `result.rates` (each `RateFields`
  vector) and `result.tunnelling_{ecb,hvb,evb}`. Non-empty values in the baseline: 27,050 equilibrium Jacobians,
  37,425 assembler rate vectors, 59,446 result rate vectors, 449 carrier-temperature fields, 221 energy-power records,
  220 assembler and 209 result tunnelling records, so every hook is exercised by the suites.
- GE1 i: `nitcad_adapt_test` links the listener. **Traced entries: 79** (every per-file entry of the assemble,
  solve, deck and adapt executables, the six `solve_*` and two `deck_*` Release entries, `adapt_goal`,
  `adapt_intelligent` and `validation`), each pinned to the P cores (0x3C3) in declaration order through a
  scratch `pin.cmd` wrapper (`start "pin" /b /wait /affinity 3C3`; a bare `""` title is mangled by MSYS). **Two
  full traced runs of the slice-0 binaries (A 13 min, B 14 min) are identical on all 79 entries** (4,437,963
  observed values, 181 MB); the files stay in scratch (FD8).
- GE1 iii: the 12 decks in `tests/decks/` (all converge; the reduced swept FinFET 61–67 s pinned, the others under
  0.7 s). The 1D junction decks use an erfc donor profile on a constant acceptor body: the deck tests'
  `diode_device` (windows `x=-1:1um`, `x=1um:3um`) turned out to give an ohmic, fully compensated bar (its current
  is linear, 7.3 mA/cm² per V), so it was not copied. GG2's comparison is the hidden `[.deck_compare]` case
  (`tests/deck/deck_compare_test.cpp`). **Two runs of the slice-0 `nitcad.exe` agree on all 38 output files.**
- **GG2 and the adaptation report (ruled 2026-11-16, owner "choose recommended": the times are excluded):** the `.nitadapt` report stores measured times (per-pass
  seconds, the cost model fitted to them, the verification and coarse seconds). Between the two identical runs the
  raw report differs in 102 bytes, against about 20 (the provenance) for a sweep file, so GG2's byte comparison
  cannot hold for it with only the provenance skipped. The comparison clears those times as it clears the
  provenance (as ruled; GG2 reads "provenance and, in a `.nitadapt` report, its measured times"); the `[adapt]`
  deck stays in GG2. Slice 1 waits for the owner's word to start.
- GE1 iv: FG4's benchmark gains `Gate tunnelling: assembly cost on the MOS capacitor` (1D and 2D); EC19 is the
  energy-balance device. Baseline, two runs pinned to CPU 0, `evaluate` median in µs (run 1 / run 2): legacy 1D
  diode 12.0/12.1, 2D MOSFET 2642/2494, 3D generated diode 753/761, Zener 12.5/11.7, electrothermal bar 212/213
  (Fermi–Dirac 333/340), energy balance 1D 396/408, 2D 5466/5849, 3D 3489/2771, gate tunnelling 1D 48.1/45.3,
  2D 171/157. The run-to-run spread reaches 26% on EC19's devices, so FG4 is judged as in Unit 40: alternating
  pinned runs, minimum per device.
- Suites: Release 226/227. `adapt_intelligent` failed AC11's timing check once under `-j 6` (a transposed solve
  0.24 ms against one factorization's 0.19 ms on a 73-node diode pass) and passes alone; it is load-sensitive,
  not changed. Debug is left to CI (owner's standing rule).

**Slice 1 (GE2) as built (2026-11-16; owner "start slice 1"; on this branch as it stands, #85 still open).**
- `src/assemble/system_data.hpp` (internal): `SystemData` holds what `DriftDiffusion`'s private section held,
  less the state: the per-node, per-edge and per-contact data, scaled tables, composition, nested types and the
  terms' Jacobian positions. `SystemState` holds the gate nodes (their bias), `psi0`/`n0`/`p0`, the metal values,
  contact biases, thermal rises, held carrier temperatures, the minority floor, and the tunnel paths, layouts,
  base triplets and pattern (R1). `DriftDiffusion` keeps `shared_ptr<const SystemData>` and an owned
  `unique_ptr<SystemState>`; its copy shares the data and copies the state. Its public signatures are
  unchanged; the former inline accessors are out of line.
- Terms take `(const SystemData&, const SystemState&, const Evaluation&)`; `Term`, `Terms`, `Composition` and
  `GateTunnelJob` are `SystemData`'s nested types, and a composition entry's `active` takes both. The helpers
  the terms call are `SystemData` members, defined where they were; those that read the state take it
  (`gate_tunnel_flux(es)`, `thermal_ohmic`, `add_gate_shift`, `assemble_heat`) and the locators take the
  pattern. `TermAccess` is removed; `DriftDiffusion` befriends `SystemData` for `of` / `state_of` alone.
- **As-built note on R1 and the positions.** GE2 lists `block_` and `EdgeTerm::ab/ba` as data, but `set_paths`
  moves every position into the new pattern (`block_`, the edges', the surface stencils', impact ionization's,
  the interfaces', the gate tunnels'). The positions stay in `SystemData`, and `set_paths` moves them in a copy of
  the data (copy on write, `make_shared<SystemData>(*data_)`), so a copy that shares the old data keeps its own
  pattern and positions. `set_paths` runs once per path freeze, so the copy is not on a hot path.
- Tests: two new cases check that a copy owns its state (R10): `band to band assemble: a copy owns its bias and
  paths` (set_bias and set_paths on a copy leave the original's residual, pattern and Jacobian values bit for
  bit; an assigned system takes the copy's) and `electrothermal assemble: a thermal bias moves the sinks of its own
  copy`. Callers outside `src/assemble` are unchanged; the two internal tests rename `TermAccess::` to
  `SystemData::`.
- **GG1:** the 79 traced entries of the slice-1 binaries are byte-identical to slice 0's run A (P-core pinned).
  The trace ran on the build before the long added lines were wrapped (whitespace only) and before the two
  tests were added.
- **GG2:** the 12 decks run with the final `nitcad.exe` agree with slice 0's run A on all 38 files
  (`[.deck_compare]`).
- **GG4:** FG4's benchmark, slice-0 and slice-1 binaries alternating, pinned to CPU 0, minimum per device. Two
  runs each of every `[.assembly_benchmark]` case: every device within 2% of slice 0 on `evaluate` and
  `residual`, or within slice 0's spread, except 1D drift-diffusion + Canali (+3.1%, spread 1.2%), so EC19's case
  was run four more times each: that device then -0.4% / -0.3%; the 3D devices +3.1% to +8.2% on `evaluate`
  against slice 0's spread of 4.2% to 9.9% (in the first two runs 3D drift-diffusion was -0.0%, 623.3 against
  623.2 us); every other device within 2% or faster. Gate tunnelling 1D / 2D: +0.1% / -0.6%. Passed under GG4's
  spread rule.
- **GG6, 9 mutants, all caught** by the `assemble.*` and `solve.*` ctest entries (control 47/47):
  M1 set_bias's changes discarded (terms read the old bias): 31 entries fail; M2 the state shared between copies
  (`shared_ptr`, copy constructor aliases it): `assemble.band_to_band` (the new copy test); M3 `set_paths` does not
  move the positions (the pattern alone relocated): `assemble.band_to_band`, `assemble.terms`, `solve.btbt`,
  `solve.rates`; M4 `set_paths` moves the positions in the shared data (`const_pointer_cast`, no copy):
  `assemble.band_to_band`; M5 SRH and Auger swapped in the table: `assemble.drift_diffusion`, `assemble.terms`,
  `solve.heavy_doping`, `solve.transient`, `solve.small_signal`; M6 electrothermal set_bias leaves the metal value:
  5 entries; M7 set_thermal_bias writes nothing: first survived (no test moved a sink after creation), caught by
  `assemble.electrothermal` after the thermal-bias copy test was added; M8 hold_temperatures keeps no
  temperatures: `assemble.energy_balance`, `solve.energy_balance`; M9 set_minority_floor a no-op:
  `assemble.drift_diffusion`. The slice-2 to slice-4 mutants of GG6's list (frequency-domain term, derived slot,
  equilibrium) belong to those slices.
- **GG5:** full Release suite (`-j 6`, `/W4 /WX`, no warnings) **227/227 pass** in 3,193 s on the final tree.
  Debug is left to CI (owner's standing rule).
- **Incident (slice 0's scratch):** a leftover first version of slice 0's `run_decks.sh` (with `start ""`,
  which MSYS turns into `\\`) was still running from 12:00 and showed "Windows cannot find '\\'" once per deck as
  the owner pressed OK. Stopped at 13:48. It could only rewrite the console logs (`.out`) of run A, which GG2 skips;
  the data files matched.

**Slice 2 (GE3) as built (2026-11-16; owner "start slice 2"; on this branch as it stands, #85 still open).**
- `include/NiTCAD/assemble/system.hpp` (renamed from `drift_diffusion.hpp`, so the header comment keeps its
  history): `assemble::System` with GE3's generic members only, plus the counts and step data every equation
  set has (`node_count`, `contact_count`, `stride`, `time_scale`, `trap_slots`, `storage_width`, `TimeStep`) and
  `bias_nonlinear` (the companion of `bias_derivative`'s state overload). `DriftDiffusion` is gone; every
  reference to it in `solve` (10 files) and the tests (31) now names `System`.
- Topic headers, each a list of free functions over `(const System&, ...)` or, for the mutators, `System&`:
  `currents.hpp` (edge, terminal, conduction and small-signal currents, `CurrentRow`, the resolution, the flux
  terms), `thermal.hpp` (`electrothermal`, `set_thermal_bias`, the sinks, `heat_scale`, contact heat, electrical
  power, Joule heat), `tunnelling.hpp` (`tunnelling`, `trace_paths`, `set_paths`, `paths`, `path_states`,
  `paths_agree`, `gate_tunnelling`, `GateTunnellingCurrents` and their currents), `energy.hpp`
  (`energy_balance`, `EnergyPowers`, `energy_powers`, hold, release, held), `charges.hpp` (trap occupancies,
  contact, gate and interface trap charges). **As-built placements GE3 does not name:** `rates` with
  `RateMechanism` and `MechanismRate` in a sixth header, `rates.hpp`; `band_edges(const System&, x)` in the
  existing `band_edges.hpp` (beside `BandEdges`, with `System` declared, not included); `set_minority_floor` /
  `minority_floor` as free functions in `system.hpp` after the class (a Newton setting, not a model). The two
  internal constants `impact_current_floor` and `matters` are `SystemData`'s (no caller outside `src/assemble`).
- **Passkey (R6):** `System::data`, `state` (const and not) and `replace_data` take a `const SystemKey&`;
  `SystemKey` is declared in the public header and defined in `src/assemble/system_data.hpp` with a private
  constructor that only `SystemData` can call (`SystemData::of`, `state_of`, `replace`). `System` no longer
  befriends `SystemData`. The public header names `SystemData`, `SystemState` and `SystemKey` only as
  declarations.
- **The bodies:** the model queries are `SystemData` members taking the state (`terminal_currents(st, x)`, ...),
  defined where they were; `src/assemble/topics.cpp` holds the free functions, each one line through the
  passkey. Three are free functions in their model's file: `terminal_current_resolution` and `energy_powers`
  call the system's own `residual` (so its trace hooks fire as before), and `set_paths` replaces the data
  (copy on write, GE2's as-built note). The expressions are unchanged: the bodies were moved by rewriting
  `data_->` and `state_->` only.
- **The frequency-domain term:** small_signal_matrix's interface pass is the term
  `Terms::InterfacesFrequency` ("interface traps, frequency domain"), moved with its expressions and order. A
  composition-table row now has a `frequency` flag; the composition keeps such terms in a list of their own
  (`frequency_terms`), so the steady list (`Terms::list`, the per-term checks, `rates`) is unchanged.
  `Term::add_frequency` adds to the complex matrix (nothing by default); the term's steady `add` is a contract
  failure. `Composition::assemble(..., frequency)` leaves out, in the frequency domain, the steady terms of a
  domain a frequency-domain term has (the interface edges), and `assemble_frequency` runs after the storage, as
  the pass did. The electrothermal table has no such term (the model refuses traps and interface charge), so
  its small-signal matrix is as before.
- The relocation loop's stand-in in `tests/solve/btbt_test.cpp` gives the same five functions as free
  functions (found by ADL), as `solve::detail::solve_relocating` now calls them.
- **GG7:** `system: the public header names no model and the old type is retired` (`tests/assemble/fields_test.cpp`):
  `drift_diffusion.hpp` does not exist; no `TermAccess` in the internal headers; the header (comments blanked)
  has no `DriftDiffusion`; in `class System`'s body every declared name (an identifier followed by
  `( ; = { , ) [`: functions, data members, nested types, their members and the parameters) contains none of
  GG7's eight words; and the parse finds the 17 members it must (so an empty parse cannot pass).
- **GG1:** the 79 traced entries of the slice-2 binaries (P-core pinned): 77 byte-identical to slice 1's
  recorded traces. The other two, `assemble.band_to_band` and `assemble.electrothermal`, hold slice 1's two copy
  tests, added after slice 1's trace was recorded; traced again with slice 1's final binaries they are
  byte-identical to slice 2's. The trace ran before the GG7 test was added (it calls no assembler).
- **GG2:** the 12 decks run with the final `nitcad.exe` agree with slice 0's run A on all 38 files
  (`[.deck_compare]`).
- **GG4:** FG4's benchmark, slice-1 and slice-2 binaries alternating, pinned to CPU 0. After two runs each,
  five 1D and 2D devices were 2.7% to 3.9% slower than slice 1's minimum and outside its two-run spread, so
  `[.assembly_benchmark]` of the solve tests ran four more times each (the deck test's gate-tunnelling cases
  were within 2%). Minimum of six per device: every device within 2% of slice 1 or within slice 1's
  six-run spread. The largest are +2.6% (1D energy balance and 2D drift-diffusion + Canali `evaluate`, 2D
  MOSFET `residual`; spreads 8.1% to 10.8%), 2D drift-diffusion `evaluate` +1.9%, and 3D drift-diffusion
  -3.0%. Passed under GG4's spread rule.
- **GG6, 11 mutants, 10 caught** by the `assemble.*` and `solve.*` ctest entries (control 47/47):
  M1 the frequency-domain interface term active in steady state (also in the steady list): 12 entries (its
  steady `add` is a contract failure); M2 the steady interface term kept in the frequency domain:
  `assemble.small_signal`, `solve.small_signal`; M3 the frequency-domain terms never run: the same two; M5 a
  topic function (`terminal_currents`) reading the first system it was called with, not its argument: 24
  entries; M6 `set_minority_floor` writing a copy of the state: `assemble.drift_diffusion`; M7 `set_paths`
  keeping the old data (positions not moved): `assemble.band_to_band`, `assemble.terms`, `solve.btbt`,
  `solve.rates`; M8 `hold_temperatures` holding in a copy: `assemble.energy_balance`, `solve.energy_balance`;
  M9 `set_thermal_bias` writing a copy: `assemble.electrothermal`; M10 `release_temperatures` a no-op: the two
  energy-balance entries; M11 a model-named member (`tunnelling`) added to `System`: `assemble.fields` (GG7).
  **M4, the frequency-domain terms added before the storage instead of after it, survives** the suites and
  the traces (`assemble.small_signal`, `solve.small_signal`, `solve.transient` byte-identical). It is
  equivalent in every traced use: s is purely imaginary there (i omega t0), so the storage adds only to
  imaginary parts, and the two orders compute the same sums (a + b = b + a; the real parts get no storage
  term). It would differ for a real s, which no caller passes. Not counted; 10 caught meets GG6's 8.
- **GG5: not completed.** The final tree builds with no warnings (`/W4 /WX`, Release). The full Release suite was
  started and stopped by the owner after 5 entries (all passed; "stop the test and push"); the branch is pushed
  for CI, whose Release and Debug jobs run the suites. Before the stop, the assemble.* and solve.* entries
  passed 47/47 (GG6's control on the final sources).

**Merge of `main` (2026-11-16, before slice 3).** PR #85 merged to `main` (9fd9fa1); `main` was merged into this
branch (59051e8). Conflicts: `CMakeLists.txt` (`nitcad_adapt_test` keeps the trace listener and #85's
`adapt_derived_test.cpp`, linked to `NiTCAD::analysis`) and this file (this section after Units 26d and 17c). #85
adds no caller of the removed types. The baselines were recorded again on the merge's binaries (P-core pinned):
80 traced entries, of which 78 are byte-identical to slice 2's; `adapt.adapt_derived` is #85's new entry, and
`assemble.fields` holds slice 2's GG7 case, added after slice 2's trace. The 12 decks agree with slice 0's run A on
all 38 files (`[.deck_compare]`), so #85 changed no GG2 output. Slice 3 is compared against these.

**Slice 3 (GE4) as built (2026-11-16; owner "start Slice 3").**
- `fields.hpp`: `Field::unknown` (bool) became `Field::kind` (`FieldKind::unknown`, `fixed`, `derived`; no field
  is `fixed` yet), and `FieldLocation` gained `edge` and `interface_edge` (one value per trap slot of the
  interface edges, in the assembler's slot order). `FieldId` lists the 21 derived quantities after the six
  unknowns. `DofMap::node_interleaved` requires node fields of kind `unknown`.
- `include/NiTCAD/assemble/derived.hpp` (public) and `src/assemble/derived.cpp`: the table, six rows in
  `DiagnosticId` order. Each row is a diagnostic, its field descriptors (kind `derived`) and an evaluator
  `(const System&, x) -> DerivedValues` (one vector per field; empty when the system lacks the field):
  `bands` (E_c, E_v, E_Fn, E_Fp per node in V_T, from `band_edges`), `rates` (SRH, Auger, radiative, impact,
  interface n/p, band-to-band n/p, net n/p in cm^-3 s^-1, from `rates`, with the sums in `rates`' order as
  `rate_fields` had them), `edge_currents` (Jn, Jp per edge), `edge_flux_terms`, `trap_occupancies` (steady)
  and `temperatures` (T, T_n, T_p in K, T0 + T0 tau with T0 the system's scaling temperature, as the solve
  layer computed them). A row's fields are evaluated together because they come from one pass (one
  `band_edges`, one `rates`). `position(row, field)` and `take(values, row, field)` find a field in a row.
- The solve layer fills through the table (`src/solve/fields.hpp`): `rate_fields`, `band_diagram(system, x,
  V_T)`, `temperatures` and `append_edge_currents`, used by `drift_diffusion_point` and the transient
  snapshots. Each call stays where its predecessor was, so the trace's order of observed values is
  unchanged. **Not through the table:** the equilibrium band diagram (`EquilibriumPoisson`, slice 4's), the
  small-signal temperature amplitudes (a response, not a state's field), and the edge currents and flux
  terms that `weighted_current` and the adjoint indicator compute with (they call `currents.hpp` directly;
  the table is for result fields).
- Tests (`tests/assemble/fields_test.cpp`): `derived fields: the table` (row order, `diagnostic(id)`, every
  field `derived` with its row's location, units, `position`, 21 distinct ids, no unknown among them) and
  `derived fields: the values are the topic functions'` (drift-diffusion, electrothermal, and energy balance
  with impact ionization, on an off-equilibrium state: every row equals its topic function bit for bit,
  the temperatures T0 + T0 tau of their slots, absent fields empty).
- **GG1:** the 80 traced entries of the slice-3 binaries are byte-identical to the merge's, except
  `assemble.fields`, whose trace is the merge's followed by the two new cases.
- **GG2:** the 12 decks agree with the merge's on all 38 files (`[.deck_compare]`).
- **GG4:** FG4's benchmark, the merge's and slice 3's binaries alternating, pinned to CPU 0, two runs each (a first
  attempt ran while a game loaded the machine, 2x slower with spreads to 55%, and was discarded). Minimum of two:
  16 of 17 devices within 2% of the merge on `evaluate` and `residual` or within the merge's spread (3D energy
  balance `evaluate` +17.7% against a 21.4% spread; slice 3 does not change `evaluate` or `residual`). The legacy
  1D diode was +2.0% / +2.2% (11.15 -> 11.37 us, 6.79 -> 6.94 us) against a two-run spread of 0.3% / 0.1%, so its
  case ran four more times each: minimum of six unchanged, the merge's six-run spread 16.1% / 16.7%, and in the four
  extra passes slice 3 was as fast or faster (best `evaluate` 12.45 against 12.64 us). Passed under GG4's spread
  rule.
- **GG6, 11 mutants, 10 caught** by the `assemble.*` and `solve.*` ctest entries (control 47/47): M1 a derived
  field reading the wrong slot (T_n from tau_p's unknowns): `assemble.fields`, `solve.energy_balance`; M2 E_c
  and E_v swapped: `assemble.fields`, `solve.heterojunction`, `solve.electrothermal`; M3 SRH in Auger's
  position: `assemble.fields`, `solve.quantum_confinement`, `solve.rates`; M4 net_n summing the hole rates:
  `solve.rates`; M5 Jn and Jp swapped: `assemble.fields` and 4 solve entries; M6 `position` always the first
  field: 10 entries; M7 `diagnostic(id)` always the first row: 25 entries; M8 the temperatures as
  T0 (1 + tau) (a rounding-level change): `assemble.fields`; M10 the solve fill taking impact ionization from
  the radiative field: `solve.rates`; M11 the hole edge currents not scaled: 4 solve entries. **M9, the
  interface rates' carriers swapped, survives** the suites and the traces (`solve.rates`, whose interface
  case is a MOS with traps and surface recombination, `solve.insulator`, `solve.mos` and `solve.transient`
  byte-identical under it). It is equivalent: the rates are read with steady traps, where each interface
  edge's electron and hole exchange are equal bit for bit, so the two vectors are the same. Not counted;
  10 caught meets GG6's 8. (M4 escapes the new values test because its state has only SRH, Auger and
  radiative, which change both carriers alike; `solve.rates` catches it.)
- **GG5:** full Release suite (`-j 6`, `/W4 /WX`, no warnings) **232/232 pass**: 227 in one run, which stopped when
  the session ended with `verify_full.tvb14` running, and the other five (`verify_full.tvb14`, `app.run_controller`,
  `app_window.analysis_window`, `app_window.self_test`, `app_window_3d`) in a second run on the same build. Debug is
  left to CI (owner's standing rule).

**Slice 4 (GE5) as built (2026-11-16; owner "start slice 4").**
- `System::create(device, scaling, models, Formulation)`: `Formulation` (`drift_diffusion`, the default, or
  `equilibrium`) is a public enum in `system.hpp`; `System::formulation()` reports it. GQ3 (a): `EquilibriumPoisson`
  (`equilibrium_poisson.hpp`, `.cpp`) is removed. Its documentation, `create_equilibrium(device, scaling, models)`
  (the create call with `Formulation::equilibrium`) and its own queries as free functions (`equilibrium_carriers`,
  `charge_neutral_potential`, `contact_potential`, `is_contact`) are in the new public topic header
  `assemble/equilibrium.hpp`. Its `band_edges`, `gate_charges` and `interface_trap_charges` are the topic headers'
  existing functions, which dispatch on the formulation; the equilibrium versions keep the old expressions and add
  no trace observation, as before.
- `src/assemble/equilibrium.cpp`: `SystemData::make_equilibrium` fills only the scaled device's tables (as the old
  `create`; the drift-diffusion model checks do not apply), a stride-1 `DofMap` (psi) and the old pattern (diagonal,
  both directions of each edge). The edges are `EdgeTerm`s with `c`, `charged`, `carriers` and `ab[0]`, `ba[0]`, so
  the separate edge type is gone. `set_bias` keeps the old rule: `check_contact_bias` with thermal equilibrium, and
  only the gates and electrodes follow the bias. The composition table (`Terms::equilibrium_table`), in order:
  equilibrium insulator rows (oxide; assign 0 and a 0 diagonal), equilibrium charge (Poisson; the slaved charge
  under Boltzmann, Fermi-Dirac and incomplete ionization, assigned), **gate** (shared), equilibrium interfaces
  (Fermi occupancy), **Poisson flux** (shared; `Terms::gate()` and `Terms::poisson_flux()` expose the two
  drift-diffusion terms), and the equilibrium contact rows (contacts) **last**. The last term assigns
  psi - psi0 on every ohmic, Schottky and electrode node, clears the row and sets its diagonal to 1. The shared
  flux adds into a Schottky node's row (it skips only ohmic and electrode rows), and this term replaces that row.
  Per entry, the order is charge, gate, interfaces, flux, then the Dirichlet rows. Every value is the old one, which
  GG1 confirms.
- `System`: `set_bias`, the steady `evaluate` and `residual` (the trace labels stay `equilibrium.evaluate.*` and
  `equilibrium.residual.*`), and the Newton hooks dispatch on the formulation. In equilibrium the hooks are
  `newton_solve`'s defaults (psi clipped, the largest |dpsi|), untraced as before. `degeneracies` is empty.
  `state_from_potential`, `stamp_contacts`, `storage`, `small_signal_matrix`, `bias_derivative`, `contact_charges`
  and a time step require the drift-diffusion formulation (NITCAD_EXPECTS).
- Solve: `solve_equilibrium` and the equilibrium-Poisson bias sweep create the system with
  `Formulation::equilibrium`. Their band diagram now goes through the derived-field table (slice 3's exception is
  gone), and the old `band_diagram(BandEdges)` overload is removed.
- Tests: 13 test files migrated (`create_equilibrium` and the free functions). New tests: `terms_test.cpp`
  `equilibrium: the composition of GE5` (sets, stride 1, term order, assignments, models only drift-diffusion reads
  change nothing), `equilibrium: each term's partials match finite differences (GG3)` (Boltzmann, Fermi-Dirac,
  incomplete ionization and both; diode 1D/2D, oxide with traps, MOS-C; every term within 1e-6; the insulator rows
  write only zeros, so they are never active) and `equilibrium: every term writes only the rows and columns it
  declares`; `schottky_test.cpp` `thermal equilibrium's anode row is psi - psi_S alone (GE5)`. Each is appended
  after its file's existing cases.
- Size: `equilibrium.cpp` 409 lines + `equilibrium.hpp` 79 against the old 295 + 153. As R2 expected, the saving is
  structural (one composition, one edge type, one solver interface), not in lines.
- **GG1:** the 80 traced entries of the slice-4 binaries against slice 3's: 78 byte-identical. `assemble.schottky`
  and `assemble.terms` are slice 3's traces exactly, followed by the new cases. The first trace pass's `deck.run`
  failed only because the mutant builds were relinking `build/release/nitcad.exe`, which that test runs; rerun alone,
  it passed and is identical.
- **GG2:** the 12 decks agree with slice 3's on all 38 files (`[.deck_compare]`).
- **GG3:** every existing FD gate passes. The equilibrium charge term's own FD check is above.
- **GG4:** FG4's benchmark, slice 3 and slice 4 binaries alternating, pinned to CPU 0. Minimum of two: 15 of 17
  devices within 2% or slice 3's spread (most faster). Two were flagged: 2D energy balance `evaluate` +33.9%
  (5.80 -> 7.77 ms, spread 2.4%; `residual` unchanged) and 1D Zener `evaluate` +6.4% (spread 5.9%). Two more
  passes each (minimum of four): **all 17 within GG4's rule**. Zener -0.6%; 2D energy balance `evaluate` +6.0%
  (4.95 -> 5.24 ms) against slice 3's four-run spread of 20.0%, `residual` +1.0%. Slice 4 changes nothing on that
  path but two formulation tests. The machine was noisy in the extra passes (several devices' spreads above 100%), so
  this passes on the spread rule, not on a quiet-machine margin.
- **GG6, 13 mutants, all caught** by the `assemble.*` and `solve.*` ctest entries: M1 the slaved charge's Fermi-Dirac
  branch using Boltzmann: `assemble.insulator`, `solve.quantum_confinement`, `solve.fermi_dirac`,
  `solve.heterojunction`; M2 equilibrium missing the interface set (its activation false; removing the row fails
  `/WX` on the unused lambda): `assemble.insulator`, `assemble.terms`, `solve.quantum_confinement`,
  `solve.insulator`; M3 the Dirichlet rows before the edge flux: `assemble.terms`, `assemble.schottky`; M4 the
  electrodes keeping the old bias after `set_bias`: `assemble.insulator`, `solve.quantum_confinement`,
  `solve.insulator`; M5 gate and charge swapped: 6 entries; M6 the interface edges' off-diagonal positions swapped:
  `solve.small_signal`; M7 the equilibrium Newton update not clipped: `solve.btbt`; M8 incomplete ionization ignored
  by the slaved charge: `solve.heterojunction`; M9 an electrode's charge without its interface half-edge flux:
  `solve.insulator`; M10 the Dirichlet rows not cleared: `assemble.schottky`; M11 `band_edges` not dispatched:
  `assemble.heterojunction`, `solve.heterojunction`; M12 `set_bias` without the thermal-equilibrium check:
  `assemble.gate`; M13 the interface edges at another node's statistics: `solve.quantum_confinement`.
- **GG5:** Release, `/W4 /WX`, ctest `-j 6` on the slice 4 tree (2026-11-17): **232/232 pass** (2,982 s). Debug is
  left to CI (owner's standing rule).
- **GG8 (owner 2026-11-16, "choose recommended for GG8"):** the χ extensibility test is its own commit after slice 4,
  not part of it, so slice 4 stays GE5 alone (GE7).

**GG8 as built (2026-11-18; owner "start GG8").**
- **The extension point** (internal, `src/assemble`; `System`'s header is unchanged). `System::create`'s
  drift-diffusion body moved into `SystemData::make_drift_diffusion(device, scaling, models, extension)`, the
  counterpart of `make_equilibrium`; `create` calls it with none and sets the zero bias as before. A
  `SystemData::CompositionExtension` (composition.hpp) holds fields and terms: `Composition::compose` puts its fields
  after the built-in ones in every node block and its terms after every term of the table, under no equation set.
  `SystemData::create_extended` (isothermal drift diffusion only) builds such a system: `System`'s constructor is its
  own, so a system created without the extension is the shell whose data and state are replaced through the passkey
  accessors. An extension's terms write their own rows, add their own pattern entries (`Term::pattern`) and find
  their positions at each add; steady evaluations only.
- **What the test showed had to change.** (i) `SystemData::electrothermal()` was `m_ == 4`, so a fourth isothermal
  field would have been taken for the heat equation; it is now a flag set from the models (`electrothermal_`).
  (ii) `System::update_size` and `apply_update` moved only ψ, n, p and the temperatures, so an added unknown would
  never have moved; the unknowns after `built_in_` are now clipped to ±max_update and measured by |dx|, as ψ. Both
  loops are skipped when there is no extension, so every existing system takes the same path as before.
- **The test** (`tests/assemble/extensibility_test.cpp`, 5 cases, `[extensibility]`). χ has its own `FieldId`
  value outside the enumerators (fields.hpp is not edited). Its row is the box-method form of d²χ/du² = ψ − ψ_ref − S
  on u = x / L (1D), held at 0 on insulator and contact nodes; the second term adds k V χ to the electron rows and
  −k V χ to the hole rows of the semiconductor nodes. (i) Block 3 → 4, the pattern gains exactly 7 entries per node
  and 2 per edge, every plain entry kept at its place, the table's terms then the two. (ii) Both terms and every
  existing term at the wider block pass the FD gate (≤ 1e-6, terms_test's method) on the 1D MOS off equilibrium; the
  cross blocks (n, p rows in χ; χ row in ψ) are present with their values. (iii) Newton's last pair above rounding
  has order ≥ 1.8 (measured 2.12 at k = 0, 2.07 at k = 1e-11). (iv) At k = 0, ψ, n, p equal the plain system's
  within 1e-12. (v) Manufactured χ* = 2 sin πu on a 1D n-bar at equilibrium, its source in χ's row and the
  generation hook (TVQ3) cancelling k V χ*: χ's and the holes' L∞ errors fall at order 2.00 over 32, 64, 128
  elements (χ 1.6e-3, 4.0e-4, 1.0e-4; p relative 9.5e-2, 2.4e-2, 5.9e-3). A fifth case starts an intrinsic bar at
  its exact equilibrium with ψ_ref = ψ_c − 200, so only χ moves, over several clipped steps.
- **k is in scaled rate units** (R0 ≈ 1e33 cm⁻³ s⁻¹ on the bar). At k = 1e-3 the minority holes moved by a factor of
  about 1e7 and the MOS Newton run failed, so the test uses k = 1e-11. ψ's and n's errors on the bar are then at
  rounding (6e-14, 9e-14) and are not gated. At k = 1e-3 they too fell at order 2.0.
- **Mutants: 9, all caught** by `nitcad_assemble_test`. M1 extension unknowns not added (`apply_update`): Newton and
  manufactured cases. M2 not measured (`update_size`): the fifth case only. On a doped device M2 is equivalent in
  practice: the minority densities' corrections carry the rounding of a large χ correction (measured on the MOS,
  update_size stayed near 1e-4 while χ moved), so Newton did not stop early; hence the intrinsic bar. M3
  `electrothermal()` back to `m_ == 4`: crash, FD case. M4 `built_in_ = m_`: crash, manufactured case. M5 extension
  terms dropped: four GG8 cases. M6 extension fields before the built-in ones: pattern case. M7 no zero bias in
  `create_extended`: Newton and manufactured cases. M8 `make_drift_diffusion` ignoring the extension: crash. M9
  `electrothermal_` never set: the electrothermal composition case, then a crash.
- **Toolset (2026-11-18).** A Visual Studio 18.12 update had replaced MSVC 14.51 with 14.52.37025. **14.52 drops
  aggregate assignments at /O1 and /O2**: `none.impact_ionization = {};` in `impact_ionization_test.cpp` left
  silicon's coefficients in place, so `assemble.impact_ionization` failed on slice 4's tree as well as on this one.
  A 20-line reproduction (a copy of a constexpr aggregate, `= {}` on a member, the type also a by-value default
  argument) prints the old values at /O1 and /O2 and zeros at /Od; 14.51.36231 prints zeros at every level. The
  owner reinstalled 14.51 and every gate below is on 14.51, as CI and slices 0–4. Q11's move to 14.52 waits for a
  fixed compiler.
- **GG1:** the 81 traced entries of this tree against slice 4's (a9cfd64, rebuilt with 14.51 in a scratch worktree),
  P-core pinned, six entries at a time. Catch2 v3 runs the cases of a file in random order, so the traces are
  compared case block by case block. **All 81 identical.** In `solve.mos` and `deck_lombardi` some cases share a
  cached solve, whose values the trace files under whichever case runs first; there the values are compared as one
  set per file (`solve.mos`: 45,050 values, `deck_lombardi`: 6,421, the same on both trees). The scratch baseline was
  first configured without `NITCAD_MKL_RUNTIME` and `NITCAD_VERIFY_FULL_TIER`, so its PARDISO cases were skipped and
  five entries seemed to differ; reconfigured as `build/release`, those five are identical too.
- **GG2:** the 12 decks, run pinned with each tree's `nitcad.exe`, agree on all 38 output files (`[.deck_compare]`, 125
  assertions).
- **GG4:** FG4's benchmark (17 devices), slice 4's and this tree's binaries alternating, pinned to CPU 0. Four
  passes: no device flagged, but only on the spread rule, as slice 4's own spread was 64–414% (passes 1 and 2 ran
  on a loaded machine), and five devices' minima were 10–12% slower (1D drift-diffusion, with Canali, 2D
  drift-diffusion `evaluate`, 2D energy balance, 3D generated diode), so the owner asked for a focused run. Those
  are in two cases (EC19's and the generated 3D mesh's, 10 devices): six more rounds of both, alternated, the order
  swapped each round. **Minimum of six: this tree faster on all 10 devices**, `evaluate` −6.4% to −13.8%,
  `residual` −4.7% to −13.0%; medians faster on 9 (the tenth, 3D drift-diffusion `evaluate`, +12.4% with its
  minimum −10.6%). Single rounds swung from −40% to +70% on one device, so the machine's noise even pinned is
  ±15–30% per run; nothing in this commit is on the `evaluate` or `residual` path of an existing system.
- **GG5:** Release, `/W4 /WX` (no warnings), MSVC 14.51, on the GG8 tree (2026-11-18): **233/233 pass**. The 226
  regular tests with `-j 6` in 919 s, then the 7 full-tier tests serially (they are memory-heavy: TVB-04 reaches
  about 20 GB) in 4,040 s (tvb04 1,356 s, tvb14 968 s, tvb20 671 s, tvb05 578 s, tvb19 466 s). A first attempt
  running the full tier inside `-j 6` was stopped after 11 passing tests. Debug is left to CI (owner's standing rule).

**Risks:** the bit-identity surface is the largest since Unit 40. Moving the bias and path state out of the system
could change which values a copy sees (GG6's first three mutants target this). Slice 2's caller migration touches
about 44 files (more if #85 adds callers). Equilibrium's variants (GE5) keep a second set of row expressions, so the
saving there is structural more than in lines. Trace and deck runs need P-core pinning (Unit 41).
