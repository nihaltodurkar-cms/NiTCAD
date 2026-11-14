<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 31 — device files, the text input deck and the batch runner (merged to `main` as `6f67510`, PR #37, 2026-10-15)

**Owner rulings (2026-10-15)**, asked before any code:
- syntax: an own sectioned `key = value` format with an in-tree parser; no TOML or JSON dependency;
- scope: device descriptions on disk, the deck turned into a device, options and a run plan, and a `nitcad`
  runner that writes sweeps and checkpoints;
- coverage: tensor grids, the 2D and 3D generators and mesh files; equilibrium and steady bias sweeps (both
  equation sets), step cutting, line search, the physics switches, electrothermal. Transient, small-signal and
  `trace_bias` are a later unit, together with their results on disk;
- units: plain numbers are cm; `um` and `nm` suffixes convert on read.

Choices made within those rulings, open to review:

| # | Choice | Why | Trade-off |
|---|---|---|---|
| DK1 | The grammar (`deck/syntax.hpp`): `[kind name]` headers, `key = value` entries, `#` comments outside double quotes, a trailing `\` continues a line; keys and names are identifiers joined by single blanks. | Line-oriented, so every error has a line; trivially diffable; no escaping rules beyond quotes. | No nesting, variables or includes; a long polygon needs continuations. |
| DK2 | A new layer `deck` above `solve` and `io` (it reads the option types of `linalg` and `assemble` through `solve`); `io` gains `device` (device files). `solve` still does not use `io`. | The deck must produce `solve::BiasOptions` and read mesh files, so it sits above both; device files belong with the other files. | One more layer in D5 (its name is open to review); `io` now depends on `physics` through `device`. |
| DK3 | Strict: an unknown section or key is an error, as are a second single-instance section and a deck with not exactly one of `[equilibrium]` and `[sweep]`. | A typo (`max_iteration`) silently ignored would run the wrong problem. | Decks written for a later build fail on an earlier one with a clear message. |
| DK4 | Syntax and meaning are separate passes. Every error is `invalid_input` "`<file>:<line>: <what>`" with the line as the context index; an error from a library (a generator, a material, `Device::create`) keeps its code and context and gets the location and section prefixed. | The user sees where, and the library's own words for what. | `Device::create`'s errors have no deck line (they concern the device as a whole) and report line 1. |
| DK5 | Only lengths take suffixes (`cm`, `um`, `nm`); every other quantity is in the unit of the C++ field it sets (V, K, cm^-3, eV, cm^-2, cm/s, s). | Owner ruling for lengths; one unit per quantity elsewhere keeps the deck and the API in step. | `1e17` must be cm^-3, not m^-3. |
| DK6 | Series: `first : last : step` with `last` included when it lies on the grid to 1e-9 of a step (and then stored as `last` exactly), a comma list, or one value held at every point; a contact not named is held at 0 V; every series longer than one value has the same length. | The common sweeps in one line; no hidden accumulation of rounding (each value is `first + k step`). | No logarithmic or adaptive series. |
| DK7 | Material parameters are overridden by their C++ member paths (`electron_mobility.mu_max`, `thermal.conductivity_W_cmK`); materials start from the built-in sets (`silicon`, `germanium`, `gaas`, `ingaas`, `sic_4h`, `algaas <x> [<share>]`, `sio2`, `insulator`). | One name per parameter, the one the code and its comments use; nothing to keep in sync. | Renaming a parameter in C++ breaks decks that set it. |
| DK8 | Regions on a tensor grid: inclusive boxes; the one region without a box is the background; later regions override earlier ones; every node must land in a region. Generated meshes take polygons (2D) or seeds (3D) per region; region *k* is the *k*-th `[region]` section, also for mesh files. A contact's nodes are its patch's semiconductor nodes (ohmic, gate) or insulator nodes (electrode), optionally those inside any of its boxes. Tensor meshes always carry their cells. | Covers the hand-built devices of the test suite without coordinates per node. | Nodes exactly on a box face belong to it; the deck author places boxes between nodes (or uses `straddle`). |
| DK9 | A restart matches the checkpoint's bias to a sweep point to 1e-12 relative (1e-12 V near zero), sweeps the points after it from the checkpoint's state with `initial_bias_V` set, and writes those points only. | Found in testing: `0 : 0.3 : 0.1` ends at 0.3 exactly, the fourth value of `0 : 0.6 : 0.1` is 0.30000000000000004; a bit-for-bit match refused a correct restart. | A sweep file of a restarted run holds the continuation, not the whole sweep. |
| DK10 | The runner: `nitcad check|run <deck>` in `src/cli/main.cpp` (not `app/`, which the layout reserves for the Win32 GUI); exit status 0 done, 1 error, 2 the sweep stopped early (outputs written), 64 usage; Ctrl+C asks the solver to stop. Outputs are written after the analysis, the checkpoint being the last completed point. | The batch runner 11.3 asks for before the GUI. `src/cli/` is a directory the layout did not list: open to review. | A crash (not a stop) loses the run; a checkpoint per point would need a per-point hook in `solve`. |
| DK11 | Device files: content kind 4 in the Unit 29 container; every material parameter under its own key (absent keys take the struct's defaults, IO2); tensor cells as their axes; the P4 certificate as issued; a loaded description goes through `Device::create`. | Bit-for-bit round trips; old files keep loading as parameters are added; no file can bypass validation. | The parameter keys are a format: never renumbered. |

Not in Unit 31: transient, small-signal and `trace_bias` decks (owner ruling); 2D structured boundary layers and
per-region mesh sizes in the deck; general 1D graph meshes (they need certificates); variables, expressions,
includes and parameter studies.

