<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 45 — generation and recombination rates and the mesh digest in results (proposed 2026-10-21; on `results/rates-digest`)

**Owner request:** VQ12 (2026-10-17) put this unit before 25b: "recombination and generation rates in results and
files; the mesh digest in sweep, transient and small-signal files", to unlock the R view and exact result–mesh
pairing. Documentation only until the owner approves.

**Why now (audited on `main` at `6d2135b`).**
- No result holds a rate. The assembler computes each mechanism only inside its residual: since Unit 40 the
  drift-diffusion rows are equation terms, one per mechanism (`terms.hpp`: SRH, Auger, radiative, interfaces,
  impact generation, local and nonlocal band-to-band tunnelling, and their electrothermal forms). The only rate in
  results is a tunnel path's (`TunnelPath::generation_cm3_s`).
- `nitcad plot` pairs a result file with a device by node count (25a·3): two meshes of equal size pass. The mesh
  digest exists (`mesh::graph_digest`, `io::mesh_digest`) and only checkpoints and P4 certificates store it.
- 25b's derived-field views and the reference's R panel need the rates per node.

**Design (proposed):**

| # | Choice | Why | Trade-off |
|---|---|---|---|
| RS1 | **Rates read from the terms.** `DriftDiffusion::rates(x)` evaluates each recombination or generation term alone on the state and takes its contribution to each node's electron and hole continuity rows, divided by the node's control volume and scaled to cm⁻³ s⁻¹. Per node: SRH, Auger, radiative and surface recombination (net, R − G of each, signed), impact-ionization and band-to-band generation (electron and hole rows separately where they differ: nonlocal tunnelling puts its electrons and holes at the path's two ends). | The stored rate is the residual's, operand for operand: no second implementation to drift from the solver. | One extra pass of those terms per stored point. |
| RS2 | **`results::RateFields`** in `BiasPoint`, `TransientSnapshot` and the small-signal operating points: the components above and the net `U = R − G` per carrier, each empty when its mechanism is off. Surface recombination is the interface edges' rate over the node's control volume (an effective volume rate at interface nodes), labelled as such. | One type for every result that has fields; empty vectors cost nothing. | A volume rate for a surface process is an approximation of its distribution, not of its total (the total is exact). |
| RS3 | **The mesh digest in `RunRecord`** (`mesh_digest`, filled by `make_run_record` from the device's mesh); not part of `input_identity` (which already covers the mesh). | Every result type carries its run record, so sweep, transient and small-signal files get it at once. | — |
| RS4 | **Files:** new optional records in the sweep, transient and small-signal files for the rates and the digest; files without them (every file written before) load as before, with no rates and digest 0 (unknown). | Old files keep loading (Unit 42's rule, PG3). | — |
| RS5 | **`nitcad plot`** checks the digest when the result has one (exact refusal), and falls back to the node count for older files; `render` gains the fields `R` (net), `srh`, `auger`, `radiative`, `surface`, `Gii`, `Gbtbt`. | The pairing VR6 asked for; the R view for 25b and the reference grid. | Seven field keys. |
| RS6 | **Out of scope:** |E| and |J| (25b, `render` reconstructs them); optical generation (no model yet); rates in checkpoints; per-edge rates. | Keeps the unit to stored rates and the digest. | — |

**Owner decisions (asked):**

| # | Question | Recommendation |
|---|---|---|
| RQ1 | Store the rates by default, or only when asked (`BiasOptions::store_rates`, deck `[solver] rates = on`)? | **By default**, with `rates = off` to save file size: a commercial viewer opens any result and shows R without rerunning. Measured cost in RG5 decides whether the default holds. |
| RQ2 | Components (as RS1), or the net rate only? | **Components and the net rate**: which mechanism dominates is the question users ask of an R map. |
| RQ3 | Surface recombination as an effective volume rate at interface nodes (RS2), or left out? | **Include, labelled** as an effective volume rate, with its exact total. |
| RQ4 | The digest refusal in `nitcad plot`: refuse a result whose digest differs, or warn? | **Refuse** (exit 1), as for a node-count mismatch: an image of the wrong mesh is worse than none. |

**Acceptance gates, fixed before any implementation:**

| # | Gate |
|---|---|
| RG1 | **Exact against the residual:** for every mechanism and every state of the existing device tests' kinds (diode, MOSFET, traps, impact, local and nonlocal BTBT, electrothermal), the volume-weighted rates reproduce that term's contribution to the continuity residual to rounding. |
| RG2 | **Physics checks:** a uniform bar at low injection gives SRH = Δn/τ; radiative B(np − n_ie²) exactly; Auger as its formula; at equilibrium every net rate is zero to rounding; the integral of impact generation equals the avalanche current's increase (Unit 19's integral). |
| RG3 | **No numerical change:** the assembly traces are byte-identical to `main`'s (`--order decl`); every suite passes unchanged; `input_identity` of the fixtures unchanged. |
| RG4 | **Files:** round trips of rates and digest bit for bit; every committed file of earlier units (`tests/io/fixtures`) loads unchanged with no rates and digest 0. |
| RG5 | **Cost measured:** time and file size per point with and without rates on the 2D MOSFET sweep and a 3D diode; recorded whatever they show. |
| RG6 | **`nitcad plot`:** a result of a different mesh with the same node count is refused by digest; an old file falls back to the node count; the R and component fields draw. |
| RG7 | **Suites and mutants:** Release and Debug at every commit, no warnings; mutants in the rate extraction (a term's sign, a missing volume division, electron and hole rows swapped, a mechanism's flag), the file records and the digest check, all caught. |
| RG8 | **No scope creep:** nothing outside RS1–RS5; ARCHITECTURE.md 6.6 and 11.3 as built after acceptance. |

**Commits (proposed):** 1) `rates()` and `RateFields`, filled by the solves (RG1–RG3, RG5's timing); 2) the digest in
`RunRecord` and the file records (RG4, RG5's sizes); 3) `nitcad plot` and the render fields (RG6).

**Risks:** the rates' extra pass slows stored points (RG5 measures it); files grow by up to eight node arrays per
point, which RQ1's off switch answers; the electrothermal terms add their rates into the heat row too (only the
continuity rows are read).

**Status:** authorized by the owner 2026-10-21 ("Approve, as recommended"): RQ1 rates stored by default
(`[solver] rates = off` to opt out), RQ2 components and the net rate, RQ3 surface recombination included and
labelled as an effective volume rate, RQ4 `nitcad plot` refuses a differing digest (exit 1). Commits 1–3 done,
**accepted by the owner (2026-10-21)**; ARCHITECTURE.md 6.6 and 11.3 as built (RG8).

**Commit 1 (the rates; RG1–RG3, RG5's timing), as built.** `DriftDiffusion::rates(x)` evaluates each generation-
recombination term alone (residual only, no Jacobian) and reads its electron and hole rows: rate = row × R0 / V in
cm⁻³ s⁻¹, positive for net recombination (one conversion for every term: `rate_scale_` = Ns / R0 and
`generation_scale_` = 1 / R0). Mechanisms by term name; local and nonlocal tunnelling add into one entry; the
electrothermal terms are read the same way. `results::RateFields` (srh, auger, radiative, impact; interface_n/p,
band_to_band_n/p; net_n/p) in `BiasPoint` and `TransientSnapshot` (small-signal operating points are `BiasPoint`s),
filled by `sweep_bias`, `trace_bias` and `solve_transient` after convergence; `BiasOptions::store_rates` (on, RQ1),
recorded in the run record only when off; deck `[solver] rates = on|off`. `rates()` reports nothing to the
assembly observer.
- **RG1, the balance:** at every node off the contacts of a converged state, the net current out of its control
  volume (from the stored edge currents) equals q U V for electrons and −q U V for holes, within 10⁻⁶ of the node's
  currents plus 5% of the run's current resolution (a node's net current below it is the rounding of large opposing
  majority fluxes: measured ≤ 9e-15 A/cm against a resolution of 1.06e-12 in a reverse-biased diode, while the
  minority carriers agree to five digits). Holds for SRH, Auger and radiative (forward diode), impact ionization
  (reverse to −4 V), interface traps and surface recombination (the small-signal tests' oxide-silicon device), local
  and nonlocal band-to-band tunnelling (nonlocal: electron and hole rates differ, each balanced) and the
  electrothermal terms.
- **RG2:** at equilibrium every rate is below 1e4 cm⁻³ s⁻¹ (against n_i/τ ≈ 1e16); forward recombination exceeds
  1e18; impact and tunnelling rates are generations (negative).
- **RG3:** assembly traces of the bias, sweep and transient tests byte-identical to `main`'s (`6d2135b`, `--order
  decl`, 852,646 lines); a run storing rates keeps the default run identity; switching them off changes no field or
  current.
- **RG5 (time):** sweeps of 7 points with and without rates, best of three: 2D diode, 14,701 nodes: 8.92 s and
  9.19 s; 3D diode, 4,961 nodes: 20.17 s and 20.24 s. The extra pass is below the runs' noise.
- **Mutation checks (RG7, Release): 10 of 10 caught** — the electron row's sign, the volume division, the rows
  swapped, Auger read as SRH, tunnelling dropped, the hole net from the electron rows (caught only once the nonlocal
  case was added: every other mechanism moves carriers in pairs), `store_rates` ignored in sweeps and in
  snapshots, the setting always recorded, the deck key ignored. (A first runner lost a test's exit status between
  two batch calls; its two survivors were rerun with a fixed runner, one then caught.) Release 22/22 and Debug 16/16,
  no warnings.

**Commit 2 (the digest and the files; RG4, RG5's sizes), as built.** `RunRecord::mesh_digest` (the device mesh's
`mesh::graph_digest`, filled by every `make_run_record`; not in `input_identity`; 0 = unknown). Files: the run record's
key 4 (the digest, written when nonzero), the points' key 18 and the snapshots' key 8 (a rates record of the ten
vectors, written when the rates were stored); absent in older files, which load with no rates and digest 0.
- **RG4:** sweep, transient and small-signal files round-trip the rates (signed zero, NaN, a denormal) and the digest
  bit for bit; a point or snapshot without rates stays without; the pre-Unit-42 fixtures load with no rates and
  digest 0; every run record (sweep, transient, by hand) carries the device mesh's digest.
- **RG5 (size), the sweeps above:** 2D diode 13.1 MB with rates against 9.0 MB without (+46%); 3D diode 4.9 against
  3.5 MB (+40%): five node arrays per point (SRH, Auger, radiative and the two nets) beside the fields, bands and
  currents. Time stays within noise (+0.1% and +4.5% this run). **For the owner (RQ1):** the default stays on as
  ruled; `rates = off` saves the 40–46%.
- **Mutation checks (RG7, Release): 7 of 7 caught** — the digest not written, not read, not filled; the point's and
  the snapshot's rates not written; the point's not read; a rates key swapped. Release 22/22 and Debug 16/16, no
  warnings.

**Commit 3 (`nitcad plot` and the render fields; RG6), as built.** `render::Field` gains `R` (net R − G),
`srh`, `auger`, `radiative`, `surface` (signed, on the signed log scale), `Gii` and `Gbtbt` (generations, positive,
log10); `field_values` refuses a rate the result does not hold, saying how to store it. `nitcad plot` refuses a
result whose run record names another mesh (RQ4: "its mesh digest is not the device's", exit 1), and falls back to
the node count for a result without a digest.
- **RG6:** a result of a mesh with the same node count (x over 3 µm instead of 2) refused by digest; the R map of a
  forward-biased diode drawn (looked at: the peak at the junction, zero at the contacts); a pre-Unit-42 sweep with
  its device drawn through the node-count fallback, its R refused ("holds no R rates"), and with another device
  refused by node count; the render fields' values, signs, scales and the missing-rate error.
- **Mutation checks (RG7, Release): 4 of 4 caught** — the digest check off, a generation not negated, SRH read as
  Auger, a missing rate not refused. Unit 45 in all: 21 of 21. Release 22/22 and Debug 16/16, no warnings.
