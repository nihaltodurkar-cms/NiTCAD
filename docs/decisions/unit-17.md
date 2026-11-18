<!-- Archived from DECISIONS.md on 2026-11-17 (branch docs/archive-17-26d), word for word; DECISIONS.md keeps its summary. -->

### Unit 17 — intelligent adaptive mesh refinement and state transfer (proposed 2026-10-31; approved with amendments 2026-10-31; on `mesh/adaptive-refinement`)

**Owner request:** "I want Intelligent Adaptive mesh refinement" (2026-10-31). ARCHITECTURE.md 6.12.4 and the unit
table's Unit 17 ("Adaptive mesh refinement and state transfer", target). The branch started from `main` `95c2254`; rebuilt on `main` `03e0bf7` (PR #80) when implementation was authorized
(Unit 54 merged). The standing conditions hold: no change to any existing mesh, device, result or run identity when
adaptation is off; criteria fixed before code; Release and Debug suites at every commit, `/W4 /WX`; hand mutants;
nothing pushed or merged without authorization.

**Owner rulings (2026-10-31):** "Approve the proposed AQ1–AQ5 and LQ1–LQ4 recommendations. Proceed with the Unit 17 /
Unit 17b split." Amendments, incorporated below before the proposal was frozen:

| # | Amendment (owner) | Where |
|---|---|---|
| OA1 | "Treat greedy benefit-per-cost allocation as a heuristic, not a globally optimal solution." | ID6 (no optimality claim), AC16 |
| OA2 | "Require independent numerical error verification, conventional AMR fallback and explicit accuracy certification." | AD11–AD13, ID1, AC17, AC20 |
| OA3 | "Add false-positive refinement/mesh overhead measurements alongside the 80% regional recall target." | AC14 |
| OA4 | "Estimate computational cost using solver timings, DOFs and Jacobian sparsity—not node count alone." | ID4, ID5, ID6, AC15 |
| OA5 | "Keep cross-run ML opt-in, with reproducible training data, model versioning and held-out device validation." | ID8, ID10, AC19 |
| OA6 | "Keep numerical accuracy independent of ML predictions." | The governing rule below; AD11–AD13; AC21 |

"Do not write implementation code yet." The design and criteria below are frozen for review; implementation starts
only on the owner's word.

**Implementation authorized (owner, 2026-11-12):** "Implement Unit 17 adaptive refinement." The design and criteria
below are implemented as frozen (no criterion changed); the slices follow in order, each a local commit. The three
proposal commits (`d363963`, `9d879be`, `a534bda`, on the old base `95c2254`) were folded into one commit on
`main` `03e0bf7`; the text is unchanged apart from this paragraph and the base named above.

**Owner definition (2026-10-31):** "Intelligent Adaptive Mesh Refinement (Intelligent AMR): adds a predictive or
learning-based decision layer that can identify where refinement is likely to be needed, estimate its benefit, and
optimize mesh allocation." The unit therefore has two layers: the **adaptive base** (items 1–5 below: indicators,
size fields, transfer, the loop; a conventional AMR) and the **intelligent decision layer** on top of it (ID1–ID9,
below the base design): prediction, benefit estimation, allocation, and learning. **The governing rule (ID1, OA6):**
the decision layer only decides *where to spend nodes*. Numerical accuracy is independent of every prediction and
every learned model: the error estimate (AD3's adjoint), the independent verification (AD11), the stopping rule (AD6)
and the accuracy certificate (AD13) are computed from solves alone and never read a prediction, a fitted parameter of
the decision layer or a `.nitlearn` model. Every result is a converged solve on a certified mesh. A poor prediction
costs time, never accuracy; repeated poor predictions hand the run to conventional AMR (AD12).

**The adaptive base (the decision layer's foundation; each part gated below).** Today the user chooses every size: h,
per-region h, boxes, s, layers, tensor lines. The adaptive loop chooses them from the solution:
1. **Goal-driven.** The user states what must be accurate (a terminal current, or V_T / g_m / swing read from the
   sweep) and a tolerance, not mesh sizes. The loop stops when the goal's change between passes, and its estimate,
   fall below the tolerance, or at a node or pass budget, which it names.
2. **Physics-aware indicators.** Where the discretization is known to lose accuracy: the potential step across an
   edge in units of V_T (Scharfetter–Gummel's accuracy), the change of log density along an edge (inversion and
   depletion edges), doping steps (junctions, before any solve), and where active models concentrate (impact
   ionization and BTBT generation, the field at a tunnelling insulator).
3. **Goal-weighted (adjoint).** An indicator only counts where the goal is sensitive to it: each node's indicator is
   weighted by the adjoint z of the goal, J^T z = ∂I/∂u (one transposed solve per goal per adaptation point). A
   junction far from the current path (a well, the substrate) is then not refined for a drain-current goal.
4. **Bias-aware.** The user lists the bias points to adapt at (or `auto`: the sweep's ends and the point of largest
   goal curvature from a first coarse sweep); marks from every point are united, so one mesh serves the whole sweep.
5. **Explains itself.** A report per pass: nodes, the goal and its estimate, which indicator marked how many nodes
   and where (regions), the stopping reason; the final mesh saved and viewable like any generated mesh.

**The questions that shape the design (owner rulings 2026-10-31: every recommendation approved):**

| # | Question | Options | Recommendation |
|---|---|---|---|
| AQ1 | How the mesh changes | (a) **Regenerate** with a size field: the existing generators already size every element through one `size_at(point)` (`generate2d.cpp:345`, `generate3d.cpp:554`, `generate_complex.cpp:785`); the indicators give a nodal h on the old mesh, graded, and the generator meets it. Every certificate (P0–P5, P4, P6) comes from the existing checked path. (b) **Local insertion** into the existing triangulation (Ruppert/Bowyer–Watson points where marked). Cheaper per pass, but P4's interface bands and the layers' exemptions must be re-proven for incremental insertion. | **(a)**: no new geometry, certificates unchanged. Cost is a full mesh per pass (2D: seconds; swept: the 2D cost; general 3D: about 130–160 s per 1e5 nodes, FC12). |
| AQ2 | Mesh kinds in this unit | (a) **Tensor** (insert lines where marked, per axis; nonlocal BTBT stays valid), **2D generated**, **swept 3D** (adapt the cross-section's size field and the axis lines). (b) Also general 3D (`generate3d`, Unit 53's complex). | **(a)**. All three are the 2D generator plus line sets; general 3D re-mesh cost is high and goes to a later unit. Structured layers (U14, UD5) are kept as given: never refined inside; an indicator saturating inside a layer is reported. |
| AQ3 | Mesh per sweep or per bias point | (a) **One mesh per sweep**: adapt at the listed points, refine-only union of marks, then the whole sweep on the final mesh. (b) A mesh per point (state transferred point to point). | **(a)**: one mesh, one run identity, curves comparable point to point; sweeps, files, plots and the app unchanged. (b) changes every downstream layer. |
| AQ4 | Goal weighting (adjoint) | (a) **In this unit**: a transposed solve in `linalg` (SparseLU and PARDISO; GMRES-ILUT falls back to SparseLU for this solve, reported) and ∂I/∂u from the weighted current (Unit 32's W1, which is linear in the state's fluxes). (b) Feature indicators only; adjoint later. | **(a)**: it is what makes the loop goal-driven rather than "refine everything steep". Gated by AC7 against feature-only. |
| AQ5 | Coarsening and anisotropy | Coarsening below the initial mesh's sizes; anisotropic (metric) refinement. | **Both out of this unit.** The initial mesh is the coarsest; anisotropy is already supplied by layers and tensor lines, which the loop adapts along their own axes (tensor lines, swept axis lines). |
| LQ1 | What "learning" learns from | (a) **Within a run** (online): the decision layer fits its own models to the run's passes and bias points (local convergence orders, error constants, cost per node) and uses them for the next decision. (b) **Across runs** (persistent): a model trained on earlier runs' adaptation histories (a versioned, checksummed `.nitlearn` file, its digest in the run identity) predicts the first mesh of a new run. (c) **Shipped pre-trained** model: needs a training corpus that does not exist yet. | **(a) and (b) in this unit, (b) opt-in**; (c) not now. (a) needs no stored state and is fully reproducible from the run. (b) is where "likely to be needed before solving" pays off most (fewer passes from the first mesh). |
| LQ2 | Model class for the cross-run predictor | (a) **Gradient-boosted regression trees**, in-tree C++ (no ML library; no Python, owner rule), deterministic training. (b) Ridge regression on engineered features. (c) A small neural network (MLP), in-tree C++. | **(a)**: handles mixed features (material, region kind, log doping, distances, bias) and their thresholds without scaling, deterministic, small (hundreds of trees × depth 4), inspectable (feature importances in the report). (b) is the baseline it must beat (ID9). (c) needs far more data than a user's run history gives. |
| LQ3 | The allocation objective | (a) **Tolerance mode**: reach the goal tolerance with the fewest nodes. (b) **Budget mode**: the smallest estimated error for a node (or wall-time) budget. | **Both**, the user picks; tolerance mode the default. |
| LQ4 | One unit or two | (a) One unit, all of it (about 9 slices, XXL). (b) **Unit 17: base + online decision layer (prediction from physics and from the sweep, benefit, allocation; LQ1 a); Unit 17b: the cross-run learned predictor (LQ1 b, LQ2).** | **(b)**: 17b needs 17's adaptation histories as its training data, and 17's gates give the baseline 17b must beat. The design of both is fixed now. |

**Design (as approved, with amendments OA1–OA6):**

| # | Choice | Evidence | Why | Trade-off |
|---|---|---|---|---|
| AD1 | **The loop** (`solve::adapt`, new): solve (the existing analysis) → indicators → mark → size field → regenerate → transfer → solve; the analysis and the device builder are passed in, so physics, solver and analysis code do not change. | 6.12.4's target loop. | Device-neutral (6.12.8). | One more layer of orchestration in `solve`. |
| AD2 | **Indicators** per edge, each normalized to a dimensionless "excess" e ≥ 0 (0 = resolved): Δψ / (k V_T) with k = 1 default; \|Δ ln n\|, \|Δ ln p\| / ln 10 per decade; \|Δ ln \|N_D − N_A\|\| for doping; generation-rate integral fraction for II/BTBT (only when the model is on). A node's indicator is the largest over its edges. | SG accuracy degrades with Δψ ≫ V_T (standard); the rest heuristic, measured by AC5–AC7. | Cheap, device-neutral, computed from fields the solver already has. | Thresholds are parameters with defaults; the defaults are fixed in slice 2 from AC4 measurements and recorded. |
| AD3 | **Goal weighting:** node weight ω_i = \|z_i\| scaled by the node's share of the goal (normalized so Σ ω_i = 1 over semiconductor nodes); marked if e_i ω_i ≥ θ max(e ω) (Dörfler-style bulk marking, θ = 0.3 default). Feature-only mode: ω ≡ 1. | Dörfler marking: standard in adaptive FEM; the adjoint weight: dual-weighted residual idea, here as a weight, not a bound. | Concentrates nodes where the goal is sensitive. | An estimate, not a bound; AC7 measures it. |
| AD4 | **Size field:** h_new = h_old / 2 at marked nodes (per pass), h_old elsewhere; then graded so h changes by at most a factor g = 1.5 per edge (fixed-point sweep), and never below `min_size` (default s / 4 where interfaces exist, else h_initial / 64). The generator queries it by exact point location on the old mesh (barycentric, linear). Tensor: a line inserted at each marked cell's midpoint along each axis where that axis' indicator marks it. | Grading: Unit 49's sliver finding. | Fed into the generators' existing `size_at`. | A new option field per generator (`size_field`); absent = bit-identical. |
| AD5 | **State transfer:** ψ, the quasi-Fermi potentials (not densities: densities interpolate badly across decades), T and carrier temperatures when present, linearly interpolated at the new nodes by exact location; interface and contact quantities recomputed. Used as Newton's initial state only: the result on the new mesh is a converged solve, so every conservation gate applies unchanged. Also exposed as `transfer_state` so a checkpoint can restart on a remeshed device (IO6's gap), still refused by digest unless transferred explicitly. | — | Fewer Newton iterations per pass; checkpoints usable after remeshing. | Interpolation across a straddled interface uses the node's own side. |
| AD6 | **Stopping:** the loop proposes to stop when \|G_k − G_{k−1}\| ≤ tol·\|G_k\| for every goal and the adjoint estimate ≤ tol; the stop is then subject to AD11's verification. Otherwise it ends at `max_passes` (default 6), `max_nodes` or a DOF/time budget, uncertified (AD13). The reason recorded. | — | The user's stated accuracy, not a mesh size. | Two consecutive small changes can be stagnation; the estimate and AD11 guard it. |
| AD7 | **Deck:** `[adapt]` section: `goal = current <contact>` (repeatable) or `vt`, `gm`, `swing`; `tolerance`; `at = <bias points>` or `auto`; `indicators = auto` or a list; `weighting = adjoint` or `none`; `policy = intelligent` or `basic`; `mode = tolerance` or `budget` with `budget_dofs` / `budget_seconds`; `verify = richardson` or `off`; `certify = required` or `report`; `max_passes`, `max_nodes`, `min_size`, `grading`, `theta`; (17b) `predictor = <file>`. Absent: no adaptation, bit-identical. | DK conventions. | The batch path reaches it without code. | — |
| AD8 | **Results, io, run record, app:** the adaptation report (per pass: nodes, goals, estimate, marks by indicator and region, time, stopping reason) in results and the run record; the final mesh saved as a generated mesh with its digest; the run identity covers the `[adapt]` settings and the final mesh digest; the app shows the passes' meshes and the report. | Unit 42 provenance rules. | Reproducible and explainable. | — |
| AD9 | **Determinism:** marking ties broken by node index; the size field and regeneration are deterministic; two runs bit-identical on one machine and thread count (6.8). | — | — | — |
| AD10 | **Existing paths:** with no `[adapt]` section and no `size_field`, every generator, device, result, file and run identity is bit-identical. | — | — | — |
| AD11 | **Independent numerical error verification (OA2):** before a result is certified, the final mesh M is refined uniformly once (every size field value halved, every tensor and axis line interval split; layers' lines halved) to M′, the analysis is solved on M′ from the transferred state, and the verified error of each goal is the Richardson estimate \|G_{M′} − G_M\| / (2^p − 1) with p the *observed* order where a third mesh exists from the loop's passes (clamped to [0.5, 3]), else p = 1 (conservative for the box method at abrupt features). This uses neither the indicators, nor the adjoint, nor any prediction: a second, independent check of AD3's estimate. The solve on M′ is reported but is not the result (the result stays on M, one mesh per sweep). `verify = richardson` is the default; `verify = off` skips it and the result is marked unverified, never certified. | Richardson extrapolation (standard); its independence from the adjoint estimate is by construction. | An adjoint estimate can be wrong where the linearization is poor (strong injection, near breakdown); a refined solve is not linearized. | Cost: one solve on about 4× (2D) or 4–8× (swept) the DOFs, measured (AC11). |
| AD12 | **Conventional AMR fallback (OA2):** the run switches from `policy = intelligent` to `policy = basic` (AD3's marking, no predictions) for every remaining pass when (a) two consecutive proposals are misses (ID1), (b) a pass lowers neither the estimate nor the goal change, or (c) AD11 rejects a proposed stop (verified error > tol): the next pass then marks by AD3 using the verified error's distribution. The switch, its trigger and its pass are recorded. | — | The intelligent layer can only speed up a conventional AMR, never replace its guarantees. | A run that falls back is no faster than basic, but no less accurate. |
| AD13 | **Explicit accuracy certification (OA2):** every adaptive result carries, per goal: the value, the adjoint estimate (AD3), the verified error (AD11) and the order used, the tolerance, the policy that produced the mesh and whether a fallback occurred, and a status: **certified** (estimate ≤ tol *and* verified error ≤ tol), **not certified** (either exceeds tol, with the margin) or **unverified** (`verify = off`, or a budget stop). The status, numbers and method are in results, the run record, the CSV table and the app; `nitcad run` exits with a distinct non-zero code when a requested `certify = required` run ends uncertified. | — | A user sees, per result, whether the stated accuracy was reached and how it was checked. | A certificate is evidence from two estimates, not a proof; the text says so. |

**The intelligent decision layer (`solve::adapt_policy`; replaces AD3's fixed-θ marking when on, `[adapt] policy =
intelligent`, the default once AC14–AC17 pass; `policy = basic` keeps AD3):**

| # | Choice | Evidence | Why | Trade-off |
|---|---|---|---|---|
| ID1 | **Propose, then verify.** The layer proposes the next size field; the base loop builds, transfers, solves and estimates. Stopping (AD6), verification (AD11) and certification (AD13) never use a prediction; only solves. A proposal that does not lower the estimate is recorded as a miss and AD3's marking is used for that pass; two consecutive misses hand the run to conventional AMR (AD12). | — | Predictions can be wrong; results cannot. | A miss costs one pass. |
| ID2 | **A priori prediction (before any solve):** a predicted size field from the device alone. Junctions: depletion width W = sqrt(2 ε (V_bi + \|V\|) / (q N_eff)) at the largest bias of the analysis, from the sampled doping (M15); the Debye length L_D = sqrt(ε V_T / (q N)) sets h ≤ L_D / 2 across each junction's predicted span. Insulator interfaces under gates: the inversion-layer scale from the gate's largest bias and oxide capacitance (classical; MLDA's λ when on). Contacts: none. Active models: II and BTBT where the predicted junction field exceeds the model's onset. | Textbook depletion and Debye scales; heuristic as a mesh rule. | Pass 0 starts near the answer instead of from uniform h: fewer passes. | A junction whose bias the analysis does not reach is over-predicted; ID5 learns it down. |
| ID3 | **Prediction across bias:** with one mesh per sweep (AQ3), the layer solves the sweep coarsely once, tracks each node's indicator over bias, and chooses the adaptation points where the indicator field changes most (the largest change in its weighted norm between neighbours), up to `max_points`; features that move with bias (pinch-off, depletion edges) are covered by the union of their predicted positions. | — | "auto" adaptation points chosen from the physics instead of the sweep's ends. | One coarse sweep of cost. |
| ID4 | **Benefit estimate per candidate:** candidates are graded patches (a marked node and its neighbours to the grading radius), tensor lines, or swept axis lines. Predicted error reduction ΔE_c = Σ_{i∈c} ω_i η_i (1 − 2^{−p_c}), η the node's indicator excess (AD2) in the goal's units through the adjoint (AD3), p_c the local convergence order. **Predicted cost (OA4), not node count alone:** C_c = ΔT_pass(DOF + ΔDOF_c, nnz + Δnnz_c) − ΔT_pass(DOF, nnz), with ΔDOF_c = added nodes × the equation set's unknowns per node (3 drift-diffusion, 4 electrothermal, 5 energy balance; insulator nodes 1), Δnnz_c the Jacobian nonzeros the candidate adds (added edges × stride², counted from the candidate's predicted local mesh: mean degree from the current mesh's region, 2D and swept separately), and ΔT_pass the cost model of ID5. | The (1 − 2^{−p}) factor: halving h with local order p (standard a priori estimate). Sparse direct cost grows faster than DOFs (fill-in), so DOFs and nonzeros, not nodes, predict time. | Benefit in the goal's units, cost in seconds: comparable across candidates and across equation sets. | p and the cost model are learned (ID5); before they are, p = 1 (the box method's order at abrupt features, measured in Units 16 and 52) and the cost model's defaults from AC11's measurements. |
| ID5 | **Online learning (LQ1 a):** after each pass the realized error change per patch (the estimate before and after, through the adjoint) gives a fitted local order p per indicator kind and region (least squares in log2, clamped to [0.5, 3]). **The cost model (OA4)** is fitted from the passes' own measured timings: T_pass = N_newton × (a · nnz(J)^b + c · nnz(L+U) + d · DOF) + T_mesh(nodes) + T_transfer, with the factorization's fill nnz(L+U) as reported by the backend (SparseLU and PARDISO report it; GMRES-ILUT reports the preconditioner's), N_newton the measured Newton iterations per point, and T_mesh, T_transfer measured; b clamped to [1, 2]. Each fit's sample count, residual and predicted-vs-measured pass time go in the report. | — | The benefit and cost estimates become the device's and machine's own, not defaults. | Few passes, few samples: fits regularized to the defaults (ridge toward p = 1, b = 1.5). The decision layer only; AD11 and AD13 never read them (OA6). |
| ID6 | **Allocation (LQ3; a heuristic, OA1):** tolerance mode: take candidates in decreasing ΔE / C until the predicted error ≤ tol / 2 (margin for misprediction) or the pass' DOF growth cap (×2.5 per pass); budget mode: the same order until the DOF or time budget is spent. Ties by node index. **No optimality is claimed:** greedy by ratio is optimal only for a fractional knapsack with fixed, additive values; here the values are predictions, overlapping patches and grading couple the candidates, and refining one region changes the error elsewhere, so the allocation is a heuristic. | Fractional knapsack's greedy optimality is the motivation only, not a property of this problem. | A simple, explainable ordering that is measured against its baseline. | Can be beaten by a better heuristic; AC16 and ID9 measure it against AD3 and record the gap. |
| ID7 | **Explanation:** per pass the report lists the top candidates (location, region, indicator, ΔE predicted vs realized, cost), the fitted p and cost model, misses (ID1), and the predicted vs actual node count and error. | — | A user can see why the mesh is where it is. | Larger reports (capped at the top 20 per pass). |
| ID8 | **Cross-run predictor (Unit 17b if LQ4 b):** training samples are (node features → realized benefit) from saved adaptation histories (every run with `[adapt]` can save its history, opt-in). Features: material, region kind, log \|N_D − N_A\|, \|∇ log N\|, distance to the nearest junction, insulator interface and contact, predicted W and L_D (ID2), the analysis' bias range, the goal kind. Model: gradient-boosted trees (LQ2), fixed hyperparameters, deterministic training (sorted samples, no randomness); `[adapt] predictor = <file>` replaces ID2's a priori field with the model's, ID1 and AD12 still guarding. **Opt-in only (OA5):** off unless `predictor` is given; never on by default, never shipped. **Reproducible training data (OA5):** an adaptation history file (versioned, checksummed, atomic, as the `io` layer's files) holds the deck text, the device and mesh digests of every pass, the provenance (Unit 42), the per-node features and realized benefits, and the feature-schema version; `nitcad learn` builds a training set from a list of history files into a manifest (each file's digest, sorted), so the same manifest always gives the same samples. | — | Learns the user's device families (a MOSFET process with varying L, a diode series). | Trained on the user's own runs only; extrapolation outside them falls back to ID2 (features out of the training range are flagged). |
| ID9 | **The layer must beat its baseline:** every intelligence gate (AC14–AC17) compares against `policy = basic` (AD3), and 17b's against 17's online layer and against ridge regression; a component that does not beat its baseline is not made the default. | — | "Intelligent" is measured, not asserted. | — |
| ID10 | **Model versioning and held-out validation (17b, OA5):** a `.nitlearn` file records its format version, the feature-schema version, the model class and hyperparameters, the training manifest's digest, the build that trained it (Unit 42 provenance), and its held-out metrics. **Held-out devices, not held-out nodes:** the training command splits by device (history file), never by node, so no device is in both sets; the metrics (rank correlation of predicted vs realized benefit, recall/precision of predicted refinement, passes saved) are measured on the held-out devices and stored in the file. A model whose schema or format version differs from the build's is refused (`invalid_input`, both versions named); the model's digest and version enter the run identity and the run record. | Device-level splitting: nodes of one device are strongly correlated, so a node split overstates accuracy (standard leakage argument). | Any prediction a run used can be traced to the exact model, data and code that made it. | Retraining after a feature change is required; old models are refused, not silently reused. |

**Acceptance criteria, fixed before any implementation:**

| # | Criterion |
|---|---|
| AC1 | Nothing changes: every existing test and pin passes; every existing deck's run identity, mesh digests and device files byte-identical. |
| AC2 | Size field honoured: on the 2D generator, swept and tensor paths, every element's size ≤ the field at its centroid (the generator's own bound), grading ≤ g measured, certificates P0–P5, P4, P6 as for any generated mesh; a constant field equal to h gives the h mesh's digest. |
| AC3 | Transfer: a linear ψ and linear quasi-Fermi potentials transfer exactly (1e-12); a diode and a MOSFET state transferred to a refined mesh converge to the same solution as a solve from scratch (1e-9 relative in currents) in fewer Newton iterations (measured, reported). |
| AC4 | Indicators: on a 1D-in-2D abrupt diode the doping and Δψ indicators mark the junction's depletion edges only; on a MOS-C in inversion the Δ ln n indicator marks the inversion layer; each indicator's edge values equal hand-computed ones on a 4-node fixture. |
| AC5 | **Efficiency gate (2D):** the 2D MOSFET (Unit 16's deck shape) from a coarse mesh, goal = I_D at the sweep's points: the adaptive mesh reaches \|I_D − I_ref\| ≤ 1% at every point above the current resolution (I_ref from a 4× uniformly refined reference) with at most 1/3 of the nodes the uniform-refinement ladder needs for the same error. Risk stated before the run: if it fails, `[!shouldfail]` with its measurement for an owner ruling; never loosened. |
| AC6 | **Diode gate:** the legacy diode's forward and reverse currents reach the 2× reference within 0.1% with fewer nodes than the existing hand-made mesh; the loop stops by AD6's goal rule, not a budget. |
| AC7 | **Goal weighting earns its place:** a MOSFET with a deep well junction off the current path; adjoint-weighted vs feature-only to the same I_D error: the weighted mesh has ≥ 25% fewer nodes, and its marks avoid the well. Risk stated: if it fails, `[!shouldfail]`, owner ruling. The adjoint estimate is within a factor 3 of the actual error at the final pass. |
| AC8 | **FinFET (swept):** from a coarse swept deck, adapt for I_D; report nodes and FC9/UC8's rule against the 2× tensor reference at every point. Measured target: UC8's rule at all 17 points within UC8's 9,367-node budget; stated before the run, not a pass/fail gate (UC8 failed 5/17 by hand). |
| AC9 | Deck, io, run record, app: `nitcad check` validates `[adapt]` with line numbers; the report and final mesh round-trip; the run identity covers settings and mesh digest; `nitcad plot` and the app show the passes. |
| AC10 | Determinism: two runs bit-identical (AD9); refusals (unknown goal, contact, indicator; goal-weighting with a backend that cannot transpose and no fallback; adapt on general 3D in this unit) with code and message. |
| AC11 | Cost: per pass, DOFs, Jacobian nonzeros, factorization fill, Newton iterations, and regeneration + transfer + solve time reported for AC5, AC6, AC8; AD11's verification solve reported separately (time and DOFs); the transposed solve ≤ the cost of one Newton iteration's linear solve. |
| AC12 | Hand mutants: at least one per design row (AD1–AD13, ID1–ID7, ID9; 17b: ID8, ID10), each survivor closed by a test or argued equivalent. |
| AC13 | Release and Debug suites at every commit, `/W4 /WX`; Debug grows by ≤ 60 s; larger cases `[.benchmark]`. |
| AC14 | **Prediction (ID2, ID3), recall and false positives (OA3):** on the diode, the MOS-C and the 2D MOSFET, measured against the region the basic loop finally refines (nodes whose final size is below the initial h): **recall** (the a priori field covers ≥ 80% of it); **false-positive refinement** reported for every device: the precision (the fraction of the predicted-refined region that the basic loop also refined; its complement is the false-positive fraction), the false-positive DOFs and nonzeros, and their cost (the time of the final pass with and without them, the latter by re-solving on the basic loop's mesh), gated: false-positive DOFs ≤ 30% of the intelligent final mesh's DOFs; **mesh overhead**: the intelligent final mesh's DOFs ≤ 1.1 × basic's at the same verified accuracy. And the intelligent loop reaches AC5's and AC6's accuracy in at least one pass fewer than `policy = basic`, summed over the three. Risk stated: `[!shouldfail]` and owner ruling if not; never loosened. |
| AC15 | **Benefit and cost estimates (ID4, ID5, OA4):** over all candidates refined in AC5's run, the predicted and realized error reductions have Spearman rank correlation ≥ 0.6; the fitted orders on the abrupt diode lie within [0.8, 1.3] (first order expected); from the third pass on, the cost model predicts each pass' measured time within ±30% on AC5, AC6 and AC8, and does so better than a node-count-only model fitted to the same passes (mean absolute relative error, both reported). |
| AC16 | **Allocation (ID6, a heuristic, OA1):** tolerance mode: the intelligent loop's final DOFs ≤ 0.8 × `policy = basic`'s at the same verified accuracy on AC5's MOSFET and AC7's well MOSFET; budget mode: at a DOF budget of 50% of basic's final DOFs, the intelligent loop's verified error is ≤ basic's at that budget. Results stated as measured comparisons, with no optimality claim. Risk stated as AC14. |
| AC17 | **The guard and the fallback (ID1, AD12):** a deliberately wrong prediction (a hand-written field refining only the substrate) still reaches the tolerance and is certified, with the misses reported, the switch to `policy = basic` recorded at the pass AD12 names, and the result within the tolerance of the reference; total time ≤ 2 × basic's. Each AD12 trigger (a, b, c) is forced once in a fixture and its switch recorded. |
| AC18 | Determinism and explanation: two runs bit-identical, report and certificate included; the report's predicted vs realized numbers equal those recomputed from the saved passes. |
| AC19 (17b) | **Cross-run learning (OA5):** a family of 8 MOSFETs (L, oxide and doping varied, generated decks); held-out by device: trained on 6, tested on the 2 held out: the learned first mesh needs ≥ 1 pass fewer than ID2's a priori field and ends with ≤ the DOFs at the same verified accuracy; gradient-boosted trees beat ridge regression on held-out benefit prediction (rank correlation); the training set built twice from the same manifest is bit-identical, and so is the trained model; a model with another schema or format version is refused; the run identity changes with the model's digest; an out-of-range device is flagged and falls back to ID2; with `predictor` absent every run is bit-identical to Unit 17's. |
| AC20 | **Certification (AD11, AD13, OA2):** on AC5, AC6, AC7 and AC8 every adaptive result carries a certificate; against the references (AC5's 4× uniform, AC6's 2×, AC8's 2× tensor) no result marked **certified** has an actual error above its tolerance; a run made to stop early (a `max_nodes` budget) is marked **unverified** or **not certified**, never certified; with `certify = required` such a run exits with the distinct code; the verified error and the adjoint estimate both reported, and their ratio to the actual error recorded. Mutants on AD11 (skipped refinement, wrong p, wrong mesh compared) and AD13 (status logic) are killed. |
| AC21 | **Accuracy independent of predictions (OA6):** for AC5's and AC6's final meshes, the estimate, verified error and certificate computed with `policy = basic`, `policy = intelligent` and (17b) a `predictor` given are bit-identical when the final mesh is the same (a fixture forcing the same mesh), and a mutant making AD11 or AD13 read any decision-layer quantity is killed by this test. |

**Slices** (each a local commit with both suites green):
1. Size fields in the 2D generator, swept and tensor paths; exact point location on the old mesh; grading
   (AD4, AD10; AC1, AC2).
2. Indicators and marking, feature-only (AD2, AD3 with ω ≡ 1; AC4).
3. State transfer, the loop with stopping, independent verification and the certificate; the diode gate
   (AD1, AD5, AD6, AD9, AD11, AD13; AC3, AC6, AC20 part).
4. Transposed solve and goal weighting (AD3; AC7).
5. The decision layer: a priori and bias prediction, benefit, the DOF/sparsity/timing cost model, online fits,
   allocation, guard and conventional fallback, explanation (ID1–ID7, ID9, AD12; AC14–AC18, AC21).
6. Deck, results, io, run record, app, certificate display and exit code (AD7, AD8, AD13; AC9, AC10).
7. Gates and cost: AC5, AC8, AC11, AC20, mutants AC12, ARCHITECTURE 6.12.4 and 11.3.

Unit 17b (approved split, LQ4 b; a separate branch and unit after 17 merges): adaptation-history files and
`nitcad learn` manifests, features and schema versioning, deterministic gradient-boosted trees and the ridge
baseline, device-level held-out validation, versioned `.nitlearn` io, deck `predictor` (opt-in), AC19 and AC21's
predictor case (ID8, ID10); 3 slices.

**Size:** Unit 17 XL (7 slices; 1, 3, 4 and 5 the bulk); Unit 17b L.

**Implementation record (from 2026-11-12).**

*Slice 1 — size fields (AD4, AD10; AC1, AC2).* `mesh/size_field.hpp`: `SizeField` (a nodal size on a background
triangulation; exact location by `orient2d` through a bucket grid, the lowest-indexed containing triangle on a shared
edge, linear interpolation clamped to the corners' range; outside every triangle the value at the nearest point of the
nearest triangle), `local_sizes` (a node's longest incident edge, the h_old of AD4), `grade_sizes`, `refine_axis`.
`MeshOptions2D::size_field` (a `shared_ptr<const SizeField>`, null by default) enters the 2D generator only through
`size_at` (`min` with the field), so a field only refines; the swept generator passes its section's field through and
adds the field's digest to `generation_digest` only when one is given (AD10: no word added without it). Two
construction points found while writing the tests, fixed before any gate run (no criterion touched):
- *Grading follows the refinement only.* Grading every node of the old mesh (AD4's "fixed-point sweep" read
  literally) spreads the s-scale sizes of every interface pair and straddle over the whole device, refining where
  nothing was marked. `grade_sizes` therefore takes seed nodes (the nodes a pass refined) and propagates g^k h from
  them only (Dijkstra on log h); with every node a seed it is the full grading (h_j ≤ g h_i on every edge, tested).
  `refine_axis` likewise grades against intervals split in the call, keeping the given axis' own ratios.
- *Axis grading factor at least 2.* A midpoint split leaves a ratio of exactly 2 beside it, so with g = 1.5 every split
  cascades over the whole axis. On axes the factor is max(g, 2) (the 2D field keeps g = 1.5). Stated in the header.
Tests (`tests/mesh/size_field_test.cpp`, 6 cases): linear reproduction (1e-12) and the outside rule; refusals;
local sizes and grading (all seeds, seeds only); axis refinement (marked, second-level grading, protected straddles,
min_interval, refusals, determinism); AC2 on the 2D generator (a constant field equal to h gives the h mesh's
`graph_digest` and points; a left-half h/4 graded field: P4 certificate kept, more nodes, concentrated where the field
is small, every triangle not exempt within the generator's bound circumradius ≤ size/√3 against the field at its
centroid, deterministic); AC2 on the swept path (same digest of the graph with a non-refining field, generation digest
names the field, a refining field adds nodes in whole axis lines).

*Slice 2 — indicators and feature-only marking (AD2, AD3 with ω ≡ 1; AC4).* `solve/adapt_indicators.hpp`:
`edge_indicators` (raw values per edge and kind), `compute_indicators` (excess e = max(0, r/t − 1), node max, dominant
kind), `unite` (AQ3's union over adaptation points), `mark_nodes` (Dörfler: e w > 0 and e w ≥ θ max). Insulator edges
carry no indicator (the potential is linear between charges there, exact for the box method). The generation indicator
is a node's share of the total G·V times the number of semiconductor nodes, t = 10.

| # | Question (put and ruled under the owner's standing instruction "always choose recommended", 2026-11-11) | Options | Recommendation and ruling |
|---|---|---|---|
| AR1 | AD2's doping indicator \|Δ ln \|N_D − N_A\|\| is near zero across an abrupt junction between equal and opposite dopings (1e17 to −1e17: Δ ln = 0), so AC4's "the doping indicator marks the junction" cannot be met as written. | (a) As written. (b) **The signed log L(N) = sign(N) ln(max(\|N\|, n_i)/n_i), indicator \|ΔL\| / ln 10**: on an edge of one sign above n_i it is AD2's formula exactly; across a junction it counts the decades down to n_i and back up. (c) The equilibrium potential step \|Δ asinh(N / 2n_i)\|. | **(b)**, ruled as recommended: it keeps AD2's meaning and units (decades of net doping) where AD2's formula works and needs no solve; (c) changes the unit. No criterion changed. |

AC4 run (`tests/solve/adapt_indicators_test.cpp`, 4 cases, 342 assertions): the four-node fixture's edge values equal
the hand-computed ones (potential, both densities, doping by AR1, generation share 4 on node 3's edges); on the abrupt
1e17 diode the doping indicator alone marks only the nodes beside the junction and the potential indicator only nodes
within W/2 + h of it (W the abrupt-junction depletion width), zero excess beyond W/2 + 2h; on a MOS-C (p 1e17, 1.5 nm
oxide, n+ poly, 1 V, the quasi-static `equilibrium_poisson` solve: drift-diffusion stalls above threshold on a MOS-C,
`bias.hpp`) the electron-density indicator alone marks only nodes within one spacing of the inversion layer (n > N_A),
the surface node included. Fixture note: a constant profile's window along its depth axis (x by default) is ignored
(`doping_profiles.hpp`), so the x-windows take `axis=y`.

*Slice 3 — transfer, the loop, verification and the certificate (AD1, AD5, AD6, AD9, AD11, AD13; AC3, AC6, AC20
part).* `mesh/locate.hpp` (`PointLocator`: simplex meshes in 2D and 3D by orientation predicates through a bucket
grid, tensor meshes by `TensorCells::stencil`; outside: the nearest node); `GeneratedMesh::section` (a swept mesh's
cross-section simplices, node p of line j being j·n2 + p; not saved in files, so no file or digest changes);
`solve/adapt_transfer.hpp` (`transfer_state`: ψ, E_Fn, E_Fp and the density offsets c_n = ln n − (ψ + E_Fn)/V_T,
c_p = ln p + (ψ + E_Fp)/V_T interpolated over the old nodes of the new node's region; n, p rebuilt from them, so a
linear ψ and linear quasi-Fermi levels transfer exactly where c is constant, one material without narrowing);
`solve/adapt.hpp` (`adapt`, `next_refinement`, the builder callback `DeviceBuilder` and `MeshRefinement`). Without
adjoint weighting (slice 4) a goal's estimate is its last change between passes.

| # | Question (ruled under the standing instruction "always choose recommended") | Options | Recommendation and ruling |
|---|---|---|---|
| AR2 | What the loop does when AD11 rejects a proposed stop. Found on AC6's first run: the reverse current's change between passes fell below tolerance at 70 nodes while its error against the 2× reference was 0.27%; AD11 rejected the stop (its verified error above the allowance; the estimate, 1.0e-12, within it), but the slice-3 loop as first written ended there, not certified. Continuing with indicator marks did not help (12 passes, the value settled 0.26% off: the indicators kept marking the same few nodes, so the missing accuracy lies where they see none; where exactly was not measured). | (a) End, not certified. (b) **Continue, the next pass marked by AD3 on the verified error's distribution** (AD12 c's wording, brought forward from slice 5 into the basic loop): per node of M, the largest over the adaptation points of \|ψ′ − ψ\|/V_T, \|ln n′ − ln n\|, \|ln p′ − ln p\| with M′'s solution interpolated at the node; uniform refinement (M′) when that marks nothing. (c) Continue with M′ itself. | **(b)**, ruled as recommended: it is AD12 (c) as frozen, reads solves only (OA6), and refines where the error is rather than everywhere; (c) multiplies 2D meshes by four per rejection. Slice 5 adds AD12's switch from the intelligent policy to this one; the basic loop has it from slice 3. No criterion changed. |

Runs (Release, `build/r3`):
- AC3: linear ψ and quasi-Fermi levels exact to 1e-12 tensor → tensor (5×4 to 9×7, not nested) and generated →
  generated (0.3 µm to 0.11 µm); the locator reproduces linear fields on tensor, 2D and swept 3D meshes. The 1e17
  diode at 0.5 V (41 → 97 nodes along x): currents equal to 1e-9, Newton 10 from scratch, 5 transferred; a lateral
  n+ p n+ MOSFET (gate 1 V, drain 0.5 V; 41×16 → 73×29): equal to 1e-9, Newton 12 from scratch, 8 transferred.
- AC6 (the legacy diode, weighted current on: the reverse current, 3.1e-9 A/cm², lies below the terminal current's
  resolution of about 4e-8): from 21 uniform nodes, 10 passes, AD11 rejecting two proposed stops (at 70 and 107
  nodes; AR2 then marked 30 and 64 nodes); stopped by the goal rule at **180 nodes** (the
  hand-made mesh: 250); I(−1 V) −3.09402e-9 against the 2× reference's −3.09688e-9 (9.2e-4, within AC6's 1e-3, a thin
  margin), I(0.6 V) 0.606567 against 0.606549 (2.9e-5); both certified; verification on 370 nodes in 0.02 s.
- AC20 part: a `max_nodes` stop and a `max_passes` stop are unverified (no verification run, verified error NaN);
  `verify = off` gives unverified after a goal stop. AD9: two runs bit-identical (node counts, goals, points,
  verified error, status). AC10 part: the refusals named in `adapt.hpp`.

*Slice 4 — the transposed solve, the goal adjoint, adjoint weighting and the estimate (AQ4, AD3, AD6; AC7).*
`linalg::LinearSolver::solve_transposed` and `can_transpose` (the equilibrated system transposed, (R A C)ᵀ y = C b,
x = R y, refined against the transposed componentwise backward error; SparseLU through Eigen's `transpose()` view,
PARDISO with iparm[11] = 2; gmres_ilut factorizes its sparse_lu fallback and reports `fell_back`, and without a
fallback cannot transpose, so `weighting = adjoint` is refused for it); `solve/adapt_adjoint.hpp` (`check_adjoint`,
`point_adjoint`, `goal_at`, `estimate_goals`); `TransferMode::exponential_fit` in `transfer_state`; the loop's
adjoint branch in `adapt.cpp`. The goal is the terminal current, or the weighted current G = T_c + Σ w_i B_i (Unit
32); its gradient is the small-signal currents' rows at s = 0 plus the weighted continuity rows, and the adjoint z
solves Jᵀ z = ∂G/∂u at the converged state.

| # | Question (ruled under the standing instruction "always choose recommended") | Options | Recommendation and ruling |
|---|---|---|---|
| AR3 | How AD3's weights and AD6's estimate are computed from the adjoint (the design fixes the dual-weighted residual, not its discretization). Measured on the legacy diode (41 to 321 uniform nodes, against the solve on M′): M's adjoint interpolated onto M′ with M's state transferred by quasi-Fermi levels gave estimates off by orders of magnitude (the transfer's residual on M′ dominated); the exponential-fit transfer with M's adjoint still missed by factors of 5 to 50 at 0.6 V; a Newton-step (nonlinear) evaluation was worse. | (a) M's adjoint interpolated to M′. (b) **M′'s own adjoint at M's state prolongated to M′ by the Scharfetter-Gummel exponential fit** (n = Σλ_j n_j e^{(ψ−ψ_j)/2V_T} / Σλ_j e^{(ψ_j−ψ)/2V_T}, holes with the signs reversed; M's recombination densities N_s kept), one factorization of M′'s Jacobian per point: D = G′(u) − G_M − z′ᵀR′(u), the estimate 2\|D\| (p = 1, as AD11 without an observed order), the weights z′_q R′_q summed onto the coarse node with the largest interpolation weight in the same region and normalized over the semiconductor nodes, ω the largest over goals and points. (c) Uniform M′ solves each pass (AD11's cost every pass). | **(b)**, ruled as recommended: it is the dual-weighted residual on M′ with one linear solve instead of (c)'s nonlinear one, and the only form measured to work: at 0.6 V within 30% of G_M′ − G_M on every mesh; at −1 V within a factor 1.6 up to 161 nodes and 10 at 321 (a change of 1e-3 of the current); at 0.3 V the uniform meshes are not asymptotic (the change flips sign from 41 to 161 nodes), within 1.3 at 321. The cost is M′'s assembly and one factorization per adaptation point and pass (AC7: 2993 to 9570 nodes, 0.07 to 0.1 s). AD11 reuses the pass's M′. AR2's difference marking stays unweighted. No criterion changed. |
| AR4 | AC7's "fine uniform reference". The test first refined the start 4× (10449 nodes), then extrapolated 2×/4×/8× (5.18618, order 1.415), against which the adjoint estimate missed AC7's factor 3 (5.8). A 16× solve (164673 nodes, 12 minutes) moved the extrapolation to 5.20577 (4×/8×/16×, order 1.165): the uniform sequence is not asymptotic, so the reference is uncertain by about 0.4%, the size of the errors measured. | (a) Solve 16× in the test (12 minutes Release, far longer Debug). (b) **Record the 4×/8×/16× extrapolation as a measured constant, with the 2× solve repeated as the guard that the deck and physics are unchanged (to 1e-5).** (c) Keep 2×/4×/8×. | **(b)**, ruled as recommended: the finest reference measured, at suite cost; (c) uses a reference known to be 0.4% coarser. The reference's own uncertainty is reported to the owner: against (c) the estimate check fails, against (b) it passes; neither 32× nor any other independent reference was run. No criterion changed. |

Runs (Release, `build/r3`):
- linalg: transposed solves match factorizing Aᵀ (real and complex, transpose not conjugate; backward error within
  1e-9 measured 2.4e-10); PARDISO transposed with MKL; gmres_ilut falls back, or refuses without a fallback.
- AQ4: the goal's gradient against central differences on the diode (terminal and weighted current at −1, 0.3 and
  0.6 V): within 1e-6 of the terms' magnitude in five cases; **the weighted current at −1 V differs by 5.8e-3** of the
  derivative (open finding: not isolated; not SRH, Auger, radiative or narrowing alone; the test bounds it at 1e-2
  in that case only).
- AC7 (a MOSFET in a p-well, 1 µm × 1.5 µm, gate 1.5 V, drain 0.5 V; I_D's reference 5.20577, AR4): from 693 nodes,
  feature-only stopped by the goal rule at 6834 nodes (I_D 5.18975, 0.31%, certified); adjoint-weighted at **2214
  nodes** (I_D 5.18151, 0.47%, certified; estimate 0.0271, AD11's verified error 0.0156, actual 0.0243, estimate /
  actual 1.12); node ratio 0.32 (AC7: at most 0.75); y-lines within 0.15 µm of the well junction 3 (the initial
  3) against feature-only's 19. 29 s.
- Open finding for the owner: AD11's verified error \|G_M′ − G_M\| / (2^p − 1) estimates M′'s error, not M's (M's
  is 2^p times it, \|G_M′ − G_M\| 2^p / (2^p − 1)); frozen, not changed.

*Slice 5a — the decision layer's prediction, guard and fallback (ID1–ID3, ID7 part, AD12; AC14, AC17, AC21).*
`solve/adapt_policy.hpp` (`Policy`, `predict_sizes`, `edge_exceeds`, `choose_points`, `PolicyOptions`); the loop's
`policy = intelligent` branch in `adapt.cpp` (the coarse sweep and point choice, the prediction's rounds before pass 0,
the miss rule, AD12's three triggers, `MeshSource` per pass, the fallback and misses in the result). The default stays
`policy = basic` (ID9: not the default until AC14–AC17 pass). Slice 5b adds ID4–ID6 and ID5's cost model.

| # | Question (ruled under the standing instruction "always choose recommended") | Options | Recommendation and ruling |
|---|---|---|---|
| AR5 | The constructions the decision layer's rows leave open. | (a) as below; (b) other thresholds. | **(a)**, ruled as recommended; no criterion changed. **Split:** slice 5 is two commits, 5a (ID1–ID3, AD12, ID7's per-pass record) and 5b (ID4–ID6, ID5, ID7's candidates). **Miss (ID1):** a proposal "lowers the estimate" when it falls by at least the miss margin (10%) of the previous pass's; the estimate compared is the largest over goals and points of estimate / allowance (the adjoint's, or the last change without weighting); a proposal that refines nothing is a miss without a solve; pass 0's mesh is judged only against ID3's coarse sweep's adjoint estimate. **Prediction (ID2):** applied before pass 0 in rounds (each round halves sizes toward the targets, at most 12, recomputed on each round's device), stopping before ID6's growth cap (2.5× the start's DOFs; the prediction is pass 0's proposal) or max_nodes; a round that would cross the cap is applied in part, to the edges of largest length / target ratio, the share found by bisection on the sorted ratios (at most 12 trial meshes), and the rounds then stop. Targets are anisotropic on tensor and swept axes (an edge's length projected on the feature's normal) and isotropic in 2D size fields. Under gates the depletion region is added (L_D / 2 within W_d, the junction's own rule; the inversion length alone covered almost nothing the basic loop refines on the MOS-C). Impact ionization and band-to-band tunnelling add nothing beyond their junction's target. **Points (ID3):** a point where a goal is at or below its current resolution is not eligible (its allowance tol·r cannot be met, AD6); the coarse sweep is pass 0's solve when the prediction refined nothing. **Hand-written prediction:** `PolicyOptions::prediction` replaces ID2 and every later proposal (AC17's wrong field). **AC17 on the diode** checks the guard; its reference is the 4× hand mesh (the 2× moves 4.7e-4 from it at −0.5 V, half the tolerance); its time bound is measured on the 2D MOSFET, the device whose substrate the criterion names. **AC13:** AC7 (28 minutes in Debug) and AC14 are hidden gates run in the Release configuration only, as the solve layer's slow gates are. |
| AR6 | AC14's devices and its measured miss. | (a) as below; (b) the generated 2D MOSFET. | **(a)**, ruled as recommended; no criterion changed. **Devices:** the diode (AC6's run), the MOS-C (substrate current, degenerate zero goal), and the 2D MOSFET on a 25 × 11 tensor mesh (I_D at Vd 0.05 V, Vg 0.5 and 1 V). The generated (isotropic) 2D MOSFET was measured with both policies to 12 passes: about 85 k nodes, I_D still moving, 47 minutes, not converged (as Unit 14's runs); it stays a probe (`[.adapt_probe]`) and is a risk for AC5 (slice 7). **Metric:** recall and precision count node–axis pairs (a node refined along x and y is two); on generated meshes one isotropic axis, refined below 0.9 of the start size. **Measured 2026-11-12** (slice 5a, before ID4–ID6): recall 0.87 / 0.83 / 0.51 (≥ 0.8 fails on the MOSFET: the prediction refines across its features, the basic loop both axes; per node it was 0.94); precision 1 / 1 / 0.996; false-positive DOFs 0 / 0 / 102 of 26460; DOFs 1116 vs 540 and 504 vs 216 (≤ 1.1× basic fails on the diode and the MOS-C), 26460 vs 31752; passes 11 vs 10, 2 vs 2, 5 vs 6 (saved 0, fails); MOSFET 76 s vs 100 s. AC14 is marked `[!shouldfail]`, re-measured after slice 5b; owner ruling then. |

*Slice 5b — benefit, cost model, online fits and allocation (ID4–ID6, ID7's candidates; AC15, AC16, AC18).* The
factor fill nnz(L + U) is reported by every backend (`FactorizationReport::fill`: SparseLU's nnzL + nnzU, PARDISO's
iparm[17], ILUT's factors) and kept per Newton solve with the Jacobian's nonzeros (`ConvergenceRecord`, not saved or
observed). `solve/adapt_policy.hpp` adds `fit_cost_model`, `fit_orders`, `order_of` and the budgets; `adapt.cpp`
the candidates, their benefits, work and predicted seconds, the allocation, the realized ratios and the refits;
`AdaptPass` the per-pass DOFs, nonzeros, fill, fitted models, predicted and listed candidates.

| # | Question (ruled under the standing instruction "always choose recommended") | Options | Recommendation and ruling |
|---|---|---|---|
| AR7 | The constructions ID4–ID6 leave open, and a conflict between ID4 and AD9. | (a) as below; (b) rank by the timing-fitted seconds (not bit-identical, AD9 and AC18 fail by construction). | **(a)**, ruled as recommended; no criterion changed. **Determinism (owner to confirm):** ID4's cost is seconds from a timing-fitted model, so ranking by it would make the mesh depend on the machine's timings, against AD9 and AC10/AC18. The allocation ranks by a deterministic work in ID5's form at its default exponent, N_newton ((W + dW)^1.5 − W^1.5), W = nnz + fill (counted, not timed), dW the candidate's added nonzeros at the mesh's mean degree and the candidate's unknowns per node, with their fill at the current ratio; the timing-fitted model predicts each candidate's and each pass's seconds (reported, gated by AC15) and enforces only a time budget (wall-clock by nature). **Adjoint only:** ID4's benefit is "in the goal's units through the adjoint"; without weighting the estimate is the last change, which a large refinement raises (measured: the tensor MOSFET's proposals took it from 2.5 to 3.5 at 19 k nodes, then the fallback refined to 51 k nodes against basic's 10.6 k), so the allocation runs with `weighting = adjoint` only, and without it the layer is slice 5a's (ID2, ID3, AD3's marks). **Candidates:** an axis interval (tensor axes, the swept axis; one line added per halving) or a node (a generated 2D size field or a swept cross-section, its size halved; 3 nodes added in 2D, 3 per line swept). An interval can be split to three levels in one pass, each level a greedy entry of benefit E (S_c / S) 2^−(l−1)p (1 − 2^−p) and 2^(l−1) times the first level's nodes (without levels every candidate with an excess was taken while the target is unreachable, θ → 0: the tensor MOSFET then needed 2× basic's DOFs); a node one level (a size field is halved at most per build). S_c is the candidate's weighted excess (ω η: intervals over their edges along the axis, with ω the larger end's; nodes their own; a swept cross-section node over its lines), S the sum over candidates, E the relative estimate (AR5). **Target and caps:** take entries in decreasing benefit over work (ties by construction order; a level after the one before it) until the predicted estimate ≤ min(1/2, E/2) (tol/2 at most; at least halving, so a proposal is not a miss by construction), the DOF growth cap (2.5× per pass) or the budget. Both caps hold on the built mesh: grading adds lines and sizes the prediction misses (measured 1.5× on the well MOSFET), so a proposal built over a cap is allocated again with the measured inflation (built over predicted added DOFs), at most three times, and the inflation is kept for the next pass. Budget mode (`dof_budget`, `seconds_budget`): the same order until the budget is spent; a pass that fits nothing, or a built mesh over the budget, stops the run (`max_nodes`, unverified); the prediction's rounds (ID2) also stop at the budget. **Realized ratios (ID5):** a taken candidate's local level, the ω-weighted mean over its place (an interval's slab; a node's disc of its local size, in its region; a swept cross-section node's cylinder) of the node level max over edges of raw / threshold for its indicator kind (the excess vanishes below the threshold, the level does not), on the next pass over this one; per halving its l-th root. Measured for the first 200 taken per pass, the first 20 listed (ID7). **Fits:** the order per (kind, region), the mean of −log2 ratio with a ridge of one sample toward 1, clamped to [0.5, 3]; the cost model per Newton iteration a nnz^b + c fill + d DOF with mesh and transfer per node, one sample a alone at b = 1.5, two a and b, three or more b on a grid of 0.05 in [1, 2] with a, c, d non-negative least squares of the relative error and a ridge (0.01 / n)(b − 1.5)² (0.05 fixed outweighed the misfit on exact data); a pass's time predicted from the earlier passes with the last pass's Newton count at this mesh's measured nonzeros and fill. The node-only model T = αN^β (AC15's comparison) is fitted beside it. **AC16 devices:** AC5's run is the generated MOSFET's (AR6), so the tensor MOSFET (AC14's) stands in, adjoint-weighted; "the same verified accuracy" is both certified at the tolerance, the errors against the 8× uniform reference reported (4× and 8× differ by 2.8% at Vg 0.5 V: the reference is not converged to 1%, and both certified results are 3% from it; for AC20, slice 7). Basic at a budget is basic with `max_nodes` at the budget's nodes. **AC15 devices:** the tensor MOSFET for AC5's run, the adjoint diode for AC6's; AC8 (the FinFET) is slice 7's. **AC18:** "bit-identical, report included" excludes the measured seconds and what is fitted to them; the report's numbers are recomputed from the saved passes (orders, benefits and realized reductions, the predicted estimate, both fitted models and the predicted seconds). |

**Slice 5b measured (2026-11-13, Release, r3):**
- AC16 tensor MOSFET (adjoint, tol 1e-2, I_D at Vg 0.5 and 1 V): basic 8 passes, 15 732 DOFs, 23 s; intelligent 4
  passes, 9 270 DOFs (0.59×, ≤ 0.8 holds), 9 s; both certified, errors against 8× 3.1% and 3.4%. Budget 7 866 DOFs:
  intelligent 7 614 DOFs, error 2.8945%; basic 3 960 DOFs, 2.8969% (holds, barely). **Passes.**
- AC16 well MOSFET (AC7's deck, adjoint, I_ref = 5.20577): basic 6 passes, 6 642 DOFs, error 0.47%; intelligent 3
  passes, 12 744 DOFs (1.92×), error 0.32%, relative estimate 0.013 at the stop against basic's 0.52: the first
  allocation, growth-capped from E = 9.5, overshoots the accuracy. Budget 3 321 DOFs: intelligent 3 255 DOFs, 6.5%;
  basic 3 240 DOFs, 7.2% (holds). **Fails in tolerance mode**; `[!shouldfail]`, owner ruling.
- AC15: rank correlation 0.97 (MOSFET, 14 candidates) and 0.997 (diode, 37) (≥ 0.6 holds); the cost model from
  the third pass within 30% (worst 5% MOSFET, 15% diode; holds) but not better than node-only (mean 3.3% vs 2.7%,
  8.0% vs 7.6%: both small and within timing noise on 2 and 12 passes; fails); the adjoint diode's orders 1.10
  (potential, 32 samples), 1.09 (electrons, 2), 0.81 (doping, 4) and 1.64 (holes, 2; outside [0.8, 1.3], fails).
  `[!shouldfail]`, owner ruling. The adjoint diode does not converge in 14 passes with either policy (its estimate
  moves between 2 and 15 on meshes a few nodes apart; basic adjoint likewise: 571, 85, 105, 2.6, 25, …); that is the
  weighted-current gradient finding of slice 4, not the allocation.
- AC18: two runs bit-identical (meshes, goals, estimates, candidates, fits, certificates), the report recomputed
  from the saved passes. **Passes.**
- AC14 is unchanged by 5b (its runs are feature-only); with adjoint weighting the intelligent loop takes 4 passes
  against basic's 8 on the tensor MOSFET (AC16's run).

**Owner rulings after slice 5b (2026-11-13; owner, "Decisions for you: choose recommended").** No criterion is
loosened or re-worded. (1) AR7's determinism construction is confirmed: the allocation ranks by the counted work
N_newton ((W + dW)^1.5 − W^1.5), and the timing-fitted model only reports (AC15) and enforces the time budget.
(2) The allocation runs with adjoint weighting only; feature-only runs keep slice 5a's layer. (3) AC14 (AR6), AC15
and the well-MOSFET AC16 stay `[!shouldfail]` with their measurements above. They are re-measured once slice 7's gates
are in, and any still unmet at the unit's close is recorded as a limitation in ARCHITECTURE.md 11.3, as for Units
46, 48, 53 and 54. (4) AR4: AC7 keeps the recorded 4×/8×/16× reference 5.20577. Its uncertainty (about 0.4%) is stated
beside AC7 and re-examined with the tensor MOSFET's unconverged 8× reference under AC20 in slice 7.

**Slice 6 (AD7, AD8, AD13; AC9, AC10).** `[adapt]` in the deck (deck/build.hpp), the runner's adaptive path
(deck/run.cpp), `results::AdaptReport` (results/adapt.hpp) filled by `solve::adapt_report`, its file
(io/adapt_report.hpp, content kind 8, `.nitadapt`), the identity stamp (`solve::stamp_adapt`), the certificates in
`nitcad run`'s output and the CSV table, exit status 3, and the app's deck vocabulary. 6b adds `nitcad plot`'s and the
app's views of the passes.

| # | Question (ruled under the standing instruction "always choose recommended") | Options | Recommendation and ruling |
|---|---|---|---|
| AR8 | The constructions AD7, AD8 and AD13 leave open. | (a) as below; (b) implement the derived goals now. | **(a)**, ruled as recommended; no criterion changed. **Goals:** `goal = current <contact>` only; `vt`, `gm` and `swing` are refused with their line ("not adapted in this unit"). The loop, its adjoint (AR3) and AD11 are built on terminal currents, and a derived goal needs its own adjoint (a functional of several points) and resolution rule: a later unit, **confirmed by the owner (2026-11-13)**. `predictor` (17b) is refused. **Defaults:** `weighting = none`, `policy = basic` (ID9: intelligent is not the default until AC14–AC17 pass, and AC14 and AC15 do not), `verify = richardson`, `certify = report`, `tolerance = 1e-3`, `max_passes = 6`; `at = auto` (absent: the same): the sweep's first and last points, or, with the intelligent policy, up to `max_points` (default 3) chosen by ID3; `at = <indices>` counts the sweep's points from 0. `mode = budget` needs the intelligent policy and `budget_dofs` and/or `budget_seconds` (the basic loop takes `max_nodes`); the budget keys are refused in tolerance mode. **Scope refusals (AC10), each with its line:** an analysis other than `[sweep]` (AQ3: one mesh per sweep), a restart (the checkpoint is the deck's mesh's), a generate3d, file or unswept extrude3d mesh (AQ2), `weighting = adjoint` with `check_adjoint`'s refusals (a backend that cannot transpose and no fallback, electrothermal, energy balance, nonlocal tunnelling) or a non-ohmic goal, a goal on a gate or electrode, a second goal on one contact, the adaptation's outputs without `[adapt]`. **Rebuild (AQ1):** the builder is `deck::build(deck, refinement)`, the whole deck rebuilt per pass: a tensor axis the loop refined is taken whole (its straddles already in it; their spans protected), a generate2d mesh and a swept cross-section take the size field, a swept mesh the refined lines; an empty refinement gives the deck's mesh (AD10). **Outputs (AD8):** `sweep`, `device`, `mesh`, `checkpoint` and `table` are the final mesh's; `adapt = <file>` the report; `adapt_meshes = <directory>` each pass's device and mesh (`pass_<k>.nitdev`, `pass_<k>.nitmesh`, written as each pass is solved through `AdaptOptions::on_mesh`): the passes' meshes the app and `nitcad plot` show. The report holds the run record, the settings, per pass the mesh digest, source, goals, estimates, marks by indicator and region, Newton counts, nonzeros, fill, the decision layer's predictions, candidates, orders and cost model, the times; the certificates; the stop, fallback and verification; its status is the worst certificate's. **Identity (AD8):** the final sweep's `input_identity` continued by FNV-1a 64 over "adapt", each setting's name and value (bit pattern) and the final mesh's digest, the settings appended as `adapt.*` (`solve::adapt_settings`); a run without `[adapt]` is unchanged (AD10). **Certification (AD13):** `nitcad run` prints each pass and each certificate, with "an estimate and an independent Richardson check, not a proof"; the CSV table gains per goal `<c> certificate`, `<c> estimate A`, `<c> verified error A`, `<c> allowed A` and `<c> order` after `substeps` (filled at the adaptation points); with `certify = required` a result that is not certified exits **3** (2 stays an early stop), its outputs written. |

**Slice 6 finding (not changed):** the deck tests' `diode_device` fixture (tests/deck/fixtures.hpp) is not a diode. Its
profiles `acceptor constant 1e17 x=-1:1um` and `donor constant 1e17 x=1um:3um` window the depth axis x, which a profile
ignores (doping_profiles.hpp: "the depth axis's entry is ignored"), so both dopants cover the whole device and it is a
compensated resistor (the anode current is linear in the bias, 1.46 mA at 0.2 V). Its tests compare a deck against the
same hand-built device, so they still hold. Slice 6's deck uses a lateral window (`axis=y`) on a 21 × 2 mesh, which
gives a junction (2.6 nA at 0.4 V). Owner to decide whether to correct the fixture (every deck test that uses it would
change its numbers).

**Slice 6b (AD8: the views).** `nitcad plot --view adapt` draws each pass's estimate over its allowance against its
DOFs (log-log, the tolerance the line at 1), from the deck's `[output] adapt` or a `.nitadapt` file; `--view passes`
tiles each pass's mesh from `adapt_meshes` (`--columns`, 4; `--size` per tile, 640×480), each pass file checked against
the report's mesh digest. Every other view of a deck with `[adapt]` draws the **final** mesh, read from its `device`
(and `mesh`) output: the results were solved on it, not on the deck's start mesh, so `plot` (and the app) would
otherwise refuse them by digest. The pass files are named by `deck::pass_device_file` / `pass_mesh_file`. The app:
`Workspace` carries the report, the passes' directory and, for a pass shown alone, its index; `workspace_from_run`
takes the final mesh's device from the run, `workspace_from_plan` from the device output or the last pass's file
(whichever has the report's final digest; an error when neither does: the report's results cannot be paired);
`workspace_of_pass` loads one pass's mesh without results; the run record gains the adaptation (status and stop, "an
estimate and an independent check, not a proof", passes and nodes, policy and weighting, fallback, verification, each
certificate) and the pass shown. The results page lists the passes (pass, mesh source, nodes, DOFs, estimate over
allowance, marked, seconds, notes) for an adaptive run only; activating a row shows that pass's mesh in the device
view with the results kept aside, and the last pass brings them back (`MainWindow::results_show_pass`). Tests: the
plot views, their refusals and a changed pass file (tests/deck/adapt_test.cpp), the workspace (tests/app/
workspace_test.cpp) and the window (tests/app/analysis_window_test.cpp). Slice 6 is committed as one commit (6a and
6b), its suites run on both together.

**Slice 7 (gates and cost: AC5, AC8, AC11, AC20, AC12; ARCHITECTURE 6.12.4 and 11.3).**

| # | Question (ruled under the standing instruction "always choose recommended") | Options | Recommendation and ruling |
|---|---|---|---|
| AR9 | The tensor MOSFET's reference (AC16's comparison, AC20's on AC5's stand-in, AR6). The 8× uniform reference of slice 5b is not converged. | (a) per-axis ladders, each axis extrapolated, recorded with a 2× guard; (b) a finer uniform ladder (32×: 590 k nodes, hours); (c) keep 8×. | **(a)**, ruled as recommended; no criterion changed. **Measured (Release, 2026-11-13), I_D (A/cm) at Vg 0.5 and 1 V.** Uniform 1×–16× (275–61 985 nodes): 0.13051, 0.13257, 0.13214, 0.12846, 0.12575 and 0.28192, 0.30445, 0.31875, 0.31947, 0.31763: not monotone, because the y error (the inversion layer: 16× is 1.25 nm) and the x error have opposite signs and different orders. Per axis: at kx = 4, ky = 32, 64, 128: 0.122849, 0.122425, 0.122309 and 0.305497, 0.304497, 0.304174 (orders 1.87, 1.63); at ky = 32, kx = 4, 8, 16, 32 (to 246 849 nodes, 14 min): 0.122849, 0.124059, 0.124385, 0.124483 and 0.305497, 0.312421, 0.314917, 0.316039 (orders 1.89 then 1.73; 1.47 then 1.15). The y step 32 → 64 at kx = 8 is within 4% (0.5 V) and 12% (1 V) of that at kx = 4: the errors add. Each axis extrapolated at its last observed order, the y correction scaled to the finer x: **0.12393 (uncertain by about 0.1%) and 0.3153 (about 0.3%: the x order still drifts)**. The 8× values are 3.7% and 1.3% high, which is slice 5b's "certified runs 3% off". The values are recorded in adapt_intelligent_test.cpp (`mosfet_ref`), with the 2× solve (0.1325677205, 0.3044490795) as the guard that the device and the physics are unchanged (as AR4). AC16 now reports its errors against them instead of solving 8× and 4× each run. |

**AC11 (measured 2026-11-13, Release, another 246 k-node solve running beside it; timings re-measured idle at
close).** The estimate's factorizations of J′ on M′ and its largest transposed solve are timed (`PointAdjoint`,
`GoalEstimate`, `AdaptPass::seconds_estimate_factorize` and `seconds_transposed`, the report row and `.nitadapt`; both
excluded from the AD9 comparison as measured seconds). The gate compares one transposed solve against one
factorization, a lower bound on a Newton iteration's linear solve (which also solves). Tensor MOSFET, adjoint: basic 8
passes (275 → 5 244 nodes; nnz 7 615 → 150 476; fill 49 890 → 2 671 164; 18–22 Newton iterations per pass), the
transposed solve 2.2–5.4% of a factorization (0.047 s against 1.32 s on the last pass); verification 21 436 nodes,
14.9 s against the final pass' 1.6 s. Intelligent 4 passes (663 → 3 090 nodes), the transposed solve 2.5–5.2%;
verification 12 900 nodes, 5.6 s. Adjoint diode, 14 passes (21 → 143 nodes), 3–13%; not verified (the max_passes stop
of slice 5b). **Passes** (the AC5 and AC8 rows with their own runs below).

**AC20 on AC5's stand-in (measured 2026-11-13; AR9's reference).** Tolerance 1e-2. Basic: certified at both points,
actual errors 0.156% and 0.431%; the verified error 0.33 and 0.87 times the actual (p 1.86 and 0.5), the adjoint
estimate 1.65 and 1.17 times. Intelligent: certified, 0.452% and 0.198%; verified 1.88 and **0.038** times the actual
(p 0.5 and 2.18), estimate 1.59 and 0.38 times. Each run stopped by a `max_nodes` budget one node above its first pass
is unverified at both points (errors 5.3–11%). **Passes.** Finding (not a failure, recorded): AD11's verified error
can underestimate the actual error by a large factor when its observed order is high (2.18 here, where the order
clamp [0.5, 3] lets it extrapolate the last change down by 2^2.18 − 1 = 3.5 and the change itself is small); the
certificate held only because the actual error was a fifth of the allowance. AC7's case now also checks that a
certified result is within its tolerance of AR4's reference.

| # | Question (ruled under the standing instruction "always choose recommended") | Options | Recommendation and ruling |
|---|---|---|---|
| AR10 | AC5's reference and its "uniform-refinement ladder needs" (AC20's "AC5's 4× uniform"). | (a) AR9's converged reference (the same continuum device), the 4× generated value reported beside it; the ladder's nodes bounded below from its h0/4 and h0/8 runs at an assumed second order; (b) as frozen: the generated mesh's own 4× reference, the ladder measured to 1%. | **(a)**, ruled as recommended; **confirmed by the owner (2026-11-13)** (it changes the reference AC5 and AC20 name, toward a more accurate one, not a looser one). The generated MOSFET is the tensor MOSFET's device (geometry, doping, contacts and gate identical; only the mesh differs), so its mesh-converged I_D is AR9's. Its 4× uniform value (5 nm, 20 295 nodes) is −13.8% and −19.3% from it (8×, 80 660 nodes: −6.1%, −11.0%), so as a reference it would certify nothing meaningful at 1%, and the ladder cannot be measured to 1% (it converges below first order between 4× and 8×; second order from 8× on would still need about 890 k nodes; the observed order about 1 M × 30). Vg 0 (above its resolution, 4.8e-9 A/cm) gets the same per-axis construction: x at ky = 32 (kx 4, 8, 16): 7.08387, 7.08502, 7.08562e-5 (order 0.94); y at kx = 4 (ky 32, 64, 128): 7.08387, 7.03040, 7.01701e-5 (order 2.0), the y step at kx = 8 equal to 0.1%: **7.0150e-5 A/cm** (about 0.05%). |

**AC5 (measured 2026-11-13, Release; generated MOSFET from h0 = 20 nm, 1 323 nodes; I_D at Vg 0, 0.5 and 1 V,
tolerance 1e-2, adjoint-weighted; against AR9 and AR10).** Uniform: h0/2 5 147 nodes (−23.6%, −27.4% at 0.5 and 1 V),
h0/4 20 295 (+12.9% at 0 V, −13.8%, −19.3%), h0/8 80 660 (+2.7%, −6.1%, −11.0%). **Basic** (12 passes, 17 s): 1 323 →
2 888 nodes, errors −15.3% and −19.9% at the last pass (θ = 0.3 marks a few nodes; growth about 10% a pass), relative
estimate 22.9, max_passes, unverified. **Intelligent** (9 passes, 40 min): the prediction's 3 303 nodes, then
proposals 7 202, 14 823, 30 104, 68 898 (+0.35%, −0.36%, −1.33%), **98 593 (+0.37%, −0.20%, −0.89%: inside 1% at every
point)**, 120 824, 129 977, 146 084 (+0.30%, −0.19%, −0.83%); relative estimate 1.05 at the stop (the estimate 1.25 times
the actual error at 1 V), so it stopped at `max_nodes` (150 000), unverified at every point. **98 593 nodes against at
least 890 k for the ladder (0.11 of it; ≤ 1/3 holds) with the intelligent policy; the basic loop does not reach 1%
within 12 passes.** The gate (`[.ac5]`, manual like `[.uc8]`: about 40 minutes) repeats the h0/4 and h0/8 runs for the
bound and the intelligent run to 100 000 nodes, and checks AC20 on it.

**AC8 (measured 2026-11-13, Release, PARDISO on 6 threads; `[.ac8]` in tests/deck/complex_finfet_test.cpp).** The
start: decks/finfet_3d_swept.deck with its axis lines about 10 nm apart (0, 20, 25, 35–50 every 5, 60, 65, 85 nm, the
gate edges straddled as the deck) and its cross-section unstructured at up to 12 nm (the structured frame unchanged):
4 718 nodes. `[adapt]`: goal the drain current, tolerance 1e-2, adjoint-weighted, `at = auto` (0 and 0.8 V for basic,
ID3's points for intelligent), max_nodes 9 367 (UC8's budget). The coarse tensor mesh (9 367) and the 2× reference
(45 375) re-solved as UC8. **Basic:** 2 passes (4 718 → 7 678 nodes; the next would exceed the budget), relative
estimate 175 then 139, stop max_nodes, unverified; UC8's rule at **9 of 17** points (fails at 0.1, 0.15 and 0.35–0.6
V: on-current 3–5% low at 0.4–0.55 V, the tensor within 0.6%), V_T and the swing within their bounds, max g_m −8.3%
outside. **Intelligent:** the prediction's 7 350 nodes and no further pass, relative 125, unverified; **6 of 17**,
V_T and max g_m outside, the swing within. UC8's hand-tuned deck (Unit 54, 8 300 nodes) met 12 of 17. **The measured
target is missed** (it was stated as a target, not a gate). The gated parts hold: every point converges, the budget
holds, both runs unverified (AC20). The cause measured: on the swept 3D mesh the adjoint estimate is 84 times the
actual error at 0 V (2.98e-10 against 3.56e-12 A) and 9 times at 0.8 V (1.29e-6 against 1.37e-7 A), so the loop sees
itself 100× from its tolerance and the budget allows one or two passes, spent where the overestimated weights point.
AC11 on AC8: per pass 6 734 and 11 550 DOFs, 181 820 and 302 032 nonzeros, fill 6.46 M and 13.1 M, 109–115 Newton
iterations over the 17 points; solve 5.8 and 10.1 s, the estimate 3.1 and 5.2 s (factorizations 1.0 and 2.1 s,
transposed solve 0.024 and 0.041 s: 2% of a factorization); no verification run (budget stops).

**Re-measured after slice 7's gates (2026-11-13, Release, idle machine; the owner ruling after slice 5b).**
`[adapt_intelligent],[adapt_goal]`: AC7 passes (adjoint 2 214 nodes against feature-only 6 834, 0.32; error 0.47%
against AR4). AC11 and AC20 pass. **AC14** unchanged (its runs are feature-only): recall 0.87 / 0.83 / 0.51, precision
1 / 1 / 0.996, DOFs 1116 vs 540 and 504 vs 216, passes saved 0: **unmet**. **AC15**: rank correlation 0.97 and 0.997
(holds); the diode's orders 1.10, 1.09, 0.81 and **1.64** (holes, 2 samples), and the cost model not better than
node-only within timing noise: **unmet**. **AC16** tolerance mode on the tensor MOSFET: 9 270 against 15 732 DOFs
(0.59), both certified, errors against AR9 0.45% and 0.43%: **holds**. **AC16 budget mode on the tensor MOSFET, against
AR9's reference: intelligent 0.6531% (7 614 DOFs) against basic's 0.6507% (3 960 DOFs) at a 7 866-DOF budget —
unmet** (a tie within AR9's 0.1% uncertainty; slice 5b's pass was against the unconverged 8× reference, 2.8945%
against 2.8969%). The case is split: tolerance mode stays a gate; budget mode is `[!shouldfail]` with this
measurement, the criterion unchanged. **AC16 on the well MOSFET** tolerance mode 12 744 against 6 642 DOFs (1.92×):
**unmet**; its budget mode holds (6.5% against 7.2%). Per the ruling, the unmet criteria at the unit's close are
limitations in ARCHITECTURE.md 11.3: AC14, AC15, AC16 (well tolerance mode; tensor budget mode), and AC8's target.

**Slice 7 suites (2026-11-13, at ca22846):** Release 217/217 (980 s, `-j 6`); Debug 198/199, the one failure
`app.live_check` (its wake count, line 53: the known flake under `-j 6`, untouched by this unit), which passes alone
three times in a row (`--repeat until-fail:3`).

**AC12 (hand mutants, 2026-11-13).** Each mutant one exact-string change, the solve, deck and mesh libraries and
`nitcad_adapt_test`, `nitcad_deck_test` and `nitcad_mesh_test` rebuilt (Release, r3), then the fast suites: every
non-hidden adapt case, the deck's adapt cases and `[size_field]`; killed if a build or a test fails. **Round 1, 37
mutants:** 23 killed; 13 survived; 1 did not build (ID5's, an unused variable under `/WX`; replaced). The survivors were
rules the end-to-end runs never pinned: AD1 (solve from scratch, not the transferred state), AD4 (halving 0.5 → 0.7;
the `min_size` skip dropped), AD9 (the tie-break reversed), AD11 (the default order 1 → 3), AD12 (the miss count not
reset), AD13 (`&&` → `||` in the status; the allowance doubled), ID1 (the margin's sign), ID3 (one point, not
`max_points`), ID6 (the target E, not E/2; the growth cap ×10), ID7 (2 listed, not 20). **Closed by tests:** the rules
are functions in solve/adapt.hpp, applied by `adapt()` in place of the inline code (no behaviour change; the fast
suites' numbers unchanged): `allowed_error`, `richardson_error`, `certificate_status`, `refined_size`, `heavier_node`,
`proposal_missed`, `MissGuard`, `allocation_target`, each tested exactly in adapt_loop_test.cpp ("the loop's
rules"), with `CostModel::per_iteration`'s form; AC6 checks that every pass after the first was transferred onto
(AD1); AC17 that `max_points` points are chosen (ID3) and that the prediction's rounds and every proposal hold the
growth cap (ID6); AC18 that the first 20 taken are listed (ID7) and that an allocation away from the cap stops at
its target (ID6; on that diode E never falls below 1, where the target's two forms differ, hence the exact test).
**Round 2:** the 13 re-aimed at the rules' new places, and ID5's (the DOF term dropped): 12 killed; ID5's and ID6's
target survived; **round 3** after their tests: both killed. **Every design row has at least one killed mutant and
none survives:** AD1 (2), AD2 (2), AD3 (2), AD4 (3), AD5, AD6 (3), AD7, AD8 (2), AD9, AD10, AD11 (4: AC20's skipped
refinement, wrong p twice, wrong mesh compared), AD12 (3), AD13 (3: AC20's status logic twice, the allowance), ID1,
ID2, ID3, ID4, ID5 (2), ID6 (2), ID7, ID9. **AC12 met; AC20's mutant clause met.** The runner first restored each file
with `git checkout`, which discarded uncommitted edits to the mutated file in round 2's first attempt; it now restores
from a copy, and each round runs on a committed tree.


---

