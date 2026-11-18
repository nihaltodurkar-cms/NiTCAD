<!-- Archived from DECISIONS.md on 2026-11-17 (branch docs/archive-17-26d), word for word; DECISIONS.md keeps its summary. -->

### Unit 17b — cross-run learning for adaptive refinement (started 2026-11-13; on `mesh/cross-run-learning`)

**Start (owner, 2026-11-13: "start Unit 17b, cross-run learning").** The design is Unit 17's frozen rows ID8, ID9
(17b's part), ID10, AC12 (ID8, ID10), AC19 and AC21's predictor case, with amendment OA5; nothing there is changed
here. What follows are the constructions those rows leave open, each put as a question with a recommendation and
ruled as recommended under the owner's standing instruction (2026-11-11, "always choose recommended"). None loosens a
criterion; LB10 narrows what one design row stores and was confirmed by the owner.

| # | Open construction | Options | Ruling |
|---|---|---|---|
| LB1 | The base. The split says "a separate branch and unit after 17 merges"; PR #81 (Unit 17) is open. | (a) `mesh/cross-run-learning` from `mesh/adaptive-refinement` at 6ffad0b (stacked); its PR to main once #81 merges (or retargeted then); (b) wait for #81. | **(a)**, ruled as recommended. Unit 17 is accepted and pushed; 17b reads its code, not main's. Any later change to #81 is merged into this branch. |
| LB2 | The sample and its label ("node features → realized benefit", ID8). The predictor's output is a first mesh, i.e. sizes. | (a) **One sample per edge** of every solved pass's mesh; its label the realized refinement ℓ = max(0, log2(L_e / h_f)), L_e the edge's length and h_f the final mesh's local size along the edge's direction at its midpoint (the extent, projected on the edge's unit direction, of the final cell containing the midpoint, by Unit 17's `PointLocator`); an edge whose midpoint lies outside every final cell is skipped. The realized benefit is the refinement the converged adaptation found worth its nodes there. (b) Per refined candidate, the adjoint estimate's realized reduction (ID5's samples): defined only where something was refined and only with adjoint weighting, and not a size. | **(a)**, ruled as recommended. ℓ is defined on every edge, with any weighting, and maps straight to a size (LB4); the per-candidate reductions stay ID5's online fit. "Benefit prediction" in AC19 and ID10 is the prediction of ℓ. |
| LB3 | The features (ID8's list) and the schema. | (a) **Schema 1**, per edge from its end nodes (the mean of the two, or the edge's own): region kind (semiconductor, insulator: a category); material (its name's index in the model's sorted material names: a category); log10 max(\|N_D − N_A\|, n_i); \|∇ log10 max(\|N\|, n_i)\| along the edge (per µm); log10 of the distance to the nearest junction point, insulator-interface point and contact node (cm; 1 cm when there is none); log10 of ID2's W and L_D of the nearest junction point and of λ of the nearest interface point (`predict_sizes`' construction; 1 cm when none); the analysis' largest junction and gate bias differences (V, ID2's); the goal kind (0: a terminal current, AR8); log10 L_e; \|cos\| between the edge and the nearest feature's normal; the nearest junction's distance over its W and the nearest interface's over its λ. 17 features. (b) ID8's list alone. | **(a)**, ruled as recommended. The additions are what the label depends on that ID8 did not list: the edge's own length and direction (an edge across a junction and one along it are refined differently, Unit 17's projected rule), and the two scaled distances (a device family varies L_D and λ; trees cannot form ratios). A change to any feature's definition is a new schema version (ID10). |
| LB4 | From a prediction to the first mesh. | (a) Per edge the target L_e 2^−ℓ̂ where ℓ̂ ≥ 0.5 (none below); per node the smallest over its edges, with that edge's direction as the normal (Unit 17's `SizeTargets`: projected on tensor and swept axes, isotropic on a generated 2D mesh); used for **the start mesh only**, in ID2's prediction rounds (at most 12, each on the refined mesh, its edges predicted afresh), ID1's guard and AD12's fallback unchanged; every later pass is Unit 17's online layer. (b) Also every later proposal (Unit 17's `prediction` hook). | **(a)**, ruled as recommended: ID8 says the model "predicts the first mesh" and replaces ID2's a priori field; ID9 measures it against 17's online layer, which (b) would remove. A new `PolicyOptions::first_prediction` hook carries it; the existing `prediction` (AC17's test hook) is unchanged. |
| LB5 | Model constants (LQ2's "fixed hyperparameters"). | (a) **Trees:** squared loss, 200 trees, depth 4, shrinkage 0.1, at least 20 samples per leaf, the start the label mean; split candidates at up to 64 cut points per feature (the training values sorted, cut points at equal ranks, duplicates removed); the best split by the largest decrease in squared error, ties to the lower feature index, then the lower cut. Categories as their index. **Ridge (the ID9 baseline):** the continuous features standardized by the training mean and deviation, categories one-hot, the intercept unpenalized, λ = 1 in standardized units, solved by Cholesky. (b) Tuned per training set. | **(a)**, ruled as recommended: fixed, deterministic, single-threaded training; tuning needs a validation split the families are too small for. |
| LB6 | Files and the manifest. | (a) **Adaptation history** `.nithist` and **model** `.nitlearn` in io's container (content kinds 9 and 10; versioned, checksummed, atomic; doubles by bit pattern); the **manifest** a text file (`nitcad-learn-manifest 1`, the schema version, then one line per history: its file digest, FNV-1a 64 of the whole file, and its path), lines sorted by digest, duplicates refused; the training set is the histories' samples in manifest order, within a history in pass then edge order. At most 20 000 samples per history, by a fixed stride (every k-th, k = ⌈n / 20 000⌉), so one large run does not swamp the family. (b) The manifest in the container too. | **(a)**, ruled as recommended: a manifest is written and read by people choosing devices; its digests make it exact. The stride is deterministic, so "the same manifest always gives the same samples" holds. |
| LB7 | The held-out split (ID10: by device). | (a) One history file is one device. `nitcad learn train <manifest> <model> [--holdout <history>...]`: the named histories held out; without `--holdout`, every fourth history in manifest order, starting with the fourth (at least one). Refused (`invalid_input`): no training or no held-out device, a held-out path not in the manifest. (b) A random split with a seed. | **(a)**, ruled as recommended (deterministic; AC19 names its 2 held-out devices). The model is trained on the training devices only and the held-out metrics stored describe that model (it is not retrained on all). |
| LB8 | "Out of the training range" (ID8). | (a) The model stores each continuous feature's training minimum and maximum and the categories seen; an edge is out of range when a continuous feature is outside [min − 5% of the range, max + 5% of it] or a category is unseen (a material name not in the model's list); **a device is out of range when more than 5% of its start mesh's edges are**: the run then uses ID2's a priori field, says so in the report and on the console, and the identity still carries the model (the run asked for it). (b) Per-edge fallback. | **(a)**, ruled as recommended: a per-edge mixture of two predictors would be neither; the device is the unit of validation (ID10). |
| LB9 | The deck. | (a) `[adapt] history = <file>` saves the adaptation history (opt-in, any policy); `predictor = <file>` loads a model (needs `policy = intelligent`, refused otherwise with its line), its digest and versions entering the identity (`adapt.predictor` appended, as AD8's settings); absent: Unit 17's runs bit-identical (AC19). (b) `history` only through the CLI. | **(a)**, ruled as recommended: the deck is the run's whole input. |
| LB10 | ID10's held-out metrics "rank correlation, recall/precision of predicted refinement, passes saved". Passes saved needs adaptive runs on the held-out devices, which training cannot do without their builders. | (a) The file stores the held-out Spearman rank correlation of ℓ̂ against ℓ (trees and ridge), and recall and precision of ℓ̂ ≥ 0.5 against ℓ ≥ 0.5 on the held-out devices' first-pass edges; passes saved is measured by AC19's runs and recorded here, not stored. (b) `nitcad learn train` runs the held-out decks' adaptations. | **(a)**, ruled as recommended; **confirmed by the owner (2026-11-13: "LB10 confirmed, choose recommended")** (it narrows what ID10 says the file stores; no criterion changes: AC19 measures passes saved either way). (b) would make training a simulation campaign. |
| LB11 | Where the history is gathered. | (a) The runner (`deck::run`), through Unit 17's `on_mesh` callback: each solved pass's features computed as it is solved, the labels at the end against the final mesh; the history holds the deck text, the run record (Unit 42 provenance), the schema version, per pass the device and mesh digests, the features and labels. (b) Inside `solve::adapt`. | **(a)**, ruled as recommended: the solve layer stays free of files and deck text (Unit 17's division). |
| LB12 | AC19's family and gate cost. | (a) 8 tensor MOSFETs as deck texts (AR6's tensor MOSFET with L_g, t_ox and N_A varied), `policy = intelligent`, adjoint-weighted, tolerance 1e-2, I_D at V_g 0.5 and 1 V; a manual Release gate `[.ac19]` like `[.ac5]`. The determinism, versioning, identity, fallback and bit-identity clauses are fast cases on small devices. | **(a)**, ruled as recommended. |

**Slices.**
1. Features and labels (LB2, LB3), the history file and its io (LB6), the deck's `history` (LB9, LB11), the manifest
   and `nitcad learn manifest` (LB6); AC19's "the training set built twice from the same manifest is bit-identical".
2. Trees and ridge (LB5), training with the device split (LB7), the `.nitlearn` file with its versions, ranges and
   held-out metrics (ID10, LB8, LB10), `nitcad learn train`; AC19's model bit-identity and version refusal.
3. The predictor in the loop (LB4, LB8), the deck's `predictor` and the identity (LB9), the report; AC19's family gate,
   out-of-range fallback and bit-identity without `predictor`; AC21's predictor case; AC12's ID8 and ID10 mutants;
   ARCHITECTURE 6.12.4 and 11.3.


**Amendment at build (LB9, 2026-11-13).** The history is an output, so it is written by `[output] adapt_history =
<file>` beside `adapt` and `adapt_meshes` (AD8's convention), not by an `[adapt] history` key; `predictor` stays an
`[adapt]` key (it is an input). Nothing else in LB9 changes.

**Implementation record.**

*Slice 1 — features, labels, histories, manifests (LB2, LB3, LB6, LB9, LB11).* `solve/adapt_policy.hpp`:
`feature_points` (ID2's junction, interface and depletion points with their spans and sizes, and the analysis's
largest junction and gate bias differences), split out of `predict_sizes` without a change to it (the fast adapt suite
and `[adapt_intelligent]` unchanged: 28 cases pass; 3 pass and 3 fail as expected). `results/learn.hpp`: schema 1's
17 features, `HistoryPass`, `AdaptHistory`, `Manifest`, `TrainingSet`. `solve/adapt_learn.hpp`: `material_key` (a
category from the material's kind, eps_r and, for a semiconductor, Eg0, Nc300 and Nv300 as exact hexadecimal
floats: materials carry no names, and two decks naming one material differently give one key), `device_digest`,
`edge_features` (nearest junction, interface and contact points by an x-sorted scan, ties to the lower index),
`realized_levels`, `HistoryRecorder`. Feature 10 is ID2's per-point largest junction difference (the diode's points
−1 to 0.6 V give 1 V, not the range 1.6 V). `io/learn_files.hpp`: the `.nithist` file (content kind 9), `file_digest`,
`make_manifest`, the text manifest, `build_training_set` with LB6's stride. `deck::Deck` keeps its text (the history
records it); `[output] adapt_history`; `deck::run` gathers the history through `on_mesh` (after `adapt_meshes`'
writer when both are given). `nitcad learn manifest <manifest> <history>...`. Tests: `[adapt_learn]` (features on
the diode, labels against a refined mesh and against itself, a run's history saved, loaded and saved again
bit-identical, the manifest the same in either input order, the training set built twice equal, the cap, and the
refusals: a duplicate, a digest that does not match, another schema, unsorted lines); the deck case (the history's
deck text, run identity and per-pass mesh digests, the CLI's manifest, its usage and a non-history refused).

**LB13 (found at build, slice 2; ruled as recommended).** ID10 says the model records "the build that trained it
(Unit 42 provenance)", and AC19 that "the trained model" built twice is bit-identical. A Unit 42 provenance holds
the run's start time and machine too, so two trainings would differ. (a) **The model keeps the build's fields only**
(version, revision, dirty, build type, compiler, Eigen version); the run's (start, wall time, OS, CPU, MKL, cache
sizes) are left empty. (b) Keep them and compare models without them. **(a)**, ruled as recommended: the build is
what ID10 names, and the file itself is then reproducible.

*Slice 2 — the model (LB5, LB7, LB10, LB13; ID9, ID10).* `solve/adapt_learn.hpp`: `cut_points`, `fit_boosted_trees`
(histogram split search over LB5's cut points, one byte per bin; the left child grown first; leaves the mean
residual times the shrinkage), `predict_trees`, `fit_ridge` (its 20-odd columns solved by an in-place Cholesky:
Eigen is private to linalg), `rank_correlation` (Spearman, average ranks), `default_held_out`, `train_predictor`
(trees and ridge on the training devices only, both measured on the held-out devices; ranges and categories over
the training samples; the manifest digest over the entries' digests, not their paths), `predict_levels`.
`results::LearnedModel`, `TreeNode`, `HeldOutMetrics`; the `.nitlearn` file (content kind 10) refusing another
format or schema with both versions named, and a tree whose children do not follow their parent or whose feature is
out of range as corrupt. `nitcad learn train <manifest> <model> [--holdout <history>]...` (a held-out path matched
as written or as the same file), printing the metrics and the importances. Tests: the parts exactly (cut points,
rank correlation with ties, the default split, trees learning a step that ridge cannot, the first split, a tie to the
lower feature, `min_leaf`, ridge recovering a linear law), training by device (sides, counts, trees beating ridge on
a step, recall and precision, trained twice equal and its files bit-identical, loaded equal and predicting the
same, another format and another schema refused by name, another content kind, the split's four refusals, a
held-out device's labels not reaching the trees), the diode's histories end to end, and the CLI (training with a
named held-out history, an unknown one, usage).


*Slice 3a — the predictor in the loop (LB4, LB8, LB9; AC19's fast clauses, AC21's predictor case).*
`PolicyOptions::first_prediction`: the start mesh's prediction rounds take it in place of ID2's field (before the test
hook `prediction`, which still drives the later proposals when set: AC21's fixture uses both); every later pass is the
layer's own. `solve::learned_targets` (LB4), `out_of_range` and `out_of_range_fraction` (LB8), `stamp_predictor` (the
identity continued over the model file's digest, format, schema and use; settings `adapt.predictor.digest_high`,
`_low`, `.format`, `.schema`, `.used`). `[adapt] predictor = <file>` is loaded at build, so `nitcad check` refuses a
model of another format or schema with its line, and needs `policy = intelligent`. `deck::run` computes the start
mesh's (the deck's own, AD10) features against the model's materials; above 5% out of range it leaves ID2's field and
records `out_of_range`. `results::AdaptReport` gains the predictor's path, digest, status and fraction (io keys 25–28,
written only with a predictor, so Unit 17's reports are unchanged and still load); `nitcad run` and `check` print
them. Tests: the range rule and the targets exactly (the 5% margin both ways, unseen and absent materials, NaN, the
5% device limit, l^ = 1, 0.4 and the 0.5 threshold); the deck case (`adapt_history` alone leaves the identity and the
deterministic report bit-identical to a run without it; a model used, its digest in the report and its file; two
models two identities; a device of another material 100% out of range, its passes' meshes those of the run without
the predictor, the identity still carrying the model with `used = 0`; the console lines; another format refused at
build); AC21 with a learned model given through `learned_targets`, the certificates bit-identical to basic's.

**AC19 (measured 2026-11-13, Release; `[.ac19]` in tests/deck/cross_run_learning_test.cpp, 18.7 min).** The family
(LB12): tensor MOSFETs from 50 nm × 20 nm uniform starts, source and drain 300 nm, I_D at V_g 0.5 and 1 V, V_d
0.05 V, intelligent, adjoint-weighted, tolerance 1e-2, `max_passes` 10. **ID2's a priori first mesh** (each run's
history saved): m0 (L_g 400 nm, t_ox 4 nm, N_A 1e17) 10 passes, 144 144 DOFs (its stops rejected twice, the basic
loop then grew it to 48 048 nodes, and the verification ran on 198 830: 14.5 min of the gate); m1 (500, 5, 2e17) 5
passes, 21 420; m2 (600, 5, 1e17) 5, 19 026; m3 (700, 6, 5e16) 5, 13 542; m4 (800, 4, 2e17) 4, 12 198; m5 (500, 3,
5e16) 4, 15 867; held out: m6 (700, 5, 1.5e17) 5, 10 476; m7 (600, 4, 7e16) 5, 15 561; every run certified, stopped
converged. **The model** (trained on m0–m5, 94 244 samples after LB6's stride; 21 637 held-out samples): held-out
rank correlation **trees 0.399, ridge 0.161** (the trees beat the ID9 baseline); first-pass recall 0.900, precision
0.722. **Learned first mesh on the held-out devices:** m6 4 passes, **12 000 DOFs** (against 5 and 10 476); m7 4
passes, 10 260 DOFs (against 5 and 15 561); both certified. **AC19's family clause is met on m7 and unmet on m6** (a
pass saved, but 1.15× the DOFs); the gate is `[!shouldfail]` with this measurement. **The other AC19 clauses hold**
(fast cases): the training set and the model built twice bit-identical; another format or schema refused; the
identity changes with the model's digest; an out-of-range device flagged and on ID2's field; without `predictor`,
runs bit-identical. **Ruled under the standing instruction, as Unit 17's unmet criteria were (2026-11-13 ruling after
slice 5b):** the m6 DOF miss is a limitation in ARCHITECTURE.md 11.3, the criterion unchanged. Reading: the learned
field covers the refinement the converged runs made (recall 0.90) and so starts nearer the final mesh (a pass fewer on
both), but at precision 0.72 it also refines where m6's own run did not, and m6's a priori run was already the
family's cheapest.

**AC12 for 17b (hand mutants, 2026-11-13).** C:\Users\disha\mut17\run17b.sh in a separate Release tree (r3): one
exact-string change each; `nitcad`, `nitcad_adapt_test` and `nitcad_deck_test` rebuilt; killed if a build or the
`[adapt_learn]`, `[adapt_policy]` or "deck adapt*" cases fail. **Round 1, 22 mutants, 18 killed.** Survivors: ID8's
stride (2 in place of 1: every capped case had a stride of at least 2), LB5's shrinkage dropped, LB5's `min_leaf` at the
split (the node gate alone held the step test), LB4's predictor on every proposal (no case looked at the later
passes). **Closed by tests:** the uncapped set holds every labelled edge and the capped one exactly LB6's stride; a
one-split tree's leaves are 0.1 of their sides' mean residual; a step at 0.9 cannot be cut with 150 a side; with a
predictor and weighting none every pass after the first is the layer's own (marks). **Round 2: all four killed.**
Every 17b row has a killed mutant and none survives: ID8 (7), LB5 (5), ID10 (5), LB8 (3), LB4 (2), LB9 (1). **AC12 met
for 17b.** (A process note: commit 7359b44 captured round 1's first mutant, applied in the working tree while the
runner ran; 4adba48 restores it, and no commit is made while the runner is active.)

**AC13 for 17b (suites, 2026-11-13, at 5cd7008).** Release: 219/219 passed (868 s). Debug: 201/201 passed (980 s).
The AC19 gate is `[.ac19]` (hidden, Release by hand); the 17b cases in the suites are `[adapt_learn]`, the AC21
predictor case and the deck cases. Unit 17 merged to `main` meanwhile (72377ad, PR #81, owner "merge pr 81"); this
branch holds 17's head (6ffad0b), so its PR to `main` carries the 17b commits only.

