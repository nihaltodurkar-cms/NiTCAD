<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 26b — the analysis workspace (started 2026-10-22; on `app/analysis-workspace`)

**Owner request:** VQ10 (live deck checking and completion, structure preview, cut lines with band diagrams, run
comparison and difference maps, mesh quality; plan page Figures 4, 7, 8 and view 25), approved after Unit 26's
acceptance ("start 26b", 2026-10-22). The owner's conditions of Units 25 and 26 hold: numerics and layer boundaries
preserved, Release and Debug suites at every commit, `/W4 /WX`, hand mutants, nothing pushed or merged without
authorization.

**What exists (audited on `main` at `e9b1d71`).** The deck's reference list of sections and keys is
`deck/build.hpp`'s comment; `deck::build` stops at the first error, `"<source>:<line>: <what>"`, and rejects an
unknown key as `unknown key '<k>' in [<kind>]`. `app::check_deck` builds on a worker (`RunController`), which is busy
while a run goes on. `render::Locator` finds the element holding a point with barycentric weights; band edges are
stored per node in every result (NaN on insulator nodes); `band_plot` and `profile_plot` take a 1D mesh only.
`Workspace` pairs results with a device by the mesh digest; `analysis` extracts the threshold. `render`'s
`Geometry` gives the linear elements of every mesh (Kuhn tetrahedra of tensor cells). Coupling areas are positive by
construction (`mesh::Mesh`, P3), so mesh quality is a matter of element shape, not of negative couplings.

**Design (taken under VQ8–VQ14's delegation, "choose what's best for Commercial TCAD"):**

| # | Choice |
|---|---|
| XD1 | **Live checking.** A `LiveChecker` with its own worker, separate from the run controller (editing goes on during a run): a check starts 400 ms after the last edit; each check carries a generation, and only the newest check's result is applied (a build cannot be cancelled mid-way, so an older one is discarded when it ends). The editor shows the error's line in the gutter and its message under the editor, as Check does. |
| XD2 | **Suggestions.** An unknown section, key or value choice gets "did you mean ..." from the vocabulary (XD3): the nearest name by edit distance, at most 2 edits or a third of the name's length. |
| XD3 | **Completion.** A deck vocabulary in `app` (`deck_vocabulary.hpp`): every section kind, its keys (with or without default), the value choices of keys that have them, and a line of help; mirroring `deck/build.hpp`'s comment, which a test reads, so the two cannot drift. Region parameter paths (`eps_r`, `electron_mobility.mu_max`, ...) are open-ended and not completed. `complete(text, line, column)` gives what fits the caret: section kinds after `[`, keys at a line's start in a section, value choices after `=`. The editor opens the list on Ctrl+Space and when a key or `[` is typed; Enter or Tab inserts. |
| XD4 | **Structure preview.** After a live check builds, the device view shows the deck's device (regions, contacts and mesh) without solving. Results stay shown while the deck's mesh digest equals theirs; otherwise the view shows the structure, marked "preview". |
| XD5 | **Cut lines.** `render::cut_line(geometry, a, b)`: the segment clipped by every element it crosses, so the samples are the crossings with element faces; a nodal field read there by the element's barycentric weights is exactly the piecewise-linear field along the line (kinks at the crossings only). Outside the mesh is a gap. `cut_plot` and `cut_band_plot` draw any field or the four band edges against the distance in µm. In the window: Cut line (two clicks; in 3D on the front-most surface point, the line running through the device), the plot in the Curves tab; also typed end points in µm. |
| XD6 | **Run comparison.** A second result (a result file, or the previous run of the same deck) as run B. Curves of A and B overlaid; for transfer curves the threshold shift ΔV_th from Unit 24's extraction. A difference map B − A of a field at a bias point both runs have (equal bias values), on the same mesh (digest) only, else refused; diverging colormap (cool-warm), symmetric about 0: potential, quasi-Fermi levels and band edges in mV; densities as Δlog10 (decades); other fields linear in their units. The record of each run side by side with the settings that differ. |
| XD7 | **Mesh quality.** Per element: in 2D the minimum angle, in 3D the minimum dihedral angle (degrees). A map on the elements (2D) or the boundary and slices (3D) and a histogram; the minimum, the median and the count below a threshold (20° in 2D, 10° in 3D) in the metrics strip. |

**Commits:** 26b·1 live checking, suggestions, completion and the structure preview (XD1–XD4); 26b·2 cut lines
(XD5); 26b·3 run comparison and mesh quality (XD6, XD7).

**Acceptance gates:**

| # | Gate | Commit |
|---|---|---|
| XG1 | Live checking: an error appears at its line after an edit without Check; a check never runs on the UI thread; a stale check's result is never applied (a deterministic test holds the worker); editing during a run checks the deck. | 26b·1 |
| XG2 | The vocabulary equals `deck/build.hpp`'s reference (every section and key it names, no other); every listed key with a sample value gives no "unknown key" error from `deck::build`; every misspelling at one edit of a key, section or choice is suggested back; completion at each caret context gives exactly the fitting candidates. | 26b·1 |
| XG3 | The structure preview's scene equals the solid and mesh scenes of the device `deck::build` makes (`nitcad plot` of the deck's device file) to the bit. | 26b·1 |
| XG4 | Cut lines: exact on a linear field; the samples are the element crossings (a tensor and a generated mesh, 2D and 3D); along a 1D device equal to `band_plot`'s points; a line leaving the mesh has a gap there; the 2D and 3D MOSFETs' band diagrams from source to drain drawn on WARP as references. | 26b·2 |
| XG5 | Comparison: B − A of a run with itself is zero; of two runs it equals the nodal difference; another mesh or a bias point not in both is refused; ΔV_th equals the difference of the extracted thresholds; a WARP reference of a difference map. | 26b·3 |
| XG6 | Mesh quality: exact angles on equilateral, right and regular elements (60°, 45°, 70.5288°); the histogram counts every element once; the generated meshes respect their generators' minimum-angle setting in 2D. | 26b·3 |
| XG7 | The window: every new command in the menus and reachable by the keyboard; a clean D3D12 debug layer and no live graphics object after the window; per-monitor DPI for the new controls. | all |
| XG8 | Suites as WG6 (Release and Debug at every commit, `/W4 /WX`, hand mutants caught); `app` still free of `physics` and `assemble`; no numerical change (the deck, solver and render references unchanged); ARCHITECTURE.md as built after acceptance. | all |

**Commit 26b·1 (live checking, suggestions, completion and the structure preview; XD1–XD4, XG1–XG3, XG7), as
built.**
- `live_check.hpp`: `LiveChecker(wake, check)` with its own `std::jthread`; `submit(text, path)` returns the check's
  generation; the worker takes only the newest submission, so edits made while a build runs coalesce into one check
  after it; a result is posted (wake) and `process()` on the UI thread applies it only when no newer submission exists,
  else counts it as `dropped()`. The check is injectable (tests hold the worker on a latch). Independent of the run
  controller: the deck is checked while a run goes on.
- `deck_vocabulary.hpp` (XD2, XD3): 15 sections, each key with its named-key arity, repeatability, a sample value,
  a line of help, its value choices (first and second word) and, for `[contact]`, `[mesh]` and `[solver]`, the
  context (the selector key's value) in which it applies; `edit_distance`, `nearest` (at most min(2, max(1, ⌊n/3⌋))
  edits, exact matches excluded), `suggestion(text, diagnostic)` from the error's own line (unknown section, key or
  choice), and `complete(text, line, column)`: section kinds after `[` (closing it unless closed or named), the keys
  of the caret's section not yet present (a key and " = ", or a named key and a space), value choices after `=`
  (waveforms with their parenthesis), nothing in a comment, on a continued line, on a header's name or before the
  first section.
- `workspace.hpp`: `workspace_from_structure(plan, deck)` (the device and geometry, no results) and `same_mesh(w,
  mesh)` (by the graph digest); `workspace_from_plan` uses the former.
- The window: an edit restarts a 400 ms timer, then the deck is submitted; a result older than the text shown is not
  applied; an error underlines its line (wavy, drawn after the Rich Edit's paint), the message with " (did you mean
  'x'?)" goes to the status bar and Find Error; a deck that builds is shown as "Device (preview, not solved)" unless
  the results shown are on the same mesh, keeping the view settings when the dimension is unchanged. Completion:
  a list under the caret (Up/Down, Enter or Tab inserts, Escape closes; Left/Right/Home/End/PgUp/PgDn close),
  opened by typing a word or `[` or `=`, and by Edit > Complete (Ctrl+Space); taking a key opens its value choices.
- **Tests** (`tests/app/`, 11 new cases): **XG1** — errors at their lines then the plan, the check on the worker
  thread; a held worker with five submissions: two results dropped, the newest applied; a check during a 1,601-point
  run; in the window: a misspelt material marked at line 11 with "(did you mean 'sio2'?)" and one check for two
  quick edits, a result older than the text not shown, `kind = omic` caught during a run; **XG2** — every vocabulary
  key with its sample value accepted by `deck::build` in every context (more than 130 probes); every word of
  `deck/build.hpp`'s reference that `deck::build` takes is listed (more than 100), parameter paths excepted; the
  menus' defaults equal the vocabulary's; every key, section and choice misspelt by one edit (a letter dropped)
  suggested back from `deck::build`'s own error; completion at each caret context; **XG3** — the preview of the 2D
  MOSFET and a generated 3D deck drawn on WARP equal `nitcad plot`'s solid and mesh PNGs of the deck to the bit;
  the preview follows a finer mesh in the view chosen; results kept on an edit that keeps the mesh; **XG7** —
  completion by the keyboard in the window (Enter, Tab, Up/Down, Escape, Ctrl+Space, a double click).
- **Finding (fixed):** `DeckPane::caret()` gave the previous line on the final empty line (Rich Edit's line index for
  a position after the last `\r`); found by the forced-completion test.
- **Finding (not changed):** `deck/build.hpp`'s reference omits keys `deck::build` takes (the `gmres_ilut` keys
  `ilut_drop`, `ilut_fill`, `restart`, and the model switches by name). The vocabulary lists what `deck::build` takes,
  so XG2's "no other" is checked against `deck::build`, not the comment.
- **Mutation checks (Release): 29 of 29 caught, one equivalent replaced.** Vocabulary, checker and workspace (18):
  the nearest name at 3 edits and completion inside a comment survived the first run and are caught by two added
  checks; three mutants that failed only by `/WX` (an unused variable) were rewritten to compile and are caught.
  Window (11): Enter's character typed after accepting, Down skipping an item, no value choices after a key, not
  opened after `=`, typing a word opens nothing, a result older than the text applied, results dropped on any edit,
  no suggestion, no live check, the view reset on every edit, Ctrl+Space not forced; three `/WX`-only mutants
  rewritten as above, of which Down past the last item is equivalent (the list box refuses an index out of range)
  and was replaced by Down skipping an item, caught by an added step on the 15-item list.
- **Suites:** Release, the app tests (43 cases, 4,423 assertions) and the window tests (11 cases) pass, no warnings
  under `/W4 /WX`; Debug, the app tests (121 s) pass, and the window tests (284 s) pass but for one self-test run
  whose WG1 gap was 248.6 ms (bound 100 ms) while the solver ran; three reruns of that self-test gave 40–51 ms (the
  26·3 record: 39 ms). The self-test makes no edit, so the live checker does not run in it; the gap is recorded as
  an outlier of the Debug build on WARP. The clean-worktree suites at this commit are recorded with 26b·2.

**Commit 26b·2 (cut lines; XD5, XG4, XG7), as built.**
- `render/cut_line.hpp`: `cut_line(geometry, a, b)` — per element (after a bounding-box test), the barycentric weights
  at a and b (Cramer, as `Locator`), affine along the line, give the interval where all are non-negative; a weight
  changing by at most 1e-10 over the line is parallel to its face (a line on a face reads both elements); the
  intervals, sorted, are joined into samples at the crossings (one where two elements meet, a `gap_before` mark
  where the line re-enters the mesh); weights below 1e-12 are zero and the rest renormalised, so a sample on a node
  reads that node exactly. `cut_values` reads a nodal field, skipping nodes of zero weight (a NaN on an insulator node
  leaves the faces away from it defined).
- `scene.hpp`: `cut_plot` and `cut_band_plot` (against "distance (µm)", a NaN point at each gap), sharing
  `profile_plot`'s and `band_plot`'s code (their output unchanged: the references pass); `add_cut_line` (the parts
  inside the mesh, a white 2 px line over a black 4 px one, on top of everything, so in 3D it shows through the
  device).
- `ViewSettings::cut` (ends in µm, in the device's dimension): the profile and bands views plot along it (needed in 2D
  and 3D), the 2D and 3D solid, mesh, field, slice and iso views draw it, any other view refuses it; 2D and 3D
  workspaces now offer the profile and bands views; `plot_arguments` gives `--cut`.
- **Deviation (owner's record): `nitcad plot --cut <x0>,<y0>[,<z0>]:<x1>,<y1>[,<z1>]`** with the same rules and
  messages. XD5 places cut lines "in the window", but VA4 has every window view drawn as `nitcad plot` draws it and
  the WG4 test compares the two images; without the option the window's cut views could not be checked that way.
  The seeds' segment parser became the shared `segment_um`. Also in both: a profile of the doping no longer needs
  results (it never used them), so a cut line on an unsolved deck plots its doping.
- The window: View > Draw cut line (Ctrl+Shift+L, a toolbar button), exclusive with Probe and Measure; two clicks on
  the device view — in 2D the point itself (`Prober::point_at` inverts the view's affine map), **in 3D the front-most
  surface node within 12 px (as Probe), not any point of the surface (a deviation from XD5's "front-most surface
  point": a node is what the view can pick reliably)** — rounded to the nanometre; the Curves tab then shows the band
  diagram along it, or the doping before a run. The inspector's "Cut line (µm)" field types or clears it (Apply). One
  cut line per window, kept while the device keeps its dimension (a run of the deck), cleared when a deck or result
  is opened. The Curves tab offers the profile in 2D and 3D; it starts on a band diagram only in 1D.
- **Tests:** `cut_line_test` (6 cases) — a linear field read to 1e-12 at every sample of tensor and generated meshes
  in 2D and 3D; the samples equal the crossings found independently by intersecting the segment with every element
  face (segment–edge in 2D, Möller–Trumbore in 3D); a horizontal line across a row of cells: 21 samples (11 faces, 10
  diagonals); lines along grid lines and through nodes: the nodes only; along a 1D device the cut band diagram equals
  `band_plot` (x to 1e-12 µm, y exactly, NaN where NaN); a U-shaped mesh: the line enters at x = 0, leaves at 1 µm,
  re-enters at 2 µm with one gap, the plot broken there, the overlay without the gap; refusals. **WARP references**
  (`mosfet2d_cut_bands`, `mosfet3d_cut_bands`, `mosfet3d_cut_line`): source to drain 0.05 µm under the oxide.
  `view_settings_test`: errors as `nitcad plot` reports them; the scene equals render's; the doping without results;
  the 2D MOSFET's band diagram, a profile of n across the oxide and a map with the line **equal to `nitcad plot`'s
  PNGs to the bit**. `run_test`: the option, the line changing the map, the refusals. The window: a typed line on the
  unsolved deck (the doping, the profile offered), kept through a run; drawn by two clicks within 0.01 µm and
  rounded to the nanometre; the band diagram in the Curves tab and the line on the map; typed, malformed, cleared;
  Probe leaving the tool; a 3D preview's line between the surface nodes the clicks pick. `probe_test`: the point
  under a pixel of a 2D view, as the window draws it (zoomed, panned) and from a turned and tilted camera.
- **Mutation checks (Release):** 27 of 27 caught — render 11 (no parallel case, covered pieces repeated, no gap,
  weights not snapped, not normalised, the box test excluding touching elements, zero weights read, equal ends
  taken, no NaN point at a gap, the overlay drawn across a gap, cut bands labelled x), app 6 (a device view without
  the line, 2D bands without a line, bands ignoring the line, `plot_arguments` repeating the start, no 2D profile
  view, the end's dimension unchecked), CLI 1 (`nitcad plot` without the line), window 9 (the Curves tab without the
  2D profile, a device view without the line, clicks ignoring the tool, ends not rounded, the Curves tab not shown,
  a run dropping the line, an unsolved line plotting the potential, the 2D inverse transposed, 3D picking nothing).
  Two survived at first and gained tests: the transposed 2D inverse (every 2D view the window builds looks straight
  down, where the swapped terms are zero — killed by the tilted camera in `probe_test`) and the end's dimension
  (killed by a cut whose start alone is right). A first run whose test filter matched no test was discarded; the
  mutation script now refuses such a filter and stops if it cannot restore a file.
- **Suites:** in the working tree, Release: render (67 cases, 109,200 assertions, the three references), app, the
  CLI cut test and the window tests pass, no warnings under `/W4 /WX`; Debug: render (67 cases, 696 s), app (123 s),
  the CLI and the 26b window tests (29 s) pass. Debug render tests on WARP must run alone: run beside other D3D
  tests they failed to create a device (0x887A0007). The clean-worktree suites at 26b·1 (3077ce4): Release 27 of 27
  (990 s); Debug 19 of 19 — render failed once beside the mutation runs (the device error above) and passed run
  alone (61 cases, 658 s); app 110 s, window 297 s. The clean-worktree suites at this commit are recorded with
  26b·3.

**Commit 26b·3 (run comparison and mesh quality; XD6, XD7, XG5, XG6, XG7), as built.**
- `render/quality.hpp`: `min_angle_deg` — a triangle's smallest angle (each from atan2(|u × v|, u · v)), a
  tetrahedron's smallest dihedral angle (π less the angle between two outward face normals, over the six face pairs);
  0 for a degenerate element (|u × v| ≤ 1e-14 size², |6V| ≤ 1e-14 size³). `element_quality` (one angle per element,
  2D and 3D only), `quality_summary` (minimum, median — the mean of the two middles for an even count — and the count
  strictly below the threshold, 20° in 2D and 10° in 3D), `quality_histogram` (18 bins of 5°, every element once).
  `scene.hpp`: `quality_view` (flat per element, 0 to the best angle — 60° in 2D, 70.53° in 3D — unless a range is
  given; 2D with the element edges, 3D on the boundary or on slices; no contacts) and `quality_plot` (the histogram as
  a step line and the threshold, 0 to 90°). `mesh_view`'s 3D edges now come from the shared `boundary_face_edges`
  (its output unchanged: the references pass).
- `render/compare.hpp`: `matching_point` (the first point of B with the same contacts, each bias within 1e-9 V);
  `difference_field` (B − A at the nodes of A's geometry: densities in Δlog₁₀ decades, NaN unless both are positive;
  the potential in mV; **the band edges and quasi-Fermi levels in meV** — XD6 says mV; they are energies, so meV is
  the unit of the same number); `difference_scale` (symmetric about 0, ±1 when everything is 0 or NaN).
- `app/compare.hpp`: `comparable` (the same mesh, by the meshes' digest), `figure_shifts` (Unit 24's figures both runs
  have, by name and curve, B − A: ΔV_th is the threshold's), `compare_records` (A's lines in order, then the lines
  only B has; `differs` marked).
- `ViewSettings`: views `quality`, `quality-histogram` (2D and 3D with elements, before a run too) and `difference`
  (offered with a run B on the same mesh). `build_scene`, `available_views` and `plot_arguments` take run B: the I-V,
  C-V and transient-current views overlay B's curve of the same contacts (by name), the series named "(A)" and
  "(B)"; the difference map is cool-warm unless another colormap is chosen, scaled by `difference_scale` unless a
  range is given, a profile in 1D. The doping's difference needs no results. Refused: no run B, another mesh, no
  bias point of B at A's biases, a contact of A that B lacks, B without the analysis the curve needs. A transient's
  snapshots are matched by bias like any point (the first equal one). Band diagrams and profiles are not overlaid.
- **Deviation (owner's record, as for 26b·2's `--cut`): `nitcad plot --compare <deck | result file>`** and the views
  `quality`, `quality-histogram`, `difference`, with the same rules and messages, so the window's new views are
  checked against `nitcad plot`'s PNGs (VA4, WG4). A result file is read with run A's device and must have its mesh
  (the result's digest, or the node count for a file without one). A deck's B outputs are loaded as A's are (the
  shared `load_outputs`); `nitcad plot` mesh quality and a doping difference need no results.
- The window: File > Compare with a result... (Ctrl+Shift+C; the device and mesh files found as Open result finds
  them), Simulation > Compare with the previous run (Ctrl+Shift+V) and Stop comparing (Ctrl+Shift+X), View > Mesh
  quality (Ctrl+Shift+Q, 2D and 3D). The previous run is the workspace with results that a run (or the structure
  preview) of the same deck replaced; Compare with the previous run takes it as run B, Stop comparing gives it back.
  Opening a deck or a result forgets both. The metrics strip shows run B's name and ΔV_th when comparing, and on
  the quality views "Min angle, median, N below T°". The Results tab adds "Run B" and "B − A" columns to the figures,
  and the run records side by side ("≠" before a line that differs). Probing a difference map reads B − A.
- **Tests:** `quality_test` (4 cases) — exact angles: equilateral 60°, right isosceles 45°, 30-60-90 30°, the regular
  tetrahedron 70.5288°, a cube corner's tetrahedron acos(1/√3), a Kuhn tetrahedron 45°, degenerate 0 (collinear, coincident,
  flat, and collinear to rounding), a sliver under 0.2°; tensor meshes 45° in 2D and 3D; generated 2D meshes at 20.7°, 25° and 30° settings: no element below
  the setting other than those the generator reports exempt, the minimum equal to its report; the 3D generator's
  minimum dihedral equal to its report; the histogram's bins and the summary (odd and even counts, the threshold
  itself not below); the map and the plot. `compare_test` (render, 3 cases) — matching points (another number of
  contacts refused, also one equal as far as A goes); a run against itself zero in every field it holds (at least
  11); two runs: the potential's difference in mV, n's in decades, the donors' log 4, the scale symmetric. **WARP reference** `mosfet2d_difference_potential` (the drain at 0.5 V less 0 V: +500 mV in
  the drain pocket, the depletion region between). `compare_test` (app, 2 cases) — the 2D MOSFET at gate work
  functions 4.1 and 4.3 eV: ΔV_th equals the difference of the two extracted thresholds and is positive; the shifts
  in A's order; the input identity marked as differing, the mesh not; a run against itself: no shift, nothing marked;
  two meshes not comparable, their figures still compared; the diode swept from the cathode: the same figure on
  another curve, no shift; an unsolved deck beside a run: B's lines after A's, marked. `view_settings_test`: the views offered (with and without
  B, another mesh), the refusals, the overlay's series, the plot arguments; the quality map, the histogram, two
  difference maps and the overlaid I-V **equal to `nitcad plot`'s PNGs to the bit**. `run_test`: B as a deck and as
  its result file give the same image; the refusals; mesh quality and the doping difference without results. The
  window: Mesh quality (the strip reads "Min angle 26.6°, median 26.6°, 0 below 20°": the MOSFET's 0.1 × 0.05 µm
  cells give atan ½), a second run with another work function, Compare with the previous run (the difference map,
  ΔV_th negative in the strip, the Results columns, the marked identity), Stop comparing (the field view back, the
  previous run kept), a result file as run B, a new deck forgetting both. The command test covers the four new
  commands' menus and unique accelerators.
- **Mutation checks (Release):** 40 of 41 caught; the 41st caught in Debug only. Render 21 (a triangle's third angle
  ignored, an angle of the wrong sign, collinear triangles not degenerate, normals not outward, the normals' angle as
  the dihedral, a face pair skipped, the 3D threshold 20°, the even median the upper middle, the threshold itself
  below, 90° in the first bin, bins rounded, a loose bias match, the contact count unchecked, A − B in decades, the
  potential in V, a zero density in B read, the scale from the largest positive, the scale not symmetric, the quality
  map's scale to 90°, every element the first's colour, the threshold line at 0), app 13 (comparable ignoring the
  mesh, A and B swapped, figures matched by name alone, differing lines not marked, lines only B has dropped, the
  difference offered on another mesh, the difference map in turbo, its scale fitted, `plot_arguments` without
  `--compare`, the overlay unnamed by run, I-V without B's curve, the histogram drawn as a map, no histogram view),
  CLI 5 (a deck B on another mesh taken, no bias-point refusal, a result file B ignored, the histogram or the doping
  difference needing results), window 11 (no previous run kept, the previous run kept after use, Stop comparing
  losing it, a new deck keeping run B, Mesh quality showing the field, views without run B, no ΔV_th in the strip, a
  difference view kept without B, the 3D threshold in 2D, A − B in the Results, differing lines not marked). Four
  survived at first and gained tests: collinear to rounding (now 0 exactly), figures matched by name alone (the
  diode swept from the cathode: the same figure on another curve), lines only B has (an unsolved deck beside a run),
  and the contact count — whose mutant reads past A's biases: a run B equal to A as far as A goes is now refused,
  caught by Debug's checked `vector` (an assertion), not in Release, where the read is unchecked.
- **Suites:** in the working tree, Release 27 of 27 (974 s); Debug 19 of 19 — `app_window` failed once in the suite
  and passed run alone (271 s): the GPU self-test's WG1 gate saw a 226 ms gap in the message loop. A first run
  failed the same way in Release (652 ms, `app_window_3d`), passing 4 of 4 after. **Open item (WG1 gap, deferred
  by the owner for investigation, 2026-10-22):** a self-test's message loop occasionally stalls past 100 ms during
  a run — on the GPU here, and in 26b·2's build as well (1,623 ms on a first GPU run, then 41 and 38 ms; 26b·3's:
  115, then 42 and 43 ms), and on WARP in 26b·1's Debug suite (248.6 ms). Not caused by this commit (26b·2's build
  shows it); the cause is not established — a first run after the adapter idles is one reading, but the WARP case
  does not fit it. WG1 is left as ruled. No warnings under `/W4 /WX`. The tests added for the four survivors came after the suites; they pass in
  Release and Debug (render quality and compare, app compare). The clean-worktree suites at 26b·2 (394c2e9): Release 27 of 27, Debug 19 of 19
  (1,687 s). The clean-worktree suites at this commit are recorded with the acceptance.

**Status:** Unit 26b built in three commits (3077ce4, 394c2e9, eb79729; XG1–XG8 met on this workstation, the WG1
gap above open). **Accepted by the owner 2026-10-22** ("accept Unit 26b"), the WG1 gap to be investigated later;
ARCHITECTURE.md 6.14, 8, 11 and 11.3 updated to "as built" (XG8). Pushed and merged to `main` as `51face0`
(PR #60, 2026-10-22, on the owner's word, the hosted CI still running). The clean-worktree suites at `eb79729`
(finished after the merge): Release 27 of 27 (980 s); Debug 18 of 19 (1,951 s) — `app_window` failed only on the
WG1 gap above (the GPU self-test's message loop stalled 278.8 ms), the open item.

**Not in Unit 26b (VA7 and later):** form editors that write the deck, project management, remote runs, plugins,
localization, an installer; completion of region parameter paths; more than two runs compared at once.
