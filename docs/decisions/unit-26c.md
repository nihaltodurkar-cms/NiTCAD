<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 26c — GUI catch-up (proposed 2026-10-22; on `app/gui-catch-up`)

**Owner request:** "GUI catch-up" (2026-10-22), chosen after Unit 26b's acceptance: bring into the window what
`nitcad plot` draws but the window does not, and what Units 25c and 26 left for the application. The owner's
answers of the same day: transient snapshots as **playback only** (the tiled filmstrip stays in `nitcad plot`); the
clip box edited by **inspector sliders** (no in-view handles); **`nitcad plot --hide <region>`** added for region
visibility. The owner's conditions of Units 25, 26 and 26b hold: numerics and layer boundaries preserved, Release
and Debug suites at every commit, `/W4 /WX`, hand mutants, nothing pushed or merged without authorization.

**What exists (audited on `app/analysis-workspace` at `195e829`).** `nitcad plot` draws arrows (`--vectors`,
`--arrows`, 300 in 2D and 400 in 3D), streamlines (`--streamlines`, `--lines` 24, `--from <contact>`, `--seeds` a
line or a 3D plane), the exploded view (`--axis`, `--gap`) and filmstrips (`--every`, `--columns`, `--frames`;
DC4: at most 12 frames on one scale fitted over all). `ViewSettings` already holds the curves' swept and sensed
contacts, the C-V frequency and the log axis, but the inspector shows no control for them. The scene tree lists the
regions without any action. `render::ViewOptions` has no notion of a hidden region; no `nitcad plot` option hides
one. Units 25a (VR7) and 25c (DC4, DC6) left animation playback and interactive clip editing "to 26"; Unit 26 did
not take them (ARCHITECTURE.md 6.14, "Not in Unit 26").

**Design:**

| # | Choice |
|---|---|
| YD1 | **Arrows and streamlines.** Inspector controls on the field view: Arrows (none, E, J, Jn, Jp) with a count; Streamlines (none, E, J, Jn, Jp) with a count and seeds from a contact or along a typed seed line (2D) or plane (3D). The same rules, defaults and messages as `nitcad plot` (the map shows the vector field's own magnitude); `ViewSettings` gains the options and `plot_arguments` gives `--vectors`, `--arrows`, `--streamlines`, `--lines`, `--from`, `--seeds`. |
| YD2 | **Exploded view.** A 3D view kind with the axis (x, y, z) and gap (µm, default 12% of the extent) in the inspector, as `--view exploded`. |
| YD3 | **Curve pickers.** On the I-V, C-V and current views the inspector shows the swept contact, the sensed contact, the C-V frequency and the log axis (the existing `ViewSettings` fields). |
| YD4 | **Region visibility.** A check box per region in the scene tree hides it. `render::ViewOptions::hidden_regions`: the solid, mesh, field, slice, iso, exploded and quality views leave out a hidden region's control-volume pieces (2D) or elements and boundary faces (3D); the faces a hidden region shared with a shown one become boundary faces; contacts on hidden regions only are left out; scales are fitted over the shown regions. `nitcad plot --hide <region>` (repeatable) with the same rules; an unknown region name is refused. |
| YD5 | **Transient playback.** On a transient run's device view: Play/Pause, a time slider over the snapshots and step keys (previous, next, first, last). The scale is fitted once over all snapshots (as DC4) and stays fixed; a frame is drawn by a timer (5 frames a second), never by a loop on the UI thread, and a frame not yet drawn is skipped rather than queued. The tiled filmstrip stays in `nitcad plot`. |
| YD6 | **Clip box sliders.** In 3D, a low and a high slider per axis over the device's extent (1,000 steps, values rounded to the nanometre) beside the typed box; moving one updates the box and the view. |

**Commits:** 26c·1 arrows, streamlines, the exploded view and the curve pickers (YD1–YD3; `app` only); 26c·2 region
visibility (YD4; `render`, `nitcad plot`, `app`); 26c·3 transient playback and the clip sliders (YD5, YD6).

**Acceptance gates:**

| # | Gate | Commit |
|---|---|---|
| YG1 | Every new view and overlay in the window — arrows and streamlines in 2D and 3D, the exploded view, the curves with chosen contacts and frequency, hidden regions — equals `nitcad plot`'s PNG for the same options to the bit (WG4's choices). | 26c·1, 26c·2 |
| YG2 | Region visibility: with no region hidden every scene is unchanged (the existing references pass); a hidden region's pieces, elements and contacts are absent and the shown regions' are unchanged; the shared faces close the shown regions; WARP references of the 2D and 3D MOSFETs with the oxide hidden. | 26c·2 |
| YG3 | Playback: every frame equals the filmstrip frame of the same snapshot on the same scale to the bit; the window's message loop meets WG1 while playing; Pause, the slider and the step keys land on the snapshot asked for. | 26c·3 |
| YG4 | Clip sliders: a slider's value and the typed box agree both ways; the view equals `nitcad plot --clip` with the slider's box. | 26c·3 |
| YG5 | The window: every new command in the menus and reachable by the keyboard; a clean D3D12 debug layer and no live graphics object after the window; per-monitor DPI for the new controls. | all |
| YG6 | Suites as WG6 (Release and Debug at every commit, `/W4 /WX`, hand mutants caught); `app` still free of `physics` and `assemble`; no numerical change; ARCHITECTURE.md as built after acceptance. | all |

**Started on the owner's word (2026-10-22, "start 26c.1 also investigate app_window").**

**Commit 26c·1 (arrows, streamlines, the exploded view, the curve pickers; YD1–YD3, YG1, YG5, YG6), as built.**
- `ViewSettings`: `vectors`, `streamlines` (E, J, Jn, Jp; empty none), `arrows`, `lines` (0: `nitcad plot`'s 300 in
  2D, 400 in 3D; 24), `from` (a contact; empty: its choice), `seeds` (a segment or, in 3D, a plane `x=<um>`); the
  view kind `exploded` with `axis`, `gap_um` (empty: 12% of the extent) and `exploded_field`. `build_scene` draws them
  with `nitcad plot`'s calls, order and messages: streamlines, then arrows, then the cut line; in 3D over their own
  magnitude the arrows on the solid device and the lines in the outline, colored by it; 3D arrows on the slices of a
  slice view, else on the exposed boundary; seeds from the named contact, else for currents the one carrying the
  most, for E the first gate or electrode, else the most biased. `plot_arguments` gives `--streamlines --lines
  --from --seeds --vectors --arrows --axis --gap` (and `--field` for an exploded view showing the field).
  `available_views` offers `exploded` in 3D, after the slices.
- **Arrows and streamlines on the field views only** (`takes_vectors`: field, slice, iso — YD1's "on the field view"):
  a setting keeps them, `build_scene` and `plot_arguments` leave them out on the other views, so the window and
  `nitcad plot` agree on every view (the exploded view in particular, where `nitcad plot` would draw them where
  the regions were before they moved).
- The inspector: Arrows (count), Streamlines (count), Seeds from (or typed, µm) — each a list and a box side by
  side — on the field views of a 2D or 3D run's results; Exploded along (gap, µm) and Show the field on the exploded
  view; Swept contact (Contact on the current view: "(all)" first), Sensed contact, Frequency (the run's frequencies)
  and Log axis on the I-V, C-V and current views. A group shows only with its views (`show_inspector_groups`; the
  inspector lays out the shown ones). Choosing arrows or streamlines on the field view shows that vector's magnitude
  (as `nitcad plot` without `--field`). Counts, seeds and the gap are read on Apply; an invalid count or gap leaves
  the setting and says why in the status bar. Every new control is a tab stop (YG5); no new command. **No probe or
  cut-line click on the exploded view** (its regions moved; `pickable`).
- **Tests:** `view_settings_test` — the exploded view's name and place in the views; `nitcad plot`'s messages (an
  unknown vector field, the counts' limits, a seed plane in 2D, a malformed segment or plane, an unknown contact, a
  1D device, the exploded view's device, axis and gap); vectors left out off the field views (the exploded scene
  unchanged); the plot arguments. **To the bit equal to `nitcad plot`'s PNGs (YG1):** 2D arrows of J over |J| and of
  E over contour levels, streamlines of E (from the gate), of J (from the contact carrying the most), of Jn from a
  named contact with arrows, from a typed segment; 3D colored arrows of J, colored lines of E, lines and arrows of J
  together, arrows of E on a slice, lines from a seed plane; the exploded view with and without a field; I-V of a
  chosen sensed contact, C-V at the second frequency, the transient current of one contact. New fixtures: a small 3D
  MOSFET (9 × 7 × 3 nodes, about a second in Debug), the diode's C-V at 1 kHz and 1 GHz. `analysis_window_test`: the
  inspector driven as a user would — groups shown and hidden per view, arrows of J (the map switches to |J|),
  counts and a refused count, lines from the gate and from a plane, the exploded view (axis, field, gap, a refused
  gap; no probe where the solid view probes), back to the field view with the lines kept; the window's images equal
  `nitcad plot`'s for its settings; tab stops; at 200% DPI the pairs side by side; the C-V pickers (frequencies,
  contacts, log, the default again).
- **Mutation checks (Release):** 31 of 31 caught — view settings 19 (Jn read as Jp, the 2D arrow and line defaults,
  colored arrows beside lines, the line color, the arrow length, 3D arrows not on the slice, the seed plane's unit,
  `--from`'s contact, the contact carrying most, E lines not from the gate, the default gap's axis, the exploded
  field ignored, `--axis`'s letter, iso without vectors, vectors kept off the field views, no exploded view, 3D
  lines not colored, the arrow limit), window 12 (the map keeping its field, a count of 0, untrimmed seeds, the
  arrows' controls on every view, a probe on the exploded view, the frequency picker, no default entry, the pairs
  overlapping, Log and Show the field swapped, swept and sensed swapped, the exploded controls always, a negative
  gap). Two survived at first and gained tests: 3D arrows and lines of one field together, and the exploded view's
  probe (the first test switched probing off by toggling it).
- **Suites (working tree, with the two changes below):** Release 27 of 27 (1,391 s), Debug 19 of 19 (2,382 s); no
  warnings under `/W4 /WX`. The two tests strengthened after the mutants pass in Release and Debug.
- **Found:** `nitcad plot <deck>` refused arrows and streamlines on the solid, mesh and quality views ("the results
  hold no field to plot"); fixed apart, on the owner's word (the next commit).

**Fix (owner, 2026-10-22, "fix it"): `nitcad plot` loads a deck's results for arrows and streamlines.** Reading a
deck, `nitcad plot` loaded the results only for the views that show them, so `--vectors` or `--streamlines` on the
solid, mesh, quality or doping views failed with "the results hold no field to plot" (a result file worked). Now
either option loads them on any view. `run_test` gains arrows on the solid view, streamlines on the mesh view,
arrows over the doping and 3D arrows on the quality view; the mutant that drops the options from the condition is
caught. No change to any image the command could draw before.

**The WG1 gap investigated (owner, 2026-10-22, "investigate app_window"; fixes on "apply them too").** The
self-test's message loop was logged message by message (each dispatched message, its window class and time, and the
wait in `GetMessage` between two ticks of the 10 ms gap timer) over six runs alone and five beside six `nitcad run`
loops. Two causes:
1. **The desktop user's input reached the self-test's window**, which opened shown and activated over the user's
   work: a title-bar drag (a modal move loop, 103 ms, the gap timer unseen meanwhile), clicks on buttons and lists
   (41–65 ms each) and streams of mouse moves over the view (each a hover probe of 1–12 ms in Debug; input is taken
   before the low-priority `WM_TIMER`, so a steady stream starves the timer: gaps up to 454 ms with the window
   working throughout). The earlier failures (1,623, 652, 279, 249 ms) came from suites run while the workstation
   was in use; that they were all input is likely, not proven.
2. **The run monitor's list repainted slowly in Debug**: 30 ms alone, 118 ms under load, with no input — each cell
   re-formatted its whole row.

Fixes: `WindowOptions::input` (default true); the self-test's window is shown without activation and disabled (it
still paints and runs its timers), so WG1 measures the window's own work. `RunPage` keeps each row's cells, formatted
when the row changes (`PointRow` gains `==`). Tests: a window made without input is disabled, a usual one not; after
a run every monitor cell equals the final rows' formatting. Mutants: 2 of 2 caught (rows kept stale, the window
taking input). Afterwards the Debug GPU self-test's longest gaps were 36.1, 38.5 and 36.6 ms, and Debug `app_window`
passed in the suite (19 of 19). WG1 stands as ruled.

**Commit 26c·2 (region visibility; YD4, YG1, YG2, YG5, YG6), as built.**
- `render::ViewOptions::hidden_regions` (region ids; an id not the device's is refused, "a hidden region is not one
  of the device's"). **Hiding goes by control volume, in every view:** a node is hidden with its region, and what is
  drawn is the shown nodes' pieces (`Geometry::owner_pieces`), so a 2D map's area and a 3D surface's enclosed volume
  are the shown regions' control volumes. Edges are drawn as the shown node's half (an edge's midpoint is where the
  two control volumes meet). The interface pieces between a shown and a hidden node become boundary faces of the
  shown side, turned toward the hidden node, so the shown regions stay closed (solid, field, clipped, opaque,
  quality, the slice and iso views' outline shell). A contact all of whose nodes are hidden is left out. Scales fit
  over the shown nodes (`render::shown_values`: NaN on hidden nodes), in the field, slice, iso, difference and
  filmstrip views (one scale over every frame's shown nodes). Isosurfaces are split into the shown pieces of each
  triangle (`IsoSurface::elements` gives each triangle's element), with normals interpolated from the corners. The
  exploded view leaves the hidden regions out and keeps the others' places. 3D quality draws the elements' pieces
  in the shown control volumes and the opened faces in their element's color, not whole elements.
- `nitcad plot --hide <region>` (repeatable; twice is once) on the solid, mesh, field, slice, iso, exploded,
  quality, difference and filmstrip views of a 2D or 3D device. Refused: an unknown name ("there is no region
  <name>"), another view or a 1D device ("hidden regions apply to the solid, mesh, field, slice, iso, exploded,
  quality, difference and filmstrip views of a 2D or 3D device"), no name (usage).
- The window: a check box per region in the scene tree (contacts and other items have none). Clearing a box hides the
  region, checking it shows it again; the status bar says which. View ▸ Show every region (Ctrl+Shift+H) checks them
  all, enabled while any is hidden. `ViewSettings::hidden` (names, in the device's order) is kept on the device views
  (`takes_hidden`: solid, mesh, field, slice, iso, exploded, quality, difference); other views leave it out, and
  `plot_arguments` gives `--hide` for each. The probe doesn't pick a hidden node; in 3D it picks the shown nodes of
  the faces the hidden regions opened. A new run resets the view settings, hidden regions with them. **Without a
  comctl32 v6 manifest the tree sends no TVN_ITEMCHANGED**, so a click or the space bar on the tree posts a message
  that reads the boxes once the tree has toggled them. `TVS_CHECKBOXES` gives every new item a box, so the window
  removes it from items that aren't regions.
- **Not hidden:** arrows and streamlines (they belong to the field, not to a region); the mesh view's edges in 3D are
  the shown halves of the domain boundary, with no edges on the opened faces.
- **Tests:** `hidden_test` (new) — nothing hidden draws as before, and an id that is not a region's is refused; 2D:
  the pieces' area is the silicon's control volumes (1.7 µm²), the gate is gone, outlines and mesh lie below the
  interface, which is drawn across the device, scales are fitted over the silicon, the filmstrip has one scale over
  shown frames, quality; 3D: the solid and field surfaces closed (volume and vector area) with either region hidden,
  opaque and clipped too, the mesh, slices (area) with their outline closed, a curved isosurface cut at the
  interface with its corners' normals, quality closed, exploded; `shown_values`. WARP references
  `mosfet2d_hidden_oxide` (potential, mesh and contours) and `mosfet3d_hidden_oxide` (n). `run_test`: `nitcad plot
  --hide` on ten 2D and 3D views, hiding changes the picture, twice is once, a filmstrip's scale is the silicon's
  (the same as `--range` of it), the refusals. `view_settings_test`: names and refusals as `nitcad plot`'s; **to the
  bit equal to `nitcad plot --hide` (YG1)** in 13 cases (2D solid, mesh, field, quality, iso with edges, difference
  over run B; 3D solid, field, exploded, quality, slice, iso, clip), the difference's scale fitted over the shown
  nodes. `probe_test`: hidden nodes not picked, the opened faces' nodes picked in 3D. `analysis_window_test`: boxes on
  the regions only, the space bar hides and shows, the window equal to `nitcad plot --hide`, no probe where the oxide
  was, kept on the mesh view, both hidden in the device's order, Show every region, a run resets, the menu entry.
- **Mutation checks (Release):** 35 of 36 caught. Render 18 (an id range off by one, no node hidden, whole edges,
  hidden edges kept, no opened faces, opened faces turned inward, the 2D interface left out, hidden contacts kept,
  the colorbar, the filmstrip and `shown_values` over every node, iso pieces unfiltered, iso normal weights, the
  exploded view keeping the region, a 3D interface into the hidden region, the solid, quality and outline shells
  not closed); `nitcad plot` 6 (difference refusing `--hide`, 1D accepting it, the field, filmstrip and difference
  fits over every node, `--hide` ignored); app 6 (difference dropping the regions, no `--hide` arguments, the field
  and difference fits, hidden nodes picked, opened faces not picked); window 5 of 6 (contacts keeping a box, Show
  every region doing nothing, the space bar ignored, the probe picking hidden nodes, the order). **Equivalent:** the
  prober hiding on every view, since every view that can be probed takes hidden regions. Four survived at first and
  gained tests: the iso normals (the isosurface was a plane), the outline shell, the filmstrip's fit, and
  `nitcad plot`'s fits (now also run against the window's to-the-bit comparisons).
- **Suites (clean worktree at `4a9d236`):** Release 27 of 27 (1,308 s), no warnings under `/W4 /WX`; Debug 18 of 19
  (1,526 s): `render` failed after 10 s, all three of its failures `D3D12CreateDevice failed (HRESULT 0x887A0007)`
  (DXGI_ERROR_DEVICE_RESET, before anything was drawn), while the other D3D suites of the same run passed; rerun
  alone in the same worktree it passed (651 s).

**Commit 26c·3 (transient playback and the clipping box's sliders; YD5, YD6, YG3, YG4, YG5, YG6), as built.**
- **Playback** on the field view of a transient run of a 2D or 3D device (`app::takes_film`), the only device view
  `nitcad plot --view filmstrip` draws. `app::film_scale` fits one scale over every snapshot as the filmstrip does
  with `--every 1` (DC4): the range if given, over the shown regions' nodes (YD4), a vector magnitude on a log scale
  floored ten decades below its largest value over all. The window fits it once and keeps it in
  `ViewSettings::film`, fitting again only when the field, scale, range or hidden regions change; a new workspace
  stops playback. With a film, `build_scene` draws the snapshot (`point`) through `render::filmstrip_scenes`, the
  filmstrip frame's own call: the same scale, and the time in the title ("ψ (V), t = 200 ps"). `plot_arguments`
  gives the snapshot's field view (`--point k`); no `nitcad plot` option draws one frame alone, so YG3 compares
  with `--frames`.
- **Controls:** in the inspector, a Snapshot slider over the snapshots, the time and "k of n" beside **Play/Pause**.
  Simulation ▸ Play / pause snapshots (Ctrl+Shift+Space), First (Alt+Home), Previous (Alt+Left), Next (Alt+Right)
  and Last snapshot (Alt+End), enabled while playback applies. Play starts from the snapshot shown, or from the
  first when the last is shown, and stops at the last. Pause, a step key, the slider, the point list and another
  tab each pause it; a step or the slider lands on the snapshot asked for. **Frames:** a 40 ms timer shows the
  snapshot the time since Play has reached (`app::playback_frame`, five a second), so a frame drawn late skips
  those it missed rather than queueing them, and the message loop is never held for more than a frame's drawing.
- **Clip sliders** on the solid and field views of a 3D device: per axis a low and a high slider side by side,
  1,000 steps over the device's extent, values rounded to the nanometre (`app::clip_value_um`, `app::clip_step`).
  Moving one makes the box over the whole device if there was none (Clipping box checked, the box typed), sets that
  side, and moves the axis's other side with it rather than letting them cross. A typed box (Apply) moves the
  sliders to its nearest steps; on an extent under 1 µm a step is finer than a nanometre, so a value, not a step,
  comes back unchanged.
- **Tests:** `view_settings_test` — playback only on a transient field view (not a sweep, 1D or the solid view),
  the refusals; the film holds every snapshot's own scale; a frame's title and kept scale; |J|'s floor and a range;
  the frame draws values below the film's floor at it; `playback_frame` (five a second, a late tick skipping, the
  last kept); the sliders' values, steps and round trips. **To the bit equal to `nitcad plot --view filmstrip
  --every 1 --frames` (YG3)** for every snapshot of a ramped 2D MOSFET (10 snapshots): ψ with contours, |J| on its
  log scale, n with the oxide hidden, ψ linear with the oxide hidden. `analysis_window_test`: the controls and
  commands on the field view, the last snapshot as its frame; First, Previous (stays at the first), Next, Last
  (stays at the last) and the slider land on their snapshot, the point list following; Play from the first runs to
  the last, every frame shown equal to its frame and none before its time, the longest message 100 ms or less
  (WG1); Play from the last starts again, Pause holds, the button plays, the slider and another tab pause; none on
  the mesh view; the menu and keys, tab stops; with the oxide hidden the scale is fitted again and the frame equals
  `--hide oxide`'s. The 3D MOSFET's sliders: six on the solid view at the whole extent, a quarter of x giving the box
  `0.5,0,0:2,1,0.8` (checked, typed) and the view equal to `nitcad plot --clip` (YG4), the high z side moved below
  the low one, a typed box moving the sliders to its steps, the field view keeping them and the mesh view without,
  side by side at 200% (YG5).
- **Mutation checks (Release):** 28 of 28 caught. App 13: the frame not floored, the film floored nine decades, over
  hidden nodes, ignoring the range, the frame titled at t = 0 and the film on every view (both caught by the
  `/WX` build, an unused variable), the frame fitting its own scale, the film in 1D, the frame rounded, past the
  last, the sliders to 0.1 nm, floored and unclamped. Window 15: Play from the last staying there, playing past
  the end, another tab or the slider leaving it playing, Previous from the first, the high side crossing the low,
  the typed box not on the sliders, a new box empty, the film not fitted again on hiding, the button's text,
  overlapping sliders, a frame every tick, the timer ignored, the controls on every view, the point list left
  behind. Five survived at first and gained tests: the floor, past the last, hiding while playing, frames before
  their time, and the point list.
- **Suites (clean worktree at `3d0e78a`):** Release 27 of 27 (1,013 s), no warnings under `/W4 /WX`; Debug 18 of
  19. Tests 1–13 passed in order (767 s; the run was stopped there), tests 14–19 run six at a time (647 s):
  `app_window` failed WG1 in the WARP self-test, the message loop's longest gap 102.1 ms (bound 100 ms; 19 of its 20
  cases passed); rerun alone it failed the same way (125.8 ms). The WARP self-test was then run twice alone in this
  worktree and twice in the one at `4c4b75b` (the WG1 fix, before 26c·2 and 26c·3): at `3d0e78a` both passed, every
  gap 36.0–42.1 ms; at `4c4b75b` both failed, one view's gap 1,095.9 ms and 215.2 ms, the other views' 36.6–38.6 ms.
  The gaps are rare single spikes on a loop that otherwise holds about 40 ms, in the build before this unit's window
  changes as in this one, so they are not this commit's. **Open item (WG1 spikes, owner, 2026-10-23, "record as
  is"):** the window takes no input since the WG1 fix, so the spikes have another cause, unknown; WG1 is left as
  ruled.
- **Owner ruling 2026-10-31 (on Unit 53's Debug suite, "opt 1"):** the debug-layer self-test (WG4, WG5:
  WARP, `--debug-layer`, Debug) failed WG1 with one 104.6 ms gap, run serially, and passed alone. A WG1 failure
  ends the self-test before the debug layer's messages and the live objects are read, so the spike hid WG5. That
  one run now takes `--gap-limit-ms 1000`: it gates WG5, and WG1 stays gated at `wg1_limit_ms()` (100 ms, 250 ms on
  the hosted runner) by every other self-test. The spikes' cause stays open; WG1 is not raised.

**Status:** Unit 26c built in three commits (c16f8e0, 4a9d236, 3d0e78a), with the `nitcad plot` fix (3e730dd) and
the WG1 fixes (4c4b75b); YG1–YG6 met on this workstation, the WG1 spikes above open. **Accepted by the owner
2026-10-23** ("accept it and merge"); ARCHITECTURE.md 6.14, 11 and 11.3 updated to "as built" (YG6). Merged to
`main` as `1be857e` (PR #61, 2026-10-23, on the owner's word, the hosted CI still queued).

**Not in Unit 26c:** in-view clip handles (the owner chose sliders), the tiled filmstrip in the window, video
export, oblique clip planes, streamtubes.
