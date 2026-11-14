<!-- Archived from DECISIONS.md on 2026-11-14 (branch docs/decisions-compress), word for word; DECISIONS.md keeps its summary. -->

### Unit 26 — the Windows application (authorized 2026-10-22; on `app/windows-application`)

**Owner request:** VQ3/VQ6 order (25a, 45, 25b, 25c, then 26), VA1–VA8 and gates WG1–WG6 (Units 25 and 26 above),
the detailed plan's start page, run monitor and results workspace (page Figures 3, 5, 6), VQ11 (bundled examples) and
VQ13 (menus open the deck at a section). On 2026-10-22 the plan page was checked against this file and the repository
and four points were raised (AQ1–AQ4); the owner ruled on them the same day and gave the go: "GO — Unit 26 is
authorized."

**What exists (audited on `main` at `e289e44`).** `solve::RunControl` (a stop token and a progress callback, 6.9);
`deck::parse`, `deck::build` (errors `"<source>:<line>: <what>"` with the line as the context index) and `deck::run`
(a stopped run is still a report with its completed part); `render` with off-screen D3D12/Direct2D drawing to images,
an orthographic and a perspective camera (`--perspective`, reference `mosfet3d_opaque_perspective`, Unit 25a), slices
(25a), the clipping box (25c); `analysis` extractors; `io` loaders. No window, no swap chain, no `app` code, no
`examples/` directory. The layer table allows `app` everything except `physics` and `assemble` (ARCHITECTURE.md 3).

**Rulings (2026-10-22):**

| # | Point | Ruling |
|---|---|---|
| AQ1 | The planned two commits do not cover the start page, run monitor and results workspace. | **Approved: three commits.** 26·1 the app logic; 26·2 the main window; 26·3 the start page, run monitor, results workspace, the four example decks and `--self-test`. |
| AQ2 | WG6 bounds scope to VA1–VA6, but VA1 and VA8 are Unit 26's. | **Approved:** WG6 covers VA1–VA6 and VA8–VA9; WG7 (menu navigation) and WG8 (every example deck tested) added; VA7's disposition stated explicitly (below). |
| AQ3 | No unit owns the examples (VQ11). | **Approved:** 26·3 owns `examples/` and the four start-page decks. |
| AQ4 | How the application is launched. | **Approved:** a separate `nitcad-app.exe`; no GUI or window code in the `nitcad` batch runner. |

Owner instructions with the rulings: preserve the architecture and the acceptance rulings; verify whether hosted CI
can open a window for `--self-test`, and if not make WG4 local-only and keep the off-screen WARP image tests in CI,
without letting the question block the work; commits scoped to their approved responsibilities; report commit hashes,
test results and deviations; no unrelated code or scope.

**VA9 (added by AQ3, AQ4).** `nitcad-app.exe`, a `/SUBSYSTEM:WINDOWS` executable of the `app` layer (L8);
`examples/` at the top level with four decks (1D pn diode, 2D MOSFET, 3D MOSFET, MOS capacitor) listed on the start
page, each run by a test.

**VA7's disposition (AQ2).** VA7 is the exclusion list and grants nothing: a graphical device editor (form editors that
write the deck, VQ13), project management, remote or cluster runs, plugins, localization and an installer are not
built in Unit 26, and no gate requires them. WG6 checks the work against VA1–VA6 and VA8–VA9 and against VA7's list.
The detailed plan's analysis workspace (live checking and completion, structure preview, cut lines, run comparison,
difference maps, mesh quality; page Figures 4, 7, 8) is Unit 26b (VQ10), not 26.

**Implementation plan (within VA1–VA9):**
- `app` logic, plain C++ under test (`include/NiTCAD/app/`, `src/app/`): the deck document (text, path, modified
  flag, load and save, a check that maps `deck::parse`/`deck::build` errors to their lines); the run controller (one
  worker `std::jthread` per run; the states idle, checking, running, stopping, done, failed; progress events queued
  on the worker and drained on the UI thread, which a caller-supplied wake function (a posted Win32 message in the
  window) prompts; Stop requests the token and the completed part is kept); the workspace (a run's report or a
  result file with its device and mesh, paired by the mesh digest as `nitcad plot` does; the run record, provenance,
  convergence history and Unit 24's figures); view settings (view, field, point, scale, range, colormap, contours,
  edges, insulators, projection, camera, slices, clip) turned into `render` scenes with `nitcad plot`'s choices; the
  menus' deck sections and keys (VQ13); the recent files and layout, saved per user.
- The window (`src/app/`, `app/main.cpp`): per-monitor DPI aware v2; the workspace of Figure 1 (menu bar, tool bar,
  scene tree, viewport tabs, metrics strip, inspector, status bar); the deck pane (Rich Edit with line numbers, errors
  at their lines); the viewport presenting `render` scenes in a DXGI swap chain on the renderer's queue; orbit, pan,
  zoom, fit, probe and measure (VA8). `render` gains what a window needs and no more: a swap-chain target that shows
  the image `draw` would return, a zoom and pan on the camera (defaults leave every image unchanged), and the screen
  position of a world point for picking.
- 26·3: the start page (new, open, recent, examples), the run monitor (a row per point, the live convergence plot,
  Stop), the results workspace (curves, Unit 24's figures, the run record), `examples/`, and `--self-test <deck>`.
- **Perspective (AQ2's proposed "camera mode with a reference image")** already exists from Unit 25a (`Projection::
  perspective`, `nitcad plot --perspective`, the WARP reference `mosfet3d_opaque_perspective`); Unit 26 only exposes it.

**Acceptance gates (final):**

| # | Gate | Commit |
|---|---|---|
| WG1 | No solve on the UI thread: during a long run the UI thread's message loop keeps running (a timer message at least every 100 ms); Stop ends the run with its completed part kept and shown. | 26·1 (logic), 26·3 (window) |
| WG2 | A deck run in the app is bit-identical to `nitcad run` on it (`deck::run`), run identity included. | 26·1 |
| WG3 | The app's logic under unit tests: the run state machine, document load and save, view settings; deck errors land on their lines. | 26·1 |
| WG4 | `--self-test` passes on WARP and on the GPU and writes the same image as `nitcad plot` for the same view. On CI only if a hosted runner can open a window: the self-test exits 77 when no window or swap chain can be made, and ctest counts that as skipped; the off-screen WARP image tests run in CI regardless. | 26·3 |
| WG5 | Per-monitor DPI re-layout, keyboard access to every command, a clean D3D12 debug layer, no leaked GPU objects. | 26·2 |
| WG6 | Suites as VG6 (Release and Debug at every commit, `/W4 /WX`, hand mutants caught); nothing outside VA1–VA6 and VA8–VA9, nothing on VA7's list; the `app` layer's includes checked against ARCHITECTURE.md 3 (no `physics`, no `assemble`); no numerical change; ARCHITECTURE.md as built after acceptance. | all |
| WG7 | The Device, Mesh, Physics and Simulation menus open the deck at their sections and list their keys for insertion (VQ13); inserting any listed key (with its default value) into a deck leaves a deck that `deck::build` accepts, with the same meaning. | 26·1 (logic), 26·2 (menus) |
| WG8 | Every deck in `examples/` runs to completion in a test, and the start page lists exactly those decks. | 26·3 |

**CI and WG4.** Whether a hosted `windows-2025` runner can create a window and a flip-model swap chain cannot be checked
from this workstation and needs a CI run, which needs a push (not authorized yet). WG4's test is therefore written to
run anywhere: where a window cannot be made it is skipped (exit 77), not failed, and the first CI run after an
authorized push answers the question. The off-screen image tests (`nitcad_render_test`, `nitcad plot`) run in CI as
before.

**Commits (approved, AQ1):** 26·1 the app logic (WG1 logic, WG2, WG3, WG7 logic); 26·2 the main window (WG5, WG7
menus); 26·3 the start page, run monitor, results workspace, `examples/`, `--self-test` (WG1 window, WG4, WG8).

**Commit 26·1 (the app logic; WG1 logic, WG2, WG3, WG7 logic, WG6), as built.** The `app` layer's first files
(`include/NiTCAD/app/`, `src/app/`; CMake target `NiTCAD::app`, linked to deck, render, analysis and io), no window:
- `document.hpp`: `DeckDocument` (open, `set_text` with the modified flag, `save` through a temporary file renamed over
  the deck, `save_as`, the title); `check_deck` (parse and build as the file, so relative paths and messages are the
  deck's); `error_line` (the line after "<source>:" in a message, the form every deck error has).
- `run_controller.hpp`: `RunController` with one `std::jthread` per check or run; the worker builds (`check_deck`) and
  runs (`deck::run` with a `RunControl` on the thread's stop token); events (built, one per Newton iteration with the
  point's wall time, finished) are queued under a mutex and the caller's wake function is called only when the queue
  was empty, so a fast solve posts one message per batch; `process()` applies them on the UI thread. States idle,
  checking, running, stopping, done, failed; a run that stops early is done with its report (`stopped_early`); the
  rows of the run monitor (one per point, its iterations, last update and residual, time) and every iteration for
  the convergence plot; `take_results` hands the plan and report to the workspace.
- `workspace.hpp`: `Workspace` (device and geometry at fixed addresses; results; the bias points `nitcad plot` lists,
  transient snapshots with their currents and rates) from a run, a deck's outputs or a result file with its device
  (and mesh) file, refused on another mesh by the mesh digest or the node count; `run_record` (mesh, points, how it
  ended, input identity, mesh digest, settings, provenance), `convergence_history`, `extracted_figures` (Unit 24: for
  a gate or electrode sweep, the threshold at maximum g_m, the subthreshold swing, the maximum g_m and the on/off
  ratio of the transfer curve of the contact carrying the most current; for an ohmic sweep, the current at its last
  bias).
- `view_settings.hpp`: `ViewSettings` (view, field, point, scale, colormap, range, contours, edges, opaque insulators,
  contacts, projection, azimuth and elevation, slices and the clip box in µm, iso levels, curve contacts and frequency,
  log axis) and `build_scene`, which makes `nitcad plot`'s calls with its choices (the ten-decade floor of vector
  magnitudes, the 0..1 scale of a zero magnitude, the swept contact) without including the command's code (AQ4);
  `available_views`, `available_fields` (the rates only when stored), `plot_arguments` (the command's options for a
  setting, slices and the box converted from µm exactly as the command converts them).
- `deck_sections.hpp` (VQ13): the sections each of the Device, Mesh, Physics and Simulation menus opens, and the keys
  with documented defaults (`[device]` 1, `[mesh]` 3 by the generator, `[models]` 16, `[solver]` 23, `[transient]` 5);
  `find_sections` on text that need not build (comments, quotes and continuations as the syntax has them);
  `insertable_keys` (not those present; a mesh's by its type); `insert_key` (after the section's last entry, or the
  section added at the end, but no second analysis section).
- `app_state.hpp`: recent files (eight, most recent first, no duplicates) and the window's layout as "key = value"
  lines, damaged lines ignored, paths in UTF-8.
- **Tests** (`tests/app/`, 22 cases): errors at their lines (syntax, meaning, untitled); open, save, save as, a failed
  save keeping the path and flag; the state machine (checks, a failed build, refusals while busy, rows equal to the
  sweep's convergence records, outputs written); **WG1** — a 1,601-point run while a hidden message-only window's
  loop pumps a 10 ms timer: 42 timer messages, the longest gap 20–23 ms (bound 100 ms), progress in many batches;
  Stop after three converged points: done, cancelled, the completed points kept and written; **WG2** — the 1D diode
  and the 2D MOSFET run by the worker and by the `nitcad` executable: every field, current, setting, input identity
  and mesh digest equal to the bit, the app's written file too; the workspace from a run, the deck's outputs and the
  result file alike, another mesh refused by digest; the run record and history; the diode's current and the
  MOSFET's threshold, swing, g_m and on/off ratio; ten views (solid, mesh, potential with contours, n with edges and
  viridis at point 1, |J|, isosurfaces, I-V on a log axis, a 1D profile, a band diagram, a diode I-V) drawn on WARP
  **equal to `nitcad plot`'s PNG for `plot_arguments` to the last bit** (max channel difference 0); **WG7** — every
  offered key of every menu inserted into five decks (tensor 1D and 2D sweeps, generated 2D and 3D meshes, a
  transient) builds with an unchanged fingerprint (device temperature and mesh digest, every model, Newton, step,
  linear, equation, current, rate and transient option): 40 to 45 insertions per deck; the recent files and layout
  round trip.
- `app_layering` (ctest, a CMake script): no `physics` or `assemble` include in `app` (WG6).
- **Finding (test fixture, not changed):** a profile's window along its depth axis is ignored
  (`device/doping_profiles.hpp`), and the depth axis defaults to x, so the deck tests' 1D diode
  (`profile = acceptor constant 1e17 x=-1:1um`, `donor ... x=1um:3um`) is uniformly compensated, not a pn junction.
  The app tests copy that diode as a generic device; the 2D MOSFET fixture puts its pockets as Gaussians along y
  with x windows, and the examples (26·3) avoid the pattern.
- **Mutation checks (Release): 15 of 16 caught, 1 equivalent and removed** — error_line ignoring the source (caught once a same-length
  source was added), save keeping the flag, one row for all points, Stop not stopping the worker, a failed run done,
  starting while running, the digest unchecked, no ten-decade floor, edges ignored, `plot_arguments` dropping
  `--edges`, no band view, insertion after the header, present keys offered again, a wrong default, duplicate recent
  files; the sensed contact restricted to ohmic contacts and to contacts other than the swept one survived as
  equivalent (a gate carries no steady current), and the condition was removed.
- **Suites:** the app tests pass in Release (2 s) and Debug (65 s), no warnings under `/W4 /WX`; the whole Release
  suite in a clean worktree at `e20f6ae`: 24/24 (850 s).

**Commit 26·2 (the main window; WG5, WG7 menus, VA1, VA2, VA4, VA8, VA9's executable), as built.**
- **`render` (what a window needs, VA4):** `WindowTarget` and `Renderer::window_target(hwnd, w, h)` (a flip-model DXGI
  swap chain on the renderer's queue), `Renderer::present(scene, target, w, h)` (the scene drawn exactly as `draw`
  draws it, then copied into the back buffer and presented, so the window shows the PNG's pixels; the swap chain
  resized when the size changes); `Camera::zoom` and `Camera::pan` (the fitted view magnified about its centre and
  moved in fractions of its half size; at 1 and 0 the arithmetic is the fitted view's to the bit, and the 20 WARP
  references pass unchanged); the 2D axes span only the part of the bounds in view (the fitted view's ticks
  unchanged); `screen_position(scene, w, h, point)` (the pixel and depth `draw` puts a point at, for picking).
  **Finding (fixed):** releasing a swap chain whose last present was still queued crashed under the D3D12 debug
  layer; a `WindowTarget` now waits for the renderer's queue before it is released.
- **App logic:** `commands.hpp` (every command with its menu, mnemonic, accelerator and its keys' text; the camera's
  keys act in the focused view, not as accelerators, so the deck editor keeps its arrows), `layout.hpp` (the panes
  at a DPI: sizes at 96 DPI scaled to the monitor's; the sides narrowed to keep the centre's minimum), `probe.hpp`
  (the node drawn under a pixel: every node in 2D, the boundary surface's nodes in 3D, the front-most within 12
  pixels; the distance between two), `ViewSettings` zoom and pan with `orbit`, `zoom_by`, `pan_by`, `fit` and
  `apply_camera` (a camera move redraws without rebuilding the scene), `workspace_from_plan` and
  `RunController::take_plan` (a deck opened is built on the worker, not the UI thread).
- **The window (`src/app/`, `nitcad_app_window`; `app/main.cpp`, `nitcad-app.exe`, a `/SUBSYSTEM:WINDOWS` executable
  separate from `nitcad`, AQ4):** per-monitor DPI aware v2, common controls 6; the workspace of the reference: the
  menu bar File, Edit, Device, Mesh, Physics, Simulation, View; a tool bar (New, Open, Save, Check, Run, Stop, Fit,
  Probe, Measure, Perspective); the scene tree (regions, contacts, mesh, views, fields, results); the tabs Device
  view (with the metrics strip: extents, nodes and elements, the field's maximum, the draw time), Curves and Deck;
  the inspector (view, field, point, scale, colormap, contours, mesh edges, contacts, opaque insulators,
  perspective, a slice plane, the clipping box, iso levels); the status bar (state, file, adapter, the cursor's
  node); splitters; the layout and recent files saved per user (`%APPDATA%\NiTCAD\app.ini`). The view
  (`viewport.cpp`) presents in its swap chain: left drag orbits (3D) or pans (2D), right or middle drag pans, the
  wheel zooms, a double click fits, a click probes or measures; keys arrows, Shift+arrows, + and -, Home. The deck
  pane (`deck_pane.cpp`): Rich Edit with a line-number gutter (error lines in red), the check's errors below (a
  double click or Enter goes to the line). Opening a deck checks it on the worker and shows its device (and the
  results its outputs hold); Run (F5) runs it on the worker (WG1), the status bar following each point, and shows
  the potential and the run's curve (the sensed contact the one carrying the most current, on a log axis). The
  Device, Mesh, Physics and Simulation menus (VQ13, WG7) list "Go to [kind]" for their sections (greyed when the
  deck has none) and "Insert into [kind]" with the insertable keys and their help.
- **Tests:** `commands_test` (unique mnemonics per menu, unique accelerators, the reference's menu bar),
  `layout_test` (the panes tile the window; at 192 DPI every pane doubles; the centre's minimum), `probe_test`
  (2D nodes found where drawn with their values; nothing out of reach; distances; in 3D the front-most of two
  nodes at one pixel; orbit, zoom, pan and fit limits; zoom and pan about the centre, both axes),
  `render/window_test` (screen positions in the layout, zoom and pan; zoom 1 is the fitted image to the bit; a
  panned view leaves the margin white; a swap chain presents field and plot scenes and follows a resize), and
  `nitcad_app_window_test` (ctest `app_window`, hidden windows on WARP): **WG5** — every command in the built
  menus and its keys in the accelerator table; a `WM_DPICHANGED` to 192 DPI re-lays the panes at the new scale
  (the children where `layout_panes` puts them, 9-point fonts at 192 DPI) and back; **WG7** — the Physics menu
  inserts `srh = on` into a new `[models]` at the caret and the deck still builds, Simulation's "Go to [sweep]"
  moves the caret to its header, a deck error (F7) lands the caret on its line; a deck opened, run (worker),
  shown, probed, zoomed, panned and fitted, the curves tab drawing I-V; **WG5** — with the D3D12 debug layer:
  no message after opening, running, zooming and switching tabs, and no Direct3D or DXGI object alive once the
  window is destroyed (the DXGI debug interface's live-object report, which sees the renderer's objects while
  the window exists).
- **Looked at:** `nitcad-app` on the 2D MOSFET deck on the RTX 5060 Ti: the regions and then the potential after
  Run, the layout as the reference (tree, tabs, inspector, metrics strip, status bar); the view drawn in 3.8 ms.
- **Mutation checks (Release): 21 of 21 caught** — zoom ignored, pan x reversed, pan y by the half width, the 2D
  axes not clamped, screen positions without the layout, the elevation and the zoom unclamped, Fit resetting the
  angles, pan y not applied (caught once a vertical pan was tested), the back-most node picked (caught once two
  nodes at one pixel were tested), no pick radius, sizes not scaled by the DPI, the centre not kept at its minimum
  (caught once a 600-pixel window was tested), a duplicate accelerator, fonts not made at the new DPI, the new DPI
  ignored, accelerators without Shift, an inserted key not shown at its line, the curve sensing the swept gate, a
  swap chain released before its present ran, Run doing nothing. Three mutants first failed to build (a `/WX`
  warning, a delimiter in the substitution) and were rewritten so they compile.
- **Suites:** the app, window and render tests pass in Release and Debug (the Debug render suite with the D3D12
  debug layer clean), no warnings under `/W4 /WX`.
- **Not in 26·2 (26·3):** the start page, the run monitor, the results page, `examples/`, `--self-test`.
- **Suites:** the whole Release suite in a clean worktree at `ea90185`: 25/25.

**Commit 26·3 (the start page, run monitor, results workspace, examples and `--self-test`; WG1 window, WG4, WG8,
VA3, VA5, VA9), as built.**
- **`examples/` (VQ11, AQ3):** four decks, each titled and described by its first two comment lines:
  `pn_diode_1d.deck` (an erfc n+ profile on a 1e16 p body, swept forward to 0.7 V: the band diagram and I-V),
  `mosfet_2d.deck` (a 2 × 1 µm n-channel MOSFET under 50 nm of oxide, Gaussian source and drain pockets, Id-Vg at
  V_D = 0.1 V, 16 points), `mosfet_3d.deck` (the same in 3D, 0.8 µm wide, a gate sweep at V_D = 0.5 V, 5 points)
  and `moscap_cv.deck` (a 10 nm n+ poly gate on 1e17 p silicon, small-signal C-V at 1 MHz). Run with `nitcad`:
  0.04 s, 0.7 s, 4.4 s and under 0.1 s (Release). Looked at: the diode's junction near 0.9 µm with the
  quasi-Fermi levels split by the bias; the MOSFET's Id-Vg from 5e-9 to 6e-2 A/cm with its subthreshold slope;
  the MOS-C falling from the oxide capacitance (3.3e-7 F/cm²) through depletion; the 3D MOSFET's inverted
  channel. **Deviation (the MOS-C):** the sweep stops at 0.25 V, the onset of inversion: from there the
  drift-diffusion steady state of this 1D device did not converge (Newton at 0.5 V, also with step cutting and a
  line search, and on a mesh graded to 0.5 nm at the oxide), its electrons coming only from generation. The deck
  says so.
- **App logic:** `examples.hpp` (`list_examples`, `examples_directory` beside the executable else the source tree's,
  `copy_example` into the user's Documents\NiTCAD\examples, a copy kept once made) and `run_monitor.hpp` (the
  monitor's columns and cells, the convergence plot of every positive update on a log axis, the progress as
  converged bias points over the planned ones).
- **The window:** tabs Start, Device view, Curves, Deck, Run monitor, Results (Ctrl+1 to Ctrl+6); the start page
  (New deck, Open deck, Open result, the recent files, the examples with their descriptions; an example opens as
  the user's copy); the run monitor (shown when Run starts: a virtual list of the rows, the progress bar, Stop, a
  status line with how the run ended, the live convergence plot in its own swap chain, made before the run starts
  and redrawn at most every quarter second and at most a twentieth of the time); the results page (Unit 24's
  figures with their curve and flags, the run record and provenance).
- **`nitcad-app --self-test <deck> --out <png> [view options]` (VA6):** opens the window (shown), checks and runs the
  deck on the worker while its message loop runs with a 10 ms timer (the longest gap must stay under 100 ms,
  WG1), shows the view in the window's swap chain, writes the same view's image at `--size` and closes; with
  `--debug-layer`, also no debug message and no live graphics object after the window. Exit 0, 1, 64 (usage) or
  77 when no window or swap chain can be made (CI, WG4).
- **Finding (fixed):** the run monitor's plot made its swap chain on its first draw, during the run; on the GPU
  that took about 130 ms of the UI thread (the self-test's WG1 measure caught it). The swap chain is now made when
  the monitor is shown, before the run starts: the longest gap is 19–32 ms on the GPU and 22–30 ms on WARP for
  the 3D example, with the plot drawn in 2–3 ms.
- **Finding (fixed, Debug):** in Debug the solver's events arrived as fast as the slower window took them, so its
  queue was never empty and the timer (and input) waited behind the posted messages: gaps of 2.2 to 3.2 s. The
  window now shows a run's progress at most every 50 ms (and on a 100 ms timer while the events pause), so taking
  an event is cheap: Debug gaps 39 ms. The self-test measures the gap while the solver runs, and reports apart
  the loading of the results and the first view when the run ends (on the UI thread: 63 ms for the 2D MOSFET and
  95 ms for the 3D one in Release, 110 ms and 809 ms in Debug).
- **Tests:** `examples_test` (**WG8**: exactly the four decks, titled and described; a deck without comments; the
  directory beside the executable; a user's copy kept; every example run by the controller to completion, its
  points all solved and its outputs written), `run_monitor_test`, and in `nitcad_app_window_test`: the start page
  lists exactly the examples (the bundled ones and a given directory); an example opens as a copy; Run shows the
  monitor, which fills while the run goes on; the results page holds the threshold and the record; **WG1** — Stop
  after three points ends the run early with the completed points kept, shown and the status saying why; every
  tab reachable from the keyboard; **WG4** (`self_test_test`) — `nitcad-app --self-test` and `nitcad plot` on the
  2D MOSFET (potential; n with viridis and contours; Id-Vg on a log axis) and the 3D MOSFET (n) give the same PNG
  to the bit on WARP, and on the RTX 5060 Ti; the self-test with the debug layer passes with nothing alive after
  the window; usage errors exit 64. A session that cannot open a window skips WG4 (exit 77): on the hosted CI
  runner this is still unverified (it needs an authorized push), and the off-screen WARP tests run there either
  way.
- **Mutation checks (Release): 14 of 14 caught** — the title not read from the deck, the description from the
  first line, the examples beside the executable ignored, a user's copy overwritten, unconverged points counted as
  done, zero updates plotted on the log axis, the point numbered from 0, an example opened in place, Run not
  showing the monitor, the monitor updated only at the end, the start page missing an example (the bundled ones;
  a given directory, caught once a test gave one), the self-test's image not the plot's size, the self-test
  ignoring the view's options (caught once a view other than the window's default was compared).
- **Observed once, not reproduced:** while two builds, a mutation run and a suite shared the machine, one GPU
  self-test image differed from `nitcad plot`'s (the files left behind were equal) and WARP self-tests saw 400 to
  880 ms gaps; on reruns the GPU and WARP comparisons were exact and the gaps 19–32 ms.
- **Suites:** the app tests (32 cases) and the window tests (8 cases, the self-test included) pass in Release and
  Debug, no warnings under `/W4 /WX`. In a clean worktree at `84c34e9`: Release 25/25, Debug 19/19 — but the
  Debug suite took 2,645 s against about 1,460 s before Unit 26, because the app and window tests ran the 3D
  MOSFET example in Debug (about 5 minutes a run: `app` 369 s, `app_window` 982 s), too long for the CI job's 60
  minutes.

**Follow-up (the 3D example's gates in the Release suite).** As the other slow gates (`solve_mosfet`, ...), the 3D
MOSFET's run (WG8) and its self-tests on WARP and the GPU (WG4) moved to hidden-tag cases run by the Release-only
ctest entries `app_examples_3d` and `app_window_3d`; the Debug suite runs the other three examples and the 2D
self-tests (WARP, the GPU, the debug layer). Debug: `app` 101 s, `app_window` 231 s; Release: `app` 3.6 s,
`app_window` 8.6 s, `app_examples_3d` 4.2 s, `app_window_3d` 10.7 s.

**Suites at every commit (clean worktrees):** `e20f6ae` (26·1) Release 24/24, Debug 18/18; `ea90185` (26·2) Release
25/25, Debug 19/19; `84c34e9` (26·3) Release 25/25, Debug 19/19. Nothing outside `render` (the window target, the
camera's zoom and pan, the 2D axes in view, screen positions), `app`, `app/`, `examples/`, tests, `CMakeLists.txt`
and this file differs from `main` (`e289e44`): no numerical change, and the 20 WARP references pass unchanged.
**Status:** Unit 26 built (WG1–WG8 met on this workstation; WG4 on the hosted CI runner unverified until an
authorized push). **Accepted by the owner 2026-10-22** ("Accept Unit 26"); ARCHITECTURE.md 6.14, 8, 11 and 11.3
updated to "as built" (WG6). Merged to `main` as `e9b1d71` (PR #59, 2026-10-22).
