<!-- Archived from DECISIONS.md on 2026-11-17 (branch docs/archive-17-26d), word for word; DECISIONS.md keeps its summary. -->

### Unit 26d — GUI quick wins (proposed 2026-11-15; on `app/gui-quick-wins`)

**Owner request (2026-11-15):** "Check the NiTCAD Interface (GUI), suggest improvements to it", then "which all are
small", then "yes write proposal". This unit covers the small items from that review. The branch starts from `main`
`dae4568`.

**Standing conditions:** the app only, with no change to numerics, results, files, run identities or `nitcad plot`
output. Win32 + D3D12 + Direct2D only. Criteria are fixed before code. Release and Debug suites run at every commit
with `/W4 /WX`, plus hand mutants. WG1 (message-loop gaps) and the self-tests stay as ruled. Nothing is pushed or
merged without authorization.

**Left out of this unit (from the same review, medium or larger):** show only the inspector controls that apply,
merge the two clipping-box inputs, image export / copy to the clipboard, TaskDialog errors, toolbar icons, tooltips
and mode grouping, right-click menus, tree/inspector consolidation, camera presets, full dark mode, and undo for view
settings. The duplicated "ψ (V)" heading also stays: it is the scene title that `nitcad plot` draws too, so changing
it would change PNG output.

**Design:**

| # | Change | Where | Behaviour |
|---|---|---|---|
| ZD1 | **Remove Apply.** The inspector's text boxes (arrow and streamline counts, seeds, exploded gap, slice at, iso levels, cut line, clipping box, and ZD3's range) take effect on Enter, and when focus leaves a box whose text changed. | `main_window.cpp` `on_inspector` (Apply branch made a function), `pre_translate` (Enter in an inspector edit), `EN_KILLFOCUS` | Same validation and messages as Apply today. Invalid text gets the status-bar message and the box shows the setting again. Each box takes effect separately (today one Apply commits all of them). The Apply button and its row are removed. |
| ZD2 | **Point labels.** The Point drop-down shows `k · <contact> = <V> V` for each contact whose bias varies over the points (the swept contact first), `k · t = <time>` for a transient (`render::time_label`), and the bare index only when nothing varies. | New `point_labels(const Workspace&)` in `view_settings` (testable without a window); `sync_inspector` uses it | Selection by index is unchanged. Volts use `{:.4g}`. |
| ZD3 | **Colour range.** A "Range" row: min and max boxes and an **Auto** checkbox. Auto (the default) clears `ViewSettings::range`. Typed values set it in the colour bar's units, exactly as `nitcad plot --range` takes them. | Inspector (field, slice, iso and exploded views; hidden on others); `on_inspector`; `sync_inspector` | Validation: min < max, both finite, and both > 0 on a log scale; otherwise a status message. Changing the field resets the range to Auto, as it does today. Inherits ZD1's Enter and focus-loss behaviour. |
| ZD4 | **Drag and drop.** Dropping one file on the window opens it the way the recent-files list does: a `.deck` with `open_deck`, any other file as a result with its sibling `.nitdev` (and `.nitmesh` if present), after `confirm_discard`. More than one file, or a folder, is refused with a status message. | `DragAcceptFiles`; `WM_DROPFILES` in `handle`; the recent-file branch of `execute` made a shared `open_path` | Refused while a run is in progress (as Open is). |
| ZD5 | **Taskbar progress.** During a run, the taskbar button shows the run controller's progress fraction (`ITaskbarList3`: normal while running, error state on failure until the next run or check, cleared on success or stop). When a run ends with the window not in the foreground, the taskbar button flashes until activated (`FlashWindowEx`). | Created on the `TaskbarButtonCreated` message; updated from `on_run_timer` and `on_run_finished` | COM is already initialized STA in `app/main.cpp`. Without a taskbar interface (e.g. the hidden test window), nothing is shown and nothing fails. |
| ZD6 | **Dark title bar.** The title bar follows Windows' app theme (`AppsUseLightTheme`) through `DWMWA_USE_IMMERSIVE_DARK_MODE`, updated on `WM_SETTINGCHANGE` "ImmersiveColorSet". | `create`, `handle` | Title bar only; the client area is unchanged. Where the attribute is unsupported, the call fails silently. |
| ZD7 | **Unit in the metrics strip.** "Max potential 1.714" becomes "Max potential 1.714 V", and the label is the field's symbol rather than its key ("Max ψ 1.714 V"). Both come from `render::field_info`'s title (the symbol, and the unit in its brackets), the text the tree and colour bar show. "Max \|B − A\|" likewise. | `update_metrics` | Text only. |
| ZD8 | **Viridis as the app's default colormap** (added by the owner 2026-11-15: "I don't like the 3D device color", answered "the real app's look", "Viridis"). New views and the Colormap drop-down start at viridis; turbo, cividis, coolwarm and soft stay selectable. This supersedes VQ5's turbo default **for the app only**: `nitcad plot` keeps turbo, so its PNGs are unchanged (ZG8). | `ViewSettings::colormap` default; `colormap_names` order (viridis first); `view_settings.cpp` difference view and equivalent command | The difference view still uses cool-warm unless the user picks a colormap: the sentinel compares to the app default rather than the literal "turbo". The view's equivalent `nitcad plot` command now includes `--colormap viridis`, so pasting it reproduces the app's image. |

**Questions (each with a recommendation):**

| # | Question | Options | Recommendation |
|---|---|---|---|
| ZQ1 | When a ZD1 box takes effect | (a) Enter and focus loss. (b) Enter only. (c) On every keystroke (redraws while typing). | **(a)**: the Windows convention for property panes. (c) would redraw 3D scenes for half-typed numbers. |
| ZQ2 | Point label when several contacts vary | (a) All varying contacts. (b) The swept contact only. | **(a)**, swept first: a 2D sweep (V_G at each V_D) is ambiguous otherwise. |
| ZQ3 | Taskbar error state after a failed run | (a) Kept until the next run or check. (b) Cleared when the window is activated. | **(a)**: a user returning to the window still sees that the run failed. |
| ZQ4 | Dark title bar | (a) Follow the system theme. (b) Always light (today). (c) A View menu choice. | **(a)**: no new setting, and it matches other Windows apps. |

**Acceptance criteria (frozen before code; tests in `tests/app`, hidden windows as today):**

| # | Criterion |
|---|---|
| ZG1 | No Apply button exists. For each ZD1 box: typing a valid value and sending Enter changes the settings and the scene; moving focus does the same; an invalid value leaves the settings unchanged, puts the message in the status bar and shows the old setting in the box. |
| ZG2 | `point_labels`: a single-contact sweep, a two-contact family (both contacts listed, swept first), a transient (times), and a run where nothing varies (index only), each against expected strings. |
| ZG3 | Range: Auto clears it; a typed min/max equals what `nitcad plot --range` gives (the scene's colour bar scale is compared); invalid or non-positive-on-log values are refused; changing the field resets to Auto. |
| ZG4 | Drop: a `.deck` opens it; a result file opens it with its `.nitdev`; two files and a folder are refused; a drop during a run is refused; unsaved edits ask first. Tested by posting `WM_DROPFILES` with a built `HDROP`. |
| ZG5 | Taskbar: the state machine (idle → running with fraction → error / none) is a pure function, tested for each run outcome. A missing taskbar interface does not fail. |
| ZG6 | The dark-mode flag is read from the theme value (a pure function, tested); the DWM call is made once at creation and once per theme change. |
| ZG7 | Metrics text includes the field's unit for every field the tree lists, and for the difference view. |
| ZG8 | No change elsewhere: every existing render and plot test passes unchanged, and app tests change only where ZD1–ZD8 change what they check. `nitcad plot` PNG digests are bit-identical. WG1 self-tests pass at their ruled bounds. |
| ZG9 | Hand mutants (at least one per ZD item, ZD8 a difference view taking viridis; dropped Enter handling, a swapped label order, a range not reset on field change, a drop accepting two files, error state cleared early, an inverted theme flag, a missing unit) are each killed by a test. |
| ZG10 | A new device view and the curves view start with viridis; the difference view starts with cool-warm and takes a chosen colormap; the equivalent command of a default field view contains `--colormap viridis` and, run through `nitcad plot`, gives the app's colour bar scale and map. |

**Slices (each one local commit):** 26d·1 ZD1 + ZD3 (the inspector); 26d·2 ZD2 + ZD7 + ZD8 (labels, text and the default colormap); 26d·3 ZD4 +
ZD5 + ZD6 (window integration); 26d·4 suites, mutants, the screenshot `docs/images/app.png` refreshed, and the
ARCHITECTURE "as built" entry.

**Owner ruling (2026-11-15):** "choose recommended": ZQ1 (a), ZQ2 (a), ZQ3 (a), ZQ4 (a). Implementation starts only
on the owner's word.

**Implementation (owner "go 26d", 2026-11-15).** The branch was rebased onto `main` `beb494f` (PR #84's MKL 2026.2)
before the first slice. Slices: 26d·1 `640f555` (ZD1, ZD3), 26d·2 `6cd8f8b` (ZD2, ZD7, ZD8), 26d·3 `dcdff44` (ZD4–ZD6,
new `app/integration`), 26d·4 this record, ARCHITECTURE's "as built" entry and `docs/images/app.png`.

Choices made in building, within the design:
- **ZD1:** a box takes effect on focus loss only when its edit's change flag is set (text typed, not text the program
  put there), so moving through the inspector does not redraw. The slice axis list commits the slice box and the
  clipping-box check commits the clipping box; an unchecked box keeps its text to be checked again. Where Apply
  silently dropped invalid slice, clipping-box and iso-level text, there are now status messages ("A slice plane is
  at one position in µm.", "A clipping box is <x0>,<y0>,<z0>:<x1>,<y1>,<z1> in µm.", "Iso levels are numbers
  separated by commas.").
- **ZD3:** with Auto on, the boxes show the scale the view was drawn on (`{:.6g}`), so one box can be typed and the
  other taken as shown; unchecking Auto keeps the drawn scale as the range; both boxes emptied mean Auto. A redraw
  while typing (playback) leaves a focused, changed box as typed.
- **ZD8:** on the difference view, viridis (the app's default) and turbo (`nitcad plot`'s, which `nitcad plot` draws
  cool-warm there) both mean cool-warm (`takes_coolwarm`), so `plot_arguments` gives a command that draws the app's
  map without changing `nitcad plot`. Choosing viridis or turbo there therefore shows cool-warm; cividis, coolwarm and
  soft are used as chosen. The curve views name no colormap (they have no colour scale). The self-test takes
  `nitcad plot`'s options, so its default colormap is turbo, keeping "same options, same image" (WG4).
- **ZD5:** the error state is for a failed run only; a check that fails (the deck does not build) is shown in the
  deck, not on the taskbar. A run stopped by Stop ends done and clears the state.
- **Tests changed (ZG8):** the 26b/26c window tests drive Enter instead of Apply; the app and window tests that
  build their own `nitcad plot` commands for the app's default settings name `--colormap viridis` (the WG4 argument
  lists, the film frames, the playback frames).

Results: Release app tests 66/66 test cases, window tests 24/24; 13 hand mutants (ZG9: Enter dropped, focus loss
ignoring the change flag, swapped label order, range kept on a field change, log range not checked, a drop
taking two files, a drop during a run, error state cleared while idle, inverted theme flag, theme applied on any
setting change, missing unit, difference taking viridis, turbo default) all killed. Full Release ctest 227/227
(58.5 min, -j 6), the render and `nitcad plot` reference tests unchanged (ZG8) and the WG1 self-tests at their ruled
bounds; the Debug suite left to CI, as for Unit 17c.

---

