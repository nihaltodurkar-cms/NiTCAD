# MOSFET Guide (Unit 26p) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A View > Guide pane that walks a learner through thirteen MOSFET lessons (basics to research level), each step driving NiTCAD's existing parameter edits, sweeps, runs and tabs and reporting the solver's own figures.

**Architecture:** Lessons are example decks whose leading comments carry the steps; a Win32-free model (`app/lesson.hpp`) parses them and formats reports; a Win32 pane (`src/app/guide_pane.*`) shows them and calls a `Host` that `MainWindow` implements with the paths it already has (`DimensionsPanel`-style `set_dimension`, `start_deck_run`, `select_tab`). The only engine change is `[sweep]` taking `$name`.

**Tech Stack:** C++23, MSVC 14.51, CMake presets, Catch2, Win32 + Direct2D.

**Spec:** `DECISIONS.md`, Unit 26p (GU1–GU6, GUC1–GUC6, GUQ1–GUQ3).

## Global Constraints

- Every deck that builds today builds to the same bits (GUC1).
- Win32 + D3D12 + Direct2D only, no new dependency, no Python.
- The models (`include/NiTCAD/app/*.hpp`, `src/app/lesson.cpp`) have no Win32 type.
- Deck edits go through the deck editor path with one undo per edit; a lesson opens as a copy (`copy_example`), never in place.
- No material value the repository does not already source; lessons use silicon and SiO2 as shipped.
- Builds `-j 6`, runs on 6 threads; targeted Release tests only; no full suite, no Debug build without the owner's word.
- No screen automation. Commit locally per task on `app/mosfet-guide`; nothing pushed.
- Figure names are the app's (`src/app/workspace.cpp` `extracted_figures`): `Threshold (max gm)`, `Subthreshold swing`, `Max transconductance`, `On/off ratio`, `Current at the last bias`.

## Review Focus

1. A lesson deck the learner has edited by hand so a `# do: set` names a parameter they deleted: the step shows the `set_dimension` error and stays (Task 6 test `guide action error keeps the step`).
2. A sweep whose runs partly fail: the report shows `—` for the failed values, not a crash or a shifted list (Task 3 test `report with a failed run`).
3. Toggling Guide while a run is busy: actions are disabled until the run ends, Back/Next still work (Task 6 test `guide actions wait for the run`).
4. The pane at 150 % and 200 % DPI and in the Aurora theme: text wraps, nothing clipped (Task 5 test `guide pane relayout at dpi`).
5. A learner opening a second lesson while the first's copy has unsaved edits: the existing `confirm_discard()` prompt runs first (Task 6 test `guide lesson switch asks to discard`).

---

### Task 1: `[sweep]` takes parameters (GU1)

**Files:**
- Modify: `src/deck/parameters.cpp:265-270` (the `takes` set and its refusal text)
- Test: `tests/deck/parameters_test.cpp`

**Interfaces:**
- Produces: decks may write `contact <name> = $Vg` and `contact drain = 0 : $Vd : 0.05` in `[sweep]`.

- [ ] **Step 1: Write the failing tests** in `tests/deck/parameters_test.cpp`:
  - `TEST_CASE("a number parameter substitutes in [sweep]")`: deck with `[parameters] Vg = 0.8`, `[sweep] contact gate = $Vg`; `deck::substituted` succeeds and the gate entry's value is `0.8`.
  - `TEST_CASE("a range in [sweep] takes parameters")`: `contact drain = 0 : $Vd : 0.05` with `Vd = 1.2` gives `0 : 1.2 : 0.05`.
  - `TEST_CASE("a missing parameter in [sweep] is refused on its line")`: `$nope` → error whose message contains the line number and `nope`.
  - Update the existing refusal case near line 224 to the new text `outside [mesh], [region], [contact], [doping] and [sweep]` and make it use `[solver]`.
- [ ] **Step 2: Run** `cmake --build --preset release -j 6 --target nitcad_deck_test && ctest --preset release -R "^deck\.parameters$"` — expected: the three new cases FAIL.
- [ ] **Step 3: Implement**: add `s.kind == "sweep"` to `takes`; refusal text `"a parameter reference outside [mesh], [region], [contact], [doping] and [sweep]"`.
- [ ] **Step 4: Run** the same command plus `-R "^deck\."` — expected: all deck tests PASS (GUC1: unchanged decks unaffected).
- [ ] **Step 5: Commit** `deck: [sweep] takes $name parameters (Unit 26p GU1)`.

### Task 2: The lesson model (GU2)

**Files:**
- Create: `include/NiTCAD/app/lesson.hpp`, `src/app/lesson.cpp`
- Modify: `CMakeLists.txt` (add `src/app/lesson.cpp` to `NiTCAD::app` beside `src/app/sweep_results.cpp:1142`; add `tests/app/lesson_test.cpp` to `nitcad_app_test` beside `tests/app/sweep_results_test.cpp`)
- Test: `tests/app/lesson_test.cpp`

**Interfaces:**
- Produces (namespace `NiTCAD::app`):
```cpp
enum class LessonTrack : unsigned char { basics, intermediate, advanced };
struct LessonAction {
    enum class Kind : unsigned char { none, set, sweep, run, edit } kind = Kind::none;
    std::string section;  // edit only, without brackets
    std::string name;     // set/sweep: parameter; edit: key
    std::string value;    // set: value; sweep: DT6 list as written; edit: value
};
struct LessonStep {
    std::string heading;
    std::string prose;                 // paragraphs joined by "\n\n", lines within a paragraph by ' '
    LessonAction action;
    std::string show;                  // "", "geometry", "fields", "curves", "results"
    std::vector<std::string> report;   // figure names, trimmed
    int line = 0;                      // the `# step:` line, 1-based
};
struct Lesson {
    std::string title;
    LessonTrack track = LessonTrack::basics;
    int order = 0;
    std::vector<LessonStep> steps;
    std::filesystem::path path;        // set by list_lessons
};
struct LessonError { int line = 0; std::string message; };  // line 0: not a lesson
[[nodiscard]] std::expected<Lesson, LessonError> parse_lesson(std::string_view deck_text);
[[nodiscard]] std::vector<Lesson> list_lessons(const std::filesystem::path& directory);  // *.deck that parse, by track then order then title
[[nodiscard]] std::string_view to_string(LessonTrack t) noexcept;  // "Basics", "Intermediate", "Advanced"
[[nodiscard]] std::expected<std::string, std::string> set_entry(std::string_view deck_text, std::string_view section,
                                                                std::string_view key, std::string_view value);
```
- Rules: only the leading comment block (lines starting `#` before the first non-comment, non-blank line) is read. Keys `lesson`, `track`, `order` before the first `# step:`. `# do:`, `# show:`, `# report:` belong to the current step; a second `# do:` in one step is an error. `set`/`sweep` names must be parameters (`read_dimensions`); `sweep` values must pass `deck::sweep_values` with the parameter's kind. `set_entry` = `find_section` + `with_value` on the key's line; error when the section or key is missing.

- [ ] **Step 1: Write the failing tests** (`tests/app/lesson_test.cpp`, tag `[lesson]`):
  - `parse a full lesson`: a two-step deck text → title, `LessonTrack::advanced`, order 10, step 1 prose with one paragraph break equals `"A b.\n\nC d."`, step 2 action `{sweep, "", "Na", "1e16, 1e17, 1e18"}`, `show == "results"`, `report == {"Threshold (max gm)", "Subthreshold swing"}`.
  - One case per action form: `set Na = 1e17`, `run`, `edit [models] field_mobility = off` (section `models`).
  - Errors with lines: unknown action `# do: dance` (its line), action before the first step, two `# do:` in one step, `sweep Lx = ...` with no parameter `Lx`, an unknown track, a non-lesson deck → line 0.
  - `list_lessons` over a temp dir with three files (advanced/1, basics/2, basics/1) and one non-lesson → order basics/1, basics/2, advanced/1.
  - `set_entry` on `[models]\nfield_mobility = on  # x` → `field_mobility = off  # x`; missing key → error.
- [ ] **Step 2: Run** `cmake --build --preset release -j 6 --target nitcad_app_test && ctest --preset release -R "^app\.lesson$"` — expected FAIL (missing header).
- [ ] **Step 3: Implement** the header and `src/app/lesson.cpp` per the rules above.
- [ ] **Step 4: Run** the same — expected PASS.
- [ ] **Step 5: Commit** `app: lesson model (Unit 26p GU2)`.

### Task 3: Reports (GU5)

**Files:**
- Modify: `include/NiTCAD/app/lesson.hpp`, `src/app/lesson.cpp`
- Test: `tests/app/lesson_test.cpp`

**Interfaces:**
- Consumes: `SweepResults`, `SweepFigures`, `sweep_figures`, `run_label` (`app/sweep_results.hpp`); `Figure`, `extracted_figures` (`app/workspace.hpp`).
- Produces:
```cpp
// One line per figure: "<name>: v1 → v2 → v3 <unit>" in the sweep's order, "—" for a run without it.
[[nodiscard]] std::string sweep_report(const SweepResults& s, const SweepFigures& f, const std::vector<std::string>& figures);
// One line per figure: "<name>: before → after <unit>"; a figure only in `after` reads "<name>: after <unit>".
[[nodiscard]] std::string run_report(const std::vector<Figure>& before, const std::vector<Figure>& after,
                                     const std::vector<std::string>& figures);
```
- Number format: 3 significant figures (`std::format("{:.3g}")`); on/off ratio as `1.2e+06`.

- [ ] **Step 1: Write the failing tests**: `sweep report lists values in order` (`Threshold (max gm): 0.31 → 0.58 → 0.92 V`), `report with a failed run` (middle run not loaded → `0.31 → — → 0.92 V`), `run report before and after` (`Subthreshold swing: 88 → 71 mV/dec`), `run report without a previous run`, `a figure no run has` → `<name>: —`.
- [ ] **Step 2: Run** `ctest --preset release -R "^app\.lesson$"` after building `nitcad_app_test` — expected FAIL.
- [ ] **Step 3: Implement** both functions.
- [ ] **Step 4: Run** — expected PASS.
- [ ] **Step 5: Commit** `app: lesson reports (Unit 26p GU5)`.

### Task 4: The guide's place in the layout and the command (GU3)

**Files:**
- Modify: `include/NiTCAD/app/layout.hpp` (`Panes` gains `Rect guide, guide_splitter`; `guide_width_96 = 320`), `src/app/layout.cpp`, `include/NiTCAD/app/app_state.hpp` (`Layout` gains `bool guide_shown = false; int guide_width = 0;`), `src/app/app_state.cpp` (save/load the two fields; an old state file without them loads the defaults), `include/NiTCAD/app/commands.hpp` (append `guide` after `tab_geometry`, comment `// Unit 26p (GU3)`), `src/app/commands.cpp` (`{Command::guide, "&View", "&Guide", {}, "F1"}` with accelerator F1)
- Test: `tests/app/layout_test.cpp`, `tests/app/app_state_test.cpp`, `tests/app/commands_test.cpp`

**Interfaces:**
- Produces: `layout_panes(...)` places `guide` at the client's right edge (width from `Layout::guide_width`, else `guide_width_96` scaled), the inspector to its left, only when `saved.guide_shown`; otherwise `guide` is empty and every other rect is what it is today.

- [ ] **Step 1: Write the failing tests**: `guide hidden leaves the panes as before` (compare all rects with `guide_shown=false` against today's values for 1600×900 @96 and @144), `guide shown takes the right edge` (`guide.right() == client_width`, `inspector.right() == guide_splitter.x`), `layout state round-trips guide fields`, `an old state file loads guide defaults`, `guide command is in View with F1`.
- [ ] **Step 2: Run** `ctest --preset release -R "^app\.(layout|app_state|commands)$"` — expected FAIL.
- [ ] **Step 3: Implement**.
- [ ] **Step 4: Run** — expected PASS.
- [ ] **Step 5: Commit** `app: guide pane layout and View > Guide (Unit 26p GU3)`.

### Task 5: The guide pane (GU3)

**Files:**
- Create: `src/app/guide_pane.hpp`, `src/app/guide_pane.cpp`
- Modify: `CMakeLists.txt` (add to `nitcad_app_window` after `src/app/sweep_dialog.cpp`; add `tests/app/guide_window_test.cpp` to `nitcad_app_window_test`)
- Test: `tests/app/guide_window_test.cpp`

**Interfaces:**
- Consumes: `Lesson`, `LessonStep`, `LessonAction`, `to_string(LessonTrack)` (Task 2).
- Produces:
```cpp
class GuidePane {
public:
    struct Host {
        virtual void guide_open(std::size_t lesson) = 0;              // the picker's choice
        virtual void guide_act(const LessonStep& step) = 0;           // the action button
    protected:
        ~Host() = default;
    };
    GuidePane(HWND parent, Host& host, int id);
    [[nodiscard]] HWND hwnd() const noexcept;
    void set_lessons(std::vector<Lesson> lessons);
    void show_lesson(std::size_t lesson);                              // step 0
    void go(int delta);                                                // Back -1 / Next +1, clamped
    void set_report(std::string text);                                 // "What happened"
    void set_error(std::string text);                                  // shown in place of the report, step kept
    void set_busy(bool busy);                                          // action disabled while a run is busy
    void set_dpi(int dpi, HFONT font, HFONT bold);
    [[nodiscard]] std::size_t lesson() const noexcept;
    [[nodiscard]] std::size_t step() const noexcept;
    [[nodiscard]] std::wstring action_label() const;                   // "Sweep Na = 1e16, 1e17, 1e18", "Run", "Set tox = 2nm", "Set field_mobility = off"; empty when none
};
```
- Controls (child windows, as `SweepDialog` does): combo of lessons grouped `Basics — 1. The gate and the channel`; static heading (bold); read-only multiline edit for prose (wraps, scrolls); action button; read-only edit for the report; Back, Next, `Step n of m` static. Tab order in that sequence; Enter on the action button activates it. Colours from the window's theme like the other panes.

- [ ] **Step 1: Write the failing tests** (tag `[guide_window]`, the existing window harness): `guide pane shows the first step`, `go clamps at both ends`, `action label per kind`, `action button calls the host with the step`, `set_error keeps the step`, `set_busy disables the action only`, `guide pane relayout at dpi` (at 144 and 192 DPI the prose edit's width equals the pane's client width minus scaled margins; no child exceeds the pane).
- [ ] **Step 2: Run** `cmake --build --preset release -j 6 --target nitcad_app_window_test && ctest --preset release -R "^app_window\.guide_window$"` — expected FAIL.
- [ ] **Step 3: Implement**.
- [ ] **Step 4: Run** — expected PASS.
- [ ] **Step 5: Commit** `app: guide pane (Unit 26p GU3)`.

### Task 6: Wiring into the main window (GU3, GU4, GU5)

**Files:**
- Modify: `src/app/main_window.hpp` (`MainWindow` also derives from `GuidePane::Host`; members `std::unique_ptr<GuidePane> guide_`, `std::vector<Lesson> lessons_`, `std::vector<Figure> guide_before_`, `const LessonStep* guide_pending_ = nullptr`), `src/app/main_window.cpp`
- Test: `tests/app/guide_window_test.cpp`

**Interfaces:**
- Consumes: `GuidePane` (Task 5); `parse_lesson`, `list_lessons`, `set_entry`, `sweep_report`, `run_report` (Tasks 2–3); `set_dimension` (`app/dimensions.hpp`); existing `confirm_discard()`, `copy_example`, `open_deck`, `geometry_edited(text)`, `start_deck_run(...)`, `select_tab(...)`, `runner_.busy()`.
- Behaviour:
  - `Command::guide` toggles `layout_.guide_shown`, creates `guide_` on first show with `list_lessons(examples dir / "lessons")`, calls `layout()`.
  - `guide_open(i)`: `confirm_discard()`, then `copy_example({lessons_[i].path, …}, documents/NiTCAD/lessons)` and `open_deck`; then `guide_->show_lesson(i)`.
  - `guide_act(step)`: `set` → `set_dimension(deck_->text(), name, value)` then `geometry_edited(*out)`; `edit` → `set_entry` then `geometry_edited`; `sweep` → `deck::sweep_values` then `start_deck_run(&pair)`; `run` → `guide_before_ = extracted_figures(*workspace_)` if any, then `start_deck_run(nullptr)`. Any error → `guide_->set_error(message)`. Then `select_tab` for `show` happens when the run finishes (or at once when there is no run).
  - On run/sweep completion (where `sweep_results_` / `workspace_` are set today): if `guide_pending_`, `guide_->set_report(sweep_report(...)` or `run_report(guide_before_, extracted_figures(*workspace_), step.report))`, select the step's `show` tab, clear `guide_pending_`.
  - `guide_->set_busy(runner_.busy())` wherever `update_command_states()` runs.

- [ ] **Step 1: Write the failing tests**: `guide toggles from View`, `guide opens a lesson as a copy` (copy path under the test's `examples_copy_dir`, original file bytes unchanged), `guide set action edits one line with undo` (Undo restores the deck), `guide sweep action starts a sweep of the named parameter`, `guide action error keeps the step`, `guide actions wait for the run`, `guide lesson switch asks to discard`, `guide report after a sweep` (a 2-value sweep of `examples/lessons/02_body_doping.deck`, Release-only tag `[.guide_run]`, report line starts `Threshold (max gm):`).
- [ ] **Step 2: Run** `ctest --preset release -R "^app_window\.guide_window$"` — expected FAIL.
- [ ] **Step 3: Implement**.
- [ ] **Step 4: Run** — expected PASS; also `ctest --preset release -R "^app_window\.(dimensions_window|sweep_window)$"` PASS (nothing broken nearby).
- [ ] **Step 5: Commit** `app: guide drives edits, sweeps and runs (Unit 26p GU4, GU5)`.

### Task 7: The thirteen lessons and their physics test (GU6, GUC4)

**Files:**
- Create: `examples/lessons/01_gate_and_channel.deck` … `13_mesh_and_models.deck` (names per GU6, two-digit order prefix)
- Create: `tests/app/lessons_test.cpp` (in `nitcad_app_test`); `CMakeLists.txt`: Release-only `add_test(NAME app_lessons COMMAND nitcad_app_test "[.lessons]")` beside `app_examples_3d`
- Test: `tests/app/lessons_test.cpp`

**Interfaces:**
- Consumes: Task 1 (`$Vg`, `$Vd` in `[sweep]`), Task 2 format, the figure names in Global Constraints.
- Every lesson starts from `examples/devices/mosfet_2d.deck`'s parametric form (2D tensor, Si/SiO2, `[parameters]` group Main), adds `Vd` and `Vg` number parameters where a lesson sweeps a bias, and keeps at most two meshes. Lesson 8 uses the vertical-field mobility already exercised by `tests/decks/lombardi_mlda_mosfet_2d.deck`; lesson 12 toggles `[models] field_mobility`.
- Each lesson's prose states the trend and the textbook relation (plain text, Unicode).

- [ ] **Step 1: Write the failing tests** (`[lessons]` fast, `[.lessons]` slow):
  - `[lessons] every lesson parses`: `list_lessons(NITCAD_EXAMPLES_DIR "/lessons").size() == 13`, tracks 4/4/5, orders 1..13.
  - `[.lessons]` one case per lesson that runs its swept step through `nitcad sweep` (`NITCAD_CLI_PATH`) into a temp dir and checks its claim from `figures.csv`: 2 V_th rises with Na (strictly); 3 V_th magnitude and SS fall as tox thins; 4 `Current at the last bias` rises with Vg; 5 SS ≥ 60 mV/dec at every value; 6 V_th rises with work function; 9 V_th falls as Lg shrinks; 10 V_th falls as Vd rises (DIBL > 0); 12 on-current with `field_mobility = on` below off. Lessons 1, 7, 8, 11, 13 check that every run converges and the reported figures exist.
- [ ] **Step 2: Run** `ctest --preset release -R "^app\.lessons$"` — expected FAIL (no lessons).
- [ ] **Step 3: Write the thirteen decks.** Each must time under 60 s per step's sweep on 6 threads; if a claim fails in its own deck, fix the deck's range or drop the claim from the prose — never ship a claim the deck does not show.
- [ ] **Step 4: Run** `ctest --preset release -R "^app\.lessons$"` then `ctest --preset release -R "^app_lessons$"` — expected PASS; record timings.
- [ ] **Step 5: Commit** `examples: thirteen MOSFET guide lessons (Unit 26p GU6)`.

### Task 8: Mutants and records (GUC6)

**Files:**
- Modify: `DECISIONS.md` (Unit 26p: per-slice as-built notes, mutant tallies), `ARCHITECTURE.md` (as-built block: guide pane, lesson format)

- [ ] **Step 1:** Run ≥ 8 mutants per slice (slice 0: Tasks 1–3; slice 1: Tasks 4–6; slice 2: Task 7), sharded over 6 cores with private logs, e.g. `[sweep]` substitution reaching `[solver]`; an action attached to the previous step; Next past the last step; a report reading the wrong figure; a lesson opened in place; `guide_shown=false` moving the inspector; busy not disabling the action; a lesson claim test comparing the wrong column. Each must be caught by the targeted tests above.
- [ ] **Step 2:** Write the as-built notes and tallies; commit `DECISIONS, ARCHITECTURE: Unit 26p as built`.
- [ ] **Step 3:** Report to the owner: branch, commits, targeted results, timings; ask before any push or full suite.
