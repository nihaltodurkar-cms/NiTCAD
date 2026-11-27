// The MOSFET guide's lessons (Unit 26p, GU2, GU5, DECISIONS.md): an example deck whose leading comment block carries
// the lesson, read here for the guide pane (src/app/guide_pane.hpp).
//
//   # lesson: Body doping                  the title
//   # track: basics                        basics | intermediate | advanced
//   # order: 2                             the lesson's place in its track
//   # step: What the doping does           a step's heading, then its prose: '# ' lines, a bare '#' between
//   # The body's acceptors ...             paragraphs
//   # do: sweep Na = 1e16, 1e17, 1e18      at most one action: set <p> = <v> | sweep <p> = <values> | run |
//   # show: results                        edit [<section>] <key> = <value>; then the tab to show
//   # report: Threshold (max gm)           and the figures to report (Unit 24's names)
//
// Only the comment lines before the deck's first other line are read; the deck itself is an ordinary deck.
#pragma once

#include <expected>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/app/sweep_results.hpp"
#include "NiTCAD/app/workspace.hpp"

namespace NiTCAD::app {

enum class LessonTrack : unsigned char { basics, intermediate, advanced };

struct LessonAction {
    enum class Kind : unsigned char { none, set, sweep, run, edit } kind = Kind::none;
    std::string section;  // edit only, without brackets
    std::string name;     // set, sweep: the parameter; edit: the key
    std::string value;    // set: the value; sweep: the values as written (DT6's syntax); edit: the value
};

struct LessonStep {
    std::string heading;
    std::string prose;                // paragraphs joined by "\n\n", a paragraph's lines by ' '
    LessonAction action;
    std::string show;                 // "", "geometry", "fields", "curves" or "results"
    std::vector<std::string> report;  // figure names
    int line = 0;                     // the `# step:` line, 1-based
};

struct Lesson {
    std::string title;
    LessonTrack track = LessonTrack::basics;
    int order = 0;
    std::vector<LessonStep> steps;
    std::filesystem::path path;  // set by list_lessons
};

struct LessonError {
    int line = 0;  // 1-based; 0: not a lesson (no `# lesson:` line)
    std::string message;
};

// The lesson a deck carries. A comment line that is not one of the keys above is prose (in a step) or ignored (before
// the first step). Errors (on their line): an unknown action or tab, an action or `show`/`report` before the
// first step, a second action in a step, a track that is not one of the three, an order that is not a whole number,
// no step, a set or sweep of a name the deck has no parameter for, sweep values not of the parameter's kind.
[[nodiscard]] std::expected<Lesson, LessonError> parse_lesson(std::string_view deck_text);

// The *.deck files of a directory that are lessons, by track, then order, then title. A missing directory has none.
[[nodiscard]] std::vector<Lesson> list_lessons(const std::filesystem::path& directory);

[[nodiscard]] std::string_view to_string(LessonTrack t) noexcept;  // "Basics", "Intermediate", "Advanced"

// The deck with `key`'s value in the first [section] replaced (one line, its spacing and comment kept). Errors: no
// such section, or no such key in it.
[[nodiscard]] std::expected<std::string, std::string> set_entry(std::string_view deck_text, std::string_view section,
                                                                std::string_view key, std::string_view value);

// GU5: one line per figure, "<name>: v1 → v2 → v3 <unit>" in the sweep's order, "—" for a run without it.
[[nodiscard]] std::string sweep_report(const SweepResults& s, const SweepFigures& f, const std::vector<std::string>& figures);

// GU5: one line per figure, "<name>: before → after <unit>"; a figure only after: "<name>: after <unit>"; neither:
// "<name>: —".
[[nodiscard]] std::string run_report(const std::vector<Figure>& before, const std::vector<Figure>& after,
                                     const std::vector<std::string>& figures);

}  // namespace NiTCAD::app
