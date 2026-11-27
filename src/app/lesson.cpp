// The MOSFET guide's lessons (lesson.hpp).
#include "NiTCAD/app/lesson.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <format>
#include <fstream>
#include <iterator>
#include <optional>
#include <system_error>
#include <tuple>

#include "NiTCAD/app/deck_sections.hpp"
#include "NiTCAD/app/dimensions.hpp"
#include "NiTCAD/deck/dimension_sweep.hpp"

namespace NiTCAD::app {

namespace {

namespace fs = std::filesystem;

std::string_view trim(std::string_view s) {
    while (!s.empty() && (s.front() == ' ' || s.front() == '\t')) s.remove_prefix(1);
    while (!s.empty() && (s.back() == ' ' || s.back() == '\t' || s.back() == '\r')) s.remove_suffix(1);
    return s;
}

std::vector<std::string_view> lines_of(std::string_view text) {
    std::vector<std::string_view> out;
    std::size_t i = 0;
    while (i <= text.size()) {
        const auto e = text.find('\n', i);
        if (e == std::string_view::npos) {
            if (i < text.size()) out.push_back(text.substr(i));
            break;
        }
        out.push_back(text.substr(i, e - i));
        i = e + 1;
    }
    return out;
}

// `name = value` split at the first '='; nullopt without one or with either side empty.
std::optional<std::pair<std::string, std::string>> assignment(std::string_view s) {
    const auto eq = s.find('=');
    if (eq == std::string_view::npos) return std::nullopt;
    const auto name = trim(s.substr(0, eq));
    const auto value = trim(s.substr(eq + 1));
    if (name.empty() || value.empty()) return std::nullopt;
    return std::pair{std::string(name), std::string(value)};
}

std::optional<LessonTrack> track_of(std::string_view s) {
    if (s == "basics") return LessonTrack::basics;
    if (s == "intermediate") return LessonTrack::intermediate;
    if (s == "advanced") return LessonTrack::advanced;
    return std::nullopt;
}

using Fail = std::unexpected<LessonError>;

Fail fail(int line, std::string message) { return Fail(LessonError{line, std::move(message)}); }

std::expected<LessonAction, std::string> action_of(std::string_view s, const std::vector<dims::Dimension>& params) {
    LessonAction a;
    if (s == "run") {
        a.kind = LessonAction::Kind::run;
        return a;
    }
    const auto sp = s.find(' ');
    const std::string_view verb = s.substr(0, sp);
    std::string_view rest = sp == std::string_view::npos ? std::string_view() : trim(s.substr(sp + 1));
    if (verb == "edit") {
        if (rest.empty() || rest.front() != '[') return std::unexpected("edit needs a [section]: '" + std::string(s) + "'");
        const auto close = rest.find(']');
        if (close == std::string_view::npos) return std::unexpected("edit needs a [section]: '" + std::string(s) + "'");
        a.section = std::string(trim(rest.substr(1, close - 1)));
        rest = trim(rest.substr(close + 1));
    } else if (verb != "set" && verb != "sweep") {
        return std::unexpected("unknown action '" + std::string(s) + "' (set, sweep, run or edit)");
    }
    const auto kv = assignment(rest);
    if (!kv || (verb == "edit" && a.section.empty())) return std::unexpected("'" + std::string(s) + "' is not <name> = <value>");
    a.name = kv->first;
    a.value = kv->second;
    if (verb == "edit") {
        a.kind = LessonAction::Kind::edit;
        return a;
    }
    const auto p = std::ranges::find(params, a.name, &dims::Dimension::name);
    if (p == params.end()) return std::unexpected("the deck has no parameter '" + a.name + "'");
    if (verb == "set") {
        a.kind = LessonAction::Kind::set;
        return a;
    }
    if (const auto v = deck::sweep_values(a.value, p->length); !v) return std::unexpected("sweep " + a.name + ": " + v.error());
    a.kind = LessonAction::Kind::sweep;
    return a;
}

}  // namespace

std::expected<Lesson, LessonError> parse_lesson(std::string_view deck_text) {
    const auto lines = lines_of(deck_text);
    // The leading comment block only (blank lines may sit in it).
    std::size_t end = 0;
    while (end < lines.size()) {
        const auto t = trim(lines[end]);
        if (!t.empty() && t.front() != '#') break;
        ++end;
    }
    Lesson lesson;
    int lesson_line = 0;
    bool action_seen = false;
    std::vector<dims::Dimension> params;
    bool params_read = false;
    std::string paragraph;  // the current step's open paragraph
    const auto close_paragraph = [&] {
        if (paragraph.empty() || lesson.steps.empty()) return;
        std::string& prose = lesson.steps.back().prose;
        if (!prose.empty()) prose += "\n\n";
        prose += paragraph;
        paragraph.clear();
    };
    for (std::size_t i = 0; i < end; ++i) {
        const int line = static_cast<int>(i) + 1;
        std::string_view t = trim(lines[i]);
        if (t.empty()) continue;
        t = trim(t.substr(1));  // after the '#'
        const auto colon = t.find(':');
        const std::string_view key = colon == std::string_view::npos ? std::string_view() : t.substr(0, colon);
        const std::string_view value = colon == std::string_view::npos ? std::string_view() : trim(t.substr(colon + 1));
        const bool in_step = !lesson.steps.empty();
        if (key == "lesson") {
            lesson.title = std::string(value);
            lesson_line = line;
        } else if (key == "track") {
            const auto tr = track_of(value);
            if (!tr) return fail(line, "track '" + std::string(value) + "' is not basics, intermediate or advanced");
            lesson.track = *tr;
        } else if (key == "order") {
            int n = 0;
            const auto r = std::from_chars(value.data(), value.data() + value.size(), n);
            if (r.ec != std::errc() || r.ptr != value.data() + value.size()) {
                return fail(line, "order '" + std::string(value) + "' is not a whole number");
            }
            lesson.order = n;
        } else if (key == "step") {
            close_paragraph();
            lesson.steps.push_back({.heading = std::string(value), .line = line});
            action_seen = false;
        } else if (key == "do" || key == "show" || key == "report") {
            if (!in_step) return fail(line, "'" + std::string(key) + ":' before the first step");
            close_paragraph();
            LessonStep& step = lesson.steps.back();
            if (key == "do") {
                if (action_seen) return fail(line, "a second action in step '" + step.heading + "'");
                if (!params_read) {
                    auto d = dims::read_dimensions(deck_text);
                    if (!d) return fail(line, d.error());
                    params = std::move(*d);
                    params_read = true;
                }
                auto a = action_of(value, params);
                if (!a) return fail(line, a.error());
                step.action = std::move(*a);
                action_seen = true;
            } else if (key == "show") {
                if (value != "geometry" && value != "fields" && value != "curves" && value != "results") {
                    return fail(line, "show '" + std::string(value) + "' is not geometry, fields, curves or results");
                }
                step.show = std::string(value);
            } else {
                std::string_view rest = value;
                while (!rest.empty()) {
                    const auto comma = rest.find(',');
                    const auto name = trim(rest.substr(0, comma));
                    if (!name.empty()) step.report.emplace_back(name);
                    rest = comma == std::string_view::npos ? std::string_view() : rest.substr(comma + 1);
                }
            }
        } else if (in_step) {
            if (t.empty()) {
                close_paragraph();
            } else {
                if (!paragraph.empty()) paragraph += ' ';
                paragraph += t;
            }
        }
    }
    close_paragraph();
    if (lesson_line == 0) return fail(0, "not a lesson (no '# lesson:' line)");
    if (lesson.steps.empty()) return fail(lesson_line, "lesson '" + lesson.title + "' has no step");
    return lesson;
}

std::vector<Lesson> list_lessons(const fs::path& directory) {
    std::vector<Lesson> out;
    std::error_code ec;
    for (const auto& e : fs::directory_iterator(directory, ec)) {
        if (!e.is_regular_file() || e.path().extension() != ".deck") continue;
        std::ifstream in(e.path(), std::ios::binary);
        const std::string text{std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>()};
        auto l = parse_lesson(text);
        if (!l) continue;
        l->path = e.path();
        out.push_back(std::move(*l));
    }
    std::ranges::sort(out, {}, [](const Lesson& l) { return std::tie(l.track, l.order, l.title); });
    return out;
}

std::string_view to_string(LessonTrack t) noexcept {
    switch (t) {
        case LessonTrack::basics: return "Basics";
        case LessonTrack::intermediate: return "Intermediate";
        case LessonTrack::advanced: return "Advanced";
    }
    return {};
}

std::expected<std::string, std::string> set_entry(std::string_view deck_text, std::string_view section, std::string_view key,
                                                  std::string_view value) {
    const auto s = find_section(deck_text, section);
    if (!s) return std::unexpected("the deck has no [" + std::string(section) + "]");
    const auto e = std::ranges::find(s->entries, key, &SectionLocation::Entry::key);
    if (e == s->entries.end()) return std::unexpected("[" + std::string(section) + "] has no '" + std::string(key) + "'");
    const auto lines = lines_of(deck_text);
    std::string out;
    out.reserve(deck_text.size() + value.size());
    for (std::size_t i = 0; i < lines.size(); ++i) {
        out += static_cast<int>(i) + 1 == e->line ? dims::with_value(lines[i], value) : std::string(lines[i]);
        if (i + 1 < lines.size() || deck_text.ends_with('\n')) out += '\n';
    }
    return out;
}

namespace {

constexpr std::string_view arrow = " \xE2\x86\x92 ", dash = "\xE2\x80\x94";  // " → ", "—"

std::string number(double v) { return std::isfinite(v) ? std::format("{:.3g}", v) : std::string(dash); }

std::string with_unit(std::string line, std::string_view unit) {
    if (!unit.empty()) line += " " + std::string(unit);
    return line;
}

const Figure* find_figure(const std::vector<Figure>& fs, std::string_view name) {
    const auto it = std::ranges::find(fs, name, &Figure::name);
    return it == fs.end() ? nullptr : &*it;
}

void append_line(std::string& out, const std::string& line) {
    if (!out.empty()) out += '\n';
    out += line;
}

}  // namespace

std::string sweep_report(const SweepResults& s, const SweepFigures& f, const std::vector<std::string>& figures) {
    std::string out;
    for (const std::string& name : figures) {
        const auto col = std::ranges::find(f.names, name);
        if (col == f.names.end()) {
            append_line(out, name + ": " + std::string(dash));
            continue;
        }
        const auto j = static_cast<std::size_t>(col - f.names.begin());
        std::string line = name + ": ";
        for (std::size_t i = 0; i < s.runs.size(); ++i) {
            if (i > 0) line += arrow;
            const auto row = std::ranges::find(f.runs, i);
            line += row == f.runs.end() ? std::string(dash) : number(f.rows[static_cast<std::size_t>(row - f.runs.begin())][j]);
        }
        append_line(out, with_unit(std::move(line), f.units[j]));
    }
    return out;
}

std::string run_report(const std::vector<Figure>& before, const std::vector<Figure>& after,
                       const std::vector<std::string>& figures) {
    std::string out;
    for (const std::string& name : figures) {
        const Figure* a = find_figure(after, name);
        if (!a) {
            append_line(out, name + ": " + std::string(dash));
            continue;
        }
        const Figure* b = find_figure(before, name);
        append_line(out, with_unit(name + ": " + (b ? number(b->value) + std::string(arrow) : std::string()) + number(a->value), a->unit));
    }
    return out;
}

}  // namespace NiTCAD::app
