#include "NiTCAD/app/dimensions.hpp"

#include "NiTCAD/app/geometry.hpp"
#include "NiTCAD/deck/syntax.hpp"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <format>
#include <optional>

namespace NiTCAD::app::dims {

namespace {

std::string_view trim(std::string_view s) {
    const auto a = s.find_first_not_of(" \t\r");
    if (a == std::string_view::npos) return {};
    const auto b = s.find_last_not_of(" \t\r");
    return s.substr(a, b - a + 1);
}

// The text's lines without their '\n' (a '\r' stays on its line).
std::vector<std::string_view> lines_of(std::string_view text) {
    std::vector<std::string_view> out;
    std::size_t at = 0;
    while (true) {
        const auto end = text.find('\n', at);
        out.push_back(text.substr(at, end == std::string_view::npos ? std::string_view::npos : end - at));
        if (end == std::string_view::npos) return out;
        at = end + 1;
    }
}

// Where a line's comment starts ('#' outside double quotes), or npos.
std::size_t comment_at(std::string_view line) {
    bool quoted = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        if (line[i] == '"') quoted = !quoted;
        else if (line[i] == '#' && !quoted) return i;
    }
    return std::string_view::npos;
}

// The [value_begin, value_end) of an entry line's value: after '=' and its blanks, before the comment's blanks.
std::pair<std::size_t, std::size_t> value_span(std::string_view line) {
    const auto eq = line.find('=');
    if (eq == std::string_view::npos) return {line.size(), line.size()};
    const auto c = std::min(comment_at(line), line.size());
    auto b = eq + 1;
    while (b < c && (line[b] == ' ' || line[b] == '\t')) ++b;
    auto e = c;
    while (e > b && (line[e - 1] == ' ' || line[e - 1] == '\t' || line[e - 1] == '\r')) --e;
    return {b, e};
}

const deck::Section* parameters_section(const deck::Deck& d) {
    const auto it = std::ranges::find(d.sections, std::string_view("parameters"), &deck::Section::kind);
    return it == d.sections.end() ? nullptr : &*it;
}

std::string_view unit_of(std::string_view literal) {
    for (const std::string_view u : {"nm", "um", "cm"}) {
        if (literal.ends_with(u)) return u;
    }
    return {};
}

std::int64_t quanta_per(std::string_view unit) { return unit == "nm" ? 10 : unit == "um" ? 10000 : 100000000; }

// q quanta written in `unit` with no trailing zeros: 220, "nm" -> "22nm"; 11 -> "1.1nm".
std::string format_in(std::int64_t q, std::string_view unit) {
    const std::int64_t per = quanta_per(unit);
    const bool negative = q < 0;
    const std::int64_t a = negative ? -q : q;
    std::string out = (negative ? "-" : "") + std::to_string(a / per);
    if (std::int64_t frac = a % per; frac != 0) {
        std::string digits = std::to_string(per + frac).substr(1);
        while (digits.back() == '0') digits.pop_back();
        out += "." + digits;
    }
    return out + std::string(unit);
}

bool bare_number(std::string_view s) {
    double v = 0.0;
    const auto [end, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
    return ec == std::errc{} && end == s.data() + s.size();
}

std::string display(const deck::Parameter& p) {
    if (!p.length) return std::format("{:g}", p.value);
    const auto q = static_cast<std::int64_t>(std::llround(p.value / geo::quantum_cm));
    const bool nm = std::abs(p.value) < 1e-4;
    std::string s = format_in(q, nm ? "nm" : "um");
    s.resize(s.size() - 2);
    return s + (nm ? " nm" : " µm");
}

// The largest 1-2-5 multiple of one quantum not above a tenth of |q|, at least one quantum.
std::int64_t ladder(std::int64_t q) {
    const double tenth = std::abs(static_cast<double>(q)) / 10.0;
    std::int64_t best = 1;
    for (std::int64_t decade = 1; decade <= 1'000'000'000'000LL; decade *= 10) {
        for (const std::int64_t m : {1, 2, 5}) {
            if (static_cast<double>(m * decade) <= tenth) best = m * decade;
        }
    }
    return best;
}

}  // namespace

std::string with_value(std::string_view line, std::string_view value) {
    const auto [b, e] = value_span(line);
    return std::string(line.substr(0, b)) + std::string(value) + std::string(line.substr(e));
}

std::expected<std::vector<Dimension>, std::string> read_dimensions(std::string_view deck_text) {
    const auto d = deck::parse(deck_text);
    if (!d) return std::unexpected(d.error().message);
    const auto resolved = deck::parameters(*d);
    if (!resolved) return std::unexpected(resolved.error().message);
    std::vector<Dimension> out;
    const deck::Section* s = parameters_section(*d);
    if (!s) return out;
    const auto lines = lines_of(deck_text);
    std::string group;
    int scanned = s->line;  // lines after the header read for "# group:" so far
    for (const deck::Entry& e : s->entries) {
        for (; scanned < e.line - 1; ++scanned) {
            const std::string_view l = trim(lines[static_cast<std::size_t>(scanned)]);
            if (l.starts_with('#')) {
                const std::string_view c = trim(l.substr(1));
                if (c.starts_with("group:")) group = std::string(trim(c.substr(6)));
            }
        }
        scanned = e.line;
        const auto p = std::ranges::find(*resolved, e.key, &deck::Parameter::name);
        if (p == resolved->end()) continue;
        Dimension r;
        r.name = e.key;
        const std::string_view raw = lines[static_cast<std::size_t>(e.line - 1)];
        const auto c = comment_at(raw);
        r.label = c == std::string_view::npos ? std::string() : std::string(trim(raw.substr(c + 1)));
        if (r.label.empty()) r.label = r.name;
        r.group = group;
        r.expr = std::string(trim(e.value));
        r.length = p->length;
        r.value = p->value;
        r.value_text = display(*p);
        r.line = e.line;
        out.push_back(std::move(r));
    }
    return out;
}

std::expected<std::string, std::string> set_dimension(std::string_view deck_text, std::string_view name, std::string_view typed) {
    const std::string who(name);
    const auto rows = read_dimensions(deck_text);
    if (!rows) return std::unexpected(who + ": " + rows.error());
    const auto row = std::ranges::find(*rows, name, &Dimension::name);
    if (row == rows->end()) return std::unexpected(who + ": no such dimension");
    std::string value(trim(typed));
    if (value.empty()) return std::unexpected(who + ": a value is needed");
    if (row->length && bare_number(value)) {
        const std::string_view u = unit_of(row->expr);
        value += u.empty() || row->expr.find('$') != std::string::npos ? (std::abs(row->value) < 1e-4 ? "nm" : "um") : std::string(u);
    }
    auto lines = lines_of(deck_text);
    const auto at = static_cast<std::size_t>(row->line - 1);
    const std::string edited = with_value(lines[at], value);
    if (edited == lines[at]) return std::string(deck_text);
    std::string out;
    out.reserve(deck_text.size() + value.size());
    for (std::size_t k = 0; k < lines.size(); ++k) {
        if (k) out += '\n';
        out += k == at ? std::string_view(edited) : lines[k];
    }
    const auto d = deck::parse(out);
    if (!d) return std::unexpected(who + ": " + d.error().message);
    if (const auto s = deck::substituted(*d); !s) return std::unexpected(who + ": " + s.error().message);
    return out;
}

std::string step_dimension(const Dimension& d, int notches) {
    if (d.expr.find('$') != std::string::npos || notches == 0) return d.expr;
    if (!d.length) {
        if (!bare_number(d.expr)) return d.expr;
        char buf[32];
        const auto r = std::to_chars(buf, buf + sizeof buf, d.value * std::pow(10.0, notches / 10.0), std::chars_format::scientific);
        return std::string(buf, r.ptr);  // shortest round trip, as a doping level is written
    }
    const std::string_view unit = unit_of(d.expr);
    if (unit.empty()) return d.expr;
    auto q = static_cast<std::int64_t>(std::llround(d.value / geo::quantum_cm));
    const bool positive = q > 0;
    for (int k = 0; k < std::abs(notches); ++k) q += (notches > 0 ? 1 : -1) * ladder(q);
    if (positive) q = std::max<std::int64_t>(q, 1);
    return format_in(q, unit);
}

namespace {

// A substituted deck length in cm: a number with nm, um or cm; unitless is cm (DK5).
std::optional<double> length_cm(std::string_view w) {
    double v = 0.0;
    const auto [end, ec] = std::from_chars(w.data(), w.data() + w.size(), v);
    if (ec != std::errc{}) return std::nullopt;
    const std::string_view unit = w.substr(static_cast<std::size_t>(end - w.data()));
    if (unit.empty() || unit == "cm") return v;
    if (unit == "um") return v * 1e-4;
    if (unit == "nm") return v * 1e-7;
    return std::nullopt;
}

}  // namespace

std::vector<double> lengths(std::string_view value) {
    std::vector<double> out;
    for (const std::string_view w : deck::value_words(value)) {
        const auto v = length_cm(w);
        if (!v) return {};
        out.push_back(*v);
    }
    return out;
}

namespace {

std::vector<double> numbers(std::string_view value) { return lengths(value); }

struct Candidate {
    Span span;
    int rank = 0;  // 0 a box, 1 a doping window, 2 a mesh piece
};

void add_box(std::vector<Candidate>& out, std::string_view value, const std::string& source) {
    const auto v = numbers(value);
    if (v.size() != 4 && v.size() != 6) return;
    const std::size_t n = v.size() / 2;
    for (std::size_t a = 0; a < n; ++a) {
        Candidate c;
        c.span = {static_cast<int>(a), v[2 * a], v[2 * a + 1], {v[0], v[2], n == 3 ? v[4] : 0.0}, source};
        c.span.at[a] = v[2 * a];
        out.push_back(std::move(c));
    }
}

// Every extent a length could show on, in a fixed order (the same for the deck at p and at p + δ).
std::vector<Candidate> candidates(const deck::Deck& d) {
    std::vector<Candidate> out;
    int profile = 0;
    for (const deck::Section& s : d.sections) {
        for (const deck::Entry& e : s.entries) {
            if ((s.kind == "region" || s.kind == "contact") && e.key == "box") {
                add_box(out, e.value, s.kind + " " + s.name);
            } else if (s.kind == "mesh" && e.key == "domain") {
                add_box(out, e.value, "mesh domain");
            } else if (s.kind == "mesh" && (e.key == "x" || e.key == "y" || e.key == "z")) {
                const int axis = e.key[0] - 'x';
                for (std::size_t at = e.value.find("linspace("); at != std::string::npos; at = e.value.find("linspace(", at + 1)) {
                    const auto close = e.value.find(')', at);
                    const std::string_view args = std::string_view(e.value).substr(at + 9, close - at - 9);
                    const auto comma1 = args.find(','), comma2 = args.find(',', comma1 + 1);
                    const auto a = length_cm(trim(args.substr(0, comma1))), b = length_cm(trim(args.substr(comma1 + 1, comma2 - comma1 - 1)));
                    if (!a || !b) continue;
                    Candidate c;
                    c.rank = 2;
                    c.span = {axis, *a, *b, {}, "mesh " + e.key};
                    c.span.at[static_cast<std::size_t>(axis)] = *a;
                    out.push_back(std::move(c));
                }
            } else if (s.kind == "doping" && e.key == "profile") {
                ++profile;
                for (const std::string_view w : deck::value_words(e.value)) {
                    if (w.size() < 3 || w[1] != '=' || (w[0] != 'x' && w[0] != 'y' && w[0] != 'z')) continue;
                    const auto colon = w.find(':');
                    if (colon == std::string_view::npos) continue;
                    const auto a = length_cm(w.substr(2, colon - 2)), b = length_cm(w.substr(colon + 1));
                    if (!a || !b) continue;
                    Candidate c;
                    c.rank = 1;
                    const int axis = w[0] - 'x';
                    c.span = {axis, *a, *b, {}, "doping " + std::to_string(profile)};
                    c.span.at[static_cast<std::size_t>(axis)] = *a;
                    out.push_back(std::move(c));
                }
            }
        }
    }
    return out;
}

}  // namespace

std::vector<Span> dimension_spans(std::string_view deck_text, std::string_view name) {
    const auto rows = read_dimensions(deck_text);
    if (!rows) return {};
    const auto row = std::ranges::find(*rows, name, &Dimension::name);
    if (row == rows->end() || !row->length) return {};
    const double delta = std::max(0.01 * std::abs(row->value), geo::quantum_cm);
    const auto moved = set_dimension(deck_text, name, std::format("{}cm", row->value + delta));
    if (!moved) return {};
    const auto d0 = deck::parse(deck_text), d1 = deck::parse(*moved);
    if (!d0 || !d1) return {};
    const auto s0 = deck::substituted(*d0), s1 = deck::substituted(*d1);
    if (!s0 || !s1) return {};
    auto c0 = candidates(*s0);
    const auto c1 = candidates(*s1);
    if (c0.size() != c1.size()) return {};
    const auto same = [&](double a, double b) { return std::abs(a - b) <= 1e-6 * delta; };
    std::array<const Candidate*, 3> best{};
    for (std::size_t k = 0; k < c0.size(); ++k) {
        const Span& a = c0[k].span;
        const Span& b = c1[k].span;
        const bool lo_moved = !same(a.lo, b.lo), hi_moved = !same(a.hi, b.hi);
        if (lo_moved == hi_moved) continue;
        if (!same(std::abs((b.hi - b.lo) - (a.hi - a.lo)), delta)) continue;
        c0[k].span.grows = (hi_moved ? b.hi - a.hi : b.lo - a.lo) > 0 ? 1 : -1;
        const Candidate*& slot = best[static_cast<std::size_t>(a.axis)];
        const double len = std::abs(a.hi - a.lo);
        if (!slot) {
            slot = &c0[k];
            continue;
        }
        const double held = std::abs(slot->span.hi - slot->span.lo);
        if (len < held && !same(len, held)) slot = &c0[k];
        else if (same(len, held) && c0[k].rank < slot->rank) slot = &c0[k];
    }
    std::vector<Span> out;
    for (const Candidate* c : best) {
        if (c) out.push_back(c->span);
    }
    return out;
}

}  // namespace NiTCAD::app::dims

namespace NiTCAD::app::dims {

std::vector<Mark> dimension_marks(std::string_view deck_text) {
    std::vector<Mark> out;
    const auto rows = read_dimensions(deck_text);
    if (!rows) return out;
    for (const Dimension& d : *rows) {
        if (!d.length) continue;
        for (const Span& s : dimension_spans(deck_text, d.name)) {
            Mark m;
            m.name = d.name;
            m.label = d.label + " = " + d.value_text;
            m.axis = s.axis;
            m.a = s.at;
            m.b = s.at;
            m.b[static_cast<std::size_t>(s.axis)] = s.hi;
            m.grows = s.grows;
            out.push_back(std::move(m));
        }
    }
    return out;
}

std::string length_literal(const Dimension& d, double cm) {
    std::string_view unit = d.expr.find('$') == std::string::npos ? unit_of(d.expr) : std::string_view();
    if (unit.empty()) unit = std::abs(cm) < 1e-4 ? "nm" : "um";
    auto q = static_cast<std::int64_t>(std::llround(cm / geo::quantum_cm));
    if (d.value > 0.0) q = std::max<std::int64_t>(q, 1);
    return format_in(q, unit);
}

std::expected<std::string, std::string> substituted_text(std::string_view deck_text) {
    const auto d = deck::parse(deck_text);
    if (!d) return std::unexpected(d.error().message);
    const auto params = deck::parameters(*d);
    if (!params) return std::unexpected(params.error().message);
    auto lines = lines_of(deck_text);
    std::vector<std::string> edited(lines.size());
    std::vector<bool> changed(lines.size(), false);
    for (const deck::Section& s : d->sections) {
        if (s.kind != "mesh" && s.kind != "region" && s.kind != "contact" && s.kind != "doping") continue;
        for (const deck::Entry& e : s.entries) {
            if (e.value.find('$') == std::string::npos) continue;
            const auto at = static_cast<std::size_t>(e.line - 1);
            if (trim(lines[at]).ends_with('\\')) return std::unexpected("line " + std::to_string(e.line) + ": a continued line with references");
            const auto v = deck::substitute(e.value, *params);
            if (!v) return std::unexpected("line " + std::to_string(e.line) + ": " + v.error());
            edited[at] = with_value(lines[at], *v);
            changed[at] = true;
        }
    }
    std::string out;
    out.reserve(deck_text.size() + 64);
    for (std::size_t k = 0; k < lines.size(); ++k) {
        if (k) out += '\n';
        out += changed[k] ? std::string_view(edited[k]) : lines[k];
    }
    return out;
}

}  // namespace NiTCAD::app::dims
