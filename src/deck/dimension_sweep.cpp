// Dimension sweeps (Unit 26n, DT6): deck/dimension_sweep.hpp.
#include "NiTCAD/deck/dimension_sweep.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <optional>

#include "values.hpp"

namespace NiTCAD::deck {

namespace {

using detail::trimmed;

base::Error invalid(std::string what) { return {base::ErrorCode::invalid_input, std::move(what), std::nullopt}; }

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

// Where an entry line's value starts and ends: after the first '=' and its blanks, up to the comment ('#' outside
// double quotes) or the line's end, trailing blanks (and a CR) left out.
std::pair<std::size_t, std::size_t> value_span(std::string_view line) {
    const std::size_t eq = line.find('=');
    std::size_t b = eq == std::string_view::npos ? line.size() : eq + 1;
    while (b < line.size() && (line[b] == ' ' || line[b] == '\t')) ++b;
    std::size_t e = line.size();
    bool quoted = false;
    for (std::size_t k = b; k < line.size(); ++k) {
        if (line[k] == '"') quoted = !quoted;
        if (line[k] == '#' && !quoted) {
            e = k;
            break;
        }
    }
    while (e > b && detail::is_blank(line[e - 1])) --e;
    return {b, e};
}

// A literal of the kind, as a sweep takes it: its value (cm for a length) and, for a length, its unit's scale.
struct Literal {
    double value = 0.0;
    double scale = 1.0;
    std::string_view unit;
};

std::expected<Literal, std::string> literal(std::string_view s, bool length) {
    s = trimmed(s);
    if (s.find('$') != std::string_view::npos) {
        return std::unexpected(std::format("'{}': a sweep takes literal values, not references", s));
    }
    const bool unit = s.ends_with("nm") || s.ends_with("um") || s.ends_with("cm");
    if (length && !unit) return std::unexpected(std::format("'{}' is not a length: give its unit (nm, um or cm)", s));
    if (!length && unit) return std::unexpected(std::format("'{}' is a length, not a number", s));
    if (!length) {
        const auto v = detail::number(s);
        if (!v) return std::unexpected(v.error());
        return Literal{.value = *v, .scale = 1.0, .unit = {}};
    }
    const auto v = detail::length(s);
    if (!v) return std::unexpected(v.error());
    const std::string_view u = s.substr(s.size() - 2);
    return Literal{.value = *v, .scale = u == "nm" ? 1e-7 : u == "um" ? 1e-4 : 1.0, .unit = u};
}

constexpr std::size_t max_values = 1000;

}  // namespace

std::expected<Override, std::string> parse_override(std::string_view text) {
    const std::size_t eq = text.find('=');
    if (eq == std::string_view::npos) return std::unexpected(std::format("'{}' is not name=value", text));
    Override o{.name = std::string(trimmed(text.substr(0, eq))), .value = std::string(trimmed(text.substr(eq + 1)))};
    if (o.name.empty()) return std::unexpected(std::format("'{}' has no parameter name", text));
    if (o.value.empty()) return std::unexpected(std::format("'{}' has no value", text));
    return o;
}

std::expected<Deck, base::Error> with_overrides(const Deck& deck, const std::vector<Override>& overrides) {
    if (overrides.empty()) return deck;
    const auto before = parameters(deck);
    if (!before) return std::unexpected(before.error());
    auto lines = lines_of(deck.text);
    std::vector<std::string> edited(lines.begin(), lines.end());
    for (const Override& o : overrides) {
        const auto p = std::ranges::find(*before, o.name, &Parameter::name);
        if (p == before->end()) {
            return std::unexpected(invalid(std::format("{}: no parameter '{}' in [parameters]", deck.source, o.name)));
        }
        const auto at = static_cast<std::size_t>(p->line - 1);
        if (at >= lines.size()) return std::unexpected(invalid(std::format("{}: parameter '{}' has no line", deck.source, o.name)));
        const std::string_view line = lines[at];
        const auto [b, e] = value_span(line);
        if (line.substr(b, e - b).ends_with('\\')) {
            return std::unexpected(invalid(std::format("{}:{}: parameter '{}' is on a continued line: it cannot be set",
                                                       deck.source, p->line, o.name)));
        }
        edited[at] = std::string(line.substr(0, b)) + o.value + std::string(line.substr(e));
    }
    std::string text;
    for (std::size_t k = 0; k < edited.size(); ++k) {
        if (k) text += '\n';
        text += edited[k];
    }
    auto out = parse(text, deck.source, deck.directory);
    if (!out) return std::unexpected(out.error());
    const auto after = parameters(*out);
    if (!after) return std::unexpected(after.error());
    for (const Override& o : overrides) {
        const auto was = std::ranges::find(*before, o.name, &Parameter::name);
        const auto now = std::ranges::find(*after, o.name, &Parameter::name);
        if (was->length != now->length) {
            return std::unexpected(invalid(std::format(
                "{}:{}: parameter '{}' is a {}: '{}' is not one", deck.source, was->line, o.name,
                was->length ? "length (nm, um or cm)" : "number (no unit)", o.value)));
        }
    }
    return out;
}

std::expected<std::vector<std::string>, std::string> sweep_values(std::string_view list, bool length) {
    list = trimmed(list);
    if (list.empty()) return std::unexpected("a sweep needs at least one value");
    std::vector<std::string> out;
    std::vector<double> values;
    if (list.find(':') != std::string_view::npos) {
        const auto parts = detail::items(list, ':');
        if (parts.size() != 3) return std::unexpected("a range is 'first : last : step'");
        const auto a = literal(parts[0], length), b = literal(parts[1], length), h = literal(parts[2], length);
        if (!a) return std::unexpected(a.error());
        if (!b) return std::unexpected(b.error());
        if (!h) return std::unexpected(h.error());
        if (h->value == 0.0 || (b->value - a->value) / h->value < 0.0) {
            return std::unexpected("the step of a range must be non-zero and point from first to last");
        }
        const double steps = std::floor((b->value - a->value) / h->value + 1e-9);
        if (steps + 1 > static_cast<double>(max_values)) return std::unexpected(std::format("more than {} values", max_values));
        for (long long k = 0; k <= static_cast<long long>(steps); ++k) {
            double v = a->value + static_cast<double>(k) * h->value;
            if (std::abs(v - b->value) <= 1e-9 * std::abs(h->value)) v = b->value;
            // In first's unit, to twelve digits: 0.1um + 0.1um is written 0.2um.
            out.push_back(std::format("{:.12g}{}", v / a->scale, a->unit));
            values.push_back(v);
        }
    } else {
        for (const std::string_view item : detail::items(list)) {
            const auto v = literal(item, length);
            if (!v) return std::unexpected(v.error());
            if (out.size() == max_values) return std::unexpected(std::format("more than {} values", max_values));
            out.emplace_back(trimmed(item));
            values.push_back(v->value);
        }
    }
    for (std::size_t i = 0; i < values.size(); ++i) {
        for (std::size_t j = 0; j < i; ++j) {
            if (values[i] == values[j]) return std::unexpected(std::format("'{}' is given twice", out[i]));
        }
    }
    return out;
}

std::filesystem::path sweep_directory(const std::filesystem::path& deck_path, std::string_view name) {
    return deck_path.parent_path() / (deck_path.stem().string() + "_sweep_" + std::string(name));
}

std::filesystem::path sweep_point_directory(const std::filesystem::path& sweep_dir, std::string_view name,
                                            std::string_view value) {
    return sweep_dir / (std::string(name) + "=" + std::string(value));
}

Outputs with_outputs_in(Outputs o, const std::filesystem::path& deck_dir, const std::filesystem::path& to) {
    const std::filesystem::path base = deck_dir.lexically_normal();
    const auto moved = [&](std::filesystem::path& p) {
        if (p.empty()) return;
        const std::filesystem::path rel = p.lexically_normal().lexically_relative(base);
        const bool below = !rel.empty() && *rel.begin() != "..";
        p = to / (below ? rel : p.filename());
    };
    for (std::filesystem::path* p : {&o.sweep, &o.device, &o.mesh, &o.checkpoint, &o.table, &o.transient,
                                     &o.small_signal, &o.adapt, &o.adapt_meshes, &o.adapt_history}) {
        moved(*p);
    }
    return o;
}

}  // namespace NiTCAD::deck
