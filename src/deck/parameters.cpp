// Named parameters: [parameters] and `$name` references (Unit 26m MT8, MQ1 (a); the rules in syntax.hpp).

#include <algorithm>
#include <cctype>
#include <charconv>
#include <format>
#include <functional>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "NiTCAD/deck/syntax.hpp"
#include "values.hpp"

namespace NiTCAD::deck {

using detail::length;
using detail::number;
using detail::trimmed;
namespace {

base::Error at(const std::string& source, int line, const std::string& what) {
    return {base::ErrorCode::invalid_input, source + ":" + std::to_string(line) + ": " + what,
            base::ErrorContext{.index = static_cast<std::size_t>(line), .value = std::nullopt}};
}

bool name_start(char c) { return std::isalpha(static_cast<unsigned char>(c)) != 0 || c == '_'; }
bool name_char(char c) { return std::isalnum(static_cast<unsigned char>(c)) != 0 || c == '_'; }
bool blank(char c) { return c == ' ' || c == '\t'; }
// What may stand before a term, and after one (an operator without spaces may follow too).
bool delimiter(char c) {
    return blank(c) || c == '\r' || c == '\n' || c == '(' || c == ')' || c == ',' || c == ':' || c == '=';
}

bool valid_name(std::string_view n) {
    if (n.empty() || !name_start(n.front())) return false;
    for (const char c : n) {
        if (!name_char(c)) return false;
    }
    return true;
}

// A term at i: a reference ($name) or a literal (a number, with nm, um or cm or none).
struct Term {
    std::size_t end = 0;
    std::string_view name;  // empty for a literal
    double value = 0.0;     // a literal's, in cm when it has a unit or none (DK5)
};

std::optional<Term> term(std::string_view s, std::size_t i) {
    if (i >= s.size()) return std::nullopt;
    if (s[i] == '$') {
        std::size_t j = i + 1;
        if (j >= s.size() || !name_start(s[j])) return std::nullopt;
        while (j < s.size() && name_char(s[j])) ++j;
        return Term{j, s.substr(i + 1, j - i - 1), 0.0};
    }
    const std::size_t from = s[i] == '+' ? i + 1 : i;  // from_chars takes '-', not '+'
    std::size_t j = s[i] == '+' || s[i] == '-' ? i + 1 : i;
    if (j >= s.size() || (std::isdigit(static_cast<unsigned char>(s[j])) == 0 && s[j] != '.')) return std::nullopt;
    double v = 0.0;
    const auto [p, ec] = std::from_chars(s.data() + from, s.data() + s.size(), v);
    if (ec != std::errc{}) return std::nullopt;
    j = static_cast<std::size_t>(p - s.data());
    double scale = 1.0;
    const std::string_view unit = s.substr(j, 2);
    if (unit == "nm") {
        scale = 1e-7;
    } else if (unit == "um") {
        scale = 1e-4;
    }
    if (unit == "nm" || unit == "um" || unit == "cm") j += 2;
    if (j < s.size() && !delimiter(s[j]) && s[j] != '+' && s[j] != '-') return std::nullopt;
    return Term{j, {}, v * scale};  // as values.hpp's length: the number times the unit's scale
}

// An operator at j (spaces on both sides or on neither): the position after it and its sign.
std::optional<std::pair<std::size_t, int>> op(std::string_view s, std::size_t j) {
    std::size_t k = j;
    while (k < s.size() && blank(s[k])) ++k;
    const bool before = k > j;
    if (k >= s.size() || (s[k] != '+' && s[k] != '-')) return std::nullopt;
    const int sign = s[k] == '+' ? 1 : -1;
    std::size_t m = ++k;
    while (m < s.size() && blank(s[m])) ++m;
    if (before != (m > k)) return std::nullopt;
    return std::pair{m, sign};
}

const Parameter* find(const std::vector<Parameter>& params, std::string_view name) {
    for (const Parameter& p : params) {
        if (p.name == name) return &p;
    }
    return nullptr;
}

// The names a value refers to, in order.
std::vector<std::string_view> references(std::string_view s) {
    std::vector<std::string_view> out;
    for (std::size_t i = 0; i < s.size(); ++i) {
        if (s[i] != '$') continue;
        if (const auto t = term(s, i); t && !t->name.empty()) out.push_back(t->name);
    }
    return out;
}

}  // namespace

std::expected<std::string, std::string> substitute(std::string_view s, const std::vector<Parameter>& params) {
    std::string out;
    out.reserve(s.size());
    std::size_t i = 0;
    while (i < s.size()) {
        const bool start = i == 0 || delimiter(s[i - 1]);
        const auto first = start ? term(s, i) : std::nullopt;
        if (!first) {
            if (s[i] == '$') return std::unexpected(std::format("'{}': a '$' not starting a reference ($name)", s));
            out += s[i++];
            continue;
        }
        std::vector<std::pair<int, Term>> terms{{1, *first}};
        std::size_t end = first->end;
        while (const auto o = op(s, end)) {
            const auto t = term(s, o->first);
            if (!t) break;
            terms.emplace_back(o->second, *t);
            end = t->end;
        }
        const bool named = std::ranges::any_of(terms, [](const auto& t) { return !t.second.name.empty(); });
        if (!named) {
            out.append(s.substr(i, end - i));  // no reference: as written
            i = end;
            continue;
        }
        for (const auto& [sign, t] : terms) {
            if (!t.name.empty() && find(params, t.name) == nullptr) {
                return std::unexpected(std::format("undefined parameter '${}'", t.name));
            }
        }
        if (terms.size() == 1) {
            out += find(params, terms.front().second.name)->text;  // a lone reference: the text as written
        } else {
            double total = 0.0;
            for (const auto& [sign, t] : terms) {
                double v = t.value;
                if (!t.name.empty()) {
                    const Parameter& p = *find(params, t.name);
                    if (!p.length) {
                        return std::unexpected(std::format("'${}' is a number ({}), not a length: only lengths add", p.name, p.text));
                    }
                    v = p.value;
                }
                total += sign * v;
            }
            out += std::format("{}cm", total);
        }
        i = end;
    }
    return out;
}

std::expected<std::vector<Parameter>, base::Error> parameters(const Deck& deck) {
    const std::string& source = deck.source;
    const Section* section = nullptr;
    for (const Section& s : deck.sections) {
        if (s.kind != "parameters") continue;
        if (section) return std::unexpected(at(source, s.line, "a second [parameters] section"));
        if (!s.name.empty()) return std::unexpected(at(source, s.line, "[parameters] takes no name"));
        section = &s;
    }
    if (!section) return std::vector<Parameter>{};
    return resolve_parameters(section->entries, source);
}

std::expected<std::vector<Parameter>, base::Error> resolve_parameters(const std::vector<Entry>& entries, const std::string& source) {
    std::vector<Parameter> resolved;
    for (std::size_t i = 0; i < entries.size(); ++i) {
        if (!valid_name(entries[i].key)) {
            return std::unexpected(at(source, entries[i].line, std::format("'{}' is not a parameter name (a letter or _, then letters, digits, _)", entries[i].key)));
        }
        for (std::size_t j = 0; j < i; ++j) {
            if (entries[j].key == entries[i].key) {
                return std::unexpected(at(source, entries[i].line, std::format("parameter '{}' given twice (first on line {})", entries[i].key, entries[j].line)));
            }
        }
    }
    // Depth first in the order written, a parameter after the ones it refers to; a reference back into the path
    // being resolved is a cycle.
    std::vector<int> state(entries.size(), 0);  // 0 not seen, 1 on the path, 2 resolved
    std::vector<std::size_t> path;
    const auto index = [&](std::string_view n) -> std::optional<std::size_t> {
        for (std::size_t j = 0; j < entries.size(); ++j) {
            if (entries[j].key == n) return j;
        }
        return std::nullopt;
    };
    std::function<std::expected<void, base::Error>(std::size_t)> visit = [&](std::size_t i) -> std::expected<void, base::Error> {
        if (state[i] == 2) return {};
        const Entry& e = entries[i];
        if (state[i] == 1) {
            std::string cycle;
            const auto from = std::ranges::find(path, i);
            for (auto it = from; it != path.end(); ++it) cycle += "$" + entries[*it].key + " -> ";
            return std::unexpected(at(source, e.line, "a cycle of parameters: " + cycle + "$" + e.key));
        }
        state[i] = 1;
        path.push_back(i);
        for (const std::string_view n : references(e.value)) {
            const auto j = index(n);
            if (!j) return std::unexpected(at(source, e.line, std::format("undefined parameter '${}'", n)));
            if (auto r = visit(*j); !r) return r;
        }
        const auto text = substitute(e.value, resolved);
        if (!text) return std::unexpected(at(source, e.line, text.error()));
        Parameter p{.name = e.key, .text = std::string(trimmed(*text)), .line = e.line};
        // A sum of lengths written out (`10nm + 20nm`, as left by the substitution or as written): its value.
        if (const auto first = term(p.text, 0); first && first->name.empty()) {
            double total = first->value;
            std::size_t end = first->end;
            int terms = 1;
            while (const auto o = op(p.text, end)) {
                const auto t = term(p.text, o->first);
                if (!t || !t->name.empty()) break;
                total += o->second * t->value;
                end = t->end;
                ++terms;
            }
            if (terms > 1 && end == p.text.size()) p.text = std::format("{}cm", total);
        }
        const std::string_view t = p.text;
        if (t.ends_with("nm") || t.ends_with("um") || t.ends_with("cm")) {
            const auto v = length(t);
            if (!v) return std::unexpected(at(source, e.line, std::format("parameter '{}': {}", e.key, v.error())));
            p.value = *v;
            p.length = true;
        } else {
            const auto v = number(t);
            if (!v) return std::unexpected(at(source, e.line, std::format("parameter '{}' = '{}' is neither a length (nm, um, cm) nor a number", e.key, e.value)));
            p.value = *v;
        }
        resolved.push_back(std::move(p));
        path.pop_back();
        state[i] = 2;
        return {};
    };
    for (std::size_t i = 0; i < entries.size(); ++i) {
        if (auto r = visit(i); !r) return std::unexpected(r.error());
    }
    return resolved;
}

std::expected<Deck, base::Error> substituted(const Deck& deck) {
    const std::string& source = deck.source;
    const auto params = parameters(deck);
    if (!params) return std::unexpected(params.error());
    Deck out;
    out.source = deck.source;
    out.directory = deck.directory;
    out.text = deck.text;
    for (const Section& s : deck.sections) {
        if (s.kind == "parameters") continue;
        Section c = s;
        const bool takes = s.kind == "mesh" || s.kind == "region" || s.kind == "contact" || s.kind == "doping" ||
                           s.kind == "sweep";  // Unit 26p (GU1): a bias as a parameter
        for (Entry& e : c.entries) {
            if (e.value.find('$') == std::string::npos) continue;
            if (!takes) {
                return std::unexpected(at(source, e.line, "a parameter reference outside [mesh], [region], [contact], [doping] and [sweep]"));
            }
            auto v = substitute(e.value, *params);
            if (!v) return std::unexpected(at(source, e.line, v.error()));
            e.value = std::move(*v);
        }
        out.sections.push_back(std::move(c));
    }
    return out;
}

std::vector<std::string_view> value_words(std::string_view s) {
    std::vector<std::string_view> out;
    std::size_t i = 0;
    while (i < s.size()) {
        if (blank(s[i]) || s[i] == ',' || s[i] == '\r' || s[i] == '\n') {
            ++i;
            continue;
        }
        std::size_t end = i;
        if (const auto first = term(s, i)) {  // a sum: its terms and operators, one word
            end = first->end;
            while (const auto o = op(s, end)) {
                const auto t = term(s, o->first);
                if (!t) break;
                end = t->end;
            }
        }
        while (end < s.size() && !blank(s[end]) && s[end] != ',' && s[end] != '\r' && s[end] != '\n') ++end;
        out.push_back(s.substr(i, end - i));
        i = end;
    }
    return out;
}

}  // namespace NiTCAD::deck
