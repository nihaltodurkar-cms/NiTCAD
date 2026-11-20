// A deck's text rewritten by minimal edits (internal to the app; Unit 26k GD3, shared with Unit 26l SD3): sections'
// lines replaced, dropped or inserted, every other line, comment, blank and line ending kept byte for byte.
#pragma once

#include <algorithm>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "NiTCAD/deck/syntax.hpp"

namespace NiTCAD::app::deck_text {

inline bool blank(char c) { return c == ' ' || c == '\t' || c == '\r'; }

inline std::string_view trim(std::string_view s) {
    while (!s.empty() && blank(s.front())) s.remove_prefix(1);
    while (!s.empty() && blank(s.back())) s.remove_suffix(1);
    return s;
}

// Tokens separated by blanks and commas.
inline std::vector<std::string_view> words(std::string_view s) {
    std::vector<std::string_view> out;
    std::size_t i = 0;
    while (i < s.size()) {
        while (i < s.size() && (blank(s[i]) || s[i] == ',')) ++i;
        std::size_t j = i;
        while (j < s.size() && !blank(s[j]) && s[j] != ',') ++j;
        if (j > i) out.push_back(s.substr(i, j - i));
        i = j;
    }
    return out;
}

inline const deck::Entry* find_entry(const deck::Section& s, std::string_view key) {
    for (const auto& e : s.entries) {
        if (e.key == key) return &e;
    }
    return nullptr;
}

inline std::optional<bool> on_off(std::string_view v) {
    v = trim(v);
    if (v == "on" || v == "yes" || v == "true") return true;
    if (v == "off" || v == "no" || v == "false") return false;
    return std::nullopt;
}

struct Lines {
    std::vector<std::string> body, eol;
    std::string nl = "\n";
};

inline Lines split(std::string_view t) {
    Lines l;
    if (t.find("\r\n") != std::string_view::npos) l.nl = "\r\n";
    std::size_t i = 0;
    while (i < t.size()) {
        const std::size_t j = t.find('\n', i);
        if (j == std::string_view::npos) {
            l.body.emplace_back(t.substr(i));
            l.eol.emplace_back();
            break;
        }
        const bool cr = j > i && t[j - 1] == '\r';
        l.body.emplace_back(t.substr(i, j - i - (cr ? 1 : 0)));
        l.eol.emplace_back(cr ? "\r\n" : "\n");
        i = j + 1;
    }
    return l;
}

inline bool continues(std::string_view s) {
    const std::size_t hash = s.find('#');
    if (hash != std::string_view::npos) s = s.substr(0, hash);
    s = trim(s);
    return !s.empty() && s.back() == '\\';
}

// 0-based line spans.
inline int entry_last(const Lines& l, int first) {
    int k = first;
    while (k + 1 < static_cast<int>(l.body.size()) && continues(l.body[static_cast<std::size_t>(k)])) ++k;
    return k;
}

inline int section_last(const Lines& l, const deck::Section& s) {
    int last = s.line - 1;
    for (const auto& e : s.entries) last = std::max(last, entry_last(l, e.line - 1));
    return last;
}

// A section's lines with entries replaced (`replace`: the new text or nullopt to drop, per entry index), new lines
// inserted after an entry (`after`, -1: after the header), and the header replaced.
struct SectionEdit {
    std::map<std::size_t, std::optional<std::string>> replace;
    std::multimap<int, std::string> after;
    std::optional<std::string> header;
};

class Rewrite {
public:
    explicit Rewrite(std::string_view text) : text_(text), l_(split(text)) {}

    [[nodiscard]] const Lines& lines() const noexcept { return l_; }
    [[nodiscard]] int size() const noexcept { return static_cast<int>(l_.body.size()); }

    [[nodiscard]] std::vector<std::string> block(const deck::Section& s) const {
        const int first = s.line - 1, last = section_last(l_, s);
        return std::vector<std::string>(l_.body.begin() + first, l_.body.begin() + last + 1);
    }

    [[nodiscard]] std::vector<std::string> rewrite(const deck::Section& s, const SectionEdit& e) const {
        std::vector<std::string> out;
        out.push_back(e.header ? *e.header : l_.body[static_cast<std::size_t>(s.line - 1)]);
        const auto put_after = [&](int k) {
            for (auto [it, end] = e.after.equal_range(k); it != end; ++it) out.push_back(it->second);
        };
        put_after(-1);
        int line = s.line;  // 0-based index of the next line to copy
        for (std::size_t k = 0; k < s.entries.size(); ++k) {
            const int first = s.entries[k].line - 1, last = entry_last(l_, first);
            for (; line < first; ++line) out.push_back(l_.body[static_cast<std::size_t>(line)]);  // blanks, comments
            if (const auto it = e.replace.find(k); it != e.replace.end()) {
                if (it->second) out.push_back(*it->second);
            } else {
                for (int j = first; j <= last; ++j) out.push_back(l_.body[static_cast<std::size_t>(j)]);
            }
            line = last + 1;
            put_after(static_cast<int>(k));
        }
        return out;
    }

    // The section's block replaced by `lines` (nothing recorded when they are its own).
    void change(const deck::Section& s, std::vector<std::string> lines) {
        const int first = s.line - 1, last = section_last(l_, s);
        if (lines != block(s)) changes_.push_back({first, last - first + 1, std::move(lines)});
    }
    // A section removed with the blank line before it.
    void drop(const deck::Section& s) {
        int first = s.line - 1;
        const int last = section_last(l_, s);
        if (first > 0 && trim(l_.body[static_cast<std::size_t>(first - 1)]).empty()) --first;
        changes_.push_back({first, last - first + 1, {}});
    }
    // Lines inserted before 0-based line `at` (the line count: at the end).
    void insert(int at, std::vector<std::string> lines) { changes_.push_back({at, 0, std::move(lines)}); }

    // The text with the changes applied in order; an insertion keeps its order among those at one place.
    [[nodiscard]] std::string apply() {
        std::stable_sort(changes_.begin(), changes_.end(), [](const Change& a, const Change& b) { return a.first < b.first; });
        std::string out;
        out.reserve(text_.size() + 256);
        int line = 0;
        const int n = size();
        const auto copy = [&](int k) { out += l_.body[static_cast<std::size_t>(k)] + l_.eol[static_cast<std::size_t>(k)]; };
        for (const Change& c : changes_) {
            for (; line < c.first && line < n; ++line) copy(line);
            if (c.first >= n && !out.empty() && !out.ends_with('\n')) out += l_.nl;  // after a last line without one
            for (const std::string& s : c.lines) out += s + l_.nl;
            line = std::max(line, c.first + c.count);
        }
        for (; line < n; ++line) copy(line);
        return out;
    }

private:
    struct Change {
        int first = 0, count = 0;  // the lines replaced (count 0: an insertion before `first`)
        std::vector<std::string> lines;
    };
    std::string_view text_;
    Lines l_;
    std::vector<Change> changes_;
};

}  // namespace NiTCAD::app::deck_text
