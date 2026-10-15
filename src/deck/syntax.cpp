#include "NiTCAD/deck/syntax.hpp"

#include <cctype>
#include <fstream>
#include <iterator>
#include <optional>
#include <utility>

namespace NiTCAD::deck {

namespace {

base::Error at(const std::string& source, int line, const std::string& what) {
    return {base::ErrorCode::invalid_input, source + ":" + std::to_string(line) + ": " + what,
            base::ErrorContext{.index = static_cast<std::size_t>(line), .value = std::nullopt}};
}

bool blank(char c) { return c == ' ' || c == '\t' || c == '\r'; }

std::string_view trim(std::string_view s) {
    while (!s.empty() && blank(s.front())) s.remove_prefix(1);
    while (!s.empty() && blank(s.back())) s.remove_suffix(1);
    return s;
}

bool identifier(std::string_view w) {
    if (w.empty() || !std::isalpha(static_cast<unsigned char>(w.front()))) return false;
    for (const char c : w) {
        if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_' && c != '.' && c != '-') {
            return false;
        }
    }
    return true;
}

// Words separated by blanks, joined by one space; nullopt if a word is not an identifier.
std::optional<std::string> words(std::string_view s) {
    std::string out;
    std::size_t i = 0;
    while (i < s.size()) {
        while (i < s.size() && blank(s[i])) ++i;
        std::size_t j = i;
        while (j < s.size() && !blank(s[j])) ++j;
        if (j > i) {
            const std::string_view w = s.substr(i, j - i);
            if (!identifier(w)) return std::nullopt;
            if (!out.empty()) out += ' ';
            out += w;
        }
        i = j;
    }
    return out;
}

// The line without its comment; nullopt for an unterminated quote.
std::optional<std::string_view> strip_comment(std::string_view line) {
    bool quoted = false;
    for (std::size_t i = 0; i < line.size(); ++i) {
        if (line[i] == '"') quoted = !quoted;
        if (line[i] == '#' && !quoted) return line.substr(0, i);
    }
    if (quoted) return std::nullopt;
    return line;
}

}  // namespace

std::expected<Deck, base::Error> parse(std::string_view text, std::string source,
                                       std::filesystem::path directory) {
    Deck deck{std::move(source), std::move(directory), {}};
    if (text.find('\0') != std::string_view::npos) {
        int line = 1;
        for (const char c : text.substr(0, text.find('\0'))) line += c == '\n' ? 1 : 0;
        return std::unexpected(at(deck.source, line, "NUL byte"));
    }
    if (text.starts_with("\xEF\xBB\xBF")) text.remove_prefix(3);  // a UTF-8 byte-order mark

    int number = 0;
    std::size_t pos = 0;
    // A logical line: the physical lines joined by trailing backslashes, without comments.
    const auto next = [&](std::string& out, int& first) -> std::expected<bool, base::Error> {
        out.clear();
        if (pos >= text.size()) return false;
        first = number + 1;
        while (pos < text.size()) {
            const std::size_t end = std::min(text.find('\n', pos), text.size());
            const std::string_view raw = text.substr(pos, end - pos);
            pos = end + 1;
            ++number;
            const auto code = strip_comment(raw);
            if (!code) return std::unexpected(at(deck.source, number, "unterminated quote"));
            const std::string_view t = trim(*code);
            if (!t.empty() && t.back() == '\\') {
                out += t.substr(0, t.size() - 1);
                out += ' ';
                if (pos >= text.size()) {
                    return std::unexpected(at(deck.source, number, "continuation at the end of the file"));
                }
                continue;
            }
            out += t;
            return true;
        }
        return true;
    };

    std::string line;
    int first = 0;
    while (true) {
        auto more = next(line, first);
        if (!more) return std::unexpected(std::move(more.error()));
        if (!*more) break;
        const std::string_view t = trim(line);
        if (t.empty()) continue;
        if (t.front() == '[') {
            if (t.back() != ']') return std::unexpected(at(deck.source, first, "a header must end with ']'"));
            const std::string_view inner = trim(t.substr(1, t.size() - 2));
            std::size_t k = 0;
            while (k < inner.size() && !blank(inner[k])) ++k;
            const std::string_view kind = inner.substr(0, k);
            const auto name = words(inner.substr(k));
            if (!identifier(kind) || !name) {
                return std::unexpected(at(deck.source, first, "malformed header '" + std::string(t) + "'"));
            }
            deck.sections.push_back({std::string(kind), *name, first, {}});
            continue;
        }
        const std::size_t eq = t.find('=');
        if (eq == std::string_view::npos) {
            return std::unexpected(at(deck.source, first, "expected 'key = value' or a [section] header"));
        }
        if (deck.sections.empty()) {
            return std::unexpected(at(deck.source, first, "an entry before the first [section] header"));
        }
        const auto key = words(t.substr(0, eq));
        if (!key || key->empty()) {
            return std::unexpected(at(deck.source, first, "malformed key '" + std::string(trim(t.substr(0, eq))) + "'"));
        }
        const std::string_view value = trim(t.substr(eq + 1));
        if (value.empty()) return std::unexpected(at(deck.source, first, "'" + *key + "' has no value"));
        deck.sections.back().entries.push_back({*key, std::string(value), first});
    }
    return deck;
}

std::expected<Deck, base::Error> read(const std::filesystem::path& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) {
        return std::unexpected(base::Error{base::ErrorCode::invalid_input,
                                           "cannot open the deck '" + path.string() + "'",
                                           std::nullopt});
    }
    const std::string text((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
    return parse(text, path.filename().string(), path.parent_path());
}

std::string unquote(std::string_view value) {
    if (value.size() >= 2 && value.front() == '"' && value.back() == '"') {
        value = value.substr(1, value.size() - 2);
    }
    return std::string(value);
}

}  // namespace NiTCAD::deck
