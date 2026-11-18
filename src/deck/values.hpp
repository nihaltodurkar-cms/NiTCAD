// Values of deck entries (private; deck/build.hpp, DECISIONS.md DK5-DK7). Each reader returns the
// value or a message without the location; the caller adds "<source>:<line>: ".
#pragma once

#include <cctype>
#include <charconv>
#include <cmath>
#include <expected>
#include <span>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

namespace NiTCAD::deck::detail {

inline bool is_blank(char c) { return c == ' ' || c == '\t' || c == '\r'; }

inline std::string_view trimmed(std::string_view s) {
    while (!s.empty() && is_blank(s.front())) s.remove_prefix(1);
    while (!s.empty() && is_blank(s.back())) s.remove_suffix(1);
    return s;
}

// Blank-separated tokens.
inline std::vector<std::string_view> tokens(std::string_view s) {
    std::vector<std::string_view> out;
    std::size_t i = 0;
    while (i < s.size()) {
        while (i < s.size() && is_blank(s[i])) ++i;
        std::size_t j = i;
        while (j < s.size() && !is_blank(s[j])) ++j;
        if (j > i) out.push_back(s.substr(i, j - i));
        i = j;
    }
    return out;
}

// Items separated by `sep`, each trimmed (empty items kept, so "1,,2" is caught by the reader).
inline std::vector<std::string_view> items(std::string_view s, char sep = ',') {
    std::vector<std::string_view> out;
    std::size_t start = 0;
    while (true) {
        const std::size_t k = s.find(sep, start);
        out.push_back(trimmed(s.substr(start, k == std::string_view::npos ? std::string_view::npos : k - start)));
        if (k == std::string_view::npos) break;
        start = k + 1;
    }
    return out;
}

// A finite real number in the C locale's form (1, -2.5, 3e-4); nothing else in the text.
inline std::expected<double, std::string> number(std::string_view s) {
    s = trimmed(s);
    double v = 0.0;
    const char* first = s.data();
    if (!s.empty() && s.front() == '+') ++first;
    const auto [end, ec] = std::from_chars(first, s.data() + s.size(), v);
    if (s.empty() || ec != std::errc{} || end != s.data() + s.size()) {
        return std::unexpected("'" + std::string(s) + "' is not a number");
    }
    if (!std::isfinite(v)) return std::unexpected("'" + std::string(s) + "' is not finite");
    return v;
}

// A length in cm: a number, optionally with the suffix cm, um or nm (DECISIONS.md DK5).
inline std::expected<double, std::string> length(std::string_view s) {
    s = trimmed(s);
    double scale = 1.0;
    if (s.ends_with("um")) {
        scale = 1e-4;
        s.remove_suffix(2);
    } else if (s.ends_with("nm")) {
        scale = 1e-7;
        s.remove_suffix(2);
    } else if (s.ends_with("cm")) {
        s.remove_suffix(2);
    }
    auto v = number(s);
    if (!v) return v;
    return *v * scale;
}

// A time in s or a frequency in Hz (DECISIONS.md Unit 61, RD1): a number, optionally followed by a
// unit, s ms us ns ps fs for a time and Hz kHz MHz GHz THz for a frequency. The unit shifts the
// text's decimal exponent before the number is read, so 2.5ns is the double 2.5e-9 is (a product
// 2.5 * 1e-9 can differ in the last bit). A unit of the other kind, or an unknown one, is an error.
enum class Dimension { time, frequency };

inline std::expected<long long, std::string> integer(std::string_view s);

namespace unit_detail {
struct Unit {
    std::string_view name;
    int exponent;
};
inline constexpr Unit times[] = {{"s", 0}, {"ms", -3}, {"us", -6}, {"ns", -9}, {"ps", -12}, {"fs", -15}};
inline constexpr Unit frequencies[] = {{"Hz", 0}, {"kHz", 3}, {"MHz", 6}, {"GHz", 9}, {"THz", 12}};
inline const Unit* find(std::span<const Unit> table, std::string_view name) {
    for (const Unit& u : table)
        if (u.name == name) return &u;
    return nullptr;
}
}  // namespace unit_detail

inline std::expected<double, std::string> quantity(std::string_view s, Dimension d) {
    s = trimmed(s);
    std::size_t cut = s.size();
    while (cut > 0 && std::isalpha(static_cast<unsigned char>(s[cut - 1]))) --cut;
    const std::string_view unit = s.substr(cut);
    if (unit.empty()) return number(s);
    std::string_view text = trimmed(s.substr(0, cut));
    using Units = std::span<const unit_detail::Unit>;
    const Units own = d == Dimension::time ? Units(unit_detail::times) : Units(unit_detail::frequencies);
    const Units other = d == Dimension::time ? Units(unit_detail::frequencies) : Units(unit_detail::times);
    const unit_detail::Unit* u = unit_detail::find(own, unit);
    if (u == nullptr) {
        const char* want = d == Dimension::time ? "a time" : "a frequency";
        if (unit_detail::find(other, unit) != nullptr) {
            return std::unexpected("'" + std::string(s) + "' is " +
                                   (d == Dimension::time ? "a frequency" : "a time") + ", not " + want);
        }
        return std::unexpected("'" + std::string(s) + "': unknown unit '" + std::string(unit) + "' for " + want +
                               (d == Dimension::time ? " (s, ms, us, ns, ps, fs)" : " (Hz, kHz, MHz, GHz, THz)"));
    }
    if (u->exponent == 0) return number(text);
    // Split off an exponent of the text's own (2.5e-3ns) and add the unit's to it.
    long long exponent = u->exponent;
    const std::size_t e = text.find_first_of("eE");
    std::string_view mantissa = text;
    if (e != std::string_view::npos) {
        std::string_view digits = text.substr(e + 1);
        if (digits.starts_with("+")) digits.remove_prefix(1);
        const auto own_exponent = integer(digits);
        if (!own_exponent) return std::unexpected("'" + std::string(s) + "' is not a number");
        exponent += *own_exponent;
        mantissa = text.substr(0, e);
    }
    if (mantissa.empty() || mantissa.find_first_of("eE") != std::string_view::npos) {
        return std::unexpected("'" + std::string(s) + "' is not a number");
    }
    auto v = number(std::string(mantissa) + "e" + std::to_string(exponent));
    if (!v) return std::unexpected("'" + std::string(s) + "' is not a number");
    return v;
}

inline std::expected<double, std::string> seconds(std::string_view s) { return quantity(s, Dimension::time); }
inline std::expected<double, std::string> hertz(std::string_view s) {
    return quantity(s, Dimension::frequency);
}

// A number with no unit, where a time or frequency unit would be a mistake (a waveform's voltage).
inline std::expected<double, std::string> voltage(std::string_view s) {
    auto v = number(s);
    if (v) return v;
    if (seconds(s) || hertz(s)) return std::unexpected("'" + std::string(trimmed(s)) + "' has a unit; a voltage is in V, with none");
    return v;
}

// An integer (no fraction or exponent).
inline std::expected<long long, std::string> integer(std::string_view s) {
    s = trimmed(s);
    long long v = 0;
    const auto [end, ec] = std::from_chars(s.data(), s.data() + s.size(), v);
    if (s.empty() || ec != std::errc{} || end != s.data() + s.size()) {
        return std::unexpected("'" + std::string(s) + "' is not an integer");
    }
    return v;
}

// on/off (also yes/no, true/false).
inline std::expected<bool, std::string> boolean(std::string_view s) {
    s = trimmed(s);
    if (s == "on" || s == "yes" || s == "true") return true;
    if (s == "off" || s == "no" || s == "false") return false;
    return std::unexpected("'" + std::string(s) + "' is not on or off");
}

// logspace(a, b, n) (Unit 43): n values from a to b evenly spaced in log10, a and b exact, both
// positive; each bound read by `read`.
template <class Read>
std::expected<std::vector<double>, std::string> logspace(std::string_view s, Read&& read) {
    if (!s.ends_with(")")) return std::unexpected("logspace( without ')'");
    const auto args = items(s.substr(9, s.size() - 10));
    if (args.size() != 3) return std::unexpected("logspace takes (first, last, count)");
    auto a = read(args[0]);
    auto b = read(args[1]);
    if (!a) return std::unexpected(a.error());
    if (!b) return std::unexpected(b.error());
    const auto n = integer(args[2]);
    if (!n) return std::unexpected(n.error());
    if (*a <= 0.0 || *b <= 0.0) return std::unexpected("logspace needs positive bounds");
    if (*n < 2 || *n > 1'000'000) return std::unexpected("logspace needs a count from 2 to a million");
    std::vector<double> out(static_cast<std::size_t>(*n));
    const double la = std::log10(*a), lb = std::log10(*b);
    for (long long k = 0; k < *n; ++k) {
        out[static_cast<std::size_t>(k)] =
            std::pow(10.0, la + (lb - la) * static_cast<double>(k) / static_cast<double>(*n - 1));
    }
    out.front() = *a;
    out.back() = *b;
    return out;
}

// Values from `a : b : step` (from a to b inclusive, a + k step; b is included when it lies on the
// grid to 1e-9 of a step), logspace(a, b, n) or a comma list; each read by `read`.
template <class Read>
std::expected<std::vector<double>, std::string> series(std::string_view s, Read&& read) {
    std::vector<double> out;
    s = trimmed(s);
    if (s.starts_with("logspace(")) return logspace(s, read);
    if (s.find(':') != std::string_view::npos) {
        const auto parts = items(s, ':');
        if (parts.size() != 3) return std::unexpected("a range is 'first : last : step'");
        auto a = read(parts[0]);
        auto b = read(parts[1]);
        auto h = read(parts[2]);
        if (!a) return std::unexpected(a.error());
        if (!b) return std::unexpected(b.error());
        if (!h) return std::unexpected(h.error());
        if (*h == 0.0 || (*b - *a) / *h < 0.0) {
            return std::unexpected("the step of a range must be non-zero and point from first to last");
        }
        const double steps = std::floor((*b - *a) / *h + 1e-9);
        if (steps > 1e6) return std::unexpected("a range of more than a million values");
        for (long long k = 0; k <= static_cast<long long>(steps); ++k) {
            const double v = *a + static_cast<double>(k) * *h;
            // The last value on the grid is `last` exactly, not its rounded multiple.
            out.push_back(std::abs(v - *b) <= 1e-9 * std::abs(*h) ? *b : v);
        }
        return out;
    }
    for (const std::string_view item : items(s)) {
        auto v = read(item);
        if (!v) return std::unexpected(v.error());
        out.push_back(*v);
    }
    return out;
}

}  // namespace NiTCAD::deck::detail
