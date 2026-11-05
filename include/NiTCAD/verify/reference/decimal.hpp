// Exact reading of a printed decimal ("-4.609e-10", "1.17", "636") for the TVB-31 parameter register (Unit 55,
// TVR-3; DECISIONS 15.9). Two results:
//   - to_double: the correctly rounded double (std::from_chars, round to nearest), which the table check compares
//     bit for bit with the shipped parameter;
//   - to_dd: the value M 10^e as a double-double, M the digit string as an integer (exact for at most 15 digits)
//     and 10^|e| exact (codata2018's exact_power_of_ten), so the only rounding is one double-double product or
//     quotient. The references evaluate in this.
// Base-free (TVD-19): standard library and verify/reference only.
#pragma once

#include <cctype>
#include <charconv>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

#include "NiTCAD/verify/reference/codata2018.hpp"
#include "NiTCAD/verify/reference/double_double.hpp"

namespace NiTCAD::verify::reference::decimal {

struct Parsed {
    bool negative = false;
    std::string digits;  // without the point or leading zeros' meaning changed
    int exponent = 0;    // value = digits * 10^exponent
};

[[nodiscard]] inline Parsed parse(std::string_view s) {
    Parsed p;
    std::size_t i = 0;
    if (i < s.size() && (s[i] == '-' || s[i] == '+')) p.negative = s[i++] == '-';
    int after_point = 0;
    bool point = false;
    for (; i < s.size() && s[i] != 'e' && s[i] != 'E'; ++i) {
        const char ch = s[i];
        if (std::isdigit(static_cast<unsigned char>(ch)) != 0) {
            p.digits += ch;
            if (point) ++after_point;
        } else if (ch == '.' && !point) {
            point = true;
        } else {
            throw std::invalid_argument("decimal: unexpected character in " + std::string(s));
        }
    }
    if (p.digits.empty()) throw std::invalid_argument("decimal: no digits in " + std::string(s));
    int e = 0;
    if (i < s.size()) {
        std::string_view rest = s.substr(i + 1);
        if (!rest.empty() && rest.front() == '+') rest.remove_prefix(1);
        const auto r = std::from_chars(rest.data(), rest.data() + rest.size(), e);
        if (r.ec != std::errc{} || r.ptr != rest.data() + rest.size()) {
            throw std::invalid_argument("decimal: bad exponent in " + std::string(s));
        }
    }
    p.exponent = e - after_point;
    return p;
}

[[nodiscard]] inline double to_double(std::string_view s) {
    double v = 0.0;
    const auto r = std::from_chars(s.data(), s.data() + s.size(), v);
    if (r.ec != std::errc{} || r.ptr != s.data() + s.size()) {
        throw std::invalid_argument("decimal: from_chars refused " + std::string(s));
    }
    return v;
}

[[nodiscard]] inline DD to_dd(std::string_view s) {
    const Parsed p = parse(s);
    std::size_t first = p.digits.find_first_not_of('0');
    if (first == std::string::npos) return DD(0.0);
    const std::string digits = p.digits.substr(first);
    if (digits.size() > 15) throw std::invalid_argument("decimal: more than 15 significant digits");
    const double m = codata2018::to_double(digits);  // an integer below 10^15: exact
    const DD v = p.exponent >= 0 ? DD(m) * codata2018::exact_power_of_ten(p.exponent)
                                 : DD(m) / codata2018::exact_power_of_ten(-p.exponent);
    return p.negative ? -v : v;
}

}  // namespace NiTCAD::verify::reference::decimal
