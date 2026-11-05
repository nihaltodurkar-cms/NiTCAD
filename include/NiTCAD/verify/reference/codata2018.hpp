// CODATA 2018 constants as NIST prints them, for TVB-01 (Unit 55, TVR-1; DECISIONS 11.8, 11.13 TB-Q1-Q3).
//
// Source (provenance PA): NIST, "Fundamental Physical Constants --- Complete Listing, 2018 CODATA adjustment",
// https://physics.nist.gov/cuu/Constants/ArchiveASCII/allascii_2018.txt (read 2026-11-05; the file's header names
// the 2018 adjustment). Each string below is the table's value column typed as printed, digit groups and "..."
// included; only the column padding is dropped. This is the one place in verify/reference that holds its own
// physical constants (DECISIONS 11.2: TVB-01's reference must not read base's).
//
// Two conversions, neither of which uses base or the compiler's reading of a literal:
//   - shifted_decimal: the printed digits with the decimal exponent moved by a whole number (an SI value in
//     F m^-1 becomes F cm^-1 by -2), as a string; exact, since only the exponent's digits change;
//   - to_dd: the printed digits as the exact rational M 10^e, M the digit string as an integer (exact in a double
//     for at most 15 digits) and 10^|e| held exactly as a double-double (10^n = 10^a 10^b with a, b <= 22, each
//     exact in a double, and their product exact in two doubles while 5^n has at most 106 bits, n <= 44), so the
//     only rounding is the one double-double quotient or product (a few u_dd, TC1).
// A printed value ending in "..." is a truncation of an exact constant (the reduced Planck constant and k in
// eV/K): its digits are not the constant's value, and TVB-01 computes those from the exact h, k and q instead.
#pragma once

#include <cctype>
#include <charconv>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

#include "NiTCAD/verify/reference/double_double.hpp"

namespace NiTCAD::verify::reference::codata2018 {

struct PrintedConstant {
    std::string_view quantity;     // the table's quantity column
    std::string_view value;        // the value column as printed
    std::string_view uncertainty;  // "(exact)" or the printed standard uncertainty
    std::string_view unit;
};

inline constexpr PrintedConstant elementary_charge{"elementary charge", "1.602 176 634 e-19", "(exact)", "C"};
inline constexpr PrintedConstant boltzmann{"Boltzmann constant", "1.380 649 e-23", "(exact)", "J K^-1"};
inline constexpr PrintedConstant boltzmann_eV{"Boltzmann constant in eV/K", "8.617 333 262... e-5", "(exact)",
                                              "eV K^-1"};
inline constexpr PrintedConstant vacuum_permittivity{"vacuum electric permittivity", "8.854 187 8128 e-12",
                                                     "0.000 000 0013 e-12", "F m^-1"};
inline constexpr PrintedConstant planck{"Planck constant", "6.626 070 15 e-34", "(exact)", "J Hz^-1"};
inline constexpr PrintedConstant reduced_planck{"reduced Planck constant", "1.054 571 817... e-34", "(exact)",
                                                "J s"};
inline constexpr PrintedConstant electron_mass{"electron mass", "9.109 383 7015 e-31", "0.000 000 0028 e-31", "kg"};

[[nodiscard]] constexpr bool is_truncated(const PrintedConstant& c) { return c.value.find("...") != c.value.npos; }

// The printed value split into its digits (no point), the number of digits before the point and the decimal
// exponent: value = digits * 10^(exponent - (size - point)), digits read as an integer.
struct Decimal {
    std::string digits;  // "1602176634"
    int point = 0;       // digits before the decimal point
    int exponent = 0;    // the printed e-exponent
};

[[nodiscard]] inline Decimal parse(std::string_view printed) {
    Decimal d;
    std::size_t i = 0;
    bool seen_point = false;
    for (; i < printed.size() && printed[i] != 'e'; ++i) {
        const char ch = printed[i];
        if (std::isdigit(static_cast<unsigned char>(ch)) != 0) {
            d.digits += ch;
            if (!seen_point) ++d.point;
        } else if (ch == '.') {
            if (!seen_point) {
                seen_point = true;
            }  // the further points of "..." mark a truncation and carry no digit
        } else if (ch != ' ') {
            throw std::invalid_argument("codata2018: unexpected character in a printed value");
        }
    }
    if (i == printed.size() || d.digits.empty()) throw std::invalid_argument("codata2018: no exponent or digits");
    const std::string_view e = printed.substr(i + 1);
    const auto r = std::from_chars(e.data(), e.data() + e.size(), d.exponent);
    if (r.ec != std::errc{} || r.ptr != e.data() + e.size()) {
        throw std::invalid_argument("codata2018: bad exponent");
    }
    return d;
}

// The printed digits with the exponent moved by `shift`, in the form "d.ddde-N" that std::from_chars reads.
[[nodiscard]] inline std::string shifted_decimal(const PrintedConstant& c, int shift = 0) {
    const Decimal d = parse(c.value);
    std::string s = d.digits.substr(0, static_cast<std::size_t>(d.point));
    if (static_cast<std::size_t>(d.point) < d.digits.size()) s += "." + d.digits.substr(static_cast<std::size_t>(d.point));
    return s + "e" + std::to_string(d.exponent + shift);
}

// The correctly rounded double of a decimal string (std::from_chars, round to nearest).
[[nodiscard]] inline double to_double(const std::string& decimal) {
    double v = 0.0;
    const auto r = std::from_chars(decimal.data(), decimal.data() + decimal.size(), v);
    if (r.ec != std::errc{} || r.ptr != decimal.data() + decimal.size()) {
        throw std::invalid_argument("codata2018: from_chars refused " + decimal);
    }
    return v;
}

// 10^n exactly as a double-double, 0 <= n <= 44.
[[nodiscard]] inline DD exact_power_of_ten(int n) {
    if (n < 0 || n > 44) throw std::invalid_argument("codata2018: 10^n exact only for 0 <= n <= 44");
    double p[23];
    p[0] = 1.0;
    for (int k = 1; k <= 22; ++k) p[k] = p[k - 1] * 10.0;  // exact: 10^k < 2^53 * 2^k
    if (n <= 22) return DD(p[n]);
    return two_prod(p[22], p[n - 22]);
}

// The printed value (shifted by `shift` decades) as a double-double: M 10^e, rounded once.
[[nodiscard]] inline DD to_dd(const PrintedConstant& c, int shift = 0) {
    const Decimal d = parse(c.value);
    if (d.digits.size() > 15) throw std::invalid_argument("codata2018: more than 15 digits");
    const double m = to_double(d.digits);  // an integer below 10^15: exact
    const int e = d.exponent + shift - static_cast<int>(d.digits.size() - static_cast<std::size_t>(d.point));
    return e >= 0 ? DD(m) * exact_power_of_ten(e) : DD(m) / exact_power_of_ten(-e);
}

}  // namespace NiTCAD::verify::reference::codata2018
