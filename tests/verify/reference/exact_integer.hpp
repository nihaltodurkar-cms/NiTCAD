// Exact big-integer arithmetic for the verify tests (Unit 55, TVR-1): checks that must not themselves round,
// such as the error-free transformations (TC1) and the comparison of a double-double with a printed decimal
// literal (TC1 constants, TC3 nodes and weights). Test-only; not part of any library.
#pragma once

#include "NiTCAD/verify/reference/double_double.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string_view>
#include <vector>

namespace exact {

// Non-negative big integer in base 2^32 (little-endian limbs).
struct Big {
    std::vector<std::uint32_t> w;

    static Big from_u64(std::uint64_t v) {
        Big b;
        while (v != 0) {
            b.w.push_back(static_cast<std::uint32_t>(v));
            v >>= 32;
        }
        return b;
    }
    void trim() {
        while (!w.empty() && w.back() == 0) w.pop_back();
    }
    void mul_small(std::uint32_t m) {
        std::uint64_t carry = 0;
        for (auto& limb : w) {
            const std::uint64_t t = static_cast<std::uint64_t>(limb) * m + carry;
            limb = static_cast<std::uint32_t>(t);
            carry = t >> 32;
        }
        if (carry != 0) w.push_back(static_cast<std::uint32_t>(carry));
    }
    void add_small(std::uint32_t a) {
        std::uint64_t carry = a;
        for (std::size_t i = 0; carry != 0; ++i) {
            if (i == w.size()) w.push_back(0);
            const std::uint64_t t = static_cast<std::uint64_t>(w[i]) + carry;
            w[i] = static_cast<std::uint32_t>(t);
            carry = t >> 32;
        }
    }
    void shl(int bits) {
        if (w.empty()) return;
        const int limbs = bits / 32;
        const int rest = bits % 32;
        w.insert(w.begin(), static_cast<std::size_t>(limbs), 0u);
        if (rest != 0) {
            std::uint32_t carry = 0;
            for (std::size_t i = static_cast<std::size_t>(limbs); i < w.size(); ++i) {
                const std::uint32_t next = w[i] >> (32 - rest);
                w[i] = (w[i] << rest) | carry;
                carry = next;
            }
            if (carry != 0) w.push_back(carry);
        }
    }
    [[nodiscard]] double to_double() const {
        double v = 0.0;
        for (std::size_t i = w.size(); i-- > 0;) v = v * 4294967296.0 + w[i];
        return v;
    }
};

inline int cmp(const Big& a, const Big& b) {
    if (a.w.size() != b.w.size()) return a.w.size() < b.w.size() ? -1 : 1;
    for (std::size_t i = a.w.size(); i-- > 0;)
        if (a.w[i] != b.w[i]) return a.w[i] < b.w[i] ? -1 : 1;
    return 0;
}
inline Big add(const Big& a, const Big& b) {
    Big r;
    std::uint64_t carry = 0;
    for (std::size_t i = 0; i < std::max(a.w.size(), b.w.size()) || carry != 0; ++i) {
        const std::uint64_t t = carry + (i < a.w.size() ? a.w[i] : 0u) + (i < b.w.size() ? b.w[i] : 0u);
        r.w.push_back(static_cast<std::uint32_t>(t));
        carry = t >> 32;
    }
    r.trim();
    return r;
}
// a - b for a >= b.
inline Big sub(const Big& a, const Big& b) {
    Big r = a;
    std::int64_t borrow = 0;
    for (std::size_t i = 0; i < r.w.size(); ++i) {
        std::int64_t t = static_cast<std::int64_t>(r.w[i]) - borrow - (i < b.w.size() ? b.w[i] : 0u);
        borrow = t < 0 ? 1 : 0;
        if (t < 0) t += 4294967296LL;
        r.w[i] = static_cast<std::uint32_t>(t);
    }
    r.trim();
    return r;
}

// Signed big integer.
struct SBig {
    Big mag;
    bool neg = false;
};
inline SBig sadd(const SBig& a, const SBig& b) {
    if (a.neg == b.neg) return {add(a.mag, b.mag), a.neg};
    const int c = cmp(a.mag, b.mag);
    if (c == 0) return {};
    return c > 0 ? SBig{sub(a.mag, b.mag), a.neg} : SBig{sub(b.mag, a.mag), b.neg};
}
inline bool operator==(const SBig& a, const SBig& b) {
    if (a.mag.w.empty() && b.mag.w.empty()) return true;
    return a.neg == b.neg && cmp(a.mag, b.mag) == 0;
}

// d = M 2^k with M an integer of at most 53 bits (M = 0 for d = 0).
struct Scaled {
    std::int64_t m = 0;
    int k = 0;
};
inline Scaled decompose(double d) {
    if (d == 0.0) return {};
    int k = 0;
    const double f = std::frexp(d, &k);
    Scaled s{static_cast<std::int64_t>(std::ldexp(f, 53)), k - 53};
    while (s.m % 2 == 0) {  // odd M: the largest k with d on the grid 2^k
        s.m /= 2;
        ++s.k;
    }
    return s;
}
// M 2^(k - base) as a signed big integer, k >= base.
inline SBig to_big(std::int64_t m, int k, int base) {
    SBig r{Big::from_u64(static_cast<std::uint64_t>(m < 0 ? -m : m)), m < 0};
    r.mag.shl(k - base);
    return r;
}
// |m1 m2| exactly.
inline Big product(std::int64_t m1, std::int64_t m2) {
    Big r = Big::from_u64(static_cast<std::uint64_t>(m1 < 0 ? -m1 : m1));
    const auto b = static_cast<std::uint64_t>(m2 < 0 ? -m2 : m2);
    // r * b by two 32-bit halves.
    Big lo = r;
    lo.mul_small(static_cast<std::uint32_t>(b));
    Big hi = r;
    hi.mul_small(static_cast<std::uint32_t>(b >> 32));
    hi.shl(32);
    return add(lo, hi);
}

// A double-double against a decimal literal ("-0.123...", digits with at most one '.' and an optional leading
// '-'), exactly: with c = C 2^E and L = N / 10^d, |c - L| <= tol_udd 2^-104 |L| iff
// 2^104 |C 10^d - N 2^-E| <= tol_udd N 2^-E (E < 0 is arranged by shifting). A zero literal requires c = 0.
// `measured_udd` is |c - L| / |L| in u_dd, approximate (for the record).
struct LiteralComparison {
    bool within = false;
    double measured_udd = 0.0;
};
inline LiteralComparison compare_with_literal(const NiTCAD::verify::reference::DD& c, std::string_view literal,
                                              std::uint32_t tol_udd) {
    bool lit_neg = false;
    if (!literal.empty() && literal.front() == '-') {
        lit_neg = true;
        literal.remove_prefix(1);
    }
    Big n;
    int d = 0;
    bool after_point = false;
    for (const char ch : literal) {
        if (ch == '.') {
            after_point = true;
            continue;
        }
        n.mul_small(10);
        n.add_small(static_cast<std::uint32_t>(ch - '0'));
        if (after_point) ++d;
    }
    n.trim();
    if (n.w.empty()) return {c.hi == 0.0 && c.lo == 0.0, (c.hi == 0.0 && c.lo == 0.0) ? 0.0 : 1e300};
    const Scaled h = decompose(c.hi);
    const Scaled l = decompose(c.lo);
    int base = std::min(c.hi == 0.0 ? 0 : h.k, c.lo == 0.0 ? 0 : l.k);
    base = std::min(base, -1);
    SBig cbig = sadd(to_big(h.m, h.k, base), to_big(l.m, l.k, base));
    for (int i = 0; i < d; ++i) cbig.mag.mul_small(10);  // C 10^d
    SBig nbig{n, lit_neg};
    nbig.mag.shl(-base);  // N 2^-E
    SBig neg_n = nbig;
    neg_n.neg = !neg_n.neg;
    Big diff = sadd(cbig, neg_n).mag;
    LiteralComparison r;
    r.measured_udd = diff.to_double() / nbig.mag.to_double() / NiTCAD::verify::reference::u_dd;
    diff.shl(104);
    Big bound = nbig.mag;
    bound.mul_small(tol_udd);
    r.within = cmp(diff, bound) <= 0;
    return r;
}

}  // namespace exact
