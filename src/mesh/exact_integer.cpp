#include "exact_integer.hpp"

#include <bit>
#include <cmath>
#include <cstddef>
#include <utility>

#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::mesh::detail {

ExactInteger::ExactInteger(std::int64_t value) {
    if (value == 0) return;
    negative_ = value < 0;
    // Two's complement: the magnitude of INT64_MIN is representable as uint64.
    const std::uint64_t m = negative_ ? ~static_cast<std::uint64_t>(value) + 1
                                      : static_cast<std::uint64_t>(value);
    limbs_ = {static_cast<std::uint32_t>(m), static_cast<std::uint32_t>(m >> 32)};
    trim(limbs_);
}

ExactInteger ExactInteger::from_shifted(std::uint64_t magnitude, bool negative, int shift) {
    NITCAD_EXPECTS(shift >= 0);
    ExactInteger r;
    if (magnitude == 0) return r;
    const auto words = static_cast<std::size_t>(shift / 32);
    const int bits = shift % 32;
    r.limbs_.assign(words, 0);
    // Up to 64 + 31 bits after the in-word shift: three limbs.
    const std::uint64_t lo = magnitude << bits;
    const std::uint64_t hi = bits == 0 ? 0 : magnitude >> (64 - bits);
    r.limbs_.push_back(static_cast<std::uint32_t>(lo));
    r.limbs_.push_back(static_cast<std::uint32_t>(lo >> 32));
    r.limbs_.push_back(static_cast<std::uint32_t>(hi));
    trim(r.limbs_);
    r.negative_ = negative;
    return r;
}

int ExactInteger::sign() const noexcept {
    if (limbs_.empty()) return 0;
    return negative_ ? -1 : 1;
}

double ExactInteger::to_double(int exp2) const {
    if (limbs_.empty()) return 0.0;
    const int bits = static_cast<int>(32 * (limbs_.size() - 1)) +
                     static_cast<int>(std::bit_width(limbs_.back()));
    // The top 64 bits, truncated; converting them rounds once, so the result is within one ulp.
    const int drop = bits > 64 ? bits - 64 : 0;
    std::uint64_t top = 0;
    for (int b = bits - 1; b >= drop; --b) {
        const std::uint32_t limb = limbs_[static_cast<std::size_t>(b / 32)];
        top = (top << 1) | ((limb >> (b % 32)) & 1u);
    }
    const double magnitude = std::ldexp(static_cast<double>(top), drop + exp2);
    return negative_ ? -magnitude : magnitude;
}

int ExactInteger::compare(const Limbs& a, const Limbs& b) noexcept {
    if (a.size() != b.size()) return a.size() < b.size() ? -1 : 1;
    for (std::size_t i = a.size(); i-- > 0;) {
        if (a[i] != b[i]) return a[i] < b[i] ? -1 : 1;
    }
    return 0;
}

ExactInteger::Limbs ExactInteger::add(const Limbs& a, const Limbs& b) {
    const Limbs& longer = a.size() >= b.size() ? a : b;
    const Limbs& shorter = a.size() >= b.size() ? b : a;
    Limbs r(longer.size() + 1, 0);
    std::uint64_t carry = 0;
    for (std::size_t i = 0; i < longer.size(); ++i) {
        const std::uint64_t s = static_cast<std::uint64_t>(longer[i]) +
                                (i < shorter.size() ? shorter[i] : 0u) + carry;
        r[i] = static_cast<std::uint32_t>(s);
        carry = s >> 32;
    }
    r[longer.size()] = static_cast<std::uint32_t>(carry);
    trim(r);
    return r;
}

ExactInteger::Limbs ExactInteger::subtract(const Limbs& a, const Limbs& b) {
    Limbs r(a.size(), 0);
    std::int64_t borrow = 0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        std::int64_t d = static_cast<std::int64_t>(a[i]) -
                         (i < b.size() ? static_cast<std::int64_t>(b[i]) : 0) - borrow;
        borrow = d < 0 ? 1 : 0;
        if (d < 0) d += std::int64_t{1} << 32;
        r[i] = static_cast<std::uint32_t>(d);
    }
    trim(r);
    return r;
}

ExactInteger::Limbs ExactInteger::multiply(const Limbs& a, const Limbs& b) {
    if (a.empty() || b.empty()) return {};
    Limbs r(a.size() + b.size(), 0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        std::uint64_t carry = 0;
        for (std::size_t j = 0; j < b.size(); ++j) {
            const std::uint64_t t = static_cast<std::uint64_t>(a[i]) * b[j] + r[i + j] + carry;
            r[i + j] = static_cast<std::uint32_t>(t);
            carry = t >> 32;
        }
        r[i + b.size()] = static_cast<std::uint32_t>(carry);
    }
    trim(r);
    return r;
}

void ExactInteger::trim(Limbs& v) noexcept {
    while (!v.empty() && v.back() == 0) v.pop_back();
}

ExactInteger ExactInteger::signed_sum(const ExactInteger& a, bool b_negative, const Limbs& b) {
    ExactInteger r;
    if (a.negative_ == b_negative) {
        r.limbs_ = add(a.limbs_, b);
        r.negative_ = a.negative_;
    } else {
        const int c = compare(a.limbs_, b);
        if (c == 0) return r;
        r.limbs_ = c > 0 ? subtract(a.limbs_, b) : subtract(b, a.limbs_);
        r.negative_ = c > 0 ? a.negative_ : b_negative;
    }
    if (r.limbs_.empty()) r.negative_ = false;
    return r;
}

ExactInteger operator+(const ExactInteger& a, const ExactInteger& b) {
    return ExactInteger::signed_sum(a, b.negative_, b.limbs_);
}

ExactInteger operator-(const ExactInteger& a, const ExactInteger& b) {
    return ExactInteger::signed_sum(a, !b.negative_ && !b.limbs_.empty(), b.limbs_);
}

ExactInteger operator*(const ExactInteger& a, const ExactInteger& b) {
    ExactInteger r;
    r.limbs_ = ExactInteger::multiply(a.limbs_, b.limbs_);
    r.negative_ = !r.limbs_.empty() && (a.negative_ != b.negative_);
    return r;
}

ExactInteger operator-(ExactInteger a) {
    if (!a.limbs_.empty()) a.negative_ = !a.negative_;
    return a;
}

int exponent_of(double x) {
    NITCAD_EXPECTS(std::isfinite(x));
    if (x == 0.0) return 1 << 20;
    int e = 0;
    (void)std::frexp(x, &e);  // x = f 2^e, 0.5 <= |f| < 1, so m = f 2^53 is an integer
    return e - 53;
}

ExactInteger scaled_integer(double x, int e0) {
    NITCAD_EXPECTS(std::isfinite(x));
    if (x == 0.0) return {};
    int e = 0;
    const double f = std::frexp(x, &e);
    const auto m = static_cast<std::int64_t>(std::ldexp(f, 53));  // exact: |f 2^53| < 2^53
    const int shift = (e - 53) - e0;
    NITCAD_EXPECTS(shift >= 0);
    const std::uint64_t magnitude = m < 0 ? static_cast<std::uint64_t>(-m)
                                          : static_cast<std::uint64_t>(m);
    return ExactInteger::from_shifted(magnitude, m < 0, shift);
}

}  // namespace NiTCAD::mesh::detail
