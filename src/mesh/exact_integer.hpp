// Signed arbitrary-precision integers: the exact stage of the geometric predicates and of the dual
// measures (Unit 16; DECISIONS.md M1, M4). Internal to the mesh layer.
//
// A finite double is m * 2^e with an integer mantissa m (|m| < 2^53). A set of doubles scaled by a
// common power of two 2^-e0 (e0 the smallest exponent among them) is a set of integers, exactly; a
// polynomial in them is then an integer, computed here without rounding, overflow or underflow.
// Only sign() is exact; to_double() rounds. Speed is not a goal: this runs only when the
// floating-point filter cannot decide, or once per element for the dual measures.
#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "NiTCAD/mesh/mesh.hpp"

namespace NiTCAD::mesh::detail {

class ExactInteger {
public:
    ExactInteger() = default;  // zero
    explicit ExactInteger(std::int64_t value);

    // (negative ? -1 : 1) * magnitude * 2^shift. Precondition (NITCAD_EXPECTS): shift >= 0.
    [[nodiscard]] static ExactInteger from_shifted(std::uint64_t magnitude, bool negative,
                                                   int shift);

    [[nodiscard]] int sign() const noexcept;
    // value * 2^exp2, rounded to a double within one unit in the last place (0 for zero; may
    // overflow to infinity or underflow to zero only if the result is outside double's range).
    [[nodiscard]] double to_double(int exp2 = 0) const;

    friend ExactInteger operator+(const ExactInteger& a, const ExactInteger& b);
    friend ExactInteger operator-(const ExactInteger& a, const ExactInteger& b);
    friend ExactInteger operator*(const ExactInteger& a, const ExactInteger& b);
    friend ExactInteger operator-(ExactInteger a);

private:
    using Limbs = std::vector<std::uint32_t>;  // magnitude, least significant first, no high zeros

    static int compare(const Limbs& a, const Limbs& b) noexcept;
    static Limbs add(const Limbs& a, const Limbs& b);
    static Limbs subtract(const Limbs& a, const Limbs& b);  // requires a >= b
    static Limbs multiply(const Limbs& a, const Limbs& b);
    static ExactInteger signed_sum(const ExactInteger& a, bool b_negative, const Limbs& b);
    static void trim(Limbs& v) noexcept;

    bool negative_ = false;  // never set for zero
    Limbs limbs_;
};

// The exponent e of x = m * 2^e with an integer mantissa |m| < 2^53; a large sentinel for 0.
// Precondition (NITCAD_EXPECTS): x finite.
[[nodiscard]] int exponent_of(double x);
// x * 2^-e0 as an exact integer. Precondition (NITCAD_EXPECTS): x finite, and x == 0 or
// exponent_of(x) >= e0.
[[nodiscard]] ExactInteger scaled_integer(double x, int e0);

// The first `dimension` coordinates of K points as exact integers, all scaled by 2^-e0 with e0 the
// smallest exponent among them (0 if every coordinate is 0); the rest are 0. A homogeneous
// polynomial of degree g in these is the polynomial in the coordinates times 2^(-g e0), so it has
// the same sign, and to_double(g * e0) recovers its value.
template <std::size_t K>
struct IntegerPoints {
    std::array<std::array<ExactInteger, 3>, K> c;
    int e0 = 0;
};

template <std::size_t K>
[[nodiscard]] IntegerPoints<K> to_integers(const std::array<Point, K>& points, int dimension) {
    const auto dim = static_cast<std::size_t>(dimension);
    int e0 = 0;
    bool any = false;
    for (const Point& p : points) {
        for (std::size_t c = 0; c < dim; ++c) {
            if (p[c] == 0.0) continue;
            const int e = exponent_of(p[c]);
            e0 = any ? (e < e0 ? e : e0) : e;
            any = true;
        }
    }
    IntegerPoints<K> r;
    r.e0 = e0;
    for (std::size_t k = 0; k < K; ++k) {
        for (std::size_t c = 0; c < dim; ++c) r.c[k][c] = scaled_integer(points[k][c], e0);
    }
    return r;
}

}  // namespace NiTCAD::mesh::detail
