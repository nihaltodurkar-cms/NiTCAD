// Exact predicates and the symbolic tie-break (Unit 16; DECISIONS.md M4). Expected signs come from
// configurations whose exact sign is known in closed form (dyadic or integer coordinates), never
// from another floating-point evaluation.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include "NiTCAD/mesh/mesh.hpp"
#include "exact_integer.hpp"
#include "predicates.hpp"

using NiTCAD::mesh::NodeId;
using NiTCAD::mesh::Point;
namespace pr = NiTCAD::mesh::predicates;
using NiTCAD::mesh::detail::ExactInteger;

namespace {

constexpr double u53 = 0x1p-53;

Point scaled(const Point& p, double s) { return {p[0] * s, p[1] * s, p[2] * s}; }
Point shifted(const Point& p, const Point& d) { return {p[0] + d[0], p[1] + d[1], p[2] + d[2]}; }

int sgn(long v) { return v > 0 ? 1 : (v < 0 ? -1 : 0); }

// Parity of a permutation: +1 even, -1 odd.
template <std::size_t K>
int parity(const std::array<std::size_t, K>& p) {
    int s = 1;
    for (std::size_t i = 0; i < K; ++i) {
        for (std::size_t j = i + 1; j < K; ++j) {
            if (p[i] > p[j]) s = -s;
        }
    }
    return s;
}

// A small deterministic generator (splitmix64), so the sampled cases are the same on every run.
struct SplitMix {
    std::uint64_t state;
    std::uint64_t next() {
        std::uint64_t z = (state += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return z ^ (z >> 31);
    }
    double uniform() { return static_cast<double>(next() >> 11) * 0x1p-53; }  // [0, 1)
    int small(int range) { return static_cast<int>(next() % static_cast<std::uint64_t>(2 * range + 1)) - range; }
};

}  // namespace

TEST_CASE("exact integers: arithmetic, sign and rounding", "[mesh][predicates]") {
    const ExactInteger a(123456789);
    const ExactInteger b(-987654321);
    CHECK((a + b).sign() == -1);
    CHECK((a - a).sign() == 0);
    CHECK((a * b).to_double() == 123456789.0 * -987654321.0);
    CHECK((-(a * b)).sign() == 1);
    // 2^200 - 1 times 2^200 + 1 = 2^400 - 1: exact through many limbs.
    const ExactInteger big = ExactInteger::from_shifted(1, false, 200);
    const ExactInteger one(1);
    const ExactInteger product = (big - one) * (big + one);
    CHECK((product - big * big + one).sign() == 0);
    CHECK(product.to_double(-400) == 1.0);
    // Scaled integers of doubles reproduce the doubles, subnormals included.
    for (const double x : {1.0, -0.1, 3e-300, -4.9e-324, 1.7e308}) {
        const int e = NiTCAD::mesh::detail::exponent_of(x);
        CHECK(NiTCAD::mesh::detail::scaled_integer(x, e).to_double(e) == x);
    }
}

TEST_CASE("exact integers: inline and heap limbs, and to_double's top 64 bits", "[mesh][predicates]") {
    // Values with a 64-bit top m (top bit set) at every shift s, plus a tail below 2^s: to_double
    // truncates to the top 64 bits, so it must give m 2^s converted once, whether the limbs fit
    // inline (24 limbs) or live on the heap (s up to 3000), and at every offset s % 32.
    std::uint64_t state = 12345;
    const auto next = [&]() {
        state += 0x9E3779B97F4A7C15ull;
        std::uint64_t z = state;
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return z ^ (z >> 31);
    };
    for (int s = 0; s <= 3000; s += (s < 128 ? 1 : 37)) {
        const std::uint64_t m = next() | (std::uint64_t{1} << 63);
        ExactInteger v = ExactInteger::from_shifted(m, false, s);
        if (s >= 64) v = v + ExactInteger::from_shifted(next(), false, s - 64);  // a tail below 2^s
        const double expected = std::ldexp(static_cast<double>(m), s - 3000);
        CHECK(v.to_double(-3000) == expected);
        CHECK((-v).to_double(-3000) == -expected);
    }
    // Small values: fewer than 64 bits, converted exactly or rounded once.
    for (const std::uint64_t m : {std::uint64_t{1}, std::uint64_t{0xFFFFFFFF}, std::uint64_t{0x1FFFFFFFFFFFFF},
                                  std::uint64_t{0x123456789ABCDEF}}) {
        CHECK(ExactInteger::from_shifted(m, false, 0).to_double() == static_cast<double>(m));
    }
    // Identities through products far beyond the inline limbs: (x + y)(x - y) = x^2 - y^2.
    ExactInteger x = ExactInteger::from_shifted(next(), false, 1500) + ExactInteger::from_shifted(next(), true, 700);
    ExactInteger y = ExactInteger::from_shifted(next(), true, 1400) + ExactInteger::from_shifted(next(), false, 3);
    CHECK(((x + y) * (x - y) - (x * x - y * y)).sign() == 0);
    CHECK((x * y + y * x - (x + x) * y).sign() == 0);
    // Back from the heap to small values: the difference of two huge equal products is zero, and
    // a huge value minus itself plus one is one.
    const ExactInteger one(1);
    CHECK((x * y - y * x + one).to_double() == 1.0);
}

TEST_CASE("orientation: signs, exact zeros and near-collinear points", "[mesh][predicates]") {
    CHECK(pr::orient2d({0, 0, 0}, {1, 0, 0}, {0, 1, 0}) == 1);
    CHECK(pr::orient2d({0, 0, 0}, {0, 1, 0}, {1, 0, 0}) == -1);
    CHECK(pr::orient2d({0.125, 0.25, 0}, {0.375, 0.75, 0}, {1.0, 2.0, 0}) == 0);
    CHECK(pr::orient3d({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}) == 1);
    CHECK(pr::orient3d({0, 0, 0}, {0, 1, 0}, {1, 0, 0}, {0, 0, 1}) == -1);
    CHECK(pr::orient3d({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0.5, 0.25, 0}) == 0);

    // q, r on the line y = x and p = (0.5 + i u, 0.5 + j u): orient(p, q, r) = (q_x - r_x)(p_x - p_y)
    // = 12 (j - i) u, so its sign is sign(j - i). Naive double evaluation gets some of these wrong.
    const Point q{12, 12, 0};
    const Point r{24, 24, 0};
    for (int i = 0; i < 32; ++i) {
        for (int j = 0; j < 32; ++j) {
            const Point p{0.5 + i * u53, 0.5 + j * u53, 0};
            INFO("i = " << i << ", j = " << j);
            CHECK(pr::orient2d(p, q, r) == sgn(j - i));
            CHECK(pr::exact::orient2d(p, q, r) == sgn(j - i));
        }
    }
    // The same in 3D: the plane z = 0 through three dyadic points, d = (x, y, k u).
    for (int k = -4; k <= 4; ++k) {
        const Point d{0.3, 0.7, k * u53};
        CHECK(pr::orient3d({0, 0, 0}, {1, 0, 0}, {0, 1, 0}, d) == sgn(k));
    }
}

TEST_CASE("orientation at extreme scales goes to the exact stage", "[mesh][predicates]") {
    // Powers of two scale exactly; the sign is scale-invariant. At 2^-960 the products underflow,
    // at 2^-530 and 2^-520 they are subnormal (rounded wrongly), at 2^900 they overflow; the filter
    // must defer in each case, and the exact stage has no range limit.
    const Point q{12, 12, 0};
    const Point r{24, 24, 0};
    for (const double s : {0x1p-960, 0x1p-530, 0x1p-520, 0x1p-500, 0x1p+500, 0x1p+900}) {
        for (int i = 0; i < 8; ++i) {
            for (int j = 0; j < 8; ++j) {
                const Point p{0.5 + i * u53, 0.5 + j * u53, 0};
                CHECK(pr::orient2d(scaled(p, s), scaled(q, s), scaled(r, s)) == sgn(j - i));
            }
        }
        const Point a = scaled({5, 0, 0}, s);
        const Point b = scaled({3, 4, 0}, s);
        const Point c = scaled({-4, 3, 0}, s);
        CHECK(pr::incircle(a, b, c, scaled({0, -5, 0}, s)) == 0);
        CHECK(pr::incircle(a, b, c, scaled({0, -5 + 0x1p-50, 0}, s)) == 1);
        CHECK(pr::incircle(a, b, c, scaled({0, -5 - 0x1p-50, 0}, s)) == -1);
    }
}

TEST_CASE("incircle: cocircular points, both orientations, translation", "[mesh][predicates]") {
    // (5, 0), (3, 4), (-4, 3) turn counterclockwise on the circle of radius 5.
    const Point a{5, 0, 0};
    const Point b{3, 4, 0};
    const Point c{-4, 3, 0};
    REQUIRE(pr::orient2d(a, b, c) == 1);
    CHECK(pr::incircle(a, b, c, {0, 5, 0}) == 0);
    CHECK(pr::incircle(a, b, c, {0, -5, 0}) == 0);
    CHECK(pr::incircle(a, b, c, {0, 0, 0}) == 1);
    CHECK(pr::incircle(a, b, c, {0, 5 - 0x1p-50, 0}) == 1);  // ulp(5) = 2^-50
    CHECK(pr::incircle(a, b, c, {0, 5 + 0x1p-50, 0}) == -1);
    CHECK(pr::incircle(b, a, c, {0, 0, 0}) == -1);  // clockwise: the sign flips
    // Translated by 1024 (exact): still cocircular, the offsets still decide.
    const Point t{1024, 1024, 0};
    CHECK(pr::incircle(shifted(a, t), shifted(b, t), shifted(c, t), shifted({0, -5, 0}, t)) == 0);
    CHECK(pr::incircle(shifted(a, t), shifted(b, t), shifted(c, t), {1024, 1029 - 0x1p-42, 0}) == 1);
    CHECK(pr::incircle(shifted(a, t), shifted(b, t), shifted(c, t), {1024, 1029 + 0x1p-42, 0}) == -1);
}

TEST_CASE("insphere: cospherical points, both orientations, translation", "[mesh][predicates]") {
    // Integer points on the sphere of radius 3: (3,0,0), (0,3,0), (0,0,3), (2,2,1), (-2,1,2), ...
    const Point a{3, 0, 0};
    const Point b{0, 3, 0};
    const Point c{0, 0, 3};
    const Point d{-2, -2, -1};
    const int o = pr::orient3d(a, b, c, d);
    REQUIRE(o != 0);
    const auto in = [&](const Point& e) {
        return o > 0 ? pr::insphere(a, b, c, d, e) : pr::insphere(b, a, c, d, e);
    };
    CHECK(in({2, 2, 1}) == 0);
    CHECK(in({-2, 1, 2}) == 0);
    CHECK(in({0, 0, 0}) == 1);
    CHECK(in({0, 0, -3 + 0x1p-51}) == 1);  // ulp(3) = 2^-51
    CHECK(in({0, 0, -3 - 0x1p-51}) == -1);
    CHECK(pr::insphere(a, b, c, d, {0, 0, 0}) == o);  // the sign follows the orientation
    const Point t{4096, -2048, 1024};
    const auto tin = [&](const Point& e) {
        const Point ta = shifted(a, t), tb = shifted(b, t), tc = shifted(c, t), td = shifted(d, t);
        return o > 0 ? pr::insphere(ta, tb, tc, td, shifted(e, t)) : pr::insphere(tb, ta, tc, td, shifted(e, t));
    };
    CHECK(tin({2, 2, 1}) == 0);
    CHECK(tin({1, 2, 2}) == 0);
}

TEST_CASE("diametral circle and sphere", "[mesh][predicates]") {
    CHECK(pr::in_diametral_circle({-1, 0, 0}, {1, 0, 0}, {0, 1, 0}) == 0);
    CHECK(pr::in_diametral_circle({-1, 0, 0}, {1, 0, 0}, {0, 1 - u53, 0}) == 1);
    CHECK(pr::in_diametral_circle({-1, 0, 0}, {1, 0, 0}, {0, 1 + 2 * u53, 0}) == -1);
    CHECK(pr::in_diametral_circle({-1, 0, 0}, {1, 0, 0}, {-1, 0, 0}) == 0);
    // The triangle's circumcircle is the unit circle in z = 0, centre at the origin.
    const Point a{1, 0, 0};
    const Point b{0, 1, 0};
    const Point c{-1, 0, 0};
    CHECK(pr::in_diametral_sphere(a, b, c, {0, 0, 1}) == 0);
    CHECK(pr::in_diametral_sphere(a, b, c, {0, -1, 0}) == 0);
    CHECK(pr::in_diametral_sphere(a, b, c, {0, 0, 1 - u53}) == 1);
    CHECK(pr::in_diametral_sphere(a, b, c, {0, 0, 1 + 2 * u53}) == -1);
    CHECK(pr::in_diametral_sphere(a, b, c, {0.5, 0.5, 0.5}) == 1);
    // An obtuse triangle: its circumcentre is outside it; the sphere still has that centre.
    // (0,0,0), (4,0,0), (2,1,0): centre (2, -1.5, 0), radius 2.5; (2, -4, 0) is on it.
    CHECK(pr::in_diametral_sphere({0, 0, 0}, {4, 0, 0}, {2, 1, 0}, {2, -4, 0}) == 0);
    CHECK(pr::in_diametral_sphere({0, 0, 0}, {4, 0, 0}, {2, 1, 0}, {2, -4 + 0x1p-51, 0}) == 1);
    CHECK(pr::exact::in_diametral_sphere({0, 0, 0}, {4, 0, 0}, {2, 1, 0}, {2, 0, 1.9}) == 1);
    CHECK(pr::in_diametral_sphere({0, 0, 0}, {4, 0, 0}, {2, 1, 0}, {2, 0, 2}) == 0);
    CHECK(pr::in_diametral_sphere({0, 0, 0}, {4, 0, 0}, {2, 1, 0}, {2, -1.5, 2.5}) == 0);
}

TEST_CASE("the filter agrees with the exact stage on near-degenerate inputs", "[mesh][predicates]") {
    SplitMix rng{20261014};
    int decided_exactly = 0;
    for (int trial = 0; trial < 2000; ++trial) {
        // Nearly cocircular / cospherical: points on a circle or sphere rounded to double, then
        // nudged by a few ulps, at a random offset.
        const double offset = rng.uniform() < 0.5 ? 0.0 : 1e3 * rng.uniform();
        std::array<Point, 5> p{};
        for (Point& x : p) {
            const double t = 6.283185307179586 * rng.uniform();
            const double z = 2.0 * rng.uniform() - 1.0;
            const double s = std::sqrt(1.0 - z * z);
            x = {offset + s * std::cos(t), offset + s * std::sin(t), offset + z};
            for (double& c : x) c += rng.small(3) * u53 * std::max(1.0, std::abs(c));
        }
        const Point flat_d{offset + std::cos(1.0), offset + std::sin(1.0), 0};
        const Point a2{p[0][0], p[0][1], 0}, b2{p[1][0], p[1][1], 0}, c2{p[2][0], p[2][1], 0};
        CHECK(pr::orient2d(a2, b2, c2) == pr::exact::orient2d(a2, b2, c2));
        CHECK(pr::incircle(a2, b2, c2, flat_d) == pr::exact::incircle(a2, b2, c2, flat_d));
        CHECK(pr::orient3d(p[0], p[1], p[2], p[3]) == pr::exact::orient3d(p[0], p[1], p[2], p[3]));
        const int s = pr::insphere(p[0], p[1], p[2], p[3], p[4]);
        CHECK(s == pr::exact::insphere(p[0], p[1], p[2], p[3], p[4]));
        CHECK(pr::in_diametral_sphere(p[0], p[1], p[2], p[3]) ==
              pr::exact::in_diametral_sphere(p[0], p[1], p[2], p[3]));
        decided_exactly += s == 0 ? 1 : 0;
    }
    // Not a check on the filter; only that the samples are near enough to degenerate to matter.
    INFO("exactly cospherical samples: " << decided_exactly);
    SUCCEED();
}

TEST_CASE("incircle_perturbed: never zero, antisymmetric, a valid tie-break", "[mesh][predicates]") {
    // Four cocircular points (a square and a Pythagorean set) with ids.
    const std::vector<std::array<Point, 4>> sets{
        {{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}},
        {{{5, 0, 0}, {3, 4, 0}, {-4, 3, 0}, {0, -5, 0}}},
    };
    const std::array<NodeId, 4> ids{7, 2, 9, 4};
    for (const auto& base : sets) {
        std::array<std::size_t, 4> perm{0, 1, 2, 3};
        std::array<Point, 4> p0{};
        std::array<NodeId, 4> i0{};
        for (std::size_t k = 0; k < 4; ++k) {
            p0[k] = base[k];
            i0[k] = ids[k];
        }
        REQUIRE(pr::incircle(p0[0], p0[1], p0[2], p0[3]) == 0);
        const int reference = pr::incircle_perturbed(p0, i0);
        CHECK(reference != 0);
        do {
            std::array<Point, 4> p{};
            std::array<NodeId, 4> id{};
            for (std::size_t k = 0; k < 4; ++k) {
                p[k] = base[perm[k]];
                id[k] = ids[perm[k]];
            }
            if (pr::orient2d(p[0], p[1], p[2]) == 0) continue;
            CHECK(pr::incircle_perturbed(p, id) == parity(perm) * reference);
        } while (std::next_permutation(perm.begin(), perm.end()));
    }
    // The documented order: the largest id's lifted value is raised most. Raising the query
    // point's lift puts it outside (coefficient -orient(a, b, c)); raising a's puts the query
    // inside (coefficient +orient(b, c, d)).
    const std::array<Point, 4> square{{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}};
    CHECK(pr::incircle_perturbed(square, {0, 1, 2, 3}) == -1);
    CHECK(pr::incircle_perturbed(square, {3, 1, 2, 0}) == 1);
    // Off a tie, the perturbed predicate is the plain one.
    CHECK(pr::incircle_perturbed({{{5, 0, 0}, {3, 4, 0}, {-4, 3, 0}, {0, 0, 0}}}, {1, 2, 3, 4}) == 1);
}

TEST_CASE("insphere_perturbed: never zero, antisymmetric over all 120 orders", "[mesh][predicates]") {
    // Five vertices of a cube (all eight are cospherical), with ids.
    const std::array<Point, 5> base{{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 1}, {1, 1, 1}}};
    const std::array<NodeId, 5> ids{3, 11, 5, 8, 1};
    std::array<std::size_t, 5> perm{0, 1, 2, 3, 4};
    std::optional<int> reference;
    int checked = 0;
    do {
        std::array<Point, 5> p{};
        std::array<NodeId, 5> id{};
        for (std::size_t k = 0; k < 5; ++k) {
            p[k] = base[perm[k]];
            id[k] = ids[perm[k]];
        }
        if (pr::orient3d(p[0], p[1], p[2], p[3]) == 0) continue;
        REQUIRE(pr::insphere(p[0], p[1], p[2], p[3], p[4]) == 0);
        const int s = pr::insphere_perturbed(p, id);
        CHECK(s != 0);
        if (!reference) reference = s * parity(perm);
        CHECK(s == parity(perm) * *reference);
        ++checked;
    } while (std::next_permutation(perm.begin(), perm.end()));
    CHECK(checked > 60);
}
