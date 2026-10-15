// The predicates' floating-point filter against its soundness argument (predicates.hpp, "Filter
// soundness"; DECISIONS.md U11), and the tie-break against the scheme it documents (M4). The exact
// value of each evaluation is the oracle: no other floating-point computation is trusted.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <span>
#include <vector>

#include "NiTCAD/mesh/mesh.hpp"
#include "exact_integer.hpp"
#include "predicates.hpp"

using NiTCAD::mesh::NodeId;
using NiTCAD::mesh::Point;
namespace pr = NiTCAD::mesh::predicates;
namespace probe = NiTCAD::mesh::predicates::probe;
using NiTCAD::mesh::detail::ExactInteger;
using NiTCAD::mesh::detail::exponent_of;
using NiTCAD::mesh::detail::scaled_integer;
using probe::Kind;

namespace {

constexpr double u53 = 0x1p-53;

struct SplitMix {
    std::uint64_t state;
    std::uint64_t next() {
        std::uint64_t z = (state += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        return z ^ (z >> 31);
    }
    double uniform() { return static_cast<double>(next() >> 11) * 0x1p-53; }
    int small(int range) {
        return static_cast<int>(next() % static_cast<std::uint64_t>(2 * range + 1)) - range;
    }
    std::size_t index(std::size_t n) { return static_cast<std::size_t>(next() % n); }
};

ExactInteger power_of_two(int k) { return ExactInteger::from_shifted(1, false, k); }
ExactInteger abs_of(const ExactInteger& x) { return x.sign() < 0 ? -x : x; }

// |value - exact| <= bound, decided in exact integers at a common power of two.
bool error_within_bound(const probe::Result& r) {
    int e = r.exponent;
    if (r.value != 0.0) e = std::min(e, exponent_of(r.value));
    if (r.bound != 0.0) e = std::min(e, exponent_of(r.bound));
    const ExactInteger exact = r.exact * power_of_two(r.exponent - e);
    const ExactInteger error = abs_of(scaled_integer(r.value, e) - exact);
    return (scaled_integer(r.bound, e) - error).sign() >= 0;
}

// Every property the soundness argument claims for one evaluation.
void check_sound(const probe::Result& r) {
    if (r.underflow) CHECK(r.filter_sign == 0);
    if (r.filter_sign == 0) return;
    CHECK(r.filter_sign == r.exact_sign);
    CHECK(error_within_bound(r));
}

int points_of(Kind k) {
    switch (k) {
        case Kind::orient2d:
        case Kind::diametral_circle: return 3;
        case Kind::insphere:
        case Kind::insphere_cofactor: return 5;
        default: return 4;
    }
}

const std::array<Kind, 8> all_kinds{Kind::orient2d,         Kind::orient3d,
                                    Kind::incircle,         Kind::insphere,
                                    Kind::diametral_circle, Kind::diametral_sphere,
                                    Kind::incircle_cofactor, Kind::insphere_cofactor};

// Points near a common circle or sphere, nudged by a few ulps, at a random offset and scale: their
// polynomials sit near zero, so the filter's decisions fall near its threshold.
std::vector<Point> near_degenerate(SplitMix& rng, int count, double scale, double offset) {
    std::vector<Point> p;
    for (int k = 0; k < count; ++k) {
        const double t = 6.283185307179586 * rng.uniform();
        const double z = 2.0 * rng.uniform() - 1.0;
        const double s = std::sqrt(1.0 - z * z);
        Point x{offset + scale * s * std::cos(t), offset + scale * s * std::sin(t),
                offset + scale * z};
        for (double& c : x) c += rng.small(4) * u53 * std::abs(c);
        p.push_back(x);
    }
    return p;
}

// Points on which the kind's polynomial vanishes exactly in real arithmetic, rounded to doubles
// and nudged by a few ulps, at a random scale and offset: the computed values sit at the noise
// level, where the filter must decide or defer correctly.
std::vector<Point> degenerate_points(Kind k, SplitMix& rng, double scale, double offset) {
    const auto rnd = [&](bool flat) {
        return Point{2 * rng.uniform() - 1, 2 * rng.uniform() - 1, flat ? 0.0 : 2 * rng.uniform() - 1};
    };
    const auto along = [](const Point& a, const Point& b, double t) {
        return Point{a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), a[2] + t * (b[2] - a[2])};
    };
    const auto on_sphere = [&](const Point& c, double r, bool flat) {
        const double t = 6.283185307179586 * rng.uniform();
        const double z = flat ? 0.0 : 2.0 * rng.uniform() - 1.0;
        const double s = std::sqrt(1.0 - z * z);
        return Point{c[0] + r * s * std::cos(t), c[1] + r * s * std::sin(t), c[2] + r * z};
    };
    std::vector<Point> p;
    switch (k) {
        case Kind::orient2d: {
            const Point a = rnd(true), b = rnd(true);
            p = {a, b, along(a, b, 2 * rng.uniform() - 0.5)};
            break;
        }
        case Kind::orient3d: {
            const Point a = rnd(false), b = rnd(false), c = rnd(false);
            p = {a, b, c, along(along(a, b, rng.uniform()), c, rng.uniform())};
            break;
        }
        case Kind::incircle:
            for (int i = 0; i < 4; ++i) p.push_back(on_sphere({0, 0, 0}, 1, true));
            break;
        case Kind::insphere:
            for (int i = 0; i < 5; ++i) p.push_back(on_sphere({0, 0, 0}, 1, false));
            break;
        case Kind::diametral_circle: {
            const Point a = rnd(true), b = rnd(true);
            const Point m = along(a, b, 0.5);
            p = {a, b, on_sphere(m, 0.5 * std::hypot(b[0] - a[0], b[1] - a[1]), true)};
            break;
        }
        case Kind::diametral_sphere: {
            const Point a = on_sphere({0, 0, 0}, 1, true), b = on_sphere({0, 0, 0}, 1, true),
                        c = on_sphere({0, 0, 0}, 1, true);
            p = {a, b, c, on_sphere({0, 0, 0}, 1, false)};  // circumcentre the origin, radius 1
            break;
        }
        case Kind::incircle_cofactor: {
            const Point a = rnd(true), b = rnd(true);
            for (int i = 0; i < 4; ++i) p.push_back(along(a, b, 2 * rng.uniform() - 0.5));
            break;
        }
        case Kind::insphere_cofactor: {
            const Point a = rnd(false), b = rnd(false), c = rnd(false);
            for (int i = 0; i < 5; ++i) p.push_back(along(along(a, b, rng.uniform()), c, rng.uniform()));
            break;
        }
    }
    for (Point& x : p) {
        for (double& c : x) {
            if (c == 0.0) continue;
            c = offset + scale * c;
            // A spread of nudges, from a few ulps to a few thousand, around the threshold 2 d u m.
            c += rng.small(3) * std::ldexp(u53, static_cast<int>(rng.index(13))) * std::abs(c);
        }
    }
    return p;
}

// The Leibniz determinant, independent of the predicates' cofactor expansion.
ExactInteger leibniz(const std::vector<std::vector<ExactInteger>>& m) {
    const std::size_t n = m.size();
    std::vector<std::size_t> perm(n);
    std::iota(perm.begin(), perm.end(), std::size_t{0});
    ExactInteger total;
    do {
        int parity = 1;
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = i + 1; j < n; ++j) {
                if (perm[i] > perm[j]) parity = -parity;
            }
        }
        ExactInteger term(parity);
        for (std::size_t i = 0; i < n; ++i) term = term * m[i][perm[i]];
        total = total + term;
    } while (std::next_permutation(perm.begin(), perm.end()));
    return total;
}

// The documented scheme, evaluated directly: rows (p, |p|^2 + eps^(r+1), 1) with r the point's
// rank in decreasing id and eps = 2^-K, scaled by 2^(K n) to integers. K is far beyond the
// cofactors' size, so the lowest rank with a nonzero cofactor decides. incircle is this
// determinant's sign, insphere its negation (predicates.cpp).
template <std::size_t N>
int lifted_oracle(const std::array<Point, N>& p, const std::array<NodeId, N>& id, int dim) {
    constexpr int K = 128;
    std::array<int, N> rank{};
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) rank[i] += id[j] > id[i] ? 1 : 0;
    }
    const int top = static_cast<int>(N);  // rank + 1 <= N
    std::vector<std::vector<ExactInteger>> m;
    for (std::size_t i = 0; i < N; ++i) {
        std::vector<ExactInteger> row;
        std::int64_t lift = 0;
        for (int c = 0; c < dim; ++c) {
            const auto x = static_cast<std::int64_t>(p[i][static_cast<std::size_t>(c)]);
            row.emplace_back(x);
            lift += x * x;
        }
        row.push_back(ExactInteger(lift) * power_of_two(K * top) + power_of_two(K * (top - rank[i] - 1)));
        row.emplace_back(std::int64_t{1});
        m.push_back(std::move(row));
    }
    const int s = leibniz(m).sign();
    return dim == 2 ? s : -s;
}

}  // namespace

TEST_CASE("the arithmetic model the filter assumes holds", "[mesh][predicates][filter]") {
    CHECK(std::fegetround() == FE_TONEAREST);
    volatile double tiny = std::numeric_limits<double>::min();
    volatile double half = 0.5;
    CHECK(tiny * half > 0.0);  // gradual underflow: no flush to zero
    CHECK(std::numeric_limits<double>::is_iec559);
}

TEST_CASE("filter depths equal the soundness argument's table", "[mesh][predicates][filter]") {
    const std::array<Point, 5> p{{{0.1, 0.2, 0.3}, {1.1, 0.4, 0.7}, {0.3, 1.9, 0.2}, {0.6, 0.5, 1.4}, {0.9, 0.8, 0.6}}};
    const auto depth = [&](Kind k, int column = 0) {
        return probe::evaluate(k, std::span<const Point>(p.data(), static_cast<std::size_t>(points_of(k))), column).depth;
    };
    CHECK(depth(Kind::orient2d) == 4);
    CHECK(depth(Kind::orient3d) == 8);
    CHECK(depth(Kind::incircle) == 11);
    CHECK(depth(Kind::insphere) == 17);
    CHECK(depth(Kind::diametral_circle) == 4);
    CHECK(depth(Kind::diametral_sphere) == 19);
    for (int c = 0; c < 4; ++c) CHECK(depth(Kind::incircle_cofactor, c) == 7);
    for (int c = 0; c < 5; ++c) CHECK(depth(Kind::insphere_cofactor, c) == 12);
}

TEST_CASE("filter: sound near its threshold, against the exact value", "[mesh][predicates][filter]") {
    SplitMix rng{424242};
    for (const Kind k : all_kinds) {
        int decided = 0;
        int deferred = 0;
        int near_threshold = 0;
        for (int trial = 0; trial < 600; ++trial) {
            const double scale = std::ldexp(1.0, rng.small(40));
            const double offset = rng.uniform() < 0.5 ? 0.0 : std::ldexp(1.0, rng.small(20)) * scale;
            const auto pts = degenerate_points(k, rng, scale, offset);
            const int column = k == Kind::incircle_cofactor ? trial % 4
                               : k == Kind::insphere_cofactor ? trial % 5
                                                              : 0;
            const probe::Result r = probe::evaluate(k, pts, column);
            INFO("kind " << static_cast<int>(k) << ", trial " << trial);
            check_sound(r);
            (r.filter_sign != 0 ? decided : deferred) += 1;
            if (r.bound > 0.0 && std::abs(r.value) > 0.25 * r.bound && std::abs(r.value) < 4.0 * r.bound) {
                ++near_threshold;
            }
        }
        INFO("kind " << static_cast<int>(k));
        CHECK(decided > 0);
        CHECK(deferred > 0);
        CHECK(near_threshold > 0);  // the threshold region itself is exercised
    }
}

TEST_CASE("filter: cancellation, overflow and underflow defer or stay sound", "[mesh][predicates][filter]") {
    SECTION("catastrophic cancellation: exact ties at large offsets") {
        // Exactly cocircular integer points far from the origin: v is pure rounding noise.
        for (const double t : {0.0, 1e6, 0x1p40, 0x1p52}) {
            const std::vector<Point> p{{t + 5, t, 0}, {t + 3, t + 4, 0}, {t - 4, t + 3, 0}, {t, t - 5, 0}};
            const auto r = probe::evaluate(Kind::incircle, p);
            CHECK(r.exact_sign == 0);
            CHECK(r.filter_sign == 0);
            check_sound(r);
        }
    }
    SECTION("overflow: the filter defers, the exact stage still decides") {
        const double big = 0x1p+300;
        const std::vector<Point> p{{5 * big, 0, 0}, {3 * big, 4 * big, 0}, {-4 * big, 3 * big, 0}, {0, -4 * big, 0}};
        const auto r = probe::evaluate(Kind::incircle, p);  // degree 4: 2^1200 overflows
        CHECK_FALSE(std::isfinite(r.value + r.magnitude));
        CHECK(r.filter_sign == 0);
        CHECK(r.exact_sign == 1);
    }
    SECTION("a product of magnitudes below the normal range sets underflow, and the filter defers") {
        // a is within 2^-540 of d while b and c are of order 1: a's lift (2^-1080) underflows, the
        // rest of the formula does not, so the bound is a normal number and only the flag defers.
        const double e = 0x1p-540;
        const std::vector<Point> p{{e, e, 0}, {1, 0.25, 0}, {-0.5, 1, 0}, {0, 0, 0}};
        const auto r = probe::evaluate(Kind::incircle, p);
        CHECK(r.underflow);
        CHECK(r.bound >= std::numeric_limits<double>::min());
        CHECK(r.filter_sign == 0);
        CHECK(r.exact_sign == pr::exact::incircle(p[0], p[1], p[2], p[3]));
    }
    SECTION("a product whose value underflows by cancellation while its magnitude does not") {
        // a and d agree to the last bit but one: the differences are 2^-552 by cancellation, their
        // magnitudes 2^-499; the value product (2^-1104) underflows, the magnitude product does not.
        const double s = 0x1p-500;
        const std::vector<Point> p{{s, s, 0}, {s * (1 + 0x1p-52), s * (1 + 0x1p-52), 0}};
        const std::vector<Point> q{p[0], {1, 0, 0}, {0, 1, 0}, p[1]};
        const auto r = probe::evaluate(Kind::incircle, q);
        CHECK(r.underflow);
        CHECK(r.filter_sign == 0);
        check_sound(r);
    }
    SECTION("a product of an exact zero and a magnitude below the normal range sets underflow") {
        // a and d share both coordinates' bits except y: dx = 0 exactly (value 0, magnitude
        // 2^-539), so dx * dx has value 0 from a zero factor and magnitude 2^-1078: only the
        // magnitude check sees it; it keeps every product's m at 0 or at least DBL_MIN.
        const double s = 0x1p-540;
        const std::vector<Point> q{{s, 1.0, 0}, {1, 0.25, 0}, {-0.5, 1, 0}, {s, 0.0, 0}};
        const auto r = probe::evaluate(Kind::incircle, q);
        CHECK(r.underflow);
        CHECK(r.filter_sign == 0);
        check_sound(r);
    }
    SECTION("a bound in the subnormal range stays sound") {
        // m just above DBL_MIN: the bound 2 d u m is subnormal, the filter may still decide.
        SplitMix rng{5};
        int subnormal_bounds = 0;
        for (int trial = 0; trial < 400; ++trial) {
            const double t = std::ldexp(1.0, -511 - static_cast<int>(rng.index(4)));
            const std::vector<Point> p{{t * (1 + rng.uniform()), t * rng.uniform(), 0},
                                       {t * rng.uniform(), t * (1 + rng.uniform()), 0},
                                       {t * rng.uniform(), t * rng.uniform(), 0}};
            const auto r = probe::evaluate(Kind::orient2d, p);
            if (r.bound > 0.0 && r.bound < std::numeric_limits<double>::min()) ++subnormal_bounds;
            check_sound(r);
        }
        CHECK(subnormal_bounds > 0);
    }
    SECTION("subnormal scales: products underflow everywhere; signs stay exact") {
        SplitMix rng{99};
        for (int trial = 0; trial < 200; ++trial) {
            const double scale = std::ldexp(1.0, -520 - static_cast<int>(rng.index(40)));
            auto pts = near_degenerate(rng, 5, scale, 0.0);
            for (const Kind k : {Kind::orient3d, Kind::insphere, Kind::diametral_sphere}) {
                const auto r = probe::evaluate(k, std::span<const Point>(pts.data(), static_cast<std::size_t>(points_of(k))));
                check_sound(r);
            }
        }
    }
}

TEST_CASE("tie-break: equals the documented lifted perturbation, evaluated directly", "[mesh][predicates][filter]") {
    // Integer points on x^2 + y^2 = 625 and on x^2 + y^2 + z^2 = 9: every 4 or 5 of them tie.
    std::vector<Point> circle;
    for (int x = -25; x <= 25; ++x) {
        for (int y = -25; y <= 25; ++y) {
            if (x * x + y * y == 625) circle.push_back({static_cast<double>(x), static_cast<double>(y), 0});
        }
    }
    std::vector<Point> sphere;
    for (int x = -3; x <= 3; ++x) {
        for (int y = -3; y <= 3; ++y) {
            for (int z = -3; z <= 3; ++z) {
                if (x * x + y * y + z * z == 9) {
                    sphere.push_back({static_cast<double>(x), static_cast<double>(y), static_cast<double>(z)});
                }
            }
        }
    }
    REQUIRE(circle.size() == 20);
    REQUIRE(sphere.size() == 30);
    SplitMix rng{2026};
    int ties2 = 0;
    int ties3 = 0;
    for (int trial = 0; trial < 600; ++trial) {
        // Translated by an integer, and half the time one point is moved off the circle.
        const double t = static_cast<double>(rng.small(500));
        std::array<Point, 4> p{};
        std::array<NodeId, 4> id{};
        for (std::size_t k = 0; k < 4; ++k) {
            p[k] = circle[rng.index(circle.size())];
            p[k][0] += t;
            p[k][1] -= t;
            id[k] = static_cast<NodeId>(rng.index(1000));
        }
        if (trial % 2 == 1) p[3][0] += 1;
        bool distinct = true;
        for (std::size_t i = 0; i < 4; ++i) {
            for (std::size_t j = i + 1; j < 4; ++j) distinct = distinct && p[i] != p[j] && id[i] != id[j];
        }
        if (!distinct || pr::orient2d(p[0], p[1], p[2]) == 0) continue;
        ties2 += pr::incircle(p[0], p[1], p[2], p[3]) == 0 ? 1 : 0;
        INFO("2D trial " << trial);
        CHECK(pr::incircle_perturbed(p, id) == lifted_oracle(p, id, 2));
    }
    for (int trial = 0; trial < 300; ++trial) {
        const double t = static_cast<double>(rng.small(500));
        std::array<Point, 5> p{};
        std::array<NodeId, 5> id{};
        for (std::size_t k = 0; k < 5; ++k) {
            p[k] = sphere[rng.index(sphere.size())];
            for (double& c : p[k]) c += t;
            id[k] = static_cast<NodeId>(rng.index(1000));
        }
        if (trial % 2 == 1) p[4][2] += 1;
        bool distinct = true;
        for (std::size_t i = 0; i < 5; ++i) {
            for (std::size_t j = i + 1; j < 5; ++j) distinct = distinct && p[i] != p[j] && id[i] != id[j];
        }
        if (!distinct || pr::orient3d(p[0], p[1], p[2], p[3]) == 0) continue;
        ties3 += pr::insphere(p[0], p[1], p[2], p[3], p[4]) == 0 ? 1 : 0;
        INFO("3D trial " << trial);
        CHECK(pr::insphere_perturbed(p, id) == lifted_oracle(p, id, 3));
    }
    CHECK(ties2 > 100);
    CHECK(ties3 > 50);
}
