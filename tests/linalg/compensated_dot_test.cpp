// linalg::compensated_dot (Unit 36, TR6): exact to the rounding of its result where the plain sum
// of the products keeps only rounding.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

#include "NiTCAD/linalg/compensated_dot.hpp"

using NiTCAD::linalg::compensated_dot;

TEST_CASE("compensated dot: a cancelling sum is exact to its rounding (TR6)", "[linalg]") {
    // 1e16 + 1 - 1e16: the plain sum loses the 1 (1e16 + 1 rounds to 1e16).
    const std::vector<double> a{1e16, 1.0, -1e16};
    const std::vector<double> b{1.0, 1.0, 1.0};
    double plain = 0.0;
    for (std::size_t q = 0; q < a.size(); ++q) plain += a[q] * b[q];
    CHECK(plain == 0.0);
    CHECK(compensated_dot(a, b) == 1.0);
    // Products whose rounding matters: (1 + 2^-30)^2 - 1 - 2^-29 = 2^-60, lost in a plain sum.
    const double x = 1.0 + std::ldexp(1.0, -30);
    const std::vector<double> c{x, -1.0, -std::ldexp(1.0, -29)};
    const std::vector<double> d{x, 1.0, 1.0};
    CHECK(compensated_dot(c, d) == std::ldexp(1.0, -60));
    // The term at `skip` is left out.
    CHECK(compensated_dot(a, b, 1) == 0.0);
    CHECK(compensated_dot(std::vector<double>{}, std::vector<double>{}) == 0.0);
}

TEST_CASE("compensated dot: equals the plain sum where nothing cancels", "[linalg]") {
    std::vector<double> a, b;
    for (int k = 1; k <= 100; ++k) {
        a.push_back(1.0 / k);
        b.push_back(static_cast<double>(k));
    }
    CHECK(std::abs(compensated_dot(a, b) / 100.0 - 1.0) <= 1e-15);
}
