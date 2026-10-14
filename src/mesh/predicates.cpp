#include "predicates.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <type_traits>

#include "NiTCAD/base/contract.hpp"
#include "exact_integer.hpp"
#include "vec_ops.hpp"

namespace NiTCAD::mesh::predicates {

namespace {

using detail::ExactInteger;

// A double with a forward error bound: |v - exact| <= gamma_d M, where M is the expression evaluated
// exactly on absolute values (m is M computed in double; on positive operands it carries the same
// relative bound). By induction: a sum's error is the operands' errors plus one rounding, so
// d = max + 1; a product's relative errors compound, (1 + gamma_a)(1 + gamma_b)(1 + u) <=
// 1 + gamma_(a+b+1), so d = a + b + 1.
// The model fl(x op y) = (x op y)(1 + delta) fails only when a product underflows (a subnormal sum
// is exact); such a product sets `underflow` and the filter defers to the exact stage.
struct Filtered {
    double v = 0.0;
    double m = 0.0;
    int d = 0;
    bool underflow = false;
};

Filtered operator+(const Filtered& a, const Filtered& b) {
    return {a.v + b.v, a.m + b.m, std::max(a.d, b.d) + 1, a.underflow || b.underflow};
}
Filtered operator-(const Filtered& a, const Filtered& b) {
    return {a.v - b.v, a.m + b.m, std::max(a.d, b.d) + 1, a.underflow || b.underflow};
}
Filtered operator*(const Filtered& a, const Filtered& b) {
    const double m = a.m * b.m;
    const bool tiny = a.m != 0.0 && b.m != 0.0 && m < std::numeric_limits<double>::min();
    return {a.v * b.v, m, a.d + b.d + 1, a.underflow || b.underflow || tiny};
}
Filtered operator-(const Filtered& a) { return {-a.v, a.m, a.d, a.underflow}; }  // exact

template <class T>
T constant(int c);
template <>
Filtered constant<Filtered>(int c) {
    const auto x = static_cast<double>(c);
    return {x, std::abs(x), 0, false};
}
template <>
ExactInteger constant<ExactInteger>(int c) {
    return ExactInteger(c);
}

using detail::cross;
using detail::dot3;
using detail::scale;
using detail::sub;
using detail::Vec;

template <class T>
T det2(const T& a, const T& b, const T& c, const T& d) {
    return a * d - b * c;
}
template <class T>
T det3(const std::array<std::array<T, 3>, 3>& m) {
    return m[0][0] * det2(m[1][1], m[1][2], m[2][1], m[2][2]) -
           m[0][1] * det2(m[1][0], m[1][2], m[2][0], m[2][2]) +
           m[0][2] * det2(m[1][0], m[1][1], m[2][0], m[2][1]);
}
template <class T>
T det4(const std::array<std::array<T, 4>, 4>& m) {
    const auto minor = [&](std::size_t skip) {
        std::array<std::array<T, 3>, 3> r;
        for (std::size_t i = 1; i < 4; ++i) {
            for (std::size_t j = 0, c = 0; j < 4; ++j) {
                if (j != skip) r[i - 1][c++] = m[i][j];
            }
        }
        return det3(r);
    };
    return m[0][0] * minor(0) - m[0][1] * minor(1) + m[0][2] * minor(2) - m[0][3] * minor(3);
}

// The sign of formula(points): the filter first (unless disabled), then exact integers.
template <std::size_t K, class Formula>
int sign_of(const std::array<Point, K>& points, int dimension, bool filter, Formula formula) {
    if (filter) {
        std::array<Vec<Filtered>, K> p;
        for (std::size_t k = 0; k < K; ++k) {
            for (std::size_t c = 0; c < 3; ++c) {
                const double x = c < static_cast<std::size_t>(dimension) ? points[k][c] : 0.0;
                p[k][c] = {x, std::abs(x), 0, false};
            }
        }
        const Filtered r = formula(p);
        if (!r.underflow && std::isfinite(r.v) && std::isfinite(r.m)) {
            // gamma_d M <= 2 d u m for d u << 1 (m within gamma_d of M).
            constexpr double u = 0x1p-53;
            const double bound = 2.0 * static_cast<double>(r.d) * u * r.m;
            if (r.v > bound) return 1;
            if (r.v < -bound) return -1;
        }
    }
    const auto exact = detail::to_integers(points, dimension);
    return formula(exact.c).sign();
}

template <class P>
using Scalar = std::remove_cvref_t<decltype(std::declval<P>()[0][0])>;

const auto orient2d_formula = [](const auto& p) {
    return det2(p[1][0] - p[0][0], p[1][1] - p[0][1], p[2][0] - p[0][0], p[2][1] - p[0][1]);
};

const auto orient3d_formula = [](const auto& p) {
    using T = Scalar<decltype(p)>;
    std::array<std::array<T, 3>, 3> m;
    for (std::size_t r = 0; r < 3; ++r) m[r] = sub(p[r + 1], p[0]);
    return det3(m);
};

// incircle with the lifted column replaced when column >= 0: a unit vector at row `column` (< 3),
// or all -1 (column == 3): the derivative of the lifted determinant with respect to that point's
// lifted value (the rows hold l_k - l_3, k < 3).
auto incircle_formula(int column) {
    return [column](const auto& p) {
        using T = Scalar<decltype(p)>;
        std::array<std::array<T, 3>, 3> m;
        for (std::size_t r = 0; r < 3; ++r) {
            const T dx = p[r][0] - p[3][0];
            const T dy = p[r][1] - p[3][1];
            m[r][0] = dx;
            m[r][1] = dy;
            if (column < 0) {
                m[r][2] = dx * dx + dy * dy;
            } else if (column == 3) {
                m[r][2] = constant<T>(-1);
            } else {
                m[r][2] = constant<T>(static_cast<int>(r) == column ? 1 : 0);
            }
        }
        return det3(m);
    };
}

// insphere = -det of the rows (p_k - p_4, |p_k - p_4|^2), k < 4; replaced columns as for incircle.
auto insphere_formula(int column) {
    return [column](const auto& p) {
        using T = Scalar<decltype(p)>;
        std::array<std::array<T, 4>, 4> m;
        for (std::size_t r = 0; r < 4; ++r) {
            const Vec<T> d = sub(p[r], p[4]);
            m[r][0] = d[0];
            m[r][1] = d[1];
            m[r][2] = d[2];
            if (column < 0) {
                m[r][3] = dot3(d, d);
            } else if (column == 4) {
                m[r][3] = constant<T>(-1);
            } else {
                m[r][3] = constant<T>(static_cast<int>(r) == column ? 1 : 0);
            }
        }
        return -det4(m);
    };
}

const auto diametral_circle_formula = [](const auto& p) {
    // p strictly inside the circle on diameter ab when (a - p) . (b - p) < 0.
    return -((p[0][0] - p[2][0]) * (p[1][0] - p[2][0]) + (p[0][1] - p[2][1]) * (p[1][1] - p[2][1]));
};

const auto diametral_sphere_formula = [](const auto& p) {
    // With u = b - a, v = c - a, n = u x v, the circumcentre is a + N / (2 |n|^2), N =
    // (|u|^2 v - |v|^2 u) x n; p = a + w is strictly inside when |w - N/(2|n|^2)| < |N|/(2|n|^2),
    // that is when w . N - |n|^2 |w|^2 > 0 (degree 6).
    const auto u = sub(p[1], p[0]);
    const auto v = sub(p[2], p[0]);
    const auto w = sub(p[3], p[0]);
    const auto n = cross(u, v);
    const auto g = sub(scale(dot3(u, u), v), scale(dot3(v, v), u));
    const auto big_n = cross(g, n);
    return dot3(w, big_n) - dot3(n, n) * dot3(w, w);
};

// The perturbed sign: the base sign, else the first nonzero lifted-value cofactor in decreasing id
// order (the point with the largest id carries the dominant perturbation).
template <std::size_t K, class Base, class Cofactor>
int perturbed(const std::array<Point, K>& p, const std::array<NodeId, K>& id, int dimension,
              Base base, Cofactor cofactor) {
    for (std::size_t i = 0; i < K; ++i) {
        for (std::size_t j = i + 1; j < K; ++j) NITCAD_EXPECTS(id[i] != id[j]);
    }
    const int s = sign_of(p, dimension, true, base);
    if (s != 0) return s;
    std::array<std::size_t, K> order;
    std::iota(order.begin(), order.end(), std::size_t{0});
    std::ranges::sort(order, [&](std::size_t a, std::size_t b) { return id[a] > id[b]; });
    for (const std::size_t q : order) {
        const int c = sign_of(p, dimension, true, cofactor(static_cast<int>(q)));
        if (c != 0) return c;
    }
    // Unreachable for a non-degenerate base simplex: the query point's cofactor is its
    // orientation, up to sign.
    NITCAD_EXPECTS(false);
    return 0;
}

int orient2d_impl(const Point& a, const Point& b, const Point& c, bool filter) {
    return sign_of(std::array<Point, 3>{a, b, c}, 2, filter, orient2d_formula);
}
int orient3d_impl(const Point& a, const Point& b, const Point& c, const Point& d, bool filter) {
    return sign_of(std::array<Point, 4>{a, b, c, d}, 3, filter, orient3d_formula);
}
int incircle_impl(const Point& a, const Point& b, const Point& c, const Point& d, bool filter) {
    return sign_of(std::array<Point, 4>{a, b, c, d}, 2, filter, incircle_formula(-1));
}
int insphere_impl(const Point& a, const Point& b, const Point& c, const Point& d, const Point& e,
                  bool filter) {
    return sign_of(std::array<Point, 5>{a, b, c, d, e}, 3, filter, insphere_formula(-1));
}
int diametral_sphere_impl(const Point& a, const Point& b, const Point& c, const Point& p,
                          bool filter) {
    return sign_of(std::array<Point, 4>{a, b, c, p}, 3, filter, diametral_sphere_formula);
}

}  // namespace

int orient2d(const Point& a, const Point& b, const Point& c) { return orient2d_impl(a, b, c, true); }
int orient3d(const Point& a, const Point& b, const Point& c, const Point& d) {
    return orient3d_impl(a, b, c, d, true);
}
int incircle(const Point& a, const Point& b, const Point& c, const Point& d) {
    return incircle_impl(a, b, c, d, true);
}
int insphere(const Point& a, const Point& b, const Point& c, const Point& d, const Point& e) {
    return insphere_impl(a, b, c, d, e, true);
}
int in_diametral_circle(const Point& a, const Point& b, const Point& p) {
    return sign_of(std::array<Point, 3>{a, b, p}, 2, true, diametral_circle_formula);
}
int in_diametral_sphere(const Point& a, const Point& b, const Point& c, const Point& p) {
    return diametral_sphere_impl(a, b, c, p, true);
}

int incircle_perturbed(const std::array<Point, 4>& p, const std::array<NodeId, 4>& id) {
    NITCAD_EXPECTS(orient2d(p[0], p[1], p[2]) != 0);
    return perturbed(p, id, 2, incircle_formula(-1), [](int q) { return incircle_formula(q); });
}

int insphere_perturbed(const std::array<Point, 5>& p, const std::array<NodeId, 5>& id) {
    NITCAD_EXPECTS(orient3d(p[0], p[1], p[2], p[3]) != 0);
    return perturbed(p, id, 3, insphere_formula(-1), [](int q) { return insphere_formula(q); });
}

namespace exact {
int orient2d(const Point& a, const Point& b, const Point& c) {
    return orient2d_impl(a, b, c, false);
}
int orient3d(const Point& a, const Point& b, const Point& c, const Point& d) {
    return orient3d_impl(a, b, c, d, false);
}
int incircle(const Point& a, const Point& b, const Point& c, const Point& d) {
    return incircle_impl(a, b, c, d, false);
}
int insphere(const Point& a, const Point& b, const Point& c, const Point& d, const Point& e) {
    return insphere_impl(a, b, c, d, e, false);
}
int in_diametral_sphere(const Point& a, const Point& b, const Point& c, const Point& p) {
    return diametral_sphere_impl(a, b, c, p, false);
}
}  // namespace exact

}  // namespace NiTCAD::mesh::predicates
