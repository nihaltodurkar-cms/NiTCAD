#include "dual_geometry.hpp"

#include <cmath>
#include <cstddef>

#include "NiTCAD/base/contract.hpp"
#include "exact_integer.hpp"
#include "vec_ops.hpp"

namespace NiTCAD::mesh::detail {

namespace {

using V = Vec<ExactInteger>;

}  // namespace

double distance(const Point& a, const Point& b) noexcept {
    return std::hypot(b[0] - a[0], b[1] - a[1], b[2] - a[2]);
}

std::array<double, 3> triangle_couplings(const std::array<Point, 3>& p) {
    const auto I = to_integers(p, 2);
    const V u = sub(I.c[1], I.c[0]);
    const V v = sub(I.c[2], I.c[0]);
    const ExactInteger twice_area = u[0] * v[1] - u[1] * v[0];  // degree 2
    NITCAD_EXPECTS(twice_area.sign() > 0);
    const double cross = twice_area.to_double(2 * I.e0);
    std::array<double, 3> piece{};
    for (std::size_t e = 0; e < 3; ++e) {
        const auto i = static_cast<std::size_t>(triangle_edges[e][0]);
        const auto j = static_cast<std::size_t>(triangle_edges[e][1]);
        const std::size_t k = 3 - i - j;
        const ExactInteger d = dot3(sub(I.c[i], I.c[k]), sub(I.c[j], I.c[k]));  // degree 2
        piece[e] = distance(p[i], p[j]) * d.to_double(2 * I.e0) / (2.0 * cross);
    }
    return piece;
}

std::array<double, 6> tetrahedron_couplings(const std::array<Point, 4>& p) {
    const auto I = to_integers(p, 3);
    const int e0 = I.e0;
    const V u = sub(I.c[1], I.c[0]);
    const V v = sub(I.c[2], I.c[0]);
    const V w = sub(I.c[3], I.c[0]);
    const ExactInteger det = dot3(u, cross(v, w));  // degree 3, six times the volume
    NITCAD_EXPECTS(det.sign() > 0);
    // The circumcentre is p0 + num / (2 det), num of degree 4.
    const V num = add(add(scale(dot3(u, u), cross(v, w)), scale(dot3(v, v), cross(w, u))),
                      scale(dot3(w, w), cross(u, v)));
    const double two_det = 2.0 * det.to_double(3 * e0);

    // Per face (opposite vertex l): its normal's length and the signed distance from it to the
    // circumcentre, positive toward l.
    std::array<double, 4> normal_length{};
    std::array<double, 4> centre_distance{};
    for (std::size_t l = 0; l < 4; ++l) {
        std::array<std::size_t, 3> f{};
        for (std::size_t q = 0, c = 0; q < 4; ++q) {
            if (q != l) f[c++] = q;
        }
        const V n = cross(sub(I.c[f[1]], I.c[f[0]]), sub(I.c[f[2]], I.c[f[0]]));  // degree 2
        const int side = dot3(sub(I.c[l], I.c[f[0]]), n).sign();
        NITCAD_EXPECTS(side != 0);
        // (centre - p_f0) . n = (num - 2 det (p_f0 - p0)) . n / (2 det), numerator of degree 6.
        const ExactInteger two = det + det;
        const ExactInteger along = dot3(sub(num, scale(two, sub(I.c[f[0]], I.c[0]))), n);
        normal_length[l] = std::sqrt(dot3(n, n).to_double(4 * e0));
        centre_distance[l] =
            static_cast<double>(side) * (along.to_double(6 * e0) / two_det) / normal_length[l];
    }

    std::array<double, 6> piece{};
    for (std::size_t e = 0; e < 6; ++e) {
        const auto i = static_cast<std::size_t>(tetrahedron_edges[e][0]);
        const auto j = static_cast<std::size_t>(tetrahedron_edges[e][1]);
        std::array<std::size_t, 2> other{};
        for (std::size_t q = 0, c = 0; q < 4; ++q) {
            if (q != i && q != j) other[c++] = q;
        }
        const double length = distance(p[i], p[j]);
        double sum = 0.0;
        for (std::size_t t = 0; t < 2; ++t) {
            const std::size_t k = other[t];      // the face ijk ...
            const std::size_t l = other[1 - t];  // ... is the one opposite l
            const ExactInteger d = dot3(sub(I.c[i], I.c[k]), sub(I.c[j], I.c[k]));
            const double h = length * d.to_double(2 * e0) / (2.0 * normal_length[l]);
            sum += h * centre_distance[l];
        }
        piece[e] = 0.5 * sum;
    }
    return piece;
}

std::array<double, 3> facet_pieces(const std::array<Point, 3>& p) {
    const auto I = to_integers(p, 3);
    const V n = cross(sub(I.c[1], I.c[0]), sub(I.c[2], I.c[0]));
    const ExactInteger nn = dot3(n, n);
    NITCAD_EXPECTS(nn.sign() > 0);
    const double normal_length = std::sqrt(nn.to_double(4 * I.e0));
    std::array<double, 3> piece{};
    for (std::size_t i = 0; i < 3; ++i) {
        const std::size_t j = (i + 1) % 3;
        const std::size_t k = (i + 2) % 3;
        const V ij = sub(I.c[i], I.c[j]);
        const V ik = sub(I.c[i], I.c[k]);
        // (L_ij^2 cot(theta_k) + L_ik^2 cot(theta_j)) / 8 with cot = dot / |n|; degree 4.
        const ExactInteger dot_k = dot3(ik, sub(I.c[j], I.c[k]));
        const ExactInteger dot_j = dot3(ij, sub(I.c[k], I.c[j]));
        const ExactInteger numerator = dot3(ij, ij) * dot_k + dot3(ik, ik) * dot_j;
        piece[i] = numerator.to_double(4 * I.e0) / (8.0 * normal_length);
    }
    return piece;
}

Point triangle_circumcentre(const std::array<Point, 3>& p) {
    const auto I = to_integers(p, 3);
    const V u = sub(I.c[1], I.c[0]);
    const V v = sub(I.c[2], I.c[0]);
    const V n = cross(u, v);
    const ExactInteger nn = dot3(n, n);
    NITCAD_EXPECTS(nn.sign() > 0);
    const V big_n = cross(sub(scale(dot3(u, u), v), scale(dot3(v, v), u)), n);  // degree 5
    const double den = 2.0 * nn.to_double(4 * I.e0);
    return {p[0][0] + big_n[0].to_double(5 * I.e0) / den,
            p[0][1] + big_n[1].to_double(5 * I.e0) / den,
            p[0][2] + big_n[2].to_double(5 * I.e0) / den};
}

Point tetrahedron_circumcentre(const std::array<Point, 4>& p) {
    const auto I = to_integers(p, 3);
    const V u = sub(I.c[1], I.c[0]);
    const V v = sub(I.c[2], I.c[0]);
    const V w = sub(I.c[3], I.c[0]);
    const ExactInteger det = dot3(u, cross(v, w));
    NITCAD_EXPECTS(det.sign() != 0);
    const V num = add(add(scale(dot3(u, u), cross(v, w)), scale(dot3(v, v), cross(w, u))),
                      scale(dot3(w, w), cross(u, v)));  // degree 4; the centre is p0 + num / (2 det)
    const double two_det = 2.0 * det.to_double(3 * I.e0);
    return {p[0][0] + num[0].to_double(4 * I.e0) / two_det,
            p[0][1] + num[1].to_double(4 * I.e0) / two_det,
            p[0][2] + num[2].to_double(4 * I.e0) / two_det};
}

}  // namespace NiTCAD::mesh::detail
