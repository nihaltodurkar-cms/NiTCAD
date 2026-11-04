// TC9 (DECISIONS Unit 55 section 11.3; readings ruled in 11.13): L0 tests of the cyclic Jacobi eigensolver.
//   - known spectra judged absolutely, |lambda^ - lambda| <= 8 n eps ||A||_2 (eps = 2^-52): diagonal matrices; the
//     integer-scaled Hilbert matrices L H_n, n = 3 .. 10, whose lambda_min and lambda_max are the extreme roots of
//     OEIS A076823's exact characteristic polynomial, found by TC5's Brent in double-double; rank-one updates
//     D + rho z z^H (roots of the secular equation, Brent in double-double) and I + z z^H; a complex unitary
//     similarity M D M^H of a known integer diagonal; Forsythe and Henrici's example (42);
//   - every eigenpair's residual ||A v - lambda v||_2 <= 8 n eps ||A||_2 and ||V^H V - I||_F <= 8 n eps, both
//     evaluated in double-double from the double results;
//   - refusals.
// Every input matrix is exact in double (integer or Gaussian-integer entries), so the references are those of
// the matrix the solver sees.
#include "NiTCAD/verify/reference/jacobi_eigen.hpp"

#include "NiTCAD/verify/reference/brent.hpp"
#include "NiTCAD/verify/reference/double_double.hpp"

#include "hilbert_charpoly_oeis.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstring>
#include <numeric>
#include <vector>

namespace ref = NiTCAD::verify::reference;
using ref::DD;
using cd = std::complex<double>;

namespace {

constexpr double eps = 0x1p-52;

DD from_decimal(const char* s) {
    const bool neg = *s == '-';
    if (neg) ++s;
    DD v(0.0);
    for (; *s; ++s) v = v * DD(10.0) + DD(static_cast<double>(*s - '0'));
    return neg ? -v : v;
}

struct DC {
    DD re, im;
};
DC mul(cd a, cd b) {
    return {DD(a.real()) * DD(b.real()) - DD(a.imag()) * DD(b.imag()),
            DD(a.real()) * DD(b.imag()) + DD(a.imag()) * DD(b.real())};
}

// max_j ||A v_j - lambda_j v_j||_2 and ||V^H V - I||_F, in double-double.
struct Quality {
    double residual = 0.0, orthogonality = 0.0;
};
Quality quality(const std::vector<cd>& a, std::size_t n, const ref::JacobiResult& r) {
    Quality q;
    for (std::size_t j = 0; j < n; ++j) {
        DD s2(0.0);
        for (std::size_t i = 0; i < n; ++i) {
            DD re(0.0), im(0.0);
            for (std::size_t k = 0; k < n; ++k) {
                const DC p = mul(a[i * n + k], r.eigenvectors[k * n + j]);
                re += p.re;
                im += p.im;
            }
            const DC l = mul(cd(r.eigenvalues[j], 0.0), r.eigenvectors[i * n + j]);
            re -= l.re;
            im -= l.im;
            s2 += re * re + im * im;
        }
        q.residual = std::max(q.residual, std::sqrt(s2.to_double()));
    }
    DD g2(0.0);
    for (std::size_t i = 0; i < n; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            DD re(i == j ? -1.0 : 0.0), im(0.0);
            for (std::size_t k = 0; k < n; ++k) {
                const DC p = mul(std::conj(r.eigenvectors[k * n + i]), r.eigenvectors[k * n + j]);
                re += p.re;
                im += p.im;
            }
            g2 += re * re + im * im;
        }
    q.orthogonality = std::sqrt(g2.to_double());
    return q;
}

// Judges one solve against a full reference spectrum (ascending); returns the worst |lambda^ - lambda| / (n eps ||A||_2).
void check_against(const char* name, const std::vector<cd>& a, std::size_t n, const std::vector<DD>& lambda) {
    const auto r = ref::jacobi_eigen(a, n);
    REQUIRE(r.status == ref::JacobiStatus::ok);
    double norm2 = 0.0;
    for (const DD& l : lambda) norm2 = std::max(norm2, std::abs(l.to_double()));
    double worst = 0.0;
    for (std::size_t j = 0; j < n; ++j)
        worst = std::max(worst, std::abs((DD(r.eigenvalues[j]) - lambda[j]).to_double()) / (n * eps * norm2));
    const Quality q = quality(a, n, r);
    INFO(name << " (n = " << n << ", " << r.sweeps << " sweeps): eigenvalues " << worst << " n eps ||A||; residual "
              << q.residual / (n * eps * norm2) << " n eps ||A||; orthogonality " << q.orthogonality / (n * eps)
              << " n eps; off " << r.off_norm);
    CHECK(worst <= 8.0);
    CHECK(q.residual <= 8.0 * n * eps * norm2);
    CHECK(q.orthogonality <= 8.0 * n * eps);
}

long long lcm_to(int m) {
    long long l = 1;
    for (int k = 2; k <= m; ++k) l = l / std::gcd(l, static_cast<long long>(k)) * k;
    return l;
}

// The extreme roots of P_n(x) = sum_k T(n, k) x^k by Brent in double-double.
struct Extremes {
    DD lmin, lmax;
    bool ok = false;
};
Extremes hilbert_extremes(int n) {
    std::vector<DD> c;
    for (int k = 0; k <= n; ++k) c.push_back(from_decimal(tc9_oeis::hilbert_charpoly[n][k]));
    const auto p = [c](const DD& x) {
        DD s(0.0);
        for (std::size_t k = c.size(); k-- > 0;) s = s * x + c[k];
        return s;
    };
    const auto sgn = [](const DD& v) { return v.hi > 0.0 ? 1 : (v.hi < 0.0 ? -1 : 0); };
    Extremes e;
    // lambda_min: the first sign change upward from 1e-20 on a log grid (200 points per decade).
    DD xa(1e-20);
    int sa = sgn(p(xa));
    bool found_min = false;
    for (int j = 1; j <= 4400 && !found_min; ++j) {
        const DD xb(std::pow(10.0, -20.0 + j / 200.0));
        const int sb = sgn(p(xb));
        if (sb != sa) {
            const auto res = ref::brent_zero(p, xa, xb, 1e-45);
            e.lmin = res.root;
            found_min = res.status == ref::BrentStatus::ok;
        }
        xa = xb;
        sa = sb;
    }
    // lambda_max: the first sign change downward from pi (lambda_max(H_n) < pi) in steps of 1e-4.
    xa = DD(3.2);
    sa = sgn(p(xa));
    bool found_max = false;
    for (int j = 1; j <= 32000 && !found_max; ++j) {
        const DD xb(3.2 - j * 1e-4);
        const int sb = sgn(p(xb));
        if (sb != sa) {
            const auto res = ref::brent_zero(p, xb, xa, 1e-45);
            e.lmax = res.root;
            found_max = res.status == ref::BrentStatus::ok;
        }
        xa = xb;
        sa = sb;
    }
    e.ok = found_min && found_max;
    return e;
}

}  // namespace

TEST_CASE("TC9: the OEIS table's structure (T(n, 0) = (-1)^n; n = 2 is 1 - 16x + 12x^2)") {
    for (int n = 1; n <= 10; ++n) CHECK(std::strcmp(tc9_oeis::hilbert_charpoly[n][0], n % 2 ? "-1" : "1") == 0);
    CHECK(std::strcmp(tc9_oeis::hilbert_charpoly[2][1], "-16") == 0);
    CHECK(std::strcmp(tc9_oeis::hilbert_charpoly[2][2], "12") == 0);
}

TEST_CASE("TC9: diagonal matrices are returned exactly, sorted, with no sweep") {
    const std::vector<double> d = {3.0, -1.0, 2.0, 0.5, -7.25};
    const std::size_t n = d.size();
    std::vector<cd> a(n * n, cd(0.0, 0.0));
    for (std::size_t i = 0; i < n; ++i) a[i * n + i] = cd(d[i], 0.0);
    const auto r = ref::jacobi_eigen(a, n);
    REQUIRE(r.status == ref::JacobiStatus::ok);
    CHECK(r.sweeps == 0);
    std::vector<double> s = d;
    std::sort(s.begin(), s.end());
    CHECK(r.eigenvalues == s);
    const Quality q = quality(a, n, r);
    CHECK(q.residual == 0.0);
    CHECK(q.orthogonality == 0.0);
}

TEST_CASE("TC9: Hilbert matrices L H_n, n = 3 .. 10: lambda_min and lambda_max within 8 n eps ||A||_2") {
    for (int n = 3; n <= 10; ++n) {
        const long long L = lcm_to(2 * n - 1);
        std::vector<cd> a(static_cast<std::size_t>(n * n));
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j) a[i * n + j] = cd(static_cast<double>(L / (i + j + 1)), 0.0);
        const Extremes e = hilbert_extremes(n);
        REQUIRE(e.ok);
        const DD lmin = DD(static_cast<double>(L)) * e.lmin, lmax = DD(static_cast<double>(L)) * e.lmax;
        const auto r = ref::jacobi_eigen(a, static_cast<std::size_t>(n));
        REQUIRE(r.status == ref::JacobiStatus::ok);
        const double tol = 8.0 * n * eps * lmax.to_double();
        const double emin = std::abs((DD(r.eigenvalues.front()) - lmin).to_double());
        const double emax = std::abs((DD(r.eigenvalues.back()) - lmax).to_double());
        const Quality q = quality(a, static_cast<std::size_t>(n), r);
        INFO("H_" << n << " (L = " << L << ", " << r.sweeps << " sweeps): lambda_min(H_n) = " << e.lmin.to_double()
                  << ", error " << emin / (n * eps * lmax.to_double()) << " n eps ||A|| (relative "
                  << emin / lmin.to_double() << ", printed only); lambda_max error "
                  << emax / (n * eps * lmax.to_double()) << " n eps ||A||; residual "
                  << q.residual / (n * eps * lmax.to_double()) << " n eps ||A||; orthogonality "
                  << q.orthogonality / (n * eps) << " n eps");
        CHECK(emin <= tol);
        CHECK(emax <= tol);
        CHECK(q.residual <= tol);
        CHECK(q.orthogonality <= 8.0 * n * eps);
    }
}

TEST_CASE("TC9: rank-one updates D + rho z z^H against the secular equation (Brent in double-double)") {
    for (std::size_t n : {std::size_t{8}, std::size_t{16}}) {
        for (double rho : {1.0, -1.0}) {
            std::vector<double> d(n);
            std::vector<cd> z(n);
            for (std::size_t k = 0; k < n; ++k) {
                d[k] = static_cast<double>(2 * k + 1);
                z[k] = cd(static_cast<double>(static_cast<int>(k % 3) - 1), static_cast<double>(k % 2 + 1));
            }
            std::vector<cd> a(n * n);
            for (std::size_t i = 0; i < n; ++i)
                for (std::size_t j = 0; j < n; ++j)
                    a[i * n + j] = rho * z[i] * std::conj(z[j]) + (i == j ? cd(d[i], 0.0) : cd(0.0, 0.0));
            double zz = 0.0;
            for (const cd& x : z) zz += std::norm(x);
            const auto f = [&](const DD& lam) {
                DD s(1.0);
                for (std::size_t k = 0; k < n; ++k) s += DD(rho * std::norm(z[k])) / (DD(d[k]) - lam);
                return s;
            };
            std::vector<DD> lambda;
            const DD delta(1e-25);
            if (rho < 0.0) lambda.push_back(ref::brent_zero(f, DD(d[0] - zz - 1.0), DD(d[0]) - delta, 1e-40).root);
            for (std::size_t k = 0; k + 1 < n; ++k)
                lambda.push_back(ref::brent_zero(f, DD(d[k]) + delta, DD(d[k + 1]) - delta, 1e-40).root);
            if (rho > 0.0) lambda.push_back(ref::brent_zero(f, DD(d[n - 1]) + delta, DD(d[n - 1] + zz + 1.0), 1e-40).root);
            REQUIRE(lambda.size() == n);
            check_against(rho > 0.0 ? "D + z z^H" : "D - z z^H", a, n, lambda);
        }
    }
}

TEST_CASE("TC9: I + z z^H (eigenvalue 1 of multiplicity n - 1, and 1 + ||z||^2)") {
    const std::size_t n = 7;
    std::vector<cd> z(n);
    for (std::size_t k = 0; k < n; ++k) z[k] = cd(static_cast<double>(k % 2) + 1.0, -static_cast<double>(k % 3));
    std::vector<cd> a(n * n);
    double zz = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        zz += std::norm(z[i]);
        for (std::size_t j = 0; j < n; ++j) a[i * n + j] = z[i] * std::conj(z[j]) + (i == j ? cd(1.0, 0.0) : cd(0.0, 0.0));
    }
    std::vector<DD> lambda(n - 1, DD(1.0));
    lambda.push_back(DD(1.0 + zz));
    check_against("I + z z^H", a, n, lambda);
}

TEST_CASE("TC9: a complex unitary similarity M D M^H of a known diagonal (n = 6 and 16)") {
    for (std::size_t n : {std::size_t{6}, std::size_t{16}}) {
        std::vector<cd> w(n);
        std::vector<double> d(n);
        double ww = 0.0;
        for (std::size_t k = 0; k < n; ++k) {
            w[k] = cd(static_cast<double>(static_cast<int>(k % 4) - 1), static_cast<double>(static_cast<int>(k % 3) - 1));
            if (k == 0) w[k] = cd(2.0, 1.0);
            ww += std::norm(w[k]);
            d[k] = static_cast<double>(static_cast<int>(3 * k) - 7);
        }
        // M = (w^H w) I - 2 w w^H: M / (w^H w) is a unitary Householder reflector, so M D M^H has eigenvalues
        // (w^H w)^2 d_k; all entries are Gaussian integers, exact in double.
        std::vector<cd> m(n * n), a(n * n, cd(0.0, 0.0));
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j) m[i * n + j] = (i == j ? cd(ww, 0.0) : cd(0.0, 0.0)) - 2.0 * w[i] * std::conj(w[j]);
        for (std::size_t i = 0; i < n; ++i)
            for (std::size_t j = 0; j < n; ++j)
                for (std::size_t k = 0; k < n; ++k) a[i * n + j] += m[i * n + k] * d[k] * std::conj(m[j * n + k]);
        double big = 0.0;
        for (const cd& x : a) big = std::max({big, std::abs(x.real()), std::abs(x.imag())});
        REQUIRE(big < 0x1p53);
        std::vector<DD> lambda;
        for (double dk : d) lambda.push_back(DD(ww * ww) * DD(dk));
        std::sort(lambda.begin(), lambda.end(), [](const DD& x, const DD& y) { return x < y; });
        check_against("M D M^H", a, n, lambda);
    }
}

TEST_CASE("TC9: Forsythe and Henrici's (42) converges with |phi| <= pi/4 (it cycles with phi = pi/2)") {
    const std::vector<cd> a = {2.0, 0.0, 1.0, 0.0, 3.0, 0.0, 1.0, 0.0, 4.0};
    const DD r2 = ref::sqrt(DD(2.0));
    check_against("FH (42)", a, 3, {DD(3.0) - r2, DD(3.0), DD(3.0) + r2});
}

TEST_CASE("TC9: invalid input is refused") {
    CHECK(ref::jacobi_eigen({}, 0).status == ref::JacobiStatus::invalid_input);
    CHECK(ref::jacobi_eigen(std::vector<cd>(17 * 17, cd(0.0, 0.0)), 17).status == ref::JacobiStatus::invalid_input);
    CHECK(ref::jacobi_eigen(std::vector<cd>(5, cd(0.0, 0.0)), 2).status == ref::JacobiStatus::invalid_input);
    CHECK(ref::jacobi_eigen({1.0, cd(0.0, 1.0), cd(0.0, 1.0), 2.0}, 2).status == ref::JacobiStatus::invalid_input);
    CHECK(ref::jacobi_eigen({cd(1.0, 1e-300), 0.0, 0.0, 2.0}, 2).status == ref::JacobiStatus::invalid_input);
    CHECK(ref::jacobi_eigen({cd(NAN, 0.0), 0.0, 0.0, 2.0}, 2).status == ref::JacobiStatus::invalid_input);
    CHECK(ref::jacobi_eigen({1.0, cd(0.0, 1.0), cd(0.0, -1.0), 2.0}, 2).status == ref::JacobiStatus::ok);
    CHECK(ref::jacobi_eigen(std::vector<cd>(16 * 16, cd(1.0, 0.0)), 16).status == ref::JacobiStatus::ok);
}
