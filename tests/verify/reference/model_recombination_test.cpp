// TVR-3 step 5 (DECISIONS Unit 55 15.9, frozen before code): L0 tests of the recombination and trap references of
// TVB-31.11-31.16 (verify/reference/model_recombination.hpp).
//   - Every partial against complex step of the reference's value, with the inputs that carry a partial (E, n1, p1)
//     moved along their given slopes; 15.9's partial rule max(1e-13 |d|, 64 eps S(d)) (MI-Q13).
//   - SRH's second method (15.3): the steady state of the capture-emission balance solved by Brent (TC5) for f, then
//     U = c_n N_t (n (1 - f) - n1 f), against the closed form.
//   - The occupancy's df/d eta = f (1 - f) against complex step; f + (1 - f) = 1.
//   - The band density's quadrature (TC3) against its closed form, to 1e-26 relative plus the quadrature's estimate.
#include "NiTCAD/verify/reference/model_recombination.hpp"

#include "NiTCAD/verify/reference/brent.hpp"
#include "NiTCAD/verify/reference/complex_step.hpp"
#include "NiTCAD/verify/reference/condition.hpp"

#include <catch2/catch_test_macros.hpp>

#include <array>
#include <cmath>
#include <complex>
#include <initializer_list>
#include <limits>
#include <utility>
#include <vector>

namespace ref = NiTCAD::verify::reference;
namespace mr = NiTCAD::verify::reference::model_recombination;
using ref::complex_d;
using ref::DD;

namespace {

constexpr double eps = std::numeric_limits<double>::epsilon();

template <std::size_t N>
std::array<DD, N> dd(const std::array<double, N>& a) {
    std::array<DD, N> r{};
    for (std::size_t i = 0; i < N; ++i) r[i] = DD(a[i]);
    return r;
}
template <std::size_t N>
std::array<complex_d, N> cx(const std::array<double, N>& a) {
    std::array<complex_d, N> r{};
    for (std::size_t i = 0; i < N; ++i) r[i] = a[i];
    return r;
}

// The partial `partial_k` of f (double-double, with its condition sums) against complex step along input i of the
// value form `value` (complex); the inputs listed in `moved` as {target, slope index} follow it along their given
// slopes (dE/dn and the like). The value form is the reference's own value where its complex step is free of
// cancellation, and an algebraically equal form where it is not (see the forms below).
template <std::size_t N, std::size_t K, class F, class V>
void check_partial(const char* what, F&& f, V&& value, const std::array<double, N>& x, std::size_t i,
                   std::size_t partial_k, std::initializer_list<std::pair<std::size_t, std::size_t>> moved) {
    auto g = [&](complex_d t) {
        std::array<complex_d, N> y = cx(x);
        const complex_d d = t - x[i];
        y[i] = t;
        for (const auto& [target, slope] : moved) y[target] = x[target] + x[slope] * d;
        return value(y);
    };
    // h scaled to the input: at the default 1e-200 the imaginary parts of these small rates underflow; the
    // truncation h^2 f''' / 6 is then ~1e-60 relative.
    const auto cs = ref::complex_step(g, x[i], 1e-30 * std::abs(x[i]));
    REQUIRE(cs.status == ref::ComplexStepStatus::ok);
    const DD d = f(dd(x))[partial_k];
    const auto S = ref::condition_sums<K>(f, dd(x));
    const double bound = std::fmax(1e-13 * std::abs(d.to_double()), 64.0 * eps * S[partial_k]);
    INFO(what << " at x[0] " << x[0] << " x[1] " << x[1] << ": analytic " << d.to_double() << " complex step "
              << cs.derivative << " bound " << bound);
    CHECK(std::abs(d.to_double() - cs.derivative) <= bound);
}

// Complex step is forward-mode differentiation of the formula as written, so a quotient N / D whose N and D move in
// proportion (a dominant capture term) cancels in Im. The value forms below take the cancellation out
// algebraically, choosing by real parts (the comparison rule of Martins et al.):
//   SRH: n p - E = n D / tau_n - n tau_p (n + n_ie) / tau_n - n n_ie - E, so R = n / tau_n - [n tau_p (n + n_ie) /
//   tau_n + n n_ie + E] / D when tau_n p dominates D, and the mirror form otherwise.
//   Trap: the smaller of f = occ / D and 1 - f = emp / D (a ratio of positive sums, with no cancellation once it is
//   the smaller); U in the hole capture form c_p ((p + p1) f - p1) when n moves (p, p1 fixed) and the electron form
//   c_n (n - (n + n1) f) when p moves.
complex_d srh_value(const std::array<complex_d, 8>& y) {
    const complex_d &n = y[0], &p = y[1], &E = y[2], &ni = y[5], &tn = y[6], &tp = y[7];
    const complex_d D = tp * (n + ni) + tn * (p + ni);
    if ((tn * p).real() >= (tp * n).real()) return n / tn - (n * tp * (n + ni) / tn + n * ni + E) / D;
    return p / tp - (p * tn * (p + ni) / tp + p * ni + E) / D;
}
complex_d trap_f(const std::array<complex_d, 8>& y) {
    const complex_d occ = y[6] * y[0] + y[7] * y[4];
    const complex_d emp = y[6] * y[2] + y[7] * y[1];
    const complex_d D = occ + emp;
    return occ.real() <= emp.real() ? occ / D : 1.0 - emp / D;
}
complex_d trap_U_hole(const std::array<complex_d, 8>& y) { return y[7] * ((y[1] + y[4]) * trap_f(y) - y[4]); }
complex_d trap_U_electron(const std::array<complex_d, 8>& y) { return y[6] * (y[0] - (y[0] + y[2]) * trap_f(y)); }

std::vector<double> densities() { return {1e-10, 1e-3, 1.0, 1e5, 1e10, 1e15, 1e18, 1e21}; }

}  // namespace

TEST_CASE("model_recombination: Scharfetter lifetime", "[verify][model_recombination]") {
    // tau0 at N = 0 and tau0 / 2 at N = N_ref, exactly; tau0 / (1 + N / N_ref) at 1e19 to the operation's rounding.
    const DD tau0(1e-5), Nref(5e16);
    CHECK(mr::scharfetter_lifetime(tau0, Nref, DD(0.0)) == tau0);
    CHECK(mr::scharfetter_lifetime(tau0, Nref, Nref) == tau0 / DD(2.0));
    CHECK(std::abs((mr::scharfetter_lifetime(tau0, Nref, DD(1e19)) * DD(201.0) - tau0).hi) <= 4.0 * ref::u_dd * 1e-5);
}

TEST_CASE("model_recombination: SRH, Auger, radiative partials against complex step", "[verify][model_recombination]") {
    const double ni = 1e10;
    for (const double n : densities()) {
        for (const double ratio : {1e-6, 0.5, 1.0, 2.0, 1e6}) {
            const double p = ratio * ni * ni / n;
            for (const double slope : {0.0, -0.3}) {
                const double E = ni * ni;
                const std::array<double, 8> x{n, p, E, slope * E / n, slope * E / p, ni, 1e-6, 3e-7};
                auto f = [](const auto& y) { return mr::srh(y); };
                check_partial<8, 3>("SRH dR/dn", f, srh_value, x, 0, 1, {{2, 3}});
                check_partial<8, 3>("SRH dR/dp", f, srh_value, x, 1, 2, {{2, 4}});
                const std::array<double, 7> xa{n, p, E, slope * E / n, slope * E / p, 2.8e-31, 9.9e-32};
                auto fa = [](const auto& y) { return mr::auger(y); };
                auto va = [](const std::array<complex_d, 7>& y) { return mr::auger(y)[0]; };
                check_partial<7, 3>("Auger dR/dn", fa, va, xa, 0, 1, {{2, 3}});
                check_partial<7, 3>("Auger dR/dp", fa, va, xa, 1, 2, {{2, 4}});
                const std::array<double, 6> xr{n, p, E, slope * E / n, slope * E / p, 1e-10};
                auto fr = [](const auto& y) { return mr::radiative(y); };
                auto vr = [](const std::array<complex_d, 6>& y) { return mr::radiative(y)[0]; };
                check_partial<6, 3>("radiative dR/dn", fr, vr, xr, 0, 1, {{2, 3}});
                check_partial<6, 3>("radiative dR/dp", fr, vr, xr, 1, 2, {{2, 4}});
            }
        }
    }
}

TEST_CASE("model_recombination: trap kinetics partials, and SRH against the balance root", "[verify][model_recombination]") {
    const double ni = 1e10;
    for (const double tau : {-20.0, 0.0, 20.0}) {
        const double n1 = ni * std::exp(tau), p1 = ni * std::exp(-tau);
        for (const double n : densities()) {
            for (const double p : {1e-5, 1e5, 1e15}) {
                for (const double g : {0.0, 0.2}) {
                    const std::array<double, 8> x{n, p, n1, g * n1 / n, p1, g * p1 / p, 1e-8, 3e-8};
                    auto f = [](const auto& y) { return mr::trap_kinetics(y); };
                    check_partial<8, 7>("trap df/dn", f, trap_f, x, 0, 2, {{2, 3}});
                    check_partial<8, 7>("trap df/dp", f, trap_f, x, 1, 3, {{4, 5}});
                    check_partial<8, 7>("trap dU/dn", f, trap_U_hole, x, 0, 5, {{2, 3}});
                    check_partial<8, 7>("trap dU/dp", f, trap_U_electron, x, 1, 6, {{4, 5}});
                }
            }
        }
    }
    // The trap rate against its electron capture form U = cn (n (1 - f) - n1 f), and f + (1 - f) = 1.
    for (const double tau : {-20.0, 0.0, 20.0}) {
        const DD n1 = DD(ni) * exp(DD(tau)), p1 = DD(ni) * exp(DD(-tau));
        for (const double n : densities()) {
            for (const double p : {1e-5, 1e5, 1e15}) {
                const auto k = mr::trap_kinetics(
                    std::array<DD, 8>{DD(n), DD(p), n1, DD(0.0), p1, DD(0.0), DD(1e-8), DD(3e-8)});
                const DD U = DD(1e-8) * (DD(n) * k[1] - n1 * k[0]);
                INFO("tau " << tau << " n " << n << " p " << p);
                // U as a difference of terms of size cn (n + n1).
                CHECK(std::abs((U - k[4]).hi) <= 32.0 * ref::u_dd * (DD(1e-8) * (DD(n) + n1)).to_double());
                CHECK(std::abs(((k[0] + k[1]) - DD(1.0)).hi) <= 4.0 * ref::u_dd);
            }
        }
    }
    // SRH's second method: the steady state of the capture-emission balance (mid-gap level, N_t = 1 cm^-3).
    for (const double n : densities()) {
        for (const double ratio : {1e-6, 0.5, 2.0, 1e6}) {
            const double p = ratio * ni * ni / n;
            const DD cn = DD(1.0) / DD(1e-6), cp = DD(1.0) / DD(3e-7), N(n), P(p), n1(ni), p1(ni);
            auto balance = [&](const DD& f) {
                return cn * (N * (DD(1.0) - f) - n1 * f) - cp * (P * f - p1 * (DD(1.0) - f));
            };
            const auto root = ref::brent_zero(balance, DD(0.0), DD(1.0), 1e-32);
            REQUIRE(root.status == ref::BrentStatus::ok);
            const DD U = cn * (N * (DD(1.0) - root.root) - n1 * root.root);
            const DD R = mr::srh(std::array<DD, 8>{N, P, DD(ni) * DD(ni), DD(0.0), DD(0.0), DD(ni), DD(1e-6), DD(3e-7)})[0];
            INFO("n " << n << " ratio " << ratio << ": U " << U.to_double() << " R " << R.to_double());
            // Brent's f is good to 6 u_dd |f| + 2t (its bound (2.18)); U's sensitivity to f is c_n (n + n1); and U is
            // a difference of terms of size c_n n.
            const double tol = (6.0 * ref::u_dd + 2e-32) * (cn * (N + n1)).to_double() +
                               16.0 * ref::u_dd * (cn * (N + n1)).to_double();
            CHECK(std::abs((U - R).hi) <= tol);
        }
    }
}

TEST_CASE("model_recombination: occupancy and band density", "[verify][model_recombination]") {
    // |x| <= 600: beyond it h f (1 - f) falls below the double range and the complex step returns 0.
    for (double x = -600.0; x <= 600.0; x += 7.3) {
        const auto cs = ref::complex_step([](complex_d t) { return mr::fermi_occupancy(t)[0]; }, x,
                                          1e-30 * std::fmax(1.0, std::abs(x)));
        REQUIRE(cs.status == ref::ComplexStepStatus::ok);
        const auto o = mr::fermi_occupancy(DD(x));
        const double d = o[2].to_double();
        CHECK(std::abs(d + cs.derivative) <= std::fmax(1e-13 * d, 64.0 * eps * (1.0 + std::abs(x)) * d));
        CHECK(std::abs(((o[0] + o[1]) - DD(1.0)).hi) <= 4.0 * ref::u_dd);
    }
    for (const std::array<double, 3> band : {std::array<double, 3>{-0.54, 0.54, 0.0},
                                            std::array<double, 3>{-0.2, 0.3, 0.25},
                                            std::array<double, 3>{0.1, 0.101, -0.4}}) {
        for (const double T : {200.0, 300.0, 500.0}) {
            const std::array<DD, 5> x{DD(1e11), DD(band[0]), DD(band[1]), DD(band[2]), DD(T)};
            for (const int k : {0, 1}) {
                const auto q = mr::band_density_quadrature(x, k);
                REQUIRE(q.ok);
                const DD c = mr::band_density_closed(x, k);
                INFO("band " << band[0] << ".." << band[1] << " E_F " << band[2] << " T " << T << " k " << k);
                CHECK(std::abs((q.value - c).hi) <= 1e-26 * std::abs(c.hi) + q.abs_error);
            }
        }
    }
}
