// L0 tests of the energy-balance boundary-layer reference (verify/reference/energy_layer.hpp; DECISIONS Unit 55
// 16.18 (a), DT-Q15 A, 16.42, TVB-13):
//   - Durand-Kerner on a quartic with known real and complex roots (1e-28), and on roots spread over ten decades
//     (1e-26 relative);
//   - the zero-field limit: the roots are +-1/L_D and +-sqrt(G / kappa) (1e-9 relative at E0 = 1e-9 V/cm);
//   - the second method (the unreduced determinant) vanishes at the roots of a hot state (1e-25);
//   - refusals.
#include "NiTCAD/verify/reference/energy_layer.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
#include <vector>

namespace ref = NiTCAD::verify::reference;
using ref::DD;

namespace {

double nearest(const std::vector<ref::ComplexDD>& roots, double re, double im) {
    double best = INFINITY;
    for (const auto& z : roots) {
        const DD dr = z.re - DD(re), di = z.im - DD(im);
        best = std::min(best, std::sqrt((dr * dr + di * di).to_double()));
    }
    return best;
}

constexpr double q = 1.602176634e-19, k = 1.380649e-23, eps = 11.7 * 8.8541878128e-14;

ref::LayerInput hot(double E0, double Th) {
    return {.q_C = q,
            .k_J_per_K = k,
            .eps_F_per_cm = eps,
            .density_cm3 = 1e17,
            .field_V_per_cm = E0,
            .mobility_cm2_Vs = 1000.0,
            .mobility_slope = -1.5,
            .temperature_K = Th,
            .lattice_K = 300.0,
            .tau_w_s = 0.35e-12};
}

}  // namespace

TEST_CASE("energy_layer: Durand-Kerner on known roots", "[energy_layer]") {
    // (s - 1)(s + 2)(s^2 - 6 s + 10) = s^4 - 5 s^3 + 2 s^2 + 22 s - 20.
    const auto r = ref::durand_kerner({DD(-20.0), DD(22.0), DD(2.0), DD(-5.0)});
    REQUIRE(r.has_value());
    CHECK(r->roots.size() == 4);
    CHECK(nearest(r->roots, 1.0, 0.0) <= 1e-28);
    CHECK(nearest(r->roots, -2.0, 0.0) <= 1e-28);
    CHECK(nearest(r->roots, 3.0, 1.0) <= 1e-28);
    CHECK(nearest(r->roots, 3.0, -1.0) <= 1e-28);
    CHECK(r->residual <= 1e-30);
    // Roots over ten decades: 1e-3, 2e6, 1e5, -1e7.
    const std::vector<double> z{1e-3, 2e6, 1e5, -1e7};
    std::vector<DD> c{DD(1.0)};  // lowest first, built by (s - z_i) products
    for (const double zi : z) {
        std::vector<DD> next(c.size() + 1, DD(0.0));
        for (std::size_t j = 0; j < c.size(); ++j) {
            next[j + 1] = next[j + 1] + c[j];
            next[j] = next[j] - DD(zi) * c[j];
        }
        c = next;
    }
    c.pop_back();  // monic
    const auto w = ref::durand_kerner(c);
    REQUIRE(w.has_value());
    for (const double zi : z) {
        CAPTURE(zi);
        CHECK(nearest(w->roots, zi, 0.0) <= 1e-26 * std::abs(zi));
    }
}

TEST_CASE("energy_layer: the zero-field limit", "[energy_layer]") {
    ref::LayerInput in = hot(1e-9, 300.0);
    const auto r = ref::energy_layer(in);
    REQUIRE(r.has_value());
    const double L_D = std::sqrt(eps * k * 300.0 / (q * q * 1e17));
    const double kappa = 2.5 * k * k / q * 1e17 * 1000.0 * 300.0, G = 1.5 * k * 1e17 / 0.35e-12;
    const double rates[2] = {1.0 / L_D, std::sqrt(G / kappa)};
    for (const double s0 : rates) {
        for (const double sign : {1.0, -1.0}) {
            double best = INFINITY;
            for (const auto& s : r->roots) best = std::min(best, std::abs(s - std::complex<double>(sign * s0, 0.0)));
            CAPTURE(s0, sign);
            CHECK(best <= 1e-9 * s0);
        }
    }
    CHECK(std::abs(r->length_cm * std::min(rates[0], rates[1]) / 2.0 - 1.0) <= 1e-9);
}

TEST_CASE("energy_layer: the determinant vanishes at the roots of a hot state", "[energy_layer]") {
    for (const auto& [E0, Th] : {std::pair{1e3, 304.0}, std::pair{1e4, 650.0}, std::pair{1e5, 2500.0}}) {
        const auto r = ref::energy_layer(hot(E0, Th));
        REQUIRE(r.has_value());
        CAPTURE(E0, r->length_cm, r->determinant_residual, r->root_residual);
        CHECK(r->determinant_residual <= 1e-25);
        CHECK(r->root_residual <= 1e-28);
        CHECK(r->length_cm > 0.0);
        CHECK(std::isfinite(r->length_cm));
    }
}

TEST_CASE("energy_layer: refusals", "[energy_layer]") {
    CHECK_FALSE(ref::durand_kerner({}).has_value());
    CHECK_FALSE(ref::durand_kerner({DD(0.0), DD(1.0)}).has_value());
    CHECK_FALSE(ref::energy_layer(hot(0.0, 300.0)).has_value());
    CHECK_FALSE(ref::energy_layer(hot(1e3, 200.0)).has_value());  // below the lattice
}
