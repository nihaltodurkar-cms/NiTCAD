// L0 tests of the thermal references (verify/reference/thermal_slab.hpp; DECISIONS Unit 55 16.3, 16.43; TVB-27,
// TVB-28):
//   - the Kirchhoff transform: u(300) = 0, the round trip (1e-29 K), du/dT = kappa (a central difference, 1e-20);
//   - the slab: its ends exactly, u linear in x;
//   - the self-heated resistor with constant kappa and mobility against its closed form
//     T = T_a + (Q / b) [L (1 - e^{-b x / kappa}) / (1 - e^{-b L / kappa}) - x], b = (3/2)(k/q) J (1e-12 K);
//     its steps-halving difference; the finite-difference second method runs in TVB-27 (CV-Q19's Debug budget);
//   - the Seebeck rise against its Gauss-Kronrod second method (1e-28 eV) and the sign; the Peltier coefficient;
//   - refusals.
#include "NiTCAD/verify/reference/thermal_slab.hpp"

#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

namespace ref = NiTCAD::verify::reference;
using ref::DD;

namespace {

constexpr double q = 1.602176634e-19, k = 1.380649e-23 / 1.602176634e-19;
const ref::Conductivity si{1.48, 1.33};

}  // namespace

TEST_CASE("thermal_slab: the Kirchhoff transform", "[thermal_slab]") {
    CHECK(ref::kirchhoff(si, DD(300.0)).hi == 0.0);
    for (const double T : {250.0, 300.0, 350.0, 400.0, 600.0}) {
        const DD u = ref::kirchhoff(si, DD(T));
        CHECK(std::abs((ref::kirchhoff_inverse(si, u) - DD(T)).to_double()) <= 1e-29 * T);
        const DD d(1e-6);
        const DD slope = (ref::kirchhoff(si, DD(T) + d) - ref::kirchhoff(si, DD(T) - d)) / (DD(2.0) * d);
        const double kappa = 1.48 * std::pow(T / 300.0, -1.33);
        CAPTURE(T);
        CHECK(std::abs(slope.to_double() / kappa - 1.0) <= 1e-12);
    }
}

TEST_CASE("thermal_slab: the slab", "[thermal_slab]") {
    const auto T = ref::slab_temperature(si, 300.0, 400.0, {0.0, 0.25, 0.5, 1.0});
    REQUIRE(T.has_value());
    CHECK((*T)[0].to_double() == 300.0);
    CHECK(std::abs((*T)[3].to_double() - 400.0) <= 1e-12);
    const DD u0 = ref::kirchhoff(si, DD(300.0)), u1 = ref::kirchhoff(si, DD(400.0));
    CHECK(std::abs((ref::kirchhoff(si, (*T)[1]) - (u0 + DD(0.25) * (u1 - u0))).to_double()) <= 1e-28);
    // kappa falls with T, so the gradient is steeper at the hot end and the profile sags: T(1/2) < 350.
    CHECK((*T)[2].to_double() < 350.0);
}

TEST_CASE("thermal_slab: the self-heated resistor's constant-coefficient limit", "[thermal_slab][.tvr4_reference_release]") {
    ref::ResistorInput in{.q_C = q,
                          .k_eV_per_K = k,
                          .length_cm = 1e-3,
                          .density_cm3 = 1e17,
                          .voltage_V = 4.0,
                          .ambient_K = 300.0,
                          .kappa = {1.48, 0.0},
                          .mu_min = 0.0,
                          .mu_max = 1400.0,
                          .N_ref = 1e17,
                          .alpha = 1.0,
                          .gamma = 0.0,
                          .steps = 128,
                          .run_second_method = false};
    std::vector<double> x;
    for (int i = 0; i <= 10; ++i) x.push_back(1e-3 * i / 10.0);
    const auto p = ref::self_heated_resistor(in, x);
    REQUIRE(p.has_value());
    const double mu = 700.0, L = 1e-3, kappa = 1.48;
    const double J = -q * 1e17 * 4.0 * mu / L;
    CHECK(std::abs(p->current_A_cm2.to_double() / J - 1.0) <= 1e-28);
    const double Q = J * J / (q * mu * 1e17), b = 1.5 * k * J;
    double worst = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        // In double-double: the bracket cancels to about 1e-3 of its terms.
        const DD beta = DD(b) / DD(kappa), X(x[i]);
        const DD expected = DD(300.0) + DD(Q) / DD(b) * (DD(L) * expm1(-beta * X) / expm1(-beta * DD(L)) - X);
        worst = std::max(worst, std::abs((p->temperature[i] - expected).to_double()));
    }
    CAPTURE(worst, p->steps_difference);
    CHECK(worst <= 1e-12);
    CHECK(p->steps_difference <= 1e-12);
    // J < 0: the electrons move toward x = L and carry (3/2) k T each, so the peak lies right of the centre.
    CHECK(p->temperature[6].to_double() > p->temperature[4].to_double());
}

TEST_CASE("thermal_slab: the Seebeck rise and the Peltier coefficient", "[thermal_slab]") {
    const auto s = ref::seebeck_rise({.k_eV_per_K = k, .density_cm3 = 1e17, .Nc300_cm3 = 2.86e19,
                                      .thermopower_exponent = -0.5, .T0 = 300.0, .T1 = 400.0});
    REQUIRE(s.has_value());
    CHECK(s->second_method <= 1e-28);
    // n = 1e17 < Nc: (5/2 + r) - ln(N / Nc) > 0, so E_Fn falls toward the hot end.
    CHECK(s->rise_eV.to_double() < 0.0);
    const double Pi = ref::peltier_coefficient(k, 1e17, 2.86e19, -0.5, 300.0).to_double();
    CHECK(std::abs(Pi / (-k * 300.0 * (2.0 - std::log(1e17 / 2.86e19))) - 1.0) <= 1e-15);
}

TEST_CASE("thermal_slab: refusals", "[thermal_slab]") {
    CHECK_FALSE(ref::slab_temperature({1.48, 1.0}, 300.0, 400.0, {0.5}).has_value());
    CHECK_FALSE(ref::slab_temperature(si, 300.0, 400.0, {1.5}).has_value());
    ref::ResistorInput in{.q_C = q, .k_eV_per_K = k, .length_cm = 1e-3, .density_cm3 = 1e17, .voltage_V = 0.0,
                          .kappa = si, .mu_max = 1400.0, .N_ref = 1e17, .alpha = 1.0};
    CHECK_FALSE(ref::self_heated_resistor(in, {0.0}).has_value());
    CHECK_FALSE(ref::seebeck_rise({.k_eV_per_K = k, .density_cm3 = 0.0, .Nc300_cm3 = 2.86e19, .T0 = 300.0,
                                   .T1 = 400.0})
                    .has_value());
}
