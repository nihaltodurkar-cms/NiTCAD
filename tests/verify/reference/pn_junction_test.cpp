// L0 tests of the pn-junction reference (verify/reference/pn_junction.hpp; DECISIONS Unit 55 14.1 (c), TVB-07).
//   - symmetric doping: u_j = 0 by symmetry, to a few u_dd of the bulk level, and psi(-x) = -psi(x);
//   - the junction condition Phi_n(u_j) = Phi_p(u_j), re-evaluated here with G(t) = expm1(t) - t (not the
//     library's series), to 1e-28 relative of Phi;
//   - the profile against the ODE itself: u'' = (q / (eps V_T)) (n - p - C) integrated outward from the junction
//     with (u_j, E_max / V_T) as its only data, by fixed-step dopri5 in double-double, on each side to the
//     depletion edge plus 3 Debye lengths; agreement to 1e-15 of V_bi (the shooting mode grows ~ e^{x / L_D});
//   - the depletion approximation with the -2 V_T correction, E_max^2 = 2 q N_eff (V_bi - 2 V_T) / eps, its known
//     asymptotic form, to 1e-3 (a sanity bound, not a reference);
//   - the collocation (I3) on the benchmark's case: two degrees (40, 48) agree to 1e-11 V, and with the first integral
//     to within the contact tail it reports;
//   - bad inputs are refused.
#include "NiTCAD/verify/reference/pn_junction.hpp"

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/verify/reference/dopri5.hpp"

namespace ref = NiTCAD::verify::reference;
using ref::DD;
using ref::u_dd;

namespace {

// Silicon-like data at 300 K (values as NiTCAD's Couderc set gives them; any positive data serve here).
ref::PnJunctionInput input(double NA, double ND) {
    return {.q_C = 1.602176634e-19,
            .eps_F_per_cm = 11.7 * 8.8541878128e-14,
            .thermal_voltage_V = 0.025851999786435535,
            .ni_cm3 = 9.673e9,
            .acceptors_cm3 = NA,
            .donors_cm3 = ND};
}

DD G_expm1(const DD& t) { return ref::expm1(t) - t; }

}  // namespace

TEST_CASE("pn junction: symmetric doping puts u_j at 0 and the profile is odd") {
    const auto j = ref::pn_junction(input(1e16, 1e16));
    REQUIRE(j.has_value());
    INFO("u_j " << j->u_j.hi << ", u_n " << j->u_n.hi);
    CHECK(std::abs(j->u_j.to_double()) <= 8.0 * u_dd * j->u_n.to_double());
    CHECK(j->built_in_V.to_double() > 0.0);
    const std::vector<double> x{-3e-5, -1e-5, -1e-6, 1e-6, 1e-5, 3e-5};
    const auto p = ref::pn_potential(*j, x);
    REQUIRE(p.has_value());
    for (std::size_t i = 0; i < 3; ++i) {
        const double a = (*p)[i].psi_V.to_double(), b = (*p)[5 - i].psi_V.to_double();
        INFO("x = " << x[i] << ": " << a << " and " << b);
        CHECK(std::abs(a + b) <= 1e-20);
    }
}

TEST_CASE("pn junction: the junction condition holds, with G evaluated another way") {
    const auto j = ref::pn_junction(input(1e16, 2e16));
    REQUIRE(j.has_value());
    const DD phin = j->n_n * G_expm1(j->u_j - j->u_n) + j->p_n * G_expm1(j->u_n - j->u_j);
    const DD phip = j->n_p * G_expm1(j->u_j - j->u_p) + j->p_p * G_expm1(j->u_p - j->u_j);
    const double rel = std::abs((phin - phip).to_double()) / phin.to_double();
    INFO("relative mismatch " << rel << ", Brent evaluations " << j->brent_evaluations);
    CHECK(rel <= 1e-28);
    CHECK(j->u_p < j->u_j);
    CHECK(j->u_j < j->u_n);
    // Neutral bulks: n_b - p_b is the net doping, n_b p_b = n_i^2.
    CHECK(std::abs((j->n_n - j->p_n).to_double() / 2e16 - 1.0) <= 1e-30);
    CHECK(std::abs((j->p_p - j->n_p).to_double() / 1e16 - 1.0) <= 1e-30);
    CHECK(std::abs((j->n_n * j->p_n).to_double() / (9.673e9 * 9.673e9) - 1.0) <= 1e-30);
}

TEST_CASE("pn junction: the profile agrees with the ODE integrated from the junction") {
    const ref::PnJunctionInput in = input(1e16, 2e16);
    const auto j = ref::pn_junction(in);
    REQUIRE(j.has_value());
    const DD c = DD(in.q_C) / (DD(in.eps_F_per_cm) * DD(in.thermal_voltage_V));  // q / (eps V_T)
    const DD ni(in.ni_cm3);
    const double vbi = j->built_in_V.to_double();
    for (const double side : {1.0, -1.0}) {
        const DD C = side > 0 ? DD(in.donors_cm3) : -DD(in.acceptors_cm3);
        const double reach = side > 0 ? (j->depletion_n_cm + DD(3.0) * j->debye_n_cm).to_double()
                                      : (j->depletion_p_cm + DD(3.0) * j->debye_p_cm).to_double();
        // y = (u, u'); x from 0 to side * reach. u' > 0 everywhere (u increases from p to n).
        const auto f = [&](const DD&, const ref::DDVector& y, ref::DDVector& k) {
            k[0] = y[1];
            k[1] = c * (ni * ref::exp(y[0]) - ni * ref::exp(-y[0]) - C);
        };
        const DD slope = j->field_max_V_per_cm / DD(in.thermal_voltage_V);
        const std::vector<double> stops{0.25, 0.5, 0.75, 1.0};
        ref::DDVector y{j->u_j, slope};
        DD x0(0.0);
        std::vector<double> xs;
        std::vector<double> us;
        for (const double s : stops) {
            const DD x1(side * reach * s);
            y = ref::dopri5_fixed(f, x0, y, x1, 500);
            x0 = x1;
            xs.push_back(x1.to_double());
            us.push_back((y[0] * DD(in.thermal_voltage_V)).to_double());
        }
        const auto p = ref::pn_potential(*j, xs);
        REQUIRE(p.has_value());
        for (std::size_t i = 0; i < xs.size(); ++i) {
            const double e = std::abs(us[i] - (*p)[i].psi_V.to_double()) / vbi;
            INFO("side " << side << ", x = " << xs[i] << " cm: ODE " << us[i] << ", first integral "
                         << (*p)[i].psi_V.to_double() << ", relative " << e);
            CHECK(e <= 1e-15);
        }
    }
}

TEST_CASE("pn junction: E_max against the depletion approximation with the -2 V_T correction") {
    // Both sides must be well depleted (u_b - u_j >> 1 on each), which a 1:100 junction is not on its heavy side.
    for (const auto& [NA, ND] : {std::pair{1e16, 2e16}, std::pair{1e16, 1e16}, std::pair{1e17, 1e17}}) {
        const ref::PnJunctionInput in = input(NA, ND);
        const auto j = ref::pn_junction(in);
        REQUIRE(j.has_value());
        const double Neff = NA * ND / (NA + ND);
        const double approx =
            std::sqrt(2.0 * in.q_C * Neff * (j->built_in_V.to_double() - 2.0 * in.thermal_voltage_V) / in.eps_F_per_cm);
        const double rel = std::abs(j->field_max_V_per_cm.to_double() / approx - 1.0);
        INFO("N_A " << NA << ", N_D " << ND << ": E_max " << j->field_max_V_per_cm.to_double() << ", approximation "
                    << approx << ", relative " << rel);
        CHECK(rel <= 1e-3);
    }
}

TEST_CASE("pn junction: Chebyshev collocation agrees with itself and with the first integral (I3)") {
    const auto j = ref::pn_junction(input(1e16, 2e16));
    REQUIRE(j.has_value());
    std::vector<double> x;
    for (int k = 0; k < 64; ++k) x.push_back(1e-4 * (2.0 * k - 63) / 63);
    const auto a = ref::pn_collocation(*j, -1e-4, 1e-4, 4, 40, x);
    const auto b = ref::pn_collocation(*j, -1e-4, 1e-4, 4, 48, x);
    const auto p = ref::pn_potential(*j, x);
    REQUIRE(a.has_value());
    INFO((b.has_value() ? std::string() : b.error()));
    REQUIRE(b.has_value());
    INFO((p.has_value() ? std::string() : p.error()));
    REQUIRE(p.has_value());
    const double VT = j->in.thermal_voltage_V;
    const double tail = VT * std::max((*p)[0].tail.to_double(), (*p)[63].tail.to_double());
    double self = 0.0, cross = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        self = std::max(self, std::abs(a->psi_V[i] - b->psi_V[i]));
        cross = std::max(cross, std::abs(a->psi_V[i] - (*p)[i].psi_V.to_double()));
    }
    INFO("degree 40 vs 48: " << self << " V; against the first integral " << cross << " V, contact tail " << tail
                             << " V; tail estimates " << a->tail_estimate_V << ", " << b->tail_estimate_V);
    CHECK(self <= 1e-11);
    CHECK(cross <= tail + a->tail_estimate_V + 1e-11);
    CHECK(cross >= 0.5 * tail);  // the difference is the truncation, which the collocation sees
}

TEST_CASE("pn junction: bad inputs are refused") {
    CHECK_FALSE(ref::pn_junction(input(0.0, 1e16)).has_value());
    CHECK_FALSE(ref::pn_junction(input(1e16, -1.0)).has_value());
    CHECK_FALSE(ref::pn_junction(input(1e16, std::nan(""))).has_value());
    const auto j = ref::pn_junction(input(1e16, 1e16));
    REQUIRE(j.has_value());
    const std::vector<double> x{0.0};
    CHECK_FALSE(ref::pn_collocation(*j, 1e-4, 2e-4, 4, 40, x).has_value());
    CHECK_FALSE(ref::pn_collocation(*j, -1e-4, 1e-4, 0, 40, x).has_value());
}
