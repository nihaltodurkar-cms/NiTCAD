// Band-gap narrowing and Auger in the device solves (ARCHITECTURE.md section 11, Unit 11).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/bandgap_narrowing.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/equilibrium.hpp"
#include "legacy_graded_mesh.hpp"

using namespace NiTCAD;

namespace {

device::Device diode(double NA, double ND) {
    mesh::Mesh m = *mesh::make_tensor_grid(legacy_graded_mesh(2e-4, 1e-4, 1e-8, 1e-6));
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (m.points()[i][0] < 1e-4) {
            acceptors[i] = NA;
        } else {
            donors[i] = ND;
        }
    }
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
}

}  // namespace

TEST_CASE("heavy doping: equilibrium carries no current where n_ie varies") {
    // 1e19 / 1e18: n_ie differs by a factor 1.8 across the junction. Without the ln(n_ie) term in
    // the Scharfetter-Gummel driving force the junction edges would carry a current of the order
    // of the one-sided fluxes (q D n / h, up to 4e8 A/cm^2 on the 1e-8 cm cells). With it every
    // edge is zero to rounding: each is held to 64 eps of its one-sided terms, bounded by
    // q mu_max V_T / h (n_a + n_b) (2 + |delta psi| / V_T) (B(-x) = B(x) + x <= 1 + |x|).
    const auto d = diode(1e19, 1e18);
    const std::vector<double> zero{0.0, 0.0};
    const auto s = solve::solve_bias(d, zero);
    REQUIRE(s.has_value());
    const double VT = base::thermal_voltage(300.0);
    const auto x = legacy_graded_mesh(2e-4, 1e-4, 1e-8, 1e-6);
    const auto& f = s->fields;
    double worst = 0.0;
    for (std::size_t k = 0; k + 1 < x.size(); ++k) {
        const double drive = 2.0 + std::abs(f.potential_V[k + 1] - f.potential_V[k]) / VT;
        const double prefactor = base::q_C * 1360.0 * VT / (x[k + 1] - x[k]) * drive;
        const double floor = 64.0 * 2.2204460492503131e-16 * prefactor;
        const double bound_n = floor * (f.n_cm3[k] + f.n_cm3[k + 1]);
        const double bound_p = floor * (f.p_cm3[k] + f.p_cm3[k + 1]);
        CAPTURE(k);
        REQUIRE(std::abs(s->edge_current_n[k]) <= bound_n);
        REQUIRE(std::abs(s->edge_current_p[k]) <= bound_p);
        worst = std::max({worst, std::abs(s->edge_current_n[k]) / bound_n,
                          std::abs(s->edge_current_p[k]) / bound_p});
    }
    UNSCOPED_INFO("worst edge current over its rounding bound " << worst);
    // Built-in potential from the effective n_ie of each side.
    const physics::Semiconductor si = physics::silicon();
    const double nie_p = physics::effective_intrinsic_density(si, 1e19, 300.0);
    const double nie_n = physics::effective_intrinsic_density(si, 1e18, 300.0);
    const double vbi = VT * std::log(1e19 * 1e18 / (nie_p * nie_n));
    const auto& phi = s->fields.potential_V;
    REQUIRE(std::abs((phi.back() - phi.front()) - vbi) <= 1e-12 * vbi);
}

TEST_CASE("heavy doping: band-gap narrowing raises the low-injection current by exp(dEg/kT)") {
    // Symmetric 1e19 / 1e19 diode at 0.5 V. The injected minority densities are
    // n_ie^2 / N e^(V/V_T) on both sides, and n_ie^2 grows by exp(dEg/kT) = 8.593
    // (dEg(1e19) = 55.6 meV). The diffusion current follows it exactly when the minority problem
    // is linear: no SRH (Auger, linear in the minority density at low injection, stays on), and
    // only the depletion width moves (V_bi drops by dEg).
    // With SRH the 0.3 V current of this junction is depletion-region recombination, which scales
    // with n_ie, not n_ie^2: measured ratio 3.14, close to sqrt(8.593) = 2.93.
    // The bias is 0.5 V because in a 1e19 device extracted currents below about 1e-8 A/cm^2 are
    // rounding noise (the majority one-sided fluxes are 5e6 A/cm^2 even on 1e-6 cm cells; measured
    // edge-to-edge spread of the total current 5e-8 to 1e-7 A/cm^2): at 0.2-0.3 V the ratio came
    // out 4.5 to 11. At 0.5 V the current is about 1e-5 A/cm^2.
    const auto d = diode(1e19, 1e19);
    const std::vector<double> bias{0.5, 0.0};
    const auto with = solve::solve_bias(d, bias, {.models = {.srh = false}});
    const auto without = solve::solve_bias(d, bias, {.models = {.srh = false, .bgn = false}});
    REQUIRE(with.has_value());
    REQUIRE(without.has_value());
    const double ratio = with->terminal_current[0] / without->terminal_current[0];
    const double expected = std::exp(physics::bandgap_narrowing_eV(physics::silicon(), 1e19) /
                                     base::thermal_voltage(300.0));
    UNSCOPED_INFO("ratio " << ratio << ", exp(dEg/kT) " << expected);
    REQUIRE(std::abs(ratio / expected - 1.0) < 0.005);  // measured 0.03%
}

TEST_CASE("heavy doping: Auger barely moves the legacy J(0.5 V) fixture (Unit 9 measurement)") {
    const auto d = diode(1e17, 1e17);
    const std::vector<double> bias{0.5, 0.0};
    const double J = solve::solve_bias(d, bias)->terminal_current[0];
    const double J_no_auger =
        solve::solve_bias(d, bias, {.models = {.auger = false}})->terminal_current[0];
    const double effect = J / J_no_auger - 1.0;
    UNSCOPED_INFO("Auger effect " << effect);
    REQUIRE(effect > 0.0);   // Auger adds recombination current
    REQUIRE(effect < 1e-6);  // measured 4.3e-7 at Unit 9
    REQUIRE(std::abs(J - 1.280e-2) / 1.280e-2 < 1e-2);
}
