// Unit 68 slice 4a (MAD11; DECISIONS.md amendment 8, owner check 8): a 1D 4H-SiC p+ n diode on
// `sic_4h_2004` breaks down within 10 % of the measured (0001) fit of Hatakeyama et al. (2004),
//     V_BD = 1940 (1e16 / N_D)^0.8 V,  N_D = 4e16 cm^-3: 639.9 V,
// drift 10 um (not punched through: about 4 um depleted at breakdown), p+ 1e19 (the owner's
// benchmark assumption), 300 K. Breakdown is where the electron ionization integral on the simulated
// field reaches 1 (multiplication diverges), the criterion of the Unit 19 gates; and the device
// with impact ionization on settles at 90 % of it.
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdio>
#include <variant>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/impact_ionization.hpp"
#include "NiTCAD/physics/materials.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "legacy_graded_mesh.hpp"

using namespace NiTCAD;

namespace {

constexpr double junction_cm = 1e-4, length_cm = 11e-4, N_D = 4e16;

physics::Semiconductor sic() {
    const auto made = physics::find_material("sic_4h_2004")->make({});
    return *physics::Semiconductor::create(std::get<physics::SemiconductorParameters>(*made));
}

device::Device diode() {
    mesh::Mesh m = *mesh::make_tensor_grid(legacy_graded_mesh(length_cm, {junction_cm}, 2e-7, 1e-6));
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) (m.points()[i][0] < junction_cm ? acceptors[i] : donors[i]) = m.points()[i][0] < junction_cm ? 1e19 : N_D;
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    return *device::Device::create({.mesh = std::move(m),
                                    .temperature_K = 300.0,
                                    .regions = {{"sic", sic()}},
                                    .node_region = std::vector<device::RegionId>(n, 0),
                                    .donors = std::move(donors),
                                    .acceptors = std::move(acceptors),
                                    .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                                                 {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
}

// The electron ionization integral from x_min (electrons enter the depletion region from the p+
// side), I_n = integral alpha_n exp(-integral (alpha_n - alpha_p)), each edge at its own field.
double ionization_integral(const device::Device& d, const std::vector<double>& psi) {
    static const physics::Semiconductor m = sic();
    const auto& en = physics::impact_ionization(m, physics::Carrier::electron);
    const auto& ep = physics::impact_ionization(m, physics::Carrier::hole);
    double phi = 0.0, I = 0.0;
    for (std::size_t i = 0; i + 1 < psi.size(); ++i) {
        const double h = d.mesh().points()[i + 1][0] - d.mesh().points()[i][0];
        const double E = std::abs(psi[i + 1] - psi[i]) / h;
        const double an = physics::impact_ionization_coefficient(en, 1.0, E).alpha;
        const double ap = physics::impact_ionization_coefficient(ep, 1.0, E).alpha;
        I += an * std::exp(-(phi + 0.5 * (an - ap) * h)) * h;
        phi += (an - ap) * h;
    }
    return I;
}

}  // namespace

TEST_CASE("sic_4h_2004: a p+ n diode breaks down within 10 % of Hatakeyama's (0001) fit", "[sic][unit68]") {
    const device::Device d = diode();
    const double target = 1940.0 * std::pow(1e16 / N_D, 0.8);
    // The field profile, reverse bias in 10 V steps (impact ionization off: below breakdown it does
    // not move the field).
    std::vector<std::vector<double>> points;
    for (int k = 0; k <= 80; ++k) points.push_back({-10.0 * k, 0.0});
    // Sub-steps on: the wide gap (n_i about 1e-8 cm^-3) needs them from the first reverse step.
    solve::BiasOptions off;
    off.step_control.max_cuts = 20;
    const auto s = solve::sweep_bias(d, points, off);
    REQUIRE(s.has_value());
    if (s->stopped) std::printf("stopped after %zu points: %s\n", s->points.size(), s->stopped->message.c_str());
    REQUIRE(!s->stopped);
    double V_bd = 0.0, before = 0.0;
    for (const results::BiasPoint& p : s->points) {
        const double I = ionization_integral(d, p.fields.potential_V), V = -p.bias_V[0];
        if (I >= 1.0) {
            V_bd = V - 10.0 * (I - 1.0) / (I - before);  // linear between the bracketing points
            break;
        }
        before = I;
    }
    std::printf("sic_4h_2004 p+ n (N_D 4e16, 10 um): breakdown %.1f V, Hatakeyama fit %.1f V (%+.1f %%)\n", V_bd, target,
                100.0 * (V_bd / target - 1.0));
    REQUIRE(V_bd > 0.0);
    CHECK(std::abs(V_bd / target - 1.0) <= 0.10);
}

// Amendments 8, 11: with impact ionization on and the plain Newton step the diode does not settle past about
// 150 V (n_i about 1e-8 cm^-3: the minority rows unresolved); with the log-density step it reaches 90 %.
TEST_CASE("sic_4h_2004: the p+ n diode with impact ionization on settles at 90 % of breakdown", "[sic][unit68]") {
    const device::Device d = diode();
    const double V_bd = 623.0;  // the ionization-integral breakdown of the test above
    solve::BiasOptions on;
    on.models.impact_ionization = true;
    on.step_control.max_cuts = 20;
    on.log_densities = true;  // amendment 11
    std::vector<std::vector<double>> up;
    for (double V = 0.0; V <= 0.9 * V_bd; V += 10.0) up.push_back({-V, 0.0});
    up.push_back({-0.9 * V_bd, 0.0});
    const auto t = solve::sweep_bias(d, up, on);
    REQUIRE(t.has_value());
    if (t->stopped) std::printf("impact on: stopped after %zu of %zu points: %s\n", t->points.size(), up.size(), t->stopped->message.c_str());
    CHECK(!t->stopped);
}
