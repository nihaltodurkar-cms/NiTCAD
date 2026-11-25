// Klaassen mobility in a solved device (Unit 66; DECISIONS.md KC6, KC7): the Unit 9 diode swept to
// 0.8 V with klaassen_equilibrium and klaassen converges at every bias; a p+ / n- diode in high
// injection shows the carrier-dependent part. Currents and Newton steps are printed for the record.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <utility>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "legacy_graded_mesh.hpp"

using namespace NiTCAD;

namespace {

// The Unit 9 diode (default_physics_test.cpp): legacy graded mesh, 1e17 / 1e17, silicon; or with
// other densities (the p side, then the n side).
device::Device diode(double NA = 1e17, double ND = 1e17) {
    mesh::Mesh m = *mesh::make_tensor_grid(legacy_graded_mesh(2e-4, 1e-4, 1e-8, 1e-6));
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const bool p_side = m.points()[i][0] < 1e-4;
        (p_side ? acceptors : donors)[i] = p_side ? NA : ND;
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

TEST_CASE("klaassen: the Unit 9 diode IV converges with klaassen_equilibrium (KC6)",
          "[klaassen]") {
    const device::Device d = diode();
    std::vector<std::vector<double>> points;
    for (int k = 0; k <= 16; ++k) points.push_back({0.05 * k, 0.0});
    const solve::BiasOptions masetti{
        .models = {.hole_mobility = assemble::HoleMobility::masetti,
                   .electron_mobility = assemble::ElectronMobility::masetti}};
    const solve::BiasOptions klaassen{
        .models = {.hole_mobility = assemble::HoleMobility::klaassen_equilibrium,
                   .electron_mobility = assemble::ElectronMobility::klaassen_equilibrium}};
    const auto m = solve::sweep_bias(d, points, masetti);
    const auto k = solve::sweep_bias(d, points, klaassen);
    REQUIRE(m.has_value());
    REQUIRE(k.has_value());
    REQUIRE_FALSE(m->stopped.has_value());
    REQUIRE_FALSE(k->stopped.has_value());
    REQUIRE(k->points.size() == points.size());
    for (std::size_t i = 0; i < points.size(); ++i) {
        const auto& pm = m->points[i];
        const auto& pk = k->points[i];
        std::printf("KC6 V=%.2f  J masetti %+.5e  klaassen_equilibrium %+.5e A/cm^2 (%+.2f%%)  "
                    "Newton %zu / %zu\n",
                    points[i][0], pm.terminal_current[0], pk.terminal_current[0],
                    pm.terminal_current[0] != 0.0
                        ? 100.0 * (pk.terminal_current[0] / pm.terminal_current[0] - 1.0)
                        : 0.0,
                    pm.convergence.iterations.size(), pk.convergence.iterations.size());
        CHECK(pk.convergence.converged);
        CHECK(std::isfinite(pk.terminal_current[0]));
    }
}

TEST_CASE("klaassen: the carrier-dependent model in a solved diode (KC7)", "[klaassen]") {
    using assemble::ElectronMobility;
    using assemble::HoleMobility;
    const solve::BiasOptions masetti{
        .models = {.hole_mobility = HoleMobility::masetti,
                   .electron_mobility = ElectronMobility::masetti}};
    const solve::BiasOptions equilibrium{
        .models = {.hole_mobility = HoleMobility::klaassen_equilibrium,
                   .electron_mobility = ElectronMobility::klaassen_equilibrium}};
    const solve::BiasOptions carriers{
        .models = {.hole_mobility = HoleMobility::klaassen,
                   .electron_mobility = ElectronMobility::klaassen}};
    const auto points = [](double top) {
        std::vector<std::vector<double>> v;
        for (int k = 0; k * 0.05 <= top + 1e-12; ++k) v.push_back({0.05 * k, 0.0});
        return v;
    };
    // Newton steps per bias: no more than Masetti's plus 2 on the Unit 9 diode.
    {
        const device::Device d = diode();
        const auto pts = points(0.8);
        const auto m = solve::sweep_bias(d, pts, masetti);
        const auto k = solve::sweep_bias(d, pts, carriers);
        REQUIRE(m.has_value());
        REQUIRE(k.has_value());
        REQUIRE_FALSE(k->stopped.has_value());
        for (std::size_t i = 0; i < pts.size(); ++i) {
            const std::size_t nm = m->points[i].convergence.iterations.size();
            const std::size_t nk = k->points[i].convergence.iterations.size();
            std::printf("KC7 Unit 9 V=%.2f  J masetti %+.5e  klaassen %+.5e  Newton %zu / %zu\n",
                        pts[i][0], m->points[i].terminal_current[0],
                        k->points[i].terminal_current[0], nm, nk);
            CHECK(k->points[i].convergence.converged);
            CHECK(nk <= nm + 2);
        }
    }
    // High injection: p+ 1e19 / n 1e15, forward to 0.9 V; the base's holes exceed 10 N_D and
    // carrier-carrier scattering lowers the current against the equilibrium-carrier model.
    {
        const device::Device d = diode(1e19, 1e15);
        const auto pts = points(0.9);
        const auto e = solve::sweep_bias(d, pts, equilibrium);
        const auto k = solve::sweep_bias(d, pts, carriers);
        REQUIRE(e.has_value());
        REQUIRE(k.has_value());
        REQUIRE_FALSE(e->stopped.has_value());
        REQUIRE_FALSE(k->stopped.has_value());
        for (std::size_t i = 0; i < pts.size(); ++i) {
            std::printf("KC7 high injection V=%.2f  J klaassen_equilibrium %+.5e  klaassen %+.5e "
                        "(%+.2f%%)  Newton %zu / %zu\n",
                        pts[i][0], e->points[i].terminal_current[0],
                        k->points[i].terminal_current[0],
                        e->points[i].terminal_current[0] != 0.0
                            ? 100.0 * (k->points[i].terminal_current[0] /
                                           e->points[i].terminal_current[0] -
                                       1.0)
                            : 0.0,
                        e->points[i].convergence.iterations.size(),
                        k->points[i].convergence.iterations.size());
            CHECK(k->points[i].convergence.converged);
        }
        const auto& last = k->points.back();
        double p_base = 0.0;
        for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
            const double xi = d.mesh().points()[i][0];
            if (xi > 1.2e-4 && xi < 1.8e-4) p_base = std::max(p_base, last.fields.p_cm3[i]);
        }
        const double ratio =
            last.terminal_current[0] / e->points.back().terminal_current[0];
        std::printf("KC7 high injection: largest base p in (1.2, 1.8) um %.3e cm^-3 (%.0f N_D); "
                    "J(0.9 V) klaassen / klaassen_equilibrium %.4f\n",
                    p_base, p_base / 1e15, ratio);
        CHECK(p_base > 10.0 * 1e15);
        CHECK(std::abs(ratio - 1.0) > 1e-3);
    }
}
