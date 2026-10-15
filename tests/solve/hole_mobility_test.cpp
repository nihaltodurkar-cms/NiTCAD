// The selectable hole mobility (Unit 34, assemble::PhysicsModels::hole_mobility): the default is
// the legacy Caughey-Thomas, bit for bit (AC1); Masetti is refused where its silicon coefficients do
// not apply (AC8) and works with the electrothermal model.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
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

bool same_bits(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (std::bit_cast<std::uint64_t>(a[i]) != std::bit_cast<std::uint64_t>(b[i])) return false;
    }
    return true;
}

// The legacy 1e17 / 1e17 diode, its p side in `p_material`, optionally with sinks at both ends.
device::Device diode(const physics::Semiconductor& p_material, bool sinks = false) {
    mesh::Mesh m = *mesh::make_tensor_grid(legacy_graded_mesh(2e-4, 1e-4, 1e-8, 1e-6));
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    std::vector<device::RegionId> region(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        if (m.points()[i][0] < 1e-4) acceptors[i] = 1e17;
        else donors[i] = 1e17;
    }
    std::vector<device::ThermalContact> thermal;
    if (sinks) {
        thermal = {{"left", "x_min", device::ThermalContactKind::isothermal, 300.0, 0.0},
                   {"right", "x_max", device::ThermalContactKind::isothermal, 300.0, 0.0}};
    }
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    return *device::Device::create({.mesh = std::move(m),
                                    .temperature_K = 300.0,
                                    .regions = {{"p", p_material}},
                                    .node_region = std::move(region),
                                    .donors = std::move(donors),
                                    .acceptors = std::move(acceptors),
                                    .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                                                 {"cathode", device::ContactKind::ohmic, std::move(cathode)}},
                                    .thermal_contacts = std::move(thermal)});
}

const std::vector<std::vector<double>> points{{0.0, 0.0}, {0.3, 0.0}, {0.5, 0.0}};

}  // namespace

TEST_CASE("hole mobility: Caughey-Thomas is the default, bit for bit (AC1)", "[hole_mobility]") {
    const auto d = diode(physics::silicon());
    solve::BiasOptions explicit_ct;
    explicit_ct.models.hole_mobility = assemble::HoleMobility::caughey_thomas;
    const auto a = solve::sweep_bias(d, points);
    const auto b = solve::sweep_bias(d, points, explicit_ct);
    REQUIRE(a.has_value());
    REQUIRE(b.has_value());
    CHECK(a->run.input_identity == b->run.input_identity);
    for (std::size_t k = 0; k < points.size(); ++k) {
        CHECK(same_bits(a->points[k].terminal_current, b->points[k].terminal_current));
        CHECK(same_bits(a->points[k].fields.p_cm3, b->points[k].fields.p_cm3));
    }
    // Masetti changes the hole current and the run identity, and is named in the settings.
    solve::BiasOptions masetti;
    masetti.models.hole_mobility = assemble::HoleMobility::masetti;
    const auto c = solve::sweep_bias(d, points, masetti);
    REQUIRE(c.has_value());
    REQUIRE_FALSE(c->stopped.has_value());
    CHECK(c->run.input_identity != a->run.input_identity);
    CHECK(c->points.back().terminal_current[0] != a->points.back().terminal_current[0]);
    const auto named = [](const results::RunRecord& r) {
        return std::ranges::any_of(r.settings, [](const auto& s) { return s.first == "models.hole_mobility"; });
    };
    CHECK_FALSE(named(a->run));
    CHECK(named(c->run));
}

TEST_CASE("hole mobility: Masetti only where its silicon coefficients apply (AC8)", "[hole_mobility]") {
    solve::BiasOptions masetti;
    masetti.models.hole_mobility = assemble::HoleMobility::masetti;
    const auto germanium = *physics::Semiconductor::create(physics::germanium_parameters);
    const auto refused = solve::sweep_bias(diode(germanium), points, masetti);
    REQUIRE_FALSE(refused.has_value());
    CHECK(refused.error().code == base::ErrorCode::invalid_input);
    CHECK(refused.error().message.find("silicon's hole mobility") != std::string::npos);
    // Silicon with another parameter changed, its hole mobility kept, is accepted.
    auto params = physics::silicon_parameters;
    params.eps_r = 11.9;
    const auto accepted = solve::sweep_bias(diode(*physics::Semiconductor::create(params)), points, masetti);
    REQUIRE(accepted.has_value());
    CHECK_FALSE(accepted->stopped.has_value());
}

TEST_CASE("hole mobility: Masetti with the electrothermal model", "[hole_mobility]") {
    solve::BiasOptions options;
    options.models.hole_mobility = assemble::HoleMobility::masetti;
    options.models.electrothermal = true;
    const auto run = solve::sweep_bias(diode(physics::silicon(), true), points, options);
    REQUIRE(run.has_value());
    REQUIRE_FALSE(run->stopped.has_value());
    // Isothermal sinks at 300 K and a small current: the temperature stays within 1 K of 300 K, and
    // the current matches the isothermal Masetti run to 1e-3.
    for (const double T : run->points.back().fields.temperature_K) CHECK(std::abs(T - 300.0) < 1.0);
    solve::BiasOptions isothermal;
    isothermal.models.hole_mobility = assemble::HoleMobility::masetti;
    const auto iso = solve::sweep_bias(diode(physics::silicon()), points, isothermal);
    REQUIRE(iso.has_value());
    const double I = run->points.back().terminal_current[0];
    const double I_iso = iso->points.back().terminal_current[0];
    CAPTURE(I, I_iso);
    CHECK(std::abs(I / I_iso - 1.0) < 1e-3);
}

TEST_CASE("hole mobility: the cost of Masetti (AC9)", "[.performance]") {
    // A 2D p-n diode, 201 x 101 nodes, swept 0 to 0.5 V in 0.1 V steps: the best of three wall
    // times with each hole mobility. The mobility is computed once per node when the system is
    // assembled, so the criterion (AC9) is a ratio within 10%.
    std::vector<double> x, y;
    for (int k = 0; k <= 200; ++k) x.push_back(2e-4 * k / 200.0);
    for (int k = 0; k <= 100; ++k) y.push_back(1e-4 * k / 100.0);
    mesh::Mesh m = *mesh::make_tensor_grid(x, y);
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (m.points()[i][0] < 1e-4) acceptors[i] = 1e17;
        else donors[i] = 1e17;
    }
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    const auto d = *device::Device::create(
        {.mesh = std::move(m), .temperature_K = 300.0, .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0), .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
    std::vector<std::vector<double>> sweep;
    for (int k = 0; k <= 5; ++k) sweep.push_back({0.1 * k, 0.0});
    struct Cost {
        double seconds;
        std::size_t iterations;  // Newton iterations over the sweep's points
    };
    const auto best = [&](assemble::HoleMobility holes) {
        solve::BiasOptions options;
        options.models.hole_mobility = holes;
        Cost c{1e300, 0};
        for (int r = 0; r < 3; ++r) {
            const auto start = std::chrono::steady_clock::now();
            const auto run = solve::sweep_bias(d, sweep, options);
            const std::chrono::duration<double> s = std::chrono::steady_clock::now() - start;
            REQUIRE(run.has_value());
            REQUIRE_FALSE(run->stopped.has_value());
            c.seconds = std::min(c.seconds, s.count());
            c.iterations = 0;
            for (const auto& p : run->points) c.iterations += p.convergence.iterations.size();
        }
        return c;
    };
    const Cost ct = best(assemble::HoleMobility::caughey_thomas);
    const Cost masetti = best(assemble::HoleMobility::masetti);
    std::printf("hole mobility cost (2D diode, %zu nodes, 6 points, best of 3): Caughey-Thomas %.3f s "
                "(%zu Newton iterations), Masetti %.3f s (%zu), ratio %.3f, per iteration %.3f\n",
                n, ct.seconds, ct.iterations, masetti.seconds, masetti.iterations,
                masetti.seconds / ct.seconds,
                (masetti.seconds / static_cast<double>(masetti.iterations)) /
                    (ct.seconds / static_cast<double>(ct.iterations)));
    CHECK(masetti.seconds / ct.seconds <= 1.10);  // AC9
}
