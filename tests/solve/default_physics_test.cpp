// The published defaults (Unit 35, DECISIONS.md S1-S5, DC9): the legacy-port gates' devices on the
// defaults, their values recorded before (legacy options) and after, and the diffusion current's
// change by the n_i^2 ratio of the two silicon band sets.
#include <catch2/catch_test_macros.hpp>

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

// The Unit 9 diode (bias_test.cpp): legacy graded mesh, 1e17 / 1e17, on `material`.
device::Device diode(const physics::Semiconductor& material) {
    mesh::Mesh m = *mesh::make_tensor_grid(legacy_graded_mesh(2e-4, 1e-4, 1e-8, 1e-6));
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) (m.points()[i][0] < 1e-4 ? acceptors : donors)[i] = 1e17;
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", material}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
}

double current(const device::Device& d, double V, const solve::BiasOptions& options = {}) {
    const std::vector<std::vector<double>> points{
        {0.0, 0.0}, {0.25 * V, 0.0}, {0.5 * V, 0.0}, {0.75 * V, 0.0}, {V, 0.0}};
    const auto s = solve::sweep_bias(d, points, options);
    REQUIRE(s.has_value());
    REQUIRE_FALSE(s->stopped.has_value());
    return s->points.back().terminal_current[0];
}

const solve::BiasOptions legacy{.models = {.hole_mobility = assemble::HoleMobility::caughey_thomas}};

}  // namespace

TEST_CASE("defaults: the Unit 9 diode, J(0.5 V) before and after (DC9)", "[defaults]") {
    // Before (legacy silicon, Caughey-Thomas holes; bias_test.cpp's gate): 1.280e-2 A/cm^2. After
    // (the defaults): 1.1399e-2, -10.9%. The diffusion current falls with n_i^2 (0.821, the next
    // case), and its hole part, the holes injected into the 1e17 n side, rises with their Masetti
    // mobility there (319 against Caughey-Thomas's 232 cm^2/(V s)): 0.821 (1 + 0.37 s) = 0.891
    // for the hole share s = 0.23 of the current.
    const double before = current(diode(physics::legacy_silicon()), 0.5, legacy);
    const double after = current(diode(physics::silicon()), 0.5);
    std::printf("Unit 9 diode J(0.5 V): legacy options %.5e, defaults %.5e A/cm^2 (%+.1f%%)\n", before,
                after, 100.0 * (after / before - 1.0));
    CHECK(std::abs(before / 1.280e-2 - 1.0) < 1e-2);  // the Unit 9 gate
    CHECK(std::abs(after / 1.1399e-2 - 1.0) < 1e-3);  // the default, pinned
}

TEST_CASE("defaults: a diffusion current changes by the n_i^2 ratio (DC9)", "[defaults]") {
    // The hole mobility held at Caughey-Thomas, so only the band set changes: in the ideal-diode
    // range the current is q n_i^2 (D_n / (N_A W_p) + D_p / (N_D W_n)) (e^(V/V_T) - 1), so the ratio
    // of the two sets' currents is (9.673e9 / 1.0674e10)^2 = 0.821 (DC9: within 2%).
    const double ni_new = physics::intrinsic_density(physics::silicon(), 300.0);
    const double ni_old = physics::intrinsic_density(physics::legacy_silicon(), 300.0);
    const double expected = (ni_new / ni_old) * (ni_new / ni_old);
    for (const double V : {0.4, 0.5}) {
        const double ratio = current(diode(physics::silicon()), V, legacy) /
                             current(diode(physics::legacy_silicon()), V, legacy);
        std::printf("n_i^2 ratio at %.1f V: measured %.4f, expected %.4f (%+.2f%%)\n", V, ratio,
                    expected, 100.0 * (ratio / expected - 1.0));
        CAPTURE(V, ratio, expected);
        CHECK(std::abs(expected - 0.821) < 0.001);
        CHECK(std::abs(ratio / expected - 1.0) < 0.02);
    }
}

TEST_CASE("defaults: the run identity follows the band forms (Unit 35)", "[defaults]") {
    // The identity digests the Passler gap and the mass fits when given (and only then, so the
    // legacy sets keep theirs).
    const std::vector<std::vector<double>> points{{0.0, 0.0}};
    const auto identity = [&](const physics::SemiconductorParameters& p) {
        return solve::make_run_record(diode(*physics::Semiconductor::create(p)), {}, points)
            .input_identity;
    };
    const auto base = identity(physics::silicon_parameters);
    CHECK(base != identity(physics::legacy_silicon_parameters));
    physics::SemiconductorParameters p = physics::silicon_parameters;
    p.passler_gap.Theta_K = 440.0;
    CHECK(identity(p) != base);
    p = physics::silicon_parameters;
    p.passler_gap.Delta = 0.5;
    CHECK(identity(p) != base);
    p = physics::silicon_parameters;
    p.dos_mass.conduction[2] *= 1.01;
    CHECK(identity(p) != base);
    p = physics::silicon_parameters;
    p.dos_mass.valence[0] *= 1.01;
    CHECK(identity(p) != base);
}
