// The weighted terminal current (Unit 32, solve::BiasOptions::weighted_current): equal to the
// contact current where both are resolved, resolved where the contact current is rounding noise,
// and without effect on anything else.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <bit>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <utility>
#include <vector>

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/trace.hpp"
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

// A 1D silicon diode: acceptors N_A for x < junction, donors N_D beyond (the legacy layout).
device::Device diode(const std::vector<double>& x, double junction, double N_A, double N_D) {
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (m.points()[i][0] < junction) acceptors[i] = N_A;
        else donors[i] = N_D;
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

std::vector<std::vector<double>> anode_points(double first, double last, double step) {
    std::vector<std::vector<double>> points;
    for (int k = 0; first + k * step <= last + 1e-12; ++k) points.push_back({first + k * step, 0.0});
    return points;
}

}  // namespace

TEST_CASE("weighted current: equal to the contact current where both are resolved", "[weighted]") {
    // The Unit 10 legacy diode, 1e17 / 1e17.
    const auto d = diode(legacy_graded_mesh(2e-4, 1e-4, 1e-8, 1e-6), 1e-4, 1e17, 1e17);
    solve::BiasOptions on;
    on.weighted_current = true;
    const auto points = anode_points(0.0, 0.6, 0.1);
    const auto sweep = solve::sweep_bias(d, points, on);
    REQUIRE(sweep.has_value());
    REQUIRE_FALSE(sweep->stopped.has_value());
    for (const results::BiasPoint& p : sweep->points) {
        REQUIRE(p.weighted_current.size() == 2);
        REQUIRE(p.weighted_current_resolution.size() == 2);
        for (std::size_t c = 0; c < 2; ++c) {
            const double tolerance =
                p.weighted_current_resolution[c] + p.terminal_current_resolution[c];
            INFO("V = " << p.bias_V[0] << ", contact " << c << ": " << p.weighted_current[c]
                        << " vs " << p.terminal_current[c] << ", tolerance " << tolerance);
            CHECK(std::abs(p.weighted_current[c] - p.terminal_current[c]) <= tolerance);
            CHECK(p.weighted_current_resolution[c] > 0.0);
        }
        // The weights of the two contacts sum to 1 at every node, so the currents sum to zero.
        CHECK(std::abs(p.weighted_current[0] + p.weighted_current[1]) <=
              1e-12 * std::abs(p.weighted_current[0]) + p.weighted_current_resolution[0]);
    }
    // The 0.5 V Unit 9 gate, J = 1.280e-2 A/cm^2 within 1%.
    CHECK(std::abs(sweep->points[5].weighted_current[0] / 1.280e-2 - 1.0) < 0.01);
}

TEST_CASE("weighted current: off by default, and on it changes nothing else", "[weighted]") {
    const auto d = diode(legacy_graded_mesh(2e-4, 1e-4, 1e-8, 1e-6), 1e-4, 1e17, 1e17);
    const auto points = anode_points(0.0, 0.3, 0.1);
    const auto off = solve::sweep_bias(d, points);
    solve::BiasOptions options;
    options.weighted_current = true;
    const auto on = solve::sweep_bias(d, points, options);
    REQUIRE(off.has_value());
    REQUIRE(on.has_value());
    CHECK(off->points.front().weighted_current.empty());
    CHECK(off->points.front().weighted_current_resolution.empty());
    CHECK(off->run.input_identity != on->run.input_identity);
    REQUIRE(off->points.size() == on->points.size());
    for (std::size_t k = 0; k < on->points.size(); ++k) {
        const auto& a = off->points[k];
        const auto& b = on->points[k];
        CHECK(same_bits(a.fields.potential_V, b.fields.potential_V));
        CHECK(same_bits(a.fields.n_cm3, b.fields.n_cm3));
        CHECK(same_bits(a.terminal_current, b.terminal_current));
        CHECK(same_bits(a.terminal_current_resolution, b.terminal_current_resolution));
        CHECK(a.convergence.iterations.size() == b.convergence.iterations.size());
    }
    // The settings name it only when on.
    const auto named = [](const results::RunRecord& r) {
        return std::ranges::any_of(r.settings, [](const auto& s) { return s.first == "weighted_current"; });
    };
    CHECK_FALSE(named(off->run));
    CHECK(named(on->run));
}

TEST_CASE("weighted current: the ideal-diode law below the contact current's resolution",
          "[weighted]") {
    // An n+ (1e19) / p (1e17) silicon diode with every recombination off and no band-gap
    // narrowing: the current is diffusion through the neutral regions, J = J0 (exp(V / V_T) - 1)
    // with J0 nearly constant (the 20 um p region's width changes by the depletion edge's shift,
    // well below 1%). At low bias the contact current, the difference of the n+ region's huge
    // fluxes, is below its resolution; the weighted current must follow the law there.
    std::vector<double> x;
    for (int k = 0; k <= 200; ++k) x.push_back(20e-4 * k / 200.0);           // p, 0.1 um
    for (int k = 1; k <= 100; ++k) x.push_back(20e-4 + 1e-4 * k / 100.0);    // n+, 10 nm
    const auto d = diode(x, 20e-4 + 1e-9, 1e17, 1e19);
    solve::BiasOptions options;
    options.weighted_current = true;
    options.models.srh = false;
    options.models.auger = false;
    options.models.bgn = false;
    const auto points = anode_points(0.0, 0.5, 0.05);
    const auto sweep = solve::sweep_bias(d, points, options);
    REQUIRE(sweep.has_value());
    REQUIRE_FALSE(sweep->stopped.has_value());
    const double V_T = 0.025852;
    // J0 from the highest bias, where both currents are resolved.
    const auto& top = sweep->points.back();
    REQUIRE(std::abs(top.terminal_current[0]) > 10.0 * top.terminal_current_resolution[0]);
    const double J0 = top.weighted_current[0] / std::expm1(top.bias_V[0] / V_T);
    int unresolved_conventionally = 0;
    for (std::size_t k = 1; k < sweep->points.size(); ++k) {
        const auto& p = sweep->points[k];
        const double V = p.bias_V[0];
        const double J = p.weighted_current[0];
        if (std::abs(p.terminal_current[0]) <= p.terminal_current_resolution[0]) ++unresolved_conventionally;
        CHECK(std::abs(J) > 100.0 * p.weighted_current_resolution[0]);
        INFO("V = " << V << ": weighted " << J << " (resolution " << p.weighted_current_resolution[0]
                    << "), contact " << p.terminal_current[0] << " (resolution "
                    << p.terminal_current_resolution[0] << ")");
        CHECK(std::abs(J / (J0 * std::expm1(V / V_T)) - 1.0) < 5e-3);
    }
    // The fixture is one where the contact current fails: most low-bias points are below its floor.
    CHECK(unresolved_conventionally >= 3);
}

TEST_CASE("weighted current: the resolution covers a loosely converged state's error",
          "[weighted]") {
    // Newton stopped at tol_update = 0.5 leaves a residual whose effect on the current is the
    // adjoint estimate: the resolution must cover the difference to a tightly converged state, and
    // be far larger than that state's own.
    const auto d = diode(legacy_graded_mesh(2e-4, 1e-4, 1e-8, 1e-6), 1e-4, 1e17, 1e17);
    const auto points = anode_points(0.0, 0.5, 0.1);
    solve::BiasOptions loose, tight;
    loose.weighted_current = tight.weighted_current = true;
    loose.newton.tol_update = 0.5;
    tight.newton.tol_update = 1e-12;
    const auto a = solve::sweep_bias(d, points, loose);
    const auto b = solve::sweep_bias(d, points, tight);
    REQUIRE(a.has_value());
    REQUIRE(b.has_value());
    REQUIRE_FALSE(a->stopped.has_value());
    REQUIRE_FALSE(b->stopped.has_value());
    int looser = 0;
    for (std::size_t k = 1; k < points.size(); ++k) {
        const auto& p = a->points[k];
        const auto& q = b->points[k];
        const double error = std::abs(p.weighted_current[0] - q.weighted_current[0]);
        CAPTURE(k, p.weighted_current[0], q.weighted_current[0], error,
                p.weighted_current_resolution[0], q.weighted_current_resolution[0]);
        CHECK(error <= p.weighted_current_resolution[0] + q.weighted_current_resolution[0]);
        // Where the loose state is really off, the estimate is its error (one more Newton step
        // nearly corrects it): measured 3e-6 relative at 0.1 V, where the current is 22 times high.
        if (error > 100.0 * q.weighted_current_resolution[0]) {
            ++looser;
            CHECK(std::abs(p.weighted_current_resolution[0] / error - 1.0) < 0.01);
        }
    }
    CHECK(looser >= 2);
}

TEST_CASE("weighted current: what it does not support", "[weighted]") {
    const auto d = diode(legacy_graded_mesh(2e-4, 1e-4, 1e-8, 1e-6), 1e-4, 1e17, 1e17);
    const std::vector<std::vector<double>> point{{0.1, 0.0}};
    solve::BiasOptions quasi_static;
    quasi_static.weighted_current = true;
    quasi_static.equations = solve::Equations::equilibrium_poisson;
    const std::vector<std::vector<double>> zero{{0.0, 0.0}};
    const auto a = solve::sweep_bias(d, zero, quasi_static);
    REQUIRE_FALSE(a.has_value());
    CHECK(a.error().message.find("drift-diffusion") != std::string::npos);
    solve::BiasOptions nonlocal;
    nonlocal.weighted_current = true;
    nonlocal.models.btbt_nonlocal = assemble::NonlocalTunnelling::kane;
    const auto b = solve::sweep_bias(d, point, nonlocal);
    REQUIRE_FALSE(b.has_value());
    // The option's own refusal, not the assembler's (this diode has no cells for the paths).
    CHECK(b.error().message.find("not carried along the mesh edges") != std::string::npos);
}

TEST_CASE("weighted current: every traced point carries it", "[weighted]") {
    const auto d = diode(legacy_graded_mesh(2e-4, 1e-4, 1e-8, 1e-6), 1e-4, 1e17, 1e17);
    solve::TraceOptions t;
    t.steady.weighted_current = true;
    t.contact = 0;
    t.start_V = {0.0, 0.0};
    t.end_V = 0.3;
    t.step_V = 0.1;
    t.max_step_V = 0.1;
    const auto run = solve::trace_bias(d, t);
    REQUIRE(run.has_value());
    REQUIRE_FALSE(run->stopped.has_value());
    REQUIRE(run->points.size() >= 3);
    for (const results::BiasPoint& p : run->points) {
        REQUIRE(p.weighted_current.size() == 2);
        CHECK(std::abs(p.weighted_current[0] - p.terminal_current[0]) <=
              p.weighted_current_resolution[0] + p.terminal_current_resolution[0]);
    }
}

namespace {

// A 1D 1e17 / 1e17 diode (junction at 1 um), optionally with germanium beyond 1.5 um behind a
// thermionic interface, isothermal sinks at both ends, or 1 ns SRH lifetimes (so that a reverse
// current, generation, is far above the contact current's resolution).
device::Device model_diode(bool hetero, bool sinks, bool short_lifetime) {
    std::vector<double> x;
    for (int k = 0; k <= 400; ++k) x.push_back(2e-4 * k / 400.0);
    if (hetero) x = *mesh::straddle_interface(x, 1.5e-4, 2e-7);
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    std::vector<device::RegionId> region(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        const double xi = m.points()[i][0];
        if (xi < 1e-4) acceptors[i] = 1e17;
        else donors[i] = 1e17;
        if (hetero && xi > 1.5e-4) region[i] = 1;
    }
    physics::SemiconductorParameters si = physics::silicon_parameters;
    if (short_lifetime) si.lifetime.tau_n0 = si.lifetime.tau_p0 = 1e-9;
    std::vector<device::Region> regions{{"silicon", *physics::Semiconductor::create(si)}};
    std::vector<device::Interface> interfaces;
    if (hetero) {
        regions.push_back({"germanium", *physics::Semiconductor::create(physics::germanium_parameters)});
        interfaces.push_back({"silicon", "germanium", device::InterfaceTransport::thermionic_emission});
    }
    std::vector<device::ThermalContact> thermal;
    if (sinks) {
        thermal = {{"left", "x_min", device::ThermalContactKind::isothermal, 300.0, 0.0},
                   {"right", "x_max", device::ThermalContactKind::isothermal, 300.0, 0.0}};
    }
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = std::move(regions),
         .node_region = std::move(region),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}},
         .interfaces = std::move(interfaces),
         .thermal_contacts = std::move(thermal)});
}

// A 2D lateral diode under an oxide with an electrode on top; the oxide interface has a fixed
// charge, a trap level and a trap band, and surface recombination.
device::Device trapped_diode() {
    std::vector<double> x, y;
    for (int k = 0; k <= 40; ++k) x.push_back(2e-4 * k / 40.0);
    for (int k = 0; k <= 20; ++k) y.push_back(0.5e-4 * k / 20.0);
    y = *mesh::straddle_interface(y, 0.4e-4, 5e-7);
    mesh::Mesh m = *mesh::make_tensor_grid(x, y);
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    std::vector<device::RegionId> region(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = m.points()[i];
        if (p[1] > 0.4e-4) {
            region[i] = 1;
            continue;
        }
        if (p[0] < 1e-4) acceptors[i] = 1e17;
        else donors[i] = 1e17;
    }
    std::vector<mesh::NodeId> left, right, top;
    for (const auto k : m.find_boundary("x_min")->nodes) {
        if (region[static_cast<std::size_t>(k)] == 0) left.push_back(k);
    }
    for (const auto k : m.find_boundary("x_max")->nodes) {
        if (region[static_cast<std::size_t>(k)] == 0) right.push_back(k);
    }
    for (const auto k : m.find_boundary("y_max")->nodes) top.push_back(k);
    device::Interface oxide{.region_a = "silicon", .region_b = "oxide", .fixed_charge_cm2 = 2e10,
                            .recombination_velocity_n_cm_s = 1e3,
                            .recombination_velocity_p_cm_s = 1e3};
    oxide.traps.levels = {{physics::TrapType::acceptor, 5e10, 0.1, 1e-15, 1e-15}};
    oxide.traps.bands = {{physics::TrapType::donor, 1e11, -0.3, 0.0, 1e-15, 1e-15}};
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}, {"oxide", physics::silicon_dioxide()}},
         .node_region = std::move(region),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(left)},
                      {"cathode", device::ContactKind::ohmic, std::move(right)},
                      {"gate", device::ContactKind::electrode, std::move(top), {},
                       {.kind = device::GateElectrode::n_poly}}},
         .interfaces = {oxide}});
}

// Whether the point's weighted current equals its contact current: within the two resolutions and
// within `relative` (1e-6; at reverse bias the contact current's own noise, its two contacts'
// mismatch, is 2e-5 of it, so 1e-4 there), the contact current resolved.
void same_current(const results::BiasPoint& p, double relative = 1e-6) {
    CAPTURE(p.weighted_current, p.weighted_current_resolution, p.terminal_current,
            p.terminal_current_resolution);
    REQUIRE(std::abs(p.terminal_current[0]) > 10.0 * p.terminal_current_resolution[0]);
    const double difference = std::abs(p.weighted_current[0] - p.terminal_current[0]);
    CHECK(difference <= p.weighted_current_resolution[0] + p.terminal_current_resolution[0]);
    CHECK(difference <= relative * std::abs(p.terminal_current[0]));
}

}  // namespace

TEST_CASE("weighted current: equal to the contact current in every model it supports",
          "[weighted]") {
    // Each model keeps the total current conserved at every node in a steady state, the identity
    // the method rests on: pair generation and recombination cancel in the total.
    struct Case {
        const char* name;
        bool hetero;
        bool sinks;
        bool short_lifetime;
        double bias;
        void (*set)(solve::BiasOptions&);
    };
    const Case cases[] = {
        {"baseline", false, false, false, 0.4, [](solve::BiasOptions&) {}},
        {"Fermi-Dirac", false, false, false, 0.4,
         [](solve::BiasOptions& o) { o.models.fermi_dirac = true; }},
        {"incomplete ionization", false, false, false, 0.4,
         [](solve::BiasOptions& o) { o.models.incomplete_ionization = true; }},
        {"field mobility", false, false, false, 0.6,
         [](solve::BiasOptions& o) { o.models.field_mobility = true; }},
        {"reverse generation", false, false, true, -8.0, [](solve::BiasOptions&) {}},
        {"impact ionization", false, false, true, -8.0,
         [](solve::BiasOptions& o) { o.models.impact_ionization = true; }},
        {"local tunnelling", false, false, true, -8.0,
         [](solve::BiasOptions& o) { o.models.btbt_local = true; }},
        {"electrothermal", false, true, false, 0.6,
         [](solve::BiasOptions& o) { o.models.electrothermal = true; }},
        {"thermionic heterojunction", true, false, false, 0.4, [](solve::BiasOptions&) {}},
    };
    for (const Case& c : cases) {
        CAPTURE(c.name);
        const auto d = model_diode(c.hetero, c.sinks, c.short_lifetime);
        solve::BiasOptions options;
        options.weighted_current = true;
        options.step_control.max_cuts = 12;
        c.set(options);
        const std::vector<std::vector<double>> points{{c.bias, 0.0}};
        const auto sweep = solve::sweep_bias(d, points, options);
        REQUIRE(sweep.has_value());
        if (sweep->stopped) FAIL(sweep->stopped->message);
        same_current(sweep->points.front(), c.bias < 0.0 ? 1e-4 : 1e-6);
    }
    // Interface traps, fixed charge and surface recombination at an oxide, with an electrode.
    solve::BiasOptions options;
    options.weighted_current = true;
    const std::vector<std::vector<double>> points{{0.0, 0.0, 0.0}, {0.4, 0.0, 0.0}};
    const auto sweep = solve::sweep_bias(trapped_diode(), points, options);
    REQUIRE(sweep.has_value());
    REQUIRE_FALSE(sweep->stopped.has_value());
    CHECK(sweep->points.back().weighted_current[2] == 0.0);
    same_current(sweep->points.back());
}
