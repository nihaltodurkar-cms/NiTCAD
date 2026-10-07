// 2D n-channel MOSFET through drift-diffusion with a gate contact (ARCHITECTURE.md section 11,
// Unit 12): the legacy tests/test_validation_2d.py MOSFET gates.
//
// Fixture (legacy _build_test_mosfet and mosfet.build_mosfet): Lg = 600 nm, Lsd = 300 nm, depth
// 200 nm; p-type body N_A = 1e17; n+ source and drain, 1e19 peak, Gaussian in depth (sigma 50 nm)
// times an erfc lateral roll-off (sigma 10 nm) at the gate edges; t_ox = 5 nm, n+ polysilicon gate;
// graded_mesh in x (focus at both gate edges) and y (focus at the surface) with nx = 90, ny = 50;
// source, gate and drain on the surface, the body contact on the bottom; the legacy default models
// (doping mobility, SRH, Auger, band-gap narrowing). Id-Vg at Vds = 0.05 V, Vg from -1 to 1 V.
//
// About 9500 nodes: each Eigen SparseLU factorization takes about 0.27 s in Release, and the sweep
// about 50 s. The tests carry the hidden tag [.mosfet] and run as their own ctest entry in the
// Release configuration only (CMakeLists.txt); in Debug they would take tens of minutes.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "legacy_graded_mesh.hpp"
#include "legacy_moscap.hpp"

using namespace NiTCAD;

namespace {

constexpr double Lg = 6e-5, Lsd = 3e-5, depth = 2e-5, N_A = 1e17, t_ox = 5e-7;
constexpr double N_peak = 1e19, sigma_y = 5e-6, sigma_lat = 1e-6;

// mosfet._sd_profile: Npeak exp(-y^2 / 2 sigma_y^2) 0.5 erfc(-s / (sqrt(2) sigma_lat)), with
// s = x_edge - x on the source side and x - x_edge on the drain side.
double sd_profile(double x, double y, double x_edge, bool source) {
    const double s = source ? x_edge - x : x - x_edge;
    return N_peak * std::exp(-(y * y) / (2.0 * sigma_y * sigma_y)) * 0.5 *
           std::erfc(-s / (std::sqrt(2.0) * sigma_lat));
}

// Contacts in the order source, drain, body, gate (the legacy add order).
device::Device mosfet(int nx = 90, int ny = 50) {
    const double L = 2 * Lsd + Lg;
    const auto x = legacy_graded_mesh(L, {Lsd, Lsd + Lg}, L / (nx * 20), L / nx, 1.15);
    const auto y = legacy_graded_mesh(depth, 0.0, depth / (ny * 20), depth / ny, 1.15);
    mesh::Mesh m = *mesh::make_tensor_grid(x, y);
    const std::size_t n = m.node_count();
    std::vector<double> donors(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = m.points()[i];
        donors[i] = sd_profile(p[0], p[1], Lsd, true) + sd_profile(p[0], p[1], Lsd + Lg, false);
    }
    std::vector<mesh::NodeId> source, drain, gate;
    for (const mesh::NodeId v : m.find_boundary("y_min")->nodes) {
        const double px = m.points()[static_cast<std::size_t>(v)][0];
        if (px <= Lsd) {
            source.push_back(v);
        } else if (px >= Lsd + Lg) {
            drain.push_back(v);
        } else {
            gate.push_back(v);
        }
    }
    auto body = m.find_boundary("y_max")->nodes;
    device::Contact g{"gate", device::ContactKind::gate, std::move(gate)};
    g.gate = {.boundary = "y_min", .oxide_thickness_cm = t_ox};
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::vector<double>(n, N_A),
         .contacts = {{"source", device::ContactKind::ohmic, std::move(source)},
                      {"drain", device::ContactKind::ohmic, std::move(drain)},
                      {"body", device::ContactKind::ohmic, std::move(body)},
                      std::move(g)}});
}

// The same MOSFET with a meshed oxide (Unit 15b): SiO2 rows on [-t_ox, -h/2] above the whole
// surface (h the first spacing of the legacy y axis, the silicon nodes moved down by h/2, the
// interface at y = 0 midway between rows), the gate an electrode on the oxide top over the channel
// (Lsd < x < Lsd + Lg), the rest of the oxide top free. Source and drain stay ohmic on the top
// silicon row (a boundary patch "surface" added to the grid), as the lumped device's on y = 0.
device::Device meshed_mosfet(int nx = 90, int ny = 50, int oxide_cells = 5) {
    const double L = 2 * Lsd + Lg;
    const auto x = legacy_graded_mesh(L, {Lsd, Lsd + Lg}, L / (nx * 20), L / nx, 1.15);
    const auto y_si = legacy_graded_mesh(depth, 0.0, depth / (ny * 20), depth / ny, 1.15);
    const double h = y_si[1];
    std::vector<double> y;
    for (int k = 0; k <= oxide_cells; ++k) y.push_back(-t_ox + (t_ox - 0.5 * h) * k / oxide_cells);
    for (const double v : y_si) y.push_back(v + 0.5 * h);
    const mesh::Mesh grid = *mesh::make_tensor_grid(x, y);
    // The top silicon row as a patch, each node with its x width.
    mesh::BoundaryPatch surface{"surface", {}, {}};
    for (std::size_t i = 0; i < grid.node_count(); ++i) {
        const auto& p = grid.points()[i];
        if (p[1] != 0.5 * h) continue;
        const auto k =
            static_cast<std::size_t>(std::lower_bound(x.begin(), x.end(), p[0]) - x.begin());
        const double left = k > 0 ? 0.5 * (x[k] - x[k - 1]) : 0.0;
        const double right = k + 1 < x.size() ? 0.5 * (x[k + 1] - x[k]) : 0.0;
        surface.nodes.push_back(static_cast<mesh::NodeId>(i));
        surface.areas.push_back(left + right);
    }
    std::vector<mesh::BoundaryPatch> patches(grid.boundary().begin(), grid.boundary().end());
    patches.push_back(surface);
    mesh::Mesh m = *mesh::Mesh::from_parts(
        2, std::vector<mesh::Point>(grid.points().begin(), grid.points().end()),
        std::vector<double>(grid.volumes().begin(), grid.volumes().end()),
        std::vector<mesh::Edge>(grid.edges().begin(), grid.edges().end()), std::move(patches));
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    std::vector<device::RegionId> region(n, 1);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = m.points()[i];
        if (p[1] < 0.0) {
            region[i] = 0;
            continue;
        }
        const double depth_y = p[1] - 0.5 * h;  // the legacy profile's y
        donors[i] =
            sd_profile(p[0], depth_y, Lsd, true) + sd_profile(p[0], depth_y, Lsd + Lg, false);
        acceptors[i] = N_A;
    }
    std::vector<mesh::NodeId> source, drain, gate;
    for (const mesh::NodeId v : m.find_boundary("surface")->nodes) {
        const double px = m.points()[static_cast<std::size_t>(v)][0];
        if (px <= Lsd) source.push_back(v);
        if (px >= Lsd + Lg) drain.push_back(v);
    }
    for (const mesh::NodeId v : m.find_boundary("y_min")->nodes) {
        const double px = m.points()[static_cast<std::size_t>(v)][0];
        if (px > Lsd && px < Lsd + Lg) gate.push_back(v);
    }
    auto body = m.find_boundary("y_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
         .node_region = std::move(region),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"source", device::ContactKind::ohmic, std::move(source)},
                      {"drain", device::ContactKind::ohmic, std::move(drain)},
                      {"body", device::ContactKind::ohmic, std::move(body)},
                      {"gate", device::ContactKind::electrode, std::move(gate)}},
         .interfaces = {{"oxide", "silicon"}}});
}

// mosfet.id_vg_sweep: equilibrium, then each Vg at the drain bias from the previous point. The
// drain current in A/cm (per unit depth).
std::vector<double> id_vg(const device::Device& d, const std::vector<double>& Vg, double Vds) {
    std::vector<std::vector<double>> points;
    for (const double v : Vg) points.push_back({0.0, Vds, 0.0, v});
    const auto sweep = solve::sweep_bias(d, points);
    REQUIRE(sweep.has_value());
    if (sweep->stopped) {
        CAPTURE(sweep->stopped->message, sweep->points.size());
        FAIL("the Id-Vg sweep stopped");
    }
    std::vector<double> Id;
    for (const auto& p : sweep->points) {
        // No current through the gate, and Kirchhoff down to the rounding floor of the extracted
        // currents (measured 5e-12 A/cm here: majority-flux cancellation in the 1e19 source and
        // drain, ARCHITECTURE.md 6.2, Unit 11).
        REQUIRE(p.terminal_current[3] == 0.0);
        const double sum = p.terminal_current[0] + p.terminal_current[1] + p.terminal_current[2];
        CAPTURE(p.bias_V[3], sum);
        REQUIRE(std::abs(sum) <= 1e-6 * std::abs(p.terminal_current[1]) + 1e-10);
        Id.push_back(p.terminal_current[1]);
    }
    return Id;
}

std::vector<double> linspace(double a, double b, int n) {
    std::vector<double> v(static_cast<std::size_t>(n));
    for (int k = 0; k < n; ++k) v[static_cast<std::size_t>(k)] = a + (b - a) * k / (n - 1);
    return v;
}

// Legacy _extract_vth_max_gm: the tangent at the largest gm, Vth = Vg - Id / gm - Vds / 2.
double vth_max_gm(const std::vector<double>& Vg, const std::vector<double>& Id, double Vds) {
    const auto gm = numpy_gradient(Id, Vg);
    const auto i =
        static_cast<std::size_t>(std::max_element(gm.begin(), gm.end()) - gm.begin());
    return Vg[i] - Id[i] / gm[i] - Vds / 2.0;
}

}  // namespace

TEST_CASE("mosfet: Id-Vg threshold, on/off ratio, swing and monotonic rise (legacy gates)",
          "[.mosfet]") {
    // One sweep serves the legacy tests test_mosfet_vth_matches_moscap_landmark (Vg -1 to 1 V in
    // 0.1 V steps) and test_id_vg_monotonic_above_threshold (0.3 to 1.5 V in 0.1 V steps).
    const auto Vg = linspace(-1.0, 1.5, 26);
    const double Vds = 0.05;
    const auto Id = id_vg(mosfet(), Vg, Vds);
    const std::vector<double> Vg_legacy(Vg.begin(), Vg.begin() + 21);
    const std::vector<double> Id_legacy(Id.begin(), Id.begin() + 21);

    const auto lm = LegacyMOSCapacitor(-N_A, t_ox).analytic_landmarks();
    const double extracted = vth_max_gm(Vg_legacy, Id_legacy, Vds);
    double swing_slope = 0.0;  // the steepest decades per volt
    for (std::size_t k = 0; k + 1 < Id_legacy.size(); ++k) {
        swing_slope = std::max(swing_slope, (std::log10(std::abs(Id[k + 1])) -
                                             std::log10(std::abs(Id[k]))) /
                                                (Vg[k + 1] - Vg[k]));
    }
    const double SS = 1000.0 / swing_slope;  // mV/decade
    CAPTURE(Id, extracted, lm.V_th, SS);
    REQUIRE(Id_legacy.back() > 1e6 * std::abs(Id_legacy.front()));  // a real off state
    REQUIRE(-0.5 < lm.V_th);
    REQUIRE(lm.V_th < 1.5);
    REQUIRE(std::abs(extracted - lm.V_th) < 0.1);
    REQUIRE(SS > 55.0);  // the thermal limit is V_T ln 10 = 59.5 mV/decade
    REQUIRE(SS < 120.0);
    // Tighter than the legacy: the swing of a long channel is V_T ln 10 (1 + C_dep / C_ox), with
    // the depletion capacitance between its values at phi_s = 2 phi_F and phi_F (depletion
    // approximation).
    const LegacyMOSCapacitor m(-N_A, t_ox);
    const auto ideal = [&](double phi_s) {
        const double W = std::sqrt(2.0 * m.eps_s * phi_s / (base::q_C * N_A));
        return 1000.0 * m.VT * std::log(10.0) * (1.0 + m.eps_s / W / m.Cox);
    };
    CAPTURE(ideal(2.0 * lm.phi_F), ideal(lm.phi_F));
    REQUIRE(SS > ideal(2.0 * lm.phi_F));
    REQUIRE(SS < ideal(lm.phi_F));
    // Monotonic above threshold, 0.3 to 1.5 V.
    for (std::size_t k = 13; k + 1 < Id.size(); ++k) REQUIRE(Id[k + 1] > Id[k]);
}

TEST_CASE("mosfet: the on current is mesh independent (legacy gate)", "[.mosfet]") {
    // Legacy test_mosfet_mesh_independence: Id(Vg = 1 V, Vds = 0.05 V) on nx = 90, ny = 50 and on
    // nx = 180, ny = 100 within 10%.
    const std::vector<double> point{0.0, 0.05, 0.0, 1.0};
    const auto coarse = solve::solve_bias(mosfet(), point);
    const auto fine = solve::solve_bias(mosfet(180, 100), point);
    REQUIRE(coarse.has_value());
    REQUIRE(fine.has_value());
    const double a = coarse->terminal_current[1], b = fine->terminal_current[1];
    CAPTURE(a, b, (a - b) / b);
    REQUIRE(std::abs(a - b) / std::abs(b) < 0.10);
}

TEST_CASE("mosfet: field mobility lowers the drain current at high drain bias (Unit 13)",
          "[.mosfet]") {
    // The legacy fixture at V_G = 1 V, drain swept to 1 V in 0.1 V steps. At V_D = 0 no current
    // flows, so that state does not depend on the mobility, and both sweeps start from the one
    // constant-mobility solve there. With field mobility the first drain step diverges under the
    // legacy damping cap of 5 V_T and converges under 1 V_T (ARCHITECTURE.md 6.2, Unit 13), so
    // that sweep uses max_update = 1. At low drain bias the lateral field is small and field
    // mobility changes little; towards 1 V the field near the drain reaches 1e4-1e5 V/cm and the
    // current falls below the constant-mobility one.
    const device::Device d = mosfet();
    const std::vector<double> gate_on{0.0, 0.0, 0.0, 1.0};
    const auto start = solve::solve_bias(d, gate_on);
    REQUIRE(start.has_value());
    std::vector<std::vector<double>> points;
    for (int k = 1; k <= 10; ++k) points.push_back({0.0, 0.1 * k, 0.0, 1.0});
    solve::BiasOptions field;
    field.models.field_mobility = true;
    field.newton.max_update = 1.0;
    const auto on = solve::sweep_bias(d, points, field, &start->fields);
    const auto off = solve::sweep_bias(d, points, {}, &start->fields);
    REQUIRE(on.has_value());
    REQUIRE(off.has_value());
    REQUIRE_FALSE(on->stopped.has_value());
    REQUIRE_FALSE(off->stopped.has_value());
    std::vector<double> ratio;
    std::vector<std::size_t> iterations;
    for (std::size_t k = 0; k < points.size(); ++k) {
        ratio.push_back(on->points[k].terminal_current[1] / off->points[k].terminal_current[1]);
        iterations.push_back(on->points[k].convergence.iterations.size());
    }
    CAPTURE(ratio, iterations, on->points.back().terminal_current[1],
            off->points.back().terminal_current[1]);
    REQUIRE(ratio.front() > 0.95);
    REQUIRE(ratio.front() < 1.0);
    for (std::size_t k = 1; k < ratio.size(); ++k) REQUIRE(ratio[k] < ratio[k - 1]);
    REQUIRE(ratio.back() < 0.9);
}

TEST_CASE("mosfet: a meshed oxide gives the lumped-oxide transfer curve (Unit 15b)", "[.mosfet]") {
    // Id-Vg at Vds = 0.05 V from -1 to 1.5 V on both devices. They differ by the oxide's fringing
    // field at the gate edges and over source and drain (the lumped oxide couples each surface
    // node vertically only) and by the meshed silicon starting h / 2 below the interface.
    const auto Vg = linspace(-1.0, 1.5, 26);
    const double Vds = 0.05;
    const auto lumped = id_vg(mosfet(), Vg, Vds);
    const auto meshed = id_vg(meshed_mosfet(), Vg, Vds);
    const std::vector<double> Vg21(Vg.begin(), Vg.begin() + 21);
    const double vth_l =
        vth_max_gm(Vg21, std::vector<double>(lumped.begin(), lumped.begin() + 21), Vds);
    const double vth_m =
        vth_max_gm(Vg21, std::vector<double>(meshed.begin(), meshed.begin() + 21), Vds);
    std::vector<double> ratio;
    for (std::size_t k = 0; k < Vg.size(); ++k) ratio.push_back(meshed[k] / lumped[k]);
    CAPTURE(lumped, meshed, ratio, vth_l, vth_m);
    // Measured: the thresholds 0.152 V apart by 0.21 mV; from Vg = -0.4 V (Id 2e-10 A/cm) up the
    // currents within 0.58% (the meshed one lower, most at 1.5 V). Below -0.4 V both are at the
    // extraction floor of about 5e-12 A/cm (sign included), so they are not compared.
    REQUIRE(std::abs(vth_m - vth_l) < 2e-3);
    for (std::size_t k = 6; k < Vg.size(); ++k) REQUIRE(std::abs(ratio[k] - 1.0) < 0.01);
}
