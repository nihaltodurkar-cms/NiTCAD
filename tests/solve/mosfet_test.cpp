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
#include "NiTCAD/mesh/tensor_grid.hpp"
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
