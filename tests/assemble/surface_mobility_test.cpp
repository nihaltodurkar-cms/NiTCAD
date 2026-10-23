// Lombardi surface mobility in the drift-diffusion assembler (Unit 46; DECISIONS.md LD1-LD8,
// LS1, LS3, LS4, LS10): E_perp of a uniform field at 0, 30, 45 and 60 degrees on 2D tensor and
// Delaunay meshes and a 3D tensor mesh, read from the edge currents' ratio to the bulk model's;
// the finite-difference Jacobian gate, isothermal and electrothermal, with and without Canali;
// the pattern unchanged with the model off; and the refusals.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numbers>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/physics/surface_mobility.hpp"

using namespace NiTCAD;
using assemble::DriftDiffusion;

namespace {

constexpr double L = 1e-4;  // [cm]

std::vector<double> uniform(double a, double b, int nodes) {
    std::vector<double> v(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) v[static_cast<std::size_t>(i)] = a + (b - a) * i / (nodes - 1);
    return v;
}

mesh::Mesh tensor(int D, int nodes) {
    const auto a = uniform(0.0, L, nodes);
    return D == 2 ? *mesh::make_tensor_grid(a, a) : *mesh::make_tensor_grid(a, a, a);
}

mesh::Mesh delaunay() {
    mesh::Geometry2D g{{{{{0, 0, 0}, {L, 0, 0}, {L, L, 0}, {0, L, 0}}, 0}},
                       {{"x_min", {0, 0, 0}, {0, L, 0}}, {"x_max", {L, 0, 0}, {L, L, 0}}}};
    mesh::MeshOptions2D o;
    o.max_size = L / 6.0;
    auto m = mesh::generate_mesh(g, o);
    REQUIRE(m.has_value());
    return std::move(m->mesh);
}

struct Block {
    physics::SemiconductorParameters material = physics::silicon_parameters;
    bool sinks = false;  // the electrothermal model's thermal contacts
};

// A uniformly doped silicon block (1e17 donors, 1e16 acceptors), ohmic contacts on x_min and x_max.
device::Device block(mesh::Mesh m, const Block& o = {}) {
    const std::size_t n = m.node_count();
    auto left = m.find_boundary("x_min")->nodes;
    auto right = m.find_boundary("x_max")->nodes;
    std::vector<device::ThermalContact> thermal;
    if (o.sinks) {
        thermal = {{"sink", "x_min", device::ThermalContactKind::isothermal, 310.0, 0.0},
                   {"package", "x_max", device::ThermalContactKind::resistance, 290.0, 1e-4}};
    }
    auto d = device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", *physics::Semiconductor::create(o.material)}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::vector<double>(n, 1e17),
         .acceptors = std::vector<double>(n, 1e16),
         .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                      {"right", device::ContactKind::ohmic, std::move(right)}},
         .thermal_contacts = std::move(thermal)});
    if (!d) FAIL(d.error().message);
    return std::move(*d);
}

assemble::PhysicsModels lombardi(bool on = true) {
    assemble::PhysicsModels m;
    m.surface_mobility = on ? assemble::SurfaceMobility::lombardi : assemble::SurfaceMobility::off;
    return m;
}

struct Noise {
    std::uint64_t state;
    double next() {
        std::uint64_t z = (state += 0x9E3779B97F4A7C15ull);
        z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
        z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
        z ^= z >> 31;
        return static_cast<double>(z >> 11) * 0x1p-53 * 2.0 - 1.0;
    }
};

// A field E0 [V/cm] at angle theta to x in the x-y plane (and, in 3D, tilted by phi out of it),
// n 1e17 and p 1e16 cm^-3, relative noise on every unknown; with four unknowns per node theta
// rises from 1 to 1.2 along x.
std::vector<double> sloped(const DriftDiffusion& dd, const device::Device& d,
                           const assemble::Scaling& s, double E0, double theta, double phi,
                           double noise, std::uint64_t seed) {
    const std::size_t m = dd.stride();
    std::vector<double> x(dd.unknowns());
    Noise z{seed};
    const double e[3] = {std::cos(phi) * std::cos(theta), std::cos(phi) * std::sin(theta),
                         std::sin(phi)};
    for (std::size_t i = 0; i < dd.node_count(); ++i) {
        const auto& p = d.mesh().points()[i];
        const double along = p[0] * e[0] + p[1] * e[1] + p[2] * e[2];
        x[m * i] = -E0 * along / s.V_T * (1.0 + noise * z.next());
        x[m * i + 1] = 1e17 / s.Ns * (1.0 + noise * z.next());
        x[m * i + 2] = 1e16 / s.Ns * (1.0 + noise * z.next());
        if (m == 4) x[4 * i + 3] = 0.2 * p[0] / L * (1.0 + noise * z.next());
    }
    return x;
}

std::vector<std::vector<double>> dense(const linalg::SparseMatrix& a) {
    const auto n = static_cast<std::size_t>(a.rows());
    std::vector<std::vector<double>> d(n, std::vector<double>(n, 0.0));
    for (std::size_t r = 0; r < n; ++r) {
        for (auto k = static_cast<std::size_t>(a.row_offsets()[r]);
             k < static_cast<std::size_t>(a.row_offsets()[r + 1]); ++k) {
            d[r][static_cast<std::size_t>(a.col_indices()[k])] = a.values()[k];
        }
    }
    return d;
}

// The finite-difference gate (as impact_ionization_test.cpp): every column, central differences,
// the worst error over the column's largest entry, each column's the smallest at steps of 1e-5,
// 1e-6 and 1e-7 max(|x|, 1).
double fd_jacobian_error(const DriftDiffusion& dd, std::vector<double> x) {
    const std::size_t n = dd.unknowns();
    std::vector<double> f(n), fp(n), fm(n);
    auto jacobian = dd.make_jacobian();
    dd.evaluate(x, f, jacobian);
    const auto j = dense(jacobian);
    double worst = 0.0;
    for (std::size_t c = 0; c < n; ++c) {
        double scale = 1e-300;
        for (std::size_t r = 0; r < n; ++r) scale = std::max(scale, std::abs(j[r][c]));
        double best = std::numeric_limits<double>::infinity();
        for (const double relative : {1e-5, 1e-6, 1e-7}) {
            const double base = x[c];
            const double step = relative * std::max(std::abs(base), 1.0);
            x[c] = base + step;
            const double up = x[c];
            dd.residual(x, fp);
            x[c] = base - step;
            const double down = x[c];
            dd.residual(x, fm);
            x[c] = base;
            double error = 0.0;
            for (std::size_t r = 0; r < n; ++r) {
                error = std::max(error, std::abs((fp[r] - fm[r]) / (up - down) - j[r][c]));
            }
            best = std::min(best, error / scale);
        }
        worst = std::max(worst, best);
    }
    return worst;
}

}  // namespace

TEST_CASE("surface mobility assemble: E_perp of a uniform field at any angle (LS3)") {
    // Uniform doping: both ends' bulk mobilities are equal, so each edge's current with Lombardi
    // over its current without is mu_L(E_perp) / mu_b, E_perp = sqrt(|E|^2 - (E . t)^2 + E_0^2).
    const physics::Semiconductor si = *physics::Semiconductor::create(physics::silicon_parameters);
    const double N = 1.1e17;
    const double mu_n = physics::caughey_thomas_mobility(si, physics::Carrier::electron, N, 300.0);
    const double mu_p = physics::masetti_mobility(physics::masetti_boron_silicon, N, 300.0);
    struct Case {
        mesh::Mesh mesh;
        double phi;
        const char* name;
    };
    std::vector<Case> cases;
    cases.push_back({tensor(2, 7), 0.0, "2D tensor"});
    cases.push_back({delaunay(), 0.0, "2D Delaunay"});
    cases.push_back({tensor(3, 5), 0.35, "3D tensor"});
    const double E0 = 2e5;
    for (Case& c : cases) {
        const device::Device d = block(std::move(c.mesh));
        const auto s = *assemble::make_scaling(d);
        const DriftDiffusion off = *DriftDiffusion::create(d, s, lombardi(false));
        const DriftDiffusion on = *DriftDiffusion::create(d, s, lombardi(true));
        for (const double degrees : {0.0, 30.0, 45.0, 60.0}) {
            const double theta = degrees * std::numbers::pi / 180.0;
            const auto x = sloped(off, d, s, E0, theta, c.phi, 0.0, 1);
            const auto j_off = off.edge_currents(x);
            const auto j_on = on.edge_currents(x);
            const double e[3] = {E0 * std::cos(c.phi) * std::cos(theta),
                                 E0 * std::cos(c.phi) * std::sin(theta), E0 * std::sin(c.phi)};
            double worst = 0.0;
            std::size_t checked = 0;
            for (std::size_t k = 0; k < d.mesh().edges().size(); ++k) {
                const mesh::Edge& edge = d.mesh().edges()[k];
                const auto& a = d.mesh().points()[static_cast<std::size_t>(edge.first)];
                const auto& b = d.mesh().points()[static_cast<std::size_t>(edge.second)];
                double along = 0.0;
                for (int q = 0; q < 3; ++q) along += e[q] * (b[q] - a[q]) / edge.length;
                const double perp = std::sqrt(E0 * E0 - along * along + 1.0);
                const double rn =
                    physics::lombardi_mobility(physics::lombardi_silicon_electrons, mu_n, N,
                                               300.0, perp).mobility / mu_n;
                const double rp = physics::lombardi_mobility(physics::lombardi_silicon_holes, mu_p,
                                                             N, 300.0, perp).mobility / mu_p;
                for (const auto [on_flux, off_flux, ratio] :
                     {std::tuple{j_on[k].first, j_off[k].first, rn},
                      std::tuple{j_on[k].second, j_off[k].second, rp}}) {
                    if (std::abs(along) < 1e-3 * E0) continue;  // no drift along the edge
                    worst = std::max(worst, std::abs(on_flux / off_flux / ratio - 1.0));
                    ++checked;
                }
            }
            CAPTURE(c.name, degrees, worst, checked);
            REQUIRE(checked > 0);
            REQUIRE(worst <= 1e-10);
        }
    }
}

TEST_CASE("surface mobility assemble: the finite-difference Jacobian gate (LS4)") {
    struct Case {
        int D;  // 2 or 3; 0: 2D Delaunay
        bool field, thermal;
        const char* name;
    };
    const Case cases[] = {{2, false, false, "2D"},
                          {2, true, false, "2D Canali"},
                          {0, false, false, "2D Delaunay"},
                          {0, true, false, "2D Delaunay Canali"},
                          {3, false, false, "3D"},
                          {3, true, false, "3D Canali"},
                          {2, false, true, "2D electrothermal"},
                          {2, true, true, "2D electrothermal Canali"},
                          {3, true, true, "3D electrothermal Canali"}};
    for (const Case& c : cases) {
        mesh::Mesh m = c.D == 0 ? delaunay() : tensor(c.D, c.D == 2 ? 6 : 4);
        const device::Device d = block(std::move(m), {.sinks = c.thermal});
        const auto s = *assemble::make_scaling(d);
        assemble::PhysicsModels models = lombardi();
        models.field_mobility = c.field;
        models.electrothermal = c.thermal;
        auto made = DriftDiffusion::create(d, s, models);
        if (!made) FAIL(c.name << ": " << made.error().message);
        const DriftDiffusion& dd = *made;
        const auto x = sloped(dd, d, s, 3e5, 0.5, c.D == 3 ? 0.3 : 0.0, 0.05, 7);
        const double error = fd_jacobian_error(dd, x);
        CAPTURE(c.name, error);
        REQUIRE(error <= 1e-6);
    }
}

// The small-signal current rows at s = 0 against central differences of terminal_currents in every
// unknown (LS10's small signal, exact on edges that leave a contact at an angle to the field).
TEST_CASE("surface mobility assemble: the current rows are the linearized terminal currents (LS10)") {
    struct Case {
        int D;
        bool field, thermal;
        const char* name;
    };
    const Case cases[] = {{2, false, false, "2D"},
                          {3, true, false, "3D Canali"},
                          {2, true, true, "2D electrothermal Canali"}};
    for (const Case& c : cases) {
        const device::Device d = block(tensor(c.D, c.D == 2 ? 6 : 4), {.sinks = c.thermal});
        const auto s = *assemble::make_scaling(d);
        assemble::PhysicsModels models = lombardi();
        models.field_mobility = c.field;
        models.electrothermal = c.thermal;
        const DriftDiffusion dd = *DriftDiffusion::create(d, s, models);
        auto x = sloped(dd, d, s, 3e5, 0.5, c.D == 3 ? 0.3 : 0.0, 0.05, 11);
        const auto rows = dd.small_signal_currents(x, 0.0);
        REQUIRE(rows.size() == 2);
        std::vector<std::vector<double>> analytic(2, std::vector<double>(dd.unknowns(), 0.0));
        for (std::size_t k = 0; k < 2; ++k) {
            for (std::size_t q = 0; q < rows[k].columns.size(); ++q) {
                analytic[k][rows[k].columns[q]] += rows[k].values[q].real();
            }
        }
        double worst = 0.0, largest = 0.0;
        for (std::size_t col = 0; col < dd.unknowns(); ++col) {
            const double base = x[col];
            const double h = 1e-7 * std::max(std::abs(base), 1.0);
            x[col] = base + h;
            const double hi = x[col];
            const auto up = dd.terminal_currents(x);
            x[col] = base - h;
            const double lo = x[col];
            const auto down = dd.terminal_currents(x);
            x[col] = base;
            for (std::size_t k = 0; k < 2; ++k) {
                largest = std::max(largest, std::abs(analytic[k][col]));
                worst = std::max(worst, std::abs((up[k] - down[k]) / (hi - lo) - analytic[k][col]));
            }
        }
        CAPTURE(c.name, worst, largest);
        REQUIRE(largest > 0.0);
        REQUIRE(worst <= 1e-6 * largest);
    }
}

TEST_CASE("surface mobility assemble: off, the pattern and residual are unchanged (LS1)") {
    const device::Device d = block(tensor(2, 6));
    const auto s = *assemble::make_scaling(d);
    const DriftDiffusion base = *DriftDiffusion::create(d, s, assemble::PhysicsModels{});
    const DriftDiffusion off = *DriftDiffusion::create(d, s, lombardi(false));
    const DriftDiffusion on = *DriftDiffusion::create(d, s, lombardi(true));
    REQUIRE(off.make_jacobian().has_same_pattern(base.make_jacobian()));
    REQUIRE(on.make_jacobian().nonzeros() > base.make_jacobian().nonzeros());
    const auto x = sloped(base, d, s, 1e5, 0.4, 0.0, 0.05, 3);
    std::vector<double> f0(base.unknowns()), f1(base.unknowns());
    base.residual(x, f0);
    off.residual(x, f1);
    REQUIRE(f0 == f1);
}

TEST_CASE("surface mobility assemble: refusals (LD6, LD8; LS10)") {
    const auto ge = physics::germanium_parameters;
    const device::Device g = block(tensor(2, 4), {.material = ge});
    const auto sg = *assemble::make_scaling(g);
    const auto refused = DriftDiffusion::create(g, sg, lombardi());
    REQUIRE_FALSE(refused.has_value());
    REQUIRE(refused.error().code == base::ErrorCode::invalid_input);
    // Silicon with another parameter changed is still silicon's mobility.
    physics::SemiconductorParameters other = physics::silicon_parameters;
    other.eps_r = 11.8;
    const device::Device o = block(tensor(2, 4), {.material = other});
    REQUIRE(DriftDiffusion::create(o, *assemble::make_scaling(o), lombardi()).has_value());
    // Impact ionization reads the fluxes' partials, which reach past its pattern.
    assemble::PhysicsModels impact = lombardi();
    impact.impact_ionization = true;
    const auto both = DriftDiffusion::create(o, *assemble::make_scaling(o), impact);
    REQUIRE_FALSE(both.has_value());
    REQUIRE(both.error().code == base::ErrorCode::invalid_input);
}
