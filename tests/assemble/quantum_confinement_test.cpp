// MLDA quantum confinement in the assemblers (Unit 47; DECISIONS.md QD1-QD7, QC3, QC4, QC6,
// QC11): the distance to the nearest interface on tensor, corner, Delaunay and lumped-gate devices;
// the densities of states times the factors with the band edges and the trap levels in place;
// the finite-difference Jacobian gate with the model on; and the refusals.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/band_edges.hpp"
#include "NiTCAD/assemble/charges.hpp"
#include "NiTCAD/assemble/equilibrium.hpp"
#include "NiTCAD/assemble/quantum_confinement.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/assemble/system.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/linalg/linear_solver.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/interface_traps.hpp"
#include "NiTCAD/physics/quantum_confinement.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/results/convergence.hpp"
#include "NiTCAD/solve/newton.hpp"
#include "scaled_device.hpp"

using namespace NiTCAD;
using assemble::System;
using assemble::create_equilibrium;

namespace {

std::vector<double> uniform(double a, double b, int nodes) {
    std::vector<double> v(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) v[static_cast<std::size_t>(i)] = a + (b - a) * i / (nodes - 1);
    return v;
}

assemble::PhysicsModels mlda(bool fermi_dirac = false) {
    assemble::PhysicsModels m;
    m.fermi_dirac = fermi_dirac;
    m.quantum_confinement = assemble::QuantumConfinement::mlda;
    return m;
}

// As insulator_test.cpp's: oxide on [-5 nm, 0) with the electrode on x_min, p-silicon (1e17) on
// (0, 30 nm] with an n+ pocket in 2D and 3D, the substrate contact on x_max; the interface at
// x = 0 between nodes at -+0.25 nm, nodes 0.5 nm apart near it, so the first nodes are well inside
// lambda. Traps, fixed charge and surface recombination on the interface.
device::Device mos(int dimension, bool traps = true, bool pocket = true) {
    auto axis = uniform(-5e-7, 3e-6, 36);
    for (double v = 0.5e-7; v < 3e-7; v += 0.5e-7) axis.push_back(v);
    std::sort(axis.begin(), axis.end());
    axis.erase(std::unique(axis.begin(), axis.end(),
                           [](double a, double b) { return std::abs(a - b) < 1e-9; }),
               axis.end());
    const auto x = *mesh::straddle_interface(axis, 0.0, 0.5e-7);
    const auto t = uniform(0.0, 3e-6, 4);
    mesh::Mesh m = dimension == 1   ? *mesh::make_tensor_grid(x)
                   : dimension == 2 ? *mesh::make_tensor_grid(x, t)
                                    : *mesh::make_tensor_grid(x, t, t);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n);
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = m.points()[i];
        const bool oxide = p[0] < 0.0;
        region[i] = oxide ? 0 : 1;
        if (oxide) continue;
        acceptors[i] = 1e17;
        if (pocket && dimension > 1 && p[1] > 2e-6 && p[0] > 2e-6) donors[i] = 1e19;
    }
    std::vector<mesh::NodeId> electrode;
    for (const mesh::NodeId v : m.find_boundary("x_min")->nodes) {
        if (dimension == 1 || m.points()[static_cast<std::size_t>(v)][1] < 1.5e-6) {
            electrode.push_back(v);
        }
    }
    auto substrate = m.find_boundary("x_max")->nodes;
    device::Interface f{"oxide", "silicon"};
    f.fixed_charge_cm2 = 3e11;
    if (traps) {
        f.traps.levels = {{.type = physics::TrapType::donor, .density_cm2 = 2e11, .energy_eV = -0.2,
                           .sigma_n_cm2 = 1e-15, .sigma_p_cm2 = 3e-16},
                          {.type = physics::TrapType::acceptor, .density_cm2 = 1e11,
                           .energy_eV = 0.25, .sigma_n_cm2 = 2e-16, .sigma_p_cm2 = 1e-15}};
        f.recombination_velocity_n_cm_s = 1e3;
        f.recombination_velocity_p_cm_s = 4e2;
    }
    device::Contact gate{"gate", device::ContactKind::electrode, std::move(electrode)};
    gate.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.4};
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
         .node_region = std::move(region),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {std::move(gate),
                      {"substrate", device::ContactKind::ohmic, std::move(substrate)}},
         .interfaces = {std::move(f)}});
}

// A p-silicon (1e17) square [0, 20 nm]^2 under oxide on two sides (x < 0 and y < 0, 2 nm thick),
// the substrate contact on x_max.
device::Device corner() {
    auto axis = uniform(-2e-7, 2e-6, 23);
    const auto a = *mesh::straddle_interface(axis, 0.0, 0.5e-7);
    mesh::Mesh m = *mesh::make_tensor_grid(a, a);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n);
    std::vector<double> acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = m.points()[i];
        const bool oxide = p[0] < 0.0 || p[1] < 0.0;
        region[i] = oxide ? 0 : 1;
        if (!oxide) acceptors[i] = 1e17;
    }
    std::vector<mesh::NodeId> substrate;
    for (const mesh::NodeId v : m.find_boundary("x_max")->nodes) {
        if (m.points()[static_cast<std::size_t>(v)][1] > 0.0) substrate.push_back(v);
    }
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
         .node_region = std::move(region),
         .donors = std::vector<double>(n, 0.0),
         .acceptors = std::move(acceptors),
         .contacts = {{"substrate", device::ContactKind::ohmic, std::move(substrate)}}});
}

// A Delaunay MOS: oxide [0, 40 nm] x [0, 3 nm] under p-silicon up to 40 nm, the interface at
// y = 3 nm with mirrored pairs at -+ s/2 (s = 0.4 nm); the gate electrode on the oxide's y_min, the
// substrate on the silicon's top.
constexpr double delaunay_s = 0.4e-7;
constexpr double delaunay_t = 3e-7;
device::Device delaunay_mos() {
    const double W = 4e-6, H = 4e-6, t = delaunay_t;
    mesh::Geometry2D g{{{{{0, 0, 0}, {W, 0, 0}, {W, t, 0}, {0, t, 0}}, 0},
                        {{{0, t, 0}, {W, t, 0}, {W, H, 0}, {0, H, 0}}, 1}},
                       {{"gate", {0, 0, 0}, {W, 0, 0}}, {"substrate", {0, H, 0}, {W, H, 0}}}};
    mesh::MeshOptions2D o;
    o.max_size = 0.6e-6;
    o.interface_resolution = delaunay_s;
    auto gm = mesh::generate_mesh(g, o);
    REQUIRE(gm.has_value());
    const std::size_t n = gm->mesh.node_count();
    std::vector<double> acceptors(n, 0.0);
    std::vector<device::RegionId> region(n);
    for (std::size_t i = 0; i < n; ++i) {
        region[i] = gm->node_region[i];
        if (region[i] == 1) acceptors[i] = 1e17;
    }
    auto gate_nodes = gm->mesh.find_boundary("gate")->nodes;
    auto substrate = gm->mesh.find_boundary("substrate")->nodes;
    device::Contact gate{"gate", device::ContactKind::electrode, std::move(gate_nodes)};
    gate.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.4};
    return *device::Device::create(
        {.mesh = std::move(gm->mesh),
         .temperature_K = 300.0,
         .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
         .node_region = std::move(region),
         .donors = std::vector<double>(n, 0.0),
         .acceptors = std::move(acceptors),
         .contacts = {std::move(gate),
                      {"substrate", device::ContactKind::ohmic, std::move(substrate)}},
         .interface_certificate = std::move(gm->certificate)});
}

// A lumped-gate MOS-C: p-silicon (1e17) on [0, 30 nm] (2D: 3 um wide) with nodes 0.2 nm apart at
// the gate face x_min, a 2 nm n+ poly gate on x_min and the substrate contact on x_max (3D: 3 um
// wide in y and z).
constexpr double lumped_h = 0.2e-7;
device::Device lumped(int dimension) {
    std::vector<double> x;
    for (int k = 0; k < 10; ++k) x.push_back(k * lumped_h);
    for (const double v : uniform(2e-7, 3e-6, 30)) x.push_back(v);
    const auto t = uniform(0.0, 3e-6, 4);
    mesh::Mesh m = dimension == 1   ? *mesh::make_tensor_grid(x)
                   : dimension == 2 ? *mesh::make_tensor_grid(x, t)
                                    : *mesh::make_tensor_grid(x, t, t);
    const std::size_t n = m.node_count();
    device::Contact gate{"gate", device::ContactKind::gate, m.find_boundary("x_min")->nodes};
    gate.gate = {.boundary = "x_min", .oxide_thickness_cm = 2e-7};
    auto substrate = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::vector<double>(n, 0.0),
         .acceptors = std::vector<double>(n, 1e17),
         .contacts = {std::move(gate),
                      {"substrate", device::ContactKind::ohmic, std::move(substrate)}}});
}

bool close(double a, double b, double rel) { return std::abs(a - b) <= rel * std::abs(b); }

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

// The house gate (as insulator_test.cpp): every column, central differences with steps
// 1e-7 max(|x|, 1), the worst error over the column's largest entry.
template <class System>
double fd_jacobian_error(const System& system, std::vector<double> x) {
    const std::size_t n = system.unknowns();
    std::vector<double> f(n), fp(n), fm(n);
    auto jacobian = system.make_jacobian();
    system.evaluate(x, f, jacobian);
    const auto j = dense(jacobian);
    double worst = 0.0;
    for (std::size_t c = 0; c < n; ++c) {
        const double base = x[c];
        const double step = 1e-7 * std::max(std::abs(base), 1.0);
        x[c] = base + step;
        const double up = x[c];
        system.residual(x, fp);
        x[c] = base - step;
        const double down = x[c];
        system.residual(x, fm);
        x[c] = base;
        double scale = 1e-300, error = 0.0;
        for (std::size_t r = 0; r < n; ++r) scale = std::max(scale, std::abs(j[r][c]));
        for (std::size_t r = 0; r < n; ++r) {
            error = std::max(error, std::abs((fp[r] - fm[r]) / (up - down) - j[r][c]));
        }
        worst = std::max(worst, error / scale);
    }
    return worst;
}

// As insulator_test.cpp's probe: the charge-neutral equilibrium state with noise, and the
// interface nodes' densities raised out of equilibrium so the trap terms are resolved.
std::vector<double> probe_state(const System& system, const System& poisson,
                                const device::Device& d, std::uint64_t seed) {
    auto x = system.state_from_potential(charge_neutral_potential(poisson));
    Noise noise{seed};
    for (std::size_t i = 0; i < system.node_count(); ++i) {
        x[3 * i] += 0.02 * noise.next();
        if (d.is_insulator(static_cast<mesh::NodeId>(i))) continue;
        x[3 * i + 1] *= 1.0 + 0.01 * noise.next();
        x[3 * i + 2] *= 1.0 + 0.01 * noise.next();
        const double px = d.mesh().points()[i][0];
        if (px > 0.0 && px < 1e-7) {
            x[3 * i + 1] = 2e-3 * (1.0 + 0.01 * noise.next());
            x[3 * i + 2] = 5e-4 * (1.0 + 0.01 * noise.next());
        }
    }
    system.stamp_contacts(x);
    return x;
}

}  // namespace

TEST_CASE("quantum confinement assemble: distances on tensor meshes (QC3)") {
    for (const int D : {1, 2, 3}) {
        const device::Device d = mos(D);
        const auto c = *assemble::mlda_factors(d, mlda());
        const double ln = physics::mlda_wavelength_cm(0.98, 300.0);
        const double lp = physics::mlda_wavelength_cm(0.49, 300.0);
        CHECK(c.wavelength_n_cm == ln);
        CHECK(c.wavelength_p_cm == lp);
        std::size_t checked = 0;
        for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
            if (d.is_insulator(static_cast<mesh::NodeId>(i))) continue;
            const double x = d.mesh().points()[i][0];
            CAPTURE(D, i, x);
            REQUIRE(close(c.distance_cm[i], x, 1e-12));
            REQUIRE(c.half_cell[i] == 0);
            REQUIRE(c.f_n[i] == physics::mlda_factor(x, ln));
            REQUIRE(c.f_p[i] == physics::mlda_factor(x, lp));
            ++checked;
        }
        CHECK(checked > 0);
    }
}

TEST_CASE("quantum confinement assemble: at a corner the nearer wall decides (QC3)") {
    const device::Device d = corner();
    const auto c = *assemble::mlda_factors(d, mlda());
    for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
        if (d.is_insulator(static_cast<mesh::NodeId>(i))) continue;
        const auto& p = d.mesh().points()[i];
        CAPTURE(i, p[0], p[1]);
        REQUIRE(close(c.distance_cm[i], std::min(p[0], p[1]), 1e-12));
    }
}

TEST_CASE("quantum confinement assemble: distances on a Delaunay mesh (QC3)") {
    const device::Device d = delaunay_mos();
    const auto c = *assemble::mlda_factors(d, mlda());
    std::size_t adjacent = 0;
    for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
        if (d.is_insulator(static_cast<mesh::NodeId>(i))) continue;
        const double depth = d.mesh().points()[i][1] - delaunay_t;
        REQUIRE(depth > 0.0);
        // Never nearer than the interface plane; exact for the mirrored pairs' silicon nodes.
        REQUIRE(c.distance_cm[i] >= depth * (1.0 - 1e-12));
        if (std::abs(depth - 0.5 * delaunay_s) <= 1e-6 * delaunay_s) {
            CAPTURE(i, depth);
            REQUIRE(close(c.distance_cm[i], depth, 1e-12));
            ++adjacent;
        }
    }
    CHECK(adjacent >= 10);
}

TEST_CASE("quantum confinement assemble: a lumped-gate node takes the half-cell mean (QC3)") {
    for (const int D : {1, 2}) {
        const device::Device d = lumped(D);
        const auto c = *assemble::mlda_factors(d, mlda());
        const auto gate = d.contacts()[0].nodes;
        for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
            const double x = d.mesh().points()[i][0];
            const bool on_gate = std::find(gate.begin(), gate.end(),
                                           static_cast<mesh::NodeId>(i)) != gate.end();
            CAPTURE(D, i, x);
            REQUIRE(c.half_cell[i] == (on_gate ? 1 : 0));
            REQUIRE(close(c.distance_cm[i], x, 1e-12));
            if (on_gate) {
                REQUIRE(c.f_n[i] == physics::mlda_mean_factor(0.5 * lumped_h, c.wavelength_n_cm));
                REQUIRE(c.f_p[i] == physics::mlda_mean_factor(0.5 * lumped_h, c.wavelength_p_cm));
                REQUIRE(c.f_n[i] > 0.0);
            }
        }
    }
}

TEST_CASE("quantum confinement assemble: N_c f_n and N_v f_p, bands and traps in place (QC4)") {
    for (const bool fd : {false, true}) {
        for (const int D : {1, 2}) {
            const device::Device d = mos(D);
            const auto scaling = *assemble::make_scaling(d);
            assemble::PhysicsModels classical;
            classical.fermi_dirac = fd;
            const auto off = *create_equilibrium(d, scaling, classical);
            const auto on = *create_equilibrium(d, scaling, mlda(fd));
            const auto c = *assemble::mlda_factors(d, mlda(fd));
            // A potential in inversion near the oxide and in the bulk, with noise.
            auto psi = charge_neutral_potential(off);
            Noise noise{3};
            for (std::size_t i = 0; i < psi.size(); ++i) {
                if (!d.is_insulator(static_cast<mesh::NodeId>(i))) psi[i] += 0.3 * noise.next();
            }
            const std::size_t n = psi.size();
            std::vector<double> n0(n), p0(n), n1(n), p1(n);
            equilibrium_carriers(off, psi, n0, p0);
            equilibrium_carriers(on, psi, n1, p1);
            const auto b0 = band_edges(off, psi);
            const auto b1 = band_edges(on, psi);
            for (std::size_t i = 0; i < n; ++i) {
                CAPTURE(fd, D, i);
                if (d.is_insulator(static_cast<mesh::NodeId>(i))) {
                    REQUIRE(n1[i] == n0[i]);
                    REQUIRE(p1[i] == p0[i]);
                    continue;
                }
                if (c.f_n[i] == 1.0 && c.f_p[i] == 1.0) {  // unchanged to the bit
                    REQUIRE(n1[i] == n0[i]);
                    REQUIRE(p1[i] == p0[i]);
                    REQUIRE(b1.conduction[i] == b0.conduction[i]);
                    continue;
                }
                REQUIRE(close(n1[i], c.f_n[i] * n0[i], 1e-14));
                REQUIRE(close(p1[i], c.f_p[i] * p0[i], 1e-14));
                REQUIRE(std::abs(b1.conduction[i] - b0.conduction[i]) <=
                        1e-14 * std::max(1.0, std::abs(b0.conduction[i])));
                REQUIRE(std::abs(b1.valence[i] - b0.valence[i]) <=
                        1e-14 * std::max(1.0, std::abs(b0.valence[i])));
            }
            // Trap occupancy depends on E_F - E_t, which MLDA does not move.
            const auto q0 = interface_trap_charges(off, psi);
            const auto q1 = interface_trap_charges(on, psi);
            REQUIRE(q0.size() == q1.size());
            for (std::size_t k = 0; k < q0.size(); ++k) {
                CAPTURE(fd, D, k, q0[k], q1[k]);
                REQUIRE(q0[k] != 0.0);
                REQUIRE(std::abs(q1[k] - q0[k]) <= 1e-12 * std::abs(q0[k]));
            }
        }
    }
}

TEST_CASE("quantum confinement assemble: drift-diffusion states carry the factors (QC4)") {
    const device::Device d = mos(2);
    const auto scaling = *assemble::make_scaling(d);
    const auto off = *System::create(d, scaling);
    const auto on = *System::create(d, scaling, mlda());
    const auto poisson = *create_equilibrium(d, scaling);
    const auto c = *assemble::mlda_factors(d, mlda());
    const auto psi = charge_neutral_potential(poisson);
    const auto x0 = off.state_from_potential(psi);
    const auto x1 = on.state_from_potential(psi);
    for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
        if (d.is_insulator(static_cast<mesh::NodeId>(i)) || is_contact(poisson, i)) continue;
        CAPTURE(i);
        REQUIRE(close(x1[3 * i + 1], c.f_n[i] * x0[3 * i + 1], 1e-14));
        REQUIRE(close(x1[3 * i + 2], c.f_p[i] * x0[3 * i + 2], 1e-14));
    }
    // The pattern is the model-off one (no new couplings).
    const auto j0 = off.make_jacobian();
    const auto j1 = on.make_jacobian();
    REQUIRE(std::vector<linalg::Index>(j0.col_indices().begin(), j0.col_indices().end()) ==
            std::vector<linalg::Index>(j1.col_indices().begin(), j1.col_indices().end()));
}

TEST_CASE("quantum confinement assemble: the finite-difference Jacobian gate (QC6)") {
    const std::vector<double> bias{-0.4, 0.1};
    const std::vector<double> equilibrium_bias{-0.4, 0.0};
    for (const int D : {1, 2, 3}) {
        const device::Device d = mos(D);
        const auto scaling = *assemble::make_scaling(d);
        for (const bool fd : {false, true}) {
            auto poisson = *create_equilibrium(d, scaling, mlda(fd));
            REQUIRE(poisson.set_bias(equilibrium_bias).has_value());
            auto psi = charge_neutral_potential(poisson);
            Noise noise{17};
            for (double& v : psi) v += 0.3 * noise.next();
            const double ep = fd_jacobian_error(poisson, psi);
            auto system = *System::create(d, scaling, mlda(fd));
            REQUIRE(system.set_bias(bias).has_value());
            const auto x = probe_state(system, poisson, d, 23);
            const double dd = fd_jacobian_error(system, x);
            CAPTURE(D, fd, ep, dd);
            REQUIRE(ep <= 1e-6);
            REQUIRE(dd <= 1e-6);
        }
    }
    // With Lombardi surface mobility as well (QC11 accepts the two together).
    for (const int D : {2, 3}) {
        const device::Device d = mos(D);
        const auto scaling = *assemble::make_scaling(d);
        auto models = mlda();
        models.surface_mobility = assemble::SurfaceMobility::lombardi;
        auto poisson = *create_equilibrium(d, scaling, models);
        REQUIRE(poisson.set_bias(equilibrium_bias).has_value());
        auto system = *System::create(d, scaling, models);
        REQUIRE(system.set_bias(bias).has_value());
        const auto x = probe_state(system, poisson, d, 23);
        const double dd = fd_jacobian_error(system, x);
        CAPTURE(D, dd);
        REQUIRE(dd <= 1e-6);
    }
    // A lumped gate and a Delaunay mesh, Boltzmann and Fermi-Dirac.
    std::vector<device::Device> others;
    others.push_back(lumped(2));
    others.push_back(delaunay_mos());
    for (const device::Device& d : others) {
        const auto scaling = *assemble::make_scaling(d);
        for (const bool fd : {false, true}) {
            auto poisson = *create_equilibrium(d, scaling, mlda(fd));
            REQUIRE(poisson.set_bias(std::vector<double>{0.8, 0.0}).has_value());
            auto system = *System::create(d, scaling, mlda(fd));
            REQUIRE(system.set_bias(std::vector<double>{0.8, 0.05}).has_value());
            auto x = system.state_from_potential(charge_neutral_potential(poisson));
            Noise noise{29};
            for (std::size_t i = 0; i < system.node_count(); ++i) {
                x[3 * i] += 0.3 * noise.next();
                if (d.is_insulator(static_cast<mesh::NodeId>(i))) continue;
                x[3 * i + 1] *= 1.0 + 0.01 * noise.next();
                x[3 * i + 2] *= 1.0 + 0.01 * noise.next();
            }
            system.stamp_contacts(x);
            const double dd = fd_jacobian_error(system, x);
            CAPTURE(fd, dd);
            REQUIRE(dd <= 1e-6);
        }
    }
}

TEST_CASE("quantum confinement assemble: refusals (QC11)") {
    const device::Device d = mos(1);
    const auto scaling = *assemble::make_scaling(d);
    const auto refused = [&](const device::Device& dev, const assemble::PhysicsModels& m,
                             const std::string& what) {
        const auto s = *assemble::make_scaling(dev);
        const auto dd = System::create(dev, s, m);
        const auto ep = create_equilibrium(dev, s, m);
        REQUIRE_FALSE(dd.has_value());
        REQUIRE_FALSE(ep.has_value());
        CAPTURE(dd.error().message);
        CHECK(dd.error().code == base::ErrorCode::invalid_input);
        CHECK(dd.error().message.find("quantum_confinement mlda: ") != std::string::npos);
        CHECK(dd.error().message.find(what) != std::string::npos);
        CHECK(ep.error().message == dd.error().message);
    };
    SECTION("electrothermal") {
        // This device has no thermal contacts, so the assembler refuses electrothermal anyway;
        // the factors' own refusal is checked directly.
        auto m = mlda();
        m.electrothermal = true;
        const auto c = assemble::mlda_factors(d, m);
        REQUIRE_FALSE(c.has_value());
        CHECK(c.error().code == base::ErrorCode::invalid_input);
        CHECK(c.error().message.find("quantum_confinement mlda: electrothermal") !=
              std::string::npos);
        CHECK_FALSE(System::create(d, scaling, m).has_value());
    }
    SECTION("masses") {
        for (const double bad : {0.0, -0.5, std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::infinity()}) {
            auto m = mlda();
            m.mlda_mass_n = bad;
            refused(d, m, "masses");
            m = mlda();
            m.mlda_mass_p = bad;
            refused(d, m, "masses");
        }
    }
    SECTION("not silicon, and silicon with another parameter changed") {
        const auto build = [](const physics::SemiconductorParameters& p) {
            mesh::Mesh m = *mesh::make_tensor_grid(
                *mesh::straddle_interface(uniform(-5e-7, 3e-6, 36), 0.0, 0.5e-7));
            const std::size_t n = m.node_count();
            std::vector<device::RegionId> region(n);
            std::vector<double> acceptors(n, 0.0);
            for (std::size_t i = 0; i < n; ++i) {
                region[i] = m.points()[i][0] < 0.0 ? 0 : 1;
                if (region[i] == 1) acceptors[i] = 1e17;
            }
            auto substrate = m.find_boundary("x_max")->nodes;
            return *device::Device::create(
                {.mesh = std::move(m),
                 .temperature_K = 300.0,
                 .regions = {{"oxide", physics::silicon_dioxide()},
                             {"body", *physics::Semiconductor::create(p)}},
                 .node_region = std::move(region),
                 .donors = std::vector<double>(n, 0.0),
                 .acceptors = std::move(acceptors),
                 .contacts = {{"substrate", device::ContactKind::ohmic, std::move(substrate)}}});
        };
        refused(build(physics::germanium_parameters), mlda(), "not silicon");
        auto other = physics::silicon_parameters;
        other.radiative_cm3_s = 1e-14;
        const device::Device ok = build(other);
        const auto s = *assemble::make_scaling(ok);
        CHECK(System::create(ok, s, mlda()).has_value());
        CHECK(create_equilibrium(ok, s, mlda()).has_value());
    }
    SECTION("nothing confines") {
        mesh::Mesh m = *mesh::make_tensor_grid(uniform(0.0, 1e-5, 20));
        const std::size_t n = m.node_count();
        auto left = m.find_boundary("x_min")->nodes;
        auto right = m.find_boundary("x_max")->nodes;
        const device::Device block = *device::Device::create(
            {.mesh = std::move(m),
             .temperature_K = 300.0,
             .regions = {{"silicon", physics::silicon()}},
             .node_region = std::vector<device::RegionId>(n, 0),
             .donors = std::vector<double>(n, 1e17),
             .acceptors = std::vector<double>(n, 0.0),
             .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                          {"right", device::ContactKind::ohmic, std::move(right)}}});
        refused(block, mlda(), "no insulator");
    }
    SECTION("a thermionic-emission interface") {
        auto shifted = physics::silicon_parameters;
        shifted.electron_affinity_eV += 0.1;
        mesh::Mesh m = *mesh::make_tensor_grid(*mesh::straddle_interface(
            *mesh::straddle_interface(uniform(-5e-7, 3e-6, 36), 0.0, 0.5e-7), 1.5e-6, 0.5e-7));
        const std::size_t n = m.node_count();
        std::vector<device::RegionId> region(n);
        std::vector<double> acceptors(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            const double x = m.points()[i][0];
            region[i] = x < 0.0 ? 0 : x < 1.5e-6 ? 1 : 2;
            if (region[i] != 0) acceptors[i] = 1e17;
        }
        auto substrate = m.find_boundary("x_max")->nodes;
        const device::Device hetero = *device::Device::create(
            {.mesh = std::move(m),
             .temperature_K = 300.0,
             .regions = {{"oxide", physics::silicon_dioxide()},
                         {"silicon", physics::silicon()},
                         {"other", *physics::Semiconductor::create(shifted)}},
             .node_region = std::move(region),
             .donors = std::vector<double>(n, 0.0),
             .acceptors = std::move(acceptors),
             .contacts = {{"substrate", device::ContactKind::ohmic, std::move(substrate)}},
             .interfaces = {{"silicon", "other", device::InterfaceTransport::thermionic_emission}}});
        refused(hetero, mlda(), "thermionic");
    }
    SECTION("accepted with Lombardi") {
        auto m = mlda();
        m.surface_mobility = assemble::SurfaceMobility::lombardi;
        CHECK(System::create(d, scaling, m).has_value());
    }
}

// Unit 62 (QP10, QV11): density gradient's refusals, the same message from both assemblers, and the
// accepted combinations, which both assemblers build (slice 2; stride 5, and 3 at equilibrium).
namespace {

assemble::PhysicsModels dg() {
    assemble::PhysicsModels m;
    m.quantum_confinement = assemble::QuantumConfinement::density_gradient;
    return m;
}

struct LumpedOptions {
    double eps_r = 3.9;
    double band_gap_eV = 0.0;  // 0: no band data of its own
    std::uint8_t tunnelling = 0;
    bool schottky = false;
};

device::Device lumped_dg(const LumpedOptions& o) {
    const auto x = uniform(0.0, 3e-6, 31);
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    device::Contact gate{"gate", device::ContactKind::gate, m.find_boundary("x_min")->nodes};
    gate.gate = {.boundary = "x_min", .oxide_thickness_cm = 2e-7, .oxide_relative_permittivity = o.eps_r};
    if (o.band_gap_eV > 0.0) {
        gate.gate.oxide_affinity_eV = 1.05;
        gate.gate.oxide_band_gap_eV = o.band_gap_eV;
        gate.gate.oxide_mass_n = 0.5;
        gate.gate.oxide_mass_p = 0.5;
    }
    gate.tunnelling = o.tunnelling;
    device::Contact substrate{"substrate", device::ContactKind::ohmic, m.find_boundary("x_max")->nodes};
    if (o.schottky) {
        substrate.kind = device::ContactKind::schottky;
        substrate.schottky = {.boundary = "x_max", .by_work_function = true, .work_function_eV = 4.6};
    }
    return *device::Device::create({.mesh = std::move(m),
                                    .temperature_K = 300.0,
                                    .regions = {{"silicon", physics::silicon()}},
                                    .node_region = std::vector<device::RegionId>(n, 0),
                                    .donors = std::vector<double>(n, 0.0),
                                    .acceptors = std::vector<double>(n, 1e17),
                                    .contacts = {std::move(gate), std::move(substrate)}});
}

// A meshed oxide of the given parameters on a p-silicon bar.
device::Device meshed_dg(const physics::InsulatorParameters& oxide) {
    mesh::Mesh m = *mesh::make_tensor_grid(*mesh::straddle_interface(uniform(-5e-7, 3e-6, 36), 0.0, 0.5e-7));
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n);
    std::vector<double> acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        region[i] = m.points()[i][0] < 0.0 ? 0 : 1;
        if (region[i] == 1) acceptors[i] = 1e17;
    }
    device::Contact gate{"gate", device::ContactKind::electrode, m.find_boundary("x_min")->nodes};
    gate.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.4};
    auto substrate = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"oxide", *physics::Insulator::create(oxide)}, {"silicon", physics::silicon()}},
         .node_region = std::move(region),
         .donors = std::vector<double>(n, 0.0),
         .acceptors = std::move(acceptors),
         .contacts = {std::move(gate), {"substrate", device::ContactKind::ohmic, std::move(substrate)}}});
}

}  // namespace

TEST_CASE("density gradient assemble: refusals and acceptances (QV11)") {
    const auto outcome = [](const device::Device& dev, const assemble::PhysicsModels& m) {
        const auto s = *assemble::make_scaling(dev);
        const auto dd = System::create(dev, s, m);
        const auto ep = create_equilibrium(dev, s, m);
        REQUIRE_FALSE(dd.has_value());
        REQUIRE_FALSE(ep.has_value());
        CHECK(dd.error().code == base::ErrorCode::invalid_input);
        CHECK(ep.error().message == dd.error().message);
        return dd.error().message;
    };
    const auto refused = [&](const device::Device& dev, const assemble::PhysicsModels& m, const std::string& what) {
        const std::string msg = outcome(dev, m);
        CAPTURE(msg);
        CHECK(msg.find("quantum_confinement density_gradient: ") != std::string::npos);
        CHECK(msg.find(what) != std::string::npos);
        CHECK(msg.find("not built yet") == std::string::npos);
    };
    const auto accepted = [&](const device::Device& dev, const assemble::PhysicsModels& m) {
        const auto direct = assemble::check_density_gradient(dev, m);
        CHECK(direct.has_value());
        const auto s = *assemble::make_scaling(dev);
        const auto dd = System::create(dev, s, m);
        const auto ep = create_equilibrium(dev, s, m);
        CHECK(dd.has_value());
        CHECK(ep.has_value());
        if (dd) CHECK(dd->stride() == 5);
        if (ep) CHECK(ep->stride() == 3);
    };
    const auto direct_refused = [](const device::Device& dev, const assemble::PhysicsModels& m,
                                   const std::string& what) {
        const auto r = assemble::check_density_gradient(dev, m);
        REQUIRE_FALSE(r.has_value());
        CAPTURE(r.error().message);
        CHECK(r.error().code == base::ErrorCode::invalid_input);
        CHECK(r.error().message.find(what) != std::string::npos);
    };
    const device::Device meshed = mos(1);

    SECTION("accepted: meshed SiO2, lumped SiO2, a lumped stack with band data, Lombardi, both statistics") {
        accepted(meshed, dg());
        accepted(mos(2), dg());
        accepted(lumped_dg({}), dg());
        accepted(lumped_dg({.eps_r = 25.0, .band_gap_eV = 5.8}), dg());
        accepted(meshed_dg(physics::silicon_dioxide_tunnelling_parameters), dg());
        auto m = dg();
        m.surface_mobility = assemble::SurfaceMobility::lombardi;
        m.fermi_dirac = true;
        m.btbt_local = true;
        accepted(mos(2), m);
    }
    SECTION("masses and gammas") {
        for (const double bad : {0.0, -0.5, std::numeric_limits<double>::quiet_NaN(),
                                 std::numeric_limits<double>::infinity()}) {
            for (double assemble::PhysicsModels::* field :
                 {&assemble::PhysicsModels::dg_mass_n, &assemble::PhysicsModels::dg_mass_p,
                  &assemble::PhysicsModels::dg_gamma_n, &assemble::PhysicsModels::dg_gamma_p}) {
                auto m = dg();
                m.*field = bad;
                refused(meshed, m, "masses and gammas");
            }
        }
    }
    SECTION("models it is refused with") {
        // The assembler may refuse some of these itself first (no thermal contact, no path), so the check is
        // also called directly.
        auto m = dg();
        m.electrothermal = true;
        direct_refused(meshed, m, "electrothermal");
        CHECK_FALSE(System::create(meshed, *assemble::make_scaling(meshed), m).has_value());
        m = dg();
        m.energy_balance = true;
        direct_refused(meshed, m, "energy_balance");
        CHECK_FALSE(System::create(meshed, *assemble::make_scaling(meshed), m).has_value());
        m = dg();
        m.btbt_nonlocal = assemble::NonlocalTunnelling::kane;
        direct_refused(meshed, m, "nonlocal band-to-band");
        CHECK_FALSE(System::create(meshed, *assemble::make_scaling(meshed), m).has_value());
    }
    SECTION("contacts it is refused with") {
        refused(lumped_dg({.tunnelling = device::all_tunnel_processes}), dg(), "gate tunnelling");
        refused(lumped_dg({.schottky = true}), dg(), "Schottky");
        refused(lumped_dg({.eps_r = 25.0}), dg(), "lumped gate's stack has no band data");
    }
    SECTION("insulators without band data, and not silicon") {
        auto other = physics::silicon_dioxide_parameters;
        other.eps_r = 7.5;
        refused(meshed_dg(other), dg(), "insulator region 'oxide' has no band data");
        mesh::Mesh m = *mesh::make_tensor_grid(uniform(0.0, 3e-6, 31));
        const std::size_t n = m.node_count();
        auto left = m.find_boundary("x_min")->nodes;
        auto right = m.find_boundary("x_max")->nodes;
        const device::Device ge = *device::Device::create(
            {.mesh = std::move(m),
             .temperature_K = 300.0,
             .regions = {{"body", *physics::Semiconductor::create(physics::germanium_parameters)}},
             .node_region = std::vector<device::RegionId>(n, 0),
             .donors = std::vector<double>(n, 1e17),
             .acceptors = std::vector<double>(n, 0.0),
             .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                          {"right", device::ContactKind::ohmic, std::move(right)}}});
        refused(ge, dg(), "not silicon");
    }
    SECTION("off and MLDA are unaffected by DG's parameters") {
        auto m = assemble::PhysicsModels{};
        m.dg_mass_n = -1.0;  // not read with the model off
        CHECK(System::create(meshed, *assemble::make_scaling(meshed), m).has_value());
        auto ml = mlda();
        ml.dg_gamma_p = 0.0;  // not read with MLDA
        CHECK(System::create(meshed, *assemble::make_scaling(meshed), ml).has_value());
    }
}

// Unit 62 slice 2 (QV1, QV3-QV5, QV9): the density-gradient rows, their Jacobian, the equilibrium
// table's slaving and symmetry.
namespace {

constexpr std::size_t dg_m = 5;

// The quantum potentials (in kT) a drift-diffusion state's Lambda rows give with Lambda = 0: the
// rows are V L + (the rest), so L_i = -f_i / J_ii (0 on the held nodes).
std::pair<std::vector<double>, std::vector<double>> dg_lambda(const System& system,
                                                              const std::vector<double>& x) {
    const std::size_t n = system.node_count();
    std::vector<double> f(system.unknowns());
    auto j = system.make_jacobian();
    system.evaluate(x, f, j);
    const auto diagonal = [&](std::size_t r) {
        for (auto k = static_cast<std::size_t>(j.row_offsets()[r]);
             k < static_cast<std::size_t>(j.row_offsets()[r + 1]); ++k) {
            if (static_cast<std::size_t>(j.col_indices()[k]) == r) return j.values()[k];
        }
        return 0.0;
    };
    std::vector<double> ln(n), lp(n);
    for (std::size_t i = 0; i < n; ++i) {
        ln[i] = -f[dg_m * i + 3] / diagonal(dg_m * i + 3);
        lp[i] = -f[dg_m * i + 4] / diagonal(dg_m * i + 4);
    }
    return {ln, lp};
}

// A silicon bar (1e16 p, [0, 2 um]) on the axis x with ohmic contacts at both ends; 2D and 3D extrude it (4
// nodes across), so every column is the 1D bar.
device::Device dg_bar(int dimension, const std::vector<double>& x) {
    const auto t = uniform(0.0, 1e-6, 4);
    mesh::Mesh m = dimension == 1   ? *mesh::make_tensor_grid(x)
                   : dimension == 2 ? *mesh::make_tensor_grid(x, t)
                                    : *mesh::make_tensor_grid(x, t, t);
    const std::size_t n = m.node_count();
    auto left = m.find_boundary("x_min")->nodes;
    auto right = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::vector<double>(n, 0.0),
         .acceptors = std::vector<double>(n, 1e16),
         .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                      {"right", device::ContactKind::ohmic, std::move(right)}}});
}

// [0, 2 um] in `cells` cells: uniform, or graded smoothly by x = 2 um (s + s^2) / 2 (ratio 3 end to
// end).
std::vector<double> dg_axis(int cells, bool graded) {
    std::vector<double> x(static_cast<std::size_t>(cells) + 1);
    for (int k = 0; k <= cells; ++k) {
        const double s = static_cast<double>(k) / cells;
        x[static_cast<std::size_t>(k)] = graded ? 2e-4 * 0.5 * (s + s * s) : 2e-4 * s;
    }
    return x;
}

// A 1D MOS wall at spacing h: p-silicon (1e17) on [0, 6 nm] with an ohmic contact at 6 nm. Meshed:
// SiO2 on [-2 nm, 0) with a metal electrode on x_min, nodes at -+h/2, -+3h/2, ... (the interface at
// the half-edge). Lumped: the silicon's first node on x = 0 under a lumped SiO2 gate, nodes h apart.
device::Device dg_wall(double h, bool lumped_gate) {
    std::vector<double> x;
    if (lumped_gate) {
        for (double v = 0.0; v < 6e-7 + 1e-3 * h; v += h) x.push_back(v);
    } else {
        for (double v = -0.5 * h; v > -2e-7 - 1e-3 * h; v -= h) x.push_back(v);
        std::reverse(x.begin(), x.end());
        for (double v = 0.5 * h; v < 6e-7 + 1e-3 * h; v += h) x.push_back(v);
    }
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n, 0);
    std::vector<double> acceptors(n, 1e17);
    for (std::size_t i = 0; i < n && !lumped_gate; ++i) {
        const bool oxide = m.points()[i][0] < 0.0;
        region[i] = oxide ? 0 : 1;
        if (oxide) acceptors[i] = 0.0;
    }
    auto gate_nodes = m.find_boundary("x_min")->nodes;
    auto substrate = m.find_boundary("x_max")->nodes;
    device::Contact gate{"gate", lumped_gate ? device::ContactKind::gate : device::ContactKind::electrode,
                         std::move(gate_nodes)};
    if (lumped_gate) {
        gate.gate = {.boundary = "x_min", .oxide_thickness_cm = 2e-7};
    } else {
        gate.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.4};
    }
    std::vector<device::Region> regions;
    if (!lumped_gate) regions.push_back({"oxide", physics::silicon_dioxide()});
    regions.push_back({"silicon", physics::silicon()});
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = std::move(regions),
         .node_region = std::move(region),
         .donors = std::vector<double>(n, 0.0),
         .acceptors = std::move(acceptors),
         .contacts = {std::move(gate), {"substrate", device::ContactKind::ohmic, std::move(substrate)}}});
}

// The state with sqrt n = rn(x) and sqrt p = rp(x) (densities 1e-3 rn^2, 2e-3 rp^2), Lambda = 0, on
// the semiconductor nodes.
template <class Rn, class Rp>
std::vector<double> dg_profile(const System& system, const device::Device& d, Rn rn, Rp rp) {
    auto x = system.state_from_potential(std::vector<double>(system.node_count(), 0.0));
    for (std::size_t i = 0; i < system.node_count(); ++i) {
        x[dg_m * i + 3] = x[dg_m * i + 4] = 0.0;
        if (d.is_insulator(static_cast<mesh::NodeId>(i))) continue;
        const double px = d.mesh().points()[i][0];
        x[dg_m * i + 1] = 1e-3 * rn(px) * rn(px);
        x[dg_m * i + 2] = 2e-3 * rp(px) * rp(px);
    }
    return x;
}

// A drift-diffusion probe with the quantum potentials resolved: the slaved state of the
// charge-neutral potential with noise (0.02 on psi, 0.3 on Lambda), 1% noise on the densities, the
// contacts stamped. The densities stay smooth from node to node: a step (as probe_state's raised
// interface densities) puts terms of 1e8 into a Lambda row whose own entry is V, and the central
// difference of that row is rounding at eps |f| / step, far above 1e-6 of its column.
std::vector<double> dg_probe(const System& system, const System& poisson, const device::Device& d,
                             std::uint64_t seed) {
    const std::size_t n = system.node_count();
    const auto psi = charge_neutral_potential(poisson);
    Noise noise{seed};
    std::vector<double> ln(n), lp(n), dpsi(psi.begin(), psi.end());
    for (std::size_t i = 0; i < n; ++i) {
        dpsi[i] += 0.02 * noise.next();
        ln[i] = 0.3 * noise.next();
        lp[i] = 0.3 * noise.next();
    }
    auto x = system.state_from_quantum_potential(dpsi, ln, lp);
    for (std::size_t i = 0; i < n; ++i) {
        if (d.is_insulator(static_cast<mesh::NodeId>(i))) continue;
        x[dg_m * i + 1] *= 1.0 + 0.01 * noise.next();
        x[dg_m * i + 2] *= 1.0 + 0.01 * noise.next();
    }
    system.stamp_contacts(x);
    return x;
}

results::ConvergenceRecord dg_newton(const System& sys, std::vector<double>& x, bool pivot_check = true);

// An equilibrium (psi, Lambda_n, Lambda_p) probe: the classical equilibrium potential with noise 1e-3
// and Lambda of order 0.3. Off a Poisson solution (the charge-neutral potential, or noise 0.02 to
// 0.1 on psi) the Poisson rows reach |f| ~ 20 to 130 and their central difference in a Lambda column
// (entry -V p or V n, small) is rounding at 1e-6 of it (measured on 3D Fermi-Dirac, incomplete
// ionization and Delaunay; the entry matches the difference at step 1e-5 to 9 digits).
std::vector<double> dg_equilibrium_probe(const System& poisson, std::uint64_t seed) {
    const auto start = charge_neutral_potential(poisson);
    std::vector<double> psi(start.begin(), start.end());
    dg_newton(poisson, psi);
    Noise noise{seed};
    std::vector<double> x(3 * psi.size());
    for (std::size_t i = 0; i < psi.size(); ++i) {
        x[3 * i] = psi[i] + 1e-3 * noise.next();
        x[3 * i + 1] = 0.3 * noise.next();
        x[3 * i + 2] = 0.3 * noise.next();
    }
    return x;
}

// min_pivot_ratio 0 where the classical system is as weakly anchored (the 1D MOS in inversion with
// only a substrate contact: 2.3e-15 without the model, 9.8e-15 with it).
results::ConvergenceRecord dg_newton(const System& sys, std::vector<double>& x, bool pivot_check) {
    linalg::SolverConfig config;
    if (!pivot_check) config.min_pivot_ratio = 0.0;
    auto solver = *linalg::LinearSolver::create(config);
    results::ConvergenceRecord record;
    solve::NewtonOptions options;
    options.tol_update = 1e-11;
    const auto ok = solve::newton_solve(sys, x, options, solver, record);
    if (!ok) FAIL(ok.error().message);
    return record;
}

// The DG equilibrium state (psi, Lambda_n, Lambda_p per node) at `bias`, from the classical one.
std::vector<double> dg_equilibrium(const device::Device& d, const assemble::PhysicsModels& m,
                                   const std::vector<double>& bias) {
    const auto scaling = *assemble::make_scaling(d);
    auto plain = m;
    plain.quantum_confinement = assemble::QuantumConfinement::off;
    auto poisson = *create_equilibrium(d, scaling, plain);
    REQUIRE(poisson.set_bias(bias).has_value());
    const auto psi0 = charge_neutral_potential(poisson);
    std::vector<double> x0(psi0.begin(), psi0.end());
    dg_newton(poisson, x0);
    auto eq = *create_equilibrium(d, scaling, m);
    REQUIRE(eq.set_bias(bias).has_value());
    std::vector<double> x(3 * x0.size(), 0.0);
    for (std::size_t i = 0; i < x0.size(); ++i) x[3 * i] = x0[i];
    dg_newton(eq, x);
    return x;
}

}  // namespace

TEST_CASE("density gradient assemble: the finite-difference Jacobian gate (QV5)") {
    // The house gate (fd_jacobian_error): the largest relative column error. The 2D and 3D MOS
    // without the n+ pocket: across its abrupt junction the minority densities jump by 1e16, and the
    // Lambda rows' terms there (1e8) put that row's central difference at rounding (reported in
    // DECISIONS.md; the entries match the differences as the step grows).
    const auto check = [](const device::Device& d, const assemble::PhysicsModels& m,
                          const std::vector<double>& eq_bias, const std::vector<double>& bias,
                          const char* what) {
        const auto scaling = *assemble::make_scaling(d);
        auto plain = m;
        plain.quantum_confinement = assemble::QuantumConfinement::off;
        auto poisson = *create_equilibrium(d, scaling, plain);
        REQUIRE(poisson.set_bias(eq_bias).has_value());
        auto dg_poisson = *create_equilibrium(d, scaling, m);
        REQUIRE(dg_poisson.set_bias(eq_bias).has_value());
        const double ep = fd_jacobian_error(dg_poisson, dg_equilibrium_probe(poisson, 31));
        auto system = *System::create(d, scaling, m);
        REQUIRE(system.set_bias(bias).has_value());
        const double dd = fd_jacobian_error(system, dg_probe(system, poisson, d, 23));
        CAPTURE(what, ep, dd);
        CHECK(ep <= 1e-6);
        CHECK(dd <= 1e-6);
    };
    const std::vector<double> eq_bias{-0.4, 0.0}, bias{-0.4, 0.1};
    for (const int D : {1, 2, 3}) {
        for (const bool fd : {false, true}) {
            auto m = dg();
            m.fermi_dirac = fd;
            CAPTURE(D, fd);
            check(mos(D, true, false), m, eq_bias, bias, "meshed oxide, traps");
        }
    }
    const device::Device mos2 = mos(2, true, false);
    auto m = dg();
    m.surface_mobility = assemble::SurfaceMobility::lombardi;
    check(mos2, m, eq_bias, bias, "Lombardi");
    m = dg();
    m.impact_ionization = true;
    m.field_mobility = true;
    check(mos2, m, eq_bias, bias, "impact ionization, field mobility");
    m = dg();
    m.incomplete_ionization = true;
    check(mos2, m, eq_bias, bias, "incomplete ionization");
    m = dg();
    m.btbt_local = true;
    check(mos2, m, eq_bias, bias, "local band-to-band");
    for (const bool fd : {false, true}) {
        m = dg();
        m.fermi_dirac = fd;
        CAPTURE(fd);
        check(lumped(1), m, {0.8, 0.0}, {0.8, 0.05}, "lumped gate 1D");
        check(lumped(2), m, {0.8, 0.0}, {0.8, 0.05}, "lumped gate 2D");
        check(delaunay_mos(), m, {0.8, 0.0}, {0.8, 0.05}, "Delaunay");
        check(meshed_dg(physics::silicon_dioxide_tunnelling_parameters), m, {0.5, 0.0}, {0.5, 0.05},
              "meshed oxide with band data");
    }
}

TEST_CASE("density gradient assemble: the Lambda rows against analytic profiles (QV3)") {
    const assemble::PhysicsModels models = dg();
    const double beta_n = physics::dg_coefficient_cm2(1.0, models.dg_mass_n, 300.0);
    const double beta_p = physics::dg_coefficient_cm2(1.0, models.dg_mass_p, 300.0);
    SECTION("interior: exp and cosh, O(h^2) on uniform and graded meshes, 2D and 3D the 1D bar") {
        // sqrt c = exp(-x / 2L): Lambda = -beta / 4L^2; sqrt c = cosh((x - 1 um) / L): -beta / L^2.
        constexpr double L_exp = 0.1e-4, L_cosh = 0.2e-4;
        const auto r_exp = [](double x) { return std::exp(-x / (2.0 * L_exp)); };
        const auto r_cosh = [](double x) { return std::cosh((x - 1e-4) / L_cosh); };
        for (const bool graded : {false, true}) {
            for (const bool cosh_profile : {false, true}) {
                const double k2 = cosh_profile ? 1.0 / (L_cosh * L_cosh) : 0.25 / (L_exp * L_exp);
                const double en = -beta_n * k2, ep = -beta_p * k2;
                std::vector<double> errors;
                for (const int cells : {20, 40, 80, 160}) {
                    const auto axis = dg_axis(cells, graded);
                    const auto profile = [&](const System& s, const device::Device& d) {
                        return cosh_profile ? dg_profile(s, d, r_cosh, r_cosh)
                                            : dg_profile(s, d, r_exp, r_exp);
                    };
                    const device::Device d = dg_bar(1, axis);
                    const auto system = *System::create(d, *assemble::make_scaling(d), models);
                    const auto [ln, lp] = dg_lambda(system, profile(system, d));
                    double worst = 0.0;
                    for (std::size_t i = 1; i + 1 < ln.size(); ++i) {
                        worst = std::max({worst, std::abs(ln[i] / en - 1.0), std::abs(lp[i] / ep - 1.0)});
                    }
                    CHECK(ln.front() == 0.0);
                    CHECK(lp.back() == 0.0);
                    errors.push_back(worst);
                    for (const int D : {2, 3}) {
                        const device::Device e = dg_bar(D, axis);
                        const auto s = *System::create(e, *assemble::make_scaling(e), models);
                        const auto [an, ap] = dg_lambda(s, profile(s, e));
                        double diff = 0.0;
                        for (std::size_t i = 0; i < an.size(); ++i) {
                            const std::size_t k = i % axis.size();
                            REQUIRE(e.mesh().points()[i][0] == d.mesh().points()[k][0]);
                            diff = std::max({diff, std::abs(an[i] - ln[k]) / std::abs(en),
                                             std::abs(ap[i] - lp[k]) / std::abs(ep)});
                        }
                        CAPTURE(graded, cosh_profile, D, cells, diff);
                        CHECK(diff <= 1e-12);
                    }
                }
                for (std::size_t k = 1; k < errors.size(); ++k) {
                    const double order = std::log2(errors[k - 1] / errors[k]);
                    CAPTURE(graded, cosh_profile, k, errors[k - 1], errors[k], order);
                    CHECK(order >= 1.9);
                    CHECK(order <= 2.1);
                }
            }
        }
    }
    SECTION("barriers: a meshed oxide and a lumped gate converge to the wall's closed form") {
        // In front of a wall with decay length x_p, sqrt c = sinh((x + x_1) / L) with coth(x_1 / L)
        // = L / x_p meets the wall's condition (sqrt c)' / sqrt c = 1 / x_p at x = 0 (the decay
        // exp(x / x_p) continued past the face) and gives Lambda = -beta / L^2 at every node.
        const physics::Semiconductor si = physics::silicon();
        const double chi = si.parameters().electron_affinity_eV;
        const double eg = physics::band_gap_eV(si, 300.0);
        const auto oxide = *assemble::detail::insulator_band_data(physics::silicon_dioxide_parameters);
        const double xpn = physics::barrier_decay_length_cm(oxide.tunnel_mass_n,
                                                            chi - oxide.electron_affinity_eV);
        const double xpp = physics::barrier_decay_length_cm(
            oxide.tunnel_mass_p, oxide.electron_affinity_eV + oxide.band_gap_eV - chi - eg);
        constexpr double L = 1e-7;
        const double x1n = L * std::atanh(xpn / L), x1p = L * std::atanh(xpp / L);
        const auto rn = [&](double x) { return std::sinh((x + x1n) / L); };
        const auto rp = [&](double x) { return std::sinh((x + x1p) / L); };
        const double en = -beta_n / (L * L), ep = -beta_p / (L * L);
        // The lumped stack's band data is SiO2's: the same wall.
        {
            const device::Device d = dg_wall(0.1e-7, true);
            const auto b = assemble::detail::gate_band_data(d.contacts()[0].gate);
            CHECK(b.electron_affinity_eV == oxide.electron_affinity_eV);
            CHECK(b.band_gap_eV == oxide.band_gap_eV);
            CHECK(b.tunnel_mass_n == oxide.tunnel_mass_n);
            CHECK(b.tunnel_mass_p == oxide.tunnel_mass_p);
        }
        // The wall node's signed errors per spacing, [case][carrier][spacing]; case 0 meshed, 1 lumped.
        // Measured (n; p alike): meshed -0.027, -0.036, -0.029, -0.019; lumped 0.86, 0.42, 0.21, 0.105
        // (first order, h / 4 x_p: the node's half box against the barrier flux).
        const double spacings[] = {0.4e-7, 0.2e-7, 0.1e-7, 0.05e-7};
        double wall_error[2][2][4] = {}, interior_error[2][4] = {};
        for (std::size_t k = 0; k < 4; ++k) {
            const double h = spacings[k];
            double interior[2] = {};
            for (const bool lumped_gate : {false, true}) {
                const device::Device d = dg_wall(h, lumped_gate);
                const auto system = *System::create(d, *assemble::make_scaling(d), models);
                const auto [ln, lp] = dg_lambda(system, dg_profile(system, d, rn, rp));
                const auto& pts = d.mesh().points();
                std::size_t wall = pts.size();
                for (std::size_t i = 0; i < pts.size(); ++i) {
                    if (pts[i][0] >= 0.0 && (wall == pts.size() || pts[i][0] < pts[wall][0])) wall = i;
                }
                const std::size_t c = lumped_gate ? 1 : 0;
                wall_error[c][0][k] = ln[wall] / en - 1.0;
                wall_error[c][1][k] = lp[wall] / ep - 1.0;
                for (std::size_t i = 0; i + 1 < pts.size(); ++i) {  // the contact node is held
                    if (i == wall || pts[i][0] < 0.0) continue;
                    interior[c] = std::max({interior[c], std::abs(ln[i] / en - 1.0), std::abs(lp[i] / ep - 1.0)});
                }
            }
            for (std::size_t c = 0; c < 2; ++c) interior_error[c][k] = interior[c];
        }
        for (std::size_t c = 0; c < 2; ++c) {
            for (std::size_t q = 0; q < 2; ++q) {
                // The rest of the silicon is the interior case: O(h^2).
                if (q == 0) {
                    for (std::size_t k = 1; k < 4; ++k) {
                        const double order = std::log2(interior_error[c][k - 1] / interior_error[c][k]);
                        CAPTURE(c, k, order);
                        CHECK(order >= 1.9);
                        CHECK(order <= 2.1);
                    }
                }
                const double* e = wall_error[c][q];
                CAPTURE(c, q, e[0], e[1], e[2], e[3]);
                // Converging once h / 2 is below x_p (0.2 nm and finer).
                CHECK(std::abs(e[2]) < std::abs(e[1]));
                CHECK(std::abs(e[3]) < std::abs(e[2]));
                // No jump: extrapolated at first order, both cases reach the closed form.
                CHECK(std::abs(2.0 * e[3] - e[2]) <= 0.02);
            }
        }
    }
}

TEST_CASE("density gradient assemble: equilibrium slaving, no currents, no net rates (QV4)") {
    struct Case {
        int which;  // 0: the meshed 1D MOS-C, 1: the lumped 1D MOS-C, 2: the 2D MOS with its n+ pocket
        double vg;
    };
    for (const Case c : {Case{0, -1.5}, Case{0, 1.0}, Case{1, -1.5}, Case{1, 1.0}, Case{2, 1.0}}) {
        for (const bool fd : {false, true}) {
            CAPTURE(c.which, c.vg, fd);
            const device::Device d = c.which == 0 ? mos(1) : c.which == 1 ? lumped(1) : mos(2);
            const auto scaling = *assemble::make_scaling(d);
            auto m = dg();
            m.fermi_dirac = fd;
            const std::vector<double> bias{c.vg, 0.0};
            const auto x = dg_equilibrium(d, m, bias);
            const std::size_t n = d.mesh().node_count();
            std::vector<double> psi(n), ln(n), lp(n);
            for (std::size_t i = 0; i < n; ++i) {
                psi[i] = x[3 * i];
                ln[i] = x[3 * i + 1];
                lp[i] = x[3 * i + 2];
            }
            auto system = *System::create(d, scaling, m);
            REQUIRE(system.set_bias(bias).has_value());
            const auto state = system.state_from_quantum_potential(psi, ln, lp);
            // QP8's relation (Boltzmann): the classical density at psi times e^-Lambda.
            if (!fd) {
                const auto classical = system.state_from_potential(psi);
                double worst = 0.0;
                for (std::size_t i = 0; i < n; ++i) {
                    if (d.is_insulator(static_cast<mesh::NodeId>(i))) continue;
                    worst = std::max({worst,
                                      std::abs(state[dg_m * i + 1] / (classical[dg_m * i + 1] * std::exp(-ln[i])) - 1.0),
                                      std::abs(state[dg_m * i + 2] / (classical[dg_m * i + 2] * std::exp(-lp[i])) - 1.0)});
                }
                CAPTURE(worst);
                CHECK(worst <= 1e-12);
            }
            // The drift-diffusion residual on the slaved state is rounding (every flux, net rate,
            // trap term and Lambda row): each row against its largest term |J_rc x_c|.
            std::vector<double> f(system.unknowns());
            auto j = system.make_jacobian();
            system.evaluate(state, f, j);
            double worst = 0.0;
            for (std::size_t r = 0; r < f.size(); ++r) {
                double scale = 1e-300;
                for (auto k = static_cast<std::size_t>(j.row_offsets()[r]);
                     k < static_cast<std::size_t>(j.row_offsets()[r + 1]); ++k) {
                    scale = std::max(
                        scale, std::abs(j.values()[k] * state[static_cast<std::size_t>(j.col_indices()[k])]));
                }
                worst = std::max(worst, std::abs(f[r]) / scale);
            }
            CAPTURE(worst);
            CHECK(worst <= 1e-10);
            // Newton from it does not move.
            auto y = state;
            const auto again = dg_newton(system, y, c.which == 2);
            REQUIRE_FALSE(again.iterations.empty());
            CAPTURE(again.iterations.front().update);
            CHECK(again.iterations.front().update <= 1e-9);
        }
    }
}

TEST_CASE("density gradient assemble: a y-uniform 2D MOS-C is the 1D one, a 3D extrusion the 2D (QV9)") {
    for (const bool fd : {false, true}) {
        for (const double vg : {-1.5, 1.5}) {
            auto m = dg();
            m.fermi_dirac = fd;
            const device::Device d1 = lumped(1), d2 = lumped(2), d3 = lumped(3);
            const auto x1 = dg_equilibrium(d1, m, {vg, 0.0});
            const auto x2 = dg_equilibrium(d2, m, {vg, 0.0});
            const auto x3 = dg_equilibrium(d3, m, {vg, 0.0});
            const std::size_t nx = d1.mesh().node_count(), n2 = d2.mesh().node_count();
            double peak = 0.0, e21 = 0.0, e32 = 0.0;
            for (std::size_t i = 0; i < nx; ++i) peak = std::max(peak, std::abs(x1[3 * i + 1]));
            const auto off = [](double a, double b) { return std::abs(a - b) / std::max(std::abs(b), 1.0); };
            for (std::size_t i = 0; i < n2; ++i) {
                const std::size_t k = i % nx;
                REQUIRE(d2.mesh().points()[i][0] == d1.mesh().points()[k][0]);
                for (std::size_t c = 0; c < 3; ++c) e21 = std::max(e21, off(x2[3 * i + c], x1[3 * k + c]));
            }
            for (std::size_t i = 0; i < d3.mesh().node_count(); ++i) {
                const std::size_t k = i % n2;
                REQUIRE(d3.mesh().points()[i][0] == d2.mesh().points()[k][0]);
                REQUIRE(d3.mesh().points()[i][1] == d2.mesh().points()[k][1]);
                for (std::size_t c = 0; c < 3; ++c) e32 = std::max(e32, off(x3[3 * i + c], x2[3 * k + c]));
            }
            CAPTURE(fd, vg, peak, e21, e32);
            CHECK(peak > 0.1);
            CHECK(e21 <= 1e-12);
            CHECK(e32 <= 1e-12);
        }
    }
}

TEST_CASE("density gradient assemble: drift-diffusion Newton is quadratic on the MOS (QV5)") {
    // The 2D MOS (with its pocket) in strong inversion from the classical equilibrium at 0.8 V,
    // Lambda = 0 (measured: 23 iterations, the last updates 6.4e-3, 2.9e-5, 5.9e-10, then 2.8e-12 at
    // rounding; without the model 6).
    const device::Device d = mos(2);
    const auto scaling = *assemble::make_scaling(d);
    auto poisson = *create_equilibrium(d, scaling);
    REQUIRE(poisson.set_bias(std::vector<double>{0.8, 0.0}).has_value());
    const auto psi0 = charge_neutral_potential(poisson);
    std::vector<double> x0(psi0.begin(), psi0.end());
    dg_newton(poisson, x0);
    auto system = *System::create(d, scaling, dg());
    REQUIRE(system.set_bias(std::vector<double>{1.0, 0.0}).has_value());
    auto x = system.state_from_potential(x0);
    const auto record = dg_newton(system, x);
    const auto& it = record.iterations;
    std::string history;
    for (const auto& r : it) {
        char buffer[32];
        std::snprintf(buffer, sizeof buffer, "%.3g ", r.update);
        history += buffer;
    }
    CAPTURE(history);
    REQUIRE(it.size() >= 3);
    // Above rounding each update is at most C times the previous squared: the order
    // log(u_{k+1} / u_k) / log(u_k / u_{k-1}) >= 1.8 where u_k < 1e-2 and u_{k+1} > 1e-10.
    std::size_t gated = 0;
    for (std::size_t k = 1; k + 1 < it.size(); ++k) {
        const double u0 = it[k - 1].update, u1 = it[k].update, u2 = it[k + 1].update;
        if (!(u1 < 1e-2 && u2 > 1e-10 && u0 > u1)) continue;
        const double order = std::log(u2 / u1) / std::log(u1 / u0);
        CAPTURE(k, order);
        CHECK(order >= 1.8);
        ++gated;
    }
    CHECK(gated >= 1);
}

TEST_CASE("density gradient assemble: the time step's Jacobian and the small-signal matrix (QP11)") {
    // The step's Jacobian against central differences (the house gate's measure, the step's
    // residual), at two rates; the small-signal matrix at s = rate (real) is that Jacobian with the
    // steady state as history (backward Euler is exact for e^(s t)), and at s = i w, without traps,
    // the steady Jacobian plus i w times the storage's.
    struct Case {
        const char* what;
        device::Device d;
        bool ionization;
    };
    std::vector<Case> cases;
    cases.push_back({"meshed 1D, traps", mos(1), false});
    cases.push_back({"meshed 2D, traps", mos(2, true, false), false});
    cases.push_back({"meshed 2D, no traps, incomplete ionization", mos(2, false, false), true});
    cases.push_back({"lumped 1D", lumped(1), false});
    for (const Case& c : cases) {
        for (const bool fd : {false, true}) {
            CAPTURE(c.what, fd);
            auto m = dg();
            m.fermi_dirac = fd;
            m.incomplete_ionization = c.ionization;
            const auto scaling = *assemble::make_scaling(c.d);
            auto plain = m;
            plain.quantum_confinement = assemble::QuantumConfinement::off;
            const bool meshed = c.d.contacts()[0].kind != device::ContactKind::gate;
            const std::vector<double> eq_bias = meshed ? std::vector<double>{-0.4, 0.0} : std::vector<double>{0.8, 0.0};
            const std::vector<double> bias = meshed ? std::vector<double>{-0.4, 0.1} : std::vector<double>{0.8, 0.05};
            auto poisson = *create_equilibrium(c.d, scaling, plain);
            REQUIRE(poisson.set_bias(eq_bias).has_value());
            auto dd = *System::create(c.d, scaling, m);
            REQUIRE(dd.set_bias(bias).has_value());
            std::vector<double> x = dg_probe(dd, poisson, c.d, 41);
            const std::vector<double> history = dd.storage(x);
            const std::vector<double> traps = assemble::trap_occupancies(dd, x);
            const std::size_t n = dd.unknowns();
            std::vector<double> f(n), fp(n), fm(n);
            for (const double factor : {1.0, 100.0}) {
                const System::TimeStep step{factor * dd.time_scale() / 0.35e-12, history, traps};
                auto jac = dd.make_jacobian();
                dd.evaluate(x, step, f, jac);
                const auto j = dense(jac);
                double worst = 0.0;
                for (std::size_t col = 0; col < n; ++col) {
                    const double base = x[col];
                    const double h = 1e-7 * std::max(std::abs(base), 1.0);
                    x[col] = base + h;
                    const double up = x[col];
                    dd.residual(x, step, fp);
                    x[col] = base - h;
                    const double down = x[col];
                    dd.residual(x, step, fm);
                    x[col] = base;
                    double scale = 1e-300, error = 0.0;
                    for (std::size_t r = 0; r < n; ++r) scale = std::max(scale, std::abs(j[r][col]));
                    for (std::size_t r = 0; r < n; ++r) {
                        error = std::max(error, std::abs((fp[r] - fm[r]) / (up - down) - j[r][col]));
                    }
                    worst = std::max(worst, error / scale);
                }
                CAPTURE(factor, worst);
                CHECK(worst <= 1e-6);
            }
            const auto jacobian_at = [&](double rate) {
                auto j = dd.make_jacobian();
                if (rate == 0.0) {
                    dd.evaluate(x, f, j);
                } else {
                    dd.evaluate(x, System::TimeStep{rate, history, traps}, f, j);
                }
                return std::vector<double>(j.values().begin(), j.values().end());
            };
            const double rate = 3.0 * dd.time_scale() / 0.35e-12;
            const std::vector<double> steady = jacobian_at(0.0), stepped = jacobian_at(rate), unit = jacobian_at(1.0);
            auto a = dd.make_small_signal_matrix();
            double scale = 0.0, worst = 0.0;
            for (const double v : stepped) scale = std::max(scale, std::abs(v));
            dd.small_signal_matrix(x, {rate, 0.0}, a);
            for (std::size_t k = 0; k < stepped.size(); ++k) {
                worst = std::max(worst, std::abs(a.values()[k] - stepped[k]));
            }
            CAPTURE(worst, scale);
            CHECK(worst <= 1e-12 * scale);
            if (!traps.empty()) continue;
            const double w = 0.7 * rate;
            dd.small_signal_matrix(x, {0.0, w}, a);
            double storage = 0.0;
            worst = 0.0;
            for (std::size_t k = 0; k < stepped.size(); ++k) {
                const std::complex<double> expected{steady[k], w * (unit[k] - steady[k])};
                worst = std::max(worst, std::abs(a.values()[k] - expected));
                storage = std::max(storage, std::abs(unit[k] - steady[k]));
            }
            CAPTURE(worst, storage);
            CHECK(storage > 0.0);
            CHECK(worst <= 1e-12 * scale);
        }
    }
}
