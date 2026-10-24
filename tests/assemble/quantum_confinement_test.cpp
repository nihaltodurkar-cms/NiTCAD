// MLDA quantum confinement in the assemblers (Unit 47; DECISIONS.md QD1-QD7, QC3, QC4, QC6,
// QC11): the distance to the nearest interface on tensor, corner, Delaunay and lumped-gate devices;
// the densities of states times the factors with the band edges and the trap levels in place;
// the finite-difference Jacobian gate with the model on; and the refusals.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/equilibrium_poisson.hpp"
#include "NiTCAD/assemble/quantum_confinement.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/interface_traps.hpp"
#include "NiTCAD/physics/quantum_confinement.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD;
using assemble::DriftDiffusion;
using assemble::EquilibriumPoisson;

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
device::Device mos(int dimension, bool traps = true) {
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
        if (dimension > 1 && p[1] > 2e-6 && p[0] > 2e-6) donors[i] = 1e19;
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
// the gate face x_min, a 2 nm n+ poly gate on x_min and the substrate contact on x_max.
constexpr double lumped_h = 0.2e-7;
device::Device lumped(int dimension) {
    std::vector<double> x;
    for (int k = 0; k < 10; ++k) x.push_back(k * lumped_h);
    for (const double v : uniform(2e-7, 3e-6, 30)) x.push_back(v);
    const auto t = uniform(0.0, 3e-6, 4);
    mesh::Mesh m = dimension == 1 ? *mesh::make_tensor_grid(x) : *mesh::make_tensor_grid(x, t);
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
std::vector<double> probe_state(const DriftDiffusion& system, const EquilibriumPoisson& poisson,
                                const device::Device& d, std::uint64_t seed) {
    auto x = system.state_from_potential(poisson.charge_neutral_potential());
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
            const auto off = *EquilibriumPoisson::create(d, scaling, classical);
            const auto on = *EquilibriumPoisson::create(d, scaling, mlda(fd));
            const auto c = *assemble::mlda_factors(d, mlda(fd));
            // A potential in inversion near the oxide and in the bulk, with noise.
            auto psi = off.charge_neutral_potential();
            Noise noise{3};
            for (std::size_t i = 0; i < psi.size(); ++i) {
                if (!d.is_insulator(static_cast<mesh::NodeId>(i))) psi[i] += 0.3 * noise.next();
            }
            const std::size_t n = psi.size();
            std::vector<double> n0(n), p0(n), n1(n), p1(n);
            off.carriers(psi, n0, p0);
            on.carriers(psi, n1, p1);
            const auto b0 = off.band_edges(psi);
            const auto b1 = on.band_edges(psi);
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
            const auto q0 = off.interface_trap_charges(psi);
            const auto q1 = on.interface_trap_charges(psi);
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
    const auto off = *DriftDiffusion::create(d, scaling);
    const auto on = *DriftDiffusion::create(d, scaling, mlda());
    const auto poisson = *EquilibriumPoisson::create(d, scaling);
    const auto c = *assemble::mlda_factors(d, mlda());
    const auto psi = poisson.charge_neutral_potential();
    const auto x0 = off.state_from_potential(psi);
    const auto x1 = on.state_from_potential(psi);
    for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
        if (d.is_insulator(static_cast<mesh::NodeId>(i)) || poisson.is_contact()[i] != 0) continue;
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
            auto poisson = *EquilibriumPoisson::create(d, scaling, mlda(fd));
            REQUIRE(poisson.set_bias(equilibrium_bias).has_value());
            auto psi = poisson.charge_neutral_potential();
            Noise noise{17};
            for (double& v : psi) v += 0.3 * noise.next();
            const double ep = fd_jacobian_error(poisson, psi);
            auto system = *DriftDiffusion::create(d, scaling, mlda(fd));
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
        auto poisson = *EquilibriumPoisson::create(d, scaling, models);
        REQUIRE(poisson.set_bias(equilibrium_bias).has_value());
        auto system = *DriftDiffusion::create(d, scaling, models);
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
            auto poisson = *EquilibriumPoisson::create(d, scaling, mlda(fd));
            REQUIRE(poisson.set_bias(std::vector<double>{0.8, 0.0}).has_value());
            auto system = *DriftDiffusion::create(d, scaling, mlda(fd));
            REQUIRE(system.set_bias(std::vector<double>{0.8, 0.05}).has_value());
            auto x = system.state_from_potential(poisson.charge_neutral_potential());
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
        const auto dd = DriftDiffusion::create(dev, s, m);
        const auto ep = EquilibriumPoisson::create(dev, s, m);
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
        CHECK_FALSE(DriftDiffusion::create(d, scaling, m).has_value());
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
        CHECK(DriftDiffusion::create(ok, s, mlda()).has_value());
        CHECK(EquilibriumPoisson::create(ok, s, mlda()).has_value());
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
        CHECK(DriftDiffusion::create(d, scaling, m).has_value());
    }
}
