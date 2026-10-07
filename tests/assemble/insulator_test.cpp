// Meshed insulators, electrodes and interface charge, traps and recombination in the assemblers
// (ARCHITECTURE.md section 11, Unit 15b): the finite-difference Jacobian gate in 1D, 2D and 3D
// under both statistics, the insulator rows, no carrier flux across a semiconductor-insulator edge,
// and the drift-diffusion interface terms reducing to the equilibrium ones (Fermi occupancy, no
// recombination) on an equilibrium state.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/equilibrium_poisson.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/interface_traps.hpp"
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

// x: oxide on [-5 nm, 0) with the electrode on x_min, p-silicon (1e17) on (0, 120 nm] with an
// n+ pocket (1e19) at the top end of the transverse axes and the substrate contact on x_max. The
// interface sits midway between the nodes at -+0.5 nm. In 2D and 3D the electrode is the part of
// x_min with y below the middle (the rest of the oxide surface is free, Neumann).
device::Device mos(int dimension, bool traps = true, double transverse_nodes = 4) {
    const auto x = *mesh::straddle_interface(uniform(-5e-7, 1.2e-5, 30), 0.0, 1e-7);
    const auto t = uniform(0.0, 3e-6, static_cast<int>(transverse_nodes));
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
        if (dimension > 1 && p[1] > 2e-6 && p[0] < 3e-6) donors[i] = 1e19;
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
        f.traps.bands = {{.type = physics::TrapType::acceptor, .density_cm2_eV = 5e11,
                          .energy_low_eV = -0.1, .energy_high_eV = 0.4}};
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

// The house gate (legacy _jacobian_probe), every column, central differences with steps
// 1e-7 max(|x|, 1).
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

// A drift-diffusion probe state: the equilibrium state of the charge-neutral guess, psi + 0.02
// noise, densities times (1 + 0.01 noise), and at the interface nodes n = 2e-3 and p = 5e-4
// (scaled): out of equilibrium, so the trap occupancy is not the Fermi one, and large enough that a
// difference step resolves the trap terms against the rest of the row (at the bulk minority
// density, about 1e-14 scaled, no step can).
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

TEST_CASE("insulator assemble: the finite-difference Jacobian gate, 1D/2D/3D, both statistics") {
    const std::vector<double> bias{-0.4, 0.1};
    const std::vector<double> equilibrium_bias{-0.4, 0.0};
    for (const int D : {1, 2, 3}) {
        const device::Device d = mos(D);
        const auto scaling = *assemble::make_scaling(d);
        for (const bool fd : {false, true}) {
            const assemble::PhysicsModels models{.fermi_dirac = fd};
            auto poisson = *EquilibriumPoisson::create(d, scaling, models);
            REQUIRE(poisson.set_bias(equilibrium_bias).has_value());
            auto psi = poisson.charge_neutral_potential();
            Noise noise{17};
            for (double& v : psi) v += 0.3 * noise.next();
            const double ep = fd_jacobian_error(poisson, psi);
            auto system = *DriftDiffusion::create(d, scaling, models);
            REQUIRE(system.set_bias(bias).has_value());
            const auto x = probe_state(system, poisson, d, 23);
            const double dd = fd_jacobian_error(system, x);
            CAPTURE(D, fd, ep, dd);
            // Measured at most 6.9e-8 (Poisson) and 7.0e-8 (drift-diffusion), in 2D and 3D; the
            // house tolerance elsewhere is 5e-5.
            REQUIRE(ep <= 1e-6);
            REQUIRE(dd <= 1e-6);
        }
    }
}

TEST_CASE("insulator assemble: the trap part of the Jacobian on its own") {
    // The same probe with and without the interface terms; their difference isolates the trap and
    // recombination columns, so a small error there cannot hide under the large flux entries.
    const device::Device with = mos(1, true), without = mos(1, false);
    const auto scaling = *assemble::make_scaling(with);
    const std::vector<double> bias{-0.4, 0.1};
    for (const bool fd : {false, true}) {
        const assemble::PhysicsModels models{.fermi_dirac = fd};
        auto a = *DriftDiffusion::create(with, scaling, models);
        auto b = *DriftDiffusion::create(without, scaling, models);
        REQUIRE(a.set_bias(bias).has_value());
        REQUIRE(b.set_bias(bias).has_value());
        const auto poisson = *EquilibriumPoisson::create(with, scaling, models);
        auto x = probe_state(a, poisson, with, 5);
        const std::size_t n = a.unknowns();
        auto ja = a.make_jacobian(), jb = b.make_jacobian();
        std::vector<double> fa(n), fb(n), up(n), down(n), up_b(n), down_b(n);
        a.evaluate(x, fa, ja);
        b.evaluate(x, fb, jb);
        const auto da = dense(ja), db = dense(jb);
        double worst = 0.0, largest = 0.0;
        for (std::size_t c = 0; c < n; ++c) {
            const double base = x[c];
            const double step = 1e-7 * std::max(std::abs(base), 1.0);
            x[c] = base + step;
            const double hi = x[c];
            a.residual(x, up);
            b.residual(x, up_b);
            x[c] = base - step;
            const double lo = x[c];
            a.residual(x, down);
            b.residual(x, down_b);
            x[c] = base;
            for (std::size_t r = 0; r < n; ++r) {
                const double analytic = da[r][c] - db[r][c];
                const double numeric =
                    ((up[r] - up_b[r]) - (down[r] - down_b[r])) / (hi - lo);
                largest = std::max(largest, std::abs(analytic));
                worst = std::max(worst, std::abs(numeric - analytic));
            }
        }
        CAPTURE(fd, worst, largest, worst / largest);
        REQUIRE(largest > 0.0);
        REQUIRE(worst <= 1e-6 * largest);  // measured 2.4e-9
    }
}

TEST_CASE("insulator assemble: insulator rows and no carrier flux across the oxide") {
    for (const int D : {1, 2}) {
        const device::Device d = mos(D);
        const auto scaling = *assemble::make_scaling(d);
        auto system = *DriftDiffusion::create(d, scaling);
        REQUIRE(system.set_bias(std::vector<double>{0.7, 0.0}).has_value());
        const auto poisson = *EquilibriumPoisson::create(d, scaling);
        auto x = probe_state(system, poisson, d, 9);
        const std::size_t n = system.unknowns();
        std::vector<double> f(n);
        auto jacobian = system.make_jacobian();
        system.evaluate(x, f, jacobian);
        const auto j = dense(jacobian);
        for (std::size_t i = 0; i < system.node_count(); ++i) {
            if (!d.is_insulator(static_cast<mesh::NodeId>(i))) continue;
            // n = p = 0 rows, decoupled from everything else.
            REQUIRE(x[3 * i + 1] == 0.0);
            REQUIRE(x[3 * i + 2] == 0.0);
            REQUIRE(f[3 * i + 1] == 0.0);
            REQUIRE(f[3 * i + 2] == 0.0);
            for (std::size_t c = 0; c < n; ++c) {
                REQUIRE(j[3 * i + 1][c] == (c == 3 * i + 1 ? 1.0 : 0.0));
                REQUIRE(j[3 * i + 2][c] == (c == 3 * i + 2 ? 1.0 : 0.0));
                // No semiconductor row reads an insulator density.
                if (!d.is_insulator(static_cast<mesh::NodeId>(c / 3))) {
                    REQUIRE(j[c][3 * i + 1] == 0.0);
                    REQUIRE(j[c][3 * i + 2] == 0.0);
                }
            }
            // An electrode node holds psi: its row is psi - psi_E.
            const bool electrode =
                d.contacts()[0].nodes.end() != std::find(d.contacts()[0].nodes.begin(),
                                                         d.contacts()[0].nodes.end(),
                                                         static_cast<mesh::NodeId>(i));
            if (electrode) {
                for (std::size_t c = 0; c < n; ++c) {
                    REQUIRE(j[3 * i][c] == (c == 3 * i ? 1.0 : 0.0));
                }
            }
        }
        // Zero electron and hole current on every edge with an insulator end; the gate carries no
        // terminal current.
        const auto currents = system.edge_currents(x);
        const auto edges = d.mesh().edges();
        for (std::size_t k = 0; k < edges.size(); ++k) {
            if (d.is_insulator(edges[k].first) || d.is_insulator(edges[k].second)) {
                REQUIRE(currents[k].first == 0.0);
                REQUIRE(currents[k].second == 0.0);
            }
        }
        REQUIRE(system.terminal_currents(x)[0] == 0.0);
        REQUIRE(system.terminal_current_resolution(x)[0] == 0.0);
        const auto bands = system.band_edges(x);
        REQUIRE(std::isnan(bands.conduction[0]));
    }
}

TEST_CASE("insulator assemble: on an equilibrium state the interface terms are the Fermi ones") {
    // At psi with n and p slaved to it (equilibrium densities, gate biased, substrate at 0 V), the
    // drift-diffusion Poisson rows equal the equilibrium ones (the SRH occupancy is the Fermi
    // function) and the interface recombination vanishes: the continuity rows hold only the
    // fluxes, which vanish too.
    for (const int D : {1, 2, 3}) {
        for (const bool fd : {false, true}) {
            const device::Device d = mos(D);
            const assemble::PhysicsModels models{.fermi_dirac = fd};
            const auto scaling = *assemble::make_scaling(d);
            auto poisson = *EquilibriumPoisson::create(d, scaling, models);
            auto system = *DriftDiffusion::create(d, scaling, models);
            const std::vector<double> bias{-0.4, 0.0};
            REQUIRE(poisson.set_bias(bias).has_value());
            REQUIRE(system.set_bias(bias).has_value());
            auto psi = poisson.charge_neutral_potential();
            Noise noise{3};
            for (std::size_t i = 0; i < psi.size(); ++i) {
                if (poisson.is_contact()[i] == 0) psi[i] += 0.5 * noise.next();
            }
            const auto x = system.state_from_potential(psi);
            std::vector<double> fe(poisson.unknowns()), fd_rows(system.unknowns());
            poisson.residual(psi, fe);
            system.residual(x, fd_rows);
            double worst = 0.0, worst_rate = 0.0, scale = 0.0;
            for (std::size_t i = 0; i < psi.size(); ++i) {
                scale = std::max(scale, std::abs(fe[i]));
                worst = std::max(worst, std::abs(fd_rows[3 * i] - fe[i]));
                worst_rate = std::max(
                    {worst_rate, std::abs(fd_rows[3 * i + 1]), std::abs(fd_rows[3 * i + 2])});
            }
            const auto qe = poisson.interface_trap_charges(psi);
            const auto qd = system.interface_trap_charges(x);
            CAPTURE(D, fd, worst, scale, worst_rate, qe, qd);
            REQUIRE(qe.size() == 1);
            REQUIRE(qe[0] != 0.0);
            // Measured: the trapped charges within 2.8e-16 relative, the Poisson rows equal to
            // rounding, the continuity rows at most 1.2e-14 (rows of magnitude about 1e2).
            REQUIRE(std::abs(qd[0] - qe[0]) <= 1e-14 * std::abs(qe[0]));
            REQUIRE(worst <= 1e-12 * scale);
            REQUIRE(worst_rate <= 1e-12);
        }
    }
}

TEST_CASE("insulator assemble: a mid-gap trap recombines as its surface recombination velocity") {
    // A mid-gap level with N sigma_n v_n = s_n and N sigma_p v_p = s_p recombines as the velocity
    // form. The trap also holds a charge, which moves the interface potential and so the interface
    // densities; with N small (and sigma large) that charge is negligible, and the continuity rows
    // of every node equal those of the velocity form at any state.
    enum class Mode { traps, velocity, none };
    const auto make = [](Mode mode) {
        const auto x = *mesh::straddle_interface(uniform(-5e-7, 1.2e-5, 30), 0.0, 1e-7);
        mesh::Mesh m = *mesh::make_tensor_grid(x);
        const std::size_t n = m.node_count();
        std::vector<device::RegionId> region(n);
        std::vector<double> acceptors(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            region[i] = m.points()[i][0] < 0.0 ? 0 : 1;
            if (region[i] == 1) acceptors[i] = 1e16;
        }
        device::Interface f{"oxide", "silicon"};
        if (mode == Mode::traps) {
            f.traps.levels = {{.type = physics::TrapType::acceptor, .density_cm2 = 2e4,
                               .energy_eV = 0.0, .sigma_n_cm2 = 3e-9, .sigma_p_cm2 = 1e-9}};
            f.traps.thermal_velocity_n_cm_s = 2e7;
            f.traps.thermal_velocity_p_cm_s = 1.5e7;
        } else if (mode == Mode::velocity) {
            f.recombination_velocity_n_cm_s = 2e10 * 3e-15 * 2e7;
            f.recombination_velocity_p_cm_s = 2e10 * 1e-15 * 1.5e7;
        }
        auto x_max = m.find_boundary("x_max")->nodes;
        return *device::Device::create(
            {.mesh = std::move(m),
             .temperature_K = 300.0,
             .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
             .node_region = std::move(region),
             .donors = std::vector<double>(n, 0.0),
             .acceptors = std::move(acceptors),
             .contacts = {{"gate", device::ContactKind::electrode, {0}},
                          {"substrate", device::ContactKind::ohmic, std::move(x_max)}},
             .interfaces = {std::move(f)}});
    };
    const device::Device a = make(Mode::traps), b = make(Mode::velocity), c = make(Mode::none);
    const auto scaling = *assemble::make_scaling(a);
    auto sa = *DriftDiffusion::create(a, scaling);
    auto sb = *DriftDiffusion::create(b, scaling);
    auto sc = *DriftDiffusion::create(c, scaling);
    const auto poisson = *EquilibriumPoisson::create(a, scaling);
    const auto x = probe_state(sa, poisson, a, 41);
    std::vector<double> fa(sa.unknowns()), fb(sb.unknowns()), fc(sc.unknowns());
    sa.residual(x, fa);
    sb.residual(x, fb);
    sc.residual(x, fc);
    // The recombination term is the continuity rows' difference from the device without it;
    // the two forms are compared on its own scale.
    double worst = 0.0, rate = 0.0, poisson_difference = 0.0;
    for (std::size_t i = 0; i < sa.node_count(); ++i) {
        for (const std::size_t r : {3 * i + 1, 3 * i + 2}) {
            worst = std::max(worst, std::abs(fa[r] - fb[r]));
            rate = std::max(rate, std::abs(fb[r] - fc[r]));
        }
        poisson_difference = std::max(poisson_difference, std::abs(fa[3 * i] - fb[3 * i]));
    }
    CAPTURE(worst, rate, worst / rate, poisson_difference);
    REQUIRE(rate > 0.0);
    REQUIRE(worst <= 1e-9 * rate);
    REQUIRE(poisson_difference > 0.0);
}
