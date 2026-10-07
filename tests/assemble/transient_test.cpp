// Time steps in the drift-diffusion assembler (ARCHITECTURE.md section 11, Unit 21): the
// finite-difference Jacobian gate of the time-step system in 1D, 2D and 3D under both statistics
// and with incomplete ionization, the trap part on its own, the steady limit of a long step, and
// the storage and charge bookkeeping that the terminal currents are built from.
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

// The Unit 15b assembler fixture: oxide on [-5 nm, 0) with the electrode on x_min, p-silicon
// (1e17) on (0, 120 nm] with an n+ pocket (1e19, 2D and 3D) and the substrate contact on x_max;
// a fixed charge, two trap levels, a trap band and surface recombination at the interface.
device::Device mos(int dimension, bool traps = true) {
    const auto x = *mesh::straddle_interface(uniform(-5e-7, 1.2e-5, 30), 0.0, 1e-7);
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

// The time-step system as a plain system, for the difference gate.
struct Stepped {
    const DriftDiffusion& dd;
    DriftDiffusion::TimeStep step;
    [[nodiscard]] std::size_t unknowns() const { return dd.unknowns(); }
    [[nodiscard]] linalg::SparseMatrix make_jacobian() const { return dd.make_jacobian(); }
    void evaluate(std::span<const double> x, std::span<double> f, linalg::SparseMatrix& j) const {
        dd.evaluate(x, step, f, j);
    }
    void residual(std::span<const double> x, std::span<double> f) const {
        dd.residual(x, step, f);
    }
};

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

// The Unit 15b probe state: equilibrium of the charge-neutral guess with noise, and at the
// interface nodes n = 2e-3 and p = 5e-4 (scaled), so the traps are out of equilibrium and their
// terms are resolved by a difference step.
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

// A step's history: the storage of a nearby state (the probe state with 5% noise on the
// densities) and trap occupancies spread over [-0.05, 1.05] (BDF2 combinations may leave [0, 1]).
struct History {
    std::vector<double> storage, traps;
};

History history(const DriftDiffusion& system, std::vector<double> x, std::uint64_t seed) {
    Noise noise{seed};
    for (std::size_t i = 0; i < system.node_count(); ++i) {
        x[3 * i + 1] *= 1.0 + 0.05 * noise.next();
        x[3 * i + 2] *= 1.0 + 0.05 * noise.next();
    }
    History h{system.storage(x), std::vector<double>(system.trap_slots())};
    for (double& f : h.traps) f = 0.5 + 0.55 * noise.next();
    return h;
}

}  // namespace

TEST_CASE("transient assemble: the finite-difference Jacobian gate of a time step") {
    // Steps of 1 ps and 1 ns: the storage term from comparable with the flux terms to small.
    const std::vector<double> bias{-0.4, 0.1};
    for (const int D : {1, 2, 3}) {
        const device::Device d = mos(D);
        const auto scaling = *assemble::make_scaling(d);
        for (const bool fd : {false, true}) {
            for (const bool ionization : {false, true}) {
                if (D == 3 && ionization) continue;  // covered in 1D and 2D
                const assemble::PhysicsModels models{.fermi_dirac = fd,
                                                     .incomplete_ionization = ionization};
                auto system = *DriftDiffusion::create(d, scaling, models);
                REQUIRE(system.set_bias(bias).has_value());
                const auto poisson = *EquilibriumPoisson::create(d, scaling, models);
                const auto x = probe_state(system, poisson, d, 23);
                const History h = history(system, x, 31);
                REQUIRE(h.traps.size() == system.trap_slots());
                REQUIRE(system.trap_slots() > 2);
                for (const double dt : {1e-12, 1e-9}) {
                    const Stepped s{system, {system.time_scale() / dt, h.storage, h.traps}};
                    const double error = fd_jacobian_error(s, x);
                    CAPTURE(D, fd, ionization, dt, error);
                    REQUIRE(error <= 1e-6);  // measured at most 1.7e-7
                }
            }
        }
    }
}

TEST_CASE("transient assemble: the trap part of a time step's Jacobian on its own") {
    // As the Unit 15b check: with and without the interface terms, so the trap columns cannot hide
    // under the flux entries. The difference step is 1e-6 max(|x|, 1), not the house 1e-7: in a
    // 1 ps step the occupancy stays near its history, the trap entries are some 1e-5 of the
    // Poisson row's terms, and the house step's rounding (|F| eps / step) reaches 3e-6 of the
    // largest entry (it doubles from 1e-7 to 1e-8, so it is rounding, not the Jacobian).
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
        const History h = history(a, x, 9);
        const std::vector<double> no_traps;
        for (const double dt : {1e-12, 1e-9, 1e-6}) {
            const DriftDiffusion::TimeStep sa{a.time_scale() / dt, h.storage, h.traps};
            const DriftDiffusion::TimeStep sb{b.time_scale() / dt, h.storage, no_traps};
            const std::size_t n = a.unknowns();
            auto ja = a.make_jacobian(), jb = b.make_jacobian();
            std::vector<double> fa(n), fb(n), up(n), down(n), up_b(n), down_b(n);
            a.evaluate(x, sa, fa, ja);
            b.evaluate(x, sb, fb, jb);
            const auto da = dense(ja), db = dense(jb);
            double worst = 0.0, largest = 0.0;
            for (std::size_t c = 0; c < n; ++c) {
                const double base = x[c];
                const double step = 1e-6 * std::max(std::abs(base), 1.0);
                x[c] = base + step;
                const double hi = x[c];
                a.residual(x, sa, up);
                b.residual(x, sb, up_b);
                x[c] = base - step;
                const double lo = x[c];
                a.residual(x, sa, down);
                b.residual(x, sb, down_b);
                x[c] = base;
                for (std::size_t r = 0; r < n; ++r) {
                    const double analytic = da[r][c] - db[r][c];
                    const double numeric =
                        ((up[r] - up_b[r]) - (down[r] - down_b[r])) / (hi - lo);
                    largest = std::max(largest, std::abs(analytic));
                    worst = std::max(worst, std::abs(numeric - analytic));
                }
            }
            CAPTURE(fd, dt, worst, largest, worst / largest);
            REQUIRE(largest > 0.0);
            REQUIRE(worst <= 1e-6 * largest);  // measured at most 3.5e-7
        }
    }
}

TEST_CASE("transient assemble: a very long step is the steady state") {
    // r -> 0: the storage term vanishes and the trap occupancy tends to the steady-state one,
    // whatever the history.
    const device::Device d = mos(2);
    const auto scaling = *assemble::make_scaling(d);
    auto system = *DriftDiffusion::create(d, scaling);
    REQUIRE(system.set_bias(std::vector<double>{-0.4, 0.1}).has_value());
    const auto poisson = *EquilibriumPoisson::create(d, scaling);
    const auto x = probe_state(system, poisson, d, 41);
    const History h = history(system, x, 43);
    const std::size_t n = system.unknowns();
    std::vector<double> steady(n), stepped(n);
    system.residual(x, steady);
    const DriftDiffusion::TimeStep step{system.time_scale() / 1e12, h.storage, h.traps};
    system.residual(x, step, stepped);
    double scale = 0.0, diff = 0.0;
    for (std::size_t k = 0; k < n; ++k) {
        scale = std::max(scale, std::abs(steady[k]));
        diff = std::max(diff, std::abs(stepped[k] - steady[k]));
    }
    CAPTURE(diff, scale);
    REQUIRE(diff <= 1e-9 * scale);
    const auto f_steady = system.trap_occupancies(x);
    const auto f_step = system.trap_occupancies(x, &step);
    double worst = 0.0;
    for (std::size_t k = 0; k < f_steady.size(); ++k) {
        REQUIRE(f_steady[k] >= 0.0);
        REQUIRE(f_steady[k] <= 1.0);
        worst = std::max(worst, std::abs(f_step[k] - f_steady[k]));
    }
    CAPTURE(worst);
    REQUIRE(worst <= 1e-9);
    const auto q_steady = system.interface_trap_charges(x);
    const auto q_step = system.interface_trap_charges(x, &step);
    REQUIRE(std::abs(q_step[0] - q_steady[0]) <= 1e-9 * std::abs(q_steady[0]));
}


TEST_CASE("transient assemble: the traps' capture imbalance is the trapped charge's change") {
    // f = (c + k occ) / (1 + k D) solves (f - c) / k = occ - f D, so the hole row's capture
    // rate_p exceeds the electron row's rate by r (Q_trap(f) - Q_trap(c)), the BDF difference of
    // the trapped charge. Only the semiconductor node's continuity rows read the occupancy, so
    // between two histories the sum of the continuity rows changes by exactly that.
    const device::Device d = mos(1);
    const auto scaling = *assemble::make_scaling(d);
    auto system = *DriftDiffusion::create(d, scaling);
    REQUIRE(system.set_bias(std::vector<double>{-0.4, 0.1}).has_value());
    const auto poisson = *EquilibriumPoisson::create(d, scaling);
    const auto x = probe_state(system, poisson, d, 51);
    const History h = history(system, x, 53);
    std::vector<double> mirrored = h.traps;
    for (double& v : mirrored) v = 1.0 - v;
    // Q_trap of the occupancies c themselves: a vanishing weight keeps f = c.
    const auto q_of = [&](const std::vector<double>& c) {
        const DriftDiffusion::TimeStep frozen{1e300, h.storage, c};
        return system.interface_trap_charges(x, &frozen)[0];
    };
    for (const double dt : {1e-12, 1e-9}) {
        const DriftDiffusion::TimeStep s1{system.time_scale() / dt, h.storage, h.traps};
        const DriftDiffusion::TimeStep s2{system.time_scale() / dt, h.storage, mirrored};
        for (const auto* s : {&s1, &s2}) {
            const auto f = system.trap_occupancies(x, s);
            for (std::size_t k = 0; k < f.size(); ++k) {  // a weighted mean of c and occ / D
                REQUIRE(f[k] >= std::min(s->traps[k], 0.0) - 1e-15);
                REQUIRE(f[k] <= std::max(s->traps[k], 1.0) + 1e-15);
            }
        }
        const std::size_t n = system.unknowns();
        std::vector<double> r1(n), r2(n);
        system.residual(x, s1, r1);
        system.residual(x, s2, r2);
        double rows = 0.0;
        for (std::size_t i = 0; i < system.node_count(); ++i) {
            rows += (r1[3 * i + 1] + r1[3 * i + 2]) - (r2[3 * i + 1] + r2[3 * i + 2]);
        }
        const double dQ1 = system.interface_trap_charges(x, &s1)[0] - q_of(h.traps);
        const double dQ2 = system.interface_trap_charges(x, &s2)[0] - q_of(mirrored);
        const double expected = (dQ1 - dQ2) * s1.rate;
        CAPTURE(dt, rows, expected);
        REQUIRE(expected != 0.0);
        REQUIRE(std::abs(rows - expected) <= 1e-6 * std::abs(expected));
    }
}

TEST_CASE("transient assemble: the contact charges and the charge inside obey Gauss's law") {
    // The Poisson rows off the Dirichlet nodes sum to the charge of the contacts (the flux leaving
    // them) plus the charge inside the device (carriers, dopants, fixed and trapped interface
    // charge), at any state; at a solution both sides vanish. The displacement current rests on
    // this identity.
    const device::Device d = mos(2);
    const auto scaling = *assemble::make_scaling(d);
    auto system = *DriftDiffusion::create(d, scaling);
    REQUIRE(system.set_bias(std::vector<double>{-0.4, 0.1}).has_value());
    const auto poisson = *EquilibriumPoisson::create(d, scaling);
    const auto x = probe_state(system, poisson, d, 61);
    const auto q = system.contact_charges(x);
    const auto g = system.gate_charges(x);
    REQUIRE(q.size() == 2);
    REQUIRE(q[0] == g[0]);  // an electrode's charge as gate_charges
    REQUIRE(g[1] == 0.0);   // which leaves the ohmic contact's out
    REQUIRE(q[1] != 0.0);
    std::vector<double> f(system.unknowns());
    system.residual(x, f);
    std::vector<char> dirichlet(system.node_count(), 0);
    for (const device::Contact& c : d.contacts()) {
        for (const mesh::NodeId v : c.nodes) dirichlet[static_cast<std::size_t>(v)] = 1;
    }
    const double volume_scale = std::pow(scaling.L_D, 2);
    double rows = 0.0, inside = 0.0, magnitude = 0.0;
    for (std::size_t i = 0; i < system.node_count(); ++i) {
        if (dirichlet[i] != 0) continue;
        rows += f[3 * i];
        if (d.is_insulator(static_cast<mesh::NodeId>(i))) continue;
        const double C = (d.donors()[i] - d.acceptors()[i]) / scaling.Ns;
        const double rho = d.mesh().volumes()[i] / volume_scale * (x[3 * i + 2] - x[3 * i + 1] + C);
        inside += rho;
        magnitude += std::abs(rho);
    }
    double area = 0.0;  // of the oxide-silicon edges
    for (const mesh::Edge& e : d.mesh().edges()) {
        if (d.is_insulator(e.first) != d.is_insulator(e.second)) area += e.coupling_area;
    }
    const double interface =
        system.interface_trap_charges(x)[0] + 3e11 * area / (scaling.Ns * volume_scale);
    const double total = q[0] + q[1] + inside + interface;
    CAPTURE(rows, total, q[0], q[1], inside, interface);
    REQUIRE(std::abs(rows - total) <= 1e-12 * (magnitude + std::abs(q[0]) + std::abs(q[1])));
}

TEST_CASE("transient assemble: occupancies outside [0, 1] keep the interface solve on its root") {
    // A BDF2 history may leave [0, 1]. With every slot's history at 1.2 (or -0.2) and a 1 fs
    // step the occupancies stay near it, and the interface charge lies outside the bounds of
    // occupancies in [0, 1]; the local solve's bracket allows for that, so psi_I is still the root
    // of Gauss's law: the Poisson rows sum to the contact charges plus the charge inside.
    const device::Device d = mos(1);
    const auto scaling = *assemble::make_scaling(d);
    auto system = *DriftDiffusion::create(d, scaling);
    REQUIRE(system.set_bias(std::vector<double>{-0.4, 0.1}).has_value());
    const auto poisson = *EquilibriumPoisson::create(d, scaling);
    const auto x = probe_state(system, poisson, d, 71);
    const History h = history(system, x, 73);
    std::vector<char> dirichlet(system.node_count(), 0);
    for (const device::Contact& c : d.contacts()) {
        for (const mesh::NodeId v : c.nodes) dirichlet[static_cast<std::size_t>(v)] = 1;
    }
    double area = 0.0;
    for (const mesh::Edge& e : d.mesh().edges()) {
        if (d.is_insulator(e.first) != d.is_insulator(e.second)) area += e.coupling_area;
    }
    for (const double c : {1.2, -0.2}) {
        const std::vector<double> traps(system.trap_slots(), c);
        const DriftDiffusion::TimeStep step{system.time_scale() / 1e-15, h.storage, traps};
        for (const double f : system.trap_occupancies(x, &step)) {
            REQUIRE(std::abs(f - c) < 1e-3);  // near the history
        }
        std::vector<double> r(system.unknowns());
        system.residual(x, step, r);
        const auto q = system.contact_charges(x, &step);
        double rows = 0.0, inside = 0.0, magnitude = 0.0;
        for (std::size_t i = 0; i < system.node_count(); ++i) {
            if (dirichlet[i] != 0) continue;
            rows += r[3 * i];
            if (d.is_insulator(static_cast<mesh::NodeId>(i))) continue;
            const double C = (d.donors()[i] - d.acceptors()[i]) / scaling.Ns;
            const double rho = d.mesh().volumes()[i] / scaling.L_D *
                               (x[3 * i + 2] - x[3 * i + 1] + C);
            inside += rho;
            magnitude += std::abs(rho);
        }
        const double interface =
            system.interface_trap_charges(x, &step)[0] + 3e11 * area / (scaling.Ns * scaling.L_D);
        const double total = q[0] + q[1] + inside + interface;
        CAPTURE(c, rows, total, interface);
        REQUIRE(std::abs(rows - total) <= 1e-12 * (magnitude + std::abs(q[0]) + std::abs(q[1])));
    }
}
