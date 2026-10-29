// Gate tunnelling in the drift-diffusion assembler (Unit 51; DECISIONS.md GD1-GD9, GC4, GC6, GC7,
// GC13): the finite-difference Jacobian gate on a 1D lumped gate and on 2D electrode stacks of one
// and two layers under both statistics, each process alone and all three; zero current at thermal
// equilibrium; the bias derivative and the small-signal current row against differences; the
// polysilicon energy limits; and the refusals of the device and the assembler.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <expected>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/equilibrium_poisson.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/generate.hpp"
#include "NiTCAD/mesh/tensor_cells.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD;
using assemble::DriftDiffusion;
using assemble::EquilibriumPoisson;

namespace {

constexpr double nm = 1e-7;
constexpr std::uint8_t all = device::all_tunnel_processes;

std::vector<double> uniform(double a, double b, int nodes) {
    std::vector<double> v(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) v[static_cast<std::size_t>(i)] = a + (b - a) * i / (nodes - 1);
    return v;
}

struct Mos {
    std::uint8_t processes = all;
    device::GateElectrode electrode = device::GateElectrode::n_poly;
    bool acceptors = true;   // p-substrate (else n)
    bool cells = true;       // the description's tensor cells (electrode stacks)
    int layers = 1;          // electrode stacks: 1 (SiO2) or 2 (SiO2 under a high-k layer)
    bool high_k_band = true;  // the high-k layer carries band data
    double eps_ox = 3.9;     // lumped gate
    bool gate_band = false;  // lumped gate: its own band data
};

// 1D: a 60 nm silicon bar, 5e17 doping, the substrate contact on x_max and a lumped gate (1.5 nm)
// on x_min.
device::Device lumped(const Mos& o) {
    mesh::Mesh m = *mesh::make_tensor_grid(uniform(0.0, 60 * nm, 31));
    const std::size_t n = m.node_count();
    std::vector<mesh::NodeId> right = m.find_boundary("x_max")->nodes;
    device::Contact gate{"gate", device::ContactKind::gate, m.find_boundary("x_min")->nodes,
                         {.boundary = "x_min", .oxide_thickness_cm = 1.5 * nm,
                          .oxide_relative_permittivity = o.eps_ox, .electrode = o.electrode,
                          .work_function_eV = 4.6}};
    if (o.gate_band) {
        gate.gate.oxide_affinity_eV = 2.4;
        gate.gate.oxide_band_gap_eV = 5.7;
        gate.gate.oxide_mass_n = 0.2;
        gate.gate.oxide_mass_p = 0.3;
    }
    gate.tunnelling = o.processes;
    auto d = device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::vector<double>(n, o.acceptors ? 0.0 : 5e17),
         .acceptors = std::vector<double>(n, o.acceptors ? 5e17 : 0.0),
         .contacts = {std::move(gate),
                      {"substrate", device::ContactKind::ohmic,
                       std::move(right)}}});
    if (!d) FAIL(d.error().message);
    return std::move(*d);
}

// A high-k insulator (Gehring Table 5.1's HfO2: eps 25, gap 5.7 eV, offset 1.5 eV to silicon) with
// test masses.
physics::Insulator high_k(bool band) {
    physics::InsulatorParameters p{.eps_r = 25.0};
    if (band) {
        p.electron_affinity_eV = 4.05 - 1.5;
        p.band_gap_eV = 5.7;
        p.tunnel_mass_n = 0.2;
        p.tunnel_mass_p = 0.3;
    }
    return *physics::Insulator::create(p);
}

// 2D: silicon y in [0, 20 nm], 5e17 doping, the substrate on y_min; SiO2 to 21.5 nm (one layer) or
// to 20.5 nm under the high-k layer to 22 nm (two); an electrode on y_max; x in [0, 20 nm].
device::Device stack(const Mos& o) {
    const std::vector<double> x = uniform(0.0, 20 * nm, 5);
    std::vector<double> y = uniform(0.0, 20 * nm, 9);
    const double top = o.layers == 1 ? 21.5 * nm : 22.0 * nm;
    for (double v = 20.5 * nm; v <= top + 1e-12; v += 0.5 * nm) y.push_back(v);
    mesh::Mesh m = *mesh::make_tensor_grid(x, y);
    auto cells = mesh::TensorCells::create(m, x, y);
    REQUIRE(cells.has_value());
    const std::size_t n = m.node_count();
    std::vector<device::Region> regions{{"silicon", physics::silicon()},
                                        {"sio2", physics::silicon_dioxide()}};
    if (o.layers == 2) regions.push_back({"high_k", high_k(o.high_k_band)});
    std::vector<device::RegionId> region(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        const double yi = m.points()[i][1];
        if (yi > 20 * nm + 1e-12) region[i] = (o.layers == 2 && yi > 20.5 * nm + 1e-12) ? 2 : 1;
    }
    std::vector<mesh::NodeId> bottom = m.find_boundary("y_min")->nodes;
    device::Contact electrode{"gate", device::ContactKind::electrode,
                              m.find_boundary("y_max")->nodes, {},
                              {.kind = o.electrode, .work_function_eV = 4.6}};
    electrode.tunnelling = o.processes;
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (region[i] == 0) (o.acceptors ? acceptors : donors)[i] = 5e17;
    }
    device::DeviceDescription desc{
        .mesh = std::move(m),
        .temperature_K = 300.0,
        .regions = std::move(regions),
        .node_region = std::move(region),
        .donors = std::move(donors),
        .acceptors = std::move(acceptors),
        .contacts = {{"substrate", device::ContactKind::ohmic, std::move(bottom)},
                     std::move(electrode)}};
    if (o.cells) desc.cells = std::move(*cells);
    auto d = device::Device::create(std::move(desc));
    if (!d) FAIL(d.error().message);
    return std::move(*d);
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

// The probe state: the charge-neutral guess at the given biases, the contacts stamped, then psi +
// 0.02 noise and n, p times (1 + 0.01 noise) off the contacts' Dirichlet nodes.
std::vector<double> probe(DriftDiffusion& system, const device::Device& d,
                          const assemble::Scaling& s, const assemble::PhysicsModels& models,
                          const std::vector<double>& bias, std::uint64_t seed, bool noise = true) {
    REQUIRE(system.set_bias(bias).has_value());
    const auto poisson = *EquilibriumPoisson::create(d, s, models);
    auto x = system.state_from_potential(poisson.charge_neutral_potential());
    system.stamp_contacts(x);
    if (!noise) return x;
    Noise r{seed};
    for (std::size_t i = 0; i < system.node_count(); ++i) {
        if (d.is_insulator(static_cast<mesh::NodeId>(i))) continue;
        x[3 * i] += 0.02 * r.next();
        x[3 * i + 1] *= 1.0 + 0.01 * r.next();
        x[3 * i + 2] *= 1.0 + 0.01 * r.next();
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

// As schottky_test.cpp: every column, central differences at steps of 1e-5, 1e-6 and 1e-7
// max(|x|, 1), the worst error over the column's largest entry.
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

// The tunnelling rows' part of the residual: with tunnelling minus without, at x.
std::vector<double> tunnel_rows(const device::Device& on, const device::Device& off,
                                const std::vector<double>& bias, std::vector<double> x,
                                const assemble::PhysicsModels& models) {
    const auto s = *assemble::make_scaling(on);
    auto a = *DriftDiffusion::create(on, s, models);
    auto b = *DriftDiffusion::create(off, s, models);
    REQUIRE(a.set_bias(bias).has_value());
    REQUIRE(b.set_bias(bias).has_value());
    std::vector<double> fa(a.unknowns()), fb(b.unknowns());
    a.residual(x, fa);
    b.residual(x, fb);
    for (std::size_t k = 0; k < fa.size(); ++k) fa[k] -= fb[k];
    return fa;
}

void require_invalid(const auto& made, const std::string& name) {
    CAPTURE(name);
    REQUIRE_FALSE(made.has_value());
    CAPTURE(made.error().message);
    REQUIRE(made.error().code == base::ErrorCode::invalid_input);
}

}  // namespace

TEST_CASE("gate tunnelling assemble: GC7 the finite-difference Jacobian, lumped gate") {
    for (const bool fermi_dirac : {false, true}) {
        for (const std::uint8_t processes : {std::uint8_t{1}, std::uint8_t{2}, std::uint8_t{4}, all}) {
            for (const bool p_type : {true, false}) {
                Mos o;
                o.processes = processes;
                o.acceptors = p_type;
                o.electrode = p_type ? device::GateElectrode::n_poly : device::GateElectrode::p_poly;
                const device::Device d = lumped(o);
                const auto s = *assemble::make_scaling(d);
                assemble::PhysicsModels models;
                models.fermi_dirac = fermi_dirac;
                auto made = DriftDiffusion::create(d, s, models);
                if (!made) FAIL(made.error().message);
                // Inversion: +2 V on an nMOS, -2 V on a pMOS (EVB needs the bands crossed).
                const double vg = p_type ? 2.0 : -2.0;
                const auto x = probe(*made, d, s, models, {vg, 0.0}, 11);
                const double error = fd_jacobian_error(*made, x);
                CAPTURE(fermi_dirac, int(processes), p_type, error);
                REQUIRE(error <= 1e-6);
            }
        }
    }
}

TEST_CASE("gate tunnelling assemble: GC7 the finite-difference Jacobian, electrode stacks") {
    for (const int layers : {1, 2}) {
        for (const bool fermi_dirac : {false, true}) {
            for (const auto electrode :
                 {device::GateElectrode::n_poly, device::GateElectrode::metal}) {
                Mos o;
                o.layers = layers;
                o.electrode = electrode;
                const device::Device d = stack(o);
                const auto s = *assemble::make_scaling(d);
                assemble::PhysicsModels models;
                models.fermi_dirac = fermi_dirac;
                auto made = DriftDiffusion::create(d, s, models);
                if (!made) FAIL(made.error().message);
                const auto x = probe(*made, d, s, models, {0.0, 2.0}, 5);
                const double error = fd_jacobian_error(*made, x);
                CAPTURE(layers, fermi_dirac, int(electrode), error);
                REQUIRE(error <= 1e-6);
            }
        }
    }
}

TEST_CASE("gate tunnelling assemble: the rows carry current and none at equilibrium (GC4)") {
    const assemble::PhysicsModels models;
    for (const int layers : {0, 1, 2}) {  // 0: lumped
        Mos on, off;
        on.layers = off.layers = std::max(layers, 1);
        off.processes = 0;
        const device::Device d_on = layers == 0 ? lumped(on) : stack(on);
        const device::Device d_off = layers == 0 ? lumped(off) : stack(off);
        const auto s = *assemble::make_scaling(d_on);
        auto system = *DriftDiffusion::create(d_on, s, models);
        const std::vector<double> zero{0.0, 0.0};
        const auto eq = probe(system, d_on, s, models, zero, 0, false);
        const auto rows = tunnel_rows(d_on, d_off, zero, eq, models);
        CAPTURE(layers);
        REQUIRE(std::ranges::all_of(rows, [](double v) { return v == 0.0; }));
        // Biased: the gate's terminal current is the tunnels' sum and the rows are nonzero.
        const std::vector<double> bias = layers == 0 ? std::vector<double>{2.0, 0.0}
                                                     : std::vector<double>{0.0, 2.0};
        const auto x = probe(system, d_on, s, models, bias, 0, false);
        const auto biased = tunnel_rows(d_on, d_off, bias, x, models);
        REQUIRE(std::ranges::any_of(biased, [](double v) { return v != 0.0; }));
        const std::size_t gate = layers == 0 ? 0 : 1;
        const auto I = system.terminal_currents(x);
        double sum = 0.0;  // T_n - T_p = -(electron rows) - (hole rows)
        for (std::size_t i = 0; i < system.node_count(); ++i) {
            sum -= biased[3 * i + 1] + biased[3 * i + 2];
        }
        auto off_system = *DriftDiffusion::create(d_off, s, models);
        REQUIRE(off_system.set_bias(bias).has_value());
        const auto I_off = off_system.terminal_currents(x);
        REQUIRE(std::abs((I[gate] - I_off[gate]) - sum) <= 1e-12 * std::abs(sum));
        REQUIRE(sum != 0.0);
    }
}

TEST_CASE("gate tunnelling assemble: the bias derivative and the small-signal row") {
    const assemble::PhysicsModels models;
    const device::Device d = lumped({});
    const auto s = *assemble::make_scaling(d);
    auto system = *DriftDiffusion::create(d, s, models);
    REQUIRE(system.bias_nonlinear());
    const auto x = probe(system, d, s, models, {1.5, 0.0}, 3);
    // dF/dV_G against central differences of the residual in the gate's bias.
    const std::vector<double> dF = system.bias_derivative(0, x);
    const double h = 1e-5;
    std::vector<double> fp(system.unknowns()), fm(system.unknowns());
    REQUIRE(system.set_bias(std::vector<double>{1.5 + h, 0.0}).has_value());
    system.residual(x, fp);
    REQUIRE(system.set_bias(std::vector<double>{1.5 - h, 0.0}).has_value());
    system.residual(x, fm);
    REQUIRE(system.set_bias(std::vector<double>{1.5, 0.0}).has_value());
    double scale = 0.0, error = 0.0;
    for (std::size_t r = 0; r < dF.size(); ++r) {
        scale = std::max(scale, std::abs(dF[r]));
        error = std::max(error, std::abs((fp[r] - fm[r]) / (2 * h) - dF[r]));
    }
    REQUIRE(error <= 1e-6 * scale);
    const auto rows = system.bias_rows(0);
    for (std::size_t r = 0; r < dF.size(); ++r) {
        if (dF[r] != 0.0) REQUIRE(std::ranges::find(rows, r) != rows.end());
    }
    // The gate's small-signal row at s = 0 against differences of its terminal current.
    const auto row = system.small_signal_currents(x, {0.0, 0.0})[0];
    double largest = 0.0;  // the row's largest change per unit relative change of its column
    for (std::size_t k = 0; k < row.columns.size(); ++k) {
        const std::size_t c = row.columns[k];
        largest = std::max(largest, std::abs(row.values[k].real()) *
                                        (c % 3 == 0 ? 1.0 : std::abs(x[c])));
    }
    for (std::size_t k = 0; k < row.columns.size(); ++k) {
        std::vector<double> y = x;
        const std::size_t c = row.columns[k];
        const double step = 1e-6 * (c % 3 == 0 ? 1.0 : std::abs(x[c]));  // densities: relative
        y[c] = x[c] + step;
        const double up = system.terminal_currents(y)[0];
        y[c] = x[c] - step;
        const double down = system.terminal_currents(y)[0];
        const double fd = (up - down) / (2 * step);
        CAPTURE(c, fd, row.values[k].real());
        const double unit = c % 3 == 0 ? 1.0 : std::abs(x[c]);
        REQUIRE(std::abs(fd - row.values[k].real()) * unit <=
                1e-5 * std::max(std::abs(fd) * unit, 1e-6 * largest));
    }
    REQUIRE(system.set_bias(std::vector<double>{1.5 + h, 0.0}).has_value());
    const double up = system.terminal_currents(x)[0];
    REQUIRE(system.set_bias(std::vector<double>{1.5 - h, 0.0}).has_value());
    const double down = system.terminal_currents(x)[0];
    REQUIRE(std::abs((up - down) / (2 * h) - row.bias.real()) <= 1e-5 * std::abs(row.bias.real()));
}

TEST_CASE("gate tunnelling assemble: GC6 EVB needs the bands crossed") {
    // EVB alone on an n-poly nMOS: zero while the silicon valence band lies below the gate's
    // conduction band (accumulation and weak inversion), positive deep in inversion.
    const assemble::PhysicsModels models;
    Mos on, off;
    on.processes = static_cast<std::uint8_t>(device::TunnelProcess::evb);
    off.processes = 0;
    const device::Device d_on = lumped(on), d_off = lumped(off);
    const auto s = *assemble::make_scaling(d_on);
    auto system = *DriftDiffusion::create(d_on, s, models);
    const auto at = [&](double vg) {
        const std::vector<double> bias{vg, 0.0};
        const auto x = probe(system, d_on, s, models, bias, 0, false);
        return tunnel_rows(d_on, d_off, bias, x, models)[3 * 0 + 2];  // the node's hole row
    };
    REQUIRE(at(-1.0) == 0.0);
    REQUIRE(at(0.05) == 0.0);
    REQUIRE(at(3.0) != 0.0);
    // GD6, GD7: the rows carry -T_n and +T_p, and T_n - T_p is the reported EVB current (both of
    // its directions with the signs that make them add). Deep inversion: the semiconductor's
    // valence electrons into the gate (the hole row only); strong accumulation: the gate's valence
    // electrons into the semiconductor's conduction band (the electron row only).
    for (const double vg : {3.0, -3.0}) {
        const std::vector<double> bias{vg, 0.0};
        const auto x = probe(system, d_on, s, models, bias, 0, false);
        const std::vector<double> rows = tunnel_rows(d_on, d_off, bias, x, models);
        REQUIRE(system.set_bias(bias).has_value());
        const double evb = system.gate_tunnelling_currents(x).evb[0];
        CAPTURE(vg, rows[1], rows[2], evb);
        REQUIRE(rows[vg > 0.0 ? 2 : 1] != 0.0);
        REQUIRE(std::abs(-(rows[1] + rows[2]) - evb) <= 1e-12 * std::abs(evb));
    }
}

TEST_CASE("gate tunnelling assemble: GC13 refusals") {
    {  // tunnelling on an ohmic contact
        mesh::Mesh m = *mesh::make_tensor_grid(uniform(0.0, 60 * nm, 5));
        const std::size_t n = m.node_count();
        device::Contact a{"a", device::ContactKind::ohmic, m.find_boundary("x_min")->nodes};
        a.tunnelling = all;
        require_invalid(
            device::Device::create({.mesh = std::move(m),
                                    .temperature_K = 300.0,
                                    .regions = {{"silicon", physics::silicon()}},
                                    .node_region = std::vector<device::RegionId>(n, 0),
                                    .donors = std::vector<double>(n, 1e17),
                                    .acceptors = std::vector<double>(n, 0.0),
                                    .contacts = {std::move(a)}}),
            "ohmic tunnelling");
    }
    const auto lumped_made = [](const Mos& o) {
        mesh::Mesh m = *mesh::make_tensor_grid(uniform(0.0, 60 * nm, 5));
        const std::size_t n = m.node_count();
        std::vector<mesh::NodeId> right = m.find_boundary("x_max")->nodes;
        device::Contact gate{"gate", device::ContactKind::gate, m.find_boundary("x_min")->nodes,
                             {.boundary = "x_min", .oxide_thickness_cm = 1.5 * nm,
                              .oxide_relative_permittivity = o.eps_ox}};
        gate.tunnelling = o.processes;
        return device::Device::create(
            {.mesh = std::move(m),
             .temperature_K = 300.0,
             .regions = {{"silicon", physics::silicon()}},
             .node_region = std::vector<device::RegionId>(n, 0),
             .donors = std::vector<double>(n, 0.0),
             .acceptors = std::vector<double>(n, 5e17),
             .contacts = {std::move(gate),
                          {"substrate", device::ContactKind::ohmic,
                           std::move(right)}}});
    };
    Mos o;
    o.processes = 8;
    require_invalid(lumped_made(o), "unknown process");
    o = {};
    o.eps_ox = 25.0;
    require_invalid(lumped_made(o), "high-k lumped stack without band data");
    o.processes = 0;
    REQUIRE(lumped_made(o).has_value());  // not tunnelling: no band data needed
    {
        Mos ok;
        ok.eps_ox = 25.0;
        ok.gate_band = true;
        const device::Device d = lumped(ok);
        REQUIRE(DriftDiffusion::create(d, *assemble::make_scaling(d), {}).has_value());
    }
    const device::Device d = lumped({});
    const auto s = *assemble::make_scaling(d);
    assemble::PhysicsModels thermal;
    thermal.electrothermal = true;
    require_invalid(DriftDiffusion::create(d, s, thermal), "electrothermal");
    {  // quantum confinement (owner ruling 2026-10-29: not validated together)
        assemble::PhysicsModels mlda;
        mlda.quantum_confinement = assemble::QuantumConfinement::mlda;
        const auto made = DriftDiffusion::create(d, s, mlda);
        require_invalid(made, "quantum confinement");
        REQUIRE(made.error().message.find("gate tunnelling with quantum confinement") !=
                std::string::npos);
    }
    Mos no_cells;
    no_cells.cells = false;
    const device::Device dc = stack(no_cells);
    require_invalid(DriftDiffusion::create(dc, *assemble::make_scaling(dc), {}),
                    "electrode without tensor cells");
    Mos no_band;
    no_band.layers = 2;
    no_band.high_k_band = false;
    const device::Device db = stack(no_band);
    require_invalid(DriftDiffusion::create(db, *assemble::make_scaling(db), {}),
                    "insulator without band data");
}
