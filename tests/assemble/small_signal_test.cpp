// Small-signal parts of the drift-diffusion assembler (ARCHITECTURE.md section 11, Unit 22): the
// matrix J + s C + T(s) is the steady Jacobian at s = 0 and a backward-Euler step's Jacobian at a
// real s (the step from the steady trap occupancy), it is holomorphic in s, and the bias
// derivative and the current rows agree with finite differences of the residual, the conduction
// currents and the contact charges, in 1D, 2D and 3D, with interface traps, a lumped gate, an
// electrode and an ohmic contact on an interface edge.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <complex>
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
using Complex = std::complex<double>;

namespace {

std::vector<double> uniform(double a, double b, int nodes) {
    std::vector<double> v(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) v[static_cast<std::size_t>(i)] = a + (b - a) * i / (nodes - 1);
    return v;
}

// The Unit 21 assembler fixture: oxide on [-5 nm, 0) with the electrode on x_min, p-silicon
// (1e17) on (0, 120 nm] with an n+ pocket (1e19, 2D and 3D) and the substrate contact on x_max;
// a fixed charge, two trap levels, a trap band and surface recombination at the interface. With
// `source` (2D) an ohmic contact on the pocket's y_max face, whose first node is the
// semiconductor end of an interface edge.
device::Device mos(int dimension, bool traps = true, bool source = false) {
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
    std::vector<mesh::NodeId> electrode, drain;
    for (const mesh::NodeId v : m.find_boundary("x_min")->nodes) {
        if (dimension == 1 || m.points()[static_cast<std::size_t>(v)][1] < 1.5e-6) {
            electrode.push_back(v);
        }
    }
    if (source) {
        for (const mesh::NodeId v : m.find_boundary("y_max")->nodes) {
            const double px = m.points()[static_cast<std::size_t>(v)][0];
            if (px > 0.0 && px < 3e-6) drain.push_back(v);
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
    std::vector<device::Contact> contacts{
        std::move(gate), {"substrate", device::ContactKind::ohmic, std::move(substrate)}};
    if (source) contacts.push_back({"source", device::ContactKind::ohmic, std::move(drain)});
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
         .node_region = std::move(region),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = std::move(contacts),
         .interfaces = {std::move(f)}});
}


// A one-cell oxide: the electrode node is the insulator end of the interface edge (Unit 15b's
// case), so its Dirichlet row must not take the interface terms. 1D, the mos() interface.
device::Device mos_one_cell() {
    std::vector<double> x{-5e-7};
    for (const double v : uniform(5e-7, 1.2e-5, 30)) x.push_back(v);
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n, 1);
    std::vector<double> acceptors(n, 1e17);
    region[0] = 0;
    acceptors[0] = 0.0;
    device::Interface f{"oxide", "silicon"};
    f.fixed_charge_cm2 = 3e11;
    f.traps.levels = {{.type = physics::TrapType::acceptor, .density_cm2 = 1e11, .energy_eV = 0.25,
                       .sigma_n_cm2 = 2e-16, .sigma_p_cm2 = 1e-15}};
    f.recombination_velocity_n_cm_s = 1e3;
    f.recombination_velocity_p_cm_s = 4e2;
    device::Contact gate{"gate", device::ContactKind::electrode, m.find_boundary("x_min")->nodes};
    gate.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = 4.4};
    auto substrate = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
         .node_region = std::move(region),
         .donors = std::vector<double>(n, 0.0),
         .acceptors = std::move(acceptors),
         .contacts = {std::move(gate),
                      {"substrate", device::ContactKind::ohmic, std::move(substrate)}},
         .interfaces = {std::move(f)}});
}

// A lumped-gate MOS capacitor (Unit 12) on 1e17 p-silicon, 1D or 2D (a gate on part of x_min).
device::Device lumped(int dimension) {
    const auto x = uniform(0.0, 1.2e-5, 40);
    mesh::Mesh m = dimension == 1 ? *mesh::make_tensor_grid(x)
                                  : *mesh::make_tensor_grid(x, uniform(0.0, 3e-6, 4));
    const std::size_t n = m.node_count();
    std::vector<mesh::NodeId> g;
    for (const mesh::NodeId v : m.find_boundary("x_min")->nodes) {
        if (dimension == 1 || m.points()[static_cast<std::size_t>(v)][1] < 2e-6) g.push_back(v);
    }
    device::Contact gate{"gate", device::ContactKind::gate, std::move(g)};
    gate.gate = {.boundary = "x_min", .oxide_thickness_cm = 5e-7, .fixed_charge_cm2 = 2e11};
    auto substrate = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"si", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::vector<double>(n, 0.0),
         .acceptors = std::vector<double>(n, 1e17),
         .contacts = {std::move(gate),
                      {"substrate", device::ContactKind::ohmic, std::move(substrate)}}});
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

template <class Scalar>
std::vector<std::vector<Scalar>> dense(const linalg::BasicSparseMatrix<Scalar>& a) {
    const auto n = static_cast<std::size_t>(a.rows());
    std::vector<std::vector<Scalar>> d(n, std::vector<Scalar>(n, Scalar{}));
    for (std::size_t r = 0; r < n; ++r) {
        for (auto k = static_cast<std::size_t>(a.row_offsets()[r]);
             k < static_cast<std::size_t>(a.row_offsets()[r + 1]); ++k) {
            d[r][static_cast<std::size_t>(a.col_indices()[k])] = a.values()[k];
        }
    }
    return d;
}

// The Unit 15b probe state: the charge-neutral guess with noise, and at the interface nodes
// n = 2e-3 and p = 5e-4 (scaled), so the traps are out of equilibrium.
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

// A case: the device, its models and bias, the probe state.
struct Case {
    const char* name;
    device::Device device;
    assemble::PhysicsModels models;
    std::vector<double> bias;
};

std::vector<Case> cases() {
    std::vector<Case> c;
    for (const int D : {1, 2, 3}) {
        c.push_back({"mos", mos(D), {}, {-0.4, 0.1}});
    }
    c.push_back({"mos fermi-dirac ionization", mos(2), {.fermi_dirac = true,
                                                       .incomplete_ionization = true},
                 {-0.4, 0.1}});
    c.push_back({"mos source", mos(2, true, true), {.field_mobility = true}, {-0.4, 0.1, 0.2}});
    c.push_back({"mos one-cell oxide", mos_one_cell(), {}, {-0.4, 0.1}});
    c.push_back({"lumped 1D", lumped(1), {}, {0.3, 0.0}});
    c.push_back({"lumped 2D", lumped(2), {.fermi_dirac = true}, {-0.2, 0.05}});
    return c;
}

// The step of rate r from the steady storage and trap occupancies at x.
struct Step {
    std::vector<double> storage, traps;
    DriftDiffusion::TimeStep step;
};

Step steady_step(const DriftDiffusion& dd, std::span<const double> x, double r) {
    Step s{dd.storage(x), dd.trap_occupancies(x), {}};
    s.step = {r, s.storage, s.traps};
    return s;
}

// The largest entry of each column (for column-relative comparisons).
template <class Scalar>
std::vector<double> column_scale(const std::vector<std::vector<Scalar>>& a) {
    std::vector<double> s(a.size(), 1e-300);
    for (const auto& row : a) {
        for (std::size_t c = 0; c < row.size(); ++c) s[c] = std::max(s[c], std::abs(row[c]));
    }
    return s;
}

}  // namespace

TEST_CASE("small-signal assemble: at s = 0 the matrix is the steady Jacobian") {
    for (Case& k : cases()) {
        const auto scaling = *assemble::make_scaling(k.device);
        auto dd = *DriftDiffusion::create(k.device, scaling, k.models);
        REQUIRE(dd.set_bias(k.bias).has_value());
        const auto poisson = *EquilibriumPoisson::create(k.device, scaling, k.models);
        const auto x = probe_state(dd, poisson, k.device, 7);
        auto j = dd.make_jacobian();
        std::vector<double> f(dd.unknowns());
        dd.evaluate(x, f, j);
        auto a = dd.make_small_signal_matrix();
        dd.small_signal_matrix(x, 0.0, a);
        const auto dj = dense(j);
        const auto da = dense(a);
        const auto scale = column_scale(dj);
        double worst = 0.0, imag = 0.0;
        for (std::size_t r = 0; r < dj.size(); ++r) {
            for (std::size_t c = 0; c < dj.size(); ++c) {
                imag = std::max(imag, std::abs(da[r][c].imag()));
                worst = std::max(worst, std::abs(da[r][c].real() - dj[r][c]) / scale[c]);
            }
        }
        CAPTURE(k.name, worst);
        REQUIRE(imag == 0.0);
        REQUIRE(worst <= 1e-14);  // measured at most 9.8e-16
    }
}

TEST_CASE("small-signal assemble: at a real s the matrix is a backward-Euler step's Jacobian") {
    for (Case& k : cases()) {
        const auto scaling = *assemble::make_scaling(k.device);
        auto dd = *DriftDiffusion::create(k.device, scaling, k.models);
        REQUIRE(dd.set_bias(k.bias).has_value());
        const auto poisson = *EquilibriumPoisson::create(k.device, scaling, k.models);
        const auto x = probe_state(dd, poisson, k.device, 11);
        for (const double dt : {1e-12, 1e-9, 1e-6}) {
            const double r = dd.time_scale() / dt;
            const Step s = steady_step(dd, x, r);
            auto j = dd.make_jacobian();
            std::vector<double> f(dd.unknowns());
            dd.evaluate(x, s.step, f, j);
            auto a = dd.make_small_signal_matrix();
            dd.small_signal_matrix(x, r, a);
            const auto dj = dense(j);
            const auto da = dense(a);
            const auto scale = column_scale(dj);
            double worst = 0.0, imag = 0.0;
            for (std::size_t row = 0; row < dj.size(); ++row) {
                for (std::size_t c = 0; c < dj.size(); ++c) {
                    imag = std::max(imag, std::abs(da[row][c].imag()));
                    worst = std::max(worst, std::abs(da[row][c].real() - dj[row][c]) / scale[c]);
                }
            }
            CAPTURE(k.name, dt, worst);
            REQUIRE(imag == 0.0);
            REQUIRE(worst <= 1e-13);  // measured at most 1.6e-15
        }
    }
}

TEST_CASE("small-signal assemble: the trap part at a real s on its own") {
    // With and without the traps, so the trap entries (some 1e-5 of the flux entries) cannot hide
    // in the comparison above.
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
        const auto x = probe_state(a, poisson, with, 5);
        for (const double dt : {1e-12, 1e-9, 1e-6}) {
            const double r = a.time_scale() / dt;
            const Step sa = steady_step(a, x, r), sb = steady_step(b, x, r);
            const std::size_t n = a.unknowns();
            auto ja = a.make_jacobian(), jb = b.make_jacobian();
            std::vector<double> f(n);
            a.evaluate(x, sa.step, f, ja);
            b.evaluate(x, sb.step, f, jb);
            auto ca = a.make_small_signal_matrix(), cb = b.make_small_signal_matrix();
            a.small_signal_matrix(x, r, ca);
            b.small_signal_matrix(x, r, cb);
            const auto dja = dense(ja), djb = dense(jb);
            const auto dca = dense(ca), dcb = dense(cb);
            double worst = 0.0, largest = 0.0;
            for (std::size_t row = 0; row < n; ++row) {
                for (std::size_t c = 0; c < n; ++c) {
                    const double step = dja[row][c] - djb[row][c];
                    const Complex small = dca[row][c] - dcb[row][c];
                    largest = std::max(largest, std::abs(step));
                    worst = std::max(worst, std::abs(small - step));
                }
            }
            CAPTURE(fd, dt, worst, largest);
            REQUIRE(largest > 0.0);
            REQUIRE(worst <= 1e-10 * largest);  // measured at most 2.1e-13
        }
    }
}

TEST_CASE("small-signal assemble: the matrix and the current rows are holomorphic in s") {
    // Cauchy-Riemann: the derivative along the real axis equals that along the imaginary axis.
    // A conjugate in place of a value would break it.
    for (Case& k : cases()) {
        const auto scaling = *assemble::make_scaling(k.device);
        auto dd = *DriftDiffusion::create(k.device, scaling, k.models);
        REQUIRE(dd.set_bias(k.bias).has_value());
        const auto poisson = *EquilibriumPoisson::create(k.device, scaling, k.models);
        const auto x = probe_state(dd, poisson, k.device, 13);
        const double r = dd.time_scale() / 1e-9;
        const Complex s0{0.3 * r, r}, h = 1e-4 * r;
        const Complex i{0.0, 1.0};
        const auto matrix = [&](Complex s) {
            auto a = dd.make_small_signal_matrix();
            dd.small_signal_matrix(x, s, a);
            return std::vector<Complex>(a.values().begin(), a.values().end());
        };
        const auto rp = matrix(s0 + h), rm = matrix(s0 - h);
        const auto ip = matrix(s0 + i * h), im = matrix(s0 - i * h);
        double worst = 0.0, largest = 0.0;
        for (std::size_t q = 0; q < rp.size(); ++q) {
            const Complex along_real = (rp[q] - rm[q]) / (2.0 * h);
            const Complex along_imag = (ip[q] - im[q]) / (2.0 * i * h);
            largest = std::max(largest, std::abs(along_real));
            worst = std::max(worst, std::abs(along_real - along_imag));
        }
        CAPTURE(k.name, worst, largest);
        REQUIRE(largest > 0.0);
        REQUIRE(worst <= 1e-6 * largest);  // measured at most 1.7e-8

        const auto rows = [&](Complex s) {
            std::vector<Complex> v;
            for (const auto& row : dd.small_signal_currents(x, s)) {
                v.insert(v.end(), row.values.begin(), row.values.end());
                v.push_back(row.bias);
            }
            return v;
        };
        const auto qp = rows(s0 + h), qm = rows(s0 - h);
        const auto jp = rows(s0 + i * h), jm = rows(s0 - i * h);
        REQUIRE(qp.size() == jp.size());
        worst = largest = 0.0;
        for (std::size_t q = 0; q < qp.size(); ++q) {
            const Complex along_real = (qp[q] - qm[q]) / (2.0 * h);
            const Complex along_imag = (jp[q] - jm[q]) / (2.0 * i * h);
            largest = std::max(largest, std::abs(along_real));
            worst = std::max(worst, std::abs(along_real - along_imag));
        }
        CAPTURE(worst, largest);
        REQUIRE(largest > 0.0);
        REQUIRE(worst <= 1e-6 * largest);  // the rows likewise
    }
}

TEST_CASE("small-signal assemble: the bias derivative is the residual's") {
    for (Case& k : cases()) {
        const auto scaling = *assemble::make_scaling(k.device);
        auto dd = *DriftDiffusion::create(k.device, scaling, k.models);
        REQUIRE(dd.set_bias(k.bias).has_value());
        const auto poisson = *EquilibriumPoisson::create(k.device, scaling, k.models);
        const auto x = probe_state(dd, poisson, k.device, 17);
        const std::size_t n = dd.unknowns();
        for (std::size_t c = 0; c < k.bias.size(); ++c) {
            const double h = 1e-3;
            std::vector<double> up_bias = k.bias, down_bias = k.bias, up(n), down(n);
            up_bias[c] += h;
            down_bias[c] -= h;
            REQUIRE(dd.set_bias(up_bias).has_value());
            dd.residual(x, up);
            REQUIRE(dd.set_bias(down_bias).has_value());
            dd.residual(x, down);
            REQUIRE(dd.set_bias(k.bias).has_value());
            const std::vector<double> d = dd.bias_derivative(c);
            double worst = 0.0, largest = 0.0;
            for (std::size_t q = 0; q < n; ++q) {
                largest = std::max(largest, std::abs(d[q]));
                worst = std::max(worst, std::abs((up[q] - down[q]) / (2.0 * h) - d[q]));
            }
            CAPTURE(k.name, c, worst, largest);
            REQUIRE(largest > 0.0);
            REQUIRE(worst <= 1e-9 * largest);  // linear in the bias; measured 3.2e-14
        }
    }
}

TEST_CASE("small-signal assemble: the current rows are the linearized total currents") {
    // At a real s = r the total current of a backward-Euler step from x (steady occupancies) is
    // conduction_currents + r contact_charges (less the history's charge): its central difference
    // in each unknown and in each bias is the row. At s = 0 it is terminal_currents.
    for (Case& k : cases()) {
        const auto scaling = *assemble::make_scaling(k.device);
        auto dd = *DriftDiffusion::create(k.device, scaling, k.models);
        REQUIRE(dd.set_bias(k.bias).has_value());
        const auto poisson = *EquilibriumPoisson::create(k.device, scaling, k.models);
        auto x = probe_state(dd, poisson, k.device, 19);
        const std::size_t contacts = k.bias.size();
        for (const double dt : {0.0, 1e-12, 1e-9}) {
            const double r = dt == 0.0 ? 0.0 : dd.time_scale() / dt;
            const Step s = steady_step(dd, x, std::max(r, 1.0));
            const auto total = [&](std::span<const double> y) {
                if (r == 0.0) return dd.terminal_currents(y);
                std::vector<double> t = dd.conduction_currents(y, s.step);
                const std::vector<double> q = dd.contact_charges(y, &s.step);
                for (std::size_t c = 0; c < t.size(); ++c) t[c] += r * q[c];
                return t;
            };
            const auto rows = dd.small_signal_currents(x, r);
            REQUIRE(rows.size() == contacts);
            std::vector<std::vector<double>> analytic(contacts,
                                                      std::vector<double>(dd.unknowns(), 0.0));
            for (std::size_t c = 0; c < contacts; ++c) {
                for (std::size_t q = 0; q < rows[c].columns.size(); ++q) {
                    REQUIRE(rows[c].values[q].imag() == 0.0);
                    analytic[c][rows[c].columns[q]] = rows[c].values[q].real();
                }
            }
            std::vector<double> worst(contacts, 0.0), largest(contacts, 0.0);
            for (std::size_t col = 0; col < dd.unknowns(); ++col) {
                const double base = x[col];
                const double h = 1e-7 * std::max(std::abs(base), 1.0);
                x[col] = base + h;
                const double hi = x[col];
                const auto up = total(x);
                x[col] = base - h;
                const double lo = x[col];
                const auto down = total(x);
                x[col] = base;
                for (std::size_t c = 0; c < contacts; ++c) {
                    const double numeric = (up[c] - down[c]) / (hi - lo);
                    largest[c] = std::max(largest[c], std::abs(analytic[c][col]));
                    worst[c] = std::max(worst[c], std::abs(numeric - analytic[c][col]));
                }
            }
            for (std::size_t c = 0; c < contacts; ++c) {
                CAPTURE(k.name, dt, c, worst[c], largest[c]);
                // At s = 0 a gate or electrode carries no current: its row is zero.
                REQUIRE((largest[c] > 0.0 || r == 0.0));
                REQUIRE(worst[c] <= 1e-6 * largest[c]);  // measured at most 2.6e-8
            }
            // The bias term: the charges at fixed x through the gate potential.
            for (std::size_t c = 0; c < contacts; ++c) {
                const double h = 1e-3;
                std::vector<double> b = k.bias;
                b[c] += h;
                REQUIRE(dd.set_bias(b).has_value());
                const auto up = total(x);
                b[c] -= 2.0 * h;
                REQUIRE(dd.set_bias(b).has_value());
                const auto down = total(x);
                REQUIRE(dd.set_bias(k.bias).has_value());
                for (std::size_t e = 0; e < contacts; ++e) {
                    const double numeric = (up[e] - down[e]) / (2.0 * h);
                    const double expected = e == c ? rows[c].bias.real() : 0.0;
                    CAPTURE(k.name, dt, c, e, numeric, expected);
                    REQUIRE(std::abs(numeric - expected) <= 1e-9 * std::max(1.0, largest[e]));
                }
            }
        }
    }
}
