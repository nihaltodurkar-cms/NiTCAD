// Gate contact terms (ARCHITECTURE.md section 11, Unit 12): electrode work functions, the scaled
// oxide coupling against the legacy kappa, the 1D rows against the legacy MOS-C rows written out by
// hand, zero boundary flux for the carriers, the FD-Jacobian gate of section 10 with gates in both
// assemblers, and the equilibrium system's gate bias.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/assemble/equilibrium_poisson.hpp"
#include "NiTCAD/assemble/gate.hpp"
#include "NiTCAD/assemble/scaling.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

using namespace NiTCAD;
using assemble::DriftDiffusion;
using assemble::EquilibriumPoisson;

namespace {

constexpr double N_A = 1e17;
constexpr double t_ox = 5e-7;
constexpr double Q_f = 2e11;

// The legacy MOS-C grading, x = L expm1(6 s) / expm1(6), with fewer nodes.
std::vector<double> moscap_axis(std::size_t nodes) {
    std::vector<double> x(nodes);
    for (std::size_t i = 0; i < nodes; ++i) {
        const double s = static_cast<double>(i) / static_cast<double>(nodes - 1);
        x[i] = 2e-4 * (std::expm1(6.0 * s) / std::expm1(6.0));
    }
    return x;
}

std::vector<double> uniform_axis(double length, int nodes) {
    std::vector<double> a(static_cast<std::size_t>(nodes));
    for (int i = 0; i < nodes; ++i) a[static_cast<std::size_t>(i)] = length * i / (nodes - 1);
    return a;
}

device::GateStack stack(const char* boundary, double Qf = Q_f) {
    return {.boundary = boundary, .oxide_thickness_cm = t_ox, .fixed_charge_cm2 = Qf};
}

// p-type MOS-C: gate (contact 0) on x_min, substrate (contact 1) on x_max; without the gate if
// with_gate is false.
device::Device moscap(mesh::Mesh m, bool with_gate = true) {
    const std::size_t n = m.node_count();
    std::vector<device::Contact> contacts;  // built before the mesh is moved into the device
    if (with_gate) {
        device::Contact g{"gate", device::ContactKind::gate, m.find_boundary("x_min")->nodes};
        g.gate = stack("x_min");
        contacts.push_back(std::move(g));
    }
    contacts.push_back(
        {"substrate", device::ContactKind::ohmic, m.find_boundary("x_max")->nodes});
    return *device::Device::create({.mesh = std::move(m),
                                    .temperature_K = 300.0,
                                    .regions = {{"silicon", physics::silicon()}},
                                    .node_region = std::vector<device::RegionId>(n, 0),
                                    .donors = std::vector<double>(n, 0.0),
                                    .acceptors = std::vector<double>(n, N_A),
                                    .contacts = std::move(contacts)});
}

// A 2D (or 3D, with z) MOSFET-like device: p-type body with n+ source and drain wells under the
// ends of the y_min face, ohmic source and drain on y_min next to a metal gate in the middle (with
// fixed charge), body contact on y_max. Contacts: source, gate, drain, body.
device::Device transistor(const std::vector<double>& z = {}) {
    const auto x = uniform_axis(3e-5, 13);
    const auto y = moscap_axis(9);
    mesh::Mesh m = z.empty() ? *mesh::make_tensor_grid(x, y) : *mesh::make_tensor_grid(x, y, z);
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = m.points()[i];
        if ((p[0] < 0.6e-5 || p[0] > 2.4e-5) && p[1] < 1e-5) donors[i] = 1e19;
    }
    std::vector<mesh::NodeId> source, gate, drain;
    for (const mesh::NodeId v : m.find_boundary("y_min")->nodes) {
        const double px = m.points()[static_cast<std::size_t>(v)][0];
        (px < 0.6e-5 ? source : px > 2.4e-5 ? drain : gate).push_back(v);
    }
    device::Contact g{"gate", device::ContactKind::gate, std::move(gate)};
    g.gate = stack("y_min");
    g.gate.electrode = device::GateElectrode::metal;
    g.gate.work_function_eV = 4.3;
    auto body = m.find_boundary("y_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::move(donors),
         .acceptors = std::vector<double>(n, N_A),
         .contacts = {{"source", device::ContactKind::ohmic, std::move(source)},
                      std::move(g),
                      {"drain", device::ContactKind::ohmic, std::move(drain)},
                      {"body", device::ContactKind::ohmic, std::move(body)}}});
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

// Section 10's gate (legacy _jacobian_probe normalization): central differences, step
// 1e-7 max(|u|, 1), every column; per column max |fd - J| over the column's largest |J| (+1e-30).
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
        double scale = 1e-30, error = 0.0;
        for (std::size_t r = 0; r < n; ++r) scale = std::max(scale, std::abs(j[r][c]));
        for (std::size_t r = 0; r < n; ++r) {
            error = std::max(error, std::abs((fp[r] - fm[r]) / (up - down) - j[r][c]));
        }
        worst = std::max(worst, error / scale);
    }
    return worst;
}

// Drift-diffusion at the given biases, at the legacy probe's state: Boltzmann carriers at the
// charge-neutral potential, then psi + 0.02 noise and n, p times (1 + 0.01 noise).
std::pair<DriftDiffusion, std::vector<double>> probe_state(const device::Device& d,
                                                           const std::vector<double>& bias,
                                                           std::uint64_t seed) {
    const auto scaling = *assemble::make_scaling(d);
    auto system = *DriftDiffusion::create(d, scaling);
    REQUIRE(system.set_bias(bias).has_value());
    const auto poisson = *EquilibriumPoisson::create(d, scaling);
    auto x = system.state_from_potential(poisson.charge_neutral_potential());
    Noise noise{seed};
    for (std::size_t i = 0; i < system.node_count(); ++i) {
        x[3 * i] += 0.02 * noise.next();
        x[3 * i + 1] *= 1.0 + 0.01 * noise.next();
        x[3 * i + 2] *= 1.0 + 0.01 * noise.next();
    }
    return {std::move(system), std::move(x)};
}

}  // namespace

TEST_CASE("gate: electrode work functions and the midgap offset") {
    const physics::Semiconductor si = physics::silicon();
    const double Eg = physics::band_gap_eV(si, 300.0);
    device::GateStack g = stack("x_min");
    REQUIRE(assemble::gate_work_function_eV(g, si, 300.0) == 4.05);  // n+ poly: chi
    REQUIRE(assemble::gate_midgap_offset_V(g, si, 300.0) == 4.05 - 4.05 - 0.5 * Eg);
    g.electrode = device::GateElectrode::p_poly;  // chi + Eg
    REQUIRE(assemble::gate_work_function_eV(g, si, 300.0) == 4.05 + Eg);
    REQUIRE(std::abs(assemble::gate_midgap_offset_V(g, si, 300.0) - 0.5 * Eg) < 1e-15);
    g.electrode = device::GateElectrode::metal;  // legacy "Al": 4.10 eV
    g.work_function_eV = 4.10;
    REQUIRE(assemble::gate_work_function_eV(g, si, 300.0) == 4.10);
    REQUIRE(std::abs(assemble::gate_midgap_offset_V(g, si, 300.0) - (0.05 - 0.5 * Eg)) < 1e-15);
    // The band gap follows the temperature (Varshni).
    g.electrode = device::GateElectrode::p_poly;
    REQUIRE(assemble::gate_work_function_eV(g, si, 400.0) ==
            4.05 + physics::band_gap_eV(si, 400.0));
}

TEST_CASE("gate: the scaled coupling is the legacy kappa times the scaled face area") {
    const physics::Semiconductor si = physics::silicon();
    const auto d = moscap(*mesh::make_tensor_grid(moscap_axis(20)));
    const auto s = *assemble::make_scaling(d);
    const double eps_s = 11.7 * base::eps0_F_per_cm, eps_ox = 3.9 * base::eps0_F_per_cm;
    const double kappa = eps_ox * s.L_D / (eps_s * t_ox);  // moscap.MOSCapacitor.kappa
    const double Cox = eps_ox / t_ox;
    for (const int D : {1, 2, 3}) {
        CAPTURE(D);
        const double area = D == 1 ? 1.0 : D == 2 ? 3e-6 : 4e-11;  // cm^(D-1)
        const assemble::GateTerm t = assemble::gate_term(stack("x_min"), si, area, D, s);
        const double scaled_area = area / std::pow(s.L_D, D - 1);
        REQUIRE(std::abs(t.coupling - kappa * scaled_area) <= 1e-15 * kappa * scaled_area);
        REQUIRE(t.offset == assemble::gate_midgap_offset_V(stack("x_min"), si, 300.0) / s.V_T);
        // The fixed charge as the legacy has it: inside -V_FB, kappa q Q_f / (C_ox V_T) per area.
        const double legacy_charge = kappa * base::q_C * Q_f / (Cox * s.V_T) * scaled_area;
        REQUIRE(std::abs(t.sheet_charge - legacy_charge) <= 1e-14 * legacy_charge);
    }
}

TEST_CASE("gate: the 1D rows are the legacy MOS-C rows written out by hand") {
    // moscap.MOSCapacitor.solve_psi, with V_FB from moscap.flatband_voltage (n+ poly, Q_f):
    //   F_0 = (psi_1 - psi_0) / h_0 + kappa (Vg - V_FB - (psi_0 - psi_b) V_T) / V_T - dV_0 rho_0
    //   F_i = (psi_{i+1} - psi_i) / h_i - (psi_i - psi_{i-1}) / h_{i-1} - dV_i rho_i
    // with rho = n - p - C. The equilibrium system uses Boltzmann carriers; the drift-diffusion
    // Poisson rows take the state's n and p.
    const auto axis = moscap_axis(40);
    const auto d = moscap(*mesh::make_tensor_grid(axis));
    const auto s = *assemble::make_scaling(d);
    const physics::Semiconductor si = physics::silicon();
    const double VT = s.V_T, Ns = s.Ns, L_D = s.L_D;
    const double nie = physics::intrinsic_density(si, 300.0) / Ns;  // no narrowing at 1e17
    const double C = -N_A / Ns;
    const double psi_b = std::asinh(C / (2.0 * nie));
    const double eps_s = 11.7 * base::eps0_F_per_cm, eps_ox = 3.9 * base::eps0_F_per_cm;
    const double kappa = eps_ox * L_D / (eps_s * t_ox), Cox = eps_ox / t_ox;
    const double Eg = physics::band_gap_eV(si, 300.0), chi = 4.05;
    const double Vfb = (chi - (chi + 0.5 * Eg - psi_b * VT)) - base::q_C * Q_f / Cox;
    const double Vg = 0.7;
    const std::vector<double> bias{Vg, 0.0};

    auto poisson = *EquilibriumPoisson::create(d, s);
    auto dd = *DriftDiffusion::create(d, s);
    REQUIRE(poisson.set_bias(bias).has_value());
    REQUIRE(dd.set_bias(bias).has_value());
    const std::size_t n = axis.size();
    Noise noise{11};
    std::vector<double> psi(n);
    for (std::size_t i = 0; i < n; ++i) {
        psi[i] = psi_b + 8.0 * std::exp(-axis[i] / 2e-5) + noise.next();
    }
    std::vector<double> x = dd.state_from_potential(psi);
    for (std::size_t i = 0; i < n; ++i) {
        psi[i] = x[3 * i];  // the ohmic contact node is stamped to its Dirichlet value
        x[3 * i + 1] *= 1.0 + 0.3 * noise.next();
        x[3 * i + 2] *= 1.0 + 0.3 * noise.next();
    }
    std::vector<double> f_eq(n), f_dd(3 * n);
    poisson.residual(psi, f_eq);
    dd.residual(x, f_dd);

    double worst_eq = 0.0, worst_dd = 0.0;
    for (std::size_t i = 0; i + 1 < n; ++i) {  // the last row is the ohmic contact's
        const double h_r = (axis[i + 1] - axis[i]) / L_D;
        const double h_l = i > 0 ? (axis[i] - axis[i - 1]) / L_D : 0.0;
        const double dV = 0.5 * (h_l + h_r);
        double flux = (psi[i + 1] - psi[i]) / h_r;
        if (i > 0) flux -= (psi[i] - psi[i - 1]) / h_l;
        if (i == 0) flux += kappa * (Vg - Vfb - (psi[0] - psi_b) * VT) / VT;
        const double rho_eq = nie * std::exp(psi[i]) - nie * std::exp(-psi[i]) - C;
        const double rho_dd = x[3 * i + 1] - x[3 * i + 2] - C;
        const double row_scale = std::abs(flux) + std::abs(dV * rho_eq) + std::abs(dV * rho_dd);
        worst_eq = std::max(worst_eq, std::abs(f_eq[i] - (flux - dV * rho_eq)) / row_scale);
        worst_dd = std::max(worst_dd, std::abs(f_dd[3 * i] - (flux - dV * rho_dd)) / row_scale);
    }
    CAPTURE(worst_eq, worst_dd);
    REQUIRE(worst_eq <= 1e-12);
    REQUIRE(worst_dd <= 1e-12);
}

TEST_CASE("gate: no carrier flux crosses the gate; only the Poisson row changes") {
    const auto axis = moscap_axis(30);
    const auto gated = moscap(*mesh::make_tensor_grid(axis));
    const auto bare = moscap(*mesh::make_tensor_grid(axis), false);
    const auto s = *assemble::make_scaling(gated);
    auto a = *DriftDiffusion::create(gated, s);
    auto b = *DriftDiffusion::create(bare, s);
    const std::vector<double> gate_bias{1.2, 0.0}, bare_bias{0.0};
    REQUIRE(a.set_bias(gate_bias).has_value());
    REQUIRE(b.set_bias(bare_bias).has_value());
    const auto [probe, x] = probe_state(gated, gate_bias, 5);
    std::vector<double> fa(a.unknowns()), fb(b.unknowns());
    a.residual(x, fa);
    b.residual(x, fb);
    const assemble::GateTerm g = assemble::gate_term(stack("x_min"), physics::silicon(), 1.0, 1, s);
    for (std::size_t k = 0; k < fa.size(); ++k) {
        CAPTURE(k);
        if (k == 0) {  // the gate node's Poisson row gains the oxide term
            const double psi_G = 1.2 / s.V_T - g.offset;
            REQUIRE(std::abs(fa[0] - fb[0] - (g.coupling * (psi_G - x[0]) + g.sheet_charge)) <=
                    1e-13 * std::abs(fa[0] - fb[0]));
        } else {
            REQUIRE(fa[k] == fb[k]);
        }
    }
    // A gate carries no current; the gate charge is the oxide flux.
    const auto currents = a.terminal_currents(x);
    REQUIRE(currents[0] == 0.0);
    const auto charges = a.gate_charges(x);
    REQUIRE(charges[1] == 0.0);
    REQUIRE(charges[0] == g.coupling * (1.2 / s.V_T - g.offset - x[0]));
}

TEST_CASE("gate: FD-Jacobian gate with gates in 1D, 2D and 3D (section 10, <= 5e-5)") {
    // Drift-diffusion: a biased MOS-C in 1D; the transistor in 2D and 3D, gate and drain biased.
    const auto [s1, x1] = probe_state(moscap(*mesh::make_tensor_grid(moscap_axis(40))),
                                      {0.8, 0.0}, 21);
    const std::vector<double> on{0.0, 1.0, 0.1, 0.0};
    const auto [s2, x2] = probe_state(transistor(), on, 22);
    const auto [s3, x3] = probe_state(transistor({0.0, 1e-5, 2.5e-5}), on, 23);
    const double e1 = fd_jacobian_error(s1, x1);
    const double e2 = fd_jacobian_error(s2, x2);
    const double e3 = fd_jacobian_error(s3, x3);
    // Equilibrium Poisson with a biased gate.
    const auto d = transistor();
    auto poisson = *EquilibriumPoisson::create(d, *assemble::make_scaling(d));
    const std::vector<double> gate_only{0.0, 1.0, 0.0, 0.0};
    REQUIRE(poisson.set_bias(gate_only).has_value());
    std::vector<double> psi = poisson.charge_neutral_potential();
    Noise noise{24};
    for (double& p : psi) p += 0.5 * noise.next();
    const double e4 = fd_jacobian_error(poisson, psi);
    UNSCOPED_INFO("worst column errors: 1D " << e1 << ", 2D " << e2 << ", 3D " << e3
                                             << ", Poisson " << e4);
    REQUIRE(e1 <= 5e-5);
    REQUIRE(e2 <= 5e-5);
    REQUIRE(e3 <= 5e-5);
    REQUIRE(e4 <= 5e-5);
}

TEST_CASE("gate: the equilibrium system takes gate biases, not ohmic ones") {
    const auto d = transistor();
    auto poisson = *EquilibriumPoisson::create(d, *assemble::make_scaling(d));
    const auto rejected = [&](std::vector<double> bias) {
        const auto r = poisson.set_bias(bias);
        REQUIRE_FALSE(r.has_value());
        REQUIRE(r.error().code == base::ErrorCode::invalid_input);
        return r.error();
    };
    rejected({0.0, 1.0, 0.0});
    REQUIRE(rejected({0.0, std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0})
                .context->index == 1);
    const auto drain = rejected({0.0, 1.0, 0.1, 0.0});
    REQUIRE(drain.context->index == 2);
    REQUIRE(drain.message.find("ohmic") != std::string::npos);

    // A gate bias of dV moves only the gate rows, by G dV / V_T.
    const auto s = *assemble::make_scaling(d);
    const std::vector<double> psi = poisson.charge_neutral_potential();
    std::vector<double> f0(psi.size()), f1(psi.size());
    poisson.residual(psi, f0);
    const std::vector<double> biased{0.0, 0.4, 0.0, 0.0};
    REQUIRE(poisson.set_bias(biased).has_value());
    poisson.residual(psi, f1);
    const auto gate = d.contacts()[1];
    const auto charges = poisson.gate_charges(psi);
    double coupling_sum = 0.0;
    for (std::size_t i = 0; i < psi.size(); ++i) {
        const bool on_gate = std::binary_search(gate.nodes.begin(), gate.nodes.end(),
                                                static_cast<mesh::NodeId>(i));
        if (!on_gate) {
            REQUIRE(f1[i] == f0[i]);
            continue;
        }
        const double G = (f1[i] - f0[i]) / (0.4 / s.V_T);
        REQUIRE(G > 0.0);
        coupling_sum += G;
    }
    REQUIRE(charges[0] == 0.0);
    REQUIRE(charges[2] == 0.0);
    REQUIRE(charges[3] == 0.0);
    REQUIRE(charges[1] != 0.0);
    // The gate couples through its whole face: the scaled areas sum to the gate width over L_D.
    const double width = 7 * 2.5e-6;  // gate nodes at x = 0.75e-5 ... 2.25e-5, cells 2.5e-6 wide
    CAPTURE(coupling_sum, width);
    const double kappa = 3.9 * s.L_D / (11.7 * t_ox);
    REQUIRE(std::abs(coupling_sum - kappa * width / s.L_D) <= 1e-12 * coupling_sum);
}
