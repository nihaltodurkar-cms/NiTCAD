// Meshed oxide, electrodes, interface charge, traps and surface recombination through the solve
// layer (ARCHITECTURE.md section 11, Unit 15b).
//
// MOS fixture as Unit 12's (legacy test_cv_physics_validation PARAMS): p-type silicon N_A = 1e17,
// t_ox = 5 nm SiO2, n+ polysilicon electrode, 300 K, the legacy MOS-C silicon mesh (L = 2 um,
// x = L expm1(6 s) / expm1(6)). The meshed device puts the oxide on [-t_ox, -h/2] (20 cells) and
// the silicon nodes at x + h/2, h the legacy first spacing, so the interface lies at x = 0 midway
// between the nodes -h/2 and h/2 (mesh::straddle_interface's convention). Each point is thermal
// equilibrium (Equations::equilibrium_poisson), the quasi-static C-V.
//
// Reference: the exact 1D solution (the first integral of Poisson's equation with Boltzmann
// statistics and uniform doping), with the fixed charge and the traps' charge at the interface:
//     V_G = V_FB0 + phi_s + Q_G / C_ox,   Q_G = -(Q_s(phi_s) + q Q_f + Q_it(phi_s)),
// Q_s the semiconductor charge, Q_it the traps' charge by the Fermi function integrated in closed
// form, V_FB0 = psi_b V_T - (chi + E_c - E_i) + phi_m.
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/mesh.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/interface_traps.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/results/solution.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/equilibrium.hpp"
#include "legacy_moscap.hpp"

using namespace NiTCAD;

namespace {

constexpr double N_A = 1e17;
constexpr double t_ox = 5e-7;
constexpr double T = 300.0;

std::vector<double> range(double first, double last, double step) {
    std::vector<double> v;
    for (int k = 0;; ++k) {
        const double x = first + k * step;
        if (x > last + 1e-9) break;
        v.push_back(x);
    }
    return v;
}

std::vector<double> legacy_axis(std::size_t nx = 1200) {
    return LegacyMOSCapacitor(-N_A, t_ox, {}, 0.0, T, 2e-4, nx).x;
}

// The meshed axis: oxide on [-t_ox, -h/2], silicon at the legacy nodes moved by h/2.
std::vector<double> meshed_axis(std::size_t nx = 1200, int oxide_cells = 20) {
    const std::vector<double> si = legacy_axis(nx);
    const double h = si[1];
    std::vector<double> x;
    for (int k = 0; k <= oxide_cells; ++k) {
        x.push_back(-t_ox + (t_ox - 0.5 * h) * k / oxide_cells);
    }
    for (const double v : si) x.push_back(v + 0.5 * h);
    return x;
}

double psi_bulk() {  // scaled neutral potential of the substrate
    return -std::asinh(N_A / (2.0 * physics::intrinsic_density(physics::silicon(), T)));
}

// E_c - E_i and E_v - E_i of silicon [eV].
std::pair<double, double> gap_edges() {
    const physics::Semiconductor si = physics::silicon();
    const double ec =
        physics::intrinsic_level_depth_eV(si, T) - si.parameters().electron_affinity_eV;
    return {ec, ec - physics::band_gap_eV(si, T)};
}

struct Stack {
    double Qf = 0.0;   // [cm^-2]
    double Dit = 0.0;  // [cm^-2 eV^-1]: acceptor-like above the flat-band Fermi level, donor-like
                       // below, so the traps hold no charge at flat band (the legacy M14 model)
    device::Electrode electrode{};
};

device::Interface oxide_interface(const Stack& s) {
    device::Interface f{"oxide", "silicon"};
    f.fixed_charge_cm2 = s.Qf;
    if (s.Dit > 0.0) {
        const auto [ec, ev] = gap_edges();
        const double E0 = psi_bulk() * base::thermal_voltage(T);  // E_F - E_i in the bulk [eV]
        f.traps.bands = {{.type = physics::TrapType::acceptor, .density_cm2_eV = s.Dit,
                          .energy_low_eV = E0, .energy_high_eV = ec},
                         {.type = physics::TrapType::donor, .density_cm2_eV = s.Dit,
                          .energy_low_eV = ev, .energy_high_eV = E0}};
    }
    return f;
}

// A MOS stack on a mesh whose x < 0 nodes are oxide: the electrode on the x_min nodes, the
// substrate contact on x_max.
device::Device meshed(mesh::Mesh m, const Stack& s = {}) {
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n);
    std::vector<double> acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const bool oxide = m.points()[i][0] < 0.0;
        region[i] = oxide ? 0 : 1;
        if (!oxide) acceptors[i] = N_A;
    }
    device::Contact gate{"gate", device::ContactKind::electrode, m.find_boundary("x_min")->nodes};
    gate.electrode = s.electrode;
    auto substrate = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = T,
         .regions = {{"oxide", physics::silicon_dioxide()}, {"silicon", physics::silicon()}},
         .node_region = std::move(region),
         .donors = std::vector<double>(n, 0.0),
         .acceptors = std::move(acceptors),
         .contacts = {std::move(gate),
                      {"substrate", device::ContactKind::ohmic, std::move(substrate)}},
         .interfaces = {oxide_interface(s)}});
}

device::Device meshed_1d(std::size_t nx = 1200, const Stack& s = {}) {
    return meshed(*mesh::make_tensor_grid(meshed_axis(nx)), s);
}

// Unit 12's lumped gate on the legacy silicon mesh.
device::Device lumped_1d(std::size_t nx = 1200, double Qf = 0.0) {
    mesh::Mesh m = *mesh::make_tensor_grid(legacy_axis(nx));
    const std::size_t n = m.node_count();
    device::Contact g{"gate", device::ContactKind::gate, m.find_boundary("x_min")->nodes};
    g.gate = {.boundary = "x_min", .oxide_thickness_cm = t_ox, .fixed_charge_cm2 = Qf};
    auto substrate = m.find_boundary("x_max")->nodes;
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = T,
         .regions = {{"silicon", physics::silicon()}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::vector<double>(n, 0.0),
         .acceptors = std::vector<double>(n, N_A),
         .contacts = {std::move(g),
                      {"substrate", device::ContactKind::ohmic, std::move(substrate)}}});
}

solve::BiasOptions quasi_static() {
    solve::BiasOptions o;
    o.equations = solve::Equations::equilibrium_poisson;
    return o;
}

std::vector<results::BiasPoint> sweep(const device::Device& d, const std::vector<double>& Vg,
                                      const solve::BiasOptions& o = quasi_static()) {
    std::vector<std::vector<double>> points;
    for (const double v : Vg) points.push_back({v, 0.0});
    auto s = solve::sweep_bias(d, points, o);
    REQUIRE(s.has_value());
    if (s->stopped) {
        CAPTURE(s->stopped->message, s->points.size());
        FAIL("the sweep stopped");
    }
    return std::move(s->points);
}

double softplus(double z) {
    return z > 0.0 ? z + std::log1p(std::exp(-z)) : std::log1p(std::exp(z));
}

// The exact 1D MOS-C.
class ExactMOS {
public:
    explicit ExactMOS(const Stack& s) : s_(s) {
        const physics::Semiconductor si = physics::silicon();
        VT_ = base::thermal_voltage(T);
        eps_ = si.parameters().eps_r * base::eps0_F_per_cm;
        Cox_ = physics::silicon_dioxide_parameters.eps_r * base::eps0_F_per_cm / t_ox;
        const double ni = physics::intrinsic_density(si, T);
        psi_b_ = psi_bulk();
        p0_ = ni * std::exp(-psi_b_);
        n0_ = ni * std::exp(psi_b_);
        const double chi = si.parameters().electron_affinity_eV;
        const double phi_m = s.electrode.kind == device::GateElectrode::n_poly ? chi
                             : s.electrode.kind == device::GateElectrode::p_poly
                                 ? chi + physics::band_gap_eV(si, T)
                                 : s.electrode.work_function_eV;
        Vfb0_ = psi_b_ * VT_ - physics::intrinsic_level_depth_eV(si, T) + phi_m;
        std::tie(ec_, ev_) = gap_edges();
    }

    [[nodiscard]] double Cox() const { return Cox_; }
    [[nodiscard]] double Vfb0() const { return Vfb0_; }

    // Semiconductor charge [C/cm^2] at surface potential phi_s [V].
    [[nodiscard]] double Qs(double phi_s) const {
        const double u = phi_s / VT_;
        const double kT = base::k_B_J_per_K * T;
        const double F = p0_ * (std::expm1(-u) + u) + n0_ * (std::expm1(u) - u);
        return -(u > 0.0 ? 1.0 : u < 0.0 ? -1.0 : 0.0) * std::sqrt(2.0 * eps_ * kT * F);
    }
    // Trapped charge [C/cm^2]: q D_it (empty donors - occupied acceptors), the Fermi function
    // integrated in closed form.
    [[nodiscard]] double Qit(double phi_s) const {
        if (s_.Dit == 0.0) return 0.0;
        const double eF = (psi_b_ + phi_s / VT_) * VT_;  // surface E_F - E_i [eV]
        const double E0 = psi_b_ * VT_;
        const auto occupied = [&](double a, double b) {  // int_a^b f dE
            return VT_ * (softplus((eF - a) / VT_) - softplus((eF - b) / VT_));
        };
        const double empty_donors = (E0 - ev_) - occupied(ev_, E0);
        return base::q_C * s_.Dit * (empty_donors - occupied(E0, ec_));
    }
    [[nodiscard]] double gate_charge(double phi_s) const {
        return -(Qs(phi_s) + base::q_C * s_.Qf + Qit(phi_s));
    }
    [[nodiscard]] double gate_voltage(double phi_s) const {
        return Vfb0_ + phi_s + gate_charge(phi_s) / Cox_;
    }
    // phi_s at V_G (bisection; V_G rises with phi_s).
    [[nodiscard]] double surface_potential(double Vg) const {
        double lo = -2.0, hi = 2.0;
        for (int k = 0; k < 200; ++k) {
            const double mid = 0.5 * (lo + hi);
            (gate_voltage(mid) < Vg ? lo : hi) = mid;
        }
        return 0.5 * (lo + hi);
    }

private:
    Stack s_;
    double VT_, eps_, Cox_, psi_b_, p0_, n0_, Vfb0_, ec_ = 0, ev_ = 0;
};

// Largest |Q_G - Q_G,exact| / C_ox [V] over a sweep.
double charge_error(const std::vector<results::BiasPoint>& points, const ExactMOS& exact) {
    double worst = 0.0;
    for (const results::BiasPoint& p : points) {
        const double Q = exact.gate_charge(exact.surface_potential(p.bias_V[0]));
        worst = std::max(worst, std::abs(p.gate_charge[0] - Q) / exact.Cox());
    }
    return worst;
}

std::size_t first_silicon(const device::Device& d) { return d.reference_node(); }

}  // namespace

TEST_CASE("insulator: the meshed MOS-C converges to the exact solution, as the lumped one") {
    // Gate charge against the exact solution over -2 to 2 V, on the legacy mesh with nx nodes and
    // its refinements; the lumped oxide of Unit 12 on the same silicon nodes for comparison. In a
    // charge-free oxide the potential is linear, which the box method reproduces exactly: the
    // meshed and lumped oxides differ by where the silicon's first control volume sits.
    const std::vector<double> Vg = range(-2.0, 2.0, 0.1);
    const ExactMOS exact({});
    std::vector<double> e_meshed, e_lumped, between;
    for (const std::size_t nx : {300, 600, 1200, 2400}) {
        const auto m = sweep(meshed_1d(nx), Vg);
        const auto l = sweep(lumped_1d(nx), Vg);
        e_meshed.push_back(charge_error(m, exact));
        e_lumped.push_back(charge_error(l, exact));
        double d = 0.0;
        for (std::size_t k = 0; k < Vg.size(); ++k) {
            d = std::max(d, std::abs(m[k].gate_charge[0] - l[k].gate_charge[0]) / exact.Cox());
        }
        between.push_back(d);
    }
    CAPTURE(e_meshed, e_lumped, between);
    // Second order in the mesh, for both. Measured [V]: meshed 5.0e-4, 1.24e-4, 3.07e-5, 7.6e-6;
    // lumped 6.6e-4, 1.7e-4, 4.3e-5, 1.08e-5; meshed against lumped 1.16e-3, 2.9e-4, 7.4e-5,
    // 1.8e-5.
    for (std::size_t k = 1; k < e_meshed.size(); ++k) {
        REQUIRE(e_meshed[k] < 0.3 * e_meshed[k - 1]);
        REQUIRE(e_lumped[k] < 0.3 * e_lumped[k - 1]);
        REQUIRE(between[k] < 0.3 * between[k - 1]);
    }
    REQUIRE(e_meshed.back() < 1e-5);
}

TEST_CASE("insulator: the fixed charge shifts the meshed C-V by -q Q_f / C_ox") {
    // With Q_f the whole curve is the exact one with Q_f, and the flat-band voltage (phi_s = 0)
    // moves by -q Q_f / C_ox.
    const std::vector<double> Vg = range(-2.0, 2.0, 0.1);
    for (const double Qf : {5e11, -1e12}) {
        const ExactMOS exact({.Qf = Qf});
        const auto points = sweep(meshed_1d(1200, {.Qf = Qf}), Vg);
        const double e = charge_error(points, exact);
        // Flat band (phi_s = 0 at the first silicon node) by secant iteration, from the exact
        // flat-band voltage.
        const double shift = -base::q_C * Qf / exact.Cox();
        const device::Device d = meshed_1d(1200, {.Qf = Qf});
        const std::size_t surface = first_silicon(d);
        const std::size_t bulk = d.contacts()[1].nodes.front();
        const auto phi_s = [&](double v) {
            const auto p = sweep(d, {v});
            return p[0].fields.potential_V[surface] - p[0].fields.potential_V[bulk];
        };
        double v0 = ExactMOS({}).Vfb0() + shift - 1e-3, v1 = v0 + 2e-3;
        double f0 = phi_s(v0), f1 = phi_s(v1);
        for (int k = 0; k < 20 && std::abs(f1) > 1e-12; ++k) {
            const double v2 = v1 - f1 * (v1 - v0) / (f1 - f0);
            v0 = v1;
            f0 = f1;
            v1 = v2;
            f1 = phi_s(v1);
        }
        const double measured = v1 - ExactMOS({}).Vfb0();
        CAPTURE(Qf, e, shift, measured, (measured - shift) / shift);
        REQUIRE(std::abs(f1) <= 1e-12);
        // The fixed charge sits at the interface (interface_edges.hpp). At flat band the silicon
        // is neutral and field-free, so the oxide carries the whole charge over exactly t_ox and
        // the discrete shift is -q Q_f / C_ox itself.
        REQUIRE(std::abs(measured - shift) <= 1e-9 * std::abs(shift));
        REQUIRE(e <= 1e-4);
    }
}

TEST_CASE("insulator: interface traps against the exact solution and the legacy stretch-out") {
    // A uniform D_it across the gap, neutral at flat band: the exact solution with the Fermi
    // function, and the legacy M14 model (moscap D_it: Q_it = -q D_it phi_s), which the Fermi
    // occupancy reduces to while the surface Fermi level is away from the band edges.
    const double Dit = 2e12;
    const std::vector<double> Vg = range(-2.0, 2.0, 0.1);
    const ExactMOS exact({.Dit = Dit});
    const device::Device d = meshed_1d(1200, {.Dit = Dit});
    const auto points = sweep(d, Vg);
    const double e = charge_error(points, exact);
    // The traps sit at the interface potential (interface_edges.hpp): second order in the mesh,
    // as without traps.
    std::vector<double> e_mesh;
    for (const std::size_t nx : {300, 600, 2400}) {
        e_mesh.push_back(charge_error(sweep(meshed_1d(nx, {.Dit = Dit}), Vg), exact));
    }
    e_mesh.insert(e_mesh.begin() + 2, e);  // 300, 600, 1200, 2400
    // The traps hold -q D_it phi_s in mid-gap, so the trapped charge reported is that.
    const std::size_t surface = first_silicon(d);
    const std::size_t bulk = d.contacts()[1].nodes.front();
    double worst_trapped = 0.0;
    for (const results::BiasPoint& p : points) {
        const double phi_s = exact.surface_potential(p.bias_V[0]);
        REQUIRE(p.interface_trap_charge.size() == 1);
        worst_trapped =
            std::max(worst_trapped, std::abs(p.interface_trap_charge[0] - exact.Qit(phi_s)));
    }
    // The legacy port with D_it, on its own mesh, its gate voltage moved by the intrinsic level as
    // in Unit 12; compared where the surface Fermi level is mid-gap (phi_s from 0 to 0.7 V: E_F
    // within 0.3 eV of E_i).
    LegacyMOSCapacitor legacy(-N_A, t_ox);
    legacy.set_interface_traps(Dit);
    const double shift = 0.5 * base::thermal_voltage(T) *
                         std::log(physics::valence_band_dos(physics::silicon(), T) /
                                  physics::conduction_band_dos(physics::silicon(), T));
    std::vector<double> legacy_Vg;
    for (const double v : Vg) legacy_Vg.push_back(v - shift);
    const LegacyCV ref = legacy.cv_sweep(legacy_Vg);
    double worst_legacy = 0.0, without_traps = 0.0;
    const auto plain = sweep(meshed_1d(1200), Vg);
    for (std::size_t k = 0; k < Vg.size(); ++k) {
        const double phi =
            points[k].fields.potential_V[surface] - points[k].fields.potential_V[bulk];
        if (ref.phi_s[k] < 0.0 || ref.phi_s[k] > 0.7) continue;
        worst_legacy = std::max(worst_legacy, std::abs(phi - ref.phi_s[k]));
        const double phi0 =
            plain[k].fields.potential_V[surface] - plain[k].fields.potential_V[bulk];
        without_traps = std::max(without_traps, std::abs(phi0 - ref.phi_s[k]));
    }
    CAPTURE(e_mesh, worst_trapped, worst_legacy, without_traps);
    for (std::size_t k = 1; k < e_mesh.size(); ++k) REQUIRE(e_mesh[k] < 0.3 * e_mesh[k - 1]);
    REQUIRE(e_mesh.back() < 2e-5);
    REQUIRE(worst_trapped < 1e-3 * base::q_C * Dit);  // within 1 mV of trapped charge
    REQUIRE(worst_legacy < 2e-3);
    REQUIRE(without_traps > 50.0 * worst_legacy);  // the stretch-out is resolved
}

TEST_CASE("insulator: Gauss's law with the meshed electrode, fixed charge and traps") {
    // Q_G + q sum_i V_i (p - n - N_A) + q Q_f + Q_it = 0 over the silicon nodes, at each point;
    // quasi-static and (below threshold) drift-diffusion. Also with a one-cell oxide: the electrode
    // node at -t_ox, the first silicon node at +t_ox, so the electrode is an end of the interface
    // edge and its charge is the half-edge flux.
    const Stack s{.Qf = 3e11, .Dit = 1e12};
    std::vector<double> one_cell{-t_ox};
    for (const double v : legacy_axis(600)) one_cell.push_back(v + t_ox);
    const device::Device thick = meshed_1d(1200, s);
    const device::Device thin = meshed(*mesh::make_tensor_grid(one_cell), s);
    std::vector<std::pair<const device::Device*, std::vector<results::BiasPoint>>> runs;
    runs.emplace_back(&thick, sweep(thick, range(-1.5, 1.5, 0.5)));
    runs.emplace_back(&thin, sweep(thin, range(-1.5, 1.5, 0.5)));
    // On this mesh the silicon half-cell at the interface is 5 nm, so in accumulation the
    // interface hole density is the node's times e^(psi_s - psi_I), about 1e8 here, and the trap
    // terms are steep in psi: the drift-diffusion Jacobian's pivot ratio falls to 1.2e-10 (6.9e-5
    // with the traps at the node, 1.5e-5 on the legacy mesh with its 0.0125 nm half-cell), close to
    // the default singularity check of 1e-11. It converges without the check.
    solve::BiasOptions unchecked;
    unchecked.linear.min_pivot_ratio = 0.0;
    runs.emplace_back(&thin, sweep(thin, range(-1.5, -0.5, 0.5), unchecked));
    for (const auto& [device, points] : runs) {
        const device::Device& d = *device;
        const auto volumes = d.mesh().volumes();
        for (const results::BiasPoint& p : points) {
            double Q = 0.0;
            for (std::size_t i = 0; i < volumes.size(); ++i) {
                if (d.is_insulator(static_cast<mesh::NodeId>(i)) || i + 1 == volumes.size()) {
                    continue;
                }
                Q += base::q_C * volumes[i] * (p.fields.p_cm3[i] - p.fields.n_cm3[i] - N_A);
            }
            const double total =
                p.gate_charge[0] + Q + base::q_C * s.Qf + p.interface_trap_charge[0];
            CAPTURE(d.mesh().node_count(), p.bias_V[0], p.gate_charge[0], Q,
                    p.interface_trap_charge[0], total);
            REQUIRE(std::abs(total) <= 1e-9 * std::abs(p.gate_charge[0]) + 1e-15);
        }
    }
}

TEST_CASE("insulator: a y-uniform 2D and 3D meshed MOS-C reproduces 1D") {
    const Stack s{.Qf = 2e11, .Dit = 5e11};
    const std::vector<double> Vg{-1.0, 0.0, 0.8};
    const auto x = meshed_axis(300);
    const auto one = sweep(meshed(*mesh::make_tensor_grid(x), s), Vg);
    const std::vector<double> y{0.0, 1e-5, 3e-5};
    const auto two = sweep(meshed(*mesh::make_tensor_grid(x, y), s), Vg);
    const auto three = sweep(meshed(*mesh::make_tensor_grid(x, y, y), s), Vg);
    const double width = 3e-5;
    for (std::size_t k = 0; k < Vg.size(); ++k) {
        CAPTURE(Vg[k], one[k].gate_charge[0], two[k].gate_charge[0] / width,
                three[k].gate_charge[0] / (width * width));
        REQUIRE(std::abs(two[k].gate_charge[0] / width - one[k].gate_charge[0]) <=
                1e-10 * std::abs(one[k].gate_charge[0]));
        REQUIRE(std::abs(three[k].gate_charge[0] / (width * width) - one[k].gate_charge[0]) <=
                1e-10 * std::abs(one[k].gate_charge[0]));
        REQUIRE(std::abs(two[k].interface_trap_charge[0] / width -
                         one[k].interface_trap_charge[0]) <=
                1e-10 * std::abs(one[k].interface_trap_charge[0]));
    }
}

TEST_CASE("insulator: drift-diffusion with the meshed oxide gives the quasi-static state") {
    // Below threshold (as Unit 12's lumped gate): the coupled system at each gate bias with the
    // substrate at 0 V is the equilibrium one, traps included.
    const Stack s{.Qf = 3e11, .Dit = 1e12};
    const device::Device d = meshed_1d(600, s);
    const std::vector<double> Vg = range(-1.5, 0.0, 0.5);
    const auto qs = sweep(d, Vg);
    const auto dd = sweep(d, Vg, {});
    for (std::size_t k = 0; k < Vg.size(); ++k) {
        double dpsi = 0.0;
        for (std::size_t i = 0; i < qs[k].fields.potential_V.size(); ++i) {
            dpsi = std::max(dpsi,
                            std::abs(qs[k].fields.potential_V[i] - dd[k].fields.potential_V[i]));
        }
        CAPTURE(Vg[k], dpsi, qs[k].gate_charge[0], dd[k].gate_charge[0]);
        REQUIRE(dpsi < 1e-9);
        REQUIRE(std::abs(qs[k].gate_charge[0] - dd[k].gate_charge[0]) <=
                1e-8 * std::abs(qs[k].gate_charge[0]));
        REQUIRE(std::abs(qs[k].interface_trap_charge[0] - dd[k].interface_trap_charge[0]) <=
                1e-8 * std::abs(qs[k].interface_trap_charge[0]) + 1e-20);
        REQUIRE(dd[k].terminal_current[0] == 0.0);
        // Insulator nodes carry no carriers and no band diagram.
        REQUIRE(dd[k].fields.n_cm3[0] == 0.0);
        REQUIRE(std::isnan(dd[k].bands.conduction_eV[0]));
    }
}

TEST_CASE("insulator: the electrode work function shifts the meshed curve rigidly") {
    // A metal electrode at phi_m moves the curve by phi_m - chi against the n+ poly one.
    const std::vector<double> Vg = range(-1.0, 1.0, 0.5);
    const auto poly = sweep(meshed_1d(600), Vg);
    const double chi = physics::silicon_parameters.electron_affinity_eV;
    const double phi_m = 4.75;
    std::vector<double> moved;
    for (const double v : Vg) moved.push_back(v + (phi_m - chi));
    const auto metal = sweep(
        meshed_1d(600, {.electrode = {.kind = device::GateElectrode::metal,
                                      .work_function_eV = phi_m}}),
        moved);
    for (std::size_t k = 0; k < Vg.size(); ++k) {
        CAPTURE(Vg[k], poly[k].gate_charge[0], metal[k].gate_charge[0]);
        REQUIRE(std::abs(poly[k].gate_charge[0] - metal[k].gate_charge[0]) <=
                1e-9 * std::abs(poly[k].gate_charge[0]) + 1e-15);
    }
}

namespace {

// A p+n diode in 2D under an oxide (Unit 15b recombination gate): emitter p+ 1e19 on x < 1 um,
// n base 1e17 to x = 21 um with the base contact on x_max and the anode on x_min (silicon nodes
// only); the silicon is w = 1 um thick in y, under 10 nm of SiO2 (the interface at y = w midway
// between nodes). The base beyond x = a = 2 um meets the oxide through an interface with surface
// recombination velocity S; the 1 um "neck" next to the junction (which holds the space-charge
// region) and the emitter meet it through none. SRH, Auger and band-gap narrowing off: holes in the
// base recombine only at the surface.
constexpr double W_E = 1e-4, a_neck = 2e-4, W_total = 21e-4, t_cap = 1e-6;
constexpr double N_emitter = 1e19, N_base = 1e17;

device::Device surface_diode(double S, double w) {
    std::vector<double> x;
    for (int k = 0; k < 16; ++k) x.push_back(k * 0.05e-4);       // 0 to 0.75 um
    for (int k = 0; k < 60; ++k) x.push_back(0.8e-4 + k * 0.01e-4);  // 0.8 to 1.39 um
    // 1.5 to 20.9 um: the neck-base boundary at a = 2 um lies midway between the nodes 1.9 and
    // 2.1 um, so the recombining boxes start exactly there.
    for (int k = 0; k < 98; ++k) x.push_back(1.5e-4 + k * 0.2e-4);
    x.push_back(W_total);
    std::vector<double> y;
    const double dy = w / 6.0;
    for (int k = 0; k < 6; ++k) y.push_back(k * dy);
    y.push_back(w - 0.5 * dy);  // the last silicon row, dy / 2 below the interface
    y.push_back(w + 0.5 * dy);  // the first oxide row
    y.push_back(w + 0.5 * dy + t_cap);
    auto grid = mesh::make_tensor_grid(x, y);
    REQUIRE(grid.has_value());
    mesh::Mesh m = std::move(*grid);
    const std::size_t n = m.node_count();
    std::vector<device::RegionId> region(n);
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto& p = m.points()[i];
        if (p[1] > w) {
            region[i] = 3;  // oxide
        } else if (p[0] < W_E) {
            region[i] = 0;  // emitter
            acceptors[i] = N_emitter;
        } else {
            region[i] = p[0] < a_neck ? 1 : 2;  // neck, base
            donors[i] = N_base;
        }
    }
    std::vector<mesh::NodeId> anode, cathode;
    for (const mesh::NodeId v : m.find_boundary("x_min")->nodes) {
        if (m.points()[static_cast<std::size_t>(v)][1] < w) anode.push_back(v);
    }
    for (const mesh::NodeId v : m.find_boundary("x_max")->nodes) {
        if (m.points()[static_cast<std::size_t>(v)][1] < w) cathode.push_back(v);
    }
    device::Interface f{"base", "oxide"};
    f.recombination_velocity_n_cm_s = S;
    f.recombination_velocity_p_cm_s = S;
    auto d = device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = T,
         .regions = {{"emitter", physics::silicon()},
                     {"neck", physics::silicon()},
                     {"base", physics::silicon()},
                     {"oxide", physics::silicon_dioxide()}},
         .node_region = std::move(region),
         .donors = std::move(donors),
         .acceptors = std::move(acceptors),
         .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                      {"cathode", device::ContactKind::ohmic, std::move(cathode)}},
         .interfaces = {std::move(f)}});
    INFO((d ? std::string() : d.error().message));
    REQUIRE(d.has_value());
    return std::move(*d);
}

// The low-injection current density [A/cm^2] of that diode. Holes: from the depletion edge
// x_n(V) a linear profile through the neck (no recombination), then in the base the lowest
// transverse mode of D p'' = 0 with -D dp/dy = S p at y = w and none at y = 0: p' ~ cos(k y),
// k tan(k w) = S / D, decaying as sinh(k (W - x)) to the ohmic contact:
//     J_p = q D_p p_n0 (e^(V/V_T) - 1) / ((a - x_n) + tanh(k (W - a)) / k).
// Electrons in the emitter: linear to the anode, J_n = q D_n n_p0 (e^(V/V_T) - 1) / W_E (the
// emitter side of the depletion region is 1e-4 of the base side, neglected). x_n from the
// depletion approximation of the one-sided junction.
double surface_diode_current(double S, double V, double w) {
    const physics::Semiconductor si = physics::silicon();
    const double VT = base::thermal_voltage(T);
    const double ni = physics::intrinsic_density(si, T);
    const double Dp =
        physics::caughey_thomas_mobility(si, physics::Carrier::hole, N_base, T) * VT;
    const double Dn =
        physics::caughey_thomas_mobility(si, physics::Carrier::electron, N_emitter, T) * VT;
    const double eps = si.parameters().eps_r * base::eps0_F_per_cm;
    const double Vbi =
        VT * (std::asinh(N_base / (2.0 * ni)) + std::asinh(N_emitter / (2.0 * ni)));
    const double xn = std::sqrt(2.0 * eps * (Vbi - V) * N_emitter /
                                (base::q_C * N_base * (N_emitter + N_base)));
    double tail = W_total - a_neck;  // tanh(k L) / k, the k -> 0 limit when S = 0
    if (S > 0.0) {
        double lo = 0.0, hi = 0.5 * 3.14159265358979323846 / w;  // k tan(k w) = S / D
        for (int it = 0; it < 200; ++it) {
            const double mid = 0.5 * (lo + hi);
            (mid * std::tan(mid * w) < S / Dp ? lo : hi) = mid;
        }
        const double k = 0.5 * (lo + hi);
        tail = std::tanh(k * (W_total - a_neck)) / k;
    }
    const double excess = std::expm1(V / VT);
    const double Jp = base::q_C * Dp * (ni * ni / N_base) * excess / ((a_neck - W_E - xn) + tail);
    const double Jn = base::q_C * Dn * (ni * ni / N_emitter) * excess / W_E;
    return Jp + Jn;
}

}  // namespace

TEST_CASE("insulator: surface recombination under the oxide against the analytic diode") {
    solve::BiasOptions o;
    o.models = {.srh = false, .auger = false, .bgn = false};
    const std::vector<double> V{0.5, 0.6};
    std::vector<double> errors;
    const double w = 1e-4;
    for (const double S : {0.0, 1e3, 3e3}) {
        const device::Device d = surface_diode(S, w);
        std::vector<std::vector<double>> points;
        for (const double v : V) points.push_back({v, 0.0});
        const auto s = solve::sweep_bias(d, points, o);
        REQUIRE(s.has_value());
        REQUIRE_FALSE(s->stopped.has_value());
        for (std::size_t k = 0; k < V.size(); ++k) {
            const double J = s->points[k].terminal_current[0] / w;  // A/cm per cm of width
            const double ref = surface_diode_current(S, V[k], w);
            CAPTURE(S, V[k], J, ref, J / ref - 1.0);
            errors.push_back(J / ref - 1.0);
            // Measured within 2.7e-4 at S = 0 and 8.3e-4 with recombination (on a 0.1 um base
            // mesh the remainder grew with S / w: 4e-3 at w = 0.25 um, S = 3e3 cm/s).
            REQUIRE(std::abs(J / ref - 1.0) < 2e-3);
        }
    }
    // Recombination raises the current several times over at S = 3e3 cm/s.
    REQUIRE(surface_diode_current(3e3, 0.6, w) > 3.0 * surface_diode_current(0.0, 0.6, w));
    CAPTURE(errors);
}

TEST_CASE("insulator: the run identity covers insulators, electrodes and interface data") {
    const std::vector<std::vector<double>> point{{0.0, 0.0}};
    const auto identity = [&](const device::Device& d, const solve::BiasOptions& o = {}) {
        return solve::make_run_record(d, o, point, nullptr).input_identity;
    };
    const auto base_id = identity(meshed_1d(300));
    REQUIRE(identity(meshed_1d(300)) == base_id);
    REQUIRE(identity(meshed_1d(300, {.Qf = 1e10})) != base_id);
    REQUIRE(identity(meshed_1d(300, {.Dit = 1e11})) != base_id);
    REQUIRE(identity(meshed_1d(300, {.electrode = {.kind = device::GateElectrode::p_poly}})) !=
            base_id);
    REQUIRE(identity(meshed_1d(300, {.electrode = {.kind = device::GateElectrode::metal,
                                                   .work_function_eV = 4.5}})) !=
            identity(meshed_1d(300, {.electrode = {.kind = device::GateElectrode::metal,
                                                   .work_function_eV = 4.6}})));
    // The trap capture data, the thermal velocity and the recombination velocities.
    const auto with = [&](auto change) {
        device::Interface f = oxide_interface({.Dit = 1e11});
        change(f);
        mesh::Mesh m = *mesh::make_tensor_grid(meshed_axis(300));
        const device::Device ref = meshed(std::move(m), {.Dit = 1e11});
        device::DeviceDescription desc{.mesh = ref.mesh(),
                                       .temperature_K = T,
                                       .regions = {ref.regions().begin(), ref.regions().end()},
                                       .node_region = {ref.node_region().begin(),
                                                       ref.node_region().end()},
                                       .donors = {ref.donors().begin(), ref.donors().end()},
                                       .acceptors = {ref.acceptors().begin(),
                                                     ref.acceptors().end()},
                                       .contacts = {ref.contacts().begin(), ref.contacts().end()},
                                       .interfaces = {f}};
        return identity(*device::Device::create(std::move(desc)));
    };
    const auto traps_id = with([](device::Interface&) {});
    REQUIRE(with([](device::Interface& f) { f.traps.bands[0].sigma_n_cm2 = 2e-15; }) != traps_id);
    REQUIRE(with([](device::Interface& f) { f.traps.thermal_velocity_n_cm_s = 2e7; }) != traps_id);
    REQUIRE(with([](device::Interface& f) { f.traps.thermal_velocity_p_cm_s = 2e7; }) != traps_id);
    REQUIRE(with([](device::Interface& f) { f.recombination_velocity_n_cm_s = 10.0; }) !=
            traps_id);
    // The insulator's permittivity.
    {
        mesh::Mesh m = *mesh::make_tensor_grid(meshed_axis(300));
        const device::Device ref = meshed(std::move(m));
        std::vector<device::Region> regions(ref.regions().begin(), ref.regions().end());
        regions[0].material = *physics::Insulator::create({.eps_r = 7.5});
        const device::Device other = *device::Device::create(
            {.mesh = ref.mesh(),
             .temperature_K = T,
             .regions = std::move(regions),
             .node_region = {ref.node_region().begin(), ref.node_region().end()},
             .donors = {ref.donors().begin(), ref.donors().end()},
             .acceptors = {ref.acceptors().begin(), ref.acceptors().end()},
             .contacts = {ref.contacts().begin(), ref.contacts().end()},
             .interfaces = {ref.interfaces().begin(), ref.interfaces().end()}});
        REQUIRE(identity(other) != base_id);
    }
}
