// MOS capacitor through the gate contact (ARCHITECTURE.md section 11, Unit 12).
//
// Fixture: the legacy tests/test_cv_physics_validation.py PARAMS: p-type substrate N_A = 1e17,
// t_ox = 5 nm, n+ polysilicon gate, 300 K, on the legacy MOS-C mesh (moscap.MOSCapacitor: L = 2 um,
// 1200 nodes graded towards the interface), gate swept from -2 to 2 V in 50 mV steps. The gate is
// contact 0 on x_min, the substrate an ohmic contact on x_max. As in the legacy, each point is
// thermal equilibrium (Equations::equilibrium_poisson): the quasi-static (low-frequency) C-V.
//
// Gates: NiTCAD reproduces the ported legacy solve (legacy_moscap.hpp) on the same mesh; the legacy
// C-V physics checks P1-P9, the temperature check and the fixed-charge shift, applied to NiTCAD's
// curve; the work-function shift of the electrode; Gauss's law; y-uniform 2D and 3D reproduce 1D;
// the drift-diffusion system with the same gate gives the same state where it converges (below
// threshold), and stalls above it (the finding recorded in ARCHITECTURE.md).
#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/analysis/curve.hpp"
#include "NiTCAD/analysis/cv.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/bandgap_narrowing.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/results/solution.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/equilibrium.hpp"
#include "legacy_moscap.hpp"

using namespace NiTCAD;

namespace {

constexpr double N_A = 1e17;
constexpr double t_ox = 5e-7;

std::vector<double> range(double first, double last, double step) {
    std::vector<double> v;
    for (int k = 0;; ++k) {
        const double x = first + k * step;
        if (x > last + 1e-9) break;
        v.push_back(x);
    }
    return v;
}

std::vector<double> legacy_axis() { return LegacyMOSCapacitor(-N_A, t_ox).x; }

device::GateStack stack(device::GateElectrode electrode = device::GateElectrode::n_poly,
                        double Qf = 0.0, double work_function_eV = 0.0) {
    return {.boundary = "x_min",
            .oxide_thickness_cm = t_ox,
            .electrode = electrode,
            .work_function_eV = work_function_eV,
            .fixed_charge_cm2 = Qf};
}

// A uniformly doped p-type MOS-C on the given mesh: gate (contact 0) on x_min, substrate (contact
// 1) on x_max.
device::Device moscap(mesh::Mesh m, device::GateStack gate = stack(), double T = 300.0) {
    const std::size_t n = m.node_count();
    auto gate_nodes = m.find_boundary("x_min")->nodes;
    auto substrate = m.find_boundary("x_max")->nodes;
    device::Contact g{"gate", device::ContactKind::gate, std::move(gate_nodes)};
    g.gate = std::move(gate);
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

device::Device moscap_1d(device::GateStack gate = stack(), double T = 300.0) {
    return moscap(*mesh::make_tensor_grid(LegacyMOSCapacitor(-N_A, t_ox, {}, 0.0, T).x),
                  std::move(gate), T);
}

struct Curve {
    std::vector<double> Vg, phi_s, Qg, C;
    std::vector<results::BiasPoint> points;
};

solve::BiasOptions quasi_static() {
    solve::BiasOptions o;
    o.equations = solve::Equations::equilibrium_poisson;
    return o;
}

// The quasi-static C-V: phi_s = potential(surface) - potential(bulk contact), the gate charge, and
// its numpy.gradient, as moscap.cv_sweep.
Curve cv(const device::Device& d, const std::vector<double>& Vg,
         const solve::BiasOptions& options = quasi_static()) {
    std::vector<std::vector<double>> points;
    for (const double v : Vg) points.push_back({v, 0.0});
    auto sweep = solve::sweep_bias(d, points, options);
    REQUIRE(sweep.has_value());
    if (sweep->stopped) {
        CAPTURE(sweep->stopped->message, sweep->points.size(),
                sweep->unfinished ? sweep->unfinished->iterations.size() : 0);
        FAIL("the C-V sweep stopped");
    }
    Curve c;
    c.Vg = Vg;
    const std::size_t bulk = d.contacts()[1].nodes.front();
    for (results::BiasPoint& p : sweep->points) {
        c.phi_s.push_back(p.fields.potential_V[0] - p.fields.potential_V[bulk]);
        c.Qg.push_back(p.gate_charge[0]);
        c.points.push_back(std::move(p));
    }
    if (Vg.size() > 1) c.C = numpy_gradient(c.Qg, Vg);  // a single point has no capacitance
    return c;
}

// Computed once: the legacy fixture's NiTCAD curve and the legacy port's.
const Curve& fixture_curve() {
    static const Curve c = cv(moscap_1d(), range(-2.0, 2.0, 0.05));
    return c;
}

const LegacyMOSCapacitor& legacy() {
    static const LegacyMOSCapacitor m(-N_A, t_ox);
    return m;
}

double median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    const std::size_t n = v.size();
    return n % 2 == 1 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

// Linear interpolation of the first crossing of f = target, rising; nullopt if none.
std::optional<double> rising_crossing(const std::vector<double>& x, const std::vector<double>& f,
                                      double target) {
    for (std::size_t i = 0; i + 1 < x.size(); ++i) {
        if (f[i] < target && target <= f[i + 1]) {
            return x[i] + (target - f[i]) * (x[i + 1] - x[i]) / (f[i + 1] - f[i]);
        }
    }
    return std::nullopt;
}

// (k T / 2) ln(Nv / Nc), 1.04 mV at 300 K: NiTCAD measures the gate potential from the intrinsic
// level, the legacy from midgap (assemble/gate.hpp), so NiTCAD's curve is the legacy one moved to
// higher gate voltage by this much.
double intrinsic_shift(double T = 300.0) {
    const physics::Semiconductor si = physics::silicon();
    return 0.5 * base::thermal_voltage(T) *
           std::log(physics::valence_band_dos(si, T) / physics::conduction_band_dos(si, T));
}

}  // namespace

TEST_CASE("mos: the gate contact reproduces the legacy MOS-C solve, moved by the intrinsic level") {
    const Curve& c = fixture_curve();
    const LegacyMOSCapacitor& m = legacy();
    const double shift = intrinsic_shift();
    std::vector<double> legacy_Vg;
    for (const double v : c.Vg) legacy_Vg.push_back(v - shift);
    const LegacyCV ref = m.cv_sweep(legacy_Vg);
    double phi_error = 0.0, q_error = 0.0, c_error = 0.0, psi_error = 0.0;
    for (std::size_t k = 0; k < c.Vg.size(); ++k) {
        phi_error = std::max(phi_error, std::abs(c.phi_s[k] - ref.phi_s[k]));
        q_error = std::max(q_error, std::abs(c.Qg[k] - ref.Qg[k]) / m.Cox);
        c_error = std::max(c_error, std::abs(c.C[k] - ref.C[k]) / m.Cox);
    }
    for (const double v : {-1.5, -0.5, 0.0, 0.5, 1.0, 2.0}) {
        const auto k = static_cast<std::size_t>(std::lround((v + 2.0) / 0.05));
        const std::vector<double> psi = m.solve_psi(v - shift);
        const auto& potential = c.points[k].fields.potential_V;
        for (std::size_t i = 0; i < psi.size(); ++i) {
            psi_error = std::max(psi_error, std::abs(potential[i] - psi[i] * m.VT));
        }
    }
    CAPTURE(shift, phi_error, q_error, c_error, psi_error);
    REQUIRE(phi_error < 1e-12);  // V
    REQUIRE(psi_error < 1e-12);  // V, every node at six biases
    REQUIRE(q_error < 1e-12);    // V (charge / C_ox)
    REQUIRE(c_error < 1e-12);    // relative to C_ox
}

TEST_CASE("mos: drift-diffusion with the gate gives the same state below threshold") {
    // Accumulation and depletion, where the drift-diffusion Newton converges: the quasi-Fermi
    // potentials are zero, so both equation sets solve the same Poisson problem with the same gate
    // row. Starting from equilibrium at 0 V the sweep goes down to -2 V.
    const device::Device d = moscap_1d();
    const auto Vg = range(-2.0, -0.3, 0.1);
    std::vector<double> down(Vg.rbegin(), Vg.rend());
    const Curve qs = cv(d, down);
    const Curve dd = cv(d, down, {});
    const physics::Semiconductor si = physics::silicon();
    const double nie = physics::effective_intrinsic_density(si, N_A, 300.0);
    const double VT = base::thermal_voltage(300.0);
    double psi_error = 0.0, q_error = 0.0, hole_error = 0.0, current = 0.0;
    for (std::size_t k = 0; k < down.size(); ++k) {
        const auto& a = dd.points[k].fields;
        const auto& b = qs.points[k].fields;
        for (std::size_t i = 0; i < a.potential_V.size(); ++i) {
            psi_error = std::max(psi_error, std::abs(a.potential_V[i] - b.potential_V[i]));
            // Majority carriers in equilibrium: p = n_ie e^(-psi / V_T).
            hole_error = std::max(
                hole_error, std::abs(a.p_cm3[i] / (nie * std::exp(-a.potential_V[i] / VT)) - 1.0));
        }
        q_error = std::max(q_error, std::abs(dd.Qg[k] - qs.Qg[k]) / legacy().Cox);
        REQUIRE(dd.points[k].terminal_current[0] == 0.0);  // a gate carries no current
        current = std::max(current, std::abs(dd.points[k].terminal_current[1]));
    }
    CAPTURE(psi_error, q_error, hole_error, current);
    REQUIRE(psi_error < 1e-12);
    REQUIRE(q_error < 1e-12);
    REQUIRE(hole_error < 1e-12);
    REQUIRE(current < 1e-12);  // A/cm^2
}

TEST_CASE("mos: drift-diffusion stalls once the inversion layer forms") {
    // The inversion layer has no ohmic contact of its own: its electrons reach the substrate only
    // through the depleted region, where n is some 1e-14 of the inversion density, and the Newton
    // corrections of those minority densities stay at the solver's round-off (ARCHITECTURE.md 6.4).
    // Recorded so that a change in this behaviour is noticed; the quasi-static solve is the
    // supported path for a MOS capacitor.
    std::vector<std::vector<double>> points;
    for (const double v : range(-0.3, 0.5, 0.1)) points.push_back({v, 0.0});
    const auto sweep = solve::sweep_bias(moscap_1d(), points);
    REQUIRE(sweep.has_value());
    REQUIRE(sweep->stopped.has_value());
    REQUIRE(sweep->stopped->code == base::ErrorCode::non_convergence);
    CAPTURE(sweep->points.size());
    REQUIRE(sweep->points.size() >= 2);  // -0.3 and -0.2 V converge
}

TEST_CASE("mos: the quasi-static sweep keeps every ohmic contact at 0 V") {
    const device::Device d = moscap_1d();
    const std::vector<std::vector<double>> points{{0.5, 0.0}, {0.5, 0.1}};
    const auto sweep = solve::sweep_bias(d, points, quasi_static());
    REQUIRE_FALSE(sweep.has_value());
    REQUIRE(sweep.error().code == base::ErrorCode::invalid_input);
    REQUIRE(sweep.error().context->index == 1);  // the point
    REQUIRE(sweep.error().context->value == 0.1);  // the contact's bias
    REQUIRE(sweep.error().message.find("bias point 1, contact 'substrate'") != std::string::npos);
    REQUIRE(sweep.error().message.find("ohmic contact at 0 V") != std::string::npos);
    // Drift-diffusion takes the same points.
    REQUIRE(solve::make_run_record(d, quasi_static(), points).input_identity !=
            solve::make_run_record(d, {}, points).input_identity);
}

TEST_CASE("mos: the quasi-static sweep starts from a potential alone") {
    const device::Device d = moscap_1d();
    const auto start = cv(d, {0.4});
    results::NodeFields potential_only{start.points[0].fields.potential_V, {}, {}};
    const std::vector<std::vector<double>> points{{0.45, 0.0}};
    const auto warm = solve::sweep_bias(d, points, quasi_static(), &potential_only);
    REQUIRE(warm.has_value());
    REQUIRE_FALSE(warm->stopped.has_value());
    const auto cold = cv(d, {0.45});
    REQUIRE(std::abs(warm->points[0].gate_charge[0] - cold.Qg[0]) <= 1e-12 * std::abs(cold.Qg[0]));
    // Drift-diffusion needs the densities too.
    const auto dd = solve::sweep_bias(d, points, {}, &potential_only);
    REQUIRE_FALSE(dd.has_value());
    REQUIRE(dd.error().message.find("positive densities") != std::string::npos);
}

TEST_CASE("mos: the run identity leaves out what the quasi-static sweep ignores") {
    const device::Device d = moscap_1d();
    const std::vector<std::vector<double>> points{{0.5, 0.0}};
    const auto id = [&](const device::Device& dev, const solve::BiasOptions& o,
                        const results::NodeFields* initial = nullptr) {
        return solve::make_run_record(dev, o, points, initial).input_identity;
    };
    solve::BiasOptions no_transport_models = quasi_static();
    no_transport_models.models.srh = false;
    no_transport_models.models.auger = false;
    no_transport_models.models.doping_mobility = false;
    REQUIRE(id(d, no_transport_models) == id(d, quasi_static()));
    solve::BiasOptions no_bgn = quasi_static();
    no_bgn.models.bgn = false;  // band-gap narrowing does enter the equilibrium
    REQUIRE(id(d, no_bgn) != id(d, quasi_static()));
    solve::BiasOptions dd_no_srh;
    dd_no_srh.models.srh = false;
    REQUIRE(id(d, dd_no_srh) != id(d, {}));
    // Initial densities: ignored by the quasi-static sweep, read by drift-diffusion.
    const auto state = cv(d, {0.4}).points[0].fields;
    results::NodeFields other = state;
    other.n_cm3[3] *= 2.0;
    REQUIRE(id(d, quasi_static(), &other) == id(d, quasi_static(), &state));
    REQUIRE(id(d, {}, &other) != id(d, {}, &state));
    // The work function field of a polysilicon gate is not read; a metal gate's is.
    device::GateStack poly = stack();
    poly.work_function_eV = 4.7;
    REQUIRE(id(moscap_1d(poly), quasi_static()) == id(d, quasi_static()));
    REQUIRE(id(moscap_1d(stack(device::GateElectrode::metal, 0.0, 4.7)), quasi_static()) !=
            id(moscap_1d(stack(device::GateElectrode::metal, 0.0, 4.6)), quasi_static()));
}

TEST_CASE("mos: the equilibrium solve reports the zero-bias gate charge") {
    for (const double Qf : {0.0, 2e11}) {
        CAPTURE(Qf);
        const device::Device d = moscap_1d(stack(device::GateElectrode::n_poly, Qf));
        const auto eq = solve::solve_equilibrium(d);
        REQUIRE(eq.has_value());
        REQUIRE(eq->gate_charge.size() == 2);
        REQUIRE(eq->gate_charge[1] == 0.0);  // the substrate is ohmic
        const auto qs = cv(d, {0.0});
        REQUIRE(std::abs(eq->gate_charge[0] - qs.Qg[0]) <= 1e-12 * std::abs(qs.Qg[0]));
        // n+ poly on p-type: 0 V is above flatband (V_FB about -1 V), so the surface is depleted
        // and the electrode charged positively against it.
        REQUIRE(eq->gate_charge[0] > 0.0);
    }
}

TEST_CASE("mos: P1 the potential satisfies the legacy discrete Poisson rows") {
    const Curve& c = fixture_curve();
    const LegacyMOSCapacitor& m = legacy();
    double worst = 0.0;
    for (const double v : {-1.5, -0.5, 0.5, 1.0}) {
        const auto k = static_cast<std::size_t>(std::lround((v + 2.0) / 0.05));
        std::vector<double> psi = c.points[k].fields.potential_V;
        for (double& p : psi) p /= m.VT;
        const std::vector<double> F = m.residual(psi, v - intrinsic_shift());
        for (const double f : F) worst = std::max(worst, std::abs(f));
    }
    CAPTURE(worst);
    REQUIRE(worst < 1e-6);  // legacy gate (interior rows); here every row, the gate row included
}

TEST_CASE("mos: P2 gate charge, fixed charge and semiconductor charge balance (Gauss)") {
    // Discrete global neutrality: the box rows sum to the flux into the substrate contact, which
    // sits in the neutral bulk 2 um away. The legacy compares a trapezoid integral within 2%.
    for (const double Qf : {0.0, 3e11}) {
        CAPTURE(Qf);
        const device::Device d = moscap_1d(stack(device::GateElectrode::n_poly, Qf));
        const Curve c = cv(d, {-1.5, 0.0, 1.5});
        const auto volumes = d.mesh().volumes();
        for (const results::BiasPoint& p : c.points) {
            double semiconductor = 0.0;
            for (std::size_t i = 0; i + 1 < volumes.size(); ++i) {  // all but the contact node
                semiconductor += base::q_C * volumes[i] *
                                 (p.fields.p_cm3[i] - p.fields.n_cm3[i] - N_A);
            }
            const double total = p.gate_charge[0] + base::q_C * Qf + semiconductor;
            CAPTURE(p.bias_V[0], p.gate_charge[0], total);
            REQUIRE(std::abs(total) < 1e-12 * std::abs(p.gate_charge[0]));
        }
    }
}

TEST_CASE("mos: P3 the gate charge matches the surface field of the profile") {
    const Curve& c = fixture_curve();
    const LegacyMOSCapacitor& m = legacy();
    for (const double v : {-1.0, 0.8}) {
        const auto k = static_cast<std::size_t>(std::lround((v + 2.0) / 0.05));
        const auto& psi = c.points[k].fields.potential_V;
        const double E_s = -(psi[1] - psi[0]) / (m.x[1] - m.x[0]);  // V/cm
        const double Q_field = -m.eps_s * E_s;                          // C/cm^2
        const double error = std::abs(Q_field + c.Qg[k]) / std::abs(c.Qg[k]);
        CAPTURE(v, error);
        REQUIRE(error < 0.05);  // legacy gate
    }
}

TEST_CASE("mos: P4 the sweep passes accumulation, flatband, depletion, inversion in order") {
    const Curve& c = fixture_curve();
    const double phiF = legacy().analytic_landmarks().phi_F;
    std::vector<std::string> sequence;
    for (const double ps : c.phi_s) {
        const std::string regime = ps < -0.2 * phiF            ? "accumulation"
                                   : std::abs(ps) <= 0.25 * phiF ? "flatband"
                                   : ps < 1.8 * phiF             ? "depletion"
                                                                 : "inversion";
        if (std::find(sequence.begin(), sequence.end(), regime) == sequence.end()) {
            sequence.push_back(regime);
        }
    }
    REQUIRE(sequence ==
            std::vector<std::string>{"accumulation", "flatband", "depletion", "inversion"});
}

TEST_CASE("mos: P5 series capacitance and doping recovered from the depletion C-V") {
    const Curve& coarse = fixture_curve();
    const LegacyMOSCapacitor& m = legacy();
    const double phiF = m.analytic_landmarks().phi_F;
    // The depletion window on the coarse sweep, re-swept in 5 mV steps (legacy).
    std::vector<double> window;
    for (std::size_t k = 0; k < coarse.Vg.size(); ++k) {
        if (0.05 * phiF < coarse.phi_s[k] && coarse.phi_s[k] < 0.95 * phiF) {
            window.push_back(coarse.Vg[k]);
        }
    }
    REQUIRE_FALSE(window.empty());
    const Curve fine = cv(moscap_1d(), range(window.front(), window.back(), 0.005));
    std::vector<double> N_points, relation;
    for (std::size_t k = 0; k < fine.Vg.size(); ++k) {
        const double ps = fine.phi_s[k];
        if (!(0.30 * phiF < ps && ps < 0.65 * phiF)) continue;
        const double c_dep = 1.0 / (1.0 / fine.C[k] - 1.0 / m.Cox);
        N_points.push_back(2.0 * ps * c_dep * c_dep / (base::q_C * m.eps_s));
        const double W = std::sqrt(2.0 * m.eps_s * ps / (base::q_C * N_A));
        const double c_series = 1.0 / (1.0 / m.Cox + W / m.eps_s);
        relation.push_back(std::abs(c_series - fine.C[k]) / fine.C[k]);
    }
    REQUIRE(N_points.size() >= 5);
    CAPTURE(N_points, relation);
    REQUIRE(std::abs(median(N_points) - N_A) / N_A < 0.20);
    for (std::size_t k = 0; k < N_points.size(); ++k) {
        REQUIRE(N_points[k] > 0.70 * N_A);
        REQUIRE(N_points[k] < 1.50 * N_A);
        if (k > 0) REQUIRE(N_points[k] - N_points[k - 1] < 1e-12);  // converges down towards N
    }
    REQUIRE(median(relation) < 0.10);
    REQUIRE(*std::max_element(relation.begin(), relation.end()) < 0.35);
}

TEST_CASE("mos: the analysis extractors read the quasi-static C-V (Unit 24)") {
    // The quasi-static capacitance from the sweep's gate charge is the curve above (numpy.gradient);
    // the accumulation capacitance approaches C_ox from below; the flat-band voltage from the
    // flat-band capacitance (Debye length) is the analytic one (as P7, 50 mV); the doping profile
    // from d(1/C_d^2)/dV in depletion (behind C_ox) approaches N_A from above (see below).
    const Curve& c = fixture_curve();
    const LegacyMOSCapacitor& m = legacy();
    const auto lm = m.analytic_landmarks();
    results::Sweep sweep;
    sweep.points = c.points;
    const auto charge = analysis::gate_charge_curve(sweep, 0, 0);
    REQUIRE(charge.has_value());
    const auto C = analysis::quasi_static_capacitance(*charge);
    REQUIRE(C.has_value());
    for (std::size_t k = 0; k < c.C.size(); ++k) {
        REQUIRE(std::abs(C->y[k] - c.C[k]) <= 1e-12 * std::abs(c.C[k]));
    }
    const auto acc = analysis::accumulation_capacitance(*C);
    REQUIRE(acc->value < m.Cox);
    REQUIRE(acc->value > 0.9 * m.Cox);
    const double C_fb = analysis::flat_band_capacitance(m.Cox, m.eps_s, N_A, 300.0);
    const auto V_fb = analysis::flat_band_voltage(*C, C_fb);
    REQUIRE(V_fb.has_value());

    const double phiF = lm.phi_F;
    std::vector<double> window;
    for (std::size_t k = 0; k < c.Vg.size(); ++k) {
        if (0.05 * phiF < c.phi_s[k] && c.phi_s[k] < 0.95 * phiF) window.push_back(c.Vg[k]);
    }
    const Curve fine = cv(moscap_1d(), range(window.front(), window.back(), 0.005));
    std::vector<double> V_dep;
    for (std::size_t k = 0; k < fine.Vg.size(); ++k) {
        if (0.30 * phiF < fine.phi_s[k] && fine.phi_s[k] < 0.65 * phiF) V_dep.push_back(fine.Vg[k]);
    }
    const auto fine_C = analysis::make_curve(fine.Vg, fine.C);
    const auto profile = analysis::doping_profile(
        *analysis::slice(*fine_C, V_dep.front(), V_dep.back()), m.eps_s, m.Cox);
    REQUIRE(profile.has_value());
    const double N_median = median(profile->doping_cm3);
    std::printf("mos C-V: C_acc / C_ox %.4f, V_FB %.4f V (analytic %.4f), N median %.4e over %zu "
                "points, from %.4e to %.4e, depth %.3e to %.3e cm\n",
                acc->value / m.Cox, V_fb->value, lm.V_FB, N_median, profile->doping_cm3.size(),
                profile->doping_cm3.front(), profile->doping_cm3.back(), profile->depth_cm.front(),
                profile->depth_cm.back());
    REQUIRE(std::abs(V_fb->value - lm.V_FB) < 0.05);
    // The 1/C^2 profile assumes a depletion edge sharp on the scale of the depletion width; here
    // the window (P5's) is in weak depletion, w = 2.8 to 4.3 Debye lengths (L_D = 12.9 nm), and the
    // profile reads high, falling towards N_A with depth (measured 1.48 to 1.27 N_A, median 1.33).
    // A junction in reverse bias recovers its doping within 1% (small_signal_test.cpp).
    for (std::size_t k = 0; k < profile->doping_cm3.size(); ++k) {
        REQUIRE(profile->doping_cm3[k] > N_A);
        if (k > 0) REQUIRE(profile->doping_cm3[k] < profile->doping_cm3[k - 1]);
    }
    REQUIRE(N_median < 1.4 * N_A);
}

TEST_CASE("mos: P6 C_min and W_max match the values derived from the parameters") {
    const Curve& c = fixture_curve();
    const auto lm = legacy().analytic_landmarks();
    const double c_min = *std::min_element(c.C.begin(), c.C.end());
    const double eps_s = legacy().eps_s;
    CAPTURE(c_min, lm.C_min);
    REQUIRE(std::abs(c_min - lm.C_min) / lm.C_min < 0.15);
    const double W = (1.0 / c_min - 1.0 / lm.C_ox) * eps_s;
    REQUIRE(std::abs(W - lm.W_max) / lm.W_max < 0.15);
}

TEST_CASE("mos: P7 flatband and threshold crossings, and the low-frequency rebound") {
    const Curve& c = fixture_curve();
    const auto lm = legacy().analytic_landmarks();
    const auto flatband = rising_crossing(c.Vg, c.phi_s, 0.0);
    REQUIRE(flatband.has_value());
    CAPTURE(*flatband, lm.V_FB);
    REQUIRE(std::abs(*flatband - lm.V_FB) < 0.05);
    const auto threshold = rising_crossing(c.Vg, c.phi_s, 2.0 * lm.phi_F);
    REQUIRE(threshold.has_value());
    CAPTURE(*threshold, lm.V_th);
    REQUIRE(std::abs(*threshold - lm.V_th) < 0.05);
    // Past C_min the capacitance rises again (the inversion layer follows the gate).
    const auto i_min = static_cast<std::size_t>(
        std::min_element(c.C.begin(), c.C.end()) - c.C.begin());
    const double slope_floor = 0.05 * c.C[i_min] / std::max(c.Vg.back() - c.Vg[i_min], 1e-9);
    bool rebound = false;
    for (std::size_t k = i_min; k + 1 < c.C.size(); ++k) {
        rebound = rebound || (c.C[k + 1] - c.C[k]) / (c.Vg[k + 1] - c.Vg[k]) > slope_floor;
    }
    REQUIRE(rebound);
}

TEST_CASE("mos: P8 quasi-static versus high-frequency behaviour in inversion") {
    const Curve& c = fixture_curve();
    const LegacyMOSCapacitor& m = legacy();
    const auto lm = m.analytic_landmarks();
    // HF approximation: the surface potential saturates at 2 phi_F.
    const auto c_hf = [&](double ps) {
        ps = std::min(std::abs(ps), 2.0 * lm.phi_F);
        const double w = std::sqrt(2.0 * m.eps_s * ps / (base::q_C * N_A));
        return 1.0 / (1.0 / m.Cox + w / m.eps_s);
    };
    std::vector<double> excess, depletion;
    for (std::size_t k = 0; k < c.Vg.size(); ++k) {
        const double hf = c_hf(c.phi_s[k]);
        if (c.Vg[k] > lm.V_th + 0.4) {
            excess.push_back((c.C[k] - hf) / hf);
        } else if (0.25 * lm.phi_F < c.phi_s[k] && c.phi_s[k] < 0.8 * lm.phi_F) {
            depletion.push_back(std::abs(c.C[k] - hf) / hf);
        }
    }
    REQUIRE_FALSE(excess.empty());
    REQUIRE(*std::min_element(excess.begin(), excess.end()) > -0.05);
    double mean = 0.0;
    for (const double e : excess) mean += e / static_cast<double>(excess.size());
    REQUIRE(mean > 0.05);
    REQUIRE_FALSE(depletion.empty());
    REQUIRE(median(depletion) < 0.15);
    REQUIRE(*std::max_element(depletion.begin(), depletion.end()) < 0.40);
}

TEST_CASE("mos: P9 the curve is converged in the voltage step and the Newton tolerance") {
    const device::Device d = moscap_1d();
    const auto coarse_v = range(-1.5, 1.5, 0.1);
    const auto fine_v = range(-1.5, 1.5, 0.05);
    const Curve coarse = cv(d, coarse_v);
    const Curve fine = cv(d, fine_v);
    std::vector<double> rel;
    for (std::size_t k = 1; k + 1 < coarse_v.size(); ++k) {
        const double interpolated = fine.C[2 * k];  // the fine grid contains the coarse points
        rel.push_back(std::abs(interpolated - coarse.C[k]) / std::abs(coarse.C[k]));
    }
    REQUIRE(median(rel) < 0.02);
    const auto trapezoid = [](const std::vector<double>& x, const std::vector<double>& f) {
        double s = 0.0;
        for (std::size_t k = 1; k < x.size(); ++k) s += 0.5 * (f[k] + f[k - 1]) * (x[k] - x[k - 1]);
        return s;
    };
    const double q_coarse = trapezoid(coarse_v, coarse.C), q_fine = trapezoid(fine_v, fine.C);
    REQUIRE(std::abs(q_fine - q_coarse) / std::abs(q_coarse) < 0.01);

    // A tighter Newton tolerance does not move the surface potential (legacy gate 1e-9 V).
    solve::BiasOptions tight = quasi_static();
    tight.newton.tol_update = 1e-12;
    tight.newton.max_iterations = 400;
    const double loose = cv(d, {1.0}).phi_s[0];
    const double strict = cv(d, {1.0}, tight).phi_s[0];
    CAPTURE(loose - strict);
    REQUIRE(std::abs(loose - strict) < 1e-9);
}

TEST_CASE("mos: C_min follows its temperature-dependent landmark") {
    std::vector<double> c_min_landmark;
    for (const double T : {275.0, 300.0, 350.0}) {
        CAPTURE(T);
        const LegacyMOSCapacitor m(-N_A, t_ox, {}, 0.0, T);
        const auto lm = m.analytic_landmarks();
        const Curve c = cv(moscap_1d(stack(), T), range(-2.0, 2.0, 0.05));
        const double c_min = *std::min_element(c.C.begin(), c.C.end());
        CAPTURE(c_min, lm.C_min);
        REQUIRE(std::abs(c_min - lm.C_min) / lm.C_min < 0.15);
        c_min_landmark.push_back(lm.C_min);
    }
    // Higher T: smaller phi_F, thinner W_max, higher C_min.
    REQUIRE(c_min_landmark[0] < c_min_landmark[1]);
    REQUIRE(c_min_landmark[1] < c_min_landmark[2]);
}

TEST_CASE("mos: fixed oxide charge shifts the curve by -q Q_f / C_ox") {
    // The legacy checks V_FB(Q_f) - V_FB(0) = -q Q_f / C_ox to 1e-6. Here the whole state shifts:
    // phi_s(V_G; Q_f) = phi_s(V_G + q Q_f / C_ox; 0), and the gate charge by -q Q_f.
    const double Qf = 1e11;
    const double Cox = legacy().Cox;
    const double shift = base::q_C * Qf / Cox;
    const std::vector<double> Vg{-1.0, 0.2, 1.2};
    std::vector<double> shifted;
    for (const double v : Vg) shifted.push_back(v + shift);
    const Curve charged = cv(moscap_1d(stack(device::GateElectrode::n_poly, Qf)), Vg);
    const Curve clean = cv(moscap_1d(), shifted);
    for (std::size_t k = 0; k < Vg.size(); ++k) {
        CAPTURE(Vg[k], charged.phi_s[k] - clean.phi_s[k]);
        REQUIRE(std::abs(charged.phi_s[k] - clean.phi_s[k]) < 1e-12);
        REQUIRE(std::abs(charged.Qg[k] - (clean.Qg[k] - base::q_C * Qf)) < 1e-12 * Cox);
    }
}

TEST_CASE("mos: the electrode work function shifts the curve rigidly") {
    // phi_m(p+) - phi_m(n+) = Eg(T); a metal at chi + 0.3 eV shifts by 0.3 V.
    const double Eg = physics::band_gap_eV(physics::silicon(), 300.0);
    const double chi = physics::silicon_parameters.electron_affinity_eV;
    const std::vector<double> Vg{-0.5, 0.4, 1.4};
    const Curve n_poly = cv(moscap_1d(), Vg);
    for (const auto& [gate, offset] :
         {std::pair{stack(device::GateElectrode::p_poly), Eg},
          std::pair{stack(device::GateElectrode::metal, 0.0, chi + 0.3), 0.3}}) {
        std::vector<double> shifted;
        for (const double v : Vg) shifted.push_back(v + offset);
        const Curve c = cv(moscap_1d(gate), shifted);
        for (std::size_t k = 0; k < Vg.size(); ++k) {
            CAPTURE(offset, Vg[k], c.phi_s[k] - n_poly.phi_s[k]);
            REQUIRE(std::abs(c.phi_s[k] - n_poly.phi_s[k]) < 1e-12);
            REQUIRE(std::abs(c.Qg[k] - n_poly.Qg[k]) < 1e-12 * legacy().Cox);
        }
    }
    // And the legacy flatband voltage of each electrode is where phi_s crosses zero.
    REQUIRE(std::abs(legacy_flatband_voltage(-N_A, t_ox, {LegacyGate::p_poly}, 0.0, 300.0) -
                     legacy_flatband_voltage(-N_A, t_ox, {}, 0.0, 300.0) - Eg) < 1e-12);
}

TEST_CASE("mos: a y-uniform 2D and 3D MOS-C reproduces 1D") {
    // Shorter mesh than the fixture (the same grading law, 300 nodes) to keep 2D and 3D small; the
    // transverse axes are non-uniform, so the gate's face areas differ from node to node.
    const LegacyMOSCapacitor m(-N_A, t_ox, {}, 0.0, 300.0, 2e-4, 300);
    const std::vector<double> y{0.0, 1e-5, 3e-5, 3.5e-5};
    const std::vector<double> z{0.0, 2e-5, 2.5e-5};
    const double width = y.back(), area = y.back() * z.back();
    const std::vector<double> Vg{-1.0, 0.0, 0.6, 1.5};
    const Curve one = cv(moscap(*mesh::make_tensor_grid(m.x)), Vg);
    const Curve two = cv(moscap(*mesh::make_tensor_grid(m.x, y)), Vg);
    const Curve three = cv(moscap(*mesh::make_tensor_grid(m.x, y, z)), Vg);
    const std::size_t nx = m.x.size();
    for (std::size_t k = 0; k < Vg.size(); ++k) {
        CAPTURE(Vg[k]);
        for (std::size_t i = 0; i < nx * y.size(); ++i) {
            REQUIRE(std::abs(two.points[k].fields.potential_V[i] -
                             one.points[k].fields.potential_V[i % nx]) < 1e-12);
        }
        for (std::size_t i = 0; i < nx * y.size() * z.size(); ++i) {
            REQUIRE(std::abs(three.points[k].fields.potential_V[i] -
                             one.points[k].fields.potential_V[i % nx]) < 1e-12);
        }
        REQUIRE(std::abs(two.Qg[k] - one.Qg[k] * width) <= 1e-12 * std::abs(one.Qg[k] * width));
        REQUIRE(std::abs(three.Qg[k] - one.Qg[k] * area) <= 1e-12 * std::abs(one.Qg[k] * area));
    }
}

