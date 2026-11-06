// TVB-13, the uniform resistor with field mobility and with energy balance (DECISIONS Unit 55 §6, 16.2 (f), 16.4;
// criteria 16.8, frozen; the layer bound 16.18 (a), DT-Q15 A; construction 16.42).
//
// Code under test: solve_bias / sweep_bias on a uniform 1D n-bar (N_D = 1e17 cm^-3) with ohmic ends, 300 K,
// Boltzmann, no band-gap narrowing, no generation-recombination, constant low-field mobilities (the materials'
// mu_max: doping_mobility off), silicon's Canali parameters.
// (i) Drift (D; fast and full tier): field_mobility on, L = 1 um, N = 64 nodes, V / L = 1e2, 1e3, 1e4, 1e5 V/cm.
//     n = n0 everywhere and E uniform, so J = q (n0 v_n(E) + p0 v_p(E)) with v(E) = mu_C(E) E, Canali's law
//     (verify/reference/model_mobility), E the contacts' potential difference over L. Judged: |J - J_ref| / J_ref
//     <= 64 eps N^2.
// (ii) Energy balance (D; full tier): energy_balance on (field_mobility off: ED4 refuses both), N = 128 and 512 nodes,
//     V / L = 1e3, 1e4 V/cm, L = 20 um. In the interior the uniform state solves the discrete equations exactly: T_n -
//     T_L = (2/3)(q/k) tau_w mu_C(E) E^2 (ED3, ED4's E_eff(T) = E). The contacts hold T_n = T_L, so there is a
//     boundary layer; its length lambda = 2 / min |Re s| over the linearized quartic's roots
//     (verify/reference/energy_layer). Judged on the nodes at distance d from the nearer contact with
//     e^{-d / lambda} <= 1e-13: max |rise - rise_ref| / rise_ref <= 64 eps N^2. Applicable (else NOT APPLICABLE)
//     while there is at least one such node. Printed, not judged: the decay rate of |rise - rise_ref| measured on
//     the solved profile near the left contact, beside min |Re s|.
// No order is computed (16.2 (f)).
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/energy_balance.hpp"
#include "NiTCAD/physics/field_mobility.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/verify/benchmark/device_data.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/benchmark/tvb32.hpp"
#include "NiTCAD/verify/reference/energy_layer.hpp"
#include "NiTCAD/verify/reference/model_mobility.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

constexpr double T = 300.0, ND = 1e17, eps_m = 0x1p-52;
constexpr double drift_length = 1e-4;
constexpr int drift_nodes = 64;
constexpr double drift_fields[4] = {1e2, 1e3, 1e4, 1e5};
// The lengths, chosen from lambda before any NiTCAD run (16.42): the judged region starts 30 lambda from each
// contact (2.9 um and 6.1 um). 1e5 V/cm is not built (16.42): its layer needs 0.81 mm a side, and at the 20 kV that
// takes, the potential's rounding (8e-11 V_T) is above Newton's tolerance and the D criterion's scale.
constexpr double eb_length[2] = {2e-3, 2e-3};
constexpr double eb_fields[2] = {1e3, 1e4};
constexpr int eb_nodes[2] = {128, 512};
const double tail = 13.0 * std::log(10.0);  // e^{-d / lambda} <= 1e-13

ReferenceInfo info() {
    ReferenceInfo r;
    r.id = "TVB-13";
    r.cls = ReferenceClass::D;
    r.independence = Independence::I1;
    r.provenance = Provenance::PD;
    r.citation = "Canali et al., IEEE Trans. Electron Devices 22, 1045 (1975) (the law, PB as 15.10); the "
                 "homogeneous energy balance and its boundary layer derived (PD, DECISIONS Unit 55 16.18 (a))";
    r.version = "DECISIONS Unit 55 16.2 (f), 16.8, 16.18 (a), 16.42";
    r.equations = {1, 2, 3};
    r.form = "J = q n0 mu_C(E) E; T_n - T_L = (2/3)(q/k) tau_w mu_C(E) E^2 beyond the layer";
    r.assumptions = {"uniform doping, ohmic ends", "constant low-field mobility", "the layer's tail measured"};
    r.derivation = "verify/reference/energy_layer.hpp";
    r.verification = "energy_layer_test: Durand-Kerner, the zero-field limit, the unreduced determinant";
    r.second_method = "the determinant of the unreduced linear system at each root";
    return r;
}

device::Device bar(int N, double length) {
    std::vector<double> x(static_cast<std::size_t>(N));
    for (int k = 0; k < N; ++k) x[static_cast<std::size_t>(k)] = length * k / (N - 1);
    x.back() = length;
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    auto a = m.find_boundary("x_min")->nodes;
    auto c = m.find_boundary("x_max")->nodes;
    return *device::Device::create({.mesh = std::move(m),
                                    .temperature_K = T,
                                    .regions = {{"silicon", physics::silicon()}},
                                    .node_region = std::vector<device::RegionId>(n, 0),
                                    .donors = std::vector<double>(n, ND),
                                    .acceptors = std::vector<double>(n, 0.0),
                                    .contacts = {{"left", device::ContactKind::ohmic, std::move(a)},
                                                 {"right", device::ContactKind::ohmic, std::move(c)}}});
}

solve::BiasOptions options(bool energy) {
    solve::BiasOptions o;
    o.models.doping_mobility = false;
    o.models.hole_mobility = assemble::HoleMobility::caughey_thomas;
    o.models.field_mobility = !energy;
    o.models.energy_balance = energy;
    o.models.srh = o.models.auger = o.models.radiative = false;
    o.models.bgn = false;
    o.models.fermi_dirac = false;
    o.newton.tol_update = 1e-12;
    o.step_control.max_cuts = 20;
    return o;
}

// The homogeneous state at field E: Canali's mobility, the rise and d mu / dT (dd, then double).
struct Homogeneous {
    double mobility, rise, slope, check;  // check: the carrier-temperature mobility at the rise, relative to mu_C
};
Homogeneous homogeneous(double mu0, double v_sat, double beta, double tau_w, double E) {
    namespace mm = reference::model_mobility;
    using reference::DD;
    const DD mu = mm::canali_mobility(DD(mu0), DD(v_sat), DD(beta), DD(E));
    const DD rise = DD(2.0) / DD(3.0) * DD(tau_w) * mu * DD(E) * DD(E) / reference::model_band::detail::k_eV<DD>();
    const mm::BalanceResult b = mm::carrier_temperature_mobility({DD(mu0), rise, DD(tau_w), DD(v_sat)}, beta);
    return {mu.to_double(), rise.to_double(), b.d_rise.to_double(), std::abs(((b.mobility - mu) / mu).to_double())};
}

}  // namespace

BenchmarkResult run_tvb13(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-13";
    r.title = "Uniform resistor: Canali drift and the homogeneous energy balance";
    r.info = info();
    const physics::Semiconductor si = physics::silicon();
    const double mu_n = si.parameters().electron_mobility.mu_max, mu_p = si.parameters().hole_mobility.mu_max;
    const physics::CanaliParameters sat_n = physics::saturation(si, physics::Carrier::electron);
    const physics::CanaliParameters sat_p = physics::saturation(si, physics::Carrier::hole);
    const double tau_n = physics::silicon_energy_relaxation.tau_n_s;  // ED5: silicon's cited value applies

    // ---- The references for (ii): the homogeneous state and its layer ----
    struct Layer {
        Homogeneous h;
        reference::LayerResult L;
    };
    std::vector<Layer> layers;
    if (tier != Tier::fast) {
        for (const double E : eb_fields) {
            const Homogeneous h = homogeneous(mu_n, sat_n.v_sat_cm_s, sat_n.beta, tau_n, E);
            const auto L = reference::energy_layer({.q_C = base::q_C,
                                                    .k_J_per_K = base::k_B_J_per_K,
                                                    .eps_F_per_cm = si.parameters().eps_r * base::eps0_F_per_cm,
                                                    .density_cm3 = ND,
                                                    .field_V_per_cm = E,
                                                    .mobility_cm2_Vs = h.mobility,
                                                    .mobility_slope = h.slope,
                                                    .temperature_K = T + h.rise,
                                                    .lattice_K = T,
                                                    .tau_w_s = tau_n});
            if (!L) {
                r.run_error = "reference: " + L.error();
                return r;
            }
            if (L->determinant_residual > 1e-20 || h.check > 1e-28) {
                r.run_error = "reference: the second method disagrees (determinant " + shortest(L->determinant_residual) +
                              ", carrier-temperature mobility " + shortest(h.check) + ")";
                return r;
            }
            std::string roots;
            for (const auto& s : L->roots) roots += " (" + shortest(s.real()) + ", " + shortest(s.imag()) + ")";
            r.notes.push_back("E " + shortest(E) + " V/cm: mu_C " + shortest(h.mobility) + ", rise " + shortest(h.rise) +
                              " K, d mu / dT " + shortest(h.slope) + "; roots [1/cm]" + roots + "; lambda " +
                              shortest(L->length_cm) + " cm, judged from d >= " + shortest(tail * L->length_cm) +
                              " cm; determinant residual " + shortest(L->determinant_residual));
            layers.push_back({h, *L});
        }
    }

    // ---- (i) Drift ----
    {
        const device::Device d = bar(drift_nodes, drift_length);
        std::vector<std::vector<double>> points;
        for (const double E : drift_fields) points.push_back({0.0, E * drift_length});
        const auto s = solve::sweep_bias(d, points, options(false));
        const double criterion = 64.0 * eps_m * drift_nodes * drift_nodes;
        const std::string crit = "relative <= 64 eps N^2 = " + shortest(criterion);
        if (!s || s->stopped || s->points.size() != points.size()) {
            r.run_error = "drift: " + (!s ? s.error().message : s->stopped ? s->stopped->message : "short sweep");
            return r;
        }
        for (std::size_t k = 0; k < points.size(); ++k) {
            const results::BiasPoint& p = s->points[k];
            const auto& f = p.fields;
            const std::size_t a = static_cast<std::size_t>(d.contacts()[0].nodes.front());
            const std::size_t b = static_cast<std::size_t>(d.contacts()[1].nodes.front());
            const double E = std::abs(f.potential_V[b] - f.potential_V[a]) / drift_length;
            namespace mm = reference::model_mobility;
            using reference::DD;
            const double vn_ref = (mm::canali_mobility(DD(mu_n), DD(sat_n.v_sat_cm_s), DD(sat_n.beta), DD(E)) * DD(E)).to_double();
            const double vp_ref = (mm::canali_mobility(DD(mu_p), DD(sat_p.v_sat_cm_s), DD(sat_p.beta), DD(E)) * DD(E)).to_double();
            const double J_ref = base::q_C * (f.n_cm3[a] * vn_ref + f.p_cm3[a] * vp_ref);
            const double J = std::abs(p.terminal_current[1]);
            const double e = std::abs(J - J_ref) / J_ref;
            ComparisonRecord c{"current, V/L = " + shortest(drift_fields[k]) + " V/cm", "A/cm^2", drift_fields[k], J,
                               J_ref, 0.0, e, crit, {}};
            c.comparison = judge(r.info, {}, e <= criterion, crit);
            r.comparisons.push_back(std::move(c));
        }
    }
    if (tier == Tier::fast) {
        apply_known_findings(r, known_findings());
        return r;
    }

    // ---- (ii) Energy balance ----
    for (const int N : eb_nodes) {
        for (std::size_t k = 0; k < std::size(eb_fields); ++k) {
            const double L = eb_length[k];
            const device::Device d = bar(N, L);
            std::vector<std::vector<double>> points;  // half decades from 1e2 V/cm: continuation
            for (double E = 1e2; E <= eb_fields[k] * (1.0 + 1e-12); E *= std::sqrt(10.0)) {
                points.push_back({0.0, E * L});
            }
            points.back()[1] = eb_fields[k] * L;
            const std::string label = "V/L = " + shortest(eb_fields[k]) + " V/cm, N = " + std::to_string(N);
            const auto s = solve::sweep_bias(d, points, options(true));
            if (!s || s->stopped || s->points.size() != points.size()) {
                r.run_error = "energy balance, " + label + ": " +
                              (!s ? s.error().message : s->stopped ? s->stopped->message : "short sweep");
                return r;
            }
            const double criterion = 64.0 * eps_m * N * N;
            const std::string crit = "relative <= 64 eps N^2 = " + shortest(criterion);
            const auto& x = d.mesh().points();
            const Layer& ly = layers[k];
            const auto& Tn = s->points.back().fields.electron_temperature_K;
            double worst = 0.0;
            int judged = 0;
            const double d_min = tail * ly.L.length_cm;
            for (std::size_t i = 0; i < Tn.size(); ++i) {
                const double dist = std::min(x[i][0], L - x[i][0]);
                if (dist < d_min) continue;
                ++judged;
                worst = std::max(worst, std::abs((Tn[i] - T) - ly.h.rise) / ly.h.rise);
            }
            // The measured tail near the left contact: ln(dev_i / dev_j) / (x_j - x_i) between the first nodes whose
            // deviation from the centre node's rise falls below 1e-3 and 1e-9 of it.
            const double centre = Tn[Tn.size() / 2] - T;
            std::size_t i3 = Tn.size(), i9 = Tn.size();
            for (std::size_t i = 1; i < Tn.size() / 2; ++i) {
                const double dev = std::abs((Tn[i] - T) - centre) / centre;
                if (i3 == Tn.size() && dev < 1e-3) i3 = i;
                if (i9 == Tn.size() && dev < 1e-9) i9 = i;
            }
            // Printed, not judged (16.42, DT-Q37): the interior's own field against V / L, and the rise against the
            // homogeneous rise at each judged node's local field (its two edges' mean).
            const auto& psi = s->points.back().fields.potential_V;
            double first = -1.0, last = -1.0, psi_first = 0.0, psi_last = 0.0, local = 0.0;
            for (std::size_t i = 1; i + 1 < Tn.size(); ++i) {
                if (std::min(x[i][0], L - x[i][0]) < d_min) continue;
                if (first < 0.0) {
                    first = x[i][0];
                    psi_first = psi[i];
                }
                last = x[i][0];
                psi_last = psi[i];
                const double El = std::abs(psi[i + 1] - psi[i - 1]) / (x[i + 1][0] - x[i - 1][0]);
                const Homogeneous hl = homogeneous(mu_n, sat_n.v_sat_cm_s, sat_n.beta, tau_n, El);
                local = std::max(local, std::abs((Tn[i] - T) - hl.rise) / hl.rise);
            }
            const double E_int = last > first ? std::abs(psi_last - psi_first) / (last - first) : 0.0;
            std::string measured = "not resolved";
            if (i3 < i9 && i9 < Tn.size()) {
                const double a3 = std::abs((Tn[i3] - T) - centre), a9 = std::abs((Tn[i9] - T) - centre);
                measured = shortest(std::log(a3 / a9) / (x[i9][0] - x[i3][0])) + " /cm";
            }
            r.notes.push_back("energy balance, " + label + ", L = " + shortest(L) + " cm: " + std::to_string(judged) +
                              " judged nodes; the tail's measured decay rate " + measured + " against min |Re s| " +
                              shortest(ly.L.slowest_rate_per_cm) + " /cm; interior field " + shortest(E_int) +
                              " V/cm ((E_int - V/L) / (V/L) " + shortest(E_int / eb_fields[k] - 1.0) +
                              "); rise against the homogeneous rise at the local field, max relative " +
                              shortest(local));
            const std::vector<Predicate> pr{
                {"judged nodes (e^{-d/lambda} <= 1e-13) >= 1", static_cast<double>(judged), 1.0, judged >= 1}};
            ComparisonRecord c{"T_n - T_L, " + label, "1", static_cast<double>(N), worst, ly.h.rise, 0.0, worst, crit,
                               {}};
            c.comparison = judge(r.info, pr, worst <= criterion, crit);
            r.comparisons.push_back(std::move(c));
        }
    }
    apply_known_findings(r, known_findings());
    return r;
}

// TVB-32 on this benchmark's devices (16.44): the drift bar (N = 64) swept from zero bias through the four fields; the
// energy-balance bar (N = 128, the coarsest) swept from zero bias by half decades to 1e4 V/cm; each equilibrium.
Tvb32Run tvb13_consistency() {
    Tvb32Run run;
    std::string error;
    for (const bool energy : {false, true}) {
        const int N = energy ? eb_nodes[0] : drift_nodes;
        const double L = energy ? eb_length[1] : drift_length;
        const device::Device d = bar(N, L);
        const std::string name = std::string("TVB-13 ") + (energy ? "energy-balance" : "drift") + " bar N = " +
                                 std::to_string(N);
        std::vector<std::vector<double>> points{{0.0, 0.0}};
        if (energy) {
            for (double E = 1e2; E <= eb_fields[1] * (1.0 + 1e-12); E *= std::sqrt(10.0)) points.push_back({0.0, E * L});
            points.back()[1] = eb_fields[1] * L;
        } else {
            for (const double E : drift_fields) points.push_back({0.0, E * L});
        }
        const solve::BiasOptions o = options(energy);
        const auto s = solve::sweep_bias(d, points, o);
        if (!s || s->stopped || s->points.size() != points.size()) {
            return {std::move(run.entries),
                    name + ": " + (!s ? s.error().message : s->stopped ? s->stopped->message : "short sweep")};
        }
        Tvb32Collector col{run.entries, name, d.contacts().size()};
        const consistency::DeviceData data = device_data(d, o.models);
        for (const auto& p : s->points) col.dc(p, d.mesh(), data, bias_label(p));
        if (!col.equilibrium(d, o.models, error)) return {std::move(run.entries), error};
    }
    return run;
}

}  // namespace NiTCAD::verify::benchmark
