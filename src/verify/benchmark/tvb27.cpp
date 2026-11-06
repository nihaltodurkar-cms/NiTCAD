// TVB-27 (lattice heat) and TVB-28 (thermoelectric identities) (DECISIONS Unit 55 §6, 16.2 (k), 16.4; criteria
// 16.8, frozen; construction 16.43).
//
// Code under test: solve_bias with the electrothermal model (Unit 23) on uniform 1D n-bars (N_D = 1e17 cm^-3),
// Boltzmann, no band-gap narrowing, no generation-recombination, Caughey-Thomas mobilities with doping.
// TVB-27 (E): T(x) in L-infinity, relative to max T_ref; p-hat in [1.8, 2.3] on N = 64 ... 512; the finest relative
// <= 1e-4; the (h, tol) factorial at N = 512 (16.2 (p)).
//   (i) The Kirchhoff slab: 10 um, one ohmic contact (no current flows), isothermal 300 K and 400 K ends. Its T
//       is the Kirchhoff closed form (verify/reference/thermal_slab (i)).
//   (ii) The self-heated resistor: 10 um, ohmic and isothermal 300 K at both ends, 4 V. Its T is the reference ODE
//       (thermal_slab (ii)), solved in dd with its finite-difference second method.
// TVB-28, on the legacy silicon set (Nc = Nc300 (T / 300)^(3/2), so S and Pi are closed forms; 16.43):
//   (iii) S (E): the open-circuit rise of the electron quasi-Fermi level across a 100 um bar with one ohmic contact
//       between 300 K and 400 K (thermal_slab (iii)); as TVB-27, p-hat and the finest relative <= 1e-4.
//   (iv) Pi (E): an isothermal 300 K bar (10 um, two ohmic contacts) at +-1 mV: the odd part of the right sink's heat
//       over the current, against -Pi_n (thermal_slab (iv)); as (iii).
//   (v) The energy residual (D): consistency check C5 (sum V I against the sinks' heat, its frozen bound, 11.6) on
//       the finest state of (ii), (iii) and both of (iv); C5's accepted limitation (the heat rows' Newton term,
//       11.13) is printed.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <expected>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/verify/benchmark/device_data.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/consistency/checks.hpp"
#include "NiTCAD/verify/reference/metrics.hpp"
#include "NiTCAD/verify/reference/thermal_slab.hpp"
#include "NiTCAD/verify/study/study.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

constexpr double ND = 1e17, tol_newton = 1e-10, tolerance = 1e-4;  // 1e-10: 16.43 (the factorial's 1e-14 failed)
constexpr int ladder[4] = {512, 256, 128, 64};  // finest first
constexpr int h_half = 1024;

ReferenceInfo info(const char* id, bool thermoelectric) {
    ReferenceInfo r;
    r.id = id;
    r.cls = ReferenceClass::E;
    r.independence = Independence::I1;
    r.provenance = Provenance::PD;
    r.citation = thermoelectric ? "Unit 23 T2, T4 (Wachutka's model as NiTCAD states it), derived to S and Pi (PD)"
                                : "the Kirchhoff transform; Unit 23 T2, T4 reduced to the uniform resistor (PD)";
    r.version = "DECISIONS Unit 55 16.2 (k), 16.8, 16.43";
    r.equations = {4};
    r.form = thermoelectric ? "E_Fn(T1) - E_Fn(T0) = -k int ((5/2 + r) - ln(N / Nc)) dT; Pi_n = T P_n"
                            : "u(T) linear; (kappa T')' = -J^2 / (q mu N) - (3/2)(k/q) J T'";
    r.assumptions = {"n = N_D throughout (the space charge neglected)", "holes neglected", "Boltzmann statistics"};
    r.derivation = "verify/reference/thermal_slab.hpp";
    r.verification = "thermal_slab_test: Kirchhoff round trip, the constant-coefficient resistor, Seebeck by GK";
    r.second_method = thermoelectric ? "Gauss-Kronrod on the Seebeck integrand"
                                     : "Kirchhoff-form finite differences, 2^14 and 2^15 cells, Richardson";
    return r;
}

struct Bar {
    double length;
    physics::SemiconductorParameters material;
    bool two_contacts;
    double T_left, T_right;
};

device::Device build(const Bar& b, int N) {
    std::vector<double> x(static_cast<std::size_t>(N));
    for (int k = 0; k < N; ++k) x[static_cast<std::size_t>(k)] = b.length * k / (N - 1);
    x.back() = b.length;
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    auto a = m.find_boundary("x_min")->nodes;
    auto c = m.find_boundary("x_max")->nodes;
    std::vector<device::Contact> contacts{{"left", device::ContactKind::ohmic, std::move(a)}};
    if (b.two_contacts) contacts.push_back({"right", device::ContactKind::ohmic, std::move(c)});
    return *device::Device::create(
        {.mesh = std::move(m),
         .temperature_K = 300.0,
         .regions = {{"silicon", *physics::Semiconductor::create(b.material)}},
         .node_region = std::vector<device::RegionId>(n, 0),
         .donors = std::vector<double>(n, ND),
         .acceptors = std::vector<double>(n, 0.0),
         .contacts = std::move(contacts),
         .thermal_contacts = {{"cold", "x_min", device::ThermalContactKind::isothermal, b.T_left, 0.0},
                              {"hot", "x_max", device::ThermalContactKind::isothermal, b.T_right, 0.0}}});
}

solve::BiasOptions options(double tol) {
    solve::BiasOptions o;
    o.models.electrothermal = true;
    o.models.electron_mobility = assemble::ElectronMobility::caughey_thomas;
    o.models.hole_mobility = assemble::HoleMobility::caughey_thomas;
    o.models.srh = o.models.auger = o.models.radiative = false;
    o.models.bgn = false;
    o.models.fermi_dirac = false;
    o.newton.tol_update = tol;
    return o;
}

// One E quantity on the ladder and the factorial's states: e(N, tol) signed relative errors.
struct Errors {
    std::vector<double> ladder;          // finest first, at tol
    double fine_tol = 0.0, half = 0.0, half_tol = 0.0;  // (512, tol/100), (1024, tol), (1024, tol/100)
};

// Judges an E quantity: the order (p_th 2) and the finest accuracy, with the factorial's verdict.
void judge_e(BenchmarkResult& r, const std::string& name, const Errors& e, const std::vector<Predicate>& pr) {
    std::vector<double> hs, es, us;
    const double u_iter = std::abs(e.ladder[0] - e.fine_tol);
    for (int k = 0; k < 4; ++k) {
        hs.push_back(1.0 / (ladder[k] - 1));
        es.push_back(std::abs(e.ladder[static_cast<std::size_t>(k)]));
        us.push_back(u_iter);
    }
    const study::OrderStudy st = study::order_study(hs, es, us);
    const study::OrderVerdict v = study::judge_order(st, {.p_th = 2.0});
    double gci = 0.0;
    if (!st.orders.empty() && st.observed > 0.0) {
        gci = reference::gci_safety_factor_three_grids * std::abs(st.error[1] - st.error[0]) /
              (std::pow(st.h[1] / st.h[0], st.observed) - 1.0);
    }
    const std::vector<double> parts{gci, u_iter};
    const Factorial f = factorial(e.ladder[0], e.fine_tol, e.half, e.half_tol, parts);
    std::string orders;
    for (const double p : st.orders) orders += " " + shortest(p);
    r.notes.push_back(name + ": errors");
    for (const double x : e.ladder) r.notes.back() += " " + shortest(x);
    r.notes.back() += "; orders" + orders + " (spread " + shortest(st.order_spread) + "); " + f.note;
    const std::string band = "p-hat in [1.8, 2.3] (p_th 2)";
    const bool holds = std::all_of(pr.begin(), pr.end(), [](const Predicate& p) { return p.holds; });
    ComparisonRecord co{name + " order in h", "1", 0.0, st.observed, 2.0, 0.0, st.observed, band, {}};
    co.comparison = !v.judged && holds ? Comparison{Outcome::not_evaluated, v.reason, {}}
                                       : judge(r.info, pr, v.pass, band);
    ComparisonRecord ca{name + " finest accuracy", "1", 0.0, std::abs(e.ladder[0]), 0.0, 0.0, std::abs(e.ladder[0]),
                        "relative <= 1e-4", {}};
    ca.comparison = judge(r.info, pr, std::abs(e.ladder[0]) <= tolerance, ca.criterion);
    for (ComparisonRecord* c : {&co, &ca}) {
        if (c->comparison.outcome != Outcome::not_applicable) block_if_interacting(*c, f);
        r.comparisons.push_back(std::move(*c));
    }
}

// The four states of a ladder plus the factorial's three, by a run function of (N, tol).
template <class F>
std::expected<Errors, std::string> errors_of(F&& run) {
    Errors e;
    for (const int N : ladder) {
        auto v = run(N, tol_newton);
        if (!v) return std::unexpected(v.error());
        e.ladder.push_back(*v);
    }
    auto a = run(ladder[0], 0.01 * tol_newton), b = run(h_half, tol_newton), c = run(h_half, 0.01 * tol_newton);
    for (auto* x : {&a, &b, &c}) {
        if (!*x) return std::unexpected(x->error() + " (factorial)");
    }
    e.fine_tol = *a;
    e.half = *b;
    e.half_tol = *c;
    return e;
}

double linf_relative(const std::vector<double>& sim, const std::vector<reference::DD>& ref) {
    double worst = 0.0, scale = 0.0;
    for (std::size_t i = 0; i < sim.size(); ++i) {
        worst = std::max(worst, std::abs(sim[i] - ref[i].to_double()));
        scale = std::max(scale, std::abs(ref[i].to_double()));
    }
    return worst / scale;
}

std::vector<double> nodes_x(const device::Device& d) {
    std::vector<double> x;
    for (const auto& p : d.mesh().points()) x.push_back(p[0]);
    return x;
}

}  // namespace

BenchmarkResult run_tvb27(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-27";
    r.title = "Lattice heat: the Kirchhoff slab and the self-heated resistor";
    r.info = info("TVB-27", false);
    if (tier == Tier::fast) {
        r.notes.push_back("full tier only (16.7)");
        return r;
    }
    const physics::SemiconductorParameters& sp = physics::silicon_parameters;
    const reference::Conductivity kappa{sp.thermal.conductivity_W_cmK, sp.thermal.conductivity_exponent};

    // (i) The slab.
    const Bar slab{1e-3, sp, false, 300.0, 400.0};
    const auto slab_errors = errors_of([&](int N, double tol) -> std::expected<double, std::string> {
        const device::Device d = build(slab, N);
        const auto p = solve::solve_bias(d, std::vector<double>{0.0}, options(tol));
        if (!p) return std::unexpected("slab, N = " + std::to_string(N) + ": " + p.error().message);
        std::vector<double> s;
        for (const double x : nodes_x(d)) s.push_back(x / slab.length);
        const auto T = reference::slab_temperature(kappa, 300.0, 400.0, s);
        if (!T) return std::unexpected("reference: " + T.error());
        return linf_relative(p->fields.temperature_K, *T);
    });
    if (!slab_errors) {
        r.run_error = slab_errors.error();
        return r;
    }
    judge_e(r, "T (Kirchhoff slab)", *slab_errors, {});

    // (ii) The self-heated resistor.
    const Bar res{1e-3, sp, true, 300.0, 300.0};
    const physics::CaugheyThomasParameters& ct = sp.electron_mobility;
    reference::ResistorInput in{.q_C = base::q_C,
                                .k_eV_per_K = base::k_B_eV_per_K,
                                .length_cm = res.length,
                                .density_cm3 = ND,
                                .voltage_V = 4.0,
                                .ambient_K = 300.0,
                                .kappa = kappa,
                                .mu_min = ct.mu_min,
                                .mu_max = ct.mu_max,
                                .N_ref = ct.N_ref,
                                .alpha = ct.alpha,
                                .gamma = ct.T_exponent};
    // The reference once, at every node of every level (the uniform grids are not nested), with its second method.
    std::vector<double> all_x;
    std::vector<std::size_t> offset;
    for (const int N : {ladder[0], ladder[1], ladder[2], ladder[3], h_half}) {
        offset.push_back(all_x.size());
        for (int k = 0; k < N; ++k) all_x.push_back(k == N - 1 ? res.length : res.length * k / (N - 1));
    }
    const auto Tref = reference::self_heated_resistor(in, all_x);
    if (!Tref) {
        r.run_error = "reference: " + Tref.error();
        return r;
    }
    if (Tref->second_method > 1e-8 || Tref->steps_difference > 1e-10) {
        r.run_error = "reference: its second method disagrees (" + shortest(Tref->second_method) + " K; steps " +
                      shortest(Tref->steps_difference) + " K)";
        return r;
    }
    double peak = 0.0;
    for (const auto& t : Tref->temperature) peak = std::max(peak, t.to_double());
    const std::string ref_note = "self-heated resistor reference: J " + shortest(Tref->current_A_cm2.to_double()) +
                                 " A/cm^2, peak " + shortest(peak) + " K, " + std::to_string(Tref->iterations) +
                                 " secant iterations; finite differences differ by " + shortest(Tref->second_method) +
                                 " K, the halved steps by " + shortest(Tref->steps_difference) + " K";
    const auto level = [&](int N) {
        const int levels[5] = {ladder[0], ladder[1], ladder[2], ladder[3], h_half};
        const std::size_t i = static_cast<std::size_t>(std::find(levels, levels + 5, N) - levels);
        return std::vector<reference::DD>(Tref->temperature.begin() + static_cast<std::ptrdiff_t>(offset[i]),
                                          Tref->temperature.begin() + static_cast<std::ptrdiff_t>(offset[i] + N));
    };
    const auto res_errors = errors_of([&](int N, double tol) -> std::expected<double, std::string> {
        const device::Device d = build(res, N);
        const auto p = solve::solve_bias(d, std::vector<double>{0.0, 4.0}, options(tol));
        if (!p) return std::unexpected("resistor, N = " + std::to_string(N) + ": " + p.error().message);
        return linf_relative(p->fields.temperature_K, level(N));
    });
    if (!res_errors) {
        r.run_error = res_errors.error();
        return r;
    }
    r.notes.push_back(ref_note);
    judge_e(r, "T (self-heated resistor)", *res_errors, {});
    apply_known_findings(r, known_findings());
    return r;
}

BenchmarkResult run_tvb28(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-28";
    r.title = "Thermoelectric identities: Seebeck, Peltier and the energy balance";
    r.info = info("TVB-28", true);
    if (tier == Tier::fast) {
        r.notes.push_back("full tier only (16.7)");
        return r;
    }
    const physics::SemiconductorParameters& sp = physics::legacy_silicon_parameters;
    const double rexp = sp.thermal.thermopower_exponent_n;
    std::vector<std::pair<std::string, consistency::CheckResult>> c5;

    // (iii) Seebeck.
    const auto seebeck = reference::seebeck_rise({.k_eV_per_K = base::k_B_eV_per_K,
                                                  .density_cm3 = ND,
                                                  .Nc300_cm3 = sp.Nc300,
                                                  .thermopower_exponent = rexp,
                                                  .T0 = 300.0,
                                                  .T1 = 400.0});
    if (!seebeck || seebeck->second_method > 1e-25) {
        r.run_error = "reference: " + (!seebeck ? seebeck.error() : "the Seebeck second method disagrees");
        return r;
    }
    const double S_ref = seebeck->rise_eV.to_double();
    const Bar sb{1e-2, sp, false, 300.0, 400.0};
    const auto s_errors = errors_of([&](int N, double tol) -> std::expected<double, std::string> {
        const device::Device d = build(sb, N);
        const auto p = solve::solve_bias(d, std::vector<double>{0.0}, options(tol));
        if (!p) return std::unexpected("Seebeck, N = " + std::to_string(N) + ": " + p.error().message);
        if (N == ladder[0] && tol == tol_newton) {
            c5.emplace_back("Seebeck bar", consistency::check_c5(*p, d.mesh(), device_data(d, options(tol).models)));
        }
        const auto& F = p->bands.electron_fermi_eV;
        return (F.back() - F.front()) / S_ref - 1.0;
    });
    if (!s_errors) {
        r.run_error = s_errors.error();
        return r;
    }
    r.notes.push_back("Seebeck reference: rise " + shortest(S_ref) + " eV (300 K to 400 K, legacy Nc300 " +
                      shortest(sp.Nc300) + "); Gauss-Kronrod differs by " + shortest(seebeck->second_method));
    judge_e(r, "S (open-circuit rise)", *s_errors, {});

    // (iv) Peltier.
    const double Pi_ref = reference::peltier_coefficient(base::k_B_eV_per_K, ND, sp.Nc300, rexp, 300.0).to_double();
    const Bar pb{1e-3, sp, true, 300.0, 300.0};
    const auto p_errors = errors_of([&](int N, double tol) -> std::expected<double, std::string> {
        const device::Device d = build(pb, N);
        const auto plus = solve::solve_bias(d, std::vector<double>{1e-3, 0.0}, options(tol));
        const auto minus = solve::solve_bias(d, std::vector<double>{-1e-3, 0.0}, options(tol));
        if (!plus || !minus) {
            return std::unexpected("Peltier, N = " + std::to_string(N) + ": " +
                                   (!plus ? plus.error().message : minus.error().message));
        }
        if (N == ladder[0] && tol == tol_newton) {
            c5.emplace_back("Peltier bar, +1 mV",
                            consistency::check_c5(*plus, d.mesh(), device_data(d, options(tol).models)));
            c5.emplace_back("Peltier bar, -1 mV",
                            consistency::check_c5(*minus, d.mesh(), device_data(d, options(tol).models)));
        }
        const double dQ = plus->thermal_contact_heat[1] - minus->thermal_contact_heat[1];
        const double dI = plus->terminal_current[1] - minus->terminal_current[1];
        return (dQ / dI) / -Pi_ref - 1.0;
    });
    if (!p_errors) {
        r.run_error = p_errors.error();
        return r;
    }
    r.notes.push_back("Peltier reference: Pi_n(300 K) " + shortest(Pi_ref) + " V");
    judge_e(r, "Pi (contact heat per current)", *p_errors, {});

    // (v) The energy residual on the self-heated resistor's finest state (TVB-27 (ii)'s device on this set).
    {
        const Bar res{1e-3, sp, true, 300.0, 300.0};
        const device::Device d = build(res, ladder[0]);
        const auto p = solve::solve_bias(d, std::vector<double>{0.0, 4.0}, options(tol_newton));
        if (!p) {
            r.run_error = "resistor (C5): " + p.error().message;
            return r;
        }
        c5.emplace(c5.begin(), "self-heated resistor, 4 V", consistency::check_c5(*p, d.mesh(),
                                                                                  device_data(d, options(tol_newton).models)));
    }
    for (auto& [label, c] : c5) {
        std::string parts;
        for (const auto& part : c.parts) {
            parts += " " + part.name + " " + shortest(part.residual) + " (<= " + shortest(part.tolerance) + ")";
        }
        r.notes.push_back("C5, " + label + ":" + parts +
                          "; accepted limitation (11.13): the heat rows' Newton term is not in the bound");
        ComparisonRecord rec{"energy residual (C5), " + label, "W/cm^2", 0.0, 0.0, 0.0, 0.0, 0.0,
                             "C5's frozen bound (11.6)", {}};
        if (!c.parts.empty()) {
            rec.q_sim = c.parts[0].residual;
            rec.metric = c.parts[0].tolerance > 0.0 ? c.parts[0].residual / c.parts[0].tolerance : 0.0;
        }
        rec.comparison = c.comparison;
        r.comparisons.push_back(std::move(rec));
    }
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
