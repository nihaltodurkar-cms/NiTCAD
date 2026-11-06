// TVB-17, the depletion capacitance of an abrupt pn junction at reverse bias (DECISIONS Unit 55 §6, 16.2 (c), 16.4;
// criteria 16.8, frozen; DT-Q4).
//
// Code under test: the drift-diffusion system, solve_small_signal and its operating points.
// Construction (TVB-07's junction):
//   - silicon on [-L, L], L = 1.5 um, N_A = 1e16 cm^-3 for x < 0 (the anode at x_min) and N_D = 2e16 cm^-3 for x > 0;
//   - ohmic contacts; no generation-recombination (SRH, Auger, radiative off); Boltzmann statistics; no band-gap
//     narrowing;
//   - uniform nodes with N even, so x = 0 is the midpoint of the middle edge, a dual-cell face.
// Biases: V = -0.5, -1, -2 V on the anode; at -2 V the p-side depletion edge (0.48 um) leaves about 25 Debye
// lengths to the contact.
// C_j = Im Y_aa / omega at 1 kHz (quasi-static: no recombination, and the dielectric relaxation time is ps).
// Reference (verify/reference/qe_junction):
//   - the QE model's C = -dQ/dV, with Q = eps E_max from TVR-2's first integral at n_i e^(V / 2 V_T), the
//     derivative taken by dual numbers;
//   - cross-checked at run time by complex step (U_ref takes their difference);
//   - inputs q, eps, V_T and n_i are NiTCAD's material data (I1, as TVB-07);
//   - the infinite-domain truncation, the reduced tail t at the contacts (pn_potential at n_i'), is also in U_ref.
// Applicability (DT-Q4), on each operating point: over the nodes where a carrier's density exceeds 1e-6 of the local
// doping, the spread of its quasi-Fermi level must be <= 1e-3 V_T. Otherwise the level is NOT APPLICABLE, with the
// measured drop.
// Criteria (E, frozen 16.8), at each V:
//   - (a) p-hat of the relative error of C_j in [1.8, 2.3], with the asymptotic check (N = 64 ... 512);
//   - (b) finest relative error <= 1e-3;
//   - U_ref <= 1e-4 (a tenth of (b)).
// U_iter (CV-Q8): the operating point's final Newton correction and the AC solve's backward error, the larger.
// The (h, tol) factorial (16.2 (p), DT-Q23 A) at N = 512: tol = 1e-12 and 1e-14, N = 512 and 1024. eps_int enters U
// (U_disc + U_iter + U_ref + |eps_int|), and above 10 % of U it blocks both E passes (DECISIONS 4.1 (c)).
// Reported: the depletion form eps / W (with -2 V_T) against the QE C (qe_junction.hpp: equal up to terms
// exponentially small in the barrier).
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <expected>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "NiTCAD/verify/benchmark/device_data.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/benchmark/tvb32.hpp"
#include "NiTCAD/verify/reference/metrics.hpp"
#include "NiTCAD/verify/reference/pn_junction.hpp"
#include "NiTCAD/verify/reference/qe_junction.hpp"
#include "NiTCAD/verify/study/study.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

constexpr double L = 1.5e-4, NA = 1e16, ND = 2e16, T = 300.0, frequency_Hz = 1e3;
constexpr double tolerance = 1e-3, tol_newton = 1e-12;

ReferenceInfo info() {
    ReferenceInfo r;
    r.id = "TVB-17";
    r.cls = ReferenceClass::E;
    r.independence = Independence::I1;
    r.provenance = Provenance::PD;
    r.citation = "Quasi-equilibrium first integral of the abrupt pn junction (TVR-2's, at n_i e^(V / 2 V_T))";
    r.version = "DECISIONS Unit 55 16.2 (c), 16.3, 16.4, 16.8; DT-Q4";
    r.equations = {1};
    r.form = "C = -dQ/dV, Q = eps E_max = sqrt(2 q eps V_T Phi_n(u_j)) at m = n_i e^(V / 2 V_T)";
    r.assumptions = {"quasi-equilibrium (predicate DT-Q4 measured)", "Boltzmann statistics", "abrupt doping",
                     "infinite domain (contact tail in U_ref)"};
    r.derivation = "QE densities n = n_i e^u, p = n_i e^(v - u) are the equilibrium ones at m after u -> u - v/2";
    r.verification = "zero bias = pn_junction; shift identity; dual vs complex step vs central difference";
    r.second_method = "complex step on Q(v) with complex Newton for u_j, at run time";
    return r;
}

std::vector<double> nodes(int N) {
    std::vector<double> x;
    for (int k = 0; k < N; ++k) x.push_back(L * (2.0 * k - (N - 1)) / (N - 1));
    return x;
}

device::Device diode(int N) {
    mesh::Mesh m = *mesh::make_tensor_grid(nodes(N));
    const std::size_t n = m.node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        if (m.points()[i][0] < 0.0) {
            acceptors[i] = NA;
        } else {
            donors[i] = ND;
        }
    }
    auto anode = m.find_boundary("x_min")->nodes;
    auto cathode = m.find_boundary("x_max")->nodes;
    return *device::Device::create({.mesh = std::move(m),
                                    .temperature_K = T,
                                    .regions = {{"silicon", physics::silicon()}},
                                    .node_region = std::vector<device::RegionId>(n, 0),
                                    .donors = std::move(donors),
                                    .acceptors = std::move(acceptors),
                                    .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                                                 {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
}

solve::BiasOptions options() {
    solve::BiasOptions o;
    o.models.srh = o.models.auger = o.models.radiative = false;
    o.models.bgn = false;
    o.models.fermi_dirac = false;
    return o;
}

struct Level {
    int N = 0;
    double h = 0.0, C = 0.0, error = 0.0, u_iter = 0.0, drop_n = 0.0, drop_p = 0.0;
};

}  // namespace

BenchmarkResult run_tvb17(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-17";
    r.title = "Depletion capacitance at reverse bias";
    r.info = info();
    if (tier == Tier::fast) {
        r.notes.push_back("full tier only (16.7)");
        return r;
    }
    const physics::Semiconductor si = physics::silicon();
    const double VT = base::thermal_voltage(T);
    const reference::PnJunctionInput in{.q_C = base::q_C,
                                        .eps_F_per_cm = si.parameters().eps_r * base::eps0_F_per_cm,
                                        .thermal_voltage_V = VT,
                                        .ni_cm3 = physics::intrinsic_density(si, T),
                                        .acceptors_cm3 = NA,
                                        .donors_cm3 = ND};
    for (const double V : {-0.5, -1.0, -2.0}) {
        const std::string at = "V = " + shortest(V) + " V";
        const auto ref = reference::qe_junction(in, V);
        const auto cs = reference::qe_capacitance_complex_step(in, V);
        if (!ref || !cs) {
            r.run_error = at + ": reference: " + (!ref ? ref.error() : cs.error());
            return r;
        }
        const double C_ref = ref->capacitance_F_per_cm2.to_double();
        // The contact tail of the infinite-domain solution, at n_i' (the shifted equilibrium).
        reference::PnJunctionInput shifted = in;
        shifted.ni_cm3 = ref->m_cm3.to_double();
        const auto eq = reference::pn_junction(shifted);
        if (!eq) {
            r.run_error = at + ": reference tail: " + eq.error();
            return r;
        }
        const std::vector<double> ends{-L, L};
        const auto tails = reference::pn_potential(*eq, ends);
        if (!tails) {
            r.run_error = at + ": reference tail: " + tails.error();
            return r;
        }
        const double tail = std::max((*tails)[0].tail.to_double(), (*tails)[1].tail.to_double());
        const double u_ref = std::abs(*cs - C_ref) / C_ref + tail;
        const double dep = (ref->depletion_form_F_per_cm2.to_double() - C_ref) / C_ref;
        r.notes.push_back(at + ": C_QE " + shortest(C_ref) + " F/cm^2, complex step differs by " +
                          shortest(std::abs(*cs - C_ref) / C_ref) + ", contact tail " + shortest(tail) +
                          "; eps / W (-2 V_T) differs by " + shortest(dep));
        // One level at N nodes and Newton tolerance tol (the factorial reruns it at tol / 100 and 2N).
        const auto run_level = [&](int N, double tol) -> std::expected<Level, std::string> {
            const device::Device d = diode(N);
            const std::vector<std::vector<double>> points{{V, 0.0}};
            solve::SmallSignalOptions o{.steady = options(), .frequencies_Hz = {frequency_Hz}};
            o.steady.newton.tol_update = tol;
            const auto s = solve::solve_small_signal(d, points, o);
            if (!s || s->stopped) {
                return std::unexpected(at + ", N = " + std::to_string(N) + ": " +
                                       (s ? s->stopped->message : s.error().message));
            }
            const auto& dc = s->points[0].dc;
            Level lv;
            lv.N = N;
            lv.h = 2.0 * L / (N - 1);
            lv.C = s->capacitance(0, 0, 0, 0);
            lv.error = std::abs(lv.C - C_ref) / C_ref;
            const auto& it = dc.convergence.iterations;
            lv.u_iter = std::max(it.empty() ? 0.0 : it.back().update, s->points[0].backward_error[0]);
            // DT-Q4: the quasi-Fermi spread where the carrier's charge matters.
            double nlo = INFINITY, nhi = -INFINITY, plo = INFINITY, phi = -INFINITY;
            for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
                const double doping = d.mesh().points()[i][0] < 0.0 ? NA : ND;
                if (dc.fields.n_cm3[i] > 1e-6 * doping) {
                    nlo = std::min(nlo, dc.bands.electron_fermi_eV[i]);
                    nhi = std::max(nhi, dc.bands.electron_fermi_eV[i]);
                }
                if (dc.fields.p_cm3[i] > 1e-6 * doping) {
                    plo = std::min(plo, dc.bands.hole_fermi_eV[i]);
                    phi = std::max(phi, dc.bands.hole_fermi_eV[i]);
                }
            }
            lv.drop_n = nhi - nlo;  // eV, i.e. V for a unit charge
            lv.drop_p = phi - plo;
            return lv;
        };
        std::vector<Level> levels;  // finest first
        bool applicable = true;
        double worst_drop = 0.0;
        for (const int N : {512, 256, 128, 64}) {
            const auto lv = run_level(N, tol_newton);
            if (!lv) {
                r.run_error = lv.error();
                return r;
            }
            worst_drop = std::max({worst_drop, lv->drop_n, lv->drop_p});
            applicable = applicable && lv->drop_n <= 1e-3 * VT && lv->drop_p <= 1e-3 * VT;
            r.notes.push_back(at + ", N = " + std::to_string(N) + ": C " + shortest(lv->C) + ", relative error " +
                              shortest(lv->error) + ", U_iter " + shortest(lv->u_iter) + ", quasi-Fermi drops " +
                              shortest(lv->drop_n) + " (n) and " + shortest(lv->drop_p) + " (p) V");
            levels.push_back(*lv);
        }
        std::vector<double> hs, es, us;
        for (const Level& l : levels) {
            hs.push_back(l.h);
            es.push_back(l.error);
            us.push_back(l.u_iter);
        }
        const study::OrderStudy st = study::order_study(hs, es, us);
        const study::OrderVerdict v = study::judge_order(st, {.p_th = 2.0});
        std::string orders;
        for (const double p : st.orders) orders += " " + shortest(p);
        r.notes.push_back(at + ": orders (finest pair first):" + orders + "; spread " + shortest(st.order_spread));
        // The (h, tol) factorial at N = 512 (16.2 (p), DT-Q23 A); h / 2 is N = 1024 (N even keeps x = 0 a face).
        // U at (h, tol): U_disc (GCI of the finest pair with p-hat, F_s = 1.25), U_iter, U_ref.
        const double u_disc = reference::gci_safety_factor_three_grids * std::abs(levels[1].error - levels[0].error) /
                              (std::pow(levels[1].h / levels[0].h, st.observed) - 1.0);
        double fac_e[3] = {0.0, 0.0, 0.0};
        {
            int k = 0;
            for (const auto& [N, tol] : {std::pair{512, 0.01 * tol_newton}, std::pair{1024, tol_newton},
                                         std::pair{1024, 0.01 * tol_newton}}) {
                const auto lv = run_level(N, tol);
                if (!lv) {
                    r.run_error = lv.error() + " (factorial)";
                    return r;
                }
                fac_e[k++] = lv->error;
            }
        }
        const std::vector<double> parts{u_disc, levels.front().u_iter, u_ref};
        const Factorial fac = factorial(levels.front().error, fac_e[0], fac_e[1], fac_e[2], parts);
        r.notes.push_back(at + ": " + fac.note + " (U_disc " + shortest(u_disc) + ")");
        const std::optional<std::string> excluded = [&]() -> std::optional<std::string> {
            if (!applicable) return std::nullopt;
            if (u_ref > 0.1 * tolerance) return "U_ref " + shortest(u_ref) + " above a tenth of the tolerance";
            return std::nullopt;
        }();
        const std::string why_na = "quasi-Fermi drop " + shortest(worst_drop) + " V above 1e-3 V_T (DT-Q4)";
        const std::string band = "p-hat in [1.8, 2.3] (p_th 2)", accuracy = "relative <= 1e-3";
        ComparisonRecord co{"C_j order, " + at, "1", V, st.observed, 2.0, u_ref, st.observed, band, {}};
        ComparisonRecord ca{"C_j finest (N = 512), " + at, "F/cm^2", V, levels.front().C, C_ref, u_ref,
                            levels.front().error, accuracy, {}};
        if (!applicable) {
            co.comparison = Comparison{Outcome::not_applicable, why_na, {}};
            ca.comparison = Comparison{Outcome::not_applicable, why_na, {}};
        } else if (excluded) {
            co.comparison = Comparison{Outcome::not_evaluated, *excluded, {}};
            ca.comparison = Comparison{Outcome::not_evaluated, *excluded, {}};
        } else {
            co.comparison =
                v.judged ? judge(r.info, {}, v.pass, band) : Comparison{Outcome::not_evaluated, v.reason, {}};
            ca.comparison = judge(r.info, {}, levels.front().error <= tolerance, accuracy);
            block_if_interacting(co, fac);
            block_if_interacting(ca, fac);
        }
        r.comparisons.push_back(std::move(co));
        r.comparisons.push_back(std::move(ca));
    }
    apply_known_findings(r, known_findings());
    return r;
}

// TVB-32 on this benchmark's device (16.44): N = 64; the small-signal run at zero bias and the three reverse biases,
// at its frequency; the equilibrium.
Tvb32Run tvb17_consistency() {
    Tvb32Run run;
    const device::Device d = diode(64);
    solve::SmallSignalOptions o{.steady = options(), .frequencies_Hz = {frequency_Hz}};
    o.steady.newton.tol_update = tol_newton;
    const auto s = solve::solve_small_signal(
        d, std::vector<std::vector<double>>{{0.0, 0.0}, {-0.5, 0.0}, {-1.0, 0.0}, {-2.0, 0.0}}, o);
    if (!s || s->stopped) return {{}, "TVB-17 diode N = 64: " + (s ? s->stopped->message : s.error().message)};
    Tvb32Collector col{run.entries, "TVB-17 diode N = 64", d.contacts().size()};
    col.ac(*s, d.mesh(), device_data(d, o.steady.models), "");
    std::string error;
    if (!col.equilibrium(d, o.steady.models, error)) return {std::move(run.entries), error};
    return run;
}

}  // namespace NiTCAD::verify::benchmark
