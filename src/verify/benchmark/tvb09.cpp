// TVB-09, the MOS capacitor on a uniform substrate (DECISIONS Unit 55 6, 14.4; criteria 14.7).
//
// Code under test: the quasi-static solve (Equations::equilibrium_poisson: Poisson with the carriers in thermal
// equilibrium, Boltzmann), the lumped gate (its Robin row, the work-function offset and the fixed oxide charge,
// assemble/gate.cpp), the gate charge it reports, and sweep_bias. Construction: p-silicon, N_A = 1e17 cm^-3, depth
// d = 0.5 um, a metal gate (phi_m = 4.7 eV) on a lumped 10 nm SiO2 (eps_r 3.9) with Q_f = 5e10 cm^-2 at x = 0, an
// ohmic substrate contact at x = d; no band-gap narrowing.
// Mesh: the inversion and accumulation layers are under 1 nm thick, so a uniform 64 ... 512-node mesh of 0.5 um
// cannot reach the asymptotic range there. The nodes are a smooth exponential grading toward the gate,
//     x_k = d (e^{a xi_k} - 1) / (e^a - 1),   xi_k = k / (N - 1),   a = 6.2,
// first spacing 0.1 nm at N = 64; refining xi uniformly keeps the box method second order (a smooth map), and h is
// the xi spacing 1 / (N - 1).
// Reference (verify/reference/mos_capacitor.hpp): Kingston-Neustadter's exact Q_s(psi_s), u_s(V_G) by Brent in
// double-double, C_LF = C_s C_ox / (C_s + C_ox) analytic and by complex step (second method, checked at run time);
// V_FB from the material data. NiTCAD's side:
// - Q_s = -(Q_G + q Q_f), Q_G the reported gate charge;
// - C_LF = dQ_G / dV_G of NiTCAD's own quasi-static charge, by Richardson-extrapolated central differences
//   (steps 4 mV and 2 mV) at the finest level; U_fd = |D(h/2) - D(h)| / 3 must be at most a tenth of the criterion;
// - V_FB = V_G - psi_s - Q_G / C_ox + psi_b - q Q_f / C_ox (the gate row's offset made explicit; psi_b the
//   substrate contact's potential), at V_G = V_FB(reference), at every level.
// Criteria (frozen 14.7): D: V_FB relative <= 64 eps. E (a): Q_s(V_G) p-hat in [1.8, 2.3] at each V_G. E (b): on the
// finest level Q_s and C_LF relative <= 1e-3. On the V_G grid V_FB + {-1, -0.5, -0.2, -0.05, 0, 0.05, 0.2, 0.5, 1,
// 1.5, 2} V (accumulation to strong inversion; threshold about V_FB + 1.3 V), a point with |Q_s| < 1e-3 of the
// grid's range is not_applicable (the predicate named): there a relative error has no scale.
// Ladder: N = 64, 128, 256, 512 at every bias (full); fast: V_FB + 1.5 V (strong inversion) only (14.7).
// Reference uncertainty: the substrate's finite depth, the infinite-depth profile's tail V_T t at x = d relative to
// the band bending, at most 1e-6 (1e-3 of the criterion, TB-Q4), else not_evaluated.
// Reference independence I1: q, eps_si, V_T, n_i, chi, E_g, N_c and N_v are NiTCAD's material data (input).
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/first_integral.hpp"
#include "NiTCAD/verify/reference/mos_capacitor.hpp"
#include "NiTCAD/verify/study/study.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

constexpr double depth_cm = 5e-5;
constexpr double acceptors_cm3 = 1e17;
constexpr double t_ox_cm = 1e-6, eps_ox_r = 3.9, phi_m_eV = 4.7, fixed_charge_cm2 = 5e10;
constexpr double grading = 6.2;
constexpr double eps = 0x1p-52;
constexpr double tolerance = 1e-3;
constexpr double fd_step_V = 4e-3;
const char* const applicability = "|Q_s| / range(Q_s) >= 1e-3";

ReferenceInfo info() {
    ReferenceInfo r;
    r.id = "TVB-09";
    r.cls = ReferenceClass::E;
    r.independence = Independence::I1;
    r.provenance = Provenance::PA;
    r.citation = "Kingston and Neustadter, J. Appl. Phys. 26, 718 (1955); Sze and Ng, Physics of Semiconductor "
                 "Devices, 3rd ed. (2007), 4.2 (the F-function and the flat-band voltage)";
    r.version = "DECISIONS Unit 55 14.4, 14.7";
    r.equations = {1, 4};
    r.form = "Q_s = -sgn(t) sqrt(2 q eps V_T Phi_b(u_s)); V_G = V_FB + V_T t - Q_s / C_ox";
    r.assumptions = {"thermal equilibrium (quasi-static)", "Boltzmann statistics", "uniform substrate",
                     "lumped oxide with fixed charge at the interface", "infinite-depth reference (truncation measured)"};
    r.derivation = "verify/reference/mos_capacitor.hpp";
    r.verification = "D: V_FB to rounding; E: Q_s orders and finest-level Q_s and C_LF";
    r.second_method = "C_LF by complex step against the analytic derivative, at run time";
    return r;
}

std::vector<double> graded(int N) {
    std::vector<double> x;
    const double scale = depth_cm / std::expm1(grading);
    for (int k = 0; k < N; ++k) x.push_back(k == N - 1 ? depth_cm : scale * std::expm1(grading * k / (N - 1)));
    return x;
}

solve::BiasOptions options() {
    solve::BiasOptions o;
    o.equations = solve::Equations::equilibrium_poisson;
    o.models.bgn = false;
    o.models.fermi_dirac = false;
    o.newton.tol_update = 1e-13;
    return o;
}

struct Solved {
    double charge;  // Q_s
    double psi_s, psi_b, gate_charge, u_iter;
    results::NodeFields fields;
};

}  // namespace

BenchmarkResult run_tvb09(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-09";
    r.title = "MOS capacitor on a uniform substrate";
    r.info = info();
    const physics::Semiconductor si = physics::silicon();
    const double T = 300.0, VT = base::thermal_voltage(T), q = base::q_C;
    const double cox = eps_ox_r * base::eps0_F_per_cm / t_ox_cm;
    const double ni = physics::intrinsic_density(si, T);
    const auto vfb = reference::flat_band_voltage({.work_function_eV = phi_m_eV,
                                                   .electron_affinity_eV = si.parameters().electron_affinity_eV,
                                                   .band_gap_eV = physics::band_gap_eV(si, T),
                                                   .conduction_dos_cm3 = physics::conduction_band_dos(si, T),
                                                   .valence_dos_cm3 = physics::valence_band_dos(si, T),
                                                   .thermal_voltage_V = VT,
                                                   .ni_cm3 = ni,
                                                   .net_doping_cm3 = -acceptors_cm3,
                                                   .fixed_charge_cm2 = fixed_charge_cm2,
                                                   .q_C = q,
                                                   .oxide_capacitance_F_per_cm2 = cox});
    if (!vfb) {
        r.run_error = "reference: " + vfb.error();
        return r;
    }
    const double VFB = vfb->to_double();
    const reference::MosInput in{.q_C = q,
                                 .eps_F_per_cm = si.parameters().eps_r * base::eps0_F_per_cm,
                                 .thermal_voltage_V = VT,
                                 .ni_cm3 = ni,
                                 .net_doping_cm3 = -acceptors_cm3,
                                 .oxide_capacitance_F_per_cm2 = cox,
                                 .flat_band_V = VFB};
    const reference::NeutralBulk bulk = reference::mos_bulk(in);
    const reference::DD debye = reference::mos_debye_cm(in);
    std::vector<double> offsets{-1.0, -0.5, -0.2, -0.05, 0.0, 0.05, 0.2, 0.5, 1.0, 1.5, 2.0};
    const std::vector<double> all_offsets = offsets;
    if (tier < Tier::full) offsets = {1.5};
    // The reference at every grid point (the range for applicability is over the whole grid, whatever the tier).
    std::vector<reference::MosPoint> ref;
    double q_min = 0.0, q_max = 0.0, u_ref = 0.0;
    for (const double dv : all_offsets) {
        const auto p = reference::mos_point(in, VFB + dv);
        if (!p) {
            r.run_error = "reference: " + p.error();
            return r;
        }
        const double Q = p->charge_C_per_cm2.to_double();
        q_min = std::min(q_min, Q);
        q_max = std::max(q_max, Q);
        if (p->complex_step_evaluated) {  // the second method
            const double a = p->capacitance.to_double(), b = p->capacitance_complex_step;
            if (!(std::abs(a - b) <= 1e-10 * a)) {
                r.run_error = "reference self-check: C_LF analytic " + shortest(a) + " and by complex step " +
                              shortest(b) + " at V_G - V_FB = " + shortest(dv);
                return r;
            }
        }
        if (dv != 0.0) {  // the substrate's finite depth
            const std::vector<double> d{depth_cm};
            const auto tail = reference::bulk_profile(bulk, p->u_s, debye, d);
            if (!tail) {
                r.run_error = "reference: " + tail.error();
                return r;
            }
            u_ref = std::max(u_ref, (*tail)[0].tail.to_double() / std::abs((p->u_s - bulk.u_b).to_double()));
        }
        ref.push_back(*p);
    }
    const double range = q_max - q_min;
    r.notes.push_back("reference: V_FB " + shortest(VFB) + " V, C_ox " + shortest(cox) + " F/cm^2, Q_s range " +
                      shortest(range) + " C/cm^2, Debye length " + shortest(debye.to_double()) +
                      " cm; truncation (relative tail at x = d) " + shortest(u_ref));
    const auto make_device = [&](int N) -> std::expected<device::Device, std::string> {
        auto m = mesh::make_tensor_grid(graded(N));
        if (!m) return std::unexpected(m.error().message);
        const std::size_t n = m->node_count();
        device::Contact gate{"gate", device::ContactKind::gate, m->find_boundary("x_min")->nodes};
        gate.gate = {.boundary = "x_min",
                     .oxide_thickness_cm = t_ox_cm,
                     .oxide_relative_permittivity = eps_ox_r,
                     .electrode = device::GateElectrode::metal,
                     .work_function_eV = phi_m_eV,
                     .fixed_charge_cm2 = fixed_charge_cm2};
        auto sub = m->find_boundary("x_max")->nodes;
        auto d = device::Device::create({.mesh = std::move(*m),
                                         .temperature_K = T,
                                         .regions = {{"silicon", si}},
                                         .node_region = std::vector<device::RegionId>(n, 0),
                                         .donors = std::vector<double>(n, 0.0),
                                         .acceptors = std::vector<double>(n, acceptors_cm3),
                                         .contacts = {std::move(gate),
                                                      {"substrate", device::ContactKind::ohmic, std::move(sub)}}});
        if (!d) return std::unexpected(d.error().message);
        return std::move(*d);
    };
    const auto read = [&](const device::Device& d, const results::BiasPoint& p) {
        const std::size_t g = static_cast<std::size_t>(d.contacts()[0].nodes.front());
        const std::size_t s = static_cast<std::size_t>(d.contacts()[1].nodes.front());
        const auto& it = p.convergence.iterations;
        return Solved{-(p.gate_charge[0] + q * fixed_charge_cm2), p.fields.potential_V[g], p.fields.potential_V[s],
                      p.gate_charge[0], it.empty() ? 0.0 : it.back().update * VT, p.fields};
    };
    // Two sweeps from flat band (where the charge-neutral start is the solution): down, then up.
    const auto solve_grid = [&](const device::Device& d, std::vector<double> dvs)
        -> std::expected<std::vector<Solved>, std::string> {
        std::vector<double> down{0.0}, up{0.0};
        for (const double v : dvs) (v < 0.0 ? down : up).push_back(v);
        std::sort(down.begin() + 1, down.end(), [](double a, double b) { return a > b; });
        std::sort(up.begin() + 1, up.end());
        std::vector<Solved> out(dvs.size());
        for (const auto* branch : {&down, &up}) {
            std::vector<std::vector<double>> points;
            for (const double v : *branch) points.push_back({VFB + v, 0.0});
            const auto sw = solve::sweep_bias(d, points, options());
            if (!sw) return std::unexpected(sw.error().message);
            if (sw->stopped) return std::unexpected("the sweep stopped: " + sw->stopped->message);
            for (std::size_t k = 0; k < branch->size(); ++k) {
                const auto at = std::find(dvs.begin(), dvs.end(), (*branch)[k]);
                if (at != dvs.end()) out[static_cast<std::size_t>(at - dvs.begin())] = read(d, sw->points[k]);
            }
        }
        return out;
    };
    const std::vector<int> ladder{64, 128, 256, 512};
    std::vector<std::vector<Solved>> levels;  // [level][offset]
    std::vector<double> vfb_sim;
    std::vector<double> c_sim(offsets.size(), 0.0), c_unc(offsets.size(), 0.0);
    for (const int N : ladder) {
        const auto d = make_device(N);
        if (!d) {
            r.run_error = "N = " + std::to_string(N) + ": " + d.error();
            return r;
        }
        std::vector<double> dvs = offsets;
        if (std::find(dvs.begin(), dvs.end(), 0.0) == dvs.end()) dvs.push_back(0.0);  // the V_FB check
        const auto solved = solve_grid(*d, dvs);
        if (!solved) {
            r.run_error = "N = " + std::to_string(N) + ": " + solved.error();
            return r;
        }
        const Solved& fb = (*solved)[static_cast<std::size_t>(std::find(dvs.begin(), dvs.end(), 0.0) - dvs.begin())];
        vfb_sim.push_back(VFB - fb.psi_s - fb.gate_charge / cox + fb.psi_b - q * fixed_charge_cm2 / cox);
        levels.emplace_back(solved->begin(), solved->begin() + static_cast<std::ptrdiff_t>(offsets.size()));
        if (N == ladder.back()) {  // C_LF of NiTCAD's quasi-static charge at the finest level
            for (std::size_t k = 0; k < offsets.size(); ++k) {
                const double V0 = VFB + offsets[k];
                double D[2];
                for (int s = 0; s < 2; ++s) {
                    const double h = fd_step_V / (1 << s);
                    double Qg[2];
                    for (int side = 0; side < 2; ++side) {
                        const std::vector<double> bias{V0 + (side == 0 ? h : -h), 0.0};
                        const auto p = solve::solve_bias(*d, bias, options(), &levels.back()[k].fields);
                        if (!p) {
                            r.run_error = "C_LF difference at V_G - V_FB = " + shortest(offsets[k]) + ": " +
                                          p.error().message;
                            return r;
                        }
                        Qg[side] = p->gate_charge[0];
                    }
                    D[s] = (Qg[0] - Qg[1]) / (2.0 * h);
                }
                c_sim[k] = (4.0 * D[1] - D[0]) / 3.0;
                c_unc[k] = std::abs(D[1] - D[0]) / 3.0;
            }
        }
        std::string s = "N = " + std::to_string(N) + ": V_FB " + shortest(vfb_sim.back()) + " V; Q_s errors";
        for (std::size_t k = 0; k < offsets.size(); ++k) {
            const double Qr = ref[static_cast<std::size_t>(
                                      std::find(all_offsets.begin(), all_offsets.end(), offsets[k]) - all_offsets.begin())]
                                  .charge_C_per_cm2.to_double();
            s += " " + shortest(offsets[k]) + ": " + shortest(Qr == 0.0 ? levels.back()[k].charge
                                                                         : (levels.back()[k].charge - Qr) / Qr);
        }
        r.notes.push_back(s);
    }
    const auto record = [&](const std::string& quantity, double sim, double refv, double metric,
                            const std::string& crit, Comparison cmp) {
        ComparisonRecord rec{quantity, "1", 1.0 / (ladder.back() - 1), sim, refv, 0.0, metric, crit, {}};
        rec.comparison = std::move(cmp);
        r.comparisons.push_back(std::move(rec));
    };
    // D: V_FB at every level.
    for (std::size_t l = 0; l < ladder.size(); ++l) {
        const double e = std::abs(vfb_sim[l] - VFB) / std::abs(VFB);
        const std::string crit = "relative <= 64 eps = " + shortest(64.0 * eps);
        record("V_FB, N = " + std::to_string(ladder[l]), vfb_sim[l], VFB, e, crit,
               judge(r.info, {}, e <= 64.0 * eps, crit));
    }
    const bool ref_ok = u_ref <= 1e-3 * tolerance;
    const std::string band = "p-hat in [1.8, 2.3]";
    std::vector<double> hs, u;
    for (auto it = ladder.rbegin(); it != ladder.rend(); ++it) hs.push_back(1.0 / (*it - 1));  // finest first
    for (std::size_t k = 0; k < offsets.size(); ++k) {
        const reference::MosPoint& p =
            ref[static_cast<std::size_t>(std::find(all_offsets.begin(), all_offsets.end(), offsets[k]) -
                                         all_offsets.begin())];
        const double Qr = p.charge_C_per_cm2.to_double(), Cr = p.capacitance.to_double();
        const std::string at = ", V_G - V_FB = " + shortest(offsets[k]) + " V";
        const std::vector<Predicate> pred{
            {applicability, range > 0.0 ? std::abs(Qr) / range : 0.0, 1e-3, range > 0.0 && std::abs(Qr) >= 1e-3 * range}};
        std::vector<double> e, ui;
        for (auto l = levels.rbegin(); l != levels.rend(); ++l) {
            e.push_back(Qr == 0.0 ? 0.0 : std::abs((*l)[k].charge - Qr) / std::abs(Qr));
            ui.push_back(Qr == 0.0 ? 0.0 : (*l)[k].u_iter * cox / std::abs(Qr));  // a potential update as charge
        }
        const study::OrderStudy st = study::order_study(hs, e, ui);
        const study::OrderVerdict v = study::judge_order(st, {.p_th = 2.0});
        if (pred[0].holds) {
            std::string orders;
            for (const double o : st.orders) orders += " " + shortest(o);
            if (st.asymptotic_evaluated) orders += ", order spread " + shortest(st.order_spread);
            r.notes.push_back("Q_s orders" + at + ":" + orders);
        }
        const Comparison order_cmp =
            !ref_ok ? missing_input("a reference within 1e-3 of the criterion (substrate depth)")
            : !pred[0].holds ? judge(r.info, pred, false, band)
            : v.judged       ? judge(r.info, pred, v.pass, band)
                             : Comparison{Outcome::not_evaluated, v.reason, pred};
        record("Q_s order" + at, st.observed, 2.0, st.observed, band, order_cmp);
        const double eq = e.front();
        const std::string crit = "relative <= 1e-3";
        record("Q_s finest" + at, levels.back()[k].charge, Qr, eq, crit,
               ref_ok ? judge(r.info, pred, eq <= tolerance, crit)
                      : missing_input("a reference within 1e-3 of the criterion (substrate depth)"));
        const double ec = std::abs(c_sim[k] - Cr) / Cr;
        r.notes.push_back("C_LF" + at + ": " + shortest(c_sim[k]) + " F/cm^2 against " + shortest(Cr) +
                          " (relative " + shortest(ec) + ", difference uncertainty " + shortest(c_unc[k] / Cr) + ")");
        const bool fd_ok = c_unc[k] <= 0.1 * tolerance * Cr;
        record("C_LF finest" + at, c_sim[k], Cr, ec, crit,
               !ref_ok  ? missing_input("a reference within 1e-3 of the criterion (substrate depth)")
               : !fd_ok ? Comparison{Outcome::not_evaluated,
                                     "the difference quotient's uncertainty " + shortest(c_unc[k] / Cr) +
                                         " exceeds a tenth of the criterion",
                                     pred}
                        : judge(r.info, pred, ec <= tolerance, crit));
    }
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
