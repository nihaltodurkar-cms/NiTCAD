// TVB-23, impact ionization: avalanche multiplication and the approach to breakdown (DECISIONS Unit 55 6, 17.2
// (h)-(l), 17.4, 17.13 (b); criteria 17.8 and the predicate limits and A-class rule of 17.14, frozen before code).
//
// Code under test: the drift-diffusion system with impact ionization (the generation alpha_n |j_n| + alpha_p |j_p|
// at each node, the field along each carrier's current), solve_bias, and trace_bias for the A part.
// Construction, fixed before any run:
//   - 1D silicon (NiTCAD's), a symmetric abrupt pn junction, N_A = N_D = 3e15 cm^-3, the p side on [-10, 0) um and the
//     n side on (0, 10] um, ohmic contacts at both ends (the anode, contact 0, on the p side); 300 K; Boltzmann, no
//     band-gap narrowing, no recombination (SRH, Auger, radiative off), constant mobility (the doping dependence off:
//     the lattice values), impact ionization with silicon's van Overstraeten-de Man coefficients; the weighted
//     terminal current (Unit 32) and the minority floor (Unit 50) on.
//     Why symmetric and this light: Q2 is measured where the reference field exceeds a tenth of its peak, and there a
//     side's majority density is about N e^-(1 + z), z = eps (E_max / 10)^2 / (2 q N V_T). A one-sided junction's
//     heavy side would put Q2 near 1e2 at every bias; at 3e15 on both sides z is 15 at the lowest judged bias (Q2
//     near 1e-7 there). Breakdown is then near 194 V.
//   - The primary current (IG-Q2 A): the TVQ3 hook's pair generation G(x) = G0 exp(-(x - x_s)^2 / (2 s^2)), s = 0.1 um,
//     equal for electrons and holes, in one neutral side: x_s = -8 um (electrons reach the junction from the p side:
//     M_n) or +8 um (holes from the n side: M_p), G0 for a generated pair current q int G dx = 2.5e-5 A/cm^2. The
//     source lies 15 s outside the widest depletion layer traced (6.5 um at M = 1000).
//   - mesh: x_k = L sinh(beta (2 k / (N - 1) - 1)) / sinh(beta), L = 10 um, beta = 2, N even: the junction lies at the
//     midpoint of the middle edge (the boxes' boundary), finest at the junction; N = 128, 256, 512, 1024 (17.8), and
//     2048 in the factorial.
// Reference (verify/reference/multiplication.hpp): I_n, phi and M_n = 1 / (1 - I_n), M_p = e^-phi(W) / (1 - I_n) on
// the QE field, with NiTCAD's q, eps, V_T, n_i and coefficients as inputs (I1). The judged biases are the reference's
// at I_n = 0.1, 0.3, 0.5, 0.7, 0.9. The simulated M is J(ionization on) / J(off) at the same bias and generation,
// from the weighted current of the anode.
// E (17.8), per injection and bias, on x = 1 - 1 / M (I_n for electrons, 1 - (1 - I_n) e^phi for holes):
//   - (a) p-hat in [1.8, 2.3] of the relative error of x on N = 128 ... 1024;
//   - (b) at N = 1024, the relative error of x <= 1e-3.
//   Applicability (17.14): Q1 = V_Omega,d / V_T <= 1e-5 (V_Omega,d = J d / (q mu_maj N) with d from the source to the
//   reference's depletion edge, J the multiplied current) and Q2 = max |p - n| / N <= 1e-5 over the nodes inside the
//   reference's field window (E >= E_max / 10), on every level (else NOT APPLICABLE); J(off) >= 1e3 times the
//   impact-ionization current floor (else NOT APPLICABLE) and >= 1e3 times its resolution (else NOT EVALUATED). A
//   comparison whose U_ref or U_iter exceeds a tenth of its tolerance is NOT EVALUATED.
//   U_ref: the reference's tolerance error and its difference from the second method (nested Kronrod, 16 panels),
//   carried to x. U_iter: the weighted currents' resolutions, carried to x: (r_on / J_on + r_off / J_off) / (M - 1).
//   The (h, tol) factorial (16.2 (p)) at N = 1024, tol = 1e-10 and 1e-12, N = 1024 and 2048, per comparison.
// A (17.14), per injection, at M_c = 10, 100, 1000: |V_sim - V_ref| <= U_num + U_ref + eps_V, V_sim where J(on) /
//   J(off) = M_c along trace_bias with the generation (interpolated in ln M between traced points; J(off) re-solved
//   at each traced bias), V_ref the reference's bias for 1 / M = 1 / M_c; eps_V = 2 Q2 (V_bi - V) + J (R_n + R_p),
//   R the neutral sides' resistances from the spec and the reference's depletion widths, Q2 at the traced point
//   past M_c; U_num the two-grid estimate on N = 512 and 1024 (F_s = 3, p = 2); U_ref the reference's error carried
//   to V by its slope. Applicability: J(off) as above. The reference's V_BR (I_n = 1) is printed beside.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <expected>
#include <numbers>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/drift_diffusion.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/impact_ionization.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/trace.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/metrics.hpp"
#include "NiTCAD/verify/reference/multiplication.hpp"
#include "NiTCAD/verify/study/study.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

constexpr double T = 300.0;
constexpr double doping_cm3 = 3e15, half_cm = 1e-3, beta = 2.0;
constexpr double source_cm = 8e-4, sigma_cm = 1e-5, pair_current_A_cm2 = 2.5e-5;
constexpr double tolerance = 1e-3, tol_newton = 1e-10;
constexpr double limit_q1 = 1e-5, limit_q2 = 1e-5, resolution_factor = 1e3;
constexpr double step_V = 4.0;
constexpr std::array<double, 5> targets{0.1, 0.3, 0.5, 0.7, 0.9};  // the reference's I_n at the judged biases
constexpr std::array<double, 3> criterion_M{10.0, 100.0, 1000.0};
constexpr std::array<int, 4> ladder{1024, 512, 256, 128};  // finest first

ReferenceInfo info() {
    ReferenceInfo r;
    r.id = "TVB-23";
    r.cls = ReferenceClass::E;
    r.independence = Independence::I1;
    r.provenance = Provenance::PD;
    r.citation = "The ionization integrals with alpha_n != alpha_p on the QE field, derived for both injections "
                 "(DECISIONS Unit 55 17.13 (b)); the form as J. M. Park, TU Wien 2004, eqs (3.26)-(3.28) (PB)";
    r.version = "DECISIONS Unit 55 17.3, 17.8, 17.13 (b), 17.14";
    r.equations = {1, 2, 3};
    r.form = "M_n = 1 / (1 - I_n), M_p = e^-phi(W) / (1 - I_n), I_n = int alpha_n e^-phi dx, phi = int (alpha_n - "
             "alpha_p) dx; alpha = gamma A exp(-gamma B / E) on the quasi-equilibrium field";
    r.assumptions = {"one-carrier injection from a neutral side (IG-Q2 A)", "the carriers' space charge negligible "
                     "(Q2, measured)", "the collected fraction unchanged by the majority return current (Q1)",
                     "300 K (gamma = 1; 17.2 (l))"};
    r.derivation = "verify/reference/multiplication.hpp";
    r.verification = "E: orders and finest-level 1 - 1/M for both injections at five biases; A: V at M = 10, 100, "
                     "1000 along the trace";
    r.second_method = "nested Kronrod-15 quadrature (16 panels per piece), at run time";
    return r;
}

std::vector<double> axis(int N) {
    std::vector<double> x;
    const double sh = std::sinh(beta);
    for (int k = 0; k < N; ++k) x.push_back(half_cm * std::sinh(beta * (2.0 * k / (N - 1) - 1.0)) / sh);
    return x;
}

std::expected<device::Device, std::string> diode(int N) {
    auto m = mesh::make_tensor_grid(axis(N));
    if (!m) return std::unexpected(m.error().message);
    const std::size_t n = m->node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) (m->points()[i][0] < 0.0 ? acceptors : donors)[i] = doping_cm3;
    auto anode = m->find_boundary("x_min")->nodes;
    auto cathode = m->find_boundary("x_max")->nodes;
    auto d = device::Device::create({.mesh = std::move(*m),
                                     .temperature_K = T,
                                     .regions = {{"silicon", physics::silicon()}},
                                     .node_region = std::vector<device::RegionId>(n, 0),
                                     .donors = std::move(donors),
                                     .acceptors = std::move(acceptors),
                                     .contacts = {{"anode", device::ContactKind::ohmic, std::move(anode)},
                                                  {"cathode", device::ContactKind::ohmic, std::move(cathode)}}});
    if (!d) return std::unexpected(d.error().message);
    return std::move(*d);
}

solve::BiasOptions options(double tol, bool impact, const assemble::ManufacturedSources* source) {
    solve::BiasOptions o;
    o.models.srh = o.models.auger = o.models.radiative = false;
    o.models.bgn = false;
    o.models.fermi_dirac = false;
    o.models.doping_mobility = false;
    o.models.hole_mobility = assemble::HoleMobility::caughey_thomas;
    o.models.impact_ionization = impact;
    o.models.manufactured = source;
    o.weighted_current = true;
    o.minority_floor = true;
    o.store_rates = false;
    o.newton.tol_update = tol;
    return o;
}

// The pair generation of one injection at the device's nodes.
assemble::ManufacturedSources source(const device::Device& d, bool electron) {
    const double G0 = pair_current_A_cm2 / (base::q_C * sigma_cm * std::sqrt(2.0 * std::numbers::pi));
    const double xs = electron ? -source_cm : source_cm;
    assemble::ManufacturedSources ms;
    for (const auto& p : d.mesh().points()) {
        const double z = (p[0] - xs) / sigma_cm;
        ms.electron_generation.push_back(G0 * std::exp(-0.5 * z * z));
    }
    ms.hole_generation = ms.electron_generation;
    return ms;
}

reference::IonizationBranches branches(const physics::ImpactIonizationCoefficients& c) {
    return {c.A_low_per_cm, c.B_low_V_per_cm, c.A_high_per_cm, c.B_high_V_per_cm, c.switch_V_per_cm};
}

// The reference at one judged bias, for both injections.
struct RefPoint {
    double V = 0.0;
    double x_n = 0.0, x_p = 0.0;      // 1 - 1/M_n, 1 - 1/M_p
    double u_n = 0.0, u_p = 0.0;      // relative uncertainty of x
    double W_p = 0.0, W_n = 0.0;      // depletion widths [cm]
    double window_p = 0.0, window_n = 0.0;  // the field window's edges [cm]
    double barrier = 0.0;             // V_bi - V
    double I_n = 0.0, phi = 0.0;
};

std::expected<RefPoint, std::string> reference_at(const reference::MultiplicationInput& in, double V) {
    const auto m = reference::multiplication(in, V);
    const auto q = reference::multiplication_quadrature(in, V, 16);
    const auto w = reference::multiplication_field_window(in, V, 0.1);
    if (!m) return std::unexpected(m.error());
    if (!q) return std::unexpected(q.error());
    if (!w) return std::unexpected(w.error());
    RefPoint r;
    r.V = V;
    r.I_n = m->ionization_integral.to_double();
    r.phi = m->phi_total.to_double();
    r.x_n = r.I_n;
    r.x_p = 1.0 - m->inverse_p.to_double();
    const double dI = m->integration_error + std::abs((q->ionization_integral - m->ionization_integral).to_double());
    const double dphi = m->integration_error + std::abs((q->phi_total - m->phi_total).to_double());
    const double e = std::exp(r.phi);
    r.u_n = dI / r.x_n;
    r.u_p = (e * dI + (1.0 - r.I_n) * e * dphi) / r.x_p;
    r.W_p = m->depletion_p_cm.to_double();
    r.W_n = m->depletion_n_cm.to_double();
    r.window_p = w->x_p_cm.to_double();
    r.window_n = w->x_n_cm.to_double();
    r.barrier = m->barrier_V.to_double();
    return r;
}

// max |p - n| / N over the nodes inside [window_p, window_n].
double q2(const device::Device& d, const results::BiasPoint& p, double lo, double hi) {
    double worst = 0.0;
    for (std::size_t i = 0; i < p.fields.n_cm3.size(); ++i) {
        const double x = d.mesh().points()[i][0];
        if (x < lo || x > hi) continue;
        worst = std::max(worst, std::abs(p.fields.p_cm3[i] - p.fields.n_cm3[i]) / doping_cm3);
    }
    return worst;
}

// One injection's values at one judged bias on one level.
struct Sample {
    double x = 0.0;          // 1 - J_off / J_on
    double error = 0.0;      // |x - x_ref| / x_ref
    double u_iter = 0.0;     // relative to x
    double J_on = 0.0, J_off = 0.0, res_off = 0.0;
    double q1 = 0.0, q2 = 0.0;
};

struct Level {
    int N = 0;
    double h = 0.0;
    std::array<std::array<Sample, 5>, 2> s{};  // [electron = 0, hole = 1][target]
    results::BiasPoint last_on[2], last_off[2];  // the states at the last judged bias (the traces start there)
};

}  // namespace

BenchmarkResult run_tvb23(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-23";
    r.title = "Impact ionization: avalanche multiplication";
    r.info = info();
    if (tier < Tier::full) return r;
    const physics::Semiconductor si = physics::silicon();
    const double VT = base::thermal_voltage(T);
    const double eps = si.parameters().eps_r * base::eps0_F_per_cm;
    const double ni = physics::intrinsic_density(si, T);
    const double gamma = physics::impact_ionization_temperature_factor(si.parameters().impact_ionization.phonon_energy_eV, T);
    const reference::MultiplicationInput in{base::q_C, eps, VT, ni, doping_cm3, doping_cm3,
                                            branches(physics::impact_ionization(si, physics::Carrier::electron)),
                                            branches(physics::impact_ionization(si, physics::Carrier::hole)), gamma};
    const double mu_n = physics::caughey_thomas_mobility(si, physics::Carrier::electron, 0.0, T);
    const double mu_p = physics::caughey_thomas_mobility(si, physics::Carrier::hole, 0.0, T);
    const auto vbr = reference::bias_at_inverse_multiplication(in, true, 0.0, -600.0, -10.0);
    if (!vbr) {
        r.run_error = "reference: V_BR: " + vbr.error();
        return r;
    }
    std::array<RefPoint, 5> ref{};
    for (std::size_t k = 0; k < targets.size(); ++k) {
        const auto V = reference::bias_at_inverse_multiplication(in, true, 1.0 - targets[k], *vbr, -1.0);
        if (!V) {
            r.run_error = "reference: the bias at I_n = " + shortest(targets[k]) + ": " + V.error();
            return r;
        }
        const auto p = reference_at(in, *V);
        if (!p) {
            r.run_error = "reference at " + shortest(*V) + " V: " + p.error();
            return r;
        }
        ref[k] = *p;
        r.notes.push_back("reference I_n " + shortest(targets[k]) + ": V " + shortest(p->V) + " V, phi " +
                          shortest(p->phi) + ", 1 - 1/M_n " + shortest(p->x_n) + ", 1 - 1/M_p " + shortest(p->x_p) +
                          ", U_ref (relative) " + shortest(p->u_n) + " / " + shortest(p->u_p) + ", W_p " +
                          shortest(p->W_p) + " cm, window [" + shortest(p->window_p) + ", " + shortest(p->window_n) +
                          "] cm");
    }
    r.notes.push_back("reference V_BR (I_n = 1) " + shortest(*vbr) + " V; mu_n " + shortest(mu_n) + ", mu_p " +
                      shortest(mu_p) + " cm^2/Vs; leakage bound q n_i^2 (D_n / (N_A d_p) + D_p / (N_D d_n)) with the "
                      "neutral lengths at the last judged bias: " +
                      shortest(base::q_C * ni * ni * VT *
                               (mu_n / (doping_cm3 * (half_cm - ref[4].W_p)) +
                                mu_p / (doping_cm3 * (half_cm - ref[4].W_n)))) +
                      " A/cm^2");

    // The ramp: steps of at most step_V to each judged bias in turn; at[k] the index of target k.
    std::vector<std::vector<double>> ramp;
    std::array<std::size_t, 5> at{};
    {
        double v = 0.0;
        for (std::size_t k = 0; k < targets.size(); ++k) {
            const int steps = static_cast<int>(std::ceil((v - ref[k].V) / step_V));
            const double from = v;
            for (int s = 1; s <= steps; ++s) {
                v = s == steps ? ref[k].V : from - step_V * s;
                ramp.push_back({v, 0.0});
            }
            at[k] = ramp.size() - 1;
        }
    }

    const auto run_level = [&](int N, double tol) -> std::expected<Level, std::string> {
        const std::string where = "N = " + std::to_string(N) + ", tol " + shortest(tol);
        const auto d = diode(N);
        if (!d) return std::unexpected(where + ": " + d.error());
        Level lv;
        lv.N = N;
        lv.h = 2.0 * half_cm / (N - 1);
        for (const bool electron : {true, false}) {
            const std::size_t c = electron ? 0 : 1;
            const assemble::ManufacturedSources ms = source(*d, electron);
            const auto on = solve::sweep_bias(*d, ramp, options(tol, true, &ms));
            const auto off = solve::sweep_bias(*d, ramp, options(tol, false, &ms));
            if (!on || on->stopped || !off || off->stopped) {
                const std::string why = !on ? on.error().message : on->stopped ? on->stopped->message
                                        : !off ? off.error().message : off->stopped->message;
                return std::unexpected(where + (electron ? ", electron" : ", hole") + " injection: " + why);
            }
            for (std::size_t k = 0; k < targets.size(); ++k) {
                const results::BiasPoint& pn = on->points[at[k]];
                const results::BiasPoint& pf = off->points[at[k]];
                Sample& s = lv.s[c][k];
                s.J_on = std::abs(pn.weighted_current[0]);
                s.J_off = std::abs(pf.weighted_current[0]);
                s.res_off = pf.weighted_current_resolution[0];
                const double M = s.J_on / s.J_off;
                s.x = 1.0 - 1.0 / M;
                const double x_ref = electron ? ref[k].x_n : ref[k].x_p;
                s.error = std::abs(s.x - x_ref) / x_ref;
                s.u_iter = (pn.weighted_current_resolution[0] / s.J_on + s.res_off / s.J_off) / (M - 1.0);
                const double d_cm = source_cm - (electron ? ref[k].W_p : ref[k].W_n);
                s.q1 = s.J_on * d_cm / (base::q_C * (electron ? mu_p : mu_n) * doping_cm3) / VT;
                s.q2 = q2(*d, pn, ref[k].window_p, ref[k].window_n);
                // The construction's own check: the collected primary current is a fraction of the pairs.
                if (!(s.J_off > 0.05 * pair_current_A_cm2 && s.J_off < pair_current_A_cm2)) {
                    return std::unexpected(where + ": the collected primary current " + shortest(s.J_off) +
                                           " A/cm^2 is not a fraction of the generated pairs' " +
                                           shortest(pair_current_A_cm2));
                }
            }
            lv.last_on[c] = on->points[at.back()];
            lv.last_off[c] = off->points[at.back()];
        }
        return lv;
    };

    std::vector<Level> levels;  // finest first
    for (const int N : ladder) {
        auto lv = run_level(N, tol_newton);
        if (!lv) {
            r.run_error = lv.error();
            return r;
        }
        for (std::size_t c = 0; c < 2; ++c) {
            std::string line = std::string(c == 0 ? "M_n" : "M_p") + " N = " + std::to_string(N) + ":";
            for (std::size_t k = 0; k < targets.size(); ++k) {
                const Sample& s = lv->s[c][k];
                line += " [I_n " + shortest(targets[k]) + ": x " + shortest(s.x) + ", error " + shortest(s.error) +
                        ", U_iter " + shortest(s.u_iter) + ", J_off " + shortest(s.J_off) + " (res " +
                        shortest(s.res_off) + "), Q1 " + shortest(s.q1) + ", Q2 " + shortest(s.q2) + "]";
            }
            r.notes.push_back(line);
        }
        levels.push_back(std::move(*lv));
    }
    // The factorial's three more levels.
    std::array<Level, 3> fac_levels;
    {
        std::size_t k = 0;
        for (const auto& [N, tol] : {std::pair{1024, 0.01 * tol_newton}, std::pair{2048, tol_newton},
                                     std::pair{2048, 0.01 * tol_newton}}) {
            auto lv = run_level(N, tol);
            if (!lv) {
                r.run_error = lv.error() + " (factorial)";
                return r;
            }
            fac_levels[k++] = std::move(*lv);
        }
    }

    const auto record = [&](const std::string& quantity, const std::string& unit, double parameter, double sim,
                            double refv, double uref, double metric, const std::string& crit) -> ComparisonRecord& {
        r.comparisons.push_back(ComparisonRecord{quantity, unit, parameter, sim, refv, uref, metric, crit, {}});
        return r.comparisons.back();
    };
    const double floor_A = assemble::DriftDiffusion::impact_current_floor;
    for (std::size_t c = 0; c < 2; ++c) {
        const std::string carrier = c == 0 ? "M_n" : "M_p";
        for (std::size_t k = 0; k < targets.size(); ++k) {
            const std::string at_k = carrier + ", I_n " + shortest(targets[k]) + " (" + shortest(ref[k].V) + " V)";
            const double u_ref = c == 0 ? ref[k].u_n : ref[k].u_p;
            const double x_ref = c == 0 ? ref[k].x_n : ref[k].x_p;
            std::vector<double> hs, es, us;
            double q1 = 0.0, q2v = 0.0, worst_res = 0.0;
            bool floor_ok = true;
            for (const Level& l : levels) {
                const Sample& s = l.s[c][k];
                hs.push_back(l.h);
                es.push_back(s.error);
                us.push_back(s.u_iter);
                q1 = std::max(q1, s.q1);
                q2v = std::max(q2v, s.q2);
                worst_res = std::max(worst_res, s.res_off / s.J_off);
                floor_ok = floor_ok && s.J_off >= resolution_factor * floor_A;
            }
            const study::OrderStudy st = study::order_study(hs, es, us);
            const study::OrderVerdict v = study::judge_order(st, {.p_th = 2.0});
            std::string orders;
            for (const double p : st.orders) orders += " " + shortest(p);
            const Sample& fine = levels.front().s[c][k];
            const double u_disc = reference::gci_safety_factor_three_grids * std::abs(es[1] - es[0]) /
                                  (std::pow(hs[1] / hs[0], st.observed) - 1.0);
            const std::vector<double> parts{u_disc, fine.u_iter, u_ref};
            const Factorial fac = factorial(fine.error, fac_levels[0].s[c][k].error, fac_levels[1].s[c][k].error,
                                            fac_levels[2].s[c][k].error, parts);
            r.notes.push_back(at_k + ": orders (finest pair first):" + orders + "; spread " +
                              shortest(st.order_spread) + "; " + fac.note + " (U_disc " + shortest(u_disc) + ")");
            // Applicability, then evaluability.
            std::optional<Comparison> blocked;
            if (q1 > limit_q1) {
                blocked = Comparison{Outcome::not_applicable, "Q1 " + shortest(q1) + " above 1e-5", {}};
            } else if (q2v > limit_q2) {
                blocked = Comparison{Outcome::not_applicable, "Q2 " + shortest(q2v) + " above 1e-5", {}};
            } else if (!floor_ok) {
                blocked = Comparison{Outcome::not_applicable, "J(off) below 1e3 times the ionization floor", {}};
            } else if (worst_res * resolution_factor > 1.0) {
                blocked = Comparison{Outcome::not_evaluated,
                                     "J(off) below 1e3 times its resolution (resolution / J " + shortest(worst_res) + ")",
                                     {}};
            } else if (u_ref > 0.1 * tolerance) {
                blocked = Comparison{Outcome::not_evaluated, "U_ref " + shortest(u_ref) + " above a tenth of the tolerance", {}};
            } else if (fine.u_iter > 0.1 * tolerance) {
                blocked = Comparison{Outcome::not_evaluated,
                                     "U_iter " + shortest(fine.u_iter) + " above a tenth of the tolerance", {}};
            }
            const std::string band = "p-hat in [1.8, 2.3] (p_th 2)";
            {
                ComparisonRecord& rec = record("1 - 1/M order, " + at_k, "1", targets[k], st.observed, 2.0, u_ref,
                                               st.observed, band);
                rec.comparison = blocked    ? *blocked
                                 : v.judged ? judge(r.info, {}, v.pass, band)
                                            : Comparison{Outcome::not_evaluated, v.reason, {}};
                block_if_interacting(rec, fac);
            }
            {
                const std::string crit = "1 - 1/M relative <= 1e-3";
                ComparisonRecord& rec = record("1 - 1/M finest (N = 1024), " + at_k, "1", targets[k], fine.x, x_ref,
                                               u_ref, fine.error, crit);
                rec.comparison = blocked ? *blocked : judge(r.info, {}, fine.error <= tolerance, crit);
                block_if_interacting(rec, fac);
            }
        }
    }

    // A: V at M_c along the trace, per injection, on N = 512 and 1024.
    ReferenceInfo a_info = r.info;
    a_info.cls = ReferenceClass::A;
    a_info.approximation_error = "eps_V = 2 Q2 (V_bi - V) + J (R_n + R_p), the carriers' space charge and the neutral "
                                 "sides' ohmic drop (17.14)";
    struct Crossing {
        double V = 0.0, J = 0.0, q2 = 0.0, res_ratio = 0.0, J_off = 0.0;
        bool found = false;
    };
    const auto trace_level = [&](const Level& lv, bool electron)
        -> std::expected<std::array<Crossing, 3>, std::string> {
        const std::size_t c = electron ? 0 : 1;
        const std::string where = std::string(electron ? "M_n" : "M_p") + " trace, N = " + std::to_string(lv.N);
        const auto d = diode(lv.N);
        if (!d) return std::unexpected(where + ": " + d.error());
        const assemble::ManufacturedSources ms = source(*d, electron);
        solve::TraceOptions o;
        o.steady = options(tol_newton, true, &ms);
        o.start_V = {ref.back().V, 0.0};
        o.end_V = *vbr - 20.0;
        o.step_V = 0.5;
        o.max_step_V = 2.0;
        o.current_limit = 2.0 * criterion_M.back() * lv.s[c].back().J_off;
        const auto t = solve::trace_bias(*d, o, &lv.last_on[c].fields);
        if (!t) return std::unexpected(where + ": " + t.error().message);
        if (t->stopped) return std::unexpected(where + ": " + t->stopped->message);
        std::vector<std::vector<double>> biases;
        for (const auto& p : t->points) biases.push_back(p.bias_V);
        const auto off = solve::sweep_bias(*d, biases, options(tol_newton, false, &ms), &lv.last_off[c].fields);
        if (!off) return std::unexpected(where + ", J(off): " + off.error().message);
        if (off->stopped) return std::unexpected(where + ", J(off): " + off->stopped->message);
        std::array<Crossing, 3> out{};
        double prev_lnM = 0.0, prev_V = 0.0;
        for (std::size_t i = 0; i < t->points.size(); ++i) {
            const auto& pn = t->points[i];
            const auto& pf = off->points[i];
            const double Jon = std::abs(pn.weighted_current[0]), Joff = std::abs(pf.weighted_current[0]);
            const double lnM = std::log(Jon / Joff), V = pn.bias_V[0];
            for (std::size_t m = 0; m < criterion_M.size(); ++m) {
                const double target = std::log(criterion_M[m]);
                if (out[m].found || i == 0 || !(lnM >= target && prev_lnM < target)) continue;
                const double f = (target - prev_lnM) / (lnM - prev_lnM);
                out[m].V = prev_V + f * (V - prev_V);
                out[m].J = criterion_M[m] * Joff;
                out[m].J_off = Joff;
                out[m].res_ratio = pf.weighted_current_resolution[0] / Joff;
                const auto w = reference::multiplication_field_window(in, out[m].V, 0.1);
                if (!w) return std::unexpected(where + ": window: " + w.error());
                out[m].q2 = q2(*d, pn, w->x_p_cm.to_double(), w->x_n_cm.to_double());
                out[m].found = true;
            }
            prev_lnM = lnM;
            prev_V = V;
        }
        r.notes.push_back(where + ": " + std::to_string(t->points.size()) + " traced points to " +
                          shortest(t->points.back().bias_V[0]) + " V, M there " +
                          shortest(std::exp(prev_lnM)));
        return out;
    };
    for (const bool electron : {true, false}) {
        const auto fine = trace_level(levels[0], electron);
        const auto coarse = trace_level(levels[1], electron);
        if (!fine || !coarse) {
            r.run_error = !fine ? fine.error() : coarse.error();
            return r;
        }
        const std::string carrier = electron ? "M_n" : "M_p";
        const double ratio = (levels[1].N - 1.0) / (levels[0].N - 1.0);
        for (std::size_t m = 0; m < criterion_M.size(); ++m) {
            const std::string at_m = carrier + " = " + shortest(criterion_M[m]);
            const auto Vr = reference::bias_at_inverse_multiplication(in, electron, 1.0 / criterion_M[m], *vbr,
                                                                       ref.back().V);
            if (!Vr) {
                r.run_error = "reference: the bias at " + at_m + ": " + Vr.error();
                return r;
            }
            const auto rp = reference_at(in, *Vr);
            const auto slope_hi = reference::multiplication(in, *Vr + 0.01);
            const auto slope_lo = reference::multiplication(in, *Vr - 0.01);
            if (!rp || !slope_hi || !slope_lo) {
                r.run_error = "reference near " + shortest(*Vr) + " V failed";
                return r;
            }
            const auto inv = [&](const reference::Multiplication& x) {
                return (electron ? x.inverse_n : x.inverse_p).to_double();
            };
            const double slope = std::abs(inv(*slope_hi) - inv(*slope_lo)) / 0.02;  // d(1/M)/dV
            const double x_ref = electron ? rp->x_n : rp->x_p;
            const double u_ref = (electron ? rp->u_n : rp->u_p) * x_ref / slope + 1e-10;
            const Crossing& f = (*fine)[m];
            const Crossing& g = (*coarse)[m];
            ComparisonRecord& rec = record("V at " + at_m + " (A), N = 1024", "V", criterion_M[m], f.V, *Vr, u_ref,
                                           std::abs(f.V - *Vr), "|V_sim - V_ref| <= U_num + U_ref + eps_V");
            if (!f.found || !g.found) {
                rec.comparison = Comparison{Outcome::not_evaluated, "the trace did not reach M = " +
                                                                        shortest(criterion_M[m]), {}};
                continue;
            }
            const double u_num = reference::gci_safety_factor_two_grids * std::abs(f.V - g.V) / (ratio * ratio - 1.0);
            const double R = (half_cm - rp->W_p) / (base::q_C * mu_p * doping_cm3) +
                             (half_cm - rp->W_n) / (base::q_C * mu_n * doping_cm3);
            const double eps_V = 2.0 * f.q2 * rp->barrier + f.J * R;
            const double bound = u_num + u_ref + eps_V;
            rec.criterion = "|V_sim - V_ref| <= U_num + U_ref + eps_V = " + shortest(bound);
            r.notes.push_back(at_m + ": V_sim " + shortest(f.V) + " (N = 512: " + shortest(g.V) + "), V_ref " +
                              shortest(*Vr) + ", U_num " + shortest(u_num) + ", U_ref " + shortest(u_ref) +
                              ", eps_V " + shortest(eps_V) + " (Q2 " + shortest(f.q2) + ", J " + shortest(f.J) +
                              " A/cm^2, R " + shortest(R) + " Ohm cm^2), J(off) " + shortest(f.J_off) +
                              " (resolution / J " + shortest(f.res_ratio) + ")");
            if (f.J_off < resolution_factor * floor_A) {
                rec.comparison = Comparison{Outcome::not_applicable, "J(off) below 1e3 times the ionization floor", {}};
            } else if (f.res_ratio * resolution_factor > 1.0) {
                rec.comparison = Comparison{Outcome::not_evaluated, "J(off) below 1e3 times its resolution", {}};
            } else {
                rec.comparison = judge(a_info, {}, std::abs(f.V - *Vr) <= bound, rec.criterion);
            }
        }
    }
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
