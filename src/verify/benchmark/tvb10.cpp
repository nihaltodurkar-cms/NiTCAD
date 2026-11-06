// TVB-10, the MOS capacitor's high-frequency C-V and deep depletion (DECISIONS Unit 55 §6, 16.2 (g), 16.4; criteria
// 16.8 (deep depletion, E) and 16.19 (HF, A), frozen; DT-Q17; construction 16.40).
//
// Code under test: TVB-09's MOS capacitor (p-silicon N_A = 1e17 cm^-3, 0.5 um, metal gate phi_m = 4.7 eV on a lumped
// 10 nm SiO2 with Q_f = 5e10 cm^-2, ohmic substrate; the same graded mesh, h the xi spacing), now with drift-diffusion
// small-signal and transient physics: SRH on (Scharfetter lifetimes), Auger and radiative off, Boltzmann, no
// band-gap narrowing, Caughey-Thomas mobility, no field mobility.
// - HF (A): the quasi-static operating point (Equations::equilibrium_poisson, as small_signal.hpp states for a MOS-C
//   in inversion), then solve_small_signal at f = 1e5, 10^5.5 and 1e6 Hz; C = Im Y_gg / w (the gate's capacitance).
//   Reference (verify/reference/mos_hf): the redistribution limit. Judged at f = 10^5.5 Hz on V_G = V_FB + {0.2, 0.5,
//   1, 1.5, 2, 2.5} V (depletion to strong inversion): |C_sim - C_ref| <= U_num + U_ref + eps_HF (16.19), eps_HF the
//   spread of C_sim over the three frequencies at the finest level; applicable while eps_HF <= 1e-3 C_ref ("on the HF
//   plateau"). U_num = U_disc (registry.hpp's u_disc on N = 512 ... 64 with N = 1024 as h/2) and the admittance
//   resolution, through the (h, tol) factorial at N = 512 (tol = 1e-12, 1e-14). U_ref = the second method's
//   difference + the quadrature error + the substrate's truncation (TVB-09's tail rule) + rounding.
// - Deep depletion (E): from the steady state at V_FB the gate is ramped linearly to V_FB + 4 V in 100 us (BDF2,
//   fixed steps of 1 us), much faster than generation fills the inversion layer. Q_s = -(Q_G + q Q_f) at V_G - V_FB =
//   1, 2, 3, 4 V (output times). Reference: mos_deep_depletion_point (no minority carriers). Applicable while the
//   electrons' charge q int n dx <= 1e-3 |Q_s| (16.2 (g)), measured. E (a): p-hat in [1.8, 2.3] on N = 512 ... 64;
//   E (b): relative <= 1e-3 at N = 512. The factorials (16.2 (p)): (h, tol) and (h, dt) at N = 512 (dt / 2 for the dt
//   factor), each must be evaluated and non-interacting (4.1 (c)). The ramp-rate departure (Q_s at a 1 ms
//   ramp) is printed, not judged.
// Reference independence I1: q, eps_si, V_T, n_i, chi, E_g, N_c and N_v are NiTCAD's material data (input).
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <expected>
#include <numbers>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "NiTCAD/solve/transient.hpp"
#include "NiTCAD/solve/waveform.hpp"
#include "NiTCAD/verify/benchmark/device_data.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/benchmark/tvb32.hpp"
#include "NiTCAD/verify/reference/first_integral.hpp"
#include "NiTCAD/verify/reference/metrics.hpp"
#include "NiTCAD/verify/reference/mos_capacitor.hpp"
#include "NiTCAD/verify/reference/mos_hf.hpp"
#include "NiTCAD/verify/study/study.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

constexpr double depth_cm = 5e-5;
constexpr double acceptors_cm3 = 1e17;
constexpr double t_ox_cm = 1e-6, eps_ox_r = 3.9, phi_m_eV = 4.7, fixed_charge_cm2 = 5e10;
constexpr double grading = 6.2;
constexpr double T = 300.0;
constexpr double tol_newton = 1e-12;
constexpr double tol_dd = 1e-10;  // deep depletion: tol / 100 = 1e-12 stays above the rounding floor at 1 us steps
constexpr int ladder[4] = {512, 256, 128, 64};  // finest first
constexpr int factorial_N = 1024;
constexpr double hf_offsets[6] = {0.2, 0.5, 1.0, 1.5, 2.0, 2.5};
constexpr double frequencies_Hz[3] = {1e5, 316227.76601683791, 1e6};
constexpr std::size_t judged_frequency = 1;
constexpr double plateau = 1e-3;          // eps_HF <= 1e-3 C_ref (16.18 (b))
constexpr double dd_offsets[4] = {1.0, 2.0, 3.0, 4.0};
constexpr double ramp_s = 1e-4, dt_dd = 1e-6;  // fixed steps (16.40)
constexpr double minority_limit = 1e-3;   // q int n dx <= 1e-3 |Q_s| (16.2 (g))
constexpr double tolerance = 1e-3;

ReferenceInfo info() {
    ReferenceInfo r;
    r.id = "TVB-10";
    r.cls = ReferenceClass::E;
    r.independence = Independence::I1;
    r.provenance = Provenance::PD;
    r.citation = "Kingston and Neustadter, J. Appl. Phys. 26, 718 (1955) (the first integral); Brews, J. Appl. Phys. "
                 "45, 1276 (1974) (the redistribution premise, abstract only, PC; the model derived, PD)";
    r.version = "DECISIONS Unit 55 16.2 (g), 16.8, 16.18 (b), 16.19, 16.40; DT-Q17";
    r.equations = {1, 2, 3, 4};
    r.form = "deep depletion: Q = -sgn(t) sqrt(2 q eps V_T N_A G(-t)), V_G = V_FB + V_T t - Q / C_ox; HF: C_s = "
             "-(dQ_s/dpsi_s) at Q_n fixed, the electrons at one quasi-Fermi level";
    r.assumptions = {"Boltzmann statistics", "uniform substrate", "lumped oxide with fixed charge at the interface",
                     "infinite-depth reference (truncation measured)",
                     "deep depletion: no minority carriers (measured: q int n dx / |Q_s|)",
                     "HF: on the plateau (measured: the spread over a decade)"};
    r.derivation = "verify/reference/mos_hf.hpp";
    r.verification = "mos_hf_test: v = 0 against KN, Q_s = Q_p + Q_n by quadrature, the depletion and inversion limits";
    r.second_method = "HF: Q_s along Q_n held fixed by Brent, Richardson central differences; deep depletion: the "
                      "charge density's quadrature";
    return r;
}

std::vector<double> graded(int N) {
    std::vector<double> x;
    const double scale = depth_cm / std::expm1(grading);
    for (int k = 0; k < N; ++k) x.push_back(k == N - 1 ? depth_cm : scale * std::expm1(grading * k / (N - 1)));
    return x;
}

std::expected<device::Device, std::string> make_device(int N) {
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
                                     .regions = {{"silicon", physics::silicon()}},
                                     .node_region = std::vector<device::RegionId>(n, 0),
                                     .donors = std::vector<double>(n, 0.0),
                                     .acceptors = std::vector<double>(n, acceptors_cm3),
                                     .contacts = {std::move(gate),
                                                  {"substrate", device::ContactKind::ohmic, std::move(sub)}}});
    if (!d) return std::unexpected(d.error().message);
    return std::move(*d);
}

solve::BiasOptions options(double tol, solve::Equations eq) {
    solve::BiasOptions o;
    o.equations = eq;
    o.models.srh = true;
    o.models.auger = o.models.radiative = false;
    o.models.bgn = false;
    o.models.fermi_dirac = false;
    o.models.field_mobility = false;
    o.newton.tol_update = tol;
    o.step_control.max_cuts = 30;
    return o;
}

// HF: C at every (V_G, f) and its resolution, at N nodes and tolerance tol.
struct HfLevel {
    std::vector<std::vector<double>> C, res;  // [bias][frequency]
};

std::expected<HfLevel, std::string> hf_level(int N, double tol, double VFB) {
    const auto d = make_device(N);
    if (!d) return std::unexpected(d.error());
    std::vector<std::vector<double>> points;
    for (const double dv : hf_offsets) points.push_back({VFB + dv, 0.0});
    solve::SmallSignalOptions o{.steady = options(tol, solve::Equations::equilibrium_poisson),
                                .frequencies_Hz = std::vector<double>(std::begin(frequencies_Hz),
                                                                      std::end(frequencies_Hz))};
    const auto s = solve::solve_small_signal(*d, points, o);
    if (!s) return std::unexpected("HF, N = " + std::to_string(N) + ": " + s.error().message);
    if (s->stopped) return std::unexpected("HF, N = " + std::to_string(N) + ": " + s->stopped->message);
    HfLevel lv;
    for (std::size_t q = 0; q < s->points.size(); ++q) {
        std::vector<double> c, r;
        for (std::size_t k = 0; k < std::size(frequencies_Hz); ++k) {
            c.push_back(s->capacitance(q, k, 0, 0));
            r.push_back(s->points[q].resolution[k][0] / (2.0 * std::numbers::pi * frequencies_Hz[k]));
        }
        lv.C.push_back(std::move(c));
        lv.res.push_back(std::move(r));
    }
    return lv;
}

// Deep depletion: Q_s and the electrons' charge at each output time, at N nodes, tolerance tol, fixed step dt and ramp time.
struct DdLevel {
    std::vector<double> Q, minority;
    std::size_t steps = 0;
};

// The deep-depletion ramp on device d (N nodes, for the messages).
std::expected<results::Transient, std::string> dd_transient(const device::Device& d, int N, double tol, double dt,
                                                            double ramp, double VFB) {
    const double top = dd_offsets[std::size(dd_offsets) - 1];
    auto gate = solve::Waveform::ramp(VFB, VFB + top, 0.0, ramp);
    if (!gate) return std::unexpected(gate.error().message);
    const std::vector<solve::Waveform> w{*gate, solve::Waveform::constant(0.0)};
    solve::TransientOptions o;
    o.steady = options(tol, solve::Equations::drift_diffusion);
    o.integrator = solve::Integrator::bdf2;
    o.t_end_s = ramp;
    o.adaptive = false;
    o.dt_initial_s = dt;
    o.dt_max_s = dt;
    o.dt_min_s = 1e-6 * dt;
    for (const double dv : dd_offsets) {
        if (dv < top) o.output_times_s.push_back(ramp * dv / top);
    }
    auto s = solve::solve_transient(d, w, o);
    if (!s) return std::unexpected("deep depletion, N = " + std::to_string(N) + ": " + s.error().message);
    if (s->stopped) return std::unexpected("deep depletion, N = " + std::to_string(N) + ": " + s->stopped->message);
    return std::move(*s);
}

std::expected<DdLevel, std::string> dd_level(int N, double tol, double dt, double ramp, double VFB) {
    const auto d = make_device(N);
    if (!d) return std::unexpected(d.error());
    const double top = dd_offsets[std::size(dd_offsets) - 1];
    const auto s = dd_transient(*d, N, tol, dt, ramp, VFB);
    if (!s) return std::unexpected(s.error());
    DdLevel lv;
    lv.steps = s->points.size();
    const auto& x = d->mesh().points();
    for (const double dv : dd_offsets) {
        const double t = ramp * dv / top;
        const auto snap = std::find_if(s->snapshots.begin(), s->snapshots.end(),
                                       [&](const auto& sn) { return std::abs(sn.time_s - t) <= 1e-12 * ramp; });
        if (snap == s->snapshots.end()) return std::unexpected("no snapshot at t = " + shortest(t));
        lv.Q.push_back(-(s->points[snap->point].contact_charge[0] + base::q_C * fixed_charge_cm2));
        double qn = 0.0;  // q int n dx, trapezoids
        for (std::size_t i = 0; i + 1 < x.size(); ++i) {
            qn += 0.5 * (snap->fields.n_cm3[i] + snap->fields.n_cm3[i + 1]) * (x[i + 1][0] - x[i][0]);
        }
        lv.minority.push_back(base::q_C * qn);
    }
    return lv;
}

std::expected<reference::DD, std::string> flat_band(const physics::Semiconductor& si) {
    const double VT = base::thermal_voltage(T), ni = physics::intrinsic_density(si, T);
    return reference::flat_band_voltage({.work_function_eV = phi_m_eV,
                                         .electron_affinity_eV = si.parameters().electron_affinity_eV,
                                         .band_gap_eV = physics::band_gap_eV(si, T),
                                         .conduction_dos_cm3 = physics::conduction_band_dos(si, T),
                                         .valence_dos_cm3 = physics::valence_band_dos(si, T),
                                         .thermal_voltage_V = VT,
                                         .ni_cm3 = ni,
                                         .net_doping_cm3 = -acceptors_cm3,
                                         .fixed_charge_cm2 = fixed_charge_cm2,
                                         .q_C = base::q_C,
                                         .oxide_capacitance_F_per_cm2 = eps_ox_r * base::eps0_F_per_cm / t_ox_cm});
}

}  // namespace

BenchmarkResult run_tvb10(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-10";
    r.title = "MOS capacitor: high-frequency C-V and deep depletion";
    r.info = info();
    if (tier == Tier::fast) {
        r.notes.push_back("full tier only (16.7)");
        return r;
    }
    const physics::Semiconductor si = physics::silicon();
    const double VT = base::thermal_voltage(T), q = base::q_C;
    const double cox = eps_ox_r * base::eps0_F_per_cm / t_ox_cm;
    const double ni = physics::intrinsic_density(si, T);
    const auto vfb = flat_band(si);
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
    const auto tail_of = [&](const reference::DD& u_s) -> std::expected<double, std::string> {
        const std::vector<double> dd{depth_cm};
        const auto t = reference::bulk_profile(bulk, u_s, debye, dd);
        if (!t) return std::unexpected(t.error());
        return (*t)[0].tail.to_double() / std::abs((u_s - bulk.u_b).to_double());
    };
    r.notes.push_back("reference: V_FB " + shortest(VFB) + " V, C_ox " + shortest(cox) + " F/cm^2, n_i " +
                      shortest(ni) + " cm^-3");
    ReferenceInfo a_info = r.info;
    a_info.cls = ReferenceClass::A;
    a_info.approximation_error = "measured: the spread of C_sim over a decade of f on the HF plateau (16.18 (b))";
    // ---- HF (A) ----
    std::vector<reference::MosHfPoint> hf_ref;
    std::vector<double> hf_uref;
    for (const double dv : hf_offsets) {
        const auto p = reference::mos_hf_point(in, VFB + dv, true);
        if (!p) {
            r.run_error = "HF reference at V_FB + " + shortest(dv) + ": " + p.error();
            return r;
        }
        const auto tail = tail_of(p->dc.u_s);
        if (!tail) {
            r.run_error = "HF reference tail: " + tail.error();
            return r;
        }
        const double C = p->capacitance_hf.to_double();
        if (!(p->second_method_difference <= 1e-15 * C)) {  // the reference's self-check (16.40)
            r.run_error = "HF reference self-check: the second method differs by " + shortest(p->second_method_difference) +
                          " F/cm^2 at V_G - V_FB = " + shortest(dv);
            return r;
        }
        hf_uref.push_back(p->second_method_difference + p->quadrature_error * C + *tail * C + 1e-15 * C);
        r.notes.push_back("HF reference, V_G - V_FB = " + shortest(dv) + " V: C_HF " + shortest(C) + ", C_LF " +
                          shortest(p->dc.capacitance.to_double()) + " F/cm^2 (the electrons following, printed), Q_n " +
                          shortest(p->minority_charge_C_per_cm2.to_double()) + " C/cm^2, dv/dt " +
                          shortest(p->dv_dt.to_double()) + "; second method differs by " +
                          shortest(p->second_method_difference) + ", truncation " + shortest(*tail));
        hf_ref.push_back(*p);
    }
    std::vector<HfLevel> hl;
    double hs[4];
    for (int k = 0; k < 4; ++k) {
        const auto s = hf_level(ladder[k], tol_newton, VFB);
        if (!s) {
            r.run_error = s.error();
            return r;
        }
        hl.push_back(*s);
        hs[k] = 1.0 / (ladder[k] - 1);
    }
    const double h_half = 1.0 / (factorial_N - 1);
    std::vector<HfLevel> hfac;  // (h, tol/100), (h/2, tol), (h/2, tol/100)
    for (const auto& [N, tol] : {std::pair{ladder[0], 0.01 * tol_newton}, std::pair{factorial_N, tol_newton},
                                 std::pair{factorial_N, 0.01 * tol_newton}}) {
        const auto s = hf_level(N, tol, VFB);
        if (!s) {
            r.run_error = s.error() + " (factorial)";
            return r;
        }
        hfac.push_back(*s);
    }
    for (std::size_t b = 0; b < std::size(hf_offsets); ++b) {
        const double dv = hf_offsets[b];
        const std::size_t k = judged_frequency;
        const double Cref = hf_ref[b].capacitance_hf.to_double();
        const auto& fineC = hl[0].C[b];
        const double spread = *std::max_element(fineC.begin(), fineC.end()) - *std::min_element(fineC.begin(), fineC.end());
        const double f[4] = {hl[0].C[b][k], hl[1].C[b][k], hl[2].C[b][k], hl[3].C[b][k]};
        const Discretization ud = u_disc(std::span<const double, 4>(hs, 4), std::span<const double, 4>(f, 4), h_half,
                                         hfac[1].C[b][k]);
        const std::vector<double> parts{ud.u ? *ud.u : 0.0, hl[0].res[b][k]};
        const Factorial fc = factorial(f[0] - Cref, hfac[0].C[b][k] - Cref, hfac[1].C[b][k] - Cref,
                                       hfac[2].C[b][k] - Cref, parts);
        const double unum = fc.ok ? fc.u : 0.0;
        const double bound = unum + hf_uref[b] + spread;
        const double metric = std::abs(f[0] - Cref);
        const std::string at = "V_G - V_FB = " + shortest(dv) + " V";
        r.notes.push_back("HF, " + at + ": C_sim " + shortest(f[0]) + " (f = 1e5, 10^5.5, 1e6 Hz: " +
                          shortest(fineC[0]) + ", " + shortest(fineC[1]) + ", " + shortest(fineC[2]) + "), C_ref " +
                          shortest(Cref) + "; eps_HF " + shortest(spread) + ", U_num " + shortest(unum) +
                          " (U_disc: " + ud.how + "; GCI orders " + shortest(ud.order) + ", " +
                          shortest(ud.order_coarse) + ", fine " + shortest(ud.order_fine) + ", resolution " +
                          shortest(hl[0].res[b][k]) + "), U_ref " + shortest(hf_uref[b]) + "; " + fc.note);
        const std::string criterion = "|C_sim - C_ref| <= U_num + U_ref + eps_HF = " + shortest(bound);
        ComparisonRecord cr{"C_HF, " + at, "F/cm^2", dv, f[0], Cref, hf_uref[b], metric, criterion, {}};
        const std::vector<Predicate> pr{{"eps_HF / C_ref <= 1e-3 (on the HF plateau)", spread / Cref, plateau,
                                         spread <= plateau * Cref}};
        cr.comparison = !(ud.u && fc.ok)
                            ? Comparison{Outcome::not_evaluated,
                                         "U_num not estimated (" + (ud.u ? std::string("the (h, tol) factorial") : ud.how) +
                                             ")",
                                         {}}
                            : judge(a_info, pr, metric <= bound, criterion);
        r.comparisons.push_back(std::move(cr));
    }

    // ---- Deep depletion (E) ----
    std::vector<reference::DeepDepletionPoint> dd_ref;
    std::vector<double> dd_uref;
    for (const double dv : dd_offsets) {
        const auto p = reference::mos_deep_depletion_point(in, VFB + dv);
        if (!p) {
            r.run_error = "deep-depletion reference at V_FB + " + shortest(dv) + ": " + p.error();
            return r;
        }
        // The truncation: beyond the depletion edge W (the depletion approximation's, sqrt(2 eps (psi - V_T) / (q N_A)))
        // the holes' profile decays as e^{-x / L_D} with an amplitude of order one, so the tail at x = d is bounded by
        // 4 e^{-(d - W) / L_D} (the equilibrium profile routine cannot take a deep-depletion band bending: it would
        // include the inversion layer).
        const double psi = p->t_s.to_double() * VT;
        const double W = std::sqrt(2.0 * in.eps_F_per_cm * std::max(psi - VT, 0.0) / (q * acceptors_cm3));
        const std::expected<double, std::string> tail = 4.0 * std::exp(-(depth_cm - W) / debye.to_double());
        const double Q = p->charge_C_per_cm2.to_double();
        dd_uref.push_back(p->second_method_difference / std::abs(Q) + *tail + 1e-15);
        r.notes.push_back("deep-depletion reference, V_G - V_FB = " + shortest(dv) + " V: Q " + shortest(Q) +
                          " C/cm^2, band bending " + shortest(p->t_s.to_double() * VT) + " V; quadrature differs by " +
                          shortest(p->second_method_difference) + ", truncation " + shortest(*tail));
        dd_ref.push_back(*p);
    }
    std::vector<DdLevel> dl;
    for (const int N : ladder) {
        const auto s = dd_level(N, tol_dd, dt_dd, ramp_s, VFB);
        if (!s) {
            r.run_error = s.error();
            return r;
        }
        dl.push_back(*s);
    }
    // The factorials at N = 512: (h, tol) and (h, dt), dt / 2 the finer step.
    DdLevel f_tol, f_h, f_h_tol, f_dt, f_h_dt;
    {
        const std::pair<DdLevel*, std::expected<DdLevel, std::string>> runs[] = {
            {&f_tol, dd_level(ladder[0], 0.01 * tol_dd, dt_dd, ramp_s, VFB)},
            {&f_h, dd_level(factorial_N, tol_dd, dt_dd, ramp_s, VFB)},
            {&f_h_tol, dd_level(factorial_N, 0.01 * tol_dd, dt_dd, ramp_s, VFB)},
            {&f_dt, dd_level(ladder[0], tol_dd, 0.5 * dt_dd, ramp_s, VFB)},
            {&f_h_dt, dd_level(factorial_N, tol_dd, 0.5 * dt_dd, ramp_s, VFB)}};
        for (const auto& [dst, s] : runs) {
            if (!s) {
                r.run_error = s.error() + " (factorial)";
                return r;
            }
            *dst = *s;
        }
    }
    const auto slow = dd_level(ladder[0], tol_dd, 10.0 * dt_dd, 10.0 * ramp_s, VFB);
    if (!slow) {
        r.run_error = slow.error() + " (1 ms ramp)";
        return r;
    }
    r.notes.push_back("deep depletion: accepted steps per run, N = 512 ... 64: " + std::to_string(dl[0].steps) + ", " +
                      std::to_string(dl[1].steps) + ", " + std::to_string(dl[2].steps) + ", " +
                      std::to_string(dl[3].steps));
    for (std::size_t b = 0; b < std::size(dd_offsets); ++b) {
        const double dv = dd_offsets[b];
        const double Qref = dd_ref[b].charge_C_per_cm2.to_double();
        const auto e = [&](const DdLevel& l) { return std::abs(l.Q[b] - Qref) / std::abs(Qref); };
        const auto se = [&](const DdLevel& l) { return (l.Q[b] - Qref) / std::abs(Qref); };
        std::vector<double> h, es, us;
        for (int k = 0; k < 4; ++k) {
            h.push_back(1.0 / (ladder[k] - 1));
            es.push_back(e(dl[static_cast<std::size_t>(k)]));
            us.push_back(std::abs(dl[0].Q[b] - f_tol.Q[b]) / std::abs(Qref));
        }
        const study::OrderStudy st = study::order_study(h, es, us);
        const study::OrderVerdict v = study::judge_order(st, {.p_th = 2.0});
        std::string orders;
        for (const double p : st.orders) orders += " " + shortest(p);
        const double u_disc3 = st.orders.empty() ? 0.0
                                                 : reference::gci_safety_factor_three_grids * std::abs(es[1] - es[0]) /
                                                       (std::pow(h[1] / h[0], st.observed) - 1.0);
        const double u_dt = std::abs(dl[0].Q[b] - f_dt.Q[b]) / std::abs(Qref);
        const std::vector<double> parts{u_disc3, us[0], u_dt, dd_uref[b]};
        const Factorial fac_tol = factorial(se(dl[0]), se(f_tol), se(f_h), se(f_h_tol), parts);
        Factorial fac_dt = factorial(se(dl[0]), se(f_dt), se(f_h), se(f_h_dt), parts);
        fac_dt.note.replace(0, std::string("factorial (h, tol)").size(), "factorial (h, dt)");
        const double minority = dl[0].minority[b] / std::abs(dl[0].Q[b]);
        const double rate = std::abs(slow->Q[b] - dl[0].Q[b]) / std::abs(Qref);
        const std::string at = "V_G - V_FB = " + shortest(dv) + " V";
        r.notes.push_back("deep depletion, " + at + ": Q_ref " + shortest(Qref) + ", errors (N = 512 ... 64)" +
                          [&] {
                              std::string s;
                              for (const double x : es) s += " " + shortest(x);
                              return s;
                          }() +
                          "; orders" + orders + ", spread " + shortest(st.order_spread) + "; electrons' charge / |Q_s| " +
                          shortest(minority) + "; U_iter " + shortest(us[0]) + ", U_dt " + shortest(u_dt) +
                          ", U_ref " + shortest(dd_uref[b]) + "; ramp-rate departure (1 ms) " + shortest(rate) +
                          "; " + fac_tol.note + "; " + fac_dt.note);
        const std::vector<Predicate> pr{
            {"q int n dx / |Q_s| <= 1e-3 (no minority carriers)", minority, minority_limit, minority <= minority_limit}};
        const std::string band = "p-hat in [1.8, 2.3] (p_th 2)", accuracy = "relative <= 1e-3";
        ComparisonRecord co{"Q_s order (deep depletion), " + at, "1", dv, st.observed, 2.0, dd_uref[b], st.observed,
                            band, {}};
        ComparisonRecord ca{"Q_s finest (deep depletion, N = 512), " + at, "C/cm^2", dv, dl[0].Q[b], Qref,
                            dd_uref[b] * std::abs(Qref), es[0], accuracy, {}};
        co.comparison = v.judged || !pr[0].holds ? judge(r.info, pr, v.pass, band)
                                                 : Comparison{Outcome::not_evaluated, v.reason, {}};
        ca.comparison = judge(r.info, pr, es[0] <= tolerance, accuracy);
        for (ComparisonRecord* c : {&co, &ca}) {
            if (c->comparison.outcome == Outcome::not_applicable) continue;
            block_if_interacting(*c, fac_tol);
            block_if_interacting(*c, fac_dt);
        }
        r.comparisons.push_back(std::move(co));
        r.comparisons.push_back(std::move(ca));
    }
    apply_known_findings(r, known_findings());
    return r;
}

// TVB-32 on this benchmark's device (16.44): N = 64; the HF run (equilibrium Poisson) at zero bias, V_FB + 0.5 V and
// V_FB + 2 V at its three frequencies; the equilibrium; the deep-depletion ramp (C6).
Tvb32Run tvb10_consistency() {
    Tvb32Run run;
    const physics::Semiconductor si = physics::silicon();
    const auto vfb = flat_band(si);
    if (!vfb) return {{}, "TVB-10: reference: " + vfb.error()};
    const double VFB = vfb->to_double();
    const auto d = make_device(64);
    if (!d) return {{}, "TVB-10: " + d.error()};
    Tvb32Collector col{run.entries, "TVB-10 MOS-C N = 64", d->contacts().size()};
    const solve::BiasOptions hf = options(1e-12, solve::Equations::equilibrium_poisson);
    const solve::SmallSignalOptions o{
        .steady = hf, .frequencies_Hz = std::vector<double>(std::begin(frequencies_Hz), std::end(frequencies_Hz))};
    const auto s = solve::solve_small_signal(
        *d, std::vector<std::vector<double>>{{0.0, 0.0}, {VFB + 0.5, 0.0}, {VFB + 2.0, 0.0}}, o);
    if (!s || s->stopped) return {{}, "TVB-10 HF, N = 64: " + (s ? s->stopped->message : s.error().message)};
    col.ac(*s, d->mesh(), device_data(*d, hf.models), "HF, ");
    std::string error;
    if (!col.equilibrium(*d, hf.models, error)) return {std::move(run.entries), error};
    const auto t = dd_transient(*d, 64, 1e-12, dt_dd, ramp_s, VFB);
    if (!t) return {std::move(run.entries), "TVB-10: " + t.error()};
    const solve::BiasOptions ddo = options(1e-12, solve::Equations::drift_diffusion);
    col.transient(*t, d->mesh(), device_data(*d, ddo.models), "deep-depletion ramp");
    return run;
}

}  // namespace NiTCAD::verify::benchmark
