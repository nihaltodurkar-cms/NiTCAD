// TVB-22, an interface trap level in a MOS capacitor: the small-signal admittance against the single-level circuit
// (DECISIONS Unit 55 6, 17.2 (e)-(g), 17.4, 17.13 (a); criteria 17.8 and the predicate limits of 17.14, frozen
// before code).
//
// Code under test: the interface edges of a meshed oxide (the level's charge at the interface potential), the
// quasi-static operating point, and solve_small_signal with the traps' dynamics T(s).
// Construction, fixed before any run:
//   - 1D: an electrode (metal, phi_m = 4.1 eV) on 10 nm of SiO2 over 1 um of n-type silicon, N_D = 1e16 cm^-3, an
//     ohmic contact at the back; 300 K; Boltzmann, no band-gap narrowing, no recombination, constant mobility (the
//     doping dependence off) mu_n = 1000, mu_p = 400 cm^2/Vs;
//   - silicon: eps_r 11.7, E_g 1.12 eV (no Varshni term), N_c 2.86e19, N_v 3.10e19, chi 4.05 eV (the spec; R-7);
//     oxide: NiTCAD's SiO2, whose eps_r the spec states as 3.9;
//   - the level: N_t = 1e10 cm^-2 at E_t - E_i = 0.20 eV, sigma_n = sigma_p = 1e-16 cm^2, v_th = 1e7 cm/s (c = 1e-9
//     cm^3/s); acceptor-like and donor-like in separate runs;
//   - the gate at V_FB - 0.20 V and V_FB - 0.30 V (depletion), V_FB the trap-free flat band of the reference;
//   - mesh: silicon nodes h_1 / 2 + L sinh(beta k / (N - 1)) / sinh(beta), beta = 4, k = 0 ... N - 1, h_1 the first
//     silicon spacing; N / 8 + 1 oxide nodes uniform on [-t_ox, -h_1 / 2]. The oxide-silicon edge has length h_1 and
//     its midpoint, where the interface lies (interface_edges.hpp), is x = 0 at every level (17.2 (e)).
// The operating point is quasi-static (Equations::equilibrium_poisson); the small-signal system is drift-diffusion.
// Frequencies: w tau = 1e-2, 1e-1, 0.3, 1, 3, 10, 100 and, for the peak, 10^-0.1, 10^-0.05, 10^0.05, 10^0.1, tau
// the reference's.
// Reference (verify/reference/trap_admittance.hpp): the operating point with the level's charge in the gate
// equation, and 1 / Y = 1 / (i w C_ox) + 1 / (i w C_s + i w C_T / (1 + i w tau)) + Z_b; cross-checked at run time
// against the direct 3 x 3 solve. Inputs: the spec above and NiTCAD's q, eps_0 and V_T.
// Reference uncertainty: the neutral length (the distance from the interface to the contact, with and without the
// depletion width), |Y(L) - Y(L - W)|, and the two methods' difference.
// Applicability (17.14), each computed by the reference from the spec: P1 <= 1e-6, P2 <= 1e-5, P3 <= 1e-6, and
// w tau_d <= 1e-2 at each frequency; else NOT APPLICABLE with the value.
// Criteria (E, 17.8), per trap type and bias:
//   - (a) p-hat in [1.8, 2.3] of the errors of G_p / w at w tau = 1 (relative to C_T / 2) and of C_p at w tau = 1e-2
//     (relative), on N = 64 ... 512; G_p / w and C_p from Y_s = 1 / (1 / Y - 1 / (i w C_ox)) with the spec's C_ox
//     (the benchmark's own removal; analysis::conductance_peak is not used, 17.2 (g));
//   - (b) at N = 512: |Y - Y_ref| <= 1e-3 w C_T at each of the seven frequencies; the peak of G_p / w, by the vertex
//     of the parabola in ln w through the five samples around w tau = 1, within 1e-3 relative of the reference's
//     by the same procedure.
//   A comparison whose reference uncertainty or admittance resolution exceeds a tenth of its tolerance is
//   NOT EVALUATED.
// The (h, tol) factorial (16.2 (p)) at N = 512 on |Y - Y_ref| / (w C_T) at w tau = 1: tol = 1e-10 and 1e-12
// (the operating point is the Poisson solve), N = 512 and 1024.
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <expected>
#include <numbers>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/insulator.hpp"
#include "NiTCAD/physics/interface_traps.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/metrics.hpp"
#include "NiTCAD/verify/reference/mos_capacitor.hpp"
#include "NiTCAD/verify/reference/trap_admittance.hpp"
#include "NiTCAD/verify/study/study.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

using Complex = std::complex<double>;

constexpr double T = 300.0;
constexpr double length_cm = 1e-4, t_ox_cm = 1e-6, beta = 4.0;
constexpr double donors_cm3 = 1e16, phi_m_eV = 4.1, eps_ox_r = 3.9;
constexpr double eps_si_r = 11.7, gap_eV = 1.12, Nc = 2.86e19, Nv = 3.10e19, chi_eV = 4.05;
constexpr double mu_n = 1000.0, mu_p = 400.0;
constexpr double N_t = 1e10, level_eV = 0.20, sigma_cm2 = 1e-16, v_th = 1e7;
// 1e-10, and 1e-12 in the factorial: 1e-14 is not reachable at N = 512 (found at the first run, before any result).
constexpr double tolerance = 1e-3, tol_newton = 1e-10;
constexpr double limit_p1 = 1e-6, limit_p2 = 1e-5, limit_p3 = 1e-6, limit_relaxation = 1e-2;
constexpr std::array<double, 2> offsets_V{-0.20, -0.30};
// w tau: the seven judged frequencies, then the four more around the peak.
constexpr std::array<double, 7> judged{1e-2, 1e-1, 0.3, 1.0, 3.0, 10.0, 100.0};
constexpr std::array<double, 5> around{-0.1, -0.05, 0.0, 0.05, 0.1};  // log10(w tau)

ReferenceInfo info() {
    ReferenceInfo r;
    r.id = "TVB-22";
    r.cls = ReferenceClass::E;
    r.independence = Independence::I1;
    r.provenance = Provenance::PD;
    r.citation = "The single-level interface-state circuit, derived from the model's rate equation (DECISIONS Unit 55 "
                 "17.13 (a)); the series C_s-R_s branch as A. Goetzberger, Interface states in Si-SiO2 interfaces "
                 "(1974), Fig. 1B (PA); the conductance peak as Nicollian and Goetzberger, BSTJ 46, 1055 (1967), 4.3";
    r.version = "DECISIONS Unit 55 17.3, 17.8, 17.14";
    r.equations = {1, 2, 3, 4};
    r.form = "1 / Y = 1 / (i w C_ox) + 1 / (i w C_s + i w C_T / (1 + i w tau)) + R_b / (1 + i w tau_d), "
             "C_T = q N_t f_0 (1 - f_0) / V_T, tau = f_0 / (c_n n_0)";
    r.assumptions = {"quasi-static operating point", "Boltzmann statistics", "majority carriers in equilibrium at the "
                     "interface (P2)", "minority carriers not exchanging with the level (P1) nor storing charge (P3)",
                     "frequency far below dielectric relaxation"};
    r.derivation = "verify/reference/trap_admittance.hpp";
    r.verification = "E: orders of G_p / w and C_p; finest-level Y at seven frequencies and the peak's place";
    r.second_method = "direct solve of the linearized surface equations (3 x 3 complex), at run time";
    return r;
}

struct Axis {
    std::vector<double> x;
    double h1;  // the first silicon spacing: the oxide-silicon edge's length
};

Axis axis(int N) {
    std::vector<double> s;
    const double sh = std::sinh(beta);
    for (int k = 0; k < N; ++k) s.push_back(length_cm * std::sinh(beta * k / (N - 1)) / sh);
    Axis a;
    a.h1 = s[1];
    const int M = N / 8;
    for (int k = 0; k <= M; ++k) a.x.push_back(-t_ox_cm + (t_ox_cm - 0.5 * a.h1) * k / M);
    for (const double v : s) a.x.push_back(v + 0.5 * a.h1);
    return a;
}

physics::SemiconductorParameters silicon() {
    physics::SemiconductorParameters p = physics::legacy_silicon_parameters;
    p.eps_r = eps_si_r;
    p.Eg0_eV = gap_eV;
    p.varshni_alpha_eV_per_K = 0.0;
    p.varshni_beta_K = 0.0;
    p.Nc300 = Nc;
    p.Nv300 = Nv;
    p.electron_affinity_eV = chi_eV;
    p.electron_mobility.mu_max = mu_n;
    p.hole_mobility.mu_max = mu_p;
    return p;
}

std::expected<device::Device, std::string> capacitor(int N, bool donor) {
    const Axis a = axis(N);
    auto m = mesh::make_tensor_grid(a.x);
    if (!m) return std::unexpected(m.error().message);
    const std::size_t n = m->node_count();
    std::vector<device::RegionId> region(n, 0);
    std::vector<double> donors(n, donors_cm3);
    for (std::size_t i = 0; i < n; ++i) {
        if (m->points()[i][0] < 0.0) {
            region[i] = 1;
            donors[i] = 0.0;
        }
    }
    const auto si = physics::Semiconductor::create(silicon());
    if (!si) return std::unexpected("material: " + si.error().message);
    device::Contact gate{"gate", device::ContactKind::electrode, m->find_boundary("x_min")->nodes};
    gate.electrode = {.kind = device::GateElectrode::metal, .work_function_eV = phi_m_eV};
    auto back = m->find_boundary("x_max")->nodes;
    device::Interface f{"silicon", "oxide"};
    f.traps.levels = {{.type = donor ? physics::TrapType::donor : physics::TrapType::acceptor,
                       .density_cm2 = N_t,
                       .energy_eV = level_eV,
                       .sigma_n_cm2 = sigma_cm2,
                       .sigma_p_cm2 = sigma_cm2}};
    f.traps.thermal_velocity_n_cm_s = v_th;
    f.traps.thermal_velocity_p_cm_s = v_th;
    auto d = device::Device::create({.mesh = std::move(*m),
                                     .temperature_K = T,
                                     .regions = {{"silicon", *si}, {"oxide", physics::silicon_dioxide()}},
                                     .node_region = std::move(region),
                                     .donors = std::move(donors),
                                     .acceptors = std::vector<double>(n, 0.0),
                                     .contacts = {std::move(gate),
                                                  {"substrate", device::ContactKind::ohmic, std::move(back)}},
                                     .interfaces = {std::move(f)}});
    if (!d) return std::unexpected(d.error().message);
    return std::move(*d);
}

solve::BiasOptions options(double tol) {
    solve::BiasOptions o;
    o.equations = solve::Equations::equilibrium_poisson;
    o.models.srh = o.models.auger = o.models.radiative = false;
    o.models.bgn = false;
    o.models.fermi_dirac = false;
    o.models.doping_mobility = false;
    o.models.hole_mobility = assemble::HoleMobility::caughey_thomas;
    o.newton.tol_update = tol;
    return o;
}

// The vertex of the parabola through (x_k, g_k), by least squares: x of the maximum.
double vertex(const std::array<double, 5>& x, const std::array<double, 5>& g) {
    // Normal equations for g = a + b x + c x^2 with centred x (the abscissae are symmetric about x[2]).
    double s2 = 0.0, s4 = 0.0, sg = 0.0, sxg = 0.0, sx2g = 0.0;
    for (std::size_t k = 0; k < 5; ++k) {
        const double d = x[k] - x[2];
        s2 += d * d;
        s4 += d * d * d * d;
        sg += g[k];
        sxg += d * g[k];
        sx2g += d * d * g[k];
    }
    const double b = sxg / s2;
    const double c = (5.0 * sx2g - s2 * sg) / (5.0 * s4 - s2 * s2);
    return x[2] - b / (2.0 * c);
}

Complex behind_oxide(Complex y, double w, double cox) { return 1.0 / (1.0 / y - 1.0 / Complex(0.0, w * cox)); }

// One level of the ladder: Y at the eleven frequencies (the seven judged first), its resolution, U_iter.
struct Level {
    int N = 0;
    double h = 0.0;
    std::vector<Complex> Y;
    std::vector<double> resolution;
    double u_iter = 0.0;
};

}  // namespace

BenchmarkResult run_tvb22(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-22";
    r.title = "Interface trap level: the small-signal admittance";
    r.info = info();
    if (tier < Tier::full) return r;
    const double VT = base::thermal_voltage(T);
    const double eps_si = eps_si_r * base::eps0_F_per_cm;
    const double cox = eps_ox_r * base::eps0_F_per_cm / t_ox_cm;
    const double ni = std::sqrt(Nc * Nv) * std::exp(-gap_eV / (2.0 * VT));
    const auto vfb = reference::flat_band_voltage({.work_function_eV = phi_m_eV,
                                                   .electron_affinity_eV = chi_eV,
                                                   .band_gap_eV = gap_eV,
                                                   .conduction_dos_cm3 = Nc,
                                                   .valence_dos_cm3 = Nv,
                                                   .thermal_voltage_V = VT,
                                                   .ni_cm3 = ni,
                                                   .net_doping_cm3 = donors_cm3,
                                                   .fixed_charge_cm2 = 0.0,
                                                   .q_C = base::q_C,
                                                   .oxide_capacitance_F_per_cm2 = cox});
    if (!vfb) {
        r.run_error = "reference: " + vfb.error();
        return r;
    }
    const double VFB = vfb->to_double();
    r.notes.push_back("reference: V_FB " + shortest(VFB) + " V, n_i " + shortest(ni) + " cm^-3, C_ox " +
                      shortest(cox) + " F/cm^2");
    const auto record = [&](const std::string& quantity, double parameter, double sim, double refv, double uref,
                            double metric, const std::string& crit, Comparison cmp) -> ComparisonRecord& {
        ComparisonRecord rec{quantity, "1", parameter, sim, refv, uref, metric, crit, {}};
        rec.comparison = std::move(cmp);
        r.comparisons.push_back(std::move(rec));
        return r.comparisons.back();
    };
    for (const bool donor : {false, true}) {
        for (const double offset : offsets_V) {
            const std::string at = std::string(donor ? "donor" : "acceptor") + ", V_FB " + shortest(offset) + " V";
            const reference::TrapAdmittanceInput in{.q_C = base::q_C,
                                                    .thermal_voltage_V = VT,
                                                    .eps_F_per_cm = eps_si,
                                                    .ni_cm3 = ni,
                                                    .net_doping_cm3 = donors_cm3,
                                                    .oxide_capacitance_F_per_cm2 = cox,
                                                    .flat_band_V = VFB,
                                                    .gate_V = VFB + offset,
                                                    .donor = donor,
                                                    .density_cm2 = N_t,
                                                    .energy_eV = level_eV,
                                                    .capture_n_cm3_s = sigma_cm2 * v_th,
                                                    .capture_p_cm3_s = sigma_cm2 * v_th,
                                                    .mobility_cm2_Vs = mu_n};
            const auto op = reference::trap_operating_point(in);
            if (!op) {
                r.run_error = at + ": reference: " + op.error();
                return r;
            }
            const double tau = op->time_constant_s, CT = op->trap_capacitance_F_per_cm2;
            // The benchmark's own removal of C_ox, checked on the circuit without a substrate: G_p / w at w tau = 1
            // must be C_T / 2.
            {
                const Complex ys = behind_oxide(reference::trap_admittance(in, *op, 1.0 / tau, 0.0), 1.0 / tau, cox);
                if (!(std::abs(ys.real() * tau / (0.5 * CT) - 1.0) <= 1e-9)) {
                    r.run_error = at + ": the removal of C_ox does not return C_T / 2 on the reference circuit";
                    return r;
                }
            }
            std::vector<double> wt(judged.begin(), judged.end());
            for (const double e : around) {
                if (e != 0.0) wt.push_back(std::pow(10.0, e));
            }
            std::vector<double> freq;
            for (const double v : wt) freq.push_back(v / (2.0 * std::numbers::pi * tau));
            r.notes.push_back(at + ": reference psi_s " + shortest(op->surface_potential_V) + " V, f_0 " +
                              shortest(op->occupancy) + ", C_s " + shortest(op->semiconductor_capacitance_F_per_cm2) +
                              ", C_T " + shortest(CT) + " F/cm^2, tau " + shortest(tau) + " s, W " +
                              shortest(op->depletion_width_cm) + " cm, E_s " + shortest(op->surface_field_V_per_cm) +
                              " V/cm; P1 " + shortest(op->p1_minority_exchange) + ", P2 " +
                              shortest(op->p2_fermi_drop) + ", P3 " + shortest(op->p3_minority_charge) +
                              ", tau_d " + shortest(op->relaxation_time_s) + " s");
            // The reference at a level's neutral length, its uncertainty and the I3 cross-check.
            struct Ref {
                std::vector<Complex> Y;
                std::vector<double> u;  // |Y(L) - Y(L - W)| + the two methods' difference
            };
            const auto reference_at = [&](int N) -> std::expected<Ref, std::string> {
                const double total = length_cm + 0.5 * axis(N).h1;
                Ref out;
                for (const double v : wt) {
                    const double w = v / tau;
                    const Complex a = reference::trap_admittance(in, *op, w, total - op->depletion_width_cm);
                    const Complex b = reference::trap_admittance(in, *op, w, total);
                    const auto c = reference::trap_admittance_direct(in, *op, w, total - op->depletion_width_cm);
                    if (!c) return std::unexpected(at + ": " + c.error());
                    out.Y.push_back(a);
                    out.u.push_back(std::abs(a - b) + std::abs(a - *c));
                }
                return out;
            };
            const auto run_level = [&](int N, double tol) -> std::expected<Level, std::string> {
                const auto d = capacitor(N, donor);
                if (!d) return std::unexpected(at + ", N = " + std::to_string(N) + ": " + d.error());
                solve::SmallSignalOptions o{.steady = options(tol), .frequencies_Hz = freq};
                const auto s = solve::solve_small_signal(*d, std::vector<std::vector<double>>{{in.gate_V, 0.0}}, o);
                if (!s || s->stopped) {
                    return std::unexpected(at + ", N = " + std::to_string(N) + ", tol " + shortest(tol) + ": " +
                                           (s ? s->stopped->message : s.error().message));
                }
                Level lv;
                lv.N = N;
                lv.h = length_cm / (N - 1);
                const auto& it = s->points[0].dc.convergence.iterations;
                lv.u_iter = it.empty() ? 0.0 : it.back().update;
                for (std::size_t k = 0; k < freq.size(); ++k) {
                    lv.Y.push_back(s->admittance(0, k, 0, 0));
                    lv.resolution.push_back(s->points[0].resolution[k][0]);
                    lv.u_iter = std::max(lv.u_iter, s->points[0].backward_error[k]);
                }
                return lv;
            };
            // The judged quantities of a level against the reference at that level's neutral length.
            struct Errors {
                std::array<double, 7> y{};  // |Y - Y_ref| / (w C_T)
                double gp = 0.0, cp = 0.0;  // G_p / w at w tau = 1 (relative to C_T / 2); C_p at w tau = 1e-2
                double peak = 0.0;          // the vertex's place, relative
            };
            const auto errors = [&](const Level& lv, const Ref& ref) {
                Errors e;
                for (std::size_t k = 0; k < judged.size(); ++k) {
                    e.y[k] = std::abs(lv.Y[k] - ref.Y[k]) / (wt[k] / tau * CT);
                }
                const auto branch = [&](const std::vector<Complex>& Y, std::size_t k) {
                    return behind_oxide(Y[k], wt[k] / tau, cox) / (wt[k] / tau);
                };
                e.gp = std::abs(branch(lv.Y, 3).real() - branch(ref.Y, 3).real()) / (0.5 * CT);
                e.cp = std::abs(branch(lv.Y, 0).imag() - branch(ref.Y, 0).imag()) / branch(ref.Y, 0).imag();
                // The five samples around w tau = 1: indices 7, 8, 3, 9, 10 (10^-0.1, 10^-0.05, 1, 10^0.05, 10^0.1).
                const std::array<std::size_t, 5> idx{7, 8, 3, 9, 10};
                std::array<double, 5> x{}, gs{}, gr{};
                for (std::size_t j = 0; j < 5; ++j) {
                    x[j] = std::log(wt[idx[j]]);
                    gs[j] = branch(lv.Y, idx[j]).real();
                    gr[j] = branch(ref.Y, idx[j]).real();
                }
                e.peak = std::abs(std::exp(vertex(x, gs) - vertex(x, gr)) - 1.0);
                return e;
            };
            std::vector<Level> levels;  // finest first
            std::vector<Errors> errs;
            Ref fine_ref;
            for (const int N : {512, 256, 128, 64}) {
                const auto lv = run_level(N, tol_newton);
                const auto ref = reference_at(N);
                if (!lv || !ref) {
                    r.run_error = !lv ? lv.error() : ref.error();
                    return r;
                }
                if (N == 512) fine_ref = *ref;
                errs.push_back(errors(*lv, *ref));
                levels.push_back(*lv);
                r.notes.push_back(at + ", N = " + std::to_string(N) + ": |Y - Y_ref| / (w C_T) at w tau = 1 " +
                                  shortest(errs.back().y[3]) + ", G_p / w error " + shortest(errs.back().gp) +
                                  ", C_p error " + shortest(errs.back().cp) + ", peak place " +
                                  shortest(errs.back().peak) + ", U_iter " + shortest(lv->u_iter));
            }
            // Applicability.
            std::string why_na;
            if (op->p1_minority_exchange > limit_p1) why_na = "P1 " + shortest(op->p1_minority_exchange) + " above 1e-6";
            if (op->p2_fermi_drop > limit_p2) why_na = "P2 " + shortest(op->p2_fermi_drop) + " above 1e-5";
            if (op->p3_minority_charge > limit_p3) why_na = "P3 " + shortest(op->p3_minority_charge) + " above 1e-6";
            // The factorial at N = 512 on the w tau = 1 error.
            std::vector<double> hs, eg, ec, us;
            for (std::size_t l = 0; l < levels.size(); ++l) {
                hs.push_back(levels[l].h);
                eg.push_back(errs[l].gp);
                ec.push_back(errs[l].cp);
                us.push_back(levels[l].u_iter);
            }
            const study::OrderStudy sg = study::order_study(hs, eg, us), sc = study::order_study(hs, ec, us);
            const study::OrderVerdict vg = study::judge_order(sg, {.p_th = 2.0}), vc = study::judge_order(sc, {.p_th = 2.0});
            const auto orders = [](const study::OrderStudy& st) {
                std::string out;
                for (const double p : st.orders) out += " " + shortest(p);
                return out + "; spread " + shortest(st.order_spread);
            };
            r.notes.push_back(at + ": G_p / w orders (finest pair first):" + orders(sg) + "; C_p orders:" + orders(sc));
            std::array<double, 3> fac_e{};
            {
                std::size_t k = 0;
                for (const auto& [N, tol] : {std::pair{512, 0.01 * tol_newton}, std::pair{1024, tol_newton},
                                             std::pair{1024, 0.01 * tol_newton}}) {
                    const auto lv = run_level(N, tol);
                    const auto ref = reference_at(N);
                    if (!lv || !ref) {
                        r.run_error = (!lv ? lv.error() : ref.error()) + " (factorial)";
                        return r;
                    }
                    fac_e[k++] = errors(*lv, *ref).y[3];
                }
            }
            const double w1 = 1.0 / tau;
            const double u_disc = reference::gci_safety_factor_three_grids * std::abs(errs[1].y[3] - errs[0].y[3]) /
                                  (std::pow(levels[1].h / levels[0].h, 2.0) - 1.0);
            const std::vector<double> parts{u_disc, levels.front().resolution[3] / (w1 * CT),
                                            fine_ref.u[3] / (w1 * CT)};
            const Factorial fac = factorial(errs.front().y[3], fac_e[0], fac_e[1], fac_e[2], parts);
            r.notes.push_back(at + ": " + fac.note + " (U_disc " + shortest(u_disc) + ")");
            const Level& fin = levels.front();
            const Errors& fe = errs.front();
            // A comparison at frequency index k with tolerance tol on a metric already divided by its scale.
            const auto judged_at = [&](ComparisonRecord& rec, bool pass, const std::string& crit, std::size_t k,
                                       double tol_abs_S) {
                const double w = wt[k] / tau;
                if (!why_na.empty()) {
                    rec.comparison = Comparison{Outcome::not_applicable, why_na, {}};
                } else if (w * op->relaxation_time_s > limit_relaxation) {
                    rec.comparison = Comparison{Outcome::not_applicable,
                                                "w tau_d " + shortest(w * op->relaxation_time_s) + " above 1e-2", {}};
                } else if (fine_ref.u[k] > 0.1 * tol_abs_S) {
                    rec.comparison = Comparison{Outcome::not_evaluated,
                                                "U_ref " + shortest(fine_ref.u[k]) + " S/cm^2 above a tenth of the tolerance", {}};
                } else if (fin.resolution[k] > 0.1 * tol_abs_S) {
                    rec.comparison = Comparison{Outcome::not_evaluated,
                                                "admittance resolution " + shortest(fin.resolution[k]) +
                                                    " S/cm^2 above a tenth of the tolerance", {}};
                } else {
                    rec.comparison = judge(r.info, {}, pass, crit);
                    block_if_interacting(rec, fac);
                }
            };
            const std::string band = "p-hat in [1.8, 2.3] (p_th 2)";
            {
                ComparisonRecord& a = record("G_p / w order (w tau = 1), " + at, offset, sg.observed, 2.0, 0.0,
                                             sg.observed, band, {});
                if (vg.judged) {
                    judged_at(a, vg.pass, band, 3, tolerance * w1 * CT);
                } else {
                    a.comparison = Comparison{Outcome::not_evaluated, vg.reason, {}};
                }
                ComparisonRecord& b = record("C_p order (w tau = 1e-2), " + at, offset, sc.observed, 2.0, 0.0,
                                             sc.observed, band, {});
                if (vc.judged) {
                    judged_at(b, vc.pass, band, 0, tolerance * wt[0] / tau * CT);
                } else {
                    b.comparison = Comparison{Outcome::not_evaluated, vc.reason, {}};
                }
            }
            for (std::size_t k = 0; k < judged.size(); ++k) {
                const std::string crit = "|Y - Y_ref| <= 1e-3 w C_T";
                ComparisonRecord& rec = record("Y finest (N = 512), w tau = " + shortest(wt[k]) + ", " + at, wt[k],
                                               std::abs(fin.Y[k]), std::abs(fine_ref.Y[k]), fine_ref.u[k],
                                               fe.y[k] / tolerance, crit, {});
                judged_at(rec, fe.y[k] <= tolerance, crit, k, tolerance * wt[k] / tau * CT);
            }
            {
                const std::string crit = "peak place relative <= 1e-3";
                ComparisonRecord& rec = record("G_p / w peak place finest (N = 512), " + at, offset, fe.peak, 0.0, 0.0,
                                               fe.peak, crit, {});
                judged_at(rec, fe.peak <= tolerance, crit, 3, tolerance * w1 * CT);
            }
        }
    }
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
