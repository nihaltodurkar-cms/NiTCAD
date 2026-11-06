// TVB-18, the diode's small-signal admittance at forward bias (DECISIONS Unit 55 §6, 16.2 (l), 16.4; A-class criteria
// 16.18 (g), 16.19, frozen; DT-Q16, DT-Q26 A, DT-Q28 B, DT-Q32 A; construction and DT-Q29 A, DT-Q30 A, DT-Q31 A,
// 16.34-16.37).
//
// Code under test: the drift-diffusion operating point and solve_small_signal on TVB-16's (ii) short diode: abrupt
// 1D silicon, N_A = 1e16 cm^-3 (x < 0, the anode at x_min) and N_D = 2e16 cm^-3, 5 um per side, ohmic contacts,
// 300 K, Boltzmann, no band-gap narrowing, no field mobility, Caughey-Thomas mobility for both carriers, SRH on
// (Scharfetter lifetimes), Auger and radiative off. The 500 um base is not built (DT-Q29 A: its series resistance
// is outside 16.18 (g), 16.34's table). Uniform nodes N = 128, 256, 512, 1024 (x = 0 a face); forward biases 0.35,
// 0.40 and 0.45 V; w tau_p = 1e-3 ... 1e3 by decades, tau_p the p base's minority lifetime (SRH's linearization).
// Y_sim is Y_11, the anode's admittance (conduction plus displacement), S/cm^2.
// Reference (verify/reference/admittance, DT-Q30 A): Y_ref = sum over the bases of Shockley's finite-base admittance
// and the edges' motion, plus i w C_j; Shockley's fixed-edge form is printed beside it.
// Judged at each (V, w) (DT-Q31 A): |Y|, arg Y and Re Y. Criterion (16.19): |q_sim - q_ref| <= U_num + U_ref + eps_Y
// for |Y| and Re Y; for the phase, a complex error of size B moves the argument by at most asin(B / |Y_ref|), so
// |arg Y_sim - arg Y_ref| <= asin(min(1, (U_num + U_ref + eps_Y) / |Y_ref|)) with U_num of |Y| (the bound of 16.19
// read on the complex value). eps_Y = |Y_ref| [max(Delta n / N) at the edges + e^v / |e^v - 1| Delta phi / V_T
// (DT-Q25 A) + J_SNS / J_diff + 3 L_D / |L_eff|] (16.18 (g)), the larger base's L_D / |L_eff|, measured on the
// N = 1024 state. Applicable (else NOT APPLICABLE): |L_eff| >= 30 L_D on both bases and eps_Y <= 0.1 |Y_ref|.
// U_num per quantity = U_disc (DT-Q26 A, DT-Q28 B: the in-range three-grid GCI on 1024-512-256 with the asymptotic
// check on 512-256-128, and the two-grid estimate on N = 2048, 1024 with F_s = 3, p_th = 2 where the order observed on
// 2048-1024-512 is within 2 +- 0.5 (DT-Q32 A); the larger; not estimated where neither) + the admittance resolution
// + |eps_int| (the (h, tol) factorial at N = 1024: N = 2048, tol = 1e-12 and 1e-14). U_ref = the Riccati second-method difference + w |C_j -
// C_j,complex step| (+ the double's rounding).
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <expected>
#include <numbers>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/recombination.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/small_signal.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/admittance.hpp"
#include "NiTCAD/verify/reference/metrics.hpp"
#include "NiTCAD/verify/reference/qe_junction.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

constexpr double NA = 1e16, ND = 2e16, T = 300.0, L = 5e-4, tol_newton = 1e-12;
constexpr double applicability = 0.1;      // eps_Y <= 0.1 |Y_ref| (16.19)
constexpr double length_ratio = 30.0;      // |L_eff| >= 30 L_D (16.18 (g))
constexpr int ladder[4] = {1024, 512, 256, 128};  // finest first
constexpr int factorial_N = 2048;
constexpr double biases[3] = {0.35, 0.40, 0.45};
constexpr double omega_tau[7] = {1e-3, 1e-2, 1e-1, 1.0, 10.0, 1e2, 1e3};

ReferenceInfo info() {
    ReferenceInfo r;
    r.id = "TVB-18";
    r.cls = ReferenceClass::A;
    r.independence = Independence::I1;
    r.provenance = Provenance::PA;
    r.citation = "Shockley, BSTJ 28 (1949) eq 4.23 (the finite base and the edges' motion derived, PD)";
    r.version = "DECISIONS Unit 55 16.2 (l), 16.3, 16.18 (g), 16.19, 16.34, 16.35, 16.37; DT-Q28 B, DT-Q29-Q31 A, DT-Q32 A";
    r.equations = {1, 2, 3, 7};
    r.form = "Y = sum q D n0 e^v / V_T kappa coth(kappa W') + edge motion + i w C_j, kappa = sqrt((1/tau + i w) / D)";
    r.assumptions = {"low injection (measured)", "quasi-equilibrium across the depletion region (measured)",
                     "no depletion recombination admittance (J_SNS / J_diff, measured)", "Boltzmann statistics",
                     "abrupt doping", "the charge-balance edge (3 L_D / |L_eff|)",
                     "no series resistance (|Y| R printed, 16.34)"};
    r.derivation = "the linearized minority diffusion D du'' = (1/tau + i w) du at the moving QE edge";
    r.verification = "admittance_test: w = 0 against d J / dV, the long- and short-base limits, Riccati agreement";
    r.second_method = "the Riccati form by Dormand-Prince in double-double, at run time; C_j by complex step";
    r.approximation_error = "measured per bias: max(Delta n / N) + e^v / |e^v - 1| Delta phi / V_T + J_SNS / J_diff "
                            "+ 3 L_D / |L_eff| (16.18 (g))";
    return r;
}

device::Device diode(int N) {
    std::vector<double> x(static_cast<std::size_t>(N));
    for (int k = 0; k < N; ++k) x[static_cast<std::size_t>(k)] = L * (2.0 * k - (N - 1)) / (N - 1);
    x.front() = -L;
    x.back() = L;
    mesh::Mesh m = *mesh::make_tensor_grid(x);
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

solve::BiasOptions options(double tol) {
    solve::BiasOptions o;
    o.models.srh = true;
    o.models.auger = o.models.radiative = false;
    o.models.bgn = false;
    o.models.fermi_dirac = false;
    o.models.field_mobility = false;
    o.models.hole_mobility = assemble::HoleMobility::caughey_thomas;
    o.models.electron_mobility = assemble::ElectronMobility::caughey_thomas;
    o.newton.tol_update = tol;
    o.step_control.max_cuts = 30;
    return o;
}

// Linear interpolation of f (on ascending x) at s; log-linear when log_scale.
double at(const std::vector<double>& x, const std::vector<double>& f, double s, bool log_scale) {
    const auto it = std::upper_bound(x.begin(), x.end(), s);
    const std::size_t k = std::clamp<std::size_t>(static_cast<std::size_t>(it - x.begin()), 1, x.size() - 1) - 1;
    const double t = (s - x[k]) / (x[k + 1] - x[k]);
    if (log_scale) return std::exp((1.0 - t) * std::log(f[k]) + t * std::log(f[k + 1]));
    return (1.0 - t) * f[k] + t * f[k + 1];
}

// One bias at one level: Y_11 and its resolution per frequency, and the measured eps terms.
struct Point {
    std::vector<std::complex<double>> Y;
    std::vector<double> res;
    double low = 0.0, phi = 0.0;
    bool solved = true;
    double stall = 0.0;
};

using References = std::vector<reference::DiodeAdmittance>;  // per bias

// Every bias at N nodes and tolerance tol, swept from equilibrium. With may_stall (the factorial's tol / 100), a
// bias whose operating point Newton cannot bring to tol is kept as not solved with its last update, and the sweep
// restarts at the next bias (DT-Q23 A); otherwise a stop is a run error.
std::expected<std::vector<Point>, std::string> solve_level(int N, double tol, const References& refs,
                                                           const std::vector<double>& freqs, double VT,
                                                           bool may_stall) {
    const device::Device d = diode(N);
    const std::size_t n = d.mesh().node_count();
    std::vector<double> x(n);
    for (std::size_t i = 0; i < n; ++i) x[i] = d.mesh().points()[i][0];
    if (!std::is_sorted(x.begin(), x.end())) return std::unexpected("nodes not in x order");
    std::vector<Point> out;
    std::size_t k0 = 0;
    while (k0 < std::size(biases)) {
        std::vector<std::vector<double>> points;
        for (std::size_t k = k0; k < std::size(biases); ++k) points.push_back({biases[k], 0.0});
        solve::SmallSignalOptions o{.steady = options(tol), .frequencies_Hz = freqs};
        const auto s = solve::solve_small_signal(d, points, o);
        if (!s) return std::unexpected("N = " + std::to_string(N) + ": " + s.error().message);
        for (std::size_t q = 0; q < s->points.size(); ++q) {
            const auto& sp = s->points[q];
            const auto& ref = refs[out.size()].dc;
            Point p;
            for (std::size_t k = 0; k < freqs.size(); ++k) {
                p.Y.push_back(s->admittance(q, k, 0, 0));
                p.res.push_back(sp.resolution[k][0]);
            }
            const auto& bp = sp.dc;
            const double xp = -ref.edge_p_cm.to_double(), xn = ref.edge_n_cm.to_double();
            p.low = std::max(std::abs(at(x, bp.fields.n_cm3, xp, true) - ref.n_p0_cm3.to_double()) / NA,
                             std::abs(at(x, bp.fields.p_cm3, xn, true) - ref.p_n0_cm3.to_double()) / ND);
            const double dn = std::abs(at(x, bp.bands.electron_fermi_eV, xn, false) -
                                       at(x, bp.bands.electron_fermi_eV, xp, false));
            const double dp =
                std::abs(at(x, bp.bands.hole_fermi_eV, xn, false) - at(x, bp.bands.hole_fermi_eV, xp, false));
            const double v = biases[out.size()] / VT;
            p.phi = std::exp(v) / std::abs(std::expm1(v)) * std::max(dn, dp) / VT;
            out.push_back(std::move(p));
            ++k0;
        }
        if (!s->stopped) break;
        if (!may_stall) return std::unexpected("N = " + std::to_string(N) + ": " + s->stopped->message);
        Point p;
        p.solved = false;
        p.Y.assign(freqs.size(), {std::nan(""), std::nan("")});
        p.res.assign(freqs.size(), std::nan(""));
        p.stall = INFINITY;
        if (s->unfinished && !s->unfinished->iterations.empty()) p.stall = s->unfinished->iterations.back().update;
        out.push_back(std::move(p));
        ++k0;
    }
    return out;
}

// U_disc by the shared rule (registry.hpp's u_disc: DT-Q26 A, DT-Q28 B, DT-Q32 A), with the orders and the reason.
std::optional<double> disc_of(const double h[4], const double f[4], double h_half, double f_half, double& order,
                              double& order2, std::string& how) {
    const Discretization d = u_disc(std::span<const double, 4>(h, 4), std::span<const double, 4>(f, 4), h_half, f_half);
    order = d.order;
    order2 = d.order_coarse;
    how = d.how;
    return d.u;
}

}  // namespace

BenchmarkResult run_tvb18(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-18";
    r.title = "Diode small-signal admittance at forward bias";
    r.info = info();
    if (tier == Tier::fast) {
        r.notes.push_back("full tier only (16.7)");
        return r;
    }
    const physics::Semiconductor si = physics::silicon();
    const double VT = base::thermal_voltage(T);
    const double ni = physics::intrinsic_density(si, T);
    const double eps = si.parameters().eps_r * base::eps0_F_per_cm;
    reference::DiodeInput in;
    in.junction = {.q_C = base::q_C, .eps_F_per_cm = eps, .thermal_voltage_V = VT, .ni_cm3 = ni,
                   .acceptors_cm3 = NA, .donors_cm3 = ND};
    const double Dn = physics::caughey_thomas_mobility(si, physics::Carrier::electron, NA, T) * VT;
    const double Dp = physics::caughey_thomas_mobility(si, physics::Carrier::hole, ND, T) * VT;
    const auto tau = [&](physics::Carrier c, double N) { return physics::scharfetter_lifetime(si, c, N); };
    in.p_side = {.length_cm = L, .minority_diffusivity_cm2_s = Dn, .tau_n_s = tau(physics::Carrier::electron, NA),
                 .tau_p_s = tau(physics::Carrier::hole, NA)};
    in.n_side = {.length_cm = L, .minority_diffusivity_cm2_s = Dp, .tau_n_s = tau(physics::Carrier::electron, ND),
                 .tau_p_s = tau(physics::Carrier::hole, ND)};
    in.srh = true;
    r.notes.push_back("spec (NiTCAD's material data, I1): n_i " + shortest(ni) + " cm^-3, D_n(N_A) " + shortest(Dn) +
                      ", D_p(N_D) " + shortest(Dp) + " cm^2/s");

    // The frequencies, from the p base's minority lifetime.
    const std::vector<double> w0{0.0};
    const auto base0 = reference::diode_admittance(in, biases[0], w0);
    if (!base0) {
        r.run_error = "reference: " + base0.error();
        return r;
    }
    const double tau_p = base0->lifetime_p_s;
    std::vector<double> omegas, freqs;
    for (const double wt : omega_tau) {
        omegas.push_back(wt / tau_p);
        freqs.push_back(wt / tau_p / (2.0 * std::numbers::pi));
    }
    const double LD = std::max(std::sqrt(eps * VT / (base::q_C * NA)), std::sqrt(eps * VT / (base::q_C * ND)));
    References refs;
    std::vector<double> cs_diff;
    for (const double V : biases) {
        auto y = reference::diode_admittance(in, V, omegas);
        const auto cs = reference::qe_capacitance_complex_step(in.junction, V);
        if (!y || !cs) {
            r.run_error = "V = " + shortest(V) + ": reference: " + (!y ? y.error() : cs.error());
            return r;
        }
        cs_diff.push_back(std::abs(*cs - y->capacitance_F_per_cm2.to_double()));
        // The depletion region's minority storage (DT-Q34 A, 16.38) and its second method.
        const auto scr2 = reference::scr_storage_second_method(in.junction, V);
        if (!scr2) {
            r.run_error = "V = " + shortest(V) + ": reference C_SCR second method: " + scr2.error();
            return r;
        }
        const double scr = y->scr_storage_F_per_cm2.to_double();
        r.notes.push_back("V = " + shortest(V) + " V: C_SCR " + shortest(scr) + " F/cm^2 (" +
                          shortest(scr / y->capacitance_F_per_cm2.to_double()) + " of C_j); second method differs by " +
                          shortest(std::abs(*scr2 - scr)) + " F/cm^2");
        refs.push_back(std::move(*y));
    }
    r.notes.push_back("tau_p " + shortest(tau_p) + " s, tau_n " + shortest(base0->lifetime_n_s) +
                      " s (the bases' minority lifetimes); L_D (larger) " + shortest(LD) + " cm");

    std::vector<std::vector<Point>> lv;
    double hs[4];
    const double h_factorial = 2.0 / (factorial_N - 1);
    for (int k = 0; k < 4; ++k) {
        const auto s = solve_level(ladder[k], tol_newton, refs, freqs, VT, false);
        if (!s) {
            r.run_error = s.error();
            return r;
        }
        lv.push_back(*s);
        hs[k] = 2.0 / (ladder[k] - 1);
    }
    std::vector<std::vector<Point>> fac;  // (h, tol/100), (h/2, tol), (h/2, tol/100)
    for (const auto& [N, tol] : {std::pair{ladder[0], 0.01 * tol_newton}, std::pair{factorial_N, tol_newton},
                                 std::pair{factorial_N, 0.01 * tol_newton}}) {
        const auto s = solve_level(N, tol, refs, freqs, VT, tol < tol_newton);
        if (!s) {
            r.run_error = s.error() + " (factorial)";
            return r;
        }
        fac.push_back(*s);
    }

    for (std::size_t b = 0; b < std::size(biases); ++b) {
        const double V = biases[b];
        const auto& ref = refs[b];
        const Point& fine = lv[0][b];
        const double sns = ref.dc.sns_A_per_cm2.to_double() / ref.dc.diffusion_A_per_cm2.to_double();
        std::string stalled;
        for (int k = 0; k < 3; ++k) {
            if (!fac[static_cast<std::size_t>(k)][b].solved) {
                stalled += std::string(k == 0 ? " (h, tol/100)" : k == 1 ? " (h/2, tol)" : " (h/2, tol/100)") +
                           " not reached, last update " + shortest(fac[static_cast<std::size_t>(k)][b].stall);
            }
        }
        for (std::size_t k = 0; k < omegas.size(); ++k) {
            const auto& rp = ref.points[k];
            const std::complex<double> Yref = rp.total;
            const double mag_ref = std::abs(Yref);
            const double leff = std::min(std::abs(rp.p.effective_length_cm), std::abs(rp.n.effective_length_cm));
            const double edge = 3.0 * LD / leff;
            const double eps_rel = fine.low + fine.phi + std::abs(sns) + edge;
            const double uref = rp.second_method_difference + omegas[k] * cs_diff[b] + 1e-15 * mag_ref;
            const std::string at_vw = "V = " + shortest(V) + " V, w tau_p = " + shortest(omega_tau[k]);
            r.notes.push_back(at_vw + ": Y_ref " + shortest(Yref.real()) + " + i " + shortest(Yref.imag()) +
                              " S/cm^2 (Shockley fixed edges " + shortest(rp.shockley.real()) + " + i " +
                              shortest(rp.shockley.imag()) + ", i w C_j share " +
                              shortest(std::abs(rp.capacitive) / mag_ref) + "), Y_sim " + shortest(fine.Y[k].real()) +
                              " + i " + shortest(fine.Y[k].imag()) + "; eps_Y/|Y|: low injection " +
                              shortest(fine.low) + ", quasi-Fermi " + shortest(fine.phi) + ", SNS " + shortest(sns) +
                              ", edge " + shortest(edge) + "; |L_eff| / L_D " + shortest(leff / LD) + "; U_ref " +
                              shortest(uref) + stalled);
            const std::vector<Predicate> pr{
                {"|L_eff| / L_D >= 30", leff / LD, length_ratio, leff >= length_ratio * LD},
                {"eps_Y / |Y_ref| <= 0.1", eps_rel, applicability, eps_rel <= applicability}};

            // One judged quantity q(Y): magnitude, phase or real part.
            enum class Q { magnitude, phase, real };
            double u_mag = 0.0;
            bool est_mag = false;
            for (const Q which : {Q::magnitude, Q::real, Q::phase}) {
                const auto q_of = [&](std::complex<double> y) {
                    return which == Q::magnitude ? std::abs(y) : which == Q::real ? y.real() : std::arg(y);
                };
                const double q_ref = q_of(Yref), q_sim = q_of(fine.Y[k]);
                const double res_q = which == Q::phase ? fine.res[k] / std::abs(fine.Y[k]) : fine.res[k];
                const double f[4] = {q_sim, q_of(lv[1][b].Y[k]), q_of(lv[2][b].Y[k]), q_of(lv[3][b].Y[k])};
                const double q_h2 = q_of(fac[1][b].Y[k]);
                double p = 0.0, p2 = 0.0;
                std::string how;
                const auto ud = disc_of(hs, f, h_factorial, q_h2, p, p2, how);
                const std::vector<double> parts{ud ? *ud : 0.0, res_q};
                const Factorial fc = factorial(q_sim - q_ref, q_of(fac[0][b].Y[k]) - q_ref, q_h2 - q_ref,
                                               q_of(fac[2][b].Y[k]) - q_ref, parts);
                const bool estimated = ud.has_value() && fc.ok;
                const double unum = fc.ok ? fc.u : 0.0;
                if (which == Q::magnitude) {
                    u_mag = unum;
                    est_mag = estimated;
                }
                const double metric = std::abs(q_sim - q_ref);
                double bound = 0.0;
                std::string name, unit, criterion;
                if (which == Q::phase) {
                    bound = std::asin(std::min(1.0, (u_mag + uref + eps_rel * mag_ref) / mag_ref));
                    name = "arg Y";
                    unit = "rad";
                    criterion = "|arg Y_sim - arg Y_ref| <= asin((U_num,|Y| + U_ref + eps_Y) / |Y_ref|) = " +
                                shortest(bound);
                } else {
                    bound = unum + uref + eps_rel * mag_ref;
                    name = which == Q::magnitude ? "|Y|" : "Re Y";
                    unit = "S/cm^2";
                    criterion = "|q_sim - q_ref| <= U_num + U_ref + eps_Y = " + shortest(bound);
                }
                r.notes.push_back(at_vw + ", " + name + ": q_ref " + shortest(q_ref) + ", q_sim " + shortest(q_sim) +
                                  ", U_num " + shortest(unum) + " (U_disc: " + how + "; GCI orders " + shortest(p) +
                                  ", " + shortest(p2) + ", resolution " + shortest(res_q) + "); " + fc.note);
                ComparisonRecord cr{name + " " + at_vw, unit, V, q_sim, q_ref, which == Q::phase ? uref / mag_ref : uref,
                                    metric, criterion, {}};
                const bool est = which == Q::phase ? est_mag && estimated : estimated;
                cr.comparison = !est ? Comparison{Outcome::not_evaluated,
                                                  "U_num not estimated (" +
                                                      (ud ? std::string("the (h, tol) factorial") : how) +
                                                      (which == Q::phase && !est_mag ? std::string(", or |Y|'s") : "") +
                                                      (stalled.empty() ? std::string() : ":" + stalled) + ")",
                                                  {}}
                                     : judge(r.info, pr, metric <= bound, cr.criterion);
                r.comparisons.push_back(std::move(cr));
            }
        }
    }
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
