// TVB-19, linear transients of a uniform n-bar (DECISIONS Unit 55 §6, 16.2 (b), (m), 16.4; criteria 16.8, frozen;
// injection by the TVQ3 hook's time profile, DT-Q3 A, 16.25; construction 16.41).
//
// Code under test: solve_transient (BDF2 and backward Euler, fixed steps) on a uniform 1D silicon bar with ohmic
// ends, 300 K, Boltzmann, no band-gap narrowing, no field mobility, Caughey-Thomas mobility for both carriers.
// Carriers are injected by the hook's generation G(x) w(t), w = 0 up to t_on (the steady state at t = 0 included),
// 1 on (t_on, t_off], 0 after; both times are restart breakpoints, so steps land on them.
// (i) Dielectric relaxation: N_D = 1e14 cm^-3, 40 um, no bias, no recombination; electrons only, a Gaussian
//     (sigma 2 um) at the centre on (0.1, 0.5] tau_d. The excess charge's integral Q(t) = -q int (n - n(0)) dx at
//     t_off + k tau_d, k = 1 ... 4; the decay rate ln(Q_1 / Q_4) / (3 tau_d) against 1 / tau_d, tau_d = eps / sigma,
//     sigma = q (mu_n n_0 + mu_p p_0) (verify/reference/linear_transient). Applicable (else NOT APPLICABLE) while
//     max (n - n(0)) / N_D <= 1e-6 (the conductivity's change, a relative model error of that size) and the excess
//     at the end nodes is <= 1e-12 of its peak (no flux at the ends).
// (ii) Haynes-Shockley: N_D = 1e15 cm^-3, 200 um, 2 V across it (E = 100 V/cm, the holes drift toward x_max), SRH
//     on (Scharfetter lifetimes); pairs, a Gaussian (sigma 3 um) at 50 um on (0.1, 2] ns. At t = 40 ns: the holes'
//     excess area int dp dx, centre and width (trapezoids on the nodes) against haynes_shockley_moments with
//     D = V_T mu_p, v = mu_p E, tau = SRH's linearization at the equilibrium (minority_diffusion). Applicable while
//     max dp / N_D <= 1e-8 (16.2 (m)) and the excess at both ends is <= 1e-12 of its peak.
// Judged (16.8): in dt at fixed h, p-hat in [1.8, 2.3] (BDF2) and [0.8, 1.3] (backward Euler), on four dt levels;
// in h at fixed dt (BDF2), p-hat in [1.8, 2.3] on four h levels; at the finest (h, dt, BDF2) the decay rate
// relative <= 1e-3, the centre <= 1e-3 of the drift length v (t - t_on), the width and area relative <= 1e-3.
// The factorials at the finest (16.2 (p)): (h, tol) and (h, dt), each must be evaluated and non-interacting (4.1 (c));
// the 2x2x2 interaction is taken as these two faces.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <expected>
#include <functional>
#include <span>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/assemble/models.hpp"
#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/recombination.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/transient.hpp"
#include "NiTCAD/solve/waveform.hpp"
#include "NiTCAD/verify/benchmark/device_data.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/benchmark/tvb32.hpp"
#include "NiTCAD/verify/reference/first_integral.hpp"
#include "NiTCAD/verify/reference/linear_transient.hpp"
#include "NiTCAD/verify/reference/metrics.hpp"
#include "NiTCAD/verify/reference/minority_diffusion.hpp"
#include "NiTCAD/verify/study/study.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

constexpr double T = 300.0, tol_newton = 1e-12, tolerance = 1e-3;
constexpr int h_ladder[4] = {2049, 1025, 513, 257};  // finest first (odd: a node at the centre)
constexpr int h_half_N = 4097;

ReferenceInfo info() {
    ReferenceInfo r;
    r.id = "TVB-19";
    r.cls = ReferenceClass::E;
    r.independence = Independence::I1;
    r.provenance = Provenance::PD;
    r.citation = "Haynes and Shockley, Phys. Rev. 81, 835 (1951) (the model, PA); the moments and the relaxation "
                 "derived (PD, verify/reference/linear_transient.hpp)";
    r.version = "DECISIONS Unit 55 16.2 (m), 16.8, 16.18 (h), 16.41";
    r.equations = {1, 2, 3};
    r.form = "Q(t) ~ e^{-t / tau_d}, tau_d = eps / sigma; the drifting, spreading, decaying Gaussian's moments";
    r.assumptions = {"linear (measured: the excess over the doping)", "no flux at the ends (measured)",
                     "uniform field (ohmic ends, uniform doping)", "Boltzmann statistics"};
    r.derivation = "verify/reference/linear_transient.hpp";
    r.verification = "linear_transient_test: the uniform-age limit, the moments' equations, the short pulse";
    r.second_method = "Dormand-Prince on the relaxation; Gauss-Kronrod on the moments' age integrals";
    return r;
}

device::Device bar(int N, double length, double ND) {
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
                                    .contacts = {{"anode", device::ContactKind::ohmic, std::move(a)},
                                                 {"cathode", device::ContactKind::ohmic, std::move(c)}}});
}

// One case (relaxation or pulse): its device, source and judged times.
struct Case {
    double length, ND, V, sigma, centre, rate, t_on, t_off, t_end;
    bool srh, pairs;
    std::vector<double> times;  // output times
};

struct Run {
    std::vector<std::vector<double>> dn, dp;  // per output time, per node: excess over t = 0
    std::vector<double> x;
};

constexpr double relaxation_ND = 1e14, pulse_ND = 1e15;

// The relaxation case at tau_d = td, and the Haynes-Shockley pulse.
Case relaxation_case(double td) {
    Case relax{.length = 4e-3, .ND = relaxation_ND, .V = 0.0, .sigma = 2e-4, .centre = 2e-3, .rate = 0.0,
               .t_on = 0.1 * td, .t_off = 0.5 * td, .t_end = 4.5 * td, .srh = false, .pairs = false, .times = {}};
    relax.rate = 1e-6 * relaxation_ND / (relax.t_off - relax.t_on) * 0.5;  // peak excess about 1e-6 N_D at most
    for (int k = 1; k <= 4; ++k) relax.times.push_back(relax.t_off + k * td);
    return relax;
}

Case pulse_case() {
    return {.length = 2e-2, .ND = pulse_ND, .V = 2.0, .sigma = 3e-4, .centre = 5e-3, .rate = 5e14, .t_on = 1e-10,
            .t_off = 2e-9, .t_end = 4e-8, .srh = true, .pairs = true, .times = {4e-8}};
}

solve::BiasOptions case_options(const Case& c, double tol) {
    solve::BiasOptions o;
    o.models.srh = c.srh;
    o.models.auger = o.models.radiative = false;
    o.models.bgn = false;
    o.models.fermi_dirac = false;
    o.models.field_mobility = false;
    o.models.hole_mobility = assemble::HoleMobility::caughey_thomas;
    o.models.electron_mobility = assemble::ElectronMobility::caughey_thomas;
    o.newton.tol_update = tol;
    return o;
}

// One case on device d at fixed step dt.
std::expected<results::Transient, std::string> solve_case(const Case& c, const device::Device& d, double dt,
                                                          solve::Integrator integ, double tol) {
    const std::size_t n = d.mesh().node_count();
    assemble::ManufacturedSources ms;
    for (std::size_t i = 0; i < n; ++i) {
        const double z = (d.mesh().points()[i][0] - c.centre) / c.sigma;
        ms.electron_generation.push_back(c.rate * std::exp(-0.5 * z * z));
    }
    ms.hole_generation = c.pairs ? ms.electron_generation : std::vector<double>(n, 0.0);
    ms.profile_times_s = {c.t_on, c.t_off};
    ms.profile_weights = {0.0, 1.0, 0.0};
    solve::TransientOptions o;
    o.steady = case_options(c, tol);
    o.steady.models.manufactured = &ms;
    o.integrator = integ;
    o.t_end_s = c.t_end;
    o.adaptive = false;
    o.dt_initial_s = dt;
    o.dt_min_s = dt * 1e-6;
    o.dt_max_s = dt;
    o.output_times_s = c.times;
    const std::vector<solve::Waveform> w{solve::Waveform::constant(c.V), solve::Waveform::constant(0.0)};
    auto s = solve::solve_transient(d, w, o);
    if (!s) return std::unexpected("N = " + std::to_string(n) + ": " + s.error().message);
    if (s->stopped) return std::unexpected("N = " + std::to_string(n) + ": " + s->stopped->message);
    return std::move(*s);
}

std::expected<Run, std::string> run_case(const Case& c, int N, double dt, solve::Integrator integ, double tol) {
    const device::Device d = bar(N, c.length, c.ND);
    const std::size_t n = d.mesh().node_count();
    Run out;
    for (std::size_t i = 0; i < n; ++i) out.x.push_back(d.mesh().points()[i][0]);
    const auto s = solve_case(c, d, dt, integ, tol);
    if (!s) return std::unexpected(s.error());
    const auto at = [&](double t) {
        return std::find_if(s->snapshots.begin(), s->snapshots.end(),
                            [&](const auto& sn) { return std::abs(sn.time_s - t) <= 1e-9 * dt; });
    };
    const auto s0 = at(0.0);
    if (s0 == s->snapshots.end()) return std::unexpected("no snapshot at t = 0");
    for (const double t : c.times) {
        const auto sn = at(t);
        if (sn == s->snapshots.end()) return std::unexpected("no snapshot at t = " + shortest(t));
        std::vector<double> a(n), b(n);
        for (std::size_t i = 0; i < n; ++i) {
            a[i] = sn->fields.n_cm3[i] - s0->fields.n_cm3[i];
            b[i] = sn->fields.p_cm3[i] - s0->fields.p_cm3[i];
        }
        out.dn.push_back(std::move(a));
        out.dp.push_back(std::move(b));
    }
    return out;
}

// Trapezoid moments of f on x: area, centre, width.
struct Moments {
    double area = 0.0, centre = 0.0, width = 0.0, peak = 0.0, ends = 0.0;
};
Moments moments(const std::vector<double>& x, const std::vector<double>& f) {
    double m0 = 0.0, m1 = 0.0, m2 = 0.0;
    Moments m;
    for (std::size_t i = 0; i + 1 < x.size(); ++i) {
        const double h = x[i + 1] - x[i];
        m0 += 0.5 * h * (f[i] + f[i + 1]);
        m1 += 0.5 * h * (x[i] * f[i] + x[i + 1] * f[i + 1]);
    }
    m.area = m0;
    m.centre = m1 / m0;
    for (std::size_t i = 0; i + 1 < x.size(); ++i) {
        const double h = x[i + 1] - x[i];
        const double a = x[i] - m.centre, b = x[i + 1] - m.centre;
        m2 += 0.5 * h * (a * a * f[i] + b * b * f[i + 1]);
    }
    m.width = std::sqrt(m2 / m0);
    for (const double v : f) m.peak = std::max(m.peak, std::abs(v));
    m.ends = std::max(std::abs(f.front()), std::abs(f.back()));
    return m;
}

// One judged quantity on the three ladders and the factorial states: q(run) and its signed error.
struct Ladders {
    std::vector<Run> dt_bdf2, dt_be, h;  // finest first; h[0] is dt_bdf2[0]
    Run tol, h2, h2_tol, dt2, h2_dt2;    // the factorials' extra states at the finest
};

}  // namespace

BenchmarkResult run_tvb19(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-19";
    r.title = "Linear transients: dielectric relaxation and the Haynes-Shockley pulse";
    r.info = info();
    if (tier == Tier::fast) {
        r.notes.push_back("full tier only (16.7)");
        return r;
    }
    const physics::Semiconductor si = physics::silicon();
    const double VT = base::thermal_voltage(T), q = base::q_C;
    const double ni = physics::intrinsic_density(si, T);
    const double eps = si.parameters().eps_r * base::eps0_F_per_cm;

    // ---- The two cases and their references ----
    const double ND1 = relaxation_ND;
    const reference::NeutralBulk b1 = reference::neutral_bulk(ni, ND1);
    const double mun1 = physics::caughey_thomas_mobility(si, physics::Carrier::electron, ND1, T);
    const double mup1 = physics::caughey_thomas_mobility(si, physics::Carrier::hole, ND1, T);
    const double sigma_over_q = mun1 * b1.n_b.to_double() + mup1 * b1.p_b.to_double();
    const auto taud = reference::dielectric_relaxation_time(eps, q, 1.0, sigma_over_q);
    if (!taud) {
        r.run_error = "reference: " + taud.error();
        return r;
    }
    const double td = taud->to_double();
    const auto relax2 = reference::relaxation_second_method(*taud, 5.0 * td);
    if (!relax2) {
        r.run_error = "reference: " + relax2.error();
        return r;
    }
    const Case relax = relaxation_case(td);

    const double ND2 = pulse_ND, L2 = 2e-2, V2 = 2.0, E2 = V2 / L2;
    const reference::NeutralBulk b2 = reference::neutral_bulk(ni, ND2);
    const double mup2 = physics::caughey_thomas_mobility(si, physics::Carrier::hole, ND2, T);
    const double inv_tau = reference::srh_linear_inverse_lifetime(
                               b2.n_b.to_double(), b2.p_b.to_double(), ni,
                               physics::scharfetter_lifetime(si, physics::Carrier::electron, ND2),
                               physics::scharfetter_lifetime(si, physics::Carrier::hole, ND2))
                               .to_double();
    const Case hs = pulse_case();
    const reference::PulseInput pin{.diffusivity_cm2_s = VT * mup2, .velocity_cm_s = mup2 * E2,
                                    .lifetime_s = 1.0 / inv_tau, .centre_cm = hs.centre, .sigma_cm = hs.sigma,
                                    .rate_cm3_s = hs.rate, .t_on_s = hs.t_on, .t_off_s = hs.t_off};
    const auto mref = reference::haynes_shockley_moments(pin, hs.times[0]);
    if (!mref) {
        r.run_error = "reference: " + mref.error();
        return r;
    }
    const double drift = pin.velocity_cm_s * (hs.times[0] - hs.t_on);
    r.notes.push_back("relaxation: tau_d " + shortest(td) + " s (mu_n " + shortest(mun1) + ", mu_p " + shortest(mup1) +
                      " cm^2/Vs); Dormand-Prince differs by " + shortest(*relax2));
    r.notes.push_back("Haynes-Shockley: D_p " + shortest(pin.diffusivity_cm2_s) + " cm^2/s, v " +
                      shortest(pin.velocity_cm_s) + " cm/s, tau " + shortest(pin.lifetime_s) + " s; at 40 ns area " +
                      shortest(mref->area_cm2.to_double()) + " cm^-2, centre " + shortest(mref->centre_cm.to_double()) +
                      " cm, width " + shortest(mref->width_cm.to_double()) + " cm; second method " +
                      shortest(mref->area_difference) + ", " + shortest(mref->centre_difference) + ", " +
                      shortest(mref->width_difference));

    // ---- The runs ----
    const auto ladders = [&](const Case& c, double dt0, Ladders& L) -> std::expected<void, std::string> {
        for (int k = 0; k < 4; ++k) {
            const double dt = dt0 * std::ldexp(1.0, k);
            auto a = run_case(c, h_ladder[0], dt, solve::Integrator::bdf2, tol_newton);
            auto b = run_case(c, h_ladder[0], dt, solve::Integrator::backward_euler, tol_newton);
            if (!a || !b) return std::unexpected(!a ? a.error() : b.error());
            L.dt_bdf2.push_back(std::move(*a));
            L.dt_be.push_back(std::move(*b));
        }
        L.h.push_back(L.dt_bdf2[0]);
        for (int k = 1; k < 4; ++k) {
            auto a = run_case(c, h_ladder[k], dt0, solve::Integrator::bdf2, tol_newton);
            if (!a) return std::unexpected(a.error());
            L.h.push_back(std::move(*a));
        }
        const auto bd = solve::Integrator::bdf2;
        auto t1 = run_case(c, h_ladder[0], dt0, bd, 0.01 * tol_newton);
        auto t2 = run_case(c, h_half_N, dt0, bd, tol_newton);
        auto t3 = run_case(c, h_half_N, dt0, bd, 0.01 * tol_newton);
        auto t4 = run_case(c, h_ladder[0], 0.5 * dt0, bd, tol_newton);
        auto t5 = run_case(c, h_half_N, 0.5 * dt0, bd, tol_newton);
        for (auto* t : {&t1, &t2, &t3, &t4, &t5}) {
            if (!*t) return std::unexpected(t->error() + " (factorial)");
        }
        L.tol = std::move(*t1);
        L.h2 = std::move(*t2);
        L.h2_tol = std::move(*t3);
        L.dt2 = std::move(*t4);
        L.h2_dt2 = std::move(*t5);
        return {};
    };
    Ladders LR, LH;
    const double dtR = td / 160.0, dtH = 5e-11;
    if (auto e = ladders(relax, dtR, LR); !e) {
        r.run_error = "relaxation, " + e.error();
        return r;
    }
    if (auto e = ladders(hs, dtH, LH); !e) {
        r.run_error = "Haynes-Shockley, " + e.error();
        return r;
    }

    // ---- Judging ----
    // A quantity: its signed relative error on a run, and the applicability predicates measured on the finest run.
    const auto judge_quantity = [&](const std::string& name, const Ladders& L, double dt0,
                                    const std::function<double(const Run&)>& err,
                                    const std::vector<Predicate>& pr) {
        const std::string band2 = "p-hat in [1.8, 2.3] (p_th 2)", band1 = "p-hat in [0.8, 1.3] (p_th 1)";
        const auto study_of = [&](const std::vector<Run>& runs, const std::vector<double>& steps) {
            std::vector<double> es, us;
            const double u_iter = std::abs(err(L.dt_bdf2[0]) - err(L.tol));
            for (const Run& x : runs) {
                es.push_back(std::abs(err(x)));
                us.push_back(u_iter);
            }
            return study::order_study(steps, es, us);
        };
        std::vector<double> dts, hs;
        for (int k = 0; k < 4; ++k) {
            dts.push_back(dt0 * std::ldexp(1.0, k));
            hs.push_back(1.0 / (h_ladder[k] - 1));
        }
        const study::OrderStudy s_b2 = study_of(L.dt_bdf2, dts), s_be = study_of(L.dt_be, dts),
                                s_h = study_of(L.h, hs);
        const study::OrderVerdict v_b2 = study::judge_order(s_b2, {.p_th = 2.0});
        const study::OrderVerdict v_be = study::judge_order(s_be, {.p_th = 1.0});
        const study::OrderVerdict v_h = study::judge_order(s_h, {.p_th = 2.0});
        const double e0 = err(L.dt_bdf2[0]);
        const auto gci = [&](const study::OrderStudy& s) {
            if (s.orders.empty() || !(s.observed > 0.0)) return 0.0;
            return reference::gci_safety_factor_three_grids * std::abs(s.error[1] - s.error[0]) /
                   (std::pow(s.h[1] / s.h[0], s.observed) - 1.0);
        };
        const std::vector<double> parts{gci(s_b2), gci(s_h), std::abs(e0 - err(L.tol))};
        const Factorial f_tol = factorial(e0, err(L.tol), err(L.h2), err(L.h2_tol), parts);
        Factorial f_dt = factorial(e0, err(L.dt2), err(L.h2), err(L.h2_dt2), parts);
        f_dt.note.replace(0, std::string("factorial (h, tol)").size(), "factorial (h, dt)");
        const auto orders = [](const study::OrderStudy& s) {
            std::string o;
            for (const double p : s.orders) o += " " + shortest(p);
            return o + " (spread " + shortest(s.order_spread) + ")";
        };
        r.notes.push_back(name + ": finest error " + shortest(e0) + "; dt orders (BDF2)" + orders(s_b2) +
                          "; dt orders (BE)" + orders(s_be) + "; h orders" + orders(s_h) + "; " + f_tol.note + "; " +
                          f_dt.note);
        const auto add = [&](const std::string& what, double metric, double target, const std::string& crit,
                             const study::OrderVerdict* v, bool met) {
            ComparisonRecord c{name + " " + what, "1", 0.0, metric, target, 0.0, metric, crit, {}};
            if (v && !v->judged && std::all_of(pr.begin(), pr.end(), [](const Predicate& p) { return p.holds; })) {
                c.comparison = Comparison{Outcome::not_evaluated, v->reason, {}};
            } else {
                c.comparison = judge(r.info, pr, v ? v->pass : met, crit);
            }
            if (c.comparison.outcome != Outcome::not_applicable) {
                block_if_interacting(c, f_tol);
                block_if_interacting(c, f_dt);
            }
            r.comparisons.push_back(std::move(c));
        };
        add("order in dt (BDF2)", s_b2.observed, 2.0, band2, &v_b2, false);
        add("order in dt (backward Euler)", s_be.observed, 1.0, band1, &v_be, false);
        add("order in h (BDF2)", s_h.observed, 2.0, band2, &v_h, false);
        add("finest accuracy", std::abs(e0), 0.0, "relative <= 1e-3", nullptr, std::abs(e0) <= tolerance);
    };

    // (i) The decay rate.
    {
        const auto rate = [&](const Run& x) {
            const auto Q = [&](std::size_t k) { return moments(x.x, x.dn[k]).area; };
            return std::log(Q(0) / Q(3)) / (3.0 * td);
        };
        const auto err = [&](const Run& x) { return rate(x) * td - 1.0; };
        double peak = 0.0, ends = 0.0;
        for (const auto& f : LR.dt_bdf2[0].dn) {
            const Moments m = moments(LR.dt_bdf2[0].x, f);
            peak = std::max(peak, m.peak);
            ends = std::max(ends, m.ends / m.peak);
        }
        const std::vector<Predicate> pr{{"max dn / N_D <= 1e-6 (linear)", peak / ND1, 1e-6, peak <= 1e-6 * ND1},
                                        {"end excess / peak <= 1e-12", ends, 1e-12, ends <= 1e-12}};
        r.notes.push_back("relaxation: rate * tau_d " + shortest(rate(LR.dt_bdf2[0]) * td) + "; max dn / N_D " +
                          shortest(peak / ND1) + ", end excess / peak " + shortest(ends));
        judge_quantity("decay rate (relaxation)", LR, dtR, err, pr);
    }
    // (ii) The pulse's moments.
    {
        const Moments fine = moments(LH.dt_bdf2[0].x, LH.dt_bdf2[0].dp[0]);
        const std::vector<Predicate> pr{
            {"max dp / N_D <= 1e-8 (low injection)", fine.peak / ND2, 1e-8, fine.peak <= 1e-8 * ND2},
            {"end excess / peak <= 1e-12", fine.ends / fine.peak, 1e-12, fine.ends <= 1e-12 * fine.peak}};
        r.notes.push_back("Haynes-Shockley at 40 ns: area " + shortest(fine.area) + ", centre " + shortest(fine.centre) +
                          ", width " + shortest(fine.width) + "; max dp / N_D " + shortest(fine.peak / ND2) +
                          ", end excess / peak " + shortest(fine.ends / fine.peak));
        const double A = mref->area_cm2.to_double(), C = mref->centre_cm.to_double(), W = mref->width_cm.to_double();
        judge_quantity("pulse area", LH, dtH,
                       [&](const Run& x) { return moments(x.x, x.dp[0]).area / A - 1.0; }, pr);
        judge_quantity("pulse centre (per drift length)", LH, dtH,
                       [&](const Run& x) { return (moments(x.x, x.dp[0]).centre - C) / drift; }, pr);
        judge_quantity("pulse width", LH, dtH,
                       [&](const Run& x) { return moments(x.x, x.dp[0]).width / W - 1.0; }, pr);
    }
    apply_known_findings(r, known_findings());
    return r;
}

// TVB-32 on this benchmark's devices (16.44): C6 on both cases at the coarsest h (N = 257) and the coarsest step of
// the dt ladder (8 dt0), relaxation with BDF2 and backward Euler, the pulse with BDF2; each device's zero-bias
// point and equilibrium.
Tvb32Run tvb19_consistency() {
    Tvb32Run run;
    const physics::Semiconductor si = physics::silicon();
    const double ni = physics::intrinsic_density(si, T), eps = si.parameters().eps_r * base::eps0_F_per_cm;
    const reference::NeutralBulk b1 = reference::neutral_bulk(ni, relaxation_ND);
    const double sigma_over_q =
        physics::caughey_thomas_mobility(si, physics::Carrier::electron, relaxation_ND, T) * b1.n_b.to_double() +
        physics::caughey_thomas_mobility(si, physics::Carrier::hole, relaxation_ND, T) * b1.p_b.to_double();
    const auto taud = reference::dielectric_relaxation_time(eps, base::q_C, 1.0, sigma_over_q);
    if (!taud) return {{}, "TVB-19: reference: " + taud.error()};
    const double td = taud->to_double();
    const int N = h_ladder[3];
    std::string error;
    struct Item {
        Case c;
        double dt;
        const char* name;
    };
    for (const Item& it : {Item{relaxation_case(td), 8.0 * td / 160.0, "relaxation"},
                           Item{pulse_case(), 8.0 * 5e-11, "Haynes-Shockley"}}) {
        const device::Device d = bar(N, it.c.length, it.c.ND);
        const std::string name = std::string("TVB-19 ") + it.name + " bar N = " + std::to_string(N);
        const solve::BiasOptions o = case_options(it.c, tol_newton);
        const consistency::DeviceData data = device_data(d, o.models);
        Tvb32Collector col{run.entries, name, d.contacts().size()};
        const auto p0 = solve::solve_bias(d, std::vector<double>{0.0, 0.0}, o);
        if (!p0) return {std::move(run.entries), name + ": " + p0.error().message};
        col.dc(*p0, d.mesh(), data, bias_label(*p0));
        if (!col.equilibrium(d, o.models, error)) return {std::move(run.entries), error};
        const bool relax = it.c.V == 0.0;
        for (const solve::Integrator integ : {solve::Integrator::bdf2, solve::Integrator::backward_euler}) {
            if (!relax && integ == solve::Integrator::backward_euler) continue;
            const auto t = solve_case(it.c, d, it.dt, integ, tol_newton);
            if (!t) return {std::move(run.entries), name + ": " + t.error()};
            col.transient(*t, d.mesh(), data,
                          std::string(integ == solve::Integrator::bdf2 ? "BDF2" : "backward Euler") + ", dt " +
                              shortest(it.dt) + " s");
        }
    }
    return run;
}

}  // namespace NiTCAD::verify::benchmark
