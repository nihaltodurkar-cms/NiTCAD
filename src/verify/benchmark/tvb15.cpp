// TVB-15, neutral-region minority diffusion with SRH (DECISIONS Unit 55 §6, 16.2 (e), 16.4; criteria 16.8 (E) and
// 16.19 (A), frozen).
//
// Code under test: the drift-diffusion system and solve_bias on a long uniform p-type bar. Construction:
//   - silicon, N_A = 1e16 cm^-3 on [0, L] with L = 1 mm (about six diffusion lengths, 164 um, so that SRH shapes
//     the profile), ohmic contacts at both ends, no bias, 300 K;
//   - SRH with Scharfetter lifetimes only (Auger, radiative and band-gap narrowing off), Boltzmann statistics,
//     Caughey-Thomas mobilities;
//   - the TVQ3 hook's localized generation G(x) = G0 exp(-(x - L/2)^2 / (2 s^2)), s = 100 um, equal for electrons and
//     holes and lumped at the nodes (CV-Q5);
//   - uniform nodes, N = 64 ... 512.
// Reference (verify/reference/minority_diffusion): the linear minority equation D u'' - u / tau = -G with
// u(0) = u(L) = 0. Its Green's function is integrated by Gauss-Kronrod in double-double, and Chebyshev collocation is
// run at run time as the second method (I3). The inputs are NiTCAD's material data at N_A (I1, as TVB-07): D = mu_n V_T
// with Caughey-Thomas mu_n, and 1 / tau the exact linearization of SRH at equilibrium (n_ie, tau_n, tau_p).
//
// DD departs from the linear model by O(Delta n / N_A): the ambipolar field and the SRH nonlinearity. Hence:
// - E (16.8), at G0 giving a reference peak of 5e-9 N_A, with the predicate max Delta n / N_A <= 1e-8 measured on
//   every solved state (else not_evaluated):
//   (a) p-hat of relative Delta n in L-infinity in [1.8, 2.3], with the asymptotic check (4 levels);
//   (b) finest (N = 512) <= 1e-4. U_ref <= 1e-5 (a tenth of (b), §4.4 E (c)), else not_evaluated.
// - A (16.18 (c), 16.19; DT-Q14), the injection ladder at N = 512: G = g0 x {1, 2, 4, 8}, g0 giving a reference
//   peak of 1e-3 N_A. Applicability: max Delta n / N_A <= 1e-1 at the largest G (else not_applicable).
//   With v(G) = (Delta n_sim - u_ref) / G, divided by the peak of u_ref / G:
//   (a) the departure's order in G from the successive differences, p_k = log2(|v_{k+2} - v_{k+1}| / |v_{k+1} - v_k|)
//       in L-infinity, judged on the smallest-G triple, in [0.8, 1.3];
//   (b) the G -> 0 limit by Richardson with that order, v_0 = v(g0) - (v(2 g0) - v(g0)) / (2^p - 1), agrees with the
//       E solve's v at the same mesh within the E comparison's U on that mesh. That U is U_disc (the GCI of the
//       finest pair, F_s = 1.25, §4.3) + U_iter + U_ref + |eps_int|.
// The (h, tol) factorial (16.2 (p), DT-Q23 A) at the E solve's N = 512: tol = 1e-12, 1e-14 and N = 512, 1023; eps_int
// enters U_E, and above 10 % of U_E it blocks the E passes (DECISIONS 4.1 (c)).
// Iterative error (CV-Q8): the final Newton correction, which for a density is relative to itself, propagated to
// Delta n: U_iter = update (1 + n0 / max Delta n). On the injection ladder the largest U_iter must be below a tenth of
// the smallest difference it orders, else both A comparisons are not_evaluated.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <expected>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/recombination.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/metrics.hpp"
#include "NiTCAD/verify/reference/minority_diffusion.hpp"
#include "NiTCAD/verify/study/study.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

constexpr double L_cm = 1e-1, NA = 1e16, sigma_cm = 1e-2, T = 300.0;
constexpr double e_peak = 5e-9, a_peak = 1e-3;          // reference peaks / N_A
constexpr double e_predicate = 1e-8, a_predicate = 1e-1;  // max Delta n / N_A
constexpr double e_tolerance = 1e-4, tol_newton = 1e-12;

ReferenceInfo info() {
    ReferenceInfo r;
    r.id = "TVB-15";
    r.cls = ReferenceClass::E;
    r.independence = Independence::I1;
    r.provenance = Provenance::PD;
    r.citation = "Linear minority-carrier diffusion with the linearized SRH lifetime (Green's function)";
    r.version = "DECISIONS Unit 55 16.2 (e), 16.4, 16.8, 16.18 (c), 16.19";
    r.equations = {1};
    r.form = "D u'' - u / tau = -G, u(0) = u(L) = 0; 1/tau = (n0 + p0) / (tau_p (n0 + n_ie) + tau_n (p0 + n_ie))";
    r.assumptions = {"low injection (predicate measured)", "quasi-neutral excess", "D = mu_n V_T (Einstein)",
                     "NiTCAD's mu_n, n_ie, tau_n, tau_p at N_A as inputs (I1)"};
    r.derivation = "first variation of SRH at equilibrium; Green's function of -D u'' + u / tau with Dirichlet ends";
    r.verification = "erfc closed form on a long bar, collocation, the derivative jump (minority_diffusion_test)";
    r.second_method = "Chebyshev collocation, 8 pieces of degree 40, in double-double, at run time";
    return r;
}

struct Inputs {
    double D = 0.0, tau = 0.0, n0 = 0.0, p0 = 0.0;
};

Inputs inputs() {
    const physics::Semiconductor si = physics::silicon();
    Inputs in;
    const double ni = physics::intrinsic_density(si, T);
    in.p0 = 0.5 * NA + std::sqrt(0.25 * NA * NA + ni * ni);
    in.n0 = ni * ni / in.p0;
    in.D = physics::caughey_thomas_mobility(si, physics::Carrier::electron, NA, T) * base::thermal_voltage(T);
    const double tn = physics::scharfetter_lifetime(si, physics::Carrier::electron, NA);
    const double tp = physics::scharfetter_lifetime(si, physics::Carrier::hole, NA);
    in.tau = 1.0 / reference::srh_linear_inverse_lifetime(in.n0, in.p0, ni, tn, tp).hi;
    return in;
}

reference::MinorityProblem problem(const Inputs& in, double G0) {
    return {.diffusivity_cm2_s = in.D, .lifetime_s = in.tau, .length_cm = L_cm, .peak_generation = G0,
            .centre_cm = 0.5 * L_cm, .width_cm = sigma_cm};
}

device::Device bar(int N) {
    std::vector<double> x(static_cast<std::size_t>(N));
    for (int i = 0; i < N; ++i) x[static_cast<std::size_t>(i)] = L_cm * i / (N - 1);
    mesh::Mesh m = *mesh::make_tensor_grid(x);
    const std::size_t n = m.node_count();
    auto left = m.find_boundary("x_min")->nodes;
    auto right = m.find_boundary("x_max")->nodes;
    return *device::Device::create({.mesh = std::move(m),
                                    .temperature_K = T,
                                    .regions = {{"silicon", physics::silicon()}},
                                    .node_region = std::vector<device::RegionId>(n, 0),
                                    .donors = std::vector<double>(n, 0.0),
                                    .acceptors = std::vector<double>(n, NA),
                                    .contacts = {{"left", device::ContactKind::ohmic, std::move(left)},
                                                 {"right", device::ContactKind::ohmic, std::move(right)}}});
}

solve::BiasOptions options(double tol = tol_newton) {
    solve::BiasOptions o;
    o.models.srh = true;
    o.models.auger = o.models.radiative = false;
    o.models.bgn = false;
    o.models.fermi_dirac = false;
    o.models.hole_mobility = assemble::HoleMobility::caughey_thomas;
    o.newton.tol_update = tol;
    return o;
}

struct Solved {
    std::vector<double> dn;  // Delta n per node
    double u_iter = 0.0;     // relative to max Delta n
    double max_dn = 0.0;
};

// Delta n = n(G) - n(0) at every node, G scaled from the unit-peak generation by G0.
std::expected<Solved, std::string> solve_level(const device::Device& d, const std::vector<double>& n_eq,
                                               const reference::MinorityProblem& unit, double G0, double n0,
                                               double tol = tol_newton) {
    const std::size_t n = d.mesh().node_count();
    assemble::ManufacturedSources ms;
    for (std::size_t i = 0; i < n; ++i) {
        const double z = (d.mesh().points()[i][0] - unit.centre_cm) / unit.width_cm;
        ms.electron_generation.push_back(G0 * std::exp(-0.5 * z * z));
    }
    ms.hole_generation = ms.electron_generation;
    solve::BiasOptions o = options(tol);
    o.models.manufactured = &ms;  // alive through the solve
    const auto s = solve::solve_bias(d, std::vector<double>{0.0, 0.0}, o);
    if (!s) return std::unexpected(s.error().message);
    Solved out;
    for (std::size_t i = 0; i < n; ++i) {
        out.dn.push_back(s->fields.n_cm3[i] - n_eq[i]);
        out.max_dn = std::max(out.max_dn, std::abs(out.dn.back()));
    }
    const auto& it = s->convergence.iterations;
    out.u_iter = (it.empty() ? 0.0 : it.back().update) * (1.0 + n0 / out.max_dn);
    return out;
}

double linf(const std::vector<double>& a) {
    double m = 0.0;
    for (const double v : a) m = std::max(m, std::abs(v));
    return m;
}

// The E error (relative Delta n in L-infinity) at N nodes with Newton tolerance tol: the factorial's extra levels.
std::expected<double, std::string> e_error(int N, double tol, const reference::MinorityProblem& unit, double G0,
                                           double u1, double n0) {
    const device::Device d = bar(N);
    const auto eq = solve::solve_bias(d, std::vector<double>{0.0, 0.0}, options(tol));
    if (!eq) return std::unexpected("N = " + std::to_string(N) + ": equilibrium: " + eq.error().message);
    const auto s = solve_level(d, eq->fields.n_cm3, unit, G0, n0, tol);
    if (!s) return std::unexpected("N = " + std::to_string(N) + ": " + s.error());
    double worst = 0.0;
    for (std::size_t i = 0; i < s->dn.size(); ++i) {
        const double x = d.mesh().points()[i][0];
        const bool end = x <= 0.0 || x >= L_cm;
        const double u = end ? 0.0 : reference::minority_excess(unit, std::clamp(x, 0.0, L_cm)).value.hi;
        worst = std::max(worst, std::abs(s->dn[i] - G0 * u));
    }
    return worst / (G0 * u1);
}

}  // namespace

BenchmarkResult run_tvb15(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-15";
    r.title = "Neutral-region minority diffusion with SRH";
    r.info = info();
    if (tier == Tier::fast) {
        r.notes.push_back("full tier only (16.7)");
        return r;
    }
    const Inputs in = inputs();
    const reference::MinorityProblem unit = problem(in, 1.0);
    const reference::MinorityValue peak = reference::minority_excess(unit, 0.5 * L_cm);
    const double u1 = peak.value.hi;  // u(L/2) per unit G0
    const double gE = e_peak * NA / u1, gA = a_peak * NA / u1;
    r.notes.push_back("D = " + shortest(in.D) + " cm^2/s, tau = " + shortest(in.tau) + " s, diffusion length " +
                      shortest(std::sqrt(in.D * in.tau)) + " cm, n0 = " + shortest(in.n0) + " cm^-3; G0 (E) = " +
                      shortest(gE) + ", g0 (A) = " + shortest(gA) + " cm^-3 s^-1");

    // I3: collocation against the Green's function at eleven points, relative to the peak.
    std::vector<double> probe;
    for (int j = 0; j <= 10; ++j) probe.push_back(L_cm * (0.05 + 0.09 * j));
    const auto col = reference::minority_collocation(unit, probe, 8, 40);
    double i3 = 0.0, quad = peak.abs_error / u1;
    bool green_ok = peak.ok;
    for (std::size_t j = 0; j < probe.size(); ++j) {
        const auto g = reference::minority_excess(unit, probe[j]);
        green_ok = green_ok && g.ok;
        quad = std::max(quad, g.abs_error / u1);
        if (col.ok) i3 = std::max(i3, std::abs((col.values[j] - g.value).hi) / u1);
    }
    const double u_ref = std::max({i3, quad, col.tail_estimate / u1});
    r.notes.push_back("I3: collocation agrees with the Green's function to " + shortest(i3) +
                      " of the peak; quadrature " + shortest(quad) + ", tail " + shortest(col.tail_estimate / u1) +
                      "; U_ref " + shortest(u_ref));
    if (!col.ok || !green_ok) {
        r.run_error = "the reference did not converge (collocation " + std::string(col.ok ? "ok" : "failed") +
                      ", quadrature " + (green_ok ? "ok" : "failed") + ")";
        return r;
    }

    struct Level {
        int N = 0;
        double h = 0.0, error = 0.0, u_iter = 0.0, max_ratio = 0.0;
        std::vector<double> u_ref;  // per node, per unit G0
        std::vector<double> dn;
        std::vector<double> n_eq;
    };
    std::vector<Level> levels;  // finest first
    bool predicate = true;
    for (const int N : {512, 256, 128, 64}) {
        const device::Device d = bar(N);
        const auto eq = solve::solve_bias(d, std::vector<double>{0.0, 0.0}, options());
        if (!eq) {
            r.run_error = "N = " + std::to_string(N) + ": equilibrium: " + eq.error().message;
            return r;
        }
        Level lv;
        lv.N = N;
        lv.h = L_cm / (N - 1);
        lv.n_eq = eq->fields.n_cm3;
        for (std::size_t i = 0; i < d.mesh().node_count(); ++i) {
            const double x = d.mesh().points()[i][0];
            const bool end = x <= 0.0 || x >= L_cm;  // u = 0 at the ohmic ends
            lv.u_ref.push_back(end ? 0.0 : reference::minority_excess(unit, std::clamp(x, 0.0, L_cm)).value.hi);
        }
        const auto s = solve_level(d, lv.n_eq, unit, gE, in.n0);
        if (!s) {
            r.run_error = "N = " + std::to_string(N) + ": " + s.error();
            return r;
        }
        double worst = 0.0;
        for (std::size_t i = 0; i < s->dn.size(); ++i) worst = std::max(worst, std::abs(s->dn[i] - gE * lv.u_ref[i]));
        lv.error = worst / (gE * u1);
        lv.u_iter = s->u_iter;
        lv.max_ratio = s->max_dn / NA;
        lv.dn = s->dn;
        predicate = predicate && lv.max_ratio <= e_predicate;
        r.notes.push_back("E N = " + std::to_string(N) + ": relative Delta n L-inf " + shortest(lv.error) +
                          ", U_iter " + shortest(lv.u_iter) + ", max Delta n / N_A " + shortest(lv.max_ratio));
        levels.push_back(std::move(lv));
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
    r.notes.push_back("E orders (finest pair first):" + orders + "; spread " + shortest(st.order_spread));
    const std::string band = "p-hat in [1.8, 2.3] (p_th 2), Delta n / N_A <= 1e-8 on every level";
    const std::string accuracy = "relative <= 1e-4, Delta n / N_A <= 1e-8, U_ref <= 1e-5";
    const auto inapplicable = [&]() -> std::optional<std::string> {
        if (!predicate) return "low injection: max Delta n / N_A above 1e-8 on a level";
        if (u_ref > 0.1 * e_tolerance) return "U_ref " + shortest(u_ref) + " above a tenth of the tolerance";
        return std::nullopt;
    }();
    ComparisonRecord co{"Delta n order (E)", "1", 4.0, st.observed, 2.0, u_ref, st.observed, band, {}};
    if (inapplicable) {
        co.comparison = Comparison{Outcome::not_evaluated, *inapplicable, {}};
    } else {
        co.comparison = v.judged ? judge(r.info, {}, v.pass, band) : Comparison{Outcome::not_evaluated, v.reason, {}};
    }
    r.comparisons.push_back(std::move(co));
    const Level& fine = levels.front();
    ComparisonRecord ca{"Delta n finest (E, N = 512)", "1", 512.0, fine.error, 0.0, u_ref, fine.error, accuracy, {}};
    ca.comparison = inapplicable ? Comparison{Outcome::not_evaluated, *inapplicable, {}}
                                 : judge(r.info, {}, fine.error <= e_tolerance, accuracy);
    r.comparisons.push_back(std::move(ca));

    // The E comparison's U on the finest mesh: GCI of the finest pair with p-hat, U_iter, U_ref.
    const double p_hat = st.observed;
    const double r_h = levels[1].h / levels[0].h;
    const double u_disc = reference::gci_safety_factor_three_grids * std::abs(levels[1].error - levels[0].error) /
                          (std::pow(r_h, p_hat) - 1.0);
    // The (h, tol) factorial at N = 512 (16.2 (p), DT-Q23 A): tol = 1e-12 and 1e-14, h / 2 at N = 1023. eps_int enters
    // U_E (DECISIONS 4.1 (c)), and above 10 % of it blocks both E passes.
    double fac_e[3] = {0.0, 0.0, 0.0};
    {
        int k = 0;
        for (const auto& [N, tol] : {std::pair{512, 0.01 * tol_newton}, std::pair{1023, tol_newton},
                                     std::pair{1023, 0.01 * tol_newton}}) {
            const auto e = e_error(N, tol, unit, gE, u1, in.n0);
            if (!e) {
                r.run_error = e.error() + " (factorial)";
                return r;
            }
            fac_e[k++] = *e;
        }
    }
    const std::vector<double> parts{u_disc, fine.u_iter, u_ref};
    const Factorial fac = factorial(fine.error, fac_e[0], fac_e[1], fac_e[2], parts);
    r.notes.push_back("E " + fac.note + " (U_disc " + shortest(u_disc) + ")");
    for (ComparisonRecord& c : r.comparisons) block_if_interacting(c, fac);
    const double u_E = fac.ok ? fac.u : u_disc + fine.u_iter + u_ref;

    // A: the injection ladder at N = 512.
    ReferenceInfo a_info = r.info;
    a_info.cls = ReferenceClass::A;
    a_info.approximation_error = "measured: the departure from the linear model on the injection ladder, O(Delta n / N_A) "
                                 "(16.18 (c), DT-Q14)";
    const device::Device d = bar(512);
    std::vector<std::vector<double>> vs;  // v(G) per node, relative to the peak of u_ref / G
    double max_ratio = 0.0, a_iter = 0.0;
    for (const double k : {1.0, 2.0, 4.0, 8.0}) {
        const double G0 = k * gA;
        const auto s = solve_level(d, fine.n_eq, unit, G0, in.n0);
        if (!s) {
            r.run_error = "A ladder G = " + shortest(G0) + ": " + s.error();
            return r;
        }
        std::vector<double> v_k;
        for (std::size_t i = 0; i < s->dn.size(); ++i) v_k.push_back((s->dn[i] - G0 * fine.u_ref[i]) / (G0 * u1));
        max_ratio = std::max(max_ratio, s->max_dn / NA);
        a_iter = std::max(a_iter, s->u_iter);
        r.notes.push_back("A G = " + shortest(k) + " g0: departure L-inf " + shortest(linf(v_k)) +
                          ", max Delta n / N_A " + shortest(s->max_dn / NA) + ", U_iter " + shortest(s->u_iter));
        vs.push_back(std::move(v_k));
    }
    std::vector<double> diff_norm;
    for (std::size_t k = 0; k + 1 < vs.size(); ++k) {
        std::vector<double> dk;
        for (std::size_t i = 0; i < vs[k].size(); ++i) dk.push_back(vs[k + 1][i] - vs[k][i]);
        diff_norm.push_back(linf(dk));
    }
    std::vector<double> a_orders;
    for (std::size_t k = 0; k + 1 < diff_norm.size(); ++k) a_orders.push_back(std::log2(diff_norm[k + 1] / diff_norm[k]));
    const double p_G = a_orders.front();
    std::vector<double> v0_minus_E;
    std::vector<double> vE;
    for (std::size_t i = 0; i < fine.dn.size(); ++i) vE.push_back((fine.dn[i] - gE * fine.u_ref[i]) / (gE * u1));
    for (std::size_t i = 0; i < vs[0].size(); ++i) {
        const double v0 = vs[0][i] - (vs[1][i] - vs[0][i]) / (std::pow(2.0, p_G) - 1.0);
        v0_minus_E.push_back(v0 - vE[i]);
    }
    const double limit_gap = linf(v0_minus_E);
    r.notes.push_back("A orders in G (smallest G first): " + shortest(a_orders[0]) + " " + shortest(a_orders[1]) +
                      "; |v_0 - v_E| " + shortest(limit_gap) + " against U_E " + shortest(u_E) + " (U_disc " +
                      shortest(u_disc) + ", U_iter " + shortest(fine.u_iter) + ", U_ref " + shortest(u_ref) +
                      ", |eps_int| " + shortest(std::abs(fac.eps_int)) + ")");
    const bool a_applies = max_ratio <= a_predicate;
    // CV-Q8 on the ladder: the iterative error below a tenth of the smallest difference it orders.
    const bool a_iterated = a_iter <= 0.1 * diff_norm.front();
    const std::string a_order_band = "departure order in G in [0.8, 1.3]; max Delta n / N_A <= 1e-1";
    const std::string a_limit =
        "G -> 0 limit agrees with the E solve within its U (U_disc + U_iter + U_ref + |eps_int|)";
    ComparisonRecord ao{"departure order in G (A)", "1", 4.0, p_G, 1.0, u_ref, p_G, a_order_band, {}};
    ComparisonRecord al{"G -> 0 limit vs the E solve (A)", "1", 512.0, limit_gap, 0.0, u_ref, limit_gap,
                       a_limit + " = " + shortest(u_E), {}};
    if (!a_applies) {
        const std::string why = "max Delta n / N_A " + shortest(max_ratio) + " above 1e-1";
        ao.comparison = Comparison{Outcome::not_applicable, why, {}};
        al.comparison = Comparison{Outcome::not_applicable, why, {}};
    } else if (!a_iterated) {
        const std::string why = "U_iter " + shortest(a_iter) + " above a tenth of the smallest difference";
        ao.comparison = Comparison{Outcome::not_evaluated, why, {}};
        al.comparison = Comparison{Outcome::not_evaluated, why, {}};
    } else if (inapplicable) {  // the limit is judged against the E comparison's U
        ao.comparison = judge(a_info, {}, p_G >= 0.8 && p_G <= 1.3, a_order_band);
        al.comparison = Comparison{Outcome::not_evaluated, *inapplicable, {}};
    } else {
        ao.comparison = judge(a_info, {}, p_G >= 0.8 && p_G <= 1.3, a_order_band);
        al.comparison = judge(a_info, {}, limit_gap <= u_E, a_limit);
    }
    r.comparisons.push_back(std::move(ao));
    r.comparisons.push_back(std::move(al));
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
