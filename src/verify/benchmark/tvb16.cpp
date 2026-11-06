// TVB-16, the ideal diode: long and short base, and the Sah-Noyce-Shockley recombination current (DECISIONS Unit 55
// §6, 16.2 (d), 16.4; A-class criteria 16.18 (f), 16.19, frozen; DT-Q5, DT-Q16, DT-Q24 A, DT-Q25 A; construction
// 16.31).
//
// Code under test: the drift-diffusion system and sweep_bias on abrupt 1D silicon diodes, N_A = 1e16 cm^-3 (x < 0, the
// anode at x_min) and N_D = 2e16 cm^-3 (x > 0), ohmic contacts, 300 K, Boltzmann statistics, no band-gap narrowing,
// no field mobility, Caughey-Thomas mobility for both carriers (constant per side: uniform doping), Auger and
// radiative off. Three configurations (16.2 (d), 16.31):
//   - (i) SRH off, 5 um per side: only the short-base law applies;
//   - (ii) short: SRH on (Scharfetter lifetimes at the local doping), the same device;
//   - (ii) long: SRH on, 500 um per side (about 3 L_n and 10 L_p), nodes x = L sinh(beta xi) / sinh(beta), beta = 7.5.
// Nodes N = 128, 256, 512, 1024 (N even: x = 0 is a face), xi = (2k - (N - 1)) / (N - 1). Biases: forward 0.2 ... 0.6 V
// by 0.1; reverse -0.5, -1, -2, -5 V in (ii) only (16.19).
// J_sim is the anode's weighted terminal current (Unit 32), with its resolution.
// Reference (verify/reference/diode): Shockley's finite-base law at the QE depletion edges plus, with SRH, the SNS
// integral over the reference's QE profile. Inputs are NiTCAD's material data (I1), printed as the spec values (R-7).
// Criterion (A, 16.19), per bias: |log10(J_sim / J_ref)| <= log10(1 + (U_num + U_ref + eps_J) / J_ref), applicable when
// eps_J <= 0.1 J_ref. eps_J / J_ref (16.18 (f), DT-Q25 A), measured on the N = 1024 state:
//   - max(Delta n / N) at the two reference edges (the minority excess against the doping);
//   - e^v / |e^v - 1| x Delta phi / V_T, Delta phi the larger quasi-Fermi drop between the edges;
//   - short base only: 3 L_D / W' (the edge's definition), the larger side.
// U_num = U_disc (Roache's three-grid GCI on N = 1024, 512, 256, F_s = 1.25, with the measured ratios;
// not estimated outside the asymptotic range: the orders of 1024-512-256 and 512-256-128 within 0.1, 4.2, DT-Q26 A) + the
// resolution + |eps_int| (the (h, tol) factorial at N = 1024: N = 2048, tol = 1e-12 and 1e-14; 16.2 (p)).
// U_ref = the SNS second-method difference + its quadrature estimate (+ the double's rounding).
// Derivatives (DT-Q24 A), between consecutive biases: the ideality n = Delta V / (V_T Delta ln J) (forward) and, in
// (ii), the reverse exponent x = Delta ln|J| / Delta ln(V_bi - V). |x_sim - x_ref| <= U_num,x + U_ref,x + eps_x, with
// eps_x = |x_ref| (eps_k + eps_k+1) / |Delta ln y| (y = J_ref or V_bi - V), U_ref,x the same with the U_ref's, and
// U_num,x the GCI of x, its share of the resolutions and the factorial's |eps_int| of x.
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
#include "NiTCAD/verify/reference/diode.hpp"
#include "NiTCAD/verify/reference/metrics.hpp"
#include "NiTCAD/verify/reference/pn_junction.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

constexpr double NA = 1e16, ND = 2e16, T = 300.0, tol_newton = 1e-12;
constexpr double applicability = 0.1;     // eps_J <= 0.1 J_ref (16.19)
constexpr double legacy_short_base = 1.321e-2;  // bias:144's analytic value at 0.5 V (DT-Q5: printed, not gated)
constexpr int ladder[4] = {1024, 512, 256, 128};  // finest first
constexpr int factorial_N = 2048;

struct Config {
    std::string name;
    double L_cm = 0.0;
    double beta = 0.0;  // 0: uniform nodes
    bool srh = false;
    bool short_base = true;
    std::vector<double> forward, reverse;
};

std::vector<Config> configs() {
    const std::vector<double> fwd{0.2, 0.3, 0.4, 0.5, 0.6}, rev{-0.5, -1.0, -2.0, -5.0};
    return {{"(i) no G-R, short base", 5e-4, 0.0, false, true, fwd, {}},
            {"(ii) SRH, short base", 5e-4, 0.0, true, true, fwd, rev},
            {"(ii) SRH, long base", 5e-2, 7.5, true, false, fwd, rev}};
}

ReferenceInfo info() {
    ReferenceInfo r;
    r.id = "TVB-16";
    r.cls = ReferenceClass::A;
    r.independence = Independence::I1;
    r.provenance = Provenance::PA;
    r.citation = "Shockley, BSTJ 28 (1949) eqs 4.23-4.24 (finite base derived, PD); Sah, Noyce and Shockley, "
                 "Proc. IRE 45 (1957) as the definition q int U dx (PD, DT-Q18)";
    r.version = "DECISIONS Unit 55 16.2 (d), 16.3, 16.18 (f), 16.19, 16.31; DT-Q24 A, DT-Q25 A";
    r.equations = {1, 2, 3, 7};
    r.form = "J = q (e^v - 1) sum D n0 / L coth(W' / L) + q int_{-x_p}^{x_n} U dx at the QE edges";
    r.assumptions = {"low injection (measured)", "quasi-equilibrium across the depletion region (measured)",
                     "Boltzmann statistics", "abrupt doping", "the charge-balance edge (3 L_D / W', short base)"};
    r.derivation = "minority diffusion D u'' = u / tau with the QE law of the junction and an ohmic contact";
    r.verification = "diode_test: short-base law, coth limits, SNS bounds and second methods";
    r.second_method = "tanh-sinh on the SNS integrals, at run time";
    r.approximation_error = "measured per bias: max(Delta n / N) at the edges + e^v / |e^v - 1| Delta phi / V_T + "
                            "3 L_D / W' (short base) (16.18 (f), DT-Q25 A)";
    return r;
}

std::vector<double> nodes(const Config& c, int N) {
    std::vector<double> x(static_cast<std::size_t>(N));
    for (int k = 0; k < N; ++k) {
        const double xi = (2.0 * k - (N - 1)) / (N - 1);
        x[static_cast<std::size_t>(k)] =
            c.beta == 0.0 ? c.L_cm * xi : c.L_cm * std::sinh(c.beta * xi) / std::sinh(c.beta);
    }
    x.front() = -c.L_cm;
    x.back() = c.L_cm;
    return x;
}

device::Device diode(const Config& c, int N) {
    mesh::Mesh m = *mesh::make_tensor_grid(nodes(c, N));
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

solve::BiasOptions options(bool srh, double tol) {
    solve::BiasOptions o;
    o.models.srh = srh;
    o.models.auger = o.models.radiative = false;
    o.models.bgn = false;
    o.models.fermi_dirac = false;
    o.models.field_mobility = false;
    o.models.hole_mobility = assemble::HoleMobility::caughey_thomas;
    o.models.electron_mobility = assemble::ElectronMobility::caughey_thomas;
    o.newton.tol_update = tol;
    o.step_control.max_cuts = 30;
    o.weighted_current = true;
    return o;
}

// One solved bias: the anode's weighted current and resolution, and the measured eps_J / J_ref terms.
struct Point {
    double V = 0.0, J = 0.0, res = 0.0;
    double low = 0.0, phi = 0.0;  // max(Delta n / N) at the edges; e^v / |e^v - 1| Delta phi / V_T
    bool solved = true;           // false: Newton did not reach tol (factorial only)
    double stall = 0.0;           // then the last Newton update
};

// Linear interpolation of f (on ascending x) at s; log-linear when log_scale.
double at(const std::vector<double>& x, const std::vector<double>& f, double s, bool log_scale) {
    const auto it = std::upper_bound(x.begin(), x.end(), s);
    const std::size_t k = std::clamp<std::size_t>(static_cast<std::size_t>(it - x.begin()), 1, x.size() - 1) - 1;
    const double t = (s - x[k]) / (x[k + 1] - x[k]);
    if (log_scale) return std::exp((1.0 - t) * std::log(f[k]) + t * std::log(f[k + 1]));
    return (1.0 - t) * f[k] + t * f[k + 1];
}

using References = std::vector<reference::DiodeCurrent>;  // per bias, forward then reverse

// Every bias of a configuration at N nodes, forward then reverse (sweeps from equilibrium). With may_stall (the
// factorial's tol / 100), a point Newton cannot bring to tol is kept as not solved, with its last update, and the
// sweep restarts at the next bias (DT-Q23 A: recorded, not loosened); otherwise a stop is a run error.
std::expected<std::vector<Point>, std::string> solve_level(const Config& c, int N, double tol, const References& refs,
                                                           double VT, bool may_stall) {
    const device::Device d = diode(c, N);
    const std::size_t n = d.mesh().node_count();
    std::vector<double> x(n);
    for (std::size_t i = 0; i < n; ++i) x[i] = d.mesh().points()[i][0];
    if (!std::is_sorted(x.begin(), x.end())) return std::unexpected("nodes not in x order");
    std::vector<Point> out;
    for (const auto* part : {&c.forward, &c.reverse}) {
        std::size_t k0 = 0;
        while (k0 < part->size()) {
            std::vector<std::vector<double>> points;
            for (std::size_t k = k0; k < part->size(); ++k) points.push_back({(*part)[k], 0.0});
            const auto s = solve::sweep_bias(d, points, options(c.srh, tol));
            if (!s) return std::unexpected(c.name + ", N = " + std::to_string(N) + ": " + s.error().message);
            for (const auto& bp : s->points) {
                const auto& ref = refs[out.size()];
                Point p;
                p.V = (*part)[k0];
                p.J = bp.weighted_current[0];
                p.res = bp.weighted_current_resolution[0];
                const double xp = -ref.edge_p_cm.to_double(), xn = ref.edge_n_cm.to_double();
                const double n_edge = at(x, bp.fields.n_cm3, xp, true), p_edge = at(x, bp.fields.p_cm3, xn, true);
                p.low = std::max(std::abs(n_edge - ref.n_p0_cm3.to_double()) / NA,
                                 std::abs(p_edge - ref.p_n0_cm3.to_double()) / ND);
                const double dn = std::abs(at(x, bp.bands.electron_fermi_eV, xn, false) -
                                           at(x, bp.bands.electron_fermi_eV, xp, false));
                const double dp = std::abs(at(x, bp.bands.hole_fermi_eV, xn, false) -
                                           at(x, bp.bands.hole_fermi_eV, xp, false));
                const double v = p.V / VT;
                p.phi = std::exp(v) / std::abs(std::expm1(v)) * std::max(dn, dp) / VT;
                out.push_back(p);
                ++k0;
            }
            if (!s->stopped) break;
            if (!may_stall) {
                return std::unexpected(c.name + ", N = " + std::to_string(N) + ": " + s->stopped->message);
            }
            Point p;
            p.V = (*part)[k0];
            p.solved = false;
            p.J = p.res = std::nan("");  // the factorial then reports "not evaluated"
            p.stall = INFINITY;
            if (s->unfinished && !s->unfinished->iterations.empty()) p.stall = s->unfinished->iterations.back().update;
            out.push_back(p);
            ++k0;
        }
    }
    return out;
}

// Roache's three-grid GCI (F_s = 1.25) with non-constant ratios, finest first; absolute. nullopt: not estimated
// (oscillatory or a non-positive order). Sets the order.
std::optional<double> gci(const double h[3], const double f[3], double& order) {
    const double e21 = f[1] - f[0], e32 = f[2] - f[1];
    if (e21 == 0.0 && e32 == 0.0) {
        order = 0.0;
        return 0.0;
    }
    if (!(e32 / e21 > 0.0)) return std::nullopt;
    const double r21 = h[1] / h[0], r32 = h[2] / h[1];
    double p = 2.0;
    for (int it = 0; it < 200; ++it) {
        const double q = std::log((std::pow(r21, p) - 1.0) / (std::pow(r32, p) - 1.0));
        p = std::abs(std::log(std::abs(e32 / e21)) + q) / std::log(r21);
    }
    order = p;
    if (!(std::isfinite(p) && p > 0.0)) return std::nullopt;
    return reference::gci_safety_factor_three_grids * std::abs(e21) / (std::pow(r21, p) - 1.0);
}

// U_disc with §4.2's asymptotic check (DT-Q26 A): the three-grid orders of levels 0-1-2 and 1-2-3 agree within 0.1,
// else U_disc is not estimated (why says which). Sets both orders.
std::optional<double> disc(const double h[4], const double f[4], double& order, double& order2, std::string& why) {
    const auto u = gci(h, f, order);
    const auto u2 = gci(h + 1, f + 1, order2);
    if (!u) {
        why = "GCI: oscillatory or a non-positive order";
        return std::nullopt;
    }
    if (!u2) {
        why = "asymptotic check: no order on the coarser triple";
        return std::nullopt;
    }
    if (std::abs(order - order2) > reference::asymptotic_order_spread) {
        why = "outside the asymptotic range (orders " + shortest(order) + ", " + shortest(order2) + ")";
        return std::nullopt;
    }
    return u;
}

}  // namespace

BenchmarkResult run_tvb16(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-16";
    r.title = "Ideal diode: long and short base, Sah-Noyce-Shockley current";
    r.info = info();
    if (tier == Tier::fast) {
        r.notes.push_back("full tier only (16.7)");
        return r;
    }
    const physics::Semiconductor si = physics::silicon();
    const double VT = base::thermal_voltage(T);
    const double ni = physics::intrinsic_density(si, T);
    const reference::PnJunctionInput junction{.q_C = base::q_C,
                                              .eps_F_per_cm = si.parameters().eps_r * base::eps0_F_per_cm,
                                              .thermal_voltage_V = VT,
                                              .ni_cm3 = ni,
                                              .acceptors_cm3 = NA,
                                              .donors_cm3 = ND};
    const double Dn = physics::caughey_thomas_mobility(si, physics::Carrier::electron, NA, T) * VT;
    const double Dp = physics::caughey_thomas_mobility(si, physics::Carrier::hole, ND, T) * VT;
    const auto tau = [&](physics::Carrier c, double N) { return physics::scharfetter_lifetime(si, c, N); };
    r.notes.push_back("spec (NiTCAD's material data, I1): n_i " + shortest(ni) + " cm^-3, D_n(N_A) " + shortest(Dn) +
                      ", D_p(N_D) " + shortest(Dp) + " cm^2/s; tau_n, tau_p " +
                      shortest(tau(physics::Carrier::electron, NA)) + ", " + shortest(tau(physics::Carrier::hole, NA)) +
                      " s (p side), " + shortest(tau(physics::Carrier::electron, ND)) + ", " +
                      shortest(tau(physics::Carrier::hole, ND)) + " s (n side)");
    const auto eq = reference::pn_junction(junction);
    if (!eq) {
        r.run_error = "reference: " + eq.error();
        return r;
    }
    const double Vbi = eq->built_in_V.to_double();

    // DT-Q5: bias:144's short-base value (1.321e-2 at 0.5 V, legacy silicon, 1e17 / 1e17, 1 um sides, W(0) edges),
    // against this reference's W(V) form on the same inputs. Printed only.
    {
        const physics::Semiconductor lg = physics::legacy_silicon();
        reference::DiodeInput li;
        li.junction = {.q_C = base::q_C, .eps_F_per_cm = 11.7 * base::eps0_F_per_cm, .thermal_voltage_V = VT,
                       .ni_cm3 = physics::intrinsic_density(lg, T), .acceptors_cm3 = 1e17, .donors_cm3 = 1e17};
        li.p_side = {.length_cm = 1e-4,
                     .minority_diffusivity_cm2_s =
                         physics::caughey_thomas_mobility(lg, physics::Carrier::electron, 1e17, T) * VT};
        li.n_side = {.length_cm = 1e-4,
                     .minority_diffusivity_cm2_s =
                         physics::caughey_thomas_mobility(lg, physics::Carrier::hole, 1e17, T) * VT};
        const auto lj = reference::diode_current(li, 0.5);
        r.notes.push_back(lj ? "DT-Q5 cross-check: the short-base law with W(V) edges on bias:144's inputs gives " +
                                   shortest(lj->total_A_per_cm2.to_double()) +
                                   " A/cm^2 at 0.5 V; bias:144's W(0) form " +
                                   shortest(legacy_short_base)
                             : "DT-Q5 cross-check: " + lj.error());
    }

    for (const Config& c : configs()) {
        reference::DiodeInput in;
        in.junction = junction;
        in.p_side = {.length_cm = c.L_cm, .minority_diffusivity_cm2_s = Dn,
                     .tau_n_s = tau(physics::Carrier::electron, NA), .tau_p_s = tau(physics::Carrier::hole, NA)};
        in.n_side = {.length_cm = c.L_cm, .minority_diffusivity_cm2_s = Dp,
                     .tau_n_s = tau(physics::Carrier::electron, ND), .tau_p_s = tau(physics::Carrier::hole, ND)};
        in.srh = c.srh;
        std::vector<double> biases = c.forward;
        biases.insert(biases.end(), c.reverse.begin(), c.reverse.end());
        References refs;
        for (const double V : biases) {
            const auto d = reference::diode_current(in, V);
            if (!d) {
                r.run_error = c.name + ", V = " + shortest(V) + ": reference: " + d.error();
                return r;
            }
            refs.push_back(*d);
        }
        const auto& r0 = refs.front();
        r.notes.push_back(c.name + ": diffusion lengths " + shortest(r0.diffusion_length_p_cm.to_double()) + " (p), " +
                          shortest(r0.diffusion_length_n_cm.to_double()) + " (n) cm; Debye lengths " +
                          shortest(r0.debye_p_cm.to_double()) + ", " + shortest(r0.debye_n_cm.to_double()) + " cm");

        // The ladder (finest first) and the factorial's three extra states.
        std::vector<std::vector<Point>> lv;
        double hs[4];
        for (int k = 0; k < 4; ++k) {
            const auto s = solve_level(c, ladder[k], tol_newton, refs, VT, false);
            if (!s) {
                r.run_error = s.error();
                return r;
            }
            lv.push_back(*s);
            hs[k] = 2.0 / (ladder[k] - 1);  // in xi
        }
        std::vector<std::vector<Point>> fac;  // (h, tol/100), (h/2, tol), (h/2, tol/100)
        for (const auto& [N, tol] : {std::pair{ladder[0], 0.01 * tol_newton}, std::pair{factorial_N, tol_newton},
                                     std::pair{factorial_N, 0.01 * tol_newton}}) {
            const auto s = solve_level(c, N, tol, refs, VT, tol < tol_newton);
            if (!s) {
                r.run_error = s.error() + " (factorial)";
                return r;
            }
            fac.push_back(*s);
        }

        const std::size_t nb = biases.size();
        std::vector<double> eps(nb), uref(nb), unum(nb);
        std::vector<bool> applies(nb), estimated(nb);
        for (std::size_t b = 0; b < nb; ++b) {
            const double V = biases[b], Jref = refs[b].total_A_per_cm2.to_double();
            const Point& fine = lv[0][b];
            const double edge =
                c.short_base ? 3.0 * std::max(refs[b].debye_p_cm.to_double() / refs[b].base_p_cm.to_double(),
                                              refs[b].debye_n_cm.to_double() / refs[b].base_n_cm.to_double())
                             : 0.0;
            eps[b] = fine.low + fine.phi + edge;
            uref[b] = (refs[b].sns_second_method_difference + refs[b].sns_quadrature_error) / std::abs(Jref) + 1e-15;
            const double f[4] = {lv[0][b].J, lv[1][b].J, lv[2][b].J, lv[3][b].J};
            double p = 0.0, p2 = 0.0;
            std::string why_d;
            const auto ud = disc(hs, f, p, p2, why_d);
            const auto e = [&](double J) { return J / Jref - 1.0; };
            const std::vector<double> parts{ud ? *ud / std::abs(Jref) : 0.0, fine.res / std::abs(Jref)};
            const Factorial fc = factorial(e(fine.J), e(fac[0][b].J), e(fac[1][b].J), e(fac[2][b].J), parts);
            std::string stalled;
            for (int k = 0; k < 3; ++k) {
                if (!fac[static_cast<std::size_t>(k)][b].solved) {
                    stalled += std::string(k == 0 ? " (h, tol/100)" : k == 1 ? " (h/2, tol)" : " (h/2, tol/100)") +
                               " not reached, last update " + shortest(fac[static_cast<std::size_t>(k)][b].stall);
                }
            }
            estimated[b] = ud.has_value() && fc.ok;
            unum[b] = fc.ok ? fc.u : 0.0;
            applies[b] = eps[b] <= applicability;
            const double metric = Jref * fine.J > 0.0 ? std::abs(std::log10(fine.J / Jref)) : INFINITY;
            const double bound = std::log10(1.0 + unum[b] + uref[b] + eps[b]);
            const std::string at_v = c.name + ", V = " + shortest(V) + " V";
            r.notes.push_back(at_v + ": J_ref " + shortest(Jref) + " (SNS " +
                              shortest(refs[b].sns_A_per_cm2.to_double()) +
                              "), J_sim " + shortest(fine.J) + "; eps_J/J: low injection " + shortest(fine.low) +
                              ", quasi-Fermi " + shortest(fine.phi) + ", edge " + shortest(edge) + "; U_num/J " +
                              shortest(unum[b]) + " (GCI orders " + shortest(p) + ", " + shortest(p2) + ", resolution/J " +
                              shortest(fine.res / std::abs(Jref)) + "), U_ref/J " + shortest(uref[b]) + "; " + fc.note +
                              stalled);
            ComparisonRecord cr{"J " + at_v, "A/cm^2", V, fine.J, Jref, uref[b] * std::abs(Jref), metric,
                                "log error <= log10(1 + (U_num + U_ref + eps_J) / J_ref) = " + shortest(bound), {}};
            const std::vector<Predicate> pr{{"eps_J / J_ref <= 0.1", eps[b], applicability, applies[b]}};
            const std::string why = "U_num not estimated (" + (ud ? std::string("the (h, tol) factorial") : why_d) +
                                    (stalled.empty() ? std::string() : ":" + stalled) + ")";
            cr.comparison = !estimated[b] ? Comparison{Outcome::not_evaluated, why, {}}
                                          : judge(r.info, pr, metric <= bound, cr.criterion);
            r.comparisons.push_back(std::move(cr));
        }

        // Derivatives between consecutive biases (DT-Q24 A): forward ideality, reverse exponent.
        const auto pairs = [&](std::size_t first, std::size_t count, bool ideality) {
            for (std::size_t b = first; b + 1 < first + count; ++b) {
                const std::size_t a = b, z = b + 1;
                const double Va = biases[a], Vz = biases[z];
                const auto x_of = [&](double Ja, double Jz) {
                    return ideality ? (Vz - Va) / (VT * (std::log(Jz) - std::log(Ja)))
                                    : (std::log(std::abs(Jz)) - std::log(std::abs(Ja))) /
                                          (std::log(Vbi - Vz) - std::log(Vbi - Va));
                };
                const double Ja = refs[a].total_A_per_cm2.to_double(), Jz = refs[z].total_A_per_cm2.to_double();
                const double x_ref = x_of(Ja, Jz);
                const double dly = ideality ? std::abs(std::log(Jz / Ja))
                                            : std::abs(std::log((Vbi - Vz) / (Vbi - Va)));
                const double share = std::abs(x_ref) / dly;
                const double eps_x = share * (eps[a] + eps[z]), uref_x = share * (uref[a] + uref[z]);
                const double x_sim = x_of(lv[0][a].J, lv[0][z].J);
                const double f[4] = {x_sim, x_of(lv[1][a].J, lv[1][z].J), x_of(lv[2][a].J, lv[2][z].J),
                                     x_of(lv[3][a].J, lv[3][z].J)};
                double p = 0.0, p2 = 0.0;
                std::string why_d;
                const auto ud = disc(hs, f, p, p2, why_d);
                const double res_x =
                    share * (lv[0][a].res / std::abs(lv[0][a].J) + lv[0][z].res / std::abs(lv[0][z].J));
                const std::vector<double> parts{ud ? *ud : 0.0, res_x};
                const Factorial fc = factorial(x_sim - x_ref, x_of(fac[0][a].J, fac[0][z].J) - x_ref,
                                               x_of(fac[1][a].J, fac[1][z].J) - x_ref,
                                               x_of(fac[2][a].J, fac[2][z].J) - x_ref, parts);
                const double bound = (fc.ok ? fc.u : 0.0) + uref_x + eps_x;
                const double metric = std::abs(x_sim - x_ref);
                const std::string name = std::string(ideality ? "ideality " : "reverse exponent ") + c.name + ", V = " +
                                         shortest(Va) + " .. " + shortest(Vz) + " V";
                r.notes.push_back(name + ": x_ref " + shortest(x_ref) + ", x_sim " + shortest(x_sim) + ", U_num " +
                                  shortest(fc.ok ? fc.u : 0.0) + " (GCI orders " + shortest(p) + ", " + shortest(p2) + "), U_ref " +
                                  shortest(uref_x) + ", eps_x " + shortest(eps_x) + "; " + fc.note);
                ComparisonRecord cr{name, "1", Va, x_sim, x_ref, uref_x, metric,
                                    "|x_sim - x_ref| <= U_num + U_ref + eps_x = " + shortest(bound), {}};
                const bool both = applies[a] && applies[z];
                const std::vector<Predicate> pr{
                    {"eps_J / J_ref <= 0.1 at both biases", std::max(eps[a], eps[z]), applicability, both}};
                cr.comparison = !(ud && fc.ok && estimated[a] && estimated[z])
                                    ? Comparison{Outcome::not_evaluated,
                                                 "U_num not estimated (" + (ud ? std::string("factorial or a bias") : why_d) + ")", {}}
                                    : judge(r.info, pr, metric <= bound, cr.criterion);
                r.comparisons.push_back(std::move(cr));
            }
        };
        pairs(0, c.forward.size(), true);
        if (!c.reverse.empty()) pairs(c.forward.size(), c.reverse.size(), false);
    }
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
