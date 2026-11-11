// TVB-08, the equilibrium abrupt heterojunction (DECISIONS Unit 55 6, 17.2 (a)-(c), 17.4; criteria 17.8, frozen
// before code).
//
// Code under test: the heterointerface in both assemblers (the band shift and the permittivity step on the edge
// joining two materials), the equilibrium solve, the drift-diffusion solve at zero bias with the default interface
// and with thermionic emission declared.
// Construction, fixed before any run:
//   - two materials on [-L, L], L = 0.5 um, the interface at x = 0; ohmic contacts at both ends; 300 K; Boltzmann
//     statistics, no band-gap narrowing, no recombination;
//   - nodes x_k = L sinh(beta xi_k) / sinh(beta), xi_k = (2k - (N - 1)) / (N - 1), beta = 3, N even: an odd map, so
//     x = 0 is the midpoint of the middle edge at every level (17.2 (b)); the spacing at the interface is 0.30 of the
//     uniform one;
//   - three junctions, each side's band data given here (the spec; R-7), with a temperature-independent gap (no
//     Varshni term) and N_c, N_v at 300 K, so no NiTCAD band function lies between the spec and either the device or
//     the reference:
//       n-N: narrow n 1e17 | wide n 1e17;   p-N: narrow p 1e17 | wide n 1e17;
//       chi-eps: silicon n 1e17 | silicon with chi + 0.25 eV and eps_r 13, n 1e17;
//     narrow: eps_r 12.9, E_g 1.424 eV, N_c 4.7e17, N_v 7.0e18, chi 4.07 eV; wide: eps_r 12.05, E_g 1.798 eV,
//     N_c 5.4e17, N_v 8.7e18, chi 3.74 eV (of the magnitude of GaAs and Al0.3Ga0.7As); silicon: eps_r 11.7, E_g
//     1.12 eV, N_c 2.86e19, N_v 3.10e19, chi 4.05 eV.
// Reference (verify/reference/heterojunction.hpp): the first integrals matched by D continuity under the electron
// affinity rule, in double-double; cross-checked at run time against Chebyshev collocation of the finite-domain
// problem (I3). Inputs: the spec above and NiTCAD's q, eps_0 and V_T (the model's constants).
// Criteria (17.8):
//   - D, at every level: the three equilibria (solve_equilibrium; solve_bias at 0 V; solve_bias at 0 V with
//     thermionic emission declared) agree, max |Delta psi| <= 64 eps max(max |psi|, V_T) + V_T (u_a + u_b), u the
//     last Newton update of each solve (C1's form, 11.13); psi is taken relative to the left contact node;
//   - E (a): psi L-infinity (relative to |V_bi|) of the zero-bias drift-diffusion solve, p-hat in [1.8, 2.3];
//   - E (b), finest level: psi L-infinity <= 1e-4 |V_bi|; the interface level within 1e-4 |V_bi|; each side's field
//     (one D gives both: one comparison) and each side's band bending (the spike and the notch) relative <= 1e-3.
//     The interface level is the permittivity-weighted mean of the two interface nodes and D the harmonic-mean
//     permittivity times the middle edge's field (exact for the two linear half-edges of the box method).
// Ladder: N = 64, 128, 256, 512. Full tier only (17.7).
// Reference uncertainty: the domain truncation V_T t at the contacts, relative to |V_bi|; above 1e-3 of the psi
// criterion the E comparisons are not_evaluated (TB-Q4).
// The (h, tol) factorial (16.2 (p)) at N = 512: tol = 1e-12 and 1e-14, N = 512 and 1024, on the psi L-infinity
// error; eps_int above 10 % of U blocks the E passes (4.1 (c)).
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <expected>
#include <string>
#include <utility>
#include <vector>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/device/device.hpp"
#include "NiTCAD/mesh/tensor_grid.hpp"
#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/solve/bias.hpp"
#include "NiTCAD/solve/equilibrium.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/heterojunction.hpp"
#include "NiTCAD/verify/reference/metrics.hpp"
#include "NiTCAD/verify/study/study.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

constexpr double L = 5e-5, beta = 3.0, T = 300.0;
constexpr double eps = 0x1p-52;
constexpr double psi_tolerance = 1e-4, level_tolerance = 1e-4, field_tolerance = 1e-3;
constexpr double tol_newton = 1e-12;

// One material of the spec.
struct Spec {
    bool silicon;  // which parameter set carries the fields the equilibrium does not read
    double eps_r, gap_eV, Nc, Nv, chi_eV;
};
constexpr Spec narrow{false, 12.9, 1.424, 4.7e17, 7.0e18, 4.07};
constexpr Spec wide{false, 12.05, 1.798, 5.4e17, 8.7e18, 3.74};
constexpr Spec silicon_a{true, 11.7, 1.12, 2.86e19, 3.10e19, 4.05};
constexpr Spec silicon_b{true, 13.0, 1.12, 2.86e19, 3.10e19, 4.30};

struct Case {
    const char* name;
    Spec left, right;
    double C_left, C_right;  // net doping, donors positive
};
constexpr Case cases[] = {{"n-N", narrow, wide, 1e17, 1e17},
                          {"p-N", narrow, wide, -1e17, 1e17},
                          {"chi-eps", silicon_a, silicon_b, 1e17, 1e17}};

physics::SemiconductorParameters parameters(const Spec& s) {
    physics::SemiconductorParameters p =
        s.silicon ? physics::legacy_silicon_parameters : physics::gallium_arsenide_parameters;
    p.eps_r = s.eps_r;
    p.Eg0_eV = s.gap_eV;
    p.varshni_alpha_eV_per_K = 0.0;
    p.varshni_beta_K = 0.0;
    p.Nc300 = s.Nc;
    p.Nv300 = s.Nv;
    p.electron_affinity_eV = s.chi_eV;
    return p;
}

reference::HeteroSide side(const Spec& s, double C) {
    return {s.eps_r * base::eps0_F_per_cm, s.chi_eV, s.gap_eV, s.Nc, s.Nv, C};
}

ReferenceInfo info() {
    ReferenceInfo r;
    r.id = "TVB-08";
    r.cls = ReferenceClass::E;
    r.independence = Independence::I1;
    r.provenance = Provenance::PB;
    r.citation = "Electron affinity rule (Anderson, Solid-State Electronics 5, 341 (1962), read as E. F. Schubert, "
                 "Band diagrams of heterostructures, eqs (17.1), (17.2)); first integrals of the 1D Poisson-Boltzmann "
                 "equation on each side matched by D continuity (derived, DECISIONS Unit 55 17.13)";
    r.version = "DECISIONS Unit 55 17.3, 17.8";
    r.equations = {1, 4};
    r.form = "eps_1 Phi_1(t_1) = eps_2 Phi_2(t_1 - V_bi / V_T), Phi_b = n_b G(t) + p_b G(-t), "
             "V_bi = [V_T ln(n_b / N_c) - chi]_2 - [V_T ln(n_b / N_c) - chi]_1";
    r.assumptions = {"thermal equilibrium", "Boltzmann statistics", "abrupt interface on a dual-cell face",
                     "no interface charge or dipole", "infinite-domain reference (truncation measured)"};
    r.derivation = "verify/reference/heterojunction.hpp";
    r.verification = "D: three equilibria agree; E: psi order and finest-level accuracy of psi, the interface level, "
                     "the fields and the band bendings";
    r.second_method = "Chebyshev collocation of the finite-domain problem with the permittivity step (I3), at run time";
    return r;
}

std::vector<double> nodes(int N) {
    std::vector<double> x;
    const double s = std::sinh(beta);
    for (int k = 0; k < N; ++k) x.push_back(L * std::sinh(beta * (2.0 * k - (N - 1)) / (N - 1)) / s);
    return x;
}

std::expected<device::Device, std::string> junction(const Case& c, int N, bool thermionic) {
    auto m = mesh::make_tensor_grid(nodes(N));
    if (!m) return std::unexpected(m.error().message);
    const std::size_t n = m->node_count();
    std::vector<double> donors(n, 0.0), acceptors(n, 0.0);
    std::vector<device::RegionId> region(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        const bool left = m->points()[i][0] < 0.0;
        region[i] = left ? 0 : 1;
        const double C = left ? c.C_left : c.C_right;
        (C > 0.0 ? donors[i] : acceptors[i]) = std::abs(C);
    }
    const auto a = physics::Semiconductor::create(parameters(c.left));
    const auto b = physics::Semiconductor::create(parameters(c.right));
    if (!a || !b) return std::unexpected("material: " + (!a ? a.error().message : b.error().message));
    auto lo = m->find_boundary("x_min")->nodes;
    auto hi = m->find_boundary("x_max")->nodes;
    std::vector<device::Interface> interfaces;
    if (thermionic) interfaces.push_back({"left", "right", device::InterfaceTransport::thermionic_emission});
    auto d = device::Device::create({.mesh = std::move(*m),
                                     .temperature_K = T,
                                     .regions = {{"left", *a}, {"right", *b}},
                                     .node_region = std::move(region),
                                     .donors = std::move(donors),
                                     .acceptors = std::move(acceptors),
                                     .contacts = {{"anode", device::ContactKind::ohmic, std::move(lo)},
                                                  {"cathode", device::ContactKind::ohmic, std::move(hi)}},
                                     .interfaces = std::move(interfaces)});
    if (!d) return std::unexpected(d.error().message);
    return std::move(*d);
}

assemble::PhysicsModels models() {
    assemble::PhysicsModels m;
    m.srh = m.auger = m.radiative = false;
    m.bgn = false;
    m.fermi_dirac = false;
    return m;
}

struct Solved {
    std::vector<double> psi;  // relative to the left contact node [V]
    double update = 0.0;      // the last Newton update (scaled)
};

std::expected<Solved, std::string> last(const std::vector<double>& potential,
                                        const results::ConvergenceRecord& convergence) {
    Solved s;
    s.psi = potential;
    for (double& v : s.psi) v -= potential.front();
    s.update = convergence.iterations.empty() ? 0.0 : convergence.iterations.back().update;
    return s;
}

std::expected<Solved, std::string> drift_diffusion(const Case& c, int N, bool thermionic, double tol) {
    const auto d = junction(c, N, thermionic);
    if (!d) return std::unexpected(d.error());
    solve::BiasOptions o;
    o.models = models();
    o.newton.tol_update = tol;
    // Unit 50's minority floor: at 1e-12 the wide-gap minorities (about 1e-5 cm^-3) never converge in the default
    // measure (found at the first run, before any result: "Newton did not converge in 100 iterations").
    o.minority_floor = true;
    const auto s = solve::solve_bias(*d, std::vector<double>{0.0, 0.0}, o);
    if (!s) return std::unexpected("solve_bias: " + s.error().message);
    return last(s->fields.potential_V, s->convergence);
}

std::expected<Solved, std::string> equilibrium(const Case& c, int N, double tol) {
    const auto d = junction(c, N, false);
    if (!d) return std::unexpected(d.error());
    solve::EquilibriumOptions o;
    o.models = models();
    o.newton.tol_update = tol;
    const auto s = solve::solve_equilibrium(*d, o);
    if (!s) return std::unexpected("solve_equilibrium: " + s.error().message);
    return last(s->fields.potential_V, s->convergence);
}

struct Level {
    int N;
    double h, psi_linf, level, field, bend_left, bend_right, u_iter, tail;
    double agree, agree_tol;  // D: the three equilibria
};

}  // namespace

BenchmarkResult run_tvb08(Tier tier) {
    BenchmarkResult r;
    r.id = "TVB-08";
    r.title = "Equilibrium abrupt heterojunction";
    r.info = info();
    if (tier < Tier::full) return r;
    const double VT = base::thermal_voltage(T);
    const auto record = [&](const std::string& quantity, double parameter, double sim, double refv, double uref,
                            double metric, const std::string& crit, Comparison cmp) -> ComparisonRecord& {
        ComparisonRecord rec{quantity, "1", parameter, sim, refv, uref, metric, crit, {}};
        rec.comparison = std::move(cmp);
        r.comparisons.push_back(std::move(rec));
        return r.comparisons.back();
    };
    for (const Case& c : cases) {
        const std::string at = c.name;
        const auto ref = reference::heterojunction(
            {base::q_C, VT, side(c.left, c.C_left), side(c.right, c.C_right)});
        if (!ref) {
            r.run_error = at + ": reference: " + ref.error();
            return r;
        }
        const double vbi = ref->built_in_V.to_double(), scale = std::abs(vbi);
        const double level_ref = ref->interface_V.to_double();
        const double D_ref = ref->displacement_C_per_cm2.to_double();
        r.notes.push_back(at + ": reference V_bi " + shortest(vbi) + " V, interface level " + shortest(level_ref) +
                          " V, D " + shortest(D_ref) + " C/cm^2, fields " +
                          shortest(ref->field_left_V_per_cm.to_double()) + " and " +
                          shortest(ref->field_right_V_per_cm.to_double()) + " V/cm, E_c - E_F at the interface " +
                          shortest(ref->ec_left_eV.to_double()) + " and " + shortest(ref->ec_right_eV.to_double()) +
                          " eV, Debye lengths " + shortest(ref->debye_left_cm.to_double()) + " and " +
                          shortest(ref->debye_right_cm.to_double()) + " cm, n_i " + shortest(ref->ni_left_cm3) +
                          " and " + shortest(ref->ni_right_cm3) + " cm^-3");
        // The I3 cross-check on the coarsest grid's nodes.
        {
            const std::vector<double> x = nodes(64);
            const auto fi = reference::hetero_potential(*ref, x);
            const auto col = reference::hetero_collocation(*ref, -L, L, 8, 48, x);
            if (!fi || !col) {
                r.run_error = at + ": reference cross-check: " + (!fi ? fi.error() : col.error());
                return r;
            }
            double worst = 0.0;
            for (std::size_t i = 0; i < x.size(); ++i) {
                worst = std::max(worst, std::abs((*fi)[i].phi_V.to_double() - col->phi_V[i]));
            }
            const double tail = VT * std::max((*fi)[0].tail.to_double(), (*fi)[x.size() - 1].tail.to_double());
            const double allowed = tail + col->tail_estimate_V + 1e-10 * std::max(scale, VT);
            r.notes.push_back(at + ": I3: collocation (8 + 8 pieces, degree 48, " +
                              std::to_string(col->newton_iterations) + " Newton steps, tail " +
                              shortest(col->tail_estimate_V) + " V) against the first integral: max |difference| " +
                              shortest(worst) + " V; truncation V_T t at the contacts " + shortest(tail) + " V");
            if (!(worst <= allowed)) {
                r.run_error = at + ": reference cross-check: collocation and first integral differ by " +
                              shortest(worst) + " V, more than " + shortest(allowed) + " V";
                return r;
            }
        }
        const double e1 = c.left.eps_r, e2 = c.right.eps_r;
        // The psi L-infinity error of the zero-bias drift-diffusion solve at (N, tol), relative to |V_bi|.
        const auto measure = [&](int N, double tol, Level* out) -> std::expected<double, std::string> {
            const auto dd = drift_diffusion(c, N, false, tol);
            if (!dd) return std::unexpected(at + ", N = " + std::to_string(N) + ": " + dd.error());
            const std::vector<double> x = nodes(N);
            const auto exact = reference::hetero_potential(*ref, x);
            if (!exact) return std::unexpected(at + ": reference: " + exact.error());
            const std::size_t n = x.size();
            double worst = 0.0;
            for (std::size_t i = 0; i < n; ++i) {
                worst = std::max(worst, std::abs(dd->psi[i] - (*exact)[i].phi_V.to_double()));
            }
            if (out != nullptr) {
                const std::size_t k = static_cast<std::size_t>(N / 2 - 1);
                const double level = (e1 * dd->psi[k] + e2 * dd->psi[k + 1]) / (e1 + e2);
                const double D = base::eps0_F_per_cm * (2.0 * e1 * e2 / (e1 + e2)) *
                                 std::abs(dd->psi[k + 1] - dd->psi[k]) / (x[k + 1] - x[k]);
                const double total = dd->psi[n - 1];
                out->N = N;
                out->h = 2.0 * L / (N - 1);
                out->psi_linf = worst / scale;
                out->level = std::abs(level - level_ref) / scale;
                out->field = std::abs(D - D_ref) / D_ref;
                out->bend_left = std::abs(level - level_ref) / std::abs(level_ref);
                out->bend_right = std::abs((total - level) - (vbi - level_ref)) / std::abs(vbi - level_ref);
                out->u_iter = dd->update * VT / scale;
                out->tail = VT * std::max((*exact)[0].tail.to_double(), (*exact)[n - 1].tail.to_double()) / scale;
                // D: the three equilibria.
                const auto eq = equilibrium(c, N, tol);
                const auto te = drift_diffusion(c, N, true, tol);
                if (!eq || !te) {
                    return std::unexpected(at + ", N = " + std::to_string(N) + ": " + (!eq ? eq.error() : te.error()));
                }
                double a = 0.0, b = 0.0, top = VT;
                for (std::size_t i = 0; i < n; ++i) {
                    a = std::max(a, std::abs(dd->psi[i] - eq->psi[i]));
                    b = std::max(b, std::abs(dd->psi[i] - te->psi[i]));
                    top = std::max(top, std::abs(dd->psi[i]));
                }
                out->agree = std::max(a, b);
                out->agree_tol = 64.0 * eps * top + VT * (dd->update + std::max(eq->update, te->update));
                r.notes.push_back(at + ", N = " + std::to_string(N) + ": psi Linf " + shortest(out->psi_linf) +
                                  ", interface level " + shortest(level) + " V, D " + shortest(D) +
                                  " C/cm^2 (relative error " + shortest(out->field) + "), V_bi " + shortest(total) +
                                  " V; equilibria differ by " + shortest(a) + " V (Poisson) and " + shortest(b) +
                                  " V (thermionic), last updates " + shortest(dd->update) + ", " +
                                  shortest(eq->update) + ", " + shortest(te->update) + "; truncation " +
                                  shortest(out->tail));
            }
            return worst / scale;
        };
        std::vector<Level> levels;  // finest first
        for (const int N : {512, 256, 128, 64}) {
            Level lv{};
            const auto e = measure(N, tol_newton, &lv);
            if (!e) {
                r.run_error = e.error();
                return r;
            }
            levels.push_back(lv);
        }
        // D, coarsest first.
        for (auto l = levels.rbegin(); l != levels.rend(); ++l) {
            const std::string crit = "max |Delta psi| <= 64 eps max(|psi|, V_T) + V_T (u_a + u_b) = " +
                                     shortest(l->agree_tol);
            record("equilibria agree, " + at + ", N = " + std::to_string(l->N), l->h, l->agree, 0.0, 0.0,
                   l->agree / l->agree_tol, crit, judge(r.info, {}, l->agree <= l->agree_tol, crit));
        }
        // E.
        double u_ref = 0.0;
        std::vector<double> hs, es, us;
        for (const Level& l : levels) {
            hs.push_back(l.h);
            es.push_back(l.psi_linf);
            us.push_back(l.u_iter);
            u_ref = std::max(u_ref, l.tail);
        }
        const bool ref_ok = u_ref <= 1e-3 * psi_tolerance;
        const study::OrderStudy st = study::order_study(hs, es, us);
        const study::OrderVerdict v = study::judge_order(st, {.p_th = 2.0});
        std::string orders;
        for (const double p : st.orders) orders += " " + shortest(p);
        r.notes.push_back(at + ": psi Linf orders (finest pair first):" + orders + "; spread " +
                          shortest(st.order_spread) + "; reference uncertainty (truncation) " + shortest(u_ref));
        // The (h, tol) factorial at N = 512.
        const double u_disc = reference::gci_safety_factor_three_grids * std::abs(levels[1].psi_linf - levels[0].psi_linf) /
                              (std::pow(levels[1].h / levels[0].h, st.observed) - 1.0);
        std::array<double, 3> fac_e{};
        {
            std::size_t k = 0;
            for (const auto& [N, tol] : {std::pair{512, 0.01 * tol_newton}, std::pair{1024, tol_newton},
                                         std::pair{1024, 0.01 * tol_newton}}) {
                const auto e = measure(N, tol, nullptr);
                if (!e) {
                    r.run_error = e.error() + " (factorial)";
                    return r;
                }
                fac_e[k++] = *e;
            }
        }
        const std::vector<double> parts{u_disc, levels.front().u_iter, u_ref};
        const Factorial fac = factorial(levels.front().psi_linf, fac_e[0], fac_e[1], fac_e[2], parts);
        r.notes.push_back(at + ": " + fac.note + " (U_disc " + shortest(u_disc) + ")");
        const std::string missing = "a reference within 1e-3 of the criterion (truncation)";
        const std::string band = "p-hat in [1.8, 2.3] (p_th 2)";
        {
            Comparison cmp = !ref_ok    ? missing_input(missing)
                             : v.judged ? judge(r.info, {}, v.pass, band)
                                        : Comparison{Outcome::not_evaluated, v.reason, {}};
            ComparisonRecord& rec =
                record("psi Linf order, " + at, hs.front(), st.observed, 2.0, u_ref, st.observed, band, std::move(cmp));
            block_if_interacting(rec, fac);
        }
        const Level& fin = levels.front();
        const auto finest = [&](const std::string& q, double metric, double tol, const std::string& of) {
            const std::string crit = "relative to " + of + " <= " + shortest(tol);
            ComparisonRecord& rec = record(q + " finest (N = 512), " + at, fin.h, metric, 0.0, u_ref, metric, crit,
                                           ref_ok ? judge(r.info, {}, metric <= tol, crit) : missing_input(missing));
            block_if_interacting(rec, fac);
        };
        finest("psi Linf", fin.psi_linf, psi_tolerance, "|V_bi|");
        finest("interface level", fin.level, level_tolerance, "|V_bi|");
        finest("fields (D), both sides", fin.field, field_tolerance, "itself");
        finest("band bending, left", fin.bend_left, field_tolerance, "itself");
        finest("band bending, right", fin.bend_right, field_tolerance, "itself");
    }
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
