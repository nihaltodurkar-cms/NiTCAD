// Physical-consistency checks C1-C5, C12 (Unit 55, TVR-1, DECISIONS 11.6). See checks.hpp.
#include "NiTCAD/verify/consistency/checks.hpp"

#include <algorithm>
#include <array>
#include <cfloat>
#include <cstdint>
#include <charconv>
#include <cmath>
#include <span>
#include <string_view>
#include <utility>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/verify/reference/double_double.hpp"

namespace NiTCAD::verify::consistency {

namespace {

using reference::DD;

std::string number(double x) {
    char buf[32];
    const auto r = std::to_chars(buf, buf + sizeof buf, x);
    return std::string(buf, r.ptr);
}

std::string at(const char* what, std::size_t i) { return std::string(what) + " " + std::to_string(i); }

// A part over many locations: tracks the location whose residual / tolerance is largest.
struct Worst {
    CheckPart part;
    double ratio = -1.0;
    bool failed = false;
    explicit Worst(std::string name) { part.name = std::move(name); }
    void add(double residual, double tolerance, std::string where) {
        // A NaN residual fails: the comparison is written so that NaN does not pass.
        if (!(residual <= tolerance)) failed = true;
        double q = tolerance > 0.0 ? residual / tolerance : (residual == 0.0 ? 0.0 : INFINITY);
        if (std::isnan(q)) q = INFINITY;
        if (q > ratio) {
            ratio = q;
            part.residual = residual;
            part.tolerance = tolerance;
            part.where = std::move(where);
        }
        part.outcome = failed ? Outcome::fail : Outcome::pass;
    }
    CheckPart done(std::string_view empty_reason) {
        if (ratio < 0.0) {  // nothing added
            part.outcome = Outcome::not_evaluated;
            part.reason = std::string(empty_reason);
        }
        return part;
    }
};

CheckPart not_evaluated_part(std::string name, std::string reason) {
    CheckPart p;
    p.name = std::move(name);
    p.outcome = Outcome::not_evaluated;
    p.reason = std::move(reason);
    return p;
}

// The overall outcome from the parts: any fail fails; else any NOT EVALUATED part makes the check NOT EVALUATED
// (a check is not passed on the parts it could evaluate alone); else pass.
CheckResult finish(std::string id, const ReferenceInfo& info, std::span<const Predicate> predicates,
                   std::vector<CheckPart> parts) {
    CheckResult r;
    r.id = std::move(id);
    r.parts = std::move(parts);
    std::string failed, missing;
    for (const auto& p : r.parts) {
        if (p.outcome == Outcome::fail) failed += (failed.empty() ? "" : ", ") + p.name;
        if (p.outcome == Outcome::not_evaluated) {
            missing += (missing.empty() ? "" : "; ") + p.name + ": " + p.reason;
        }
    }
    if (!failed.empty()) {
        r.comparison = judge(info, predicates, false, failed);
    } else if (!missing.empty()) {
        r.comparison = judge(info, predicates, true, "parts evaluated");
        if (r.comparison.outcome == Outcome::pass) r.comparison = missing_input(missing);
    } else {
        r.comparison = judge(info, predicates, true, r.id);
    }
    return r;
}

// Applicability first: if a predicate is false (or one is missing), the parts are not evaluated.
bool applicable(const ReferenceInfo& info, std::span<const Predicate> predicates, CheckResult& out) {
    const Comparison c = judge(info, predicates, true, "applicability");
    if (c.outcome == Outcome::not_applicable || c.outcome == Outcome::not_evaluated) {
        out.comparison = c;
        return false;
    }
    return true;
}

Predicate zero_bias(const results::BiasPoint& point) {
    double m = 0.0;
    for (double v : point.bias_V) m = std::max(m, std::abs(v));
    return {"all biases 0", m, 0.0, m == 0.0};
}

Predicate one_thermal_temperature(const results::BiasPoint& point) {
    double spread = 0.0;
    if (!point.thermal_bias_K.empty()) {
        const auto [lo, hi] = std::minmax_element(point.thermal_bias_K.begin(), point.thermal_bias_K.end());
        spread = *hi - *lo;
    }
    return {"thermal contacts at one temperature", spread, 0.0, spread == 0.0};
}

// The smallest nonzero terminal resolution (gates carry 0), or 0 when there is none.
double smallest_resolution(const results::BiasPoint& point) {
    double m = 0.0;
    for (double r : point.terminal_current_resolution) {
        if (r > 0.0 && (m == 0.0 || r < m)) m = r;
    }
    return m;
}

bool sized(const std::vector<double>& v, std::size_t n) { return v.size() == n; }

// |F_n - F_p| <= 64 eps max(|E_c|, |E_v|, V_T) + V_T u_last at each semiconductor node, u_last the result's last
// Newton update (owner ruling 2026-11-04 after measurement, DECISIONS 11.13: at a drift-diffusion point the
// split carries the Newton stopping error; 0 for the equilibrium solve, which is held to rounding).
CheckPart fermi_part(const results::BandDiagram& b, std::size_t nn, const DeviceData& d, double u_last) {
    if (!sized(b.conduction_eV, nn) || !sized(b.valence_eV, nn) || !sized(b.electron_fermi_eV, nn) ||
        !sized(b.hole_fermi_eV, nn)) {
        return not_evaluated_part("quasi-Fermi levels", "band diagram not stored");
    }
    if (!(std::isfinite(u_last) && u_last >= 0.0)) {
        return not_evaluated_part("quasi-Fermi levels", "Newton update not finite and >= 0");
    }
    const double vt = base::k_B_eV_per_K * d.temperature_K;
    Worst w("quasi-Fermi levels");
    for (std::size_t i = 0; i < nn; ++i) {
        if (d.insulator[i] != 0) continue;
        const double scale = std::max({std::abs(b.conduction_eV[i]), std::abs(b.valence_eV[i]), vt});
        w.add(std::abs(b.electron_fermi_eV[i] - b.hole_fermi_eV[i]),
              rounding_factor * DBL_EPSILON * scale + vt * u_last, at("node", i));
    }
    return w.done("no semiconductor nodes");
}

}  // namespace

CheckResult check_c1(const results::EquilibriumResult& eq, const mesh::Mesh& m, const DeviceData& d) {
    const ReferenceInfo info = check_info("C1");
    // The equilibrium solve is at zero bias by construction and stores no currents (its carriers follow psi
    // with one Fermi level), so only the quasi-Fermi part applies, held to rounding.
    const std::array<Predicate, 2> pred{Predicate{"all biases 0", 0.0, 0.0, true},
                                        Predicate{"thermal contacts at one temperature", 0.0, 0.0, true}};
    std::vector<CheckPart> parts;
    parts.push_back(fermi_part(eq.bands, m.node_count(), d, 0.0));
    return finish("C1", info, pred, std::move(parts));
}

ReferenceInfo check_info(std::string_view id) {
    ReferenceInfo r;
    r.id = std::string(id);
    r.cls = ReferenceClass::D;
    r.independence = Independence::I2;
    r.provenance = Provenance::PD;
    r.citation = "DECISIONS Unit 55 5.1 and 11.6, row " + std::string(id) + " (restated in result quantities)";
    r.version = "1";
    r.derivation = "DECISIONS 11.6 row " + std::string(id);
    r.verification = "defect-injection test, tests/verify/consistency";
    r.second_method = "not required (class D)";
    if (id == "C1") {
        r.equations = {1, 2, 3, 15};
        r.form = "I_k <= res_k; |I_edge| <= min res; |F_n - F_p| <= 64 eps max(|E_c|, |E_v|, V_T)";
        r.assumptions = {"no optical or prescribed generation (none exists in NiTCAD)", "steady state"};
        r.applicability = {"all biases 0", "thermal contacts at one temperature"};
    } else if (id == "C2") {
        r.equations = {2, 7, 8};
        r.form = "q sum_i |R_k,i| V_i <= min nonzero res";
        r.assumptions = {"as C1"};
        r.applicability = {"all biases 0", "thermal contacts at one temperature"};
    } else if (id == "C3") {
        r.equations = {1};
        r.form = "sum_f C_f (psi_j - psi_i) + q V (p - n + N_D - N_A) + stored charges = 0 per cell or group";
        r.assumptions = {"dopants fully ionized (Boltzmann or Fermi-Dirac statistics)",
                         "the interface model only on charged interfaces' insulator-semiconductor edges"};
        r.applicability = {"converged state"};
    } else if (id == "C4") {
        r.equations = {15};
        r.form = "|sum_k I_k| <= sum_k res_k";
    } else if (id == "C5") {
        r.equations = {9, 15};
        r.form = "|sum_c V_c I_c - sum_k H_k| <= E_max sum_c res_c + 64 eps (sum_c |V_c I_c| + sum_k |H_k|)";
        r.assumptions = {"steady state (a BiasPoint)",
                         "global balance only: the per-node heat is not in the results (listed gap)",
                         "no Schottky contact, gate tunnelling or energy balance (refused with the model)"};
        r.applicability = {"electrothermal on", "thermal_contact_heat present"};
    } else if (id == "C6") {
        r.equations = {1, 15};
        r.form = "beta h I_disp = Q_k - a1 Q_k-1 + a2 Q_k-2 per step (BE, variable-step BDF2); gate conduction 0";
        r.assumptions = {"the step's formula from its stored order and steps (BE after a start or breakpoint)",
                         "bookkeeping only: the solver forms I_disp from the same stored contact charges",
                         "the trapezoidal integral against the summed local error is a diagnostic, not judged"};
        r.applicability = {"transient with a gate or electrode"};
    } else if (id == "C12") {
        r.equations = {2, 4, 9, 10};
        r.form = "n, p > 0; T > 0; occupancy in [0, 1]";
    }
    return r;
}

CheckResult check_c1(const results::BiasPoint& point, const mesh::Mesh& m, const DeviceData& d) {
    const ReferenceInfo info = check_info("C1");
    const std::array<Predicate, 2> pred{zero_bias(point), one_thermal_temperature(point)};
    CheckResult out;
    out.id = "C1";
    if (!applicable(info, pred, out)) return out;
    std::vector<CheckPart> parts;
    const std::size_t nc = d.contacts.size(), nn = m.node_count(), ne = m.edges().size();

    if (point.terminal_current.size() != nc || point.terminal_current_resolution.size() != nc) {
        parts.push_back(not_evaluated_part("terminal currents", "terminal currents or resolutions not stored"));
    } else {
        Worst w("terminal currents");
        for (std::size_t k = 0; k < nc; ++k) {
            w.add(std::abs(point.terminal_current[k]), point.terminal_current_resolution[k], "contact " + d.contacts[k].name);
        }
        parts.push_back(w.done("no contacts"));
    }

    const double edge_tol = smallest_resolution(point);
    if (!sized(point.edge_current_n, ne) || !sized(point.edge_current_p, ne)) {
        parts.push_back(not_evaluated_part("edge currents", "edge currents not stored"));
    } else if (edge_tol == 0.0) {
        parts.push_back(not_evaluated_part("edge currents", "no nonzero terminal resolution"));
    } else {
        Worst w("edge currents");
        for (std::size_t e = 0; e < ne; ++e) {
            w.add(std::abs(point.edge_current_n[e]), edge_tol, at("electron edge", e));
            w.add(std::abs(point.edge_current_p[e]), edge_tol, at("hole edge", e));
        }
        parts.push_back(w.done("no edges"));
    }

    if (point.convergence.iterations.empty()) {
        parts.push_back(not_evaluated_part("quasi-Fermi levels", "no Newton record (the stopping term)"));
    } else {
        parts.push_back(fermi_part(point.bands, nn, d, point.convergence.iterations.back().update));
    }
    return finish("C1", info, pred, std::move(parts));
}

CheckResult check_c2(const results::BiasPoint& point, const mesh::Mesh& m, const DeviceData& d) {
    (void)d;
    const ReferenceInfo info = check_info("C2");
    const std::array<Predicate, 2> pred{zero_bias(point), one_thermal_temperature(point)};
    CheckResult out;
    out.id = "C2";
    if (!applicable(info, pred, out)) return out;
    const auto& r = point.rates;
    const std::array<std::pair<const char*, const std::vector<double>*>, 8> mech{{
        {"srh", &r.srh},
        {"auger", &r.auger},
        {"radiative", &r.radiative},
        {"impact", &r.impact},
        {"interface_n", &r.interface_n},
        {"interface_p", &r.interface_p},
        {"band_to_band_n", &r.band_to_band_n},
        {"band_to_band_p", &r.band_to_band_p},
    }};
    const double tol = smallest_resolution(point);
    const auto vol = m.volumes();
    std::vector<CheckPart> parts;
    bool any = false;
    for (const auto& [name, v] : mech) {
        if (v->empty()) continue;
        any = true;
        if (v->size() != m.node_count()) {
            parts.push_back(not_evaluated_part(name, "rate vector not one per node"));
            continue;
        }
        if (tol == 0.0) {
            parts.push_back(not_evaluated_part(name, "no nonzero terminal resolution"));
            continue;
        }
        DD sum(0.0);
        for (std::size_t i = 0; i < v->size(); ++i) sum += reference::two_prod(std::abs((*v)[i]), vol[i]);
        Worst w(name);
        w.add(base::q_C * sum.to_double(), tol, "integral over the device");
        parts.push_back(w.done(""));
    }
    if (!any) parts.push_back(not_evaluated_part("rates", "no rates stored"));
    return finish("C2", info, pred, std::move(parts));
}

CheckResult check_c4(const results::BiasPoint& point, const mesh::Mesh& m, const DeviceData& d) {
    (void)m;
    const ReferenceInfo info = check_info("C4");
    CheckResult out;
    out.id = "C4";
    if (!applicable(info, {}, out)) return out;
    const std::size_t nc = d.contacts.size();
    std::vector<CheckPart> parts;
    auto kcl = [&parts, nc](const char* name, const std::vector<double>& current, const std::vector<double>& res) {
        if (current.size() != nc || res.size() != nc) {
            parts.push_back(not_evaluated_part(name, "currents or resolutions not stored"));
            return;
        }
        DD sum(0.0), tol(0.0);
        for (std::size_t k = 0; k < nc; ++k) {
            sum += DD(current[k]);
            tol += DD(res[k]);
        }
        Worst w(name);
        w.add(std::abs(sum.to_double()), tol.to_double(), "sum over contacts");
        parts.push_back(w.done(""));
    };
    kcl("terminal currents", point.terminal_current, point.terminal_current_resolution);
    if (!point.weighted_current.empty()) {
        kcl("weighted currents", point.weighted_current, point.weighted_current_resolution);
    }
    return finish("C4", info, {}, std::move(parts));
}

// The steady energy balance (owner rulings 2026-11-04, DECISIONS 11.13): the power the contacts deliver against
// the heat leaving through the thermal contacts, to Unit 23's U8 scale restated in result quantities.
CheckResult check_c5(const results::BiasPoint& point, const mesh::Mesh& m, const DeviceData& d) {
    const ReferenceInfo info = check_info("C5");
    CheckResult out;
    out.id = "C5";
    const std::array<Predicate, 2> predicates{
        Predicate{"electrothermal on", d.electrothermal ? 1.0 : 0.0, 1.0, d.electrothermal},
        Predicate{"thermal_contact_heat present", static_cast<double>(point.thermal_contact_heat.size()), 1.0,
                  !point.thermal_contact_heat.empty()}};
    if (!applicable(info, predicates, out)) return out;
    const char* name = "global energy balance";
    const std::size_t nc = d.contacts.size();
    const std::size_t nn = m.node_count();
    std::string missing;
    if (!sized(point.bias_V, nc) || !sized(point.terminal_current, nc) ||
        !sized(point.terminal_current_resolution, nc)) {
        missing = "biases, currents or resolutions not stored";
    } else if (point.thermal_contact_heat.size() != d.thermal_contacts) {
        missing = "thermal_contact_heat not one per thermal contact";
    } else if (!sized(point.bands.conduction_eV, nn) || !sized(point.bands.valence_eV, nn)) {
        missing = "bands not stored";
    } else if (!sized(point.fields.temperature_K, nn)) {
        missing = "lattice temperature not stored";
    } else if (!sized(d.thermopower_n, nn) || !sized(d.thermopower_p, nn)) {
        missing = "thermopower exponents not in the device data";
    }
    if (!missing.empty()) return finish("C5", info, predicates, {not_evaluated_part(name, missing)});

    DD power(0.0), heat(0.0), magnitude(0.0), resolution(0.0);
    double e_max = 0.0;
    for (std::size_t c = 0; c < nc; ++c) {
        const DD vi = DD(point.bias_V[c]) * DD(point.terminal_current[c]);
        power += vi;
        magnitude += DD(std::abs(vi.to_double()));
        resolution += DD(point.terminal_current_resolution[c]);
        e_max = std::max(e_max, std::abs(point.bias_V[c]));
    }
    for (double h : point.thermal_contact_heat) {
        heat += DD(h);
        magnitude += DD(std::abs(h));
    }
    double t_max = 0.0, r_max = 0.0;
    for (std::size_t i = 0; i < nn; ++i) {
        t_max = std::max(t_max, point.fields.temperature_K[i]);
        if (d.insulator[i] != 0) continue;
        e_max = std::max({e_max, std::abs(point.bands.conduction_eV[i]), std::abs(point.bands.valence_eV[i])});
        r_max = std::max({r_max, std::abs(d.thermopower_n[i] + 2.5), std::abs(d.thermopower_p[i] + 2.5)});
    }
    e_max += r_max * base::k_B_eV_per_K * t_max;
    const double residual = std::abs((power - heat).to_double());
    const double tolerance =
        e_max * resolution.to_double() + rounding_factor * DBL_EPSILON * magnitude.to_double();
    Worst w(name);
    w.add(residual, tolerance, "sum over contacts");
    return finish("C5", info, predicates, {w.done("")});
}

// Transient charge bookkeeping (owner rulings 2026-11-04, DECISIONS 11.13). The BDF coefficients are derived here
// from the formulas' definitions, not read from the solver: backward Euler y_k - y_k-1 = h f_k; variable-step BDF2
// from the quadratic through (t_k-2, t_k-1, t_k) differentiated at t_k, scaled so the y_k coefficient is 1.
CheckResult check_c6(const results::Transient& run, const mesh::Mesh& m, const DeviceData& d) {
    (void)m;
    const ReferenceInfo info = check_info("C6");
    CheckResult out;
    out.id = "C6";
    double gates = 0.0;
    for (const auto& c : d.contacts) {
        if (c.kind == ContactKind::gate || c.kind == ContactKind::electrode) gates += 1.0;
    }
    const std::array<Predicate, 1> predicates{
        Predicate{"transient with a gate or electrode", gates, 1.0, gates > 0.0}};
    if (!applicable(info, predicates, out)) return out;

    const char* identity = "displacement identity";
    const std::size_t nc = d.contacts.size();
    const auto& pts = run.points;
    if (pts.size() < 2) return finish("C6", info, predicates, {not_evaluated_part(identity, "fewer than two time points")});
    for (std::size_t k = 0; k < pts.size(); ++k) {
        const auto& p = pts[k];
        if (!sized(p.contact_charge, nc) || !sized(p.displacement_current, nc) || !sized(p.conduction_current, nc) ||
            !sized(p.terminal_current, nc)) {
            return finish("C6", info, predicates,
                          {not_evaluated_part(identity, "charges or currents not one per contact at " + at("point", k))});
        }
    }

    struct Coefficients {
        double beta, a1, a2;
    };
    Worst w(identity);
    std::vector<DD> lhs_sum(nc, DD(0.0)), rhs_sum(nc, DD(0.0));
    for (std::size_t k = 1; k < pts.size(); ++k) {
        const auto& p = pts[k];
        const double h = p.step_s;
        Coefficients f{};
        if (p.order == 1) {
            f = {1.0, 1.0, 0.0};
        } else if (p.order == 2 && k >= 2 && pts[k - 1].step_s > 0.0) {
            const double r = h / pts[k - 1].step_s;  // the step ratio
            const double den = 1.0 + 2.0 * r;
            f = {(1.0 + r) / den, (1.0 + r) * (1.0 + r) / den, r * r / den};
        } else {
            return finish("C6", info, predicates,
                          {not_evaluated_part(identity, "order " + std::to_string(p.order) +
                                                            " without its earlier steps at " + at("point", k))});
        }
        if (!(h > 0.0)) {
            return finish("C6", info, predicates, {not_evaluated_part(identity, "step not positive at " + at("point", k))});
        }
        for (std::size_t c = 0; c < nc; ++c) {
            const DD lhs = DD(f.beta) * DD(h) * DD(p.displacement_current[c]);
            DD rhs = DD(p.contact_charge[c]) - DD(f.a1) * DD(pts[k - 1].contact_charge[c]);
            double operands = std::abs(p.contact_charge[c]) + std::abs(f.a1 * pts[k - 1].contact_charge[c]) +
                              std::abs(lhs.to_double());
            if (p.order == 2) {
                rhs += DD(f.a2) * DD(pts[k - 2].contact_charge[c]);
                operands += std::abs(f.a2 * pts[k - 2].contact_charge[c]);
            }
            lhs_sum[c] += lhs;
            rhs_sum[c] += rhs;
            w.add(std::abs((lhs - rhs).to_double()), rounding_factor * DBL_EPSILON * operands,
                  at("point", k) + ", contact " + d.contacts[c].name);
        }
    }
    std::vector<CheckPart> parts{w.done("")};

    // Gates and electrodes without gate tunnelling: no conduction current; terminal = displacement, exactly.
    Worst g("gate and electrode currents");
    bool any_gate = false;
    for (std::size_t c = 0; c < nc; ++c) {
        const auto& k = d.contacts[c];
        if ((k.kind != ContactKind::gate && k.kind != ContactKind::electrode) || k.tunnelling) continue;
        any_gate = true;
        for (std::size_t i = 0; i < pts.size(); ++i) {
            const auto& p = pts[i];
            const double r = std::abs(p.conduction_current[c]) + std::abs(p.terminal_current[c] - p.displacement_current[c]);
            g.add(r, 0.0, at("point", i) + ", contact " + k.name);
        }
    }
    if (any_gate) parts.push_back(g.done(""));

    // Diagnostics (not judged): the cumulative identity and the trapezoidal integral per contact.
    double rtol = std::nan("");
    for (const auto& [name, value] : run.run.settings) {
        if (name == "transient.rtol") rtol = value;
    }
    for (std::size_t c = 0; c < nc; ++c) {
        const std::string& name = d.contacts[c].name;
        out.diagnostics.push_back({"cumulative identity", (lhs_sum[c] - rhs_sum[c]).to_double(),
                                   std::abs(rhs_sum[c].to_double()), "contact " + name,
                                   "sum_k beta h I_disp - sum_k (Q_k - a1 Q_k-1 + a2 Q_k-2)"});
        DD trapezoid(0.0);
        double budget = 0.0;
        std::size_t unestimated = 0;
        for (std::size_t k = 1; k < pts.size(); ++k) {
            trapezoid += DD(0.5) * DD(pts[k].step_s) *
                         (DD(pts[k].displacement_current[c]) + DD(pts[k - 1].displacement_current[c]));
            if (std::isnan(pts[k].error_ratio)) {
                ++unestimated;
            } else {
                budget += pts[k].error_ratio * rtol * std::abs(pts[k].contact_charge[c] - pts[k - 1].contact_charge[c]);
            }
        }
        const DD dq = DD(pts.back().contact_charge[c]) - DD(pts.front().contact_charge[c]);
        out.diagnostics.push_back(
            {"trapezoidal integral", (trapezoid - dq).to_double(), budget, "contact " + name,
             "integral of I_disp dt - (Q_N - Q_0); scale sum error_ratio rtol |dQ_k| over " +
                 std::to_string(pts.size() - 1 - unestimated) + " estimated steps (" + std::to_string(unestimated) +
                 " without an estimate)" + (std::isnan(rtol) ? "; rtol not in the run record" : "")});
    }
    CheckResult r = finish("C6", info, predicates, std::move(parts));
    r.diagnostics = std::move(out.diagnostics);
    return r;
}

CheckResult check_c12(const results::NodeFields& f, const results::BiasPoint* point, const mesh::Mesh& m,
                      const DeviceData& d) {
    const ReferenceInfo info = check_info("C12");
    CheckResult out;
    out.id = "C12";
    if (!applicable(info, {}, out)) return out;
    const std::size_t nn = m.node_count();
    std::vector<CheckPart> parts;
    if (f.n_cm3.size() != nn || f.p_cm3.size() != nn) {
        parts.push_back(not_evaluated_part("densities", "densities not stored"));
    } else {
        // Exact: a density is positive or it is not. The residual is 1 at a violation, 0 otherwise.
        Worst w("densities");
        for (std::size_t i = 0; i < nn; ++i) {
            if (d.insulator[i] != 0) continue;
            w.add(f.n_cm3[i] > 0.0 ? 0.0 : 1.0, 0.0, at("electron node", i) + " n = " + number(f.n_cm3[i]));
            w.add(f.p_cm3[i] > 0.0 ? 0.0 : 1.0, 0.0, at("hole node", i) + " p = " + number(f.p_cm3[i]));
        }
        parts.push_back(w.done("no semiconductor nodes"));
    }
    {
        Worst w("temperatures");
        auto scan = [&w](const char* what, const std::vector<double>& v) {
            for (std::size_t i = 0; i < v.size(); ++i) {
                w.add(v[i] > 0.0 ? 0.0 : 1.0, 0.0, std::string(what) + " " + std::to_string(i) + " T = " + number(v[i]));
            }
        };
        scan("lattice node", f.temperature_K);
        scan("electron node", f.electron_temperature_K);
        scan("hole node", f.hole_temperature_K);
        if (point != nullptr) scan("thermal contact", point->thermal_bias_K);
        w.add(d.temperature_K > 0.0 ? 0.0 : 1.0, 0.0, "device temperature T = " + number(d.temperature_K));
        parts.push_back(w.done(""));
    }
    if (d.interface_traps) {
        parts.push_back(not_evaluated_part("trap occupancies", "occupancies not stored in the results"));
    }
    return finish("C12", info, {}, std::move(parts));
}

namespace {

// Disjoint sets over the nodes (path halving, union by index): the groups of 11.13's C3 ruling.
struct Groups {
    std::vector<std::size_t> parent;
    explicit Groups(std::size_t n) : parent(n) {
        for (std::size_t i = 0; i < n; ++i) parent[i] = i;
    }
    std::size_t find(std::size_t i) {
        while (parent[i] != i) i = parent[i] = parent[parent[i]];
        return i;
    }
    void unite(std::size_t a, std::size_t b) {
        a = find(a);
        b = find(b);
        if (a != b) parent[std::max(a, b)] = std::min(a, b);
    }
};

// The declared interface of an edge whose flux the interface model splits (an insulator-semiconductor edge
// between the regions of a charged interface), or -1: the same edges the assembler marks (scaled_device.cpp).
std::int32_t charged_interface(const DeviceData& d, std::size_t a, std::size_t b) {
    if (d.insulator[a] == d.insulator[b]) return -1;
    const std::int32_t ra = d.region[a], rb = d.region[b];
    if (ra == rb) return -1;
    for (std::size_t k = 0; k < d.interfaces.size(); ++k) {
        const InterfaceData& f = d.interfaces[k];
        if ((f.region_a == ra && f.region_b == rb) || (f.region_a == rb && f.region_b == ra)) {
            return f.charged ? static_cast<std::int32_t>(k) : -1;
        }
    }
    return -1;
}

// One cell or group's Gauss balance: the residual sum (flux in across its boundary + its space charge + the
// stored charges inside it), its rounding scale and its Newton term.
struct Balance {
    DD residual{0.0};
    double scale = 0.0;
    double newton = 0.0;
    bool grouped = false;   // a gate, an electrode or a charged interface is in it
    std::int32_t contact = -1;  // an ohmic or Schottky contact it belongs to (skipped), else -1
    std::size_t first = 0;      // its lowest node, for the location
    std::string members;        // what it holds, for the location
};

struct Inputs {
    const results::NodeFields& fields;
    const std::vector<double>& gate_charge;
    const std::vector<double>& interface_trap_charge;
    const results::ConvergenceRecord& convergence;
};

CheckResult c3(const Inputs& in, const mesh::Mesh& m, const DeviceData& d) {
    const ReferenceInfo info = check_info("C3");
    const std::array<Predicate, 1> pred{
        Predicate{"converged state", in.convergence.converged ? 1.0 : 0.0, 1.0, in.convergence.converged}};
    CheckResult out;
    out.id = "C3";
    if (!applicable(info, pred, out)) return out;
    const std::size_t nn = m.node_count(), nc = d.contacts.size();
    const auto& f = in.fields;
    std::vector<CheckPart> parts;
    const char* const cells = "cells";
    if (d.incomplete_ionization) {
        parts.push_back(not_evaluated_part(cells, "incomplete ionization: the ionized densities are not stored"));
        return finish("C3", info, pred, std::move(parts));
    }
    if (!sized(f.potential_V, nn) || !sized(f.n_cm3, nn) || !sized(f.p_cm3, nn)) {
        parts.push_back(not_evaluated_part(cells, "potential or densities not stored"));
        return finish("C3", info, pred, std::move(parts));
    }
    if (in.convergence.iterations.empty()) {
        parts.push_back(not_evaluated_part(cells, "no Newton record (the stopping term)"));
        return finish("C3", info, pred, std::move(parts));
    }
    const double u = in.convergence.iterations.back().update;
    if (!(std::isfinite(u) && u >= 0.0)) {
        parts.push_back(not_evaluated_part(cells, "Newton update not finite and >= 0"));
        return finish("C3", info, pred, std::move(parts));
    }
    const double vt = base::k_B_eV_per_K * d.temperature_K;  // V_T in V
    const auto vol = m.volumes();
    const auto edges = m.edges();
    const auto& psi = f.potential_V;

    // The groups: each contact's nodes, each charged edge's ends, each charged interface's edges together.
    Groups g(nn);
    std::vector<std::int32_t> node_contact(nn, -1);
    for (std::size_t k = 0; k < nc; ++k) {
        const auto& nodes = d.contacts[k].nodes;
        for (std::size_t j = 0; j < nodes.size(); ++j) {
            const auto i = static_cast<std::size_t>(nodes[j]);
            node_contact[i] = static_cast<std::int32_t>(k);
            g.unite(i, static_cast<std::size_t>(nodes[0]));
        }
    }
    std::vector<std::int32_t> edge_interface(edges.size(), -1);
    std::vector<std::int64_t> interface_node(d.interfaces.size(), -1);
    std::vector<double> interface_area(d.interfaces.size(), 0.0);
    for (std::size_t e = 0; e < edges.size(); ++e) {
        const auto a = static_cast<std::size_t>(edges[e].first), b = static_cast<std::size_t>(edges[e].second);
        const std::int32_t k = charged_interface(d, a, b);
        if (k < 0) continue;
        edge_interface[e] = k;
        interface_area[static_cast<std::size_t>(k)] += edges[e].coupling_area;
        g.unite(a, b);
        auto& first = interface_node[static_cast<std::size_t>(k)];
        if (first < 0) first = static_cast<std::int64_t>(a);
        g.unite(a, static_cast<std::size_t>(first));
    }

    // Balances, one per root.
    std::vector<std::int64_t> slot(nn, -1);
    std::vector<Balance> bal;
    std::vector<std::size_t> owner(nn);
    for (std::size_t i = 0; i < nn; ++i) {
        const std::size_t r = g.find(i);
        if (slot[r] < 0) {
            slot[r] = static_cast<std::int64_t>(bal.size());
            bal.emplace_back();
            bal.back().first = i;
        }
        owner[i] = static_cast<std::size_t>(slot[r]);
    }
    const double q = base::q_C;
    for (std::size_t i = 0; i < nn; ++i) {
        Balance& b = bal[owner[i]];
        const std::int32_t k = node_contact[i];
        if (k >= 0) {
            const ContactKind kind = d.contacts[static_cast<std::size_t>(k)].kind;
            if (kind == ContactKind::ohmic || kind == ContactKind::schottky) b.contact = k;
            else b.grouped = true;
        }
        if (d.insulator[i] != 0) continue;
        // q V (p - n + N_D - N_A), summed exactly in double-double.
        const double V = vol[i];
        DD rho = DD(f.p_cm3[i]) - DD(f.n_cm3[i]) + DD(d.donors[i]) - DD(d.acceptors[i]);
        b.residual += rho * DD(q * V);
        b.scale += q * V * (f.n_cm3[i] + f.p_cm3[i] + d.donors[i] + d.acceptors[i]);
        b.newton += q * V * (f.n_cm3[i] + f.p_cm3[i]);
    }
    for (std::size_t e = 0; e < edges.size(); ++e) {
        const auto a = static_cast<std::size_t>(edges[e].first), b = static_cast<std::size_t>(edges[e].second);
        const double ea = d.relative_permittivity[a], eb = d.relative_permittivity[b];
        const double eps = ea == eb ? ea : 2.0 * ea * eb / (ea + eb);  // the discretization's harmonic mean
        const double C = base::eps0_F_per_cm * eps * edges[e].coupling_area / edges[e].length;
        double mag = C * (std::abs(psi[a]) + std::abs(psi[b]));
        if (edge_interface[e] >= 0) {
            // The two half-edges (conductance 2 eps_end c) at an interface potential between the ends.
            mag += 2.0 * base::eps0_F_per_cm * std::max(ea, eb) * edges[e].coupling_area / edges[e].length *
                   (std::abs(psi[a]) + std::abs(psi[b]));
        }
        Balance& ba = bal[owner[a]];
        Balance& bb = bal[owner[b]];
        ba.scale += mag;
        bb.scale += mag;
        if (owner[a] == owner[b]) continue;  // inside one group: cancels in its sum
        const DD flux = reference::two_prod(C, psi[b]) - reference::two_prod(C, psi[a]);  // into a
        ba.residual += flux;
        bb.residual -= flux;
        ba.newton += 2.0 * C * vt;
        bb.newton += 2.0 * C * vt;
    }
    // The stored charges.
    bool gate_missing = false, trap_missing = false;
    for (std::size_t k = 0; k < nc; ++k) {
        const ContactData& c = d.contacts[k];
        if ((c.kind != ContactKind::gate && c.kind != ContactKind::electrode) || c.nodes.empty()) continue;
        Balance& b = bal[owner[static_cast<std::size_t>(c.nodes[0])]];
        b.members += (b.members.empty() ? "" : ", ") + std::string(c.kind == ContactKind::gate ? "gate " : "electrode ") + c.name;
        if (in.gate_charge.size() != nc) {
            gate_missing = true;
            continue;
        }
        const double Q = in.gate_charge[k];
        b.residual += DD(Q);
        b.scale += std::abs(Q);
        if (c.kind == ContactKind::gate) {
            DD area(0.0);
            for (std::size_t j = 0; j < c.nodes.size() && j < c.gate_areas.size(); ++j) {
                const double A = c.gate_areas[j];
                area += DD(A);
                const double cox = c.oxide_capacitance_F_cm2 * A;
                b.scale += 2.0 * cox * std::abs(psi[static_cast<std::size_t>(c.nodes[j])]);
                b.newton += cox * vt;
            }
            const double S = q * c.gate_fixed_charge_cm2 * area.to_double();
            b.residual += DD(S);
            b.scale += std::abs(S);
        }
    }
    for (std::size_t k = 0; k < d.interfaces.size(); ++k) {
        if (interface_node[k] < 0) continue;
        Balance& b = bal[owner[static_cast<std::size_t>(interface_node[k])]];
        b.grouped = true;
        b.members += (b.members.empty() ? "" : ", ") + std::string("interface ") + std::to_string(k);
        const InterfaceData& fi = d.interfaces[k];
        const double S = q * fi.fixed_charge_cm2 * interface_area[k];
        b.residual += DD(S);
        b.scale += std::abs(S);
        if (in.interface_trap_charge.size() != d.interfaces.size()) {
            if (fi.trap_density_cm2 > 0.0) trap_missing = true;
            continue;
        }
        const double Qt = in.interface_trap_charge[k];
        b.residual += DD(Qt);
        // Occupancy moves at most 2 u per trap for a scaled update u (the interface level moves at most u, the
        // densities by a factor e^(+-u)): the Newton term, and the traps' own rounding.
        const double Nt = q * fi.trap_density_cm2 * interface_area[k];
        b.scale += std::abs(Qt) + Nt;
        b.newton += 2.0 * Nt;
    }

    // Parts: the cells alone, the groups, and global neutrality (the sum of every evaluated residual, with each
    // skipped contact's charge taken as its cells' Gauss flux).
    Worst wc(cells), wg("groups");
    DD total(0.0);
    double total_tol = 0.0;
    std::size_t groups = 0;
    std::string contact_charges;
    for (const Balance& b : bal) {
        const double tol = rounding_factor * DBL_EPSILON * b.scale + b.newton * u;
        if (b.contact >= 0) {
            // The contact's charge balances its cells: minus their residual (reported, not judged).
            contact_charges += (contact_charges.empty() ? "" : ", ") + d.contacts[static_cast<std::size_t>(b.contact)].name +
                               " " + number(-b.residual.to_double());
            continue;
        }
        total += b.residual;
        total_tol += tol;
        const double r = std::abs(b.residual.to_double());
        if (b.grouped) {
            ++groups;
            wg.add(r, tol, "group at node " + std::to_string(b.first) + " (" + b.members + ")");
        } else {
            wc.add(r, tol, at("node", b.first));
        }
    }
    parts.push_back(wc.done("no cell restated alone"));
    if (gate_missing) parts.push_back(not_evaluated_part("groups", "gate charges not stored"));
    else if (trap_missing) parts.push_back(not_evaluated_part("groups", "interface trap charges not stored"));
    else if (groups > 0) parts.push_back(wg.done(""));
    if (gate_missing || trap_missing) {
        // A stored charge is missing from the sum as from its group.
        parts.push_back(not_evaluated_part("global neutrality (sum of the other parts' residuals)",
                                           gate_missing ? "gate charges not stored" : "interface trap charges not stored"));
    } else {
        Worst w("global neutrality (sum of the other parts' residuals)");
        w.add(std::abs(total.to_double()), total_tol,
              "device; contact charges from their cells' flux: " + (contact_charges.empty() ? "none" : contact_charges));
        parts.push_back(w.done(""));
    }
    return finish("C3", info, pred, std::move(parts));
}

}  // namespace

CheckResult check_c3(const results::BiasPoint& point, const mesh::Mesh& m, const DeviceData& d) {
    return c3({point.fields, point.gate_charge, point.interface_trap_charge, point.convergence}, m, d);
}

CheckResult check_c3(const results::EquilibriumResult& eq, const mesh::Mesh& m, const DeviceData& d) {
    return c3({eq.fields, eq.gate_charge, eq.interface_trap_charge, eq.convergence}, m, d);
}

std::vector<std::string> validate(const DeviceData& d, const mesh::Mesh& m) {
    std::vector<std::string> p;
    const std::size_t n = m.node_count();
    if (!(std::isfinite(d.temperature_K))) p.push_back("temperature not finite");
    if (d.donors.size() != n) p.push_back("donors not one per node");
    if (d.acceptors.size() != n) p.push_back("acceptors not one per node");
    if (d.relative_permittivity.size() != n) p.push_back("permittivity not one per node");
    if (d.insulator.size() != n) p.push_back("insulator flags not one per node");
    if (d.region.size() != n) p.push_back("regions not one per node");
    if (!d.thermopower_n.empty() && d.thermopower_n.size() != n) p.push_back("thermopower (n) not one per node");
    if (!d.thermopower_p.empty() && d.thermopower_p.size() != n) p.push_back("thermopower (p) not one per node");
    for (const auto& c : d.contacts) {
        if (c.name.empty()) p.push_back("contact name empty");
        for (auto node : c.nodes) {
            if (node < 0 || static_cast<std::size_t>(node) >= n) p.push_back("contact " + c.name + " node out of range");
        }
        if (c.kind == ContactKind::gate) {
            if (c.gate_areas.size() != c.nodes.size()) p.push_back("gate " + c.name + " areas not one per node");
            if (!(std::isfinite(c.oxide_capacitance_F_cm2) && c.oxide_capacitance_F_cm2 > 0.0)) {
                p.push_back("gate " + c.name + " oxide capacitance not finite and positive");
            }
            if (!std::isfinite(c.gate_fixed_charge_cm2)) p.push_back("gate " + c.name + " fixed charge not finite");
        }
    }
    for (const auto& f : d.interfaces) {
        if (!std::isfinite(f.fixed_charge_cm2) || !(std::isfinite(f.trap_density_cm2) && f.trap_density_cm2 >= 0.0)) {
            p.push_back("interface charge or trap density not finite");
        }
    }
    return p;
}

}  // namespace NiTCAD::verify::consistency
