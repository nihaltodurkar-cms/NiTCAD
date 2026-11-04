// Physical-consistency checks C1, C2, C4, C12 (Unit 55, TVR-1, DECISIONS 11.6). See checks.hpp.
#include "NiTCAD/verify/consistency/checks.hpp"

#include <algorithm>
#include <array>
#include <cfloat>
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
    } else if (id == "C4") {
        r.equations = {15};
        r.form = "|sum_k I_k| <= sum_k res_k";
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

std::vector<std::string> validate(const DeviceData& d, const mesh::Mesh& m) {
    std::vector<std::string> p;
    const std::size_t n = m.node_count();
    if (!(std::isfinite(d.temperature_K))) p.push_back("temperature not finite");
    if (d.donors.size() != n) p.push_back("donors not one per node");
    if (d.acceptors.size() != n) p.push_back("acceptors not one per node");
    if (d.relative_permittivity.size() != n) p.push_back("permittivity not one per node");
    if (d.insulator.size() != n) p.push_back("insulator flags not one per node");
    for (const auto& c : d.contacts) {
        if (c.name.empty()) p.push_back("contact name empty");
        for (auto node : c.nodes) {
            if (node < 0 || static_cast<std::size_t>(node) >= n) p.push_back("contact " + c.name + " node out of range");
        }
    }
    return p;
}

}  // namespace NiTCAD::verify::consistency
