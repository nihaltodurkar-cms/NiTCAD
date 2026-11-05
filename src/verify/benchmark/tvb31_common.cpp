// Shared machinery of the TVB-31 entries (Unit 55, TVR-3; tvb31_common.hpp).
#include "tvb31_common.hpp"

#include <bit>
#include <utility>

#include "NiTCAD/physics/energy_balance.hpp"
#include "NiTCAD/physics/mobility.hpp"
#include "NiTCAD/physics/surface_mobility.hpp"
#include <charconv>
#include <cstdint>
#include <stdexcept>

namespace NiTCAD::verify::benchmark::tvb31 {

namespace {

bool starts_with(std::string_view s, std::string_view p) { return s.substr(0, p.size()) == p; }

// "name[k]" -> k in [0, 3], or -1.
int index_of(std::string_view field, std::string_view name) {
    if (!starts_with(field, name) || field.size() != name.size() + 3 || field[name.size()] != '[' ||
        field.back() != ']') {
        return -1;
    }
    const char d = field[name.size() + 1];
    return d >= '0' && d <= '3' ? d - '0' : -1;
}

Provenance provenance_of(reg::Grade g) {
    switch (g) {
        case reg::Grade::PA: return Provenance::PA;
        case reg::Grade::PB: return Provenance::PB;
        case reg::Grade::PC: return Provenance::PC;
        case reg::Grade::PD: return Provenance::PD;
        case reg::Grade::L: return Provenance::L;
        case reg::Grade::C: return Provenance::PD;  // pinned through I0 below
    }
    return Provenance::PC;
}

std::string_view grade_name(reg::Grade g) {
    switch (g) {
        case reg::Grade::PA: return "PA";
        case reg::Grade::PB: return "PB";
        case reg::Grade::PC: return "PC";
        case reg::Grade::PD: return "PD";
        case reg::Grade::L: return "L";
        case reg::Grade::C: return "C";
    }
    return "?";
}

// Shipped values outside the semiconductor sets (the Masetti and Lombardi sets, the energy-balance constants).
std::optional<double> other_shipped_value(std::string_view binding) {
    const auto masetti = [](const physics::MasettiParameters& m, std::string_view g) -> std::optional<double> {
        if (g == "mu_min1") return m.mu_min1;
        if (g == "mu_min2") return m.mu_min2;
        if (g == "mu_max") return m.mu_max;
        if (g == "mu_1") return m.mu_1;
        if (g == "P_c") return m.P_c;
        if (g == "C_r") return m.C_r;
        if (g == "C_s") return m.C_s;
        if (g == "alpha") return m.alpha;
        if (g == "beta") return m.beta;
        if (g == "T_exponent") return m.T_exponent;
        return std::nullopt;
    };
    const auto lombardi = [](const physics::LombardiParameters& l, std::string_view g) -> std::optional<double> {
        if (g == "B") return l.B;
        if (g == "C") return l.C;
        if (g == "lambda") return l.lambda;
        if (g == "delta") return l.delta;
        return std::nullopt;
    };
    const auto tail = [binding](std::string_view prefix) -> std::optional<std::string_view> {
        if (!starts_with(binding, prefix)) return std::nullopt;
        return binding.substr(prefix.size());
    };
    if (const auto g = tail("masetti_boron_silicon.")) return masetti(physics::masetti_boron_silicon, *g);
    if (const auto g = tail("masetti_phosphorus_silicon.")) return masetti(physics::masetti_phosphorus_silicon, *g);
    if (const auto g = tail("lombardi_silicon_electrons.")) return lombardi(physics::lombardi_silicon_electrons, *g);
    if (const auto g = tail("lombardi_silicon_holes.")) return lombardi(physics::lombardi_silicon_holes, *g);
    if (binding == "energy_relaxation.electron") return physics::silicon_energy_relaxation.tau_n_s;
    if (binding == "energy_relaxation.hole") return physics::silicon_energy_relaxation.tau_p_s;
    if (binding == "energy_balance_regularization") return physics::energy_balance_regularization;
    return std::nullopt;
}

}  // namespace

double* bound_field(physics::SemiconductorParameters& p, std::string_view binding, std::string_view set) {
    const std::string prefix = std::string(set) + ".";
    if (!starts_with(binding, prefix)) return nullptr;
    const std::string_view f = binding.substr(prefix.size());
    if (f == "eps_r") return &p.eps_r;
    if (f == "Eg0_eV") return &p.Eg0_eV;
    if (f == "varshni_alpha_eV_per_K") return &p.varshni_alpha_eV_per_K;
    if (f == "varshni_beta_K") return &p.varshni_beta_K;
    if (f == "Nc300") return &p.Nc300;
    if (f == "Nv300") return &p.Nv300;
    if (f == "electron_affinity_eV") return &p.electron_affinity_eV;
    if (f == "passler_gap.alpha_eV_per_K") return &p.passler_gap.alpha_eV_per_K;
    if (f == "passler_gap.Theta_K") return &p.passler_gap.Theta_K;
    if (f == "passler_gap.Delta") return &p.passler_gap.Delta;
    if (const int k = index_of(f, "dos_mass.conduction"); k >= 0) return &p.dos_mass.conduction[k];
    if (const int k = index_of(f, "dos_mass.valence"); k >= 0) return &p.dos_mass.valence[k];
    if (f == "bandgap_narrowing.E0_eV") return &p.bandgap_narrowing.E0_eV;
    if (f == "bandgap_narrowing.N0") return &p.bandgap_narrowing.N0;
    for (const auto& [name, c] : {std::pair{"electron_mobility.", &p.electron_mobility}, std::pair{"hole_mobility.", &p.hole_mobility}}) {
        if (!starts_with(f, name)) continue;
        const std::string_view g = f.substr(std::string_view(name).size());
        if (g == "mu_min") return &c->mu_min;
        if (g == "mu_max") return &c->mu_max;
        if (g == "N_ref") return &c->N_ref;
        if (g == "alpha") return &c->alpha;
        if (g == "T_exponent") return &c->T_exponent;
    }
    for (const auto& [name, c] :
         {std::pair{"electron_saturation.", &p.electron_saturation}, std::pair{"hole_saturation.", &p.hole_saturation}}) {
        if (!starts_with(f, name)) continue;
        const std::string_view g = f.substr(std::string_view(name).size());
        if (g == "v_sat_cm_s") return &c->v_sat_cm_s;
        if (g == "beta") return &c->beta;
    }
    return nullptr;
}

std::optional<double> shipped_value(std::string_view binding) {
    physics::SemiconductorParameters legacy = physics::legacy_silicon_parameters;
    physics::SemiconductorParameters si = physics::silicon_parameters;
    if (const double* v = bound_field(legacy, binding, "legacy_silicon")) return *v;
    if (const double* v = bound_field(si, binding, "silicon")) return *v;
    if (const auto v = other_shipped_value(binding)) return v;
    return std::nullopt;
}

physics::SemiconductorParameters register_parameters(physics::SemiconductorParameters base, std::string_view set,
                                                     std::initializer_list<std::string_view> entries) {
    for (const reg::Row& row : reg::rows) {
        if (row.kind != reg::Kind::printed || row.source_value.empty()) continue;
        if (std::find(entries.begin(), entries.end(), row.entry) == entries.end()) continue;
        // A silicon binding also sets the default set's copy of a legacy field the default set inherits.
        double* f = bound_field(base, row.binding, set);
        if (f == nullptr && set == "silicon") f = bound_field(base, row.binding, "legacy_silicon");
        if (f != nullptr) *f = ref::decimal::to_double(row.source_value);
    }
    return base;
}

DD register_dd(std::string_view binding) {
    for (const reg::Row& row : reg::rows) {
        if (row.binding == binding && !row.source_value.empty()) return ref::decimal::to_dd(row.source_value);
    }
    throw std::invalid_argument("register: no printed value for " + std::string(binding));
}

namespace {

double half_unit(std::string_view printed) { return 0.5 * std::pow(10.0, ref::decimal::parse(printed).exponent); }

// A law_300K row: the shipped value against coefficient 300^exponent, to the change of the law's value for half a
// unit in the last printed digit of the coefficient and of the exponent (summed). A mismatch is a register finding.
ComparisonRecord law_check(const ReferenceInfo& info, const reg::Row& row, double shipped, const std::string& quantity) {
    const auto law = std::find_if(reg::laws.begin(), reg::laws.end(),
                                  [&](const reg::Law& l) { return l.binding == row.binding; });
    ComparisonRecord c{quantity, std::string(row.unit), 300.0, shipped, 0.0, 0.0, 0.0, "", {}};
    if (law == reg::laws.end()) {
        c.criterion = "law of the law_300K row known";
        c.comparison = missing_input("the law of " + std::string(row.binding));
        return c;
    }
    const DD A = ref::decimal::to_dd(law->coefficient);
    const DD e = ref::decimal::to_dd(law->exponent);
    const DD v = A * pow(DD(300.0), e);
    const double value = v.to_double();
    const double tol = std::abs(value / A.to_double()) * half_unit(law->coefficient) +
                       std::abs(value * std::log(300.0)) * half_unit(law->exponent);
    c.q_ref = value;
    c.metric = std::abs(shipped - value) / tol;
    c.criterion = "within the law " + std::string(law->coefficient) + " T^" + std::string(law->exponent) +
                  " at 300 K to its printed digits (+-" + shortest(tol) + ")";
    // The register's own 300 K value must round the law (a check of the register, not of NiTCAD).
    if (std::abs(ref::decimal::to_double(row.source_value) - value) > half_unit(row.source_value)) {
        c.comparison.outcome = Outcome::not_evaluated;
        c.comparison.reason = "the register's law_300K value " + std::string(row.source_value) + " does not round the law";
        return c;
    }
    const bool met = c.metric <= 1.0;
    c.comparison = judge(info, {}, met, c.criterion);
    if (c.comparison.outcome == Outcome::fail) {
        c.comparison.outcome = Outcome::finding;
        c.comparison.reason = "register mismatch (15.9; grade " + std::string(grade_name(row.grade)) + "): shipped " +
                              shortest(shipped) + ", law at 300 K " + shortest(value) + " (" + std::string(row.source) +
                              ")";
    }
    return c;
}

}  // namespace

void table_checks(BenchmarkResult& r, std::string_view entry) {
    for (const reg::Row& row : reg::rows) {
        if (row.entry != entry) continue;
        ReferenceInfo info = r.info;
        info.cls = ReferenceClass::M;
        info.provenance = provenance_of(row.grade);
        info.independence = row.grade == reg::Grade::C ? Independence::I0 : Independence::I1;
        info.citation = std::string(row.source);
        info.applicability.clear();
        info.second_method.clear();
        const std::string quantity = "parameter " + std::string(row.binding);
        const std::optional<double> shipped = shipped_value(row.binding);
        if (!shipped) {
            ComparisonRecord c{quantity, std::string(row.unit), 0.0, 0.0, 0.0, 0.0, 0.0, "binding resolved", {}};
            c.comparison = missing_input("the register binding " + std::string(row.binding));
            r.comparisons.push_back(std::move(c));
            continue;
        }
        if (row.kind == reg::Kind::law_300K) {
            r.comparisons.push_back(law_check(info, row, *shipped, quantity));
            continue;
        }
        const bool printed = !row.source_value.empty();
        const double source = printed ? ref::decimal::to_double(row.source_value) : 0.0;
        const bool equal = printed && std::bit_cast<std::uint64_t>(source) == std::bit_cast<std::uint64_t>(*shipped);
        const std::string criterion =
            printed ? "bit-equal to the correctly rounded source decimal " + std::string(row.source_value)
                    : std::string("bit-equal to a source decimal (none read: grade L)");
        ComparisonRecord c{quantity, std::string(row.unit), 0.0, *shipped, source, 0.0, equal ? 0.0 : 1.0, criterion,
                           {}};
        c.comparison = judge(info, {}, equal, criterion);
        if (row.grade == reg::Grade::C && !equal) {
            // Tag C: pinned to the ruling's value (15.9); I0 would make the miss a finding, but a broken pin fails.
            c.comparison.outcome = Outcome::fail;
            c.comparison.reason = "pinned choice changed (tag C): shipped " + shortest(*shipped) + ", ruling " +
                                  std::string(row.source_value) + " (" + std::string(row.source) + ")";
        } else if (c.comparison.outcome == Outcome::fail) {
            // 15.9: a register mismatch is a finding, recorded with its source line.
            c.comparison.outcome = Outcome::finding;
            c.comparison.reason = "register mismatch (15.9; grade " + std::string(grade_name(row.grade)) +
                                  "): shipped " + shortest(*shipped) + ", source " + std::string(row.source_value) +
                                  " (" + std::string(row.source) + ")";
        }
        r.comparisons.push_back(std::move(c));
    }
}

std::vector<double> temperature_grid(double lo, double hi, std::initializer_list<double> beyond) {
    std::vector<double> g;
    const int first = static_cast<int>(std::ceil(64.0 * std::log10(lo) - 1e-9));
    const int last = static_cast<int>(std::floor(64.0 * std::log10(hi) + 1e-9));
    for (int i = first; i <= last; ++i) g.push_back(std::pow(10.0, i / 64.0));
    if (g.empty() || g.front() != lo) g.insert(g.begin(), lo);
    if (g.back() != hi) g.push_back(hi);
    if (lo < 300.0 && 300.0 < hi) g.push_back(300.0);
    for (const double b : beyond) g.push_back(b);
    return g;
}

std::vector<double> log_grid(double lo, double hi, int per) {
    std::vector<double> g;
    const int first = static_cast<int>(std::ceil(per * std::log10(lo) - 1e-9));
    const int last = static_cast<int>(std::floor(per * std::log10(hi) + 1e-9));
    for (int i = first; i <= last; ++i) g.push_back(std::pow(10.0, static_cast<double>(i) / per));
    if (g.empty() || g.front() != lo) g.insert(g.begin(), lo);
    if (g.back() != hi) g.push_back(hi);
    return g;
}

int per_decade(Tier) { return 64; }

Predicate in_range(std::string name, double value, double lo, double hi) {
    return Predicate{std::move(name), value, value < lo ? lo : hi, value >= lo && value <= hi};
}

}  // namespace NiTCAD::verify::benchmark::tvb31
