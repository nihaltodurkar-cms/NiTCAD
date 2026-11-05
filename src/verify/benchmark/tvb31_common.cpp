// Shared machinery of the TVB-31 entries (Unit 55, TVR-3; tvb31_common.hpp).
#include "tvb31_common.hpp"

#include <bit>
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
    return nullptr;
}

std::optional<double> shipped_value(std::string_view binding) {
    physics::SemiconductorParameters legacy = physics::legacy_silicon_parameters;
    physics::SemiconductorParameters si = physics::silicon_parameters;
    if (const double* v = bound_field(legacy, binding, "legacy_silicon")) return *v;
    if (const double* v = bound_field(si, binding, "silicon")) return *v;
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

Predicate in_range(std::string name, double value, double lo, double hi) {
    return Predicate{std::move(name), value, value < lo ? lo : hi, value >= lo && value <= hi};
}

}  // namespace NiTCAD::verify::benchmark::tvb31
