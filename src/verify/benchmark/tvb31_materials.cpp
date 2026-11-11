// TVB-31.30-31.31, the material tables (Unit 55, TVR-3 step 9; DECISIONS 15.3, 15.9, 15.19).
//
// Code under test: the shipped parameter sets physics::silicon_dioxide_parameters and
// silicon_dioxide_tunnelling_parameters (insulator.hpp); germanium_parameters, gallium_arsenide_parameters,
// indium_gallium_arsenide_parameters, silicon_carbide_4h_parameters and algaas_parameters(x) (semiconductor.hpp).
// Reference: the register rows of verify/reference/model_parameters.hpp, typed from the sources as printed.
//   31.30 the SiO2 rows (15.9's table rule; MI-Q26).
//   31.31 the rows of the four constant sets; AlGaAs's composition laws at x = 0, 0.05, ..., 0.45 against Ioffe's
//         printed laws, to half a unit in the last printed digit of each coefficient (MI-Q23); and the coverage of
//         MI-Q25: a field of a constant set with no row, not 0 and not the legacy silicon value is reported
//         not_evaluated (missing input).
#include <array>
#include <bit>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "NiTCAD/physics/semiconductor.hpp"
#include "NiTCAD/verify/benchmark/registry.hpp"
#include "NiTCAD/verify/reference/decimal.hpp"
#include "NiTCAD/verify/reference/model_parameters.hpp"
#include "tvb31_common.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

namespace reg = verify::reference::model_parameters;
namespace dec = verify::reference::decimal;
using tvb31::DD;
using tvb31::in_range;

double half_unit(std::string_view printed) { return 0.5 * std::pow(10.0, dec::parse(printed).exponent); }

// The shipped value of an AlGaAs law's binding at x.
double algaas_value(std::string_view binding, double x) {
    physics::SemiconductorParameters p = physics::algaas_parameters(x).value();
    if (binding == "algaas.band_gap_300K") {
        return physics::band_gap_eV(physics::Semiconductor::create(p).value(), 300.0);
    }
    const double* f = tvb31::bound_field(p, binding, "algaas");
    return f != nullptr ? *f : NAN;
}

// The law and its tolerance at x (MI-Q23): the change for half a unit in each coefficient's last printed digit.
std::pair<DD, double> law_at(const reg::CompositionLaw& law, double x) {
    const DD X(x);
    if (law.form == reg::LawForm::polynomial) {
        DD v(0.0), xk(1.0);
        double tol = 0.0;
        for (const std::string_view c : law.c) {
            if (c.empty()) continue;
            v = v + dec::to_dd(c) * xk;
            tol += half_unit(c) * std::abs(xk.hi);
            xk = xk * X;
        }
        return {v, tol};
    }
    const DD A = dec::to_dd(law.scale);
    const DD u = dec::to_dd(law.c[0]) + dec::to_dd(law.c[1]) * X;
    const DD v = A * u * sqrt(u);
    const double tol = std::abs(v.hi) * half_unit(law.scale) / std::abs(A.hi) +
                       1.5 * std::abs((A * sqrt(u)).hi) * (half_unit(law.c[0]) + half_unit(law.c[1]) * std::abs(x));
    return {v, tol};
}

void composition_checks(BenchmarkResult& r) {
    for (const reg::CompositionLaw& law : reg::composition_laws) {
        ReferenceInfo info = r.info;
        info.provenance = Provenance::PC;
        info.citation = std::string(law.source);
        info.applicability = {"x within the law's printed range"};
        const double x_max = dec::to_double(law.x_max);
        for (int i = 0; i <= 9; ++i) {
            const double x = i / 20.0;
            const double sim = algaas_value(law.binding, x);
            const auto [q, tol] = law_at(law, x);
            const double diff = std::abs((DD(sim) - q).hi);
            const std::string criterion = "within the printed law (" + std::string(law.source) +
                                          ") to its printed digits (+-" + shortest(tol) + ") (15.19, MI-Q23)";
            ComparisonRecord c{"composition " + std::string(law.binding), std::string(law.unit), x, sim,
                               q.to_double(), 0x1p-80 * std::abs(q.hi), diff / tol, criterion, {}};
            const std::vector<Predicate> pr{in_range("x within the law's printed range", x, 0.0, x_max)};
            c.comparison = judge(info, pr, std::isfinite(sim) && diff <= tol, criterion);
            r.comparisons.push_back(std::move(c));
        }
    }
}

bool has_row(std::string_view binding) {
    for (const reg::Row& row : reg::rows) {
        if (row.binding == binding) return true;
    }
    return false;
}

// MI-Q25: every field of a constant set has a row, is 0, or is bit-equal to the legacy silicon value.
void coverage_checks(BenchmarkResult& r) {
    physics::SemiconductorParameters si = physics::legacy_silicon_parameters;
    for (const auto& [set, params] : tvb31::non_silicon_sets()) {
        physics::SemiconductorParameters p = params;
        for (const std::string_view field : tvb31::field_names()) {
            const std::string binding = std::string(set) + "." + std::string(field);
            const double* v = tvb31::bound_field(p, binding, set);
            const double* s = tvb31::bound_field(si, "legacy_silicon." + std::string(field), "legacy_silicon");
            if (v == nullptr || s == nullptr) {
                ComparisonRecord c{"coverage " + binding, "", 0.0, 0.0, 0.0, 0.0, 0.0, "field resolved", {}};
                c.comparison = missing_input("the field " + binding + " (resolver)");
                r.comparisons.push_back(std::move(c));
                continue;
            }
            if (has_row(binding) || *v == 0.0 ||
                std::bit_cast<std::uint64_t>(*v) == std::bit_cast<std::uint64_t>(*s)) {
                continue;
            }
            ComparisonRecord c{"coverage " + binding, "", 0.0, *v, 0.0, 0.0, 0.0,
                               "a register row, 0, or the legacy silicon value (15.19, MI-Q25)", {}};
            c.comparison = missing_input("a register row for " + binding + " (MI-Q25)");
            r.comparisons.push_back(std::move(c));
        }
    }
}

}  // namespace

BenchmarkResult run_tvb31_30(Tier) {
    BenchmarkResult r;
    r.id = "TVB-31.30";
    r.title = "SiO2 data";
    r.info = tvb31::entry_info(r.id, Provenance::PD, "the register's SiO2 rows (15.19, MI-Q26)",
                               "each shipped SiO2 parameter bit-equal to the correctly rounded printed decimal",
                               {}, {1, 9, 12}, "verify/reference/model_parameters.hpp",
                               "15.9's table rule; tag C pinned");
    tvb31::table_checks(r, "31.30");
    apply_known_findings(r, known_findings());
    return r;
}

BenchmarkResult run_tvb31_31(Tier) {
    BenchmarkResult r;
    r.id = "TVB-31.31";
    r.title = "Non-silicon material sets";
    r.info = tvb31::entry_info(r.id, Provenance::PD,
                               "the register's 31.31 rows and AlGaAs composition laws, typed from the Ioffe NSM "
                               "archive (15.19, MI-Q23, MI-Q24)",
                               "each shipped parameter bit-equal to the printed decimal; algaas_parameters(x) within "
                               "the printed laws to their printed digits",
                               {}, {1, 5, 6, 7, 9}, "verify/reference/model_parameters.hpp",
                               "15.9's table rule; MI-Q23's law rule; MI-Q25's coverage");
    tvb31::table_checks(r, "31.31");
    composition_checks(r);
    coverage_checks(r);
    apply_known_findings(r, known_findings());
    return r;
}

}  // namespace NiTCAD::verify::benchmark
