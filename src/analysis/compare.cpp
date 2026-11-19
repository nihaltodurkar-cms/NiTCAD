// compare_sweeps (compare.hpp; DECISIONS.md Unit 60, XC1-XC2).
#include "NiTCAD/analysis/compare.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <limits>
#include <span>

namespace NiTCAD::analysis {

namespace {

// Accumulates one quantity: the largest difference and the largest ratio to the allowance.
class Tally {
public:
    explicit Tally(std::string name) { q_.name = std::move(name); }

    void add(double a, double b, double allowance, std::size_t point, std::size_t index) {
        constexpr double inf = std::numeric_limits<double>::infinity();
        const bool both_nan = std::isnan(a) && std::isnan(b);  // insulator nodes' bands; never here today
        const double d = std::abs(a - b);
        double ratio = 0.0;
        if (a == b || both_nan) {
            ratio = 0.0;
        } else if (std::isnan(d)) {
            ratio = inf;  // a NaN on one side only
        } else {
            q_.max_difference = std::max(q_.max_difference, d);
            ratio = allowance > 0.0 ? d / allowance : inf;
        }
        if (ratio > q_.max_ratio) {
            q_.max_ratio = ratio;
            q_.point = point;
            q_.index = index;
        }
        ++q_.compared;
    }

    [[nodiscard]] bool used() const noexcept { return q_.compared > 0; }
    [[nodiscard]] QuantityComparison take() && { return std::move(q_); }

private:
    QuantityComparison q_;
};

double relative_allowance(double a, double b, double relative, double floor) {
    return relative * std::max(std::abs(a), std::abs(b)) + floor;
}

double at_or_zero(std::span<const double> v, std::size_t i) { return i < v.size() ? v[i] : 0.0; }

double largest(std::span<const double> v) {
    double m = 0.0;
    for (const double x : v)
        if (std::isfinite(x)) m = std::max(m, std::abs(x));
    return m;
}

// Why two sweeps cannot be compared value by value; empty when they can.
std::string structure_mismatch(const results::Sweep& a, const results::Sweep& b) {
    if (a.points.size() != b.points.size())
        return std::format("point count {} against {}", a.points.size(), b.points.size());
    if (a.stopped.has_value() != b.stopped.has_value())
        return a.stopped ? "the first run stopped early, the second did not"
                         : "the second run stopped early, the first did not";
    for (std::size_t k = 0; k < a.points.size(); ++k) {
        const auto& p = a.points[k];
        const auto& q = b.points[k];
        if (p.bias_V.size() != q.bias_V.size())
            return std::format("point {}: contact count {} against {}", k, p.bias_V.size(), q.bias_V.size());
        if (p.bias_V != q.bias_V) return std::format("point {}: the bias differs", k);
        if (p.fields.potential_V.size() != q.fields.potential_V.size())
            return std::format("point {}: node count {} against {}", k, p.fields.potential_V.size(),
                               q.fields.potential_V.size());
        const auto sizes_differ = [](const std::vector<double>& x, const std::vector<double>& y) {
            return x.size() != y.size();
        };
        const struct {
            const std::vector<double>& x;
            const std::vector<double>& y;
            const char* name;
        } fields[] = {
            {p.fields.n_cm3, q.fields.n_cm3, "n_cm3"},
            {p.fields.p_cm3, q.fields.p_cm3, "p_cm3"},
            {p.fields.temperature_K, q.fields.temperature_K, "temperature_K"},
            {p.fields.electron_temperature_K, q.fields.electron_temperature_K, "electron_temperature_K"},
            {p.fields.hole_temperature_K, q.fields.hole_temperature_K, "hole_temperature_K"},
            {p.fields.quantum_potential_n_eV, q.fields.quantum_potential_n_eV, "quantum_potential_n_eV"},
            {p.fields.quantum_potential_p_eV, q.fields.quantum_potential_p_eV, "quantum_potential_p_eV"},
            {p.terminal_current, q.terminal_current, "terminal_current"},
            {p.terminal_current_resolution, q.terminal_current_resolution, "terminal_current_resolution"},
            {p.weighted_current, q.weighted_current, "weighted_current"},
            {p.weighted_current_resolution, q.weighted_current_resolution, "weighted_current_resolution"},
            {p.gate_charge, q.gate_charge, "gate_charge"},
            {p.interface_trap_charge, q.interface_trap_charge, "interface_trap_charge"},
        };
        for (const auto& f : fields)
            if (sizes_differ(f.x, f.y))
                return std::format("point {}: {} has {} values against {}", k, f.name, f.x.size(), f.y.size());
    }
    return {};
}

}  // namespace

SweepComparison compare_sweeps(const results::Sweep& a, const results::Sweep& b,
                               const CompareTolerances& t) {
    SweepComparison out;
    out.identity_equal = a.run.input_identity == b.run.input_identity;
    out.mesh_digest_equal = a.run.mesh_digest == b.run.mesh_digest;
    out.reason = structure_mismatch(a, b);
    if (!out.reason.empty()) {
        out.outcome = CompareOutcome::not_comparable;
        return out;
    }

    Tally current("terminal_current"), weighted("weighted_current"), gate("gate_charge"),
        traps("interface_trap_charge"), psi("potential_V"), n("n_cm3"), p("p_cm3"), lattice("temperature_K"),
        electron_t("electron_temperature_K"), hole_t("hole_temperature_K"),
        lambda_n("quantum_potential_n_eV"), lambda_p("quantum_potential_p_eV");

    for (std::size_t k = 0; k < a.points.size(); ++k) {
        const auto& x = a.points[k];
        const auto& y = b.points[k];
        for (std::size_t c = 0; c < x.terminal_current.size(); ++c) {
            const double floor = std::max(at_or_zero(x.terminal_current_resolution, c),
                                          at_or_zero(y.terminal_current_resolution, c));
            current.add(x.terminal_current[c], y.terminal_current[c],
                        relative_allowance(x.terminal_current[c], y.terminal_current[c], t.relative, floor), k, c);
        }
        for (std::size_t c = 0; c < x.weighted_current.size(); ++c) {
            const double floor = std::max(at_or_zero(x.weighted_current_resolution, c),
                                          at_or_zero(y.weighted_current_resolution, c));
            weighted.add(x.weighted_current[c], y.weighted_current[c],
                         relative_allowance(x.weighted_current[c], y.weighted_current[c], t.relative, floor), k, c);
        }
        for (std::size_t c = 0; c < x.gate_charge.size(); ++c)
            gate.add(x.gate_charge[c], y.gate_charge[c],
                     relative_allowance(x.gate_charge[c], y.gate_charge[c], t.relative, t.charge_floor), k, c);
        for (std::size_t c = 0; c < x.interface_trap_charge.size(); ++c)
            traps.add(x.interface_trap_charge[c], y.interface_trap_charge[c],
                      relative_allowance(x.interface_trap_charge[c], y.interface_trap_charge[c], t.relative,
                                         t.charge_floor),
                      k, c);

        const auto& fx = x.fields;
        const auto& fy = y.fields;
        for (std::size_t i = 0; i < fx.potential_V.size(); ++i)
            psi.add(fx.potential_V[i], fy.potential_V[i], t.potential_V, k, i);
        const double dens_floor =
            t.density_floor * std::max({largest(fx.n_cm3), largest(fx.p_cm3), largest(fy.n_cm3), largest(fy.p_cm3)});
        for (std::size_t i = 0; i < fx.n_cm3.size(); ++i)
            n.add(fx.n_cm3[i], fy.n_cm3[i], relative_allowance(fx.n_cm3[i], fy.n_cm3[i], t.relative, dens_floor), k, i);
        for (std::size_t i = 0; i < fx.p_cm3.size(); ++i)
            p.add(fx.p_cm3[i], fy.p_cm3[i], relative_allowance(fx.p_cm3[i], fy.p_cm3[i], t.relative, dens_floor), k, i);
        for (std::size_t i = 0; i < fx.temperature_K.size(); ++i)
            lattice.add(fx.temperature_K[i], fy.temperature_K[i], t.temperature_K, k, i);
        for (std::size_t i = 0; i < fx.electron_temperature_K.size(); ++i)
            electron_t.add(fx.electron_temperature_K[i], fy.electron_temperature_K[i], t.temperature_K, k, i);
        for (std::size_t i = 0; i < fx.hole_temperature_K.size(); ++i)
            hole_t.add(fx.hole_temperature_K[i], fy.hole_temperature_K[i], t.temperature_K, k, i);
        // Unit 62: the quantum potentials (eV) at the potential's tolerance (V).
        for (std::size_t i = 0; i < fx.quantum_potential_n_eV.size(); ++i)
            lambda_n.add(fx.quantum_potential_n_eV[i], fy.quantum_potential_n_eV[i], t.potential_V, k, i);
        for (std::size_t i = 0; i < fx.quantum_potential_p_eV.size(); ++i)
            lambda_p.add(fx.quantum_potential_p_eV[i], fy.quantum_potential_p_eV[i], t.potential_V, k, i);
    }

    for (Tally* q : {&current, &weighted, &gate, &traps, &psi, &n, &p, &lattice, &electron_t, &hole_t, &lambda_n, &lambda_p})
        if (q->used()) out.quantities.push_back(std::move(*q).take());
    out.outcome = std::ranges::all_of(out.quantities, &QuantityComparison::pass) ? CompareOutcome::match
                                                                                  : CompareOutcome::differ;
    return out;
}

std::string_view to_string(CompareOutcome outcome) noexcept {
    switch (outcome) {
    case CompareOutcome::match: return "match";
    case CompareOutcome::differ: return "differ";
    case CompareOutcome::not_comparable: return "not comparable";
    }
    return "unknown";
}

}  // namespace NiTCAD::analysis
