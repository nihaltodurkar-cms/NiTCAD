// Refinement indicators and marking (Unit 17, slice 2; DECISIONS.md Unit 17 AD2, AD3, AR1).
#include "NiTCAD/solve/adapt_indicators.hpp"

#include <algorithm>
#include <cmath>
#include <string>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::solve {
namespace {

using base::ErrorCode;

base::Error invalid(std::string message) {
    return {ErrorCode::invalid_input, std::move(message), std::nullopt};
}

constexpr double ln10 = 2.302585092994046;

bool positive_finite(double x) { return std::isfinite(x) && x > 0.0; }

// The generation part of a node's impact and band-to-band rates [cm^-3 s^-1] (rates are R - G).
double generation_at(const results::RateFields& r, std::size_t i) {
    double g = 0.0;
    if (i < r.impact.size()) g += std::max(0.0, -r.impact[i]);
    double btbt = 0.0;
    if (i < r.band_to_band_n.size()) btbt = std::max(btbt, -r.band_to_band_n[i]);
    if (i < r.band_to_band_p.size()) btbt = std::max(btbt, -r.band_to_band_p[i]);
    return g + btbt;
}

}  // namespace

std::expected<EdgeIndicators, base::Error> edge_indicators(const device::Device& device,
                                                           const results::NodeFields* state,
                                                           const results::RateFields* rates,
                                                           const IndicatorOptions& o) {
    if (!positive_finite(o.potential_vt) || !positive_finite(o.density_decades) || !positive_finite(o.doping_decades) ||
        !positive_finite(o.generation_share)) {
        return std::unexpected(invalid("indicator thresholds must be finite and positive"));
    }
    const mesh::Mesh& m = device.mesh();
    const std::size_t nn = m.node_count();
    if (state != nullptr) {
        if (state->potential_V.size() != nn || state->n_cm3.size() != nn || state->p_cm3.size() != nn ||
            (!state->temperature_K.empty() && state->temperature_K.size() != nn)) {
            return std::unexpected(invalid("indicator state fields must be one per mesh node"));
        }
    }
    const auto edges = m.edges();
    EdgeIndicators out;
    for (auto& v : out.raw) v.assign(edges.size(), 0.0);
    const auto on = [&](Indicator k) { return o.enabled[static_cast<std::size_t>(k)]; };

    // Per node: the signed log of net doping (AR1), and the generation share.
    std::vector<char> semi(nn, 0);
    std::size_t n_semi = 0;
    std::vector<double> L(nn, 0.0);
    const double T = device.temperature_K();
    for (std::size_t i = 0; i < nn; ++i) {
        const auto id = static_cast<mesh::NodeId>(i);
        if (device.is_insulator(id)) continue;
        semi[i] = 1;
        ++n_semi;
        const double ni = physics::intrinsic_density(device.material(id), T);
        const double N = device.net_doping(id);
        const double mag = std::log(std::max(std::abs(N), ni) / ni);
        L[i] = N < 0.0 ? -mag : mag;
    }
    std::vector<double> share;
    if (on(Indicator::generation) && rates != nullptr &&
        (!rates->impact.empty() || !rates->band_to_band_n.empty() || !rates->band_to_band_p.empty())) {
        const auto vol = m.volumes();
        double total = 0.0;
        share.assign(nn, 0.0);
        for (std::size_t i = 0; i < nn; ++i) {
            if (semi[i] == 0) continue;
            share[i] = generation_at(*rates, i) * vol[i];
            total += share[i];
        }
        if (total > 0.0) {
            for (double& s : share) s = s / total * static_cast<double>(n_semi);
        } else {
            share.clear();
        }
    }

    for (std::size_t e = 0; e < edges.size(); ++e) {
        const auto a = static_cast<std::size_t>(edges[e].first), b = static_cast<std::size_t>(edges[e].second);
        if (semi[a] == 0 || semi[b] == 0) continue;
        if (on(Indicator::doping)) out.raw[static_cast<std::size_t>(Indicator::doping)][e] = std::abs(L[a] - L[b]) / ln10;
        if (!share.empty()) out.raw[static_cast<std::size_t>(Indicator::generation)][e] = std::max(share[a], share[b]);
        if (state == nullptr) continue;
        if (on(Indicator::potential)) {
            const double Te = state->temperature_K.empty() ? T : 0.5 * (state->temperature_K[a] + state->temperature_K[b]);
            out.raw[static_cast<std::size_t>(Indicator::potential)][e] =
                std::abs(state->potential_V[a] - state->potential_V[b]) / base::thermal_voltage(Te);
        }
        const auto decades = [](double x, double y) {
            return positive_finite(x) && positive_finite(y) ? std::abs(std::log(x) - std::log(y)) / ln10 : 0.0;
        };
        if (on(Indicator::electron_density)) {
            out.raw[static_cast<std::size_t>(Indicator::electron_density)][e] = decades(state->n_cm3[a], state->n_cm3[b]);
        }
        if (on(Indicator::hole_density)) {
            out.raw[static_cast<std::size_t>(Indicator::hole_density)][e] = decades(state->p_cm3[a], state->p_cm3[b]);
        }
    }
    return out;
}

std::expected<IndicatorField, base::Error> compute_indicators(const device::Device& device,
                                                            const results::NodeFields* state,
                                                            const results::RateFields* rates,
                                                            const IndicatorOptions& o) {
    auto raw = edge_indicators(device, state, rates, o);
    if (!raw) return std::unexpected(raw.error());
    const mesh::Mesh& m = device.mesh();
    const auto edges = m.edges();
    const std::array<double, indicator_count> threshold{o.potential_vt, o.density_decades, o.density_decades,
                                                        o.doping_decades, o.generation_share};
    IndicatorField f;
    f.excess.assign(m.node_count(), 0.0);
    f.dominant.assign(m.node_count(), Indicator::potential);
    std::array<std::vector<char>, indicator_count> touched;
    for (auto& t : touched) t.assign(m.node_count(), 0);
    for (std::size_t e = 0; e < edges.size(); ++e) {
        for (std::size_t k = 0; k < indicator_count; ++k) {
            const double x = std::max(0.0, raw->raw[k][e] / threshold[k] - 1.0);
            if (!(x > 0.0)) continue;
            for (const mesh::NodeId v : {edges[e].first, edges[e].second}) {
                const auto i = static_cast<std::size_t>(v);
                touched[k][i] = 1;
                // Kinds go in increasing order: a later kind replaces only when strictly larger.
                if (x > f.excess[i]) {
                    f.excess[i] = x;
                    f.dominant[i] = static_cast<Indicator>(k);
                }
            }
        }
    }
    for (std::size_t k = 0; k < indicator_count; ++k) {
        f.nodes_with_excess[k] = static_cast<std::size_t>(std::ranges::count(touched[k], 1));
    }
    return f;
}

std::expected<IndicatorField, base::Error> unite(std::span<const IndicatorField> fields) {
    if (fields.empty()) return std::unexpected(invalid("no indicator fields to unite"));
    IndicatorField out = fields[0];
    for (std::size_t j = 1; j < fields.size(); ++j) {
        const IndicatorField& f = fields[j];
        if (f.excess.size() != out.excess.size() || f.dominant.size() != out.dominant.size()) {
            return std::unexpected(invalid("indicator fields of different sizes"));
        }
        for (std::size_t i = 0; i < f.excess.size(); ++i) {
            if (f.excess[i] > out.excess[i]) {
                out.excess[i] = f.excess[i];
                out.dominant[i] = f.dominant[i];
            }
        }
        for (std::size_t k = 0; k < indicator_count; ++k) {
            out.nodes_with_excess[k] = std::max(out.nodes_with_excess[k], f.nodes_with_excess[k]);
        }
    }
    return out;
}

std::expected<Marking, base::Error> mark_nodes(const device::Device& device, const IndicatorField& field,
                                               std::span<const double> weights, double theta) {
    if (!(theta > 0.0 && theta <= 1.0)) return std::unexpected(invalid("theta must be in (0, 1]"));
    const std::size_t nn = device.mesh().node_count();
    if (field.excess.size() != nn || field.dominant.size() != nn) {
        return std::unexpected(invalid("an indicator field must have one entry per device node"));
    }
    if (!weights.empty()) {
        if (weights.size() != nn) return std::unexpected(invalid("goal weights must be one per node"));
        for (const double w : weights) {
            if (!(std::isfinite(w) && w >= 0.0)) return std::unexpected(invalid("goal weights must be finite and >= 0"));
        }
    }
    const auto ew = [&](std::size_t i) { return field.excess[i] * (weights.empty() ? 1.0 : weights[i]); };
    double top = 0.0;
    for (std::size_t i = 0; i < nn; ++i) top = std::max(top, ew(i));
    Marking mk;
    mk.marked.assign(nn, 0);
    mk.by_region.assign(device.regions().size(), 0);
    mk.threshold = theta * top;
    if (!(top > 0.0)) return mk;
    const auto region = device.node_region();
    for (std::size_t i = 0; i < nn; ++i) {
        const double v = ew(i);
        if (v > 0.0 && v >= mk.threshold) {
            mk.marked[i] = 1;
            ++mk.count;
            ++mk.by_kind[static_cast<std::size_t>(field.dominant[i])];
            ++mk.by_region[static_cast<std::size_t>(region[i])];
        }
    }
    return mk;
}

}  // namespace NiTCAD::solve
