// State transfer between meshes (Unit 17, slice 3; DECISIONS.md Unit 17 AD5).
#include "NiTCAD/solve/adapt_transfer.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "NiTCAD/base/constants.hpp"

namespace NiTCAD::solve {
namespace {

base::Error invalid(std::string message) { return {base::ErrorCode::invalid_input, std::move(message), std::nullopt}; }

double safe_exp(double x) { return std::exp(std::clamp(x, -700.0, 700.0)); }

}  // namespace

std::expected<results::NodeFields, base::Error> transfer_state(const device::Device& from,
                                                               const mesh::PointLocator& where,
                                                               const results::NodeFields& s,
                                                               const results::BandDiagram& bands,
                                                               const device::Device& to, TransferReport* report,
                                                               TransferMode mode) {
    const std::size_t no = from.mesh().node_count();
    const auto sized = [&](const std::vector<double>& v) { return v.size() == no; };
    if (!sized(s.potential_V) || !sized(s.n_cm3) || !sized(s.p_cm3) || !sized(bands.electron_fermi_eV) ||
        !sized(bands.hole_fermi_eV)) {
        return std::unexpected(invalid("a transferred state needs its fields and bands one per old node"));
    }
    const bool thermal = !s.temperature_K.empty();
    const bool carriers_hot = !s.electron_temperature_K.empty() || !s.hole_temperature_K.empty();
    if ((thermal && !sized(s.temperature_K)) || (!s.electron_temperature_K.empty() && !sized(s.electron_temperature_K)) ||
        (!s.hole_temperature_K.empty() && !sized(s.hole_temperature_K))) {
        return std::unexpected(invalid("a transferred state's temperatures must be one per old node"));
    }
    if (where.node_count() != no) return std::unexpected(invalid("the locator is not over the old mesh"));
    if (from.regions().size() != to.regions().size()) {
        return std::unexpected(invalid("transfer needs the same regions on both meshes"));
    }
    for (std::size_t r = 0; r < from.regions().size(); ++r) {
        if (device::is_insulator(from.regions()[r]) != device::is_insulator(to.regions()[r])) {
            return std::unexpected(invalid("transfer needs the same regions on both meshes"));
        }
    }

    // Per old node: the quantities to interpolate.
    const double T0 = from.temperature_K();
    const auto vt_old = [&](std::size_t k) { return base::thermal_voltage(thermal ? s.temperature_K[k] : T0); };
    std::vector<double> cn(no, 0.0), cp(no, 0.0);
    for (std::size_t k = 0; k < no; ++k) {
        if (from.is_insulator(static_cast<mesh::NodeId>(k))) continue;
        const double vt = vt_old(k);
        cn[k] = std::log(s.n_cm3[k]) - (s.potential_V[k] + bands.electron_fermi_eV[k]) / vt;
        cp[k] = std::log(s.p_cm3[k]) + (s.potential_V[k] + bands.hole_fermi_eV[k]) / vt;
    }
    const auto old_region = from.node_region();
    std::vector<std::vector<std::size_t>> by_region(from.regions().size());
    for (std::size_t k = 0; k < no; ++k) by_region[static_cast<std::size_t>(old_region[k])].push_back(k);

    const auto pts = to.mesh().points();
    const auto new_region = to.node_region();
    const std::size_t nn = pts.size();
    results::NodeFields out;
    out.potential_V.assign(nn, 0.0);
    out.n_cm3.assign(nn, 0.0);
    out.p_cm3.assign(nn, 0.0);
    if (thermal) out.temperature_K.assign(nn, T0);
    if (!s.electron_temperature_K.empty()) out.electron_temperature_K.assign(nn, T0);
    if (!s.hole_temperature_K.empty()) out.hole_temperature_K.assign(nn, T0);
    TransferReport rep;
    rep.nodes = nn;
    for (std::size_t i = 0; i < nn; ++i) {
        const auto r = static_cast<std::size_t>(new_region[i]);
        if (by_region[r].empty()) return std::unexpected(invalid("a region of the new mesh has no node on the old one"));
        const mesh::Interpolant loc = where.locate(pts[i]);
        rep.outside += loc.inside ? 0 : 1;
        std::array<std::size_t, 8> node{};
        std::array<double, 8> w{};
        int m = 0;
        double sum = 0.0;
        for (int j = 0; j < loc.count; ++j) {
            const auto k = static_cast<std::size_t>(loc.nodes[static_cast<std::size_t>(j)]);
            if (static_cast<std::size_t>(old_region[k]) != r) continue;
            node[static_cast<std::size_t>(m)] = k;
            w[static_cast<std::size_t>(m)] = loc.weights[static_cast<std::size_t>(j)];
            sum += w[static_cast<std::size_t>(m)];
            ++m;
        }
        if (!(sum > 0.0)) {
            // The nearest old node of the region (lowest index on ties).
            ++rep.other_region;
            std::size_t best = by_region[r][0];
            double bd = std::numeric_limits<double>::infinity();
            const auto opts = from.mesh().points();
            for (const std::size_t k : by_region[r]) {
                const double dx = opts[k][0] - pts[i][0], dy = opts[k][1] - pts[i][1], dz = opts[k][2] - pts[i][2];
                const double d = dx * dx + dy * dy + dz * dz;
                if (d < bd) {
                    bd = d;
                    best = k;
                }
            }
            node[0] = best;
            w[0] = 1.0;
            m = 1;
            sum = 1.0;
        }
        const auto interp = [&](auto&& value) {
            double v = 0.0;
            for (int j = 0; j < m; ++j) v += w[static_cast<std::size_t>(j)] * value(node[static_cast<std::size_t>(j)]);
            return v / sum;
        };
        const double psi = interp([&](std::size_t k) { return s.potential_V[k]; });
        out.potential_V[i] = psi;
        if (thermal) out.temperature_K[i] = interp([&](std::size_t k) { return s.temperature_K[k]; });
        if (carriers_hot) {
            if (!s.electron_temperature_K.empty()) {
                out.electron_temperature_K[i] = interp([&](std::size_t k) { return s.electron_temperature_K[k]; });
            }
            if (!s.hole_temperature_K.empty()) {
                out.hole_temperature_K[i] = interp([&](std::size_t k) { return s.hole_temperature_K[k]; });
            }
        }
        if (to.is_insulator(static_cast<mesh::NodeId>(i))) continue;
        const double vt = base::thermal_voltage(thermal ? out.temperature_K[i] : to.temperature_K());
        if (mode == TransferMode::exponential_fit) {
            const double nn_ = interp([&](std::size_t k) { return s.n_cm3[k] * std::exp((psi - s.potential_V[k]) / (2.0 * vt)); });
            const double dn = interp([&](std::size_t k) { return std::exp((s.potential_V[k] - psi) / (2.0 * vt)); });
            const double np_ = interp([&](std::size_t k) { return s.p_cm3[k] * std::exp((s.potential_V[k] - psi) / (2.0 * vt)); });
            const double dp = interp([&](std::size_t k) { return std::exp((psi - s.potential_V[k]) / (2.0 * vt)); });
            out.n_cm3[i] = nn_ / dn;
            out.p_cm3[i] = np_ / dp;
            continue;
        }
        const double fn = interp([&](std::size_t k) { return bands.electron_fermi_eV[k]; });
        const double fp = interp([&](std::size_t k) { return bands.hole_fermi_eV[k]; });
        const double c_n = interp([&](std::size_t k) { return cn[k]; });
        const double c_p = interp([&](std::size_t k) { return cp[k]; });
        out.n_cm3[i] = safe_exp(c_n + (psi + fn) / vt);
        out.p_cm3[i] = safe_exp(c_p - (psi + fp) / vt);
    }
    if (report != nullptr) *report = rep;
    return out;
}

}  // namespace NiTCAD::solve
