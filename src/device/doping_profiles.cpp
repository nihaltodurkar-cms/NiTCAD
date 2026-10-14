#include "NiTCAD/device/doping_profiles.hpp"

#include <cmath>
#include <cstddef>
#include <optional>
#include <string>

namespace NiTCAD::device {

namespace {

base::Error invalid(std::string message, std::size_t index) {
    return {base::ErrorCode::invalid_input, std::move(message),
            base::ErrorContext{.index = index, .value = std::nullopt}};
}

std::optional<base::Error> check_profile(const DopingProfile& p, int dimension, std::size_t k) {
    if (!(std::isfinite(p.peak_cm3) && p.peak_cm3 >= 0.0)) return invalid("doping profile: peak not finite and >= 0", k);
    if (p.axis < 0 || p.axis >= dimension) return invalid("doping profile: axis outside the mesh dimension", k);
    if (p.shape != DepthShape::constant && !(std::isfinite(p.width) && p.width > 0.0)) {
        return invalid("doping profile: width not finite and positive", k);
    }
    if (!std::isfinite(p.center) || !std::isfinite(p.lateral) || p.lateral < 0.0) {
        return invalid("doping profile: center or lateral not finite, or lateral < 0", k);
    }
    for (int b = 0; b < dimension; ++b) {
        const auto i = static_cast<std::size_t>(b);
        if (std::isnan(p.lo[i]) || std::isnan(p.hi[i]) || p.lo[i] > p.hi[i]) {
            return invalid("doping profile: a window with lo > hi or a NaN bound", k);
        }
    }
    return std::nullopt;
}

double window(double x, double lo, double hi, double lateral) {
    if (lateral == 0.0) {
        const double below = x < lo ? 0.0 : (x == lo ? 0.5 : 1.0);
        const double above = x > hi ? 0.0 : (x == hi ? 0.5 : 1.0);
        return below * above;
    }
    const double a = std::isinf(lo) ? 2.0 : std::erfc((lo - x) / lateral);
    const double b = std::isinf(hi) ? 0.0 : std::erfc((hi - x) / lateral);
    return 0.5 * (a - b);
}

}  // namespace

double evaluate(const DopingProfile& p, const mesh::Point& x) {
    const auto axis = static_cast<std::size_t>(p.axis);
    const double d = x[axis];
    double depth = 1.0;
    if (p.shape == DepthShape::gaussian) {
        const double u = (d - p.center) / p.width;
        depth = std::exp(-0.5 * u * u);
    } else if (p.shape == DepthShape::erfc) {
        depth = 0.5 * std::erfc((d - p.center) / p.width);
    }
    double lateral = 1.0;
    for (std::size_t b = 0; b < 3; ++b) {
        if (b == axis || (std::isinf(p.lo[b]) && std::isinf(p.hi[b]))) continue;
        lateral *= window(x[b], p.lo[b], p.hi[b], p.lateral);
    }
    return p.peak_cm3 * depth * lateral;
}

std::expected<Doping, base::Error> sample(std::span<const DopingProfile> profiles, const mesh::Mesh& m,
                                          std::span<const RegionId> node_region, std::span<const Region> regions) {
    for (std::size_t k = 0; k < profiles.size(); ++k) {
        if (auto e = check_profile(profiles[k], m.dimension(), k)) return std::unexpected(*e);
    }
    const std::size_t n = m.node_count();
    if (node_region.size() != n) return std::unexpected(invalid("node_region needs one region per node", n));
    Doping out{std::vector<double>(n, 0.0), std::vector<double>(n, 0.0)};
    for (std::size_t i = 0; i < n; ++i) {
        const RegionId r = node_region[i];
        if (r < 0 || static_cast<std::size_t>(r) >= regions.size()) {
            return std::unexpected(invalid("node region out of range", i));
        }
        // Profiles beyond the dimension are evaluated with zero coordinates there.
        if (is_insulator(regions[static_cast<std::size_t>(r)])) continue;
        for (const DopingProfile& p : profiles) {
            const double v = evaluate(p, m.points()[i]);
            (p.dopant == Dopant::donor ? out.donors : out.acceptors)[i] += v;
        }
    }
    return out;
}

}  // namespace NiTCAD::device
