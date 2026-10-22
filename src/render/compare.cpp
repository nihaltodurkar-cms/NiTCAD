// Run comparison (compare.hpp).
#include "NiTCAD/render/compare.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace NiTCAD::render {

std::optional<std::size_t> matching_point(const results::BiasPoint& a, std::span<const results::BiasPoint> b) noexcept {
    for (std::size_t k = 0; k < b.size(); ++k) {
        const auto& v = b[k].bias_V;
        if (v.size() != a.bias_V.size()) continue;
        bool equal = true;
        for (std::size_t c = 0; c < v.size() && equal; ++c) equal = std::abs(v[c] - a.bias_V[c]) <= 1e-9;
        if (equal) return k;
    }
    return std::nullopt;
}

std::expected<DifferenceField, base::Error> difference_field(const device::Device& device_a,
                                                             const device::Device& device_b, const Geometry& geometry,
                                                             const results::BiasPoint& a, const results::BiasPoint& b,
                                                             Field f) {
    auto va = field_values(device_a, geometry, a, f);
    if (!va) return std::unexpected(std::move(va.error()));
    auto vb = field_values(device_b, geometry, b, f);
    if (!vb) return std::unexpected(std::move(vb.error()));
    const std::string_view key = field_info(f).key;
    DifferenceField d;
    d.values.resize(va->size());
    const bool density = f == Field::electrons || f == Field::holes || f == Field::donors || f == Field::acceptors;
    const bool energy = f == Field::conduction_band || f == Field::valence_band || f == Field::electron_fermi ||
                        f == Field::hole_fermi;
    for (std::size_t i = 0; i < d.values.size(); ++i) {
        const double x = (*va)[i], y = (*vb)[i];
        if (density) {
            d.values[i] = x > 0.0 && y > 0.0 ? std::log10(y) - std::log10(x) : std::numeric_limits<double>::quiet_NaN();
        } else if (f == Field::potential || energy) {
            d.values[i] = (y - x) * 1e3;
        } else {
            d.values[i] = y - x;
        }
    }
    const std::string delta = "\xCE\x94";  // Δ
    if (density) d.title = delta + "log\xE2\x82\x81\xE2\x82\x80 " + std::string(key) + " (decades)";
    else if (f == Field::potential) d.title = delta + "\xCF\x88 (mV)";
    else if (energy) d.title = delta + std::string(key) + " (meV)";
    else if (f == Field::net_doping) d.title = delta + "(ND \xE2\x88\x92 NA) (cm\xE2\x81\xBB\xC2\xB3)";
    else d.title = delta + std::string(field_info(f).title);
    return d;
}

ColorScale difference_scale(std::span<const double> values) {
    double m = 0.0;
    for (const double v : values) {
        if (std::isfinite(v)) m = std::max(m, std::abs(v));
    }
    if (!(m > 0.0)) m = 1.0;
    return *ColorScale::create(ScaleKind::linear, -m, m);
}

}  // namespace NiTCAD::render
