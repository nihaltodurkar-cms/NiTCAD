#include "NiTCAD/physics/heavy_ion.hpp"

#include <cmath>
#include <cstddef>
#include <numbers>
#include <string>
#include <utility>

#include "NiTCAD/base/constants.hpp"

namespace NiTCAD::physics {

namespace {

std::unexpected<base::Error> refuse(std::string message, double value) {
    return std::unexpected(base::Error{base::ErrorCode::invalid_input,
                                       "heavy-ion track: " + std::move(message),
                                       base::ErrorContext{.index = std::nullopt, .value = value}});
}

}  // namespace

double let_from_mev_cm2_mg(double mev_cm2_mg) noexcept {
    // pairs per cm = x [MeV cm^2/mg] rho [mg/cm^3] / E_pair [MeV]; times q.
    const double rho_mg_cm3 = silicon_density_g_cm3 * 1e3;
    const double pair_MeV = silicon_pair_energy_eV * 1e-6;
    return mev_cm2_mg * rho_mg_cm3 / pair_MeV * base::q_C;
}

std::expected<HeavyIonTrack, base::Error> check_heavy_ion_track(HeavyIonTrack t, int dimension) {
    if (dimension < 1 || dimension > 3) return refuse("dimension must be 1, 2 or 3", dimension);
    const auto d = static_cast<std::size_t>(dimension);
    for (std::size_t k = 0; k < 3; ++k) {
        if (!std::isfinite(t.start_cm[k]) || !std::isfinite(t.direction[k])) {
            return refuse("start and direction must be finite", static_cast<double>(k));
        }
        if (k >= d && (t.start_cm[k] != 0.0 || t.direction[k] != 0.0)) {
            return refuse("start or direction has a component beyond the mesh's dimension",
                          static_cast<double>(k));
        }
    }
    for (const auto& [name, v] : {std::pair{"length", t.length_cm}, std::pair{"LET", t.let_C_per_cm},
                                  std::pair{"t0", t.t0_s}, std::pair{"sigma", t.sigma_s}}) {
        if (!(std::isfinite(v) && v > 0.0)) return refuse(std::string(name) + " must be finite and positive", v);
    }
    for (const auto& [name, v] : {std::pair{"radius", t.radius_cm}, std::pair{"depth", t.depth_cm},
                                  std::pair{"area", t.area_cm2}}) {
        if (!(std::isfinite(v) && v >= 0.0)) return refuse(std::string(name) + " must be finite and not negative", v);
    }
    const double norm = std::hypot(t.direction[0], t.direction[1], t.direction[2]);
    if (!(norm > 0.0)) return refuse("the direction must not be zero", norm);
    for (double& c : t.direction) c /= norm;
    if (t.t0_s < 4.0 * t.sigma_s) {
        return refuse("t0 must be at least 4 sigma (the start at t = 0 has no generation)", t.t0_s);
    }
    if (dimension == 1) {
        if (!(t.area_cm2 > 0.0)) return refuse("a 1D track needs a positive area", t.area_cm2);
        if (t.depth_cm != 0.0) return refuse("depth is for 2D tracks only", t.depth_cm);
        return t;
    }
    if (t.area_cm2 != 0.0) return refuse("area is for 1D tracks only", t.area_cm2);
    if (!(t.radius_cm > 0.0)) return refuse("the radius must be positive", t.radius_cm);
    if (dimension == 2 && !(t.depth_cm > 0.0)) return refuse("a 2D track needs a positive depth", t.depth_cm);
    if (dimension == 3 && t.depth_cm != 0.0) return refuse("depth is for 2D tracks only", t.depth_cm);
    return t;
}

double heavy_ion_pairs(const HeavyIonTrack& t, int dimension, const std::array<double, 3>& r) noexcept {
    const std::array<double, 3> v{r[0] - t.start_cm[0], r[1] - t.start_cm[1], r[2] - t.start_cm[2]};
    const double l = v[0] * t.direction[0] + v[1] * t.direction[1] + v[2] * t.direction[2];
    if (!(l >= 0.0 && l <= t.length_cm)) return 0.0;
    const double per_length = t.let_C_per_cm / base::q_C;  // pairs per cm of track
    if (dimension == 1) return per_length / t.area_cm2;
    const double px = v[0] - l * t.direction[0];
    const double py = v[1] - l * t.direction[1];
    const double pz = v[2] - l * t.direction[2];
    const double u2 = (px * px + py * py + pz * pz) / (t.radius_cm * t.radius_cm);
    const double w = t.radius_cm;
    const double norm = dimension == 2 ? std::sqrt(std::numbers::pi) * w * t.depth_cm
                                       : std::numbers::pi * w * w;
    return per_length * std::exp(-u2) / norm;
}

double heavy_ion_time_fraction(const HeavyIonTrack& t, double ta, double tb) noexcept {
    const double a = (ta - t.t0_s) / t.sigma_s;
    const double b = (tb - t.t0_s) / t.sigma_s;
    if (a >= 0.0) return 0.5 * (std::erfc(a) - std::erfc(b));    // both after the centre
    if (b <= 0.0) return 0.5 * (std::erfc(-b) - std::erfc(-a));  // both before it
    return 1.0 - 0.5 * (std::erfc(b) + std::erfc(-a));
}

double heavy_ion_time_profile(const HeavyIonTrack& t, double ts) noexcept {
    const double u = (ts - t.t0_s) / t.sigma_s;
    return std::exp(-u * u) / (t.sigma_s * std::sqrt(std::numbers::pi));
}

}  // namespace NiTCAD::physics
