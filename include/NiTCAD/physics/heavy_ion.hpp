// Heavy-ion (single-event) charge deposition (Unit 65; DECISIONS.md HD1-HD10).
//
// A track is a segment from `start_cm` along the unit vector `direction` of length `length_cm`.
// At a point r its axial coordinate is l = (r - r0) . d and its radial distance
// rho = |r - r0 - l d|. Electron-hole pairs are generated at
//
//     G(r, t) = (LET / q) R(rho) T(t)  for 0 <= l <= L, and 0 elsewhere,
//
// LET the charge deposited per track length [C/cm], R the radial profile normalised over the
// track's cross-section and T the time profile normalised over t:
//   R, Gaussian: 3D exp(-rho^2 / w^2) / (pi w^2); 2D exp(-rho^2 / w^2) / (sqrt(pi) w D), D the
//   device depth the 2D section stands for (HQ1 a); 1D (a track along +x or -x) 1 / A, A the
//   track's cross-section;
//   T(t) = exp(-((t - t0) / s)^2) / (s sqrt(pi)).
// A 1D mesh is per unit cross-section and a 2D mesh per unit depth (mesh.hpp), so the charge a
// track deposits in the device's units (C / cm^(3-D)) is LET L / A in 1D, LET L / D in 2D and
// LET L in 3D.
#pragma once

#include <array>
#include <expected>

#include "NiTCAD/base/error.hpp"

namespace NiTCAD::physics {

struct HeavyIonTrack {
    std::array<double, 3> start_cm{};      // r0; coordinates beyond the dimension are 0
    std::array<double, 3> direction{};     // d; check_heavy_ion_track normalises it
    double length_cm = 0.0;                // L
    double let_C_per_cm = 0.0;             // LET
    double radius_cm = 0.0;                // w (2D and 3D)
    double depth_cm = 0.0;                 // D (2D only)
    double area_cm2 = 0.0;                 // A (1D only)
    double t0_s = 0.0;                     // the pulse's centre
    double sigma_s = 0.0;                  // s, its characteristic width
};

// Silicon's pair energy [eV] and density [g/cm^3] (HD: the MeV cm^2/mg conversion).
inline constexpr double silicon_pair_energy_eV = 3.6;
inline constexpr double silicon_density_g_cm3 = 2.329;

// LET [C/cm] of a LET given as an energy loss per areal density [MeV cm^2/mg], through silicon's
// pair energy and density: LET = x 1e6 eV/MeV (rho 1e3 mg/cm^3) / E_pair q.
[[nodiscard]] double let_from_mev_cm2_mg(double mev_cm2_mg) noexcept;

// The track's inputs in `dimension` (1, 2 or 3), checked; the direction is returned normalised.
// invalid_input when: a value is not finite; length, LET, t0 or sigma not positive; the direction
// zero or with a component beyond the dimension; t0 < 4 sigma (HD6: the steady state at t = 0 has
// no generation); 2D and 3D without a positive radius; 2D without a positive depth; 1D without a
// positive area or not along +x or -x; an area outside 1D or a depth outside 2D.
[[nodiscard]] std::expected<HeavyIonTrack, base::Error> check_heavy_ion_track(HeavyIonTrack track,
                                                                             int dimension);

// (LET / q) R(rho) [cm^-3 per unit of T], 0 outside 0 <= l <= L. Preconditions: the track as
// check_heavy_ion_track returns it for `dimension`.
[[nodiscard]] double heavy_ion_pairs(const HeavyIonTrack& track, int dimension,
                                     const std::array<double, 3>& point_cm) noexcept;

// The integral of T over [t_a, t_b]: (erf((t_b - t0) / s) - erf((t_a - t0) / s)) / 2, computed
// with erfc on the pulse's far side so the tails keep their relative accuracy.
[[nodiscard]] double heavy_ion_time_fraction(const HeavyIonTrack& track, double t_a_s,
                                             double t_b_s) noexcept;

// T(t) [1/s].
[[nodiscard]] double heavy_ion_time_profile(const HeavyIonTrack& track, double t_s) noexcept;

}  // namespace NiTCAD::physics
