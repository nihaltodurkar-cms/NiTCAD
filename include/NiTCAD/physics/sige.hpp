// Si(1-x)Ge(x) parameter sets (Unit 68, slice 1, MAD7; DECISIONS.md amendment 5).
//
// Source M7, the Ioffe NSM SiGe pages (F. Schaffler's chapter of Levinshtein, Rumyantsev, Shur,
// eds., Properties of Advanced Semiconductor Materials, Wiley 2001), with the owner's rulings
// S1-S11 and F1-F4 of 2026-11-26. Every set is anchored to NiTCAD's silicon: at x = 0 the band,
// mobility and recombination values are silicon's exactly.
//
// Band edges, all at 300 K, against relaxed silicon. A film of Ge fraction x grown on a relaxed
// substrate of fraction x_s (x_s = x: relaxed):
//     a(x) = 5.431 + 0.20 x + 0.027 x^2 [A];  eps_par = (a(x_s) - a(x)) / a(x);
//     eps_perp = -2 (C12 / C11) eps_par, C12 / C11 = s / (1 - s), s = 0.278 - 0.005 x (Poisson);
//     weighted-average valence offset (Rieger-Vogl): dEv_av = (0.47 - 0.06 x) x (composition only);
//     spin-orbit splitting D0 = 0.044 + 0.246 x;
//     valence splitting dE = 2 b (eps_perp - eps_par), heavy hole D0 / 3 - dE / 2, light hole
//     -D0 / 6 + dE / 4 + (D0^2 + D0 dE + 9 dE^2 / 4)^(1/2) / 2, the top their maximum;
//     conduction: the Delta gap's hydrostatic shift (Xd + Xu / 3 - a)(2 eps_par + eps_perp), the
//     Delta2 / Delta4 splits 2 Xu (eps_perp - eps_par) / 3 and -Xu (eps_perp - eps_par) / 3, the
//     lowest valley taken;
//     deformation potentials (Van de Walle and Martin 1986, theory) linear in x.
// The relaxed gap is Eg_silicon + (-0.41 x + 0.008 x^2) (Ioffe's 300 K fit with Si's 1.12 eV
// replaced by silicon's own gap, ruling F3a).
#pragma once

#include <expected>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::physics {

// The largest Ge fraction (ruling S3: Ioffe's mobility lines hold to 0.3).
inline constexpr double sige_max_fraction = 0.30;

// A film's band edges against relaxed silicon at 300 K [eV].
struct SiGeBands {
    double in_plane_strain;          // eps_par
    double average_valence_offset;   // dEv_av (Rieger-Vogl), the quantity MAC9 checks (S10)
    double valence_edge_offset;      // top of the film's valence band minus silicon's
    double gap_change;               // the film's gap minus silicon's
    double conduction_edge_offset;   // valence_edge_offset + gap_change
};

// Errors: invalid_input when x or x_substrate lies outside [0, sige_max_fraction].
[[nodiscard]] std::expected<SiGeBands, base::Error> sige_bands(double x, double x_substrate);

// `sige <x>` (relaxed) and `sige_strained <x>` (pseudomorphic on (001) silicon). Errors: as
// sige_bands.
[[nodiscard]] std::expected<SemiconductorParameters, base::Error> sige_parameters(double x, bool strained);

// `si_strained <x_substrate>`: silicon pseudomorphic on relaxed SiGe; silicon's set with the
// strained band edges (S11: silicon's mobility). Errors: as sige_bands.
[[nodiscard]] std::expected<SemiconductorParameters, base::Error> strained_silicon_parameters(double x_substrate);

}  // namespace NiTCAD::physics
