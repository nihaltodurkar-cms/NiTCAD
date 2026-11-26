#include "NiTCAD/physics/sige.hpp"

#include <algorithm>
#include <cmath>
#include <format>
#include <optional>

namespace NiTCAD::physics {
namespace {

double lattice_A(double x) { return 5.431 + 0.20 * x + 0.027 * x * x; }
double spin_orbit_eV(double x) { return 0.044 + 0.246 * x; }
// Van de Walle and Martin (1986) theory values for Si and Ge, linear in x (the Ioffe page's advice).
double lerp(double si, double ge, double x) { return si + (ge - si) * x; }

// The heavy- and light-hole energies above the weighted average valence band [eV].
double valence_top_eV(double x, double eps_par, double eps_perp) {
    const double d0 = spin_orbit_eV(x);
    const double dE = 2.0 * lerp(-2.35, -2.55, x) * (eps_perp - eps_par);
    const double heavy = d0 / 3.0 - dE / 2.0;
    const double light = -d0 / 6.0 + dE / 4.0 + 0.5 * std::sqrt(d0 * d0 + d0 * dE + 2.25 * dE * dE);
    return std::max(heavy, light);
}

// Sasso-Rinaldi's composition interpolation at 1e16 cm^-3 and 300 K (owner-supplied, ruling F2b),
// as a ratio to its x = 0 value: v(x) / v(0), xn = x / 0.3.
double saturation_ratio(double x, double v_si, double v_sige, double c_v, double alpha) {
    const double xn = x / 0.30;
    const double v = v_si * (1.0 - xn) + v_sige * xn + std::pow(1.0 - xn, alpha) * std::pow(xn, alpha) * c_v;
    return v / v_si;
}

std::optional<base::Error> check_fraction(const char* what, double x) {
    if (std::isfinite(x) && x >= 0.0 && x <= sige_max_fraction) return std::nullopt;
    return base::Error{base::ErrorCode::invalid_input,
                       std::format("{} {} lies outside [0, {}] (the SiGe mobility data stop at 0.3)", what, x,
                                   sige_max_fraction)};
}

}  // namespace

std::expected<SiGeBands, base::Error> sige_bands(double x, double x_substrate) {
    if (auto e = check_fraction("the Ge fraction", x)) return std::unexpected(std::move(*e));
    if (auto e = check_fraction("the substrate Ge fraction", x_substrate)) return std::unexpected(std::move(*e));
    const double eps_par = (lattice_A(x_substrate) - lattice_A(x)) / lattice_A(x);
    const double s = 0.278 - 0.005 * x;
    const double eps_perp = -2.0 * s / (1.0 - s) * eps_par;
    const double average = (0.47 - 0.06 * x) * x;
    const double top = valence_top_eV(x, eps_par, eps_perp);
    const double hydrostatic = lerp(1.72, 1.31, x) * (2.0 * eps_par + eps_perp);
    const double xu = lerp(9.16, 9.42, x);
    const double valley = std::min(2.0 / 3.0 * xu * (eps_perp - eps_par), -1.0 / 3.0 * xu * (eps_perp - eps_par));
    // The relaxed film's conduction edge sits D0 / 3 + Eg above its average valence band.
    const double relaxed_gap_change = -0.41 * x + 0.008 * x * x;
    const double gap_change = relaxed_gap_change + spin_orbit_eV(x) / 3.0 + hydrostatic + valley - top;
    const double valence = average + top - spin_orbit_eV(0.0) / 3.0;
    return SiGeBands{eps_par, average, valence, gap_change, valence + gap_change};
}

std::expected<SemiconductorParameters, base::Error> sige_parameters(double x, bool strained) {
    auto bands = sige_bands(x, strained ? 0.0 : x);
    if (!bands) return std::unexpected(std::move(bands.error()));
    SemiconductorParameters p = silicon_parameters;
    p.eps_r = 11.7 + 4.5 * x;
    p.Eg0_eV += bands->gap_change;
    p.electron_affinity_eV -= bands->conduction_edge_offset;
    // F3a: the Varshni change of each end member from 300 K (Si 4.73e-4 eV/K, 636 K; Ge 4.774e-4
    // eV/K, 235 K, Palankovski and Selberherr as the owner gave them).
    p.alloy_gap = {.fraction = x, .alpha_a_eV_per_K = 4.73e-4, .beta_a_K = 636.0, .alpha_b_eV_per_K = 4.774e-4,
                   .beta_b_K = 235.0};
    // S1, S2: silicon's Caughey-Thomas set; mu_max on Ioffe's lines scaled to silicon at x = 0.
    p.electron_mobility.mu_max *= (1396.0 - 4315.0 * x) / 1396.0;
    p.hole_mobility.mu_max *= (450.0 - 865.0 * x) / 450.0;
    // F1b and the x = 0 anchor: silicon's conduction mass fit; the valence fit scaled by the hole
    // DOS mass (0.81 - 0.47 x) m0, relative to its x = 0 value, to the 3/2.
    const double hole_dos = std::pow((0.81 - 0.47 * x) / 0.81, 1.5);
    for (double& c : p.dos_mass.valence) c *= hole_dos;
    // S6: band-gap narrowing (1 - 0.5 x) silicon's (owner's fit); F2b saturation velocities; F4a
    // acceptor (boron) 45 - 34.6 x meV, donor (phosphorus) 45 - 33 x meV.
    p.bandgap_narrowing.E0_eV *= 1.0 - 0.5 * x;
    p.electron_saturation.v_sat_cm_s *= saturation_ratio(x, 9.954e6, 6.093e6, -5.262e6, 0.74665);
    p.hole_saturation.v_sat_cm_s *= saturation_ratio(x, 1.1412e7, 8.1459e6, -2.22e6, 0.8223);
    p.ionization.acceptor_eV = 0.045 - 0.0346 * x;
    p.ionization.donor_eV = 0.045 - 0.033 * x;
    // No SiGe data in hand: impact ionization and band-to-band tunnelling off (all zero), and no
    // thermal conductivity (electrothermal refused) until the owner gives the temperature exponent
    // of Ioffe's k = 0.046 + 0.084 x. Heat capacity from Ioffe's specific heat, molar mass and
    // density.
    p.impact_ionization = {};
    p.band_to_band = {};
    const double molar_g = 28.0855 * (1.0 - x) + 72.630 * x;
    const double density_g_cm3 = 2.329 + 3.493 * x - 0.499 * x * x;
    p.thermal.conductivity_W_cmK = 0.0;
    p.thermal.conductivity_exponent = 0.0;
    p.thermal.heat_capacity_J_cm3K = (19.6 + 2.9 * x) / molar_g * density_g_cm3;
    return p;
}

std::expected<SemiconductorParameters, base::Error> strained_silicon_parameters(double x_substrate) {
    auto bands = sige_bands(0.0, x_substrate);
    if (!bands) return std::unexpected(std::move(bands.error()));
    SemiconductorParameters p = silicon_parameters;
    p.Eg0_eV += bands->gap_change;
    p.electron_affinity_eV -= bands->conduction_edge_offset;
    return p;
}

}  // namespace NiTCAD::physics
