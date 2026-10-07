// Interface traps and surface recombination at a semiconductor-insulator interface (Unit 15b; not
// in the legacy, whose only interface-trap model is moscap's lumped capacitance q D_it, M14).
//
// A trap level at energy E_t in the gap, with sheet density N_t [cm^-2] and capture cross-sections
// sigma_n and sigma_p, exchanges carriers with the semiconductor node at the interface. Its
// electron occupancy is the steady-state Shockley-Read-Hall one,
//
//     f = (c_n n + c_p p1) / (c_n (n + n1) + c_p (p + p1)),   c = sigma v_th,
//
// and it recombines carriers at the rate (per unit area)
//
//     U = N_t c_n c_p (n p - n1 p1) / (c_n (n + n1) + c_p (p + p1)).
//
// n1 and p1 are the densities with the Fermi level at E_t: n1 = gamma_n n_ie e^tau and
// p1 = gamma_p n_ie e^-tau, tau = (E_t - E_i) / kT, with gamma the degeneracy factor of the node's
// own density (1 under Boltzmann statistics; physics::fermi_dirac_degeneracy). At equilibrium
// n = gamma_n n_ie e^eta and p = gamma_p n_ie e^-eta with eta = (E_F - E_i) / kT, so
// n1 / n = p / p1 = e^(tau - eta) and f reduces exactly to the Fermi function
// 1 / (1 + e^(tau - eta)), and U to 0, under either statistics. n1 p1 is then the equilibrium
// product of statistics.hpp.
//
// Charge: a donor-like trap is neutral when occupied and +q when empty, N_t (1 - f); an acceptor-
// like trap is neutral when empty and -q when occupied, -N_t f.
//
// A trap band of uniform density D_it [cm^-2 eV^-1] between two energies is integrated by
// composite 6-point Gauss-Legendre quadrature on panels no wider than kT (trap_band_levels). The
// occupancy is smooth on the scale kT (its poles lie pi kT off the real energy axis), so the rule
// converges fast: on a 1.08 eV band at 300 K (252 levels) the equilibrium occupied density is
// within 7.1e-16 of the band's density of the closed-form integral (physics tests).
//
// The plain surface recombination velocity form, for an interface given s_n and s_p instead of
// traps, is SRH with a mid-gap level per unit area: U = (n p - E) / ((n + n_ie) / s_p +
// (p + n_ie) / s_n), physics::srh_recombination with tau_n = 1 / s_n and tau_p = 1 / s_p.
#pragma once

#include <expected>
#include <vector>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::physics {

enum class TrapType : unsigned char { donor, acceptor };

// One trap level. The energy is measured from the intrinsic level of the semiconductor side.
struct TrapLevel {
    TrapType type = TrapType::acceptor;
    double density_cm2 = 0.0;       // N_t [cm^-2], finite and >= 0
    double energy_eV = 0.0;         // E_t - E_i, inside the gap
    double sigma_n_cm2 = 1e-15;     // electron capture cross-section, finite and positive
    double sigma_p_cm2 = 1e-15;     // hole capture cross-section, finite and positive

    bool operator==(const TrapLevel&) const = default;
};

// A band of traps of uniform density between two energies (from the intrinsic level).
struct TrapBand {
    TrapType type = TrapType::acceptor;
    double density_cm2_eV = 0.0;    // D_it [cm^-2 eV^-1], finite and >= 0
    double energy_low_eV = 0.0;     // lower edge, inside the gap
    double energy_high_eV = 0.0;    // upper edge, inside the gap and above the lower one
    double sigma_n_cm2 = 1e-15;
    double sigma_p_cm2 = 1e-15;

    bool operator==(const TrapBand&) const = default;
};

struct InterfaceTraps {
    std::vector<TrapLevel> levels;
    std::vector<TrapBand> bands;
    double thermal_velocity_n_cm_s = 1e7;  // v_th of electrons, finite and positive
    double thermal_velocity_p_cm_s = 1e7;  // v_th of holes, finite and positive

    bool operator==(const InterfaceTraps&) const = default;
};

// Errors (invalid_input; the message names the quantity, the context value is it, and the index the
// level or band): a density not finite and >= 0, a cross-section or thermal velocity not finite and
// positive, an energy not finite or outside the gap of `m` at temperature_K (E_v - E_i to
// E_c - E_i, intrinsic-level offsets from intrinsic_level_depth_eV), or a band whose upper edge is
// not above its lower one. Preconditions as physics::band_gap_eV.
[[nodiscard]] std::expected<void, base::Error> check_interface_traps(const InterfaceTraps& traps,
                                                                     const Semiconductor& m,
                                                                     double temperature_K);

// The quadrature levels of a trap band at temperature_K: per panel six levels at the Gauss-Legendre
// points, each with density D_it times its weight (the levels' densities sum to D_it times the
// band width). Precondition (NITCAD_EXPECTS): the band passes check_interface_traps and
// temperature_K is finite and positive.
[[nodiscard]] std::vector<TrapLevel> trap_band_levels(const TrapBand& band, double temperature_K);

// The equilibrium occupancy of a level, from x = tau - eta = (E_t - E_F) / kT.
struct FermiOccupancy {
    double occupied;  // f = 1 / (1 + e^x)
    double empty;     // 1 - f = 1 / (1 + e^-x), evaluated directly
    double d_eta;     // df / d eta = f (1 - f)
};

[[nodiscard]] FermiOccupancy fermi_occupancy(double x) noexcept;

// The steady-state occupancy and the recombination rate per trap density of one level at the node
// densities n and p, with partials. n1 and p1 as above, dn1_dn = d n1 / dn (n1 d ln gamma_n / dn
// under Fermi-Dirac, else 0) and dp1_dp likewise; cn = sigma_n v_th,n and cp = sigma_p v_th,p.
// Homogeneous of degree zero (f) and one (rate) in the densities. A per-node kernel with no checks;
// finite for n, p >= 0, positive n1, p1, cn, cp.
struct TrapKinetics {
    double occupied, empty;    // f and 1 - f
    double d_occupied_dn, d_occupied_dp;
    double rate;               // U / N_t = cn cp (n p - n1 p1) / D
    double d_rate_dn, d_rate_dp;
};

[[nodiscard]] constexpr TrapKinetics trap_kinetics(double n, double p, double n1, double dn1_dn,
                                                   double p1, double dp1_dp, double cn,
                                                   double cp) noexcept {
    const double D = cn * (n + n1) + cp * (p + p1);
    const double occ = cn * n + cp * p1;
    const double emp = cn * n1 + cp * p;
    const double dD_dn = cn * (1.0 + dn1_dn);
    const double dD_dp = cp * (1.0 + dp1_dp);
    const double excess = n * p - n1 * p1;
    const double k = cn * cp;
    return {
        occ / D,
        emp / D,
        (cn * D - occ * dD_dn) / (D * D),
        (cp * dp1_dp * D - occ * dD_dp) / (D * D),
        k * excess / D,
        k * ((p - dn1_dn * p1) * D - excess * dD_dn) / (D * D),
        k * ((n - n1 * dp1_dp) * D - excess * dD_dp) / (D * D),
    };
}

}  // namespace NiTCAD::physics
