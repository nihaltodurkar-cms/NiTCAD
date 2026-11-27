#include "NiTCAD/physics/gallium_oxide.hpp"

#include <cmath>
#include <limits>

#include "temperature.hpp"

namespace NiTCAD::physics {

SemiconductorParameters ga2o3_beta_parameters() {
    const auto dos_300 = [](double mass) { return detail::dos_constant * std::pow(mass, 1.5) * std::pow(300.0, 1.5); };
    SemiconductorParameters p{
        .eps_r = 10.87,  // E || [010] (Fiedler 2019)
        .Eg0_eV = 4.43,  // the 300 K indirect gap, held at every temperature
        .varshni_alpha_eV_per_K = 0.0,
        .varshni_beta_K = 0.0,
        .Nc300 = dos_300(0.28),  // one valley
        .Nv300 = dos_300(3.4),   // placeholder (owner item 8)
        .electron_affinity_eV = 4.00,
        .electron_mobility = {.mu_min = 20.0, .mu_max = 155.0, .N_ref = 1.0e18, .alpha = 0.8, .T_exponent = 0.0},
        // Constant placeholder: mu_min = mu_max (N_ref and alpha then change nothing).
        .hole_mobility = {.mu_min = 10.0, .mu_max = 10.0, .N_ref = 1.0e17, .alpha = 1.0, .T_exponent = 0.0},
        .lifetime = {.tau_n0 = 1.2e-8, .tau_p0 = 1.2e-8, .N_ref = std::numeric_limits<double>::max()},
        .auger = {.Cn = 0.0, .Cp = 0.0},
        .bandgap_narrowing = no_bandgap_narrowing,
        .electron_saturation = {.v_sat_cm_s = 1.5e7, .beta = 0.0, .F_c_V_per_cm = 1.54e5, .gamma = 2.47},
        .hole_saturation = {.v_sat_cm_s = 0.0, .beta = 0.0},  // not given (owner ruling 1)
        .radiative_cm3_s = 0.0,
        .ionization = {.donor_eV = 0.030, .acceptor_eV = 1.10, .donor_degeneracy = 2.0, .acceptor_degeneracy = 4.0},
        .richardson = {.electron = 0.0, .hole = 0.0},
    };
    // Electrons only, one law over all fields (fitted 1-8 MV/cm); no temperature factor.
    p.impact_ionization = {.electron = {.A_low_per_cm = 2.16e6, .B_low_V_per_cm = 1.77e7, .A_high_per_cm = 2.16e6,
                                        .B_high_V_per_cm = 1.77e7, .switch_V_per_cm = 1.0},
                           .hole = {},
                           .phonon_energy_eV = 0.0};
    p.thermal = {.conductivity_W_cmK = 0.228, .conductivity_exponent = 0.0, .heat_capacity_J_cm3K = 0.49 * 5.95};
    return p;
}

}  // namespace NiTCAD::physics
