// Insulator parameter set (Unit 15b): a region with no carriers, in which only Poisson's equation
// is solved. Its parameters are the relative permittivity and, for the electrothermal model (Unit
// 23), its thermal data (conduction and heat storage only). The legacy never meshes an insulator
// (its gates are lumped oxides, device/contact.hpp); its SiO2 permittivity, EPS_OX_R = 3.9, is
// carried as the built-in set. Any other insulator is given by its permittivity.
//
// An Insulator is validated once, at construction.
#pragma once

#include <expected>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::physics {

struct InsulatorParameters {
    double eps_r;  // relative permittivity
    // Heat conduction and storage (Unit 23; the thermopower exponents are not read: no carriers).
    // Conductivity 0 (the default): no thermal data.
    ThermalParameters thermal{};
    // Band data for gate tunnelling (Unit 51, DECISIONS.md GD5): the electron affinity chi_ox, the
    // band gap and the tunnelling masses [m0] of electrons, holes and valence-band electrons (EVB;
    // 0 takes the electron mass). All 0 (the default): no band data, and a tunnelling path through
    // the insulator is refused.
    double electron_affinity_eV = 0.0;
    double band_gap_eV = 0.0;
    double tunnel_mass_n = 0.0;
    double tunnel_mass_p = 0.0;
    double tunnel_mass_evb = 0.0;

    bool operator==(const InsulatorParameters&) const = default;
};

// Whether p carries gate-tunnelling band data (a positive gap and both masses).
[[nodiscard]] constexpr bool has_band_data(const InsulatorParameters& p) noexcept {
    return p.band_gap_eV > 0.0 && p.tunnel_mass_n > 0.0 && p.tunnel_mass_p > 0.0;
}

// The tunnelling mass of valence-band electrons: tunnel_mass_evb, or the electron mass when 0.
[[nodiscard]] constexpr double evb_mass(const InsulatorParameters& p) noexcept {
    return p.tunnel_mass_evb > 0.0 ? p.tunnel_mass_evb : p.tunnel_mass_n;
}

// SiO2 (legacy moscap.EPS_OX_R). Thermal data (Unit 23): a conductivity of 0.014 W/(cm K),
// temperature independent (thermal oxide and fused silica, about 1.4 W/(m K) at 300 K); rho c =
// 2.20 g/cm^3 times 0.74 J/(g K) (the heat capacity of SiO2 near 44.4 J/(mol K) at 298 K, NIST-JANAF
// tables, over 60.08 g/mol). No band data: a device without gate tunnelling is unchanged (GC1).
inline constexpr InsulatorParameters silicon_dioxide_parameters{
    .eps_r = 3.9,
    .thermal = {.conductivity_W_cmK = 0.014, .conductivity_exponent = 0.0,
                .heat_capacity_J_cm3K = 1.628, .thermopower_exponent_n = -0.5,
                .thermopower_exponent_p = -0.5}};

// SiO2 with its gate-tunnelling band data (Unit 51, GD5), used where a device tunnels: a gap of
// 9.0 eV and a conduction-band offset of 3.0 eV to silicon (chi = 4.05 eV), so chi_ox = 1.05 eV
// (Gehring, TU Wien 2003, Table 5.1, from Leroy 2003); tunnelling masses 0.50 m0 for electrons and
// holes, the values usually assumed (owner ruling 2026-10-28).
inline constexpr InsulatorParameters silicon_dioxide_tunnelling_parameters{
    .eps_r = silicon_dioxide_parameters.eps_r,
    .thermal = silicon_dioxide_parameters.thermal,
    .electron_affinity_eV = 1.05,
    .band_gap_eV = 9.0,
    .tunnel_mass_n = 0.50,
    .tunnel_mass_p = 0.50,
    .tunnel_mass_evb = 0.0};

class Insulator {
public:
    // Errors: invalid_input if eps_r is not finite and positive, or the thermal parameters are
    // rejected by check_thermal_parameters (the context value is the offending one).
    [[nodiscard]] static std::expected<Insulator, base::Error> create(
        const InsulatorParameters& parameters);

    [[nodiscard]] const InsulatorParameters& parameters() const noexcept { return parameters_; }

private:
    explicit Insulator(const InsulatorParameters& p) noexcept : parameters_(p) {}

    InsulatorParameters parameters_;
};

// Insulator::create(silicon_dioxide_parameters), which cannot fail.
[[nodiscard]] Insulator silicon_dioxide();

}  // namespace NiTCAD::physics
