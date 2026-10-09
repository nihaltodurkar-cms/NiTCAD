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

    bool operator==(const InsulatorParameters&) const = default;
};

// SiO2 (legacy moscap.EPS_OX_R). Thermal data (Unit 23): a conductivity of 0.014 W/(cm K),
// temperature independent (thermal oxide and fused silica, about 1.4 W/(m K) at 300 K); rho c =
// 2.20 g/cm^3 times 0.74 J/(g K) (the heat capacity of SiO2 near 44.4 J/(mol K) at 298 K, NIST-JANAF
// tables, over 60.08 g/mol).
inline constexpr InsulatorParameters silicon_dioxide_parameters{
    .eps_r = 3.9,
    .thermal = {.conductivity_W_cmK = 0.014, .conductivity_exponent = 0.0,
                .heat_capacity_J_cm3K = 1.628, .thermopower_exponent_n = -0.5,
                .thermopower_exponent_p = -0.5}};

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
