// Insulator parameter set (Unit 15b): a region with no carriers, in which only Poisson's equation
// is solved. Its one parameter is the relative permittivity. The legacy never meshes an insulator
// (its gates are lumped oxides, device/contact.hpp); its SiO2 permittivity, EPS_OX_R = 3.9, is
// carried as the built-in set. Any other insulator is given by its permittivity.
//
// An Insulator is validated once, at construction.
#pragma once

#include <expected>

#include "NiTCAD/base/error.hpp"

namespace NiTCAD::physics {

struct InsulatorParameters {
    double eps_r;  // relative permittivity

    bool operator==(const InsulatorParameters&) const = default;
};

// SiO2 (legacy moscap.EPS_OX_R).
inline constexpr InsulatorParameters silicon_dioxide_parameters{.eps_r = 3.9};

class Insulator {
public:
    // Errors: invalid_input if eps_r is not finite and positive (the context value is it).
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
