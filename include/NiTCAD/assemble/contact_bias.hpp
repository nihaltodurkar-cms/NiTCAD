// The rule for a set of contact biases, shared by both assemblers' set_bias and by the solve layer,
// which checks every bias point before it solves anything.
#pragma once

#include <expected>
#include <span>

#include "NiTCAD/base/error.hpp"
#include "NiTCAD/device/contact.hpp"

namespace NiTCAD::assemble {

// bias_V holds one bias in V per contact, in the order of `kinds`. Errors (invalid_input): the
// sizes differ; a bias is not finite; with thermal_equilibrium, an ohmic contact's bias is not 0
// (thermal equilibrium has one Fermi level, that of the ohmic contacts; only gates may be biased).
// The context index is the contact and the value its bias.
[[nodiscard]] std::expected<void, base::Error> check_contact_bias(
    std::span<const device::ContactKind> kinds, std::span<const double> bias_V,
    bool thermal_equilibrium);

}  // namespace NiTCAD::assemble
