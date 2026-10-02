#include "NiTCAD/assemble/contact_bias.hpp"

#include <cmath>
#include <cstddef>

namespace NiTCAD::assemble {

std::expected<void, base::Error> check_contact_bias(std::span<const device::ContactKind> kinds,
                                                    std::span<const double> bias_V,
                                                    bool thermal_equilibrium) {
    if (bias_V.size() != kinds.size()) {
        return std::unexpected(base::Error{base::ErrorCode::invalid_input,
                                           "one bias per contact is required", std::nullopt});
    }
    for (std::size_t c = 0; c < bias_V.size(); ++c) {
        const base::ErrorContext at{.index = c, .value = bias_V[c]};
        if (!std::isfinite(bias_V[c])) {
            return std::unexpected(
                base::Error{base::ErrorCode::invalid_input, "contact bias is not finite", at});
        }
        if (thermal_equilibrium && kinds[c] == device::ContactKind::ohmic && bias_V[c] != 0.0) {
            return std::unexpected(base::Error{
                base::ErrorCode::invalid_input,
                "thermal equilibrium needs every ohmic contact at 0 V", at});
        }
    }
    return {};
}

}  // namespace NiTCAD::assemble
