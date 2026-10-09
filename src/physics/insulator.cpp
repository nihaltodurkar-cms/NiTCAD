#include "NiTCAD/physics/insulator.hpp"

#include <cmath>
#include <utility>

#include "NiTCAD/base/contract.hpp"
#include "NiTCAD/physics/thermal.hpp"

namespace NiTCAD::physics {

std::expected<Insulator, base::Error> Insulator::create(const InsulatorParameters& parameters) {
    if (!(std::isfinite(parameters.eps_r) && parameters.eps_r > 0.0)) {
        return std::unexpected(base::Error{
            base::ErrorCode::invalid_input, "insulator parameter eps_r must be finite and positive",
            base::ErrorContext{.index = std::nullopt, .value = parameters.eps_r}});
    }
    if (auto ok = check_thermal_parameters(parameters.thermal, "insulator"); !ok) {
        return std::unexpected(std::move(ok.error()));
    }
    return Insulator{parameters};
}

Insulator silicon_dioxide() {
    auto m = Insulator::create(silicon_dioxide_parameters);
    NITCAD_EXPECTS(m.has_value());
    return *m;
}

}  // namespace NiTCAD::physics
