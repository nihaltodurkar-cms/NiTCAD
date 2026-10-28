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
    // Band data: each value finite and >= 0; given at all (any nonzero), the gap and both masses
    // positive and the affinity finite.
    const double band[] = {parameters.electron_affinity_eV, parameters.band_gap_eV,
                           parameters.tunnel_mass_n, parameters.tunnel_mass_p,
                           parameters.tunnel_mass_evb};
    bool any = false;
    for (double v : band) {
        if (!(std::isfinite(v) && v >= 0.0)) {
            return std::unexpected(base::Error{
                base::ErrorCode::invalid_input,
                "insulator band data (affinity, band gap, tunnelling masses) must be finite and >= 0",
                base::ErrorContext{.index = std::nullopt, .value = v}});
        }
        any = any || v != 0.0;
    }
    if (any && !has_band_data(parameters)) {
        return std::unexpected(base::Error{
            base::ErrorCode::invalid_input,
            "insulator band data needs a positive band gap and positive electron and hole tunnelling masses",
            base::ErrorContext{.index = std::nullopt, .value = parameters.band_gap_eV}});
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
