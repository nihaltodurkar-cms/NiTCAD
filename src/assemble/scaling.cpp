#include "NiTCAD/assemble/scaling.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/physics/semiconductor.hpp"

namespace NiTCAD::assemble {

std::expected<Scaling, base::Error> make_scaling(const device::Device& device,
                                                 std::optional<double> Ns_override) {
    if (Ns_override && !(std::isfinite(*Ns_override) && *Ns_override > 0.0)) {
        return std::unexpected(base::Error{
            base::ErrorCode::invalid_input, "Ns override must be finite and positive",
            base::ErrorContext{.index = std::nullopt, .value = *Ns_override}});
    }
    const double T = device.temperature_K();
    const physics::Semiconductor& reference = device.material(0);
    Scaling s{};
    s.temperature_K = T;
    s.V_T = base::thermal_voltage(T);
    s.eps_F_per_cm = reference.parameters().eps_r * base::eps0_F_per_cm;
    s.n_i = physics::intrinsic_density(reference, T);
    if (Ns_override) {
        s.Ns = *Ns_override;
    } else {
        double max_doping = 0.0;
        for (std::size_t i = 0; i < device.mesh().node_count(); ++i) {
            max_doping = std::max(max_doping,
                                  std::abs(device.net_doping(static_cast<mesh::NodeId>(i))));
        }
        s.Ns = std::max(max_doping, s.n_i);
    }
    s.L_D = std::sqrt(s.eps_F_per_cm * s.V_T / (base::q_C * s.Ns));
    s.D0 = 1.0;
    s.J0 = base::q_C * s.D0 * s.Ns / s.L_D;
    s.R0 = s.D0 * s.Ns / (s.L_D * s.L_D);
    return s;
}

}  // namespace NiTCAD::assemble
