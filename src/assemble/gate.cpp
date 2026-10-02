#include "NiTCAD/assemble/gate.hpp"

#include <cmath>

#include "NiTCAD/base/constants.hpp"
#include "NiTCAD/base/contract.hpp"

namespace NiTCAD::assemble {

double gate_work_function_eV(const device::GateStack& gate, const physics::Semiconductor& m,
                             double temperature_K) {
    const double chi = m.parameters().electron_affinity_eV;
    switch (gate.electrode) {
        case device::GateElectrode::n_poly:
            return chi;
        case device::GateElectrode::p_poly:
            return chi + physics::band_gap_eV(m, temperature_K);
        case device::GateElectrode::metal:
            return gate.work_function_eV;
    }
    NITCAD_EXPECTS(false);  // a validated device has no other electrode
    return 0.0;
}

double gate_midgap_offset_V(const device::GateStack& gate, const physics::Semiconductor& m,
                            double temperature_K) {
    return gate_work_function_eV(gate, m, temperature_K) - m.parameters().electron_affinity_eV -
           0.5 * physics::band_gap_eV(m, temperature_K);
}

GateTerm gate_term(const device::GateStack& gate, const physics::Semiconductor& m, double area_cm,
                   int dimension, const Scaling& scaling) {
    const double area = area_cm / std::pow(scaling.L_D, dimension - 1);
    const double eps_ratio =
        gate.oxide_relative_permittivity * base::eps0_F_per_cm / scaling.eps_F_per_cm;
    return {eps_ratio * area / (gate.oxide_thickness_cm / scaling.L_D),
            gate_midgap_offset_V(gate, m, scaling.temperature_K) / scaling.V_T,
            gate.fixed_charge_cm2 * area / (scaling.Ns * scaling.L_D)};
}

}  // namespace NiTCAD::assemble
