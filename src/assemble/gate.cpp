#include "NiTCAD/assemble/gate.hpp"

#include <cmath>
#include <utility>

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

double gate_intrinsic_offset_V(const device::GateStack& gate, const physics::Semiconductor& m,
                               double temperature_K) {
    // E_c - E_i = Eg / 2 + (k T / 2) ln(Nc / Nv) [eV].
    const double conduction_to_intrinsic =
        0.5 * physics::band_gap_eV(m, temperature_K) +
        0.5 * base::thermal_voltage(temperature_K) *
            std::log(physics::conduction_band_dos(m, temperature_K) /
                     physics::valence_band_dos(m, temperature_K));
    return gate_work_function_eV(gate, m, temperature_K) - m.parameters().electron_affinity_eV -
           conduction_to_intrinsic;
}

double electrode_work_function_eV(const device::Electrode& electrode, double temperature_K) {
    const physics::Semiconductor poly = physics::silicon();
    const double chi = poly.parameters().electron_affinity_eV;
    switch (electrode.kind) {
        case device::GateElectrode::n_poly:
            return chi;
        case device::GateElectrode::p_poly:
            return chi + physics::band_gap_eV(poly, temperature_K);
        case device::GateElectrode::metal:
            return electrode.work_function_eV;
    }
    NITCAD_EXPECTS(false);  // a validated device has no other electrode
    return 0.0;
}

GateTerm gate_term(const device::GateStack& gate, const physics::Semiconductor& m, double area_cm,
                   int dimension, const Scaling& scaling) {
    const double area = area_cm / std::pow(scaling.L_D, dimension - 1);
    return {permittivity_ratio(gate.oxide_relative_permittivity, scaling) * area /
                (gate.oxide_thickness_cm / scaling.L_D),
            gate_intrinsic_offset_V(gate, m, scaling.temperature_K) / scaling.V_T,
            gate.fixed_charge_cm2 * area / (scaling.Ns * scaling.L_D)};
}

GateNodes::GateNodes(std::vector<std::int32_t> contact, std::vector<GateTerm> term, double V_T)
    : contact_(std::move(contact)), term_(std::move(term)), V_T_(V_T) {
    NITCAD_EXPECTS(contact_.size() == term_.size() && V_T > 0.0);
    psi_gate_.assign(contact_.size(), 0.0);
    for (std::size_t i = 0; i < contact_.size(); ++i) {
        if (contact_[i] >= 0) psi_gate_[i] = -term_[i].offset;  // zero bias
    }
}

void GateNodes::set_bias(std::span<const double> bias_V) {
    for (std::size_t i = 0; i < contact_.size(); ++i) {
        if (contact_[i] < 0) continue;
        const auto c = static_cast<std::size_t>(contact_[i]);
        NITCAD_EXPECTS(c < bias_V.size());
        psi_gate_[i] = bias_V[c] / V_T_ - term_[i].offset;
    }
}

std::vector<double> GateNodes::charges(std::span<const double> state, std::size_t stride,
                                       std::size_t contacts) const {
    NITCAD_EXPECTS(state.size() == stride * contact_.size());
    std::vector<double> charge(contacts, 0.0);
    for (std::size_t i = 0; i < contact_.size(); ++i) {
        if (contact_[i] < 0) continue;
        const auto c = static_cast<std::size_t>(contact_[i]);
        NITCAD_EXPECTS(c < contacts);
        charge[c] += term_[i].coupling * (psi_gate_[i] - state[stride * i]);
    }
    return charge;
}

}  // namespace NiTCAD::assemble
