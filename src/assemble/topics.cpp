// The free functions of the topic headers (Unit 57, DECISIONS.md GE3): each reaches the system's
// data and state through the passkey accessors and calls the SystemData member of its model,
// defined in that model's source file. Those that call the system's own residual
// (terminal_current_resolution, energy_powers) and set_paths, which replaces the data, are defined
// there too.
#include "NiTCAD/assemble/band_edges.hpp"
#include "NiTCAD/assemble/charges.hpp"
#include "NiTCAD/assemble/currents.hpp"
#include "NiTCAD/assemble/energy.hpp"
#include "NiTCAD/assemble/rates.hpp"
#include "NiTCAD/assemble/system.hpp"
#include "NiTCAD/assemble/thermal.hpp"
#include "NiTCAD/assemble/tunnelling.hpp"
#include "system_data.hpp"

namespace NiTCAD::assemble {

namespace {

const SystemData& data(const System& s) noexcept { return SystemData::of(s); }
const SystemState& state(const System& s) noexcept { return SystemData::state_of(s); }

}  // namespace

// --- system.hpp ------------------------------------------------------------------------------

void set_minority_floor(System& system, bool on) noexcept {
    SystemData::state_of(system).minority_floor_ = on;
}
bool minority_floor(const System& system) noexcept { return state(system).minority_floor_; }
void set_log_densities(System& system, bool on) noexcept {
    SystemData::state_of(system).log_densities_ = on;
}
bool log_densities(const System& system) noexcept { return state(system).log_densities_; }

// --- currents.hpp ----------------------------------------------------------------------------

std::vector<std::pair<double, double>> edge_currents(const System& system,
                                                     std::span<const double> x) {
    return data(system).edge_currents(state(system), x);
}
std::vector<double> terminal_currents(const System& system, std::span<const double> x) {
    return data(system).terminal_currents(state(system), x);
}
std::vector<double> edge_flux_terms(const System& system, std::span<const double> x) {
    return data(system).edge_flux_terms(state(system), x);
}
std::vector<double> conduction_currents(const System& system, std::span<const double> x,
                                        const System::TimeStep& step) {
    return data(system).conduction_currents(state(system), x, step);
}
std::vector<CurrentRow> small_signal_currents(const System& system, std::span<const double> x,
                                              std::complex<double> s) {
    return data(system).small_signal_currents(state(system), x, s);
}

// --- thermal.hpp -----------------------------------------------------------------------------

bool electrothermal(const System& system) noexcept { return data(system).electrothermal(); }
std::expected<void, base::Error> set_thermal_bias(System& system,
                                                  std::span<const double> temperature_K) {
    return data(system).set_thermal_bias(SystemData::state_of(system), temperature_K);
}
bool heat_sinks_complete(const System& system) noexcept { return data(system).sinks_complete_; }
bool reaches_heat_sink(const System& system, std::size_t node) {
    return data(system).reaches_heat_sink(node);
}
double heat_scale(const System& system) noexcept { return data(system).heat_scale_; }
std::vector<double> thermal_contact_heat(const System& system, std::span<const double> x,
                                         const System::TimeStep* step) {
    return data(system).thermal_contact_heat(state(system), x, step);
}
double electrical_power(const System& system, std::span<const double> x) {
    return data(system).electrical_power(state(system), x);
}
std::vector<std::pair<double, double>> edge_joule_heat(const System& system,
                                                       std::span<const double> x) {
    return data(system).edge_joule_heat(state(system), x);
}

// --- tunnelling.hpp --------------------------------------------------------------------------

bool tunnelling(const System& system) noexcept { return data(system).tunnelling(); }
TunnelPaths trace_paths(const System& system, std::span<const double> x) {
    return data(system).trace_paths(state(system), x);
}
const TunnelPaths& paths(const System& system) noexcept { return state(system).paths_; }
std::vector<TunnelPathState> path_states(const System& system, std::span<const double> x) {
    return data(system).path_states(state(system), x);
}
bool paths_agree(const System& system, const TunnelPaths& a, const TunnelPaths& b,
                 std::span<const double> x) {
    return data(system).paths_agree(state(system), a, b, x);
}
bool gate_tunnelling(const System& system) noexcept { return data(system).gate_tunnelling(); }
GateTunnellingCurrents gate_tunnelling_currents(const System& system,
                                                std::span<const double> x) {
    return data(system).gate_tunnelling_currents(state(system), x);
}

// --- energy.hpp ------------------------------------------------------------------------------

bool energy_balance(const System& system) noexcept { return data(system).energy_balance(); }
void hold_temperatures(System& system, std::span<const double> x) {
    data(system).hold_temperatures(SystemData::state_of(system), x);
}
void release_temperatures(System& system) noexcept { SystemData::state_of(system).held_.clear(); }
bool temperatures_held(const System& system) noexcept { return !state(system).held_.empty(); }

// --- charges.hpp -----------------------------------------------------------------------------

std::vector<double> trap_occupancies(const System& system, std::span<const double> x,
                                     const System::TimeStep* step) {
    return data(system).trap_occupancies(state(system), x, step);
}
std::vector<double> contact_charges(const System& system, std::span<const double> x,
                                    const System::TimeStep* step) {
    return data(system).contact_charges(state(system), x, step);
}
std::vector<double> gate_charges(const System& system, std::span<const double> x,
                                 const System::TimeStep* step) {
    return data(system).gate_charges(state(system), x, step);
}
std::vector<double> interface_trap_charges(const System& system, std::span<const double> x,
                                           const System::TimeStep* step) {
    return data(system).interface_trap_charges(state(system), x, step);
}

// --- rates.hpp, band_edges.hpp -----------------------------------------------------------------

std::vector<MechanismRate> rates(const System& system, std::span<const double> x) {
    return data(system).rates(state(system), x);
}
BandEdges band_edges(const System& system, std::span<const double> x) {
    return data(system).band_edges(state(system), x);
}

}  // namespace NiTCAD::assemble
