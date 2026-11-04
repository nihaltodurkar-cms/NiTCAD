// Copies a device into the consistency layer's plain data. See device_data.hpp.
#include "NiTCAD/verify/benchmark/device_data.hpp"

namespace NiTCAD::verify::benchmark {

consistency::DeviceData device_data(const device::Device& device) {
    consistency::DeviceData d;
    d.temperature_K = device.temperature_K();
    d.donors.assign(device.donors().begin(), device.donors().end());
    d.acceptors.assign(device.acceptors().begin(), device.acceptors().end());
    const std::size_t n = device.mesh().node_count();
    d.relative_permittivity.resize(n);
    d.insulator.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        const auto node = static_cast<mesh::NodeId>(i);
        d.relative_permittivity[i] = device.relative_permittivity(node);
        d.insulator[i] = device.is_insulator(node) ? 1 : 0;
    }
    for (const auto& c : device.contacts()) {
        consistency::ContactData k;
        k.name = c.name;
        switch (c.kind) {
            case device::ContactKind::ohmic: k.kind = consistency::ContactKind::ohmic; break;
            case device::ContactKind::gate: k.kind = consistency::ContactKind::gate; break;
            case device::ContactKind::electrode: k.kind = consistency::ContactKind::electrode; break;
            case device::ContactKind::schottky: k.kind = consistency::ContactKind::schottky; break;
        }
        k.nodes = c.nodes;
        d.contacts.push_back(std::move(k));
    }
    for (const auto& f : device.interfaces()) {
        if (!f.traps.levels.empty() || !f.traps.bands.empty()) d.interface_traps = true;
        if (f.fixed_charge_cm2 != 0.0) d.interface_fixed_charge = true;
    }
    d.thermal_contacts = device.thermal_contacts().size();
    return d;
}

}  // namespace NiTCAD::verify::benchmark
