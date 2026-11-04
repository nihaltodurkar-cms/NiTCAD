// Copies a device into the consistency layer's plain data. See device_data.hpp.
#include "NiTCAD/verify/benchmark/device_data.hpp"

#include <cstdint>
#include <variant>

#include "NiTCAD/base/constants.hpp"

namespace NiTCAD::verify::benchmark {

namespace {

std::int32_t region_id(const device::Device& device, const std::string& name) {
    const auto regions = device.regions();
    for (std::size_t r = 0; r < regions.size(); ++r) {
        if (regions[r].name == name) return static_cast<std::int32_t>(r);
    }
    return -1;
}

}  // namespace

consistency::DeviceData device_data(const device::Device& device) {
    consistency::DeviceData d;
    d.temperature_K = device.temperature_K();
    d.donors.assign(device.donors().begin(), device.donors().end());
    d.acceptors.assign(device.acceptors().begin(), device.acceptors().end());
    d.region.assign(device.node_region().begin(), device.node_region().end());
    const mesh::Mesh& m = device.mesh();
    const std::size_t n = m.node_count();
    d.relative_permittivity.resize(n);
    d.insulator.resize(n);
    d.thermopower_n.assign(n, 0.0);
    d.thermopower_p.assign(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const auto node = static_cast<mesh::NodeId>(i);
        d.relative_permittivity[i] = device.relative_permittivity(node);
        d.insulator[i] = device.is_insulator(node) ? 1 : 0;
        const auto& material = device.regions()[static_cast<std::size_t>(device.node_region()[i])].material;
        if (const auto* s = std::get_if<physics::Semiconductor>(&material)) {
            d.thermopower_n[i] = s->parameters().thermal.thermopower_exponent_n;
            d.thermopower_p[i] = s->parameters().thermal.thermopower_exponent_p;
        }
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
        if (c.kind == device::ContactKind::gate) {
            // The face area of each node on the gate's patch (both lists increase), as the gate term reads it.
            const mesh::BoundaryPatch* patch = m.find_boundary(c.gate.boundary);
            if (patch != nullptr) {
                std::size_t j = 0;
                for (const mesh::NodeId v : c.nodes) {
                    while (j < patch->nodes.size() && patch->nodes[j] < v) ++j;
                    k.gate_areas.push_back(j < patch->nodes.size() && patch->nodes[j] == v ? patch->areas[j] : 0.0);
                }
            }
            k.oxide_capacitance_F_cm2 =
                base::eps0_F_per_cm * c.gate.oxide_relative_permittivity / c.gate.oxide_thickness_cm;
            k.gate_fixed_charge_cm2 = c.gate.fixed_charge_cm2;
        }
        d.contacts.push_back(std::move(k));
    }
    for (const auto& f : device.interfaces()) {
        if (!f.traps.levels.empty() || !f.traps.bands.empty()) d.interface_traps = true;
        if (f.fixed_charge_cm2 != 0.0) d.interface_fixed_charge = true;
        consistency::InterfaceData g;
        g.region_a = region_id(device, f.region_a);
        g.region_b = region_id(device, f.region_b);
        g.charged = device::has_interface_charge_or_recombination(f);
        g.fixed_charge_cm2 = f.fixed_charge_cm2;
        for (const auto& t : f.traps.levels) g.trap_density_cm2 += t.density_cm2;
        for (const auto& b : f.traps.bands) g.trap_density_cm2 += b.density_cm2_eV * (b.energy_high_eV - b.energy_low_eV);
        d.interfaces.push_back(g);
    }
    d.thermal_contacts = device.thermal_contacts().size();
    return d;
}

consistency::DeviceData device_data(const device::Device& device, const assemble::PhysicsModels& models) {
    consistency::DeviceData d = device_data(device);
    d.incomplete_ionization = models.incomplete_ionization;
    d.electrothermal = models.electrothermal;
    return d;
}

}  // namespace NiTCAD::verify::benchmark
