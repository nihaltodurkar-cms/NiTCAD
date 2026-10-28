#include "NiTCAD/io/device_file.hpp"

#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "NiTCAD/mesh/tensor_cells.hpp"
#include "container.hpp"
#include "record.hpp"

namespace NiTCAD::io {

namespace {

using detail::corrupt;
using detail::Reader;
using detail::Writer;

namespace device_key {
constexpr std::uint16_t mesh = 1, temperature = 2, region = 3, node_region = 4, donors = 5,
                        acceptors = 6, contact = 7, interface = 8, cells = 9, thermal_contact = 10,
                        certificate = 11;
}
namespace region_key {
constexpr std::uint16_t name = 1, kind = 2, parameters = 3;  // kind 0 semiconductor, 1 insulator
}
namespace contact_key {
constexpr std::uint16_t name = 1, kind = 2, nodes = 3, gate = 4, electrode = 5, schottky = 6;
constexpr std::uint16_t gate_boundary = 1, gate_thickness = 2, gate_permittivity = 3,
                        gate_electrode = 4, gate_work_function = 5, gate_fixed_charge = 6;
constexpr std::uint16_t electrode_kind = 1, electrode_work_function = 2;
constexpr std::uint16_t schottky_boundary = 1, schottky_by_work_function = 2,
                        schottky_work_function = 3, schottky_barrier = 4, schottky_richardson_n = 5,
                        schottky_richardson_p = 6;
}  // namespace contact_key
namespace interface_key {
constexpr std::uint16_t region_a = 1, region_b = 2, transport = 3, fixed_charge = 4, traps = 5,
                        velocity_n = 6, velocity_p = 7;
constexpr std::uint16_t level = 1, band = 2, thermal_velocity_n = 3, thermal_velocity_p = 4;
constexpr std::uint16_t trap_type = 1, trap_density = 2, trap_energy = 3, trap_energy_high = 4,
                        trap_sigma_n = 5, trap_sigma_p = 6;
}  // namespace interface_key
namespace cells_key {
constexpr std::uint16_t x = 1, y = 2, z = 3;
}
namespace thermal_key {
constexpr std::uint16_t name = 1, boundary = 2, kind = 3, temperature = 4, resistance = 5;
}

// Every thermal parameter with its key, for a writer or a reader (`f(key, member)`).
template <class T, class F>
void thermal_fields(T& t, F&& f) {
    f(1, t.conductivity_W_cmK);
    f(2, t.conductivity_exponent);
    f(3, t.heat_capacity_J_cm3K);
    f(4, t.thermopower_exponent_n);
    f(5, t.thermopower_exponent_p);
}

// Every semiconductor parameter with its key. The keys are the format: never renumber one; a new
// parameter takes a new key (thermal data are a nested record, key 100).
template <class P, class F>
void semiconductor_fields(P& p, F&& f) {
    f(1, p.eps_r);
    f(2, p.Eg0_eV);
    f(3, p.varshni_alpha_eV_per_K);
    f(4, p.varshni_beta_K);
    f(5, p.Nc300);
    f(6, p.Nv300);
    f(7, p.electron_affinity_eV);
    const auto mobility = [&f](int base, auto& m) {
        f(base, m.mu_min);
        f(base + 1, m.mu_max);
        f(base + 2, m.N_ref);
        f(base + 3, m.alpha);
        f(base + 4, m.T_exponent);
    };
    mobility(10, p.electron_mobility);
    mobility(15, p.hole_mobility);
    f(20, p.lifetime.tau_n0);
    f(21, p.lifetime.tau_p0);
    f(22, p.lifetime.N_ref);
    f(23, p.auger.Cn);
    f(24, p.auger.Cp);
    f(25, p.bandgap_narrowing.E0_eV);
    f(26, p.bandgap_narrowing.N0);
    f(27, p.electron_saturation.v_sat_cm_s);
    f(28, p.electron_saturation.beta);
    f(29, p.hole_saturation.v_sat_cm_s);
    f(30, p.hole_saturation.beta);
    f(31, p.radiative_cm3_s);
    f(32, p.ionization.donor_eV);
    f(33, p.ionization.acceptor_eV);
    f(34, p.ionization.donor_degeneracy);
    f(35, p.ionization.acceptor_degeneracy);
    f(36, p.richardson.electron);
    f(37, p.richardson.hole);
    const auto impact = [&f](int base, auto& c) {
        f(base, c.A_low_per_cm);
        f(base + 1, c.B_low_V_per_cm);
        f(base + 2, c.A_high_per_cm);
        f(base + 3, c.B_high_V_per_cm);
        f(base + 4, c.switch_V_per_cm);
    };
    impact(40, p.impact_ionization.electron);
    impact(45, p.impact_ionization.hole);
    f(50, p.impact_ionization.phonon_energy_eV);
    f(51, p.band_to_band.A_per_cm3_s);
    f(52, p.band_to_band.B_V_per_cm);
    f(53, p.band_to_band.direct_gap);
    f(54, p.band_to_band.electron_mass);
    f(55, p.band_to_band.hole_mass);
    // Unit 35: the Passler gap and the density-of-states mass fits (absent: legacy forms).
    f(56, p.passler_gap.alpha_eV_per_K);
    f(57, p.passler_gap.Theta_K);
    f(58, p.passler_gap.Delta);
    for (int k = 0; k < 4; ++k) {
        f(60 + k, p.dos_mass.conduction[k]);
        f(64 + k, p.dos_mass.valence[k]);
    }
}

constexpr std::uint16_t thermal_record_key = 100;

// A writer visitor: doubles by bit pattern, bools as 0 / 1.
struct Put {
    Writer& w;
    void operator()(int key, double v) const { w.f64(static_cast<std::uint16_t>(key), v); }
    void operator()(int key, bool v) const { w.u64(static_cast<std::uint16_t>(key), v ? 1 : 0); }
};

// A reader visitor: a present key overwrites the member, an absent one leaves its default (IO2).
struct Get {
    const Reader& r;
    std::optional<base::Error>& error;
    void operator()(int k, double& v) const {
        const auto key = static_cast<std::uint16_t>(k);
        if (error || !r.has(key)) return;
        if (auto x = r.f64(key)) v = *x; else error = std::move(x.error());
    }
    void operator()(int k, bool& v) const {
        const auto key = static_cast<std::uint16_t>(k);
        if (error || !r.has(key)) return;
        if (auto x = r.u64(key)) v = *x != 0; else error = std::move(x.error());
    }
};

Writer write_thermal(const physics::ThermalParameters& t) {
    Writer w;
    thermal_fields(t, Put{w});
    return w;
}

std::expected<physics::ThermalParameters, base::Error> read_thermal(const Reader& owner) {
    physics::ThermalParameters t{};
    if (!owner.has(thermal_record_key)) return t;
    NITCAD_IO_TRY(r, owner.record(thermal_record_key));
    std::optional<base::Error> error;
    thermal_fields(t, Get{r, error});
    if (error) return std::unexpected(std::move(*error));
    return t;
}

Writer write_region(const device::Region& region) {
    Writer w;
    w.text(region_key::name, region.name);
    Writer p;
    if (const auto* s = std::get_if<physics::Semiconductor>(&region.material)) {
        w.u64(region_key::kind, 0);
        semiconductor_fields(s->parameters(), Put{p});
        p.record(thermal_record_key, write_thermal(s->parameters().thermal));
    } else {
        const auto& i = std::get<physics::Insulator>(region.material).parameters();
        w.u64(region_key::kind, 1);
        p.f64(1, i.eps_r);
        p.record(thermal_record_key, write_thermal(i.thermal));
    }
    w.record(region_key::parameters, p);
    return w;
}

std::expected<device::Region, base::Error> read_region(const Reader& r) {
    NITCAD_IO_TRY(name, r.text(region_key::name));
    NITCAD_IO_TRY(kind, r.u64(region_key::kind));
    NITCAD_IO_TRY(p, r.record(region_key::parameters));
    NITCAD_IO_TRY(thermal, read_thermal(p));
    std::optional<base::Error> error;
    if (kind == 0) {
        physics::SemiconductorParameters s{};
        semiconductor_fields(s, Get{p, error});
        if (error) return std::unexpected(std::move(*error));
        s.thermal = thermal;
        NITCAD_IO_TRY(material, physics::Semiconductor::create(s));
        return device::Region{std::move(name), std::move(material)};
    }
    if (kind == 1) {
        physics::InsulatorParameters i{};
        NITCAD_IO_TRY(eps_r, p.f64(1));
        i.eps_r = eps_r;
        i.thermal = thermal;
        NITCAD_IO_TRY(material, physics::Insulator::create(i));
        return device::Region{std::move(name), std::move(material)};
    }
    return std::unexpected(corrupt("unknown region kind"));
}

Writer write_contact(const device::Contact& c) {
    Writer w;
    w.text(contact_key::name, c.name);
    w.u64(contact_key::kind, static_cast<std::uint64_t>(c.kind));
    w.i64s(contact_key::nodes, std::vector<std::int64_t>(c.nodes.begin(), c.nodes.end()));
    if (c.kind == device::ContactKind::gate) {
        Writer g;
        g.text(contact_key::gate_boundary, c.gate.boundary);
        g.f64(contact_key::gate_thickness, c.gate.oxide_thickness_cm);
        g.f64(contact_key::gate_permittivity, c.gate.oxide_relative_permittivity);
        g.u64(contact_key::gate_electrode, static_cast<std::uint64_t>(c.gate.electrode));
        g.f64(contact_key::gate_work_function, c.gate.work_function_eV);
        g.f64(contact_key::gate_fixed_charge, c.gate.fixed_charge_cm2);
        w.record(contact_key::gate, g);
    }
    if (c.kind == device::ContactKind::electrode) {
        Writer e;
        e.u64(contact_key::electrode_kind, static_cast<std::uint64_t>(c.electrode.kind));
        e.f64(contact_key::electrode_work_function, c.electrode.work_function_eV);
        w.record(contact_key::electrode, e);
    }
    if (c.kind == device::ContactKind::schottky) {  // Unit 48: a new optional record
        Writer s;
        s.text(contact_key::schottky_boundary, c.schottky.boundary);
        s.u64(contact_key::schottky_by_work_function, c.schottky.by_work_function ? 1 : 0);
        s.f64(contact_key::schottky_work_function, c.schottky.work_function_eV);
        s.f64(contact_key::schottky_barrier, c.schottky.barrier_eV);
        s.f64(contact_key::schottky_richardson_n, c.schottky.richardson_n);
        s.f64(contact_key::schottky_richardson_p, c.schottky.richardson_p);
        w.record(contact_key::schottky, s);
    }
    return w;
}

// A small enumeration value, or corrupt (Device::create names an out-of-range kind; this keeps a
// value that does not fit the enumeration's type from wrapping into a valid one).
template <class E>
std::expected<E, base::Error> enumerator(const Reader& r, std::uint16_t key, const char* what) {
    NITCAD_IO_TRY(v, r.u64(key));
    if (v > std::numeric_limits<std::underlying_type_t<E>>::max()) {
        return std::unexpected(corrupt(std::string("unknown ") + what));
    }
    return static_cast<E>(v);
}

std::expected<device::Contact, base::Error> read_contact(const Reader& r) {
    device::Contact c;
    NITCAD_IO_TRY(name, r.text(contact_key::name));
    NITCAD_IO_TRY(kind, enumerator<device::ContactKind>(r, contact_key::kind, "contact kind"));
    NITCAD_IO_TRY(nodes, r.i64s(contact_key::nodes));
    c.name = std::move(name);
    c.kind = kind;
    for (const std::int64_t v : nodes) {
        if (v < 0 || v > std::numeric_limits<mesh::NodeId>::max()) {
            return std::unexpected(corrupt("contact node out of range"));
        }
        c.nodes.push_back(static_cast<mesh::NodeId>(v));
    }
    if (r.has(contact_key::gate)) {
        NITCAD_IO_TRY(g, r.record(contact_key::gate));
        NITCAD_IO_TRY(boundary, g.text(contact_key::gate_boundary));
        NITCAD_IO_TRY(thickness, g.f64(contact_key::gate_thickness));
        NITCAD_IO_TRY(permittivity, g.f64(contact_key::gate_permittivity));
        NITCAD_IO_TRY(electrode, enumerator<device::GateElectrode>(g, contact_key::gate_electrode,
                                                                   "gate electrode"));
        NITCAD_IO_TRY(work_function, g.f64(contact_key::gate_work_function));
        NITCAD_IO_TRY(fixed_charge, g.f64(contact_key::gate_fixed_charge));
        c.gate = {std::move(boundary), thickness, permittivity, electrode, work_function,
                  fixed_charge};
    }
    if (r.has(contact_key::electrode)) {
        NITCAD_IO_TRY(e, r.record(contact_key::electrode));
        NITCAD_IO_TRY(ekind, enumerator<device::GateElectrode>(e, contact_key::electrode_kind,
                                                               "electrode kind"));
        NITCAD_IO_TRY(work_function, e.f64(contact_key::electrode_work_function));
        c.electrode = {ekind, work_function};
    }
    if (r.has(contact_key::schottky)) {
        NITCAD_IO_TRY(s, r.record(contact_key::schottky));
        NITCAD_IO_TRY(boundary, s.text(contact_key::schottky_boundary));
        NITCAD_IO_TRY(by_work_function, s.u64(contact_key::schottky_by_work_function));
        if (by_work_function > 1) return std::unexpected(corrupt("unknown Schottky barrier kind"));
        NITCAD_IO_TRY(work_function, s.f64(contact_key::schottky_work_function));
        NITCAD_IO_TRY(barrier, s.f64(contact_key::schottky_barrier));
        NITCAD_IO_TRY(richardson_n, s.f64(contact_key::schottky_richardson_n));
        NITCAD_IO_TRY(richardson_p, s.f64(contact_key::schottky_richardson_p));
        c.schottky = {std::move(boundary), by_work_function == 1, work_function, barrier,
                      richardson_n, richardson_p};
    }
    return c;
}

std::expected<physics::TrapType, base::Error> trap_type(const Reader& r) {
    NITCAD_IO_TRY(v, r.u64(interface_key::trap_type));
    if (v > static_cast<std::uint64_t>(physics::TrapType::acceptor)) {
        return std::unexpected(corrupt("unknown trap type"));
    }
    return static_cast<physics::TrapType>(v);
}

Writer write_interface(const device::Interface& f) {
    Writer w;
    w.text(interface_key::region_a, f.region_a);
    w.text(interface_key::region_b, f.region_b);
    w.u64(interface_key::transport, static_cast<std::uint64_t>(f.transport));
    w.f64(interface_key::fixed_charge, f.fixed_charge_cm2);
    Writer t;
    for (const physics::TrapLevel& l : f.traps.levels) {
        Writer lw;
        lw.u64(interface_key::trap_type, static_cast<std::uint64_t>(l.type));
        lw.f64(interface_key::trap_density, l.density_cm2);
        lw.f64(interface_key::trap_energy, l.energy_eV);
        lw.f64(interface_key::trap_sigma_n, l.sigma_n_cm2);
        lw.f64(interface_key::trap_sigma_p, l.sigma_p_cm2);
        t.record(interface_key::level, lw);
    }
    for (const physics::TrapBand& b : f.traps.bands) {
        Writer bw;
        bw.u64(interface_key::trap_type, static_cast<std::uint64_t>(b.type));
        bw.f64(interface_key::trap_density, b.density_cm2_eV);
        bw.f64(interface_key::trap_energy, b.energy_low_eV);
        bw.f64(interface_key::trap_energy_high, b.energy_high_eV);
        bw.f64(interface_key::trap_sigma_n, b.sigma_n_cm2);
        bw.f64(interface_key::trap_sigma_p, b.sigma_p_cm2);
        t.record(interface_key::band, bw);
    }
    t.f64(interface_key::thermal_velocity_n, f.traps.thermal_velocity_n_cm_s);
    t.f64(interface_key::thermal_velocity_p, f.traps.thermal_velocity_p_cm_s);
    w.record(interface_key::traps, t);
    w.f64(interface_key::velocity_n, f.recombination_velocity_n_cm_s);
    w.f64(interface_key::velocity_p, f.recombination_velocity_p_cm_s);
    return w;
}

std::expected<device::Interface, base::Error> read_interface(const Reader& r) {
    device::Interface f;
    NITCAD_IO_TRY(a, r.text(interface_key::region_a));
    NITCAD_IO_TRY(b, r.text(interface_key::region_b));
    NITCAD_IO_TRY(transport, enumerator<device::InterfaceTransport>(r, interface_key::transport,
                                                                    "interface transport"));
    NITCAD_IO_TRY(fixed_charge, r.f64(interface_key::fixed_charge));
    NITCAD_IO_TRY(velocity_n, r.f64(interface_key::velocity_n));
    NITCAD_IO_TRY(velocity_p, r.f64(interface_key::velocity_p));
    f.region_a = std::move(a);
    f.region_b = std::move(b);
    f.transport = transport;
    f.fixed_charge_cm2 = fixed_charge;
    f.recombination_velocity_n_cm_s = velocity_n;
    f.recombination_velocity_p_cm_s = velocity_p;
    NITCAD_IO_TRY(t, r.record(interface_key::traps));
    NITCAD_IO_TRY(levels, t.records(interface_key::level));
    for (const Reader& l : levels) {
        NITCAD_IO_TRY(type, trap_type(l));
        NITCAD_IO_TRY(density, l.f64(interface_key::trap_density));
        NITCAD_IO_TRY(energy, l.f64(interface_key::trap_energy));
        NITCAD_IO_TRY(sigma_n, l.f64(interface_key::trap_sigma_n));
        NITCAD_IO_TRY(sigma_p, l.f64(interface_key::trap_sigma_p));
        f.traps.levels.push_back({type, density, energy, sigma_n, sigma_p});
    }
    NITCAD_IO_TRY(bands, t.records(interface_key::band));
    for (const Reader& bd : bands) {
        NITCAD_IO_TRY(type, trap_type(bd));
        NITCAD_IO_TRY(density, bd.f64(interface_key::trap_density));
        NITCAD_IO_TRY(low, bd.f64(interface_key::trap_energy));
        NITCAD_IO_TRY(high, bd.f64(interface_key::trap_energy_high));
        NITCAD_IO_TRY(sigma_n, bd.f64(interface_key::trap_sigma_n));
        NITCAD_IO_TRY(sigma_p, bd.f64(interface_key::trap_sigma_p));
        f.traps.bands.push_back({type, density, low, high, sigma_n, sigma_p});
    }
    NITCAD_IO_TRY(vn, t.f64(interface_key::thermal_velocity_n));
    NITCAD_IO_TRY(vp, t.f64(interface_key::thermal_velocity_p));
    f.traps.thermal_velocity_n_cm_s = vn;
    f.traps.thermal_velocity_p_cm_s = vp;
    return f;
}

Writer write_thermal_contact(const device::ThermalContact& t) {
    Writer w;
    w.text(thermal_key::name, t.name);
    w.text(thermal_key::boundary, t.boundary);
    w.u64(thermal_key::kind, static_cast<std::uint64_t>(t.kind));
    w.f64(thermal_key::temperature, t.temperature_K);
    w.f64(thermal_key::resistance, t.resistance_K_cm2_W);
    return w;
}

std::expected<device::ThermalContact, base::Error> read_thermal_contact(const Reader& r) {
    NITCAD_IO_TRY(name, r.text(thermal_key::name));
    NITCAD_IO_TRY(boundary, r.text(thermal_key::boundary));
    NITCAD_IO_TRY(kind, enumerator<device::ThermalContactKind>(r, thermal_key::kind,
                                                               "thermal contact kind"));
    NITCAD_IO_TRY(temperature, r.f64(thermal_key::temperature));
    NITCAD_IO_TRY(resistance, r.f64(thermal_key::resistance));
    return device::ThermalContact{std::move(name), std::move(boundary), kind, temperature,
                                  resistance};
}

}  // namespace

std::expected<void, base::Error> save_device(const std::filesystem::path& path,
                                             const device::Device& d) {
    Writer w;
    w.record(device_key::mesh, detail::write_mesh(d.mesh()));
    w.f64(device_key::temperature, d.temperature_K());
    for (const device::Region& r : d.regions()) w.record(device_key::region, write_region(r));
    w.i64s(device_key::node_region,
           std::vector<std::int64_t>(d.node_region().begin(), d.node_region().end()));
    w.f64s(device_key::donors, d.donors());
    w.f64s(device_key::acceptors, d.acceptors());
    for (const device::Contact& c : d.contacts()) w.record(device_key::contact, write_contact(c));
    for (const device::Interface& f : d.interfaces()) {
        w.record(device_key::interface, write_interface(f));
    }
    if (const mesh::TensorCells* cells = d.cells()) {
        Writer c;
        constexpr std::uint16_t axis_key[3] = {cells_key::x, cells_key::y, cells_key::z};
        for (int a = 0; a < cells->dimension(); ++a) c.f64s(axis_key[a], cells->axis(a));
        w.record(device_key::cells, c);
    }
    for (const device::ThermalContact& t : d.thermal_contacts()) {
        w.record(device_key::thermal_contact, write_thermal_contact(t));
    }
    if (const mesh::InterfaceCertificate* c = d.interface_certificate()) {
        w.record(device_key::certificate, detail::write_certificate(*c));
    }
    return detail::write_file(path, Content::device, w);
}

std::expected<device::Device, base::Error> load_device(const std::filesystem::path& path) {
    NITCAD_IO_TRY(bytes, detail::read_file(path, Content::device));
    NITCAD_IO_TRY(r, Reader::parse(bytes));
    NITCAD_IO_TRY(mesh_record, r.record(device_key::mesh));
    NITCAD_IO_TRY(m, detail::read_mesh(mesh_record));
    NITCAD_IO_TRY(temperature, r.f64(device_key::temperature));
    device::DeviceDescription d{.mesh = std::move(m), .temperature_K = temperature};
    NITCAD_IO_TRY(regions, r.records(device_key::region));
    for (const Reader& region : regions) {
        NITCAD_IO_TRY(g, read_region(region));
        d.regions.push_back(std::move(g));
    }
    NITCAD_IO_TRY(node_region, r.i64s(device_key::node_region));
    for (const std::int64_t v : node_region) {
        if (v < 0 || v > std::numeric_limits<device::RegionId>::max()) {
            return std::unexpected(corrupt("node region out of range"));
        }
        d.node_region.push_back(static_cast<device::RegionId>(v));
    }
    NITCAD_IO_TRY(donors, r.f64s(device_key::donors));
    NITCAD_IO_TRY(acceptors, r.f64s(device_key::acceptors));
    d.donors = std::move(donors);
    d.acceptors = std::move(acceptors);
    NITCAD_IO_TRY(contacts, r.records(device_key::contact));
    for (const Reader& c : contacts) {
        NITCAD_IO_TRY(contact, read_contact(c));
        d.contacts.push_back(std::move(contact));
    }
    NITCAD_IO_TRY(interfaces, r.records(device_key::interface));
    for (const Reader& f : interfaces) {
        NITCAD_IO_TRY(interface, read_interface(f));
        d.interfaces.push_back(std::move(interface));
    }
    if (r.has(device_key::cells)) {
        NITCAD_IO_TRY(c, r.record(device_key::cells));
        NITCAD_IO_TRY(x, c.f64s_or_empty(cells_key::x));
        NITCAD_IO_TRY(y, c.f64s_or_empty(cells_key::y));
        NITCAD_IO_TRY(z, c.f64s_or_empty(cells_key::z));
        NITCAD_IO_TRY(cells, mesh::TensorCells::create(d.mesh, x, y, z));
        d.cells = std::move(cells);
    }
    NITCAD_IO_TRY(thermal, r.records(device_key::thermal_contact));
    for (const Reader& t : thermal) {
        NITCAD_IO_TRY(contact, read_thermal_contact(t));
        d.thermal_contacts.push_back(std::move(contact));
    }
    if (r.has(device_key::certificate)) {
        NITCAD_IO_TRY(c, r.record(device_key::certificate));
        NITCAD_IO_TRY(certificate, detail::read_certificate(c));
        d.interface_certificate = certificate;
    }
    return device::Device::create(std::move(d));
}

}  // namespace NiTCAD::io
