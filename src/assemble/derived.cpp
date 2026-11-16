// The derived-field table (Unit 57 slice 3, DECISIONS.md GE4): each evaluator calls its topic
// function and moves the result into one vector per field, so the values are the function's.
#include "NiTCAD/assemble/derived.hpp"

#include <array>
#include <utility>

#include "NiTCAD/assemble/band_edges.hpp"
#include "NiTCAD/assemble/charges.hpp"
#include "NiTCAD/assemble/currents.hpp"
#include "NiTCAD/assemble/energy.hpp"
#include "NiTCAD/assemble/rates.hpp"
#include "NiTCAD/assemble/system.hpp"
#include "NiTCAD/assemble/thermal.hpp"
#include "NiTCAD/base/contract.hpp"
#include "system_data.hpp"

namespace NiTCAD::assemble {

namespace {

constexpr Field derived(FieldId id, std::string_view name, std::string_view units,
                        FieldLocation location, FieldDomain domain) {
    return {id, name, units, location, FieldKind::derived, domain};
}

constexpr auto node = FieldLocation::node;
constexpr auto semi = FieldDomain::semiconductor;
constexpr auto all = FieldDomain::all;

constexpr std::array band_fields{
    derived(FieldId::conduction_band, "conduction band edge", "V_T", node, semi),
    derived(FieldId::valence_band, "valence band edge", "V_T", node, semi),
    derived(FieldId::electron_fermi, "electron quasi-Fermi level", "V_T", node, semi),
    derived(FieldId::hole_fermi, "hole quasi-Fermi level", "V_T", node, semi),
};

constexpr std::string_view rate_units = "cm^-3 s^-1";
constexpr std::array rate_fields{
    derived(FieldId::srh_rate, "SRH recombination", rate_units, node, semi),
    derived(FieldId::auger_rate, "Auger recombination", rate_units, node, semi),
    derived(FieldId::radiative_rate, "radiative recombination", rate_units, node, semi),
    derived(FieldId::impact_rate, "impact ionization", rate_units, node, semi),
    derived(FieldId::interface_rate_n, "interface rate, electrons", rate_units, node, semi),
    derived(FieldId::interface_rate_p, "interface rate, holes", rate_units, node, semi),
    derived(FieldId::band_to_band_rate_n, "band-to-band rate, electrons", rate_units, node, semi),
    derived(FieldId::band_to_band_rate_p, "band-to-band rate, holes", rate_units, node, semi),
    derived(FieldId::net_rate_n, "net rate, electrons", rate_units, node, semi),
    derived(FieldId::net_rate_p, "net rate, holes", rate_units, node, semi),
};

constexpr std::array current_fields{
    derived(FieldId::electron_edge_current, "electron edge current", "J0", FieldLocation::edge,
            all),
    derived(FieldId::hole_edge_current, "hole edge current", "J0", FieldLocation::edge, all),
};

constexpr std::array flux_fields{
    derived(FieldId::edge_flux_terms, "edge flux terms", "J0", FieldLocation::edge, all),
};

constexpr std::array trap_fields{
    derived(FieldId::trap_occupancy, "trap occupancy", "1", FieldLocation::interface_edge, all),
};

constexpr std::array temperature_fields{
    derived(FieldId::lattice_temperature, "lattice temperature", "K", node, all),
    derived(FieldId::electron_temperature_K, "electron temperature", "K", node, semi),
    derived(FieldId::hole_temperature_K, "hole temperature", "K", node, semi),
};

DerivedValues bands(const System& system, std::span<const double> x) {
    BandEdges b = band_edges(system, x);
    DerivedValues v;
    v.reserve(band_fields.size());
    v.push_back(std::move(b.conduction));
    v.push_back(std::move(b.valence));
    v.push_back(std::move(b.electron_fermi));
    v.push_back(std::move(b.hole_fermi));
    return v;
}

// As each mechanism's rows see it, and the sums in rates()'s order (net_n, net_p).
DerivedValues mechanism_rates(const System& system, std::span<const double> x) {
    using M = RateMechanism;
    const std::size_t nodes = system.node_count();
    DerivedValues v(rate_fields.size());
    std::vector<double>& net_n = v[8];
    std::vector<double>& net_p = v[9];
    net_n.assign(nodes, 0.0);
    net_p.assign(nodes, 0.0);
    for (MechanismRate& m : rates(system, x)) {
        for (std::size_t i = 0; i < nodes; ++i) {
            net_n[i] += m.electron[i];
            net_p[i] += m.hole[i];
        }
        switch (m.mechanism) {
            case M::srh: v[0] = std::move(m.electron); break;
            case M::auger: v[1] = std::move(m.electron); break;
            case M::radiative: v[2] = std::move(m.electron); break;
            case M::impact: v[3] = std::move(m.electron); break;
            case M::interfaces:
                v[4] = std::move(m.electron);
                v[5] = std::move(m.hole);
                break;
            case M::band_to_band:
                v[6] = std::move(m.electron);
                v[7] = std::move(m.hole);
                break;
        }
    }
    return v;
}

DerivedValues currents(const System& system, std::span<const double> x) {
    DerivedValues v(current_fields.size());
    for (const auto& [jn, jp] : edge_currents(system, x)) {
        v[0].push_back(jn);
        v[1].push_back(jp);
    }
    return v;
}

DerivedValues flux_terms(const System& system, std::span<const double> x) {
    DerivedValues v;
    v.push_back(edge_flux_terms(system, x));
    return v;
}

DerivedValues occupancies(const System& system, std::span<const double> x) {
    DerivedValues v;
    v.push_back(trap_occupancies(system, x));
    return v;
}

// T = T0 + T0 tau per node (T0 the scaling temperature; for the carriers, the lattice's, which
// energy balance holds at T0, Unit 52).
DerivedValues temperatures(const System& system, std::span<const double> x) {
    DerivedValues v(temperature_fields.size());
    const double T0 = SystemData::of(system).temperature_K_;
    const DofMap& dofs = system.dofs();
    const std::size_t block = dofs.block();
    const auto fill = [&](std::vector<double>& out, FieldId rise) {
        const std::size_t k = *dofs.slot(rise);
        out.reserve(system.node_count());
        for (std::size_t i = 0; i < system.node_count(); ++i) {
            out.push_back(T0 + T0 * x[block * i + k]);
        }
    };
    if (electrothermal(system)) fill(v[0], FieldId::temperature);
    if (energy_balance(system)) {
        fill(v[1], FieldId::electron_temperature);
        fill(v[2], FieldId::hole_temperature);
    }
    return v;
}

constexpr std::array table{
    Diagnostic{DiagnosticId::bands, "bands", band_fields, &bands},
    Diagnostic{DiagnosticId::rates, "rates", rate_fields, &mechanism_rates},
    Diagnostic{DiagnosticId::edge_currents, "edge currents", current_fields, &currents},
    Diagnostic{DiagnosticId::edge_flux_terms, "edge flux terms", flux_fields, &flux_terms},
    Diagnostic{DiagnosticId::trap_occupancies, "trap occupancies", trap_fields, &occupancies},
    Diagnostic{DiagnosticId::temperatures, "temperatures", temperature_fields, &temperatures},
};

}  // namespace

std::span<const Diagnostic> diagnostics() noexcept { return table; }

const Diagnostic& diagnostic(DiagnosticId id) noexcept {
    return table[static_cast<std::size_t>(id)];
}

std::size_t position(const Diagnostic& diagnostic, FieldId field) {
    for (std::size_t k = 0; k < diagnostic.fields.size(); ++k) {
        if (diagnostic.fields[k].id == field) return k;
    }
    NITCAD_EXPECTS(false && "the diagnostic has no such field");
    return diagnostic.fields.size();
}

std::vector<double> take(DerivedValues& values, const Diagnostic& diagnostic, FieldId field) {
    NITCAD_EXPECTS(values.size() == diagnostic.fields.size());
    return std::move(values[position(diagnostic, field)]);
}

}  // namespace NiTCAD::assemble
