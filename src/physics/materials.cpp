#include "NiTCAD/physics/materials.hpp"

#include <array>

namespace NiTCAD::physics {
namespace {

template <const SemiconductorParameters& P>
std::expected<MaterialParameters, base::Error> fixed_semiconductor(std::span<const double>) {
    return P;
}

template <const InsulatorParameters& P>
std::expected<MaterialParameters, base::Error> fixed_insulator(std::span<const double>) {
    return P;
}

// `insulator`: the generic insulator, its eps_r given in the region (0 until then, refused by
// Insulator::create if left unset).
constexpr InsulatorParameters generic_insulator{.eps_r = 0.0};

std::expected<MaterialParameters, base::Error> algaas(std::span<const double> a) {
    auto p = a.size() == 2 ? algaas_parameters(a[0], a[1]) : algaas_parameters(a[0]);
    if (!p) return std::unexpected(std::move(p.error()));
    return *p;
}

// The sources named in each set's header comment (MAD3: no new reading for existing sets).
constexpr std::array<Citation, 2> silicon_sources{{
    {"R. Couderc, M. Amara, M. Lemiti, J. Appl. Phys. 115, 093705 (2014)", "band set (Unit 35)"},
    {"legacy materials.py SILICON", "every other field"},
}};
constexpr std::array<Citation, 1> legacy_silicon_sources{{
    {"legacy materials.py SILICON", "every field (Unit 15)"},
}};
constexpr std::array<Citation, 2> legacy_set_sources{{
    {"legacy materials.py (Adachi handbook band values)", "band, mobility, Auger, lifetimes (Unit 15)"},
    {"Ioffe NSM archive", "radiative coefficient"},
}};
constexpr std::array<Citation, 2> sic_sources{{
    {"legacy materials.py (Baliga)", "eps_r, chi, Eg(300 K), mobility end points"},
    {"Ioffe NSM archive", "radiative coefficient, heat capacity"},
}};
constexpr std::array<Citation, 1> algaas_sources{{
    {"legacy materials.algaas", "linear in x, conduction share of the gap step"},
}};
constexpr std::array<Citation, 1> sio2_sources{{
    {"legacy moscap.EPS_OX_R; NIST-JANAF", "eps_r 3.9; heat capacity (Unit 23)"},
}};

constexpr std::array<MaterialEntry, 9> registry{{
    {"silicon", MaterialKind::semiconductor, 0, 0, "", "silicon, Couderc band set (default)", silicon_sources,
     fixed_semiconductor<silicon_parameters>},
    {"silicon_legacy", MaterialKind::semiconductor, 0, 0, "", "silicon, the legacy set", legacy_silicon_sources,
     fixed_semiconductor<legacy_silicon_parameters>},
    {"germanium", MaterialKind::semiconductor, 0, 0, "", "germanium", legacy_set_sources,
     fixed_semiconductor<germanium_parameters>},
    {"gaas", MaterialKind::semiconductor, 0, 0, "", "GaAs, Gamma-valley gap", legacy_set_sources,
     fixed_semiconductor<gallium_arsenide_parameters>},
    {"ingaas", MaterialKind::semiconductor, 0, 0, "", "In0.53Ga0.47As, lattice-matched to InP", legacy_set_sources,
     fixed_semiconductor<indium_gallium_arsenide_parameters>},
    {"sic_4h", MaterialKind::semiconductor, 0, 0, "", "4H-SiC", sic_sources,
     fixed_semiconductor<silicon_carbide_4h_parameters>},
    {"algaas", MaterialKind::semiconductor, 1, 2, "algaas takes the Al fraction and optionally the conduction share",
     "Al(x)Ga(1-x)As, 0 <= x <= 0.45", algaas_sources, algaas},
    {"sio2", MaterialKind::insulator, 0, 0, "", "SiO2", sio2_sources, fixed_insulator<silicon_dioxide_parameters>},
    {"insulator", MaterialKind::insulator, 0, 0, "", "generic insulator, eps_r given in the region", {},
     fixed_insulator<generic_insulator>},
}};

}  // namespace

std::span<const MaterialEntry> material_registry() noexcept { return registry; }

const MaterialEntry* find_material(std::string_view name) noexcept {
    for (const MaterialEntry& e : registry) {
        if (e.name == name) return &e;
    }
    return nullptr;
}

}  // namespace NiTCAD::physics
