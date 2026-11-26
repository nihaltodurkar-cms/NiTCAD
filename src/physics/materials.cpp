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

// --- material parameters by member path (DK7) --------------------------------------------------

#define NITCAD_THERMAL_PARAMETERS(X)                                                           \
    X(thermal.conductivity_W_cmK) X(thermal.conductivity_exponent) X(thermal.heat_capacity_J_cm3K) \
    X(thermal.thermopower_exponent_n) X(thermal.thermopower_exponent_p)

#define NITCAD_SEMICONDUCTOR_PARAMETERS(X)                                                     \
    X(eps_r) X(Eg0_eV) X(varshni_alpha_eV_per_K) X(varshni_beta_K) X(Nc300) X(Nv300)           \
    X(electron_affinity_eV) X(electron_mobility.mu_min) X(electron_mobility.mu_max)            \
    X(electron_mobility.N_ref) X(electron_mobility.alpha) X(electron_mobility.T_exponent)      \
    X(hole_mobility.mu_min) X(hole_mobility.mu_max) X(hole_mobility.N_ref)                     \
    X(hole_mobility.alpha) X(hole_mobility.T_exponent) X(lifetime.tau_n0) X(lifetime.tau_p0)   \
    X(lifetime.N_ref) X(auger.Cn) X(auger.Cp) X(bandgap_narrowing.E0_eV)                       \
    X(bandgap_narrowing.N0) X(electron_saturation.v_sat_cm_s) X(electron_saturation.beta)      \
    X(hole_saturation.v_sat_cm_s) X(hole_saturation.beta) X(radiative_cm3_s)                   \
    X(ionization.donor_eV) X(ionization.acceptor_eV) X(ionization.donor_degeneracy)            \
    X(ionization.acceptor_degeneracy) X(richardson.electron) X(richardson.hole)                \
    X(impact_ionization.electron.A_low_per_cm) X(impact_ionization.electron.B_low_V_per_cm)    \
    X(impact_ionization.electron.A_high_per_cm) X(impact_ionization.electron.B_high_V_per_cm)  \
    X(impact_ionization.electron.switch_V_per_cm) X(impact_ionization.hole.A_low_per_cm)       \
    X(impact_ionization.hole.B_low_V_per_cm) X(impact_ionization.hole.A_high_per_cm)           \
    X(impact_ionization.hole.B_high_V_per_cm) X(impact_ionization.hole.switch_V_per_cm)        \
    X(impact_ionization.phonon_energy_eV) X(band_to_band.A_per_cm3_s)                          \
    X(band_to_band.B_V_per_cm) X(band_to_band.electron_mass) X(band_to_band.hole_mass)         \
    X(passler_gap.alpha_eV_per_K) X(passler_gap.Theta_K) X(passler_gap.Delta)                  \
    X(energy_relaxation.tau_n_s) X(energy_relaxation.tau_p_s)                                  \
    NITCAD_THERMAL_PARAMETERS(X)

// Unit 51 (GD5, GD10): an insulator's gate-tunnelling band data after its permittivity and heat.
#define NITCAD_INSULATOR_PARAMETERS(X)                                                         \
    X(eps_r) NITCAD_THERMAL_PARAMETERS(X) X(electron_affinity_eV) X(band_gap_eV) X(tunnel_mass_n) \
    X(tunnel_mass_p) X(tunnel_mass_evb)

}  // namespace

double* semiconductor_parameter(SemiconductorParameters& p, std::string_view path) noexcept {
#define NITCAD_X(member) \
    if (path == #member) return &p.member;
    NITCAD_SEMICONDUCTOR_PARAMETERS(NITCAD_X)
#undef NITCAD_X
    return nullptr;
}

double* insulator_parameter(InsulatorParameters& p, std::string_view path) noexcept {
#define NITCAD_X(member) \
    if (path == #member) return &p.member;
    NITCAD_INSULATOR_PARAMETERS(NITCAD_X)
#undef NITCAD_X
    return nullptr;
}

std::vector<std::pair<std::string_view, double>> parameter_values(const MaterialParameters& m) {
    std::vector<std::pair<std::string_view, double>> out;
#define NITCAD_X(member) out.emplace_back(#member, p.member);
    if (const auto* s = std::get_if<SemiconductorParameters>(&m)) {
        const SemiconductorParameters& p = *s;
        NITCAD_SEMICONDUCTOR_PARAMETERS(NITCAD_X)
    } else {
        const InsulatorParameters& p = std::get<InsulatorParameters>(m);
        NITCAD_INSULATOR_PARAMETERS(NITCAD_X)
    }
#undef NITCAD_X
    return out;
}

std::span<const MaterialEntry> material_registry() noexcept { return registry; }

const MaterialEntry* find_material(std::string_view name) noexcept {
    for (const MaterialEntry& e : registry) {
        if (e.name == name) return &e;
    }
    return nullptr;
}

}  // namespace NiTCAD::physics
