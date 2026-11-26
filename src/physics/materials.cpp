#include "NiTCAD/physics/materials.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>

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

// `metal` (Unit 68, MAD6): not solved, its nodes held by an electrode; eps_r 1 is never read.
constexpr InsulatorParameters metal_region{.eps_r = 1.0, .metal = true};
constexpr std::array<Citation, 1> metal_sources{{
    {"NiTCAD (MAD6)", "a metal region is not solved: no parameter is read; the work function is the region's"},
}};
constexpr std::array<ParameterCitation, 11> metal_parameter_sources{{
    {"eps_r", 0, "default"},
    {"thermal.conductivity_W_cmK", 0, "default"},
    {"thermal.conductivity_exponent", 0, "default"},
    {"thermal.heat_capacity_J_cm3K", 0, "default"},
    {"thermal.thermopower_exponent_n", 0, "default"},
    {"thermal.thermopower_exponent_p", 0, "default"},
    {"electron_affinity_eV", 0, "default"},
    {"band_gap_eV", 0, "default"},
    {"tunnel_mass_n", 0, "default"},
    {"tunnel_mass_p", 0, "default"},
    {"tunnel_mass_evb", 0, "default"},
}};

constexpr std::array<MaterialEntry, 10> registry{{
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
    {"metal", MaterialKind::metal, 0, 0, "",
     "metal region held by one electrode contact; work_function or metal = <preset> in the region", metal_sources,
     fixed_insulator<metal_region>, metal_parameter_sources, true},
}};

// --- metal presets (MAD6, MAD14; amendment 4) ---------------------------------------------------

constexpr std::string_view michaelson = "H. B. Michaelson, J. Appl. Phys. 48, 4729 (1977)";
constexpr std::string_view michaelson_where = "Table III, polycrystalline preferred value";
constexpr std::string_view greczynski = "G. Greczynski, arXiv:2405.10919 (2024)";
constexpr std::string_view poly = "degenerate polysilicon: Fermi level at the silicon band edge";
constexpr std::string_view sze = "S. M. Sze, K. K. Ng, Physics of Semiconductor Devices, 3rd ed. (Wiley 2007), ch. 3, p. 180";

constexpr std::array<MetalPreset, 19> presets{{
    {"al", MetalPresetKind::metal, 4.28, {michaelson, michaelson_where}},
    {"au", MetalPresetKind::metal, 5.1, {michaelson, michaelson_where}},
    {"co", MetalPresetKind::metal, 5.0, {michaelson, michaelson_where}},
    {"cr", MetalPresetKind::metal, 4.5, {michaelson, michaelson_where}},
    {"cu", MetalPresetKind::metal, 4.65, {michaelson, michaelson_where}},
    {"mo", MetalPresetKind::metal, 4.6, {michaelson, michaelson_where}},
    {"ni", MetalPresetKind::metal, 5.15, {michaelson, michaelson_where}},
    {"pd", MetalPresetKind::metal, 5.12, {michaelson, michaelson_where}},
    {"pt", MetalPresetKind::metal, 5.65, {michaelson, michaelson_where}},
    {"ru", MetalPresetKind::metal, 4.71, {michaelson, michaelson_where}},
    {"ti", MetalPresetKind::metal, 4.33, {michaelson, michaelson_where}},
    {"w", MetalPresetKind::metal, 4.55, {michaelson, michaelson_where}},
    // UPS on sputtered films; air-exposed films read 4.9 (TiN) and 4.7 (TaN) +- 0.3 eV (M21).
    {"tin", MetalPresetKind::metal, 4.70, {greczynski, "UPS work function, mean of three TiN films (4.69-4.71 eV)"}},
    {"tan", MetalPresetKind::metal, 4.41, {greczynski, "UPS work function of TaN"}},
    {"n_poly", MetalPresetKind::n_poly, 0.0, {poly, "phi_m = chi of silicon"}},
    {"p_poly", MetalPresetKind::p_poly, 0.0, {poly, "phi_m = chi + Eg(T) of silicon"}},
    // Slice 5c (M18): the highest reported barrier of each system, as the table gives it.
    {"nisi", MetalPresetKind::silicide, 0.0, {sze, "Table 4, NiSi (orthorhombic) on n-Si"}, 0.75},
    {"tisi2", MetalPresetKind::silicide, 0.0, {sze, "Table 4, TiSi2 (orthorhombic, C54) on n-Si"}, 0.60},
    {"cosi2", MetalPresetKind::silicide, 0.0, {sze, "Table 4, CoSi2 (cubic) on n-Si"}, 0.64},
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

std::string_view nearest_name(std::string_view word, std::span<const std::string_view> names) {
    const auto distance = [](std::string_view a, std::string_view b) {
        std::vector<std::size_t> row(b.size() + 1);
        for (std::size_t j = 0; j <= b.size(); ++j) row[j] = j;
        for (std::size_t i = 1; i <= a.size(); ++i) {
            std::size_t diagonal = row[0];
            row[0] = i;
            for (std::size_t j = 1; j <= b.size(); ++j) {
                const std::size_t up = row[j];
                row[j] = std::min({row[j] + 1, row[j - 1] + 1, diagonal + (a[i - 1] == b[j - 1] ? 0u : 1u)});
                diagonal = up;
            }
        }
        return row[b.size()];
    };
    const std::size_t limit = std::max<std::size_t>(1, word.size() / 3);
    std::string_view best;
    std::size_t best_distance = limit + 1;
    // A tie goes to a name the word is a dropped-letter form of ('nii': 'nisi' over 'ni').
    const auto dropped_from = [&](std::string_view n) {
        std::size_t i = 0;
        for (const char c : n) i += i < word.size() && word[i] == c;
        return i == word.size() && n.size() > word.size();
    };
    for (const std::string_view n : names) {
        const std::size_t d = distance(word, n);
        if (d < best_distance || (d == best_distance && !best.empty() && dropped_from(n) && !dropped_from(best))) {
            best = n;
            best_distance = d;
        }
    }
    return best;
}

std::uint64_t parameter_hash(const MaterialParameters& m) {
    std::uint64_t h = 14695981039346656037ull;
    const auto bytes = [&h](const void* data, std::size_t n) {
        const auto* b = static_cast<const unsigned char*>(data);
        for (std::size_t i = 0; i < n; ++i) h = (h ^ b[i]) * 1099511628211ull;
    };
    const auto value = [&bytes](double v) {
        const auto bits = std::bit_cast<std::uint64_t>(v);
        bytes(&bits, sizeof bits);
    };
    const unsigned char kind = m.index() == 0 ? 0 : 1;
    bytes(&kind, 1);
    for (const auto& [path, v] : parameter_values(m)) {
        bytes(path.data(), path.size());
        value(v);
    }
    if (const auto* s = std::get_if<SemiconductorParameters>(&m)) {
        const unsigned char direct = s->band_to_band.direct_gap ? 1 : 0;
        bytes(&direct, 1);
        for (const double c : s->dos_mass.conduction) value(c);
        for (const double c : s->dos_mass.valence) value(c);
    }
    return h;
}

std::span<const MaterialEntry> material_registry() noexcept { return registry; }

std::span<const MetalPreset> metal_presets() noexcept { return presets; }

const MetalPreset* find_metal_preset(std::string_view name) noexcept {
    for (const MetalPreset& p : presets) {
        if (p.name == name) return &p;
    }
    return nullptr;
}

const MaterialEntry* find_material(std::string_view name) noexcept {
    for (const MaterialEntry& e : registry) {
        if (e.name == name) return &e;
    }
    return nullptr;
}

}  // namespace NiTCAD::physics
