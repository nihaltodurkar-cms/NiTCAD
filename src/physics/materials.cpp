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

// --- dielectrics with band data (MAD5, MAD13; slice 5b) -----------------------------------------
// chi = 4.05 eV (silicon) - the conduction-band offset to silicon. The hole tunnelling mass, and any
// value no source in hand gives, is a general approximate value (owner ruling 2026-11-26), cited
// "approximate"; the thermal set is SiO2's, approximate for every dielectric.

constexpr double chi_si = 4.05;
constexpr double mass_general = 0.50;  // m0, the value SiO2 already assumes (owner ruling 2026-10-28)

constexpr InsulatorParameters dielectric(double eps_r, double chi, double gap, double mass_n) {
    return {.eps_r = eps_r, .thermal = silicon_dioxide_parameters.thermal, .electron_affinity_eV = chi,
            .band_gap_eV = gap, .tunnel_mass_n = mass_n, .tunnel_mass_p = mass_general, .tunnel_mass_evb = 0.0};
}

struct Cite {
    std::size_t source;
    std::string_view basis;
};

// eps_r, chi, gap and electron mass cited per material; thermal and hole mass approximate (source
// `approx`); tunnel_mass_evb 0 means the electron mass (insulator.hpp).
constexpr std::array<ParameterCitation, 11> dielectric_cites(Cite eps, Cite chi, Cite gap, Cite mass_n, std::size_t approx) {
    return {{{"eps_r", eps.source, eps.basis},
             {"thermal.conductivity_W_cmK", approx, "approximate"},
             {"thermal.conductivity_exponent", approx, "approximate"},
             {"thermal.heat_capacity_J_cm3K", approx, "approximate"},
             {"thermal.thermopower_exponent_n", approx, "approximate"},
             {"thermal.thermopower_exponent_p", approx, "approximate"},
             {"electron_affinity_eV", chi.source, chi.basis},
             {"band_gap_eV", gap.source, gap.basis},
             {"tunnel_mass_n", mass_n.source, mass_n.basis},
             {"tunnel_mass_p", approx, "approximate"},
             {"tunnel_mass_evb", approx, "default"}}};
}

constexpr Citation robertson{"J. Robertson, Eur. Phys. J. Appl. Phys. 28, 265 (2004)",
                             "Table 2: K, gap, conduction-band offset to Si; chi = 4.05 eV - offset"};
constexpr Citation salmani{"M. Salmani-Jelodar et al., arXiv:1502.06178 (2015)",
                           "Table I: electron tunnelling mass (secondary, from Yeo 2002, Robertson 2000)"};
constexpr Citation approx_note{"NiTCAD, owner ruling 2026-11-26 (Unit 68 slice 5b)",
                               "general approximate values where no source in hand gives one: hole tunnelling "
                               "mass 0.50 m0, SiO2's thermal set, and the values the material names"};

constexpr std::array<Citation, 3> robertson_salmani_sources{{robertson, salmani, approx_note}};
constexpr std::array<Citation, 2> robertson_sources{{robertson, approx_note}};

constexpr InsulatorParameters hfo2 = dielectric(25.0, chi_si - 1.4, 5.8, 0.11);
constexpr InsulatorParameters al2o3 = dielectric(9.0, chi_si - 2.8, 8.8, 0.35);
constexpr InsulatorParameters si3n4 = dielectric(7.0, chi_si - 2.4, 5.3, 0.40);
constexpr InsulatorParameters la2o3 = dielectric(30.0, chi_si - 2.3, 6.0, 0.27);
constexpr auto robertson_salmani_cites = dielectric_cites({0, "measured"}, {0, "derived"}, {0, "measured"}, {1, "secondary"}, 2);

constexpr InsulatorParameters zro2 = dielectric(25.0, chi_si - 1.5, 5.8, mass_general);
constexpr InsulatorParameters tio2 = dielectric(80.0, chi_si - 0.0, 3.5, mass_general);
constexpr InsulatorParameters ta2o5 = dielectric(22.0, chi_si - 0.35, 4.4, mass_general);
constexpr auto robertson_cites = dielectric_cites({0, "measured"}, {0, "derived"}, {0, "measured"}, {1, "approximate"}, 1);

// MgO (MAD13): amorphous MgO on (100)Si by internal photoemission, the conduction band 3.37 eV
// above the silicon valence band top (offset 3.37 - 1.12 = 2.25 eV), gap 6.1 eV; eps_r approximate.
constexpr std::array<Citation, 2> mgo_sources{{
    {"V. V. Afanas'ev et al., Appl. Phys. Lett. 96, 052103 (2010), doi:10.1063/1.3294328",
     "abstract: a-MgO gap 6.1 eV, conduction band 3.37 eV above the Si valence band top"},
    {approx_note.source, "eps_r 9.8, electron tunnelling mass 0.50 m0, hole mass, thermal set"},
}};
constexpr InsulatorParameters mgo = dielectric(9.8, chi_si - (3.37 - 1.12), 6.1, mass_general);
constexpr auto mgo_cites = dielectric_cites({1, "approximate"}, {0, "derived"}, {0, "measured"}, {1, "approximate"}, 1);

// AlN as a gate or passivation insulator (MAD5): gap 6.07 eV, chi 1.83-1.94 eV (the midpoint).
constexpr std::array<Citation, 2> aln_sources{{
    {"arXiv:2407.17360 (2024), Si/AlN", "gap 6.07 eV, electron affinity 1.83-1.94 eV (midpoint taken)"},
    {approx_note.source, "eps_r 8.5, electron tunnelling mass 0.50 m0, hole mass, thermal set"},
}};
constexpr InsulatorParameters aln = dielectric(8.5, 1.885, 6.07, mass_general);
constexpr auto aln_cites = dielectric_cites({1, "approximate"}, {0, "derived"}, {0, "measured"}, {1, "approximate"}, 1);

// Sapphire (MAD13): crystalline Al2O3, its own set; gap and offset as al2o3 (M5), eps_r along c
// approximate (c-plane substrates).
constexpr std::array<Citation, 3> sapphire_sources{{
    robertson, salmani, {approx_note.source, "eps_r 11.6 along c, hole mass, thermal set"}}};
constexpr InsulatorParameters sapphire = dielectric(11.6, chi_si - 2.8, 8.8, 0.35);
constexpr auto sapphire_cites = dielectric_cites({2, "approximate"}, {0, "secondary"}, {0, "secondary"}, {1, "secondary"}, 2);

// h-BN (MAD13): gap 5.955 eV (indirect, M16); out-of-plane eps_r and chi approximate.
constexpr std::array<Citation, 2> hbn_sources{{
    {"G. Cassabois, P. Valvin, B. Gil, arXiv:1512.02962 (Nat. Photon. 10, 262, 2016)", "indirect gap 5.955 eV"},
    {approx_note.source, "eps_r 3.76 out of plane, chi 2.0 eV, electron tunnelling mass 0.50 m0, hole mass, thermal set"},
}};
constexpr InsulatorParameters hbn = dielectric(3.76, 2.0, 5.955, mass_general);
constexpr auto hbn_cites = dielectric_cites({1, "approximate"}, {1, "approximate"}, {0, "measured"}, {1, "approximate"}, 1);

// SiON (MAD13): x the Si3N4 share, every value linear between SiO2 (its tunnelling set) and si3n4.
constexpr std::array<Citation, 4> sion_sources{{
    {"NiTCAD (MAD13)", "linear in x between sio2 (x = 0, its Unit 51 tunnelling set) and si3n4 (x = 1)"},
    robertson, salmani, approx_note}};
constexpr auto sion_cites = dielectric_cites({0, "interpolated"}, {0, "interpolated"}, {0, "interpolated"}, {0, "interpolated"}, 3);

std::expected<MaterialParameters, base::Error> sion(std::span<const double> a) {
    const double x = a[0];
    if (!(x >= 0.0 && x <= 1.0)) {
        return std::unexpected(base::Error{base::ErrorCode::invalid_input, "sion: the Si3N4 share x must lie in [0, 1]"});
    }
    const InsulatorParameters& o = silicon_dioxide_tunnelling_parameters;
    const auto mix = [x](double a0, double a1) { return (1.0 - x) * a0 + x * a1; };
    InsulatorParameters p = o;
    p.eps_r = mix(o.eps_r, si3n4.eps_r);
    p.electron_affinity_eV = mix(o.electron_affinity_eV, si3n4.electron_affinity_eV);
    p.band_gap_eV = mix(o.band_gap_eV, si3n4.band_gap_eV);
    p.tunnel_mass_n = mix(o.tunnel_mass_n, si3n4.tunnel_mass_n);
    p.tunnel_mass_p = mix(o.tunnel_mass_p, si3n4.tunnel_mass_p);
    return p;
}

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

constexpr std::array<MaterialEntry, 22> registry{{
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
    {"hfo2", MaterialKind::insulator, 0, 0, "", "HfO2, amorphous gate dielectric", robertson_salmani_sources,
     fixed_insulator<hfo2>, robertson_salmani_cites, true},
    {"al2o3", MaterialKind::insulator, 0, 0, "", "Al2O3, amorphous gate dielectric", robertson_salmani_sources,
     fixed_insulator<al2o3>, robertson_salmani_cites, true},
    {"si3n4", MaterialKind::insulator, 0, 0, "", "Si3N4", robertson_salmani_sources, fixed_insulator<si3n4>,
     robertson_salmani_cites, true},
    {"zro2", MaterialKind::insulator, 0, 0, "", "ZrO2", robertson_sources, fixed_insulator<zro2>, robertson_cites, true},
    {"la2o3", MaterialKind::insulator, 0, 0, "", "La2O3", robertson_salmani_sources, fixed_insulator<la2o3>,
     robertson_salmani_cites, true},
    {"tio2", MaterialKind::insulator, 0, 0, "", "TiO2 (no conduction-band offset to Si)", robertson_sources,
     fixed_insulator<tio2>, robertson_cites, true},
    {"ta2o5", MaterialKind::insulator, 0, 0, "", "Ta2O5", robertson_sources, fixed_insulator<ta2o5>, robertson_cites, true},
    {"mgo", MaterialKind::insulator, 0, 0, "", "MgO, amorphous on Si", mgo_sources, fixed_insulator<mgo>, mgo_cites, true},
    {"aln_dielectric", MaterialKind::insulator, 0, 0, "", "AlN as a gate or passivation insulator", aln_sources,
     fixed_insulator<aln>, aln_cites, true},
    {"sapphire", MaterialKind::insulator, 0, 0, "", "sapphire, crystalline Al2O3 substrate (eps_r along c)",
     sapphire_sources, fixed_insulator<sapphire>, sapphire_cites, true},
    {"hbn", MaterialKind::insulator, 0, 0, "", "h-BN, gate dielectric for layered channels (eps_r out of plane)",
     hbn_sources, fixed_insulator<hbn>, hbn_cites, true},
    {"sion", MaterialKind::insulator, 1, 1, "sion takes the Si3N4 share x, 0 <= x <= 1",
     "SiON, linear between SiO2 (x = 0) and Si3N4 (x = 1)", sion_sources, sion, sion_cites, true},
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
