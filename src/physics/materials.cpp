#include "NiTCAD/physics/materials.hpp"

#include "NiTCAD/physics/iii_v.hpp"
#include "NiTCAD/physics/gallium_oxide.hpp"
#include "NiTCAD/physics/nitrides.hpp"
#include "NiTCAD/physics/sige.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>

namespace NiTCAD::physics {
namespace {
// --- material parameters by member path (DK7) --------------------------------------------------

#define NITCAD_THERMAL_PARAMETERS(X)                                                           \
    X(thermal.conductivity_W_cmK) X(thermal.conductivity_exponent) X(thermal.heat_capacity_J_cm3K) \
    X(thermal.thermopower_exponent_n) X(thermal.thermopower_exponent_p)

// Unit 68 slice 2: the Caughey-Thomas exponents added to a mobility set, and a whole set.
#define NITCAD_CT_EXPONENTS(X, s) X(s.mu_min_T_exponent) X(s.N_ref_T_exponent) X(s.alpha_T_exponent)
#define NITCAD_CT(X, s) X(s.mu_min) X(s.mu_max) X(s.N_ref) X(s.alpha) X(s.T_exponent) NITCAD_CT_EXPONENTS(X, s)

#define NITCAD_SEMICONDUCTOR_PARAMETERS(X)                                                   \
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
    X(energy_relaxation.tau_n_s) X(energy_relaxation.tau_p_s) X(alloy_gap.fraction)             \
    X(alloy_gap.alpha_a_eV_per_K) X(alloy_gap.beta_a_K) X(alloy_gap.alpha_b_eV_per_K)          \
    X(alloy_gap.beta_b_K) NITCAD_CT_EXPONENTS(X, electron_mobility)                            \
    NITCAD_CT_EXPONENTS(X, hole_mobility) X(alloy_mobility.fraction)                           \
    X(alloy_mobility.C_mu_electron) NITCAD_CT(X, alloy_mobility.electron_b)                    \
    NITCAD_CT(X, alloy_mobility.hole_b) NITCAD_THERMAL_PARAMETERS(X) X(alloy_mobility.C_mu_hole)  \
    X(electron_saturation.temperature_A) X(hole_saturation.temperature_A) X(electron_saturation.F_c_V_per_cm) \
    X(electron_saturation.gamma) X(hole_saturation.F_c_V_per_cm) X(hole_saturation.gamma)

// Unit 51 (GD5, GD10): an insulator's gate-tunnelling band data after its permittivity and heat.
#define NITCAD_INSULATOR_PARAMETERS(X)                                                         \
    X(eps_r) NITCAD_THERMAL_PARAMETERS(X) X(electron_affinity_eV) X(band_gap_eV) X(tunnel_mass_n) \
    X(tunnel_mass_p) X(tunnel_mass_evb)


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
// chi = 4.05 eV (silicon) - the conduction-band offset to silicon. The hole tunnelling masses, and
// the electron masses no source in hand gives, are the owner's values (2026-11-26); every unsourced value is cited
// "approximate" (owner ruling 2026-11-26); the thermal set is SiO2's, approximate for every dielectric.

constexpr double chi_si = 4.05;

// k in W/cm K and heat capacity in J/cm3 K: the owner's thin-film starting values near 300 K
// (2026-11-26); the conductivity exponent (0) and thermopower exponents (-0.5) stay SiO2's.
constexpr InsulatorParameters dielectric(double eps_r, double chi, double gap, double mass_n, double mass_p, double k,
                                         double heat_capacity) {
    ThermalParameters t = silicon_dioxide_parameters.thermal;
    t.conductivity_W_cmK = k;
    t.heat_capacity_J_cm3K = heat_capacity;
    return {.eps_r = eps_r, .thermal = t, .electron_affinity_eV = chi, .band_gap_eV = gap, .tunnel_mass_n = mass_n,
            .tunnel_mass_p = mass_p, .tunnel_mass_evb = 0.0};
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
                               "general approximate values where no source in hand gives one, and the tunnelling "
                               "masses, thermal conductivities and heat capacities the owner supplied (2026-11-26, "
                               "from thin-film reports not read here); the conductivity and thermopower exponents are "
                               "SiO2's; and the values each material names"};

constexpr std::array<Citation, 3> robertson_salmani_sources{{robertson, salmani, approx_note}};
constexpr std::array<Citation, 2> robertson_sources{{robertson, approx_note}};

constexpr InsulatorParameters hfo2 = dielectric(25.0, chi_si - 1.4, 5.8, 0.11, 0.60, 0.010, 1.16);
constexpr InsulatorParameters si3n4 = dielectric(7.0, chi_si - 2.4, 5.3, 0.40, 0.35, 0.012, 2.30);
constexpr InsulatorParameters la2o3 = dielectric(30.0, chi_si - 2.3, 6.0, 0.27, 0.50, 0.010, 2.15);
constexpr auto robertson_salmani_cites = dielectric_cites({0, "measured"}, {0, "derived"}, {0, "measured"}, {1, "secondary"}, 2);

// Al2O3 as deposited (ALD) gate dielectric: gap 6.4 eV and offset 1.6 eV (chi 2.45 eV), the ALD
// row of Robertson & Wallace 2015 as the owner gave it (2026-11-26; not read here). Robertson 2004's
// 8.8 / 2.8 eV is marked "not ALD" and is sapphire's below.
constexpr std::array<Citation, 4> al2o3_sources{{
    robertson, salmani, approx_note,
    {"J. Robertson, R. M. Wallace, Mater. Sci. Eng. R 88, 1 (2015)", "ALD Al2O3: gap 6.4 eV, offset 1.6 eV (owner-supplied)"}}};
constexpr InsulatorParameters al2o3 = dielectric(9.0, chi_si - 1.6, 6.4, 0.35, 1.30, 0.015, 2.50);
constexpr auto al2o3_cites = dielectric_cites({0, "measured"}, {3, "secondary"}, {3, "secondary"}, {1, "secondary"}, 2);

constexpr InsulatorParameters zro2 = dielectric(25.0, chi_si - 1.5, 5.8, 0.28, 0.60, 0.005, 2.70);
constexpr InsulatorParameters tio2 = dielectric(80.0, chi_si - 0.0, 3.5, 0.50, 3.80, 0.025, 3.03);
constexpr InsulatorParameters ta2o5 = dielectric(22.0, chi_si - 0.35, 4.4, 0.25, 1.00, 0.0012, 2.52);
constexpr auto robertson_cites = dielectric_cites({0, "measured"}, {0, "derived"}, {0, "measured"}, {1, "approximate"}, 1);

// MgO (MAD13): amorphous MgO on (100)Si by internal photoemission, the conduction band 3.37 eV
// above the silicon valence band top (offset 3.37 - 1.12 = 2.25 eV), gap 6.1 eV; eps_r approximate.
// chi 1.80 eV is Si-referenced (4.05 - 2.25), an interface-derived effective affinity, not MgO's
// vacuum-referenced surface affinity (-0.40 eV, Rohlfing et al., PRL 91, 256802, 2003, via Heinisch
// et al. 2012), which does not set the Si/MgO barrier.
constexpr std::array<Citation, 2> mgo_sources{{
    {"V. V. Afanas'ev et al., Appl. Phys. Lett. 96, 052103 (2010), doi:10.1063/1.3294328",
     "abstract: a-MgO gap 6.1 eV, conduction band 3.37 eV above the Si valence band top"},
    {approx_note.source, "eps_r 9.8, tunnelling masses (owner-supplied), thermal set"},
}};
constexpr InsulatorParameters mgo = dielectric(9.8, chi_si - (3.37 - 1.12), 6.1, 0.35, 3.00, 0.10, 3.33);
constexpr auto mgo_cites = dielectric_cites({1, "approximate"}, {0, "derived"}, {0, "measured"}, {1, "approximate"}, 1);

// AlN as a gate or passivation insulator (MAD5): gap 6.07 eV, chi 1.83-1.94 eV (midpoint, 1.89).
constexpr std::array<Citation, 2> aln_sources{{
    {"arXiv:2407.17360 (2024), Si/AlN", "gap 6.07 eV, electron affinity 1.83-1.94 eV (midpoint 1.89 taken)"},
    {approx_note.source, "eps_r 8.5, tunnelling masses (owner-supplied), thermal set"},
}};
constexpr InsulatorParameters aln = dielectric(8.5, 1.89, 6.07, 0.33, 3.50, 0.10, 2.40);
constexpr auto aln_cites = dielectric_cites({1, "approximate"}, {0, "derived"}, {0, "measured"}, {1, "approximate"}, 1);

// Sapphire (MAD13): crystalline Al2O3, its own set: gap 8.8 eV and offset 2.8 eV, Robertson 2004's
// Al2O3 row marked "not ALD"; the electron mass is alumina's (Salmani), a proxy; eps_r along c
// approximate (c-plane substrates); thermal conductivity along c (owner-supplied, 2026-11-26).
constexpr std::array<Citation, 3> sapphire_sources{{
    {robertson.source, "Table 2, Al2O3 (not ALD): gap 8.8 eV, offset 2.8 eV; chi = 4.05 eV - offset"},
    salmani,
    {approx_note.source, "eps_r 11.6 along c, hole tunnelling mass, thermal set along c (owner-supplied)"}}};
constexpr InsulatorParameters sapphire = dielectric(11.6, chi_si - 2.8, 8.8, 0.35, 0.36, 0.30, 3.04);
constexpr auto sapphire_cites = dielectric_cites({2, "approximate"}, {0, "derived"}, {0, "measured"}, {1, "secondary"}, 2);

// h-BN (MAD13): gap 5.955 eV (indirect, M16); out-of-plane eps_r and chi approximate.
constexpr std::array<Citation, 2> hbn_sources{{
    {"G. Cassabois, P. Valvin, B. Gil, arXiv:1512.02962 (Nat. Photon. 10, 262, 2016)", "indirect gap 5.955 eV"},
    {approx_note.source, "eps_r 3.76 out of plane, chi 2.0 eV, tunnelling masses (owner-supplied), thermal set"},
}};
constexpr InsulatorParameters hbn = dielectric(3.76, 2.0, 5.955, 2.21, 1.33, 0.033, 1.80);  // k out of plane
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

// --- SiGe (MAD7; slice 1, amendment 5) -----------------------------------------------------------

#define NITCAD_PATH(member) std::string_view{#member},
constexpr std::string_view semiconductor_paths[] = {NITCAD_SEMICONDUCTOR_PARAMETERS(NITCAD_PATH)};
#undef NITCAD_PATH

// One citation per semiconductor path: `fallback`, overridden by the rules in order (a rule names a
// path, or a group with its trailing dot).
template <std::size_t N>
constexpr auto semiconductor_cites(Cite fallback, const std::array<ParameterCitation, N>& rules) {
    std::array<ParameterCitation, std::size(semiconductor_paths)> out{};
    for (std::size_t i = 0; i < out.size(); ++i) {
        const std::string_view path = semiconductor_paths[i];
        out[i] = {path, fallback.source, fallback.basis};
        for (const ParameterCitation& r : rules) {
            if (path == r.path || (r.path.ends_with('.') && path.starts_with(r.path))) out[i] = {path, r.source, r.basis};
        }
    }
    return out;
}

constexpr Citation ioffe_sige{"Ioffe NSM archive, SiGe (F. Schaffler, in M. Levinshtein et al., eds., Properties of "
                              "Advanced Semiconductor Materials, Wiley 2001)",
                              "eps_r; gap fit; Rieger-Vogl offset, spin-orbit splitting, deformation potentials; lattice "
                              "constant, Poisson ratio; hole DOS mass; mobility lines; specific heat and density"};
constexpr Citation silicon_inherited{"NiTCAD silicon (Couderc et al. 2014; legacy materials.py)",
                                     "inherited: the x = 0 anchor (rulings S1, S5, S8, S11, F1b)"};
constexpr Citation sige_rulings{"owner rulings 2026-11-26 (Unit 68 amendment 5)",
                                "BGN (1 - 0.5 x) silicon's (owner's fit); saturation velocity ratios of Sasso and Rinaldi "
                                "2010 (owner-supplied, unread); ionization energies linear between Si and Ge; Varshni end "
                                "points of Palankovski and Selberherr"};
constexpr Citation sige_none{"NiTCAD: no SiGe data in hand",
                             "off: impact ionization, band-to-band tunnelling, thermal conductivity (its temperature "
                             "exponent awaits the owner)"};
constexpr std::array<Citation, 4> sige_sources{{ioffe_sige, silicon_inherited, sige_rulings, sige_none}};

constexpr auto sige_cites = semiconductor_cites(
    {1, "secondary"},
    std::array<ParameterCitation, 17>{{{"eps_r", 0, "fit"},
                                       {"Eg0_eV", 0, "derived"},
                                       {"electron_affinity_eV", 0, "derived"},
                                       {"electron_mobility.mu_max", 0, "derived"},
                                       {"hole_mobility.mu_max", 0, "derived"},
                                       {"bandgap_narrowing.E0_eV", 2, "approximate"},
                                       {"electron_saturation.v_sat_cm_s", 2, "secondary"},
                                       {"hole_saturation.v_sat_cm_s", 2, "secondary"},
                                       {"ionization.donor_eV", 2, "interpolated"},
                                       {"ionization.acceptor_eV", 2, "interpolated"},
                                       {"impact_ionization.", 3, "default"},
                                       {"band_to_band.", 3, "default"},
                                       {"alloy_gap.", 2, "secondary"},
                                       {"alloy_gap.fraction", 0, "derived"},
                                       {"thermal.conductivity_W_cmK", 3, "default"},
                                       {"thermal.conductivity_exponent", 3, "default"},
                                       {"thermal.heat_capacity_J_cm3K", 0, "derived"}}});
constexpr auto si_strained_cites = semiconductor_cites(
    {1, "secondary"}, std::array<ParameterCitation, 3>{{{"Eg0_eV", 0, "derived"},
                                                       {"electron_affinity_eV", 0, "derived"},
                                                       {"alloy_gap.", 3, "default"}}});

std::expected<MaterialParameters, base::Error> sige(std::span<const double> a) {
    auto p = sige_parameters(a[0], false);
    if (!p) return std::unexpected(std::move(p.error()));
    return *p;
}
std::expected<MaterialParameters, base::Error> sige_strained(std::span<const double> a) {
    auto p = sige_parameters(a[0], true);
    if (!p) return std::unexpected(std::move(p.error()));
    return *p;
}
std::expected<MaterialParameters, base::Error> si_strained(std::span<const double> a) {
    auto p = strained_silicon_parameters(a[0]);
    if (!p) return std::unexpected(std::move(p.error()));
    return *p;
}

// --- nitrides (MAD8; slice 2, amendment 6) -------------------------------------------------------

constexpr Citation levinshtein_nitrides{
    "M. Levinshtein, S. Rumyantsev, M. Shur, eds., Properties of Advanced Semiconductor Materials (Wiley 2001), ch. 1-3",
    "eps_r; GaN gap and Varshni; electron and hole DOS masses; radiative coefficient; 300 K thermal conductivity; "
    "specific heat and density"};
constexpr Citation piprek_aln{"J. Piprek, ed., Nitride Semiconductor Devices: Principles and Simulation (Wiley-VCH 2007), "
                              "Section 2",
                              "AlN gap Eg(0) 6.10 eV and Varshni alpha 2.63e-3 eV/K, beta 2082 K"};
constexpr Citation nitride_owner{"owner values and rulings 2026-11-27 (Unit 68 amendment 6)",
                                 "InN gap and Varshni; Caughey-Thomas mobility and its exponents (AlN alpha +0.18 "
                                 "provisional); constant hole mobility; saturation (AlN beta 1, not 0.45); SRH; Auger; "
                                 "dopant levels (InN donor fully ionized fallback); GaN impact ionization; thermal "
                                 "exponents; Q_c alignment, bowing, alloy mobility"};
constexpr Citation nitride_none{"NiTCAD: no nitride data in hand",
                                "off: band-gap narrowing, band-to-band tunnelling; AlN, InN and alloy impact ionization; "
                                "alloy thermal conductivity"};
constexpr std::array<Citation, 4> nitride_sources{{levinshtein_nitrides, piprek_aln, nitride_owner, nitride_none}};

// Rules every nitride shares, then each set's own (later rules win).
template <std::size_t N>
constexpr auto nitride_rules(const std::array<ParameterCitation, N>& own) {
    constexpr std::array<ParameterCitation, 12> common{{{"electron_affinity_eV", 2, "derived"},
                                                        {"Nc300", 0, "derived"},
                                                        {"Nv300", 0, "derived"},
                                                        {"electron_mobility.", 2, "secondary"},
                                                        {"hole_mobility.", 2, "secondary"},
                                                        {"lifetime.", 2, "secondary"},
                                                        {"auger.", 2, "secondary"},
                                                        {"electron_saturation.", 2, "secondary"},
                                                        {"hole_saturation.", 2, "secondary"},
                                                        {"ionization.", 2, "secondary"},
                                                        {"bandgap_narrowing.", 3, "default"},
                                                        {"band_to_band.", 3, "default"}}};
    std::array<ParameterCitation, 12 + N> out{};
    for (std::size_t i = 0; i < 12; ++i) out[i] = common[i];
    for (std::size_t i = 0; i < N; ++i) out[12 + i] = own[i];
    return out;
}
constexpr auto gan_cites = semiconductor_cites(
    {0, "secondary"}, nitride_rules(std::array<ParameterCitation, 4>{{{"electron_affinity_eV", 2, "secondary"},
                                                                      {"impact_ionization.", 2, "secondary"},
                                                                      {"thermal.conductivity_exponent", 2, "secondary"},
                                                                      {"thermal.heat_capacity_J_cm3K", 0, "derived"}}}));
constexpr auto aln_semiconductor_cites = semiconductor_cites(
    {0, "secondary"}, nitride_rules(std::array<ParameterCitation, 7>{{{"Eg0_eV", 1, "secondary"},
                                                                      {"varshni_alpha_eV_per_K", 1, "secondary"},
                                                                      {"varshni_beta_K", 1, "secondary"},
                                                                      {"electron_saturation.beta", 2, "approximate"},
                                                                      {"impact_ionization.", 3, "default"},
                                                                      {"thermal.conductivity_exponent", 2, "secondary"},
                                                                      {"thermal.heat_capacity_J_cm3K", 0, "derived"}}}));
constexpr auto inn_cites = semiconductor_cites(
    {0, "secondary"}, nitride_rules(std::array<ParameterCitation, 7>{{{"Eg0_eV", 2, "secondary"},
                                                                      {"varshni_alpha_eV_per_K", 2, "secondary"},
                                                                      {"varshni_beta_K", 2, "secondary"},
                                                                      {"ionization.donor_eV", 2, "approximate"},
                                                                      {"impact_ionization.", 3, "default"},
                                                                      {"thermal.conductivity_exponent", 2, "secondary"},
                                                                      {"thermal.heat_capacity_J_cm3K", 0, "derived"}}}));
constexpr auto nitride_alloy_cites = semiconductor_cites(
    {0, "interpolated"},
    nitride_rules(std::array<ParameterCitation, 9>{{{"Eg0_eV", 2, "derived"},
                                                    {"alloy_gap.", 0, "derived"},
                                                    {"alloy_mobility.", 2, "secondary"},
                                                    {"alloy_mobility.fraction", 0, "derived"},
                                                    {"Nc300", 2, "derived"},
                                                    {"impact_ionization.", 3, "default"},
                                                    {"thermal.conductivity_W_cmK", 3, "default"},
                                                    {"thermal.conductivity_exponent", 3, "default"},
                                                    {"thermal.heat_capacity_J_cm3K", 0, "interpolated"}}}));

template <SemiconductorParameters (*make)()>
std::expected<MaterialParameters, base::Error> nitride_binary(std::span<const double>) {
    return make();
}
template <std::expected<SemiconductorParameters, base::Error> (*make)(double)>
std::expected<MaterialParameters, base::Error> nitride_alloy(std::span<const double> a) {
    auto p = make(a[0]);
    if (!p) return std::unexpected(std::move(p.error()));
    return *p;
}

// --- III-V (MAD10; slice 3, amendment 7) ---------------------------------------------------------

constexpr Citation vurgaftman{"I. Vurgaftman, J. R. Meyer, L. R. Ram-Mohan, J. Appl. Phys. 89, 5815 (2001) (M1), "
                              "as the owner gave the values 2026-11-27",
                              "gaps and Varshni (GaP, AlAs: X valley), electron and X-valley masses, Luttinger "
                              "gamma_1 and gamma_2, valence-band offsets, InGaAs and InAlAs bowing"};
constexpr Citation sotoodeh{"M. Sotoodeh, A. H. Khalid, A. A. Rezazadeh, J. Appl. Phys. 87, 2890 (2000), Table I, "
                            "as the owner gave it",
                            "Caughey-Thomas mobility of GaAs, InP, InAs, GaP and AlAs; GaSb and InSb from a "
                            "temperature-extended fit (InSb exponents provisional)"};
constexpr Citation palankovski_quay{"V. Palankovski (thesis, TU Wien 2000), Tables 3.30-3.32; R. Quay et al., Mater. "
                                    "Sci. Semicond. Process. 3, 149 (2000), Table 1",
                                    "saturation velocities and their temperature coefficient A; alloy C_mu"};
constexpr Citation iii_v_owner{"owner values and rulings 2026-11-27 (Unit 68 amendment 7)",
                               "permittivity; SRH lifetimes (representative or placeholder); Auger; radiative; thermal "
                               "data; GaSb and InSb saturation; alloy thermal fit C_k; rulings: VBO alignment, full "
                               "ionization, beta 2"};
constexpr Citation iii_v_none{"NiTCAD: no value given",
                              "off: band-gap narrowing, impact ionization, Kane band-to-band pair; InSb and alloy "
                              "lifetimes (SRH refused); GaP and AlAs Auger and radiative"};
constexpr std::array<Citation, 5> iii_v_sources{{vurgaftman, sotoodeh, palankovski_quay, iii_v_owner, iii_v_none}};

// Rules every III-V set shares, then each set's own (later rules win).
template <std::size_t N>
constexpr auto iii_v_rules(const std::array<ParameterCitation, N>& own) {
    constexpr std::array<ParameterCitation, 23> common{{{"Eg0_eV", 0, "secondary"},
                                                        {"varshni_alpha_eV_per_K", 0, "secondary"},
                                                        {"varshni_beta_K", 0, "secondary"},
                                                        {"Nc300", 0, "derived"},
                                                        {"Nv300", 0, "derived"},
                                                        {"electron_affinity_eV", 0, "derived"},
                                                        {"electron_mobility.", 1, "secondary"},
                                                        {"hole_mobility.", 1, "secondary"},
                                                        {"electron_saturation.", 2, "secondary"},
                                                        {"hole_saturation.", 2, "secondary"},
                                                        {"electron_saturation.beta", 4, "default"},
                                                        {"hole_saturation.beta", 4, "default"},
                                                        {"lifetime.", 3, "approximate"},
                                                        {"lifetime.N_ref", 4, "default"},
                                                        {"ionization.", 4, "default"},
                                                        {"richardson.", 4, "default"},
                                                        {"bandgap_narrowing.", 4, "default"},
                                                        {"impact_ionization.", 4, "default"},
                                                        {"band_to_band.", 4, "default"},
                                                        {"thermal.thermopower_exponent_n", 4, "default"},
                                                        {"thermal.thermopower_exponent_p", 4, "default"},
                                                        {"thermal.conductivity_exponent", 3, "approximate"},
                                                        {"thermal.heat_capacity_J_cm3K", 3, "derived"}}};
    std::array<ParameterCitation, 23 + N> out{};
    for (std::size_t i = 0; i < 23; ++i) out[i] = common[i];
    for (std::size_t i = 0; i < N; ++i) out[23 + i] = own[i];
    return out;
}
constexpr auto iii_v_tunnelling_cites = semiconductor_cites(
    {3, "secondary"}, iii_v_rules(std::array<ParameterCitation, 2>{{{"band_to_band.electron_mass", 0, "derived"},
                                                                     {"band_to_band.hole_mass", 0, "derived"}}}));
constexpr auto gaas_2001_cites = semiconductor_cites(
    {3, "secondary"},
    iii_v_rules(std::array<ParameterCitation, 5>{{{"band_to_band.electron_mass", 0, "derived"},
                                                  {"band_to_band.hole_mass", 0, "derived"},
                                                  {"lifetime.tau_n0", 3, "secondary"},
                                                  {"lifetime.tau_p0", 3, "secondary"},
                                                  {"thermal.conductivity_exponent", 3, "secondary"}}}));
constexpr auto iii_v_indirect_cites = semiconductor_cites(
    {3, "secondary"}, iii_v_rules(std::array<ParameterCitation, 2>{{{"auger.", 4, "default"},
                                                                     {"radiative_cm3_s", 4, "default"}}}));
constexpr auto gasb_cites = semiconductor_cites(
    {3, "secondary"},
    iii_v_rules(std::array<ParameterCitation, 4>{{{"electron_saturation.", 3, "approximate"},
                                                  {"hole_saturation.", 3, "approximate"},
                                                  {"electron_saturation.beta", 4, "default"},
                                                  {"thermal.conductivity_exponent", 3, "secondary"}}}));
constexpr auto insb_cites = semiconductor_cites(
    {3, "secondary"},
    iii_v_rules(std::array<ParameterCitation, 8>{{{"electron_saturation.", 3, "approximate"},
                                                  {"hole_saturation.", 3, "approximate"},
                                                  {"electron_mobility.T_exponent", 3, "approximate"},
                                                  {"electron_mobility.N_ref_T_exponent", 3, "approximate"},
                                                  {"hole_mobility.T_exponent", 3, "approximate"},
                                                  {"hole_mobility.N_ref_T_exponent", 3, "approximate"},
                                                  {"lifetime.tau_n0", 4, "default"},
                                                  {"lifetime.tau_p0", 4, "default"}}}));
constexpr auto iii_v_alloy_cites = semiconductor_cites(
    {0, "interpolated"},
    iii_v_rules(std::array<ParameterCitation, 12>{{{"Eg0_eV", 0, "derived"},
                                                   {"alloy_gap.", 0, "derived"},
                                                   {"alloy_mobility.", 1, "secondary"},
                                                   {"alloy_mobility.C_mu_electron", 2, "secondary"},
                                                   {"alloy_mobility.C_mu_hole", 2, "secondary"},
                                                   {"alloy_mobility.fraction", 0, "derived"},
                                                   {"lifetime.tau_n0", 4, "default"},
                                                   {"lifetime.tau_p0", 4, "default"},
                                                   {"band_to_band.electron_mass", 0, "derived"},
                                                   {"band_to_band.hole_mass", 0, "derived"},
                                                   {"thermal.conductivity_W_cmK", 3, "fit"},
                                                   {"thermal.conductivity_exponent", 3, "approximate"}}}));

template <SemiconductorParameters (*make)()>
std::expected<MaterialParameters, base::Error> iii_v_binary(std::span<const double>) {
    return make();
}

// --- beta-Ga2O3 (MAD11; slice 4b, amendment 9) -------------------------------------------------

// The owner's values of 2026-11-27, each from its source; [010] throughout (owner ruling 2).
constexpr std::array<Citation, 12> ga2o3_sources{{
    {"C. Fiedler et al., ECS J. Solid State Sci. Technol. 8, Q3083 (2019)", "eps_r perpendicular to (010), 10.87"},
    {"T. Onuma et al., Jpn. J. Appl. Phys. 54, 112601 (2015); M. Higashiwaki et al., AAPPS Bull. 32, 3 (2022)",
     "300 K indirect gap 4.43 eV, held at every temperature (no bulk temperature law in hand)"},
    {"Y. Mohamed et al., Appl. Phys. Lett. 101, 132106 (2012)", "electron affinity 4.00 eV"},
    {"M. Higashiwaki et al., AAPPS Bull. 32, 3 (2022)", "electron DOS mass 0.28, one valley"},
    {"Li et al., Nanomaterials 16, 961 (2026), Table 1", "hole DOS mass 3.4, a minority-carrier placeholder"},
    {"Jang, Atmaca, Cha, Micromachines 13, 1185 (2022)",
     "electron Caughey-Thomas mobility (300 K, no temperature exponents); SRH lifetimes 1.2e-8 s"},
    {"RSC Advances (2026), doi:10.1039/D6RA02684B, Table 1", "hole mobility 10, constant, a placeholder"},
    {"K. Ghosh, U. Singisetti, J. Appl. Phys. 122, 035702 (2017), Table I, y = [010]",
     "electron Barnes law: v_sat 1.5e7 cm/s, F_c 1.54e5 V/cm, gamma 2.47 (fitted to about 0.5 MV/cm)"},
    {"K. Ghosh, U. Singisetti, J. Appl. Phys. 124, 085707 (2018), Table I, y = [010]",
     "electron impact ionization A 2.16e6 /cm, B 1.77e7 V/cm (fitted 1-8 MV/cm); holes none"},
    {"A. T. Neal et al., Appl. Phys. Lett. 113, 062101 (2018)", "Si donor 30 meV; Mg acceptor 1.10 eV above E_v"},
    {"Crystals 15, 909 (2025)", "[010] thermal conductivity 0.228 W/(cm K); c_p 0.49 J/(g K) times rho 5.95 g/cm^3"},
    {"NiTCAD: no value given (owner rulings 2026-11-27)",
     "off: gap temperature law, Auger, radiative, band-gap narrowing, hole saturation law, hole impact ionization, "
     "impact-ionization temperature factor, band-to-band tunnelling, thermal-conductivity exponent"},
}};
constexpr auto ga2o3_cites = semiconductor_cites(
    {11, "default"}, std::array<ParameterCitation, 14>{{{"eps_r", 0, "measured"},
                                                        {"Eg0_eV", 1, "approximate"},
                                                        {"electron_affinity_eV", 2, "measured"},
                                                        {"Nc300", 3, "derived"},
                                                        {"Nv300", 4, "approximate"},
                                                        {"electron_mobility.", 5, "fit"},
                                                        {"lifetime.", 5, "secondary"},
                                                        {"hole_mobility.", 6, "approximate"},
                                                        {"electron_saturation.", 7, "fit"},
                                                        {"impact_ionization.electron.", 8, "fit"},
                                                        {"impact_ionization.electron.switch_V_per_cm", 11, "default"},
                                                        {"ionization.", 9, "measured"},
                                                        {"thermal.conductivity_W_cmK", 10, "measured"},
                                                        {"thermal.heat_capacity_J_cm3K", 10, "derived"}}});

constexpr std::array<MaterialEntry, 40> registry{{
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
    {"al2o3", MaterialKind::insulator, 0, 0, "", "Al2O3, amorphous (ALD) gate dielectric", al2o3_sources,
     fixed_insulator<al2o3>, al2o3_cites, true},
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
    {"sige", MaterialKind::semiconductor, 1, 1, "sige takes the Ge fraction x, 0 <= x <= 0.3",
     "Si(1-x)Ge(x), relaxed, 0 <= x <= 0.3", sige_sources, sige, sige_cites, true},
    {"sige_strained", MaterialKind::semiconductor, 1, 1, "sige_strained takes the Ge fraction x, 0 <= x <= 0.3",
     "Si(1-x)Ge(x) pseudomorphic on (001) silicon, 0 <= x <= 0.3", sige_sources, sige_strained, sige_cites, true},
    {"si_strained", MaterialKind::semiconductor, 1, 1, "si_strained takes the substrate's Ge fraction, 0 <= x <= 0.3",
     "silicon pseudomorphic on relaxed Si(1-x)Ge(x)", sige_sources, si_strained, si_strained_cites, true},
    {"gan", MaterialKind::semiconductor, 0, 0, "", "GaN, wurtzite, isotropic", nitride_sources,
     nitride_binary<gan_parameters>, gan_cites, true},
    {"aln", MaterialKind::semiconductor, 0, 0, "", "AlN, wurtzite, isotropic", nitride_sources,
     nitride_binary<aln_parameters>, aln_semiconductor_cites, true},
    {"inn", MaterialKind::semiconductor, 0, 0, "", "InN, wurtzite, isotropic (refused above about 494 K)",
     nitride_sources, nitride_binary<inn_parameters>, inn_cites, true},
    {"algan", MaterialKind::semiconductor, 1, 1, "algan takes the Al fraction x, 0 <= x <= 1",
     "Al(x)Ga(1-x)N between GaN (x = 0) and AlN (x = 1)", nitride_sources, nitride_alloy<algan_parameters>,
     nitride_alloy_cites, true},
    {"ingan", MaterialKind::semiconductor, 1, 1, "ingan takes the In fraction x, 0 <= x <= 1",
     "In(x)Ga(1-x)N between GaN (x = 0) and InN (x = 1)", nitride_sources, nitride_alloy<ingan_parameters>,
     nitride_alloy_cites, true},
    {"gaas_2001", MaterialKind::semiconductor, 0, 0, "",
     "GaAs on Vurgaftman 2001 bands, Sotoodeh mobility (the alloys' GaAs end)", iii_v_sources,
     iii_v_binary<gaas_2001_parameters>, gaas_2001_cites, true},
    {"inp", MaterialKind::semiconductor, 0, 0, "", "InP", iii_v_sources, iii_v_binary<inp_parameters>,
     iii_v_tunnelling_cites, true},
    {"inas", MaterialKind::semiconductor, 0, 0, "", "InAs", iii_v_sources, iii_v_binary<inas_parameters>,
     iii_v_tunnelling_cites, true},
    {"gap", MaterialKind::semiconductor, 0, 0, "", "GaP, indirect (X-valley conduction minimum)", iii_v_sources,
     iii_v_binary<gap_parameters>, iii_v_indirect_cites, true},
    {"alas", MaterialKind::semiconductor, 0, 0, "", "AlAs, indirect (X-valley conduction minimum)", iii_v_sources,
     iii_v_binary<alas_parameters>, iii_v_indirect_cites, true},
    {"gasb", MaterialKind::semiconductor, 0, 0, "", "GaSb", iii_v_sources, iii_v_binary<gasb_parameters>, gasb_cites,
     true},
    {"insb", MaterialKind::semiconductor, 0, 0, "", "InSb (no SRH lifetimes: srh = off, or give them)", iii_v_sources,
     iii_v_binary<insb_parameters>, insb_cites, true},
    // Unit 68 slice 3 (MAD10): `ingaas` with no fraction stays the legacy entry above; with one, this.
    {"ingaas", MaterialKind::semiconductor, 1, 1,
     "ingaas takes the In fraction x, 0 <= x <= 1, or none (the legacy In0.53 set)",
     "In(x)Ga(1-x)As between gaas_2001 (x = 0) and InAs (x = 1) (no SRH lifetimes)", iii_v_sources,
     nitride_alloy<ingaas_parameters>, iii_v_alloy_cites, true},
    {"inalas", MaterialKind::semiconductor, 1, 1, "inalas takes the In fraction x, 0 <= x <= 1",
     "In(x)Al(1-x)As between AlAs (x = 0) and InAs (x = 1) (no SRH lifetimes)", iii_v_sources,
     nitride_alloy<inalas_parameters>, iii_v_alloy_cites, true},
    {"ga2o3_beta", MaterialKind::semiconductor, 0, 0, "", "beta-Ga2O3 along [010], n-type only (300 K gap; Barnes electron law)",
     ga2o3_sources, nitride_binary<ga2o3_beta_parameters>, ga2o3_cites, true, true},
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
    // Paths added after slice 0 enter only when nonzero, so a set's hash does not move when the
    // parameter list grows (MAC1): the alloy gap (slice 1), the alloy mobility and the
    // Caughey-Thomas exponents (slice 2), the saturation velocity's temperature law (slice 3), the Barnes law
    // (slice 4b).
    const auto added_later = [](std::string_view p) {
        return p.starts_with("alloy_gap.") || p.starts_with("alloy_mobility.") || p.ends_with(".mu_min_T_exponent") ||
               p.ends_with(".N_ref_T_exponent") || p.ends_with(".alpha_T_exponent") || p.ends_with(".temperature_A") ||
               p.ends_with(".F_c_V_per_cm") || p.ends_with(".gamma");
    };
    for (const auto& [path, v] : parameter_values(m)) {
        if (v == 0.0 && added_later(path)) continue;
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

const MaterialEntry* find_material(std::string_view name, std::size_t arguments) noexcept {
    const MaterialEntry* last = nullptr;
    for (const MaterialEntry& e : registry) {
        if (e.name != name) continue;
        if (arguments >= e.arguments_min && arguments <= e.arguments_max) return &e;
        last = &e;
    }
    return last;
}

}  // namespace NiTCAD::physics
